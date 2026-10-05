"""Step 2: the turbulent two-pitch split step against flat grids, 10 km MMF link.

THE QUESTION. In vacuum (vacuum_check.py) the Schmidt two-pitch chain at
2048 px holds the received field as well as a flat 8192 px grid. Does it
hold the TURBULENT statistics too? The open risk is the scatter margin c of
the blurred extents (Schmidt (2010), DOI 10.1117/3.866274, Ch. 9,
Eqs. (9.84) and (9.85), printed p. 173): light scattered past the extent
wraps or meets the absorber, and the scatter-angle spectrum has a power-law
wing (validation/receiver_cone_clip/).

THE CHAIN. One code for every route, so only the pitch differs: Ch. 8,
Eq. (8.18), printed p. 139, with the linear pitch rule Eq. (8.8), printed
p. 136, and the screen operator T of Ch. 9, Eq. (9.3), printed p. 150, on
the planes of the production screen plan (the equal-Rytov terrestrial
planner, `turbulent_grid`). A flat grid is the chain with dx1 = dxn (every
m_i = 1, so the quadratic phases vanish). Each screen is a production
`ScreenFactory` screen on the pitch of its own plane, at the plan r0 and
L0 = 25 m. The absorber is the book super-Gaussian on every plane
(Listing 8.1, printed p. 142). The production runner at its own 2048 px grid
runs too, as a cross-check of this code.

THE RECORD. One npz for each route: the receive crop (a square of half-side
1.5 x the aperture radius, complex64, centre pixel at the middle) of each
trial and of one vacuum run, the pitch, the launched power and the wall
time. Every figure of merit (bucket power, point intensity, MMF coupling,
tip-tilt) is read from the crops after the run (`--summary`).

Run it from the repository root (the GPU run needs the cupy venv):

    python -m validation.terrestrial_sampling.turbulent_chain --trials 20 --backend cupy
    python -m validation.terrestrial_sampling.turbulent_chain --summary
"""

import argparse
import dataclasses
import math
import os
import time

import numpy as np

from olb.waveoptics.field import Begin
from olb.waveoptics.grid import GridSpec
from olb.waveoptics.mmf import mmf_coupling_efficiency
from olb.waveoptics.priority import boost_process_priority
from olb.waveoptics.propagators import set_fft_backend, xp
from olb.waveoptics.schmidt.fresnel import super_gaussian_absorber
from olb.waveoptics.schmidt.turbulence import max_partial_step
from olb.waveoptics.turbulence.run import propagate_turbulent_scenario
from olb.waveoptics.turbulence.sampling import PRESETS, turbulent_grid
from olb.waveoptics.turbulence.screens import ScreenFactory
from validation.terrestrial_mmf_tiptilt.run import (APERTURE_M, CORE_RADIUS_M,
                                                    LAM, NA, PATH_M, WAIST_M,
                                                    build_scenario)

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
K = 2 * math.pi / LAM
L0 = 25.0
SEED = 20261005
CROP_R = 1.5 * APERTURE_M / 2
F_MMF = math.pi * (APERTURE_M / 2) * CORE_RADIUS_M / (1.12 * LAM)

# name: (N, dx1, dxn); None pitches take the production side at that N.
ROUTES = {
    "flat2048": (2048, None, None),
    "flat4096": (4096, None, None),
    "flat8192": (8192, None, None),
    "tp2048_c4": (2048, 2.0e-3, 3.0e-3),
    "tp2048_c2": (2048, 2.0e-3, 4.0e-3),
    "tp4096_c4": (4096, 1.0e-3, 3.0e-3),
}


def production(cn2):
    scn, geo = build_scenario(cn2)
    grid, plan = turbulent_grid(scn, geo, preset="standard", L0_m=L0)[:2]
    return scn, geo, grid, plan


def planes(plan, n, dx1, dxn):
    """The chain planes: 0, the screen planes, L, split to the step cap."""
    cap = 0.999 * max_partial_step(dx1, dxn, n, LAM)   # Ch. 9, Eq. (9.89)
    base = np.concatenate(([0.0], plan.z_m, [PATH_M]))
    z, screen = [0.0], [None]
    for i in range(base.size - 1):
        k = max(1, math.ceil((base[i + 1] - base[i]) / cap))
        for s in range(1, k + 1):
            z.append(base[i] + (base[i + 1] - base[i]) * s / k)
            last = s == k and i < plan.z_m.size
            screen.append(i if last else None)
    return np.array(z), screen


def run_route(name, cn2, n_trials, start):
    _, _, grid, plan = production(cn2)
    n, dx1, dxn = ROUTES[name]
    if dx1 is None:
        dx1 = dxn = grid.size_m / n
    a = xp()
    z, screen = planes(plan, n, dx1, dxn)
    delta = dx1 + (dxn - dx1) * z / PATH_M                # Eq. (8.8)
    dz = np.diff(z)
    m = [float(x) for x in delta[1:] / delta[:-1]]
    dz = [float(x) for x in dz]
    idx = a.arange(n, dtype=a.float32) - n // 2
    idx2 = idx[None, :] ** 2 + idx[:, None] ** 2           # centred
    idx2_f = a.fft.ifftshift(idx2)                         # FFT layout
    absorber = a.asarray(super_gaussian_absorber(n).astype(np.float32))
    facs = {}
    for j, s in enumerate(screen):
        if s is not None:
            p = round(float(delta[j]), 12)
            if p not in facs:
                facs[p] = ScreenFactory(n, p, L0_m=L0, dtype=np.float32)

    r2 = idx2 * dx1 ** 2
    U0 = (math.sqrt(2 / math.pi) / WAIST_M) * a.exp(-r2 / WAIST_M ** 2)
    U0 = (U0 * (r2 <= (APERTURE_M / 2) ** 2)).astype(a.complex64)
    p_launch = float(a.sum(a.abs(U0) ** 2)) * dx1 ** 2
    h = math.ceil(CROP_R / dxn) + 1
    c = n // 2

    def trial(k):
        U = U0 * a.exp(1j * (K / 2 * (1 - m[0]) / dz[0] * r2)).astype(
            a.complex64)                                   # entry Q
        for i in range(len(dz)):
            df = 1.0 / (n * delta[i])
            q2 = a.exp(-1j * (2 * math.pi ** 2 * dz[i] / m[i] / K * df ** 2)
                       * idx2_f).astype(a.complex64)       # Eq. (6.12)
            U = a.fft.ifft2(q2 * a.fft.fft2(a.fft.ifftshift(U / m[i])))
            U = a.fft.fftshift(U)
            s = screen[i + 1]
            if s is not None and k is not None:
                rng = np.random.default_rng([SEED, k, s])
                phi = facs[round(float(delta[i + 1]), 12)].make(
                    float(plan.r0_m[s]), rng)
                U = U * a.exp(1j * a.asarray(phi, dtype=a.float32))
            U = U * absorber
        U = U * a.exp(1j * (K / 2 * (m[-1] - 1) / (m[-1] * dz[-1])
                            * idx2 * dxn ** 2)).astype(a.complex64)   # exit Q
        crop = U[c - h:c + h + 1, c - h:c + h + 1]
        return crop.get() if hasattr(crop, "get") else np.asarray(crop)

    vac = trial(None)
    crops, t0 = [], time.perf_counter()
    for k in range(start, start + n_trials):
        crops.append(trial(k).astype(np.complex64))
    wall = (time.perf_counter() - t0) / n_trials
    return dict(crops=np.array(crops), vac=vac.astype(np.complex64),
                dxn=dxn, dx1=dx1, n=n, p_launch=p_launch, planes=z.size,
                s_per_trial=wall)


def run_runner(cn2, n_trials, start, backend):
    """The production runner at its own grid, the cross-check."""
    scn, geo, grid, plan = production(cn2)
    t0 = time.perf_counter()
    res = propagate_turbulent_scenario(
        scn, geo, n_trials=n_trials, seed=SEED, grid=grid, plan=plan,
        L0_m=L0, start_index=start, fft_backend=backend)
    wall = (time.perf_counter() - t0) / n_trials
    return dict(collected=np.array([t.collected_power for t in res.trials]),
                mmf=np.array([t.mmf_eta for t in res.trials]),
                s_per_trial=wall)


# ---------------------------------------------------------------- the reads

def bucket(crop, dx):
    h = crop.shape[-1] // 2
    i = np.arange(-h, h + 1) * dx
    disc = (i[None, :] ** 2 + i[:, None] ** 2) <= (APERTURE_M / 2) ** 2
    return disc


def mmf_eta(crop, dx):
    """MMF coupling of one crop: clip, zero-pad for a fine focal pixel."""
    disc = bucket(crop, dx)
    npad = 1 << max(8, math.ceil(math.log2(0.75 / dx)))   # >= 2.5 um focal px
    E = np.zeros((npad, npad), complex)
    h, c = crop.shape[0] // 2, npad // 2
    E[c - h:c + h + 1, c - h:c + h + 1] = crop * disc
    F = Begin(npad * dx, LAM, npad)
    F.field = E
    return mmf_coupling_efficiency(F, APERTURE_M, CORE_RADIUS_M, F_MMF,
                                   numerical_aperture=NA)


def summary():
    print(f"{'route':10s} {'n':>4s} {'s/trial':>8s} {'pen dB':>8s} "
          f"{'s2_P':>7s} {'s2_I pt':>8s} {'MMF mean':>9s} {'MMF p5':>7s} "
          f"{'loss p5 dB':>10s}")
    for name in ROUTES:
        path = os.path.join(DATA, f"{name}.npz")
        if not os.path.exists(path):
            continue
        d = np.load(path)
        dx, crops, vac = float(d["dxn"]), d["crops"], d["vac"]
        disc = bucket(vac, dx)
        pv = (np.abs(vac[disc]) ** 2).sum()
        p = (np.abs(crops[:, disc]) ** 2).sum(axis=1) / pv
        h = vac.shape[0] // 2
        pt = np.abs(crops[:, h, h]) ** 2 / abs(vac[h, h]) ** 2
        eta = np.array([mmf_eta(cr, dx) for cr in crops])
        loss = -10 * np.log10(p * eta)
        print(f"{name:10s} {len(p):4d} {float(d['s_per_trial']):8.2f} "
              f"{-10 * np.log10(p.mean()):8.3f} {p.var() / p.mean() ** 2:7.3f} "
              f"{pt.var() / pt.mean() ** 2:8.3f} {eta.mean():9.4f} "
              f"{np.percentile(eta, 5):7.4f} {np.percentile(loss, 95):10.3f}")
    path = os.path.join(DATA, "runner2048.npz")
    if os.path.exists(path):
        d = np.load(path)
        p = d["collected"]
        print(f"{'runner2048':10s} {len(p):4d} {float(d['s_per_trial']):8.2f} "
              f"{'(raw)':>8s} {p.var() / p.mean() ** 2:7.3f} {'':8s} "
              f"{d['mmf'].mean():9.4f} {np.percentile(d['mmf'], 5):7.4f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=20)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--cn2", type=float, default=1e-14)
    ap.add_argument("--backend", default="numpy")
    ap.add_argument("--routes", nargs="*", default=list(ROUTES) + ["runner2048"])
    ap.add_argument("--max-8192", type=int, default=5,
                    help="trial cap for the costly flat8192 reference")
    ap.add_argument("--summary", action="store_true")
    args = ap.parse_args()
    if args.summary:
        summary()
        return
    os.makedirs(DATA, exist_ok=True)
    boost_process_priority()          # an ssh launch on bigfraw is throttled
    set_fft_backend(args.backend)
    for name in args.routes:
        t = time.perf_counter()
        if name == "runner2048":
            out = run_runner(args.cn2, args.trials, args.start, args.backend)
        else:
            n = args.trials if name != "flat8192" else min(args.trials,
                                                           args.max_8192)
            out = run_route(name, args.cn2, n, args.start)
        np.savez(os.path.join(DATA, f"{name}.npz"), **out)
        print(f"{name}: {out['s_per_trial']:.2f} s/trial, "
              f"{time.perf_counter() - t:.0f} s total", flush=True)


if __name__ == "__main__":
    main()

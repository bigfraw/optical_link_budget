r"""Experiment (b): a grid that follows the beam.

RECOMMENDATION (read first): BURY the beam-sized screens (step 2) and the
per-plane-pitch chain (step 3) for the CURRENT wired scenarios; the flat grid
already wins. Step 1 CONFIRMS the space slab fills the grid (the plane-wave
input covers >99% of the interior at every plane), so beam-sizing saves nothing
there. Step 2 saves a little screen-generation time on the terrestrial launch
end but shifts sigma2_I past the 10% kill line because the small screens starve
the low-frequency tilt. Step 3 (gap S-14) reproduces the reference correctly on
the same grid but buys no pixel-operations on either wired case: the space slab
holds a constant extent, and the terrestrial graded-pitch variant that DOES cut
pixels drifts mean power and sigma2_I past the kill line. The per-plane-pitch
win is real only for a long co-moving space path, which the runner already
avoids by simulating the ~20 km slab alone with a plane-wave input.

WHAT THIS SCRIPT MEASURES. Three steps, all report-only, no production change.

STEP 1 (space). It VERIFIES the claim that the space slab gains nothing from
beam-sized screens, because the input plane wave fills the grid. It propagates
one turbulent snapshot of the 30 deg standard slab and reports the irradiance
support (the 90% encircled-energy radius and the fill fraction) at every screen
plane.

STEP 2 (terrestrial). It prototypes beam-extent-sized screens: at each plane it
generates an m x m screen at the SAME pixel pitch that covers the beam support
plus a guard, embeds it in the full grid with zero phase outside, and measures
the speed and the error against the full-screen reference (the same statistics
as experiment (a)). This saves SCREEN GENERATION only, not the FFT cost.

STEP 3 (gap S-14). It prototypes a per-plane-pitch chain built on the schmidt
angular-spectrum kernel, with the linear pitch schedule of
sampling.partial_grid_spacing and the plane geometry checked by constraints 1 to
4 (check_sampling). The screens live on per-plane pitches (screen_r0 stays
valid; the pitch enters phase_screen directly). It compares the total
pixel-operations and the wall time against the one-flat-grid Forvard reference,
at equal accuracy, for the terrestrial case AND the space slab.

KILL CRITERION (stated before the run): a variant DIES if it cannot hold the
mean power inside 0.1 dB AND sigma2_I inside 10% of the reference while saving
time (pixel-operations or wall time).

Sources:
- Schmidt, DOI 10.1117/3.866274, Ch. 8 (the partial-propagation chain,
  Eq. (8.18), printed p. 139; the pitch rule Eq. (8.8), printed p. 136; the
  step cap Eq. (8.24), printed p. 144) and Ch. 7 (the four constraints).
- Fried, DOI 10.1364/JOSA.56.001372 (r0). Andrews and Phillips,
  DOI 10.1117/3.626196, Ch. 7 (the long-term beam radius).

Run from the repository root:
    python -m validation.waveoptics_speed.beam_grid_experiment
"""

import json
import math
import os
import platform
import time
import warnings

import numpy as np

from olb.beam import free_space_radius, virtual_waist
from olb.geometry import CircularOrbit, HorizontalPath
from olb.scenario import (Channel, SpaceScenario, TerrestrialChannel,
                          TerrestrialScenario)
from olb.terminal import SMF, Terminal, Transmitter
from olb.turbulence.profiles import DEFAULT_HS, default_cn2_profile
from olb.waveoptics.field import Begin, Field, Power
from olb.waveoptics.propagators import Forvard
from olb.waveoptics.run import (_clip, _launch_aperture, _normalised_gauss)
from olb.waveoptics.schmidt.fresnel import angular_spectrum
from olb.waveoptics.schmidt.sampling import check_sampling, partial_grid_spacing
from olb.waveoptics.smf import coupling_efficiency
from olb.waveoptics.sources import GaussBeam
from olb.waveoptics.threader import Threader
from olb.waveoptics.turbulence.run import _screen_seed
from olb.waveoptics.turbulence.sampling import PRESETS
from olb.waveoptics.turbulence.screens import ScreenFactory, phase_screen, Screen
from olb.waveoptics.turbulence.splitstep import (_apply_mask, _substeps,
                                                 super_gaussian_boundary)

# Reuse the case setup, trial runner, and timing of experiment (a).
from validation.waveoptics_speed.coarse_screen_experiment import (
    LAM, SEED, full_builder, run_trials, setup_case, stats, time_stack,
    warm_aotools)

HERE = os.path.dirname(__file__)
N_TRIALS = 200


# ---------------------------------------------------------------------------
# the cases (standard preset, to give step 3 a meatier grid)
# ---------------------------------------------------------------------------

def terrestrial_case():
    scn = TerrestrialScenario(
        near=Terminal(aperture_m=0.10, wavelength_m=LAM,
                      transmitter=Transmitter(waist_m=0.02)),
        far=Terminal(aperture_m=0.20, wavelength_m=LAM, detector=SMF()),
        channel=TerrestrialChannel(path_length_m=2000.0, cn2=5e-15))
    return dict(name="terrestrial 2km standard", preset="standard",
                scenario=scn, geometry=HorizontalPath(2000.0),
                hs=None, cn2=None, is_space=False)


def space_case():
    ground = Terminal(aperture_m=0.50, wavelength_m=LAM,
                      transmitter=Transmitter(waist_m=0.05))
    scn = SpaceScenario(ground=ground,
                        space=Terminal(aperture_m=0.30, wavelength_m=LAM),
                        direction="downlink", channel=Channel(altitude_m=600e3))
    hs = DEFAULT_HS
    cn2 = default_cn2_profile(scn.channel.site, hs)
    return dict(name="space downlink 30deg standard", preset="standard",
                scenario=scn, geometry=CircularOrbit(altitude_m=600e3,
                                                     elevation_deg=[30.0]),
                hs=hs, cn2=cn2, is_space=True)


# ===========================================================================
# STEP 1: the space slab fills the grid
# ===========================================================================

def step1_space_support(setup, seed=SEED, trial=0):
    """Report the irradiance support at each screen plane of a space slab."""
    grid, plan, mask = setup["grid"], setup["plan"], setup["mask"]
    preset = setup["preset"]
    n, dx = grid.n, grid.pixel_m
    fac = ScreenFactory(n, dx)
    stack = [fac.make(plan.r0_m[j],
                      np.random.default_rng(_screen_seed(seed, trial, j)))
             for j in range(plan.r0_m.size)]

    xx = (np.arange(n) - n // 2) * dx
    R = np.hypot(xx[:, None], xx[None, :])
    r_flat = (1.0 - preset.boundary_width_frac) * grid.size_m / 2
    interior = R <= r_flat
    Rin = R[interior]
    order = np.argsort(Rin)
    Rsort = Rin[order]

    def support(F):
        I = np.abs(F.field) ** 2
        Iin = I[interior]
        c = np.cumsum(Iin[order])
        c = c / c[-1]
        r90 = float(Rsort[np.searchsorted(c, 0.9)])
        fill = float((Iin > 0.5 * Iin.mean()).mean())
        return r90 / r_flat, fill

    F = Begin(grid.size_m, LAM, n)
    max_step = n * dx ** 2 / LAM
    rows = []
    r90, fill = support(F)
    rows.append({"plane": "input", "z_km": 0.0, "r90_over_interior": r90,
                 "fill_fraction": fill})
    here = 0.0
    for j, (zi, scr) in enumerate(zip(plan.z_m, stack)):
        for dz in _substeps(zi - here, max_step):
            F = Forvard(F, dz)
            F = _apply_mask(F, mask)
        r90, fill = support(F)
        rows.append({"plane": f"screen {j}", "z_km": float(zi / 1e3),
                     "r90_over_interior": r90, "fill_fraction": fill})
        F = Screen(F, scr)
        F = _apply_mask(F, mask)
        here = zi
    for dz in _substeps(plan.z_total_m - here, max_step):
        F = Forvard(F, dz)
        F = _apply_mask(F, mask)
    r90, fill = support(F)
    rows.append({"plane": "rx", "z_km": float(plan.z_total_m / 1e3),
                 "r90_over_interior": r90, "fill_fraction": fill})
    return rows


# ===========================================================================
# STEP 2: beam-extent-sized screens (terrestrial)
# ===========================================================================

def _beam_radius(scn, z, r0_total, z_total):
    """Estimate the beam support radius at plane z, with a scatter cone.

    Free-space radius from olb.beam (it reads the Transmitter divergence), plus
    the turbulence scatter cone (lambda/r0_total) z. A guard is added by the
    caller. See Andrews and Phillips, DOI 10.1117/3.626196, Ch. 7.
    """
    t = scn.tx_terminal.transmitter
    w_free = float(free_space_radius(t.waist_m, z, t.divergence_rad, LAM))
    scatter = (LAM / r0_total) * z
    tx_ap = _launch_aperture(scn.tx_terminal)[0]
    return max(w_free, tx_ap / 2.0) + scatter


def beam_sized_builder(setup, guard=2.0):
    """Build a screen stack where each screen covers the beam plus a guard.

    Each screen has m x m pixels at the SAME pitch dx, embedded centred in the
    full n x n grid with zero phase outside. The m x m physical side is smaller,
    so the subharmonics see a higher fundamental frequency: the beam-sized
    screen carries LESS large-scale tilt. That is the error this step measures.
    """
    grid, plan = setup["grid"], setup["plan"]
    scn = setup["scenario"]
    n, dx = grid.n, grid.pixel_m
    r0 = plan.r0_m
    m_pix = []
    for j in range(r0.size):
        r_beam = guard * _beam_radius(scn, plan.z_m[j], plan.r0_total_m,
                                      plan.z_total_m)
        m = int(2 * math.ceil((2 * r_beam / dx) / 2))       # even
        m_pix.append(int(min(max(m, 16), n)))
    facs = {m: ScreenFactory(m, dx) for m in set(m_pix)}

    def build(entropy, k):
        out = []
        for j in range(r0.size):
            m = m_pix[j]
            rng = np.random.default_rng(_screen_seed(entropy, k, j))
            small = facs[m].make(r0[j], rng)
            if m == n:
                out.append(small)
            else:
                big = np.zeros((n, n))
                off = (n - m) // 2
                big[off:off + m, off:off + m] = small
                out.append(big)
        return out
    return build, m_pix


# ===========================================================================
# STEP 3: the per-plane-pitch chain (gap S-14)
# ===========================================================================

def _flat_hops(plan, n, dx):
    """Count the Forvard hops of the flat reference split step."""
    max_step = n * dx ** 2 / LAM
    z = np.concatenate(([0.0], plan.z_m, [plan.z_total_m]))
    gaps = np.diff(z)
    return int(sum(max(1, int(math.ceil(g / max_step))) if g > 0 else 0
                   for g in gaps))


def _chain_planes(plan, n, dx1, dxn):
    """Build the plane list, pitches, and screen flags of the chain.

    The base planes are the input, the screen planes, and the receiver. Each
    segment is subdivided so every step obeys the Schmidt step cap
    min(pitch)^2 N / lambda (Ch. 8, Eq. (8.24), printed p. 144). The pitch of
    every plane is the linear schedule of partial_grid_spacing.
    """
    z_total = plan.z_total_m
    base = np.concatenate(([0.0], plan.z_m, [z_total]))
    is_screen_base = [False] + [True] * plan.z_m.size + [False]

    def pitch(z):
        return float(partial_grid_spacing(dx1, dxn, z / z_total))

    planes, screens_r0idx = [0.0], [None]
    for i in range(base.size - 1):
        za, zb = base[i], base[i + 1]
        pa, pb = pitch(za), pitch(zb)
        max_step = min(pa, pb) ** 2 * n / LAM
        nsub = max(1, int(math.ceil((zb - za) / max_step))) if zb > za else 1
        for s in range(1, nsub + 1):
            zc = za + (zb - za) * s / nsub
            planes.append(float(zc))
            # only the true screen base planes carry a screen
            screens_r0idx.append((i) if (s == nsub and is_screen_base[i + 1]
                                         and i < plan.z_m.size) else None)
    # map screen planes to r0 index (0..n_screens-1)
    idx = []
    counter = 0
    for flag in screens_r0idx:
        if flag is not None:
            idx.append(counter)
            counter += 1
        else:
            idx.append(None)
    pitches = [pitch(z) for z in planes]
    return np.array(planes), pitches, idx


def chain_builder(setup, n_proto, dx1, dxn):
    """Return (build_stack, planes, pitches, screen_idx, absorber).

    build_stack(entropy, k) -> list of screens, one per plane (None where a
    plane carries no screen), each on that plane's pitch.
    """
    plan = setup["plan"]
    planes, pitches, screen_idx = _chain_planes(plan, n_proto, dx1, dxn)
    r0 = plan.r0_m
    # one factory per unique screen-plane pitch.
    facs = {}
    for p, si in zip(pitches, screen_idx):
        if si is not None and p not in facs:
            facs[p] = ScreenFactory(n_proto, p)

    def build(entropy, k):
        stack = []
        for p, si in zip(pitches, screen_idx):
            if si is None:
                stack.append(None)
            else:
                rng = np.random.default_rng(_screen_seed(entropy, k, si))
                stack.append(facs[p].make(r0[si], rng))
        return stack
    return build, planes, pitches, screen_idx


def _chain_start(setup, n_proto, dx1):
    """The chain start field array (space plane wave, or terrestrial launch)."""
    if setup["is_space"]:
        return np.ones((n_proto, n_proto), dtype=complex)
    scn = setup["scenario"]
    tx = scn.tx_terminal
    t = tx.transmitter
    w_v, offset = virtual_waist(t.waist_m, t.divergence_rad, LAM)
    F0 = _normalised_gauss(GaussBeam(Begin(n_proto * dx1, LAM, n_proto), w_v))
    if offset > 0:
        from olb.waveoptics.propagators import GForvard
        F0 = GForvard(F0, offset)
    F0 = _clip(F0, *_launch_aperture(tx))
    return F0.field.astype(complex)


def as_chain(U0, planes, pitches, stack, absorber):
    """Propagate U0 along the per-plane-pitch angular-spectrum chain."""
    U = U0
    for i in range(len(planes) - 1):
        dz = planes[i + 1] - planes[i]
        if dz > 0:
            U = angular_spectrum(U, LAM, pitches[i], dz, dx2=pitches[i + 1])
        U = U * absorber
        if stack[i + 1] is not None:
            U = U * np.exp(1j * stack[i + 1])
    return U


def run_chain_trials(setup, n_proto, dx1, dxn, n_trials, threader, seed=SEED):
    """Run the chain and return powers, etas, and geometry diagnostics."""
    scn, plan = setup["scenario"], setup["plan"]
    preset = setup["preset"]
    rx, is_space = setup["rx"], setup["is_space"]
    build, planes, pitches, screen_idx = chain_builder(setup, n_proto, dx1, dxn)
    absorber = super_gaussian_boundary(n_proto, preset.boundary_width_frac)
    U0 = _chain_start(setup, n_proto, dx1)
    side_n = n_proto * dxn

    # the space vacuum baseline on the SAME chain (flat screens).
    if is_space:
        flat_stack = [None] * len(planes)
        U_vac = as_chain(U0, planes, pitches, flat_stack, absorber)
        F_vac = Begin(side_n, LAM, n_proto)
        F_vac.field = U_vac
        p_reference = Power(_clip(F_vac, rx.aperture_m, rx.obscuration_ratio))
    else:
        p_reference = float((np.abs(U0) ** 2).sum() * dx1 ** 2)

    def run_one(k):
        stack = build(seed, k)
        U = as_chain(U0, planes, pitches, stack, absorber)
        F = Begin(side_n, LAM, n_proto)
        F.field = U
        collected = _clip(F, rx.aperture_m, rx.obscuration_ratio)
        power = float(Power(collected) / p_reference)
        eta = (float(coupling_efficiency(collected, rx.aperture_m))
               if isinstance(rx.detector, SMF) else None)
        return power, eta

    results = threader.map(run_one, range(n_trials))
    powers = np.array([r[0] for r in results], dtype=float)
    etas = np.array([r[1] for r in results if r[1] is not None], dtype=float)

    n_steps = len(planes) - 1
    diag = {"n_proto": n_proto, "dx1_mm": dx1 * 1e3, "dxn_mm": dxn * 1e3,
            "n_planes": len(planes), "n_steps": n_steps,
            "chain_ffts": 2 * n_steps,
            "side_launch_m": n_proto * dx1, "side_rx_m": side_n}
    # constraint check on the worst step (min pitch, longest sub-step).
    dz = np.diff(planes)
    worst = int(np.argmax(dz)) if dz.size else 0
    D1 = n_proto * min(pitches)
    rules = check_sampling(D1, D1, pitches[worst], pitches[worst + 1], n_proto,
                           LAM, float(dz[worst]))
    diag["worst_step_constraints"] = {r.name.split(":")[0]: bool(r.satisfied)
                                      for r in rules}
    return powers, etas, diag


def time_chain_prop(setup, n_proto, dx1, dxn, reps=5, seed=SEED):
    """Median wall time of ONE chain propagation (screens prebuilt)."""
    build, planes, pitches, screen_idx = chain_builder(setup, n_proto, dx1, dxn)
    absorber = super_gaussian_boundary(n_proto, setup["preset"].boundary_width_frac)
    U0 = _chain_start(setup, n_proto, dx1)
    stack = build(seed, 0)
    as_chain(U0, planes, pitches, stack, absorber)      # warm
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter()
        as_chain(U0, planes, pitches, stack, absorber)
        ts.append(time.perf_counter() - t0)
    return float(np.median(ts))


def time_flat_prop(setup, reps=5, seed=SEED):
    """Median wall time of ONE flat Forvard split-step propagation."""
    grid, plan, mask = setup["grid"], setup["plan"], setup["mask"]
    n, dx = grid.n, grid.pixel_m
    fac = ScreenFactory(n, dx)
    stack = [fac.make(plan.r0_m[j],
                      np.random.default_rng(_screen_seed(seed, 0, j)))
             for j in range(plan.r0_m.size)]
    from olb.waveoptics.turbulence.splitstep import split_step
    F_start = (Begin(grid.size_m, LAM, n) if setup["is_space"]
               else setup["start"])
    split_step(F_start, plan.z_m, stack, plan.z_total_m, boundary=mask)  # warm
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter()
        split_step(F_start, plan.z_m, stack, plan.z_total_m, boundary=mask)
        ts.append(time.perf_counter() - t0)
    return float(np.median(ts))


# ===========================================================================
# main
# ===========================================================================

def _cmp(ref, mean, sigma2):
    """Return (mean_db, sigma2_pct)."""
    return (10.0 * np.log10(mean / ref["mean_power"]),
            100.0 * (sigma2 - ref["sigma2_I"]) / ref["sigma2_I"])


def run_step2(setup, ref, threader):
    """Beam-sized screens on the terrestrial case."""
    build, m_pix = beam_sized_builder(setup)
    powers, etas = run_trials(setup, build, N_TRIALS, threader)
    s = stats(powers, etas)
    s["build_s"] = time_stack(build)
    mean_db, sig_pct = _cmp(ref, s["mean_power"], s["sigma2_I"])
    s["mean_db"], s["sigma2_pct"] = float(mean_db), float(sig_pct)
    s["m_pix"] = m_pix
    s["n"] = setup["grid"].n
    s["speedup_build"] = float(ref["build_s"] / s["build_s"])
    s["saves_time"] = bool(s["build_s"] < ref["build_s"])
    s["killed"] = bool(not (abs(mean_db) <= 0.1 and abs(sig_pct) <= 10.0
                            and s["saves_time"]))
    print(f'  m_pix per screen: {m_pix} (grid n={setup["grid"].n})')
    print(f'  beam-sized: mean {s["mean_power"]:.5f} ({mean_db:+.2f} dB)  '
          f'sigma2_I {s["sigma2_I"]:.5f} ({sig_pct:+.1f}%)  '
          f'eta {s["smf_eta"]:.5f}  build {s["build_s"]*1e3:.1f} ms '
          f'({s["speedup_build"]:.2f}x)  killed={s["killed"]}')
    return s


def run_step3(setup, ref, threader, variants):
    """The per-plane-pitch chain, one or more (n_proto, dx1, dxn) variants."""
    grid, plan = setup["grid"], setup["plan"]
    dx_ref = grid.pixel_m
    flat_hops = _flat_hops(plan, grid.n, dx_ref)
    flat_ffts = 2 * flat_hops
    flat_ops = flat_ffts * grid.n ** 2 * math.log2(grid.n)
    flat_time = time_flat_prop(setup)
    print(f'  flat reference: n={grid.n}, hops={flat_hops}, '
          f'ffts={flat_ffts}, ops={flat_ops:.3e}, prop {flat_time*1e3:.1f} ms')
    out = {"flat": {"n": grid.n, "hops": flat_hops, "ffts": flat_ffts,
                    "ops": float(flat_ops), "prop_s": flat_time}}
    rows = []
    for name, (n_proto, f1, fn) in variants.items():
        dx1, dxn = dx_ref * f1, dx_ref * fn
        powers, etas, diag = run_chain_trials(setup, n_proto, dx1, dxn,
                                              N_TRIALS, threader)
        s = stats(powers, etas)
        mean_db, sig_pct = _cmp(ref, s["mean_power"], s["sigma2_I"])
        prop_s = time_chain_prop(setup, n_proto, dx1, dxn)
        ops = diag["chain_ffts"] * n_proto ** 2 * math.log2(n_proto)
        ops_ratio = ops / flat_ops
        time_ratio = prop_s / flat_time
        saves = bool(ops < flat_ops or prop_s < flat_time)
        killed = bool(not (abs(mean_db) <= 0.1 and abs(sig_pct) <= 10.0
                           and saves))
        row = dict(name=name, mean_power=s["mean_power"], smf_eta=s["smf_eta"],
                   sigma2_I=s["sigma2_I"], mean_db=float(mean_db),
                   sigma2_pct=float(sig_pct), ops=float(ops),
                   ops_ratio=float(ops_ratio), prop_s=prop_s,
                   time_ratio=float(time_ratio), saves_time=saves,
                   killed=killed, **diag)
        rows.append(row)
        print(f'  {name:16s} n={n_proto} steps={diag["n_steps"]:2d} '
              f'ffts={diag["chain_ffts"]:2d}  mean {s["mean_power"]:.5f} '
              f'({mean_db:+.2f} dB)  sigma2 {s["sigma2_I"]:.5f} '
              f'({sig_pct:+.1f}%)  ops {ops_ratio:.2f}x  '
              f'time {time_ratio:.2f}x  killed={killed}')
        print(f'                   constraints {diag["worst_step_constraints"]}')
    out["variants"] = rows
    return out


def main():
    t_start = time.time()
    warm_aotools()
    threader = Threader()

    terr = terrestrial_case()
    space = space_case()
    terr_setup = setup_case(terr)
    space_setup = setup_case(space)

    # references (flat, full-resolution) for the accuracy comparison.
    print("building flat references (200 trials each)")
    terr_ref = _reference(terr_setup, threader)
    space_ref = _reference(space_setup, threader)
    print(f'  terrestrial ref: mean {terr_ref["mean_power"]:.5f}  '
          f'sigma2_I {terr_ref["sigma2_I"]:.5f}  eta {terr_ref["smf_eta"]:.5f}')
    print(f'  space ref:       mean {space_ref["mean_power"]:.5f}  '
          f'sigma2_I {space_ref["sigma2_I"]:.5f} +/- '
          f'{space_ref["se_sigma2_I"]:.5f}')

    print("\nSTEP 1: the space slab fills the grid (one snapshot, 30 deg std)")
    step1 = step1_space_support(space_setup)
    print("  plane        z[km]   r90/interior   fill_fraction")
    for r in step1:
        print(f'  {r["plane"]:12s} {r["z_km"]:6.1f}   '
              f'{r["r90_over_interior"]:9.3f}      {r["fill_fraction"]:.4f}')

    print("\nSTEP 2: beam-sized screens (terrestrial 2km standard)")
    step2 = run_step2(terr_setup, terr_ref, threader)

    print("\nSTEP 3: per-plane-pitch chain (gap S-14)")
    print(" terrestrial:")
    terr_variants = {
        "match n256 1:1": (terr_setup["grid"].n, 1.0, 1.0),
        "graded n128 1:4": (128, 1.0, 4.0),
        "graded n192 1:2": (192, 1.0, 2.0),
    }
    step3_terr = run_step3(terr_setup, terr_ref, threader, terr_variants)
    print(" space slab:")
    space_variants = {
        "match n512 1:1": (space_setup["grid"].n, 1.0, 1.0),
        "half n256 1:1": (256, 1.0, 1.0),
    }
    step3_space = run_step3(space_setup, space_ref, threader, space_variants)

    out = {
        "environment": {"numpy": np.__version__,
                        "platform": platform.platform(),
                        "cores": os.cpu_count()},
        "n_trials": N_TRIALS,
        "references": {"terrestrial": terr_ref, "space": space_ref},
        "step1_space_support": step1,
        "step2_beam_sized_screens": step2,
        "step3_terrestrial": step3_terr,
        "step3_space": step3_space,
    }
    path = os.path.join(HERE, "beam_grid_experiment_results.json")
    with open(path, "w") as fp:
        json.dump(out, fp, indent=1)

    print("\n=== SUMMARY ===")
    print(f'  step1: space fill fraction stays '
          f'>= {min(r["fill_fraction"] for r in step1):.3f} at every plane '
          f'-> the plane wave fills the grid; beam-sizing saves nothing.')
    print(f'  step2 (terrestrial beam-sized): killed={step2["killed"]} '
          f'(sigma2 {step2["sigma2_pct"]:+.1f}%, mean {step2["mean_db"]:+.2f} '
          f'dB, build {step2["speedup_build"]:.2f}x)')
    print("  step3 terrestrial:")
    for v in step3_terr["variants"]:
        print(f'    {v["name"]:16s} killed={v["killed"]}  ops {v["ops_ratio"]:.2f}x '
              f'time {v["time_ratio"]:.2f}x  sigma2 {v["sigma2_pct"]:+.1f}%  '
              f'mean {v["mean_db"]:+.2f} dB')
    print("  step3 space:")
    for v in step3_space["variants"]:
        print(f'    {v["name"]:16s} killed={v["killed"]}  ops {v["ops_ratio"]:.2f}x '
              f'time {v["time_ratio"]:.2f}x  sigma2 {v["sigma2_pct"]:+.1f}%  '
              f'mean {v["mean_db"]:+.2f} dB')
    print(f'\nwrote {path}')
    print(f'(elapsed {time.time() - t_start:.0f} s)')


def _reference(setup, threader):
    """The flat full-resolution reference stats."""
    build = full_builder(setup["grid"], setup["plan"])
    powers, etas = run_trials(setup, build, N_TRIALS, threader)
    s = stats(powers, etas)
    s["build_s"] = time_stack(build)
    return s


if __name__ == "__main__":
    main()

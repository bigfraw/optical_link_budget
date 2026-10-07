"""Validity checks of the 10 km MMF smoke campaigns (smoke.py), read post hoc.

For each Cn2 it reads the WIDE stored field (radius 1.5 m) and measures:

1. The beam: the second-moment radius and the centroid wander of each trial,
   against the free-space W, the strong-fluctuation long-term W_LT
   (Andrews and Phillips, DOI 10.1117/3.626196, Ch. 7, Eq. (57), printed
   p. 242) and the beam wander (Ch. 6, Eq. (93), printed p. 203).
2. The point scintillation index near the axis, against the strong Gaussian
   index (Ch. 8, DOI 10.1117/3.626196).
3. The sampling: the speckle (intensity correlation) width in pixels, the
   coherence radius rho_0 of a Gaussian beam (Ch. 6, Table, printed p. 193)
   in pixels, and the wrapped phase step per pixel (the slope route warns at
   2.8 rad).
4. The tip-tilt: the gradient tilt (mean wrapped phase gradient over the
   aperture, Andrews G-tilt, Ch. 6, Eq. (82), printed p. 200) and the
   centroid tilt (intensity-weighted gradient, the quad-cell reading) of
   sub-apertures of D = 1, 2 and 4 inch tiled inside r < 0.5 m, against the
   aperture angle-of-arrival variance at the Gaussian-beam r0 and L0 = 25 m
   (Ch. 6, Eq. (83) and (84), printed p. 201).
5. The fibre: the collected power, the MMF coupling with no tip-tilt and
   with a perfect tip-tilt (run.summarise).

It writes figures/ and prints one table for each Cn2.

    python -m validation.terrestrial_mmf_tiptilt.analyse_smoke

It compares the sized 2048 px grid with a 4096 px grid of the same side
(smoke.py --n 4096). Each value carries a bootstrap SE over the trials.
"""

import math
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from olb.turbulence.andrews.beam import beam_params, effective_beam_params
from olb.turbulence.andrews.scintillation import (rytov_variance,
                                                  scintillation_index)
from olb.turbulence.andrews.structure import coherence_radius
from olb.turbulence.andrews.wander import beam_wander_variance
from olb.turbulence.angle_of_arrival import aperture_arrival_angle_variance
from olb.terminal import TipTilt
from validation.terrestrial_mmf_tiptilt.run import (APERTURE_M, CN2, HERE, LAM,
                                                    PATH_M, WAIST_M)
from validation.terrestrial_sampling.turbulent_chain import mmf_eta
from validation.terrestrial_mmf_tiptilt.smoke import smoke_campaign

FIG = os.path.join(HERE, "figures")
K = 2 * math.pi / LAM
L0 = 25.0
INCH = 0.0254
D_SUB = (INCH, 2 * INCH, 4 * INCH)
R_TILE = 0.5


def analytic(cn2):
    bp = beam_params(WAIST_M, LAM, PATH_M)
    s2r = float(rytov_variance(LAM, PATH_M, cn2))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        w_lt = float(effective_beam_params(bp, s2r).w)
        rho0 = float(coherence_radius(LAM, PATH_M, cn2, wave="gaussian",
                                      beam=bp))
        si = float(scintillation_index(LAM, PATH_M, cn2, wave="gaussian",
                                       regime="strong", beam=bp,
                                       tracked=False))
        wander = math.sqrt(beam_wander_variance(
            WAIST_M, LAM, PATH_M, cn2, spectrum="exponential", L0=L0))
        aoa = {D: math.sqrt(aperture_arrival_angle_variance(
            D, 2.1 * rho0, LAM, L0=L0)) for D in D_SUB}
    return dict(W=float(bp.w), W_LT=w_lt, s2r=s2r, rho0=rho0, si=si,
                wander=wander, wander_diff=wander_diffractive(cn2), aoa=aoa)


def wander_diffractive(cn2):
    """Per-axis rms beam wander (m) with the DIFFRACTIVE beam size W(z).

    <r_c^2> = 8 pi^2 int_0^L (L-z)^2 int k^3 Phi_n(k) exp(-k^2 W(z)^2) dk dz,
    the form of Andrews and Phillips, DOI 10.1117/3.626196, Ch. 6, Eq. (93),
    printed p. 203, with the von Karman Phi_n = 0.033 Cn2 (k^2+k0^2)^(-11/6),
    k0 = 2 pi/L0 (Ch. 3, Eq. (20), printed p. 64). `beam_wander_variance`
    puts the GEOMETRIC W(z) in the filter, which reads 2 to 3x high for this
    strongly diffracting launch (Lambda0 ~ 240).
    """
    from scipy.integrate import quad
    zr = math.pi * WAIST_M ** 2 / LAM
    k0 = 2 * math.pi / L0
    zs = np.geomspace(1e-2, PATH_M, 300)
    inner = [quad(lambda q, w=WAIST_M * math.sqrt(1 + (z / zr) ** 2):
                  q ** 3 * 0.033 * (q * q + k0 * k0) ** (-11 / 6)
                  * math.exp(-q * q * w * w), 0, np.inf, limit=400)[0]
             for z in zs]
    v = 8 * math.pi ** 2 * cn2 * np.trapezoid((PATH_M - zs) ** 2 * inner, zs)
    return math.sqrt(v / 2)


def grads(U, dx):
    """Wrapped phase gradients (rad/m), on the pixel grid (edge-padded)."""
    gx = np.angle(U[:, 1:] * np.conj(U[:, :-1])) / dx
    gy = np.angle(U[1:, :] * np.conj(U[:-1, :])) / dx
    return np.pad(gx, ((0, 0), (0, 1)), mode="edge"), \
        np.pad(gy, ((0, 1), (0, 0)), mode="edge")


def sub_tilts(U, dx, X, Y, D):
    """G-tilt and centroid tilt (rad) of each sub-aperture inside R_TILE."""
    gx, gy = grads(U, dx)
    I = np.abs(U) ** 2
    c = np.arange(-R_TILE, R_TILE + 1e-9, D)
    N0, h = U.shape[0] // 2, math.ceil(D / 2 / dx) + 1
    out = []
    for cx in c:
        for cy in c:
            if math.hypot(cx, cy) > R_TILE - D / 2:
                continue
            j, i = N0 + round(cx / dx), N0 + round(cy / dx)
            win = np.s_[i - h:i + h + 1, j - h:j + h + 1]
            m = (X[win] - cx) ** 2 + (Y[win] - cy) ** 2 <= (D / 2) ** 2
            w, ax_, ay_ = I[win][m], gx[win][m], gy[win][m]
            out.append([ax_.mean(), ay_.mean(),
                        (w * ax_).sum() / w.sum(),
                        (w * ay_).sum() / w.sum()])
    return np.array(out) / K            # (n_sub, 4): Gx, Gy, Cx, Cy


def corr_width(I, X, Y, dx):
    """1/e half-width (m) of the intensity autocovariance inside r < R_TILE."""
    m = X ** 2 + Y ** 2 <= R_TILE ** 2
    d = np.where(m, I - I[m].mean(), 0.0)
    a = np.fft.fftshift(np.abs(np.fft.ifft2(np.abs(np.fft.fft2(d,
                        s=(2 * d.shape[0],) * 2)) ** 2)))
    c = a.shape[0] // 2
    row = (a[c, c:] + a[c:, c]) / (2 * a[c, c])
    j = np.argmax(row < 1 / math.e)
    return dx * (j - 1 + (row[j - 1] - 1 / math.e) / (row[j - 1] - row[j]))


def boot(fn, n, B=1000, seed=1):
    """The value of fn(all trials) and its bootstrap SE over the trials."""
    rng = np.random.default_rng(seed)
    idx = np.arange(n)
    reps = [fn(rng.choice(idx, n)) for _ in range(B)]
    return float(fn(idx)), float(np.std(reps))


def own_tt_eta(U, dx):
    """MMF eta of the 1 inch aperture, no tilt and with the G-tilt removed."""
    h = math.ceil(1.5 * APERTURE_M / 2 / dx) + 1
    U = np.pad(U, h)                    # a compact patch crop can be smaller
    c = U.shape[0] // 2
    crop = U[c - h:c + h + 1, c - h:c + h + 1]
    x = np.arange(-h, h + 1) * dx
    X, Y = np.meshgrid(x, x)
    m = X ** 2 + Y ** 2 <= (APERTURE_M / 2) ** 2
    gx, gy = grads(crop, dx)
    return mmf_eta(crop, dx), mmf_eta(
        crop * np.exp(-1j * (gx[m].mean() * X + gy[m].mean() * Y)), dx)


def analyse(cn2, n_px):
    camp = smoke_campaign(cn2, 10, "cupy", n_px)
    n = camp.n_stored
    a = analytic(cn2)
    F0 = camp.field(0)
    dx, N = F0.dx, F0.N
    x = (np.arange(N) - N // 2) * dx
    X, Y = np.meshgrid(x, x)
    r = np.hypot(X, Y)
    m_tile = r < R_TILE

    cent, corr, step, dark, pt, own = [], [], [], [], [], []
    tilts = {D: [] for D in D_SUB}
    mean_I = np.zeros_like(r)
    field0 = None
    for i in range(n):
        U = camp.field(i).field.astype(np.complex128)
        field0 = U if i == 0 else field0
        I = np.abs(U) ** 2
        mean_I += I / n
        P = I.sum()
        cent.append(((X * I).sum() / P, (Y * I).sum() / P))
        corr.append(corr_width(I, X, Y, dx))
        gx, gy = grads(U, dx)
        step.append(np.abs(np.r_[gx[m_tile], gy[m_tile]]) * dx)
        dark.append(np.mean(I[m_tile] < 0.05 * I[m_tile].mean()))
        pt.append(I[r < 0.1])
        own.append(own_tt_eta(U, dx))
        for D in D_SUB:
            tilts[D].append(sub_tilts(U, dx, X, Y, D))
    cent, pt, own = np.array(cent), np.array(pt), np.array(own)
    corr, dark = np.array(corr), np.array(dark)
    step = np.concatenate(step)
    tilts = {D: np.array(t) for D, t in tilts.items()}
    det = camp.scenario.rx_terminal.detector
    p = np.array([t.collected_power for t in
                  camp.load(n, fields=False).trials])
    eta_pkg_tt = camp.recouple_compensated([TipTilt()], det)

    def tilt_rms(D, ix):
        t = tilts[D][ix]
        t = t - t.mean(axis=0)
        return math.sqrt((t[:, :, :2] ** 2).sum(axis=0).mean() / (len(ix) - 1))

    def si(ix):
        q = pt[ix]
        return q.var() / q.mean() ** 2

    def dbm(v):
        return lambda ix: -10 * math.log10(v[ix].mean())

    st = {
        "dx_mm": (1e3 * dx, 0.0),
        "ap_px": (APERTURE_M / dx, 0.0),
        "wander_mm": boot(lambda ix: 1e3 * cent[ix].std(axis=0, ddof=1).mean(), n),
        "si_point": boot(si, n),
        "speckle_mm": boot(lambda ix: 1e3 * corr[ix].mean(), n),
        "speckle_px": boot(lambda ix: corr[ix].mean() / dx, n),
        "step_p99_rad": (float(np.percentile(step, 99)), 0.0),
        "step>2.8_pct": (100 * float(np.mean(step > 2.8)), 0.0),
        "dark_pct": boot(lambda ix: 100 * dark[ix].mean(), n),
        "collected_db": boot(dbm(p), n),
        "bucket_s2": boot(lambda ix: p[ix].var() / p[ix].mean() ** 2, n),
        "eta_open": boot(lambda ix: own[ix, 0].mean(), n),
        "eta_tt_G": boot(lambda ix: own[ix, 1].mean(), n),
        "eta_tt_pkg": boot(lambda ix: eta_pkg_tt[ix].mean(), n),
        "tt_gain_db": boot(lambda ix: 10 * math.log10(
            own[ix, 1].mean() / own[ix, 0].mean()), n),
        "fibre_open_db": boot(dbm(p * own[:, 0]), n),
        "fibre_tt_db": boot(dbm(p * own[:, 1]), n),
    }
    for D in D_SUB:
        st[f"tilt{D / INCH:.0f}in_urad"] = boot(
            lambda ix, D=D: 1e6 * tilt_rms(D, ix), n)
    return dict(a=a, x=x, field=field0, mean_I=mean_I, dx=dx, X=X, Y=Y,
                st=st, n=n, n_px=camp.grid.n)


def reference(k, a):
    """The analytic value of one metric, or None."""
    if k.startswith("tilt"):
        return 1e6 * a["aoa"][INCH * float(k[4:k.index("in")])]
    return {"wander_mm": 1e3 * a["wander_diff"], "si_point": a["si"],
            "speckle_mm": 1e3 * a["rho0"], "collected_db": 35.64}.get(k)


def table(res):
    for cn2 in CN2:
        cols = res[cn2]
        a = cols[0]["a"]
        print(f"\n=== Cn2 {cn2:.0e}  sigma_R^2 {a['s2r']:.2f}  "
              f"(value +/- bootstrap SE over {cols[0]['n']} trials)")
        print(f"  {'metric':15s}" + "".join(f"{d['n_px']:>18d} px"
                                            for d in cols) + f"{'analytic':>11s}")
        for k in cols[0]["st"]:
            row = ""
            for d in cols:
                v, e = d["st"][k]
                row += (f"{v:12.3f} +/-{e:6.3f}" if e else f"{v:21.3f}")
            ref = reference(k, a)
            print(f"  {k:15s}{row}" + (f"{ref:11.2f}" if ref is not None else ""))


def figures(res, tag):
    os.makedirs(FIG, exist_ok=True)
    fig, ax = plt.subplots(len(res), 4, figsize=(17, 4.2 * len(res)))
    for i, (cn2, d) in enumerate(res.items()):
        x, U = d["x"], d["field"]
        I = np.abs(U) ** 2
        e = [x[0], x[-1], x[0], x[-1]]
        ax[i, 0].imshow(10 * np.log10(I / I.max() + 1e-6), extent=e,
                        origin="lower", vmin=-30, cmap="inferno")
        for rad, sty in ((d["a"]["W"], "c--"), (d["a"]["W_LT"], "w:")):
            t = np.linspace(0, 2 * np.pi, 200)
            ax[i, 0].plot(rad * np.cos(t), rad * np.sin(t), sty, lw=1)
        ax[i, 0].set_title(f"Cn2 {cn2:.0e}, trial 0, |U|^2 dB\n"
                           "cyan: vacuum W, white: W_LT")
        z = np.abs(x) <= 0.06
        sub = np.ix_(z, z)
        ez = [x[z][0], x[z][-1], x[z][0], x[z][-1]]
        t = np.linspace(0, 2 * np.pi, 100)
        for j, (img, cm, ttl) in enumerate((
                (I[sub], "inferno", "|U|^2, +/-60 mm"),
                (np.angle(U[sub]), "twilight", "wrapped phase, +/-60 mm"))):
            ax[i, j + 1].imshow(img, extent=ez, origin="lower", cmap=cm,
                                interpolation="nearest")
            ax[i, j + 1].plot(APERTURE_M / 2 * np.cos(t),
                              APERTURE_M / 2 * np.sin(t), "c-", lw=1.5)
            ax[i, j + 1].set_title(f"{ttl} (pixel {1e3 * d['dx']:.1f} mm)\n"
                                   "cyan: 1 inch aperture")
        rr = np.hypot(d["X"], d["Y"]).ravel()
        b = np.linspace(0, 1.5, 40)
        prof = np.histogram(rr, b, weights=d["mean_I"].ravel())[0] / \
            np.histogram(rr, b)[0]
        rc = 0.5 * (b[1:] + b[:-1])
        ax[i, 3].semilogy(rc, prof / prof[0], "k.-", label="field mean (10)")
        for rad, lab in ((d["a"]["W"], "vacuum W"), (d["a"]["W_LT"], "W_LT")):
            ax[i, 3].semilogy(rc, np.exp(-2 * (rc ** 2 - rc[0] ** 2) / rad ** 2),
                              label=f"{lab} {rad:.2f} m")
        ax[i, 3].set_ylim(1e-3, 2)
        ax[i, 3].set_xlabel("r (m)")
        ax[i, 3].legend(fontsize=8)
        ax[i, 3].set_title("mean irradiance profile")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, f"smoke_fields_{tag}.png"), dpi=90)
    plt.close(fig)


def tilt_figure(res):
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    Ds = np.geomspace(0.02, 0.12, 30)
    Dm = np.array(D_SUB)
    for c, cn2 in zip(("C0", "C1", "C2"), CN2):
        a = res[cn2][0]["a"]
        ax.loglog(Ds * 1e3, [1e6 * math.sqrt(aperture_arrival_angle_variance(
            D, 2.1 * a["rho0"], LAM, L0=L0)) for D in Ds], c + "-",
            label=f"{cn2:.0e} Andrews AoA")
        for d, mk, off in zip(res[cn2], "os", (0.97, 1.03)):
            v = [d["st"][f"tilt{D / INCH:.0f}in_urad"] for D in Dm]
            ax.errorbar(Dm * 1e3 * off, [q[0] for q in v], [q[1] for q in v],
                        fmt=c + mk, ms=5, capsize=2,
                        label=f"field {d['n_px']} px")
    ax.set_xlabel("sub-aperture D (mm)")
    ax.set_ylabel("G-tilt rms per axis (urad)")
    ax.set_title("o: 2048 px, s: 4096 px (bootstrap SE)")
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "smoke_tilt_grids.png"), dpi=110)


def main():
    res = {cn2: [analyse(cn2, n_px) for n_px in (None, 4096)] for cn2 in CN2}
    table(res)
    for j, tag in enumerate(("n2048", "n4096")):
        figures({c: res[c][j] for c in CN2}, tag)
    tilt_figure(res)


if __name__ == "__main__":
    main()

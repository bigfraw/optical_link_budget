"""Does a STACK of phase screens hold the statistics of ONE screen? Phase only.

THE QUESTION (owner, 2026-09-04). The FFT-plus-subharmonic screen generator is
close to Kolmogorov, but not faithful: it misses the low-frequency (tip and
tilt) band below the grid fundamental (backlog 2-N2, validation/screens/).
When a split-step plan divides the SAME turbulence among MORE screens, does the
stack lose MORE of that band than one screen of the same composite r0? If it
does, a screen-count sweep carries a generator bias, and a many-screen plan
reads a little LESS effective turbulence and a little LESS fibre-coupling loss.

THE TEST. No propagation. Each configuration is a list of per-screen r0 values
taken from the plans of the tail-convergence study
(validation/tail_convergence/). The script draws K stacks with the production
`ScreenFactory` on the pinned 1024 px grid of that study, SUMS the screens of a
stack into one phase map, and measures three things against the Kolmogorov
theory of the COMPOSITE r0:

  1. The phase structure function D(r) at six separations, averaged over the
     whole grid (both axes). Theory: D(r) = 6.88 (r/r0)^(5/3), Fried,
     DOI 10.1364/JOSA.56.001372.
  2. Delta1, the piston-removed phase variance over a D = 0.7 m aperture.
     Theory: 1.0299 (D/r0)^(5/3), Noll, DOI 10.1364/JOSA.66.000207, Table IV.
  3. Delta3, the tip-tilt-removed phase variance over the same aperture.
     Theory: 0.134 (D/r0)^(5/3), Noll, the same table.

Each number is a RATIO to its theory. A ratio of 1.0 is a faithful screen. The
25 disjoint 0.7 m apertures that tile the grid interior give the aperture
statistics of each draw; the standard error comes from the spread of the
per-draw means over the K draws, so the correlation between the apertures of
one draw does not fake a small bar.

THE READING. Delta3 near 1.0 with Delta1 under 1.0 says the missing power is
ALL tip and tilt. A Delta1 ratio that does not change with the screen count
says the deficit is a fixed FRACTION of each screen (Kolmogorov is scale-free,
so the fraction depends on D over the grid side and on the subharmonic depth,
not on r0), and stacking is then unbiased.

The composite r0 of a stack: r0c = (SUM r0_i^(-5/3))^(-3/5), Andrews and
Phillips, DOI 10.1117/3.626196, Ch. 12, Eq. (23).

THE OUTER SCALE (`--L0`, owner-requested 2026-09-04). The default reference is
Kolmogorov with L0 = inf, which the production screens also claim (`L0_m =
inf`). But a grid of side L with three subharmonic levels holds scales up to
about 27 L and nothing beyond, so against an INFINITE outer scale a "deficit"
is built in. With `--L0 <m>` the screens are drawn with that von Karman outer
scale AND the reference becomes the von Karman theory at the same L0:

  D(r) = 2 [B(0) - B(r)], B from the closed form of Assemat and Wilson,
         DOI 10.1364/OE.14.000988, Eq. (5)
         (validation.screens.helpers.vk_covariance_closed);
  Delta1 = INT 2 pi f Phi(f) [1 - (2 J1(pi D f) / (pi D f))^2] df, the
         piston filter of Noll, DOI 10.1364/JOSA.66.000207, Eq. (8);
  Delta3 = Delta1 - 2 <a2^2>, the per-axis Z-tilt variance from
         validation.screens.helpers.tilt_filter_variance (the same Eq. (8)),

with Phi(f) the von Karman phase PSD of Schmidt, DOI 10.1117/3.866274, Ch. 9,
Eq. (9.50), printed p. 161 (olb.waveoptics.schmidt.turbulence). The script
checks that route against the Kolmogorov closed forms at a large L0 before it
draws a screen. So a ratio of 1.0 at a finite L0 says the generator is
faithful to the spectrum it was asked for, and any residual is a REAL
generator deficit. See backlog 2-P5.

VALIDATION ONLY. It reads the production generator and it changes no olb
module. Run it from the repository root:

    python -m validation.screen_stacking.screen_stacking [--draws 100] [--L0 25]
"""

import argparse
import json
import os
import time
import warnings

import numpy as np
from scipy.special import j1

from olb.waveoptics.schmidt.turbulence import von_karman_phase_psd
from olb.waveoptics.turbulence.screens import ScreenFactory
from validation.screens.helpers import (F_GRID, _trapz, tilt_filter_variance,
                                        vk_covariance_closed)
from validation.tail_convergence.tail_convergence import (build_cases,
                                                          scenario_and_geometry)

HERE = os.path.dirname(os.path.abspath(__file__))

APERTURE_M = 0.7
SEPARATIONS_M = (0.05, 0.1, 0.2, 0.35, 0.7, 1.4)
SEED = 1234

# Noll, DOI 10.1364/JOSA.66.000207, Table IV: the residual phase variance after
# the first J Zernike modes are removed, in units of (D/r0)^(5/3).
NOLL_DELTA1 = 1.0299
NOLL_DELTA3 = 0.134
# Fried, DOI 10.1364/JOSA.56.001372: the phase structure function constant.
FRIED_D = 6.88


def composite_r0(r0s):
    """Give the composite r0 of a stack (Andrews Ch. 12, Eq. (23))."""
    return float(np.sum(np.asarray(r0s, dtype=float) ** (-5 / 3)) ** (-3 / 5))


def reference(r0c, L0_m, seps_m, aperture_m):
    """Give the theory of D(r), Delta1 and Delta3 for one composite r0.

    An infinite L0 takes the Kolmogorov closed forms (Fried; Noll, Table IV).
    A finite L0 integrates the von Karman phase PSD through the Noll piston
    and tilt filters (Noll, DOI 10.1364/JOSA.66.000207, Eq. (8)) on the
    F_GRID of validation.screens.helpers, and takes D(r) from the closed-form
    covariance (Assemat and Wilson, DOI 10.1364/OE.14.000988, Eq. (5)).

    Returns:
        The triple (D(r) array [rad^2], Delta1 [rad^2], Delta3 [rad^2]).
    """
    seps_m = np.asarray(seps_m, dtype=float)
    dr = (aperture_m / r0c) ** (5 / 3)
    if not np.isfinite(L0_m):
        return (FRIED_D * (seps_m / r0c) ** (5 / 3),
                NOLL_DELTA1 * dr, NOLL_DELTA3 * dr)
    psd = lambda f: von_karman_phase_psd(f, r0c, L0_m)
    b = vk_covariance_closed(np.concatenate(([0.0], seps_m)), r0c, L0_m)
    d_struct = 2.0 * (b[0] - b[1:])
    x = np.pi * aperture_m * F_GRID
    piston = (2.0 * j1(x) / x) ** 2
    delta1 = float(_trapz(2.0 * np.pi * F_GRID * psd(F_GRID) * (1.0 - piston),
                          F_GRID))
    delta3 = delta1 - 2.0 * tilt_filter_variance(psd, aperture_m)
    return d_struct, delta1, delta3


def check_reference(aperture_m):
    """Assert the finite-L0 route meets the Kolmogorov forms at a large L0.

    The von Karman forms approach Kolmogorov SLOWLY, as a power (r/L0)^(1/3)
    of the ratio (the small-argument expansion of the covariance of Assemat
    and Wilson, DOI 10.1364/OE.14.000988, Eq. (5)); at L0 = 1e4 m the
    structure function still reads 0.96 of Kolmogorov at r = 0.2 m. At
    L0 = 1e7 m every correction is under 1 percent for r <= 0.7 m, and F_GRID
    (F_MIN = 1e-8) still resolves that outer scale, so every ratio is held to
    1.5 percent.
    """
    r0, L0 = 0.1, 1.0e7
    seps = np.array([0.05, 0.2])
    d_vk, d1_vk, d3_vk = reference(r0, L0, seps, aperture_m)
    d_k, d1_k, d3_k = reference(r0, np.inf, seps, aperture_m)
    r_d = d_vk / d_k
    assert np.all(np.abs(r_d - 1.0) < 0.015), r_d
    assert abs(d3_vk / d3_k - 1.0) < 0.015, (d3_vk, d3_k)
    assert abs(d1_vk / d1_k - 1.0) < 0.015, (d1_vk, d1_k)
    return r_d, d1_vk / d1_k, d3_vk / d3_k


def configurations():
    """Give the named per-screen r0 lists, from the tail-convergence plans."""
    scn, geom = scenario_and_geometry(30.0)
    c = build_cases(scn, geom, ["pin05", "pin09", "pin25", "gnd09x4"])
    grid = c["pin09"]["grid"]
    cfg = {
        "ground: 1 screen (pin09 bottom)": c["pin09"]["plan"].r0_m[-1:],
        "ground: 4 sub-screens (gnd09x4)": c["gnd09x4"]["plan"].r0_m[-4:],
        "whole plan: 5 screens (pin05)": c["pin05"]["plan"].r0_m,
        "whole plan: 9 screens (pin09)": c["pin09"]["plan"].r0_m,
        "whole plan: 25 screens (pin25)": c["pin25"]["plan"].r0_m,
    }
    return grid, {k: np.asarray(v, dtype=float) for k, v in cfg.items()}


def aperture_masks(n, pixel_m, aperture_m):
    """Give the disjoint aperture discs that tile the grid, and the fit basis."""
    yy, xx = np.mgrid[:n, :n]
    ap_px = int(round(aperture_m / pixel_m))
    n_ap = n // ap_px
    out = []
    for i in range(n_ap):
        for j in range(n_ap):
            cy, cx = (i + 0.5) * ap_px, (j + 0.5) * ap_px
            m = ((yy - cy) ** 2 + (xx - cx) ** 2) * pixel_m ** 2 \
                <= (aperture_m / 2) ** 2
            x = (xx[m] - xx[m].mean()) * pixel_m
            y = (yy[m] - yy[m].mean()) * pixel_m
            basis = np.stack([np.ones_like(x), x, y], axis=1)
            out.append((m, basis))
    return out


def measure(fac, r0s, masks, seps_px, n_draws, seed):
    """Measure one configuration. Give the per-draw arrays."""
    rng = np.random.default_rng(seed)
    d_struct = np.zeros((n_draws, len(seps_px)))
    delta1 = np.zeros(n_draws)
    delta3 = np.zeros(n_draws)
    for k in range(n_draws):
        phi = np.sum(np.stack(fac.make_stack(r0s, rng)), axis=0)
        for i, s in enumerate(seps_px):
            d_struct[k, i] = (np.mean((phi[:, s:] - phi[:, :-s]) ** 2)
                              + np.mean((phi[s:, :] - phi[:-s, :]) ** 2)) / 2
        v1, v3 = [], []
        for m, basis in masks:
            v = phi[m]
            coef, *_ = np.linalg.lstsq(basis, v, rcond=None)
            v1.append(v.var())
            v3.append((v - basis @ coef).var())
        delta1[k] = np.mean(v1)
        delta3[k] = np.mean(v3)
    return d_struct, delta1, delta3


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--draws", type=int, default=100,
                    help="the stacks drawn for each configuration")
    ap.add_argument("--L0", type=float, default=np.inf,
                    help="the von Karman outer scale of the screens AND of "
                         "the reference [m]; the default is infinite")
    args = ap.parse_args()
    L0 = float(args.L0)
    tag = "" if not np.isfinite(L0) else f"_L0{L0:g}"

    log_path = os.path.join(HERE, f"screen_stacking{tag}.log")
    with open(log_path, "w", encoding="utf-8"):
        pass

    def say(text=""):
        print(text, flush=True)
        with open(log_path, "a", encoding="utf-8") as fh:
            fh.write(text + "\n")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        grid, cfg = configurations()
    n, dx = int(grid.n), float(grid.pixel_m)
    fac = ScreenFactory(n, dx, L0_m=L0)
    masks = aperture_masks(n, dx, APERTURE_M)
    seps_m = np.asarray(SEPARATIONS_M)
    seps_px = np.round(seps_m / dx).astype(int)
    r_d, r_1, r_3 = check_reference(APERTURE_M)

    say("The phase-screen STACKING test (phase only, no propagation)")
    say(f"grid {n} px, {dx * 1e3:.2f} mm pixel; {len(masks)} apertures of "
        f"{APERTURE_M} m per draw; {args.draws} draws per configuration")
    say(f"outer scale L0 = {L0:g} m (the screens AND the reference)")
    say(f"reference self-check at L0 = 1e4 m against Kolmogorov: D(r) "
        f"{r_d[0]:.4f} / {r_d[1]:.4f}, Delta1 {r_1:.4f}, Delta3 {r_3:.4f}")
    say("Every number is a RATIO to the theory of the composite r0 at this L0 "
        "(Fried / Noll when L0 is infinite, von Karman otherwise). +- is the "
        "standard error over the draws.")
    say()
    hdr = (f"{'configuration':<34s}{'r0c[cm]':>8s}"
           + "".join(f"{f'D({s:g})':>9s}" for s in seps_m)
           + f"{'Delta1':>14s}{'Delta3':>14s}")
    say(hdr)
    say("-" * len(hdr))

    rows = []
    for name, r0s in cfg.items():
        t0 = time.time()
        r0c = composite_r0(r0s)
        d_s, d1, d3 = measure(fac, r0s, masks, seps_px, args.draws, SEED)
        th_d, th_1, th_3 = reference(r0c, L0, seps_m, APERTURE_M)
        ratio_d = d_s.mean(axis=0) / th_d
        se_d = d_s.std(axis=0, ddof=1) / np.sqrt(args.draws) / th_d
        r1 = d1.mean() / th_1
        e1 = d1.std(ddof=1) / np.sqrt(args.draws) / th_1
        r3 = d3.mean() / th_3
        e3 = d3.std(ddof=1) / np.sqrt(args.draws) / th_3
        say(f"{name:<34s}{r0c * 100:8.2f}"
            + "".join(f"{v:9.3f}" for v in ratio_d)
            + f"{r1:8.3f}+-{e1:5.3f}{r3:8.3f}+-{e3:5.3f}"
            f"   ({time.time() - t0:.0f} s)")
        rows.append({"configuration": name, "n_screens": int(r0s.size),
                     "r0_screens_m": r0s.tolist(), "r0_composite_m": r0c,
                     "structure_ratio": ratio_d.tolist(),
                     "structure_ratio_se": se_d.tolist(),
                     "delta1_ratio": float(r1), "delta1_ratio_se": float(e1),
                     "delta3_ratio": float(r3), "delta3_ratio_se": float(e3)})
    say()

    # The verdict: the stack against the single screen, and the count series.
    by = {r["configuration"]: r for r in rows}
    g1 = by["ground: 1 screen (pin09 bottom)"]
    g4 = by["ground: 4 sub-screens (gnd09x4)"]
    p5 = by["whole plan: 5 screens (pin05)"]
    p25 = by["whole plan: 25 screens (pin25)"]

    def _cmp(a, b, label):
        d = b["delta1_ratio"] - a["delta1_ratio"]
        bar = float(np.hypot(a["delta1_ratio_se"], b["delta1_ratio_se"]))
        say(f"  {label}: Delta1 ratio {a['delta1_ratio']:.3f} -> "
            f"{b['delta1_ratio']:.3f}, a change of {d:+.3f} against a bar of "
            f"{bar:.3f} ({abs(d) / bar:.1f} sigma).")
        return {"delta": float(d), "bar": bar}

    say("THE VERDICT")
    verdict = {
        "ground_1_to_4": _cmp(g1, g4, "the ground layer, 1 -> 4 screens"),
        "plan_5_to_25": _cmp(p5, p25, "the whole plan, 5 -> 25 screens"),
        "single_screen_deficit": {
            "delta1_ratio_range": [min(r["delta1_ratio"] for r in rows),
                                   max(r["delta1_ratio"] for r in rows)],
            "delta3_ratio_range": [min(r["delta3_ratio"] for r in rows),
                                   max(r["delta3_ratio"] for r in rows)]},
    }
    say(f"  every configuration: Delta1 ratio "
        f"{verdict['single_screen_deficit']['delta1_ratio_range'][0]:.3f} to "
        f"{verdict['single_screen_deficit']['delta1_ratio_range'][1]:.3f}, "
        f"Delta3 ratio "
        f"{verdict['single_screen_deficit']['delta3_ratio_range'][0]:.3f} to "
        f"{verdict['single_screen_deficit']['delta3_ratio_range'][1]:.3f}. "
        + ("Against an INFINITE outer scale the missing power is tip and "
           "tilt beyond the subharmonic reach (backlog 2-P5, 2-N2)."
           if not np.isfinite(L0) else
           f"The reference is von Karman at L0 = {L0:g} m, so any shortfall "
           "is a REAL generator deficit."))

    out = {"study": "screen_stacking", "n_draws": args.draws, "seed": SEED,
           "L0_m": None if not np.isfinite(L0) else L0,
           "grid_n": n, "pixel_m": dx, "aperture_m": APERTURE_M,
           "separations_m": list(SEPARATIONS_M), "rows": rows,
           "verdict": verdict}
    json_path = os.path.join(HERE, f"screen_stacking{tag}_results.json")
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    say(f"wrote {json_path}")


if __name__ == '__main__':
    main()

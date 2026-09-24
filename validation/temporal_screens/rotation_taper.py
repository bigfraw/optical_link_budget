"""Taper and detrend variants of the Fourier rotation of a cropped screen.

WHY. The rotated-strip time axis (`olb.waveoptics.turbulence.temporal`) makes
one frame from a thin strip: it crops a square patch of side n(1 + 2*margin),
it turns the patch back with `rotate_fourier` (the exact three-shear route of
Unser, Thevenaz and Yaroslavsky, DOI 10.1109/83.469963), and it keeps the
central n by n. The crop pixels drive the cost of a frame, so a SHORT margin is
worth money. Two errors set the margin:

  1. CORNER CLIPPING. Pure geometry. The kept square reads its source in the
     ROTATED square, whose axis-aligned footprint is n(|cos t| + |sin t|). A
     margin under 0.5*(footprint/n - 1) cuts real data away, and no window
     repairs that. It is 0.05 n at 6.5 deg and 0.21 n at 45 deg.
  2. SEAM RINGING. Each shear is a phase ramp on the 1-D transform of a line
     (the Fourier shift theorem; Schmidt, DOI 10.1117/3.866274, Ch. 2), so the
     FFT reads each line as PERIODIC. The step between the two ends of a crop
     rings inward with a sinc tail, which falls as 1/distance.

THE VARIANTS. Each one attacks the step at the ends of the lines.

  V0  the plain `rotate_fourier` call (the control).
  V1  DETREND. Fit the mean and the tilt plane of the crop, rotate the
      residual, and add the plane back ANALYTICALLY. A plane g.q + c maps to
      (R^T g).p + c under out(p) = patch(R p), so the plane needs no FFT. Most
      of the power of a von Karman screen is at the largest scales, so the
      plane holds most of the end step.
  V2  TUKEY TAPER. Multiply the crop by a square (Chebyshev) cosine-tapered
      window that is 1 over the axis-aligned FOOTPRINT of the rotated square
      and rolls to 0 at the edge of the crop, then rotate. The kept square
      reads only the flat part, so the taper does not touch the answer; it
      only removes the end step. It needs a margin over the geometric one.
  V3  V1 then V2.
  V4  PER-SHEAR TAPER. Taper the ends of every line with the same 1-D cosine
      rule BEFORE each of the three shear transforms.
  V5  REFLECT PAD. Mirror the crop out to twice its side, rotate, and cut the
      centre. The mirror makes the line continuous at the old edge, so the step
      goes away without a change of the data.

THE TEST. One 1024 px von Karman screen from `ScreenFactory` (r0 = 10 cm,
L0 = 25 m, dx = 1 cm), 4 seeds, n = 256, theta = 6.5 and 45 deg, margins
0.0625, 0.125, 0.25 and 0.5 n. The REFERENCE of a variant is the same rotation
of the FULL periodic 1024 px screen, central n by n: that field IS periodic, so
it carries no seam.

THE MEASURES. The rms error of the kept centre; the rms error against the
distance from the edge of the kept square; the rms error over the outer 16 px
ring (the worst part, and the part a window hurts first); and the ratio of the
phase structure function D(lag) of the kept centre to the reference at the lags
1, 4 and 16 px, on both axes. A window that flattens the field near the flat
top would show up as a D ratio under 1 even when the rms error looks small, so
the two measures are read together.

Run:  python -m validation.temporal_screens.rotation_taper
"""

import os
import time

import numpy as np

from olb.waveoptics.turbulence.screens import ScreenFactory
from olb.waveoptics.turbulence.temporal import _shear, rotate_fourier

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
FIGS = os.path.join(HERE, "figures")

N_BIG = 1024        # the full periodic screen, in px.
N_FRAME = 256       # the kept frame side n, in px.
DX = 0.01           # the pixel pitch, in m.
R0 = 0.10           # the Fried parameter of the screen, in m.
L0 = 25.0           # the outer scale, in m (the site operating value).
SEEDS = (3, 11, 29, 47)
ANGLES = (6.5, 45.0)
MARGINS = (0.0625, 0.125, 0.25, 0.5)
LAGS = (1, 4, 16)
RING = 16           # the width of the outer ring, in px.

EDGE_RMS_LIMIT = 1e-2      # rad. The owner target of the edge ring.
D_RATIO_LIMIT = 0.01       # 1 percent on |D ratio - 1|.


# ---------------------------------------------------------------- helpers ---

def centre(a, m):
    """Cut the central m by m square about the rotation centre.

    The three shears turn about the pixel `m // 2`, so the cut uses the same
    centre and not the middle of the array.
    """
    k = a.shape[0] // 2
    return a[k - m // 2:k + m // 2, k - m // 2:k + m // 2]


def footprint(n, theta_rad):
    """Give the axis-aligned footprint of the rotated n by n square, in px."""
    return n * (abs(np.cos(theta_rad)) + abs(np.sin(theta_rad)))


def min_geometric_margin(theta_rad):
    """Give the margin, in units of n, that the corner geometry alone needs."""
    return 0.5 * (footprint(1.0, theta_rad) - 1.0)


def _taper_1d(m, flat_px):
    """Give the 1-D cosine (Tukey) profile of a line of m px.

    It is 1 where |c| <= flat_px / 2, with c the coordinate about the rotation
    centre `m // 2`, and it falls to 0 with a raised cosine at the ends of the
    line. `flat_px` is the footprint of the rotated square, so the kept data
    never leaves the flat top.
    """
    c = np.abs(np.arange(m) - m // 2)
    half = 0.5 * float(flat_px)
    edge = 0.5 * m
    if half >= edge:
        raise ValueError("the margin is too small: the footprint fills the "
                         "crop, so no taper fits.")
    t = np.clip((c - half) / (edge - half), 0.0, 1.0)
    return 0.5 * (1.0 + np.cos(np.pi * t))


def _plane_fit(a):
    """Give the least-squares mean and tilt of `a` about the rotation centre.

    The coordinates x and y, once their own mean is removed, are orthogonal to
    each other and to the constant on a rectangular grid, so the fit is three
    one-line sums and it needs no matrix solve.

    Returns:
        (c0, gx, gy), the plane c0 + gx*x + gy*y with x and y the coordinates
        about the pixel `m // 2`.
    """
    m = a.shape[0]
    c = (np.arange(m) - m // 2).astype(float)
    cc = c - c.mean()
    denom = float(np.sum(cc ** 2)) * m
    gx = float(np.sum(a * cc[None, :])) / denom
    gy = float(np.sum(a * cc[:, None])) / denom
    c0 = float(a.mean()) - gx * c.mean() - gy * c.mean()
    return c0, gx, gy


def _plane(m, c0, gx, gy):
    """Raster the plane c0 + gx*x + gy*y about the pixel `m // 2`."""
    c = (np.arange(m) - m // 2).astype(float)
    return c0 + gx * c[None, :] + gy * c[:, None]


def _three_shear(patch, theta_rad, flat_px=None):
    """Rotate with the three shears, with an OPTIONAL taper before each shear.

    With `flat_px` None this is `rotate_fourier` itself, line for line, for an
    angle inside [-45, 45] deg (no quarter turn is needed there). The
    `__main__` gate asserts that the two agree bit for bit.

    Args:
        patch:     a real square array.
        theta_rad: the angle, in rad, inside [-45, 45] deg.
        flat_px:   the flat width of the 1-D taper, in px, or None for no
                   taper (V4 against V0).

    Returns:
        A real array of the shape of `patch`.
    """
    m = patch.shape[0]
    c = (np.arange(m) - m // 2).astype(patch.dtype)
    t = np.tan(0.5 * theta_rad)
    s = np.sin(theta_rad)
    w = None if flat_px is None else _taper_1d(m, flat_px)
    out = patch
    for shift, axis in ((t * c, 1), (-s * c, 0), (t * c, 1)):
        if w is not None:
            out = out * (w[None, :] if axis == 1 else w[:, None])
        out = _shear(out, shift, axis)
    return out


def _tukey_2d(m, flat_px):
    """Give the square (Chebyshev) 2-D Tukey window of the crop.

    It is the outer maximum of the 1-D profile on the two axes, so its flat top
    is the axis-aligned SQUARE of side `flat_px`, which holds the whole rotated
    square that the kept centre reads.
    """
    w = _taper_1d(m, flat_px)
    return np.minimum(w[None, :], w[:, None])


# --------------------------------------------------------------- variants ---

def variant(name, crop, theta_rad, n):
    """Rotate one crop with one variant and give the kept n by n centre.

    Args:
        name:      "V0" to "V5".
        crop:      the square crop, of side m.
        theta_rad: the angle of the rotation, in rad.
        n:         the side of the kept frame, in px.

    Returns:
        The kept n by n array, or None when the margin cannot hold the
        variant (a taper needs room past the rotation footprint).
    """
    m = crop.shape[0]
    flat = footprint(n, theta_rad)
    if name in ("V2", "V3", "V4") and flat >= m:
        return None                      # no room for the roll-off.
    if name == "V0":
        return centre(rotate_fourier(crop, theta_rad), n)
    if name == "V1":
        c0, gx, gy = _plane_fit(crop)
        res = crop - _plane(m, c0, gx, gy)
        # out(p) = patch(R p), so a plane gradient g maps to R^T g.
        cs, sn = np.cos(theta_rad), np.sin(theta_rad)
        back = _plane(m, c0, gx * cs + gy * sn, -gx * sn + gy * cs)
        return centre(rotate_fourier(res, theta_rad) + back, n)
    if name == "V2":
        return centre(rotate_fourier(crop * _tukey_2d(m, flat), theta_rad), n)
    if name == "V3":
        c0, gx, gy = _plane_fit(crop)
        res = (crop - _plane(m, c0, gx, gy)) * _tukey_2d(m, flat)
        cs, sn = np.cos(theta_rad), np.sin(theta_rad)
        back = _plane(m, c0, gx * cs + gy * sn, -gx * sn + gy * cs)
        return centre(rotate_fourier(res, theta_rad) + back, n)
    if name == "V4":
        return centre(_three_shear(crop, theta_rad, flat_px=flat), n)
    if name == "V5":
        pad = m // 2
        wide = np.pad(crop, pad, mode="reflect")
        return centre(rotate_fourier(wide, theta_rad), n)
    raise ValueError(f"unknown variant {name!r}")


VARIANTS = ("V0", "V1", "V2", "V3", "V4", "V5")


# ---------------------------------------------------------------- measures ---

def structure_ratios(a, ref, lags=LAGS, mask=None):
    """Give max |D(lag)/D_ref(lag) - 1| over the two axes, for each lag.

    D(lag) is the mean square of the phase difference at that lag, the
    structure function of the screen (Andrews and Phillips,
    DOI 10.1117/3.626196, Ch. 3).

    `mask` keeps only the pairs whose TWO ends are True. The outer-ring mask
    is the test of the caveat of a window: a taper flattens the field near its
    flat top, and that shows up in D before it shows up in the rms error.
    """
    out = []
    for k in lags:
        worst = 0.0
        for axis in (0, 1):
            sl_a = (slice(k, None), slice(None))
            sl_b = (slice(None, -k), slice(None))
            if axis == 1:
                sl_a, sl_b = sl_a[::-1], sl_b[::-1]
            if mask is None:
                sel = Ellipsis
            else:
                sel = mask[sl_a] & mask[sl_b]
            d = float(np.mean(((a[sl_a] - a[sl_b]) ** 2)[sel]))
            d_ref = float(np.mean(((ref[sl_a] - ref[sl_b]) ** 2)[sel]))
            worst = max(worst, abs(d / d_ref - 1.0))
        out.append(worst)
    return out


def edge_distance(n):
    """Give the distance of every pixel from the edge of the n by n square."""
    yy, xx = np.indices((n, n))
    return np.minimum.reduce([xx, yy, n - 1 - xx, n - 1 - yy])


def measure(err, dist, n):
    """Give (rms, edge-ring rms, the profile against the distance)."""
    rms = float(np.sqrt(np.mean(err ** 2)))
    ring = float(np.sqrt(np.mean(err[dist < RING] ** 2)))
    prof = np.array([np.sqrt(np.mean(err[dist == k] ** 2))
                     for k in range(n // 2)])
    return rms, ring, prof


# ------------------------------------------------------------------- study ---

def run():
    """Run the whole sweep and give the table, the profiles and the maps."""
    n = N_FRAME
    dist = edge_distance(n)
    # key -> lists over the seeds; the table averages the SQUARES.
    acc = {}
    for seed in SEEDS:
        screen = ScreenFactory(N_BIG, DX, L0_m=L0).make(
            R0, np.random.default_rng(seed))
        for deg in ANGLES:
            th = np.deg2rad(deg)
            ref = centre(rotate_fourier(screen, th), n)
            for mg in MARGINS:
                m = int(round(n * (1.0 + 2.0 * mg)))
                m += m % 2
                crop = centre(screen, m)
                for name in VARIANTS:
                    got = variant(name, crop, th, n)
                    key = (name, deg, mg)
                    if got is None:
                        acc.setdefault(key, None)
                        continue
                    err = np.abs(got - ref)
                    rms, ring, prof = measure(err, dist, n)
                    ratios = structure_ratios(got, ref)
                    ring_ratios = structure_ratios(got, ref,
                                                   mask=dist < 2 * RING)
                    row = acc.setdefault(key, [])
                    row.append((rms, ring, prof, ratios, ring_ratios))
    return acc


def table_rows(acc):
    """Turn the accumulator into one row for each variant, angle and margin."""
    rows = []
    for name in VARIANTS:
        for deg in ANGLES:
            for mg in MARGINS:
                got = acc.get((name, deg, mg))
                if not got:
                    rows.append(dict(variant=name, angle=deg, margin=mg,
                                     rms=np.nan, ring=np.nan,
                                     d=[np.nan] * len(LAGS),
                                     d_ring=[np.nan] * len(LAGS), ok=False,
                                     note="no room for the taper"))
                    continue
                rms = float(np.sqrt(np.mean([g[0] ** 2 for g in got])))
                ring = float(np.sqrt(np.mean([g[1] ** 2 for g in got])))
                d = [float(np.max([g[3][i] for g in got]))
                     for i in range(len(LAGS))]
                d_ring = [float(np.max([g[4][i] for g in got]))
                          for i in range(len(LAGS))]
                ok = (ring < EDGE_RMS_LIMIT
                      and all(x < D_RATIO_LIMIT for x in d)
                      and all(x < D_RATIO_LIMIT for x in d_ring))
                rows.append(dict(variant=name, angle=deg, margin=mg, rms=rms,
                                 ring=ring, d=d, d_ring=d_ring, ok=ok,
                                 note=""))
    return rows


def write_csv(rows, path):
    """Write the table as a CSV."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    head = ["variant", "angle_deg", "margin_n", "rms_rad", "edge_ring_rms_rad"]
    head += [f"d{k}px_ratio_err" for k in LAGS]
    head += [f"d{k}px_ratio_err_ring" for k in LAGS] + ["pass", "note"]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(",".join(head) + "\n")
        for r in rows:
            fh.write(",".join([r["variant"], f"{r['angle']:g}",
                               f"{r['margin']:g}", f"{r['rms']:.4e}",
                               f"{r['ring']:.4e}"]
                              + [f"{x:.4e}" for x in r["d"]]
                              + [f"{x:.4e}" for x in r["d_ring"]]
                              + [str(int(r["ok"])), r["note"]]) + "\n")


def figures(acc, best):
    """Make the profile figure and the error-map grid of the best variant."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    os.makedirs(FIGS, exist_ok=True)
    n = N_FRAME

    # Figure 1: the rms error against the distance from the edge.
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0),
                             constrained_layout=True)
    styles = {0.125: "-", 0.25: "--"}
    colours = dict(zip(VARIANTS, plt.rcParams["axes.prop_cycle"].by_key()
                       ["color"]))
    for ax, deg in zip(axes, ANGLES):
        for name in VARIANTS:
            for mg, ls in styles.items():
                got = acc.get((name, deg, mg))
                if not got:
                    continue
                prof = np.sqrt(np.mean([g[2] ** 2 for g in got], axis=0))
                ax.semilogy(np.arange(prof.size), prof, ls,
                            color=colours[name], lw=1.3,
                            label=f"{name}, margin {mg:g} n")
        ax.axhline(EDGE_RMS_LIMIT, color="k", ls=":", lw=0.8)
        ax.set_title(f"theta {deg:g} deg (geometry alone needs "
                     f"{min_geometric_margin(np.deg2rad(deg)):.2f} n)",
                     fontsize=10)
        ax.set_xlabel("distance from the edge of the kept n x n square (px)")
        ax.set_ylabel("rms error (rad)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=6, ncol=2)
    fig.suptitle("Taper and detrend variants of the Fourier rotation: the "
                 f"error into the interior (n = {n} px, {len(SEEDS)} seeds, "
                 f"r0 = {R0 * 100:.0f} cm, L0 = {L0:.0f} m)")
    fig.savefig(os.path.join(FIGS, "rotation_taper_profiles.png"), dpi=130)
    plt.close(fig)

    # Figure 2: the error maps of the best variant, seed 0.
    screen = ScreenFactory(N_BIG, DX, L0_m=L0).make(
        R0, np.random.default_rng(SEEDS[0]))
    fig, axes = plt.subplots(len(ANGLES), len(MARGINS),
                             figsize=(15.0, 7.4), constrained_layout=True)
    im = None
    for i, deg in enumerate(ANGLES):
        th = np.deg2rad(deg)
        ref = centre(rotate_fourier(screen, th), n)
        for j, mg in enumerate(MARGINS):
            ax = axes[i, j]
            m = int(round(n * (1.0 + 2.0 * mg)))
            m += m % 2
            got = variant(best, centre(screen, m), th, n)
            ax.set_xticks([])
            ax.set_yticks([])
            if got is None:
                ax.set_title(f"theta {deg:g} deg, margin {mg:g} n\n"
                             "no room for the taper", fontsize=9)
                continue
            err = np.abs(got - ref)
            im = ax.imshow(np.maximum(err, 1e-6), norm=LogNorm(1e-4, 1e0),
                           cmap="magma", origin="lower")
            ax.set_title(f"theta {deg:g} deg, margin {mg:g} n\ncrop {m} px, "
                         f"rms {np.sqrt(np.mean(err ** 2)):.2e} rad",
                         fontsize=9)
    if im is not None:
        fig.colorbar(im, ax=axes, shrink=0.8,
                     label="|rotated crop - reference| (rad), log")
    fig.suptitle(f"The best variant ({best}): the error of the kept centre "
                 "against the rotation of the full periodic screen")
    fig.savefig(os.path.join(FIGS, "rotation_taper_maps.png"), dpi=130)
    plt.close(fig)


def timing():
    """Time one frame of V0 and of the winner at n = 512, margin 0.25 n."""
    n, mg, th = 512, 0.25, np.deg2rad(45.0)
    m = int(round(n * (1.0 + 2.0 * mg)))
    m += m % 2
    crop = ScreenFactory(m, DX, L0_m=L0).make(R0, np.random.default_rng(5))
    out = {}
    for name in VARIANTS:
        if variant(name, crop, th, n) is None:
            continue
        t0 = time.perf_counter()
        for _ in range(3):
            variant(name, crop, th, n)
        out[name] = (time.perf_counter() - t0) / 3.0
    return m, out


if __name__ == '__main__':
    # ---- 1. the local three-shear route IS `rotate_fourier` ----
    probe = np.random.default_rng(1).standard_normal((64, 64))
    for deg in (6.5, -20.0, 45.0):
        a = _three_shear(probe, np.deg2rad(deg))
        b = rotate_fourier(probe, np.deg2rad(deg))
        assert np.array_equal(a, b), (deg, float(np.abs(a - b).max()))

    # ---- 2. the analytic plane rotation of V1 is exact on a pure plane ----
    m_p = 64
    plane = _plane(m_p, 2.0, 0.03, -0.07)
    th_p = np.deg2rad(6.5)
    c0, gx, gy = _plane_fit(plane)
    assert abs(c0 - 2.0) < 1e-9 and abs(gx - 0.03) < 1e-12, (c0, gx, gy)
    cs, sn = np.cos(th_p), np.sin(th_p)
    want = _plane(m_p, 2.0, 0.03 * cs - 0.07 * sn, -0.03 * sn - 0.07 * cs)
    got = variant("V1", plane, th_p, 32)
    assert np.allclose(got, centre(want, 32), atol=1e-10), \
        float(np.abs(got - centre(want, 32)).max())

    # ---- 3. the sweep ----
    acc = run()
    rows = table_rows(acc)
    csv_path = os.path.join(DATA, "rotation_taper.csv")
    write_csv(rows, csv_path)

    # The winner: the variant with the lowest edge-ring rms at margin 0.25 n,
    # taken as the worst of the two angles.
    def worst_ring(name, mg):
        vals = [r["ring"] for r in rows
                if r["variant"] == name and r["margin"] == mg]
        return np.nan if any(np.isnan(v) for v in vals) else max(vals)

    ranked = [(worst_ring(v, 0.25), v) for v in VARIANTS if v != "V0"]
    ranked = [p for p in ranked if not np.isnan(p[0])]
    ranked.sort()
    best = ranked[0][1]

    # ---- 4. the winner beats the control at margin 0.25 n ----
    assert ranked[0][0] < worst_ring("V0", 0.25), (ranked, worst_ring("V0",
                                                                     0.25))

    figures(acc, best)
    m_t, times = timing()

    print(f"screen {N_BIG} px, frame n = {N_FRAME} px, dx = {DX * 100:.0f} cm,"
          f" r0 = {R0 * 100:.0f} cm, L0 = {L0:.0f} m, {len(SEEDS)} seeds")
    print("dD<k> is |D(k px) ratio - 1| over the whole frame, and dR<k> is the"
          " same over the outer 32 px band.")
    print(f"{'var':>4} {'deg':>5} {'marg':>6} {'rms':>10} {'ring':>10} "
          + " ".join(f"{'dD' + str(k):>8}" for k in LAGS) + " "
          + " ".join(f"{'dR' + str(k):>8}" for k in LAGS) + "  pass")
    for r in rows:
        if np.isnan(r["rms"]):
            print(f"{r['variant']:>4} {r['angle']:5g} {r['margin']:6g} "
                  f"{'-':>10} {'-':>10} "
                  + " ".join(f"{'-':>8}" for _ in LAGS * 2)
                  + f"   {r['note']}")
            continue
        print(f"{r['variant']:>4} {r['angle']:5g} {r['margin']:6g} "
              f"{r['rms']:10.3e} {r['ring']:10.3e} "
              + " ".join(f"{x:8.2e}" for x in r["d"]) + " "
              + " ".join(f"{x:8.2e}" for x in r["d_ring"])
              + f"  {'yes' if r['ok'] else 'no'}")
    print()
    print("smallest margin that passes (edge ring < "
          f"{EDGE_RMS_LIMIT:g} rad AND all |D ratio - 1| < "
          f"{D_RATIO_LIMIT * 100:g} percent):")
    for name in VARIANTS:
        for deg in ANGLES:
            good = [r["margin"] for r in rows
                    if r["variant"] == name and r["angle"] == deg and r["ok"]]
            print(f"  {name}  {deg:5g} deg  "
                  + (f"{min(good):g} n" if good else "none of the tested"))
    print()
    print(f"best variant: {best}")
    print(f"one frame at n = 512, margin 0.25 n (crop {m_t} px):")
    for name, t in times.items():
        print(f"  {name}  {t * 1e3:8.1f} ms  "
              f"({t / times['V0']:.2f} x V0)")
    print(f"table  {csv_path}")
    print(f"figures {FIGS}")
    print("self-check passed")

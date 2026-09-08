"""Read the stored terrestrial campaigns ONE time, and write a per-trial table.

This is step 1 and step 2 of backlog 1-9, the data half. The twelve
terrestrial fidelity-2 campaigns of `validation/terrestrial_campaigns/` store
the complex receive-plane field of every trial on a 5 cm disc. The fit of the
power distributions needs about twenty scalars for each trial, not the field.
So this script runs on the machine that holds the campaigns, it reads each
stored field one time, and it writes one small `.npz` table for each cell.
`fit_distributions.py` then reads those tables anywhere.

What the table holds, for each trial (see COLUMNS):

- the BUCKET power inside the 10 cm and the 5 cm receive aperture, in grid
  units (`Campaign.recollect`);
- the SINGLE-MODE-FIBRE coupling efficiency of the 10 cm and the 5 cm
  aperture, three ways: UNTRACKED at the nominal focal plane (the in-run
  quantity), TRACKED at the received-curvature focus (`defocus_m =
  curvature_focus_shift`, S. A. Self, DOI 10.1364/AO.22.000658), and with the
  TIP-TILT REMOVED by the perfect modal corrector of `olb.waveoptics.
  compensation` (the first 3 Noll modes, R. J. Noll,
  DOI 10.1364/JOSA.66.000207, sensed from the wrapped-gradient slopes of the
  field, which is the terrestrial sensing rule of the runner);
- the Noll tilt coefficients a2, a3 of the 10 cm and the 5 cm aperture, in
  radians, and the largest phase step per pixel (the slope aliasing test,
  `olb.waveoptics.compensation.max_abs_step`);
- the point irradiance: the mean |E|^2 inside an on-axis disc of radius
  POINT_RADIUS_M (at least the centre pixel), in grid units;
- the MULTIMODE light-bucket coupling of the 10 cm aperture into a 100 um
  core of numerical aperture 0.22 behind a 0.25 m focal length, so the optic
  cone (NA_optic = 0.20) sits inside the fibre acceptance and the core is the
  only gate (`olb.waveoptics.mmf.mmf_coupling_efficiency`).

The in-run `collected_power` and `smf_eta` of the trial are copied too, so
the post-hoc read is cross-checked against the run.

Run it from the repository root on the machine that holds the campaigns:

    python -m validation.fibre_fade_models.extract_trials --workers auto

It reads the campaign roots of `validation/terrestrial_campaigns/
run_campaigns.py` (the `OLB_TERRESTRIAL_CAMPAIGNS_ROOT` override applies).
The tables go to `validation/fibre_fade_models/trials/`.
"""

import argparse
import json
import os
import sys
import time

import numpy as np

from olb.models.coupling.terrestrial import curvature_focus_shift
from olb.terminal import MMF, SMF
from olb.waveoptics.compensation import (ApertureModes, circle, max_abs_step,
                                         wrapped_gradient)
from olb.waveoptics.turbulence.campaign import Campaign
from olb.waveoptics.turbulence.run import SLOPE_MIN_MODES, _PostTail
from validation.terrestrial_campaigns import run_campaigns as rc

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "trials")

# The receive apertures of the study. The 10 cm aperture is the clip of the
# campaign; the 5 cm aperture is a post-hoc crop inside the stored disc.
APERTURES_M = (0.10, 0.05)

# The point-irradiance disc. The 1-6 certification used an 8 mm disc (0.14 of
# the Fresnel scale at 2 km). A 10 km cell has a 3.8 to 6.4 mm pixel, so the
# disc holds one to four pixels there: a coarse point.
POINT_RADIUS_M = 4e-3

# The multimode light bucket of the 10 cm aperture. A 100 um core (radius
# 50 um) of NA 0.22 behind f = 0.25 m: NA_optic = 0.05 / 0.25 = 0.20 < 0.22,
# so the optic cone sits inside the fibre acceptance and the NA gate passes
# every ray. The diffraction spot radius lam f / (pi R) = 2.5 um sits far
# inside the core, so the coupling is the spot-walk against the core.
MMF_CORE_RADIUS_M = 50e-6
MMF_NA = 0.22
MMF_FOCAL_LENGTH_M = 0.25

COLUMNS = (
    "collected_power_stored",   # the in-run 10 cm bucket power
    "smf_eta_stored",           # the in-run 10 cm untracked SMF eta
    "bucket_10cm",              # the post-hoc 10 cm bucket power
    "bucket_5cm",
    "smf10_untracked",          # SMF eta, 10 cm, defocus 0
    "smf10_tracked",            # SMF eta, 10 cm, defocus = dz_curv
    "smf10_tiptilt",            # SMF eta, 10 cm, tilt removed, defocus 0
    "smf10_tiptilt_tracked",    # SMF eta, 10 cm, tilt removed, tracked
    "smf5_untracked",
    "smf5_tracked",
    "smf5_tiptilt",
    "a2_10cm",                  # Noll x tilt of the 10 cm aperture, rad
    "a3_10cm",                  # Noll y tilt of the 10 cm aperture, rad
    "a2_5cm",
    "a3_5cm",
    "max_step_10cm",            # the largest phase step per pixel, rad
    "point_irradiance",         # the mean |E|^2 in the on-axis disc
    "centre_irradiance",        # the |E|^2 of the axis pixel
    "mmf10",                    # MMF eta, 10 cm, the light bucket
)


class _Mask:
    """The point-irradiance pixels of one crop."""

    def __init__(self, ref, radius_m):
        Y, X = ref.mgrid_cartesian
        dist_sq = X ** 2 + Y ** 2
        keep = dist_sq <= radius_m ** 2
        if not keep.any():
            keep = dist_sq == dist_sq.min()
        self.point = keep
        self.centre = dist_sq == dist_sq.min()


class ExtractTrials:
    """The per-trial callable of the extraction (see the module docstring).

    The masks, the fibre modes, the defocus phases and the modal bases are
    built ONE time for each block, in the record context, and they are reused
    for every trial of the block. This is the `_CoupleTrials` pattern of
    `olb.waveoptics.turbulence.campaign`.
    """

    def __init__(self, dz_curv_m, obscuration_ratio=0.0):
        # dz_curv_m: the received-curvature focus shift of each aperture, a
        # dict aperture_m -> m.
        self.dz_curv_m = dict(dz_curv_m)
        self.obscuration_ratio = float(obscuration_ratio)

    def _build(self, rec):
        ctx = rec.context
        if "tails" in ctx:
            return ctx
        patch, lam = rec.patch, rec.lam
        tails, modes = {}, {}
        for D in APERTURES_M:
            dz = self.dz_curv_m[D]
            tails[(D, "bucket")] = _PostTail(patch, D, self.obscuration_ratio,
                                             lam, detector=None)
            tails[(D, "untracked")] = _PostTail(
                patch, D, self.obscuration_ratio, lam,
                detector=SMF(optimal_focus=True))
            tails[(D, "tracked")] = _PostTail(
                patch, D, self.obscuration_ratio, lam,
                detector=SMF(optimal_focus=True, defocus_m=dz))
            side = tails[(D, "bucket")].side
            modes[D] = ApertureModes(
                SLOPE_MIN_MODES, side,
                circle(side, D / patch.pixel_m, self.obscuration_ratio))
        tails["mmf"] = _PostTail(
            patch, APERTURES_M[0], self.obscuration_ratio, lam,
            detector=MMF(core_radius_m=MMF_CORE_RADIUS_M,
                         focal_length_m=MMF_FOCAL_LENGTH_M,
                         numerical_aperture=MMF_NA))
        ctx["tails"] = tails
        ctx["modes"] = modes
        ctx["point"] = _Mask(tails[(APERTURES_M[0], "bucket")]._ref,
                             POINT_RADIUS_M)
        return ctx

    def __call__(self, rec):
        ctx = self._build(rec)
        tails, modes, point = ctx["tails"], ctx["modes"], ctx["point"]
        E = rec.array
        out = {
            "collected_power_stored": float(rec.scalars["collected_power"]),
            "smf_eta_stored": float(rec.scalars["smf_eta"]),
        }
        tilt_removed = {}
        for D in APERTURES_M:
            tag = "10cm" if D == 0.10 else "5cm"
            bucket = tails[(D, "bucket")]
            out[f"bucket_{tag}"] = bucket.power(E)
            out[f"smf{tag[:-2]}_untracked"] = tails[(D, "untracked")].eta(E)
            out[f"smf{tag[:-2]}_tracked"] = tails[(D, "tracked")].eta(E)
            # The tilt from the wrapped-gradient slopes of the CLIPPED field,
            # the terrestrial sensing rule of the runner. The fit takes
            # SLOPE_MIN_MODES modes and keeps the tilt pair only, as the
            # runner does for a TipTilt stage.
            m = modes[D]
            sx, sy = wrapped_gradient(bucket.clip(E), mask=m.mask)
            coeffs = m.estimate_from_slopes(sx, sy)
            out[f"a2_{tag}"] = float(coeffs[1])
            out[f"a3_{tag}"] = float(coeffs[2])
            if D == APERTURES_M[0]:
                out["max_step_10cm"] = max_abs_step(sx, sy)
            keep = np.zeros_like(coeffs)
            keep[1:3] = coeffs[1:3]
            E_tt = m.apply(E, keep, sign=-1)
            tilt_removed[D] = E_tt
            out[f"smf{tag[:-2]}_tiptilt"] = tails[(D, "untracked")].eta(E_tt)
        out["smf10_tiptilt_tracked"] = tails[(APERTURES_M[0], "tracked")].eta(
            tilt_removed[APERTURES_M[0]])
        I = np.abs(E) ** 2
        out["point_irradiance"] = float(I[point.point].mean())
        out["centre_irradiance"] = float(I[point.centre].mean())
        out["mmf10"] = float(tails["mmf"].eta(E))
        return np.array([out[c] for c in COLUMNS], dtype=np.float64)


# ---------------------------------------------------------------------------
# The cells
# ---------------------------------------------------------------------------

def default_cells():
    """Give the (path_m, cn2, preset, suffix) of the twelve full campaigns.

    The four `standard` cells past 2 km ran under the two memory-cut opt-ins
    (`fft_backend="scipy"`, `screen_generator="olb-lean"`), so their roots
    carry the `_scipy_lean` suffix. See the terrestrial_campaigns README.
    """
    cells = []
    for path_m in rc.PATHS_M:
        for cn2 in rc.CN2:
            for preset in rc.PRESET_NAMES:
                lean = preset == "standard" and path_m > 2e3
                cells.append((path_m, cn2, preset, lean))
    return cells


def open_cell(path_m, cn2, preset, lean):
    """Reopen one stored campaign, and give it with its cell record."""
    fft_backend = "scipy" if lean else rc.FFT_BACKEND
    generator = "olb-lean" if lean else rc.SCREEN_GENERATOR
    scn, geom = rc.build_scenario(path_m, cn2, "collimated")
    root = os.path.join(rc.campaigns_root(),
                        rc.cell_tag(path_m, cn2, preset, "collimated", False,
                                    rc.settings_suffix(fft_backend,
                                                       generator)))
    camp = Campaign(scn, geom, root, seed=rc.SEED, preset=preset,
                    block_size=50, patch_radius_m=rc.PATCH_RADIUS_M,
                    L0_m=rc.L0_M, precision=rc.PRECISION,
                    fft_backend=fft_backend, screen_generator=generator)
    cell_path = os.path.join(root, "cell.json")
    with open(cell_path) as fh:
        cell = json.load(fh)
    return camp, cell, scn


def focus_shifts(path_m, cn2):
    """Give the received-curvature focus shift of each aperture, in m."""
    out = {}
    for D in APERTURES_M:
        scn, _ = rc.build_scenario(path_m, cn2, "collimated", rx_aperture_m=D)
        out[D] = float(curvature_focus_shift(scn))
    return out


def extract_cell(spec, n_trials, workers, say):
    """Extract one cell, and write its table."""
    path_m, cn2, preset, lean = spec
    camp, cell, scn = open_cell(*spec)
    tag = rc.cell_tag(path_m, cn2, preset, "collimated", False, "")
    dz = focus_shifts(path_m, cn2)
    n = camp.n_stored if n_trials is None else min(n_trials, camp.n_stored)
    say(f"{tag}: {camp.n_stored} stored trials, {n} to read, "
        f"patch {camp.patch.indices.size} px at {camp.patch.pixel_m * 1e3:.2f} mm, "
        f"dz_curv 10 cm {dz[0.10] * 1e6:.1f} um, 5 cm {dz[0.05] * 1e6:.1f} um")
    t0 = time.perf_counter()
    values = camp.map_trials(ExtractTrials(dz), n_trials=n, workers=workers,
                             screen_phase=False)
    wall = time.perf_counter() - t0
    values = np.asarray(values, dtype=np.float64)
    # The cross-check of the post-hoc read against the run. The order of a
    # crop sum differs from the full-grid sum, so a value moves at the float
    # rounding level only.
    k = {c: i for i, c in enumerate(COLUMNS)}
    d_power = np.abs(values[:, k["bucket_10cm"]]
                     / values[:, k["collected_power_stored"]] - 1.0).max()
    d_eta = np.abs(values[:, k["smf10_untracked"]]
                   / values[:, k["smf_eta_stored"]] - 1.0).max()
    steps = values[:, k["max_step_10cm"]]
    say(f"{tag}: {n} trials in {wall:.1f} s ({wall / n * 1e3:.0f} ms/trial); "
        f"post-hoc vs stored: power {d_power:.2e}, eta {d_eta:.2e}; "
        f"slope step > 2.8 rad in {(steps > 2.8).mean() * 100:.1f} % of trials")
    meta = {
        "cell": cell,
        "tag": tag,
        "root": camp.root_dir,
        "n_trials": int(n),
        "wall_s": float(wall),
        "columns": list(COLUMNS),
        "apertures_m": list(APERTURES_M),
        "dz_curv_m": {f"{D:.2f}": dz[D] for D in APERTURES_M},
        "point_radius_m": POINT_RADIUS_M,
        "mmf": {"core_radius_m": MMF_CORE_RADIUS_M, "na": MMF_NA,
                "focal_length_m": MMF_FOCAL_LENGTH_M},
        "slope_min_modes": int(SLOPE_MIN_MODES),
        "crosscheck": {"power_max_rel": float(d_power),
                       "eta_max_rel": float(d_eta)},
        "scenario": repr(scn),
    }
    os.makedirs(OUT_DIR, exist_ok=True)
    np.savez_compressed(os.path.join(OUT_DIR, f"trials_{tag}.npz"),
                        values=values, columns=np.array(COLUMNS),
                        meta=json.dumps(meta))
    with open(os.path.join(OUT_DIR, f"trials_{tag}.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    return values


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cells", nargs="*", default=None,
                    help="cell tokens path:cn2:preset, e.g. 2km:3e-15:rapid "
                         "(default: all twelve)")
    ap.add_argument("--trials", type=int, default=None,
                    help="read only the first N trials (a smoke test)")
    ap.add_argument("--workers", default=None,
                    help="None, an int, or 'auto' (see Campaign.map_trials)")
    ap.add_argument("--log", default=os.path.join(HERE, "extract_trials.log"))
    args = ap.parse_args(argv)

    workers = args.workers
    if workers is not None and workers != "auto":
        workers = int(workers)
    cells = default_cells()
    if args.cells:
        wanted = set(args.cells)
        cells = [c for c in cells
                 if rc.cell_token(c[0], c[1], c[2]) in wanted]
        if not cells:
            raise SystemExit(f"no cell matches {sorted(wanted)}")

    log = open(args.log, "a")

    def say(text):
        line = f"[{time.strftime('%H:%M:%S')}] {text}"
        print(line, flush=True)
        log.write(line + "\n")
        log.flush()

    say(f"extract_trials: {len(cells)} cells, workers={workers}, "
        f"trials={args.trials}, python {sys.version.split()[0]}")
    t0 = time.perf_counter()
    for spec in cells:
        extract_cell(spec, args.trials, workers, say)
    say(f"ALL DONE in {(time.perf_counter() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()

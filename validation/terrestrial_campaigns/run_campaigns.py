"""Store a backbone set of terrestrial fidelity-2 campaigns, and report the cost.

WHAT IT DOES. It runs six terrestrial (horizontal) fidelity-2 cases through the
production `Campaign` store (olb.waveoptics.turbulence.campaign): three path
lengths (2, 5 and 10 km) crossed with two turbulence strengths (Cn2 = 3e-15 and
1e-14 m^-2/3), each one at the `rapid` and at the `standard` preset. Each trial
stores the collected power of the 10 cm receive aperture, the single-mode-fibre
coupling efficiency of that aperture, and the complex receive-plane field patch
of radius 5 cm. So a smaller aperture, another detector, another focal length
and another defocus are a POST-HOC `recouple` or `recollect`, with no new
propagation.

WHAT IT IS NOT. This module is a RUNNER. It stores the campaigns and it reports
what they cost. It does NO distribution analysis: the lognormal against
gamma-gamma question, the index split, the fade quantiles, the rapid-against-
standard verdict and the comparison with the analytic terrestrial Terms are all
DEFERRED (owner decision). See the README.

WHY THE DATASET. The stored trials serve four open questions:
- backlog 1-8, gate (b): is the terrestrial fade lognormal over this band?
- the rapid preset: is it usable for a terrestrial link?
- backlog 2-N2: how does a receive aperture that holds much of the beam behave?
- the single-mode-fibre coupling distribution of a horizontal link.

NO TIP-TILT CORRECTION. Fidelity 2 models no adaptive optics and no tracking
(backlog 2-AO), so every fibre number here is UNTRACKED and it holds the full
beam wander.

Sources:
- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.
  (2005), DOI 10.1117/3.626196. Ch. 4, Eqs. (7) and (8), printed p. 87 (the
  free-space Gaussian beam radius and the irradiance profile); Ch. 6 (the
  coherence radius rho_0 and the Fried parameter r_0 = 2.1 rho_0); Ch. 8,
  Eq. (20), printed p. 264 (the plane-wave Rytov variance).
- Shaklan and Roddier, Appl. Opt. 27 (1988) 2334, DOI 10.1364/AO.27.002334.
  The single-mode-fibre coupling parameter a, which `optimal_focus` sets to
  1.12.
- S. A. Self, Appl. Opt. 22 (1983) 658, DOI 10.1364/AO.22.000658. The focus
  shift of a diverging Gaussian input, which `curvature_focus_shift` gives.
- Schmidt, Numerical Simulation of Optical Wave Propagation, (2010),
  DOI 10.1117/3.866274, Ch. 9. The split-step method of the trials.

Run it from the repository root:

    python -m validation.terrestrial_campaigns.run_campaigns --dry-run
    python -m validation.terrestrial_campaigns.run_campaigns --smoke --workers 8
    python -m validation.terrestrial_campaigns.run_campaigns --workers 8 --block-size 50
    python -m validation.terrestrial_campaigns.run_campaigns --cells 2km:3e-15:rapid
    python -m validation.terrestrial_campaigns.run_campaigns --launch diverged
    python -m validation.terrestrial_campaigns.run_campaigns \
        --fft-backend scipy --screen-generator olb-lean --workers auto

THE TWO SPEED OPT-INS. `--fft-backend scipy` and `--screen-generator olb-lean`
make a trial faster (validation/memory_cut/). They agree with the default
settings at the rounding level of single precision, about 6e-7 in the
collected power and in the coupling efficiency, but each one ENTERS the
campaign fingerprint. So a campaign that takes them gets its OWN root, through
the suffix `_scipy`, `_lean` or `_scipy_lean`, and it never mixes with a store
that the default settings made.
"""

import argparse
import csv
import json
import math
import os
import subprocess
import sys
import threading
import time
import warnings

import numpy as np

from olb.beam import gaussz, virtual_waist
from olb.geometry import HorizontalPath
from olb.models.coupling.terrestrial import (_smf_optics,
                                             curvature_focus_shift)
from olb.scenario import TerrestrialChannel, TerrestrialScenario
from olb.terminal import SMF, Terminal, Transmitter
from olb.turbulence.andrews.scintillation import rytov_variance
from olb.turbulence.andrews.structure import coherence_radius, fried_parameter
from olb.waveoptics.turbulence import Campaign
from olb.waveoptics.turbulence.sampling import turbulent_grid

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# The cell table. Three paths x two Cn2 values x two presets = twelve roots,
# of which six run at one launch setting.
# ---------------------------------------------------------------------------

PATHS_M = (2e3, 5e3, 10e3)
CN2 = (3e-15, 1e-14)
PRESET_NAMES = ("rapid", "standard")

LAM = 1550e-9
WAIST_M = 5e-3
TX_APERTURE_M = 0.10
RX_APERTURE_M = 0.10

# The receive diameters of the study. Only the LARGEST one enters the
# propagation: the smaller one is a post-hoc crop of the stored field patch.
APERTURES_M = (0.05, 0.10)

# The stored field disc. It equals the clip radius of the 10 cm aperture.
PATCH_RADIUS_M = 0.05

# The outer scale of the screens. The owner fixed it at 25 m on 2026-09-05
# (backlog 2-P5): L0 = inf is grid-dependent and it is pessimistic on the
# fibre fade tail.
L0_M = 25.0

SEED = 20260906
PRECISION = "single"

# The DEFAULT settings of the two speed opt-ins. "numpy" and "olb" are the
# settings of record: the eight finished cells hold them. The opt-ins
# "scipy" and "olb-lean" agree with them at the rounding level of single
# precision (about 6e-7 in the collected power and in the coupling
# efficiency, validation/memory_cut/), but they ENTER the campaign
# fingerprint, so they need their own roots.
FFT_BACKEND = "numpy"
SCREEN_GENERATOR = "olb"

# The optional diverged launch. The diffraction divergence of a 5 mm waist is
# about 1.0e-4 rad, so this is a beam opened by about a factor of two. Same
# value as validation/lognormal_certification.
DIVERGENCE_RAD = 2.0e-4

# The trial count that the smoke projection reports against.
PROJECT_TRIALS = 2000

# The default trial and block counts of a smoke run and of a full run.
SMOKE_TRIALS, SMOKE_BLOCK = 16, 2
FULL_TRIALS, FULL_BLOCK = 2000, 50

ENV_ROOT = "OLB_TERRESTRIAL_CAMPAIGNS_ROOT"


# ---------------------------------------------------------------------------
# The cells
# ---------------------------------------------------------------------------

def _path_tag(path_m):
    """Give the short path label of a path length, for example "2km"."""
    return f"{path_m / 1e3:g}km"


def cell_token(path_m, cn2, preset):
    """Give the command-line token of one cell, for example 2km:3e-15:rapid."""
    return f"{_path_tag(path_m)}:{cn2:.0e}:{preset}"


def all_cells():
    """Give the (path_m, cn2, preset) triples of the whole table."""
    return [(p, c, q) for p in PATHS_M for c in CN2 for q in PRESET_NAMES]


def parse_cells(tokens):
    """Turn command-line cell tokens into (path_m, cn2, preset) triples.

    Args:
        tokens: a list of tokens such as "2km:3e-15:rapid", or None for the
                whole table.

    Returns:
        A list of (path_m, cn2, preset) triples, in the given order.

    Raises:
        ValueError: a token does not name a cell of the table.
    """
    table = all_cells()
    if not tokens:
        return table
    known = {cell_token(*c): c for c in table}
    out = []
    for token in tokens:
        if token not in known:
            raise ValueError(
                f"unknown cell {token!r}. The cells are: "
                + ", ".join(sorted(known)))
        out.append(known[token])
    return out


def settings_suffix(fft_backend=FFT_BACKEND, screen_generator=SCREEN_GENERATOR):
    """Give the root suffix of the two speed opt-ins.

    Each opt-in ENTERS the campaign fingerprint, so a campaign that takes one
    of them cannot reopen a store that a default run made. The suffix keeps
    the two stores apart: "_scipy", "_lean", or "_scipy_lean". The defaults
    give an empty suffix, so every old root keeps its name.

    Args:
        fft_backend:      "numpy" (the default) or "scipy".
        screen_generator: "olb" (the default) or "olb-lean".

    Returns:
        The suffix string.
    """
    return (("_scipy" if fft_backend != FFT_BACKEND else "")
            + ("_lean" if screen_generator != SCREEN_GENERATOR else ""))


def cell_tag(path_m, cn2, preset, launch, smoke, suffix=""):
    """Give the directory name of one campaign."""
    return (f"L{_path_tag(path_m)}_cn2{cn2:.0e}_{preset}"
            + ("_diverged" if launch == "diverged" else "")
            + suffix
            + ("_smoke" if smoke else ""))


def campaigns_root():
    """Give the parent directory of every campaign of this study."""
    return os.environ.get(ENV_ROOT) or os.path.join(HERE, "campaigns")


def build_scenario(path_m, cn2, launch, rx_aperture_m=RX_APERTURE_M):
    """Build one terrestrial case: a small launch, a fibre-coupled receiver.

    Args:
        path_m:        the horizontal path length, in m.
        cn2:           the constant Cn2 of the path, in m^-2/3.
        launch:        "collimated" or "diverged".
        rx_aperture_m: the receive aperture diameter, in m.

    Returns:
        The pair (TerrestrialScenario, HorizontalPath).
    """
    divergence = DIVERGENCE_RAD if launch == "diverged" else None
    near = Terminal(aperture_m=TX_APERTURE_M, wavelength_m=LAM,
                    transmitter=Transmitter(waist_m=WAIST_M,
                                            divergence_rad=divergence))
    far = Terminal(aperture_m=rx_aperture_m, wavelength_m=LAM,
                   detector=SMF(optimal_focus=True, sensitivity_dbm=-40))
    channel = TerrestrialChannel(path_length_m=path_m,
                                 attenuation_db_per_km=0.0, cn2=cn2)
    return (TerrestrialScenario(near=near, far=far, channel=channel),
            HorizontalPath(path_m))


# ---------------------------------------------------------------------------
# The cell record
# ---------------------------------------------------------------------------

def beam_radius_m(path_m, launch):
    """Give the VACUUM received beam radius w(L) of the launch, in m.

    A deliberately diverged transmitter is an ordinary Gaussian beam from a
    virtual waist w_v at the distance d behind the aperture (olb.beam), so the
    free-space radius at the range L is w(L) = gaussz(w_v, d + L). Source:
    Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 4, Eqs. (7)
    and (8), printed p. 87. Copied from
    validation/lognormal_certification/lognormal_certification.py.
    """
    divergence = DIVERGENCE_RAD if launch == "diverged" else None
    w_v, offset = virtual_waist(WAIST_M, divergence, LAM)
    return float(gaussz(w_v, offset + path_m, LAM))


def fill_fraction(diameter_m, w_L):
    """Give the vacuum power fraction that a diameter D catches, eta_fill.

    A Gaussian beam of 1/e^2 radius w carries the irradiance
    I(r) = (2/(pi w^2)) exp(-2 r^2 / w^2). The radial integral to r = D/2 is

        eta_fill = 1 - exp(-2 (D/2)^2 / w^2).

    Source of the irradiance profile: Andrews and Phillips, 2nd ed. (2005),
    DOI 10.1117/3.626196, Ch. 4, Eq. (8), printed p. 87. Copied from
    validation/lognormal_certification/lognormal_certification.py.
    """
    return float(1.0 - np.exp(-2.0 * (diameter_m / 2.0) ** 2 / w_L ** 2))


def aperture_records(path_m, cn2, launch, w_L):
    """Give the per-aperture record of one cell.

    Each entry holds the beam-fill fraction, the received-curvature focus shift
    and the single-mode-fibre optics that `optimal_focus` resolves for that
    aperture.

    Args:
        path_m: the path length, in m.
        cn2:    the Cn2 of the path, in m^-2/3.
        launch: "collimated" or "diverged".
        w_L:    the vacuum received beam radius, in m.

    Returns:
        A list of dicts, one for each diameter of APERTURES_M.
    """
    out = []
    for a in APERTURES_M:
        scn, _ = build_scenario(path_m, cn2, launch, rx_aperture_m=a)
        f, w_m = _smf_optics(scn.rx_terminal.detector, a, LAM)
        out.append({
            "aperture_m": float(a),
            "eta_fill": fill_fraction(a, w_L),
            # The true focus of the received diverging beam sits BEYOND the
            # focal plane. Set detector.defocus_m to this value for a tracked
            # (aligned) coupler. Source: S. A. Self, Appl. Opt. 22 (1983) 658,
            # DOI 10.1364/AO.22.000658.
            "curvature_focus_shift_m": float(curvature_focus_shift(scn)),
            # optimal_focus picks f = pi*(D/2)*w_m/(lambda*1.12). Source:
            # Shaklan and Roddier, DOI 10.1364/AO.27.002334; olb holds the
            # rule in olb.models.coupling._common and
            # olb.models.coupling.terrestrial._smf_optics.
            "smf_focal_length_m": float(f),
            "smf_mode_field_radius_m": float(w_m),
        })
    return out


def cell_record(path_m, cn2, preset, launch, camp, report, sizer_warnings):
    """Build the cell.json content of one campaign.

    Args:
        path_m:         the path length, in m.
        cn2:            the Cn2 of the path, in m^-2/3.
        preset:         the quality preset name.
        launch:         "collimated" or "diverged".
        camp:           the Campaign object.
        report:         the SamplingReport of a fresh sizing, or None.
        sizer_warnings: the warning texts of the campaign construction.

    Returns:
        A JSON-ready dict.
    """
    w_L = beam_radius_m(path_m, launch)
    # The plane-wave Rytov variance sigma_R^2 = 1.23 Cn2 k^(7/6) L^(11/6).
    # Source: Andrews and Phillips, DOI 10.1117/3.626196, Ch. 8, Eq. (20),
    # printed p. 264.
    sigma2_r = float(rytov_variance(LAM, path_m, cn2, wave="plane"))
    # rho_0 is the 1/e coherence radius, and r_0 = 2.1 rho_0 is the Fried
    # parameter. Source: the same book, Ch. 6, Eqs. (56) and (64), printed
    # pp. 193 and 194.
    rho0 = float(coherence_radius(LAM, path_m, cn2, wave="plane"))
    return {
        "cell": cell_token(path_m, cn2, preset),
        "path_length_m": float(path_m),
        "cn2_m_m23": float(cn2),
        "preset": preset,
        "launch": launch,
        "divergence_rad": (DIVERGENCE_RAD if launch == "diverged" else None),
        "wavelength_m": LAM,
        "waist_m": WAIST_M,
        "sigma2_R_plane": sigma2_r,
        "rho0_m": rho0,
        "r0_m": float(fried_parameter(rho0)),
        "beam_radius_m": w_L,
        "apertures": aperture_records(path_m, cn2, launch, w_L),
        "grid": {"n": int(camp.grid.n), "side_m": float(camp.grid.size_m),
                 "pixel_m": float(camp.grid.size_m / camp.grid.n),
                 "scaled": bool(camp.grid.scaled)},
        "screens": {"n": int(camp.plan.z_m.size),
                    "z_m": camp.plan.z_m.tolist(),
                    "sigma2_r": camp.plan.sigma2_r.tolist(),
                    "r0_m": camp.plan.r0_m.tolist(),
                    "r0_total_m": float(camp.plan.r0_total_m)},
        "sizer": {
            "warnings": list(sizer_warnings),
            "step_over_limit_max": (None if report is None
                                    else float(report.step_over_limit_max)),
            "pixels_per_r0": (None if report is None
                              else float(report.pixels_per_r0)),
            "grid_margin": (None if report is None
                            else float(report.grid_margin)),
            "fresnel_pixels_min": (None if report is None
                                   else float(report.fresnel_pixels_min)),
            "n_clamped": (None if report is None else bool(report.n_clamped)),
            "clamp_factor": (None if report is None
                             else float(report.clamp_factor)),
            "feature_m": (None if report is None else float(report.feature_m)),
            "feature_pixels": (None if report is None
                               else float(report.feature_pixels)),
        },
        "campaign": {
            "root": camp.root_dir,
            "seed": int(camp.seed),
            "block_size": int(camp.block_size),
            "L0_m": L0_M,
            "precision": camp.precision,
            "patch_radius_m": float(camp.patch_radius_m),
            "patch_n": int(camp.patch.n),
            "patch_pixels": int(camp.patch.indices.size),
            "patch_pixel_m": float(camp.patch.pixel_m),
            "fft_backend": camp.fft_backend,
            "screen_generator": camp.screen_generator,
            "fingerprint": camp.fingerprint,
        },
    }


def make_campaign(path_m, cn2, preset, launch, block_size, smoke,
                  fft_backend=FFT_BACKEND, screen_generator=SCREEN_GENERATOR):
    """Size one campaign, and give it with its sizing record.

    The construction only SIZES the grid and plans the screens. It runs no
    trial. A campaign that already sits on disk takes its grid and its plan
    from the manifest, so it never re-sizes.

    Args:
        path_m:     the path length, in m.
        cn2:        the Cn2 of the path, in m^-2/3.
        preset:     the quality preset name.
        launch:     "collimated" or "diverged".
        block_size: the trials in one block.
        smoke:      True gives the campaign its own `_smoke` root.
        fft_backend:      "numpy" or the "scipy" opt-in.
        screen_generator: "olb" or the "olb-lean" opt-in.

    Returns:
        The tuple (Campaign, SamplingReport or None, the warning texts).
    """
    scn, geom = build_scenario(path_m, cn2, launch)
    root = os.path.join(campaigns_root(),
                        cell_tag(path_m, cn2, preset, launch, smoke,
                                 settings_suffix(fft_backend,
                                                 screen_generator)))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        # A fresh sizing gives the SamplingReport, which the Campaign does not
        # keep. It is cheap: it plans the screens and it runs no trial.
        try:
            _, _, report = turbulent_grid(scn, geom, preset=preset, L0_m=L0_M)
        except Exception:               # pragma: no cover - a sizer failure
            report = None
        camp = Campaign(scn, geom, root, seed=SEED, preset=preset,
                        block_size=block_size,
                        patch_radius_m=PATCH_RADIUS_M, L0_m=L0_M,
                        precision=PRECISION, fft_backend=fft_backend,
                        screen_generator=screen_generator)
    texts = sorted({str(w.message) for w in caught})
    return camp, report, texts


def write_cell_json(camp, record):
    """Write cell.json next to the campaign root."""
    path = os.path.join(camp.root_dir, "cell.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(record, fh, indent=1)
    return path


# ---------------------------------------------------------------------------
# The cost model of the dry run
# ---------------------------------------------------------------------------

def worker_bytes(grid_n):
    """Give the working bytes of ONE trial on one worker.

    The split step holds the field and, under the lazy screen stack, two
    screens at a time. In single precision the field is complex64 (8 B for
    each pixel) and a screen is float32 (4 B for each pixel). So

        bytes = n^2 * (8 + 2 * 4).

    This is a FLOOR, not a measurement: the FFT plan, the boundary mask and
    python itself add more. The smoke run measures the true working set.

    Args:
        grid_n: the pixel count of one grid side.

    Returns:
        The bytes as an int.
    """
    n2 = int(grid_n) ** 2
    return n2 * 8 + 2 * n2 * 4


def dry_run_rows(specs):
    """Give the dry-run table rows of the sized campaigns.

    Args:
        specs: a list of dicts from `size_all`.

    Returns:
        A list of string lists, the header first.
    """
    head = ["cell", "preset", "sigma_R^2", "rho0 cm", "w(L) cm", "n px",
            "side m", "px mm", "screens", "MB/worker", "clamp", "px/feat",
            "warn", "root"]
    rows = [head]
    for s in specs:
        rec, camp = s["record"], s["camp"]
        rows.append([
            f"{_path_tag(s['path_m'])}:{s['cn2']:.0e}",
            s["preset"],
            f"{rec['sigma2_R_plane']:.3f}",
            f"{rec['rho0_m'] * 100:.2f}",
            f"{rec['beam_radius_m'] * 100:.1f}",
            f"{camp.grid.n:d}",
            f"{camp.grid.size_m:.3f}",
            f"{rec['grid']['pixel_m'] * 1e3:.2f}",
            f"{camp.plan.z_m.size:d}",
            f"{worker_bytes(camp.grid.n) / 2 ** 20:.1f}",
            (f"{rec['sizer']['clamp_factor']:.1f}x"
             if rec['sizer']['clamp_factor'] is not None else "-"),
            (f"{rec['sizer']['feature_pixels']:.1f}"
             if rec['sizer']['feature_pixels'] is not None else "-"),
            f"{len(rec['sizer']['warnings']):d}",
            # The directory name carries the opt-in suffix, so the dry run
            # shows WHICH store each cell writes into.
            os.path.basename(camp.root_dir),
        ])
    return rows


def dry_run_warning_lines(specs):
    """Give one line for each sizer warning, so the dry run shows the text.

    The table gives the count only. A clamped grid says WHICH sampling rules
    it breaks, and that text is the flag a reader needs before a run.
    """
    lines = []
    for s in specs:
        for w in s["record"]["sizer"]["warnings"]:
            lines.append(f"  {_path_tag(s['path_m'])}:{s['cn2']:.0e}:"
                         f"{s['preset']}: {w}")
    return lines


def print_table(rows, say):
    """Print a table with one space-padded column for each field."""
    width = [max(len(r[i]) for r in rows) for i in range(len(rows[0]))]
    for k, row in enumerate(rows):
        say("  " + "  ".join(v.rjust(width[i]) if i else v.ljust(width[i])
                             for i, v in enumerate(row)))
        if k == 0:
            say("  " + "  ".join("-" * w for w in width))


# ---------------------------------------------------------------------------
# The memory poll
# ---------------------------------------------------------------------------

def _python_stats():
    """Give (count, total_mb, max_single_mb) of every python.exe.

    It reads `tasklist`, as
    validation/campaign_resources/campaign_resources.py `_python_processes`
    does. This function adds the LARGEST single working set, which the shared
    one does not report, because a pool worker is the value that sizes the
    machine. Off Windows it gives NaN.
    """
    if not sys.platform.startswith("win"):
        return np.nan, np.nan, np.nan
    try:
        out = subprocess.run(
            ["tasklist", "/FI", "IMAGENAME eq python.exe", "/FO", "CSV", "/NH"],
            capture_output=True, text=True, timeout=10).stdout
    except Exception:
        return np.nan, np.nan, np.nan
    n, total, biggest = 0, 0.0, 0.0
    for row in csv.reader(out.splitlines()):
        if len(row) >= 5 and row[0].lower() == "python.exe":
            n += 1
            digits = "".join(ch for ch in row[4] if ch.isdigit())
            mb = float(digits or 0) / 1024
            total += mb
            biggest = max(biggest, mb)
    return n, total, biggest


class PythonPoller:
    """A background poll of the python.exe working sets.

    Use it as a context manager. It keeps the peak of the process count, of
    the summed working set and of the largest single working set.
    """

    def __init__(self, period_s=2.0):
        self.period_s = float(period_s)
        self.peak_count = 0.0
        self.peak_total_mb = 0.0
        self.peak_single_mb = 0.0
        self.n_samples = 0
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._loop, daemon=True)

    def __enter__(self):
        self._thread.start()
        return self

    def __exit__(self, *exc):
        self._stop.set()
        self._thread.join(timeout=self.period_s * 3)

    def _loop(self):
        while not self._stop.wait(self.period_s):
            n, total, biggest = _python_stats()
            if not np.isfinite(total):
                continue
            self.n_samples += 1
            self.peak_count = max(self.peak_count, n)
            self.peak_total_mb = max(self.peak_total_mb, total)
            self.peak_single_mb = max(self.peak_single_mb, biggest)

    def summary(self):
        """Give the peaks as a JSON-ready dict."""
        return {"n_samples": int(self.n_samples),
                "peak_python_count": float(self.peak_count),
                "peak_python_total_mb": float(self.peak_total_mb),
                "peak_python_single_mb": float(self.peak_single_mb)}


# ---------------------------------------------------------------------------
# The logger
# ---------------------------------------------------------------------------

def make_say(*paths):
    """Give a print function that also appends to each named log file."""
    def say(text=""):
        print(text, flush=True)
        for path in paths:
            with open(path, "a", encoding="utf-8") as fh:
                fh.write(text + "\n")
    return say


# ---------------------------------------------------------------------------
# The work
# ---------------------------------------------------------------------------

def size_all(cells, launch, block_size, smoke, order,
             fft_backend=FFT_BACKEND, screen_generator=SCREEN_GENERATOR):
    """Size every cell, and give the specs in the run order.

    Args:
        cells:      the (path_m, cn2, preset) triples.
        launch:     "collimated" or "diverged".
        block_size: the trials in one block.
        smoke:      True gives the `_smoke` roots.
        order:      "cheap-first" sorts by grid.n^2 * n_screens; "table" keeps
                    the given order.
        fft_backend:      "numpy" or the "scipy" opt-in.
        screen_generator: "olb" or the "olb-lean" opt-in.

    Returns:
        A list of dicts with the keys path_m, cn2, preset, camp, record, cost.
    """
    specs = []
    for path_m, cn2, preset in cells:
        camp, report, texts = make_campaign(path_m, cn2, preset, launch,
                                            block_size, smoke, fft_backend,
                                            screen_generator)
        record = cell_record(path_m, cn2, preset, launch, camp, report, texts)
        specs.append({"path_m": path_m, "cn2": cn2, "preset": preset,
                      "camp": camp, "record": record,
                      "cost": int(camp.grid.n) ** 2 * int(camp.plan.z_m.size)})
    if order == "cheap-first":
        specs.sort(key=lambda s: s["cost"])
    return specs


def smoke_cross_check(camp):
    """Check the post-hoc recouple against the stored in-run scalars.

    The plan expects the two to agree bit for bit, because `recouple` runs the
    SAME coupling code on the SAME stored field. Single precision can move the
    last bits, so the test takes a relative tolerance.

    Args:
        camp: the Campaign, with at least one stored block.

    Returns:
        A JSON-ready dict of the check.
    """
    res = camp.load(fields=True)
    stored = np.array([t.smf_eta for t in res.trials], dtype=float)
    back = camp.recouple(SMF(optimal_focus=True), aperture_m=RX_APERTURE_M)
    rel = float(np.max(np.abs(back / stored - 1.0)))
    small = min(APERTURES_M)
    power_small = camp.recollect(aperture_m=small)
    eta_small = camp.recouple(SMF(optimal_focus=True), aperture_m=small)
    ok = (rel < 1e-6
          and np.all(np.isfinite(power_small)) and np.all(power_small > 0.0)
          and np.all(np.isfinite(eta_small))
          and np.all((eta_small > 0.0) & (eta_small < 1.0)))
    return {"n_trials": int(stored.size),
            "max_rel_diff_smf_eta": rel,
            "mean_smf_eta_10cm": float(stored.mean()),
            "mean_smf_eta_5cm": float(eta_small.mean()),
            "min_smf_eta_5cm": float(eta_small.min()),
            "max_smf_eta_5cm": float(eta_small.max()),
            "min_collected_power_5cm": float(power_small.min()),
            "passed": bool(ok)}


def run_smoke(specs, args, say):
    """Run the smoke pass: a few trials for each cell, with the cost.

    Args:
        specs: the sized cell specs.
        args:  the parsed arguments.
        say:   the log function.

    Returns:
        A list of JSON-ready result dicts, one for each cell.
    """
    out = []
    for spec in specs:
        camp = spec["camp"]
        token = spec["record"]["cell"]
        say(f"--- {token} ---")
        say(f"  root        : {camp.root_dir}")
        say(f"  grid        : {camp.grid.n} px, {camp.grid.size_m:.3f} m, "
            f"{camp.plan.z_m.size} screens")
        say(f"  on disk     : {camp.n_stored} trials")
        for text in spec["record"]["sizer"]["warnings"]:
            say(f"  sizer warns : {text}")
        n_before = camp.n_stored
        t0 = time.time()
        with PythonPoller(period_s=args.sample_s) as poller:
            stored = camp.run(args.n_trials, workers=args.workers,
                              progress=True)
        wall = time.time() - t0
        n_new = int(stored) - int(n_before)
        # The rate counts the NEW trials only: a resumed campaign already
        # holds its old blocks.
        s_per_trial = wall / n_new if n_new > 0 else float("nan")
        hours = PROJECT_TRIALS * s_per_trial / 3600.0
        check = smoke_cross_check(camp)
        mem = poller.summary()
        say(f"  ran         : {n_new} new trials in {wall:.1f} s "
            f"({s_per_trial:.3f} s/trial)")
        say(f"  projection  : {PROJECT_TRIALS} trials at {args.workers} "
            f"workers = {hours:.2f} h")
        say(f"  peak python : {mem['peak_python_total_mb']:.0f} MB total, "
            f"{mem['peak_python_single_mb']:.0f} MB largest, "
            f"{mem['peak_python_count']:.0f} processes")
        say(f"  recouple    : max relative difference "
            f"{check['max_rel_diff_smf_eta']:.2e}, "
            f"passed = {check['passed']}")
        say()
        row = dict(spec["record"])
        row.update({"n_new_trials": n_new, "n_stored": int(stored),
                    "wall_s": wall, "s_per_trial": s_per_trial,
                    "projected_hours_2000": hours,
                    "workers": args.workers, "block_size": args.block_size,
                    "cross_check": check})
        row.update(mem)
        out.append(row)
    return out


def run_full(specs, args, say):
    """Run the full pass: the trials of each cell, with a per-cell log.

    Args:
        specs: the sized cell specs.
        args:  the parsed arguments.
        say:   the shared log function.
    """
    n_blocks = -(-args.n_trials // args.block_size)
    # "auto" sizes itself against the block count, so the warning is for an
    # explicit worker count only.
    if isinstance(args.workers, int) and n_blocks < args.workers:
        say(f"WARNING: only {n_blocks} blocks for {args.workers} workers; the "
            f"pool can use at most {n_blocks} processes. Lower --block-size.")
    for spec in specs:
        camp = spec["camp"]
        token = spec["record"]["cell"]
        tag = cell_tag(spec["path_m"], spec["cn2"], spec["preset"],
                       args.launch, False,
                       settings_suffix(args.fft_backend,
                                       args.screen_generator))
        cell_say = make_say(*args.log_paths,
                            os.path.join(HERE, f"run_{tag}.log"))
        cell_say(f"--- {token} ---")
        cell_say(f"  root        : {camp.root_dir}")
        cell_say(f"  grid        : {camp.grid.n} px, {camp.grid.size_m:.3f} m, "
                 f"{camp.plan.z_m.size} screens, L0 = {L0_M} m")
        cell_say(f"  on disk     : {camp.n_stored} trials before this call")
        for text in spec["record"]["sizer"]["warnings"]:
            cell_say(f"  sizer warns : {text}")
        n_before = camp.n_stored
        t0 = time.time()
        stored = camp.run(args.n_trials, workers=args.workers, progress=True)
        wall = time.time() - t0
        n_new = int(stored) - int(n_before)
        assert stored >= args.n_trials, (stored, args.n_trials)
        res = camp.load(args.n_trials, fields=False)
        power = np.array([t.collected_power for t in res.trials], dtype=float)
        eta = np.array([t.smf_eta for t in res.trials], dtype=float)
        # Loss is POSITIVE dB, so the mean loss is -10 log10 of the mean
        # fraction. These two lines are a SANITY check only, not an analysis.
        cell_say(f"  ran         : {n_new} new trials in {wall:.1f} s "
                 f"({wall / max(n_new, 1):.3f} s/trial), {stored} on disk")
        cell_say(f"  mean loss   : collected power "
                 f"{-10.0 * math.log10(power.mean()):.2f} dB, "
                 f"SMF coupling {-10.0 * math.log10(eta.mean()):.2f} dB")
        cell_say()


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def parse_workers(text):
    """Turn the --workers argument into an int or the string "auto".

    Args:
        text: the command-line text.

    Returns:
        An int, or the string "auto".

    Raises:
        argparse.ArgumentTypeError: the text is neither.
    """
    if text == "auto":
        return "auto"
    try:
        return int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"--workers takes an integer or 'auto', not {text!r}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cells", nargs="*", default=None,
                    help="the cells to run, as path:cn2:preset tokens, for "
                         "example 2km:3e-15:rapid. The default is all six "
                         "cells of one preset pair.")
    ap.add_argument("--n-trials", type=int, default=None,
                    help=f"the trials for each cell (default {FULL_TRIALS}, "
                         f"or {SMOKE_TRIALS} with --smoke)")
    ap.add_argument("--workers", type=parse_workers, default=8,
                    help="the process-pool size of Campaign.run: an integer, "
                         "or 'auto' to let the campaign size the pool")
    ap.add_argument("--fft-backend", choices=("numpy", "scipy"),
                    default=FFT_BACKEND,
                    help=f"the Forvard transform backend (default "
                         f"{FFT_BACKEND}). 'scipy' is a speed OPT-IN, and it "
                         f"gives the root the suffix _scipy.")
    ap.add_argument("--screen-generator", choices=("olb", "olb-lean"),
                    default=SCREEN_GENERATOR,
                    help=f"the phase-screen generator (default "
                         f"{SCREEN_GENERATOR}). 'olb-lean' is a speed OPT-IN, "
                         f"and it gives the root the suffix _lean.")
    ap.add_argument("--block-size", type=int, default=None,
                    help=f"the trials in one block (default {FULL_BLOCK}, or "
                         f"{SMOKE_BLOCK} with --smoke)")
    ap.add_argument("--launch", choices=("collimated", "diverged"),
                    default="collimated",
                    help="the transmit beam. 'diverged' opens it to "
                         f"{DIVERGENCE_RAD} rad (a spot check).")
    ap.add_argument("--order", choices=("cheap-first", "table"),
                    default="cheap-first",
                    help="cheap-first runs the small grids first")
    ap.add_argument("--dry-run", action="store_true",
                    help="size every cell, print the table, and run nothing")
    ap.add_argument("--smoke", action="store_true",
                    help="a few trials for each cell, in _smoke roots, with "
                         "the cost and the recouple cross-check")
    ap.add_argument("--sample-s", type=float, default=2.0,
                    help="the memory poll period of the smoke run, in s")
    args = ap.parse_args(argv)
    if args.n_trials is None:
        args.n_trials = SMOKE_TRIALS if args.smoke else FULL_TRIALS
    if args.block_size is None:
        args.block_size = SMOKE_BLOCK if args.smoke else FULL_BLOCK

    log_name = "smoke.log" if args.smoke else "run_campaigns.log"
    args.log_paths = () if args.dry_run else (os.path.join(HERE, log_name),)
    say = make_say(*args.log_paths)

    cells = parse_cells(args.cells)
    say(f"terrestrial fidelity-2 campaigns, {time.strftime('%Y-%m-%d %H:%M:%S')}")
    say(f"cells         : {len(cells)}  ({args.order} order)")
    say(f"launch        : {args.launch}")
    say(f"store         : {campaigns_root()}")
    say(f"fixed         : lambda {LAM * 1e9:.0f} nm, waist {WAIST_M * 1e3:.0f} mm, "
        f"rx {RX_APERTURE_M * 100:.0f} cm, patch {PATCH_RADIUS_M * 100:.0f} cm, "
        f"L0 {L0_M:.0f} m, seed {SEED}, {PRECISION} precision")
    say(f"speed opt-ins : fft {args.fft_backend}, screens "
        f"{args.screen_generator}, root suffix "
        f"{settings_suffix(args.fft_backend, args.screen_generator) or '(none)'}")
    say()

    specs = size_all(cells, args.launch, args.block_size, args.smoke,
                     args.order, args.fft_backend, args.screen_generator)
    for spec in specs:
        write_cell_json(spec["camp"], spec["record"])

    print_table(dry_run_rows(specs), say)
    warn_lines = dry_run_warning_lines(specs)
    if warn_lines:
        say("")
        say("sizer warnings (a clamped grid names the rules it breaks):")
        for line in warn_lines:
            say(line)
    say()
    if args.dry_run:
        say("dry run: nothing ran. Drop --dry-run to store the trials.")
        return

    if args.smoke:
        results = run_smoke(specs, args, say)
        path = os.path.join(HERE, "smoke_results.json")
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(results, fh, indent=2)
        say("smoke summary:")
        rows = [["cell", "preset", "n px", "screens", "s/trial",
                 "h per 2000", "peak MB", "check"]]
        for r in results:
            rows.append([r["cell"].rsplit(":", 1)[0], r["preset"],
                         f"{r['grid']['n']:d}", f"{r['screens']['n']:d}",
                         f"{r['s_per_trial']:.3f}",
                         f"{r['projected_hours_2000']:.2f}",
                         f"{r['peak_python_total_mb']:.0f}",
                         "pass" if r["cross_check"]["passed"] else "FAIL"])
        print_table(rows, say)
        say()
        say(f"outputs       : {path}")
    else:
        run_full(specs, args, say)
        say("done. Every cell holds its trials in "
            f"{campaigns_root()} and its record in cell.json.")


if __name__ == "__main__":
    main()

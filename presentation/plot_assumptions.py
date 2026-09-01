"""Slide 7: the assumption / constraint / violation guardrail.

Every olb model states the regime it is valid in. A physics function OWNS its
assumptions (an @assumes decorator, a DOI-cited Constraint), a Term factory
traces them, and Budget.check() flags a scenario that breaks one. This slide
shows FOUR real breaks, pulled LIVE from the budgets below -- nothing on the
slide is typed by hand:

  SPACE
    * downlink, FIDELITY 1: weak-turbulence Rytov theory breaks (sigma_I^2 > 1)
    * uplink,   FIDELITY 0: the extended-Marechal Strehl breaks (sigma^2 > 1 rad^2)
  TERRESTRIAL
    * FIDELITY 0: a hard-clipped launch inside the Rayleigh range breaks the
      far-field truncation efficiency
  BOTH
    * a central obscuration breaks the filled-circular-aperture averaging filter
      (Gap 8: the book gives no annular form)

Two of the four are FIDELITY-0 limits (the analytic rung), which is the point:
the cheap closed-form model tells you, in its own voice, when to climb the
fidelity ladder.

Run from the repo root:
    python -m presentation.plot_assumptions
"""
import warnings

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from presentation.common import (use_style, figpath, WAVELENGTH_M, RED, GREY,
                                  FID_COLORS)
from olb import (SpaceScenario, TerrestrialScenario, Site, Channel,
                 TerrestrialChannel, CircularOrbit, HorizontalPath, Terminal,
                 Transmitter, Aperture, downlink_budget, terrestrial_budget)
from olb.models.coupling._common import (smf_eta_max_from_a, _smf_large_residual,
                                         SMF_OPTIMAL_A, SMF_SMALL_RESIDUAL_LIMIT)
from olb.turbulence.ao import NOLL_PISTON

WL = WAVELENGTH_M
STRONG_SITE = Site(cn2_ground=5e-13)     # a bad-seeing night, to break the weak forms
ALT_M = 500e3
ELEV_DEG = 20.0

BLUE = FID_COLORS[0]      # fidelity-0 / space accent
GREEN = FID_COLORS[2]     # terrestrial accent
INK = "#1a1a1a"
FAINT = "#f4f6f9"


# ----------------------------------------------------------------------------
# Build the three broken budgets and pull the real numbers.
# ----------------------------------------------------------------------------
def _first_violation(term, needle):
    """Return the first violation on a Term that contains `needle` (else '')."""
    a = term.assumptions
    if a is None:
        return ""
    for v in a.violations:
        if needle.lower() in v.lower():
            return v
    return ""


def build_cards():
    """Run the budgets once and distil each break into a card dict."""
    cards = {}

    # --- SPACE, downlink, fidelity 1: weak-fluctuation Rytov limit (#3) ------
    gnd = Terminal(aperture_m=0.3, obscuration_ratio=0.3, wavelength_m=WL,
                   pointing_jitter_rad=2e-6, detector=Aperture(sensitivity_dbm=-45.0))
    sat = Terminal(aperture_m=0.10, wavelength_m=WL, pointing_jitter_rad=1e-6,
                   transmitter=Transmitter(waist_m=0.04, power_dbm=30.0))
    dl = SpaceScenario(ground=gnd, space=sat, direction="downlink",
                       channel=Channel(site=STRONG_SITE, altitude_m=ALT_M))
    dlb = downlink_budget(dl, CircularOrbit(ALT_M, elevation_deg=ELEV_DEG),
                          fidelity=1)
    scint = next(t for t in dlb.terms if t.category == "coupling")
    s2I = float(scint.meta["sigma2_I"])
    cards["weak"] = dict(
        accent=BLUE, fidelity=1, link="downlink",
        title="Weak-turbulence fade",
        assume="Lognormal Rytov theory: the fade is a weak perturbation.",
        limit="sigma_I^2 = sigma_R^2 < 1",
        cite="Andrews & Phillips, Ch. 8  ·  DOI 10.1117/3.626196",
        scenario=f"Cn2(0)=5e-13 m^-2/3, {ELEV_DEG:.0f}° elevation, 0.3 m bucket",
        fired=f"sigma_I^2 = {s2I:.2f}  ≥ 1  → weak Rytov theory does not hold",
        source="plane_wave_scintillation_index",
    )

    # --- SPACE, fidelity 0: fibre coupling — TWO regimes (#9) ----------------
    # Not a break: the SMF coupling Term switches equation at the residual
    # boundary, so it stays valid past sigma^2 = 1 rad^2 instead of breaking.
    # The numbers come from the real coupling kernel (coupling/_common.py).
    eta_max = float(smf_eta_max_from_a(SMF_OPTIMAL_A))          # 0.8145
    s2x = float(SMF_SMALL_RESIDUAL_LIMIT)                       # 1.0 rad^2
    r_marechal = float(np.exp(-s2x))                            # Marechal /eta_max
    r_dikmelik = float(_smf_large_residual(s2x, eta_max) / eta_max)  # D-D /eta_max
    cards["coupling"] = dict(
        accent=BLUE, fidelity=0, link="downlink SMF  (shared w/ terrestrial)",
        title="Single-mode fibre coupling",
        assume="Gaussian modal overlap; the residual phase variance sets the loss.",
        case_small=("σ² < 1", "η = η_max · e^(−σ²)"),
        case_large=("σ² ≥ 1", f"η = η_max · [1 + σ²/{NOLL_PISTON:.2f}]^(−6/5)"),
        cite="Maréchal (σ²<1)  ·  Dikmelik-Davidson 2005, DOI 10.1364/AO.44.004946",
        crossover=f"switches at σ² = 1 rad²   ·   "
                  f"{r_marechal:.2f}·η_max → {r_dikmelik:.2f}·η_max",
    )

    # --- TERRESTRIAL, fidelity 0: far-field truncation limit (#1) ------------
    near = Terminal(aperture_m=0.06, wavelength_m=WL, pointing_jitter_rad=5e-6,
                    transmitter=Transmitter(waist_m=0.05, power_dbm=30.0))
    far = Terminal(aperture_m=0.06, obscuration_ratio=0.2, wavelength_m=WL,
                   detector=Aperture(sensitivity_dbm=-40.0))
    L = 1000.0
    tscn = TerrestrialScenario(
        near=near, far=far,
        channel=TerrestrialChannel(site=Site(), path_length_m=L,
                                   attenuation_db_per_km=0.0, cn2=1e-14))
    tb = terrestrial_budget(tscn, HorizontalPath(L))
    trunc = next(t for t in tb.terms if t.name == "transmit Gaussian efficiency")
    alpha = float(trunc.meta["alpha"])
    clip_db = float(trunc.mean_db)
    zr = float(np.pi * 0.05 ** 2 / WL)
    cards["farfield"] = dict(
        accent=GREEN, fidelity=0, link="terrestrial",
        title="Far-field truncated beam",
        assume="Far-field (Fraunhofer) Gaussian into a circular aperture.",
        limit="range  >  Rayleigh range z_R",
        cite="near-field flag  ·  gaussian_efficiency_term",
        scenario=f"hard clip α={alpha:.2f} ({clip_db:.1f} dB), path {L/1e3:.1f} km",
        fired=f"range {L/1e3:.1f} km  <  z_R = {zr/1e3:.2f} km  → far-field η invalid",
        source="tx_gaussian_efficiency_term",
    )

    # --- BOTH: central obscuration (#11) -------------------------------------
    obsc_dl = _first_violation(scint, "obscuration")
    obsc_terr = _first_violation(
        next(t for t in tb.terms if t.name == "scintillation"), "obscuration")
    cards["obscuration"] = dict(
        dl_ratio=0.30, terr_ratio=0.20,
        dl_txt=obsc_dl, terr_txt=obsc_terr,
    )
    return cards


# ----------------------------------------------------------------------------
# Draw one slide.
# ----------------------------------------------------------------------------
def _round(ax, x, y, w, h, fc, ec="none", lw=0, r=0.018, z=1, alpha=1.0):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
        mutation_aspect=h / w if w else 1, fc=fc, ec=ec, lw=lw, zorder=z,
        alpha=alpha, transform=ax.transAxes, clip_on=False))


def draw_card(ax, x, y, w, h, c):
    """Draw one break card top-down with a descending cursor (no overlaps)."""
    accent = c["accent"]
    _round(ax, x, y, w, h, "white", ec="#dfe3ea", lw=1.1, z=2)
    _round(ax, x, y, 0.006, h, accent, z=3, r=0.006)          # accent spine

    pad = 0.018
    tx = x + pad + 0.006
    innerw = w - 2 * pad - 0.006
    cy = y + h - 0.030                                         # row-centre cursor

    # fidelity chip + link
    ax.text(tx, cy, f"FIDELITY {c['fidelity']}", transform=ax.transAxes,
            fontsize=8.5, weight="bold", color="white", va="center", ha="left",
            bbox=dict(boxstyle="round,pad=0.32", fc=accent, ec="none"))
    ax.text(tx + 0.094, cy, c["link"], transform=ax.transAxes, fontsize=9.5,
            color=GREY, va="center", ha="left", style="italic")

    # title
    cy -= 0.046
    ax.text(tx, cy, c["title"], transform=ax.transAxes, fontsize=13,
            weight="bold", color=INK, va="center", ha="left")

    # ASSUMES label + text
    cy -= 0.040
    ax.text(tx, cy, "ASSUMES", transform=ax.transAxes, fontsize=7.8,
            weight="bold", color=accent, va="center", ha="left")
    cy -= 0.030
    n = c["assume"].count("\n") + 1
    ax.text(tx, cy, c["assume"], transform=ax.transAxes, fontsize=9.6,
            color=INK, va="top", ha="left")
    cy -= 0.029 * n + 0.008

    # LIMIT band (the DOI-cited constraint) + citation under it
    bandh = 0.036
    _round(ax, tx, cy - bandh, innerw, bandh, FAINT, z=3, r=0.010)
    ax.text(tx + 0.012, cy - bandh / 2, c["limit"], transform=ax.transAxes,
            fontsize=10.5, family="monospace", weight="bold", color=INK,
            va="center", ha="left")
    cy -= bandh + 0.022
    ax.text(tx, cy, c["cite"], transform=ax.transAxes, fontsize=8.2,
            color=GREY, va="center", ha="left")

    # SCENARIO knob
    cy -= 0.032
    ax.text(tx, cy, "SCENARIO", transform=ax.transAxes, fontsize=7.8,
            weight="bold", color=GREY, va="center", ha="left")
    ax.text(tx + 0.058, cy, c["scenario"], transform=ax.transAxes, fontsize=9.2,
            color=GREY, va="center", ha="left")

    # FIRED flag (red band)
    fh = 0.048
    cy -= 0.022 + fh
    _round(ax, tx, cy, innerw, fh, "#fdecec", z=3, r=0.012)
    ax.text(tx + 0.012, cy + fh / 2, "✕", transform=ax.transAxes, fontsize=12,
            weight="bold", color=RED, va="center", ha="left")
    ax.text(tx + 0.034, cy + fh / 2, c["fired"], transform=ax.transAxes,
            fontsize=9.4, weight="bold", color=RED, va="center", ha="left")


def draw_coupling_card(ax, x, y, w, h, c):
    """Draw the two-regime fibre-coupling card (a handled boundary, not a break)."""
    accent = c["accent"]
    _round(ax, x, y, w, h, "white", ec="#dfe3ea", lw=1.1, z=2)
    _round(ax, x, y, 0.006, h, accent, z=3, r=0.006)

    pad = 0.018
    tx = x + pad + 0.006
    innerw = w - 2 * pad - 0.006
    cy = y + h - 0.030

    # chip + link
    ax.text(tx, cy, f"FIDELITY {c['fidelity']}", transform=ax.transAxes,
            fontsize=8.5, weight="bold", color="white", va="center", ha="left",
            bbox=dict(boxstyle="round,pad=0.32", fc=accent, ec="none"))
    ax.text(tx + 0.094, cy, c["link"], transform=ax.transAxes, fontsize=9.5,
            color=GREY, va="center", ha="left", style="italic")

    # title + assume
    cy -= 0.044
    ax.text(tx, cy, c["title"], transform=ax.transAxes, fontsize=13,
            weight="bold", color=INK, va="center", ha="left")
    cy -= 0.038
    ax.text(tx, cy, "ASSUMES", transform=ax.transAxes, fontsize=7.8,
            weight="bold", color=accent, va="center", ha="left")
    cy -= 0.027
    ax.text(tx, cy, c["assume"], transform=ax.transAxes, fontsize=9.4,
            color=INK, va="top", ha="left")
    cy -= 0.034

    # two regime bands, each: a case tag + the mono equation (one shared cite)
    for tag, eqn in (c["case_small"], c["case_large"]):
        bandh = 0.036
        _round(ax, tx, cy - bandh, innerw, bandh, FAINT, z=3, r=0.010)
        ax.text(tx + 0.012, cy - bandh / 2, tag, transform=ax.transAxes,
                fontsize=9.5, family="monospace", weight="bold", color=accent,
                va="center", ha="left")
        ax.text(tx + 0.066, cy - bandh / 2, eqn, transform=ax.transAxes,
                fontsize=10, family="monospace", weight="bold", color=INK,
                va="center", ha="left")
        cy -= bandh + 0.014
    ax.text(tx + 0.012, cy, c["cite"], transform=ax.transAxes, fontsize=8.0,
            color=GREY, va="center", ha="left")
    cy -= 0.020

    # crossover strip (handled, not a violation) -> neutral accent band
    fh = 0.044
    cy -= 0.006 + fh
    _round(ax, tx, cy, innerw, fh, "#eef4ef", z=3, r=0.012)
    ax.text(tx + 0.012, cy + fh / 2, "✓", transform=ax.transAxes, fontsize=12,
            weight="bold", color=GREEN, va="center", ha="left")
    ax.text(tx + 0.034, cy + fh / 2, c["crossover"], transform=ax.transAxes,
            fontsize=9.2, weight="bold", color="#2f6b46", va="center", ha="left")


def make_slide(cards, out="7_assumptions_guardrail.png"):
    use_style()
    fig = plt.figure(figsize=(13.33, 7.5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # ---- title band ----
    ax.text(0.038, 0.945, "The model tells you when you have left its regime",
            fontsize=20, weight="bold", color=INK, va="center", ha="left")
    ax.text(0.038, 0.900,
            "Each Term owns DOI-cited assumptions. Three scenarios break one and "
            "get flagged; a fourth crosses the boundary by switching equation — "
            "all computed live.",
            fontsize=11, color=GREY, va="center", ha="left")

    # ---- column headers ----
    ax.text(0.038, 0.836, "SPACE", fontsize=12.5, weight="bold", color=BLUE,
            va="center", ha="left")
    ax.add_line(plt.Line2D([0.095, 0.487], [0.836, 0.836], color="#dfe3ea",
                           lw=1.2, transform=ax.transAxes))
    ax.text(0.517, 0.836, "TERRESTRIAL", fontsize=12.5, weight="bold",
            color=GREEN, va="center", ha="left")
    ax.add_line(plt.Line2D([0.622, 0.962], [0.836, 0.836], color="#dfe3ea",
                           lw=1.2, transform=ax.transAxes))

    # ---- cards ----
    colw = 0.449
    ch = 0.360
    top_y, bot_y = 0.470, 0.070
    # SPACE: two stacked cards (a break, then a handled boundary)
    draw_card(ax, 0.038, top_y, colw, ch, cards["weak"])
    draw_coupling_card(ax, 0.038, bot_y, colw, ch, cards["coupling"])
    # TERRESTRIAL: one card (top), aligned with the space top card
    draw_card(ax, 0.513, top_y, colw, ch, cards["farfield"])

    # ---- shared obscuration strip (fires on BOTH), lower-right ----
    o = cards["obscuration"]
    x, y, w, h = 0.513, bot_y, colw, ch
    _round(ax, x, y, w, h, "white", ec="#dfe3ea", lw=1.1, z=2)
    _round(ax, x, y, w, 0.006, "#8a8f98", z=3, r=0.006)
    pad = 0.022
    ax.text(x + pad, y + h - 0.034, "FIRES ON BOTH LINKS", transform=ax.transAxes,
            fontsize=8.5, weight="bold", color="white", va="center", ha="left",
            bbox=dict(boxstyle="round,pad=0.32", fc="#8a8f98", ec="none"))
    ax.text(x + pad, y + h - 0.086, "Central obscuration", transform=ax.transAxes,
            fontsize=13.5, weight="bold", color=INK, va="center", ha="left")
    ax.text(x + pad, y + h - 0.132,
            "The filled circular-aperture averaging filter has no annular\n"
            "form. No obscured-aperture filter exists in the book (Gap 8).",
            transform=ax.transAxes, fontsize=9.8, color=INK, va="top", ha="left")
    # two mini flags
    for i, (lbl, ratio) in enumerate((("downlink rx", o["dl_ratio"]),
                                      ("terrestrial far", o["terr_ratio"]))):
        yy = y + 0.088 - i * 0.044
        _round(ax, x + pad, yy - 0.017, w - 2 * pad, 0.034, "#fdecec", z=3, r=0.010)
        ax.text(x + pad + 0.010, yy, "✕", transform=ax.transAxes, fontsize=10.5,
                weight="bold", color=RED, va="center", ha="left")
        ax.text(x + pad + 0.030, yy,
                f"{lbl}: obscuration ratio {ratio:.2f}  →  filter not modelled",
                transform=ax.transAxes, fontsize=9.4, weight="bold", color=RED,
                va="center", ha="left")

    # ---- footer ----
    ax.text(0.038, 0.032,
            "The analytic rung names its own breaking point — it flags a broken "
            "regime, or switches to the valid equation. Either way the budget "
            "stays honest about when to climb the fidelity ladder.",
            fontsize=9.6, color=GREY, style="italic", va="center", ha="left")

    path = figpath(out)
    fig.savefig(path)
    plt.close(fig)
    return path


# ----------------------------------------------------------------------------
# Console report (the "dedicated script" half).
# ----------------------------------------------------------------------------
def report(cards):
    print("=" * 74)
    print("ASSUMPTION / CONSTRAINT / VIOLATION  — four live breaks")
    print("=" * 74)
    for key in ("weak", "farfield"):
        c = cards[key]
        print(f"\n[{c['link'].upper():>26}  fidelity {c['fidelity']}]  {c['title']}")
        print(f"   assumes : {c['assume'].replace(chr(10),' ')}")
        print(f"   limit   : {c['limit']}   ({c['cite']})")
        print(f"   scenario: {c['scenario']}")
        print(f"   BROKEN  : {c['fired']}")

    c = cards["coupling"]
    print(f"\n[{c['link'].upper():>26}  fidelity {c['fidelity']}]  {c['title']}")
    print(f"   assumes : {c['assume']}")
    print(f"   small σ²: {c['case_small'][1]}")
    print(f"   large σ²: {c['case_large'][1]}")
    print(f"   cite    : {c['cite']}")
    print(f"   HANDLED : {c['crossover']}")

    o = cards["obscuration"]
    print("\n[            BOTH LINKS]  Central obscuration (Gap 8)")
    print(f"   downlink   : {o['dl_txt'][:120]}")
    print(f"   terrestrial: {o['terr_txt'][:120]}")
    print("=" * 74)


def main():
    import sys
    try:
        sys.stdout.reconfigure(encoding="utf-8")   # the report prints unicode
    except Exception:
        pass
    warnings.simplefilter("ignore")           # check() warns; we render instead
    cards = build_cards()
    report(cards)
    path = make_slide(cards)
    print(f"\nsaved slide -> {path}")


if __name__ == "__main__":
    main()

"""The turbulent grid sizer, and the screen-placement planner.

The vacuum sizer olb.waveoptics.grid.GridSpec.for_scenario is not sufficient
for a turbulent path. Turbulence does two things to the sampling. It spreads
the beam, so the grid needs a wider side. It also adds coherence structure at
the Fried scale r0, so the grid needs a finer pixel. This module gives the
grid, the screen positions, and an honest report of what the grid achieves.

The module WARNS. It does not raise on a sampling problem, because an honest
warning is better than a silent bad answer. It follows the pattern of
GridSpec.for_scenario.

THE PLAN IS A LIST OF SLABS. Each screen carries the integrated Cn2 of one
slab of the path. The screen sits at the centre of that slab (terrestrial), or
at the Cn2-weighted centre of the layers that it holds (space).

THE BOUNDARY MASK IS ALWAYS ON. The subharmonic content of a screen is not
periodic on the grid. The Forvard propagator IS periodic. So the split step
needs the absorbing mask of olb.waveoptics.turbulence.splitstep. The mask is
exactly 1.0 inside the radius (1 - boundary_width_frac) of the half-side. This
sizer keeps the receive aperture and the scattered light inside that untouched
interior.

Sources:
- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.,
  DOI 10.1117/3.626196. Ch. 7, Eq. (57): the long-term beam radius W_LT.
  Ch. 8, Eq. (20): the plane-wave Rytov variance. Ch. 12, Eqs. (14), (36) and
  (38): the slant secant and the path weight of the downlink index.
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB, DOI 10.1117/3.866274. The scattering cone is Ch. 9, Eqs. (9.84) and
  (9.85), printed p. 173, at c = 2. The pixel rules are Sec. 9.4, printed
  p. 172. The per-screen cap is Listing 9.5, lines 37 and 38, printed p. 175.
  The absorbing boundary is Ch. 8, Eq. (8.1), printed p. 134. The range limit
  z_max = N dx^2 / lambda is Ch. 7, Eq. (7.59), printed p. 127, and Ch. 9,
  Eq. (9.89), printed p. 174. NOTE: this sizer evaluates NONE of the three
  turbulent geometry constraints, Eqs. (9.86) to (9.88), printed pp. 173 and
  174, and it matches no moment of the layer rule, Eq. (9.65), printed p. 164.
  See docs/schmidt-crosscheck.md, gaps S-21 and S-22.
- Martin and Flatte, Intensity images and statistics from numerical simulation
  of wave propagation in 3-D random media, DOI 10.1364/AO.27.002111. The
  pixel-per-coherence-length rule of a split-step simulation.
- Fried, DOI 10.1364/JOSA.56.001372. The Fried parameter r0.
"""

import warnings
from dataclasses import dataclass

import numpy as np

from ...turbulence.andrews.beam import (beam_params, effective_beam_params,
                                        wavenumber)
from ...turbulence.andrews.paths import sec_zeta
from ...turbulence.andrews.scintillation import rytov_variance
from ...turbulence.profiles import DEFAULT_HS, default_cn2_profile
from ..grid import N_MIN, GridSpec, _features, beam_magnification, forvard_max_z
from .screens import screen_r0

# The number of pixels across the smallest hard edge. The vacuum sizer asks for
# 16. A turbulent grid is much wider than a vacuum grid, so 16 pixels across a
# small aperture makes an impossible pixel count. 8 is the practical value.
# Pass a manual GridSpec to the runner for a finer edge.
PIXELS_PER_FEATURE = 8

# The largest screen count that the planner builds. A path that asks for more
# gets the cap and a warning.
MAX_SCREENS = 500


@dataclass(frozen=True)
class QualityPreset:
    """One named set of sampling rules.

    Attributes:
        name:                 the name of the preset.
        pixels_per_r0:        the pixel pitch obeys dx <= r0_total /
                              pixels_per_r0. See Martin and Flatte,
                              DOI 10.1364/AO.27.002111. Schmidt,
                              DOI 10.1117/3.866274, Sec. 9.4, printed p. 172,
                              gives the same rule in prose, from Johnston and
                              Lane: the phase step between two adjacent
                              samples stays below pi more than 99.7% of the
                              time. With Ch. 9, Eq. (9.44), printed p. 160,
                              that reads dx <= 0.332 r0, which is 3.01 pixels
                              per r0. The standard preset value of 3 lands on
                              it.
        guard:                the ratio of the grid half-side to the beam
                              radius. It has the same meaning as the guard of
                              GridSpec.for_scenario.
        n_max:                the largest pixel count.
        sigma2_r_screen_max:  the largest plane-wave Rytov contribution of one
                              screen. A screen that is stronger than this value
                              breaks the thin-screen approximation. The book
                              rule is rmax = 0.1 of Schmidt,
                              DOI 10.1117/3.866274, Listing 9.5, lines 37 and
                              38, printed p. 175, which the book credits to
                              Martin and Flatte. THE TWO NUMBERS ARE NOT THE
                              SAME QUANTITY. The book caps the LOG-AMPLITUDE
                              variance sigma_chi^2 of Ch. 9, Eqs. (9.64) and
                              (9.74), printed pp. 163 and 165. This field caps
                              the PLANE-WAVE RYTOV variance sigma_R^2, and
                              sigma_R^2 = 4 sigma_chi^2 (the self-check of
                              olb.waveoptics.schmidt.turbulence measures
                              3.9994). So rmax = 0.1 is a cap of 0.4 on this
                              field, and 0.05 / 0.10 / 0.25 are 8x / 4x / 1.6x
                              STRICTER than the book. olb is conservative here,
                              and it is not wrong.
        min_screens:          the smallest screen count. A weak path passes
                              sigma2_r_screen_max with one screen, but one
                              screen gives phase only and no scintillation, so
                              a floor is needed. TO REVISE: the integers
                              15/9/5 have NO derivation and NO DOI, unlike every
                              other preset field. Schmidt gives NO screen-count
                              floor either: Eq. (9.90), printed p. 174, is a
                              sampling floor only, and the 11 planes of
                              Sec. 9.5.2, printed p. 177, come with no formula.
                              The principled replacement is the layer moment
                              rule, Eq. (9.65), printed p. 164, which fixes the
                              screen positions and strengths together and gives
                              a real floor of 4. See the WP7 gate verdict in
                              docs/schmidt-crosscheck.md, and the _merge_layers
                              fallback note. (CLAUDE.md open items.)
        fresnel_weight_min:   the Rytov share above which a screen must obey
                              the Fresnel-scale pixel rule. A screen that
                              carries less than this share of the total Rytov
                              variance is exempt. THE EXEMPTION IS AN olb RULE.
                              Schmidt, DOI 10.1117/3.866274, Sec. 9.4, printed
                              p. 172, applies the rule to every step. The
                              exemption saves real time, so it stays. See
                              docs/schmidt-crosscheck.md, gap S-25.
        boundary_width_frac:  the width of the absorbing band, as a fraction of
                              the half-side. It goes to
                              splitstep.super_gaussian_boundary. The book
                              parameterises the boundary differently: Schmidt,
                              DOI 10.1117/3.866274, Listing 9.7, line 19,
                              printed p. 179, gives one half-width of 0.47 N
                              pixels from the centre. The two forms do not
                              convert. See gap S-15.
    """

    name: str
    pixels_per_r0: float
    guard: float
    n_max: int
    sigma2_r_screen_max: float
    min_screens: int
    fresnel_weight_min: float
    boundary_width_frac: float


PRESETS = {
    "reference": QualityPreset("reference", 4, 4, 4096, 0.05, 15, 0.005, 0.125),
    "standard": QualityPreset("standard", 3, 3, 2048, 0.10, 9, 0.02, 0.125),
    "rapid": QualityPreset("rapid", 2, 2, 1024, 0.25, 5, 0.05, 0.10),
}


@dataclass(frozen=True)
class ScreenPlan:
    """Where the screens sit, and what each screen carries.

    Attributes:
        z_m:            the distance of each screen from the INPUT plane, in m.
                        The values go up.
        cn2_int_m13:    the integrated Cn2 of each screen, in m^(1/3). The
                        value carries the slant factor already.
        r0_m:           the Fried parameter of each screen, in m.
        sigma2_r:       the plane-wave Rytov contribution of each screen.
        z_total_m:      the length of the gridded path, in m.
        r0_total_m:     the composite Fried parameter of the whole path, in m.
        direction:      "terrestrial" or "down". The space plan always
                        propagates DOWN the atmosphere. An uplink reads the
                        result through reciprocity. See
                        olb.waveoptics.turbulence.run.
    """

    z_m: np.ndarray
    cn2_int_m13: np.ndarray
    r0_m: np.ndarray
    sigma2_r: np.ndarray
    z_total_m: float
    r0_total_m: float
    direction: str


@dataclass(frozen=True)
class SamplingReport:
    """What the grid ACHIEVES, against what the preset asks for.

    Attributes:
        pixels_per_r0:       the achieved r0_total / dx.
        grid_margin:         the untouched interior half-side, divided by the
                             beam radius plus the scattering cone. A value of
                             1.0 means the light just fits. The preset guard is
                             the target.
        fresnel_pixels_min:  the smallest achieved pixel count across a
                             REQUIRED Fresnel scale sqrt(lambda z). It is
                             infinite when no screen passes
                             fresnel_weight_min. A value of 2 or more is good.
        step_over_limit_max: the largest planned gap between two screens,
                             divided by forvard_max_z. A value of 1.0 or less
                             is good. The split-step engine cuts a longer gap
                             into sub-steps. The limit is the step cap of
                             Schmidt, DOI 10.1117/3.866274, Ch. 8, Eq. (8.24),
                             printed p. 144, repeated as Ch. 9, Eq. (9.89),
                             printed p. 174. THE ROUTE DIFFERS: the book SETS
                             the plane count from that cap, Eq. (9.90). This
                             planner sets the count from the Cn2 profile and
                             only REPORTS the ratio. See gap S-10.
        sigma2_r_screen_max: the largest per-screen Rytov contribution that the
                             plan holds. It is a plane-wave Rytov variance, so
                             it is 4 times the book quantity. See
                             QualityPreset.sigma2_r_screen_max.
        n_clamped:           True means the pixel count hit n_max.
        warnings:            a tuple of the warning texts that the sizer sent.
    """

    pixels_per_r0: float
    grid_margin: float
    fresnel_pixels_min: float
    step_over_limit_max: float
    sigma2_r_screen_max: float
    n_clamped: bool
    warnings: tuple


def _screen_rytov(k, cn2_int, z_to_rx):
    """Give the plane-wave Rytov contribution of one slab.

        d(sigma_R^2) = 2.25 k^(7/6) (INT Cn2 dz) (z_to_rx)^(5/6)

    The plane-wave Rytov variance is the path integral
    sigma_R^2 = 2.25 k^(7/6) INT Cn2(z) (L - z)^(5/6) dz. A constant Cn2 gives
    the familiar 1.23 Cn2 k^(7/6) L^(11/6), because 2.25 * (6/11) = 1.23. See
    Andrews and Phillips, DOI 10.1117/3.626196, Ch. 8, Eq. (20), and the same
    constant in Ch. 12, Eqs. (36) and (38).

    THE BOOK USES THE OTHER VARIANCE. Schmidt, DOI 10.1117/3.866274, Ch. 9,
    Eqs. (9.63) and (9.73), printed pp. 163 and 165, give the same path weight
    with the constant 0.563, which is the LOG-AMPLITUDE variance
    sigma_chi^2. The ratio is 2.25 / 0.563 = 3.997, so the value that this
    function returns is sigma_R^2 = 4 sigma_chi^2. Do not compare it directly
    with the book cap rmax = 0.1 of Listing 9.5, printed p. 175.

    Args:
        k:        the wavenumber, in rad/m.
        cn2_int:  the integrated Cn2 of the slab, in m^(1/3).
        z_to_rx:  the distance from the slab to the RECEIVER, in m.

    Returns:
        The Rytov contribution. The type follows the input.
    """
    return (2.25 * k ** (7.0 / 6.0) * np.asarray(cn2_int, dtype=float)
            * np.asarray(z_to_rx, dtype=float) ** (5.0 / 6.0))


def _composite_r0(r0_m):
    """Add the screen Fried parameters: r0 = (SUM r0_i^(-5/3))^(-3/5).

    The phase variances of independent screens add, and the variance goes as
    r0^(-5/3). See Fried, DOI 10.1364/JOSA.56.001372, and Andrews and Phillips,
    DOI 10.1117/3.626196, Ch. 12, Eq. (23). The same sum is Schmidt,
    DOI 10.1117/3.866274, Ch. 9, Eq. (9.71), printed p. 165. It is the
    PLANE-wave composite; the spherical one is Eq. (9.72) on the same page.
    """
    return float(np.sum(np.asarray(r0_m, dtype=float) ** (-5.0 / 3.0))
                 ** (-3.0 / 5.0))


def _merge_layers(weights, cap, min_groups):
    """Group adjacent layers so that each group stays under the Rytov cap.

    The function only MERGES. It does not split a layer, because a Cn2 profile
    gives no sub-layer structure. A profile whose single layer is stronger than
    the cap keeps that layer, and the caller warns.

    Args:
        weights:    the Rytov contribution of each layer, in path order.
        cap:        the largest Rytov contribution of one group.
        min_groups: the smallest group count. A merge that gives fewer groups
                    than this value is refused, and each layer keeps its own
                    screen.

    Returns:
        A list of lists of layer indices.
    """
    groups, current, acc = [], [], 0.0
    for i, w in enumerate(weights):
        if current and acc + w > cap:
            groups.append(current)
            current, acc = [], 0.0
        current.append(i)
        acc += w
    if current:
        groups.append(current)
    if len(groups) < min_groups:
        # TO REVISE (CLAUDE.md open items). This fallback couples the screen
        # count to the profile sampling: it keeps ONE screen per layer, so a
        # weak space slab gets len(DEFAULT_HS) = 20 screens, and a finely
        # sampled profile would give hundreds. It should instead clamp UP to
        # EXACTLY min_groups contiguous Cn2-weighted groups, and the caller
        # should WARN only when len(weights) < min_groups (the model cannot
        # split one layer). Fix this together with the min_screens
        # justification, then re-run the three turbulent examples.
        # The book rule for a layering is Schmidt, DOI 10.1117/3.866274,
        # Ch. 9, Eq. (9.65), printed p. 164: the layered Cn2 must match the
        # continuous profile for the moments 0 <= m <= 7. That rule fixes the
        # positions and the strengths together, it decouples the screen count
        # from the profile sampling, and it gives a floor of 4 screens. This
        # merge matches moment 0 only. See gap S-22 and the WP7 gate verdict
        # in docs/schmidt-crosscheck.md.
        return [[i] for i in range(len(weights))]
    return groups


def _plan_terrestrial(scenario, geometry, preset, lam):
    """Build the screen plan and the physical extent of a horizontal path."""
    p = preset
    k = wavenumber(lam)
    tx, rx = scenario.tx_terminal, scenario.rx_terminal
    t = tx.transmitter
    if t is None:
        raise ValueError('turbulent_grid: the transmit terminal has no '
                         'Transmitter')
    # The sizer takes the LONGEST range of the geometry, exactly as
    # GridSpec.for_scenario does. The runner takes one range only.
    z_total = float(np.max(np.asarray(geometry.slant_range_m, dtype=float)))
    cn2 = float(scenario.channel.cn2)

    sigma2_total = float(rytov_variance(lam, z_total, cn2, wave='plane'))

    # The screen count. Start from the mean-share estimate, then raise it until
    # the STRONGEST screen (the one farthest from the receiver) obeys the cap.
    n_s = max(p.min_screens,
              int(np.ceil(sigma2_total / p.sigma2_r_screen_max)))
    for _ in range(20):
        z = (np.arange(n_s) + 0.5) * z_total / n_s
        cn2_int = np.full(n_s, cn2 * z_total / n_s)
        s2 = _screen_rytov(k, cn2_int, z_total - z)
        if s2.max() <= p.sigma2_r_screen_max or n_s >= MAX_SCREENS:
            break
        n_s = min(MAX_SCREENS,
                  int(np.ceil(n_s * s2.max() / p.sigma2_r_screen_max)))

    r0 = screen_r0(cn2_int, lam)
    r0_total = _composite_r0(r0)

    # The beam radius at the receiver. beam_magnification reads the Transmitter
    # divergence, so a deliberately diverged beam gives its own free-space
    # radius. The Andrews pair gives the TURBULENT broadening factor on top of
    # it. See Andrews and Phillips, DOI 10.1117/3.626196, Ch. 7, Eq. (57).
    w_free = t.waist_m * beam_magnification(scenario, z_total)
    bp = beam_params(t.waist_m, lam, z_total)
    spread = float(effective_beam_params(bp, sigma2_total).w / bp.w)
    w_lt = w_free * spread

    tx_aperture = t.aperture_m if t.aperture_m is not None else tx.aperture_m
    tx_obscuration = (t.obscuration_ratio if t.obscuration_ratio is not None
                      else tx.obscuration_ratio)
    r_beam = max(w_lt, tx_aperture / 2, rx.aperture_m / 2)
    feature = min([t.waist_m] + _features(tx_aperture, tx_obscuration)
                  + _features(rx.aperture_m, rx.obscuration_ratio))

    plan = ScreenPlan(z_m=z, cn2_int_m13=cn2_int, r0_m=r0, sigma2_r=s2,
                      z_total_m=z_total, r0_total_m=r0_total,
                      direction="terrestrial")
    return plan, r_beam, feature


def _plan_space(scenario, geometry, preset, lam, hs, cn2_profile, warns):
    """Build the screen plan and the physical extent of the atmosphere slab.

    The gridded path is the DOWNLINK slab only. The satellite is outside the
    atmosphere, so the vacuum part of the path carries no turbulence. An uplink
    reads the same slab through reciprocity. See
    olb.waveoptics.turbulence.run.
    """
    p = preset
    k = wavenumber(lam)
    ground = scenario.ground
    hs = DEFAULT_HS if hs is None else np.asarray(hs, dtype=float)
    if cn2_profile is None:
        cn2_profile = default_cn2_profile(scenario.channel.site, hs)
    cn2_profile = np.asarray(cn2_profile, dtype=float)

    # The sizer takes the LOWEST elevation of the geometry, because that is the
    # longest slant path and the worst sampling case.
    elevation = float(np.min(np.asarray(geometry.elevation_deg, dtype=float)))
    sec = float(sec_zeta(elevation))

    # The integrated Cn2 of each layer, on the slant path. This is the
    # _cn2_layers pattern of olb.models.coupling.fast, times the airmass.
    # See Andrews and Phillips, DOI 10.1117/3.626196, Ch. 12, Eq. (14).
    layer_cn2 = cn2_profile * np.gradient(hs) * sec

    h_top = float(hs[-1])
    z_total = h_top * sec
    # Order the layers from the TOP of the slab to the ground: the field enters
    # at the top. The distance of a layer to the ground receiver is h*sec.
    order = np.argsort(-hs)
    h_ord = hs[order]
    cn2_ord = layer_cn2[order]
    z_ord = (h_top - h_ord) * sec
    w_ord = _screen_rytov(k, cn2_ord, np.maximum(h_ord * sec, 1e-9))

    groups = _merge_layers(w_ord, p.sigma2_r_screen_max, p.min_screens)
    z = np.array([float(np.average(z_ord[g], weights=cn2_ord[g]))
                  if cn2_ord[g].sum() > 0 else float(np.mean(z_ord[g]))
                  for g in groups])
    cn2_int = np.array([float(cn2_ord[g].sum()) for g in groups])
    s2 = np.array([float(w_ord[g].sum()) for g in groups])
    r0 = screen_r0(cn2_int, lam)
    r0_total = _composite_r0(r0)

    if w_ord.max() > p.sigma2_r_screen_max:
        warns.append(
            f"turbulent_grid: one Cn2 LAYER carries a Rytov contribution of "
            f"{w_ord.max():.3g}, past the preset cap of "
            f"{p.sigma2_r_screen_max:.3g}. The planner only merges layers; it "
            f"does not split one. Give a finer hs grid.")

    d_ground = ground.aperture_m
    t = ground.transmitter
    if t is not None and t.aperture_m is not None:
        d_ground = max(d_ground, t.aperture_m)
    r_beam = d_ground / 2
    feature = min(_features(d_ground, ground.obscuration_ratio))

    plan = ScreenPlan(z_m=z, cn2_int_m13=cn2_int, r0_m=r0, sigma2_r=s2,
                      z_total_m=z_total, r0_total_m=r0_total, direction="down")
    return plan, r_beam, feature


def turbulent_grid(scenario, geometry, *, preset="standard", hs=None,
                   cn2_profile=None, L0_m=np.inf):
    """Size a turbulent split-step grid, and plan the screens.

    THE EXTENT RULE. The grid holds the beam AND the light that the turbulence
    scatters out of it:

        side = [guard * 2 * r_beam + 2 * (lambda / r0_total) * z] / (1 - b)

    The first part is the vacuum extent rule of GridSpec.for_scenario. The
    second part is the scattering cone: turbulence scatters light through the
    angle lambda/r0, and that light must stay off the edge of the periodic
    grid. The added term is c lambda z / r0 with c = 2, from Schmidt,
    DOI 10.1117/3.866274, Ch. 9, Eqs. (9.84) and (9.85), printed p. 173. The
    book states that c = 2 holds 97% of the light and c = 4 holds 99%, and its
    own Listing 9.6, line 2, printed p. 177, uses c = 2. THE ROUTE DIFFERS:
    the book adds the blur to the extents D1' and D2' and then feeds
    constraints 1 to 3, Eqs. (9.86) to (9.88), printed pp. 173 and 174. This
    sizer adds it to the grid SIDE and checks no constraint. See
    docs/schmidt-crosscheck.md, gap S-21. The divisor (1 - b) makes room for
    the absorbing band of the boundary mask, where b is boundary_width_frac.

    THE PIXEL RULE. The pixel obeys three limits, and the smallest wins:

        dx <= r0_total / pixels_per_r0    the coherence structure
        dx <= sqrt(lambda z_i) / 2        the Fresnel scale of screen i
        dx <= feature / (P / 2)           the hard edges, P = PIXELS_PER_FEATURE

    The first limit comes from Schmidt, DOI 10.1117/3.866274, Sec. 9.4, printed
    p. 172, and from Martin and Flatte, DOI 10.1364/AO.27.002111. The book
    states it in prose, from Johnston and Lane, and with Ch. 9, Eq. (9.44),
    printed p. 160, it reads dx <= 0.332 r0, that is 3.01 pixels per r0. The
    second limit keeps the irradiance correlation width sampled; the width is
    the Fresnel scale of the distance from the screen to the receiver (Andrews
    and Phillips, DOI 10.1117/3.626196, Ch. 8). It is the SAME rule as the
    sqrt(lambda z)/2 pitch cap of Schmidt, Sec. 9.4, printed p. 172, which the
    book also credits to Johnston and Lane. Only a screen that carries more
    than fresnel_weight_min of the total Rytov variance must obey it. A weak
    screen close to the receiver is exempt, because it adds almost no
    scintillation. That exemption is an olb rule; the book gives none.

    THE PIXEL COUNT. n is the next power of two of side/dx, as Schmidt,
    DOI 10.1117/3.866274, Listing 7.2, line 13, printed p. 128, does, inside
    the interval [256, n_max]. The clamp has no book source. A clamp does NOT
    shrink the side, because the extent is physics. The pixel grows instead,
    and the report says so.

    Args:
        scenario:    a SpaceScenario or a TerrestrialScenario.
        geometry:    an object with slant_range_m (terrestrial) or
                     elevation_deg (space). The sizer takes the worst case: the
                     longest range, or the lowest elevation.
        preset:      the name of a preset in PRESETS, or a QualityPreset.
        hs:          the height grid of the Cn2 profile, in m. None takes
                     DEFAULT_HS. Space only.
        cn2_profile: the zenith Cn2 profile on hs. None takes the site profile.
                     Space only.
        L0_m:        the outer scale, in m. The SIZER does not read it. The
                     runner passes it to phase_screen. It stays here so that
                     one call site holds all the turbulence options.

    Returns:
        A tuple (GridSpec, ScreenPlan, SamplingReport).

    Raises:
        ValueError:  the preset name is unknown, or a terrestrial transmit
                     terminal has no Transmitter.

    Warns:
        UserWarning: the grid misses a sampling rule. The report holds the same
                     texts, and it holds the ACHIEVED numbers.
    """
    if isinstance(preset, str):
        if preset not in PRESETS:
            raise ValueError(f"turbulent_grid: unknown preset {preset!r}. Use "
                             f"one of {sorted(PRESETS)}.")
        preset = PRESETS[preset]
    p = preset
    lam = scenario.tx_terminal.wavelength_m
    warns = []

    if hasattr(scenario, "direction"):
        plan, r_beam, feature = _plan_space(
            scenario, geometry, p, lam, hs, cn2_profile, warns)
    else:
        plan, r_beam, feature = _plan_terrestrial(scenario, geometry, p, lam)

    # ---- the extent ----
    scatter = (lam / plan.r0_total_m) * plan.z_total_m
    side_core = p.guard * 2 * r_beam + 2 * scatter
    side = side_core / (1.0 - p.boundary_width_frac)

    # ---- the pixel ----
    z_to_rx = plan.z_total_m - plan.z_m
    total_s2 = float(plan.sigma2_r.sum())
    share = (plan.sigma2_r / total_s2 if total_s2 > 0
             else np.zeros_like(plan.sigma2_r))
    needs_fresnel = share > p.fresnel_weight_min
    dx_wanted = min(plan.r0_total_m / p.pixels_per_r0,
                    feature / (PIXELS_PER_FEATURE / 2))
    if needs_fresnel.any():
        dx_wanted = min(dx_wanted,
                        float(np.sqrt(lam * z_to_rx[needs_fresnel]).min()) / 2)

    n_wanted = int(2 ** np.ceil(np.log2(side / dx_wanted)))
    n = int(min(max(n_wanted, N_MIN), p.n_max))
    n_clamped = n_wanted > p.n_max
    grid = GridSpec(size_m=side, n=n, scaled=False)

    # ---- the achieved numbers ----
    dx = grid.pixel_m
    achieved_r0 = plan.r0_total_m / dx
    grid_margin = ((1.0 - p.boundary_width_frac) * side / 2) / (r_beam + scatter)
    fresnel_pixels = (float((np.sqrt(lam * z_to_rx[needs_fresnel]) / dx).min())
                      if needs_fresnel.any() else float('inf'))
    gaps = np.diff(np.concatenate(([0.0], plan.z_m, [plan.z_total_m])))
    step_ratio = (float(gaps.max() / forvard_max_z(grid, lam))
                  if gaps.size else 0.0)

    if n_clamped:
        warns.append(
            f"turbulent_grid: the pixel count wants {n_wanted}, but n_max is "
            f"{p.n_max}. The grid keeps its side and takes a coarse pixel.")
    if achieved_r0 < p.pixels_per_r0:
        warns.append(
            f"turbulent_grid: the grid gives {achieved_r0:.2f} pixels per r0, "
            f"and the {p.name} preset asks for {p.pixels_per_r0}. The screens "
            f"lose the small-scale phase. Schmidt, DOI 10.1117/3.866274, "
            f"Ch. 9.")
    if fresnel_pixels < 2.0:
        warns.append(
            f"turbulent_grid: the grid gives {fresnel_pixels:.2f} pixels per "
            f"Fresnel scale. 2 or more are necessary for the scintillation. "
            f"Andrews and Phillips, DOI 10.1117/3.626196, Ch. 8.")
    if plan.sigma2_r.max() > p.sigma2_r_screen_max * (1.0 + 1e-9):
        warns.append(
            f"turbulent_grid: the strongest screen carries a Rytov "
            f"contribution of {plan.sigma2_r.max():.3g}, past the preset cap "
            f"of {p.sigma2_r_screen_max:.3g}. The thin-screen approximation is "
            f"weak here.")
    if plan.z_m.size >= MAX_SCREENS:
        warns.append(
            f"turbulent_grid: the plan hit the screen cap of {MAX_SCREENS}. "
            f"The path is very strong. Use a shorter path, or accept the "
            f"error.")

    report = SamplingReport(
        pixels_per_r0=float(achieved_r0), grid_margin=float(grid_margin),
        fresnel_pixels_min=fresnel_pixels, step_over_limit_max=step_ratio,
        sigma2_r_screen_max=float(plan.sigma2_r.max()) if plan.z_m.size else 0.0,
        n_clamped=bool(n_clamped), warnings=tuple(warns))
    for text in warns:
        warnings.warn(text)
    return grid, plan, report


if __name__ == '__main__':
    import time

    from ...geometry import CircularOrbit, HorizontalPath
    from ...scenario import (Channel, SpaceScenario, TerrestrialChannel,
                             TerrestrialScenario)
    from ...terminal import Terminal, Transmitter
    from ...turbulence.ao import plane_wave_fried_parameter_profile

    t_start = time.time()
    lam = 1550e-9

    def terrestrial(path_m, cn2, waist=0.02, tx_ap=0.10, rx_ap=0.20):
        """Build one horizontal case."""
        return TerrestrialScenario(
            near=Terminal(aperture_m=tx_ap, wavelength_m=lam,
                          transmitter=Transmitter(waist_m=waist)),
            far=Terminal(aperture_m=rx_ap, wavelength_m=lam),
            channel=TerrestrialChannel(path_length_m=path_m, cn2=cn2))

    # ---- 1. a hand-checked terrestrial case ----
    # 2 km, Cn2 = 5e-15. The numbers below repeat the rules by hand.
    L1, cn2_1 = 2000.0, 5e-15
    scn1 = terrestrial(L1, cn2_1)
    with warnings.catch_warnings(record=True) as caught1:
        warnings.simplefilter("always")
        g1, plan1, rep1 = turbulent_grid(scn1, HorizontalPath(L1))
    assert not caught1, [str(w.message) for w in caught1]
    p_std = PRESETS["standard"]

    k1 = 2 * np.pi / lam
    s2_hand = 1.23 * cn2_1 * k1 ** (7 / 6) * L1 ** (11 / 6)
    assert abs(plan1.sigma2_r.sum() / s2_hand - 1.0) < 0.02, \
        (plan1.sigma2_r.sum(), s2_hand)
    r0_hand = (0.423 * k1 ** 2 * cn2_1 * L1) ** (-3 / 5)
    assert abs(plan1.r0_total_m / r0_hand - 1.0) < 1e-9, (plan1.r0_total_m,
                                                          r0_hand)
    # The screens sit at the slab centres, and they share the path equally.
    assert plan1.z_m.size >= p_std.min_screens, plan1.z_m.size
    assert abs(plan1.cn2_int_m13.sum() - cn2_1 * L1) < 1e-18
    assert abs(plan1.z_m[0] - 0.5 * L1 / plan1.z_m.size) < 1e-9
    assert plan1.direction == "terrestrial"
    assert plan1.z_total_m == L1
    # Each screen obeys the Rytov cap, and the plan is well sampled.
    assert plan1.sigma2_r.max() <= p_std.sigma2_r_screen_max
    assert rep1.pixels_per_r0 >= p_std.pixels_per_r0, rep1.pixels_per_r0
    assert rep1.grid_margin >= 1.0, rep1.grid_margin
    assert rep1.fresnel_pixels_min >= 2.0, rep1.fresnel_pixels_min
    assert not rep1.n_clamped and rep1.warnings == ()
    assert g1.n >= N_MIN and not g1.scaled
    # The receive aperture stays inside the untouched interior of the mask.
    assert (rx_half := (1 - p_std.boundary_width_frac) * g1.size_m / 2) > 0.10, \
        rx_half

    # ---- 2. a space case: r0_total against the analytic profile ----
    hs = DEFAULT_HS
    orbit30 = CircularOrbit(altitude_m=600e3, elevation_deg=[30.0])
    space = SpaceScenario(
        ground=Terminal(aperture_m=0.50, wavelength_m=lam,
                        transmitter=Transmitter(waist_m=0.05)),
        space=Terminal(aperture_m=0.30, wavelength_m=lam),
        direction="downlink", channel=Channel(altitude_m=600e3))
    with warnings.catch_warnings(record=True) as caught2:
        warnings.simplefilter("always")
        g2, plan2, rep2 = turbulent_grid(space, orbit30)
    cn2_prof = default_cn2_profile(space.channel.site, hs)
    r0_analytic = plane_wave_fried_parameter_profile(cn2_prof, hs, lam, 30.0)
    # The plan r0 uses the Fried constant 0.423, and the Andrews chain uses
    # 0.4240. The layer sum uses np.gradient, and the analytic value uses
    # np.trapz. The two differences are small, and they partly cancel.
    assert abs(plan2.r0_total_m / r0_analytic - 1.0) < 0.01, \
        (plan2.r0_total_m, r0_analytic)
    assert plan2.direction == "down"
    assert np.all(np.diff(plan2.z_m) > 0), plan2.z_m
    assert abs(plan2.z_total_m - hs[-1] * 2.0) < 1e-6      # sec(60 deg) = 2
    # The top screen is at the top of the slab, and no screen is past the end.
    assert plan2.z_m[0] >= 0.0 and plan2.z_m[-1] <= plan2.z_total_m

    # ---- 3. the presets are monotone ----
    # A strong 3 km path, so that no preset hits the 256-pixel floor.
    scn3 = terrestrial(3000.0, 5e-14, waist=0.05, tx_ap=0.25, rx_ap=0.30)
    trio = {}
    for name in ("reference", "standard", "rapid"):
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            trio[name] = turbulent_grid(scn3, HorizontalPath(3000.0),
                                        preset=name)
    n_ref, n_std, n_rap = (trio[x][0].n for x in
                           ("reference", "standard", "rapid"))
    s_ref, s_std, s_rap = (trio[x][1].z_m.size for x in
                           ("reference", "standard", "rapid"))
    assert n_ref >= n_std >= n_rap, (n_ref, n_std, n_rap)
    assert s_ref >= s_std >= s_rap, (s_ref, s_std, s_rap)

    # ---- 4. an impossible case warns, and it does not raise ----
    # Cn2 = 1e-12 over 5 km gives r0 near 2 mm. No pixel count holds it.
    scn4 = terrestrial(5000.0, 1e-12)
    with warnings.catch_warnings(record=True) as caught4:
        warnings.simplefilter("always")
        g4, plan4, rep4 = turbulent_grid(scn4, HorizontalPath(5000.0))
    assert caught4, "an impossible case must warn"
    assert rep4.n_clamped and g4.n == PRESETS["standard"].n_max
    assert rep4.warnings, rep4
    assert any("n_max" in w for w in rep4.warnings), rep4.warnings
    assert rep4.pixels_per_r0 < PRESETS["standard"].pixels_per_r0
    # The report stays honest: the numbers are the ACHIEVED ones.
    assert abs(rep4.pixels_per_r0 - plan4.r0_total_m / g4.pixel_m) < 1e-9

    # ---- 5. a weak screen near the receiver is exempt ----
    # The lowest Cn2 layer of the space case sits 2 m from the ground receiver.
    # Its Fresnel scale is 1.8 mm, so the Fresnel rule would ask for a 0.9 mm
    # pixel on a 3 m grid. That is over 3000 pixels. The layer carries far less
    # than fresnel_weight_min of the Rytov variance, so the plan ignores it.
    z_to_rx2 = plan2.z_total_m - plan2.z_m
    share2 = plan2.sigma2_r / plan2.sigma2_r.sum()
    near_i = int(np.argmin(z_to_rx2))
    dx_if_forced = float(np.sqrt(lam * z_to_rx2[near_i]) / 2)
    n_if_forced = g2.size_m / dx_if_forced
    assert share2[near_i] < PRESETS["standard"].fresnel_weight_min, share2[near_i]
    assert n_if_forced > 4 * g2.n, (n_if_forced, g2.n)
    assert g2.pixel_m > dx_if_forced, (g2.pixel_m, dx_if_forced)
    # The screens that DO pass the threshold are still sampled.
    assert rep2.fresnel_pixels_min >= 2.0, rep2.fresnel_pixels_min

    # ---- the printed tables ----
    print("case 1, terrestrial, 2 km, Cn2 = 5e-15, standard preset:")
    print(f"  grid side               {g1.size_m:11.4f} m")
    print(f"  pixels per side         {g1.n:11d}")
    print(f"  pixel pitch             {g1.pixel_m * 1e3:11.4f} mm")
    print(f"  r0 total                {plan1.r0_total_m * 1e2:11.3f} cm")
    print(f"  sigma2_R total          {plan1.sigma2_r.sum():11.4f}")
    print(f"  screens                 {plan1.z_m.size:11d}")
    print(f"  pixels per r0           {rep1.pixels_per_r0:11.2f}")
    print(f"  grid margin             {rep1.grid_margin:11.2f}")
    print(f"  Fresnel pixels, min     {rep1.fresnel_pixels_min:11.2f}")
    print(f"  step / Forvard limit    {rep1.step_over_limit_max:11.3f}")
    print("")
    head = f"  {'i':>3}{'z [m]':>10}{'r0 [cm]':>10}{'sigma2_r':>11}{'share':>9}"
    print("  the screen table:")
    print(head)
    sh1 = plan1.sigma2_r / plan1.sigma2_r.sum()
    for i, (zz, rr, ss, hh) in enumerate(zip(plan1.z_m, plan1.r0_m,
                                             plan1.sigma2_r, sh1)):
        print(f"  {i:>3}{zz:>10.1f}{rr * 1e2:>10.2f}{ss:>11.4f}{hh:>9.4f}")
    print("")
    print("case 2, space, 30 deg elevation, 600 km, standard preset:")
    print(f"  slab length             {plan2.z_total_m * 1e-3:11.3f} km")
    print(f"  grid side               {g2.size_m:11.4f} m")
    print(f"  pixels per side         {g2.n:11d}")
    print(f"  pixel pitch             {g2.pixel_m * 1e3:11.4f} mm")
    print(f"  screens                 {plan2.z_m.size:11d}")
    print(f"  r0, plan                {plan2.r0_total_m * 1e2:11.3f} cm")
    print(f"  r0, analytic profile    {r0_analytic * 1e2:11.3f} cm")
    print(f"  ratio                   {plan2.r0_total_m / r0_analytic:11.5f}")
    print(f"  sigma2_R total          {plan2.sigma2_r.sum():11.4f}")
    print(f"  pixels per r0           {rep2.pixels_per_r0:11.2f}")
    print(f"  Fresnel pixels, min     {rep2.fresnel_pixels_min:11.2f}")
    print(f"  step / Forvard limit    {rep2.step_over_limit_max:11.3f}")
    print(f"  nearest screen share    {share2[near_i]:11.5f} "
          f"(exempt below {PRESETS['standard'].fresnel_weight_min})")
    print(f"  its forced pixel count  {n_if_forced:11.0f}")
    print("")
    print("case 3, preset monotonicity, 3 km, Cn2 = 5e-14:")
    print(f"  {'preset':<12}{'n':>8}{'screens':>10}{'pixels/r0':>12}"
          f"{'side [m]':>11}")
    for name in ("reference", "standard", "rapid"):
        gg, pp, rr = trio[name]
        print(f"  {name:<12}{gg.n:>8}{pp.z_m.size:>10}"
              f"{rr.pixels_per_r0:>12.2f}{gg.size_m:>11.3f}")
    print("")
    print("case 4, an impossible path, 5 km, Cn2 = 1e-12:")
    print(f"  r0 total                {plan4.r0_total_m * 1e3:11.4f} mm")
    print(f"  grid side               {g4.size_m:11.3f} m")
    print(f"  pixels per side         {g4.n:11d} (clamped)")
    print(f"  pixels per r0           {rep4.pixels_per_r0:11.3f}")
    for w in rep4.warnings:
        # A decimal point is a full stop too, so cut on the length.
        print(f"  warning: {w[16:88]}...")
    print("")
    print(f"(elapsed {time.time() - t_start:.1f} s)")
    print("self-check passed")

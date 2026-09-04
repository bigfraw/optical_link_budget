"""One end-to-end field propagation for an olb scenario.

The module puts the pieces together: it launches the transmit beam, it clips
the beam at the launch aperture, it propagates the beam over the range, it
clips the beam at the receive aperture, and it couples the beam into a
single-mode fibre. It gives the losses of each step in positive dB.

This is the fidelity-2 no-turbulence validator. It builds NO Term and it
touches NO budget. Compare its numbers against the analytic fidelity-0 Terms
(olb.models.geometric and olb.models.gaussian_efficiency). The two agree in the
far field with a light truncation. They disagree in the near field with a hard
truncation, because the analytic transmit efficiency is a far-field quantity.
That disagreement is the reason for this layer.

Sources:
- Goodman, Introduction to Fourier Optics, ISBN 978-0974707723. The scalar
  diffraction of a truncated beam.
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB, DOI 10.1117/3.866274. The propagator choice and the sampling limits.
- Siegman, Lasers, ISBN 978-0935702118. The Gaussian beam and the ABCD route.
"""

from dataclasses import dataclass

import numpy as np

from ..beam import virtual_waist
from ..terminal import SMF
from .field import Begin, Field, Normal, Power
from .grid import GridSpec, beam_magnification
from .lenses import Convert, Lens, LensFresnel
from .propagators import Fresnel, GForvard
from .smf import coupling_efficiency
from .sources import CircAperture, CircScreen, GaussBeam

# Below this clipped power fraction the launch aperture takes almost nothing
# from the beam. Then the field stays a pure Gaussian, and the exact ABCD route
# (GForvard) applies. Above it the field carries the aperture edge, so the
# Fresnel convolution applies. See Schmidt, DOI 10.1117/3.866274.
PURE_GAUSS_CLIP = 1e-6

# The optimal single-mode-fibre coupling parameter a = pi*(D/2)*w_m/(lambda*f).
# It gives the eta_max peak 0.8145. Source: Shaklan and Roddier, Appl. Opt. 27
# (1988) 2334, DOI 10.1364/AO.27.002334.
SMF_OPTIMAL_A = 1.12

# The default fibre mode field RADIUS, in m: SMF-28 at 1550 nm (mode field
# diameter 10.4 um). It is the same default as
# olb.models.coupling.terrestrial.SMF28_MODE_FIELD_RADIUS_M.
SMF28_MODE_FIELD_RADIUS_M = 5.2e-6


@dataclass(frozen=True)
class WaveResult:
    """The result of one field propagation.

    Attributes:
        stages:            a list of (label, Field) pairs. The labels are
                           "launch", "after tx clip", "at rx plane" and
                           "after rx clip".
        grid:              the GridSpec that the propagation used.
        tx_truncation_db:  the power that the launch aperture takes, in
                           positive dB.
        geometric_loss_db: the power that the receive aperture does not
                           collect, in positive dB.
        smf_coupling_db:   the single-mode-fibre coupling loss, in positive dB.
                           None if the receive terminal has no SMF detector.
        propagator:        the name of the propagator that ran, a string:
                           "GForvard", "Fresnel" or "LensFresnel".
    """

    stages: list
    grid: GridSpec
    tx_truncation_db: float
    geometric_loss_db: float
    smf_coupling_db: float
    propagator: str


def _loss_db(power_out, power_in):
    """Give the loss of one step in positive dB."""
    # The + 0.0 turns the negative zero of a lossless step into a plain zero.
    return float(-10 * np.log10(power_out / power_in)) + 0.0


def _normalised_gauss(Fin):
    """Scale a Gaussian field to a power of 1.0 and keep the ABCD bookkeeping.

    Normal() scales the amplitude array, but it does not touch the Gaussian
    amplitude _A that GForvard reads. So this helper updates _A too. Then the
    ABCD route keeps the same power as the array.
    """
    scale = 1.0 / np.sqrt(Power(Fin))
    Fout = Normal(Fin)
    Fout._A = Fin._A * scale
    return Fout


def _launch_aperture(tx):
    """Give the launch aperture diameter and obscuration of a transmitter.

    The bistatic override of olb.models.gaussian_efficiency applies: the
    Transmitter values win when they are set. If not, the owning Terminal
    values apply.
    """
    t = tx.transmitter
    aperture_m = t.aperture_m if t.aperture_m is not None else tx.aperture_m
    obscuration = (t.obscuration_ratio if t.obscuration_ratio is not None
                   else tx.obscuration_ratio)
    return aperture_m, obscuration


def _clip(Fin, aperture_m, obscuration_ratio):
    """Clip a field at a circular aperture with an optional obscuration."""
    Fout = CircAperture(Fin, aperture_m / 2)
    if obscuration_ratio > 0:
        Fout = CircScreen(Fout, obscuration_ratio * aperture_m / 2)
    return Fout


def _smf_focal_length(detector, aperture_m, lam):
    """Give the focal length of the single-mode-fibre coupling optic, in m.

    An explicit SMF.focal_length_m wins. Else SMF.optimal_focus derives the
    focal length from the a = 1.12 coupling parameter,
    f = pi*(D/2)*w_m/(lambda*1.12), with the SMF-28 mode field radius as the
    default w_m. Source: Shaklan and Roddier, Appl. Opt. 27 (1988) 2334,
    DOI 10.1364/AO.27.002334. This is the SAME rule as
    olb.models.coupling.terrestrial._smf_optics. The one-way dependency
    (waveoptics does not import models) is the reason for the second copy.

    Args:
        detector:   the SMF detector.
        aperture_m: the receive aperture diameter, in m.
        lam:        the wavelength, in m.

    Returns:
        The focal length in m, or None when the detector sets neither an
        explicit focal length nor optimal_focus.
    """
    if detector.focal_length_m is not None:
        return float(detector.focal_length_m)
    if detector.optimal_focus:
        w_m = detector.mode_field_radius_m
        if w_m is None:
            w_m = SMF28_MODE_FIELD_RADIUS_M
        return np.pi * (aperture_m / 2.0) * w_m / (lam * SMF_OPTIMAL_A)
    return None


def _smf_eta(detector, collected, aperture_m, lam):
    """Give the single-mode-fibre coupling efficiency of a collected field.

    The overlap reads the detector defocus: the fibre tip sits at
    z = f + SMF.defocus_m. A non-zero defocus needs a focal length, so the
    function resolves it with _smf_focal_length. See olb.waveoptics.smf.

    Args:
        detector:   the SMF detector.
        collected:  the clipped receive-plane Field.
        aperture_m: the receive aperture diameter, in m.
        lam:        the wavelength, in m.

    Returns:
        The coupling efficiency, a float between 0 and 1.

    Raises:
        ValueError: the detector sets a defocus but no focal length.
    """
    f_smf = _smf_focal_length(detector, aperture_m, lam)
    if detector.defocus_m != 0.0 and f_smf is None:
        raise ValueError(
            "SMF.defocus_m needs a focal length to make the defocus phase. "
            "Set SMF.focal_length_m, or set SMF.optimal_focus=True.")
    return float(coupling_efficiency(collected, aperture_m,
                                     defocus_m=detector.defocus_m,
                                     focal_length_m=f_smf))


def propagate_scenario(scenario, geometry, grid=None):
    """Propagate the transmit beam of a scenario to the receive aperture.

    The steps are: launch, launch-aperture clip, free-space propagation,
    receive-aperture clip, and the fibre coupling.

    The propagation step selects one of three routes. An almost untouched
    Gaussian takes the exact ABCD route (GForvard) at any range. A clipped
    field on a flat grid takes the Fresnel convolution. A clipped field on a
    SCALED grid takes the co-moving lens recipe (LensFresnel), which is the
    route for a long space link. GridSpec.for_scenario() selects the grid
    route, and grid.scaled records it.

    Args:
        scenario: a SpaceScenario or a TerrestrialScenario. The transmit
                  terminal needs a Transmitter.
        geometry: an object with slant_range_m. The range must be ONE value. A
                  sweep over the elevation belongs to the caller.
        grid:     an optional GridSpec. None derives the grid with
                  GridSpec.for_scenario().

    Returns:
        A WaveResult.

    Raises:
        ValueError: the geometry gives more than one range.
    """
    if grid is None:
        grid = GridSpec.for_scenario(scenario, geometry)

    range_m = np.asarray(geometry.slant_range_m, dtype=float)
    if range_m.size != 1:
        raise ValueError(
            f"propagate_scenario: the geometry gives {range_m.size} ranges. "
            "Give one range, and loop in the caller.")
    z = float(range_m.reshape(-1)[0])

    tx = scenario.tx_terminal
    rx = scenario.rx_terminal
    t = tx.transmitter
    lam = tx.wavelength_m

    # ---- the launch ----
    # A deliberately diverged beam starts at a virtual waist behind the
    # aperture. See olb.beam. The beam then has the asked-for radius in the
    # aperture plane.
    w_v, offset = virtual_waist(t.waist_m, t.divergence_rad, lam)
    F = GaussBeam(Begin(grid.size_m, lam, grid.n), w_v)
    F = _normalised_gauss(F)
    if offset > 0:
        F = GForvard(F, offset)
    launch = F
    power_launch = Power(launch)

    # ---- the launch-aperture clip ----
    tx_aperture_m, tx_obscuration = _launch_aperture(tx)
    clipped = _clip(launch, tx_aperture_m, tx_obscuration)
    power_clipped = Power(clipped)
    tx_truncation_db = _loss_db(power_clipped, power_launch)

    # ---- the propagation ----
    # An almost untouched Gaussian keeps the exact ABCD route, at any range.
    # A clipped field carries the aperture edge, so it needs a numerical
    # propagator: the Fresnel convolution on a flat grid, or the co-moving
    # lens recipe on a scaled grid. See olb.waveoptics.lenses and Schmidt,
    # DOI 10.1117/3.866274, Ch. 7.
    if 1.0 - power_clipped / power_launch < PURE_GAUSS_CLIP:
        at_rx = GForvard(launch, z)
        propagator = "GForvard"
    elif grid.scaled:
        # m comes from the same call that sized the grid, so the grid grows
        # by exactly the factor that the beam grows by.
        m = beam_magnification(scenario, z)
        f_lens = z / (m - 1)
        at_rx = Convert(LensFresnel(Lens(clipped, f_lens), -f_lens, z))
        propagator = "LensFresnel"
    else:
        at_rx = Fresnel(clipped, z)
        propagator = "Fresnel"
    power_at_rx = Power(at_rx)

    # ---- the receive-aperture clip ----
    collected = _clip(at_rx, rx.aperture_m, rx.obscuration_ratio)
    geometric_loss_db = _loss_db(Power(collected), power_at_rx)

    # ---- the fibre coupling ----
    smf_coupling_db = None
    if isinstance(rx.detector, SMF):
        # The overlap reads SMF.defocus_m: the fibre tip sits at z = f + dz.
        eta = _smf_eta(rx.detector, collected, rx.aperture_m, collected.lam)
        smf_coupling_db = _loss_db(eta, 1.0)

    stages = [("launch", launch), ("after tx clip", clipped),
              ("at rx plane", at_rx), ("after rx clip", collected)]
    return WaveResult(stages=stages, grid=grid,
                      tx_truncation_db=tx_truncation_db,
                      geometric_loss_db=geometric_loss_db,
                      smf_coupling_db=smf_coupling_db,
                      propagator=propagator)


if __name__ == '__main__':
    from ..geometry import CircularOrbit, HorizontalPath
    from ..models.gaussian_efficiency import tx_efficiency_loss_db
    from ..models.geometric import geometric_loss_db
    from ..scenario import (Channel, SpaceScenario, TerrestrialChannel,
                            TerrestrialScenario)
    from ..terminal import Terminal, Transmitter

    lam = 1550e-9

    def case(waist_m, tx_aperture_m, rx_aperture_m, range_m, grid,
             divergence_rad=None):
        """Run one case and give the fidelity-2 and the fidelity-0 numbers."""
        scn = TerrestrialScenario(
            near=Terminal(aperture_m=tx_aperture_m, wavelength_m=lam,
                          transmitter=Transmitter(waist_m=waist_m,
                                                  divergence_rad=divergence_rad)),
            far=Terminal(aperture_m=rx_aperture_m, wavelength_m=lam),
            channel=TerrestrialChannel(path_length_m=range_m))
        res = propagate_scenario(scn, HorizontalPath(range_m), grid=grid)
        analytic_tx = float(tx_efficiency_loss_db(tx_aperture_m, waist_m))
        analytic_geo = float(geometric_loss_db(range_m, waist_m, rx_aperture_m,
                                               wavelength=lam,
                                               divergence_rad=divergence_rad))
        return res, analytic_tx, analytic_geo

    # ---- case 1: far field, light truncation ----
    # waist 20 mm, tx aperture 100 mm (alpha = 2.5), rx aperture 200 mm,
    # range 30 km. The Rayleigh range is about 0.81 km, so 30 km is far field.
    zR = np.pi * 0.02 ** 2 / lam
    assert 30e3 > 20 * zR, zR
    far, far_tx, far_geo = case(0.02, 0.10, 0.20, 30e3, GridSpec(4.0, 2048))
    assert far.propagator == "Fresnel", far.propagator
    assert abs(far.geometric_loss_db - far_geo) < 0.1, \
        (far.geometric_loss_db, far_geo)
    assert abs(far.tx_truncation_db - far_tx) < 0.05, \
        (far.tx_truncation_db, far_tx)
    # No SMF detector, so there is no coupling number.
    assert far.smf_coupling_db is None
    # The stages carry the four labelled fields.
    assert [label for label, _ in far.stages] == [
        "launch", "after tx clip", "at rx plane", "after rx clip"]
    assert all(isinstance(F, Field) for _, F in far.stages)
    assert abs(Power(far.stages[0][1]) - 1.0) < 1e-9      # the launch is normal

    # ---- case 2: near field, hard truncation ----
    # waist 120 mm, tx aperture 150 mm (alpha = 0.625), range 1 km. The
    # Rayleigh range is about 29 km, so the receiver sits deep inside it. The
    # far-field analytic total does NOT hold here.
    zR_near = np.pi * 0.12 ** 2 / lam
    assert 1e3 < zR_near / 10, zR_near
    near, near_tx, near_geo = case(0.12, 0.15, 0.20, 1e3, GridSpec(1.0, 1024))
    fidelity2_total = near.tx_truncation_db + near.geometric_loss_db
    analytic_total = near_tx + near_geo
    assert abs(fidelity2_total - analytic_total) > 0.1, \
        (fidelity2_total, analytic_total)

    # ---- case 3: a diverged beam, almost no truncation ----
    # waist 50 mm, tx aperture 500 mm (alpha = 5), so the clip takes almost
    # nothing. The field stays a pure Gaussian, and the exact ABCD route runs.
    # A 3x diverged beam spreads more, so it costs more geometric loss.
    theta = 3 * lam / (np.pi * 0.05)
    div, div_tx, div_geo = case(0.05, 0.50, 0.20, 20e3, GridSpec(4.0, 1024),
                                divergence_rad=theta)
    assert div.propagator == "GForvard", div.propagator
    assert abs(div.geometric_loss_db - div_geo) < 0.1, \
        (div.geometric_loss_db, div_geo)
    assert div.tx_truncation_db < 1e-6
    collimated = case(0.05, 0.50, 0.20, 20e3, GridSpec(4.0, 1024))[0]
    assert div.geometric_loss_db > collimated.geometric_loss_db

    # ---- the physical limits ----
    for res in (far, near, div):
        assert res.tx_truncation_db >= 0.0 and res.geometric_loss_db >= 0.0
        assert res.grid.n >= 256

    # ---- an SMF detector gives a coupling loss ----
    smf_scn = TerrestrialScenario(
        near=Terminal(aperture_m=0.10, wavelength_m=lam,
                      transmitter=Transmitter(waist_m=0.02)),
        far=Terminal(aperture_m=0.20, wavelength_m=lam, detector=SMF()),
        channel=TerrestrialChannel(path_length_m=30e3))
    smf_res = propagate_scenario(smf_scn, HorizontalPath(30e3),
                                 grid=GridSpec(4.0, 2048))
    # A flat wavefront over the receive pupil cannot beat the Ruilier maximum
    # of 0.8145, which is 0.891 dB. See olb.waveoptics.smf.
    assert smf_res.smf_coupling_db > 0.89, smf_res.smf_coupling_db

    # ---- case 4: a space link on the co-moving grid ----
    # An uplink to a 600 km orbit at zenith. The launch waist is 50 mm, the
    # launch aperture is 100 mm (alpha = 1.0, a real truncation) with a 0.3
    # central obscuration, and the receive aperture is 500 mm. The automatic
    # sizer must select the scaled route, and the propagation must select
    # LensFresnel.
    space_scn = SpaceScenario(
        ground=Terminal(aperture_m=0.10, obscuration_ratio=0.3,
                        wavelength_m=lam,
                        transmitter=Transmitter(waist_m=0.05)),
        space=Terminal(aperture_m=0.50, wavelength_m=lam),
        direction="uplink", channel=Channel(altitude_m=600e3))
    orbit = CircularOrbit(altitude_m=600e3, elevation_deg=[90.0])
    z_space = float(np.max(orbit.slant_range_m))
    space = propagate_scenario(space_scn, orbit)      # the automatic grid
    assert space.grid.scaled, "a 600 km link must take the co-moving route"
    assert space.propagator == "LensFresnel", space.propagator

    space_tx = float(tx_efficiency_loss_db(0.10, 0.05, obscuration_ratio=0.3))
    # The receive terminal has no obscuration, so the launch obscuration is
    # charged one time only, by the transmit efficiency.
    space_geo = float(geometric_loss_db(z_space, 0.05, 0.50, wavelength=lam))

    # ONLY THE TOTAL COMPARES. The fidelity-0 pair splits the loss in a way
    # that the field does not: tx_efficiency_loss_db is an on-axis far-field
    # GAIN ratio, and geometric_loss_db is the power fraction of the
    # UNtruncated Gaussian in the receive aperture. The fidelity-2 pair is
    # plain power bookkeeping at each plane. The product of the fidelity-0
    # pair is the collected power fraction in the far field with a small
    # receiver, so the TOTALS agree. See the printed table below.
    p_launch = Power(space.stages[0][1])
    p_collected = Power(space.stages[3][1])
    space_total = _loss_db(p_collected, p_launch)
    assert abs(space_total - (space_tx + space_geo)) < 0.1, \
        (space_total, space_tx + space_geo)
    # The grid holds the far-field power: the co-moving grid keeps 95 percent
    # of the clipped power. The rest is the diffraction tail past the grid
    # edge, which the 500 mm receiver would never collect.
    keep = Power(space.stages[2][1]) / Power(space.stages[1][1])
    assert keep > 0.9, keep

    # ---- an array geometry is a caller loop, not a case ----
    try:
        propagate_scenario(smf_scn, HorizontalPath([1e3, 2e3]),
                           grid=GridSpec(4.0, 256))
        raise AssertionError("an array range must raise ValueError")
    except ValueError:
        pass

    header = (f"{'case':<26}{'fid-2 tx':>10}{'fid-0 tx':>10}"
              f"{'fid-2 geo':>11}{'fid-0 geo':>11}{'fid-2 tot':>11}"
              f"{'fid-0 tot':>11}")
    print(header)
    print("-" * len(header))
    for name, res, a_tx, a_geo in (
            ("far field, alpha=2.5", far, far_tx, far_geo),
            ("near field, alpha=0.625", near, near_tx, near_geo),
            ("far field, 3x diverged", div, div_tx, div_geo)):
        print(f"{name:<26}{res.tx_truncation_db:>10.3f}{a_tx:>10.3f}"
              f"{res.geometric_loss_db:>11.3f}{a_geo:>11.3f}"
              f"{res.tx_truncation_db + res.geometric_loss_db:>11.3f}"
              f"{a_tx + a_geo:>11.3f}")
    print("")
    print("all the numbers are losses in positive dB.")
    print(f"far-field propagator      {far.propagator}")
    print(f"near-field propagator     {near.propagator}")
    print(f"diverged propagator       {div.propagator}")
    print(f"SMF coupling loss, 30 km  {smf_res.smf_coupling_db:.3f} dB")
    print("")
    print("space uplink, 600 km at zenith, the co-moving grid:")
    print(f"  range                   {z_space * 1e-3:11.1f} km")
    print(f"  propagator              {space.propagator:>11}")
    print(f"  grid, launch side       {space.grid.size_m:11.3f} m "
          f"({space.grid.n} pixels)")
    print(f"  grid, receive side      "
          f"{space.grid.size_m * beam_magnification(space_scn, z_space):11.3f} m")
    print(f"  power kept on the grid  {keep:11.4f}")
    print(f"  fid-2 launch truncation {space.tx_truncation_db:11.3f} dB")
    print(f"  fid-0 launch truncation {space_tx:11.3f} dB")
    print(f"  fid-2 geometric spread  {space.geometric_loss_db:11.3f} dB")
    print(f"  fid-0 geometric spread  {space_geo:11.3f} dB")
    print(f"  fid-2 TOTAL, collected  {space_total:11.3f} dB")
    print(f"  fid-0 TOTAL             {space_tx + space_geo:11.3f} dB")
    print(f"  difference              {space_total - (space_tx + space_geo):11.3f} dB")
    print("  The split differs, because the two fidelities cut the loss in")
    print("  two places. The TOTAL is the number that compares.")
    print("self-check passed")

"""The grid specification for a fidelity-2 field propagation.

A field propagation needs two numbers: the physical side of the square grid
and the number of pixels along that side. The two numbers fight each other. A
wide grid holds the far-field beam, but it makes the pixels coarse. A fine
pixel resolves the launch aperture, but it makes the grid small.

GridSpec holds the two numbers. Build it directly for manual work. Call
GridSpec.for_scenario() to derive it from an olb scenario and geometry.

TWO ROUTES. A short link takes a FLAT grid: one side for the whole path. A
long space link makes the beam grow by a factor of 100 or more, so a flat grid
must hold the far-field beam AND resolve the launch aperture. No practical
pixel count does both. Such a link takes the SCALED (co-moving) route of
olb.waveoptics.lenses: the grid starts at the launch plane, and it grows with
the beam by the magnification m = w(z)/w(0). The `scaled` attribute records
the choice, and propagate_scenario() reads it to select the propagator.

The class WARNS. It does not raise. A link that NEITHER route samples well
gets an honest warning, because that is better than a silent bad answer.

Sources:
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB, DOI 10.1117/3.866274. The angular-spectrum range limit
  z_max = N * dx^2 / lambda is constraint 4, Ch. 7, Eq. (7.59), printed
  p. 127, at m = 1. The same limit is the step cap of Ch. 8, Eq. (8.24),
  printed p. 144. The pixel count goes up to the next power of two, as
  Listing 7.1, line 11, printed p. 124, does.
  THE EXTENT RULE AND THE PIXEL-PER-FEATURE RULE OF THIS MODULE ARE NOT THE
  BOOK'S. The book sizes the grid from the illuminated diameter and the
  region of interest, Ch. 7, Eq. (7.18), printed p. 120, and it lets the
  wrapped light come up to the edge of that region. This module keeps a
  fixed margin around the beam instead. The book also gives no fixed
  pixels-per-feature number: it picks 50 points across the aperture in
  Listing 7.1, printed p. 124, and 30 points in the Ch. 8 example, printed
  p. 144. See docs/schmidt-crosscheck.md, gaps S-07 and S-16.
"""

import warnings
from dataclasses import dataclass

import numpy as np

from ..beam import free_space_radius

# The floor on the pixel count. It is a convenience of this module, not
# physics. Schmidt, DOI 10.1117/3.866274, gives no floor: the three worked
# examples use N = 128 (Ch. 7, printed p. 123), N = 512 (Ch. 7, printed
# p. 128) and N = 128 (Ch. 8, printed p. 144).
N_MIN = 256


def beam_magnification(scenario, z):
    """Give the grid magnification m = w(z)/w(0) of the transmit beam.

    The scaled route of olb.waveoptics.lenses grows the grid by this factor,
    so the grid holds the same number of beam radii at each end of the link.
    A deliberately diverged beam brings its own w(z), because
    olb.beam.free_space_radius reads the Transmitter divergence.

    The sizer and propagate_scenario() both call this function, so the two
    always use the same m.

    Args:
        scenario: a SpaceScenario or a TerrestrialScenario.
        z:        the range, in m.

    Returns:
        The magnification, a float. It is 1.0 at z = 0.
    """
    tx = scenario.tx_terminal
    t = tx.transmitter
    return float(free_space_radius(t.waist_m, z, t.divergence_rad,
                                   tx.wavelength_m)) / t.waist_m


def _features(aperture_m, obscuration_ratio):
    """List the hard edges of one aperture, as radii in m.

    The pixel must resolve the smallest edge. A central obscuration is an
    edge too, and it is smaller than the aperture.
    """
    radii = [aperture_m / 2]
    if obscuration_ratio > 0:
        radii.append(obscuration_ratio * aperture_m / 2)
    return radii


@dataclass(frozen=True)
class GridSpec:
    """The square simulation grid.

    Attributes:
        size_m: the physical side of the grid, in m. It is the side at the
                LAUNCH plane when scaled is True.
        n:      the number of pixels along one side.
        scaled: True selects the scaled (co-moving) route. The grid then
                grows with the beam, and the side at the receive plane is
                size_m times the magnification. False keeps a flat grid.
    """

    size_m: float
    n: int
    scaled: bool = False

    @property
    def pixel_m(self):
        """The distance between two pixels, in m."""
        return self.size_m / self.n

    @classmethod
    def for_scenario(cls, scenario, geometry, guard=4.0, pixels_per_feature=16,
                     n_max=4096):
        """Derive a grid from a scenario and a geometry.

        The function tries the FLAT route first, and it takes the SCALED
        (co-moving) route when the flat route cannot resolve the apertures.

        The FLAT EXTENT rule: size = guard * 2 * r_max. r_max is the largest of
        the free-space beam radius at the launch plane, the free-space beam
        radius at the longest range, the transmit aperture radius, and the
        receive aperture radius. The guard keeps the beam away from the grid
        edge, because the FFT propagators are periodic. THE GUARD IS AN olb
        RULE, not a book rule. Schmidt, DOI 10.1117/3.866274, Ch. 7,
        Eqs. (7.18) and (7.20), printed p. 120, size the grid from the
        illuminated diameter D_illum and the region of interest D2, and they
        let the wrapped light come up to the edge of D2. The observation
        extent D2 never enters this rule. See docs/schmidt-crosscheck.md,
        gaps S-07 and S-16.

        The SCALED EXTENT rule: the grid starts at the launch plane, so
        size = guard * 2 * max(launch radii) only. The launch radii are the
        transmit waist and the transmit aperture radius. The grid then grows
        with the beam by the magnification m = w(z)/w(0). See
        olb.waveoptics.lenses. THE BOOK GIVES NO EQUATION FOR A CO-MOVING
        GRID: Schmidt, DOI 10.1117/3.866274, Ch. 6, text, printed p. 87, names
        the Coles and Rubio angular-grid method and does not develop it. The
        book's own answer to the same problem is the scaling parameter m of
        Ch. 6, Eq. (6.65), printed p. 100, on a FLAT grid. See gap S-13.

        The RESOLUTION rule: the smallest feature gets pixels_per_feature
        pixels across it. The features are the transmit waist and the hard
        edges of the two apertures (each aperture radius, and each central
        obscuration radius). The scaled route measures a receive feature at
        the LAUNCH plane, so it divides that feature by m. The pixel count
        goes up to the next power of two, as Schmidt, DOI 10.1117/3.866274,
        Listing 7.1, line 11, printed p. 124, does, and it stays in the
        interval [256, n_max]. The book gives no pixels-per-feature equation.
        It picks 50 points across the aperture in Listing 7.1, printed p. 124,
        and 30 points in the Ch. 8 example, printed p. 144, so the default of
        16 is coarser than both.

        The transmit aperture obeys the bistatic rule of
        olb.models.gaussian_efficiency: the Transmitter aperture_m and
        obscuration_ratio win when they are set. If not, the owning Terminal
        values apply.

        Args:
            scenario:           a SpaceScenario or a TerrestrialScenario.
            geometry:           an object with slant_range_m (a scalar or an
                                array). The function takes the largest range.
            guard:              the ratio of the grid half-side to r_max.
            pixels_per_feature: the number of pixels across the smallest
                                feature.
            n_max:              the largest pixel count.

        Returns:
            A GridSpec. The scaled attribute says which route it holds.

        Warns:
            UserWarning: NEITHER route resolves the smallest feature, or the
                range of a flat grid is longer than forvard_max_z().
        """
        tx = scenario.tx_terminal
        rx = scenario.rx_terminal
        t = tx.transmitter
        if t is None:
            raise ValueError('GridSpec.for_scenario: the transmit terminal '
                             'has no Transmitter')
        lam = tx.wavelength_m
        z_max = float(np.max(np.asarray(geometry.slant_range_m, dtype=float)))

        # The bistatic override. See olb.models.gaussian_efficiency.
        tx_aperture_m = t.aperture_m if t.aperture_m is not None else tx.aperture_m
        tx_obscuration = (t.obscuration_ratio if t.obscuration_ratio is not None
                          else tx.obscuration_ratio)

        w0 = t.waist_m
        w_end = float(free_space_radius(w0, z_max, t.divergence_rad, lam))
        launch_feature = min([w0] + _features(tx_aperture_m, tx_obscuration))
        rx_feature = min(_features(rx.aperture_m, rx.obscuration_ratio))

        def build(size_m, feature, scaled):
            """Round the pixel count up, clamp it, and count the pixels won."""
            pixel_wanted = feature / (pixels_per_feature / 2)
            n_wanted = int(2 ** np.ceil(np.log2(size_m / pixel_wanted)))
            grid = cls(size_m=size_m, n=int(min(max(n_wanted, N_MIN), n_max)),
                       scaled=scaled)
            return grid, 2 * feature / grid.pixel_m

        # ---- the flat route ----
        r_max = max(w0, w_end, tx_aperture_m / 2, rx.aperture_m / 2)
        flat, flat_pixels = build(guard * 2 * r_max,
                                  min(launch_feature, rx_feature), False)
        if flat_pixels >= pixels_per_feature:
            z_limit = forvard_max_z(flat, lam)
            if z_max > z_limit:
                warnings.warn(
                    f"GridSpec: the range ({z_max:.3g} m) is past the sampling "
                    f"limit N*dx^2/lambda = {z_limit:.3g} m. Schmidt, "
                    f"DOI 10.1117/3.866274. The spectral propagator (Forvard) "
                    f"aliases here. Use Fresnel, or accept the error.")
            return flat

        # ---- the scaled (co-moving) route ----
        # The receive features come back to the launch plane through m, so the
        # one pixel serves both planes. The recipe needs a lens of the focal
        # length z/(m - 1), so a beam that does not grow (m = 1) has no scaled
        # route. See olb.waveoptics.lenses.
        m = beam_magnification(scenario, z_max)
        scaled, scaled_pixels = build(
            guard * 2 * max(w0, tx_aperture_m / 2),
            min(launch_feature, rx_feature / m), True)
        if m > 1.0 + 1e-9 and scaled_pixels >= pixels_per_feature:
            return scaled

        warnings.warn(
            f"GridSpec: NEITHER route resolves the smallest feature under the "
            f"n_max clamp ({n_max}). The flat grid gives {flat_pixels:.1f} "
            f"pixels, and the co-moving grid gives {scaled_pixels:.1f}. "
            f"{pixels_per_feature} were asked for. The aperture edges are "
            f"coarse. Raise n_max, or accept the error.")
        return flat


def forvard_max_z(grid, wavelength_m):
    """Calculate the longest range that the grid samples well.

    z_max = N * dx^2 / lambda. Past this range the quadratic phase of the
    transfer function turns faster than one sample, so the spectral propagator
    aliases.

    This IS constraint 4 of Schmidt (2010), DOI 10.1117/3.866274, Ch. 7,
    Eq. (7.59), printed p. 127: N >= lambda z / (dx1 dx2). A flat grid has
    dx1 = dx2 = dx, so the rule inverts to z <= N dx^2 / lambda. The same
    expression is the partial-propagation step cap of Ch. 8, Eq. (8.24),
    printed p. 144, with min(dx1, dxn) = dx, and the turbulent form of Ch. 9,
    Eq. (9.89), printed p. 174. The constant is DERIVED, not a guess.
    olb.waveoptics.schmidt.sampling.angular_spectrum_max_z reproduces it.

    Args:
        grid:          a GridSpec.
        wavelength_m:  the wavelength, in m.

    Returns:
        The largest well-sampled range, in m.
    """
    return grid.n * grid.pixel_m ** 2 / wavelength_m


if __name__ == '__main__':
    from ..geometry import CircularOrbit, HorizontalPath
    from ..scenario import (Channel, SpaceScenario, TerrestrialChannel,
                            TerrestrialScenario)
    from ..terminal import Terminal, Transmitter

    lam = 1550e-9

    # ---- the manual route ----
    g = GridSpec(size_m=4.0, n=1024)
    assert abs(g.pixel_m - 4.0 / 1024) < 1e-18
    assert abs(forvard_max_z(g, lam) - 1024 * g.pixel_m ** 2 / lam) < 1e-9

    # ---- a terrestrial reference case ----
    # waist 10 mm, tx aperture 100 mm, rx aperture 12.7 mm, range 10 km.
    terr = TerrestrialScenario(
        near=Terminal(aperture_m=0.1, wavelength_m=lam,
                      transmitter=Transmitter(waist_m=0.01)),
        far=Terminal(aperture_m=12.7e-3, wavelength_m=lam),
        channel=TerrestrialChannel(path_length_m=10e3))
    path = HorizontalPath(10e3)
    # n_max=8192, because a 12.7 mm receive aperture on a 4 m grid is a small
    # feature. The default 4096 clamp gives only 13 pixels across it.
    gt = GridSpec.for_scenario(terr, path, n_max=8192)

    w_end = float(free_space_radius(0.01, 10e3, None, lam))
    assert abs(gt.size_m - 8 * w_end) < 1e-9          # guard 4, r_max = w(z)
    # The Gaussian edge intensity is small at both ends of the path.
    edge = gt.size_m / 2
    for w in (0.01, w_end):
        assert np.exp(-2 * edge ** 2 / w ** 2) < 1e-10
    # The receive aperture gets 16 pixels or more.
    pixels_on_rx = 12.7e-3 / gt.pixel_m
    assert pixels_on_rx >= 16, pixels_on_rx

    # ---- a space case takes the co-moving route ----
    # waist 50 mm, launch aperture 100 mm with a 0.3 obscuration, receive
    # aperture 500 mm, 600 km at zenith. A flat grid must hold w(z) = 5.9 m
    # AND resolve a 50 mm waist, so it fails. The scaled route sizes the
    # LAUNCH plane only, and it resolves both ends.
    orbit = CircularOrbit(altitude_m=600e3, elevation_deg=[90.0])
    z600 = float(np.max(orbit.slant_range_m))
    space = SpaceScenario(
        ground=Terminal(aperture_m=0.1, obscuration_ratio=0.3, wavelength_m=lam,
                        transmitter=Transmitter(waist_m=0.05)),
        space=Terminal(aperture_m=0.5, wavelength_m=lam),
        direction="uplink", channel=Channel(altitude_m=600e3))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("error")              # any warning is a failure
        gs = GridSpec.for_scenario(space, orbit)
    assert gs.scaled, "a 600 km link must take the co-moving route"
    assert abs(gs.size_m - 8 * 0.05) < 1e-12        # guard 4 on the launch
    m600 = beam_magnification(space, z600)
    assert m600 > 100, m600
    # The receive aperture gets 16 pixels or more, after the magnification.
    pixels_on_rx_space = 0.5 / (gs.size_m * m600 / gs.n)
    assert pixels_on_rx_space >= 16, pixels_on_rx_space
    # The launch obscuration is a hard edge too, and it stays resolved.
    assert 2 * (0.3 * 0.05) / gs.pixel_m >= 16

    # ---- one deliberately impossible case still warns ----
    # A 300 mm launch aperture with a 50 mm receive aperture at 600 km. The
    # scaled route must resolve a 25 mm receive radius through m = 118, which
    # is 0.2 mm at the launch plane, on a 2.4 m grid. No n_max holds that.
    hard = SpaceScenario(
        ground=Terminal(aperture_m=0.3, wavelength_m=lam,
                        transmitter=Transmitter(waist_m=0.05)),
        space=Terminal(aperture_m=0.05, wavelength_m=lam),
        direction="uplink", channel=Channel(altitude_m=600e3))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        gh = GridSpec.for_scenario(hard, orbit)
    assert len(caught) == 1, [str(w.message) for w in caught]
    assert "NEITHER route" in str(caught[0].message)
    assert gh.n == 4096 and not gh.scaled           # the clamp holds

    # A scenario with no Transmitter is not a propagation case.
    try:
        GridSpec.for_scenario(
            TerrestrialScenario(near=Terminal(aperture_m=0.1),
                                far=Terminal(aperture_m=0.1)), path)
        raise AssertionError("a missing Transmitter must raise ValueError")
    except ValueError:
        pass

    print("terrestrial reference case, 10 km:")
    print(f"  beam radius at range  {w_end * 1e3:9.2f} mm")
    print(f"  grid side             {gt.size_m:9.3f} m")
    print(f"  pixels per side       {gt.n:9d}")
    print(f"  pixel pitch           {gt.pixel_m * 1e3:9.4f} mm")
    print(f"  pixels across rx      {pixels_on_rx:9.1f}")
    print(f"  Forvard max range     {forvard_max_z(gt, lam):9.1f} m")
    print("")
    print("space case, 600 km, the co-moving route:")
    print(f"  magnification m       {m600:9.1f}")
    print(f"  grid side at launch   {gs.size_m:9.3f} m")
    print(f"  grid side at receiver {gs.size_m * m600:9.3f} m")
    print(f"  pixels per side       {gs.n:9d}")
    print(f"  pixel at launch       {gs.pixel_m * 1e3:9.4f} mm")
    print(f"  pixel at receiver     {gs.pixel_m * m600 * 1e3:9.4f} mm")
    print(f"  pixels across rx      {pixels_on_rx_space:9.1f}")
    print("")
    print("space case, 600 km, a 50 mm receiver: no route works:")
    print(f"  grid side             {gh.size_m:9.3f} m")
    print(f"  pixels per side       {gh.n:9d}")
    for w in caught:
        print(f"  warning: {str(w.message).split('.')[0]}.")
    print("self-check passed")

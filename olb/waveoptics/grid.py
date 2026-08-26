"""The grid specification for a fidelity-2 field propagation.

A field propagation needs two numbers: the physical side of the square grid
and the number of pixels along that side. The two numbers fight each other. A
wide grid holds the far-field beam, but it makes the pixels coarse. A fine
pixel resolves the launch aperture, but it makes the grid small.

GridSpec holds the two numbers. Build it directly for manual work. Call
GridSpec.for_scenario() to derive it from an olb scenario and geometry.

The class WARNS. It does not raise. A long space link cannot be sampled well on
a grid of a practical size, so an honest warning is better than a silent bad
answer.

Sources:
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB, DOI 10.1117/3.866274, Ch. 6 and Ch. 7. The grid extent rule, the
  pixel-per-feature rule, and the angular-spectrum range limit
  z_max = N * dx^2 / lambda.
"""

import warnings
from dataclasses import dataclass

import numpy as np

from ..beam import free_space_radius

N_MIN = 256


@dataclass(frozen=True)
class GridSpec:
    """The square simulation grid.

    Attributes:
        size_m: the physical side of the grid, in m.
        n:      the number of pixels along one side.
    """

    size_m: float
    n: int

    @property
    def pixel_m(self):
        """The distance between two pixels, in m."""
        return self.size_m / self.n

    @classmethod
    def for_scenario(cls, scenario, geometry, guard=4.0, pixels_per_feature=16,
                     n_max=4096):
        """Derive a grid from a scenario and a geometry.

        The EXTENT rule: size = guard * 2 * r_max. r_max is the largest of the
        free-space beam radius at the launch plane, the free-space beam radius
        at the longest range, the transmit aperture radius, and the receive
        aperture radius. The guard keeps the beam away from the grid edge,
        because the FFT propagators are periodic. See Schmidt,
        DOI 10.1117/3.866274, Ch. 6.

        The RESOLUTION rule: the smallest feature gets pixels_per_feature pixels
        across it. The features are the transmit waist, the transmit aperture
        radius, the receive aperture radius, and the receive obscuration radius.
        The pixel count goes up to the next power of two, and it stays in the
        interval [256, n_max].

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
            A GridSpec.

        Warns:
            UserWarning: the range is longer than forvard_max_z(), or the n_max
                clamp leaves the smallest feature with too few pixels.
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

        w0 = t.waist_m
        w_end = float(free_space_radius(w0, z_max, t.divergence_rad, lam))
        r_max = max(w0, w_end, tx_aperture_m / 2, rx.aperture_m / 2)
        size_m = guard * 2 * r_max

        features = [w0, tx_aperture_m / 2, rx.aperture_m / 2]
        rx_obscuration_m = rx.obscuration_ratio * rx.aperture_m / 2
        if rx_obscuration_m > 0:
            features.append(rx_obscuration_m)
        feature = min(features)

        pixel_wanted = feature / (pixels_per_feature / 2)
        n_wanted = int(2 ** np.ceil(np.log2(size_m / pixel_wanted)))
        n = int(min(max(n_wanted, N_MIN), n_max))
        grid = cls(size_m=size_m, n=n)

        # The honest warnings.
        pixels_on_feature = 2 * feature / grid.pixel_m
        if pixels_on_feature < pixels_per_feature:
            warnings.warn(
                f"GridSpec: the n_max clamp ({n_max}) leaves only "
                f"{pixels_on_feature:.1f} pixels across the smallest feature "
                f"({feature * 1e3:.3f} mm radius). {pixels_per_feature} were "
                f"asked for. The aperture edges are coarse.")
        z_limit = forvard_max_z(grid, lam)
        if z_max > z_limit:
            warnings.warn(
                f"GridSpec: the range ({z_max:.3g} m) is past the sampling "
                f"limit N*dx^2/lambda = {z_limit:.3g} m. Schmidt, "
                f"DOI 10.1117/3.866274. The spectral propagator (Forvard) "
                f"aliases here. Use Fresnel, or accept the error.")
        return grid


def forvard_max_z(grid, wavelength_m):
    """Calculate the longest range that the grid samples well.

    z_max = N * dx^2 / lambda. Past this range the quadratic phase of the
    transfer function turns faster than one sample, so the spectral propagator
    aliases. See Schmidt, DOI 10.1117/3.866274, Ch. 6.

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

    # ---- a space case warns ----
    space = SpaceScenario(
        ground=Terminal(aperture_m=0.3, wavelength_m=lam,
                        transmitter=Transmitter(waist_m=0.05)),
        space=Terminal(aperture_m=0.05, wavelength_m=lam),
        direction="uplink", channel=Channel(altitude_m=600e3))
    orbit = CircularOrbit(altitude_m=600e3, elevation_deg=[90.0])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        gs = GridSpec.for_scenario(space, orbit)
    assert len(caught) >= 1, "a 600 km link must warn"
    assert gs.n == 4096                                # the clamp holds

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
    print("space case, 600 km:")
    print(f"  grid side             {gs.size_m:9.3f} m")
    print(f"  pixels per side       {gs.n:9d}")
    print(f"  pixel pitch           {gs.pixel_m * 1e3:9.4f} mm")
    print(f"  Forvard max range     {forvard_max_z(gs, lam):9.3g} m")
    for w in caught:
        print(f"  warning: {str(w.message).split('.')[0]}.")
    print("self-check passed")

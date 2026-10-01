"""Tabulate the velocity-aberration loss of a corner-cube return.

A corner cube returns the beam along the incoming ray in its own frame. The
satellite moves, so in the ground frame the return lobe leaves 2 v_perp / c
off the station (the point-ahead angle; Degnan, DOI 10.1029/GD025p0133). The
lobe of an unobscured cube of aperture D is the Airy pattern, so the station
reads the fraction [2 J1(x) / x]^2 of the on-axis peak, x = pi D theta / lambda
(Born and Wolf, DOI 10.1017/CBO9781139644181, Sec. 8.5.2). This is
`olb.links.retro_space.retro_velocity_aberration_term`.

The script gives:
1. the loss against x (where the station sits in the lobe);
2. the loss of typical cubes at typical angles;
3. the best cube size: the return to the station scales as
   D^2 (the capture) * D^2 (the on-axis gain) * [2 J1(x) / x]^2, which is
   proportional to [x J1(x)]^2 and peaks where d(x J1)/dx = x J0(x) = 0, the
   first zero of J0, x = 2.405;
4. the loss over one day of real ISS passes;
5. the loss of a 5 cm and a 35 mm cube from the flat CircularOrbit angle against the
   overhead spherical angle, 20 to 90 deg, 500 and 1500 km (backlog 0-P18).

Run from the repository root:

    python -m validation.retro_velocity_aberration.airy_aberration
"""

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import j1, jn_zeros

from olb.geometry import _C, CircularOrbit, Satellite, TLEPass
from olb.links.retro_space import retro_velocity_aberration_term
from olb.scenario import Channel, SpaceScenario
from olb.terminal import Terminal, Transmitter
from validation.point_ahead_geometry.tle_point_ahead import ISS, STATION

LAM = 1550e-9


def airy_loss_db(x):
    """Give -10 log10([2 J1(x) / x]^2), with the x -> 0 limit 0 dB."""
    x = np.asarray(x, dtype=float)
    xs = np.where(x > 0, x, 1.0)
    frac = np.where(x > 0, (2 * j1(xs) / xs) ** 2, 1.0)
    with np.errstate(divide="ignore"):
        return -10 * np.log10(frac)


def retro(D, altitude_m=500e3):
    """Give a retro scenario with a cube of aperture D."""
    return SpaceScenario(
        ground=Terminal(aperture_m=0.4, wavelength_m=LAM,
                        transmitter=Transmitter(waist_m=0.06)),
        space=Terminal(aperture_m=D, wavelength_m=LAM),
        direction="retro", channel=Channel(altitude_m=altitude_m))


def main():
    print("1. the loss against x = pi D theta / lambda")
    for x, where in ((0.0, "the lobe centre"), (1.0, "about half-way out"),
                     (2.405, "the best cube size"), (3.8317, "the first null"),
                     (5.1356, "the first bright ring")):
        print(f"   x = {x:6.3f}  {airy_loss_db(x) + 0.0:7.2f} dB  ({where})")

    print(f"\n2. typical cubes, lambda = {LAM * 1e9:.0f} nm")
    geom = CircularOrbit(altitude_m=500e3, elevation_deg=[30.0])
    theta_geom = float(np.ravel(geom.point_ahead_rad)[0])
    print(f"   {'D':>6} {'lambda/D':>10} {'x @ 50 urad':>12} {'loss':>8}"
          f"  {'x @ 500 km 30 deg':>18} {'loss':>8}")
    for D in (0.01, 0.02, 0.05, 0.10):
        x50 = np.pi * D * 50e-6 / LAM
        term = retro_velocity_aberration_term(retro(D), geom)
        print(f"   {D * 100:4.0f} cm {LAM / D * 1e6:7.1f} urad {x50:12.2f} "
              f"{airy_loss_db(x50):6.2f} dB  {float(np.ravel(term.meta['x'])[0]):18.2f} "
              f"{float(np.ravel(term.mean_db)[0]):6.2f} dB")
    print(f"   (the 500 km, 30 deg CircularOrbit angle is "
          f"{theta_geom * 1e6:.1f} urad)")

    print("\n3. the best cube size")
    best = minimize_scalar(lambda x: -(x * j1(x)) ** 2, bounds=(0.5, 3.8),
                           method="bounded").x
    print(f"   the numeric maximum of [x J1(x)]^2 is at x = {best:.4f}; the "
          f"first zero of J0 is {jn_zeros(0, 1)[0]:.4f}")
    k = best / np.pi
    print(f"   so D_opt = {k:.3f} lambda / theta: "
          + ", ".join(f"{k * LAM / th * 100:.1f} cm at {th * 1e6:.0f} urad"
                      for th in (20e-6, 30e-6, 50e-6)))

    print("\n4. one day of ISS passes, elevation above 10 deg")
    iss = TLEPass.from_window(*ISS, **STATION, start_utc=(2014, 1, 21, 0, 0, 0),
                              duration_s=86400, step_s=30.0)
    vis = iss.elevation_deg > 10.0
    for D in (0.01, 0.02, 0.05):
        loss = np.asarray(retro_velocity_aberration_term(
            retro(D, 420e3), iss).mean_db)[vis]
        print(f"   {D * 100:3.0f} cm cube: {loss.min():6.2f} to "
              f"{loss.max():6.2f} dB (median {np.median(loss):.2f} dB)")

    print("\n5. the flat CircularOrbit angle against the overhead spherical "
          "angle (backlog 0-P18)")
    for D, h in ((D, h) for D in (0.05, 0.035) for h in (500e3, 1500e3)):
        print(f"   {D * 1e3:.0f} mm, {h / 1e3:4.0f} km: {'el':>3} {'th flat':>8} {'th sph':>7}"
              f" {'loss flat':>10} {'loss sph':>9} {'delta':>7}")
        for el in range(20, 91, 10):
            # The old flat-Earth form 2 v_orb sin(el) / c, against the
            # CircularOrbit overhead pass on a sphere, 2 v_orb cos(eta) / c
            # (Degnan, DOI 10.1029/GD025p0133; see ../point_ahead_geometry/).
            flat = 2 * Satellite(h).orbital_speed * np.sin(np.radians(el)) / _C
            sph = float(np.ravel(CircularOrbit(h, [el]).point_ahead_rad)[0])
            lf, ls = (airy_loss_db(np.pi * D * t / LAM) for t in (flat, sph))
            print(f"   {'':14} {el:3d} {flat * 1e6:8.1f} {sph * 1e6:7.1f}"
                  f" {lf:10.1f} {ls:9.1f} {ls - lf:+7.1f}")


if __name__ == "__main__":
    main()

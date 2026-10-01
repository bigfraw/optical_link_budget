"""The LAGEOS-1 pass of 2026-10-01 over Yarragadee: the point-ahead angle and the single-cube Airy loss.

This script records the numbers that start the array-speckle study (see
README.md). It propagates the CelesTrak TLE of NORAD 8820 (LAGEOS-1), finds
the time at which the TLE matches the RC700 mount readout, and applies
`retro_velocity_aberration_term` to ONE unspoiled 3.81 cm cube. LAGEOS has
426 SPOILED cubes, so the Airy loss is NOT the LAGEOS loss; the script shows
why the single-cube model does not apply (the offset is past the first Airy
null at 532 nm and sits near it at 1064 nm).

Run from the repository root:

    python -m validation.retro_array_speckle.lageos_pass
"""

from types import SimpleNamespace as NS

import numpy as np
from skyfield.api import load, wgs84, EarthSatellite

from olb.geometry import TLEPass
from olb.links.retro_space import retro_velocity_aberration_term

# CelesTrak GP, fetched 2026-10-01 05:20 UTC (epoch 2026-09-30 12:56 UTC).
L1 = "1 08820U 76039A   26273.53930453  .00000007  00000+0  00000+0 0  9990"
L2 = "2 08820 109.8213 215.3928 0044588 279.7080 260.7759  6.38664837919866"
# Yarragadee, MOBLAS-5 (ILRS 7090).
STATION = dict(lat_deg=-29.0464, lon_deg=115.3467, alt_m=244.0)
# The RC700 mount readout (J2000).
RA_DEG = (18 + 17 / 60 + 54.8 / 3600) * 15
DEC_DEG = -(49 + 33 / 60 + 28 / 3600)
D_CUBE = 0.0381                                   # one LAGEOS cube [m]
ARCSEC = np.pi / 180 / 3600


def main():
    ts = load.timescale()
    p = TLEPass.from_window(L1, L2, **STATION, start_utc=(2026, 10, 1, 5, 0, 0),
                            duration_s=3900, step_s=10.0, name="LAGEOS-1")
    topo = (EarthSatellite(L1, L2, "LAGEOS-1", ts)
            - wgs84.latlon(STATION["lat_deg"], STATION["lon_deg"],
                           STATION["alt_m"]))
    ra, dec, _ = topo.at(p.times).radec()
    ra0, de0 = np.deg2rad(RA_DEG), np.deg2rad(DEC_DEG)
    cos_sep = (np.sin(dec.radians) * np.sin(de0)
               + np.cos(dec.radians) * np.cos(de0) * np.cos(ra.radians - ra0))
    sep = np.rad2deg(np.arccos(np.clip(cos_sep, -1, 1)))
    i = int(np.argmin(sep))
    assert sep[i] < 0.5, "the TLE does not match the mount readout"
    print(f"readout match: {p.times[i].utc_strftime()}, "
          f"{sep[i] * 60:.1f} arcmin, el {p.elevation_deg[i]:.1f} deg, "
          f"range {p.slant_range_m[i] / 1e3:.0f} km")

    up = p.elevation_deg > 0
    for lam in (532e-9, 1064e-9):
        cube = NS(space=NS(aperture_m=D_CUBE, wavelength_m=lam))
        loss = np.atleast_1d(retro_velocity_aberration_term(cube, p).mean_db)
        x = np.pi * D_CUBE * p.point_ahead_rad[i] / lam
        print(f"\n{lam * 1e9:.0f} nm, one unspoiled {D_CUBE * 100:.2f} cm cube "
              f"(first Airy null at x = 3.83)")
        print(f"  at the readout: PAA {p.point_ahead_rad[i] * 1e6:.2f} urad "
              f"({p.point_ahead_rad[i] / ARCSEC:.2f} arcsec), x = {x:.2f}, "
              f"loss {loss[i]:.2f} dB")
        print(f"  over the pass: PAA {p.point_ahead_rad[up].min() / ARCSEC:.2f}"
              f" to {p.point_ahead_rad[up].max() / ARCSEC:.2f} arcsec, loss "
              f"{loss[up].min():.2f} to {loss[up].max():.2f} dB")


if __name__ == "__main__":
    main()

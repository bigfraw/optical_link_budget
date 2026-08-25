'''
The Chapter 12 slant path: a real SpaceScenario swept across elevation.

Andrews and Phillips, "Laser Beam Propagation through Random Media", 2nd ed.,
SPIE Press (2005), DOI: 10.1117/3.626196, Ch. 12:
    Eq. (15), printed p. 490   the Kolmogorov spectrum, the ONLY one Ch. 12 uses
    Eq. (38), printed p. 495   downlink point scintillation index
    Eq. (39), printed p. 496   downlink aperture-averaged index, WEAK theory
    Eq. (50), printed p. 502   uplink beam-wander displacement
    Eq. (54), printed p. 503   uplink UNTRACKED on-axis index
    Eqs. (57), (58), printed p. 504   uplink TRACKED on-axis index
    Eq. (23), printed p. 492   the slant Fried parameter
    Eqs. (29), (30), printed p. 493   the isoplanatic angle
The Cn2 profile is the Hufnagel-Valley 5/7 model of Ch. 12, Eq. (3), printed
p. 481, through olb.turbulence.profiles.get_c2n.

WHAT THE TABLES SHOW.
  1. The downlink index falls with the aperture, because the aperture averages
     the irradiance. Eq. (39) is WEAK theory, so read the marked column.
  2. The uplink UNTRACKED index is more than an order of magnitude above the
     tracked one for a small transmit beam. The whole difference is beam wander,
     which is the effect that Ch. 12.6.3 exists to model. Tracking removes the
     WANDER term, not the Rytov term.
  3. The isoplanatic angle shrinks at a low elevation, because the slant path
     crosses more turbulence. The uplink coherence radius rho_0 is the radius AT
     THE SATELLITE (Ch. 12, Eqs. (24) to (27)). The book states below Eq. (27),
     printed p. 492, that it is many times larger than a satellite, so it is
     metres, not centimetres. It is NOT the ground-referred Fried parameter.

TWO LIMITS THE PACKAGE STATES. The slant forms use a plane-parallel atmosphere
with no Earth-curvature correction, and the book limits the weak-fluctuation
slant results to a zenith angle no larger than ZENITH_LIMIT_DEG. The table marks
the rows that break it. And Ch. 12 REFUSES an inner scale or an outer scale on
every slant form: it uses the Kolmogorov spectrum only. The last table prints
the outer-scale PROFILE of Ch. 12, Eq. (68), printed p. 510, and says where it
may be used.

Run from the repo root:
    python -m examples.andrews.slant_paths
'''

from olb import Channel, CircularOrbit, Site, SpaceScenario, Terminal, Transmitter
from olb.turbulence.andrews import (OUTER_SCALE_MODELS, ZENITH_LIMIT_DEG,
                                    beam_params, downlink_scintillation_index,
                                    isoplanatic_angle, outer_scale_profile,
                                    point_ahead_angle, sec_zeta,
                                    uplink_coherence_radius,
                                    uplink_scintillation_index)
from olb.turbulence.profiles import DEFAULT_HS, get_c2n

# --- configuration ----------------------------------------------------------
WAVELENGTH_M = 1550e-9
ALTITUDE_M = 600e3
GROUND_APERTURE_M = 0.7
SPACE_APERTURE_M = 0.08
UPLINK_WAIST_M = 0.10
ORBIT_SPEED_M_S = 7562.0        # circular LEO at 600 km
ELEVATIONS_DEG = (15.0, 20.0, 30.0, 45.0, 60.0, 90.0)


def scenarios():
    '''Build one downlink and one uplink SpaceScenario over the same channel.'''
    channel = Channel(site=Site(cn2_ground=1.7e-14), altitude_m=ALTITUDE_M)
    ground = Terminal(aperture_m=GROUND_APERTURE_M, wavelength_m=WAVELENGTH_M,
                      transmitter=Transmitter(waist_m=UPLINK_WAIST_M))
    space = Terminal(aperture_m=SPACE_APERTURE_M, wavelength_m=WAVELENGTH_M,
                     transmitter=Transmitter(waist_m=0.035))
    return (SpaceScenario(ground=ground, space=space, direction="downlink",
                          channel=channel),
            SpaceScenario(ground=ground, space=space, direction="uplink",
                          channel=channel))


def mark(elevation_deg):
    '''Flag an elevation that breaks the book's own weak-theory zenith limit.'''
    zenith = 90.0 - elevation_deg
    return "*" if zenith > ZENITH_LIMIT_DEG else " "


def print_downlink(scenario, hs, cn2):
    '''Print the downlink point index and the aperture-averaged index.'''
    D = scenario.rx_terminal.aperture_m
    print(f"DOWNLINK, rx aperture D={D*100:.0f} cm, "
          f"lambda={WAVELENGTH_M*1e9:.0f} nm, H-V 5/7")
    print(f"  {'elev':>6} {'sec(zeta)':>10} | {'Eq. (38) point':>15} "
          f"{'Eq. (39) with D':>16} {'ratio':>8}")
    print("  " + "-" * 62)
    for e in ELEVATIONS_DEG:
        point = float(downlink_scintillation_index(hs, cn2, WAVELENGTH_M, e))
        aperture = float(downlink_scintillation_index(hs, cn2, WAVELENGTH_M, e,
                                                      D=D))
        print(f"  {e:>5.0f}{mark(e)} {float(sec_zeta(e)):>10.4f} | "
              f"{point:>15.5f} {aperture:>16.5f} {aperture/point:>8.4f}")
    print(f"  * zenith angle above {ZENITH_LIMIT_DEG:.0f} deg: outside the "
          f"book's weak-theory slant limit.\n")


def print_uplink(scenario, hs, cn2):
    '''Print the tracked and the untracked uplink index, r0 and theta0.'''
    w0 = scenario.tx_terminal.transmitter.waist_m
    print(f"UPLINK, transmit waist W0={w0*100:.0f} cm, orbit "
          f"{ALTITUDE_M/1e3:.0f} km")
    print(f"  {'elev':>6} {'range [km]':>11} | {'tracked':>10} "
          f"{'untracked':>11} {'ratio':>8} | {'rho0 sat [m]':>13} "
          f"{'theta0 [urad]':>14}")
    print("  " + "-" * 82)
    for e in ELEVATIONS_DEG:
        geom = CircularOrbit(ALTITUDE_M, elevation_deg=e)
        beam = beam_params(w0, WAVELENGTH_M, geom.slant_range_m)
        kw = dict(altitude_m=ALTITUDE_M)
        tracked = float(uplink_scintillation_index(hs, cn2, WAVELENGTH_M, e,
                                                   beam, tracked=True, **kw))
        untracked = float(uplink_scintillation_index(hs, cn2, WAVELENGTH_M, e,
                                                     beam, tracked=False, **kw))
        r0 = float(uplink_coherence_radius(hs, cn2, WAVELENGTH_M, e, beam,
                                           **kw))
        theta0 = float(isoplanatic_angle(hs, cn2, WAVELENGTH_M, e))
        print(f"  {e:>5.0f}{mark(e)} {geom.slant_range_m/1e3:>11.1f} | "
              f"{tracked:>10.5f} {untracked:>11.5f} {untracked/tracked:>8.2f} "
              f"| {r0:>13.2f} {theta0*1e6:>14.3f}")
    print(f"  point-ahead angle at {ORBIT_SPEED_M_S:.0f} m/s: "
          f"{float(point_ahead_angle(ORBIT_SPEED_M_S))*1e6:.2f} urad "
          f"(Ch. 12, Eq. (61)).")
    print("  The untracked index carries the beam wander. Tracking removes that "
          "term only.")
    print("  rho0 is measured AT THE SATELLITE, so it is metres. See Conflict "
          "C-02.\n")


def print_outer_scale(heights_m):
    '''Print the outer-scale profile models, and state where they may be used.'''
    print("outer-scale profile L0(h), Ch. 12, Eq. (68), printed p. 510")
    print(f"  {'h [m]':>9} | " + " ".join(f"{n:>10}" for n in OUTER_SCALE_MODELS))
    print("  " + "-" * (11 + 11 * len(OUTER_SCALE_MODELS)))
    for h in heights_m:
        row = " ".join(f"{float(outer_scale_profile(h, n)):>10.4f}"
                       for n in OUTER_SCALE_MODELS)
        print(f"  {h:>9.0f} | {row}")
    print("  REFUSED on every Ch. 12 slant index: Ch. 12, Eq. (15), printed "
          "p. 490, uses the\n  Kolmogorov spectrum only. For a two-scale index "
          "use a single homogeneous path\n  through "
          "andrews.scintillation.weak_two_scale_index.\n")


if __name__ == '__main__':
    downlink, uplink = scenarios()
    hs = DEFAULT_HS
    cn2 = get_c2n(hs, 21.0, downlink.channel.site.cn2_ground)
    print_downlink(downlink, hs, cn2)
    print_uplink(uplink, hs, cn2)
    print_outer_scale([0.0, 1000.0, 5000.0, 8500.0, 20000.0])

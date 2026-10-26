'''
Cross-check of the "Uplink Pointing Budget - Bias & Jitter" widget against olb.

The widget gives the mean irradiance at a LEO satellite for a Gaussian uplink
beam with a static pointing bias and a random pointing jitter. NO turbulence.
This script builds the same link in olb and runs uplink_budget(fidelity=1,
turbulence=False). With turbulence off, the fidelity-1 budget is deterministic:
the geometric Term, the pointing-jitter Term and (reported apart) the slant
extinction Term.

The launch is a 5 cm waist in a 30 cm aperture (no truncation Term, as in the
widget), with the divergence set DIRECTLY
(Transmitter.divergence_rad), so the beam radius at the satellite is
theta * R. The satellite aperture has an area of 1 m^2, so the geometric loss is
the loss "per m^2" and P * 10^(-loss/10) is the irradiance in W/m^2.

olb has no Term for a static pointing BIAS. The script adds it with the
closed form of the widget: the jitter broadens the mean beam to
theta_eff^2 = theta^2 + 4 sigma^2, and the bias b costs
exp(-2 b^2 / theta_eff^2) (Toyoshima et al., JOSA A 19, 567 (2002),
DOI 10.1364/JOSAA.19.000567).

Widget inputs (2026-10-06): P = 12 W, lambda = 1550 nm, R = 1193 km,
sigma = 3 urad per axis, bias = 9.8 urad (the frame of the screenshot),
throughput 0.60, threshold 25 uW/m^2.

RESULT. At theta = 60 urad half-angle the pointing losses agree (0.04 dB jitter,
0.23 dB bias), but the widget irradiance (22.73 uW/m^2, "space loss 78.8 dB")
is wrong: the correct value is about 840 uW/m^2 (space loss 39.06 dB per m^2).

Run from the repository root:
    python -m validation.widget_pointing_xcheck.widget_pointing_xcheck
'''

import numpy as np
from scipy.optimize import brentq

from olb.geometry import CircularOrbit
from olb.scenario import SpaceScenario, Channel
from olb.terminal import Terminal, Transmitter
from olb.links.uplink import uplink_budget

P_W = 12.0
LAM = 1550e-9
RANGE_M = 1193e3
ALTITUDE_M = 500e3        # assumed; only the extinction (elevation) reads it
SIGMA_RAD = 3e-6
BIAS_RAD = 9.8e-6
THROUGHPUT = 0.60
THRESHOLD_W_M2 = 25e-6
WAIST_M = 0.05
WIDGET_60URAD = dict(jitter_db=0.04, bias_db=0.23, space_db=78.8,
                     irradiance_w_m2=22.73e-6)


def db(x):
    '''Loss in dB of a power ratio x (positive for x < 1).'''
    return -10 * np.log10(x)


def run(theta_rad, geom):
    '''Return the budget numbers for one divergence half-angle.'''
    scen = SpaceScenario(
        ground=Terminal(aperture_m=6 * WAIST_M, wavelength_m=LAM,
                        pointing_jitter_rad=SIGMA_RAD,
                        transmitter=Transmitter(waist_m=WAIST_M,
                                                divergence_rad=theta_rad)),
        space=Terminal(aperture_m=2 * np.sqrt(1 / np.pi), wavelength_m=LAM),
        direction="uplink", channel=Channel(altitude_m=ALTITUDE_M))
    terms = {t.name: float(t.mean_db) for t in
             uplink_budget(scen, geom, fidelity=1, turbulence=False).terms}
    theta_eff2 = theta_rad**2 + 4 * SIGMA_RAD**2
    bias_db = db(np.exp(-2 * BIAS_RAD**2 / theta_eff2))
    loss = (terms["geometric spreading"] + db(THROUGHPUT)
            + terms["pointing jitter"] + bias_db)
    return dict(space_db=terms["geometric spreading"],
                jitter_db=terms["pointing jitter"], bias_db=bias_db,
                ext_db=terms["atmospheric extinction"],
                other=sorted(set(terms) - {"geometric spreading",
                                           "pointing jitter",
                                           "atmospheric extinction"}),
                irradiance=P_W * 10 ** (-loss / 10))


if __name__ == '__main__':
    el = brentq(lambda e: CircularOrbit(ALTITUDE_M, e).slant_range_m - RANGE_M,
                1, 90)
    geom = CircularOrbit(ALTITUDE_M, el)
    print(f"R = {geom.slant_range_m / 1e3:.1f} km (h = {ALTITUDE_M / 1e3:.0f} km, "
          f"elevation {el:.2f} deg), P = {P_W} W, sigma = {SIGMA_RAD * 1e6} urad, "
          f"bias = {BIAS_RAD * 1e6:.1f} urad\n")
    print(f"{'theta half':>10} {'w(R) m':>7} {'space dB':>8} {'jitter':>7} "
          f"{'bias':>6} {'I uW/m2':>8} {'margin dB':>9} {'+ext uW/m2':>10}")
    res = {}
    for theta_urad in [60]:
        theta = theta_urad * 1e-6
        r = res[theta_urad] = run(theta, geom)
        assert not r["other"], r["other"]   # untruncated launch, as in the widget
        # Hand check: on-axis Gaussian irradiance 2P/(pi w^2), w = theta * R.
        w = theta * geom.slant_range_m
        hand = (P_W * 2 / (np.pi * w**2) * THROUGHPUT
                * theta**2 / (theta**2 + 4 * SIGMA_RAD**2)
                * np.exp(-2 * BIAS_RAD**2 / (theta**2 + 4 * SIGMA_RAD**2)))
        assert abs(db(r["irradiance"] / hand)) < 0.02, (r["irradiance"], hand)
        print(f"{theta_urad:>7} ur {w:7.1f} {r['space_db']:8.2f} "
              f"{r['jitter_db']:7.3f} {r['bias_db']:6.3f} "
              f"{r['irradiance'] * 1e6:8.1f} "
              f"{-db(r['irradiance'] / THRESHOLD_W_M2):9.2f} "
              f"{r['irradiance'] * 10 ** (-r['ext_db'] / 10) * 1e6:10.1f}")

    r, w = res[60], WIDGET_60URAD
    print("\nwidget at 60 urad half-angle:")
    print(f"  jitter {w['jitter_db']} dB (olb {r['jitter_db']:.3f}), bias "
          f"{w['bias_db']} dB (olb {r['bias_db']:.3f})")
    print(f"  space loss {w['space_db']} dB (olb {r['space_db']:.2f} dB per m^2)")
    print(f"  irradiance {w['irradiance_w_m2'] * 1e6} uW/m^2 (olb "
          f"{r['irradiance'] * 1e6:.1f} uW/m^2, widget low by "
          f"{db(w['irradiance_w_m2'] / r['irradiance']):.1f} dB)")
    assert abs(r["jitter_db"] - w["jitter_db"]) < 0.01
    assert abs(r["bias_db"] - w["bias_db"]) < 0.01
    print("self-check passed")

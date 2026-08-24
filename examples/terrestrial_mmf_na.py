'''
Terrestrial multimode-fibre coupling: the numerical-aperture angular gate.

This script shows the etendue trade that the multimode-fibre numerical aperture
(NA) sets. A multimode fibre is a light bucket: it collects the focal-spot power
that lands inside the hard core disk. A shorter focal length shrinks the spot, so
the spot sits deeper in the core and tolerates more tip-tilt walk-off. BUT the
spot size and the focusing cone are locked by the diffraction invariant

    w_s * NA_optic = lambda / pi ,   NA_optic = (D/2) / f .

So a shorter focal length steepens the cone. When the cone half-angle NA_optic
passes the fibre NA, the fibre does not guide the steep rays, and the coupled
power falls by min(1, (NA/NA_optic)^2). See docs/physics.md section 6c.

Two sweeps:
  1. Fixed optics, sweep the fibre NA. The gate turns on below NA_optic.
  2. Fixed fibre NA, sweep the focal length. A shorter f wins on the walk-off but
     loses to the gate, so a best focal length exists.

Run from the repo root:
    python -m examples.terrestrial_mmf_na
'''

import warnings

import numpy as np

from olb import (TerrestrialScenario, TerrestrialChannel, HorizontalPath,
                 Terminal, Transmitter)
from olb.terminal import MMF
from olb.models.coupling import mmf_coupling_term

warnings.simplefilter("ignore")

# --- configuration ----------------------------------------------------------
WAVELENGTH_M = 1550e-9
PATH_LENGTH_M = 3e3          # horizontal path [m]
CN2 = 1e-15                  # constant Cn2 [m^-2/3]
TX_WAIST_M = 0.02            # transmit Gaussian waist [m]
NEAR_APERTURE_M = 0.30       # transmit (near) telescope [m]
FAR_APERTURE_M = 0.20        # receive (far) telescope [m]
CORE_RADIUS_M = 25e-6        # multimode-fibre core radius [m]
FOCAL_LENGTH_M = 0.50        # receive focusing optic [m] -> NA_optic = 0.20
RX_JITTER_RAD = 5e-6         # receive mechanical jitter, per-axis 1-sigma [rad]


def make_scenario(*, focal_length_m=FOCAL_LENGTH_M, numerical_aperture=None):
    '''Build a terrestrial MMF scenario: near transmits, far receives.'''
    return TerrestrialScenario(
        near=Terminal(aperture_m=NEAR_APERTURE_M, wavelength_m=WAVELENGTH_M,
                      transmitter=Transmitter(waist_m=TX_WAIST_M)),
        far=Terminal(aperture_m=FAR_APERTURE_M, wavelength_m=WAVELENGTH_M,
                     pointing_jitter_rad=RX_JITTER_RAD,
                     detector=MMF(core_radius_m=CORE_RADIUS_M,
                                  focal_length_m=focal_length_m,
                                  numerical_aperture=numerical_aperture)),
        channel=TerrestrialChannel(path_length_m=PATH_LENGTH_M,
                                   attenuation_db_per_km=0.0, cn2=CN2))


def evaluate(**kw):
    '''Return the MMF coupling breakdown of one scenario.'''
    scn = make_scenario(**kw)
    hp = HorizontalPath(scn.channel.path_length_m)
    t = mmf_coupling_term(scn, hp)
    m = t.meta
    # The dB losses split additively: the static loss folds in the NA gate, so the
    # pure spot-in-core loss is the static loss minus the gate. Then
    # static_pure + gate + walkoff = mean.
    return {
        "na_optic": m["na_optic"],
        "spot_um": m["spot_radius_m"] * 1e6,
        "static_db": m["static_loss_db"] - m["na_gate_loss_db"],
        "walkoff_db": m["walkoff_mean_db"],
        "gate_db": m["na_gate_loss_db"],
        "mean_db": float(t.mean_db),
        "q99_db": float(t.quantile_db(0.99)),
    }


def print_header(r):
    '''Print the fixed optics: spot size and the focusing cone NA_optic.'''
    print(f"Terrestrial MMF coupling: L={PATH_LENGTH_M/1e3:.1f} km, "
          f"far D={FAR_APERTURE_M*100:.0f} cm, core={CORE_RADIUS_M*1e6:.0f} um, "
          f"f={FOCAL_LENGTH_M*100:.0f} cm")
    print(f"  focal spot radius={r['spot_um']:.2f} um (<< core: deep in the "
          f"bucket), cone NA_optic={r['na_optic']:.3f}, rx jitter="
          f"{RX_JITTER_RAD*1e6:.1f} urad\n")


def sweep_na(na_values):
    '''Sweep the fibre NA at the fixed optics. The gate turns on below NA_optic.'''
    print("sweep fibre NA (fixed optics):")
    print(f"  {'NA':>6} | {'spot':>7} {'walkoff':>8} {'NA gate':>8} | "
          f"{'mean dB':>8} {'99% dB':>8}")
    print("  " + "-" * 58)
    for na in na_values:
        r = evaluate(numerical_aperture=na)
        print(f"  {na:>6.3f} | {r['static_db']:>7.3f} {r['walkoff_db']:>8.3f} "
              f"{r['gate_db']:>8.3f} | {r['mean_db']:>8.3f} {r['q99_db']:>8.3f}")
    print("  NA >= NA_optic: gate off (0 dB). NA < NA_optic: gate = "
          "-10log10((NA/NA_optic)^2).\n")


def sweep_focal(focal_values, numerical_aperture):
    '''
    Sweep the focal length at a fixed fibre NA.

    A shorter focal length shrinks the spot, so the walk-off loss falls. But it
    steepens the cone, so the NA gate grows. A best focal length balances the two.
    '''
    print(f"sweep focal length (fibre NA={numerical_aperture:.2f}):")
    print(f"  {'f cm':>6} | {'NA_opt':>7} {'spot um':>8} {'walkoff':>8} "
          f"{'NA gate':>8} | {'mean dB':>8}")
    print("  " + "-" * 62)
    best = None
    for f in focal_values:
        r = evaluate(focal_length_m=f, numerical_aperture=numerical_aperture)
        mark = ""
        if best is None or r["mean_db"] < best[1]:
            best = (f, r["mean_db"])
        print(f"  {f*100:>6.1f} | {r['na_optic']:>7.3f} {r['spot_um']:>8.2f} "
              f"{r['walkoff_db']:>8.3f} {r['gate_db']:>8.3f} | {r['mean_db']:>8.3f}")
    print(f"  best focal length: {best[0]*100:.0f} cm "
          f"(mean coupling loss {best[1]:.3f} dB)\n")


if __name__ == '__main__':
    print_header(evaluate())

    # 1. The fibre NA gate. At NA_optic=0.20, a fibre NA below 0.20 gates the cone.
    #    NA=0.20 gives 6.02 dB (factor 0.25); NA>=0.20 gives 0 dB.
    sweep_na([0.10, 0.14, 0.20, 0.30])

    # 2. The etendue trade. A shorter focal length wins on the walk-off but loses
    #    to the NA gate, so a best focal length exists for a common NA=0.20 fibre.
    sweep_focal([0.30, 0.50, 0.80, 1.20, 2.00], numerical_aperture=0.20)

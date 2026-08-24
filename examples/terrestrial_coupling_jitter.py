'''
Terrestrial fibre-coupling efficiency: beam wander, tx jitter, and rx jitter.

This script evaluates the coupling efficiency of a terrestrial (horizontal-path)
single-mode-fibre receiver. It separates the loss into the three pointing
mechanisms that the user asked for, and it keeps the free-space loss apart from
the fibre loss.

Three mechanisms, two groups:

  Free-space (the power that reaches the far aperture):
    - Tx jitter. The transmit beam has a mechanical pointing jitter. The beam
      centre moves at the far aperture, so the collected power falls. This is the
      transmit pointing loss (olb.models.pointing). It is a free-space loss, not
      a fibre-coupling loss.

  Fibre coupling (the fraction of the collected power that enters the fibre):
    - Beam wander. The turbulence moves the arriving beam, so the wavefront tilts
      at the receiver. The tilt moves the focal spot on the fibre tip. This is
      the walk-off, contribution A (Dios beam-wander arrival tilt).
    - Rx jitter. The receive terminal has a mechanical pointing jitter. It also
      moves the focal spot on the fibre tip. This is the walk-off, contribution B.
    - The static mode match plus the higher-order residual set the coupling floor
      (olb.models.coupling.terrestrial_smf_coupling_term, higher-order only,
      because the walk-off owns the tip-tilt).

The walk-off loss is exponential in dB, and its mean is linear in the tilt
variance. The tilt variance is the sum of the beam-wander part and the rx-jitter
part. So the mean walk-off loss splits into the two parts without a second run.

Note: the receive-aperture angle-of-arrival term (the Andrews aperture tilt,
contribution C) is not included yet. See docs/andrews-crosscheck.md batch 2.

Run from the repo root:
    python -m examples.terrestrial_coupling_jitter
'''

import warnings

import numpy as np

from olb import (TerrestrialScenario, TerrestrialChannel, HorizontalPath,
                 Terminal, Transmitter, SMF)
from olb.models.pointing import pointing_loss_term
from olb.models.coupling import terrestrial_smf_coupling_term, smf_walkoff_term

warnings.simplefilter("ignore")

# --- configuration ----------------------------------------------------------
WAVELENGTH_M = 1550e-9
PATH_LENGTH_M = 3e3          # horizontal path [m]
CN2 = 1e-15                  # constant Cn2 [m^-2/3]
TX_WAIST_M = 0.02            # transmit Gaussian waist [m]
TX_DIVERGENCE_RAD = None     # transmit divergence; None = collimated
NEAR_APERTURE_M = 0.30       # transmit (near) telescope [m]
FAR_APERTURE_M = 0.20        # receive (far) telescope [m]
# The receive optic is at optimal focus: the model derives the focal length from
# the mode field radius and the aperture at a=1.12 (SMF optimal_focus=True).
MODE_FIELD_RADIUS_M = 5.2e-6   # single-mode-fibre mode field radius [m] (SMF-28)
TX_JITTER_RAD = 5e-6         # transmit mechanical jitter, per-axis 1-sigma [rad]
RX_JITTER_RAD = 5e-6         # receive mechanical jitter, per-axis 1-sigma [rad]

_K = 20.0 / np.log(10.0)     # dB per (displacement^2 / w_eff^2); see walk-off


def make_scenario(*, cn2=CN2, tx_jitter=TX_JITTER_RAD, rx_jitter=RX_JITTER_RAD):
    '''Build a terrestrial SMF scenario: near transmits, far receives.'''
    return TerrestrialScenario(
        near=Terminal(aperture_m=NEAR_APERTURE_M, wavelength_m=WAVELENGTH_M,
                      pointing_jitter_rad=tx_jitter,
                      transmitter=Transmitter(waist_m=TX_WAIST_M,
                                              divergence_rad=TX_DIVERGENCE_RAD)),
        far=Terminal(aperture_m=FAR_APERTURE_M, wavelength_m=WAVELENGTH_M,
                     pointing_jitter_rad=rx_jitter,
                     detector=SMF(mode_field_radius_m=MODE_FIELD_RADIUS_M,
                                  optimal_focus=True)),
        channel=TerrestrialChannel(path_length_m=PATH_LENGTH_M,
                                   attenuation_db_per_km=0.0, cn2=cn2))


def eta_from_loss_db(loss_db):
    '''Coupling efficiency eta from a positive-dB loss.'''
    return 10.0 ** (-np.asarray(loss_db) / 10.0)


def evaluate(scenario):
    '''
    Return the loss breakdown of one terrestrial SMF scenario.

    Split the walk-off mean into the beam-wander part and the rx-jitter part.
    The walk-off mean is linear in the tilt variance, so the split uses the
    variance meta of one walk-off Term.

    Returns:
        dict
            The mean loss [dB] of each mechanism, and the 99% fade [dB] of the
            stochastic terms.
    '''
    hp = HorizontalPath(scenario.channel.path_length_m)

    # Free-space: the transmit jitter moves the beam at the far aperture.
    tx = pointing_loss_term(scenario, hp)

    # Fibre floor: static mode match + higher-order residual (tip-tilt removed,
    # because the walk-off owns the tip-tilt).
    floor = terrestrial_smf_coupling_term(scenario, hp, drop_tiptilt=True)

    # Fibre walk-off: the received tip-tilt moves the focal spot on the fibre.
    wo = smf_walkoff_term(scenario, hp)
    f = wo.meta["focal_length_m"]
    w_eff = wo.meta["w_eff_m"]
    mean_wander = _K * f ** 2 * wo.meta["sigma2_wander"] / w_eff ** 2
    mean_jitter = _K * f ** 2 * wo.meta["sigma2_jitter"] / w_eff ** 2
    # The two parts sum to the walk-off mean (a linear split).
    assert np.isclose(mean_wander + mean_jitter, wo.mean_db), \
        (mean_wander, mean_jitter, wo.mean_db)

    return {
        "tx_jitter_mean": float(tx.mean_db),
        "tx_jitter_q99": float(tx.quantile_db(0.99)) if tx.stochastic else float(tx.mean_db),
        "floor_mean": float(floor.mean_db),
        "wander_mean": float(mean_wander),
        "rx_jitter_mean": float(mean_jitter),
        "walkoff_mean": float(wo.mean_db),
        "walkoff_q99": float(wo.quantile_db(0.99)),
        "w_eff_um": w_eff * 1e6,
        "spot_um": wo.meta["spot_radius_m"] * 1e6,
        "focal_cm": f * 100.0,
    }


def _row(name, group, loss, extra=""):
    eta = eta_from_loss_db(loss)
    print(f"  {name:<34}{group:<12}{loss:>8.3f}   {eta * 100:>6.2f}   {extra}")


def print_breakdown(scenario):
    '''Print the full loss breakdown at the nominal configuration.'''
    r = evaluate(scenario)
    print(f"Terrestrial SMF coupling: L={PATH_LENGTH_M/1e3:.1f} km, "
          f"Cn2={scenario.channel.cn2:.0e}, tx waist={TX_WAIST_M*100:.1f} cm, "
          f"far D={FAR_APERTURE_M*100:.0f} cm")
    print(f"  focal length={r['focal_cm']:.1f} cm (optimal focus), mode radius="
          f"{MODE_FIELD_RADIUS_M*1e6:.2f} um, spot radius={r['spot_um']:.2f} um, "
          f"w_eff={r['w_eff_um']:.2f} um")
    print(f"  tx jitter={TX_JITTER_RAD*1e6:.1f} urad, "
          f"rx jitter={RX_JITTER_RAD*1e6:.1f} urad\n")
    print(f"  {'mechanism':<34}{'group':<12}{'loss dB':>8}   {'eta %':>6}")
    print("  " + "-" * 62)
    _row("tx mechanical jitter", "free-space", r["tx_jitter_mean"],
         f"99% fade {r['tx_jitter_q99']:.2f} dB")
    _row("mode match + higher-order", "fibre", r["floor_mean"], "mean-only")
    _row("beam wander (turbulence tilt)", "fibre", r["wander_mean"], "walk-off A")
    _row("rx mechanical jitter", "fibre", r["rx_jitter_mean"], "walk-off B")
    print("  " + "-" * 62)
    fibre_total = r["floor_mean"] + r["walkoff_mean"]
    total = fibre_total + r["tx_jitter_mean"]
    _row("FIBRE coupling total", "fibre", fibre_total,
         f"99% walk-off {r['walkoff_q99']:.2f} dB")
    _row("TOTAL (incl. free-space)", "both", total)


def print_sweep(label, values, key, unit, **kw):
    '''Print a one-variable sweep of the loss breakdown.'''
    print(f"\nsweep {label}:")
    print(f"  {label:>10} | {'tx_jit':>7} {'floor':>7} {'wander':>7} "
          f"{'rx_jit':>7} | {'fibre_eta':>9}")
    print("  " + "-" * 60)
    for v in values:
        r = evaluate(make_scenario(**{key: v}, **kw))
        fibre_total = r["floor_mean"] + r["walkoff_mean"]
        print(f"  {v * unit:>10.3g} | {r['tx_jitter_mean']:>7.3f} "
              f"{r['floor_mean']:>7.3f} {r['wander_mean']:>7.3f} "
              f"{r['rx_jitter_mean']:>7.3f} | {eta_from_loss_db(fibre_total) * 100:>8.2f}%")


if __name__ == '__main__':
    print_breakdown(make_scenario())

    # Beam wander grows with the turbulence. The tx and rx jitter do not.
    print_sweep("Cn2", [1e-16, 3e-16, 1e-15, 3e-15, 1e-14], "cn2", 1.0)

    # The rx jitter walk-off grows with the receive jitter. The beam wander and
    # the tx jitter do not change.
    print_sweep("rx_jit_ur", [0.0, 2e-6, 5e-6, 10e-6, 20e-6], "rx_jitter", 1e6)

    # The tx jitter free-space loss grows with the transmit jitter. The fibre
    # terms do not change.
    print_sweep("tx_jit_ur", [0.0, 2e-6, 5e-6, 10e-6, 20e-6], "tx_jitter", 1e6)

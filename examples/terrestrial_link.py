'''
A terrestrial (horizontal-path) link at fidelity 0: aperture vs single-mode fibre.

This example builds the terrestrial budget for a ground-to-ground horizontal
path. It sweeps the link distance and it shows the budget for two receivers: a
power-in-bucket aperture and a single-mode fibre (SMF). Both run at FIDELITY 0.

A terrestrial link uses a TerrestrialScenario: near = the transmit end, far =
the receive end. The channel is a TerrestrialChannel (path length, horizontal
extinction, constant Cn2). The geometry is a HorizontalPath (the range).

Two receiver front ends, two very different budgets:

  Aperture (bucket): the budget is the three deterministic Terms (geometric
  spreading, horizontal extinction, pointing jitter) PLUS the horizontal
  Gaussian-beam scintillation Term. That Term is a real analytic lognormal
  turbulence fade: the point index is the on-axis Gaussian beam-wave index
  sigma2_I(0, L) (Dios et al. 2004), averaged over the receive aperture with the
  Andrews weak Kolmogorov factor (Andrews and Phillips 2005, Ch. 10). Every Term
  has a closed-form fade, so the aperture budget reports an analytic 99% fade
  margin. A larger receive aperture averages more scintillation, so the fade
  shrinks: the aperture-averaging win.

  SMF (single-mode fibre): the budget adds the fibre-coupling loss. That Term is
  FIDELITY 0: it is the MEAN coupling loss from the horizontal Gaussian-beam r0
  and the compensation stack (tip-tilt, AO). It is an effective-r0,
  weak-turbulence approximation (the Noll / Dikmelik-Davidson forms evaluated at
  the Gaussian-beam r0). It models NO fade. So it LOCKS the budget to fidelity 0,
  and the budget then refuses a fade margin: it would be misleading to add the
  pointing fade to a coupling MEAN.

FIDELITY 1 (the statistical coupling fade) is NOT available for a terrestrial
link. The downlink gets it from FAST, but FAST is a far-field / plane-wave-source
engine; a near-field finite Gaussian beam needs a split-step beam-propagation
model, which is not built yet.

Run from the repo root:
    python -m examples.terrestrial_link
'''

import numpy as np

from olb import (TerrestrialScenario, TerrestrialChannel, Site, HorizontalPath,
                 Terminal, Transmitter, Aperture, SMF, TipTilt, AO,
                 terrestrial_budget)

WAVELENGTH_M = 1550e-9
CN2 = 1e-14                # good site: keeps the weak-turbulence regime over the sweep
ATTEN_DB_PER_KM = 0.0       # clear-air horizontal extinction
LAUNCH_DBM = 30.0           # 1 W launch
SENSITIVITY_DBM = -40.0


def _near():
    '''The transmit (near) terminal: a 0.2 m director, a 2 cm waist, 1 W.'''
    return Terminal(aperture_m=0.05, wavelength_m=WAVELENGTH_M,
                    pointing_jitter_rad=5e-6,
                    transmitter=Transmitter(waist_m=0.05/2*0.8, power_dbm=LAUNCH_DBM))


def _budget(path_length_m, detector, compensation=None):
    '''Build the terrestrial budget for one distance and one receiver.'''
    far = Terminal(aperture_m=0.05, wavelength_m=WAVELENGTH_M,
                   detector=detector, compensation=compensation or [])
    channel = TerrestrialChannel(site=Site(), path_length_m=path_length_m,
                                 attenuation_db_per_km=ATTEN_DB_PER_KM, cn2=CN2)
    scenario = TerrestrialScenario(near=_near(), far=far, channel=channel)
    return terrestrial_budget(scenario, HorizontalPath(path_length_m))


def _coupling_term(budget):
    '''Return the SMF coupling Term of a budget, or None.'''
    return next((t for t in budget.terms if t.category == "coupling"), None)


def main():
    # --- 1. Itemised budgets at 3 km ---------------------------------------
    L = 3e3
    ap = _budget(L, Aperture(sensitivity_dbm=SENSITIVITY_DBM))
    smf = _budget(L, SMF(sensitivity_dbm=SENSITIVITY_DBM))

    print("=" * 66)
    print(f"Terrestrial link, {L / 1e3:.0f} km horizontal path, Cn2={CN2:.0e}")
    print("=" * 66)

    print("\nAPERTURE (bucket) receiver -- fidelity 0")
    print("-" * 66)
    print(ap.to_frame()[["name", "category", "mean_db"]].to_string(index=False))
    ap_fade = float(ap.fade_margin_db(0.99))
    scint = next(t for t in ap.terms if t.name == "scintillation")
    # The pointing-only fade (scintillation off) shows how much the turbulence
    # deepens the fade.
    ap_noscint = terrestrial_budget(ap.scenario, HorizontalPath(L),
                                    scintillation=False)
    print(f"\n  mean loss        {float(ap.total_loss_db()):6.2f} dB")
    print(f"  99% analytic fade {ap_fade:6.2f} dB   (fade available: "
          f"{ap.provides_fade})")
    print(f"  99% link margin  {70.0 - ap_fade:6.2f} dB")
    print(f"  turbulence: scintillation IS in this budget "
          f"(sigma2_I={scint.meta['sigma2_I']:.4f}, A={scint.meta['aperture_averaging_factor']:.3f}, "
          f"sigma2_P={scint.meta['sigma2_P']:.4f}).")
    print(f"  it deepens the 99% fade from "
          f"{float(ap_noscint.fade_margin_db(0.99)):.2f} dB (pointing only) to "
          f"{ap_fade:.2f} dB.")

    print("\nSINGLE-MODE FIBRE (SMF) receiver -- fidelity 0 (mean-only)")
    print("-" * 66)
    print(smf.to_frame()[["name", "category", "mean_db"]].to_string(index=False))
    coup = _coupling_term(smf)
    print(f"\n  mean loss        {float(smf.total_loss_db()):6.2f} dB")
    print(f"  fibre coupling   {float(coup.mean_db):6.2f} dB   "
          f"(eta={coup.meta['eta']:.3f}, r0={coup.meta['r0_m'] * 100:.1f} cm)")
    print(f"  mean link margin {70.0 - float(smf.total_loss_db()):6.2f} dB")
    print(f"  fade available:  {smf.provides_fade}   -> the budget refuses a "
          "fade margin:")
    try:
        smf.fade_margin_db(0.99)
    except ValueError as e:
        print(f"    {str(e).split('.')[0]}.")
    print("  why it is locked (always-on caveats on the coupling Term):")
    for v in coup.assumptions.violations:
        print(f"    - {v.split(':')[0]}: {v.split(':', 1)[1].strip().split('.')[0]}.")

    # --- 2. Distance sweep --------------------------------------------------
    # The SMF coupling Term reads a scalar path length, so sweep in a loop (one
    # budget per distance), not with an array geometry.
    print("\n" + "=" * 66)
    print("Distance sweep (mean loss; aperture also gives a 99% fade)")
    print("-" * 66)
    print(f"{'range':>7} | {'aperture mean':>13} {'aperture 99%':>12} | "
          f"{'SMF mean':>9} {'SMF coupling':>12} {'SMF fade':>9}")
    for km in (1.0, 2.0, 3.0, 5.0):
        Lm = km * 1e3
        b_ap = _budget(Lm, Aperture(sensitivity_dbm=SENSITIVITY_DBM))
        b_smf = _budget(Lm, SMF(sensitivity_dbm=SENSITIVITY_DBM))
        c = _coupling_term(b_smf)
        print(f"{km:5.1f}km | {float(b_ap.total_loss_db()):11.2f}dB "
              f"{float(b_ap.fade_margin_db(0.99)):10.2f}dB | "
              f"{float(b_smf.total_loss_db()):7.2f}dB {float(c.mean_db):10.2f}dB "
              f"{'locked':>9}")

    # --- 3. Adaptive optics buys back the coupling --------------------------
    # Tip-tilt removes the first three Zernikes; AO(200) removes many more. Each
    # lowers the residual wavefront, so the fibre couples more of the field.
    print("\n" + "=" * 66)
    print(f"SMF coupling loss at {L / 1e3:.0f} km vs compensation (fidelity 0)")
    print("-" * 66)
    for comp, label in ((None, "none (uncorrected)"), ([TipTilt()], "tip-tilt"),
                        ([TipTilt(), AO(200)], "tip-tilt + AO(200)")):
        c = _coupling_term(_budget(L, SMF(sensitivity_dbm=SENSITIVITY_DBM), comp))
        print(f"  {label:20s} coupling {float(c.mean_db):6.2f} dB   "
              f"eta={c.meta['eta']:.3f}")

    # --- 4. The aperture-averaging win --------------------------------------
    # A power-in-bucket receiver averages the scintillation over its aperture. A
    # larger aperture averages more, so the flux index and the 99% fade shrink.
    print("\n" + "=" * 66)
    print(f"Aperture-averaging win: scintillation fade vs receive aperture "
          f"at {L / 1e3:.0f} km")
    print("-" * 66)
    print(f"{'D (cm)':>7} | {'A':>7} {'sigma2_P':>10} {'99% fade':>10}")
    for D_m in (0.025, 0.05, 0.1, 0.2):
        far = Terminal(aperture_m=D_m, wavelength_m=WAVELENGTH_M,
                       detector=Aperture(sensitivity_dbm=SENSITIVITY_DBM))
        channel = TerrestrialChannel(site=Site(), path_length_m=L,
                                     attenuation_db_per_km=ATTEN_DB_PER_KM, cn2=CN2)
        b = terrestrial_budget(TerrestrialScenario(near=_near(), far=far,
                                                   channel=channel),
                               HorizontalPath(L))
        s = next(t for t in b.terms if t.name == "scintillation")
        print(f"{D_m * 100:7.1f} | {s.meta['aperture_averaging_factor']:7.3f} "
              f"{s.meta['sigma2_P']:10.4f} {float(s.quantile_db(0.99)):9.3f}dB")

    print("\n" + "=" * 66)
    print("Fidelity 1 (the statistical coupling fade) is not built for terrestrial:\n"
          "FAST is a far-field / plane-wave-source engine, but a near-field finite\n"
          "Gaussian beam needs a split-step beam-propagation model. Until then, the\n"
          "SMF terrestrial budget reports the MEAN coupling loss only.")


if __name__ == "__main__":
    main()

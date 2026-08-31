'''
Non-focal-plane (defocused) detector sensing for a terrestrial link.

This script validates the defocus model in olb.models.coupling.terrestrial and
the bidirectional wrapper in olb.links.bidirectional. It is pure analytic, so it
runs WITHOUT aotools. The one optional fidelity-2 cross-check is guarded, so a
missing aotools does not fail the run.

The physics is geometric optics. The detector sits at z = f + dz, with dz the
defocus distance and f the coupling focal length. Two effects follow:

  - the focused spot grows to w_det = gaussz(w_s, dz_eff), w_s =
    lambda*f/(pi*(D/2)) the diffraction spot radius. At large |dz_eff| this tends
    to the geometric blur (D/2)*|dz_eff|/f (Andrews and Phillips 2005, Ch. 4,
    DOI 10.1117/3.626196). The received beam is a DIVERGING Gaussian, so its true
    focus is at z = f + dz_curv with dz_curv = f^2/(R_rx - f) (S. A. Self,
    Appl. Opt. 22, 658 (1983), DOI 10.1364/AO.22.000658). That curvature defocus
    is ALWAYS charged, so the spot argument is dz_eff = dz - dz_curv;
  - the spot centre moves by d_spot = (f+dz)*theta, with theta the arrival tilt
    (ray-optics chief-ray of a thin lens). At focus (dz=0) the lever is f; off
    focus the longer lever arm (f+dz) moves the spot more.

This script checks each of these limits, then it shows the bidirectional wrapper
tie both the transmit divergence and the receive defocus to one dz.

Run from the repository root:
    python validation/defocus/defocus_sensing.py
'''

import os
import sys
import warnings

import numpy as np

# Allow a direct "python validation/defocus/defocus_sensing.py" run: add the
# repository root to the import path (the script directory is two levels down).
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from olb.beam import gaussz
from olb.scenario import TerrestrialScenario, TerrestrialChannel
from olb.geometry import HorizontalPath
from olb.terminal import Terminal, Transmitter, MMF
from olb.models.coupling import terrestrial_mmf_coupling_term
from olb.models.coupling.terrestrial import _spot_offset_sigma, _mmf_focal_length
from olb.links.bidirectional import defocused_terminal, bidirectional_terrestrial

LAM = 1550e-9
D_RX = 0.025                 # far aperture diameter [m]
A_CORE = 25e-6            # multimode-fibre core radius [m]
L_PATH = 1e3             # horizontal path length [m]


def _mmf_scenario(dz, *, jitter, cn2, focal_length):
    '''Build a terrestrial scenario with a defocused MMF far receiver.'''
    detector = MMF(core_radius_m=A_CORE, focal_length_m=focal_length, defocus_m=dz,
                   sensitivity_dbm=-38)
    return TerrestrialScenario(
        near=Terminal(aperture_m=0.3, wavelength_m=LAM,
                      transmitter=Transmitter(waist_m=0.02, power_dbm=30)),
        far=Terminal(aperture_m=D_RX, wavelength_m=LAM, pointing_jitter_rad=jitter,
                     detector=detector),
        channel=TerrestrialChannel(path_length_m=L_PATH, attenuation_db_per_km=0.5,
                                   cn2=cn2))


def _mmf_term(dz, *, jitter=1e-6, cn2=1e-16, focal_length):
    hpath = HorizontalPath(L_PATH)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return terrestrial_mmf_coupling_term(
            _mmf_scenario(dz, jitter=jitter, cn2=cn2,
                          focal_length=focal_length), hpath)


def main():
    # Optimal-focus focal length so the focal spot matches the core (a=1.12).
    f = _mmf_focal_length(MMF(core_radius_m=A_CORE, optimal_focus=True), D_RX, LAM)
    w_s = LAM * f / (np.pi * (D_RX / 2.0))
    print(f"MMF receiver: D={D_RX} m, core radius={A_CORE * 1e6:.0f} um, "
          f"f={f:.2f} m, focal spot w_s={w_s * 1e6:.1f} um")
    print()

    # (a) Sweep the defocus outward. The spot grows, so the mean loss and the 99%
    #     fade both grow with |dz|. A low jitter keeps the static spill in charge.
    print("(a) defocus sweep (low jitter, so the static spill leads):")
    print(f"    {'dz [mm]':>8} {'w_det [um]':>11} {'mean [dB]':>10} {'99% [dB]':>9}")
    dz_grid = np.array([0.0, 0.5, 1.0, 2.0, 4.0, 8.0]) * 1e-3
    means, fades, w_dets = [], [], []
    for dz in dz_grid:
        t = _mmf_term(dz, jitter=1e-6, cn2=1e-16, focal_length=f)
        q99 = t.quantile_db(0.99)
        means.append(t.mean_db)
        fades.append(float(q99))
        w_dets.append(t.meta["spot_radius_detector_m"])
        print(f"    {dz * 1e3:>8.2f} {t.meta['spot_radius_detector_m'] * 1e6:>11.1f} "
              f"{t.mean_db:>10.3f} {float(q99):>9.3f}")
    assert all(np.diff(means) > 0.0), means
    assert all(np.diff(fades) > 0.0), fades
    assert all(np.diff(w_dets) > 0.0), w_dets
    print("    PASS: the spot, the mean loss, and the 99% fade grow with |dz|.")
    print()

    # (b) w_det against gaussz, and against the large-dz geometric blur. The spot
    #     grows from the TRUE focus, so the argument is dz_eff = dz - dz_curv, with
    #     dz_curv the received-curvature focus shift (always charged; see
    #     olb.models.coupling.terrestrial and the resolution section of
    #     validation/defocus/fidelity2_mmf_coupling_gap.md).
    print("(b) spot radius: model vs gaussz vs geometric blur:")
    dz_big = 0.05
    t_big = _mmf_term(dz_big, jitter=1e-9, cn2=1e-18, focal_length=f)
    w_model = t_big.meta["spot_radius_detector_m"]
    dz_eff = t_big.meta["dz_eff_m"]
    w_gauss = gaussz(w_s, dz_eff, LAM)
    w_blur = (D_RX / 2.0) * abs(dz_eff) / f
    print(f"    dz={dz_big * 1e3:.0f} mm (dz_curv={t_big.meta['curvature_defocus_m'] * 1e3:.3f} mm, "
          f"dz_eff={dz_eff * 1e3:.3f} mm):")
    print(f"    model {w_model * 1e6:.1f} um   gaussz {w_gauss * 1e6:.1f} um   "
          f"blur (D/2)dz_eff/f {w_blur * 1e6:.1f} um")
    assert np.isclose(w_model, w_gauss)
    assert np.isclose(w_model, w_blur, rtol=0.02), (w_model, w_blur)
    print("    PASS: w_det = gaussz(w_s, dz_eff), and tends to the geometric blur.")
    print()

    # (c) The chief-ray tilt lever. For a received arrival tilt theta, the spot
    #     moves by d_spot = (f+dz)*theta. Confirm it against _spot_offset_sigma:
    #     the lever is f at focus and grows to (f+dz) off focus.
    print("(c) chief-ray tilt lever: d_spot = (f+dz)*theta:")
    d_tilt_focus = _spot_offset_sigma(f, 0.0, 2.0 * (1e-6) ** 2)
    d_tilt_off = _spot_offset_sigma(f, dz_big, 2.0 * (1e-6) ** 2)
    print(f"    theta=1 urad  ->  focus lever {d_tilt_focus / 1e-6:>8.2f} m   "
          f"dz={dz_big * 1e3:.0f} mm lever {d_tilt_off / 1e-6:>8.2f} m")
    assert np.isclose(d_tilt_focus, f * 1e-6)
    assert np.isclose(d_tilt_off, (f + dz_big) * 1e-6)
    print("    PASS: the tilt lever is f at focus and (f+dz) off focus.")
    print()

    # (d) The bidirectional wrapper: one dz widens the received beam AND defocuses
    #     the coupling. Report the geometric spreading loss increase.
    print("(d) bidirectional wrapper: one dz drives divergence AND defocus:")
    near = Terminal(aperture_m=0.2, wavelength_m=LAM, pointing_jitter_rad=2e-6,
                    transmitter=Transmitter(waist_m=0.05, power_dbm=30),
                    detector=MMF(core_radius_m=A_CORE, optimal_focus=True,
                                 sensitivity_dbm=-38))
    far = Terminal(aperture_m=0.2, wavelength_m=LAM, pointing_jitter_rad=2e-6,
                   transmitter=Transmitter(waist_m=0.05, power_dbm=30),
                   detector=MMF(core_radius_m=A_CORE, optimal_focus=True,
                                sensitivity_dbm=-38))
    chan = TerrestrialChannel(path_length_m=L_PATH, attenuation_db_per_km=0.5,
                              cn2=1e-16)
    geom = HorizontalPath(L_PATH)

    def _geo(budget):
        return next(t for t in budget.terms if t.name == "geometric spreading").mean_db

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        geos = []
        for dz in (0.0, 3e-3, 6e-3):
            b = bidirectional_terrestrial(near, far, chan, geom,
                                          near_defocus_m=dz, far_defocus_m=dz)
            div = b.forward.scenario.tx_terminal.transmitter.divergence_rad
            div_urad = 0.0 if div is None else div * 1e6
            defocus = b.forward.scenario.rx_terminal.detector.defocus_m
            geos.append(_geo(b.forward))
            print(f"    dz={dz * 1e3:>4.1f} mm  tx divergence {div_urad:>6.2f} urad  "
                  f"rx defocus {defocus * 1e3:>4.1f} mm  geo loss {_geo(b.forward):>6.3f} dB")
    assert geos[0] < geos[1] < geos[2], geos
    dt = defocused_terminal(near, 3e-3)
    assert dt.transmitter.divergence_rad is not None and dt.detector.defocus_m == 3e-3
    print("    PASS: a larger dz widens the beam (more geometric loss) and defocuses "
          "the coupling.")
    print()

    # (e) fidelity-2 cross-check of the defocus TREND. The wave-optics layer now
    #     reads detector.defocus_m (the defocus quadratic phase in the field).
    #     Compare in a QUIET regime (weak turbulence, tiny jitter), so the coupling
    #     loss is STATIC-SPILL dominated in both models, and defocus grows the spot
    #     and the loss in both.
    #
    #     IMPORTANT (regime dependence and a model difference): the defocus effect
    #     is NOT always a loss increase. When WALK-OFF dominates (a small spot and
    #     a large tilt), a defocus grows the spot and makes it LESS walk-off
    #     sensitive, so the loss can DROP. And the analytic MMF Term folds the
    #     receive MECHANICAL jitter into the walk-off, but the fidelity-2 mmf_eta
    #     does NOT (it carries the turbulence tilt only). So the two models agree on
    #     the TREND only in the quiet regime here; they differ where the mechanical
    #     jitter walk-off dominates.
    print("(e) fidelity-2 cross-check of the defocus trend (quiet, spill regime):")
    try:
        from olb.models.waveoptics import run_fidelity2
        from olb.links.terrestrial import terrestrial_budget

        def _f2_mmf_loss(dz):
            s = _mmf_scenario(dz, jitter=1e-6, cn2=1e-15,
                              focal_length=f)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                wave = run_fidelity2(s, HorizontalPath(L_PATH), preset="rapid",
                                     n_trials=24, seed=3)
                bud = terrestrial_budget(s, HorizontalPath(L_PATH), fidelity=2,
                                         wave=wave)
            mmf = next(t for t in bud.terms if t.name == "receive coupling (MMF)")
            return mmf.mean_db

        f2_focus = _f2_mmf_loss(0.0)
        f2_defoc = _f2_mmf_loss(1e-3)
        an_focus = _mmf_term(0.0, jitter=1e-6, cn2=1e-15, focal_length=f).mean_db
        an_defoc = _mmf_term(1e-3, jitter=1e-6, cn2=1e-15, focal_length=f).mean_db
        print(f"    fidelity-2 MMF loss:  focus {f2_focus:6.2f} dB  "
              f"defocus 1 mm {f2_defoc:6.2f} dB  ({f2_defoc - f2_focus:+.2f} dB)")
        print(f"    analytic  MMF loss:   focus {an_focus:6.2f} dB  "
              f"defocus 1 mm {an_defoc:6.2f} dB  ({an_defoc - an_focus:+.2f} dB)")
        print(f"    NOTE: the LEVELS differ by ~{f2_focus - an_focus:.1f} dB. BOTH "
              f"models now charge the received-beam WAVEFRONT")
        print(f"          CURVATURE: the received beam is a diverging Gaussian, so "
              f"its true focus sits dz_curv BEYOND the lens")
        print(f"          geometric f, and the analytic Terms read "
              f"dz_eff = defocus_m - dz_curv. What remains is the")
        print(f"          Airy-versus-Gaussian light-bucket spot-shape gap (2-W1), "
              f"about 1 to 1.5 dB. See the resolution")
        print(f"          section of validation/defocus/fidelity2_mmf_coupling_gap.md.")
        assert f2_defoc > f2_focus, (f2_defoc, f2_focus)
        assert an_defoc > an_focus, (an_defoc, an_focus)
        print("    PASS: in the quiet regime both the fidelity-2 and the analytic "
              "MMF coupling loss grow with defocus (direction only).")
    except ImportError:
        print("    aotools not installed; skipping the fidelity-2 cross-check.")
    print()

    print("defocus_sensing PASS")


if __name__ == '__main__':
    main()

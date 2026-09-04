# Non-focal-plane (defocused) detector sensing

This folder validates the defocus model for a terrestrial link. The detector sits
at `z = f + dz`, with `dz` the defocus distance and `f` the coupling focal length.
`dz = 0` puts the detector at focus, so the model reproduces the focal-plane
behaviour exactly.

The received beam is a DIVERGING Gaussian, so the true focus of the coupling
optic sits BEYOND the focal plane, at `z = f + dz_curv` with
`dz_curv = f^2/(R_rx - f)` (S. A. Self, Appl. Opt. 22, 658 (1983),
DOI 10.1364/AO.22.000658). The coupling Terms ALWAYS charge that curvature, at
the actual fibre plane, so the spot growth reads
`dz_eff = defocus_m - dz_curv`. `optimal_focus` stays a focal-LENGTH rule and
never moves the detector; `olb.models.coupling.curvature_focus_shift(scenario)`
gives the `defocus_m` of a tracked (aligned) coupler. See the RESOLUTION appendix
of [fidelity2_mmf_coupling_gap.md](fidelity2_mmf_coupling_gap.md).

The physics is geometric optics:

- **Spot growth.** The focused spot grows to `w_det = gaussz(w_s, dz_eff)`, with
  `w_s = lambda*f/(pi*(D/2))` the diffraction spot radius. At large `dz` this
  tends to the geometric blur `(D/2)*(dz/f)`. Source: Andrews and Phillips,
  2nd ed. (2005), DOI 10.1117/3.626196, Ch. 4.
- **Spot displacement.** The spot centre moves by `d_spot = (f+dz)*theta`, with
  `theta` the arrival tilt (ray-optics chief-ray of a thin lens). At focus
  (`dz = 0`) the lever is `f`; off focus the longer lever arm `(f+dz)` moves the
  spot more.

## Script

| File | Purpose |
| --- | --- |
| [defocus_sensing.py](defocus_sensing.py) | Pure-analytic checks of the defocus model and the bidirectional wrapper. It sweeps `dz` for a multimode-fibre receiver, cross-checks `w_det` against `gaussz` and the geometric blur, confirms the chief-ray tilt lever, and demonstrates the bidirectional wrapper. One optional fidelity-2 cross-check is guarded, so a missing aotools does not fail the run. |

Run it from the repository root:

    python validation/defocus/defocus_sensing.py

## Notes

- A single-mode fibre carries an extra modal quadratic-phase mismatch off focus.
  Its MEAN penalty is now MODELLED, in the coupling Term, by the closed form
  `eta(a, c) = 2 a^2 |(1 - exp(-(a^2 - i c)))/(a^2 - i c)|^2` with
  `c = pi*dz_eff*(D/2)^2/(lambda*f^2)` (Shaklan and Roddier,
  DOI 10.1364/AO.27.002334; Ruilier and Cassaing, JOSA A 18 (2001) 143,
  DOI 10.1364/JOSAA.18.000143). What stays geometric is the walk-off
  DISPLACEMENT response of `terrestrial_smf_walkoff_term`, so that fade is
  OPTIMISTIC off focus and it carries a loud assumptions flag. Use a multimode
  fibre (a light bucket), fidelity 2, or the full modal model above.
- The fidelity-2 wave-optics layer reads `MMF.defocus_m`, so it cross-checks the
  multimode defocus trend. It reads `SMF.defocus_m` too (backlog 2-W2, DONE
  2026-09-04), so the single-mode closed form has a field reference: the
  `olb/waveoptics/smf.py` self-check matches `smf_eta_defocused(a=1.12, c)` to
  four decimals at c = 0, 1, 2 and 4.

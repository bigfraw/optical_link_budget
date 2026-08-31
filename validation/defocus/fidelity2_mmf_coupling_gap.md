# Report: the fidelity-2 vs analytic MMF coupling-loss gap

Author: Claude Code session, 2026-08-31. For review by fable.

## Summary

A terrestrial MMF (light-bucket) receiver reads a much larger coupling loss at
fidelity 2 than the analytic (fidelity-0) MMF Term: about 7 to 8 dB more, at the
same scenario. This report finds the CAUSE. It is NOT primarily the Airy-vs-
Gaussian spot-shape difference (the "2-W1 reference gap"), as first claimed. The
dominant cause is the WAVEFRONT CURVATURE of the received beam: the received beam
is a diverging Gaussian, so its true focus is not at the coupling lens geometric
focal plane f. The fidelity-2 model observes at f; the analytic model assumes a
flat received wavefront (best focus). The gap is a focus-plane modelling
mismatch, not a bug.

## Scenario

Terrestrial horizontal link (`olb.scenario.TerrestrialScenario`):
- wavelength 1550 nm, path length L = 5 km, weak turbulence Cn2 = 1e-15.
- launch: Gaussian waist w0 = 0.02 m, COLLIMATED (no divergence set).
- receive aperture D = 0.2 m.
- detector: MMF, core radius a = 25 um, focal length f = 4.524 m (the
  "optimal_focus" value: f = pi*(D/2)*a/(lambda*1.12), the a=1.12 spot-to-core
  match).

## Observation

| quantity | value |
|---|---|
| fidelity-2 vacuum MMF eta at f | 0.19 to 0.20 -> 7.0 to 7.2 dB loss |
| analytic static MMF eta (Gaussian spot) | 0.92 -> 0.37 dB loss |
| gap | about 6.6 to 7.8 dB |

The gap is present in the VACUUM (no-turbulence) run, so it is NOT a turbulence
effect. Refining the grid (rapid preset, dx_focal 8.64 um, core radius 2.89 focal
pixels -> standard preset, dx_focal 5.96 um, core radius 4.20 pixels) moves the
eta only from 0.189 to 0.201. So it is NOT a coarse focal-sampling artefact.

## Root cause: received-beam wavefront curvature

The received beam is a diverging Gaussian. Its parameters at L = 5 km:
- Rayleigh range zR = pi*w0^2/lambda = 810 m.
- 1/e^2 radius w(L) = 0.125 m (so the D = 0.2 m aperture is moderately filled).
- phase-front radius of curvature R(L) = L*(1 + (zR/L)^2) = 5131 m.

The wavefront SAG across the aperture is (D/2)^2/(2R) = 0.63 waves. That is a
strong DEFOCUS aberration in the pupil. A lens of focal length f focuses a beam
of input curvature R at the plane z' = 1/(1/f + 1/R) = 4.520 m, so the TRUE focus
is 3.99 mm short of the geometric focal plane f = 4.524 m.

The fidelity-2 model puts the detector at f and observes there, so it sees a beam
that is defocused by ~4 mm. The analytic MMF Term uses the flat-wavefront focal
spot w_s = lambda*f/(pi*(D/2)) = 22.3 um, which is the BEST-FOCUS spot of a
collimated (flat) input. It ignores the received-beam curvature.

## Proof

Scan an applied defocus on the SAME fidelity-2 vacuum field (move the observation
plane away from f), and find where the coupling peaks:

```
defocus_m = -8.0 mm -> eta 0.198  (7.03 dB)
defocus_m = -6.0 mm -> eta 0.531  (2.75 dB)
defocus_m = -4.0 mm -> eta 0.719  (1.43 dB)   <-- peak, = the predicted -3.99 mm
defocus_m = -2.0 mm -> eta 0.535  (2.72 dB)
defocus_m =  0.0 mm -> eta 0.201  (6.96 dB)   <-- the geometric focal plane f
defocus_m = +2.0 mm -> eta 0.039  (14.06 dB)
```

The coupling peaks at defocus_m = -4 mm, exactly the curvature-predicted focus
shift. So the loss at the geometric focal plane is an inherent defocus from the
received-beam curvature, NOT the spot shape.

## Decomposition of the ~7 dB

- about 5.5 dB: the received-beam curvature. The true focus is ~4 mm from f, and
  the fidelity-2 detector is at f.
- about 1.0 dB: the residual gap at true best focus (eta 0.72 -> 1.43 dB vs the
  analytic 0.37 dB). THIS is the real Airy-vs-Gaussian light-bucket spot-shape
  difference (the diffraction pattern of a truncated pupil has slow rings that a
  Gaussian model omits). It is the 2-W1 fibre-coupling reference gap, and it is
  small.

## Why the two models differ (neither is a bug)

- The analytic MMF Term (`olb.models.coupling.terrestrial.terrestrial_mmf_coupling_term`)
  assumes a FLAT received wavefront and the best-focus spot w_s = lambda*f/(pi*(D/2)).
  It is OPTIMISTIC: it assumes the coupler is at the true image plane.
- The fidelity-2 model (`olb.waveoptics.turbulence.run` -> `olb.waveoptics.mmf.mmf_coupling_efficiency`)
  uses the ACTUAL curved received field and observes at the geometric focal plane
  f. It is more realistic for a fixed coupler, but it charges the full
  curvature defocus.
- `optimal_focus` derives f from the a=1.12 spot-to-core match assuming a flat
  wavefront. It does NOT place the detector at the true image plane of a curved
  input.

For a SPACE downlink from a distant source the received beam is nearly planar (R
huge), so this gap vanishes. It is specific to a finite-path (terrestrial or
near-field) link, where the received beam carries real curvature.

## Questions for fable (the modelling decision)

1. What is the intended convention for fibre coupling of a curved received beam:
   BEST-FOCUS (the coupler tracks the true image plane) or FIXED geometric f?
2. If best-focus is intended: the fidelity-2 run should observe at the curvature-
   shifted plane. The received-beam curvature is computable (see
   `olb.beam.launch_curvature` and the Gaussian R(z)), so the run could set an
   internal focus offset, or `optimal_focus` could target the true image plane.
   The analytic Term is then already correct (it assumes best focus).
3. If fixed-f is intended: the analytic MMF Term is OPTIMISTIC for a finite-path
   link. It should add the curvature defocus (the received R(z) gives the focus
   shift, which feeds the existing defocus model, `defocus_m` in
   `olb.waveoptics.mmf` and the analytic `w_det` growth).
4. Does the SMF coupling carry the same issue? The analytic SMF coupling also
   assumes a flat wavefront (the a-parameter mode match); a curved received beam
   would lower the modal overlap too. Worth checking with the same method.

## How to reproduce

From the repository root:

```
python validation/defocus/defocus_sensing.py     # section (f) shows the gap
```

The diagnostic that isolates the vacuum coupling and the defocus scan:
build the scenario above, run `olb.models.waveoptics.run_fidelity2`, take the
vacuum receive field `bundle.vacuum.stages[3][1]`, and call
`olb.waveoptics.mmf.mmf_coupling_efficiency(field, D, a, f, defocus_m=...)` over a
range of defocus_m. The peak locates the true focus.

## Files

- `olb/waveoptics/mmf.py` — `focal_intensity`, `mmf_coupling_efficiency`
  (now take `defocus_m`).
- `olb/waveoptics/turbulence/run.py` — the per-trial MMF branch.
- `olb/models/coupling/terrestrial.py` — the analytic MMF Term (flat-wavefront
  spot assumption).
- `olb/links/terrestrial.py` — `_terrestrial_fidelity2_terms` (routes MMF at
  fidelity 2).
- `validation/defocus/defocus_sensing.py` — section (e).

---

# RESOLUTION (2026-08-31)

The text above stays as written, as the record of the investigation. This section
appends what changed. One number in it is WRONG, and the correction is item 1.

## 1. The fidelity-2 `defocus_m` sign was INVERTED

`olb.waveoptics.mmf.focal_intensity` applied the pupil phase
`exp(+i*pi*defocus_m*rho^2/(lambda*f^2))`. In the phase convention of this port a
diverging beam carries `+i*k*r^2/2R` (see `olb.waveoptics.propagators.GForvard`)
and a lens applies `-i*k*r^2/2f` (see `olb.waveoptics.lenses.Lens`), so the plane
`z = f + dz` is the pupil factor `exp(-i*pi*dz*rho^2/(lambda*f^2))`. The knob was
therefore sign-inverted against its own docstring, and every `defocus_m` scan came
out MIRRORED. The sign is now minus.

So the "Proof" scan above reads the right MAGNITUDE but the wrong SIGN. The true
focus of a DIVERGING received beam is **+3.99 mm BEYOND f**, not 3.99 mm short of
it. The `-4 mm` peak in that table was the artifact. A thin lens images a diverging
input beyond its focal plane:

    dz_curv = f^2 / (R_rx - f)   > 0        (S. A. Self, Appl. Opt. 22, 658 (1983),
                                             DOI 10.1364/AO.22.000658)

An independent 1-D Fresnel (Hankel) quadrature of the truncated curved pupil, run
on the PHYSICAL axis with no olb code in the loop, confirms it:

```
detector at z = f + dz (direct quadrature, report scenario)
   dz = -4.0 mm  eta = 0.038  (14.20 dB)
   dz =  0.0 mm  eta = 0.196  ( 7.08 dB)     <-- the geometric focal plane f
   dz = +4.0 mm  eta = 0.709  ( 1.49 dB)     <-- peak, = dz_curv = +3.99 mm
```

The port's `defocus_m` scan now peaks at `+4 mm` too, and it tracks the quadrature
to better than 0.1 dB across the whole +/-8 mm sweep.

Note that the old report's line "A lens of focal length f focuses a beam of input
curvature R at the plane z' = 1/(1/f + 1/R) = 4.520 m" used the CONVERGING-input
sign for R. With the diverging input the imaging equation is
`1/z' = 1/f - 1/R_rx`, giving `z' = 4.528 m = f + 3.99 mm`.

## 2. The decided convention: the curvature defocus is ALWAYS charged

Answering question 3 of the list above: **fixed-f is intended**, and the analytic
Terms were indeed optimistic. The convention now wired into
`olb/models/coupling/terrestrial.py`:

- the received phase-front radius is `R_rx = phase_front_radius(w0, L, ...)`, a new
  pure-physics helper in `olb/beam.py` (Andrews and Phillips 2005, Ch. 4,
  Eqs. (7)-(8), DOI 10.1117/3.626196);
- the true focus sits at `z = f + dz_curv`, `dz_curv = f^2/(R_rx - f)`;
- the detector sits at `z = f + defocus_m`, so its distance from the TRUE focus is
  **`dz_eff = defocus_m - dz_curv`**, and every spot-size and aberration term uses
  `dz_eff`;
- the chief-ray tilt lever `(f+dz)*theta` keeps the PHYSICAL `dz`: the detector
  position, not the focus position, sets it;
- `optimal_focus` keeps its meaning (a focal-LENGTH rule) and never moves the
  detector;
- a user who wants a tracked (aligned-at-true-focus) coupler sets
  `detector.defocus_m = curvature_focus_shift(scenario)`, a new public helper.

For a SPACE link `R_rx` is enormous, so `dz_curv` is about zero and nothing
changes. The downlink and space Terms are untouched.

## 3. The closed forms now wired

Report scenario: 1550 nm, L = 5 km, collimated w0 = 0.02 m, D = 0.2 m, core 25 um,
f = 4.5242 m. Then `R_rx = 5131.5 m` and `dz_curv = +3.992 mm`.

**MMF (light bucket).** The spot growth reads `dz_eff`:
`w_det = gaussz(w_s, dz_eff)`, `eta = 1 - exp(-2*a_core^2/w_det^2)`.

| fibre plane | dz_eff | w_det | eta | loss |
|---|---|---|---|---|
| at f (`defocus_m=0`) | -3.99 mm | 91.0 um | 0.140 | **8.54 dB** |
| at the true focus (`defocus_m=dz_curv`) | 0.00 mm | 22.3 um | 0.919 | 0.37 dB |

The fidelity-2 field reads **7.08 dB** at f, against the analytic 8.54 dB. The
1.5 dB residual is the known 2-W1 Airy-versus-Gaussian light-bucket gap (the
truncated pupil makes an Airy pattern whose slow rings a Gaussian spot model
omits), the same effect the "Decomposition" section above measured as about 1 dB
at true focus. It is NOT chased here.

**SMF (modal overlap).** The flat-wavefront `eta_max(a) = 2*((1-exp(-a^2))/a)^2` is
replaced by the defocus-aberrated closed form, new in
`olb/models/coupling/_common.py`:

    eta(a, c) = 2*a^2 * |(1 - exp(-(a^2 - i*c))) / (a^2 - i*c)|^2,
    c = pi * dz_eff * (D/2)^2 / (lambda * f^2)

It reduces EXACTLY to `eta_max(a)` at `c = 0` (asserted in the module self-check),
and it depends on `|c|` only. Sources: Shaklan and Roddier, DOI
10.1364/AO.27.002334 (the `a` parameter); Ruilier and Cassaing, JOSA A 18, 143
(2001), DOI 10.1364/JOSAA.18.000143 (aberrated single-mode coupling).

At the report scenario with `a = 1.12` and the fibre at f: `c = -3.95 rad`, so
`eta = 0.215` = **6.68 dB**, of which **5.79 dB** is the curvature penalty on the
0.8145 flat-wavefront value. The turbulence residual (extended Marechal /
Dikmelik-Davidson) multiplies this aberrated eta exactly as it multiplied
`eta_max` before, and the `turbulence=False` branch charges it too, because the
curvature is static optics, not turbulence.

Answering question 4 of the list above: yes, the SMF carries the same issue, and
it is now modelled.

## 4. What stays open

- The **~1.5 dB Airy-versus-Gaussian gap (2-W1)** at true focus. The analytic
  Terms use a Gaussian spot; the field uses the real truncated-pupil Airy pattern.
  Owner-gated, not chased.
- The SMF **walk-off** Term (`terrestrial_smf_walkoff_term`) now grows its spot
  over `dz_eff`, but its DISPLACEMENT response stays geometric (a two-Gaussian
  overlap). The MEAN modal penalty is modelled in the coupling Term; the walk-off
  response is not. Its loud flag says exactly that, and still fires on a
  deliberate `defocus_m`.
- An SMF detector with **no resolvable coupling optics** (no focal length) has no
  `a` and no `c`, so it keeps the plain `eta_max` and carries a flag that the
  curvature penalty is NOT modelled there.
- `olb/links/bidirectional.py`: a POSITIVE `dz` (a converging launch) is outside
  the fidelity-0 divergence model, and one `dz` drives launch AND receive, so a
  deliberately diverged monostatic terminal pays `|dz| + dz_curv` of receive
  defocus. Documented in the module and `defocused_terminal` docstrings.

## Files changed in the resolution

- `olb/waveoptics/mmf.py` — the `defocus_m` sign fix, plus a self-check that a
  diverging pupil couples best at a POSITIVE `defocus_m`.
- `olb/beam.py` — new `phase_front_radius`.
- `olb/models/coupling/_common.py` — new `smf_eta_defocused(a, c)`.
- `olb/models/coupling/terrestrial.py` — `_received_curvature`, the public
  `curvature_focus_shift`, `_flag_curvature`, and the `dz_eff` wiring in all
  three terrestrial coupling Terms.
- `olb/links/terrestrial.py` — the fidelity-2 MMF self-check now defocuses AWAY
  from the true focus (a negative `defocus_m`).
- `olb/links/bidirectional.py` — docstring limits only.

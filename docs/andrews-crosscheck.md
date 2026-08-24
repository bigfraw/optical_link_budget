# Andrews cross-check tracker

This file tracks the sections of Andrews and Phillips that we must check our
equations against, or bring into our analytical modules. The user flags the
sections in batches. Each entry keeps a status.

Source book: Andrews and Phillips, *Laser Beam Propagation through Random
Media*, 2nd ed. (SPIE Press, 2005). DOI: 10.1117/3.626196. The equation numbers
are per chapter, so "Ch. 6, Eq. (108)" means equation (108) of Chapter 6. The
end of each chapter collates its equations. Use that collation as the quick
index.

## Status keys

- **flagged** — the user marked the section. We have not checked it yet.
- **checked** — we compared the section to the code. See the note for the
  result.
- **incorporated** — the section is now in a module. See the note for the
  module and the `docs/physics.md` section.

## Policy: strong fluctuation regime

For the strong fluctuation regime, we use numerical simulation, not an analytic
extrapolation of weak-fluctuation theory. So each flagged section that gives a
weak-fluctuation form must also record its regime limit (for example the
Rytov-variance bound). When a scenario passes that limit, the analytic Term is
not valid, and the budget must fall back to a numerical path. See also
`olb/assumptions.py`, which each Term uses to declare its turbulence regime.

## Batch 1 — terrestrial and slant-path Gaussian beam

| Section | Equations | Topic | Maps to | Status |
| --- | --- | --- | --- | --- |
| 6.8.1 | Ch. 6, Eq. (108); Eq. (40), (45) from 6.3.1 | Mean irradiance. Turbulence-induced beam spread for an arbitrary refractive-index structure function. | Beam spread / effective long-term beam radius. `olb/beam.py`, `olb/turbulence/gaussian_fried.py`, `docs/physics.md` §1 and §5e. | flagged |
| 6.6.1 | Ch. 6, Eq. (88) | Beam wander. | Beam-wander variance that folds into the coupled-flux wander term and the Dios off-axis model. `olb/turbulence/uplink_flux.py`, `olb/turbulence/beam_wave_scintillation.py`, `docs/physics.md` §5c and §5d. | flagged |
| 6.8 (general) | slant-path extension | Extension to slant paths for an arbitrary Cn2. | Slant-path generalisation of the Gaussian-beam turbulence Terms. Ties to the CLAUDE.md "Next task" (curvature past the collimated case). | flagged |
| 6.7 | temporal spectra | Temporal spectra of the beam parameters. | The planned temporal-vs-snapshot option. See the temporal-statistics side-step. | flagged |
| 8.2 | scintillation index | Scintillation index for a tracked and an untracked Gaussian beam. Find the restrictions of weak-fluctuation theory for this case. | Scintillation index Terms; the weak/strong regime limit that sets when we switch to a numerical path. `olb/turbulence/plane_wave_scintillation.py`, `olb/turbulence/beam_wave_scintillation.py`, `docs/physics.md` §5b and §5d. | flagged |

## Batch 2 — Gaussian-beam angle of arrival / aperture tip-tilt

| Section | Equations | Topic | Maps to | Status |
| --- | --- | --- | --- | --- |
| (owner to specify) | (owner to specify) | Aperture angle-of-arrival "corrugation" tip-tilt of a Gaussian beam (the classic plane-wave form ~0.182*(D/r0)^(5/3)*(lambda/D)^2). | The second, smaller received tip-tilt term. `olb/turbulence/angle_of_arrival.py` `aperture_arrival_angle_variance` (a stub that raises NotImplementedError). The working received tip-tilt is the beam-wander term only. | **DEFERRED** — owner to specify the explicit Andrews form. Do not guess the coefficient. |

### Notes for batch 2

- The working received tip-tilt used by the coupling Terms is the beam-wander
  arrival tilt, `wander_arrival_angle_variance` (Dios et al. 2004, DOI
  10.1364/AO.43.003866). It reuses the `beam_wander_variance` kernel.
- The aperture angle-of-arrival "corrugation" term is a separate, smaller
  contribution. It is DEFERRED. The stub `aperture_arrival_angle_variance` in
  `olb/turbulence/angle_of_arrival.py` raises NotImplementedError. Fill it with
  the explicit Andrews and Phillips coefficient when the owner specifies it.

### Notes for batch 1

- Section 6.8.1 gives the mean irradiance for an arbitrary structure function.
  Eq. (40) and Eq. (45) come from 6.3.1 and give the base Gaussian-beam wave
  parameters that Eq. (108) uses. So check Eq. (40) and Eq. (45) first.
- Section 6.8 extends the beam-spread and mean-irradiance forms to a slant path
  with an arbitrary Cn2. This is useful for the slant-path generalisation.
- Section 8.2 needs the weak-fluctuation limit written down next to the
  equation. When a terrestrial or slant scenario is in strong turbulence, the
  analytic scintillation index is not valid, and we use a numerical path (see
  the policy above).

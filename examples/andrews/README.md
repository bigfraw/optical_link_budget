# examples/andrews — the Andrews and Phillips foundation layer

Nine runnable scripts. Each one prints a labelled table of numbers that the
package computed. No script invents a number.

Every equation comes from one book:

> L. C. Andrews and R. L. Phillips, *Laser Beam Propagation through Random
> Media*, 2nd ed., SPIE Press (2005). DOI: 10.1117/3.626196

Each module docstring names the chapter, the equation number and the printed
page of the physics that it shows. The package itself is
`olb/turbulence/andrews/`. Its documentation is `docs/physics.md` Section 5h,
and the equation-by-equation audit is `docs/andrews-crosscheck.md`.

## How to run

From the repository root:

    python -m examples.andrews.spectra_and_scales

The `-m` form is the house idiom of `examples/`. It puts the repository root on
the import path, so `import olb` works. Run each script the same way, with its
own module name.

## The scripts

| Script | What it prints |
| --- | --- |
| `spectra_and_scales.py` | The five refractive-index spectra at one Cn2 over four decades of wavenumber, each divided by Kolmogorov, then the inner-scale and the outer-scale sweeps. Ch. 3, Eqs. (18) to (23). |
| `beam_and_coherence.py` | The Gaussian-beam parameters across the input curvature f0 (focused, collimated, diverged), the effective strong-turbulence pair against sigma_R^2, the coherence radius and the Fried parameter for the three wave types, and the CURVATURE-GENERAL Fried parameter. Ch. 4, Ch. 6, Ch. 7. |
| `scintillation_regimes.py` | The scintillation index from weak fluctuation to saturation for the plane wave, the spherical wave and the Gaussian beam, with the two regime boundaries marked; the inner-scale effect on the weak index; the tracked and the untracked Gaussian beam off axis. Ch. 8, Ch. 9, Ch. 12. |
| `distributions_and_fades.py` | Lognormal, gamma-gamma and K at ONE matched scintillation index: the fade depth at each quantile, the probability of fade, the fade rate and the mean fade time, and the numerical proof of the identity Pr(fade) = <n> <t>. Ch. 9, Ch. 11, Ch. 12. |
| `wander_two_routes.py` | Beam wander by the Andrews Ch. 6 route and by the Dios/Belmonte kernel route, on a terrestrial case and an uplink case. Prints the measured ratio 3.5024 = 7.25/2.07 and the long-term-radius residual 1.7512. Conflict C-01. |
| `aperture_averaging.py` | The book's soft-aperture chain against the legacy Churnside fits across the aperture, for the plane wave and the spherical wave, plus the Gaussian-beam all-regime form. Prints the 5 % to 13 % Churnside optimism. Ch. 10, Sec. 10.3. |
| `temporal_statistics.py` | The temporal irradiance spectrum (weak, strong, and strong with an aperture), the quasi-frequency and its band dependence, the Greenwood frequency and tau0 across elevation, and the fade numbers they feed. Ch. 8, Ch. 9, Ch. 11, Ch. 14. |
| `slant_paths.py` | A real `SpaceScenario` swept across elevation: the downlink point and aperture-averaged index, the uplink tracked and untracked index, the uplink coherence radius, the isoplanatic angle, and the outer-scale profile. Ch. 12. |
| `downlink_budget_models.py` | The capstone. One full downlink budget at two elevations with `model="lognormal"`, `"gamma_gamma"` and `"auto"`. Prints the mean loss, the 99 % fade, which model the selector chose, and the aperture-averaging caveat. |

## Wiring status: what is LIVE and what is NOT

The foundation layer holds more physics than the link budgets read today. This
list is honest about that split.

### Live in the package

- **The gamma-gamma downlink Term.** `downlink_scintillation_term(...,
  model="gamma_gamma")` and `model="auto"` compose
  `andrews.scintillation`, `andrews.distributions` and the Term adapter
  `olb/models/fade.py`. The `"auto"` selector switches at the house limit
  sigma_I^2 = 0.25.
- **The delegations.** Each of these keeps its old name and signature, and its
  body now calls the foundation layer: `plane_wave_scintillation.
  coherence_radius`, `plane_wave_scintillation.aperture_averaged_index_andrews`,
  `gaussian_fried.plane_wave_coherence_radius`,
  `gaussian_fried.plane_wave_fried_parameter`,
  `gaussian_fried.output_beam_params`, `gaussian_fried.effective_beam_params`,
  `gaussian_fried.rytov_std`, `ao.plane_wave_fried_parameter_profile`, and
  `angle_of_arrival.aperture_arrival_angle_variance`.
- **The Term adapter.** `olb/models/fade.py` turns any irradiance model into the
  three Term faces. Any new distribution needs no new decibel code.

### Available, but NOT wired into a budget

- **`downlink_budget` has no `model` keyword.** It always asks for the lognormal
  Term. `downlink_budget_models.py` therefore calls it with
  `scintillation=False` and adds the chosen Term itself. The owner must decide
  the default, because the gamma-gamma Term drops the aperture averaging and so
  it moves the budget total by several dB at a low elevation.
- **The tracked uplink index.** `andrews.paths.uplink_scintillation_index(...,
  tracked=True)` gives the scintillation floor that a tracked uplink keeps
  (Ch. 12, Eqs. (57) to (60)). `uplink_budget` does NOT read it. The
  beacon-plus-adaptive-optics uplink budget is still phase-only and flags itself
  `NO SCINTILLATION`. This is olb gap 2.
- **The inner scale and the outer scale.** `andrews.spectra`,
  `andrews.structure` and `andrews.scintillation.weak_two_scale_index` all take
  `l0` and `L0`. No Term passes either one today. Every Ch. 12 slant form
  REFUSES both, because Chapter 12 uses the Kolmogorov spectrum only.
- **The temporal module.** The irradiance temporal spectrum, the quasi-frequency,
  the Greenwood frequency, the coherence time, the fade rate and the mean fade
  time all work, and no budget reads any of them. A caller must set the
  frequency band, because the quasi-frequency has no upper limit of its own.
- **The Andrews beam-wander route.** `andrews.wander` measures only. The uplink
  and terrestrial budgets use the Dios/Belmonte kernel constant 2.07. See
  Conflict C-01 and the adjudicated position in `wander_two_routes.py`.
- **The Andrews aperture-averaging chain.** `andrews.aperture.averaged_index` is
  reachable through `plane_wave_scintillation.aperture_averaged_index_andrews`,
  but the downlink lognormal Term still uses the numerical Airy-filter integral,
  and the terrestrial Term still uses the Churnside fit.
- **The K distribution and the lognormal-Rician PDF.** Both are built. No Term
  uses either.
- **The curvature-general Fried parameter.** It is CLOSED at the physics layer:
  `andrews.beam.beam_params(w0, lambda, z, f0)` feeds
  `andrews.structure.coherence_radius(..., wave="gaussian", beam=...)`, which
  `beam_and_coherence.py` shows. It is NOT closed at the model layer: the
  single-path `gaussian_fried.gaussian_fried_parameter` keeps its collimated
  signature, and the terrestrial fibre-coupling call site passes no `f0`. See
  the CLAUDE.md "Next task" and `docs/physics.md` Section 5e.

### Documented refusals that the scripts show

The package refuses a form that the book does not print, and it names the
citation. The scripts print those refusals in place instead of hiding them:

- the Gaussian-beam Rytov variance for a CONVERGENT beam (Theta0 < 1);
- the Gaussian-beam aperture form when Omega_G < Lambda;
- the K distribution below sigma_I^2 = 1;
- an inner scale or an outer scale on any Ch. 12 slant form, and on the temporal
  spectrum in any regime;
- an aperture-averaged downlink index in the moderate-to-strong regime;
- the annular (centrally obscured) receive aperture, which needs another source.

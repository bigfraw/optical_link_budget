# API reference: results and budget entry points

This page documents the result layer and the four per-link budget functions of
the `olb` package. It gives the exact signatures, keyword arguments, and
defaults from the source.

Sign convention (from `olb.units`): loss is positive dB, gain is negative dB.
The budget adds the value of each Term directly. A `+3` dB Term adds 3 dB of
loss. A `-3` dB Term is a gain and removes 3 dB of loss.

---

## `Term` (`olb/results.py`)

A `Term` is one line of the link budget. It gives three views of the same
contribution. Each view answers a different question.

### The three faces

- `term.mean_db` — the deterministic value or the expected loss `E[loss]`. It is
  a float or a numpy array. A deterministic Term, such as geometric spreading,
  sets only this face.
- `term.quantile_db(p)` — the analytic loss at availability `p`. It returns the
  `p`-quantile of the loss of this Term. It returns `None` for a Term that has
  only Monte Carlo samples and no closed form. That `None` value tells the
  budget to use Monte Carlo, not the analytic sum.
- `term.sample_db(n, rng)` — `n` Monte Carlo draws of the contribution. The
  shape is `(n, *base_shape)`. A deterministic Term broadcasts its mean into the
  samples.

### `stochastic`

`term.stochastic` is a property. It is `True` when the Term has a sampler, and
`False` when it has no sampler. A Term with only a mean is not stochastic.

### `quantile_db(p)` fallback rules

- The Term returns `quantile(p)` when it has a closed-form `quantile`.
- The Term returns `mean_db` when it has no sampler. A deterministic Term is a
  constant across availability.
- The Term returns `None` when it has a sampler but no quantile. This is a
  Monte-Carlo-only Term.

### The mean-only Term

A Term sets `mean_only=True` when it models only the mean of a quantity that
really fluctuates. One example is a fibre-coupling loss from the mean residual
wavefront. A mean-only Term is not the same as a deterministic Term. A
deterministic Term (geometric spreading) has a mean that is its whole
distribution. A mean-only Term has no trustworthy fade. It therefore locks the
whole budget to fidelity 0 (see `provides_fade` below).

---

## `Budget` (`olb/results.py`)

A `Budget` is a list of Terms with the optional top-line values (transmit power,
receive sensitivity).

### `Budget(terms=None, tx_power_dbm=None, rx_sensitivity_dbm=None, scenario=None)`

The constructor stores the Terms and the scenario. It reads the top-line values
from the scenario terminals when the caller gives none. The transmit power comes
from the transmit terminal `Transmitter.power_dbm`. The receive sensitivity
comes from the receive terminal `Detector.sensitivity_dbm`. Both are optional.

### `add(term)`

Append a Term. Return `self`, so the calls chain.

### `total_loss_db()`

Return the sum of the deterministic (mean) contributions in dB. This is the
deterministic view. It ignores the fade.

### `to_frame()`

Return the itemised budget as a pandas DataFrame. There is one row per Term. The
columns are `name`, `category`, `mean_db`, `stochastic`, and `note`. Use this
frame to read the budget line by line.

### `assumptions_frame()`

Return the model assumptions as a pandas DataFrame. There is one row per Term.
The columns are `name`, `beam_type`, `regime`, `spectrum`, `validity`, `ok`, and
`violations`. The `ok` column is `False` when the scenario breaks an assumption.
The `violations` column gives the reasons.

Use `to_frame()` to read the loss values. Use `assumptions_frame()` to check the
model validity. They answer different questions.

### `check(warn=True)`

Find the Terms whose assumptions the scenario breaks. Return a list of
`(term_name, reason)` pairs. Issue a warning for each violation when `warn` is
`True`.

`check()` also tests cross-term spectrum consistency. Every turbulence-bearing
Term models the same atmosphere, so the Terms must agree on the spectrum. A
budget that mixes, for example, a Kolmogorov analytic Term with a von Karman
Term is inconsistent. `check()` then adds one `("budget", reason)` pair. Terms
with no turbulence spectrum (geometric or pointing) do not constrain this.

### `mean_only_terms()`

Return the Terms that model only a mean (fidelity-0, no fade). See the mean-only
Term above.

### `provides_fade`

`provides_fade` is a property. It is `False` when any Term is mean-only. A fade
adds the fade of every Term. A mean-only Term has no fade, so the total tail is
understated. In that case the budget reports the mean only and refuses the fade.
This is the mean-only lock.

### `fade_margin_db(availability)`

Return the analytic loss at the given availability. This is the sum of the
`p`-quantile loss of each Term. This bound ignores the independence of the
Terms. It is an upper bound, because it adds the worst case of every Term
together. Use `monte_carlo()` for the joint distribution.

`fade_margin_db()` raises `ValueError` in two cases:

- Any Term is mean-only (fidelity-0). The budget then has no trustworthy fade,
  so it refuses. Read `total_loss_db()` for the mean, or use a statistical
  (fidelity-1) model.
- Any Term is Monte-Carlo-only (no closed-form quantile). Evaluate those Terms
  with `monte_carlo()`.

### `monte_carlo(n, rng=None, availabilities=(0.99,))`

Draw `n` joint samples of the total loss and summarise the distribution. Monte
Carlo asks every Term for samples. The code samples every Term with the same
`rng` and sums the samples per draw. This keeps the correlations inside a Term
(for example the coupled-flux wander and scintillation). Independent Terms
combine correctly by construction.

Parameters:

- `n` — the number of Monte Carlo draws.
- `rng` — a seeded `numpy.random.Generator` for reproducibility. The default is
  a fresh generator.
- `availabilities` — an iterable of availabilities at which to report the fade
  (loss) level. The default is `(0.99,)`.

Return a dict with these keys:

- `total_loss_samples` — the loss samples, shape `(n, *base_shape)`.
- `mean_loss_db` — the mean loss `E[loss]`.
- `fade_db` — a dict `{availability: loss level}`, or `None` when a mean-only
  Term suppresses the fade.
- `fade_available` — `True` when the budget gives a fade, else `False`.
- `received_dbm` — the received power, or `None` when no transmit power is set.
- `margin_db` — a dict `{availability: margin}`, or `None`.

When a mean-only Term is present, `monte_carlo()` reports the mean, warns, and
suppresses the fade and the fade-based margin. This is the mean-only lock in the
Monte Carlo path.

---

## `Assumptions` (`olb/assumptions.py`)

Each model attaches an `Assumptions` record to its Term. The record states the
regime the model is valid in.

### `Assumptions(beam_type, turbulence_regime, spectrum, validity="", violations=[])`

The fields are `beam_type`, `turbulence_regime`, `spectrum`, `validity` (the
numeric limit in words), and `violations` (the reasons the scenario breaks the
model).

Use the string constants so every Term uses the same words: `BEAM_PLANE_WAVE`,
`BEAM_SPHERICAL_WAVE`, `BEAM_GAUSSIAN`, `BEAM_NA`; `REGIME_WEAK`,
`REGIME_MODERATE`, `REGIME_STRONG`, `REGIME_NA`; `SPECTRUM_KOLMOGOROV`,
`SPECTRUM_VON_KARMAN`, `SPECTRUM_NA`.

### `ok`

`ok` is a property. It is `True` when the scenario breaks no assumption. It is
`False` when the `violations` list is not empty.

### `flag(reason)`

Add one reason that the scenario breaks an assumption. Return `self`.

---

## `units` (`olb/units.py`)

The module re-exports `todB`, `fromdB`, `todBm`, and `fromdBm` from
`my_analysis_modules`. It adds two helper functions for the olb sign convention.

### `loss_db(transmission)`

Return the loss in dB (positive) for a linear power transmission fraction in
`(0, 1]`. A transmission of `1.0` gives 0 dB. The formula is
`-10 * log10(transmission)`.

### `combine_db(*terms_db)`

Sum independent dB contributions (losses positive, gains negative). Return the
total dB. It broadcasts the inputs, then sums along the first axis.

---

## The four budget entry points

Every function returns a `Budget` with the scenario set. Each function assembles
a fixed set of Terms.

| Function | Module | Link | Keyword defaults |
|---|---|---|---|
| `uplink_budget` | `olb/links/uplink.py` | Ground-to-space uplink | `turbulence=True, tau_zenith=None, n_samples=3000, cn2_profile=None` |
| `downlink_budget` | `olb/links/downlink.py` | Space-to-ground downlink | `tau_zenith=None, scintillation=True, turbulence=True, n_samples=2000, smf_fidelity="fast", fast_params=None` |
| `retro_space_budget` | `olb/links/retro_space.py` | Retroreflected ground-to-space | `turbulence=True, tau_zenith=None, n_samples=3000, cn2_profile=None, retro_loss_db=0.0, smf_fidelity="fast", fast_params=None` |
| `terrestrial_budget` | `olb/links/terrestrial.py` | Horizontal ground-to-ground | `scintillation=True, turbulence=True` |

`tau_zenith=None` selects `extinction.DEFAULT_TAU_ZENITH`, which is `0.05`
(the near-IR clear-sky zenith optical depth). All keyword arguments after the
first two are keyword-only.

### `uplink_budget(scenario, geometry, *, turbulence=True, tau_zenith=None, n_samples=3000, cn2_profile=None)`

Assemble the uplink budget. The Terms are the geometric loss, the atmospheric
loss, and, when `turbulence` is `True`, the turbulence physics. The
turbulence Term depends on the pre-compensation source on the scenario
(`SpaceScenario.precompensation`, see `api-terminal-scenario.md`).

- No source (`None`): the uplink is uncorrected. The turbulence Term is the
  coupled-flux Monte Carlo (`uplink_turbulence_term`, beam wander plus
  scintillation). It is a Monte-Carlo-only Term (`quantile=None`), so the budget
  must use `monte_carlo()`. This Term also carries the tracking jitter.
- `DownlinkBeacon` with an `AO` stage: the uplink is pre-compensated. The
  coupled-flux Term is REPLACED by two adding analytic wavefront Terms: the AO
  fitting error (`uplink_fitting_term`, category `fitting`) and the point-ahead
  anisoplanatism (`uplink_point_ahead_term`, category `anisoplanatism`). Both
  are mean-only, so the budget then locks to fidelity 0.

  > **MAJOR LIMITATION — no scintillation.** The two pre-compensation Terms
  > model the PHASE only. The replaced coupled-flux Term carried the
  > scintillation, so the pre-compensated budget MISSES the scintillation and
  > understates the deep fade. Adaptive optics corrects the phase, not the
  > amplitude, so a real corrected uplink still scintillates. Do NOT trust the
  > corrected uplink fade until a scintillation Term is added. Both Terms flag
  > this, so `Budget.check()` warns.
- `DownlinkBeacon` with only a tip-tilt stage: no order above the tilt is
  corrected, so the uplink stays uncorrected (coupled flux).
- `LaserGuideStar`: not modelled yet. `uplink_budget` raises
  `NotImplementedError`.

Other rules:

- `turbulence` — add the turbulence Term when `True`.
- Pointing jitter folds into the coupled-flux turbulence Term. A standalone
  pointing-loss Term is added only when that Term is absent: `turbulence` is
  `False`, or the pre-compensation Terms replace it. The jitter is never lost
  and never double-counted.
- The transmit Gaussian-efficiency Term is opt-in. It fires only when the
  transmit terminal has a `Transmitter` and the launch aperture truncates the
  beam by more than `TX_TRUNCATION_MIN_DB` (`1e-2` dB).
- `cn2_profile=None` builds a default zenith Cn2 profile, so the budget runs
  without the `fast` package.

The budget-building Terms:

- `uplink_turbulence_term(scenario, geometry, n_samples=3000, n_apertures=1,
  hs=None, cn2_profile=None)` builds the uncorrected turbulence Term. It reads
  the transmit waist, the divergence, the wavelength, and the site Cn2. The
  divergence enters the beam broadening and the scintillation index.
- `uplink_point_ahead_term(scenario, geometry, hs=None, cn2_profile=None,
  max_order="auto")` builds the point-ahead anisoplanatism Term. It is the
  decorrelation residual of the corrected Zernike orders across the point-ahead
  angle (see `physics.md` section 5g). `max_order="auto"` reads the AO order
  from the transmit terminal: an `AO(n_modes)` stage sets the highest corrected
  radial order; no AO stage gives the infinite-order upper bound. The phase
  variance becomes a loss with the extended Marechal approximation. The Term is
  mean-only.
- `uplink_fitting_term(scenario, geometry, hs=None, cn2_profile=None)` builds the
  AO fitting-error Term. It is the Noll residual of the uncorrected high orders
  (see `physics.md` section 5f). An empty compensation stack gives the total
  uncorrected phase variance. The Term is mean-only. It models the phase only,
  not the scintillation.

Examples: `examples/uplink_sim.py`, `examples/uplink_divergence.py`,
`examples/build_a_link.py`.

### `downlink_budget(scenario, geometry, *, tau_zenith=None, scintillation=True, turbulence=True, n_samples=2000, smf_fidelity="fast", fast_params=None)`

Assemble the downlink budget. The Terms are the geometric loss, the atmospheric
loss, the pointing loss, and one turbulence effect. The downlink keeps its
standalone pointing Term, unlike the uplink.

- `scintillation` — add the lognormal downlink scintillation Term when `True`.
  Every downlink Term has a closed-form quantile, so the downlink budget
  supports the analytic fade.
- `turbulence` — the master turbulence switch. When `False`, drop every
  turbulence quantity. The budget keeps the deterministic Terms (geometric,
  atmospheric, pointing) and any static coupling loss, but it drops the
  scintillation and the turbulence part of the receive-coupling Term.
- The receive terminal is opt-in. When the receive terminal has a detector, the
  receive-coupling Term owns the receive-side turbulence physics. It replaces
  the standalone scintillation Term. An `Aperture` detector reproduces the plain
  scintillation, so the total is unchanged. An `SMF` detector adds the
  fibre-coupling loss and the coupling fade.
- `smf_fidelity` — the SMF coupling model. `"fast"` (the default) is the
  fidelity-1 true modal overlap; it needs the `fast-aosim` package. `"mean"` is
  the analytic mean-only model with no fade. A mean-only coupling Term locks the
  budget to fidelity 0.
- `n_samples` — the FAST Monte Carlo draws for the SMF fidelity-1 coupling. It
  is ignored for an `Aperture` detector and for `smf_fidelity="mean"`.

`downlink_scintillation_term(scenario, geometry, *, model="lognormal",
aperture_average=True, hs=None, cn2_profile=None)` builds the scintillation Term
on its own. The `model` argument selects the physics through an auto-select
dispatch:

- `"lognormal"` — the analytic weak-fluctuation plane-wave Term. It is the only
  implemented model.
- `"gamma_gamma"` — a reserved slot for the moderate-to-strong regime. It raises
  `NotImplementedError`.
- `"montecarlo"` — a reserved slot for the phase-screen model. It raises
  `NotImplementedError`.
- `"auto"` — the selector layer. It returns the lognormal Term now, and warns
  when `sigma2_I` exceeds the weak-fluctuation limit.

An unknown `model` name raises `ValueError`.

Examples: `examples/downlink_terminal.py`, `examples/build_a_link.py`.

### `retro_space_budget(scenario, geometry, *, turbulence=True, tau_zenith=None, n_samples=3000, cn2_profile=None, retro_loss_db=0.0, smf_fidelity="fast", fast_params=None)`

Assemble the retroreflected ground-to-space budget as a retransmission. The
retroreflector re-transmits the beam. The budget is an up-leg transmission
followed by a down-leg transmission, with the retro aperture as the hinge.

The Terms are the up-leg Terms (geometric, atmospheric, opt-in launch
truncation, and turbulence when `turbulence` is `True`), the down-leg Terms
(geometric, atmospheric, top-hat correction, and the receive-side scintillation
or coupling Term), and the fixed retro-reflection Term. The up-leg Term names
carry an `"uplink "` prefix. The down-leg Term names carry a `"downlink "`
prefix.

- `turbulence` — add the up-leg coupled-flux turbulence Term when `True`. The
  up-leg jitter folds into the turbulence Term, exactly as the uplink does. A
  standalone up-leg pointing Term is added only when `turbulence` is `False`.
- `retro_loss_db` — the fixed loss of the retroreflection in dB. The default is
  `0.0`.
- `smf_fidelity` and `fast_params` — the return-leg receive coupling follows the
  downlink rule. An `SMF` ground detector adds the fibre-coupling loss.

This is the space model only. It assumes a long slant range, a fully diverged
return, and independent turbulence on the two legs. Do not use it for a short
terrestrial retro link.

`retro_budget` is a backward-compatible alias of `retro_space_budget`, kept in
`olb/links/__init__.py` (there is no `retro.py` file). Prefer
`retro_space_budget` in new code.

Example: `examples/retro_link.py`.

### `terrestrial_budget(scenario, geometry, *, scintillation=True, turbulence=True)`

Assemble the terrestrial (horizontal-path) budget. The Terms are the geometric
spreading, the horizontal Beer-Lambert extinction, the pointing jitter, an
opt-in launch truncation, and one receive-side turbulence effect. The
receive-side effect follows the far-terminal detector:

- No detector, or an `Aperture` (bucket) detector: the horizontal Gaussian-beam
  scintillation Term (`terrestrial_scintillation_term`). It is a real analytic
  fade.
- An `SMF` detector: the fidelity-0 mean-only fibre-coupling Term
  (`terrestrial_smf_coupling_term`) replaces the scintillation Term. When the SMF
  sets the coupling optics (`focal_length_m` and `mode_field_radius_m`, or
  `optimal_focus=True`), the budget also adds the receive tip-tilt walk-off fade
  Term (`terrestrial_smf_walkoff_term`). The walk-off then owns the tip-tilt, so the coupling
  Term keeps the higher-order residual only (`drop_tiptilt=True`), and the
  tip-tilt is not counted two times. The coupling Term is mean-only, so it locks
  the budget to fidelity 0 and the budget refuses a fade margin. The walk-off
  Term carries a real fade, but the mean-only lock still holds.
- An `MMF` (light-bucket) detector: the multimode-fibre coupling Term
  (`terrestrial_mmf_coupling_term`) replaces the scintillation Term. It is the encircled
  energy of the focal spot inside the hard core, offset by the received tip-tilt
  (a flat-top acceptance, not a mode overlap). It is not mean-only, so an MMF
  budget keeps its fade.

Flags:

- `scintillation` — add the horizontal Gaussian-beam scintillation Term for an
  `Aperture` or no-detector receiver when `True` (the default). Set it to
  `False` to keep only the deterministic Terms, for example to sweep an array
  path length. The scintillation Term is scalar-only, so it does not broadcast.
- `turbulence` — the master turbulence switch. When `False`, drop every
  turbulence quantity: no scintillation Term, and the fibre-coupling Terms keep
  only their static parts. The SMF coupling Term becomes the static mode-match
  loss, the MMF Term keeps its spot-overfill loss, and the walk-off Term keeps
  only the receive mechanical jitter (the beam-wander tilt drops). The
  deterministic Terms and the transmit pointing jitter stay. So a coupling budget
  with angular jitter still runs, only without turbulence.

The receive-side Terms:

- `terrestrial_scintillation_term(scenario, geometry, *, n_grid=400)` builds the
  horizontal Gaussian-beam scintillation Term. It is a full analytic lognormal
  fade with all three faces. It uses the on-axis Gaussian beam-wave scintillation
  index (Dios et al., Applied Optics 43 (2004) 3866, Eq. 16) and the Andrews
  weak-turbulence aperture-averaging factor. It raises `ValueError` when the near
  terminal has no `Transmitter`.
- `terrestrial_smf_coupling_term(scenario, geometry, *, n_grid=64,
  drop_tiptilt=False, turbulence=True)` builds the mean-only single-mode-fibre
  coupling loss for a horizontal Gaussian beam. `drop_tiptilt=True` removes the
  tip-tilt from the residual, so the walk-off Term can own it. See `physics.md`
  section 6c.
- `terrestrial_smf_walkoff_term(scenario, geometry, *, n_grid=64, turbulence=True)` builds the
  receive tip-tilt walk-off fade. The received tip-tilt (beam wander plus the
  receive mechanical jitter) moves the focal spot on the fibre tip by `f*theta`.
  The fade is exponential in dB. It needs the coupling optics
  (`focal_length_m` and `mode_field_radius_m`, or `optimal_focus=True`), else it
  raises `ValueError`. See `physics.md` section 6c.
- `terrestrial_mmf_coupling_term(scenario, geometry, *, n_grid=64, turbulence=True)` builds the
  multimode-fibre coupling Term: the static spot-in-core overfill loss plus the
  walk-off fade. It needs a focal length (`focal_length_m` or
  `optimal_focus=True`), else it raises `ValueError`. See `physics.md` section 6c.

Examples: `examples/terrestrial_link.py`,
`examples/terrestrial_coupling_jitter.py`.

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
The columns are `name`, `beam_type`, `regime`, `spectrum`, `validity`,
`provenance`, `n_constraints`, `ok`, and `violations`. The `ok` column is `False`
when the scenario breaks an assumption. The `violations` column gives the
reasons. The `provenance` column lists the traced physics source names (empty for
a hand-built or untraced record), and `n_constraints` counts the function-owned
constraints the Term inherited.

Use `to_frame()` to read the loss values. Use `assumptions_frame()` to check the
model validity. They answer different questions.

### `constraints_frame()`

Return the function-owned constraints as a pandas DataFrame, one row per Term and
constraint. The columns are `name` (the Term), `source` (the physics function
that owns the constraint, `<module>.<qualname>`), `kind` (the assumption axis
slug), `statement` (the ASD-STE100 limit), `doi`, and `where` (the printed
citation). This unfolds the `(source, Constraint)` pairs that a traced Term
carries, so a reader sees which function owns each validity limit. A Term with no
traced constraints contributes no row.

### `check(warn=True)`

Find the Terms whose assumptions the scenario breaks. Return a list of
`(term_name, reason)` pairs. Issue a warning for each violation when `warn` is
`True`. A reason that came from a traced physics check or a source-tagged factory
flag carries a `[<source>]` prefix, so the reader sees which function or factory
raised it.

`check()` also tests cross-term spectrum consistency. Every turbulence-bearing
Term models the same atmosphere, so the Terms must agree on the spectrum. A
budget that mixes, for example, a Kolmogorov analytic Term with a von Karman
Term is inconsistent. `check()` then adds one `("budget", reason)` pair. Terms
with no turbulence spectrum (geometric or pointing) do not constrain this.

`check()` also runs an untraced-Term guard. A `turbulence` or `coupling` Term
reads physics functions that now own their assumptions, so its record must carry
traced `provenance`. An EMPTY provenance means the factory did not open the
collection context, so the Term's assumptions are unverified; `check()` reports
it. A legitimately untraced Term (a wave-optics simulation, the external FAST
part) self-declares a `"untraced: ..."` provenance and passes the guard.

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

Each Term carries an `Assumptions` record. The record states the regime the model
is valid in. A factory builds it from the trace of the physics functions that
ran (see below), or by hand for a Term with no traced physics.

### `Assumptions(beam_type, turbulence_regime, spectrum, validity="", violations=[], constraints=[], provenance=[])`

The fields are `beam_type`, `turbulence_regime`, `spectrum`, `validity` (the
numeric limit in words), `violations` (the reasons the scenario breaks the
model), `constraints` (the traced `(source, Constraint)` pairs), and `provenance`
(the traced physics source names). The last two default to empty, so the 23
hand-built records break nothing.

Use the string constants so every Term uses the same words: `BEAM_PLANE_WAVE`,
`BEAM_SPHERICAL_WAVE`, `BEAM_GAUSSIAN`, `BEAM_NA`; `REGIME_WEAK`,
`REGIME_MODERATE`, `REGIME_STRONG`, `REGIME_NA`; `SPECTRUM_KOLMOGOROV`,
`SPECTRUM_VON_KARMAN`, `SPECTRUM_NA`.

### `ok`

`ok` is a property. It is `True` when the scenario breaks no assumption. It is
`False` when the `violations` list is not empty.

### `flag(reason, source=None)`

Add one reason that the scenario breaks an assumption. Return `self`. A
scenario-level fact that the physics never sees (a central obscuration, the
extended-Marechal limit, NO SCINTILLATION) passes `source`
(`source="factory:links.downlink"`), so the violation carries the same
`[source] reason` prefix as a traced check.

## The function-owned assumption mechanism (`olb/assumptions.py`)

A physics function states its own validity through the `@assumes(...)` decorator.
This is the mechanism the Term factories in `olb/turbulence/**`, `olb/links/`, and
`olb/models/` use; see [architecture.md](architecture.md) Section 5 for the
design.

### `Constraint(kind, statement, doi, where="", check=None)`

A frozen record of one validity limit. `kind` is one slug from `KINDS` (an unknown
kind raises); `statement` is one ASD-STE100 sentence; `doi` and `where` cite the
source; `check(args, result) -> Optional[str]` optionally tests the run. A check
returns one reason string when the scenario breaks the limit, or `None`. A check
never warns and never raises.

### `@assumes(*constraints, beam_type=..., turbulence_regime=..., spectrum=...)`

Decorate a public physics function. The decorator stores a `FuncAssumptions`
record on `wrapper.__assumptions__` and, inside a collection context only,
registers the record and runs each constraint check. `module_assumptions(...)`
returns a decorator that carries module-wide defaults; `assumes` is the
no-default form.

### `trace_assumptions()`

A context manager. A Term factory opens `with trace_assumptions() as trace:`
around its physics calls; every decorated function that runs registers to the
`trace`. `trace.merge(beam_type=..., turbulence_regime=..., spectrum=...,
validity=...)` folds the trace into one `Assumptions` record. Outside a context
the decorator adds no work and the numeric output is unchanged.

### `merge_assumptions(*records, validity="")`

Recompose finished `Assumptions` records into one (the retro link folds the
uplink and downlink records) without a trace of its own.

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
| `uplink_budget` | `olb/links/uplink.py` | Ground-to-space uplink | `fidelity=1, turbulence=True, tau_zenith=None, n_samples=3000, cn2_profile=None, wave=None` |
| `downlink_budget` | `olb/links/downlink.py` | Space-to-ground downlink | `fidelity=1, tau_zenith=None, scintillation=True, turbulence=True, n_samples=2000, fast_params=None, scint_model="lognormal", wave=None` |
| `retro_space_budget` | `olb/links/retro_space.py` | Retroreflected ground-to-space | `fidelity=1, turbulence=True, tau_zenith=None, n_samples=3000, cn2_profile=None, retro_loss_db=0.0, fast_params=None` |
| `terrestrial_budget` | `olb/links/terrestrial.py` | Horizontal ground-to-ground | `fidelity=0, scintillation=True, turbulence=True, wave=None` |

`tau_zenith=None` selects `extinction.DEFAULT_TAU_ZENITH`, which is `0.05`
(the near-IR clear-sky zenith optical depth). All keyword arguments after the
first two are keyword-only.

### The `fidelity` ladder

Every budget takes ONE whole-path `fidelity` argument. It sets the model of the
turbulence physics for the whole link. It replaces the old per-component knobs
(`smf_fidelity`, `precomp_fidelity`, the `scint_model="montecarlo"` value, and
`wave_result`).

- `fidelity=0` — analytic. Closed-form Terms. An SMF receiver gets the mean-only
  fibre-coupling Term.
- `fidelity=1` — statistical. FAST modal-overlap coupling and the coupled-flux
  Monte Carlo uplink. It carries a real fade.
- `fidelity=2` — wave optics. The turbulence physics is a field simulation. It
  appears as TWO Terms: a deterministic vacuum-optics Term (the full
  no-turbulence loss from launch to detector, from `propagate_scenario`) and a
  stochastic turbulence Term. Only the analytic extinction and pointing Terms
  stay at fidelity 2. A terrestrial `MMF` receiver adds ONE more Term, the
  light-bucket core coupling (see below). Fidelity 2 needs a precomputed `wave` bundle (a
  `Fidelity2Bundle` from `olb.models.waveoptics.run_fidelity2`); the budget never
  runs the simulation itself.

#### The fidelity-2 coupling faces

When a receiver has a fibre or a light-bucket detector, the fidelity-2
turbulence penalty is a coupling-category face of the wave-optics record. Two
factories in `olb.models.waveoptics` build it from a turbulent run:

- `waveoptics_smf_coupling_term(result, **kwargs)` reduces the per-trial
  single-mode-fibre efficiency `smf_eta`.
- `waveoptics_mmf_coupling_term(result, **kwargs)` reduces the per-trial
  multimode-fibre (light-bucket) efficiency `mmf_eta`.

Each factory reduces the per-trial efficiency to the three Term faces (an
empirical mean, an empirical quantile, and a resampling sampler), category
`coupling`. Each per-trial efficiency is the ABSOLUTE coupling efficiency, so
the loss `-10*log10(eta)` already holds the static floor: `smf_eta` holds the
static mode-match floor, and `mmf_eta` holds the static encircled-energy floor.
No extra floor is added. Both Terms carry a real fade (the turbulent tilt walks
the focused spot off the fixed on-axis core). The receive mechanical jitter is
not in these Terms; it is a separate analytic Term in the budget.

`mmf_eta` also holds the NON-FOCAL-PLANE detector. The runner reads
`MMF.defocus_m`, so the field is focused to the plane `z = f + defocus_m`
(a quadratic pupil phase; see `api-waveoptics.md` section 4a). At the focal plane
(`defocus_m = 0`) this is the plain focal-plane coupling.

`waveoptics_smf_coupling_term` is the fidelity-2 companion of the fidelity-1
FAST `smf_fast_term`. `waveoptics_mmf_coupling_term` is the fidelity-2 companion
of a fidelity-1 FAST MMF Term that does NOT exist: the light bucket has no
analytic and no FAST model, so this Term is the only statistical MMF coupling
model in olb. Both are re-exported from `olb.models.coupling`, so a coupling
Term is discoverable there whatever its fidelity.

Not every fidelity fits every link. Each budget section below gives the exact
mapping and the cases that raise.

### `uplink_budget(scenario, geometry, *, fidelity=1, turbulence=True, tau_zenith=None, n_samples=3000, cn2_profile=None, wave=None)`

Assemble the uplink budget at a chosen `fidelity`. The Terms are the geometric
loss, the atmospheric loss, and, when `turbulence` is `True`, the turbulence
physics. The turbulence route depends on BOTH the `fidelity` and the
pre-compensation source on the scenario (`SpaceScenario.precompensation`, see
`api-terminal-scenario.md`).

UNCORRECTED (no source, or a tip-tilt-only `DownlinkBeacon`):

- `fidelity=0`: raises `ValueError`. There is no analytic mean-only model for an
  uncorrected uplink (beam wander plus scintillation has no closed form).
- `fidelity=1` (the default): the coupled-flux Monte Carlo Term
  (`uplink_turbulence_term`, beam wander plus scintillation). It is a
  Monte-Carlo-only Term (`quantile=None`), so the budget must use
  `monte_carlo()`. This Term also carries the tracking jitter, so there is no
  standalone pointing Term. The coupled-flux fade reads only the waist `w0`, so it
  is obscuration-blind: a set launch `obscuration_ratio` gives the same fade. The
  Term flags this in `budget.check()` and grades the severity by the obscuration
  radius against the waist. The mean loss stays correct through the
  launch-truncation Term. Use fidelity 2 for the fade past a small obscuration
  (see `docs/physics.md` Section 5c).
- `fidelity=2`: two wave-optics Terms. A deterministic vacuum-optics Term (launch
  truncation plus geometric spread plus satellite-aperture capture) and a
  stochastic turbulence Term from the reciprocity overlap `eta_turb` (Shapiro,
  DOI 10.1364/JOSA.61.000492). They replace the geometric, launch-truncation, and
  coupled-flux Terms. The reciprocity overlap holds no jitter, so the standalone
  pointing Term stays. It needs the precomputed `wave` bundle.

PRE-COMPENSATED (`DownlinkBeacon` with an `AO` stage):

- `fidelity=0`: the AO error budget. Two adding analytic wavefront Terms: the AO
  fitting error (`uplink_fitting_term`, category `fitting`) and the point-ahead
  anisoplanatism (`uplink_point_ahead_term`, category `anisoplanatism`). Both are
  mean-only, so the budget then locks to fidelity 0.

  > **Fidelity-0 LIMITATION — no scintillation, no fade.** The two analytic
  > pre-compensation Terms model the PHASE only, and both are mean-only. The
  > replaced coupled-flux Term carried the scintillation, so the fidelity-0
  > pre-compensated budget has no scintillation and no fade of any kind. This is
  > a recorded DECISION (2026-08-27, backlog 0-W1): no trustworthy analytic form
  > exists for the scintillation of a pre-compensated beam. The fidelity-1 FAST
  > Term below is the model of record. Both analytic Terms flag this, and they
  > also flag a residual past the extended-Marechal limit
  > (sigma2 > 1 rad^2, T. S. Ross, DOI 10.1364/AO.48.001812), so
  > `Budget.check()` warns. The budget still returns: the geometric, extinction,
  > and pointing Terms stay exact, and `turbulence=False` gives the
  > geometric-only budget.

- `fidelity=1` (the default, the model of record): ONE Monte-Carlo Term from
  `uplink_fast_term` (category `turbulence`). FAST computes the residual phase of
  the adaptive optics with the point-ahead decorrelation (`DTHETA` from
  `geometry.point_ahead_rad`), plus the uncorrected log-amplitude, and gives the
  flux at the satellite by reciprocity (Shapiro, DOI 10.1364/JOSA.61.000492;
  Farley, DOI 10.1364/OE.458659). The Term carries the scintillation AND a real
  fade, so `fade_margin_db()` works. It is the pure turbulence penalty: the
  standalone pointing Term still fires and carries the mechanical jitter. Needs
  `fast-aosim`; without it the `ImportError` names the fallback.
- `fidelity=2`: raises `ValueError`. The reciprocity screens carry no
  adaptive-optics correction or point-ahead decorrelation, so wave optics does
  not model a pre-compensated uplink. Use `fidelity=1` (FAST).

`LaserGuideStar`: not modelled yet. `uplink_budget` raises `NotImplementedError`.

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
  mean-only. It flags a residual past the Marechal limit (sigma2 > 1 rad^2,
  T. S. Ross, DOI 10.1364/AO.48.001812), where the exponential form overstates
  the loss.
- `uplink_fitting_term(scenario, geometry, hs=None, cn2_profile=None)` builds the
  AO fitting-error Term. It is the Noll residual of the uncorrected high orders
  (see `physics.md` section 5f). An empty compensation stack gives the total
  uncorrected phase variance. The Term is mean-only. It models the phase only,
  not the scintillation, and it carries the same Marechal-limit flag.
- `uplink_fast_term(scenario, geometry, *, hs=None, cn2_profile=None,
  n_samples=1000, fast_params=None)` (in `olb.models.fast`) builds the
  fidelity-1 pre-compensated turbulence Term. FAST overlaps the ground-pupil
  field with the adaptive-optics residual phase (point-ahead decorrelation
  included) and a log-normal log-amplitude; by reciprocity that overlap is the
  uplink on-axis flux at the satellite. The Term is the pure turbulence
  penalty (no static loss) with an empirical mean, quantile, and sampler.
  Scalar elevation only. Needs `fast-aosim`. It requires a `Transmitter`
  waist and an `AO` stage on the ground terminal, and it refuses a
  non-uplink scenario.

Examples: `examples/uplink_sim.py`, `examples/build_a_link.py`,
`validation/coupling_checks/uplink_divergence.py`.

### `downlink_budget(scenario, geometry, *, fidelity=1, tau_zenith=None, scintillation=True, turbulence=True, n_samples=2000, fast_params=None, scint_model="lognormal", wave=None)`

Assemble the downlink budget at a chosen `fidelity`. The Terms are the geometric
loss, the atmospheric loss, the pointing loss, and one turbulence effect. The
downlink keeps its standalone pointing Term, unlike the uplink.

The `fidelity` maps to the receive-side turbulence model:

- `fidelity=0`: an `SMF` detector gets the mean-only analytic fibre-coupling Term
  (`downlink_coupling_term(smf_fidelity="mean")`). An `Aperture` or no detector
  gets the analytic scintillation Term (`scint_model`).
- `fidelity=1` (the default): an `SMF` detector gets the FAST modal-overlap Term
  (`downlink_coupling_term(smf_fidelity="fast")`); it needs the `fast-aosim`
  package. An `Aperture` or no detector uses the SAME analytic scintillation Term
  as fidelity 0. Fidelity 0 and 1 COINCIDE for an aperture: the closed-form
  lognormal or gamma-gamma is the model of record and already carries a fade.
  Only the SMF coupling model changes between the two tiers.
- `fidelity=2`: two wave-optics Terms. A deterministic vacuum-optics Term
  (geometric spread plus aperture capture plus vacuum fibre coupling over the
  full slant range) and a stochastic turbulence Term (the slab penalty). They
  replace the geometric and the scintillation or coupling Terms. Only the
  analytic extinction and pointing Terms stay. It needs the precomputed `wave`
  bundle.

An `MMF` (light-bucket) receive detector is a special case. At fidelity 0 and
fidelity 1 the downlink receive-coupling raises `NotImplementedError`: olb has no
analytic and no FAST MMF coupling model, because the encircled energy of the
focal spot on a fixed core needs the field. At fidelity 2 the downlink budget
routes an `MMF` receiver. It builds THREE Terms: the deterministic vacuum-optics
Term (geometry and truncation only, NO coupling), the aperture-power
scintillation Term (`collected_power`), and the MMF coupling Term
(`waveoptics_mmf_coupling_term`, the absolute core-capture from `mmf_eta`). The
MMF coupling holds the static encircled-energy floor, so NO vacuum coupling
baseline is subtracted (this is the difference from the SMF composite, which
does subtract one).

Other rules:

- `scintillation` — add the analytic scintillation Term for an aperture or
  no-detector receiver at fidelity 0/1 when `True`. Every fidelity-0/1 downlink
  Term has a closed-form quantile, so the downlink budget supports the analytic
  fade.
- `turbulence` — the master turbulence switch for fidelity 0/1. When `False`,
  drop every turbulence quantity. The budget keeps the deterministic Terms
  (geometric, atmospheric, pointing) and any static coupling loss, but it drops
  the scintillation and the turbulence part of the receive-coupling Term.
- The receive terminal is opt-in. When the receive terminal has a detector, the
  receive-coupling Term owns the receive-side turbulence physics. It replaces
  the standalone scintillation Term. An `Aperture` detector reproduces the plain
  scintillation, so the total is unchanged. An `SMF` detector adds the
  fibre-coupling loss and the coupling fade.
- `n_samples` — the FAST Monte Carlo draws for the SMF fidelity-1 coupling. It
  is ignored for an `Aperture` detector and at fidelity 0.
- `scint_model` — the analytic scintillation MODEL for an APERTURE receiver:
  `"lognormal"` (the default), `"gamma_gamma"`, or `"auto"`. It is NOT a fidelity
  axis. It applies at fidelity 0/1 only, and it selects the analytic aperture
  physics (see `downlink_scintillation_term` below).

`downlink_budget` builds its scintillation Term with `scint_model="lognormal"`
and `aperture_average=True`. That is the safe default for a normal elevation. For
a low elevation, build the Term with `model="auto"` and add it to your own
`Budget`, because the gamma-gamma Term drops the aperture averaging and so it
changes the total by several dB.

`downlink_scintillation_term(scenario, geometry, *, model="lognormal",
aperture_average=True, hs=None, cn2_profile=None)` builds the analytic
scintillation Term on its own (fidelity 0/1, aperture). The `model` argument
selects the physics through an auto-select dispatch:

- `"lognormal"` — the analytic weak-fluctuation plane-wave Term.
- `"gamma_gamma"` — the analytic gamma-gamma Term. It holds at every fluctuation
  strength (Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 12,
  Eq. (40), printed p. 497). It models a POINT receiver, because the book gives
  no aperture-averaged downlink index in that regime, and its `Assumptions`
  record flags that. It takes a scalar elevation only; an elevation array raises
  `NotImplementedError`.
- `"auto"` — the selector layer. It reads the point `sigma2_I` and returns the
  lognormal Term below `LOGNORMAL_PDF_LIMIT = 0.25`, or the gamma-gamma Term
  at or above it. For an elevation array that breaks the limit it keeps the
  lognormal Term and warns.

An unknown `model` name raises `ValueError`. The wave-optics downlink is NOT a
`model` name: it is the whole-path `fidelity=2` route of `downlink_budget` (two
Terms). The old `model="montecarlo"` value is gone.

Examples: `examples/downlink_terminal.py`, `examples/build_a_link.py`.

### `retro_space_budget(scenario, geometry, *, fidelity=1, turbulence=True, tau_zenith=None, n_samples=3000, cn2_profile=None, retro_loss_db=0.0, fast_params=None)`

Assemble the retroreflected ground-to-space budget as a retransmission. The
retroreflector re-transmits the beam. The budget is an up-leg transmission
followed by a down-leg transmission, with the retro aperture as the hinge.

The Terms are the up-leg Terms (geometric, atmospheric, opt-in launch
truncation, and turbulence when `turbulence` is `True`), the down-leg Terms
(geometric, atmospheric, top-hat correction, and the receive-side scintillation
or coupling Term), and the fixed retro-reflection Term. The up-leg Term names
carry an `"uplink "` prefix. The down-leg Term names carry a `"downlink "`
prefix.

- `fidelity` — the DOWN-leg receive-coupling model: `1` (the default, FAST modal
  overlap) or `0` (analytic mean-only). The UP-leg turbulence stays the
  coupled-flux Monte Carlo at either value, because there is no analytic
  mean-only uncorrected uplink model. So the up-leg is fidelity 1 regardless.
  `fidelity=2` (wave optics) is NOT supported and raises `ValueError`: the folded
  double pass shares its screens (the two legs are correlated), which needs its
  own design. A `fidelity` other than 0, 1, or 2 also raises.
- `turbulence` — add the up-leg coupled-flux turbulence Term when `True`. The
  up-leg jitter folds into the turbulence Term, exactly as the uplink does. A
  standalone up-leg pointing Term is added only when `turbulence` is `False`.
- `retro_loss_db` — the fixed loss of the retroreflection in dB. The default is
  `0.0`.
- `fast_params` — extra FAST parameters for the fidelity-1 down-leg coupling. The
  return-leg receive coupling follows the downlink rule. An `SMF` ground detector
  adds the fibre-coupling loss.

This is the space model only. It assumes a long slant range, a fully diverged
return, and independent turbulence on the two legs. Do not use it for a short
terrestrial retro link.

`retro_budget` is a backward-compatible alias of `retro_space_budget`, kept in
`olb/links/__init__.py` (there is no `retro.py` file). Prefer
`retro_space_budget` in new code.

Example: `examples/retro_link.py`.

### `terrestrial_budget(scenario, geometry, *, fidelity=0, scintillation=True, turbulence=True, wave=None)`

Assemble the terrestrial (horizontal-path) budget at a chosen `fidelity`.

The `fidelity` maps to the turbulence model:

- `fidelity=0` (the default, analytic). The deterministic Terms (geometric
  spreading, horizontal extinction, pointing jitter) are exact. The receive-side
  turbulence effect is the mean SMF coupling and walk-off, the MMF coupling, or
  the scintillation Term (see the detector rule below).
- `fidelity=1`: raises `ValueError`. It is unavailable for a terrestrial link.
  FAST is a far-field plane-wave-source model; a near-field finite Gaussian beam
  needs the split-step model of fidelity 2 (backlog 1-1).
- `fidelity=2`: the wave-optics Terms. A deterministic vacuum-optics Term (launch
  truncation plus geometric spread plus aperture capture plus vacuum fibre
  coupling) and a stochastic turbulence Term (the fade). They replace the
  geometric, launch-truncation, scintillation, and coupling Terms. Only the
  analytic extinction and pointing Terms stay. An `MMF` receiver gets ONE more
  Term, the light-bucket core coupling
  (`waveoptics_mmf_coupling_term`): it is the ABSOLUTE core capture relative to
  the COLLECTED power, so it does not double-count the aperture capture, and it
  already holds the detector defocus. An
  `Aperture` receiver gets the aperture-power penalty only. It needs the
  precomputed `wave` bundle. A `fidelity` other than 0, 1, or 2 raises
  `ValueError`.

At fidelity 0 the Terms are the geometric spreading, the horizontal
Beer-Lambert extinction, the pointing jitter, an opt-in launch truncation, and
one receive-side turbulence effect. The receive-side effect follows the
far-terminal detector:

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
  energy of the spot inside the hard core, displaced by the received tip-tilt
  (a flat-top acceptance, not a mode overlap). It is not mean-only, so an MMF
  budget keeps its fade.

The received CURVATURE is always charged. A terrestrial received beam is a
diverging Gaussian, so the true focus of the coupling optic is BEYOND its focal
plane, at `z = f + dz_curv` (`dz_curv = f^2/(R_rx - f)`, S. A. Self, Appl. Opt.
22, 658 (1983), DOI 10.1364/AO.22.000658). Every fidelity-0 terrestrial coupling
Term evaluates the detector at `dz_eff = defocus_m - dz_curv`, in both the
turbulent and the `turbulence=False` branch, because the curvature is static
optics. `optimal_focus` stays a focal-LENGTH rule and never moves the detector.
`olb.models.coupling.curvature_focus_shift(scenario)` returns `dz_curv`, so a
tracked (aligned) coupler is `detector.defocus_m = curvature_focus_shift(...)`.
A scenario with no launch beam charges no curvature and flags itself OPTIMISTIC.
See `physics.md` section 6a.

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
  The mean coupling of the Term is the defocus-aberrated closed form
  `smf_eta_defocused(a, c)` (`olb/models/coupling/_common.py`), so the mean
  received-curvature penalty is MODELLED, not only flagged.
- `terrestrial_smf_walkoff_term(scenario, geometry, *, n_grid=64, turbulence=True)` builds the
  receive tip-tilt walk-off fade. The received tip-tilt (beam wander plus the
  receive mechanical jitter) moves the spot on the fibre tip by
  `(f + defocus_m)*theta`. The spot radius is `gaussz(w_s, dz_eff)`. The fade is
  exponential in dB. It needs the coupling optics
  (`focal_length_m` and `mode_field_radius_m`, or `optimal_focus=True`), else it
  raises `ValueError`. The walk-off DISPLACEMENT response is GEOMETRIC ONLY, so
  the Term flags itself when `defocus_m` is not zero. See `physics.md` section 6c.
- `terrestrial_mmf_coupling_term(scenario, geometry, *, n_grid=64, turbulence=True)` builds the
  multimode-fibre coupling Term: the static spot-in-core overfill loss plus the
  walk-off fade. It needs a focal length (`focal_length_m` or
  `optimal_focus=True`), else it raises `ValueError`. See `physics.md` section 6c.
- `curvature_focus_shift(scenario)` (in `olb.models.coupling`) returns the
  received-curvature focus shift `dz_curv` in m of the receive optics. It raises
  `ValueError` when the receive detector is not a fibre, or when the focal length
  cannot be resolved.

Examples: `examples/terrestrial_link.py`,
`validation/coupling_checks/terrestrial_coupling_jitter.py`.

### The bidirectional terrestrial wrapper (`olb/links/bidirectional.py`)

A monostatic terminal uses ONE collimator to transmit and to receive, so ONE
fibre-plane defocus `dz` drives both sides: it diverges the launched beam AND it
moves the detector off the focal plane. This wrapper ties the two to one `dz`.

- `defocused_terminal(terminal, dz_m, *, focal_length_m=None)` returns a NEW
  `Terminal` (the input is not mutated). An `SMF` or `MMF` detector gets
  `defocus_m = dz_m`. A `Transmitter` gets the divergence
  `theta(dz) = sqrt(theta_diff^2 + (W0*|dz|/f^2)^2)`, with
  `theta_diff = lambda/(pi*W0)` and `f` the collimator focal length (Andrews and
  Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 4). `dz = 0` keeps the beam
  collimated. `focal_length_m` None reads `f` from the detector optics; a
  transmit terminal with no resolvable `f` raises `ValueError`.
- `bidirectional_terrestrial(near, far, channel, geometry, *,
  near_defocus_m=0.0, far_defocus_m=0.0, **budget_kwargs)` returns the
  `BidirectionalBudget` namedtuple `(forward, reverse)`: the near->far budget and
  the far->near budget, each from `terrestrial_budget`. The two share the one
  `TerrestrialChannel`. The extra keywords go straight to `terrestrial_budget`.

TWO LIMITS of this fidelity-0 wrapper:

1. Only the DIVERGING side is modelled. `theta(dz)` reads `|dz|` only, and a
   `Transmitter` cannot hold a converging beam. So `dz > 0` (a converging launch)
   is OUTSIDE the model: the wrapper gives it the divergence of the mirror-image
   diverging launch. Use `dz < 0`, or a fidelity-2 field model.
2. One `dz` drives BOTH directions. The received beam is itself a diverging
   Gaussian, so its true focus is already `dz_curv` beyond the focal plane. A
   deliberately diverged monostatic terminal therefore pays `|dz| + dz_curv` of
   receive defocus, and the coupling Terms now CHARGE it. There is no free best
   focus for a monostatic terminal.

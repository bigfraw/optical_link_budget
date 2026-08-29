# Function-owned assumptions with automatic Term inheritance

## Context

olb flags a broken model assumption today only when a Term FACTORY remembers
to (a) hand-build an `Assumptions` record and (b) compute the violation test
itself. The 120 low-level physics functions in `olb/turbulence/**` return bare
floats; their validity limits live only in docstrings. The catalogue
(2026-08-29 exploration) found: 87/120 functions state at least one assumption
in prose, only ~33 enforce anything at runtime, 67 are prose-only, and the
taxonomy has ~21 assumption KINDS — far beyond the 3 fields (beam_type,
regime, spectrum) that `olb/assumptions.py` carries. Zero `Assumptions`
records exist in the physics layer; all 23 sit at factory level. Live
examples of the gap: the terrestrial SMF walk-off Term declares REGIME_WEAK
and never flags; `ZENITH_LIMIT_DEG = 60` is exported and enforced NOWHERE.

GOAL: the assumption belongs to the FUNCTION. A decorator attaches a
machine-readable record and optional runtime checks to each physics function.
A factory opens a collection context; every decorated function that runs
inside registers itself and its violations. The Term inherits the union
automatically — a forgotten dependency becomes impossible.

Decisions fixed with the owner (2026-08-29):
- Mechanism: decorator on the function (`@assumes(...)`).
- Collection: automatic call tracing through `contextvars`.
- First-pass scope: `olb/turbulence/**` + the Term factories in `olb/links/`
  and `olb/models/`. `olb/waveoptics/` waits (numerical-sampling assumptions,
  a different family).

Execution model: the owner is low on Fable usage. This plan contains
PASTE-READY OPUS SUBAGENT PROMPTS. Step one of execution copies this plan
into the repo as `docs/assumptions-refactor-plan.md` so every subagent can
read the spec from disk.

## Design spec

### D1. `Constraint` record (new, in `olb/assumptions.py`)

```python
@dataclass(frozen=True)
class Constraint:
    kind: str                 # one slug from KINDS; __post_init__ raises on unknown
    statement: str            # ASD-STE100 prose: "The Rytov variance obeys sigma_R^2 < 1."
    doi: str                  # "10.1117/3.626196"
    where: str = ""           # "Ch. 5, Eq. (15), printed p. 140"
    check: Optional[Callable] = field(default=None, compare=False)
```

`KINDS` = frozenset of ~21 slugs from the catalogue: `beam-type`, `regime`,
`spectrum`, `pdf-shape`, `launch-curvature`, `receiver`, `obscuration`,
`aperture-order`, `on-axis`, `tracking`, `path-homogeneity`, `path-weight`,
`geometry`, `zenith`, `field-region`, `variance-convention`,
`tilt-convention`, `approximation`, `frozen-flow`, `isoplanatism`,
`not-built`, `conflict`. No enum class, no subclasses. Shared physics
constraints (ZENITH, HOMOGENEOUS_PATH, ...) are module-level `Constraint`
instances in the module that OWNS the physics; `assumptions.py` stays
physics-free.

### D2. Check signature

`check(args: dict, result) -> Optional[str]`. `args` = bound arguments with
defaults applied (Signature cached at decoration time); `result` = the return
value, so a predicate can test a quantity the function computed (the
sigma2_R-in-the-result case). Return None (ok) or one ASD-STE100 reason
string. Vector inputs: test `np.all`, report the worst value. Checks NEVER
warn and NEVER raise; warnings stay factory/user-level. Template:
`andrews.scintillation.rytov_weak` — a regime check returns a reason only on
the "hard" tier; "soft" stays a factory warning (binary violation model,
no severities).

### D3. `@assumes` and module defaults

```python
def module_assumptions(*, beam_type=None, turbulence_regime=None,
                       spectrum=None, constraints=()): ...   # returns a decorator
assumes = module_assumptions()   # the no-default form, exported
```

PER-FUNCTION DECORATION IS THE NORM: `@assumes(CONSTRAINT_A, beam_type=...,
turbulence_regime=...)`. Shared constraints are defined ONCE as module-level
`Constraint` instances and passed explicitly to each function that carries
them — no duplication, no inheritance. `module_assumptions(...)` is OPTIONAL
sugar, allowed ONLY when the statement is true of EVERY public function in
the file (e.g. homogeneous path in andrews/scintillation.py, frozen flow in
temporal.py). The union semantics cannot REMOVE a module default, so one
exception in a file forbids the module form for that constraint — when in
doubt, put it on the function. Function constraints = union(module defaults,
per-function); headline keywords override module defaults. The decorator
stores a frozen `FuncAssumptions(source, beam_type, turbulence_regime,
spectrum, constraints)` on `wrapper.__assumptions__`
(`source = f"{module}.{qualname}"`) — coverage is machine-auditable with no
context. `functools.wraps` preserved.

### D4. Zero overhead without a context

```python
_ACTIVE: ContextVar[Optional[Trace]] = ContextVar("olb_assumptions", default=None)
# wrapper: if _ACTIVE.get() is None: return func(*a, **kw)   # one get, nothing else
# else: run func, trace.record(record), run each constraint check,
#       trace.violate(source, constraint, reason) on a non-None reason
```

Numeric output byte-identical outside a context. `contextvars` is
thread-local: a ThreadPoolExecutor worker does NOT inherit the caller's
context — document in the module docstring (the untraced guard in D7 is the
safety net).

### D5. Collection context and merge

`trace_assumptions()` = contextmanager; sets/resets the ContextVar; yields a
`Trace` with `records: dict[source, FuncAssumptions]` (dedup by source) and
`violations: list[str]` formatted `"[<source>] <reason>"`, deduped exactly.
`trace.merge(beam_type=..., turbulence_regime=..., spectrum=...,
validity="") -> Assumptions`:

1. Headline fields: explicit factory kwarg wins; else the unanimous
   non-`*_NA` traced value; a real conflict joins with "/" AND appends a
   `[merge]` violation telling the factory to state the field.
2. Constraints: union as `(source, Constraint)` pairs, dedup by
   `(source, kind, statement)`.
3. Violations: union of prefixed strings.
4. `validity`: the factory prose, unchanged in role.

Nesting: inner context shadows outer (natural ContextVar behavior;
document). Also provide `merge_assumptions(*records) -> Assumptions` for
`links/retro_space.py`, which recomposes finished Terms.

### D6. Extended `Assumptions` (backward compatible)

Append two defaulted fields to the existing dataclass:
`constraints: list = []` (of `(source, Constraint)`) and
`provenance: list = []` (traced source names). `violations` STAYS a list of
plain strings so `Budget.check()` and `assumptions_frame()` consumers keep
working; provenance is carried in the string prefix. `flag(reason,
source=None)` gains an optional keyword (`source="factory:links.downlink"`)
— this is how scenario-level facts the physics never sees (obscuration_ratio,
Marechal, NO SCINTILLATION) coexist with function-owned entries. `.ok`
unchanged. The 23 existing hand-built records break nothing.

### D7. `olb/results.py` changes

- `assumptions_frame()`: add `provenance` (joined) and `n_constraints`
  columns; keep every existing column.
- New `constraints_frame()`: one row per (term, constraint): name, source,
  kind, statement, doi, where.
- `Budget.check()`: NEW untraced-Term guard — a `turbulence`- or
  `coupling`-category Term whose assumptions record has EMPTY provenance is
  reported ("the factory did not open the collection context"). Legitimately
  untraced Terms self-declare: `provenance=["untraced: wave-optics
  simulation"]` (models/waveoptics.py), `["untraced: fast-aosim"]` (external
  FAST part).

### D8. Constant centralisation

`WEAK_FLUCTUATION_LIMIT = 0.25` is defined twice
(`plane_wave_scintillation.py:58`, `uplink_flux.py:~63`); both are the house
PDF rule canonical as `LOGNORMAL_PDF_LIMIT = 0.25` in
`andrews/scintillation.py`. Both modules re-alias:
`WEAK_FLUCTUATION_LIMIT = LOGNORMAL_PDF_LIMIT` (import name kept for
links/downlink, models/fast, models/waveoptics). Keep the "do not change to
1.0" comment only at the canonical site. Do NOT merge with
`WEAK_REGIME_LIMIT = 1.0` — the PDF-shape axis (sigma2_I) and the regime
axis (sigma2_R) stay separate; the `kind` slugs encode that.

## Work packages

Dependency graph: **WP0 → (WP1a..e, WP2a..c parallel) → (WP3a..d parallel)
→ WP4.** Disjoint file sets within a tier; physics annotation changes
nothing observable until a factory opens a context, so WP1/WP2 land in any
order.

| WP | Files | One line |
|----|-------|----------|
| WP0 | assumptions.py, results.py, plane_wave_scintillation.py + uplink_flux.py (alias only) | The core mechanism + untraced guard + constant |
| WP1a | andrews/spectra.py, andrews/beam.py | Foundation others delegate to |
| WP1b | andrews/scintillation.py, andrews/structure.py | THE worked example; rytov_weak checks |
| WP1c | andrews/aperture.py, andrews/distributions.py | Receiver/obscuration kinds; pdf-shape checks |
| WP1d | andrews/paths.py, andrews/wander.py | ZENITH first enforced; tracking; C-01 conflict tag |
| WP1e | andrews/temporal.py | Frozen flow; no Term consumer, zero risk |
| WP2a | plane_wave_scintillation.py, beam_wave_scintillation.py, gaussian_fried.py | Parent layer, slant path |
| WP2b | uplink_flux.py, coupled_flux.py, angle_of_arrival.py | Dios kernels; the only warns in the layer |
| WP2c | ao.py, anisoplanatism.py, profiles.py | Noll convention; Stone; profiles mostly NA |
| WP3a | links/downlink.py, models/fade.py, models/waveoptics.py | Worked factory example + untraced self-declares |
| WP3b | links/uplink.py, links/retro_space.py | merge_assumptions for retro |
| WP3c | links/terrestrial.py, models/coupling/terrestrial.py, models/coupling/_common.py | Closes the walk-off gap |
| WP3d | models/fast.py, models/coupling/downlink.py | FAST: traced olb part + untraced declare |
| WP4 | docs/, README, CLAUDE.md, coverage audit | /update skill + mechanical coverage count |

No change: `models/{extinction,geometric,pointing,gaussian_efficiency}.py`
(REGIME_NA, no turbulence calls).

Migration rule for every WP3 factory: (a) context around the physics calls
only; (b) DELETE a factory `flag()` only where a traced check now produces an
equivalent violation, and lock the replacement in the self-check (same
trigger still yields `not ok`); (c) KEEP every `warnings.warn` verbatim —
user-visible warnings must not change; (d) scenario-level flags stay factory
`flag(..., source=...)`; (e) `meta` keys unchanged.

## Paste-ready Opus subagent prompts

Preamble to prepend to EVERY prompt below (the shared rules):

```
You work in D:\repos\optical_link_budget (Python package `olb`, optical
link budgets). First read docs/assumptions-refactor-plan.md IN FULL — it is
the spec. Follow the ponytail rule: the laziest solution that works, no
speculative abstraction, reuse the shared kernels. Hard conventions:
ALL docstrings, comments, and commit messages use ASD-STE100 Simplified
Technical English. Every equation and every Constraint cites its source by
DOI (Andrews and Phillips adds chapter, equation number, printed page from
DOI 10.1117/3.626196). Loss is positive dB. Package-relative imports.
Every touched module keeps/extends its `if __name__ == '__main__':`
self-check and `python -m olb.<module>` must pass from the repo root.
Decorate PUBLIC functions only. A pure delegator gets NO decorator (the
traced inner call carries the record); a wrapper that ADDS an assumption
decorates only its delta. Checks never warn and never raise. Do not change
any numeric output: one representative call per module must return the
identical value outside any context.
```

### Prompt WP0 (run first, alone)

```
Implement the core assumptions mechanism, spec sections D1-D8 of
docs/assumptions-refactor-plan.md.

1. olb/assumptions.py: add KINDS, Constraint (frozen dataclass, __post_init__
   raises on an unknown kind), FuncAssumptions, module_assumptions()/assumes,
   the ContextVar wiring with zero-overhead bypass, trace_assumptions(),
   Trace (records dedup by source; violations as "[source] reason" strings,
   exact dedup), Trace.merge(...), merge_assumptions(*records), and extend
   Assumptions with `constraints` and `provenance` (defaulted, backward
   compatible) and flag(reason, source=None).
2. olb/results.py: assumptions_frame() gains provenance and n_constraints
   columns; add constraints_frame(); Budget.check() gains the untraced guard
   for turbulence/coupling Terms with empty provenance (D7).
3. Centralise the 0.25 constant (D8): olb/turbulence/plane_wave_scintillation.py
   and olb/turbulence/uplink_flux.py alias WEAK_FLUCTUATION_LIMIT =
   LOGNORMAL_PDF_LIMIT imported from olb/turbulence/andrews/scintillation.py.
   Keep the house-rule comment at the canonical site only.
4. NEW self-check in assumptions.py proving: (a) identical return value with
   and without a context; (b) in-context registration, a failing check
   appends a source-prefixed violation, no warning emitted; (c) dedup on
   repeated calls; (d) headline-conflict merge produces the [merge] flag;
   (e) unknown kind raises; (f) nested contexts shadow; (g) a
   threading.Thread worker does not see the caller's trace.
5. Extend the results.py self-check: the two new columns, constraints_frame,
   and the untraced guard (fake turbulence Term, empty provenance, check()
   reports it).

Acceptance: python -m olb.assumptions, python -m olb.results,
python -m olb.turbulence.plane_wave_scintillation,
python -m olb.turbulence.uplink_flux all pass; only ONE definition of the
0.25 literal remains (grep).
```

### Prompts WP1a-WP1e and WP2a-WP2c (parallel after WP0)

Common body (fill in the FILES line per package):

```
Annotate the physics functions in FILES with the @assumes decorator, spec
sections D1-D4 of docs/assumptions-refactor-plan.md. The worked example is
section "Worked example" below (WP1b lands it; if it is already merged,
imitate olb/turbulence/andrews/scintillation.py directly).

Method per module: (1) read the module docstring and each function docstring;
every stated assumption becomes a Constraint with its kind slug, one
ASD-STE100 statement, and the DOI + chapter/eq/printed page already present
in the prose; (2) define shared Constraint instances once at module level and
pass them explicitly to each @assumes(...) that carries them; use
module_assumptions(...) defaults ONLY for a statement true of EVERY public
function in the file — one exception forbids the module form, and when in
doubt decorate the function; (3) numeric limits that the prose states become check callables —
reuse the existing module constants (WEAK_REGIME_LIMIT, LOGNORMAL_PDF_LIMIT,
UPLINK_SIGMA2X_LIMIT, ZENITH_LIMIT_DEG, ...), never re-derive a threshold;
(4) an existing raise/NotImplementedError stays a raise — additionally record
it as kind="not-built" or the matching kind; (5) cross-source conflict notes
(C-01 etc.) become kind="conflict" constraints citing both sources.

Extend the module self-check with three blocks: (1) value parity — one
representative call returns the identical float outside any context;
(2) inside trace_assumptions() the expected sources and kinds register;
(3) one deliberately out-of-range call yields a violation containing the
source prefix, and warnings.catch_warnings(record=True) proves the physics
layer emitted nothing.
```

Package-specific FILES lines and notes:

- **WP1a** — `olb/turbulence/andrews/spectra.py`, `olb/turbulence/andrews/beam.py`.
  Spectra: one spectrum kind per builder; keep the existing `_refuse`/`_need`
  raises. Beam: launch-curvature (f0 sign selects the branch) and paraxial
  approximation kinds.
- **WP1b** — `olb/turbulence/andrews/scintillation.py`, `olb/turbulence/andrews/structure.py`.
  THIS PACKAGE LANDS THE WORKED EXAMPLE (see plan section below). Build the
  regime check on rytov_weak ("hard" tier only). The Gaussian second weak
  condition (sigma_R^2 * Lambda^(5/6) < 1, today "the caller must test")
  becomes a real check. structure.py: tilt-convention (G-tilt not Noll),
  the sqrt(L/k) << D condition of angle_of_arrival_variance (today
  explicitly ungated) becomes a check, coherence_radius outer-scale
  restriction.
- **WP1c** — `olb/turbulence/andrews/aperture.py`, `olb/turbulence/andrews/distributions.py`.
  Aperture: receiver kind (soft Gaussian aperture, NOT the Airy filter —
  conflict C-06), obscuration (no annular model, gap G-108), aperture-order
  (Omega_G >= Lambda stays a raise + constraint). Distributions: pdf-shape
  check on the lognormal faces against LOGNORMAL_PDF_LIMIT (needs sigma2 in
  args); k_params keeps its sigma_I^2 > 1 raise + a regime constraint;
  gamma-gamma records all-strength validity.
- **WP1d** — `olb/turbulence/andrews/paths.py`, `olb/turbulence/andrews/wander.py`.
  paths: ZENITH_LIMIT_DEG becomes an enforced check (kind="zenith") on every
  slant function — its FIRST enforcement anywhere; geometry (plane-parallel,
  no Earth curvature) module default; downlink strong branch point-receiver
  constraint; uplink tracked=True carries a constraint stating it is
  OPTIMISTIC and does not model a pre-compensated uplink (GAP 2). wander:
  variance-convention (radial, never /2), approximation (geometrical
  optics, lambda drops out), weak-fluctuation, C-01 conflict tag (Andrews
  7.25 vs Dios/Belmonte 2.07), the silent negative-clip in
  short_term_beam_radius becomes a check that reports the broken weak limit.
- **WP1e** — `olb/turbulence/andrews/temporal.py`.
  frozen-flow kind on everything; quasi_frequency records the band-dependence
  as a constraint; keep the three existing raises; no Term consumer exists,
  so this package cannot break a budget.
- **WP2a** — `olb/turbulence/plane_wave_scintillation.py`, `olb/turbulence/beam_wave_scintillation.py`, `olb/turbulence/gaussian_fried.py`.
  plane_wave_scintillation_index: its RESULT is sigma2_R — carry the regime
  check on the result; the Airy aperture filter gets receiver + obscuration
  constraints (uniform circular, unobscured); Churnside fits get their own
  regime/inner-scale constraints (explicitly not the Andrews fit).
  beam_wave: collimated-launch constraint, the sigma_chi^2 ~ 0.6 Dios
  reliability bound as a check. gaussian_fried: collimated (Theta0=1)
  constraints; spherical_wave_fried_parameter carries a path-homogeneity
  constraint ("do not use for an uplink"); the profile form keeps its
  explicit weak-turbulence block as constraints.
- **WP2b** — `olb/turbulence/uplink_flux.py`, `olb/turbulence/coupled_flux.py`, `olb/turbulence/angle_of_arrival.py`.
  uplink_flux._flux_result: keep BOTH existing warnings verbatim; add the
  traced check at the hard tier (sigma_x^2 >= the Dios bound), plus
  obscuration constraint (pure unclipped Gaussian launch — the index is
  blind to obscuration). coupled_flux: path-weight (transmitter-referred,
  "do not flip"), variance-convention (radial), C-01 conflict tag on
  beam_wander_variance. angle_of_arrival: tilt-convention (G-tilt), the
  ungated sqrt(L/k) << D condition becomes a check.
- **WP2c** — `olb/turbulence/ao.py`, `olb/turbulence/anisoplanatism.py`, `olb/turbulence/profiles.py`.
  ao.apply_compensation: tilt-convention (Noll — explicitly the OTHER
  convention from structure.py; record both sides of C-04); the
  extended-Marechal validity as a constraint. anisoplanatism: isoplanatism
  kind, Stone finite-aperture assumptions, the classic form records its
  up-to-10x overestimate. profiles: mostly NA data plumbing — annotate only
  real assumptions (HV57/Bufton model citations), no checks.

### Prompts WP3a-WP3d (parallel after ALL of WP1+WP2)

Common body:

```
Wire the assumption trace into the Term factories in FILES, spec sections
D5-D7 of docs/assumptions-refactor-plan.md. The worked factory example is
links/downlink.py (WP3a lands it; if merged, imitate it).

Pattern: open `with trace_assumptions() as trace:` around the PHYSICS CALLS
only (not Term assembly); build the record with trace.merge(beam_type=...,
turbulence_regime=..., spectrum=..., validity=<the existing prose>); keep
scenario-level facts the physics never sees as
assumptions.flag(..., source="factory:<module>").

Migration rules: DELETE a factory flag/violation computation ONLY where a
traced check now produces an equivalent violation, and lock the replacement
in the self-check (the same trigger case still yields `not
term.assumptions.ok`). KEEP every existing warnings.warn verbatim. Keep all
Term.meta keys unchanged. Newly enforced constraints (zenith, Gaussian
second condition) may ADD violations to existing cases: assert MEMBERSHIP of
expected violations, never exact counts.

Extend each module self-check with, per wired factory: (1)
term.assumptions.provenance is non-empty and names the expected physics
sources; (2) every pre-existing violation trigger still yields not-ok;
(3) Budget.check() on the module demo budget reports no untraced-guard entry.
```

Package-specific FILES lines and notes:

- **WP3a** — `olb/links/downlink.py`, `olb/models/fade.py`, `olb/models/waveoptics.py`.
  Lands the worked factory example in _lognormal_term/_gamma_gamma_term/
  _auto_select. fade.irradiance_fade_term: UNCHANGED except self-check —
  factories pass the merged record explicitly (auto-merge inside the adapter
  is speculative; do not build it). waveoptics.py: no trace (physics is the
  sim); its Terms self-declare provenance=["untraced: wave-optics
  simulation"] so the untraced guard stays quiet; the constant import is
  already an alias after WP0. Keep the t20-not-ok / t60-ok / obscuration
  self-check assertions.
- **WP3b** — `olb/links/uplink.py`, `olb/links/retro_space.py`.
  Uplink: the Dios sigma_x^2 hard gate moves to the traced uplink_flux
  check; Marechal and NO-SCINTILLATION/NO-FADE stay factory flags with
  source=. retro_space recomposes finished Terms: use
  merge_assumptions(up.assumptions, down.assumptions), no trace of its own.
- **WP3c** — `olb/links/terrestrial.py`, `olb/models/coupling/terrestrial.py`, `olb/models/coupling/_common.py`.
  The three-tier soft/hard warning ladder STAYS factory logic (soft is not a
  violation); only the hard-tier flag and the pdf-shape flag migrate to
  traced checks. THIS PACKAGE CLOSES THE KNOWN GAP: the SMF walk-off Term
  (models/coupling/terrestrial.py, ~line 424) declares REGIME_WEAK and never
  flags — its tilt-variance physics calls (angle_of_arrival,
  gaussian_fried) are now traced, so the weak-limit violation arrives
  automatically; assert that in the self-check with a strong-turbulence
  case. _common.py static coupling stays NA/untraced-exempt (deterministic,
  not a turbulence Term? it IS category "coupling" — give it
  provenance=["untraced: static optics"] to satisfy the guard).
- **WP3d** — `olb/models/fast.py`, `olb/models/coupling/downlink.py`.
  FAST Terms: trace the olb-side physics (plane_wave_scintillation_index,
  scintillation_index, ao chain) AND append "untraced: fast-aosim" to
  provenance for the external Monte Carlo part. coupling/downlink: trace the
  ao chain; the mean-only tier keeps mean_only semantics untouched.

### Prompt WP4 (last, alone)

```
Final sweep for the assumptions refactor, docs/assumptions-refactor-plan.md.

1. Add a coverage audit to the olb/assumptions.py self-check: walk the
   in-scope modules (the 16 physics modules of olb/turbulence/**), count
   functions bearing __assumptions__, assert the count against a stated
   floor (~85; delegators and pure-algebra helpers are exempt) so a dropped
   decorator is caught mechanically.
2. Run the repo /update skill: reconcile docs/architecture.md,
   docs/physics.md, docs/api-budget.md, docs/andrews-crosscheck.md (link the
   conflict-tag constraints C-01..C-07 to their kind="conflict" records),
   README.md, and the CLAUDE.md architecture/current-state notes with the
   new mechanism.
3. Run every in-scope module self-check (loop python -m over the touched
   modules) and the examples/ scripts that touch links; confirm
   Budget.check() output on the example budgets is warning-compatible with
   before (same warnings, plus provenance-prefixed violation strings).
```

## Worked example (annotation pattern, lands in WP1b)

```python
# top of olb/turbulence/andrews/scintillation.py
from ...assumptions import Constraint, module_assumptions, SPECTRUM_KOLMOGOROV

assumes = module_assumptions(
    spectrum=SPECTRUM_KOLMOGOROV,
    constraints=(
        Constraint("path-homogeneity",
                   "One path length L and one scalar Cn2. No profile integral.",
                   "10.1117/3.626196", "Ch. 8, Eq. (4), printed p. 261"),
    ),
)

def _weak_regime_check(args, result):
    '''Return a reason when the Rytov gate reads "hard". No warning here.'''
    label = rytov_weak(float(np.max(result)))
    if label == 'hard':
        return (f"sigma_R^2 = {float(np.max(result)):.3f} >= "
                f"{WEAK_REGIME_LIMIT}; the Rytov weak theory does not hold.")
    return None

WEAK_REGIME_CONSTRAINT = Constraint(
    "regime", "Weak fluctuation: sigma_R^2 < 1.",
    "10.1117/3.626196", "Ch. 8, text below Eq. (23), printed pp. 264-265",
    check=_weak_regime_check)

@assumes(WEAK_REGIME_CONSTRAINT, turbulence_regime=REGIME_WEAK)
def rytov_variance(wavelength, z, cn2, *, wave='plane'):
    ...unchanged body...
```

And the factory pattern (lands in WP3a, `links/downlink._lognormal_term`):

```python
with trace_assumptions() as trace:
    sigma2_I = plane_wave_scintillation_index(elev, wavelength, hs, cn2_profile)
    if aperture_average:
        sigma2_P = aperture_averaged_scintillation_index(D, elev, wavelength,
                                                         hs, cn2_profile)
assumptions = trace.merge(
    beam_type=BEAM_PLANE_WAVE, turbulence_regime=REGIME_WEAK,
    spectrum=SPECTRUM_KOLMOGOROV, validity="...unchanged prose...")
if aperture_average and obscuration > 0.0:
    assumptions.flag("The receive aperture has a central obscuration ...",
                     source="factory:links.downlink")
```

## Verification

- WP0: `python -m olb.assumptions`, `python -m olb.results`, plus the two
  alias modules; grep shows one 0.25 definition.
- Every WP1/WP2 module: `python -m olb.turbulence...<module>` with the three
  new self-check blocks (value parity / registration / violation with no
  warning).
- Every WP3 module: `python -m olb.links...` with provenance + regression +
  guard-quiet assertions.
- WP4: the coverage audit, the full self-check loop, examples/ link scripts,
  and warning-compatibility of `Budget.check()` on the demo budgets.

## Risks and owner decisions

1. **Violation strings become load-bearing.** The `[source]` prefix changes
   the text of migrated violations; the known consumers (results.py joins,
   self-check asserts) are updated in WP3. Grep for other string-matching
   before merge.
2. **Newly enforced constraints flip `ok` in existing cases.** The zenith
   check and the Gaussian second condition will newly flag low-elevation /
   focused-beam cases the moment WP3 lands. This is the POINT of the
   refactor (catching missed flags), and WP3 asserts membership not counts —
   but the owner should expect, for example, a 25-30 deg downlink to read
   not-ok where it read ok before. Sign-off happens by reviewing the WP3
   self-check output.
3. **Threading.** Workers do not inherit contextvars; the untraced guard is
   the net. Documented rule: open the context in the thread that makes the
   physics calls, or pass `contextvars.copy_context()`. No code now
   (waveoptics out of scope).
4. **Delegator noise.** "Pure delegator gets no decorator" rule; the WP4
   coverage audit exempts them.
5. **Check closures.** A Constraint with a check lives in the module whose
   function it checks; only check-less constraints may be shared across
   modules (avoids import-order coupling).
6. Deferred (YAGNI): soft-tier labels in the machine record; structured
   violation objects; auto-merge inside fade.irradiance_fade_term; the
   waveoptics numerical-sampling assumption family.

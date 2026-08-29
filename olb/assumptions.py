'''
Model assumptions for a Term, and the check that flags a scenario that breaks them.

Each analytic or Monte Carlo model is valid only in a regime. This module gives
a small record that states the regime. A model attaches an Assumptions record to
its Term. The model adds a reason to `violations` when the scenario breaks an
assumption. A Budget then reports which terms are misrepresentative.

Three constraints matter for the optical propagation models:
- beam_type: the wavefront the model assumes (plane wave, spherical wave, or
  Gaussian beam).
- turbulence_regime: the fluctuation strength the model assumes (weak, moderate,
  or strong). The regime is tied to a numeric bound on the scintillation index.
- spectrum: the turbulence spectrum the model assumes (for example Kolmogorov
  with no inner scale and no outer scale).

Use the string constants below so every term uses the same words.

Function-owned assumptions
--------------------------
A physics function states its own validity through the `@assumes` decorator. The
decorator attaches a machine-readable `FuncAssumptions` record to the function
and, INSIDE a collection context only, registers the record and runs each
`Constraint` check. A Term factory opens the context with `trace_assumptions()`
around its physics calls; every decorated function that runs registers itself, so
the Term inherits the union automatically and a forgotten dependency is
impossible. `Trace.merge(...)` turns the trace into one `Assumptions` record.

The context is a `contextvars.ContextVar`. Outside a context the decorator does
ONE `get` and calls the function, so the numeric output is byte-identical. A
`ContextVar` is thread-local: a `ThreadPoolExecutor` or `threading.Thread` worker
does NOT inherit the caller's context, so open the context in the thread that
makes the physics calls (the untraced guard in `Budget.check()` is the net).
'''

import functools
import inspect
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
from typing import Callable, Optional

# Beam type (wavefront the model assumes).
BEAM_PLANE_WAVE = "plane wave"
BEAM_SPHERICAL_WAVE = "spherical wave"
BEAM_GAUSSIAN = "Gaussian beam"
BEAM_NA = "not applicable"

# Turbulence regime (fluctuation strength the model assumes).
REGIME_WEAK = "weak"
REGIME_MODERATE = "moderate"
REGIME_STRONG = "strong"
REGIME_NA = "not applicable"

# Turbulence spectrum. The five named models are the ones that
# olb.turbulence.andrews.spectra builds. Source: Andrews and Phillips, Laser
# Beam Propagation through Random Media, 2nd ed. (2005), DOI 10.1117/3.626196,
# Ch. 3, Eqs. (18) to (23), printed pp. 67 to 69.
SPECTRUM_KOLMOGOROV = "Kolmogorov, no inner or outer scale"      # Eq. (18)
SPECTRUM_TATARSKII = "Tatarskii, finite inner scale"             # Eq. (19)
SPECTRUM_VON_KARMAN = "von Karman, finite inner or outer scale"  # Eq. (20)
SPECTRUM_EXPONENTIAL = "exponential outer-scale cut-off"         # Eq. (21)
SPECTRUM_MODIFIED = "modified atmospheric, inner-scale bump"     # Eqs. (22), (23)
SPECTRUM_NA = "not applicable"

# The assumption KINDS. Each Constraint carries one slug from this set, so a
# reader can group the validity limits by axis. The catalogue (2026-08-29) found
# these ~21 recurring kinds across olb/turbulence/**. No enum, no subclass: a
# slug is a plain string, and __post_init__ rejects an unknown one.
KINDS = frozenset({
    "beam-type",            # the wavefront the model assumes
    "regime",               # the fluctuation strength (Rytov axis)
    "spectrum",             # the turbulence spectrum
    "pdf-shape",            # the fade PDF the model trusts (lognormal axis)
    "launch-curvature",     # the launch phase-front curvature branch
    "receiver",             # the receive-aperture shape the filter assumes
    "obscuration",          # a central obscuration the model does or does not carry
    "aperture-order",       # a size ordering of aperture versus beam
    "on-axis",              # an on-axis (no wander) evaluation
    "tracking",             # a tilt-removal (tracked) assumption
    "path-homogeneity",     # one L and one scalar Cn2 versus a profile integral
    "path-weight",          # the reference end of a path-weighted integral
    "geometry",             # plane-parallel path, no Earth curvature
    "zenith",               # a zenith-angle (elevation) bound
    "field-region",         # a near-field or far-field restriction
    "variance-convention",  # radial versus one-axis variance
    "tilt-convention",      # G-tilt versus Noll (Zernike) tilt
    "approximation",        # a stated mathematical approximation
    "frozen-flow",          # the Taylor frozen-flow hypothesis
    "isoplanatism",         # an angular-decorrelation (anisoplanatic) assumption
    "not-built",            # a branch that raises instead of returning a value
    "conflict",             # a cross-source disagreement noted in the code
})


@dataclass(frozen=True)
class Constraint:
    '''One validity limit of a physics function.

    A Constraint is a machine-readable statement of one assumption. It carries
    the KIND slug, one ASD-STE100 statement, the source DOI, and the printed
    location. An optional `check(args, result) -> Optional[str]` tests the
    scenario at run time: it returns None when the scenario is inside the limit,
    or one ASD-STE100 reason string when the scenario breaks it. A check NEVER
    warns and NEVER raises (warnings stay factory-level).
    '''
    kind: str                 # one slug from KINDS
    statement: str            # ASD-STE100 prose
    doi: str                  # source DOI, e.g. "10.1117/3.626196"
    where: str = ""           # e.g. "Ch. 5, Eq. (15), printed p. 140"
    check: Optional[Callable] = field(default=None, compare=False)

    def __post_init__(self):
        if self.kind not in KINDS:
            raise ValueError(
                f"unknown Constraint kind {self.kind!r}; use one of "
                f"{sorted(KINDS)}."
            )


@dataclass(frozen=True)
class FuncAssumptions:
    '''The assumption record the `@assumes` decorator stores on a function.

    `source` is "<module>.<qualname>", so coverage is machine-auditable with no
    context. The headline fields are None when the function does not constrain
    that axis. `constraints` is the union of the module defaults and the
    per-function constraints.
    '''
    source: str
    beam_type: Optional[str] = None
    turbulence_regime: Optional[str] = None
    spectrum: Optional[str] = None
    constraints: tuple = ()


@dataclass
class Assumptions:
    '''The regime one model is valid in, and the reasons a scenario breaks it.'''
    beam_type: str
    turbulence_regime: str
    spectrum: str
    validity: str = ""                       # the numeric limit, in words
    violations: list = field(default_factory=list)   # reasons the scenario breaks the model
    constraints: list = field(default_factory=list)  # (source, Constraint) pairs, from the trace
    provenance: list = field(default_factory=list)   # the traced physics source names

    @property
    def ok(self) -> bool:
        '''Return True when the scenario breaks no assumption.'''
        return len(self.violations) == 0

    def flag(self, reason, source=None):
        '''Add one reason that the scenario breaks an assumption. Return self.

        A scenario-level fact that the physics never sees (a central obscuration,
        the extended-Marechal limit, NO SCINTILLATION) passes `source` so the
        violation carries the same "[source] reason" prefix as a traced check.
        '''
        entry = reason if source is None else f"[{source}] {reason}"
        self.violations.append(entry)
        return self


# ----------------------------------------------------------------------------
# The collection context.
# ----------------------------------------------------------------------------

# The active trace. None outside a context, so the decorator does ONE get and
# returns. A ContextVar is thread-local (a worker thread does not inherit it).
_ACTIVE: ContextVar = ContextVar("olb_assumptions", default=None)


class Trace:
    '''The record a `trace_assumptions()` context collects.

    `records` maps each physics source name to its FuncAssumptions (dedup by
    source: the same function that runs many times registers once). `violations`
    is a list of "[source] reason" strings, deduped exactly.
    '''

    def __init__(self):
        self.records = {}       # source -> FuncAssumptions
        self.violations = []    # "[source] reason", exact dedup

    def record(self, rec: FuncAssumptions):
        '''Register one function record (dedup by source).'''
        self.records[rec.source] = rec

    def violate(self, source, reason):
        '''Append one "[source] reason" violation (exact dedup).'''
        entry = f"[{source}] {reason}"
        if entry not in self.violations:
            self.violations.append(entry)

    def merge(self, *, beam_type=None, turbulence_regime=None, spectrum=None,
              validity="") -> "Assumptions":
        '''Fold the trace into one Assumptions record.

        A headline field takes the explicit factory kwarg when given; otherwise
        the unanimous non-NA traced value; a real disagreement joins the values
        with "/" AND appends a "[merge]" violation that tells the factory to
        state the field. The constraints are the union of the traced
        (source, Constraint) pairs, deduped by (source, kind, statement). The
        violations are the traced ones plus any "[merge]" note. `validity` is the
        factory prose, unchanged in role.
        '''
        merge_violations = []
        records = list(self.records.values())
        beam = _resolve_headline(
            beam_type, [r.beam_type for r in records], BEAM_NA,
            "beam_type", merge_violations)
        regime = _resolve_headline(
            turbulence_regime, [r.turbulence_regime for r in records], REGIME_NA,
            "turbulence_regime", merge_violations)
        spec = _resolve_headline(
            spectrum, [r.spectrum for r in records], SPECTRUM_NA,
            "spectrum", merge_violations)

        pairs = _union_constraints(records)
        return Assumptions(
            beam_type=beam,
            turbulence_regime=regime,
            spectrum=spec,
            validity=validity,
            violations=list(self.violations) + merge_violations,
            constraints=pairs,
            provenance=list(self.records.keys()),
        )


@contextmanager
def trace_assumptions():
    '''Open a collection context and yield the Trace.

    Every decorated physics function that runs inside registers itself and its
    violations. A nested context shadows the outer one (natural ContextVar
    behaviour): a function that runs in the inner context registers to the inner
    Trace only.
    '''
    trace = Trace()
    token = _ACTIVE.set(trace)
    try:
        yield trace
    finally:
        _ACTIVE.reset(token)


def module_assumptions(*, beam_type=None, turbulence_regime=None, spectrum=None,
                       constraints=()):
    '''Return an `@assumes` decorator factory that carries these module defaults.

    The returned `assumes(*constraints, beam_type=..., turbulence_regime=...,
    spectrum=...)` decorates one function. The function record is the UNION of
    the module defaults and the per-function constraints; a per-function headline
    keyword overrides the module default for that axis. Use a module default ONLY
    for a statement true of EVERY public function in the file; one exception
    forbids the module form for that axis (the union cannot REMOVE a default).
    '''
    module_beam = beam_type
    module_regime = turbulence_regime
    module_spectrum = spectrum
    module_constraints = tuple(constraints)

    def assumes(*func_constraints, beam_type=None, turbulence_regime=None,
                spectrum=None):
        headline_beam = beam_type if beam_type is not None else module_beam
        headline_regime = (turbulence_regime if turbulence_regime is not None
                           else module_regime)
        headline_spectrum = spectrum if spectrum is not None else module_spectrum
        merged = _dedup_constraints(module_constraints + tuple(func_constraints))
        checks = [c for c in merged if c.check is not None]

        def decorate(func):
            source = f"{func.__module__}.{func.__qualname__}"
            record = FuncAssumptions(
                source=source,
                beam_type=headline_beam,
                turbulence_regime=headline_regime,
                spectrum=headline_spectrum,
                constraints=merged,
            )
            sig = inspect.signature(func) if checks else None

            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                trace = _ACTIVE.get()
                if trace is None:                       # zero-overhead bypass
                    return func(*args, **kwargs)
                result = func(*args, **kwargs)
                trace.record(record)
                if checks:
                    bound = sig.bind(*args, **kwargs)
                    bound.apply_defaults()
                    bound_args = bound.arguments
                    for c in checks:
                        reason = c.check(bound_args, result)
                        if reason:
                            trace.violate(source, reason)
                return result

            wrapper.__assumptions__ = record
            return wrapper

        return decorate

    return assumes


# The no-default form, exported for a per-function decoration.
assumes = module_assumptions()


def merge_assumptions(*records, validity="") -> "Assumptions":
    '''Recompose finished Assumptions records into one (for links/retro_space.py).

    A Term that is built from finished Terms (the retro link folds the uplink and
    the downlink) recomposes their records without a trace of its own. The rule
    matches Trace.merge: a headline field is the unanimous non-NA value or the
    "/"-joined disagreement with a "[merge]" note; the constraints, the
    provenance, and the violations are the unions.
    '''
    records = [r for r in records if r is not None]
    merge_violations = []
    beam = _resolve_headline(
        None, [r.beam_type for r in records], BEAM_NA, "beam_type",
        merge_violations)
    regime = _resolve_headline(
        None, [r.turbulence_regime for r in records], REGIME_NA,
        "turbulence_regime", merge_violations)
    spec = _resolve_headline(
        None, [r.spectrum for r in records], SPECTRUM_NA, "spectrum",
        merge_violations)

    constraints = []
    seen = set()
    for r in records:
        for source, c in getattr(r, "constraints", ()):
            key = (source, c.kind, c.statement)
            if key not in seen:
                seen.add(key)
                constraints.append((source, c))

    provenance = []
    for r in records:
        for source in getattr(r, "provenance", ()):
            if source not in provenance:
                provenance.append(source)

    violations = []
    for r in records:
        for v in r.violations:
            if v not in violations:
                violations.append(v)
    violations += merge_violations

    joined_validity = validity or " | ".join(
        v for v in dict.fromkeys(r.validity for r in records if r.validity))
    return Assumptions(
        beam_type=beam,
        turbulence_regime=regime,
        spectrum=spec,
        validity=joined_validity,
        violations=violations,
        constraints=constraints,
        provenance=provenance,
    )


# ----------------------------------------------------------------------------
# Helpers.
# ----------------------------------------------------------------------------

def _dedup_constraints(constraints):
    '''Return the constraints deduped by (kind, statement), order preserved.'''
    out = []
    seen = set()
    for c in constraints:
        key = (c.kind, c.statement)
        if key not in seen:
            seen.add(key)
            out.append(c)
    return tuple(out)


def _resolve_headline(explicit, values, na_value, field_name, violations):
    '''Resolve one headline axis (see Trace.merge for the rule).'''
    if explicit is not None:
        return explicit
    seen = []
    for v in values:
        if v is not None and v != na_value and v not in seen:
            seen.append(v)
    if not seen:
        return na_value
    if len(seen) == 1:
        return seen[0]
    joined = "/".join(sorted(seen))
    violations.append(
        f"[merge] the traced {field_name} values disagree ({joined}); the "
        f"factory must state {field_name} explicitly."
    )
    return joined


def _union_constraints(records):
    '''Return the (source, Constraint) union, deduped by (source, kind, statement).'''
    pairs = []
    seen = set()
    for rec in records:
        for c in rec.constraints:
            key = (rec.source, c.kind, c.statement)
            if key not in seen:
                seen.add(key)
                pairs.append((rec.source, c))
    return pairs


if __name__ == '__main__':
    import threading
    import warnings

    # A throwaway constraint whose check fires when the result exceeds a bound.
    def _over_one(args, result):
        if float(result) > 1.0:
            return f"the value {float(result):.3f} exceeds the limit 1.0."
        return None

    _LIMIT = Constraint("regime", "The value obeys value <= 1.",
                        "10.0000/test", "test suite", check=_over_one)

    demo = module_assumptions(spectrum=SPECTRUM_KOLMOGOROV)

    @demo(_LIMIT, beam_type=BEAM_PLANE_WAVE, turbulence_regime=REGIME_WEAK)
    def scaled(x):
        '''Return x unchanged (a stand-in physics function).'''
        return x

    @demo(beam_type=BEAM_GAUSSIAN, turbulence_regime=REGIME_WEAK)
    def other(x):
        '''A second function with a different beam type, for the merge conflict.'''
        return x

    # (a) Identical return value with and without a context.
    outside = scaled(0.5)
    with trace_assumptions():
        inside = scaled(0.5)
    assert outside == inside == 0.5, "the return value must not change in a context"

    # (b) In-context registration; a failing check appends a source-prefixed
    #     violation; the physics layer emits NO warning.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with trace_assumptions() as trace:
            scaled(2.0)                 # 2.0 > 1.0 -> the check fires
    assert f"{__name__}.scaled" in trace.records, "the source must register"
    assert any(v.startswith(f"[{__name__}.scaled]") for v in trace.violations), \
        "the violation must carry the source prefix"
    assert len(caught) == 0, "a check must not warn"

    # (c) Dedup on repeated calls: same source once, same violation once.
    with trace_assumptions() as trace:
        scaled(2.0)
        scaled(2.0)
        scaled(3.0)                     # a different value, but the same reason text differs
    assert list(trace.records) == [f"{__name__}.scaled"], "the source dedups"
    # Two distinct reason strings (2.000 and 3.000) but each appears once.
    assert len(trace.violations) == len(set(trace.violations)) == 2, \
        "violations dedup exactly"

    # (d) A headline conflict on merge produces the [merge] flag.
    with trace_assumptions() as trace:
        scaled(0.5)                     # plane wave
        other(0.5)                      # Gaussian beam
    merged = trace.merge(validity="demo")
    assert "/" in merged.beam_type, "a disagreement joins with /"
    assert any(v.startswith("[merge]") and "beam_type" in v
               for v in merged.violations), "the merge appends a [merge] flag"
    assert merged.spectrum == SPECTRUM_KOLMOGOROV, "the unanimous module default holds"
    assert not merged.ok

    # A clean merge: unanimous fields, no [merge] flag.
    with trace_assumptions() as trace:
        scaled(0.5)
    clean = trace.merge(turbulence_regime=REGIME_WEAK, validity="ok demo")
    assert clean.beam_type == BEAM_PLANE_WAVE and clean.turbulence_regime == REGIME_WEAK
    assert clean.spectrum == SPECTRUM_KOLMOGOROV
    assert clean.provenance == [f"{__name__}.scaled"]
    assert len(clean.constraints) == 1 and clean.constraints[0][1].kind == "regime"
    assert clean.ok

    # (e) An unknown kind raises.
    try:
        Constraint("bogus-kind", "nope", "10.0/x")
    except ValueError as e:
        assert "unknown Constraint kind" in str(e)
    else:
        raise AssertionError("an unknown kind must raise")

    # (f) Nested contexts shadow: a call in the inner context registers to the
    #     inner Trace only.
    with trace_assumptions() as outer:
        with trace_assumptions() as inner:
            scaled(0.5)
        assert f"{__name__}.scaled" in inner.records
        assert f"{__name__}.scaled" not in outer.records, "the inner context shadows"
        scaled(0.5)
        assert f"{__name__}.scaled" in outer.records, "the outer context resumes"

    # (g) A worker thread does not see the caller's trace (ContextVar is
    #     thread-local), so the physics runs untraced there.
    seen = {}

    def _worker():
        seen["active"] = _ACTIVE.get()
        scaled(2.0)                     # runs, but registers nowhere

    with trace_assumptions() as trace:
        t = threading.Thread(target=_worker)
        t.start()
        t.join()
    assert seen["active"] is None, "a worker thread must not inherit the context"
    assert not trace.violations, "the worker's call must not reach the caller trace"

    # merge_assumptions recomposes finished records (the retro-link path). The
    # uplink record (plane wave, one traced source) folds with a downlink record
    # (Gaussian), so the beam types disagree and a [merge] flag appears.
    with trace_assumptions() as ta:
        scaled(0.5)
    up = ta.merge(validity="up")
    down = Assumptions(BEAM_GAUSSIAN, REGIME_WEAK, SPECTRUM_KOLMOGOROV,
                       validity="down")
    recomposed = merge_assumptions(up, down)
    assert "/" in recomposed.beam_type
    assert any(v.startswith("[merge]") for v in recomposed.violations)
    assert recomposed.provenance == [f"{__name__}.scaled"]
    assert recomposed.validity == "up | down"

    # Coverage audit: walk olb.turbulence and count the DISTINCT functions that
    # (a) carry __assumptions__ and (b) are defined in the module the walk visits
    # (f.__module__ == that module, to skip a re-export). The measured count is
    # 88 across 18 turbulence modules (2026-08-29). Assert a floor of 85, because
    # a pure delegator and a pure-algebra helper are deliberately exempt (they get
    # no decorator, per the refactor plan). A dropped decorator drops the count
    # below the floor and this assertion catches it with no context.
    import importlib
    import pkgutil

    import olb.turbulence as _turb_pkg

    _covered = set()
    for _mi in pkgutil.walk_packages(_turb_pkg.__path__, _turb_pkg.__name__ + "."):
        _mod = importlib.import_module(_mi.name)
        for _obj in vars(_mod).values():
            _rec = getattr(_obj, "__assumptions__", None)
            if _rec is None:
                continue
            if getattr(_obj, "__module__", None) != _mi.name:
                continue                       # a re-export, counted at its home
            _covered.add(_rec.source)
    _n_covered = len(_covered)
    assert _n_covered >= 85, (
        f"assumptions coverage dropped to {_n_covered} decorated turbulence "
        f"functions (floor 85); a physics function lost its @assumes decorator."
    )
    print(f"assumptions coverage: {_n_covered} decorated turbulence functions.")

    print("assumptions.py self-check passed.")

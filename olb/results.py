'''
The result objects that every model makes and the budget uses.

A Term is one line of the link budget. It gives three views of the same
contribution. The choice of analytic or Monte Carlo is the choice of view:

    term.mean_db            deterministic value / expected loss  (+ = loss)
    term.quantile_db(p)     analytic loss at availability p, if a closed form
                            exists; None for terms that only have samples
    term.sample_db(n, rng)  n Monte Carlo draws of the contribution

A deterministic term, such as geometric loss, sets only mean_db. olb broadcasts
it into samples and uses a constant quantile. A statistical term with a closed
form, such as log-normal scintillation, gives quantile_db and a sampler. A term
with only a Monte Carlo model, such as the coupled-flux beam wander and
scintillation, gives only a sampler, and quantile_db returns None. This value
tells the budget to use Monte Carlo, not the analytic sum.

A Budget is a list of Terms with the optional top-line values (tx power, rx
sensitivity). It reports the deterministic total, an itemised table, an analytic
fade margin, or a full Monte Carlo of the joint distribution.
'''

import warnings
from dataclasses import dataclass, field
from typing import Callable, Optional, Union

import numpy as np
import pandas as pd

from .assumptions import Assumptions, SPECTRUM_NA

# A loss or gain value in dB: a scalar, or a per-geometry array (for example one
# value per elevation angle). Loss is positive, gain negative.
NumDB = Union[float, np.ndarray]


@dataclass
class Term:
    '''One contribution to the link budget (loss positive, gain negative, dB).'''
    name: str
    category: str                       # geometric | atmospheric | turbulence | pointing | system
    mean_db: NumDB = 0.0                 # deterministic value / E[loss], scalar or per-geometry array
    sampler: Optional[Callable[[int, np.random.Generator], np.ndarray]] = None
    quantile: Optional[Callable[[float], NumDB]] = None
    note: str = ""
    meta: dict = field(default_factory=dict)
    assumptions: Optional[Assumptions] = None   # the regime the model is valid in
    mean_only: bool = False   # fidelity-0: models only E[loss] of a fluctuating quantity, no fade

    @property
    def stochastic(self) -> bool:
        return self.sampler is not None

    def sample_db(self, n, rng):
        '''n Monte Carlo draws, shape (n, *base_shape). Deterministic terms broadcast.'''
        if self.sampler is not None:
            return np.asarray(self.sampler(n, rng))
        base = np.shape(self.mean_db)
        return np.broadcast_to(np.asarray(self.mean_db), (n, *base)).copy()

    def quantile_db(self, p) -> Optional[NumDB]:
        '''
        Loss at availability p (the p-quantile of the loss of this term).
        Return None for Monte-Carlo-only terms that have no closed form. This
        value tells the budget that the analytic sum cannot include this term.
        '''
        if self.quantile is not None:
            return self.quantile(p)
        if self.sampler is None:
            return self.mean_db          # deterministic: constant across availability
        return None


class EmpiricalSampler:
    '''The two stochastic faces of a Term, built from a finite set of loss samples.

    A Monte-Carlo Term (FAST, wave optics) holds a set of per-realisation loss
    samples (dB, positive = loss) and needs the same two faces from them: a
    sampler that draws any number of values, and an empirical quantile. This
    class gives both, and it OWNS the tail adequacy, so no caller re-derives it:

        sampler = EmpiricalSampler(loss_db)
        Term(..., mean_db=sampler.mean_db, sampler=sampler, quantile=sampler.quantile)

    - the object is callable, so it is the Term sampler: it bootstrap-resamples
      the loss samples (draw with replacement);
    - `quantile(p)` is the empirical loss at availability p;
    - a p-availability quantile reads only n*(1-p) samples in the tail. Below
      `TAIL_MIN` samples there the quantile is not a design number, so
      `undersampled(p)` reports it and `quantile(p)` warns. This is a PROPERTY of
      the sampler, not of any caller: `term.sampler.undersampled(p)` answers it.
    '''

    # The fewest samples that must sit deeper than the availability for the
    # quantile to be a trustworthy design number. An olb convention, not a book
    # value: below about ten tail samples the empirical quantile is noise.
    TAIL_MIN = 10.0

    def __init__(self, loss_db):
        self.loss_db = np.asarray(loss_db, dtype=float).ravel()
        if self.loss_db.size == 0:
            raise ValueError("EmpiricalSampler needs at least one loss sample.")

    @property
    def n(self) -> int:
        '''The number of loss samples.'''
        return int(self.loss_db.size)

    @property
    def mean_db(self) -> float:
        '''The empirical mean loss (the Term mean face).'''
        return float(self.loss_db.mean())

    def tail_count(self, p) -> float:
        '''The number of samples deeper than availability p, n*(1-p).'''
        return self.n * (1.0 - float(p))

    def undersampled(self, p) -> bool:
        '''True when too few samples sit in the p tail to trust the quantile.'''
        return self.tail_count(p) < self.TAIL_MIN

    def trustworthy_n(self, p) -> int:
        '''The sample count that puts about TAIL_MIN samples in the p tail.'''
        return int(round(self.TAIL_MIN / (1.0 - float(p))))

    def quantile(self, p) -> float:
        '''The empirical loss at availability p (the Term quantile face).'''
        if self.undersampled(p):
            warnings.warn(
                f"the {float(p):g} quantile reads only about "
                f"{self.tail_count(p):.1f} of {self.n} samples in the tail, so it "
                f"is UNDER-SAMPLED and not a design number. Raise the sample count "
                f"(about {self.trustworthy_n(p)} for a trustworthy quantile)."
            )
        return float(np.quantile(self.loss_db, p))

    def __call__(self, n, rng) -> np.ndarray:
        '''Bootstrap resample: n draws from the loss samples, with replacement.'''
        return rng.choice(self.loss_db, size=n, replace=True)


class Budget:
    '''A collection of Terms that form one link budget.'''

    def __init__(self, terms=None, tx_power_dbm=None, rx_sensitivity_dbm=None,
                 scenario=None):
        self.terms = list(terms) if terms else []
        self.scenario = scenario
        # Use the scenario terminals' top-line values when the caller gives none.
        # The launch power lives on the transmit terminal's Transmitter; the
        # sensitivity lives on the receive terminal's Detector. Both are optional.
        tx = getattr(scenario, "tx_terminal", None)
        rx = getattr(scenario, "rx_terminal", None)
        transmitter = getattr(tx, "transmitter", None)
        detector = getattr(rx, "detector", None)
        self.tx_power_dbm = tx_power_dbm if tx_power_dbm is not None \
            else getattr(transmitter, "power_dbm", None)
        self.rx_sensitivity_dbm = rx_sensitivity_dbm if rx_sensitivity_dbm is not None \
            else getattr(detector, "sensitivity_dbm", None)

    def add(self, term: Term):
        '''Append a Term and return self (chainable).'''
        self.terms.append(term)
        return self

    # --- fidelity -----------------------------------------------------------

    def mean_only_terms(self):
        '''
        Return the Terms that model only a mean (fidelity-0, no fade).

        A mean-only Term gives the expected loss of a quantity that really
        fluctuates (for example a fibre-coupling loss computed from the mean
        residual wavefront). It is NOT the same as a deterministic Term such as
        geometric spreading, whose mean IS its whole distribution. A mean-only
        Term therefore has no trustworthy fade, and it locks the whole budget to
        fidelity 0.
        '''
        return [t for t in self.terms if t.mean_only]

    @property
    def provides_fade(self) -> bool:
        '''
        False when any Term is mean-only, so a budget fade would mislead.

        A fade margin (or a Monte Carlo fade) adds the fade of every Term. When
        one Term is a mean-only fidelity-0 approximation, its fade is missing, so
        the total tail is understated. In that case the budget reports the mean
        only and refuses the fade.
        '''
        return not any(t.mean_only for t in self.terms)

    # --- deterministic view -------------------------------------------------

    def total_loss_db(self):
        '''Sum of the deterministic (mean) contributions [dB].'''
        return sum(np.asarray(t.mean_db) for t in self.terms)

    def to_frame(self):
        '''Itemised budget as a DataFrame (one row per Term).'''
        rows = [{
            "name": t.name,
            "category": t.category,
            "mean_db": t.mean_db,
            "stochastic": t.stochastic,
            "note": t.note,
        } for t in self.terms]
        return pd.DataFrame(rows, columns=["name", "category", "mean_db",
                                           "stochastic", "note"])

    # --- model assumptions --------------------------------------------------

    def assumptions_frame(self):
        '''
        Table of the model assumptions, one row per Term.

        Each row states the beam type, the turbulence regime, the spectrum, and
        the validity limit. The `ok` column is False when the scenario breaks an
        assumption. The `violations` column gives the reasons.

        Returns:
            pandas.DataFrame
        '''
        rows = []
        for t in self.terms:
            a = t.assumptions
            if a is None:
                rows.append({"name": t.name, "beam_type": "", "regime": "",
                             "spectrum": "", "validity": "", "ok": True,
                             "violations": ""})
            else:
                rows.append({
                    "name": t.name,
                    "beam_type": a.beam_type,
                    "regime": a.turbulence_regime,
                    "spectrum": a.spectrum,
                    "validity": a.validity,
                    "ok": a.ok,
                    "violations": "; ".join(a.violations),
                })
        return pd.DataFrame(rows, columns=["name", "beam_type", "regime",
                                           "spectrum", "validity", "ok",
                                           "violations"])

    def check(self, warn=True):
        '''
        Find the terms whose assumptions the scenario breaks.

        Parameters:
            warn : bool
                Issue a warning for each broken assumption when True.

        Returns:
            list of (str, str)
                One (term name, reason) pair for each violation. A mixed-spectrum
                budget adds one ("budget", reason) pair (see below).
        '''
        found = []
        for t in self.terms:
            a = t.assumptions
            if a is None:
                continue
            for reason in a.violations:
                found.append((t.name, reason))
                if warn:
                    warnings.warn(
                        f"{t.name}: the scenario breaks a model assumption. "
                        f"{reason} The analytic result can be misrepresentative."
                    )

        # Cross-term spectrum consistency. Every turbulence-bearing term models
        # the SAME atmosphere, so they must agree on the spectrum. A budget that
        # mixes, say, a Kolmogorov analytic term with a von Karman FAST term is
        # physically inconsistent. Terms with SPECTRUM_NA (no turbulence spectrum,
        # such as geometric or pointing) do not constrain this.
        spectra = {t.assumptions.spectrum for t in self.terms
                   if t.assumptions is not None
                   and t.assumptions.spectrum not in ("", SPECTRUM_NA)}
        if len(spectra) > 1:
            reason = ("the budget mixes turbulence spectra "
                      f"({', '.join(sorted(spectra))}); the terms model the same "
                      "atmosphere and must assume one spectrum.")
            found.append(("budget", reason))
            if warn:
                warnings.warn(f"budget: {reason}")
        return found

    # --- analytic fade ------------------------------------------------------

    def fade_margin_db(self, availability):
        '''
        Analytic loss at the given availability. This is the sum of the
        p-quantile loss of each term. This bound ignores the independence of the
        terms. It is an upper bound because it adds the worst case of every term
        together. Use monte_carlo() for the joint distribution.

        Raises:
            ValueError
                If any term is mean-only (fidelity-0): the budget then has no
                trustworthy fade, so it refuses rather than return a misleading
                number. Read total_loss_db() for the mean.
            ValueError
                If any term is Monte-Carlo-only (no closed-form quantile).
                Evaluate those terms with monte_carlo().
        '''
        mean_only = self.mean_only_terms()
        if mean_only:
            names = ", ".join(repr(t.name) for t in mean_only)
            raise ValueError(
                f"Budget has mean-only (fidelity-0) term(s) {names}: they model "
                "only the expected loss, not the fade. A fade margin would add the "
                "other terms' fades to a fibre-coupling MEAN and misrepresent the "
                "tail. Read total_loss_db() for the mean loss, or use a statistical "
                "(fidelity-1) coupling model to get the coupling fade."
            )
        total = 0.0
        for t in self.terms:
            q = t.quantile_db(availability)
            if q is None:
                raise ValueError(
                    f"Term {t.name!r} has no closed-form quantile; use "
                    "monte_carlo() to evaluate this budget."
                )
            total = total + np.asarray(q)
        return total

    # --- Monte Carlo --------------------------------------------------------

    def monte_carlo(self, n, rng=None, availabilities=(0.99,)):
        '''
        Draw n joint samples of the total loss and summarise the distribution.

        The code samples every term with the same rng and sums the samples per
        draw. This keeps the correlations inside a term (for example coupled-flux
        wander and scintillation). Independent terms combine correctly by
        construction.

        Parameters:
            n : int
                Number of Monte Carlo draws.
            rng : numpy.random.Generator, optional
                Seeded generator for reproducibility.
            availabilities : iterable of float
                Availabilities at which to report the fade (loss) level.

        Returns:
            dict
                total_loss_samples : ndarray, shape (n, *base_shape)
                mean_loss_db       : ndarray, E[loss]
                fade_db            : dict {availability: loss level}
                received_dbm       : ndarray or None (if tx power set)
                margin_db          : dict {availability: margin} or None
        '''
        rng = np.random.default_rng() if rng is None else rng
        samples = sum(t.sample_db(n, rng) for t in self.terms)
        samples = np.asarray(samples)

        # A mean-only (fidelity-0) term broadcasts a constant instead of a fade,
        # so the joint tail is understated. Report the mean, but refuse the fade
        # and the fade-based margin rather than show a misleading number.
        provides_fade = self.provides_fade
        if provides_fade:
            fade = {a: np.percentile(samples, a * 100.0, axis=0)
                    for a in availabilities}
        else:
            names = ", ".join(repr(t.name) for t in self.mean_only_terms())
            warnings.warn(
                f"budget has mean-only (fidelity-0) term(s) {names}; the fade and "
                "the fade margin are suppressed. Only the mean is reported. Use a "
                "statistical (fidelity-1) coupling model for the fade."
            )
            fade = None

        received = None
        margin = None
        if self.tx_power_dbm is not None:
            received = self.tx_power_dbm - samples
            if self.rx_sensitivity_dbm is not None and fade is not None:
                margin = {a: (self.tx_power_dbm - fade[a]) - self.rx_sensitivity_dbm
                          for a in availabilities}

        return {
            "total_loss_samples": samples,
            "mean_loss_db": samples.mean(axis=0),
            "fade_db": fade,
            "fade_available": provides_fade,
            "received_dbm": received,
            "margin_db": margin,
        }

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
from typing import Callable, Optional

import numpy as np
import pandas as pd

from .assumptions import Assumptions


@dataclass
class Term:
    '''One contribution to the link budget (loss positive, gain negative, dB).'''
    name: str
    category: str                       # geometric | atmospheric | turbulence | pointing | system
    mean_db: object = 0.0               # float or np.ndarray: deterministic value / E[loss]
    sampler: Optional[Callable[[int, np.random.Generator], np.ndarray]] = None
    quantile: Optional[Callable[[float], object]] = None
    note: str = ""
    meta: dict = field(default_factory=dict)
    assumptions: Optional[Assumptions] = None   # the regime the model is valid in

    @property
    def stochastic(self) -> bool:
        return self.sampler is not None

    def sample_db(self, n, rng):
        '''n Monte Carlo draws, shape (n, *base_shape). Deterministic terms broadcast.'''
        if self.sampler is not None:
            return np.asarray(self.sampler(n, rng))
        base = np.shape(self.mean_db)
        return np.broadcast_to(np.asarray(self.mean_db), (n, *base)).copy()

    def quantile_db(self, p):
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
                One (term name, reason) pair for each violation.
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
                If any term is Monte-Carlo-only (no closed-form quantile).
                Evaluate those terms with monte_carlo().
        '''
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

        fade = {a: np.percentile(samples, a * 100.0, axis=0)
                for a in availabilities}

        received = None
        margin = None
        if self.tx_power_dbm is not None:
            received = self.tx_power_dbm - samples
            if self.rx_sensitivity_dbm is not None:
                margin = {a: (self.tx_power_dbm - fade[a]) - self.rx_sensitivity_dbm
                          for a in availabilities}

        return {
            "total_loss_samples": samples,
            "mean_loss_db": samples.mean(axis=0),
            "fade_db": fade,
            "received_dbm": received,
            "margin_db": margin,
        }

'''
Turn an irradiance distribution model into a link-budget Term.

olb/turbulence/andrews/distributions.py gives irradiance models on the
NORMALISED irradiance I, with E[I] = 1. It holds no decibels. This module is the
one place that converts such a model into the three faces of a Term.

The conversion. The received power is proportional to I, and the loss in dB is
    loss_db = -10 log10(I).
So the three faces of the Term are:
    mean_db     = -(10 / ln 10) E[ln I]
    quantile(p) = -10 log10( I(1 - p) )
    sampler     = -10 log10( draws of I )

Two points about the mean face. First, the budget adds Terms in dB, so the mean
of the dB loss is the correct quantity to add, and that is -(10/ln10) E[ln I],
not -10 log10(E[I]). Second, E[ln I] is negative for every unit-mean model, so
mean_db is a positive loss. That agrees with the olb convention: loss is
positive dB.

Availability. The Term contract asks quantile(p) for the loss at availability p.
The link is at or better than that loss a fraction p of the time. A deeper loss
is a SMALLER irradiance, so the loss at availability p reads the (1 - p)
quantile of the irradiance. That is the 1 - p in the code.
'''

import numpy as np

from ..results import Term

_LN10 = np.log(10.0)


def irradiance_fade_term(name, category, *, mean_log, quantile, rvs, note="",
                         assumptions=None, meta=None) -> Term:
    '''
    Build a Term from a normalised-irradiance model with E[I] = 1.

    Parameters:
        name : str
            The Term name, for example "scintillation".
        category : str
            The Term category: geometric | atmospheric | turbulence | pointing |
            system | coupling.
        mean_log : float or ndarray
            E[ln I] of the model. Get it from the model's mean_log face, for
            example olb.turbulence.andrews.lognormal_mean_log(sigma_l2).
        quantile : callable
            f(p) -> the p-quantile of I. Bind the model parameters before you
            pass it, for example
            lambda p: lognormal_quantile(p, sigma_l2).
        rvs : callable
            f(n, rng) -> n draws of I. Bind the model parameters the same way.
        note : str
            One short line for the itemised budget table.
        assumptions : Assumptions, optional
            The regime the model is valid in. See olb.assumptions.
        meta : dict, optional
            Extra values to carry on the Term, for example the scintillation
            index and the model name.

    Returns:
        Term
            With mean_db, quantile(p) and sampler(n, rng) all set. The Term is
            NOT mean-only, because all three faces are real.
    '''
    mean_db = -(10.0 / _LN10) * np.asarray(mean_log, dtype=float)
    if np.ndim(mean_db) == 0:
        mean_db = float(mean_db)

    def _quantile(p):
        # A loss at availability p is the (1 - p) quantile of the irradiance.
        return -10.0 * np.log10(quantile(1.0 - p))

    def _sampler(n, rng):
        return -10.0 * np.log10(rvs(n, rng))

    return Term(
        name=name,
        category=category,
        mean_db=mean_db,
        sampler=_sampler,
        quantile=_quantile,
        note=note,
        meta=dict(meta) if meta else {},
        assumptions=assumptions,
    )


if __name__ == '__main__':
    from ..turbulence.andrews import (gamma_gamma_mean_log, gamma_gamma_params,
                                      gamma_gamma_quantile, gamma_gamma_rvs,
                                      lognormal_mean_log, lognormal_params,
                                      lognormal_quantile, lognormal_rvs)

    # === physics ============================================================

    # A model with no fluctuation gives no loss and no fade.
    flat = irradiance_fade_term(
        "flat", "turbulence",
        mean_log=0.0,
        quantile=lambda p: np.ones_like(np.asarray(p, dtype=float)),
        rvs=lambda n, rng: np.ones(n),
    )
    print(f"[physics] zero-fluctuation mean err = {abs(flat.mean_db):.3e}")
    print(f"[physics] zero-fluctuation q99  err = "
          f"{abs(flat.quantile_db(0.99)):.3e}")
    assert flat.mean_db == 0.0 and flat.quantile_db(0.99) == 0.0

    # The fade at high availability is deeper than the mean loss, and the fade
    # at low availability is shallower. That is the sign convention.
    sl2 = lognormal_params(0.2)
    ln_term = irradiance_fade_term(
        "lognormal", "turbulence",
        mean_log=lognormal_mean_log(sl2),
        quantile=lambda p: lognormal_quantile(p, sl2),
        rvs=lambda n, rng: lognormal_rvs(n, sl2, rng),
        note="test",
        meta={"model": "lognormal"},
    )
    q99, q50, q01 = (ln_term.quantile_db(0.99), ln_term.quantile_db(0.5),
                     ln_term.quantile_db(0.01))
    print(f"[physics] loss order q01<q50<q99: {q01:.4f} < {q50:.4f} < {q99:.4f}")
    assert q01 < q50 < q99
    assert ln_term.mean_db > 0.0 and not ln_term.mean_only

    # The sampled mean of the dB loss matches mean_db.
    rng = np.random.default_rng(7)
    draws = ln_term.sample_db(400_000, rng)
    print(f"[physics] sampled mean          err = "
          f"{abs(draws.mean() - ln_term.mean_db):.4f} dB")
    assert abs(draws.mean() - ln_term.mean_db) < 0.01

    # The sampled 99% loss matches the analytic quantile.
    p99 = float(np.percentile(draws, 99.0))
    print(f"[physics] sampled 99% fade      err = {abs(p99 - q99):.4f} dB")
    assert abs(p99 - q99) < 0.02

    # The gamma-gamma model works through the same adapter.
    s = 0.5 * np.log(1.0 + 0.2)
    aa, bb = gamma_gamma_params(s, s)
    gg_term = irradiance_fade_term(
        "gamma-gamma", "turbulence",
        mean_log=gamma_gamma_mean_log(aa, bb),
        quantile=lambda p: gamma_gamma_quantile(p, aa, bb),
        rvs=lambda n, rng: gamma_gamma_rvs(n, aa, bb, rng),
    )
    gg_draws = gg_term.sample_db(400_000, rng)
    print(f"[physics] gamma-gamma mean      err = "
          f"{abs(gg_draws.mean() - gg_term.mean_db):.4f} dB")
    assert abs(gg_draws.mean() - gg_term.mean_db) < 0.01

    # === reduction ==========================================================
    # Byte parity with the RETIRED inline lognormal faces. olb/links/downlink.py
    # and olb/links/terrestrial.py once built the three dB faces inline; they now
    # route through this adapter (backlog I-2 / crosscheck TL-01..04). Rebuild the
    # retired formula here and assert the downlink Term (which now uses the
    # adapter) reproduces it, so a future change to the adapter cannot silently
    # move the budget numbers. The mean and the sampler are BYTE identical; the
    # quantile matches to machine precision (~1e-16 dB), because the adapter takes
    # -10 log10(exp(x)) where the retired code took -10 x / ln10 directly.
    # Scenario: a 0.7 m ground aperture at 1550 nm, a 600 km orbit, 30 deg.
    from scipy.stats import norm
    from ..scenario import SpaceScenario, Channel
    from ..geometry import CircularOrbit
    from ..terminal import Terminal, Transmitter
    from ..links.downlink import downlink_scintillation_term
    from ..turbulence.profiles import DEFAULT_HS, default_cn2_profile

    lam = 1550e-9
    space = Terminal(aperture_m=0.05, wavelength_m=lam,
                     transmitter=Transmitter(waist_m=0.035))
    ground = Terminal(aperture_m=0.7, wavelength_m=lam)
    scenario = SpaceScenario(ground=ground, space=space, direction="downlink",
                             channel=Channel(altitude_m=600e3))
    cn2 = default_cn2_profile(scenario.channel.site, DEFAULT_HS)
    term = downlink_scintillation_term(scenario, CircularOrbit(600e3, 30.0),
                                       cn2_profile=cn2)

    # The retired inline faces, verbatim from the old downlink/terrestrial code.
    sigma2_P = term.meta["sigma2_P"]
    sl2 = np.log(1.0 + sigma2_P)
    sl = np.sqrt(sl2)
    old_mean = (5.0 / _LN10) * sl2
    old_quant = lambda p: -10.0 / _LN10 * (-sl2 / 2.0 + sl * norm.ppf(1.0 - p))
    old_sampler = lambda n, rng: -10.0 * np.log10(
        rng.lognormal(mean=-sl2 / 2.0, sigma=sl, size=n))

    d_mean = abs(term.mean_db - old_mean)
    print(f"[reduce ] downlink Term vs retired inline mean err = {d_mean:.3e} dB "
          f"({term.mean_db:.6f} dB)")
    assert d_mean == 0.0

    for p in (0.01, 0.99):
        d_q = abs(term.quantile_db(p) - old_quant(p))
        print(f"[reduce ] quantile({p}) err = {d_q:.3e} dB "
              f"({term.quantile_db(p):.6f} dB)")
        assert d_q < 1e-12

    # The sampler agrees draw for draw with the retired one at the same seed.
    a = term.sample_db(50_000, np.random.default_rng(1234))
    b = old_sampler(50_000, np.random.default_rng(1234))
    print(f"[reduce ] sampler err = {np.max(np.abs(a - b)):.3e} dB "
          f"(n=50000, same seed)")
    assert np.max(np.abs(a - b)) == 0.0
    assert a.shape == (50_000,)

    # WP3a: the adapter carries a TRACED assumptions record through unchanged.
    # links/downlink._gamma_gamma_term hands this adapter a merged, provenance-
    # bearing record (the factory owns the merge; the adapter builds none of its
    # own). Prove the pass-through keeps the provenance and the constraints.
    from ..assumptions import (Assumptions, Constraint, BEAM_PLANE_WAVE,
                               REGIME_STRONG, SPECTRUM_KOLMOGOROV)
    traced_rec = Assumptions(
        beam_type=BEAM_PLANE_WAVE, turbulence_regime=REGIME_STRONG,
        spectrum=SPECTRUM_KOLMOGOROV, validity="traced demo",
        constraints=[("olb.demo.phys", Constraint(
            "regime", "Weak fluctuation: sigma_R^2 < 1.", "10.1117/3.626196",
            "Ch. 8, printed p. 264"))],
        provenance=["olb.demo.phys"])
    carried = irradiance_fade_term(
        "scintillation", "turbulence", mean_log=gamma_gamma_mean_log(aa, bb),
        quantile=lambda p: gamma_gamma_quantile(p, aa, bb),
        rvs=lambda n, rng: gamma_gamma_rvs(n, aa, bb, rng),
        assumptions=traced_rec)
    assert carried.assumptions is traced_rec
    assert carried.assumptions.provenance == ["olb.demo.phys"]
    assert len(carried.assumptions.constraints) == 1
    print(f"[reduce ] traced record carried through: provenance "
          f"{carried.assumptions.provenance}")

    print("self-check passed")

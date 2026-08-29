'''
Refractive-index power spectrum models of Andrews and Phillips.

This module gives the three-dimensional spatial power spectral density Phi_n of
the refractive-index fluctuations. Every other Andrews module builds its integral
on one of these five models.

Source of every equation:
    L. C. Andrews and R. L. Phillips, "Laser Beam Propagation through Random
    Media", 2nd ed., SPIE Press (2005). DOI: 10.1117/3.626196
Chapter 3, Sec. 3.3, printed pp. 66 to 72. Each function names its section, its
equation number, and its printed page.

The five models:

- `kolmogorov` (Eq. (18)) is the plain inertial-range power law. It has no inner
  scale and no outer scale.
- `tatarskii` (Eq. (19)) adds a Gaussian cut at the inner scale.
- `von_karman` (Eq. (20)) adds an outer scale. With an inner scale it becomes the
  "modified von Karman" spectrum, which is the lower expression of Eq. (20).
- `exponential` (Eq. (21)) is a second outer-scale model. Its outer-scale
  wavenumber carries a scaling constant C0 that the application chooses.
- `modified_atmospheric` (Eq. (22), with the alternative outer-scale form
  Eq. (23)) is the Andrews-Hill spectrum. It is the only model in the list that
  has the high-wavenumber bump. The book states the bump changes the
  scintillation index.

UNITS. `kappa` is the scalar spatial wavenumber [rad/m]. `cn2` is the
refractive-index structure constant [m^-2/3]. `l0` and `L0` are the inner scale
and the outer scale [m]. The result Phi_n has units of m^3.

OUTER-SCALE WAVENUMBER. The book is explicit that the choice is not unique. It
prints k0 = 2*pi/L0 at Eq. (20) and again at Eq. (22), with "or sometimes
k0 = 1/L0". It prints k0 = 4*pi/L0 at Eq. (23), "where, in some cases, we might
instead define k0 = 2*pi/L0 or k0 = 8*pi/L0". The Ch. 9 scintillation model uses
C0 = 8*pi (text below Eq. (21), printed p. 68). Each function below states the
constant it uses and lets the caller change it.

This module holds physics only. It returns no decibels. Each builder declares
its spectrum model through the @assumes decorator (see olb.assumptions), so a
Term that runs the builder inside a trace_assumptions() context inherits the
spectrum limit automatically.
'''

import numpy as np

from ...assumptions import (
    Constraint, assumes,
    SPECTRUM_KOLMOGOROV, SPECTRUM_TATARSKII, SPECTRUM_VON_KARMAN,
    SPECTRUM_EXPONENTIAL, SPECTRUM_MODIFIED,
)

# The Kolmogorov constant of the refractive-index spectrum. Source: Andrews and
# Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 3, Eq. (18), printed p. 67.
KOLMOGOROV_CONSTANT = 0.033

# Inner-scale wavenumber constants. Source: Ch. 3, Eq. (19), printed p. 67
# (km = 5.92/l0) and Ch. 3, Eq. (22), printed p. 69 (kl = 3.3/l0).
TATARSKII_KM = 5.92
MODIFIED_KL = 3.3

# Default outer-scale scaling constant, k0 = C0/L0. Source: Ch. 3, Eq. (20),
# printed p. 68, and Eq. (22), printed p. 69.
VON_KARMAN_C0 = 2.0 * np.pi

# Outer-scale scaling constant of the exponential spectrum. Source: Ch. 3, text
# below Eq. (21), printed p. 68: C0 = 4*pi approximates the von Karman spectrum,
# and the Ch. 9 scintillation model uses C0 = 8*pi.
EXPONENTIAL_C0 = 4.0 * np.pi

# Outer-scale scaling constant of the alternative modified form, Eq. (23).
# Source: Ch. 3, Eq. (23), printed p. 69.
MODIFIED_EQ23_C0 = 4.0 * np.pi


# ---------------------------------------------------------------------------
# Function-owned validity limits (see olb.assumptions).
#
# Each builder carries ONE spectrum Constraint that names the model and states
# which scales the model holds. The `_refuse`/`_need` raises still enforce a
# missing or a forbidden scale; the Constraint records the same limit for a
# reader and for Budget.check(). No numeric run-time check applies here: a
# spectrum is a modelling choice, not a regime bound.
# ---------------------------------------------------------------------------
_DOI = "10.1117/3.626196"

_KOLMOGOROV_SPECTRUM = Constraint(
    "spectrum",
    "This is the Kolmogorov spectrum. It has no inner scale and no outer "
    "scale. It holds only over the inertial subrange 1/L0 << kappa << 1/l0. An "
    "extension to all wavenumbers makes some integrals diverge.",
    _DOI, "Ch. 3, Eq. (18), printed p. 67")

_TATARSKII_SPECTRUM = Constraint(
    "spectrum",
    "This is the Tatarskii spectrum. It adds a Gaussian cut at the inner scale "
    "l0. It has no outer scale, so it keeps the singularity at kappa = 0.",
    _DOI, "Ch. 3, Eq. (19), printed p. 67")

_VON_KARMAN_SPECTRUM = Constraint(
    "spectrum",
    "This is the von Karman spectrum. It carries the outer scale L0. With l0 "
    "set it becomes the modified von Karman spectrum and carries the inner "
    "scale too. The outer-scale constant k0 = C0/L0 belongs to the "
    "application.",
    _DOI, "Ch. 3, Eq. (20), printed p. 68")

_EXPONENTIAL_SPECTRUM = Constraint(
    "spectrum",
    "This is the exponential outer-scale spectrum. It carries the outer scale "
    "L0 and cuts the low wavenumbers. It has no inner scale.",
    _DOI, "Ch. 3, Eq. (21), printed p. 68")

_MODIFIED_SPECTRUM = Constraint(
    "spectrum",
    "This is the modified atmospheric (Andrews-Hill) spectrum. It carries the "
    "high-wavenumber bump and needs the inner scale l0. The outer scale L0 is "
    "optional.",
    _DOI, "Ch. 3, Eqs. (22) and (23), printed p. 69")


def _base(kappa, cn2):
    '''Return the two inputs as float arrays.'''
    return (np.asarray(kappa, dtype=float), np.asarray(cn2, dtype=float))


def _need(value, name, func):
    '''Refuse a missing scale.'''
    if value is None:
        raise ValueError(f'{func} needs {name}')
    return float(value)


def _refuse(l0, L0, func, allowed):
    '''Refuse a scale that the named model does not carry.'''
    if l0 is not None and 'l0' not in allowed:
        raise ValueError(f'{func} has no inner scale; l0 must be None')
    if L0 is not None and 'L0' not in allowed:
        raise ValueError(f'{func} has no outer scale; L0 must be None')


@assumes(_KOLMOGOROV_SPECTRUM, spectrum=SPECTRUM_KOLMOGOROV)
def kolmogorov(kappa, cn2, l0=None, L0=None):
    '''
    Return the Kolmogorov spectrum Phi_n(kappa).

    formula:
        Phi_n = 0.033 Cn2 kappa^(-11/3),   1/L0 << kappa << 1/l0
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 3,
    Eq. (18), printed p. 67. Ch. 6, Eq. (28), printed p. 187, restates it.

    VALIDITY. The book states the model holds only over the inertial subrange.
    To use it over all wavenumbers, the caller must accept L0 = infinity and
    l0 = 0. The book warns that this extension makes some integrals diverge
    (Ch. 3, text below Eq. (18), printed p. 67).

    This model has no inner scale and no outer scale. A value for l0 or L0 is
    refused, not ignored.
    '''
    _refuse(l0, L0, 'kolmogorov', allowed=())
    kappa, cn2 = _base(kappa, cn2)
    return KOLMOGOROV_CONSTANT * cn2 * kappa ** (-11.0 / 3.0)


@assumes(_TATARSKII_SPECTRUM, spectrum=SPECTRUM_TATARSKII)
def tatarskii(kappa, cn2, l0=None, L0=None):
    '''
    Return the Tatarskii spectrum Phi_n(kappa) with a finite inner scale.

    formula:
        Phi_n = 0.033 Cn2 kappa^(-11/3) exp(-kappa^2 / km^2),
        km = 5.92 / l0,   kappa >> 1/L0
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 3,
    Eq. (19), printed p. 67. Ch. 6, Eq. (29), printed p. 187, restates it.

    The constant 5.92 makes the structure function take its quadratic form below
    the inner scale (Ch. 3, footnote 3, printed p. 67).

    This model has no outer scale. It keeps the kappa = 0 singularity of the
    Kolmogorov model, so its covariance function does not exist (Ch. 3, text
    below Eq. (19), printed p. 67).
    '''
    _refuse(None, L0, 'tatarskii', allowed=('l0',))
    l0 = _need(l0, 'l0', 'tatarskii')
    kappa, cn2 = _base(kappa, cn2)
    km = TATARSKII_KM / l0
    return (KOLMOGOROV_CONSTANT * cn2 * kappa ** (-11.0 / 3.0)
            * np.exp(-(kappa / km) ** 2))


@assumes(_VON_KARMAN_SPECTRUM, spectrum=SPECTRUM_VON_KARMAN)
def von_karman(kappa, cn2, l0=None, L0=None, *, c0=VON_KARMAN_C0):
    '''
    Return the von Karman spectrum Phi_n(kappa).

    With `l0` set to None the result is the plain von Karman spectrum, which
    carries the outer scale only. With `l0` set the result is the modified von
    Karman spectrum, which carries both scales.

    formula:
        Phi_n = 0.033 Cn2 / (kappa^2 + k0^2)^(11/6),          l0 is None
        Phi_n = 0.033 Cn2 exp(-kappa^2/km^2)
                / (kappa^2 + k0^2)^(11/6),                    l0 is set
        km = 5.92 / l0,   k0 = C0 / L0
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 3,
    Eq. (20), printed p. 68. Ch. 6, Eq. (28), printed p. 187, and Ch. 8, Eq. (25),
    printed p. 265, restate it.

    OUTER-SCALE CONSTANT. The book prints "k0 = 2*pi/L0 (or sometimes
    k0 = 1/L0)" at Eq. (20), printed p. 68. This function uses C0 = 2*pi by
    default. Note that the book uses C0 = 4*pi at Ch. 3, Eq. (23), printed p. 69,
    for the alternative outer-scale form of the MODIFIED spectrum, and C0 = 8*pi
    in the Ch. 9 scintillation model. So the constant belongs to the application,
    not to the spectrum.

    In the inertial subrange k0 << kappa << km the result goes to the Kolmogorov
    spectrum (Ch. 3, text below Eq. (20), printed p. 68).
    '''
    L0 = _need(L0, 'L0', 'von_karman')
    kappa, cn2 = _base(kappa, cn2)
    k0 = c0 / L0
    out = KOLMOGOROV_CONSTANT * cn2 / (kappa ** 2 + k0 ** 2) ** (11.0 / 6.0)
    if l0 is None:
        return out
    km = TATARSKII_KM / float(l0)
    return out * np.exp(-(kappa / km) ** 2)


@assumes(_EXPONENTIAL_SPECTRUM, spectrum=SPECTRUM_EXPONENTIAL)
def exponential(kappa, cn2, l0=None, L0=None, *, c0=EXPONENTIAL_C0):
    '''
    Return the exponential outer-scale spectrum Phi_n(kappa).

    formula:
        Phi_n = 0.033 Cn2 kappa^(-11/3) [1 - exp(-kappa^2 / k0^2)],
        k0 = C0 / L0,   0 <= kappa << 1/l0
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 3,
    Eq. (21), printed p. 68.

    OUTER-SCALE CONSTANT. The book states the scaling constant depends on the
    application: C0 = 4*pi approximates the von Karman spectrum, and the Ch. 9
    scintillation model uses C0 = 8*pi (Ch. 3, text below Eq. (21), printed
    p. 68). This function uses C0 = 4*pi by default.

    This model has no inner scale. A value for l0 is refused, not ignored.
    '''
    _refuse(l0, None, 'exponential', allowed=('L0',))
    L0 = _need(L0, 'L0', 'exponential')
    kappa, cn2 = _base(kappa, cn2)
    k0 = c0 / L0
    return (KOLMOGOROV_CONSTANT * cn2 * kappa ** (-11.0 / 3.0)
            * (1.0 - np.exp(-(kappa / k0) ** 2)))


@assumes(_MODIFIED_SPECTRUM, spectrum=SPECTRUM_MODIFIED)
def modified_atmospheric(kappa, cn2, l0=None, L0=None, *, outer='karman',
                         c0=None):
    '''
    Return the modified atmospheric (Andrews-Hill) spectrum Phi_n(kappa).

    This is the only model in this module that carries the high-wavenumber bump.
    The book states the bump appears in the measured temperature spectrum, and
    that it puts a matching bump into the structure function and the
    scintillation index (Ch. 3, Sec. 3.3.3, printed pp. 68 to 70).

    Parameters:
        kappa : float or numpy.ndarray
            Scalar spatial wavenumber [rad/m].
        cn2 : float or numpy.ndarray
            Refractive-index structure constant [m^-2/3].
        l0 : float
            Inner scale [m]. Required.
        L0 : float, optional
            Outer scale [m]. None gives an infinite outer scale.
        outer : str
            "karman" gives the Eq. (22) outer-scale form, which divides by
            (kappa^2 + k0^2)^(11/6) with k0 = 2*pi/L0. "exponential" gives the
            Eq. (23) form, which multiplies by [1 - exp(-kappa^2/k0^2)] with
            k0 = 4*pi/L0.
        c0 : float, optional
            Outer-scale scaling constant, k0 = C0/L0. None takes the default of
            the chosen form.

    formula (Eq. (22), outer="karman"):
        Phi_n = 0.033 Cn2 [1 + 1.802 (kappa/kl) - 0.254 (kappa/kl)^(7/6)]
                exp(-kappa^2/kl^2) / (kappa^2 + k0^2)^(11/6),
        kl = 3.3 / l0,   k0 = 2*pi/L0
    formula (Eq. (23), outer="exponential"):
        Phi_n = 0.033 Cn2 [1 + 1.802 (kappa/kl) - 0.254 (kappa/kl)^(7/6)]
                [1 - exp(-kappa^2/k0^2)] exp(-kappa^2/kl^2) kappa^(-11/3),
        kl = 3.3 / l0,   k0 = 4*pi/L0
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 3,
    Eqs. (22) and (23), printed p. 69. Ch. 6, Eq. (31), printed p. 187, restates
    Eq. (22). Ch. 9, Eqs. (3), (5) and (6), printed pp. 327 and 328, write the
    same two filters as f(kappa l0) and g(kappa L0), and there the book uses
    kappa_0 = 8*pi/L0.
    '''
    l0 = _need(l0, 'l0', 'modified_atmospheric')
    kappa, cn2 = _base(kappa, cn2)
    kl = MODIFIED_KL / l0
    ratio = kappa / kl
    bump = 1.0 + 1.802 * ratio - 0.254 * ratio ** (7.0 / 6.0)
    core = KOLMOGOROV_CONSTANT * cn2 * bump * np.exp(-ratio ** 2)

    if outer == 'karman':
        k0 = 0.0 if L0 is None else (VON_KARMAN_C0 if c0 is None else c0) / L0
        return core / (kappa ** 2 + k0 ** 2) ** (11.0 / 6.0)
    if outer == 'exponential':
        core = core * kappa ** (-11.0 / 3.0)
        if L0 is None:
            return core
        k0 = (MODIFIED_EQ23_C0 if c0 is None else c0) / L0
        return core * (1.0 - np.exp(-(kappa / k0) ** 2))
    raise ValueError(f'outer must be "karman" or "exponential", not {outer!r}')


# The five models, keyed by a short name. A plain dict, not a registry.
SPECTRA = {
    'kolmogorov': kolmogorov,
    'tatarskii': tatarskii,
    'von_karman': von_karman,
    'exponential': exponential,
    'modified': modified_atmospheric,
}


if __name__ == '__main__':
    # ---------------- physics self-checks ----------------
    cn2_ref = 1e-14
    l0_ref = 5e-3
    L0_ref = 10.0

    # Every model gives a positive spectrum over the inertial range.
    inertial = np.logspace(0.0, 2.0, 50)     # 1 to 100 rad/m
    for name, func in SPECTRA.items():
        kwargs = {}
        if name in ('tatarskii', 'modified'):
            kwargs['l0'] = l0_ref
        if name in ('von_karman', 'exponential'):
            kwargs['L0'] = L0_ref
        value = func(inertial, cn2_ref, **kwargs)
        assert np.all(value > 0.0), name
        assert np.all(np.isfinite(value)), name

    # The spectrum falls as the wavenumber grows.
    kol = kolmogorov(inertial, cn2_ref)
    assert np.all(np.diff(kol) < 0.0)

    # The Kolmogorov slope is -11/3 on a log-log plot.
    slope = np.polyfit(np.log(inertial), np.log(kol), 1)[0]
    assert abs(slope + 11.0 / 3.0) < 1e-10, slope

    # The Tatarskii cut kills the spectrum above km = 5.92/l0.
    km = TATARSKII_KM / l0_ref
    assert tatarskii(3.0 * km, cn2_ref, l0=l0_ref) < 1e-3 * kolmogorov(
        3.0 * km, cn2_ref)

    # The outer scale flattens the von Karman spectrum below k0 = 2*pi/L0.
    k0 = VON_KARMAN_C0 / L0_ref
    low = np.array([0.01 * k0, 0.02 * k0])
    vk_low = von_karman(low, cn2_ref, L0=L0_ref)
    assert abs(vk_low[1] / vk_low[0] - 1.0) < 1e-3, vk_low

    # The exponential model also flattens, and it goes to zero at kappa = 0.
    assert exponential(1e-9, cn2_ref, L0=L0_ref) < exponential(
        1.0, cn2_ref, L0=L0_ref)

    # The modified spectrum carries the bump. Scale it by the Kolmogorov law and
    # look near kappa*l0 ~ 1, where Ch. 3, Fig. 3.6, printed p. 70, shows the
    # rise above 1.
    kb = np.logspace(-1.0, 1.0, 400) / l0_ref
    scaled = modified_atmospheric(kb, cn2_ref, l0=l0_ref) / kolmogorov(kb,
                                                                      cn2_ref)
    assert scaled.max() > 1.1, scaled.max()
    bump_at = kb[np.argmax(scaled)] * l0_ref
    assert 0.2 < bump_at < 2.0, bump_at

    # A model refuses a scale it does not carry.
    for call in ((kolmogorov, {'l0': l0_ref}), (kolmogorov, {'L0': L0_ref}),
                 (tatarskii, {'l0': l0_ref, 'L0': L0_ref}),
                 (exponential, {'l0': l0_ref, 'L0': L0_ref})):
        try:
            call[0](inertial, cn2_ref, **call[1])
        except ValueError:
            pass
        else:
            raise AssertionError(f'{call[0].__name__} must refuse {call[1]}')

    # ---------------- REDUCTION checks ----------------
    # 1. kolmogorov equals the literal 0.033 Cn2 kappa^(-11/3) that the parent
    # module `olb.turbulence.plane_wave_scintillation` writes inline inside
    # `_scintillation_integral`. That module holds the only other copy in olb.
    from .. import plane_wave_scintillation as pws

    parent_shape = 0.033 * cn2_ref * pws._KAPPA ** (-11.0 / 3.0)
    mine_shape = kolmogorov(pws._KAPPA, cn2_ref)
    err_kol = float(np.max(np.abs(mine_shape / parent_shape - 1.0)))
    assert err_kol < 1e-12, err_kol
    print(f'REDUCTION kolmogorov vs the inline parent copy : '
          f'max rel err = {err_kol:.3e}  (target 1e-12)')

    # 2. von_karman with a vanishing inner scale and an unbounded outer scale
    # goes to the Kolmogorov spectrum over the inertial range.
    tiny_l0, huge_L0 = 1e-9, 1e9
    vk = von_karman(inertial, cn2_ref, l0=tiny_l0, L0=huge_L0)
    err_vk = float(np.max(np.abs(vk / kol - 1.0)))
    assert err_vk < 1e-6, err_vk
    print(f'REDUCTION von_karman(l0->0, L0->inf) -> kolmogorov : '
          f'max rel err = {err_vk:.3e}  (target 1e-6)')

    # 3. The same limit for the other two scale models.
    tat = tatarskii(inertial, cn2_ref, l0=tiny_l0)
    err_tat = float(np.max(np.abs(tat / kol - 1.0)))
    expo = exponential(inertial, cn2_ref, L0=huge_L0)
    err_exp = float(np.max(np.abs(expo / kol - 1.0)))
    assert err_tat < 1e-6 and err_exp < 1e-6, (err_tat, err_exp)
    print(f'REDUCTION tatarskii(l0->0) = {err_tat:.3e}   '
          f'exponential(L0->inf) = {err_exp:.3e}  (target 1e-6)')

    # 4. The modified spectrum keeps the Kolmogorov law far below the bump.
    tiny = np.logspace(-4.0, -3.0, 20) / l0_ref
    mod = modified_atmospheric(tiny, cn2_ref, l0=l0_ref)
    err_mod = float(np.max(np.abs(mod / kolmogorov(tiny, cn2_ref) - 1.0)))
    print(f'REDUCTION modified -> kolmogorov at kappa l0 ~ 1e-3 : '
          f'max rel err = {err_mod:.3e}  (no target, the bump term is linear '
          f'in kappa l0)')

    # ---------------- assumption-trace checks ----------------
    import warnings

    from ...assumptions import trace_assumptions

    # (1) Value parity: the return value is byte-identical with and without a
    # collection context.
    kol_plain = float(kolmogorov(10.0, cn2_ref))
    with trace_assumptions():
        kol_traced = float(kolmogorov(10.0, cn2_ref))
    assert kol_plain == kol_traced, (kol_plain, kol_traced)

    # (2) Registration: every builder registers its source and its spectrum
    # Constraint, and its headline spectrum.
    with trace_assumptions() as tr:
        for name, func in SPECTRA.items():
            kwargs = {}
            if name in ('tatarskii', 'modified'):
                kwargs['l0'] = l0_ref
            if name in ('von_karman', 'exponential'):
                kwargs['L0'] = L0_ref
            func(inertial, cn2_ref, **kwargs)
    sources = set(tr.records)
    for fn in ('kolmogorov', 'tatarskii', 'von_karman', 'exponential',
               'modified_atmospheric'):
        assert any(s.endswith('.' + fn) for s in sources), fn
    kinds = {c.kind for r in tr.records.values() for c in r.constraints}
    assert kinds == {'spectrum'}, kinds
    specs = {r.spectrum for r in tr.records.values()}
    assert {SPECTRUM_KOLMOGOROV, SPECTRUM_TATARSKII, SPECTRUM_VON_KARMAN,
            SPECTRUM_EXPONENTIAL, SPECTRUM_MODIFIED} == specs, specs

    # (3) The module carries no numeric run-time check (a spectrum is a
    # modelling choice, not a regime bound), so a refused scale still RAISES,
    # and the physics layer emits no warning inside a context.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        with trace_assumptions():
            kolmogorov(inertial, cn2_ref)
            try:
                kolmogorov(inertial, cn2_ref, l0=l0_ref)
            except ValueError:
                pass
            else:
                raise AssertionError('kolmogorov must still refuse l0 in a context')
    assert len(caught) == 0, 'the physics layer must not warn'
    print(f'assumption trace : {len(sources)} sources, kinds {sorted(kinds)}')

    print(f'Phi_n(kappa=10, Cn2=1e-14) = {kolmogorov(10.0, cn2_ref):.4e} m^3')
    print('self-check passed')

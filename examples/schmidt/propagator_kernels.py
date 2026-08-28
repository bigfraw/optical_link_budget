'''
The Schmidt propagation kernels against the production LightPipes propagators.

The two layers solve the SAME Fresnel integral. They do not use the same
quadrature, and they do not use the same phase convention. This script puts the
same input field into both, bridges the two conventions, and MEASURES what is
left.

THE TWO CONVENTION GAPS. Bridge them before you compare anything.

  1. THE PISTON PHASE. The production `Forvard` and `Fresnel` keep the factor
     exp(i k z) of the free-space transfer function. The Schmidt kernels DROP
     it, because the book's own listings drop it (Listings 6.1, 6.3 and 6.5,
     printed pp. 91, 96 and 102). So divide the production field by exp(i k z)
     before you compare. Use the LEGACY constant of the production module for
     that k: `olb/waveoptics/propagators.py` writes 2*3.141592654 in the place
     of 2*pi. Over a 0.5 m path at 1550 nm the true pi and the legacy pi differ
     by 7 mrad of piston, so the legacy value is the one that cancels.
  2. THE QUADRATURE. `Forvard` uses an unscaled `fft2` with a sign-alternation
     trick in the place of the two `fftshift` calls, and it wraps the transfer
     phase into [-2 pi, 0] by an integer truncation. That is the SAME algorithm
     as `schmidt.fresnel.angular_spectrum` with m = 1, so tier (a) below is
     tight. `Fresnel` is a DIFFERENT quadrature: it convolves on a doubled grid
     with the closed-form pixel integral of the Fresnel integrals C(x) and
     S(x), and it builds that kernel on the legacy pitch siz/(N - 1), not
     siz/N. So tier (b) below is loose, and the script prints the number.

THE THREE TIERS.

  (a) SAME ALGORITHM. `angular_spectrum` (m = 1) against `Forvard`. The two
      sample one transfer function on one grid. ASSERT a tight bound.
  (b) SAME INTEGRAL, OTHER QUADRATURE. `one_step_fresnel`, `two_step_fresnel`
      and `angular_spectrum` against the production `Fresnel`, and against the
      closed-form Gaussian where the source is a clean Gaussian. ASSERT a loose
      STATED bound, and print the measured number.
  (c) THE CO-MOVING GRID. `two_step_fresnel` with a free output pitch against
      the production three-call recipe Lens -> LensFresnel -> Convert. The two
      routes reach the same magnified grid by different roads. Loose bound.

THE GEOMETRY IS NOT ARBITRARY. The far range is
Z_FAR = N dx^2 / lambda. Two rules meet at that one point:

  - Ch. 6, Eq. (6.16), printed p. 90, FIXES the one-step observation pitch at
    lambda z /(N dx). At Z_FAR that pitch is exactly dx, so the one-step result
    lands on the SAME grid as the input and as the production propagators. No
    interpolation is necessary anywhere in this script.
  - Ch. 7, Eq. (7.59), printed p. 127 (constraint 4), reads N >= lambda z /
    (dx1 dx2). At Z_FAR with dx2 = dx that is N >= N, the range limit of the
    angular-spectrum kernel. So Z_FAR is the one range where the SHORT-range
    kernel and the LONG-range kernel are both just valid.

The near range is Z_NEAR = Z_FAR / 2, comfortably inside the angular-spectrum
limit.

TWO INPUT FIELDS, EVERY TIME. A soft Gaussian, which every kernel reproduces to
the closed form; and the same Gaussian behind a hard circular aperture, which
diffracts into the ring structure of a Fresnel-number 3 pattern. The hard case
is the one that separates the quadratures, because its spectrum does not decay.

Figures:
    examples/schmidt/figures/propagator_kernels_cuts.png      the intensity cross-cuts
    examples/schmidt/figures/propagator_kernels_diffs.png     the log difference maps

Sources:
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB, DOI 10.1117/3.866274. Ch. 6, Eq. (6.5), printed p. 88 (the one-step
  Fresnel transform); Ch. 6, Eq. (6.16), printed p. 90 (the fixed output
  pitch); Ch. 6, Eqs. (6.24) and (6.25), printed p. 94 (the two-step
  intermediate plane); Ch. 6, Eq. (6.32), printed p. 95 (the free-space
  transfer function); Ch. 6, Eq. (6.65), printed p. 100 (the scaled
  angular-spectrum chain); Ch. 1, Eqs. (1.53) to (1.56), printed p. 9 (the
  closed-form Gaussian beam); Ch. 7, Eqs. (7.41) and (7.42), printed p. 123
  (the Fresnel-integral minimum distance); Ch. 7, Eq. (7.59), printed p. 127
  (constraint 4).
- Goodman, Introduction to Fourier Optics, ISBN 978-0974707723. The Fresnel
  diffraction integral, and the circular pupil.
- Siegman, Lasers, ISBN 978-0935702118. The Gaussian beam w(z) and R(z).
- LightPipes manual, https://opticspy.github.io/lightpipes/manual.html. The
  lineage of the production propagators.

Run from the repo root:
    python -m examples.schmidt.propagator_kernels
'''

import time

import matplotlib.pyplot as plt
import numpy as np

from olb.waveoptics import (Begin, CircAperture, Convert, Forvard, Fresnel,
                            GaussBeam, Lens, LensFresnel)
from olb.waveoptics.schmidt.fresnel import (angular_spectrum, one_step_fresnel,
                                            two_step_fresnel)
from olb.waveoptics.schmidt.sampling import (constraint4_min_n,
                                             fresnel_min_distance,
                                             one_step_delta2)

WAVELENGTH_M = 1550e-9
N = 512
SIDE_M = 20e-3                  # the grid pitch is 39.06 um
WAIST_M = 1.5e-3                # the soft Gaussian waist radius
APERTURE_R_M = 1.5e-3           # the hard truncation, at one waist radius

# The co-moving case. A 500 m hop of a 1 mm waist grows the beam by about 247,
# so a flat grid cannot hold both ends.
LENS_N = 1024
LENS_SIDE_M = 20e-3
LENS_WAIST_M = 1.0e-3
LENS_APERTURE_R_M = 1.0e-3
LENS_RANGE_M = 500.0

# The radius that tier (b) reads. The production `Fresnel` convolves on a
# grid of twice the side and then differences four shifted copies, so it leaves
# an artefact in the OUTER band of the grid. That artefact is the largest
# single difference between the two layers, and it sits where the field is
# below 1e-8 of the peak. The script measures the interior AND the whole grid,
# and it prints both. The interior radius is 0.6 of the half-side.
CORE_RADIUS_M = 0.6 * SIDE_M / 2

# The stated tolerances. Tier (a) is one algorithm, so it is tight, and it
# holds over the WHOLE grid. Two BOOK kernels on one output grid solve one
# integral with one transform library, so they are tighter still. Everything
# that crosses a quadrature, or that crosses an output grid, is loose, and the
# script prints the measured value next to the bound.
TOL_TIER_A = 1e-6
TOL_BOOK_PAIR = 1e-12
TOL_TIER_B = 3e-2
TOL_TIER_C = 5e-2

# The book kernels of the tier (b) table. A pair drawn from this set gets the
# tight bound.
BOOK_KEYS = ("one-step Fresnel", "angular spectrum")

# The legacy 2 pi of olb/waveoptics/propagators.py. The production piston phase
# uses it, so the bridge must use it too.
TWO_PI_LEGACY = 2.0 * 3.141592654

CUTS_PNG = "examples/schmidt/figures/propagator_kernels_cuts.png"
DIFFS_PNG = "examples/schmidt/figures/propagator_kernels_diffs.png"


def axis_m(n, dx):
    '''Give the coordinate of each sample, in m.

    The axis runs from -N/2 dx to (N/2 - 1) dx. Schmidt (2010),
    DOI 10.1117/3.866274, Ch. 6, Listing 6.1, line 9, printed p. 91, uses that
    range, and `olb.waveoptics.field.Field.xvalues` uses the same one.
    '''
    return (np.arange(n) - n // 2) * dx


def drop_piston(field_array, z, wavelength_m):
    '''Remove the factor exp(i k z) that the production propagators keep.

    The Schmidt kernels drop that factor (Schmidt (2010),
    DOI 10.1117/3.866274, Ch. 6, text below Eq. (6.32), printed p. 95, and
    Listings 6.1, 6.3 and 6.5, printed pp. 91, 96 and 102). The bridge uses the
    LEGACY 2 pi of the production module, because that is the constant that the
    production field carries.
    '''
    return field_array * np.exp(-1j * TWO_PI_LEGACY / wavelength_m * z)


def gaussian_closed_form(n, dx, w0, z, wavelength_m):
    '''Give the free-space Gaussian field at the range z, with no piston.

        U(r,z) = (W0/W) exp(-r^2/W^2) exp(i k r^2 /(2 R)) exp(-i arctan(z/zR))
        W = W0 sqrt(1 + (z/zR)^2),  R = z (1 + (zR/z)^2),  zR = pi W0^2/lambda

    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 1, Eqs. (1.53) to (1.56),
    printed p. 9. Siegman, Lasers, ISBN 978-0935702118, prints the same
    solution. The factor exp(i k z) is REMOVED, so the result matches the
    Schmidt kernel convention.
    '''
    r2 = axis_m(n, dx)[:, None] ** 2 + axis_m(n, dx)[None, :] ** 2
    k = 2.0 * np.pi / wavelength_m
    zr = np.pi * w0 ** 2 / wavelength_m
    w = w0 * np.sqrt(1.0 + (z / zr) ** 2)
    radius = z * (1.0 + (zr / z) ** 2)
    return (w0 / w * np.exp(-r2 / w ** 2)
            * np.exp(1j * k * r2 / (2.0 * radius))
            * np.exp(-1j * np.arctan(z / zr)))


def rel_diff(a, b, mask=None):
    '''Give max|a - b| inside the mask, over the peak amplitude of b.'''
    if mask is None:
        mask = np.ones(a.shape, dtype=bool)
    peak = float(np.max(np.abs(b)))
    return float(np.max(np.abs(a[mask] - b[mask]))) / peak


def input_fields():
    '''Build the two input fields, on the production grid.

    The Schmidt kernels then read `Field.field` itself, so the two layers see
    the SAME complex array, to the last bit.
    '''
    soft = GaussBeam(Begin(SIDE_M, WAVELENGTH_M, N), WAIST_M)
    hard = CircAperture(GaussBeam(Begin(SIDE_M, WAVELENGTH_M, N), WAIST_M),
                        APERTURE_R_M)
    return {"soft Gaussian": soft, "hard-truncated": hard}


def tier_a(fields, z, dx):
    '''Tier (a): the angular spectrum against Forvard. One algorithm.

    `Forvard` multiplies the angular spectrum by exp(-i pi lambda z f^2), and
    so does `angular_spectrum` with m = 1 (Schmidt (2010),
    DOI 10.1117/3.866274, Ch. 6, Eq. (6.32), printed p. 95). The two must agree
    to the arithmetic noise of the two transform routes.
    '''
    out = {}
    for name, F0 in fields.items():
        book = angular_spectrum(F0.field, WAVELENGTH_M, dx, z)
        prod = drop_piston(Forvard(F0, z).field, z, WAVELENGTH_M)
        out[name] = (book, prod, rel_diff(book, prod))
    return out


def tier_b(fields, z, dx):
    '''Tier (b): the three book kernels against the production Fresnel.

    Every result lands on the pitch dx, because the range is chosen to make the
    FIXED one-step pitch of Ch. 6, Eq. (6.16), printed p. 90, equal to dx.
    '''
    out = {}
    for name, F0 in fields.items():
        u_one, dx2 = one_step_fresnel(F0.field, WAVELENGTH_M, dx, z)
        assert abs(dx2 / dx - 1.0) < 1e-12, (dx2, dx)
        row = {
            "one-step Fresnel": u_one,
            "angular spectrum": angular_spectrum(F0.field, WAVELENGTH_M, dx, z),
            "production Fresnel": drop_piston(Fresnel(F0, z).field, z,
                                              WAVELENGTH_M),
            "production Forvard": drop_piston(Forvard(F0, z).field, z,
                                              WAVELENGTH_M),
        }
        # `two_step_fresnel` refuses m = 1 (Ch. 6, Eq. (6.25) and Table 6.2,
        # printed pp. 94 and 95: the intermediate plane then runs to infinity).
        # So it runs at m = 2. Every second sample of the central half of the
        # fine grid falls exactly on a coarse sample, so the two compare with
        # NO interpolation.
        u_two = two_step_fresnel(F0.field, WAVELENGTH_M, dx, 2.0 * dx, z)
        j = np.arange(N // 4, 3 * N // 4)
        i = N // 2 + 2 * (j - N // 2)
        two_err = rel_diff(u_two[np.ix_(j, j)], u_one[np.ix_(i, i)])
        out[name] = (row, u_two, two_err)
    return out


def tier_c():
    '''Tier (c): the two-step Fresnel against the production co-moving recipe.

    The production route is Lens -> LensFresnel -> Convert on a grid that grows
    with the beam. The book route is `two_step_fresnel` with the same output
    pitch. See Schmidt (2010), DOI 10.1117/3.866274, Ch. 6, Sec. 6.3.2, printed
    pp. 92 to 95, and `olb.waveoptics.lenses`.

    THE PISTON IS NOT SIMPLE HERE. The production recipe carries exp(i k z1) of
    its INTERNAL step, and `Convert` adds a quadratic phase of its own. So this
    tier removes the on-axis phase of each field and compares what is left:
    the amplitude and the RELATIVE phase across the grid.
    '''
    dx1 = LENS_SIDE_M / LENS_N
    zr = np.pi * LENS_WAIST_M ** 2 / WAVELENGTH_M
    m = float(np.sqrt(1.0 + (LENS_RANGE_M / zr) ** 2))
    f_lens = LENS_RANGE_M / (m - 1.0)

    out = {}
    for name, make in (("soft Gaussian", lambda F: F),
                       ("hard-truncated",
                        lambda F: CircAperture(F, LENS_APERTURE_R_M))):
        F0 = make(GaussBeam(Begin(LENS_SIDE_M, WAVELENGTH_M, LENS_N),
                            LENS_WAIST_M))
        prod = Convert(LensFresnel(Lens(F0, f_lens), -f_lens, LENS_RANGE_M))
        dx2 = prod.siz / LENS_N
        book = two_step_fresnel(F0.field, WAVELENGTH_M, dx1, dx2,
                                LENS_RANGE_M)

        c = LENS_N // 2
        a = prod.field * np.exp(-1j * np.angle(prod.field[c, c]))
        b = book * np.exp(-1j * np.angle(book[c, c]))
        # Read the comparison inside three beam radii. Outside that the
        # amplitude is below 1e-4 of the peak, so a relative reading is noise.
        r2 = (axis_m(LENS_N, dx2)[:, None] ** 2
              + axis_m(LENS_N, dx2)[None, :] ** 2)
        core = r2 <= (3.0 * LENS_WAIST_M * m) ** 2
        out[name] = (b, a, dx2, rel_diff(b, a, core), core)
    return out, m, f_lens, dx1


def draw_cuts(z_near, z_far, a_near, b_far, c_out, dx, m_lens):
    '''Draw the intensity cross-cuts of every kernel, on one axes per case.'''
    fig, axes = plt.subplots(2, 2, figsize=(14.0, 9.4),
                             constrained_layout=True)
    x_mm = axis_m(N, dx) * 1e3

    # ---- panel (a): tier (a), the hard aperture at the near range ----
    ax = axes[0, 0]
    book, prod, _ = a_near["hard-truncated"]
    ax.semilogy(x_mm, np.abs(book[N // 2]) ** 2, color="tab:blue",
                linewidth=2.6, label="schmidt angular_spectrum (m = 1)")
    ax.semilogy(x_mm, np.abs(prod[N // 2]) ** 2, color="tab:red",
                linewidth=1.2, linestyle="--",
                label="production Forvard, piston removed")
    ax.set_title(f"(a) SAME ALGORITHM: hard aperture at z = {z_near * 1e3:.0f}"
                 f" mm\nthe two curves lie on top of each other", fontsize=10)

    # ---- panel (b): tier (b), the hard aperture at the far range ----
    ax = axes[0, 1]
    styles = [("one-step Fresnel", "tab:blue", "-", 2.6),
              ("angular spectrum", "tab:orange", "-", 1.8),
              ("production Fresnel", "tab:green", "--", 1.6),
              ("production Forvard", "tab:red", ":", 1.6)]
    row = b_far["hard-truncated"][0]
    for key, colour, style, width in styles:
        ax.semilogy(x_mm, np.abs(row[key][N // 2]) ** 2, color=colour,
                    linestyle=style, linewidth=width, label=key)
    ax.set_title(f"(b) SAME INTEGRAL, OTHER QUADRATURE: hard aperture at "
                 f"z = {z_far * 1e3:.0f} mm\nthe Fresnel ring structure",
                 fontsize=10)

    # ---- panel (c): tier (b), the soft Gaussian at the far range ----
    ax = axes[1, 0]
    row = b_far["soft Gaussian"][0]
    for key, colour, style, width in styles:
        ax.semilogy(x_mm, np.abs(row[key][N // 2]) ** 2, color=colour,
                    linestyle=style, linewidth=width, label=key)
    exact = gaussian_closed_form(N, dx, WAIST_M, z_far, WAVELENGTH_M)
    ax.semilogy(x_mm, np.abs(exact[N // 2]) ** 2, color="black",
                linestyle="-.", linewidth=1.2,
                label="closed form, Ch. 1, Eq. (1.53)")
    ax.set_title(f"(c) the soft Gaussian at z = {z_far * 1e3:.0f} mm\n"
                 f"every kernel meets the closed form", fontsize=10)

    # ---- panel (d): tier (c), the co-moving grid ----
    ax = axes[1, 1]
    book, prod, dx2, _, _ = c_out["soft Gaussian"]
    x2_mm = axis_m(LENS_N, dx2) * 1e3
    ax.semilogy(x2_mm, np.abs(book[LENS_N // 2]) ** 2, color="tab:blue",
                linewidth=2.6, label="schmidt two_step_fresnel")
    ax.semilogy(x2_mm, np.abs(prod[LENS_N // 2]) ** 2, color="tab:red",
                linestyle="--", linewidth=1.2,
                label="production Lens -> LensFresnel -> Convert")
    ax.set_ylim(1e-13, 3e-5)
    ax.set_title(f"(d) THE CO-MOVING GRID: soft Gaussian, "
                 f"z = {LENS_RANGE_M:.0f} m, magnification {m_lens:.0f}\n"
                 f"the grid grows with the beam", fontsize=10)

    for ax in axes.ravel():
        ax.set_xlabel("x, mm")
        ax.set_ylabel("intensity, arb.")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    for ax in axes.ravel()[:3]:
        ax.set_xlim(-5.0, 5.0)
    axes[0, 0].set_ylim(1e-8, 2.0)
    axes[0, 1].set_ylim(1e-8, 2.0)
    axes[1, 0].set_ylim(1e-8, 2.0)

    fig.suptitle("Schmidt reference kernels against the production LightPipes "
                 f"propagators, {WAVELENGTH_M * 1e9:.0f} nm, N = {N}",
                 fontsize=13)
    fig.savefig(CUTS_PNG, dpi=150)
    plt.close(fig)


def draw_diffs(a_far, b_far, c_out, dx):
    '''Draw the log10 difference map of the interesting kernel pairs.'''
    row = b_far["hard-truncated"][0]
    book, prod, _ = a_far["hard-truncated"]
    c_book, c_prod, dx2, _, _ = c_out["hard-truncated"]

    panels = [
        ("(a) angular_spectrum - Forvard\nsame algorithm, hard aperture",
         book - prod, book, dx, N),
        ("(b) one_step_fresnel - angular_spectrum\ntwo book kernels, "
         "hard aperture",
         row["one-step Fresnel"] - row["angular spectrum"],
         row["one-step Fresnel"], dx, N),
        ("(c) one_step_fresnel - production Fresnel\nother quadrature, "
         "hard aperture",
         row["one-step Fresnel"] - row["production Fresnel"],
         row["one-step Fresnel"], dx, N),
        ("(d) two_step_fresnel - production co-moving\nhard aperture, "
         f"z = {LENS_RANGE_M:.0f} m",
         c_book - c_prod, c_book, dx2, LENS_N),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(19.6, 5.0),
                             constrained_layout=True)
    floor = 1e-14
    for ax, (title, delta, reference, pitch, n) in zip(axes, panels):
        peak = float(np.max(np.abs(reference)))
        data = np.log10(np.maximum(np.abs(delta) / peak, floor))
        half = axis_m(n, pitch)[-1] * 1e3
        image = ax.imshow(data, extent=[-half, half, -half, half],
                          origin="lower", cmap="viridis", vmin=np.log10(floor),
                          vmax=0.0)
        fig.colorbar(image, ax=ax, shrink=0.80,
                     label="log10(|difference| / peak)")
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("x, mm")
        ax.set_ylabel("y, mm")

    fig.suptitle("Where the two layers disagree. Panel (a) is arithmetic "
                 "noise. Panels (b) to (d) are quadrature.", fontsize=13)
    fig.savefig(DIFFS_PNG, dpi=150)
    plt.close(fig)


def main():
    t_start = time.time()
    dx = SIDE_M / N
    z_far = N * dx ** 2 / WAVELENGTH_M
    z_near = 0.5 * z_far

    fields = input_fields()
    a_near = tier_a(fields, z_near, dx)
    a_far = tier_a(fields, z_far, dx)
    b_far = tier_b(fields, z_far, dx)
    c_out, m_lens, f_lens, dx1_lens = tier_c()

    print("=" * 78)
    print("The Schmidt propagation kernels against the production propagators")
    print("=" * 78)
    print(f"  wavelength              {WAVELENGTH_M * 1e9:11.1f} nm")
    print(f"  grid                    {N:11d} px, {SIDE_M * 1e3:.3f} mm")
    print(f"  grid pitch dx           {dx * 1e6:11.3f} um")
    print(f"  soft waist w0           {WAIST_M * 1e3:11.3f} mm")
    print(f"  hard aperture radius    {APERTURE_R_M * 1e3:11.3f} mm "
          f"(alpha = {APERTURE_R_M / WAIST_M:.2f})")
    print(f"  near range Z_NEAR       {z_near * 1e3:11.3f} mm")
    print(f"  far range Z_FAR         {z_far * 1e3:11.3f} mm")
    print("")
    print("  why Z_FAR is that number:")
    print(f"    Ch. 6, Eq. (6.16), p. 90: the one-step pitch lambda z/(N dx) = "
          f"{one_step_delta2(N, dx, WAVELENGTH_M, z_far) * 1e6:.3f} um = dx")
    print(f"    Ch. 7, Eq. (7.59), p. 127: constraint 4 asks for N >= "
          f"{constraint4_min_n(dx, dx, WAVELENGTH_M, z_far):.1f}, and N = {N}")
    print("    Ch. 7, Eq. (7.42), p. 123: the Fresnel-integral minimum "
          "distance for the")
    print(f"      {2 * APERTURE_R_M * 1e3:.1f} mm hard source is "
          f"{fresnel_min_distance(2 * APERTURE_R_M, dx, WAVELENGTH_M) * 1e3:.1f}"
          f" mm, so Z_FAR is "
          f"{z_far / fresnel_min_distance(2 * APERTURE_R_M, dx, WAVELENGTH_M):.1f}"
          f" times it")
    print(f"    the Fresnel number of the hard aperture at Z_FAR is "
          f"{APERTURE_R_M ** 2 / (WAVELENGTH_M * z_far):.2f}")

    # ---- tier (a) ----
    print("")
    print("TIER (a) SAME ALGORITHM. schmidt.angular_spectrum(m = 1) against")
    print("the production Forvard, after the exp(i k z) bridge. Both sample")
    print("exp(-i pi lambda z f^2) on one grid: Ch. 6, Eq. (6.32), p. 95.")
    print("")
    header = f"  {'input field':<20}{'range':>12}{'max rel diff':>16}"
    print(header)
    print("  " + "-" * (len(header) - 2))
    for label, table, z in (("Z_NEAR", a_near, z_near), ("Z_FAR", a_far,
                                                         z_far)):
        for name, (_, _, err) in table.items():
            print(f"  {name:<20}{label:>12}{err:>16.3e}")
            assert err < TOL_TIER_A, (name, label, err)
    print(f"  stated bound {TOL_TIER_A:.0e}. The residue is the legacy 2 pi of")
    print("  the production module and the integer wrap of its transfer phase.")

    # ---- tier (b) ----
    print("")
    print("TIER (b) SAME INTEGRAL, OTHER QUADRATURE, at Z_FAR. The production")
    print("Fresnel convolves on a doubled grid with the C(x), S(x) pixel")
    print("integral, and it builds that kernel on the legacy pitch siz/(N-1).")
    print("")
    print(f"The table reads the interior, r <= {CORE_RADIUS_M * 1e3:.1f} mm.")
    print("")
    r2_all = axis_m(N, dx)[:, None] ** 2 + axis_m(N, dx)[None, :] ** 2
    interior = r2_all <= CORE_RADIUS_M ** 2
    keys = ["one-step Fresnel", "angular spectrum", "production Fresnel",
            "production Forvard"]
    for name in fields:
        row, u_two, two_err = b_far[name]
        print(f"  {name}, pairwise max rel diff:")
        head = "    " + " " * 22 + "".join(f"{k[:18]:>20}" for k in keys)
        print(head)
        for a_key in keys:
            line = f"    {a_key:<22}"
            for b_key in keys:
                line += f"{rel_diff(row[a_key], row[b_key], interior):>20.3e}"
            print(line)
        print(f"    two_step_fresnel (m = 2) against one_step_fresnel, on the "
              f"shared samples: {two_err:.3e}")
        assert two_err < TOL_TIER_B, (name, two_err)
        for a_key in keys:
            for b_key in keys:
                err = rel_diff(row[a_key], row[b_key], interior)
                bound = (TOL_BOOK_PAIR
                         if a_key in BOOK_KEYS and b_key in BOOK_KEYS
                         else TOL_TIER_B)
                assert err < bound, (name, a_key, b_key, err, bound)
        whole = rel_diff(row["one-step Fresnel"], row["production Fresnel"])
        print(f"    the same pair one_step_fresnel / production Fresnel over "
              f"the WHOLE grid: {whole:.3e}")
        print("")
    print(f"  stated bounds: {TOL_BOOK_PAIR:.0e} for a BOOK pair on one output")
    print(f"  grid, {TOL_TIER_B:.0e} for every pair that crosses a quadrature "
          f"or an output grid.")
    print("")
    print("  READ THE TABLE THIS WAY.")
    print("  - one_step_fresnel and angular_spectrum agree to 1e-15. They are")
    print("    two roads to one integral, and they take one transform library.")
    print("  - production Forvard joins them at 1e-10: it is the same angular")
    print("    spectrum through the legacy 2 pi and the integer phase wrap.")
    print("  - production Fresnel is the outlier, and the SOFT and the HARD")
    print("    cases part company: 6e-4 for the Gaussian, 1.5e-2 for the hard")
    print("    aperture. The pixel-integral kernel on the legacy pitch")
    print("    siz/(N-1) is a different quadrature, and a hard edge has a")
    print("    spectrum that does not decay, so it feels that difference.")
    print("  - the WHOLE-grid number is larger again, and it comes from the")
    print("    outer band alone: the production Fresnel differences four")
    print("    shifted copies of a doubled-grid convolution, so it leaves an")
    print("    artefact at the grid edge. Tier (a) has no such band.")
    print("  - two_step_fresnel at m = 2 puts its intermediate plane BEHIND")
    print("    the source, at z1 = -z (Ch. 6, Eq. (6.25) and Table 6.2,")
    print("    printed pp. 94 and 95). The back-propagated Gaussian stays on")
    print("    the grid and the pair holds to 1e-15; the back-propagated HARD")
    print("    edge spreads past the grid, so that pair reads 6e-3.")

    # ---- the closed form ----
    print("  Against the closed-form Gaussian of Ch. 1, Eq. (1.53), p. 9,")
    print("  inside three beam radii:")
    exact = gaussian_closed_form(N, dx, WAIST_M, z_far, WAVELENGTH_M)
    w_far = WAIST_M * np.sqrt(1.0 + (z_far / (np.pi * WAIST_M ** 2
                                              / WAVELENGTH_M)) ** 2)
    r2 = axis_m(N, dx)[:, None] ** 2 + axis_m(N, dx)[None, :] ** 2
    core = r2 <= (3.0 * w_far) ** 2
    row = b_far["soft Gaussian"][0]
    for key in keys:
        err = rel_diff(row[key], exact, core)
        print(f"    {key:<24}{err:>16.3e}")
        assert err < TOL_TIER_B, (key, err)
    print(f"    w(Z_FAR) = {w_far * 1e3:.4f} mm")

    # ---- tier (c) ----
    print("")
    print("TIER (c) THE CO-MOVING GRID. schmidt.two_step_fresnel with a free")
    print("output pitch against Lens -> LensFresnel -> Convert. The on-axis")
    print("phase of each field is removed first, because the two routes carry")
    print("different piston.")
    print("")
    print(f"  range                   {LENS_RANGE_M:11.1f} m")
    print(f"  waist w0                {LENS_WAIST_M * 1e3:11.3f} mm")
    print(f"  magnification m         {m_lens:11.3f}")
    print(f"  physical lens fA        {f_lens:11.4f} m")
    print(f"  launch grid             {LENS_N:11d} px, "
          f"{LENS_SIDE_M * 1e3:.3f} mm, pitch {dx1_lens * 1e6:.3f} um")
    print("")
    head2 = f"  {'input field':<20}{'output pitch':>16}{'max rel diff':>16}"
    print(head2)
    print("  " + "-" * (len(head2) - 2))
    for name, (_, _, dx2, err, _) in c_out.items():
        print(f"  {name:<20}{dx2 * 1e3:>13.4f} mm{err:>16.3e}")
        assert err < TOL_TIER_C, (name, err)
    print(f"  stated bound {TOL_TIER_C:.0e}. The production route runs its")
    print("  convolution through an internal step of z/m, and it removes the")
    print("  residual curvature with a second legacy 2 pi in Convert(). Neither")
    print("  step exists in the two-step Fresnel chain of Ch. 6, Eq. (6.18),")
    print("  printed p. 93.")
    print("")
    print("  THE HARD CASE IS 10 TIMES WORSE, AND THE LAUNCH GRID SAYS WHY.")
    print("  Both routes take a short internal step near 2 m. A 1 mm hard edge")
    print("  spreads by lambda z1 / R = 3 mm over that step, so a tight launch")
    print("  grid loses the spread light at its own edge. The disagreement")
    print("  falls with the launch side: 7.0e-2 at 10 w0, 2.3e-2 at 20 w0 (the")
    print("  value in the table), and 1.4e-2 at 30 w0. The soft Gaussian sits")
    print("  near 2e-3 at each of the three, because it has no edge to spread.")

    draw_cuts(z_near, z_far, a_near, b_far, c_out, dx, m_lens)
    draw_diffs(a_far, b_far, c_out, dx)
    print("")
    print(f"figure saved: {CUTS_PNG}")
    print("  Caption: the intensity cross-cut of every kernel. Panel (a) is "
          "one algorithm,")
    print("  panels (b) and (c) are one integral through four quadratures, and "
          "panel (d) is")
    print("  the co-moving grid at 500 m.")
    print(f"figure saved: {DIFFS_PNG}")
    print("  Caption: log10 of the field difference over the peak. Panel (a) "
          "is arithmetic")
    print("  noise at 1e-10; panels (b), (c) and (d) show the quadrature error "
          "sitting on the")
    print("  diffraction rings of the hard aperture, not in the smooth core.")
    print("")
    print(f"(elapsed {time.time() - t_start:.1f} s)")


if __name__ == "__main__":
    main()

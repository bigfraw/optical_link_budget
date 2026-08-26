'''
The grid artefact of the FFT propagators, made visible on purpose.

A fidelity-2 field propagation runs on a square grid of a limited side. The FFT
propagators treat that grid as PERIODIC. So the grid behaves like a square
waveguide: the light that leaves one edge comes back at the opposite edge. A
grid that is too small therefore folds the beam tail back onto the beam, and
the answer is wrong.

The script runs one well-behaved Gaussian link three times:

  (a) a proper grid, 8 times the final beam radius, with Forvard. The result is
      clean.
  (b) the SAME link on a grid of only 2 times the final beam radius, with
      Forvard. The wrap-around artefact appears.
  (c) the same small grid with Fresnel. That propagator convolves on a doubled
      grid, so its edge behaviour is different, but the small grid still cuts
      the beam tail.

The fourth panel gives the analytic ABCD result (GForvard) on the small grid.
That route has no FFT and no grid artefact, so it is the reference.

The figure goes to `examples/waveoptics/grid_artefacts.png`.

Sources:
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB, DOI 10.1117/3.866274, Ch. 6 and Ch. 7. The periodic artefact of the
  spectral propagator and the grid-extent rule.

Run from the repo root:
    python -m examples.waveoptics.grid_artefacts
'''

import matplotlib.pyplot as plt
import numpy as np

from olb.waveoptics import (Begin, CircAperture, Forvard, Fresnel, GaussBeam,
                            GForvard, Intensity, Power)

WAVELENGTH_M = 1550e-9
WAIST_M = 5e-3              # the waist radius at the launch plane
RANGE_M = 200.0
N = 512

PNG = "examples/waveoptics/grid_artefacts.png"


def analytic_radius(w0, z, lam):
    '''Give the free-space Gaussian radius at the range z.

    w(z) = w0 * sqrt(1 + (z/zR)^2), with zR = pi*w0^2/lambda. See Siegman,
    Lasers, ISBN 978-0935702118.
    '''
    zR = np.pi * w0 ** 2 / lam
    return w0 * np.sqrt(1 + (z / zR) ** 2), zR


def bucket_power(field, radius_m):
    '''Give the power fraction inside a circle of the given radius.'''
    return float(Power(CircAperture(field, radius_m)) / Power(field))


def run(grid_size_m, propagator):
    '''Launch the beam on one grid and propagate it with one propagator.'''
    F0 = GaussBeam(Begin(grid_size_m, WAVELENGTH_M, N), WAIST_M)
    return propagator(F0, RANGE_M)


def draw(panels, w_z):
    '''Draw the four log-intensity panels on one shared colour scale.'''
    peak = max(Intensity(field).max() for _, field in panels)
    floor = 1e-4
    fig, axes = plt.subplots(1, len(panels), figsize=(4.6 * len(panels), 4.6),
                             constrained_layout=True)
    for ax, (title, field) in zip(np.atleast_1d(axes).ravel(), panels):
        half_mm = field.siz / 2 * 1e3
        extent = [-half_mm, half_mm, -half_mm, half_mm]
        data = np.log10(np.maximum(Intensity(field) / peak, floor))
        image = ax.imshow(data, extent=extent, origin="lower", cmap="magma",
                          vmin=np.log10(floor), vmax=0.0)
        ax.add_patch(plt.Circle((0.0, 0.0), w_z * 1e3, fill=False,
                                linestyle="--", linewidth=1.2,
                                edgecolor="white"))
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("x, mm")
        ax.set_ylabel("y, mm")
        fig.colorbar(image, ax=ax, shrink=0.80, label="log10(I / I_peak)")

    fig.suptitle("The FFT grid artefact: one Gaussian link, "
                 f"waist {WAIST_M * 1e3:.0f} mm, range {RANGE_M:.0f} m, "
                 f"{WAVELENGTH_M * 1e9:.0f} nm. The dashed circle is w(z).",
                 fontsize=13)
    fig.savefig(PNG, dpi=150)
    return fig


def main():
    w_z, zR = analytic_radius(WAIST_M, RANGE_M, WAVELENGTH_M)
    big = 8 * w_z               # the grid-extent rule of Schmidt, Ch. 6
    # 2 w(z) is the documented failure case of the propagators self-check.
    small = 2.0 * w_z           # deliberately too small

    clean = run(big, Forvard)
    wrapped = run(small, Forvard)
    convolved = run(small, Fresnel)
    exact_big = run(big, GForvard)
    exact_small = run(small, GForvard)

    print("=" * 72)
    print("The periodic grid artefact of the FFT propagators")
    print("=" * 72)
    print(f"  wavelength              {WAVELENGTH_M * 1e9:9.1f} nm")
    print(f"  waist radius w0         {WAIST_M * 1e3:9.3f} mm")
    print(f"  range z                 {RANGE_M:9.1f} m")
    print(f"  Rayleigh range zR       {zR:9.3f} m")
    print(f"  analytic w(z)           {w_z * 1e3:9.4f} mm")
    print(f"  proper grid side        {big * 1e3:9.2f} mm "
          f"({big / w_z:.2f} w(z))")
    print(f"  small grid side         {small * 1e3:9.2f} mm "
          f"({small / w_z:.2f} w(z))")
    print(f"  pixels per side         {N:9d}")
    print("")
    print("Bucket power inside the radius w(z). The reference of each row is")
    print("the analytic ABCD route (GForvard) on the SAME grid.")
    print("")
    header = (f"{'run':<38}{'bucket':>10}{'reference':>11}"
              f"{'rel. error':>12}")
    print(header)
    print("-" * len(header))
    rows = [("proper grid, Forvard", clean, exact_big),
            ("small grid, Forvard", wrapped, exact_small),
            ("small grid, Fresnel", convolved, exact_small),
            ("small grid, GForvard (analytic)", exact_small, exact_small)]
    for name, field, reference in rows:
        b = bucket_power(field, w_z)
        r = bucket_power(reference, w_z)
        print(f"{name:<38}{b:>10.6f}{r:>11.6f}{abs(b - r) / r:>11.2%}")
    print("")
    print(f"free-space Gaussian value, 1 - exp(-2) = {1 - np.exp(-2):.6f}")
    print("  The small grid cuts the beam tail, so each small-grid bucket sits")
    print("  above that value. The Forvard row adds the wrap-around error on")
    print("  top of the cut. See Schmidt, DOI 10.1117/3.866274, Ch. 6.")

    panels = [
        (f"(a) proper grid ({big / w_z:.1f} w(z)), Forvard:\n"
         "the beam stays away from the edge", clean),
        (f"(b) grid too small ({small / w_z:.2f} w(z)), Forvard:\n"
         "the FFT wraps the beam at the edge", wrapped),
        (f"(c) grid too small ({small / w_z:.2f} w(z)), Fresnel:\n"
         "no wrap, but the grid cuts the tail", convolved),
        (f"(d) grid too small ({small / w_z:.2f} w(z)), GForvard:\n"
         "the analytic route, no grid artefact", exact_small),
    ]
    draw(panels, w_z)
    print("")
    print(f"figure saved: {PNG}")
    plt.show()


if __name__ == "__main__":
    main()

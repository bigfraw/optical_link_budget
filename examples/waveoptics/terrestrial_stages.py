'''
The fidelity-2 field propagation, stage by stage, on a near-field link.

The script runs the no-turbulence field validator of `olb.waveoptics` on a
horizontal link that BREAKS the analytic fidelity-0 total. The transmit waist
is 120 mm and the launch aperture is 150 mm, so the truncation ratio alpha is
0.625. The range is 1 km, but the Rayleigh range of that waist is about 29 km.
The receiver therefore sits deep inside the near field, where the far-field
analytic transmit efficiency does not hold.

The script does three things:

  1. It propagates the scenario with `propagate_scenario` and it draws each
     stage of the propagation.
  2. It prints the fidelity-0 analytic numbers next to the fidelity-2 numbers,
     with the difference.
  3. It adds a retroreflected return leg with the raw primitives: a clip at a
     corner-cube aperture, a Fresnel propagation back over the same range, and
     a clip at the original launch aperture.

The figure goes to `examples/waveoptics/figures/terrestrial_stages.png`.

The layer builds NO Term and it changes NO budget. See the README.

Run from the repo root:
    python -m examples.waveoptics.terrestrial_stages
'''

import matplotlib.pyplot as plt
import numpy as np

from olb import HorizontalPath, Terminal, Transmitter, SMF
from olb.models.gaussian_efficiency import tx_efficiency_loss_db
from olb.models.geometric import geometric_loss_db
from olb.scenario import TerrestrialChannel, TerrestrialScenario
from olb.waveoptics import (CircAperture, Fresnel, GridSpec, Intensity, Power,
                            propagate_scenario)

WAVELENGTH_M = 1550e-9
WAIST_M = 0.12              # the transmit waist at the launch plane
TX_APERTURE_M = 0.15        # alpha = (D/2)/w0 = 0.625, a hard truncation
RX_APERTURE_M = 0.20
RANGE_M = 1e3
CORNER_CUBE_M = 63.5e-3     # a 2.5 inch hollow corner cube

# The manual grid of the run.py self-check. The automatic grid gives the same
# extent rule, but this value keeps the script fast and it keeps the numbers
# equal to the self-check.
GRID = GridSpec(size_m=1.0, n=1024)

PNG = "examples/waveoptics/figures /terrestrial_stages.png"


def build_scenario():
    '''Build the near-field terrestrial scenario with an SMF receiver.'''
    near = Terminal(aperture_m=TX_APERTURE_M, wavelength_m=WAVELENGTH_M,
                    transmitter=Transmitter(waist_m=WAIST_M))
    far = Terminal(aperture_m=RX_APERTURE_M, wavelength_m=WAVELENGTH_M,
                   detector=SMF())
    channel = TerrestrialChannel(path_length_m=RANGE_M)
    return TerrestrialScenario(near=near, far=far, channel=channel)


def loss_db(power, reference):
    '''Give a positive-dB loss of one power against a reference power.'''
    return float(-10 * np.log10(power / reference)) + 0.0


def retro_leg(field_at_rx):
    '''Retroreflect the received field and bring it back to the launch plane.

    Three primitives do the whole return leg: a clip at the corner-cube
    aperture, a Fresnel propagation back over the range, and a clip at the
    original launch aperture.
    '''
    on_cube = CircAperture(field_at_rx, CORNER_CUBE_M / 2)
    back = Fresnel(on_cube, RANGE_M)
    collected = CircAperture(back, TX_APERTURE_M / 2)
    return on_cube, back, collected


def draw(panels, power_launch):
    '''Draw one intensity panel for each stage of the two legs.

    Each panel shows the full grid, so the grid edge stays visible. The title
    holds the stage name and the power of that stage in positive dB against the
    launch power. A dashed white circle marks the aperture of the stage.
    '''
    n_col = 4
    n_row = int(np.ceil(len(panels) / n_col))
    fig, axes = plt.subplots(n_row, n_col, figsize=(4.2 * n_col, 4.0 * n_row),
                             constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()

    for ax, (label, field, aperture_m, log_scale) in zip(axes, panels):
        half_mm = field.siz / 2 * 1e3
        extent = [-half_mm, half_mm, -half_mm, half_mm]
        I = Intensity(field)
        if log_scale:
            floor = 1e-6
            data = np.log10(np.maximum(I / I.max(), floor))
            cmap, vmin, vmax = "magma", np.log10(floor), 0.0
            bar_label = "log10(I / I_max)"
        else:
            data = I
            cmap, vmin, vmax = "inferno", 0.0, None
            bar_label = "intensity, W/m^2"
        image = ax.imshow(data, extent=extent, origin="lower", cmap=cmap,
                          vmin=vmin, vmax=vmax)
        if aperture_m is not None:
            ax.add_patch(plt.Circle((0.0, 0.0), aperture_m / 2 * 1e3,
                                    fill=False, linestyle="--", linewidth=1.2,
                                    edgecolor="white"))
        ax.set_title(f"{label}\n{loss_db(Power(field), power_launch):.2f} dB "
                     "below launch", fontsize=10)
        ax.set_xlabel("x, mm")
        ax.set_ylabel("y, mm")
        fig.colorbar(image, ax=ax, shrink=0.82, label=bar_label)

    for ax in axes[len(panels):]:
        ax.axis("off")

    fig.suptitle("Fidelity-2 field propagation, near-field terrestrial link: "
                 f"waist {WAIST_M * 1e3:.0f} mm, launch aperture "
                 f"{TX_APERTURE_M * 1e3:.0f} mm, range "
                 f"{RANGE_M / 1e3:.0f} km, {WAVELENGTH_M * 1e9:.0f} nm",
                 fontsize=13)
    fig.savefig(PNG, dpi=150)
    return fig


def main():
    scenario = build_scenario()
    result = propagate_scenario(scenario, HorizontalPath(RANGE_M), grid=GRID)
    stages = dict(result.stages)
    power_launch = Power(stages["launch"])

    # ---- 1. the two fidelities, side by side ----
    analytic_tx = float(tx_efficiency_loss_db(TX_APERTURE_M, WAIST_M))
    analytic_geo = float(geometric_loss_db(RANGE_M, WAIST_M, RX_APERTURE_M,
                                           wavelength=WAVELENGTH_M))
    rows = [("launch truncation", result.tx_truncation_db, analytic_tx),
            ("geometric spread", result.geometric_loss_db, analytic_geo),
            ("truncation + spread",
             result.tx_truncation_db + result.geometric_loss_db,
             analytic_tx + analytic_geo)]

    print("=" * 68)
    print("Near-field terrestrial link, fidelity 2 against fidelity 0")
    print("=" * 68)
    print(f"  wavelength            {WAVELENGTH_M * 1e9:9.1f} nm")
    print(f"  transmit waist        {WAIST_M * 1e3:9.1f} mm")
    print(f"  launch aperture       {TX_APERTURE_M * 1e3:9.1f} mm "
          f"(alpha = {(TX_APERTURE_M / 2) / WAIST_M:.3f})")
    print(f"  receive aperture      {RX_APERTURE_M * 1e3:9.1f} mm")
    print(f"  range                 {RANGE_M:9.1f} m")
    print(f"  Rayleigh range        {np.pi * WAIST_M ** 2 / WAVELENGTH_M:9.1f} m")
    print(f"  grid side             {GRID.size_m:9.3f} m")
    print(f"  pixels per side       {GRID.n:9d}")
    print(f"  propagator            {result.propagator:>9}")
    print("")
    print(f"{'quantity':<22}{'fidelity 2':>12}{'fidelity 0':>12}"
          f"{'difference':>12}")
    print("-" * 58)
    for name, fid2, fid0 in rows:
        print(f"{name:<22}{fid2:>11.3f}dB{fid0:>11.3f}dB{fid2 - fid0:>11.3f}dB")
    print("")
    print(f"  fibre coupling loss   {result.smf_coupling_db:9.3f} dB "
          "(single-mode fibre, fidelity 2)")
    print("  all the numbers are losses in positive dB.")
    print("  The analytic transmit efficiency is a FAR-FIELD form. The receiver")
    print("  sits inside the Rayleigh range here, so the two totals disagree.")

    # ---- 2. the retro return leg, three primitives ----
    on_cube, back, collected = retro_leg(stages["at rx plane"])
    power_cube = Power(on_cube)
    power_back = Power(collected)

    print("")
    print("=" * 68)
    print(f"Retroreflected return leg, {CORNER_CUBE_M * 1e3:.1f} mm corner cube")
    print("=" * 68)
    print(f"  clip at the corner cube      {loss_db(power_cube, Power(stages['at rx plane'])):9.3f} dB")
    print(f"  return leg, cube to aperture {loss_db(power_back, power_cube):9.3f} dB")
    print(f"  round trip, launch to launch {loss_db(power_back, power_launch):9.3f} dB")
    print("  The return leg starts at the corner cube, so it is a NEW source of")
    print("  the cube diameter. It spreads much more than the outbound beam.")

    # ---- 3. the figure ----
    panels = [
        ("1. launch", stages["launch"], TX_APERTURE_M, False),
        ("2. after the launch clip", stages["after tx clip"], TX_APERTURE_M,
         False),
        ("3. at the receive plane", stages["at rx plane"], RX_APERTURE_M,
         False),
        ("4. after the receive clip", stages["after rx clip"], RX_APERTURE_M,
         False),
        ("4b. receive clip, log scale", stages["after rx clip"], RX_APERTURE_M,
         True),
        ("5. retro: on the corner cube", on_cube, CORNER_CUBE_M, False),
        ("6. retro: back at the launch plane", back, TX_APERTURE_M, True),
        ("7. retro: collected by the launch aperture", collected,
         TX_APERTURE_M, True),
    ]
    draw(panels, power_launch)
    print("")
    print(f"figure saved: {PNG}")
    plt.show()


if __name__ == "__main__":
    main()

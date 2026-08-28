'''
The fidelity-2 field propagation on a 600 km space link, against fidelity 0.

The script launches a truncated Gaussian beam through an obscured circular
aperture, and it propagates the field to a satellite at 600 km. The launch
waist is 50 mm and the launch aperture is 100 mm, so the truncation ratio
alpha is 1.0. A 0.3 central obscuration blocks the middle of the aperture,
which is the Cassegrain secondary of a real beam director.

THE CO-MOVING GRID. A flat grid cannot do this link. The beam radius grows
from 50 mm to 5.9 m, a factor of 118. A flat grid must hold the 5.9 m beam AND
resolve the 15 mm obscuration, so it needs more than 40000 pixels per side.
`GridSpec.for_scenario` therefore selects the SCALED route: the grid starts at
the launch plane, and it grows with the beam. `propagate_scenario` reads the
`scaled` attribute, and it runs the three-call lens recipe of
`olb.waveoptics.lenses`. See Schmidt, DOI 10.1117/3.866274, Ch. 7.

THE RESULT. The receiver is in the absolute far field, so the fidelity-2 total
and the fidelity-0 analytic total agree to about 0.01 dB. The SPLIT between
the two losses does NOT agree, because the two fidelities cut the loss in two
different places. The script prints both, and it says why.

The figure goes to `examples/waveoptics/figures/space_farfield.png`.

The layer builds NO Term and it changes NO budget. See the README.

Run from the repo root:
    python -m examples.waveoptics.space_farfield
'''

import matplotlib.pyplot as plt
import numpy as np

from olb import CircularOrbit, Terminal, Transmitter
from olb.models.gaussian_efficiency import tx_efficiency_loss_db
from olb.models.geometric import geometric_loss_db
from olb.scenario import Channel, SpaceScenario
from olb.waveoptics import (Intensity, Power, beam_magnification,
                            propagate_scenario)

WAVELENGTH_M = 1550e-9
WAIST_M = 0.05              # the transmit waist at the launch plane
TX_APERTURE_M = 0.10        # alpha = (D/2)/w0 = 1.0, a real truncation
TX_OBSCURATION = 0.3        # the linear central-obscuration ratio
RX_APERTURE_M = 0.50        # the satellite receive telescope
ALTITUDE_M = 600e3

PNG = "examples/waveoptics/figures/space_farfield.png"


def build_scenario():
    '''Build the uplink scenario: an obscured 100 mm launch to a 600 km orbit.'''
    ground = Terminal(aperture_m=TX_APERTURE_M,
                      obscuration_ratio=TX_OBSCURATION,
                      wavelength_m=WAVELENGTH_M,
                      transmitter=Transmitter(waist_m=WAIST_M))
    space = Terminal(aperture_m=RX_APERTURE_M, wavelength_m=WAVELENGTH_M)
    return SpaceScenario(ground=ground, space=space, direction="uplink",
                         channel=Channel(altitude_m=ALTITUDE_M))


def loss_db(power, reference):
    '''Give a positive-dB loss of one power against a reference power.'''
    return float(-10 * np.log10(power / reference)) + 0.0


def draw(stages, losses, beam_radius_m):
    '''Draw the launch plane, the obscured annulus, and the far-field pattern.

    The first two panels sit at the launch plane, so their axes are in mm. The
    last two panels sit at the receive plane on the co-moving grid, so their
    axes are in m. Both take a log scale, because the obscuration puts deep
    interference rings in the far field.

    The receive aperture is 0.5 m on a 47 m grid, so no map panel can show it.
    The fourth panel is therefore a radial cut. It puts the receive aperture,
    the beam radius, and the analytic Gaussian on one axis.
    '''
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 9.6),
                             constrained_layout=True)
    axes = axes.ravel()
    launch, clipped, at_rx = stages
    panels = [
        (axes[0], "1. launch, before the aperture", launch, 1e3, "mm",
         [(TX_APERTURE_M / 2, "aperture")], False),
        (axes[1], "2. after the aperture and the obscuration", clipped, 1e3,
         "mm", [(TX_APERTURE_M / 2, "aperture"),
                (TX_OBSCURATION * TX_APERTURE_M / 2, "obscuration")], False),
        (axes[2], "3. at the receive plane, 600 km", at_rx, 1.0, "m",
         [(beam_radius_m, "beam radius w(z)")], True),
    ]

    for ax, label, field, unit, axis_unit, circles, log in panels:
        half = field.siz / 2 * unit
        extent = [-half, half, -half, half]
        I = Intensity(field)
        if log:
            floor = 1e-5
            data = np.log10(np.maximum(I / I.max(), floor))
            cmap, vmin, vmax = "magma", np.log10(floor), 0.0
            bar_label = "log10(I / I_max)"
        else:
            data = I
            cmap, vmin, vmax = "inferno", 0.0, None
            bar_label = "intensity, W/m^2"
        image = ax.imshow(data, extent=extent, origin="lower", cmap=cmap,
                          vmin=vmin, vmax=vmax)
        for radius_m, name in circles:
            ax.add_patch(plt.Circle((0.0, 0.0), radius_m * unit, fill=False,
                                    linestyle="--", linewidth=1.4,
                                    edgecolor="white", label=name))
        ax.legend(loc="upper right", fontsize=8, framealpha=0.4,
                  labelcolor="white")
        ax.set_title(f"{label}\n{losses[label]}", fontsize=10)
        ax.set_xlabel(f"x, {axis_unit}")
        ax.set_ylabel(f"y, {axis_unit}")
        fig.colorbar(image, ax=ax, shrink=0.82, label=bar_label)

    _radial_cut(axes[3], at_rx, beam_radius_m, losses)

    fig.suptitle("Fidelity-2 field propagation on a co-moving grid, 600 km "
                 f"uplink\nwaist {WAIST_M * 1e3:.0f} mm, launch aperture "
                 f"{TX_APERTURE_M * 1e3:.0f} mm, obscuration "
                 f"{TX_OBSCURATION:g}, receive aperture "
                 f"{RX_APERTURE_M * 1e3:.0f} mm, "
                 f"{WAVELENGTH_M * 1e9:.0f} nm", fontsize=12)
    fig.savefig(PNG, dpi=150)
    return fig


def _radial_cut(ax, at_rx, beam_radius_m, losses):
    '''Draw the far-field radial profile against the analytic Gaussian.

    The propagated field carries the truncation rings and the obscuration
    rings. The analytic fidelity-0 model has neither: it is the UNtruncated
    Gaussian exp(-2 r^2/w(z)^2). The two differ on the axis, which is the
    reason that the loss SPLIT of the two fidelities differs.
    '''
    label = "4. radial cut at the receive plane"
    centre = at_rx.N // 2
    r_m = at_rx.xvalues[centre:]
    profile = Intensity(at_rx)[centre, centre:]
    profile = profile / profile[0]
    gauss = np.exp(-2 * r_m ** 2 / beam_radius_m ** 2)

    ax.semilogy(r_m, np.maximum(profile, 1e-6), color="tab:red",
                label="fidelity 2, the propagated field")
    ax.semilogy(r_m, gauss, color="tab:blue", linestyle="--",
                label="fidelity 0, the untruncated Gaussian")
    ax.axvline(RX_APERTURE_M / 2, color="black", linestyle=":",
               label=f"receive aperture radius, {RX_APERTURE_M / 2:.2f} m")
    ax.axvline(beam_radius_m, color="grey", linestyle=":",
               label=f"beam radius w(z), {beam_radius_m:.2f} m")
    ax.set_xlim(0.0, at_rx.siz / 2)
    ax.set_ylim(1e-6, 2.0)
    ax.set_xlabel("radius, m")
    ax.set_ylabel("I(r) / I(0)")
    ax.set_title(f"{label}\n{losses[label]}", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)


def main():
    scenario = build_scenario()
    orbit = CircularOrbit(altitude_m=ALTITUDE_M, elevation_deg=[90.0])
    range_m = float(np.max(orbit.slant_range_m))

    result = propagate_scenario(scenario, orbit)      # the automatic grid
    stages = dict(result.stages)
    power_launch = Power(stages["launch"])
    power_clipped = Power(stages["after tx clip"])
    power_at_rx = Power(stages["at rx plane"])
    power_collected = Power(stages["after rx clip"])

    magnification = beam_magnification(scenario, range_m)
    rayleigh_m = np.pi * WAIST_M ** 2 / WAVELENGTH_M
    beam_radius_m = WAIST_M * magnification

    # ---- the grid ----
    print("=" * 72)
    print("600 km uplink, fidelity 2 on the co-moving grid")
    print("=" * 72)
    print(f"  wavelength              {WAVELENGTH_M * 1e9:11.1f} nm")
    print(f"  transmit waist          {WAIST_M * 1e3:11.1f} mm")
    print(f"  launch aperture         {TX_APERTURE_M * 1e3:11.1f} mm "
          f"(alpha = {(TX_APERTURE_M / 2) / WAIST_M:.3f})")
    print(f"  launch obscuration      {TX_OBSCURATION:11.2f}")
    print(f"  receive aperture        {RX_APERTURE_M * 1e3:11.1f} mm")
    print(f"  range at zenith         {range_m * 1e-3:11.1f} km")
    print(f"  Rayleigh range          {rayleigh_m * 1e-3:11.3f} km "
          f"({range_m / rayleigh_m:.0f} Rayleigh ranges to the receiver)")
    print(f"  beam radius at range    {beam_radius_m:11.3f} m")
    print("")
    print(f"  grid route              "
          f"{'scaled (co-moving)' if result.grid.scaled else 'flat':>11}")
    print(f"  propagator              {result.propagator:>11}")
    print(f"  pixels per side         {result.grid.n:11d}")
    print(f"  grid side at launch     {result.grid.size_m:11.4f} m")
    print(f"  grid side at receiver   {result.grid.size_m * magnification:11.3f} m")
    print(f"  magnification m         {magnification:11.1f}")
    print(f"  pixel at launch         {result.grid.pixel_m * 1e3:11.4f} mm")
    print(f"  pixel at receiver       "
          f"{result.grid.pixel_m * magnification * 1e3:11.3f} mm")
    print(f"  power kept on the grid  {power_at_rx / power_clipped:11.4f}  "
          "(the rest is the far tail past the grid edge)")

    # ---- the two fidelities ----
    analytic_tx = float(tx_efficiency_loss_db(TX_APERTURE_M, WAIST_M,
                                              obscuration_ratio=TX_OBSCURATION))
    # The receive terminal carries no obscuration, so the launch obscuration is
    # charged one time only, by the transmit efficiency above.
    analytic_geo = float(geometric_loss_db(range_m, WAIST_M, RX_APERTURE_M,
                                           wavelength=WAVELENGTH_M))
    total_fid2 = loss_db(power_collected, power_launch)
    total_fid0 = analytic_tx + analytic_geo
    rows = [("launch truncation", result.tx_truncation_db, analytic_tx),
            ("geometric spread", result.geometric_loss_db, analytic_geo),
            ("TOTAL, launch to fibre", total_fid2, total_fid0)]

    print("")
    print(f"{'quantity':<26}{'fidelity 2':>12}{'fidelity 0':>12}"
          f"{'difference':>12}")
    print("-" * 62)
    for name, fid2, fid0 in rows:
        print(f"{name:<26}{fid2:>11.3f}dB{fid0:>11.3f}dB{fid2 - fid0:>11.3f}dB")
    print("  all the numbers are losses in positive dB.")
    print("")
    print("  The SPLIT differs, because the two fidelities cut the loss in two")
    print("  different places. The fidelity-0 launch truncation is an on-axis")
    print("  far-field GAIN ratio, and the fidelity-0 geometric spread is the")
    print("  power fraction of the UNtruncated Gaussian in the receive")
    print("  aperture. The fidelity-2 numbers are plain power bookkeeping at")
    print("  each plane. Compare the TOTAL, not the two parts.")

    # ---- the figure ----
    collected_db = loss_db(power_collected, power_launch)
    losses = {
        "1. launch, before the aperture": "reference power, 0.00 dB",
        "2. after the aperture and the obscuration":
            f"{loss_db(power_clipped, power_launch):.2f} dB below launch",
        "3. at the receive plane, 600 km":
            f"beam radius {beam_radius_m:.2f} m, grid side "
            f"{result.grid.size_m * magnification:.1f} m",
        "4. radial cut at the receive plane":
            f"the receive aperture collects {collected_db:.2f} dB below launch",
    }
    draw((stages["launch"], stages["after tx clip"], stages["at rx plane"]),
         losses, beam_radius_m)
    print("")
    print(f"figure saved: {PNG}")
    print("")
    print("  The two fidelities MATCH here, to about 0.01 dB, because the")
    print("  receiver sits in the absolute far field: 600 km is more than 100")
    print("  Rayleigh ranges. The near-field case is where they disagree. Run")
    print("  `python -m examples.waveoptics.terrestrial_stages` for that one:")
    print("  a 1 km link inside a 29 km Rayleigh range, where the two totals")
    print("  differ by more than 8 dB.")
    plt.show()


if __name__ == "__main__":
    main()

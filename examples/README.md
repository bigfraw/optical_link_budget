# olb examples

These scripts are the curated user-facing examples. Each one builds a link
budget and prints the result. Run each script from the repository root as a
module:

    python -m examples.uplink_sim

Read [build_a_link.py](build_a_link.py) first. It shows the shape of the
package: hardware lives on a `Terminal`, the channel holds no hardware, and the
direction sets the transmit and the receive roles. Then read
[custom_budget.py](custom_budget.py) if you must assemble a `Budget` by hand
from the Term factories.

Loss is positive dB. Gain is negative dB.

## The scripts

| File | Link type | Fidelity | Purpose |
| --- | --- | --- | --- |
| [build_a_link.py](build_a_link.py) | uplink + downlink | 0, 1 | Start here. It builds a bistatic station: separate transmit and receive apertures, so one Terminal for each role. |
| [uplink_sim.py](uplink_sim.py) | uplink | 1 | The canonical ground-to-satellite uplink. It prints the itemised budget, the assumptions, the Monte Carlo fade, and an elevation sweep. |
| [downlink_terminal.py](downlink_terminal.py) | downlink | 0, 1 | How the receive front end changes the downlink: a bucket aperture, a fibre with no correction, and a fibre with tip-tilt or adaptive optics. |
| [smf_fidelity_benchmark.py](smf_fidelity_benchmark.py) | downlink | 0 vs 1 | A Term-level benchmark of the two single-mode-fibre coupling fidelities: the cheap analytic mean against the FAST statistical model. It needs `fast-aosim`. |
| [retro_link.py](retro_link.py) | retro | 1 | The retroreflected space link. One ground station transmits the up-leg and receives the return. The budget carries both legs. It needs `fast-aosim`. |
| [terrestrial_link.py](terrestrial_link.py) | terrestrial | 0 | The horizontal ground-to-ground link. It compares a bucket aperture against a single-mode fibre, then sweeps the distance, the compensation, and the receive aperture. |
| [custom_budget.py](custom_budget.py) | any | 1 | How to assemble a `Budget` by hand from the Term factories, when a pre-built budget in `olb/links` does not fit. |

## The fidelity tiers

- **Fidelity 0** is the analytic mean. A Term gives a closed-form mean loss, and
  sometimes a closed-form fade.
- **Fidelity 1** is the statistical model. A Term gives samples, so the Budget
  gives a Monte Carlo fade and a margin.
- **Fidelity 2** is wave optics. The code propagates a real complex field on a
  grid.

Fidelity 0 and fidelity 1 live in these top-level examples. Fidelity 2 lives in
[waveoptics/](waveoptics/). The demonstrations of the Andrews and Phillips
physics layer live in [andrews/](andrews/).

Fidelity 1 does NOT exist for the terrestrial fibre link. The FAST engine is a
far-field, plane-wave-source model, but a near-field finite Gaussian beam needs
a split-step beam-propagation model. So the terrestrial fibre budget gives the
mean coupling loss only. See [terrestrial_link.py](terrestrial_link.py).

## The validation scripts

The owner's cross-check and validation scripts live in
[../validation/](../validation/). They are not curated examples.

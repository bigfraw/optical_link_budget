# olb example scripts

The `examples/` directory holds the curated runnable scripts. Each script builds
a link budget and prints the result. The `validation/` directory holds the
owner's cross-check scripts. This guide describes each one.

## How to run

Run every script from the repository root. Each script uses package-relative
imports, so start it as a module:

    python -m examples.uplink_sim

## Prerequisites

The core scripts need only `olb` and its shared kernels from
`my_analysis_modules`. Set the `MY_ANALYSIS_MODULES` environment variable, or
place that repository at `D:\repos\my_analysis_modules`. Two scripts use the
FAST fidelity-1 single-mode-fibre coupling. That path needs `fast-aosim`.
Install it with `pip install fast-aosim`. Loss is positive dB. Gain is negative
dB.

---

## [uplink_sim.py](../examples/uplink_sim.py)

The canonical ground-to-satellite uplink. It is the best first read. The script
builds one `SpaceScenario`, assembles the uplink budget, and evaluates it.

- API: `uplink_budget(scenario, geom, n_samples=...)`; a `SpaceScenario` with
  `direction="uplink"`; a ground `Terminal` with a `Transmitter`; a space
  `Terminal` with an `Aperture` detector; `CircularOrbit` geometry.
- The script prints the itemised budget, the model assumptions, and the broken
  assumptions. It then runs `budget.monte_carlo(...)` for the fade and the
  margin. The uplink turbulence term is a Monte Carlo model, so the budget uses
  `monte_carlo()`, not the analytic fade.
- Output: an itemised table, an assumptions table, a Monte Carlo mean and fade,
  and an elevation sweep of the 99% margin.
- Run: `python -m examples.uplink_sim`

## [downlink_terminal.py](../examples/downlink_terminal.py)

The satellite transmits. The ground station receives. The script shows how the
receive `Terminal` changes the downlink budget across five front ends.

- API: `downlink_budget(scenario, geom)`; a `SpaceScenario` with
  `direction="downlink"`; receive detectors `Aperture` and `SMF`; a
  `compensation` stack of `TipTilt` and `AO`.
- It compares a power-in-bucket aperture, a single-mode fibre with no
  correction, and a fibre with tip-tilt or adaptive optics. The fibre couples
  only the matched field, so wavefront correction recovers the coupling loss.
- Output: an itemised table per front end, plus a summary table of the coupling
  loss, the total loss, and the 99% fade.
- Run: `python -m examples.downlink_terminal`

## [custom_budget.py](../examples/custom_budget.py)

The lower-level path. It composes the four models into one `Budget` by hand,
instead of a link factory. Use it when a pre-built budget in `olb/links` does
not fit: a link that no factory covers, a Term swap, or a what-if.

- API: `Budget([...], scenario=...)` built from `geometric_loss_term`,
  `slant_extinction_term`, `pointing_loss_term`, and `uplink_turbulence_term`.
  A FAST-free Cn2 profile comes from `get_c2n`, so the script runs anywhere.
- A hand-built Budget does not protect you from a double count. The script sets
  the ground jitter to zero, because it keeps the pointing Term apart from the
  coupled-flux turbulence Term.
- It asserts that the table has all four terms and that the margin is finite.
- Output: an itemised table, the broken assumptions, a Monte Carlo mean, a 99%
  fade, and a 99% margin.
- Run: `python -m examples.custom_budget`

## [build_a_link.py](../examples/build_a_link.py)

The bistatic pattern. The transmit aperture and the receive aperture differ, so
each role is its own `Terminal`. The script builds four terminals and runs the
link both ways.

- API: two `SpaceScenario` objects (one uplink, one downlink) over one shared
  `Channel`; four `Terminal` objects; `uplink_budget` and `downlink_budget`;
  `downlink_budget(..., smf_fidelity="fast")`.
- An aperture is a `Terminal` parameter, so a bistatic station needs a separate
  transmit and receive `Terminal`. Each direction wires in its own pair.
- Output: an itemised table and a Monte Carlo fade for the uplink and the
  downlink. The downlink path needs `fast-aosim`.
- Run: `python -m examples.build_a_link`

## [retro_link.py](../examples/retro_link.py)

The retroreflected link. One ground station transmits the up-leg and receives
the return. The satellite is a passive retroreflector. The budget carries both
legs.

- API: `retro_space_budget(scenario, geom, n_samples=..., smf_fidelity="fast")`;
  a `SpaceScenario` with `direction="retro"`; one bistatic ground `Terminal`
  whose `Transmitter` carries its own `aperture_m` (the beam director) separate
  from the receive telescope.
- The retro direction makes the transmit terminal and the receive terminal the
  same object. The script sweeps three elevations, checks the assumptions, and
  prints the weak-fluctuation regime status per leg. The up-leg and the down-leg
  use different metrics and limits.
- Output: an itemised table and a Monte Carlo fade per elevation, an assumptions
  report, and the per-leg regime status. The path needs `fast-aosim`.
- Run: `python -m examples.retro_link`

## [smf_fidelity_benchmark.py](../examples/smf_fidelity_benchmark.py)

A term-level benchmark. It compares the two single-mode-fibre coupling
fidelities on a no-AO downlink over an elevation sweep.

- API: `downlink_coupling_term(scenario, geom, smf_fidelity=...)` called directly, not
  through a budget; `smf_fidelity="mean"` (the cheap analytic mean, no fade)
  versus `smf_fidelity="fast"` (the FAST fidelity-1 statistical model).
- The FAST reference gives the mean, the quantile, and the deep-fade tail. The
  cheap mean-only model has no fade, so the 99% margin can never come from it.
- Output: a table per elevation with the mean-only loss, the FAST mean, the FAST
  99% loss, the plane-wave scintillation index, and the regime flag. The FAST
  path needs `fast-aosim`.
- Run: `python -m examples.smf_fidelity_benchmark`

## [terrestrial_link.py](../examples/terrestrial_link.py)

The horizontal ground-to-ground family. It is the most involved script. It
builds the terrestrial budget for an aperture and a fibre receiver, then sweeps
distance, compensation, and receive aperture.

- API: `terrestrial_budget(scenario, HorizontalPath(...))`; a
  `TerrestrialScenario` with `near` and `far` terminals; a `TerrestrialChannel`
  with `path_length_m`, `attenuation_db_per_km`, and `cn2`. This family has no
  direction.
- The aperture receiver gets an analytic 99% fade from the horizontal
  Gaussian-beam scintillation term. The fibre receiver is fidelity 0 (mean-only
  coupling), so its budget refuses a fade margin. The script shows the
  aperture-averaging win: a larger receive aperture shrinks the fade.
- Output: two itemised budgets at 3 km, a distance sweep, a compensation sweep
  of the fibre coupling loss, and an aperture-averaging table.
- Run: `python -m examples.terrestrial_link`

---

## The validation scripts ([validation/](../validation/))

The `validation/` directory at the repository root holds the owner's cross-check
and validation scripts. They are not curated user examples. They can be
specific, they can overlap, and they can be rough. Run each one as a module, for
example `python -m validation.uplink_divergence`. For the folder guide, see
[validation/README.md](../validation/README.md).

## [uplink_divergence.py](../validation/uplink_divergence.py)

A trade study. The ground station widens (diverges) the transmit beam on
purpose. The script shows the trade across four divergence values.

- API: `uplink_budget(...)`; `Transmitter(waist_m=..., divergence_rad=...)`;
  `dataclasses.replace` to vary only the divergence per case.
- A wider beam raises the geometric spreading loss but lowers the pointing and
  turbulence loss. The 99% fade tightens most, so a moderate divergence can
  improve the 99% margin. The divergence is the far-field half-angle. It cannot
  be smaller than the diffraction limit.
- Output: an itemised table per case, a Monte Carlo margin per case, and a
  summary table. The script names the divergence with the best 99% margin.
- Run: `python -m validation.uplink_divergence`

## [terrestrial_coupling_jitter.py](../validation/terrestrial_coupling_jitter.py)

A terrestrial single-mode-fibre receiver. It splits the coupling loss into three
pointing mechanisms and keeps the free-space loss apart from the fibre loss:

- The transmit jitter is a free-space loss. The beam misses the far aperture.
- The beam wander is a fibre-coupling loss. The turbulence tilts the arriving
  wavefront, so the focal spot walks off the fibre (the walk-off, contribution
  A).
- The receive jitter is a fibre-coupling loss. It also moves the focal spot
  (the walk-off, contribution B).

- API: `terrestrial_budget`, `terrestrial_smf_walkoff_term`, `terrestrial_smf_coupling_term`,
  and `pointing_loss_term`. The receive terminal is `SMF(optimal_focus=True)`, so
  the focal length comes from the mode field radius and the aperture at `a=1.12`.
- Output: a loss breakdown by mechanism, then three sweeps (over Cn2, the receive
  jitter, and the transmit jitter) that show each mechanism scales on its own.
- Run: `python -m validation.terrestrial_coupling_jitter`

## [mmf_coupling_validation.py](../validation/mmf_coupling_validation.py)

A validation of the multimode-fibre (light-bucket) coupling. It plots the coupled
power against the incident angle for the correct encircled-energy model and for
the old, wrong Gaussian roll-off, at two spot sizes.

- API: `olb.models.coupling.terrestrial._mmf_encircled_efficiency`. The correct model has a
  flat top: a small spot loses almost nothing until it nears the core edge, where
  it collects about half the power (about 3 dB). The old model wrongly lost power
  from zero angle.
- Output: a two-panel PNG (`mmf_coupling_vs_angle.png`, or a path you pass), plus
  the coupled power at half the edge angle and at the edge.
- Run: `python -m validation.mmf_coupling_validation [out.png]`

## [terrestrial_mmf_na.py](../validation/terrestrial_mmf_na.py)

A terrestrial multimode-fibre link that shows the numerical-aperture angular gate.
The focusing cone `NA_optic = (D/2)/f` must stay within the fibre NA, or the fibre
does not guide the steep rays. The spot size and the cone are locked by the
diffraction invariant `w_s * NA_optic = lambda/pi`, so a shorter focal length
tolerates more walk-off but pays a larger gate.

- API: `terrestrial_mmf_coupling_term` with `MMF(numerical_aperture=...)`. The loss splits
  additively: spot-in-core + NA gate + walk-off = mean.
- Output: an NA sweep (the gate turns on below `NA_optic`), then a focal-length
  sweep that shows the etendue trade and a best focal length.
- Run: `python -m validation.terrestrial_mmf_na`

---

## The Andrews foundation suite ([examples/andrews/](../examples/andrews/))

The `examples/andrews/` directory holds a separate suite. Each script
demonstrates the Andrews and Phillips foundation layer
(`olb/turbulence/andrews/`) one topic at a time: the spectra and the scales,
the Gaussian-beam parameters, the scintillation regimes, the irradiance
distributions and fades, the two beam-wander routes, the aperture averaging,
the temporal statistics, the slant paths, and the downlink distribution
selection. Most scripts print book values, not a link margin.

Every equation cites its chapter, its equation number, and its printed page from
Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196. For the per-script
guide and the wired-versus-available status, see the suite README,
[examples/andrews/README.md](../examples/andrews/README.md). Run each script as a
module, for example `python -m examples.andrews.scintillation_regimes`.

---

## The wave-optics suite ([examples/waveoptics/](../examples/waveoptics/))

The `examples/waveoptics/` directory holds six scripts for the fidelity-2 field
propagation layer (`olb/waveoptics/`). Each script propagates a real complex
field on a square grid, prints a table of numbers, and saves a figure next to the
script. The first three have NO turbulence. The last three add the turbulent
split step of `olb/waveoptics/turbulence/`.

The vacuum scripts:

- `terrestrial_stages.py` — the stage-by-stage propagation of a NEAR-FIELD
  terrestrial link. It prints the fidelity-2 numbers against the fidelity-0
  analytic Terms. The two totals disagree by 8.28 dB, because the analytic
  transmit efficiency is a far-field form. The script then adds a retroreflected
  return leg by hand, from three primitives.
- `space_farfield.py` — the FAR-FIELD check on a space link: a 50 mm waist through
  a 100 mm launch aperture with a 0.3 central obscuration, to a 500 mm receive
  aperture at 600 km. A flat grid cannot hold this link, so
  `GridSpec.for_scenario` selects the co-moving route and `propagate_scenario`
  runs `LensFresnel`. The two TOTALS agree to 0.011 dB. The SPLIT does not agree,
  because the fidelity-0 transmit efficiency is an on-axis far-field gain ratio
  and the fidelity-2 number is a power ratio. Compare the total only.
- `grid_artefacts.py` — the deliberate failure. It shows the FFT wrap-around
  artefact of a grid that is too small, against the analytic ABCD route
  (`GForvard`).

The turbulent scripts. Each one runs for about four to five minutes on a desktop,
each one prints the SAMPLING REPORT of its grid and the per-trial wall times, and
each one carries the fixed seed `SEED = 20260826`, so a second run repeats the
first one exactly. Each one saves the figure and opens NO window, because a
blocking window would hold the terminal for minutes.

- `turbulent_terrestrial.py` — a 2 km horizontal link at `Cn2 = 3e-15`
  (`sigma2_R = 0.21`, firmly weak; the script ASSERTS that, because every
  analytic target here is a weak-fluctuation form). It runs 120 snapshots THREE
  times on the same grid, the same screens and the same seeds, and it changes
  only the receive aperture: a 3-pixel pinhole, a 30 mm sampling bucket, and the
  100 mm budget aperture with its single-mode fibre. Headline: the pinhole index
  and the 30 mm bucket index agree with the Dios on-axis form and the Andrews
  weak aperture-averaging factor; the 100 mm bucket does NOT, because it holds 78
  percent of the beam and the split step conserves power; and the fidelity-0
  fibre-coupling Term reads about 2.5 dB MORE loss than the field.
- `turbulent_downlink.py` — a 600 km downlink into a 500 mm obscured fibre
  receiver, at 30, 60 and 90 degrees, 70 snapshots each, `rapid` preset.
  Headline: the aperture scintillation index agrees with the fidelity-0
  plane-wave integral at every elevation, and the fibre coupling does not agree
  with the fidelity-1 FAST Term. The field reads 0.7 dB less loss at 30 degrees
  and 2.9 dB less at the zenith. The script prints the static mode-match floor of
  each model, so the turbulence part can be read alone, and it names the
  candidate causes without picking one.
- `turbulent_uplink_reciprocity.py` — a 600 km uplink at the zenith and at 30
  degrees, 100 snapshots each, `rapid` preset. The satellite is outside the grid,
  so the uplink flux comes from the reciprocity overlap of the propagated
  downlink field with the ground transmit mode (Shapiro,
  DOI 10.1364/JOSA.61.000492). Headline: that loss goes against the Dios
  coupled-flux Monte Carlo of `olb.turbulence.uplink_flux`, and the MEANS agree
  inside 1 dB at both elevations. The 30-degree row is a REPORT, not a test,
  because the coupled-flux model already says `weak_fluctuation_valid = False`
  there. The TAILS are reported, not tested: a field Monte Carlo reaches deeper
  than a parametric lognormal. Both terminals carry zero pointing jitter.

No budget consumes the layer. Run each script as a module, for example
`python -m examples.waveoptics.terrestrial_stages`. For the per-script guide and
the status, see the suite README,
[examples/waveoptics/README.md](../examples/waveoptics/README.md).

---

The four link families map to their example: uplink -> `uplink_sim.py`,
downlink -> `downlink_terminal.py`, retro -> `retro_link.py`, terrestrial ->
`terrestrial_link.py`.

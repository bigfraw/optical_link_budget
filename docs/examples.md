# olb example scripts

The `examples/` directory holds the curated runnable scripts. Each script builds
a link budget and prints the result. The `validation/` directory holds the
owner's cross-check scripts. This guide describes each one.

## How to run

Run every script from the repository root. Each script uses package-relative
imports, so start it as a module:

    python -m examples.uplink_sim

## Prerequisites

The core scripts need only `olb`, which is self-contained (no sibling
repository). Two scripts use the FAST fidelity-1 single-mode-fibre coupling.
That path needs `fast-aosim`.
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
  `downlink_budget(..., fidelity=1)`; `run_fidelity2` plus
  `uplink_budget/downlink_budget(..., fidelity=2, wave=...)` for the wave-optics
  tier (200 trials, threaded).
- An aperture is a `Terminal` parameter, so a bistatic station needs a separate
  transmit and receive `Terminal`. Each direction wires in its own pair.
- Output: an itemised table and a Monte Carlo fade for the uplink and the
  downlink at each fidelity rung. The downlink path needs `fast-aosim`, and the
  fidelity-2 rung runs two Monte Carlo propagations.
- Run: `python -m examples.build_a_link`

## [retro_link.py](../examples/retro_link.py)

The retroreflected link. One ground station transmits the up-leg and receives
the return. The satellite is a passive retroreflector. The budget carries both
legs.

- API: `retro_space_budget(scenario, geom, n_samples=..., fidelity=1)`;
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
  versus `smf_fidelity="fast"` (the FAST fidelity-1 statistical model). The
  coupling Term keeps its `smf_fidelity` parameter: `smf_fidelity="mean"` is the
  fidelity-0 model and `smf_fidelity="fast"` is the fidelity-1 model of the Term.
  At the budget level these are now selected with the whole-path `fidelity=0` and
  `fidelity=1`.
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
example `python -m validation.coupling_checks.uplink_divergence`. For the folder guide, see
[validation/README.md](../validation/README.md).

## [coupling_checks/uplink_divergence.py](../validation/coupling_checks/uplink_divergence.py)

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
- Run: `python -m validation.coupling_checks.uplink_divergence`

## [coupling_checks/terrestrial_coupling_jitter.py](../validation/coupling_checks/terrestrial_coupling_jitter.py)

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
- Run: `python -m validation.coupling_checks.terrestrial_coupling_jitter`

## [coupling_checks/mmf_coupling_validation.py](../validation/coupling_checks/mmf_coupling_validation.py)

A validation of the multimode-fibre (light-bucket) coupling. It plots the coupled
power against the incident angle for the correct encircled-energy model and for
the old, wrong Gaussian roll-off, at two spot sizes.

- API: `olb.models.coupling.terrestrial._mmf_encircled_efficiency`. The correct model has a
  flat top: a small spot loses almost nothing until it nears the core edge, where
  it collects about half the power (about 3 dB). The old model wrongly lost power
  from zero angle.
- Output: a two-panel PNG (`coupling_checks/figures/mmf_coupling_vs_angle.png`, or a path you pass), plus
  the coupled power at half the edge angle and at the edge.
- Run: `python -m validation.coupling_checks.mmf_coupling_validation [out.png]`

## [coupling_checks/terrestrial_mmf_na.py](../validation/coupling_checks/terrestrial_mmf_na.py)

A terrestrial multimode-fibre link that shows the numerical-aperture angular gate.
The focusing cone `NA_optic = (D/2)/f` must stay within the fibre NA, or the fibre
does not guide the steep rays. The spot size and the cone are locked by the
diffraction invariant `w_s * NA_optic = lambda/pi`, so a shorter focal length
tolerates more walk-off but pays a larger gate.

- API: `terrestrial_mmf_coupling_term` with `MMF(numerical_aperture=...)`. The loss splits
  additively: spot-in-core + NA gate + walk-off = mean.
- Output: an NA sweep (the gate turns on below `NA_optic`), then a focal-length
  sweep that shows the etendue trade and a best focal length.
- Run: `python -m validation.coupling_checks.terrestrial_mmf_na`

## [defocus/defocus_sensing.py](../validation/defocus/defocus_sensing.py)

A non-focal-plane (defocused) detector on a terrestrial link. The detector sits
at `z = f + defocus_m`, and the received diverging beam puts the TRUE focus at
`dz_curv` BEYOND the focal plane, so the coupling reads
`dz_eff = defocus_m - dz_curv`.

- API: `terrestrial_mmf_coupling_term`, `terrestrial_smf_walkoff_term`,
  `olb.models.coupling.curvature_focus_shift`, and the
  `olb.links.bidirectional` wrapper.
- Output: a `dz` sweep for a multimode receiver, the spot radius against
  `gaussz` and the geometric blur, the chief-ray tilt lever, and the
  bidirectional demonstration. One fidelity-2 cross-check is guarded, so a
  missing `aotools` does not fail it.
- Run: `python validation/defocus/defocus_sensing.py`
- Write-up: `validation/defocus/fidelity2_mmf_coupling_gap.md`.

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

The `examples/waveoptics/` directory holds eleven scripts for the fidelity-2 field
propagation layer (`olb/waveoptics/`). Each script propagates a real complex
field on a square grid, prints a table of numbers, and saves a figure next to the
script. The first three have NO turbulence. The next three add the turbulent
split step of `olb/waveoptics/turbulence/`. The seventh wires the layer into the
three link budgets. Two draw the focused spot on a multimode-fibre core, one
bins the focal spot onto a tracking-camera pixel grid, and one stores a campaign
of trials on disk.

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

The turbulent scripts. The two space scripts run in about one minute on a
desktop, and the terrestrial one in about three (work package 7 cut the space
screen count from 20 to 5). Each one prints the SAMPLING REPORT of its grid and
the per-trial wall times, and
each one carries the fixed seed `SEED = 20260826`, so a second run repeats the
first one exactly. Each one saves the figure and opens NO window, because a
blocking window would hold the terminal for minutes.

- `turbulent_terrestrial.py` — a 2 km horizontal link at `Cn2 = 3e-15`
  (`sigma2_R = 0.21`, firmly weak; the script ASSERTS that, because every
  analytic target here is a weak-fluctuation form). It runs 300 snapshots THREE
  times on the same grid, the same screens and the same seeds, and it changes
  only the receive aperture: a 3-pixel pinhole, a 30 mm sampling bucket, and the
  100 mm budget aperture with its single-mode fibre. Headline: the pinhole index
  and the 30 mm bucket index agree with the Dios on-axis form and the Andrews
  weak aperture-averaging factor; the 100 mm bucket does NOT, because it holds 78
  percent of the beam and the split step conserves power; and the fidelity-0
  fibre-coupling Term reads about 2.4 dB MORE loss than the field (4.61 dB
  against 2.25 dB). The fidelity-0 Term keeps the printed 4.61 dB. This script
  uses a bare `SMF()` with no coupling optics, so the received-curvature defocus
  charge (see `physics.md` section 6a) does NOT fire, and the Term flags the
  curvature as NOT modelled. That charge applies only to an `SMF` with a focal
  length or `optimal_focus`. The 2.25 dB field value is the Monte Carlo result
  of the default `olb` screen generator and 300 snapshots; the earlier 2.32 dB
  record used 120 snapshots. The horizontal planner takes no `Cn2` layer list,
  so work package 7 did not move this script: it keeps its 9 screens.
- `turbulent_downlink.py` — a 600 km downlink into a 500 mm obscured fibre
  receiver, at 30, 60 and 90 degrees, 70 snapshots each, `rapid` preset, 5
  screens. Headline: the aperture scintillation index agrees with the fidelity-0
  plane-wave integral at every elevation (the ratios are 1.01, 1.19 and 1.28,
  against a 17 percent Monte Carlo error), and the fibre coupling does not agree
  with the fidelity-1 FAST Term. The field reads 2.7 dB less loss at 30 degrees
  and 3.9 dB less at the zenith; on the turbulence part alone, 1.8 to 3.0 dB
  less. The script prints the static mode-match floor of
  each model, so the turbulence part can be read alone, and it names the
  candidate causes without picking one.
- `turbulent_uplink_reciprocity.py` — a 600 km uplink at the zenith and at 30
  degrees, 200 snapshots each, `rapid` preset, 5 screens. The satellite is
  outside the grid,
  so the uplink flux comes from the reciprocity overlap of the propagated
  downlink field with the ground transmit mode (Shapiro,
  DOI 10.1364/JOSA.61.000492). Headline: that loss goes against the Dios
  coupled-flux Monte Carlo of `olb.turbulence.uplink_flux`, and the MEANS agree
  to 0.19 dB at the zenith and 1.05 dB at 30 degrees. The 30-degree row is a
  REPORT, not a test,
  because the coupled-flux model already says `weak_fluctuation_valid = False`
  there. The TAILS are reported, not tested: a field Monte Carlo reaches deeper
  than a parametric lognormal. Both terminals carry zero pointing jitter.

The budget-wiring script:

- `budget_wiring.py` — the whole-path fidelity-2 wiring of the three link
  budgets. Fidelity is a whole-path choice. At `fidelity=2` the entire path is a
  field simulation, and the budget shows TWO wave-optics Terms: a DETERMINISTIC
  vacuum-optics Term (the full no-turbulence loss from launch to detector: launch
  truncation, geometric spread, aperture capture, and vacuum fibre coupling) and
  a STOCHASTIC turbulence Term (the fade). Together they replace the analytic
  geometric, launch-truncation, scintillation, and coupling Terms. Only the
  analytic extinction (molecular absorption) and pointing (mechanical jitter)
  Terms stay. The budget never runs the split step. The caller runs both
  propagations one time with `run_fidelity2(scenario, geometry, ...)`, which
  gives a `Fidelity2Bundle`, then passes the bundle to the budget. The script
  wires all three links: terrestrial single-mode fibre, an uncorrected uplink,
  and a downlink aperture. For the terrestrial link, `fidelity=2` unlocks the
  fade margin that the fidelity-0 mean-only coupling Term refuses.
  - API: `run_fidelity2(scenario, geometry, n_trials=..., seed=..., preset=...,
    threader=...)`, then `terrestrial_budget(..., fidelity=2, wave=bundle)`,
    `uplink_budget(..., fidelity=2, wave=bundle)`, and
    `downlink_budget(..., fidelity=2, wave=bundle)`.
  - It uses the RAPID preset and 200 snapshots. It runs in about 75 seconds on a
    desktop. It carries the fixed seed `SEED = 20260828`.
  - Output: per link, the fidelity-0/1 default for reference, then the
    fidelity-2 vacuum-optics loss, the turbulence mean, the 90% fade, and the
    budget total.
  - Run: `python -m examples.waveoptics.budget_wiring`

The multimode-fibre snapshot scripts. Each one draws the focused spot on a
multimode-fibre (light-bucket) core, turbulent against no turbulence. Each one
runs ONE snapshot of each field (a picture, not a statistics run), and each one
builds the fidelity-2 `waveoptics_mmf_coupling_term`. The two links sit in
OPPOSITE corners of the turbulence. Both use the shared helper
`olb.waveoptics.mmf.focal_intensity`.

- `mmf_core_psf.py` — a 600 km downlink into a ground light-bucket receiver. The
  big aperture gives a large `D/r0`, so the turbulence broadens the focused spot
  and the spot spills past the core (a real loss). The still-atmosphere spot fits
  inside the core. The figure goes to
  `examples/waveoptics/figures/mmf_core_psf.png`.
  - Run: `python -m examples.waveoptics.mmf_core_psf`
- `mmf_core_psf_terrestrial.py` — the terrestrial sibling: a 5 km horizontal link
  at `Cn2 = 5e-15`, into a small 25 mm receiver. The aperture is smaller than one
  coherence cell (`D/r0` about 0.55), so the focused spot stays compact and the
  big core holds it. The turbulence loss lives in the scintillation and the
  wander, not in the spot that spills the core. The figure goes to
  `examples/waveoptics/figures/mmf_core_psf_terrestrial.png`.
  - Run: `python -m examples.waveoptics.mmf_core_psf_terrestrial`

The camera tracking script:

- `camera_tracking.py` — the fidelity-2 focal spot on a tracking camera. A
  600 km downlink at 30 degrees into a 0.7 m ground telescope with a `Camera`
  detector. It propagates a handful of turbulent snapshots
  (`propagate_turbulent_field`), clips each one at the ground aperture, and bins
  the focal spot onto the camera pixels with
  `olb.waveoptics.camera.camera_image`. For each snapshot it prints the
  centroid (in pixels and in microradians, through the plate scale
  theta = x/f), the second-moment spot radius from `spot_metrics`, and the
  fraction of the collected power on the sensor. One still-atmosphere row gives
  the instrument floor. The script builds NO budget change and NO Term: the
  `Camera` is a diagnostic front end. The figure goes to
  `examples/waveoptics/figures/camera_tracking.png`.
  - Run: `python -m examples.waveoptics.camera_tracking`

The campaign script:

- `campaign_demo.py` — a set of turbulent trials on disk, in blocks, read as a
  fidelity-2 BUDGET. It builds 1000 downlink snapshots (600 km, 30 deg, rapid
  preset) as ten blocks of 100, with `workers=4` (ONE warm process pool, and
  each block runs serially inside its process), into
  `examples/waveoptics/_campaign_demo2/`. A second run of the script computes
  NOTHING: the blocks sit on disk. The script then shows the canonical flow: ONE
  scenario and ONE orbit serve the campaign AND the budgets, and the campaign
  goes straight into the `wave` slot. `downlink_budget(scenario, orbit,
  fidelity=2, wave=campaign)` gives 45.60 dB (extinction 0.43 dB + geometric
  spreading 31.48 dB + wave-optics turbulence 13.68 dB, with no pointing
  jitter), and `multi_detector_budgets(scenario, orbit, arms, fidelity=2,
  wave=campaign)` gives 47.15 dB for the `SMF` arm (splitter 1.55 dB) and
  37.23 dB for the `MMF` light-bucket arm (splitter 5.23 dB). The last section
  is DIAGNOSTIC, outside the budget flow: `campaign.recouple(SMF(),
  aperture_m=0.20)` couples the SAME stored fields into a smaller receive
  aperture with no new propagation (mean eta 0.27681). The first run takes about
  50 s, the second about 18 s.
  - API: `Campaign(scenario, geometry, root, seed=..., preset=...,
    block_size=..., sizing_aperture_m=...)`, then `campaign.run(1000, workers=4)`
    and `wave=campaign` on a budget. `campaign.recouple(detector,
    aperture_m=...)` is the diagnostic face.
  - Run: `python examples/waveoptics/campaign_demo.py`

The fidelity-0 and fidelity-1 defaults are unchanged: a budget consumes the
wave-optics layer only when the caller sets `fidelity=2` and gives it a bundle.
Run each script as a module, for example
`python -m examples.waveoptics.terrestrial_stages`. For the per-script guide and
the status, see the suite README,
[examples/waveoptics/README.md](../examples/waveoptics/README.md).

---

## The Schmidt foundation suite ([examples/schmidt/](../examples/schmidt/))

The `examples/schmidt/` directory holds three scripts for the Schmidt numerical
foundation layer (`olb/waveoptics/schmidt/`). Each script puts the book method
against the production wave-optics code, prints a table with a citation on every
row, and saves its figures next to the script. No script changes an `olb`
module.

- `propagator_kernels.py` — the book kernels against the production
  propagators, in three tiers. It bridges the piston phase and the quadrature
  first. Headline: the same algorithm agrees to 1e-10; the one-step and
  two-step Fresnel kernels against the production `Fresnel` agree to 6e-4 in
  the interior for a soft Gaussian and 1.5e-2 for a hard truncation; and the
  two-step kernel against the co-moving `Lens -> LensFresnel -> Convert` recipe
  agrees to 1.7e-3 and 2.3e-2 at a magnification of 247.
- `sampling_and_edges.py` — a gallery of deliberate sampling failures, each one
  paired with the grid that obeys the rule, then the rule checker on the real
  production grids. It also plots the two absorber shapes on one axes.
- `screens_and_turbulence.py` — the screen generators against Eq. (9.44). The
  book subharmonic generator reaches 0.88 to 0.93 of theory over
  `r/r0 = 0.3` to 1.6, and the `aotools` generator of the production layer
  reads 1 to 3 percent above it. The script also proves the factor-4 bridge
  between the two per-screen variance conventions from the live code: 3.9994.

Every equation cites its chapter, its equation number and its printed page from
Schmidt (2010), DOI 10.1117/3.866274. The layer is validation only: no budget,
no Term and no sizer reads it. For the per-script guide, the measured numbers
and the wiring status, see the suite README,
[examples/schmidt/README.md](../examples/schmidt/README.md). Run each script as
a module, for example `python -m examples.schmidt.propagator_kernels`.

---

The four link families map to their example: uplink -> `uplink_sim.py`,
downlink -> `downlink_terminal.py`, retro -> `retro_link.py`, terrestrial ->
`terrestrial_link.py`.

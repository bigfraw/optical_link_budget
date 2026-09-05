# examples/waveoptics — the fidelity-2 field propagation layer

Eleven runnable scripts, in six groups. Each one propagates a real complex field
on a square grid with `olb.waveoptics`, it prints a labelled table of numbers,
and it saves its figures in `examples/waveoptics/figures/`.

- Three vacuum scripts have NO turbulence. They show the vacuum layer at its
  two limits, and they show how the grid breaks.
- Three turbulent scripts add the turbulent split step of
  `olb/waveoptics/turbulence/`, and they compare its statistics with the
  analytic fidelity-0 and fidelity-1 models of the budgets.
- One budget-wiring demo shows the layer inside the three budgets at
  `fidelity=2`.
- Two multimode-fibre demos draw the focused spot on a light-bucket core.
- One camera demo bins the focal spot onto a tracking-camera pixel grid.
- One campaign demo stores the trials on disk and reads them as a budget.

The package is `olb/waveoptics/`. It is a trimmed port of LightPipes
(BSD-3-Clause, see `olb/waveoptics/LIGHTPIPES_LICENSE.txt`), plus two modules
that read an olb scenario: `grid.py` (GridSpec) and `run.py`
(`propagate_scenario`). The turbulent sub-package adds `screens.py`,
`splitstep.py`, `sampling.py`, `run.py`
(`propagate_turbulent_scenario`, and `propagate_turbulent_field` for one
snapshot as a complex field) and `campaign.py` (`Campaign`, the store of trials
on disk).

The physics comes from three sources:

> J. W. Goodman, *Introduction to Fourier Optics*, ISBN 978-0974707723.
> J. D. Schmidt, *Numerical Simulation of Optical Wave Propagation with
> Examples in MATLAB*, SPIE Press (2010). DOI: 10.1117/3.866274.
> A. E. Siegman, *Lasers*, ISBN 978-0935702118.

## How to run

From the repository root:

    python -m examples.waveoptics.terrestrial_stages

The `-m` form is the house idiom of `examples/`. It puts the repository root on
the import path, so `import olb` works. Run each script the same way, with its
own module name. The three vacuum scripts save the figure and then they call
`plt.show()`. A machine with no display shows no window, but the figure file is
already on the disk. Every script that runs a Monte Carlo saves the figure and
opens NO window at all, because a blocking window would hold the terminal.

Every Monte Carlo script keeps its trials in a CAMPAIGN
(`olb.waveoptics.turbulence.Campaign`), the store of trials on disk. Each script
has its own directory under `examples/waveoptics/_campaigns/<script>/<case>/`,
and that tree is git-ignored. The script calls
`campaign.run(n, workers=4, progress=True)`: the blocks run on ONE warm process
pool, and the run prints a progress line for each block. So:

- A SECOND run of a script computes NO trial. It reads the blocks from the disk.
  Only the single-snapshot field pictures propagate again, and they take one to
  five seconds.
- The script reads the record back with `campaign.load(fields=False)`, which
  gives the same trial record that a direct run gives. The statistics code is
  the same code as before.
- The campaign also stores the receive-plane FIELD of each trial, on a disc of
  the receive-aperture radius. `campaign.recollect(aperture_m=...)` gives the
  collected power of a smaller aperture, and `campaign.recouple(detector, ...)`
  gives the coupling efficiency of another detector. Both are post hoc: they
  need no new propagation.

Every seed is explicit. A campaign seeds trial k off `(seed, k)`, so a run of n
trials equals the concatenation of its blocks, trial for trial. Each script
prints how many phase screens it created (one fresh screen stack for each
trial): 2,700 for the terrestrial script (9 screens x 300 trials), 1,055 for the
downlink and 2,000 for the uplink.

The three turbulent scripts and the two multimode-fibre scripts still call
`propagate_turbulent_field` for their PICTURES (`*_field.png`, and the two focal
spots), because a picture needs the wide field and the campaign keeps the
receive disc only. `camera_tracking.py` no longer calls it: it rebuilds the
stored field of each trial on the full grid.

## The vacuum scripts

| Script | What it prints and draws |
| --- | --- |
| `terrestrial_stages.py` | The stage-by-stage propagation of a NEAR-FIELD terrestrial link (waist 120 mm, launch aperture 150 mm, range 1 km, 1550 nm, single-mode-fibre receiver). It prints the fidelity-2 numbers against the fidelity-0 analytic Terms: the two totals disagree by 8.28 dB, because the analytic transmit efficiency is a far-field form. It then adds a retroreflected return leg with three primitives (a clip at a 63.5 mm corner cube, a Fresnel propagation back, a clip at the launch aperture). Figure: `terrestrial_stages.png`, eight intensity panels with physical axes in mm and the aperture of each stage. |
| `space_farfield.py` | The FAR-FIELD check on a space link (waist 50 mm, launch aperture 100 mm with a 0.3 central obscuration, receive aperture 500 mm, 600 km at zenith, 1550 nm). A flat grid cannot hold this link, so `GridSpec.for_scenario` selects the co-moving route and `propagate_scenario` runs `LensFresnel`. It prints the grid, and the fidelity-2 numbers against the fidelity-0 analytic Terms: the two TOTALS agree to 0.011 dB, but the split does not. Figure: `space_farfield.png`, the launch plane, the obscured annulus, the far-field map in log10, and a radial cut against the untruncated Gaussian. |
| `grid_artefacts.py` | The deliberate failure. One Gaussian link (waist 5 mm, range 200 m, 1550 nm) on a proper grid and on a grid of only 2 times the final beam radius. The FFT propagator treats the grid as periodic, so the small grid folds the beam tail back onto the beam. It prints the bucket power of each run against the analytic ABCD route, with the relative error. Figure: `grid_artefacts.png`, four log-intensity panels on one colour scale. |

## The turbulent scripts

Each of the three puts one turbulent split-step Monte Carlo against the
analytic model that a budget already uses. Each one prints the SAMPLING REPORT
of its grid, so the reader sees what the preset achieved, and each one prints
the mean, the smallest and the largest per-trial wall time.

Each turbulent script also saves a SECOND figure, `*_field.png`: the amplitude
and the phase of the received complex field of one snapshot (trial 0 of the
run), with the receive aperture drawn as a dashed ring. It comes from
`propagate_turbulent_field`, a diagnostic entry point that returns the field
without extending the scalar trial record.

| Script | What it prints and draws |
| --- | --- |
| `turbulent_terrestrial.py` | A 2 km horizontal link at Cn2 = 3e-15 (the plane-wave Rytov variance is 0.21, firmly weak; the script ASSERTS that, because every analytic target here is a weak-fluctuation form). It runs ONE campaign of 300 snapshots on the `standard` preset, sized for the 100 mm budget aperture, and it reads the two smaller apertures POST HOC with `campaign.recollect`: a 3-pixel pinhole (8.81 mm) and a 30 mm sampling bucket. So three receive apertures come from one Monte Carlo. The grid is 256 px over 0.7522 m (pixel 2.938 mm), with 9 screens, r0 = 10.66 cm and no sampling warning. The pinhole index reads 0.908 of the Dios on-axis form, and the 30 mm bucket index reads 0.891 of the Andrews aperture-averaging factor; both PASS the loose band. The 100 mm bucket reads 0.268 and it does NOT, and the printed capture fraction says why: it holds 78 percent of the beam, and the split step conserves power. The fibre-coupling Term reads 2.4 dB more loss than the field does (4.61 dB against 2.23 dB). The horizontal planner takes no Cn2 layer list, so work package 7 did not move this script: it keeps its 9 screens, and the run makes 2,700 of them. The first run takes about 97 s at four workers; the second takes 3.2 s and it computes no trial. Figure: `turbulent_terrestrial.png`, the two bucket fades against their lognormals, and the fibre-coupling histogram against the mean-only Term. |
| `turbulent_downlink.py` | A 600 km downlink into a 500 mm obscured fibre receiver, at 30, 60 and 90 degrees. Each elevation is its OWN campaign of 70 snapshots on the `rapid` preset, with 5 screens after work package 7; a one-trial campaign gives the still-atmosphere floor. The grid is 256 px over 2.220 m at 30 degrees (r0 = 12.42 cm), 1.572 m at 60 degrees (r0 = 17.27 cm) and 1.477 m at the zenith (r0 = 18.83 cm), with no sampling warning. The aperture scintillation index sits near the fidelity-0 plane-wave integral at every elevation: the ratios are 1.17, 1.30 and 1.41, against a 17 percent Monte Carlo error, so the first two PASS and the zenith row reads CHECK. The fidelity-2 fibre coupling reads 12.54 dB at 30 degrees, 10.60 dB at 60 degrees and 9.97 dB at the zenith, over a static mode-match floor of 1.95 dB; the turbulence part alone is 10.59, 8.65 and 8.02 dB, and the 99% fade is 30.06, 24.48 and 31.21 dB. The FAST (fidelity-1) comparison did not run here, because `fast-aosim` was not installed (see below). The run makes 1,055 screens; the first run takes 35 s and the second 1.5 s. Figure: `turbulent_downlink.png`, the index against elevation on the analytic curve, and the coupling loss with error bars against the FAST mean and the FAST 99% fade. |
| `turbulent_uplink_reciprocity.py` | A 600 km uplink at the zenith and at 30 degrees. Each elevation is its own campaign of 200 snapshots on the `rapid` preset, with 5 screens. The satellite is outside the grid, so the uplink flux comes from the reciprocity overlap of the propagated downlink field with the ground transmit mode (Shapiro, DOI 10.1364/JOSA.61.000492); the campaign stores that overlap as `eta_turb`. That loss goes against the Dios coupled-flux Monte Carlo of `olb.turbulence.uplink_flux` (3,000 draws), which the uplink Term calls, and which is cheap enough to run each time. The MEANS are 0.95 dB apart at the zenith (field 5.66 dB, flux 6.61 dB) and 1.54 dB apart at 30 degrees (field 9.87 dB, flux 11.41 dB), and the field reads the smaller loss of the two. The coupled-flux model now reports `weak_fluctuation_valid = True` on BOTH rows, so both are inside the regime of its own Rytov model. The TAILS are reported, not tested: a field Monte Carlo reaches deeper than a parametric lognormal (99% fade: field 21.15 dB against flux 18.36 dB at the zenith, 27.72 dB against 28.22 dB at 30 degrees). Both terminals carry zero pointing jitter, because the coupled-flux model folds a jitter into the same wander variance and the overlap does not. The run makes 2,000 screens; the first run takes 69 s and the second 1.9 s. Figure: `turbulent_uplink_reciprocity.png`, the two loss distributions on one axis for each elevation. |

### The FAST comparison of the downlink script

`turbulent_downlink.py` compares its fibre coupling with the fidelity-1 FAST
Term. `fast-aosim` was NOT installed for the run above, so the FAST columns did
not run and the FAST numbers of this README could not be refreshed. The earlier
run read 2.7 dB less loss at 30 degrees and 3.9 dB less at the zenith, and 1.8
to 3.0 dB less on the turbulence part alone. Do not take those numbers as
current.

The current finding on that gap comes from `validation/waveoptics_vs_fast/`, not
from this example: the gap is an OUTER-SCALE artifact. At the physical
L0 = 25 m, uncorrected and like for like, FAST and the field AGREE to about
0.3 dB from 20 to 90 degrees. The old gap reappears only at the grid-dependent
L0 = inf, because FAST is more sensitive to L0.

## What the three turbulent scripts show together

The layer reproduces the analytic models where the analytic assumptions hold,
and it departs from them where they do not:

- the scintillation indices agree, at every elevation and on the horizontal
  path, whenever the receive aperture SAMPLES the beam;
- the aperture-averaged form fails when the aperture nearly HOLDS the beam,
  because the split step conserves power and the analytic factor does not know
  about that;
- on the horizontal path the fidelity-0 fibre-coupling Term reads about 2.4 dB
  MORE loss than the field. On the space downlink the FAST comparison is
  matched to about 0.3 dB once the outer scale is matched (see the note above).

The terrestrial coupling gap is the reason no fidelity-2 Term is wired by
default. It moves budget numbers, so the default is an owner decision.

## What the three vacuum scripts show together

`terrestrial_stages.py` and `space_farfield.py` show the layer at work, at the
two limits: a near-field link where the analytic total FAILS, and a far-field
link where it holds to 0.011 dB. `grid_artefacts.py` shows how the layer
breaks. Read the third one before you trust a number from the other two on a
new link:

- The grid side must hold about 8 times the largest beam radius on the path.
  `GridSpec.for_scenario` applies that rule with a guard factor of 4.
- The range must stay below `forvard_max_z(grid, wavelength) = N*dx^2/lambda`.
  `GridSpec.for_scenario` WARNS above that limit for a flat grid. It does not
  raise.
- A clipped field is not a pure Gaussian, so the exact ABCD route (`GForvard`)
  refuses it. `propagate_scenario` therefore selects `Fresnel` after a hard
  launch clip, and `GForvard` after an almost clean one.

## The co-moving route

A space link makes the beam grow by a factor of 100 or more. A flat grid must
then hold the far-field beam AND resolve the launch aperture, which no
practical pixel count does. `GridSpec.for_scenario` therefore tries the flat
grid first, and it falls back to the SCALED (co-moving) grid: the grid starts
at the launch plane, and it grows with the beam by the magnification
m = w(z)/w(0). `GridSpec.scaled` records the choice, and `propagate_scenario`
reads it to run the three-call lens recipe of `olb/waveoptics/lenses.py`
(`Lens`, `LensFresnel`, `Convert`). See Schmidt, DOI 10.1117/3.866274, Ch. 7.
`GridSpec.for_scenario` warns only when NEITHER route resolves the apertures.

## The budget-wiring demo

| script | what it shows |
| --- | --- |
| `budget_wiring.py` | The wave-optics layer wired into the three budgets at `fidelity=2`. Each link keeps its trials in its own `Campaign` under `_campaigns/budget_wiring/`, and a Campaign IS a fidelity-2 wave record: it goes straight into the `wave` slot of the budget. So the calls are `Campaign(scenario, geometry, root, seed=20260828, preset="rapid", block_size=50)`, then `campaign.run(200, workers=4, progress=True)`, then `terrestrial_budget/uplink_budget/downlink_budget(fidelity=2, wave=campaign)`. Each fidelity-2 budget shows TWO Terms: a deterministic geometric loss and a stochastic turbulence Term. The terrestrial SMF link gives a vacuum-optics (wave) loss of 7.996 dB, a turbulence mean of 3.646 dB, a 90% fade of 8.456 dB, and a total of 13.142 dB, against a fidelity-0 mean-only coupling reference of 8.125 dB; fidelity 2 UNLOCKS the fade margin that the mean-only Term refuses. The uncorrected uplink at 60 degrees gives an analytic geometric loss of 37.691 dB, a turbulence mean of 9.386 dB, a 90% fade of 18.646 dB, and a total of 48.927 dB. The downlink aperture at 30 degrees gives 48.568 dB analytic geometric, a turbulence mean of 0.039 dB, a 90% fade of 1.373 dB, and a total of 49.042 dB, against a fidelity-1 scintillation reference of 0.1162 dB. All three run on a 256 px grid with 5 screens (`rapid`) and 200 trials in four blocks of 50. The cold run takes 112 s in total (34, 39 and 39 s); the second run takes about 1 s and it recomputes no block. The stores hold 8.5 MB (terrestrial), 1.4 MB (downlink) and 140 kB (uplink). Each default budget (fidelity 0/1) is unchanged. |

## The multimode-fibre coupling demo

| script | what it shows |
| --- | --- |
| `mmf_core_psf.py` | The fidelity-2 multimode-fibre (light-bucket) coupling Term, and the focused spot on the core. A 600 km downlink into a ground light bucket at 30 degrees (r0 = 12.42 cm, D/r0 = 5.63) runs a 60-snapshot campaign on the `rapid` preset (256 px over 2.664 m, 5 screens). It builds `olb.models.waveoptics.waveoptics_mmf_coupling_term` from the stored trials and it prints the three Term faces: the mean coupling loss 0.771 dB, the 99% fade 1.758 dB, and the static encircled-energy floor 0.141 dB (a second campaign of one trial on flat screens). So the turbulence part of the mean is 0.630 dB. This is the ONLY statistical MMF model in olb; there is no analytic or FAST sibling. It then focuses one turbulent snapshot and one still snapshot on the SAME grid, and it draws the focal-plane intensity on the core for each, with the core edge as a dashed circle. With a 3.69 m focal length, a 50 um core and NA 0.20, the still spot holds 0.968 of the power (0.14 dB) and the turbulent snapshot holds 0.874 (0.59 dB). So the light bucket is FORGIVING here too: the broadened spot mostly stays inside the core. The panel titles give the capture fraction, so the picture and the numbers agree. It uses the shared focal helper `olb.waveoptics.mmf.focal_intensity`, so it duplicates no FFT. The first run takes 12 s and the second 2 s. Figure: `mmf_core_psf.png`. |
| `mmf_core_psf_terrestrial.py` | The TERRESTRIAL sibling of `mmf_core_psf.py`, in the OPPOSITE corner. A 5 km horizontal link at Cn2 = 5e-15 into a SMALL 25 mm receiver with a 50 um-core light bucket. The regime points two ways: the aperture is smaller than one coherence cell (r0 = 4.52 cm, so D/r0 = 0.55, a mild pupil phase), but the scintillation is strong (sigma_R^2 = 1.90). So the loss is FADE-dominated, not a spot that spills the core. This is a spot PICTURE, not a statistics run: the Term record is a campaign of ONE trial on the `rapid` preset (512 px over 0.994 m, 14 screens), and the still-atmosphere floor is a second one-trial campaign, so it draws NO fade distribution. It builds the same `waveoptics_mmf_coupling_term` from that single trial (which proves the Term is buildable for a terrestrial MMF) and it prints the single-snapshot coupling loss 0.287 dB against the static floor 0.256 dB, so the turbulence part of the coupling is 0.031 dB. The two focal spots (one turbulent, one still) look ALIKE, and both stay inside the big core (0.943 still against 0.936 turbulent). That is the physics result: the light bucket is forgiving when D/r0 < 1. The script does NOT shrink the core to manufacture a contrast, and it prints faithfully whether the turbulent snapshot spills the core. The first run takes 9 s and the second 5 s (the pictures). Figure: `mmf_core_psf_terrestrial.png`. |

## The camera tracking demo

| script | what it shows |
| --- | --- |
| `camera_tracking.py` | The fidelity-2 focal spot on a tracking camera: pixels, centroid and jitter. A 600 km downlink at 30 degrees into a 0.7 m ground telescope with a `Camera` detector. It runs five turbulent snapshots as a CAMPAIGN on a 3-times zoomed grid (768 px over 7.993 m, 5 screens, r0 = 12.42 cm, D/r0 = 5.63), it REBUILDS the stored field of each trial on the full grid, it clips each one at the ground aperture, and it bins the focal spot onto the camera pixels with `olb.waveoptics.camera.camera_image`. It no longer calls `propagate_turbulent_field`. For each snapshot it prints the measured centroid (in pixels and in microradians on the sky, through the plate scale theta = x/f), the second-moment spot radius, and the fraction of the collected power on the sensor. With a 21.28 m focal length the still spot is 30.0 um (2.0 px) and the plate scale is 0.705 urad/px. One STILL-atmosphere row (a second one-trial campaign) gives the instrument floor: its centroid sits 0.035 px off the axis and its rms radius is 10.86 px. The five turbulent centroids run from -8.9 px to 9.2 px, the rms radius grows to 17.0 to 18.7 px, and all the power stays on the sensor (1.0000). The centroid scatter is 3.80 urad in x and 4.53 urad in y. The script builds NO budget change and NO Term: the `Camera` is a diagnostic front end (see `olb.terminal.Camera`; the power-to-pixel-brightness model is backlog 2-W3). The first run takes 7 s and the second 0.8 s. Figure: `camera_tracking.png`. |

## The campaign demo

| script | what it shows |
| --- | --- |
| `campaign_demo.py` | The CANONICAL minimal flow: one scenario, one geometry, one campaign, then the budgets. It builds 1000 downlink snapshots (600 km, 30 degrees, `rapid` preset) as ten blocks of 100 on a warm pool of four processes, then it reads the store two ways: `downlink_budget(scenario, orbit, fidelity=2, wave=campaign)`, and `multi_detector_budgets(scenario, orbit, arms, fidelity=2, wave=campaign)` for a two-arm split. A second run computes NOTHING. It ends with one DIAGNOSTIC section: `campaign.recouple(SMF(), aperture_m=0.20)` couples the SAME stored fields into a smaller receive aperture, with no new propagation. Run it with `python examples/waveoptics/campaign_demo.py`. |

## Status: wave optics is WIRED at `fidelity=2` (both the turbulent and vacuum layers)

The `olb/waveoptics/` package is complete and each module holds a self-check
(`python -m olb.waveoptics.run`, `python -m olb.waveoptics.grid`,
`python -m olb.waveoptics.turbulence.run`, and so on).
It agrees with the fidelity-0 analytic Terms in the far field with a light
truncation, to 0.02 dB. On the 600 km space link of `space_farfield.py`, with
a hard truncation and a central obscuration, the two TOTALS agree to 0.011 dB.

Compare the TOTAL, not the two parts. The two fidelities cut the loss in two
different places: the fidelity-0 launch truncation is an on-axis FAR-FIELD
gain ratio, and the fidelity-0 geometric spread is the power fraction of the
UNtruncated Gaussian in the receive aperture. The fidelity-2 numbers are plain
power bookkeeping at each plane. The product of the fidelity-0 pair is the
collected power fraction, so the totals compare, but the split does not.

**BOTH layers are now wired at `fidelity=2`** (2026-08-28,
`olb/models/waveoptics.py` plus the four budgets). A fidelity-2 budget shows the
DETERMINISTIC geometric Term and the STOCHASTIC turbulence Term (from the split
step); only the analytic extinction and pointing Terms stay. The budget never
runs the sim: the caller gives it a `wave` record.

There are TWO ways to make that record, and both are supported:

- `olb.models.waveoptics.run_fidelity2` is the ONE-CALL API. It runs the trials
  in memory and it gives a `Fidelity2Bundle`.
- A `Campaign` IS a fidelity-2 wave record too. It stores the trials on disk, it
  resumes, and it goes straight into the `wave` slot. Every example in this
  suite now uses the campaign store, because a second run then computes nothing.
  `campaign_demo.py` is the canonical minimal flow.

Every default budget (fidelity 0/1) is unchanged. STILL owner-gated:

- whether wave optics ever becomes the DEFAULT. The terrestrial script shows how
  much that would move a total: the scintillation stays put, and the horizontal
  fibre coupling moves by about 2.4 dB (the field reads the smaller loss). On the
  space downlink, `validation/waveoptics_vs_fast/` shows that the FAST gap
  closes to about 0.3 dB once the outer scale matches. Until the owner resolves
  that reference-model gap, `fidelity=2` is opt-in — you A/B the field against
  the fidelity-0/1 incumbents without perturbing a published total.
- an AUTOMATIC fidelity selector is the next owner-requested step.

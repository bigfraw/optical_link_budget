# examples/waveoptics — the fidelity-2 field propagation layer

Six runnable scripts, in two groups. Each one propagates a real complex field
on a square grid with `olb.waveoptics`, it prints a labelled table of numbers,
and it saves its figures next to the script.

- The first three scripts have NO turbulence. They show the vacuum layer at its
  two limits, and they show how the grid breaks.
- The last three scripts add the turbulent split step of
  `olb/waveoptics/turbulence/`, and they compare its statistics with the
  analytic fidelity-0 and fidelity-1 models of the budgets.

The package is `olb/waveoptics/`. It is a trimmed port of LightPipes
(BSD-3-Clause, see `olb/waveoptics/LIGHTPIPES_LICENSE.txt`), plus two modules
that read an olb scenario: `grid.py` (GridSpec) and `run.py`
(`propagate_scenario`). The turbulent sub-package adds `screens.py`,
`splitstep.py`, `sampling.py` and `run.py`
(`propagate_turbulent_scenario`, and `propagate_turbulent_field` for one
snapshot as a complex field). `threader.py` (`Threader`) runs the independent
trials across threads.

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
already on the disk. The three turbulent scripts save the figure and they open
NO window at all, because each one runs for minutes and a blocking window would
hold the terminal.

Each turbulent script runs the trials across threads with a `Threader`, and it
prints a progress line for each block of trials. The trials are independent
snapshots, and the FFT of the split step releases the GIL, so the threads give
a real speed-up: the terrestrial script drops from about six minutes to about
one on a desktop. Every seed is explicit, so a second run — threaded or not —
repeats the first one exactly. Each script also prints how many phase screens
it created (one fresh screen stack for each trial): about 3,200 for the
terrestrial script, 4,200 for the downlink, and 4,000 for the uplink.

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
the mean, the smallest and the largest per-trial wall time. (A per-trial time
holds the thread contention, so the sum of the per-trial times runs past the
real wall time of a threaded run.)

Each turbulent script also saves a SECOND figure, `*_field.png`: the amplitude
and the phase of the received complex field of one snapshot (trial 0 of the
run), with the receive aperture drawn as a dashed ring. It comes from
`propagate_turbulent_field`, a diagnostic entry point that returns the field
without extending the scalar trial record.

| Script | What it prints and draws |
| --- | --- |
| `turbulent_terrestrial.py` | A 2 km horizontal link at Cn2 = 3e-15 (the plane-wave Rytov variance is 0.21, firmly weak; the script ASSERTS that, because every analytic target here is a weak-fluctuation form). It runs 120 snapshots THREE times on the same screens and the same seeds, and it changes only the receive aperture: a 3-pixel pinhole, a 30 mm sampling bucket, and the 100 mm budget aperture with its single-mode fibre. The pinhole index and the 30 mm bucket index agree with the Dios on-axis form and the Andrews aperture-averaging factor. The 100 mm bucket does NOT, and the printed capture fraction says why: it holds 78 percent of the beam, and the split step conserves power. The fibre-coupling Term reads 2.3 dB more loss than the field does (4.61 dB against 2.32 dB). Work package 7 did not move this script: the horizontal planner takes no Cn2 layer list, so it keeps its 9 screens. Figure: `turbulent_terrestrial.png`, the two bucket fades against their lognormals, and the fibre-coupling histogram against the mean-only Term. |
| `turbulent_downlink.py` | A 600 km downlink into a 500 mm obscured fibre receiver, at 30, 60 and 90 degrees, 70 snapshots each, 5 screens after work package 7. The aperture scintillation index agrees with the fidelity-0 plane-wave integral at every elevation: the ratios are 1.01, 1.19 and 1.28, against a 17 percent Monte Carlo error. The fibre coupling does not agree with the fidelity-1 FAST Term: the field reads 2.7 dB less loss at 30 degrees and 3.9 dB less at the zenith, and 1.8 to 3.0 dB less on the turbulence part alone. The script prints the static mode-match floor of each model, so the turbulence part can be read alone, and it names the candidate causes without picking one. Figure: `turbulent_downlink.png`, the index against elevation on the analytic curve, and the coupling loss with error bars against the FAST mean and the FAST 99% fade. |
| `turbulent_uplink_reciprocity.py` | A 600 km uplink at the zenith and at 30 degrees, 200 snapshots each, 5 screens after work package 7. The satellite is outside the grid, so the uplink flux comes from the reciprocity overlap of the propagated downlink field with the ground transmit mode (Shapiro, DOI 10.1364/JOSA.61.000492). That loss goes against the Dios coupled-flux Monte Carlo of `olb.turbulence.uplink_flux`, which the uplink Term calls. The MEANS agree to 0.19 dB at the zenith and 1.05 dB at 30 degrees, and the field reads the smaller loss of the two. The 30-degree row is a REPORT, not a test: the coupled-flux model already says `weak_fluctuation_valid = False` there, so it is outside the regime of its own Rytov model. The TAILS are reported, not tested: a field Monte Carlo reaches deeper than a parametric lognormal. Both terminals carry zero pointing jitter, because the coupled-flux model folds a jitter into the same wander variance and the overlap does not. Figure: `turbulent_uplink_reciprocity.png`, the two loss distributions on one axis for each elevation. |

## What the three turbulent scripts show together

The layer reproduces the analytic models where the analytic assumptions hold,
and it departs from them where they do not:

- the scintillation indices agree, at every elevation and on the horizontal
  path, whenever the receive aperture SAMPLES the beam;
- the aperture-averaged form fails when the aperture nearly HOLDS the beam,
  because the split step conserves power and the analytic factor does not know
  about that;
- the fibre-coupling models disagree by 2 to 4 dB, and the field reads LESS
  loss than both the fidelity-0 Term and the fidelity-1 FAST Term.

That last gap is the reason no fidelity-2 Term is wired. It moves budget
numbers, so the default is an owner decision.

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

## Status: built and self-checked, but NOT wired

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

**NO budget consumes it.** The layer builds no Term, and it changes no budget
number. A fidelity-2 Term is an owner-gated later step, for two reasons:

- the vacuum part is the no-turbulence validator that flags the near-field and
  far-field limits of the analytic Terms;
- a fidelity-2 Term would move the totals of an existing budget, so the owner
  must decide the default. The three turbulent scripts show HOW MUCH it would
  move them: the scintillation stays put, and the fibre coupling moves by 1 to
  3 dB.

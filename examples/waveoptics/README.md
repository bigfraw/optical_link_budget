# examples/waveoptics — the fidelity-2 field propagation layer

Three runnable scripts. Each one propagates a real complex field on a square
grid with `olb.waveoptics`, it prints a labelled table of numbers, and it saves
a figure next to the script.

The package is `olb/waveoptics/`. It is a trimmed port of LightPipes
(BSD-3-Clause, see `olb/waveoptics/LIGHTPIPES_LICENSE.txt`), plus two modules
that read an olb scenario: `grid.py` (GridSpec) and `run.py`
(`propagate_scenario`).

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
own module name. Each script saves its figure and then it calls `plt.show()`.
A machine with no display shows no window, but the figure file is already on
the disk.

## The scripts

| Script | What it prints and draws |
| --- | --- |
| `terrestrial_stages.py` | The stage-by-stage propagation of a NEAR-FIELD terrestrial link (waist 120 mm, launch aperture 150 mm, range 1 km, 1550 nm, single-mode-fibre receiver). It prints the fidelity-2 numbers against the fidelity-0 analytic Terms: the two totals disagree by 8.28 dB, because the analytic transmit efficiency is a far-field form. It then adds a retroreflected return leg with three primitives (a clip at a 63.5 mm corner cube, a Fresnel propagation back, a clip at the launch aperture). Figure: `terrestrial_stages.png`, eight intensity panels with physical axes in mm and the aperture of each stage. |
| `space_farfield.py` | The FAR-FIELD check on a space link (waist 50 mm, launch aperture 100 mm with a 0.3 central obscuration, receive aperture 500 mm, 600 km at zenith, 1550 nm). A flat grid cannot hold this link, so `GridSpec.for_scenario` selects the co-moving route and `propagate_scenario` runs `LensFresnel`. It prints the grid, and the fidelity-2 numbers against the fidelity-0 analytic Terms: the two TOTALS agree to 0.011 dB, but the split does not. Figure: `space_farfield.png`, the launch plane, the obscured annulus, the far-field map in log10, and a radial cut against the untruncated Gaussian. |
| `grid_artefacts.py` | The deliberate failure. One Gaussian link (waist 5 mm, range 200 m, 1550 nm) on a proper grid and on a grid of only 2 times the final beam radius. The FFT propagator treats the grid as periodic, so the small grid folds the beam tail back onto the beam. It prints the bucket power of each run against the analytic ABCD route, with the relative error. Figure: `grid_artefacts.png`, four log-intensity panels on one colour scale. |

## What the three scripts show together

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
(`python -m olb.waveoptics.run`, `python -m olb.waveoptics.grid`, and so on).
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

- the layer has NO turbulence today. It is the no-turbulence validator that
  flags the near-field and far-field limits of the analytic Terms;
- a fidelity-2 Term would move the totals of an existing budget, so the owner
  must decide the default.

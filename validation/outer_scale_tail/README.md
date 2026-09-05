# The fidelity-2 SMF outer-scale fade-tail study (2-P5, 2-I3, 0-W4)

This study answers two questions with one matched-seed design:

1. **2-P5.** Does the assumed outer scale `L0` bias the fidelity-2 single-mode-
   fibre (SMF) fade tail of a downlink, and by how much? The default screens run
   `L0 = inf`, but a finite grid holds no cell larger than about 27 times the
   grid side, so `L0 = inf` is not physical: the SMF loss then depends on the
   grid size. This study turns the earlier "order of 2 dB at p5" ESTIMATE into a
   measured number.
2. **2-I3.** Is the shipped `rapid` preset a safe DEFAULT? It compares `rapid`
   (5 screens, a coarse grid) against a well-resolved reference at the same
   outer scales.

It is VALIDATION ONLY: it reads the production layer and changes no `olb` module.

## The method: a matched-seed L0 pair, crossed with the preset

The study builds a 2 x 2 of `config` x `L0`:

- **config** `ref` = the well-resolved reference (the standard preset pinned to
  15 screens, 1024 px) and `rapid` = the shipped rapid preset (5 screens, its
  own coarser grid).
- **L0** = `inf` (the grid-dependent Kolmogorov limit) and `25 m` (a finite von
  Karman outer scale).

WITHIN a config the two L0 campaigns share ONE grid, ONE screen plan and ONE
seed, so only the outer-scale filter of the screens changes. That is a clean
matched-seed pair: a per-quantile difference is the outer scale alone. ACROSS
configs the grids and screen counts differ, so `rapid` against `ref` is a
QUANTILE comparison, not matched-seed.

NO CORRECTION: fidelity 2 applies no tip-tilt removal and no adaptive optics
(backlog 2-AO), so the SMF power is the raw uncorrected atmosphere. The SMF tail
is where the outer scale bites most, because the fibre pays the full received
tilt.

The scenario is the hero downlink: 1550 nm, a 700 mm ground SMF terminal, a
100 mm space terminal at 500 km, `cn2_ground = 1.7e-14`, wind 21 m/s.

## The runs and the results

Two elevations, 1000 trials for each case.

**2-P5, the outer-scale effect (reference config, matched-seed, `inf -> 25 m`):**

| quantity | 30 deg | 20 deg |
| --- | --- | --- |
| SMF p5 delta | -2.49 dB (3.0 sigma) | -2.83 dB (2.4 sigma) |
| SMF p1 delta | -3.13 dB (1.3 sigma) | -2.83 dB (1.0 sigma) |
| point p5 delta | +0.06 dB (not resolved) | -0.08 dB (not resolved) |

The finite outer scale gives LESS fade (the safe direction), by about 2.5 to
2.8 dB at p5, ROBUST across elevation. The POINT (centre-pixel) fade does not
move, which confirms the mechanism is the fibre TILT, not scintillation. So the
`L0 = inf` default is about 2.5 to 2.8 dB PESSIMISTIC on the SMF p5 tail.

**2-I3, rapid against the reference at the physical `L0 = 25 m`** (rapid minus
reference, in dB, as (30 deg / 20 deg) pairs):

| quantile | delta (30 deg / 20 deg) |
| --- | --- |
| mean | +0.09 / -0.02 |
| p50 | -0.25 / -0.25 |
| p10 | +0.22 / +0.86 |
| p5 | +0.33 / +0.19 |
| p1 | +0.03 / +0.08 |

At `L0 = 25 m` rapid and the reference are NEAR-IDENTICAL: the mean, p50, p5 and
p1 agree inside about 0.3 dB at both elevations, although rapid uses 5 screens
and 512 px against the reference 15 screens and 1024 px. The one wrinkle is p10
at 20 deg, where rapid reads +0.86 dB more loss (about 1.6 sigma, the safe
direction); at 30 deg even p10 agrees. So rapid IS supported as the default at
the 25 m outer scale on this scenario, with a footnote to watch p10 at low
elevation. This is TWO catalogue points (30 and 20 deg, downlink, SMF); the
strong-Cn2 and terrestrial rows are still needed before the blanket switch.

The `L0 = inf` arm is NOT the operating point (it is grid-dependent). It is kept
only as the baseline that MEASURES the outer-scale bias.

## The decision

Owner decision (2026-09-05): run fidelity-2 sims at a FIXED `L0 = 25 m`, not
`inf`. This is the working site outer scale until 0-W4 threads an explicit `L0`
to the analytic tilt Terms as well.

## Run it

From the repository root:

    python -m validation.outer_scale_tail.outer_scale_tail
    python -m validation.outer_scale_tail.outer_scale_tail --elevation 20
    python -m validation.outer_scale_tail.outer_scale_tail --L0 inf 25 --configs ref rapid
    python -m validation.outer_scale_tail.outer_scale_tail --analyse-only

It writes a results JSON and a run log for each elevation (`_el30`, `_el20`),
and the figures to `figures/`. The campaign blocks go to `campaigns/` (gitignored).

## Sources

- Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196. Ch. 8 the Rytov
  variance; Ch. 12 the Hufnagel-Valley profile.
- Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.50), the von Karman phase
  spectrum with the outer-scale term.
- Fried, DOI 10.1364/JOSA.56.001372. The Fried parameter r0.

# olb validation scripts

This folder holds the owner's cross-check and validation scripts. They are not
curated user examples. Each one checks one model against another model, or
against a known result.

These scripts can be specific, they can overlap, and they can be rough. Read
[../examples/](../examples/) first if you want the curated set.

Each study has its own subfolder. A subfolder keeps its scripts, its results
JSON files and its run logs at its top level, and its figures in a `figures/`
subfolder. A script that draws a figure writes it to `figures/`, and it makes
that folder if the folder is absent.

Run each script from the repository root as a module:

    python -m validation.coupling_checks.uplink_divergence

## coupling_checks/

Four small trade studies and coupling checks. Each one is independent.

| File | Purpose |
| --- | --- |
| [coupling_checks/uplink_divergence.py](coupling_checks/uplink_divergence.py) | A trade study. It widens the uplink transmit beam on purpose, then it finds the divergence with the best 99% margin. |
| [coupling_checks/terrestrial_coupling_jitter.py](coupling_checks/terrestrial_coupling_jitter.py) | It splits the terrestrial single-mode-fibre coupling loss into three pointing mechanisms, then it sweeps each one. |
| [coupling_checks/terrestrial_mmf_na.py](coupling_checks/terrestrial_mmf_na.py) | It shows the numerical-aperture angular gate of a terrestrial multimode-fibre link, then it sweeps the focal length. |
| [coupling_checks/mmf_coupling_validation.py](coupling_checks/mmf_coupling_validation.py) | It plots the multimode-fibre coupled power against the incident angle, for the correct encircled-energy model and for the old, wrong Gaussian roll-off. It writes `figures/mmf_coupling_vs_angle.png`. |

## uplink_sigma2i/

The uplink fidelity-1 against fidelity-2 scintillation-index investigation.
These scripts measure the fidelity-1 Dios coupled-flux uplink against the
fidelity-2 wave-optics field solve. The full write-up is
[uplink_sigma2i/UPLINK_SIGMA2I_INVESTIGATION.md](uplink_sigma2i/UPLINK_SIGMA2I_INVESTIGATION.md)
(RESOLVED 2026-08-28: fidelity 1 over-predicts sigma2_I by 2x to 7x for a FILLED
launch, because the unsaturated Dios off-axis term runs past its validity at the
beam edge; fidelity 2 is the trustworthy leg).

| File | Purpose |
| --- | --- |
| [uplink_sigma2i/uplink_farfield_reciprocity.py](uplink_sigma2i/uplink_farfield_reciprocity.py) | Mode-matched fidelity-1 versus fidelity-2 uplink scintillation through the reciprocity far-field map. It measures EACH ingredient of the fidelity-1 model on its own (on-axis sigma2_I, wander variance, long/short-term widths, beam-frame index), not only the headline number. |
| [uplink_sigma2i/uplink_obscuration_dios_vs_waveoptics.py](uplink_sigma2i/uplink_obscuration_dios_vs_waveoptics.py) | How far the fidelity-1 Dios uplink can be trusted through a centrally obscured (annular) launch pupil. Dios reads the launch through one number (the waist w0), so its sigma2_I is flat in the obscuration ratio; the sweep straddles the point where the obscuration blocks the beam core. It writes `figures/uplink_obscuration*.png`. |
| [uplink_sigma2i/uplink_obscuration_farfield.py](uplink_sigma2i/uplink_obscuration_farfield.py) | The fidelity-2 VACUUM far-field spot at the satellite for each obscuration radius. It is the picture behind the mean-loss curve: as the obscuration grows past the waist, the Gaussian core is blocked and the surviving ring paints a broad Airy-like pattern. It writes `figures/uplink_obscuration_farfield.png`. |
| [uplink_sigma2i/dios_fig5_replication.py](uplink_sigma2i/dios_fig5_replication.py) | Replicate Dios et al. 2004, Fig. 5 (DOI 10.1364/AO.43.003866): uplink log-amplitude variance against transmit waist, GEO at 0.84 um, 90 and 30 deg. Puts both olb legs on the paper's case (`--fid1` for the analytic leg only). If the vendored coupled-flux curve overlays the paper's line, the port is faithful. |
| [uplink_sigma2i/dios_fig5_plot.py](uplink_sigma2i/dios_fig5_plot.py) | Plot the Dios Fig. 5 replication (fidelity-1 curves, fidelity-2 points). Run after `dios_fig5_replication.py`. It writes `figures/dios_fig5_replication.png`. |

## defocus/

The non-focal-plane (defocused) detector study, and the fidelity-0 against
fidelity-2 terrestrial coupling gap that it resolved. The detector sits at
`z = f + defocus_m`; the received beam is a diverging Gaussian, so its TRUE focus
sits at `dz_curv = f^2/(R_rx - f)` BEYOND the focal plane (S. A. Self, Appl. Opt.
22, 658 (1983), DOI 10.1364/AO.22.000658). The coupling Terms always charge that
curvature. See [defocus/README.md](defocus/README.md) and the report
[defocus/fidelity2_mmf_coupling_gap.md](defocus/fidelity2_mmf_coupling_gap.md),
whose RESOLUTION appendix records the `defocus_m` sign fix, the always-charged
curvature convention, and the new closed forms.

| File | Purpose |
| --- | --- |
| [defocus/defocus_sensing.py](defocus/defocus_sensing.py) | Pure-analytic checks of the defocus model and the bidirectional wrapper: a `dz` sweep for a multimode-fibre receiver, the two lateral-sensitivity limits, the spot radius against `gaussz` and the geometric blur, and the chief-ray lever. |
| [defocus/fidelity2_mmf_coupling_gap.md](defocus/fidelity2_mmf_coupling_gap.md) | The write-up: why the fidelity-0 MMF Term read about 7 dB more loss than the field, and how the received-curvature defocus closed most of it (about 1.2 dB left, the known 2-W1 Airy-versus-Gaussian spot-shape gap). |

## screens/

The phase-screen low-frequency (tip/tilt) study. A screen on a finite grid
holds no power below its grid fundamental; that missing band is the tip and the
tilt. The study measures how much tilt each screen route holds against the
analytic value, and it settles two open rows in `docs/schmidt-crosscheck.md`.
It is VALIDATION ONLY: it reads the production layer and changes no olb module.
See [screens/FINDINGS.md](screens/FINDINGS.md) for the write-up. The measured
tables are in `screens/data/`.

| File | Purpose |
| --- | --- |
| [screens/helpers.py](screens/helpers.py) | The analytic truths and the shared estimators. |
| [screens/oversize_crop.py](screens/oversize_crop.py) | Arm 1: the Fourier and the oversize-and-crop screens. |
| [screens/infinite_screen_stats.py](screens/infinite_screen_stats.py) | Arm 2: the spatial statistics of the extruded screens. |
| [screens/extrusion_stationarity.py](screens/extrusion_stationarity.py) | Arm 3: the drift test of the extrusion. |

## lognormal_certification/

Certify the cheap analytic aperture-averaged lognormal power draw against the
fidelity-2 split-step Monte Carlo (backlog 1-6). It measures the POINT index
`sigma2_I`, the aperture-averaging filter `A` and the deep-fade quantiles apart,
so an INDEX error, a FILTER error and a SHAPE error do not mix. QUICK-MODE
reading (2026-09-01): the lognormal FAMILY holds (the shape leg agrees inside
0.13 dB at the 5 % fade over `D/rho_0` = 0.2 to 7.9); the analytic point index is
10 to 20 % HIGH; and the FILTER is the fault, because it OVER-AVERAGES by 1.4 to
2.9 times over `D/rho_0` = 1 to 8. The absolute impact stays under 0.30 dB at the
5 % fade. The D = 40 cm collimated column is BEAM-FILLING-LIMITED
(`eta_fill = 0.87`) and must not be read as a filter error. The `--full` tail run
is still to run.
See [lognormal_certification/README.md](lognormal_certification/README.md).

| File | Purpose |
| --- | --- |
| [lognormal_certification/lognormal_certification.py](lognormal_certification/lognormal_certification.py) | The certification run: a `D/rho_0` sweep from a point-like aperture to strong averaging, for a collimated and a diverged launch. ONE propagation for each trial serves the whole aperture sweep (`propagate_turbulent_field`), so every aperture reads the same atmosphere; a matched-seed check against `propagate_turbulent_scenario` proves the two agree bit for bit. It reports the point index, the effective averaging factor, the beam-fill fraction and the absolute fade spread. It writes a results JSON, a run log and `figures/lognormal_certification.png`. `--full` raises the trial count for the 1 % fade tail. |

## tail_convergence/

The fidelity-2 single-mode-fibre fade-tail convergence study (backlog 2-I2T),
with the large-campaign measurements of backlog 2-N6. THE QUESTION: does the
deep SMF fade tail (p10, p5, p1) of a fidelity-2 downlink CONVERGE as the
near-ground `Cn2` is resolved with more, thinner phase screens? The mean is
already validated flat (`docs/schmidt-crosscheck.md`, work package 7), so the
tail is the open risk: it sets the link availability margin. The study PINS the
grid and it moves the screens only, because the shipped sizer refines the grid
with the screen count. Each case is one resumable `Campaign`.
See [tail_convergence/README.md](tail_convergence/README.md).

| File | Purpose |
| --- | --- |
| [tail_convergence/tail_convergence.py](tail_convergence/tail_convergence.py) | The study itself. Six cases: the shipped default, four pinned-grid screen counts (9, 15, 25, 40), and a near-ground refinement that splits the bottom screen into four equal-`Cn2` sub-screens. It reports p50 / p10 / p5 / p1 of the composite SMF loss with bootstrap intervals, the fade depth, the aperture and the point scintillation index, and the 2-N6 campaign numbers (wall time, seconds per trial, disk bytes, load memory, and the growth of the tail estimate with the trial count). It writes a results JSON, a run log and three figures. |

## outer_scale_tail/

The fidelity-2 SMF outer-scale fade-tail study (backlog 2-P5, 2-I3, 0-W4). It
answers two questions with one matched-seed 2 x 2 of `config` (a well-resolved
reference against the shipped `rapid` preset) by `L0` (`inf` against `25 m`).
MEASURED (2026-09-05, 30 and 20 deg, 1000 trials each): the finite outer scale
gives 2.5 dB (30 deg, 3.0 sigma) to 2.8 dB (20 deg, 2.4 sigma) LESS SMF p5 fade
than `L0 = inf`, and the point fade does not move, so the bias is the fibre TILT
and the `L0 = inf` default is that much PESSIMISTIC. At the physical `L0 = 25 m`,
`rapid` tracks the reference inside about 0.3 dB at the mean, p50, p5 and p1
(both elevations), with one small p10 wrinkle at 20 deg (+0.86 dB, safe), so
rapid is a defensible default on this scenario. Owner decision: run fidelity 2 at
a FIXED `L0 = 25 m`. See [outer_scale_tail/README.md](outer_scale_tail/README.md).

| File | Purpose |
| --- | --- |
| [outer_scale_tail/outer_scale_tail.py](outer_scale_tail/outer_scale_tail.py) | The study. Two configs (`ref`, `rapid`) crossed with the outer-scale values, a matched-seed L0 pair within each config on ONE grid and plan. It reports p50 / p10 / p5 / p1 of the composite SMF loss and the point fade with bootstrap intervals, the within-config outer-scale delta, the rapid-against-reference comparison, and a rapid-as-default verdict. It writes a results JSON and a run log for each elevation, and figures to `figures/`. |

## waveoptics_vs_fast/

The fidelity-2 field against FAST and the analytic model: the space-downlink SMF
coupling-loss gap (backlog 2-W1, 2-AO). It runs all three models UNCORRECTED
(NOAO) and like-for-like at a matched outer scale, and it reports the gap
FAST-minus-field against elevation. MEASURED (2026-09-05): the 0.7-2.9 dB gap of
the older informal comparison is an OUTER-SCALE artifact. At the physical
`L0 = 25 m` FAST and the field AGREE (gap -0.34 to +0.12 dB, 20 to 90 deg, all
within 1 sigma of zero); the gap only reopens at the grid-dependent `L0 = inf`
(+0.35 to +1.17 dB), because FAST is more outer-scale-sensitive. The analytic
term stays 1 to 2.5 dB optimistic. This certifies the uncorrected rung only.
See [waveoptics_vs_fast/README.md](waveoptics_vs_fast/README.md).

| File | Purpose |
| --- | --- |
| [waveoptics_vs_fast/waveoptics_vs_fast.py](waveoptics_vs_fast/waveoptics_vs_fast.py) | The study. FAST (`smf_fast_term`, NOAO), the field (a `Campaign` process pool), and the analytic term, per elevation, at a matched `L0`. An NPXLS convergence guard pins the FAST grid first. `--L0`, `--field-mode` (process / thread / serial), `--workers` and `--block-size` (the effective process count is `min(workers, ceil(n_trials/block_size))`). It writes a results JSON and a run log tagged by outer scale and field mode, and figures to `figures/`. |

## screen_stacking/

The phase-screen STACKING test, phase only. Does a stack of N screens hold the
statistics of one screen of the same composite `r0`? It is the generator half
of the tail-convergence study, and it answers the owner's hypothesis that a
many-screen plan loses more of the tip-tilt band than one screen. VERDICT
(2026-09-04): not at L0 = inf, where every count misses the SAME 16 to 25
percent of the aperture phase variance (all tip and tilt, `Delta3` = 1.00) —
and that deficit is the OUTER SCALE, not the generator: the grid holds scales
to 27 x its side (95 m) and the screens match a von Karman L0 = 95 m theory
exactly. At L0 = 25 m, judged against the von Karman theory, one screen reads
1.00 +-0.03 and a 5- to 25-screen plan 0.97 to 0.94 (a mild few-percent
stacking drift). The production `L0 = inf` default therefore claims an outer
scale it cannot deliver, worth a MEASURED 2.5 to 2.8 dB at p5 of the SMF fade
(30 and 20 deg; see `outer_scale_tail/`, backlog 2-P5). See
[screen_stacking/README.md](screen_stacking/README.md).

| File | Purpose |
| --- | --- |
| [screen_stacking/screen_stacking.py](screen_stacking/screen_stacking.py) | Draws 100 stacks for each of five per-screen `r0` lists (the ground layer as 1 or 4 screens; the whole plan as 5, 9 or 25) on the pinned 1024 px grid, sums them, and measures the structure function and the Noll `Delta1` / `Delta3` aperture variances as ratios to theory with standard errors. |

## vacuum_loss/

The fidelity-2 vacuum (no-turbulence) geometric loss against the analytic
geometric Term. It shows that a terrestrial far-field link agrees to about
0.15 dB, and that the full-path space solve is grid-noise-limited (the loss
scatters +/- 1 to 4 dB and does not converge at a practical grid size). That
measurement is why a space fidelity-2 budget takes the ANALYTIC geometric Term
by default.

| File | Purpose |
| --- | --- |
| [vacuum_loss/vacuum_loss_validation.py](vacuum_loss/vacuum_loss_validation.py) | The cross-check itself. It writes `vacuum_loss_results.json`. |

## waveoptics_speed/

The fidelity-2 speed campaign (P0 to P4; see
`docs/waveoptics-efficiency-plan.md`). Each is a measurement script with one
results JSON and one run log; none touches production code.

| File | Purpose |
| --- | --- |
| [waveoptics_speed/profile_baseline.py](waveoptics_speed/profile_baseline.py) | P0: where one turbulent trial spends its time. It gates the rest of the plan (screen generation is 80 to 84% of a trial). |
| [waveoptics_speed/screen_generator_check.py](waveoptics_speed/screen_generator_check.py) | P1: the fast, cached olb `ScreenFactory` against the aotools baseline (structure function, speed, accuracy). |
| [waveoptics_speed/generator_validation.py](waveoptics_speed/generator_validation.py) | A broad validity pass on the olb generator across geometries, presets, the outer scale, and the FADE TAIL. Verdict: a trustworthy drop-in. |
| [waveoptics_speed/coarse_screen_experiment.py](waveoptics_speed/coarse_screen_experiment.py) | P2 experiment (a): coarse screens plus interpolation. BURIED (loses the Fresnel-scale phase that builds scintillation). |
| [waveoptics_speed/beam_grid_experiment.py](waveoptics_speed/beam_grid_experiment.py) | P2 experiment (b): a grid that follows the beam. BURIED for the wired scenarios (the flat grid already wins). |
| [waveoptics_speed/scaling_study.py](waveoptics_speed/scaling_study.py) | P3: how trials scale across workers (threads, processes, batched split step). Processes beat threads; threads saturate at 8 to 16 workers. |
| [waveoptics_speed/make_plots.py](waveoptics_speed/make_plots.py) | Draw one PNG per speed task from its results JSON, into `figures/`. Skips a task whose JSON is absent. |

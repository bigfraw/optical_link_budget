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
| [screens/n_columns_sweep.py](screens/n_columns_sweep.py) | Arm 4 (2026-09-08): what `n_columns` buys on the fixed float64 aotools kernel. It sweeps `n_columns` 2 to 32 at L0 2.56 and 25 m, N 128 and 512, 8 seeds, plus a record-length study. See FINDINGS Q6. |
| [screens/subharmonic_start.py](screens/subharmonic_start.py) | Arm 5 (2026-09-08): a SUBHARMONIC initial frame against the stock plain start. It removes the optical spin-up at row 0 (the piston-removed variance 0.90 to 1.03 and the 1 m Z-tilt 1.04 to 1.08 of theory, against 0.32 and 0.71 for the plain start). See FINDINGS Q7. |
| [screens/kept_start.py](screens/kept_start.py) | Arm 6 (2026-09-08, 32 seeds, 2 SE): is the subharmonic start KEPT on a 0.2 L0 frame? No, and it is not carried either: the plain and the subharmonic start read the SAME extrusion-axis excess (+0.19 to +0.27 at 0.5 L0, 2 SE 0.08, the Q6 value), so the recursion imposes its own axis statistics whatever the first frame held. The subharmonic start buys the first frame only. See FINDINGS Q8. |
| [screens/PLAN_n_columns.md](screens/PLAN_n_columns.md) | The plan of arm 4, the aotools code review (PR 111, issues 107 and 109), and the parked Phase 4 (the time-axis speed evaluation, NOT started). |

**The extrusion in the LEO point-ahead regime (2026-09-08).** The regime of
record: a point-ahead angle of at most 10 arcsec, so the uplink and the
downlink beams sit 0.97 m apart at 20 km; a ground aperture of at most 1 m;
a realistic outer scale of 20 to 50 m; the usual sampling constraints, so a
grid of more than 512 px with a side of 4 to 7 m. That is a frame side of
0.1 to 0.35 outer scales, the cell of the sweep (N 512, L0 25 m, side 0.2 L0)
where NO `n_columns` passes the extrusion-axis band. The reading splits by
scale:

- The PAA anisoplanatism is served. D(r) is isotropic inside 5 percent from
  5 cm to 1 m at `n_columns = 2`. Eight columns do NOT help there: they clean
  the 2.5 m end and roughen the 5 cm end by 15 percent, so keep 2 columns for
  the aperture scale.
- The absolute level is the real error, and it is a spin-up and frame-width
  effect, not `n_columns`: the row-lag variance reads 0.91 of the theory and
  the Z-tilt variance over 1 m reads 0.84 to 0.88, so a pointing-fade depth
  from this screen is about 0.6 dB optimistic at every `n_columns`.
- The extrusion-axis over-correlation (+0.1 to +0.23 in rho) lives at lags of
  0.1 to 1 L0, so 2 to 50 m of translation. A snapshot study never sees it.
  Under frozen flow it is the tilt and piston spectrum below v / (0.1 L0):
  below about 100 Hz for the top layer at 250 m/s of apparent slew, below
  about 2 Hz for the ground layer, so a temporal tilt model reads too slow in
  the tracking band and its fade durations come out too benign.

Not measured yet: the absolute D(1 m) on both axes at L0 20 and 50 m (the
sweep stored correlation coefficients for the 25 m cell). That one read closes
the PAA question with a number.

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

## fibre_fade_models/

The terrestrial power-distribution study (backlog 1-9, steps 1 to 3 DONE
2026-09-08; the certification of record is physics.md Section 9m). It reads
the twelve 2-TC terrestrial campaigns post hoc and fits every family of
`andrews/distributions.py` to the bucket power, the fibre-coupled power and
the point irradiance, with the tilt split out. VERDICTS: the fidelity-0
analytic lognormal holds a bucket to `sigma_R^2 = 0.7` and a refit
gamma-gamma above it; NO free route holds a 10 cm fibre (the shipped chain
under-reads the p5 fade by 1.3 to 11 dB, the tilt is the tail), and its
fitted family of record is the lognormal-Rician; the received tilt is the
aperture angle of arrival at the Gaussian r0 reduced by the outer scale, and
the walk-off Term reads half its variance.
See [fibre_fade_models/README.md](fibre_fade_models/README.md).

| File | Purpose |
| --- | --- |
| [fibre_fade_models/extract_trials.py](fibre_fade_models/extract_trials.py) | The one-pass read of the stored campaigns on bigfraw (`Campaign.map_trials`): the bucket power, the SMF coupling three ways (untracked, tracked focus, tilt removed), the Noll tilt pair, the slope aliasing test, the point irradiance and an MMF coupling, for each trial. |
| [fibre_fade_models/fit_distributions.py](fibre_fade_models/fit_distributions.py) | The fits (moment, maximum likelihood, and the free analytic routes), the verdict tables, `fit_results.json` and the three exceedance figures. |
| [fibre_fade_models/r_rule.py](fibre_fade_models/r_rule.py) | The free-fibre rule search. With `sigma_z^2` pinned to the bucket index and the lognormal-Rician `r` fitted alone, the TILT-REMOVED fibre follows `r = 2.3 / sigma2_HO` (leave-one-cell-out 12 of 12, one launch); the untracked fibre has no rule in `r` alone. It motivates the PROPOSED composite free route (the higher-order lognormal-Rician times the walk-off fade). It writes `r_rule.log`. |

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

## terrestrial_campaigns/

The terrestrial fidelity-2 backbone dataset. It stores twelve resumable
`Campaign` sets: three path lengths (2, 5 and 10 km) crossed with two
turbulence strengths (`Cn2` = 3e-15 and 1e-14), each one at the `rapid` and at
the `standard` preset, at the fixed `L0 = 25 m`. It is a RUNNER only: it stores
the trials and it reports the cost. Every analysis is DEFERRED (the fade
distribution, the index split, the fade quantiles, the rapid-against-standard
verdict, the comparison with the analytic terrestrial Terms). The dataset
serves backlog 1-8 gate (b), the rapid-preset question, backlog 2-N2 beam
filling, and the single-mode-fibre coupling distribution of a horizontal link.
STATUS: the run is DONE (2026-09-06). Twelve campaigns of 2000 trials, about
0.9 GB, stay on bigfraw; the four `standard` cells past 2 km ran with the two
speed opt-ins, under the `_scipy_lean` roots. The run logs, the `cell.json`
records and the block census are committed under
`terrestrial_campaigns/records/`.
See [terrestrial_campaigns/README.md](terrestrial_campaigns/README.md).

| File | Purpose |
| --- | --- |
| [terrestrial_campaigns/run_campaigns.py](terrestrial_campaigns/run_campaigns.py) | The runner. `--dry-run` sizes every cell and prints the grid, the screen count and the memory of each one; `--smoke` runs a few trials for each cell and reports the seconds for each trial, the projected hours, the peak working set and a post-hoc `recouple` cross-check; the plain call stores the trials. Each campaign root gets a `cell.json` with the Rytov variance, `rho_0`, `w(L)`, the beam-fill fraction and the curvature focus shift of each aperture, the grid, the screen plan and every sizer warning. |
| [terrestrial_campaigns/waist_bias_check.py](terrestrial_campaigns/waist_bias_check.py) | The grid-bias check. It propagates the launch field in VACUUM on each cell grid, with no screens, and it reports the pixels across the 5 mm waist, the second-moment beam radius against `olb.beam.gaussz`, and the 10 cm bucket power against the analytic Gaussian capture, with the difference in dB. It needs no stored trials. It writes `waist_bias_check.log` and `waist_bias_check_results.json`. |
| [terrestrial_campaigns/optin_crosscheck.py](terrestrial_campaigns/optin_crosscheck.py) | The opt-in cross-check. It opens the default-settings `L5km_cn23e-15_standard` campaign and the opt-in `L5km_cn23e-15_standard_scipy_lean` campaign, and it compares them trial by trial: the maximum and the median relative difference of the collected power and of `smf_eta`, plus the mean loss and the p5 and p1 fades of both. It needs BOTH stores, so it runs on bigfraw; `--dry-run` prints the two roots. It writes `optin_crosscheck.log` and `optin_crosscheck_results.json`. |

## terrestrial_screen_count/

The terrestrial SCREEN-COUNT convergence sweep (backlog 2-TC1). At the Schmidt
per-screen cap the 10 km / `Cn2` = 1e-14 standard cell asks for 35 screens, and
that cap is a thin-screen VALIDITY rule (Schmidt, DOI 10.1117/3.866274, Listing
9.5, printed p. 175), not a convergence result. The study holds that cell's GRID
fixed and it OVERRIDES the screen count with a caller plan at n = 5, 10 and 15
(the owner stopped it before the planned 20), against the 35-screen backbone
campaign as the reference (it is reopened, never rerun). It compares the collected power, the centre-pixel irradiance, the
fibre coupling and a 5 cm bucket: the index and the p10 / p5 / p1 fades, each
with a bootstrap interval. THE 10 km CELL IS DONE (2026-09-06): the counts 5,
10 and 15 ran at 1000 trials each, and the owner stopped the sweep before the
20-screen count. The verdict is BY RECEIVER KIND. An SMF receiver shows NO
detectable screen-count effect: the p5 and p1 deltas carry no sign pattern and
they sit inside the noise, because the fibre fade in saturation is tilt and
low-order phase, which few screens already carry. A BUCKET receiver shows a
small, consistently OPTIMISTIC bias: the 10 cm bucket p5 reads 0.46 to 0.62 dB
less fade, with the same sign at every count. THE 5 km CELL IS LAUNCHED
(2026-09-06 on bigfraw, 12 workers, the counts 5 / 10 / 15 / 20 at 2000 trials
each). Its results are PENDING. That cell sits at the `min_screens` floor, so it
asks the other half of the question: is the floor enough?
The worker plateau after the memory cut is measured here too.
See [terrestrial_screen_count/README.md](terrestrial_screen_count/README.md).

| File | Purpose |
| --- | --- |
| [terrestrial_screen_count/screen_count_sweep.py](terrestrial_screen_count/screen_count_sweep.py) | The sweep. It takes ANY cell (`--path-km`, `--cn2`, `--counts`, `--n-trials`, `--tolerance-db`), and the file names carry the cell tag. `--dry-run` sizes every count and prints the grid, the per-screen `sigma2_r` maximum and the projected cost. The plain call stores the trials cheapest first (resumable), then analyses; `--run-only` stores the trials and skips the analysis (so a staged run does not pay the analysis every time); `--analyse-only` reads what is stored. It first runs a REFERENCE-COUNT check: an override plan at the reference count must give bit-identical trials, which proves the caller plan is the only change. It writes a results JSON, a run log and a figure, and it prints ONE table with the reference row first, a delta against the reference for each quantity, and TWO verdicts for each count: CONVERGED / NOT CONVERGED against the bootstrap noise, and a pass or fail against a dB tolerance. |

## receiver_cone_clip/

The CONVERGING absorbing boundary (backlog 2-P3, route (b)). Only the light
inside the back-projected cone of the receive aperture can land on it, so a
mask that follows that cone should hold the aperture field and remove most of
the grid extent. The test holds the GRID, the screen plan and the screens
fixed, and it changes ONLY the mask. VERDICT (2026-09-06, TWO rules and two
cells): NO. In VACUUM the ray-geometry claim is right at both path lengths: the
cone removes 93 to 99.6 percent of the launched power and the aperture field is
untouched from `c = 1.5` up. Under TURBULENCE neither rule holds the field. The
MAPPED rule (the ray map `z/L` on the scatter term too) reaches a field RMS of
1.2e-2 at `c = 5`, against a 1e-3 limit. The UNMAPPED rule (the ray map on the
aperture term only, the scatter allowance un-mapped and read from the whole
path) is ten times better, 1.24e-3 at `c = 5` at 5 km and 5.6e-3 at 10 km, and
it still fails at every factor. The residue is the power-law TAIL of the
scattered light, not the cone shape: the field RMS falls as about `c^-4.2`
(5 km) and `c^-2.1` (10 km), so the field rule needs about `c = 5.3` and
`c = 11.5`, whose implied sides are 1.24 x and 3.25 x LARGER than the sizer
side. The budget-visible metrics alone (0.01 dB and 1e-3 of efficiency, all a
Term reads) pass at `c = 2` at 5 km (a 1.71 x smaller side, but the same 1024
px) and at `c = 5` at 10 km (a side 1.5 x LARGER than the sizer). So the saving
is not real at 10 km, the cell that motivated the route. See
[receiver_cone_clip/README.md](receiver_cone_clip/README.md).

| File | Purpose |
| --- | --- |
| [receiver_cone_clip/receiver_cone_clip.py](receiver_cone_clip/receiver_cone_clip.py) | The test. It rebuilds the production `split_step` loop by hand and asserts `numpy.array_equal` against `propagate_turbulent_field`, then reruns the same screens with the converging mask. `--cone-rule` (`mapped` is the first form, `unmapped` the default), `--guard` (`on` keeps the `3 w(z)` near-field guard, `off` is the default because the guard covers the whole grid at 5 km), `--trials`, `--factors`, `--preset`, `--seed`, `--path-km`, `--cn2`. Every output file carries the rule, the guard and the cell, and the run prints the c-versus-metric table for turbulence and for vacuum, the smallest passing factor under the FIELD rule and under the BUDGET-VISIBLE rule, and the grid side and pixel count that each implies. |

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
| [waveoptics_speed/fair_scaling_rerun.py](waveoptics_speed/fair_scaling_rerun.py) | The FAIR rerun of P3 (2026-09-04): it pins the BLAS thread pool before numpy imports and it measures a warm process pool in steady state. VERDICT: threads and processes TIE on wall time for ONE run (the Windows pool spawn costs 2.5 to 4.4 s); processes win 1.15x to 1.7x only when the pool stays warm across many blocks, which is `Campaign`. The 8-to-16-worker plateau is the machine (memory bandwidth, the hybrid P/E cores), not the GIL. |
| [waveoptics_speed/make_plots.py](waveoptics_speed/make_plots.py) | Draw one PNG per speed task from its results JSON, into `figures/`. Skips a task whose JSON is absent. |

## memory_cut/

The MEMORY cut of one fidelity-2 trial, and the two speed OPT-INS (2026-09-06).
It has two halves. The first half checks the bit-identical cut: the lazy screen
stack (`split_step` reads its screens one at a time), the `Forvard` transfer-
function cache, and the pool sizer of `olb/waveoptics/resources.py`
(`Campaign.run(workers="auto")`). The second half measures the two OPT-INS that
are NOT defaults, because each one moves every seeded number at the rounding
level: `screen_generator="olb-lean"` (the same physics and the same random
stream through fewer full-grid passes) and `fft_backend="scipy"`. The
measurements are from a 4-core container, so the RATIOS carry over and the
seconds do not. The worker plateau that goes with them was measured later on
bigfraw; it is in `terrestrial_screen_count/`.
See [memory_cut/README.md](memory_cut/README.md).

| File | Purpose |
| --- | --- |
| [memory_cut/memory_cut_check.py](memory_cut/memory_cut_check.py) | The bit-identical half. It proves the lazy screens, the `Forvard` cache and the pool sizer, and it times one serial trial with the cache off and on (single 9 screens 3.61 to 2.58 s; double 15 screens 8.09 to 6.64 s). It writes `memory_cut_check.json`. |
| [memory_cut/screen_generator_lean.py](memory_cut/screen_generator_lean.py) | The opt-in half. The lean generator against the default one (the same draw to 1e-7 in float32 and 2e-16 to 4e-16 in float64, the fitted `r0` inside the standard error), the cost of one screen, the raw transform time of the two FFT backends, and one trial in each of the four combinations (1.38x together). It writes `screen_generator_lean.json`. |

## fast_stone_pointahead/

The FAST against Stone point-ahead anisoplanatism study (backlog 1-5). An
uplink terminal senses turbulence on a downlink beacon and applies the
conjugate phase; the satellite moves during the round trip, so the correction
decorrelates over the point-ahead angle and a residual phase variance stays.
olb holds TWO models of that residual — fidelity 1 (FAST, `uplink_fast_term`;
Farley et al., DOI 10.1364/OE.458659) and fidelity 0 (the Stone 1994 modal law,
`uplink_point_ahead_term`; DOI 10.1364/JOSAA.11.000347) — and this study
compares them at matched conditions. VERDICT (2026-09-02): MATCH, both routes
validated (physics.md Section 9j). With the servo off the PAOLA filter reduces
exactly to `2 - 2cos(delta_r . kappa)` (to 9e-16), and at MATCHED mode sets the
two agree to about 5 % across the full sweep; the Term-level factor is the mode
set (the FAST mask keeps piston and tilts, so its analytic partner is Stone
`remove='none'`).
See [fast_stone_pointahead/README.md](fast_stone_pointahead/README.md).

| File | Purpose |
| --- | --- |
| [fast_stone_pointahead/fast_stone_pointahead.py](fast_stone_pointahead/fast_stone_pointahead.py) | The comparison. It drives the FAST PAOLA aniso-servo residual and the Stone modal decorrelation residual on the same `Cn2` profile, aperture and corrected order, with the servo and the sensor effects off, then it sweeps the point-ahead angle and the corrected order. It also compares the fitting sides (FAST `sim.fitting_error` against the Noll residual) and the full Terms, and it decomposes the Term-level gap into the mode set, the FAST auto-grid truncation and the Marechal-versus-Monte-Carlo mapping. |

## waveoptics_ao/

The perfect-AO validation of the fidelity-2 layer (backlog 2-AO). The layer
removes the first N Noll modes of the wavefront over the receive aperture
(`olb/waveoptics/compensation/`); this study measures that chain on the hero
0.7 m SMF downlink (V0 to V5, on the cupy backend, 2026-09-07, L0 = 25 m;
physics.md Section 9l). VERDICTS: the modal chain matches the Noll residual law
inside 2 % from J = 3 up (J = 1 reads the outer scale); the summed-screen source
is TRUSTED to 20 deg; the post-hoc route equals the in-run route to 1.4e-07; the
SMF p5 fade improves by 9.4 / 22.1 / 24.3 dB at 30 deg for TipTilt / AO(10) /
AO(21), and the bucket power does not move; the corrected field and the tracked
fidelity-1 FAST Term agree on the mean to -0.5 to +0.1 dB. THE SOURCE RULE
(V5): screens on a space link, slopes on a terrestrial one — a summed screen
weights every plane equally while the arriving terrestrial tilt carries the
`(1 - z/L)` path lever, so a screen-sensed terrestrial AO is wrong.
See [waveoptics_ao/README.md](waveoptics_ao/README.md).

| File | Purpose |
| --- | --- |
| [waveoptics_ao/waveoptics_ao.py](waveoptics_ao/waveoptics_ao.py) | The study. It runs (or reopens) the AO campaigns and answers the five questions: is the summed screen phase a valid space sensing source (V0); does the modal chain obey Noll's residual law (V2); does the correction move the SMF fade in the safe direction and how far is it from the tracked FAST Term (V1); does the wrapped-gradient slope source agree with the summed-screen source (V3), survive a terrestrial link (V4), and which source is right on a terrestrial path and by how much (V5). |

## anisoplanatism_screens/

The screen plan against the point-ahead anisoplanatism (backlog 2-P4, the
gates). Two PHASE-ONLY studies, on a laptop CPU, with no campaign and no
propagation. Study 1 asks whether the production screen plan reproduces the
continuous Stone variance (Stone and others, DOI 10.1364/JOSAA.11.000347), and
whether a plan cut by the anisoplanatic weight does better. Study 2 asks whether
one screen drawn once and read through two integer-pixel windows gives that
variance back. THE CONVENTION: the TILT STAYS IN (`remove="piston"`, owner
decision 2026-09-11), because the terminal senses the downlink beacon tilt and
the steering mirror adds the point-ahead offset geometrically. VERDICTS: NO
PLANNER CHANGE — the production nine-screen plan reads 0.987 to 0.995 of the
continuous Stone sum over the 135 owner-regime cells at `L0 = inf` AND at
`L0 = 25 m`, and the equal-anisoplanatic-weight families are worse; THE WINDOW
RULE HOLDS at both outer scales — 52 of 54 cells sit inside the 0.95 to 1.05
band or inside their own 2-sigma bar; THE OUTER-SCALE KERNEL IS VALIDATED — the
measured variance and the von Karman Stone reference both fall by the SAME 1.7
percent between the two scales. The tilt multiplies the variance by 2.4, and the
25 m outer scale removes 1.6 percent. See
[anisoplanatism_screens/README.md](anisoplanatism_screens/README.md).

| File | Purpose |
| --- | --- |
| [anisoplanatism_screens/common.py](anisoplanatism_screens/common.py) | The shared case. It holds the hero uplink scenario, the plan builders (the production plan, the equal-Rytov override, the ground split, and the two equal-anisoplanatic-weight families) and the log helpers. |
| [anisoplanatism_screens/anisoplanatism_screens.py](anisoplanatism_screens/anisoplanatism_screens.py) | Study 1. It evaluates the delta-layer Stone sum on each screen plan and compares it with the continuous integral, over elevations, apertures, angles and corrected orders. It takes `--remove` and `--L0`, and it writes the log, the results JSON, `production_table[_L0<m>].md` and `figures/discrete_vs_continuous[_L0<m>].png`. |
| [anisoplanatism_screens/phase_only_shift.py](anisoplanatism_screens/phase_only_shift.py) | Study 2. It draws each screen one time on an oversize grid with the production `ScreenFactory`, reads the beacon window and the shifted uplink window, sums the difference over the plan, and fits the Noll bands over many apertures. It takes `--L0`, and it writes the log, the results JSON and `figures/phase_only_shift[_L0<m>].png`. |

## waveoptics_pointahead/

The point-ahead validation of the fidelity-2 uplink (backlog 2-P4). The runner
makes one more propagation pass for each point-ahead angle, through a laterally
shifted window of the SAME screens, so the drop of the reciprocity overlap
between the beacon direction and the point-ahead direction IS the point-ahead
anisoplanatism. Nine campaigns of the hero 0.7 m uplink, 4400 trials, on the
cupy backend, `L0 = 25 m`, 2026-09-11; physics.md Section 9n. VERDICTS: p0, the
record holds together — the stored zero angle equals the beacon overlap bit for
bit, the post-hoc read of the stored planes matches an in-run campaign to
9.0e-07, and the regeneration route is BIT IDENTICAL on the campaign backend;
p1, at 30 deg the field and the fidelity-1 FAST Term agree inside 0.5 dB at
every AO cell (the FAST run-to-run spread is 0.2 to 0.3 dB at 1000 draws), and
at 20 deg FAST reads 0.7 to 1.5 dB ABOVE the field, which is EXPECTED, because
FAST propagates no field and holds no saturation; the Stone column overstates
everywhere (the extended Marechal saturates past 1 rad^2) and it is a REPORT,
not a gate; p2, the p5 fade penalty at the geometry angle is 8 to 9 dB at 30 deg
and 11 to 12 dB at 20 deg for AO(10)/AO(21), against 2.7 to 3.9 dB on the mean
(tip-tilt 2.4 dB, uncorrected zero), so a link sized on the MEAN is under-sized;
p3, 24 of 24 screen-count readings are flat from 5 screens up. THE MODEL OF
RECORD STAYS FIDELITY 1 (FAST); whether the fidelity-2 Term takes that role is
an open owner decision. The run takes under 30 minutes of GPU wall time and it
writes 5 GB. See
[waveoptics_pointahead/README.md](waveoptics_pointahead/README.md).

| File | Purpose |
| --- | --- |
| [waveoptics_pointahead/waveoptics_pointahead.py](waveoptics_pointahead/waveoptics_pointahead.py) | The driver. `--study p0 p1 p2 p3` runs (or reopens) the point-ahead campaigns and answers the four questions: do the record identities hold (p0); how large is the penalty against FAST and against Stone (p1); what does the point ahead do to the fade, not only to the mean (p2); is the penalty converged in the screen count (p3). It derives the three corrected stacks from the uncorrected campaign through `Campaign.recouple_point_ahead`, and it writes a log, a results JSON and two figures. |

## gtilt_sensing/

The G-tilt against wrapped-slopes tilt-sensing study. For a terrestrial
fidelity-2 TipTilt correction, does the far-field intensity centroid (the
G-tilt, the intensity-weighted mean phase gradient; Tyler 1994,
DOI 10.1364/JOSAA.11.000358) sense the tilt better than the shipped
wrapped-gradient slopes, which alias where the local phase step passes pi? It
reopens each stored terrestrial 2-TC campaign and runs the correction post hoc
through both routes. VERDICT (2026-09-08): the wrapped slopes WIN; the G-tilt
is kept as an opt-in source only, not a default. See
[gtilt_sensing/README.md](gtilt_sensing/README.md).

| File | Purpose |
| --- | --- |
| [gtilt_sensing/gtilt_vs_slopes.py](gtilt_sensing/gtilt_vs_slopes.py) | The comparison. It reopens the terrestrial campaigns and runs three post-hoc coupling passes (uncorrected, slope-sensed TipTilt, G-tilt-sensed TipTilt) on the SAME stored fields, so the only change is the sensing source, and it reports the SMF fade each way. |

## precision/

The single-against-double precision check of the fidelity-2 layer. It runs the
SAME turbulent trials two times — once in double precision (complex128,
float64 screens, the pre-2026-09-05 record) and once in single precision
(complex64, float32 screens, the default since 2026-09-05) — on one seed, so
trial k sees the same atmosphere, and it measures the difference. Single
precision exists because a `Campaign` is memory-bandwidth bound, so half the
bytes for each element is the one change with real leverage. MEASURED: the two
agree to parts per million. See [precision/README.md](precision/README.md).

| File | Purpose |
| --- | --- |
| [precision/precision_check.py](precision/precision_check.py) | The check. It runs matched-seed double and single trials of a downlink and reports the maximum relative difference of the per-trial scalars. It boosts its own process priority (it is a direct runner call over ssh; see the priority note). |

## posthoc_speed/

The post-hoc read-back of a stored campaign on the CROP (2026-09-07). A stored
trial holds the patch DISC pixels only; the old read scattered them into the
FULL grid and swept the zero padding (14.3x the useful pixels at 1024 px) and
rebuilt the masks and the fibre mode for every trial. The read now works on the
square CROP that just holds the disc, under ONE rule: PUPIL-plane quantities on
the crop, FOCAL-plane quantities on the padded grid (an MMF or a Camera
focuses). MEASURED: the pupil routes agree to 3e-15 (the padded MMF route is
bit-identical) and one trial is 2.4x to 12.5x faster; a process pool loses on a
light read (the Windows spawn) and wins 1.28x on a compensated read.
See [posthoc_speed/README.md](posthoc_speed/README.md).

| File | Purpose |
| --- | --- |
| [posthoc_speed/posthoc_speed.py](posthoc_speed/posthoc_speed.py) | The benchmark. It reads a stored campaign both ways (the crop against the full-grid comparison route, `compact=True` against `compact=False`) and reports the agreement of each read-back quantity and the per-trial speed-up, for the light reads and for a compensated read, serial and pooled. It writes `posthoc_speed.json`. |

## gpu_fft/

The CUDA (cupy) FFT backend (backlog 2-N8), an explicit opt-in built in three
milestones (2026-09-07). A fidelity-2 trial is FFTs, and the 12-worker CPU pool
on bigfraw is memory-bandwidth bound, so a device with several times the host
memory bandwidth is the one lever left on the plateau. The white noise stays a
numpy PCG64 host draw, so the same seed gives the same atmosphere, and a device
trial agrees with a host trial at the float32 rounding level. RESULT: one GPU
stream beats the 12-worker CPU pool by about 4x at 1024 and 2048 px. It needs
the `gpu` extra and the `olb-gpu-venv`. See [gpu_fft/README.md](gpu_fft/README.md).

| File | Purpose |
| --- | --- |
| [gpu_fft/gpu_microbench.py](gpu_fft/gpu_microbench.py) | The pre-code look: raw `fft2` at 512 to 4096 px (numpy, scipy, cupy), one real CPU trial of the hero downlink under the numpy and the scipy backends (a Forvard hook counts the plan's transforms), and a SYNTHETIC device trial that replays the hops and the screen transforms on the device. It gates the design before any package change. |
| [gpu_fft/cupy_backend_check.py](gpu_fft/cupy_backend_check.py) | Milestone 1, the kernel certification: does `ScreenFactory.make` give the SAME screen on the device and the host for one seed (about 1e-7 rms), does `split_step` give the same receive field for the same start field and screens (about 1e-6 rms), and the speed of one 1024 px split step under the three backends. It writes `cupy_backend_check.json`. |
| [gpu_fft/cupy_campaign_check.py](gpu_fft/cupy_campaign_check.py) | Milestone 2, the backend end to end through the ONE knob (`fft_backend="cupy"` on the runner, `Campaign`, `run_waveoptics` and `run_fidelity2`): the trials agree with the host to the float32 level (1024 and 2048 px, runner and campaign), and the speed. It writes `cupy_campaign_check.json`. |
| [gpu_fft/cupy_trial_profile.py](gpu_fft/cupy_trial_profile.py) | Milestone 3, the profile and the fixes: it wraps the stages of the REAL trial loop with timers (so the measured code is the shipped code) to find the host cost that milestone 2 left, and it records the before and the after of moving the Forvard factors, the tail and the setup onto the device (18.7x on one serial 1024 px trial). It writes `cupy_trial_profile.json`. |

## campaign_resources/

The large-campaign resource monitor. `campaign_resources.py` runs one
fidelity-2 space downlink through the production `Campaign` store and records
the CPU and the RAM of the machine while the trials run, to VERIFY that the
process pool does the parallel work. The case: a 700 mm ground aperture with a
30 % central obscuration, an SMF detector, a 500 km orbit at 30 deg, the
standard preset, L0 = 25 m, seed 20260905; the default is 4000 trials in blocks
of 50 on 16 workers. Its README also carries the ssh-launch rules for a desktop
run. See [campaign_resources/README.md](campaign_resources/README.md).

| File | Purpose |
| --- | --- |
| [campaign_resources/campaign_resources.py](campaign_resources/campaign_resources.py) | The monitored run. `--workers <N>` / `--workers auto` / `--threads` select the parallelism, `--precision single` the element size. It samples the per-core utility and the working set over the run and reports whether the pool saturates the machine, alongside the trials-per-second throughput. |

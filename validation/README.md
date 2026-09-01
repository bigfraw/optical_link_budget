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
| [waveoptics_speed/cache_check.py](waveoptics_speed/cache_check.py) | P4: the opt-in disk cache. A hit returns in about a millisecond; a grow computes only the new blocks. |
| [waveoptics_speed/make_plots.py](waveoptics_speed/make_plots.py) | Draw one PNG per speed task from its results JSON, into `figures/`. Skips a task whose JSON is absent. |

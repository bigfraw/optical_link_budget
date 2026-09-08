# The terrestrial power distributions, bucket and fibre (backlog 1-9)

STATUS (2026-09-08): steps 1 to 3 of backlog 1-9 are DONE. The twelve
terrestrial fidelity-2 campaigns of `validation/terrestrial_campaigns/`
(2 / 5 / 10 km, `Cn2` 3e-15 / 1e-14, `rapid` / `standard`, 2000 trials each,
`L0 = 25 m`, single precision) are read post hoc, every candidate family is
fitted to the BUCKET power, to the FIBRE-coupled power and to the POINT
irradiance, and each receiver case carries a verdict. Step 4 (the wiring of
`terrestrial_budget(fidelity=1)`) WAITS for an owner decision; Section 6
says what the decision is. Step 5 is the entry `docs/physics.md` Section 9m.

Two scripts:

- `extract_trials.py` runs on the machine that holds the campaigns (bigfraw).
  It reads every stored receive field ONE time through `Campaign.map_trials`
  and it writes one small table for each cell under `trials/` (about 20
  scalars for each trial: the 10 cm and 5 cm bucket power, the SMF coupling
  three ways, the Noll tilt pair, the slope aliasing test, the point
  irradiance, an MMF light-bucket coupling). The post-hoc read reproduces the
  in-run `collected_power` and `smf_eta` to 3e-7 on every cell. The twelve
  tables (3.3 MB) sit on bigfraw and on the laptop; they are git-ignored.
- `fit_distributions.py` runs anywhere. It fits the families, it prints the
  verdict tables to `fit_distributions.log`, it writes `fit_results.json`
  and the three figures under `figures/`. About 10 minutes.

## 1. The question

Backlog 1-9 asks for a `terrestrial_budget(fidelity=1)` that gives a real
fade for BOTH receiver kinds from a distribution that fidelity 2 has
certified, at the cost of a draw and not a field solve. Three sub-questions:

1. Which family holds for the BUCKET power, and does the fidelity-0 analytic
   lognormal (`terrestrial_scintillation_term`, `sigma2_P = A sigma2_I`) hold
   with its OWN parameters past the one weak path of 1-6?
2. Which family holds for the FIBRE-coupled power, and does any FREE route
   (the shipped fidelity-0 chain, or the lognormal-Rician of the plan) hold?
3. How much of the fibre fade is the TILT, and does the tilt that the walk-off
   Term reads match the field?

## 2. What is measured

Every power is normalised by its own sample mean, so a fade is the loss
below the mean, positive dB, and the fade at the exceedance q is the loss the
link exceeds a fraction q of the time. The 10, 5 and 1 percent fades are
reported with a 300-sample bootstrap 16 to 84 percent band.

The families, from `olb/turbulence/andrews/distributions.py` (Andrews and
Phillips, 2nd ed. (2005), DOI 10.1117/3.626196):

| family | parameters | how they are set |
| --- | --- | --- |
| lognormal (analytic) | `sigma2_P` | the fidelity-0 Term (FREE) |
| lognormal (refit index) | `sigma2_P` | the measured normalised variance (the 1-8 route) |
| lognormal (MLE) | mean and variance of `ln P` | maximum likelihood |
| gamma-gamma (MLE) | alpha, beta | maximum likelihood on the PDF |
| K (MLE) | alpha | maximum likelihood |
| negative exponential | none | the speckle limit (FREE) |
| lognormal-Rician (MLE) | r, `sigma_z^2` | a binned maximum likelihood on the Ch. 9 Eq. (133) PDF |
| lognormal-Rician (Strehl map, olb heuristic) | r, `sigma_z^2` | `r = S/(1-S)`, `S = exp(-sigma_phi^2)` from the Noll residual (extended Marechal, T. S. Ross, DOI 10.1364/AO.48.001812), `sigma_z^2 = ln(1 + sigma2_P)` (FREE). This map is olb's OWN heuristic from the 2026-09-05 plan, NOT a result of the book: Andrews and Phillips state at printed p. 369 that no map from atmospheric conditions to r and `sigma_z^2` is known. It identifies the Marechal core with the Rician coherent phasor and the halo with the Gaussian random part. |
| shipped F0 chain | none | the mean-only higher-order coupling plus the exponential walk-off fade of `olb/models/coupling/terrestrial.py` (FREE) |
| composite | none | the analytic lognormal bucket times that walk-off fade (FREE) |

THE VERDICT RULE. A family HOLDS on a case when it agrees with the empirical
fade inside 0.5 dB at 5 percent AND inside 1.0 dB at 1 percent. The 0.5 dB is
the "notable" rule of 1-6. A case reads FREE when a family with analytic
parameters holds, REFIT when only a fitted family holds, and NONE otherwise.

THE TILT SPLIT. The tilt of each trial is sensed from the wrapped-gradient
slopes of the clipped field (`olb.waveoptics.compensation`, the terrestrial
sensing rule of the runner), fitted with 21 Noll modes, and the tilt pair is
removed from the field before a second coupling. A cell whose phase step per
pixel passes 2.8 rad in more than 10 percent of its trials is ALIASED: its
tilt-removed cases and its tilt table are not results. That is every 5 km /
1e-14 and every 10 km cell (34 to 97 percent of the trials). The bucket, the
untracked fibre and the MMF do not read the slope sensor, so they stand there.
TWO CAVEATS FOLLOW (owner reading, 2026-09-08). First, an aliased
"tilt removed" case has NOT had all its tilt removed (the sensor reads the
tilt too small), so its fade is the untracked fade minus a part of the tilt,
not the higher-order fade: do not trust it. Second, the same steep phase that
defeats the sensor means the FIELD of those cells is under-sampled on the
clamped grid (the sizer wanted 8192 px and got 2048), and the untracked
coupled loss there leans OPTIMISTIC: the `standard` preset reads 0.2 dB
(10 km / 3e-15) and 0.7 dB (10 km / 1e-14) more fade at p5 than `rapid`, and
the 2-TC1 screen sweep pointed the same way for a bucket. Under 1 dB at the
grid we have, but not converged.

## 3. The bucket: the analytic lognormal holds while `sigma_R^2 <= 1`

Bucket 10 cm, the analytic lognormal against the empirical fade (delta at
5 and 1 percent, dB; a negative value is an OPTIMISTIC model):

| cell | `sigma_R^2` | index measured | index analytic | delta p5 | delta p1 | route | best refit family |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2 km, 3e-15 | 0.21 | 0.015 | 0.009 | -0.2 | -0.3 | FREE | any |
| 2 km, 1e-14 | 0.71 | 0.048 | 0.031 | -0.5 | -0.6 | FREE | gamma-gamma |
| 5 km, 3e-15 | 1.14 | 0.15 | 0.12 | -0.5 / -0.6 | -0.8 / -0.9 | FREE (standard) / REFIT (rapid) | gamma-gamma |
| 5 km, 1e-14 | 3.81 | 0.43 / 0.46 | 0.41 | -1.0 / -0.8 | -1.8 / -1.3 | REFIT | gamma-gamma, lognormal-Rician |
| 10 km, 3e-15 | 4.07 | 0.68 / 0.75 | 0.75 | -0.3 / -0.5 | -0.4 / -0.9 | FREE | lognormal (MLE) |
| 10 km, 1e-14 | 13.6 | 0.85 / 1.13 | 2.51 | +4.6 / +2.7 | +6.7 / +4.2 | NONE / REFIT | lognormal (refit) |

Four readings:

- The 1-6 verdict extends to `sigma_R^2 = 0.7`: the lognormal SHAPE holds
  and the analytic index reads about 1.5 times low, which costs 0.2 to 0.6 dB
  at 5 percent. At `sigma_R^2 = 1.1` the analytic route is on the 0.5 dB
  edge.
- From `sigma_R^2 = 3.8` the lognormal shape itself goes: the refit
  lognormal misses the 1 percent fade by 1.3 to 1.6 dB, and the gamma-gamma
  fitted by maximum likelihood holds (0.0 / 0.0 and 0.2 / 0.4 dB). This is
  the book's statement that the gamma-gamma is the model of every regime.
- The 10 km / 3e-15 FREE reading is luck, not physics: the analytic index
  (0.75) sits at the measured value because the Dios on-axis index and the
  Churnside filter miss in opposite directions there.
- The saturated cell (`sigma_R^2 = 13.6`) is not a clean test. The two
  presets disagree on the index (0.85 against 1.13), the skew of `ln P` is
  POSITIVE (+0.4 / +0.1, no other cell is), and the 2-TC1 sweep found that
  cell screen-count sensitive for a bucket. The analytic index (2.5) is two
  to three times too high there, so the fidelity-0 Term is PESSIMISTIC in
  saturation. Do not read a family verdict from it.

The 5 cm bucket reads the same way with a wider margin (FREE in 8 of 12; the
index is larger, the filter error smaller). The POINT irradiance follows the
analytic Dios on-axis index to 3 percent up to `sigma_R^2 = 1.1` (0.073 /
0.074, 0.249 / 0.248, 0.45 / 0.43) and the gamma-gamma MLE holds it in 11 of
12 cells, the lognormal only while weak. The MMF light bucket (a 100 um core
of NA 0.22 behind f = 0.25 m) is a bucket here: its index equals the 10 cm
bucket index to three digits in every cell, because the spot walk (a few
micrometres) never reaches the core edge.

## 4. The fibre: no free route holds, and the reason is the tilt

The SMF 10 cm untracked coupled power against the shipped fidelity-0 chain
(the mean-only higher-order coupling Term plus the exponential walk-off
Term), delta at 5 / 1 percent, dB, and the best fitted family:

| cell | index | skew `ln P` | p5 | p1 | shipped F0 chain | lognormal-Rician (Strehl map) | lognormal-Rician (MLE) | gamma-gamma (MLE) | route |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2 km, 3e-15 | 0.08 | -0.7 | 2.4 | 4.0 | -1.3 / -2.1 | +6.6 / +11 | 0.0 / -0.5 | 0.0 / -0.6 | REFIT |
| 2 km, 1e-14 | 0.34 | -1.3 | 6.6 | 10.3 | -2.8 / -3.8 | +6.1 / +9.4 | +0.1 / +0.8 | -0.3 / -0.9 | REFIT |
| 5 km, 3e-15 | 0.36 | -1.8 | 6.4 | 11.1 | -4.3 / -7.3 | +6.2 / +8.6 | -0.1 / -0.7 | -0.1 / -1.7 | REFIT |
| 5 km, 1e-14 | 1.65 | -0.9 | 16.3 | 23.3 | -8.7 / -10.4 | -1.9 / -1.8 | -1.5 / -1.4 | +0.1 / +1.8 | NONE / REFIT (K) |
| 10 km, 3e-15 | 1.60 | -0.8 | 14.3 | 20.8 | -10.8 / -14.7 | +0.8 / +1.5 | -0.2 / +0.4 | -0.2 / +0.2 | REFIT |
| 10 km, 1e-14 | 4.0 | -0.4 | 17.6 | 25.1 | -4.9 / -3.7 | +0.5 / +0.2 | +0.8 / +0.6 | +2.1 / +3.3 | FREE (speckle) |

(the `rapid` row of each cell; the `standard` rows agree inside the bars.)

Five readings:

- THE SHIPPED CHAIN UNDER-READS THE FIBRE FADE EVERYWHERE, by 1.3 dB at
  5 percent on the weakest cell and by 8 to 11 dB on the strong ones. It has
  no higher-order fade (the coupling Term is mean-only) and its walk-off tilt
  is too small (Section 5). The composite of the lognormal bucket times the
  walk-off reads the same way. So there is NO free route for a 10 cm fibre.
- THE TILT IS THE HEAVY TAIL. The skew of `ln P` is -0.7 to -1.8 untracked
  and -0.2 to -0.5 with the tilt removed. Removing the tilt takes the 5
  percent fade from 2.4 to 1.8 dB (2 km, 3e-15), from 6.6 to 3.5 dB (2 km,
  1e-14) and from 6.4 to 3.9 dB (5 km, 3e-15): 0.6, 3.1 and 2.5 dB of tilt
  at 5 percent, 1.5, 5.2 and 5.4 dB at 1 percent, on the three unaliased
  cells. With the tilt removed EVERY fitted family holds on those cells,
  and even the Strehl-map lognormal-Rician comes inside 0.8 / 1.4 dB. The 5 cm
  fibre carries far less tilt (0.1 to 0.5 dB at 5 percent), so its fade is a
  gamma-gamma to 0.1 dB in 9 of 12 cells and the composite free route holds
  on the two weakest.
- THE FITTED FAMILY OF RECORD IS THE LOGNORMAL-RICIAN. Fitted by maximum
  likelihood it holds the untracked 10 cm fibre in 9 of 12 cells (the two
  5 km / 1e-14 cells and the `rapid` saturated cell fail by 0.8 to 1.5 dB at
  5 percent). The gamma-gamma holds in 5 of 12: it misses the 1 percent
  tail by 0.9 to 2.0 dB at `sigma_R^2 = 0.7` to `1.1`, where the two-scale
  amplitude model has no room for a phase tail. A two-parameter fit is the
  floor: the one-parameter refit lognormal fails on every untracked 10 cm
  case past the weakest cell.
- THE STREHL-MAP LOGNORMAL-RICIAN (olb's heuristic) IS WRONG UNTIL SATURATION. With
  `S = exp(-sigma_phi^2)` from the FULL Noll residual the coherent part is
  too small at D/r0 of 1 to 2 (the 2 km and 5 km cells): the family reads 6
  to 12 dB pessimistic. The plan's Section 3.1 guessed this: a 5 to 10 cm
  aperture on a few km holds one or two phase cells, so a Gaussian random
  part is the wrong limit. At `sigma_R^2 = 13.6` (D/r0 of 5) the map lands
  (+0.5 / +0.2 dB), the negative exponential is 4 to 5 dB away, and the
  fibre is at the speckle limit that the plan predicted from the downlink
  numbers. Between the two, at `sigma_R^2 = 4`, the map is 1 to 2 dB out.
- THE TRACKED FOCUS CHANGES THE MEAN, NOT THE FAMILY. The received-curvature
  focus shift (S. A. Self, DOI 10.1364/AO.22.000658) moves the coupling
  mean, and on the 2 km cells the untracked-and-tracked case reads a
  slightly heavier tail (skew -1.6 to -2.4) that no family holds (NONE).
  Everywhere else the tracked and the untracked rows agree inside the bars.

## 5. The tilt: the aperture angle of arrival at the Gaussian-beam r0, with the outer scale

The per-axis tilt angle variance measured from the Noll pair, against the
analytic tilts, on the eight unaliased tilt tables (the four 2 km cells, the
two 5 km / 3e-15 cells, at 10 cm and 5 cm):

| tilt model | measured / model, 10 cm | measured / model, 5 cm |
| --- | --- | --- |
| the beam-wander arrival tilt that `terrestrial_smf_walkoff_term` reads (`olb.turbulence.angle_of_arrival.wander_arrival_angle_variance`) | 1.65 to 2.1 | 2.1 to 2.6 |
| the Noll Zernike tilt `0.182 (D/r0)^(5/3) (lam/D)^2` at the PLANE-wave r0 | 0.28 | 0.28 |
| the Noll Zernike tilt at the GAUSSIAN-beam r0 of the coupling Term (`gaussian_fried_parameter_profile`) | 0.74 (0.73 to 0.76) | 0.74 (0.73 to 0.76) |
| the same, times the von Karman outer-scale factor of Andrews Ch. 6, Eq. (83) at `L0 = 25 m` (0.763 at 10 cm, 0.812 at 5 cm; `olb.turbulence.andrews.structure.angle_of_arrival_variance`) | 0.97 | 0.91 |

Three readings:

- THE RECEIVED TILT IS AN APERTURE TILT, NOT A BEAM-FRAME WANDER. It grows
  from 10 cm to 5 cm by 1.26, which is the `D^(-1/3)` law of the angle of
  arrival to three digits; the wander model does not move with D at all.
  Its size is the Noll tilt at the Gaussian-beam r0 (the plane-wave r0 reads
  3.6 times too high), reduced by the outer scale: with the book's von
  Karman factor at the operating `L0 = 25 m` the prediction lands inside
  3 percent at 10 cm and 9 percent at 5 cm, on every unaliased cell. The
  constant 0.74 is therefore the outer scale, the same effect that moves the
  downlink SMF tail by 2.5 to 2.8 dB (2-P5).
- THE WALK-OFF TERM READS HALF THE TILT VARIANCE, so its exponential fade
  (linear in the variance) is about half the true tilt fade in dB. That is
  the first half of the shipped-chain deficit of Section 4. The Andrews
  Ch. 6 aperture route already sits in `andrews/structure.py` and in
  `angle_of_arrival.aperture_arrival_angle_variance`, unused by the
  terrestrial Term (backlog Conflict C-01). Re-pointing the Term at the
  aperture tilt with the Gaussian r0 and the site `L0` is a one-line physics
  change that moves every terrestrial fibre budget (an owner decision).
- Past `sigma_R^2 = 3.8` the slope sensor aliases, and the measured ratios
  fall (0.63 to 0.68 at 5 km / 1e-14 and 10 km / 3e-15, 0.40 to 0.54 at
  10 km / 1e-14) because the wrapped gradient loses tilt. Those rows are
  not a test of the model.

## 6. The verdict and the decision for step 4

| receiver | free route | fitted family that holds | where nothing holds |
| --- | --- | --- | --- |
| bucket (and the MMF light bucket) | the fidelity-0 analytic lognormal, while `sigma_R^2 <= 0.7` (inside 0.6 dB at 5 percent) and at the 0.5 dB edge at 1.1 | gamma-gamma (MLE) in every cell to `sigma_R^2 = 4`; the refit lognormal to 1.1 | the saturated 10 km / 1e-14 cell, which is also grid-limited |
| point irradiance | the analytic lognormal, while weak | gamma-gamma (MLE), 11 of 12 | the saturated `standard` cell, by 0.5 / 1.0 dB |
| SMF 5 cm (D/r0 to 1) | the composite lognormal times walk-off on the two weakest cells only | gamma-gamma (MLE) to 0.1 dB in 9 of 12 | the saturated `rapid` cell |
| SMF 10 cm (D/r0 1 to 5) | none, except the speckle limit in saturation | lognormal-Rician (MLE) in 9 of 12; the gamma-gamma in 5 | 5 km / 1e-14 and the saturated `rapid` cell |
| SMF, tilt removed | none (the Strehl-map lognormal-Rician is 0.8 to 2 dB pessimistic) | every family, on the unaliased cells | the aliased cells are not a result |

So the terrestrial fidelity-1 rung of 1-9 CANNOT be a free analytic draw for
a fibre, and for a bucket it is the fidelity-0 Term that already exists up
to `sigma_R^2` of about 1. The rung that the data supports is the CALIBRATED
route of 1-8, generalised: a short `Campaign` for the scenario (a few
hundred trials; the family parameters of the bucket and the 5 cm fibre are
stable at 2000), a maximum-likelihood fit of the family of record (the
gamma-gamma for a bucket, the lognormal-Rician for a fibre), and a Term
built through `olb/models/fade.py` from the fitted parameters, with the
measured mean coupling as its mean face. Three things the owner must decide
before that is built:

1. Whether a fidelity-1 rung that needs a short field campaign is the rung
   wanted (the 1-8 proposal said yes; it costs minutes on bigfraw and
   seconds on the GPU backend for each scenario).
2. The lognormal-Rician needs a CDF, a quantile and a sampler in
   `andrews/distributions.py` (it holds the PDF only, 0-W7). This study built
   the numeric CDF in the script; the Term needs it in the package.
3. Whether to re-point `terrestrial_smf_walkoff_term` at the aperture
   angle-of-arrival tilt (the Gaussian r0, the site `L0`), which changes
   every terrestrial fibre budget at fidelity 0. It fixes half the deficit;
   the other half (the higher-order phase fade of 1.8 to 3.9 dB at 5 percent
   on the unaliased cells) has no analytic Term today.

## 6a. A rule for the Rice factor: it exists for the higher-order part only (2026-09-08)

The fitted lognormal-Rician holds the fibre, so the owner asked whether a
RULE for r against D/r0 could make it a free route. `r_rule.py` tests that.

THE IDENTIFIABILITY PROBLEM FIRST. At a small index the two-parameter fit
cannot tell r from `sigma_z^2`: on the weak cells it pushed `sigma_z^2` to
zero and put the real bucket scintillation (index 0.015 to 0.05) into r. So
the script PINS `sigma_z^2 = ln(1 + sigma2_P)` to the bucket index and fits
r alone. The measured index and the ANALYTIC index of the fidelity-0 Term
give the same r inside a few percent, so the analytic Term can feed
`sigma_z^2`. In words: `sigma_z^2` is the variance of the log of the
lognormal envelope z in `I = z |A + X|^2`, the amplitude-scintillation
factor, the same lognormal the bucket uses; r carries the phase.

THE RULE, a power law `r = a x^b` on the pinned r, judged on each cell with
the law fitted on the OTHER cells (leave one cell out), pass rule as above:

| case | x = D/r0 (Gaussian r0) | x = Noll residual `sigma2_phi` | scatter in ln r | leave-one-out |
| --- | --- | --- | --- | --- |
| tilt removed (12 cases, 3 unaliased cells x 2 apertures x 2 presets) | `r = 18 (D/r0)^(-1.7)` | `r = 2.3 / sigma2_HO^1.03` | 0.41 (a factor 1.5) | 12 of 12 hold; worst 0.45 / 0.9 dB |
| untracked (24 cases, r > 0.05 on 18) | `r = 3.0 (D/r0)^(-2.9)` | `r = 3.1 / sigma2_full^1.7` | 0.53 (a factor 1.7) | 8 of 24 hold; misses of 3 to 6 dB |

Three readings:

- WITH THE TILT REMOVED a one-parameter rule works, and the `sigma2_HO` form
  is the Strehl map's small-residual limit (`r ~ 1/sigma2`) times 2.3: the
  higher-order phase makes the overlap fluctuate about 2.3 times less than a
  Gaussian phasor cloud of the same lost power would, the few-mode argument
  of Section 4. It is tested to `sigma2_HO = 0.15` and D/r0 = 1.1 only, on
  one launch (5 mm collimated), one wavelength, the optimal coupling
  parameter and `L0 = 25 m`. The hope that the Gaussian-beam r0 of
  `gaussian_fried_parameter_profile` self-normalises another waist is an
  ASSUMPTION until a second waist is run.
- THE UNTRACKED FIBRE HAS NO RULE IN r ALONE. The tilt is one bounded mode
  with its own exponential-in-dB statistic; stuffing it into the Rice factor
  of a Gaussian-cloud model distorts the whole curve, and the r that the fit
  wants falls steeply and scatters. The fitted family still DESCRIBES the
  untracked fade (9 of 12, Section 4); it is the PREDICTION of r that fails.
- THE ROUTE THIS POINTS TO is a free COMPOSITE for the fibre:
  `P/<P> = LR_HO(r = 2.3 / sigma2_HO, sigma_z^2 = ln(1 + sigma2_P)) x
  walk-off(tilt)`, with the walk-off fade olb already has but fed the RIGHT
  tilt variance (the aperture angle of arrival at the Gaussian r0 times the
  outer-scale factor of Section 5, measured to 3 percent), not the
  beam-wander one. Every piece is analytic and the 2.3 is the one fitted
  constant. NOT BUILT (owner pause, 2026-09-08). Its open assumption is the
  independence of the tilt factor from the higher-order factor (tilt and
  defocus are correlated modes); the test is to add the composite to
  `fit_distributions.py`, once with the modelled tilt and once with the
  MEASURED tilt variance, so a miss can be blamed on the tilt model or on
  the independence separately.

## 7. Files

| file | purpose |
| --- | --- |
| `extract_trials.py` | the one-pass read of the campaigns, on bigfraw |
| `fit_distributions.py` | the fits, the verdicts, the figures |
| `r_rule.py` | the pinned-`sigma_z^2` fit of r, the power-law rules, the leave-one-cell-out test (Section 6a); writes `r_rule.log` and `r_rule.json`, git-ignored |
| `extract_trials.log` | the extraction record (12 cells, 32 minutes, 8 workers); git-ignored, on disk |
| `fit_distributions.log` | every table of every cell; git-ignored, regenerated by the fit script |
| `fit_results.json` | every number of the log, and the fitted parameters; git-ignored, regenerated |
| `trials/trials_<cell>.npz` and `.json` | the per-trial tables and their metadata; git-ignored, on bigfraw and the laptop |
| `figures/cdf_bucket_10cm.png` | the bucket exceedance curves against the families |
| `figures/cdf_smf_10cm_untracked.png` | the untracked 10 cm fibre |
| `figures/cdf_smf_10cm_tilt_removed.png` | the same with the tilt removed |

## 8. What was not done

- A FOCUSED launch (blocked by 0-P17: olb cannot express a converging beam).
- The Zernike Monte Carlo of the plan (Section 3.2 of the 2026-09-05 note,
  now Section 9 below): the field data decided the families without it.
- The temporal extension (Section 10 below).
- The saturated cell needs a finer grid before its bucket verdict counts.

## 9. The candidate models (the 2026-09-05 plan, kept for the record)

The plan proposed four routes for the fibre fade. The sweep judged them:

- 3.1, lognormal times Rician with `r = S/(1-S)` from the Noll residual: the
  FAMILY is right (it is the fitted family of record), the PARAMETER MAP is
  wrong below saturation, as the plan's Section 6 suspected (one or two
  phase cells across a 5 to 10 cm aperture, not a Gaussian random part).
- 3.2, the Zernike Monte Carlo with no propagation: not built; a candidate
  for the parameter map at D/r0 of 1 to 3, where the map fails.
- 3.3, the coupled-power variance from the fourth-order coherence: not
  started.
- 3.4, the receiver-kind split: done (Sections 3 and 4).

## 10. The temporal extension (not started)

The fade DEPTH work above has a temporal twin. A first guess, "convolve the
bucket temporal spectrum with the tilt spectrum", is WRONG: the coupling is
a nonlinear function of ALL the modes, each with its own temporal spectrum
(Conan, Rousset and Madec 1995, DOI 10.1364/JOSAA.12.001559: the tilt falls
as `f^(-2/3)` at low frequency, every higher order is flat to a knee near
`0.3 (n + 1) v / D` and rolls off as `f^(-17/3)`), so the higher orders put
their power at HIGHER frequencies than the tilt. The product of two
independent processes gives the convolution of their spectra, so
`P_bucket * eta` is that product, and the gap is the `eta` spectrum itself:
either the frozen-flow axis of `olb/waveoptics/turbulence/temporal.py`
(the truth), or a Zernike temporal Monte Carlo drawn from the Conan spectra
(the cheap twin of 3.2).

## 11. Conventions

- ASD-STE100 in every docstring, comment and commit message.
- Every equation cites a DOI.
- Loss is positive dB.
- Run as `python -m validation.fibre_fade_models.<script>` from the root.

# Past-reference audit of the olb documentation (2026-09-06)

Scope: docstrings and comments in `olb/**`, the reference docs (`README.md`,
`docs/physics.md`, `docs/architecture.md`, `docs/getting-started.md`,
`docs/examples.md`, `docs/api-*.md`), and the `examples/**/README.md` files.
Excluded on purpose: `docs/backlog.md`, the two crosscheck trackers, the two
plan docs, `validation/`, and `CLAUDE.md`.

Categories:
- A = dependency / vendoring history (my_analysis_modules, cache.py, fso_spot_size, TN-2)
- B = dated decision or status stamp ("(2026-08-27)", "owner decision", "DONE", "WIRED")
- C = a prior bug or prior behaviour described ("before the fix", "the old behaviour", "unchanged", "bit for bit the old record")
- D = prior name / rename / removed thing ("legacy", "retired name", "there is no X")
- E = backlog / work-package / gap / conflict id in prose (WP3a, 2-W2, Gap 3, C-05, TL-01, I-2)
- F = "now" / "still" / "today" relative-time phrasing
- G = a whole status-narrative section

Actions: DELETE, REWRITE (present tense only), MOVE (to backlog), KEEP (with reason).

## 1. Totals

| Area | Entries | Dominant | Worst files |
|---|---|---|---|
| olb root + models + links | 100 | F, E | links/terrestrial.py 17, links/uplink.py 13, links/downlink.py 12, models/coupling/terrestrial.py 12, terminal.py 8 |
| olb/turbulence + olb/waveoptics | 79 (8 KEEP) | C, D | waveoptics/turbulence/sampling.py 12, waveoptics/turbulence/run.py 12, turbulence/coupled_flux.py 7, andrews/paths.py 6, andrews/structure.py 5 |
| Reference docs + example READMEs | 75 | E, B, G | docs/physics.md 23, docs/architecture.md 13, examples/waveoptics/README.md 10 |
| docs/api-*.md | 38 | B | api-waveoptics.md 24, api-budget.md 12, api-terminal-scenario.md 2 |
| **Total** | **292** | | |

Clean files: `olb/__init__.py`, `units.py`, `beam.py`, `geometry.py`, `multidetector.py`, `models/__init__.py`, `models/geometric.py`, `models/pointing.py`, `models/splitter.py`, `models/coupling/__init__.py`, `docs/README.md`, `examples/README.md`.

## 2. FIX FIRST: notes that are not just stale but WRONG

1. `olb/waveoptics/turbulence/screens.py:15-16, 224-227` say aotools is the default generator and a required dependency. The default is "olb" (`run.py:333`); aotools is an optional lazy import.
2. `olb/waveoptics/__init__.py:30`, `olb/waveoptics/run.py:8-9`, `olb/waveoptics/turbulence/__init__.py:5-6`, `olb/waveoptics/turbulence/run.py:7` all say the layer "builds NO Term and changes NO budget". `olb/models/waveoptics.py` builds Terms from both runners.
3. `docs/api-budget.md:257-258`: "re-exports todB ... from my_analysis_modules". They are vendored in `olb/units.py`.
4. `olb/links/uplink.py:315-316` ships `my_analysis_modules.satellite.SatellitePass.point_ahead_angle()` inside Term validity text, so it appears in `assumptions_frame()`. The geometry lives in `olb.geometry`.
5. `olb/models/coupling/terrestrial.py:556`: "is replaced by the DEFOCUS-ABERRATED closed form" ships in Term validity text.
6. `examples/andrews/README.md:95-101` says the terrestrial fibre call site passes no f0; `docs/physics.md:639` says it does. It also points at a CLAUDE.md "Next task" section that does not exist.
7. `docs/api-waveoptics.md:589-603` quotes an "estimated 2 dB at p5"; the measured value is 2.5 to 2.8 dB (`validation/outer_scale_tail/`).
8. `docs/getting-started.md:16-18`: the my_analysis_modules sentence, plus a pointer to a README "Dependency" section that is actually called "Dependencies".
9. `olb/links/terrestrial.py:165`: "See memory dios-scintillation-convergence / pointing-jitter-into-beta" names private assistant memory notes no reader can open.

## 3. The headline complaint: my_analysis_modules and other vendoring history (A)

| Location | Text | Action |
|---|---|---|
| docs/architecture.md:496 | heading "Self-contained: the vendored physics (formerly my_analysis_modules)" | REWRITE heading |
| docs/architecture.md:498-500 | "olb no longer depends on my_analysis_modules. It once borrowed ... through olb/_deps.py" | DELETE |
| docs/architecture.md:510-512 | "_deps.py is deleted. Each vendored copy is verbatim ... cross-validated bit-for-bit against the original" | DELETE first sentence and "against the original" |
| docs/getting-started.md:16-18 | "kernels it once borrowed from my_analysis_modules are now vendored" | DELETE |
| docs/api-budget.md:257-258 | re-export from my_analysis_modules | REWRITE (wrong) |
| docs/physics.md:499, 506 | "the vendored kernels", "The kernels are vendored in" | REWRITE ("the kernels are in") |
| olb/turbulence/coupled_flux.py:5-7, 19, 502-504 | "borrowed from the my_analysis_modules kernel ... vendored here VERBATIM"; "nothing from my_analysis_modules"; "cross-check against the my_analysis_modules original is in the commit" | REWRITE / DELETE |
| olb/turbulence/andrews/paths.py:826-828 | "The kernel ... in my_analysis_modules/coupled_flux.py returns a DIFFERENT quantity" | REWRITE (point at olb.turbulence.coupled_flux) |
| olb/turbulence/ao.py:192-196 | "NEW HOME ... used the Fried 1966 constant 0.423 before" | REWRITE |
| olb/turbulence/gaussian_fried.py:30-34 | "NEW HOME. ... now live in the Andrews foundation package ... call the new home" | REWRITE ("delegate to andrews.beam") |
| olb/waveoptics/smf.py:8-12 | "transcribes two helpers from the shared kernel repository (my_analysis_modules ...)" | DELETE (keep the LightPipes attribution at 24-25) |
| olb/waveoptics/turbulence/campaign.py:17-21 | "HISTORY. This module REPLACED the P4 scalar cache (cache.py, retired 2026-09-04)" | DELETE |
| olb/waveoptics/turbulence/fingerprint.py:21-23, 35-37 | "HISTORY. The key was born in the retired cache.py ..." | DELETE / REWRITE |
| docs/api-waveoptics.md:1120-1122 | "This key came from the P4 scalar cache (cache.py) ... RETIRED on 2026-09-04" | DELETE |
| olb/models/extinction.py:24, 44 | "comes from fso_spot_size.airmass in the sibling TN-2 analysis repo"; "borrowed from fso_spot_size.airmass" | REWRITE |
| olb/models/gaussian_efficiency.py:23-26, 220-222 | "corrected antenna-gain form ... validated against tn2_kepler test_gauss_prop"; "The stale form (with 2/alpha^2) gives eta > 1" | REWRITE |
| olb/links/uplink.py:315-316 | my_analysis_modules.satellite in shipped validity text | REWRITE |

## 4. Whole status-narrative sections (G): owner to decide

| Location | Section | Suggested |
|---|---|---|
| README.md:131-275 | "## Roadmap" (mermaid status trees, "Next / planned", node id NT6 in prose) | MOVE to backlog, leave a pointer |
| README.md:200-223 | "Fidelity ladder — the fibre-coupling instance" mixes description with "is WIRED (2026-08-28)" | REWRITE, strip status |
| docs/physics.md:1658-2111 | "## 9. Measured validity": every entry "Measured (YYYY-MM-DD)" | KEEP as a measurement record; PRUNE backlog ids and "earlier run was wrong" sentences |
| docs/architecture.md:419-458 | "Newly enforced constraints, and the open follow-ups" (refactor change-log; "CLOSED (2026-09-04)" twice; "old name WEAK_FLUCTUATION_LIMIT is fully retired") | MOVE; keep a short list of enforced constraints |
| docs/api-waveoptics.md:21-29, 493-500 | "Status: the layer IS wired ... (2026-08-28) ... remaining owner gate" | REWRITE as a scope note; MOVE the owner gate |
| docs/api-waveoptics.md:1197 | heading "boost: the process priority boost (2026-09-05)" | drop the date |
| examples/andrews/README.md:40-101 | "## Wiring status: what is LIVE and what is NOT" | owner to decide (one entry is stale, see 2.6) |
| examples/waveoptics/README.md:105-118 | "The FAST comparison of the downlink script": superseded numbers marked "Do not take those numbers as current" | REWRITE to one sentence pointing at validation/waveoptics_vs_fast |
| examples/waveoptics/README.md:191-231 | "## Status: wave optics is WIRED at fidelity=2" | MOVE / REWRITE |
| examples/schmidt/README.md:78-91 | "## Wiring status / LIVE: nothing"; WP6 retrofit narrative | REWRITE to "validation only" |

## 5. Recurring stamps (one grep-and-strip pass each)

- **Work-package prefixes** `WP3a`, `WP3b`, `WP3c`, `WP3d`, `WP6`, `WP7`, `P0`, `P4`: olb/models/fade.py:208, fast.py:677/765, waveoptics.py:1119, coupling/downlink.py:359, coupling/terrestrial.py:1230, links/uplink.py:890, links/downlink.py:683/870/942, links/terrestrial.py:648, links/retro_space.py:275/293, waveoptics/turbulence/sampling.py:133/148/184, screens.py:203, docs/physics.md:1585/1644/2041, docs/api-waveoptics.md:720, examples/waveoptics/README.md:101/102, examples/schmidt/README.md:87/106.
- **Backlog / crosscheck ids** `I-1`, `I-2`, `I-4`, `1-1`, `1-2`, `0-W1`, `0-W3`, `0-W4`, `0-P11`, `2-W1`, `2-W2`, `2-W3`, `2-I2`, `2-I2T`, `2-N6`, `2-P5`, `Gap 2`, `Gap 3`, `S-06/09/14/15`, `C-01`, `C-04`, `C-05`, `TL-01..05`, `NT6`: olb/scenario.py:187/276, sweep.py:96/112, models/fade.py:154, coupling/_common.py:130, coupling/terrestrial.py:516/1197/1265/1393, links/downlink.py:90/97/258, links/terrestrial.py:92/176/408/598, andrews/paths.py:727, waveoptics/turbulence/sampling.py:468/600/729/908, schmidt/__init__.py:28, docs/physics.md:639/646/798/1096/1332/1347/2041/2099, docs/architecture.md:50, docs/api-budget.md:418/646, docs/api-waveoptics.md:589/653/1531, README.md:225, examples/andrews/README.md:34/73/86, examples/schmidt/README.md:99/111/113/126/131. Two live doc pointers may KEEP (links/downlink.py:428, links/terrestrial.py:180) if andrews-crosscheck.md stays authoritative.
- **Dates** 2026-08-25/27/28/29, 2026-09-04/05: olb/assumptions.py:67/542, coupled_flux.py:411/468/570, uplink_flux.py:227, andrews/paths.py:170/715, coupling/terrestrial.py:272/771/1250, links/uplink.py:326/586, waveoptics/priority.py:6/16, campaign.py:769, fingerprint.py:21, turbulence/run.py:264/1285, README.md:211/220, docs/physics.md:639/798/1352/1794/1913/2099, docs/architecture.md:154/446/456, docs/api-budget.md:418, docs/api-waveoptics.md:21/589/809-825/961/1084/1104/1120/1197/1211/1216/1531, examples/andrews/README.md:73, examples/waveoptics/README.md:207.
- **"now" / "still" / "today" / "unchanged" / "the old behaviour" / "bit for bit the old record"**: the largest bucket. olb/terminal.py:135/147/152/203/221/398/401/415/447; links/terrestrial.py:234/578/654/662/667/792/821/870/902; links/uplink.py:101/149/813/897/906/1092; links/downlink.py:39/709/733/948; models/fast.py:297/569; coupling/terrestrial.py:745; bidirectional.py:41/112; results.py:299; waveoptics/turbulence/run.py:334/343/528/552/555/1134; mmf.py:186; andrews/aperture.py:564, beam.py:272, scintillation.py:954, structure.py:721/731/752/764; beam_wave_scintillation.py:54; plane_wave_scintillation.py:54; uplink_flux.py:376; docs/examples.md:118/344/428/468/472; docs/api-waveoptics.md:197/247/256/486/838/842/844/1041; docs/api-budget.md:120/359/730/780; docs/api-terminal-scenario.md:134/221; docs/physics.md:848/895/1056; docs/architecture.md:98/129/437/440. Fix: state the present default.
- **"legacy array planner"**: waveoptics/turbulence/sampling.py:468/600/732/930/982, turbulence/run.py:515/755, campaign.py:258, models/waveoptics.py:131/825, docs/api-waveoptics.md:659. The array planner is still supported, so call it "the discrete-array planner".
- **"retired inline faces" migration regressions**: models/fade.py:154-200 (5 sites incl. a print string), links/terrestrial.py:598-618 (incl. a print string). The assertions stay; rename the baseline "the reference closed form".
- **Prior-bug narratives** (C, DELETE): coupled_flux.py:411-412 ("a parenthesis closed too early"), 468-471 ("An earlier patch removed this weight ... back (2026-08-25)"); gaussian_fried.py:495 ("the old (1 - xi) weight"); uplink_flux.py:42/227/466 (old vertical-grid call); andrews/paths.py:724 ("An earlier docstring called it the floor ... wrong"); docs/physics.md:1784 (defocus SIGN fault), 562-566 (flag migrated from the old axis); docs/api-waveoptics.md:256 ("the sign is now RIGHT-WAY-ROUND"), 1084-1090 (uplink patch too small before 2026-09-05), 1211 (validation scripts patched the initializer by hand); docs/architecture.md:421-423 ("NOT a regression"), 443-445; links/downlink.py:424 ("NOT the factor-of-4 conflation"); coupling/terrestrial.py:304 ("the old factory patch"); api-budget.md:296/593 (old knobs `smf_fidelity`, `model="montecarlo"` gone).
- **"TILT DEFINITION - THE OWNER MADE THIS CHOICE"**: turbulence/angle_of_arrival.py:135, andrews/structure.py:483. State the convention.
- **Removed-thing statements** (D): links/downlink.py:587/638 ("there is no 'no detector' special case"); docs/physics.md:1182 and architecture.md:443 (retired `WEAK_FLUCTUATION_LIMIT`); physics.md:877-889 ("Old home / New home" delegation table, "Nothing above changed"); api-budget.md:630 ("there is no retro.py file").

## 6. KEEP (live alias or load-bearing lineage)

| Location | Why |
|---|---|
| olb/links/__init__.py:11 and docs/api-budget.md:630-632 | `retro_budget` alias is defined and exported |
| olb/waveoptics/lenses.py:76-79, propagators.py:113/216/226 | "legacy" C++ LightPipes constants justify live literals |
| olb/waveoptics/smf.py:24-25 | BSD-3 LightPipes attribution |
| olb/waveoptics/mmf.py:192-193 | explains a live unused argument |
| olb/waveoptics/schmidt/turbulence.py:38-39 | numpy trapz/trapezoid shim |
| olb/waveoptics/turbulence/campaign.py:369-371 | explains the live manifest `defaults` fallback |
| olb/waveoptics/turbulence/run.py:943-947 | aotools 1.0.7 scipy workaround |
| docs/architecture.md:108-109, docs/api-waveoptics.md:519 | temporal.py stub status is a live API fact |
| docs/examples.md:189-197 | "the old Gaussian roll-off" is what the figure plots (rename "naive") |
| examples/schmidt/README.md:133-149 | book errata, not repo history |
| olb/models/gaussian_efficiency.py:229/279 | "TN-2 launch" names an example scenario; relabel if TN-2 is not a live reference |

## 7. Per-file counts

| File | Entries | Dominant |
|---|---|---|
| olb/links/terrestrial.py | 17 | F/E |
| olb/links/uplink.py | 13 | C/F |
| olb/links/downlink.py | 12 | E/F |
| olb/models/coupling/terrestrial.py | 12 | C/E |
| olb/waveoptics/turbulence/sampling.py | 12 | D/E |
| olb/waveoptics/turbulence/run.py | 12 | C |
| olb/terminal.py | 8 | F |
| olb/turbulence/coupled_flux.py | 7 | A/C |
| olb/models/fade.py | 6 | C |
| olb/models/fast.py | 6 | C/E |
| olb/turbulence/andrews/paths.py | 6 | B |
| olb/turbulence/andrews/structure.py | 5 | D |
| olb/links/retro_space.py | 4 | E |
| olb/turbulence/uplink_flux.py | 4 | C |
| olb/waveoptics/turbulence/campaign.py | 4 | A/B |
| olb/models/waveoptics.py, coupling/downlink.py, gaussian_efficiency.py, waveoptics/turbulence/screens.py | 3 each | |
| olb/assumptions.py, scenario.py, sweep.py, models/extinction.py, links/bidirectional.py, turbulence/gaussian_fried.py, waveoptics/mmf.py, priority.py, smf.py, fingerprint.py | 2 each | |
| 13 further package files | 1 each | |
| docs/physics.md | 23 | E (+ section 9) |
| docs/api-waveoptics.md | 24 | B |
| docs/architecture.md | 13 | B/A |
| docs/api-budget.md | 12 | C/D |
| examples/waveoptics/README.md | 10 | G/C |
| docs/examples.md | 9 | F |
| examples/schmidt/README.md | 7 | E/G |
| README.md | 6 | B/G |
| examples/andrews/README.md | 6 | G/E |
| docs/api-terminal-scenario.md | 2 | F |
| docs/getting-started.md | 1 | A |

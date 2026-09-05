# Wave optics against FAST: the space-downlink SMF gap (2-W1, 2-AO)

This study measures backlog 2-W1: does the fidelity-2 wave-optics field read LESS
single-mode-fibre (SMF) coupling loss than the fidelity-1 FAST model on the space
downlink, and by how much? It reads the production layer and changes no `olb`
module.

## The method

Three models, per elevation, ALL UNCORRECTED (no adaptive optics, no tip-tilt).
That is the only fair comparison, because the fidelity-2 layer applies no
correction (backlog 2-AO). The hero downlink terminal has an empty compensation
stack, so FAST runs `AO_MODE = NOAO` on its own.

- **FAST** (fidelity 1): `smf_fast_term`. `mean_db` is the mean per-sample loss.
- **field** (fidelity 2): the split-step field through a `Campaign` PROCESS pool.
  The per-trial loss is `-10 log10(collected_power * smf_eta)` (the same
  composite as `validation/tail_convergence/`).
- **analytic** (fidelity 0): `downlink_coupling_term(smf_fidelity="mean")`.

PARITY (or the comparison is a confound): the SAME Cn2, the SAME outer scale
`L0` on the field AND FAST, mean-of-dB on both, and `defocus_m = 0` (a 500 km
downlink is far field, so the received curvature is about zero). The analytic
term is L0-agnostic, so it is an approximate reference, not a full parity match.
An NPXLS convergence guard pins the FAST grid first, because FAST can undersample
the low-order tilt in NOAO mode.

## THE FINDING: the gap is an OUTER-SCALE artifact

Measured 2026-09-05, 400 trials per case, 20/30/60/90 deg. The gap is
FAST minus field (positive = field reads LESS loss, the 2-W1 claim):

| elev | gap at L0 = 25 m | gap at L0 = inf |
| --- | --- | --- |
| 20 deg | -0.34 +-0.36 | +0.48 +-0.37 |
| 30 deg | -0.21 +-0.37 | +0.35 +-0.38 |
| 60 deg | +0.07 +-0.38 | +1.17 +-0.42 |
| 90 deg | +0.12 +-0.38 | +0.71 +-0.44 |

At the PHYSICAL `L0 = 25 m` (the operating choice, backlog 2-P5) FAST and the
field AGREE: the gap is -0.34 to +0.12 dB, all within 1 sigma of zero, 0 of 4 in
the backlog's 0.7-2.9 dB claim band. At `L0 = inf` the gap REOPENS to +0.35 to
+1.17 dB (2 of 4 in the band). So the old "field reads 0.7-2.9 dB less than FAST"
was largely an OUTER-SCALE / parity artifact, not a wave-optics error.

WHY: FAST is MORE outer-scale-sensitive than the field. Going inf -> 25 m, FAST
drops about 1.6 dB at 20 deg while the field drops about 0.8 dB. And FAST's NOAO
tilt is grid-defined at inf: the NPXLS guard pinned 512 at inf against 128 at
25 m. So the tilt the fibre pays is set by the outer scale, and once both models
see the same physical `L0` they converge. This is the same mechanism as backlog
2-P5 (the outer scale sets the fibre tilt).

The analytic fidelity-0 term stays about 1 to 2.5 dB OPTIMISTIC against both FAST
and the field (it is L0-agnostic).

CAVEAT: this certifies the UNCORRECTED rung only (backlog 2-AO). Fidelity 2
applies no correction, so this is NOT a reference for an AO-corrected link.

## Run it

From the repository root:

    python -m validation.waveoptics_vs_fast.waveoptics_vs_fast
    python -m validation.waveoptics_vs_fast.waveoptics_vs_fast --L0 inf
    python -m validation.waveoptics_vs_fast.waveoptics_vs_fast --field-mode thread

Outputs are tagged by outer scale and field mode
(`waveoptics_vs_fast_L0<value>_<mode>_results.json` and `.log`), and the figures
go to `figures/`. The field campaign blocks go to `campaigns/` (gitignored).

## Parallelism note (how to use the machine)

The field runs through a `Campaign` (backlog 2-N4a). Its `--field-mode`:

- `process` (default) -- a warm `ProcessPoolExecutor` of `--workers` processes,
  one block per process, no GIL. This SATURATES the machine.
- `thread` -- bare `run_fidelity2` with a `Threader`. The threads share memory
  but the GIL caps them at about 0.35 efficiency, so a thread run does NOT
  saturate.
- `serial` -- one trial at a time.

IMPORTANT: the effective process count is `min(workers, number_of_blocks)`, and
`number_of_blocks = ceil(n_trials / block_size)`. So `--block-size` must be small
enough to make at least `--workers` blocks, or the extra workers sit idle. For
`n_trials` trials on `W` workers, use `--block-size <= n_trials / W`. See
`docs/api-waveoptics.md` Section 9g.

## Sources

- O. J. D. Farley and others, Opt. Express 30(13), 23050 (2022),
  DOI 10.1364/OE.458659. The FAST method (fidelity 1).
- Y. Dikmelik and F. M. Davidson, Appl. Opt. 44(23), 4946 (2005),
  DOI 10.1364/AO.44.004946. The analytic uncorrected SMF coupling curve.
- C. Ruilier, Proc. SPIE 3350, 319 (1998), DOI 10.1117/12.317094. The 0.8145
  mode-match limit.
- Schmidt (2010), DOI 10.1117/3.866274, Ch. 9. The split-step field solve.
- Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196. Ch. 12 the
  Hufnagel-Valley profile; Ch. 8 the Rytov variance.

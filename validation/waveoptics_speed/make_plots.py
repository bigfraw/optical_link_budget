r"""Plot the fidelity-2 speed-plan results (P0 to P3) for a visual read.

The script reads the results JSON of each speed task in this folder and it
writes one PNG for each. It touches no production code and it runs no
simulation; it only draws the numbers the tasks already measured. A task whose
JSON is absent is skipped, so the script runs before P3 finishes and again
after.

Figures:
- profile_baseline.png     (P0) where one trial spends its time.
- screen_generator.png     (P1) the fast generator: speed and accuracy.
- grid_experiments.png     (P2) the two buried grid ideas, against their kill
                                 lines.
- generator_validation.png the broad drop-in check, including the fade tail.
- scaling_study.png        (P3) the parallel-scaling curves, when present.

Run from the repository root:
    python -m validation.waveoptics_speed.make_plots
"""

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(__file__)
DPI = 130
KILL = "#c1121f"          # the kill line / bad
OK = "#2a9d8f"            # good / olb
REF = "#264653"           # reference / aotools
ACCENT = "#e76f43"        # a third series


def _load(name):
    """Load one results JSON, or None when it is absent."""
    path = os.path.join(HERE, name)
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def _save(fig, name):
    """Save a figure and report the path."""
    path = os.path.join(HERE, name)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {os.path.relpath(path)}")


def _short(name):
    """Shorten a case label for an axis tick."""
    return (name.replace("space ", "").replace("downlink", "down")
            .replace("uplink", "up").replace(" 30deg", " 30")
            .replace("terrestrial", "terr").replace(" standard", " std")
            .replace(" rapid", " rap").replace("2km", "2km"))


# ---------------------------------------------------------------------------
# P0 - the profiling baseline
# ---------------------------------------------------------------------------
def plot_p0(d):
    cases = d["cases"]
    names = [_short(c["name"]) for c in cases]
    y = range(len(cases))
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))

    # A: the share stacked bar
    scr = [c["share_pct"]["screen_gen"] for c in cases]
    prop = [c["share_pct"]["propagation"] for c in cases]
    scal = [c["share_pct"]["scalar_reads"] for c in cases]
    ax[0].barh(y, scr, color=ACCENT, label="screen generation")
    ax[0].barh(y, prop, left=scr, color=REF, label="split step")
    ax[0].barh(y, scal, left=[s + p for s, p in zip(scr, prop)],
               color=OK, label="scalar reads")
    ax[0].set_yticks(list(y))
    ax[0].set_yticklabels(names, fontsize=8)
    ax[0].set_xlabel("share of one trial [%]")
    ax[0].set_title("P0: where one trial spends its time")
    ax[0].legend(fontsize=8, loc="lower right")
    for i, s in enumerate(scr):
        ax[0].text(s / 2, i, f"{s:.0f}%", va="center", ha="center",
                   color="white", fontsize=8)

    # B: base FFT vs subharmonic within screen gen (seconds)
    base = [c["timing_s"]["screen_base_fft"] for c in cases]
    sub = [c["timing_s"]["screen_subharmonic"] for c in cases]
    ax[1].barh(y, base, color=REF, label="base FFT draw")
    ax[1].barh(y, sub, left=base, color=KILL, label="subharmonic add")
    ax[1].set_yticks(list(y))
    ax[1].set_yticklabels(names, fontsize=8)
    ax[1].set_xlabel("screen-generation time [s]")
    ax[1].set_title("P0: subharmonics dominate screen gen (~7x the base)")
    ax[1].legend(fontsize=8, loc="lower right")

    # C: raw fft2 microbench
    fb = d["fft_microbench"]
    ns = sorted(int(k) for k in fb)
    npy = [fb[str(n)]["numpy_fft2_s"] * 1e3 for n in ns]
    scp = [fb[str(n)]["scipy_fft2_w1_s"] * 1e3 for n in ns]
    ax[2].loglog(ns, npy, "o-", color=REF, label="numpy fft2")
    ax[2].loglog(ns, scp, "s--", color=ACCENT, label="scipy fft2 (workers=1)")
    ax[2].set_xlabel("grid n")
    ax[2].set_ylabel("one fft2 [ms]")
    ax[2].set_title("P0: raw fft2 cost")
    ax[2].set_xticks(ns)
    ax[2].set_xticklabels(ns)
    ax[2].xaxis.set_minor_formatter(plt.NullFormatter())  # no overlap
    ax[2].grid(True, which="both", alpha=0.3)
    ax[2].legend(fontsize=8)

    fig.suptitle("P0 - fidelity-2 profiling baseline", fontweight="bold")
    _save(fig, "profile_baseline.png")


# ---------------------------------------------------------------------------
# P1 - the fast screen generator
# ---------------------------------------------------------------------------
def plot_p1(d):
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))

    # A: speed-up vs n
    sp = sorted(d["speed_per_screen"], key=lambda r: r["n"])
    ns = [r["n"] for r in sp]
    x = range(len(ns))
    make = [r["speedup_make"] for r in sp]
    stack = [r["speedup_stack"] for r in sp]
    w = 0.38
    ax[0].bar([i - w / 2 for i in x], make, w, color=OK, label="make")
    ax[0].bar([i + w / 2 for i in x], stack, w, color=REF,
              label="make_stack")
    ax[0].axhline(1.0, color="grey", ls=":", lw=1)
    ax[0].set_xticks(list(x))
    ax[0].set_xticklabels(ns)
    ax[0].set_xlabel("grid n")
    ax[0].set_ylabel("speed-up over aotools (x)")
    ax[0].set_title("P1: per-screen speed-up")
    ax[0].legend(fontsize=8)
    for i, v in enumerate(stack):
        ax[0].text(i + w / 2, v + 0.2, f"{v:.0f}x", ha="center", fontsize=8)

    # B: structure function ratio vs theory, with the [0.85, 1.02] band
    sf = d["structure_function"]
    ax[1].axhspan(0.85, 1.02, color=OK, alpha=0.15, label="accept band")
    ax[1].axhline(1.0, color="grey", ls=":", lw=1)
    ax[1].plot(sf["r_over_r0"], sf["ratio"], "o-", color=REF,
               label="olb D_phi / theory")
    ax[1].set_ylim(0.7, 1.1)
    ax[1].set_xlabel("r / r0")
    ax[1].set_ylabel("D_phi ratio to Fried theory")
    ax[1].set_title("P1: screen structure function")
    ax[1].legend(fontsize=8, loc="lower left")

    # C: statistical equivalence (200 trials) - sigma2_I with bars + wall time
    eq = d["statistical_equivalence"]
    labels = ["aotools", "olb"]
    s2 = [eq["aotools"]["sigma2_I"], eq["olb"]["sigma2_I"]]
    se = [eq["aotools"]["se_sigma2_I"], eq["olb"]["se_sigma2_I"]]
    ax[2].bar(labels, s2, yerr=se, capsize=6, color=[REF, OK])
    ax[2].axhline(eq["comparison"]["sigma2_I_analytic"], color=KILL, ls="--",
                  label="analytic")
    ax[2].set_ylabel("aperture sigma2_I (200 trials)")
    ap = eq["comparison"]["sigma2_I_sigmas_apart"]
    mp = eq["comparison"]["mean_power_sigmas_apart"]
    ax[2].set_title(f"P1: equivalence (sigma2_I {ap:.2f} sig, "
                    f"mean {mp:.2f} sig)")
    ax[2].legend(fontsize=8)
    # annotate the wall-time win
    wt_a, wt_o = eq["aotools"]["wall_s"], eq["olb"]["wall_s"]
    ax[2].text(0.5, max(s2) * 0.5,
               f"wall: {wt_a:.1f}s -> {wt_o:.1f}s\n({wt_a / wt_o:.1f}x)",
               ha="center", fontsize=8,
               bbox=dict(boxstyle="round", fc="white", ec="grey"))

    fig.suptitle("P1 - fast phase-screen generator (ScreenFactory)",
                 fontweight="bold")
    _save(fig, "screen_generator.png")


# ---------------------------------------------------------------------------
# P2 - coarse screens and the beam-following grid (both buried)
# ---------------------------------------------------------------------------
def plot_p2(coarse, beam):
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))

    # A: coarse screens - speedup vs sigma2 error, with the +/-5% kill band
    ax[0].axhspan(-5, 5, color=OK, alpha=0.15, label="+/-5% accept band")
    ax[0].axhline(0, color="grey", lw=1)
    case_colors = [REF, ACCENT]
    for ci, case in enumerate(coarse["cases"]):
        cfgs = case["configs"]
        for key, r in cfgs.items():
            m = "o" if "fft" in key else "^"
            ax[0].scatter(r["speedup"], r["sigma2_pct"], marker=m,
                          color=case_colors[ci % 2], s=55, alpha=0.85)
        # label the case once
        ax[0].scatter([], [], color=case_colors[ci % 2],
                      label=_short(case["name"]))
    ax[0].scatter([], [], marker="o", color="grey", label="FFT pad")
    ax[0].scatter([], [], marker="^", color="grey", label="bicubic")
    ax[0].set_xlabel("build speed-up (x)")
    ax[0].set_ylabel("aperture sigma2_I error [%]")
    ax[0].set_title("P2a: coarse screens - every saver leaves the band -> BURY")
    ax[0].legend(fontsize=8)

    # B: beam grid step 3 variants - sigma2 error with +/-10% kill band
    variants = []
    for grp in ("step3_terrestrial", "step3_space"):
        for v in beam.get(grp, {}).get("variants", []):
            variants.append((v["name"], v["sigma2_pct"], v.get("time_ratio"),
                             v.get("ops_ratio"), v["killed"]))
    # add the step-2 beam-sized-screen point
    s2 = beam.get("step2_beam_sized_screens")
    if s2:
        variants.append(("terr beam-sized", s2["sigma2_pct"],
                         1.0 / s2["speedup_build"], None, s2["killed"]))
    names = [v[0] for v in variants]
    errs = [v[1] for v in variants]
    killed = [v[4] for v in variants]
    y = range(len(names))
    colors = [KILL if k else OK for k in killed]
    ax[1].axvspan(-10, 10, color=OK, alpha=0.15, label="+/-10% accept band")
    ax[1].axvline(0, color="grey", lw=1)
    ax[1].barh(list(y), errs, color=colors)
    ax[1].set_yticks(list(y))
    ax[1].set_yticklabels(names, fontsize=8)
    ax[1].set_xlabel("aperture sigma2_I error [%]")
    ax[1].set_title("P2b: beam-following grid - savers break accuracy -> BURY")
    ax[1].legend(fontsize=8, loc="lower left")

    fig.suptitle("P2 - grid cost experiments (both buried)",
                 fontweight="bold")
    _save(fig, "grid_experiments.png")


# ---------------------------------------------------------------------------
# generator validation - the broad drop-in check
# ---------------------------------------------------------------------------
def plot_validation(d):
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))

    # A: equivalence sigmas apart per case (mean and sigma2)
    eq = d["equivalence"]
    labels = [_short(c["label"]) for c in eq]
    mp = [c["mean_power_cmp"]["sigmas_apart"] for c in eq]
    s2 = [c["sigma2_I_cmp"]["sigmas_apart"] for c in eq]
    x = range(len(eq))
    w = 0.38
    ax[0, 0].bar([i - w / 2 for i in x], mp, w, color=REF, label="mean power")
    ax[0, 0].bar([i + w / 2 for i in x], s2, w, color=OK, label="sigma2_I")
    ax[0, 0].axhline(2.0, color=KILL, ls="--", label="2 sigma line")
    ax[0, 0].set_xticks(list(x))
    ax[0, 0].set_xticklabels(labels, fontsize=8, rotation=15, ha="right")
    ax[0, 0].set_ylabel("|aotools - olb| [sigma]")
    ax[0, 0].set_title("equivalence: all agree well inside 2 sigma")
    ax[0, 0].legend(fontsize=8)

    # B: the fade tail - collected power vs quantile, both generators
    ft = d["fade_tail"]
    q = ft["quantiles"]
    xa = range(len(q))
    ax[0, 1].errorbar([i - 0.08 for i in xa], ft["aotools"]["quantile"],
                      yerr=ft["aotools"]["se_quantile"], fmt="o", color=REF,
                      capsize=4, label="aotools")
    ax[0, 1].errorbar([i + 0.08 for i in xa], ft["olb"]["quantile"],
                      yerr=ft["olb"]["se_quantile"], fmt="s", color=OK,
                      capsize=4, label="olb")
    ax[0, 1].set_xticks(list(xa))
    ax[0, 1].set_xticklabels([f"{p * 100:g}%" for p in q])
    ax[0, 1].set_xlabel("fade quantile")
    ax[0, 1].set_ylabel("collected power (vacuum = 1)")
    ax[0, 1].set_title("fade tail: the margin-setting low quantiles agree")
    ax[0, 1].legend(fontsize=8)

    # C: converged sigma2_I vs analytic
    cv = d["converged_vs_analytic"]
    labs = ["aotools", "olb"]
    vals = [cv["aotools"]["sigma2_I"], cv["olb"]["sigma2_I"]]
    ses = [cv["aotools"]["se_sigma2_I"], cv["olb"]["se_sigma2_I"]]
    ax[1, 0].bar(labs, vals, yerr=ses, capsize=6, color=[REF, OK])
    ax[1, 0].axhline(cv["sigma2_I_analytic"], color=KILL, ls="--",
                     label="analytic")
    ax[1, 0].set_ylabel("converged sigma2_I (2000 trials)")
    ax[1, 0].set_title(f"converged index vs analytic "
                       f"({cv['gen_cmp']['sigmas_apart']:.2f} sig apart)")
    ax[1, 0].legend(fontsize=8)

    # D: outer scale structure function
    os_ = d["outer_scale"]
    r = os_["r_over_r0"]
    ax[1, 1].plot(r, os_["d_theory"], "k--", label="Kolmogorov theory")
    ax[1, 1].plot(r, os_["d_olb_infinite"], "^-", color="grey",
                  label="olb, L0=inf")
    ax[1, 1].plot(r, os_["d_aotools_finite"], "o-", color=REF,
                  label="aotools, L0=25m")
    ax[1, 1].plot(r, os_["d_olb_finite"], "s-", color=OK, label="olb, L0=25m")
    ax[1, 1].set_xlabel("r / r0")
    ax[1, 1].set_ylabel("D_phi [rad^2]")
    ax[1, 1].set_title("outer scale: olb tracks aotools within 3.5%")
    ax[1, 1].legend(fontsize=8)

    fig.suptitle("Generator validation - olb is a trustworthy drop-in, "
                 "fade tail safe", fontweight="bold")
    _save(fig, "generator_validation.png")


# ---------------------------------------------------------------------------
# P3 - the parallel scaling study (drawn when present)
# ---------------------------------------------------------------------------
def plot_p3(d):
    threads = d["threads"]
    procs = d["processes"]
    batched = d["batched"]
    rec = d["recommendation"]
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))

    # A: trials/s vs workers - threads vs processes, downlink standard, both gens
    hero = "space downlink 30deg standard"
    for gen, col in (("aotools", REF), ("olb", OK)):
        key = f"{hero} | {gen}"
        tw = [r["workers"] for r in threads[key]]
        tr = [r["trials_per_s"] for r in threads[key]]
        pw = [r["workers"] for r in procs[key]]
        pr = [r["trials_per_s"] for r in procs[key]]
        ax[0, 0].plot(tw, tr, "o-", color=col, label=f"{gen} threads")
        ax[0, 0].plot(pw, pr, "s--", color=col, label=f"{gen} processes")
    ax[0, 0].set_xlabel("workers")
    ax[0, 0].set_ylabel("trials / second")
    ax[0, 0].set_title(f"scaling: {_short(hero)} (processes > threads)")
    ax[0, 0].legend(fontsize=8)
    ax[0, 0].grid(True, alpha=0.3)

    # B: parallel efficiency vs workers - all three cases, threads (olb)
    for name, col in (("terrestrial 2km standard", REF),
                      ("space downlink 30deg rapid", ACCENT),
                      ("space downlink 30deg standard", OK)):
        rows = threads[f"{name} | olb"]
        w = [r["workers"] for r in rows]
        e = [r["efficiency"] for r in rows]
        ax[0, 1].plot(w, e, "o-", color=col, label=_short(name))
    ax[0, 1].axhline(0.5, color="grey", ls=":", lw=1)
    ax[0, 1].set_xlabel("workers")
    ax[0, 1].set_ylabel("parallel efficiency (olb, threads)")
    ax[0, 1].set_title("efficiency: threads saturate by ~8-16 workers")
    ax[0, 1].legend(fontsize=8)
    ax[0, 1].grid(True, alpha=0.3)

    # C: batched FFT - trials/s vs B, olb, three backends
    rows = batched["olb"]["rows"]
    B = [r["B"] for r in rows]
    for bk, col, mk in (("numpy", REF, "o"), ("scipy_w1", ACCENT, "^"),
                        ("scipy_w32", OK, "s")):
        y = [r[bk + "_trials_per_s"] for r in rows]
        ax[1, 0].plot(B, y, mk + "-", color=col, label=bk)
    ax[1, 0].set_xlabel("batch size B")
    ax[1, 0].set_ylabel("trials / second (olb)")
    ax[1, 0].set_title("batched FFT plateaus below the process pool")
    ax[1, 0].set_xscale("log", base=2)
    ax[1, 0].set_xticks(B)
    ax[1, 0].set_xticklabels(B)
    ax[1, 0].xaxis.set_minor_formatter(plt.NullFormatter())
    ax[1, 0].legend(fontsize=8)
    ax[1, 0].grid(True, alpha=0.3)

    # D: the recommendation - best speed-up over one worker, per case|gen
    keys = sorted(rec)
    labels = [_short(k) for k in keys]
    ups = [rec[k]["speedup_over_1"] for k in keys]
    modes = [rec[k]["mode"] for k in keys]
    counts = [rec[k]["count"] for k in keys]
    y = range(len(keys))
    colors = [OK if m == "processes" else ACCENT for m in modes]
    ax[1, 1].barh(list(y), ups, color=colors)
    ax[1, 1].set_yticks(list(y))
    ax[1, 1].set_yticklabels(labels, fontsize=7)
    ax[1, 1].set_xlabel("best speed-up over 1 worker (x)")
    ax[1, 1].set_title("recommended mode per case (green=processes, "
                       "orange=threads)")
    for i, (u, m, c) in enumerate(zip(ups, modes, counts)):
        ax[1, 1].text(u + 0.1, i, f"{m[:4]} x{c}", va="center", fontsize=7)
    ax[1, 1].set_xlim(0, max(ups) * 1.25)

    fig.suptitle("P3 - parallel scaling (seed contract holds bit-for-bit "
                 "on every route)", fontweight="bold")
    _save(fig, "scaling_study.png")


# ---------------------------------------------------------------------------
def main():
    p0 = _load("profile_baseline_results.json")
    if p0:
        plot_p0(p0)
    p1 = _load("screen_generator_check_results.json")
    if p1:
        plot_p1(p1)
    coarse = _load("coarse_screen_experiment_results.json")
    beam = _load("beam_grid_experiment_results.json")
    if coarse and beam:
        plot_p2(coarse, beam)
    val = _load("generator_validation_results.json")
    if val:
        plot_validation(val)
    p3 = _load("scaling_study_results.json")
    if p3:
        plot_p3(p3)
    else:
        print("scaling_study_results.json absent; P3 plot skipped (re-run "
              "after P3 finishes).")
    print("done")


if __name__ == "__main__":
    main()

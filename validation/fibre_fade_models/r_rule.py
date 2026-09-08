"""Fit the Rice factor r alone, with sigma_z^2 pinned, and test a rule for r.

A follow-up of `fit_distributions.py` (backlog 1-9, 2026-09-08). The
two-parameter lognormal-Rician fit cannot separate r from sigma_z^2 at a
small index, so this script PINS sigma_z^2 = ln(1 + sigma2_P) to the bucket
index (the measured one and the analytic one of the fidelity-0 Term; the two
give the same r inside a few percent) and fits r alone on each fibre case.
It then fits a power law r = a x^b, with x the Gaussian-beam D/r0 or the Noll
residual phase variance sigma2_phi, and it judges the law on each cell with
the law fitted on the OTHER cells (leave one cell out). The pass rule is the
one of `fit_distributions.py`: inside 0.5 dB at 5 percent and 1.0 dB at
1 percent. Family: Andrews and Phillips, 2nd ed. (2005), Ch. 9, Eq. (133),
DOI 10.1117/3.626196.

Run it from the repository root after `fit_distributions.py`:

    python -m validation.fibre_fade_models.r_rule

It writes `r_rule.json` and `r_rule.log` (both git-ignored).
"""
import builtins
import json
import os
import sys

import numpy as np
from scipy.optimize import minimize_scalar

from olb.turbulence.andrews.distributions import lognormal_rician_pdf
from validation.fibre_fade_models import fit_distributions as fd

HERE = os.path.dirname(os.path.abspath(__file__))
_log = open(os.path.join(HERE, "r_rule.log"), "w")


def print(*args, **kw):          # noqa: A001 - the log twin of print
    text = " ".join(str(a) for a in args)
    builtins.print(text, flush=True)
    _log.write(text + chr(10))

tables = fd.load_tables()
rows = []
for tag, (cols, meta) in tables.items():
    cell = meta["cell"]
    L, cn2 = float(cell["path_length_m"]), float(cell["cn2_m_m23"])
    dz = {float(k): v for k, v in meta["dz_curv_m"].items()}
    alias = float((cols["max_step_10cm"] > fd.STEP_WARN_RAD).mean()) > fd.ALIAS_FRAC
    for D, tagD in ((0.10, "10"), (0.05, "5")):
        feeds = fd.analytic_feeds(L, cn2, D, dz[D])
        b = cols[f"bucket_{tagD}cm"]
        idx_b = float(b.var() / b.mean() ** 2)
        sz2_meas = float(np.log1p(idx_b))
        sz2_an = float(np.log1p(feeds["scint"]["sigma2_P"]))
        c = feeds["coupling"]
        dr0 = D / c["r0_gauss_m"]
        for variant in ("untracked", "tiptilt"):
            if variant == "tiptilt" and alias:
                continue
            x = b * cols[f"smf{tagD}_{variant}"]
            x = x / x.mean()
            emp = {q: float(np.quantile(-10 * np.log10(x), 1 - q)) for q in (0.05, 0.01)}
            lx = np.log(x)
            edges = np.linspace(lx.min() - 1e-9, lx.max() + 1e-9, 51)
            counts, _ = np.histogram(lx, edges)
            centres = np.exp(0.5 * (edges[1:] + edges[:-1]))
            widths = np.diff(np.exp(edges))
            keep = counts > 0

            def nll(t, sz2):
                r = float(np.exp(t))
                with np.errstate(all="ignore"):
                    p = lognormal_rician_pdf(centres[keep], r, sz2) * widths[keep]
                if not np.all(np.isfinite(p)):
                    return 1e30
                return float(-(counts[keep] * np.log(np.maximum(p, 1e-300))).sum())

            out = {"tag": tag, "D": D, "variant": variant, "dr0": dr0,
                   "s2_full": c["sigma2_res_full"], "s2_ho": c["sigma2_res_ho"],
                   "idx_bucket": idx_b, "p5": emp[0.05], "p1": emp[0.01]}
            for name, sz2 in (("meas", sz2_meas), ("an", sz2_an)):
                res = minimize_scalar(lambda t: nll(t, sz2), bounds=(-7, 9), method="bounded")
                r = float(np.exp(res.x))
                f = fd.lognormal_rician_fades(r, sz2)
                out[f"r_{name}"] = r
                out[f"d5_{name}"] = f[0.05] - emp[0.05]
                out[f"d1_{name}"] = f[0.01] - emp[0.01]
            rows.append(out)
            print(f"{tag:<26s} D={D:.2f} {variant:<9s} D/r0g={dr0:4.2f} s2full={out['s2_full']:5.2f} s2ho={out['s2_ho']:5.2f} "
                  f"idxB={idx_b:.3f} | r(sz2 meas)={out['r_meas']:7.2f} d5={out['d5_meas']:+.2f} d1={out['d1_meas']:+.2f} "
                  f"| r(sz2 an)={out['r_an']:7.2f} d5={out['d5_an']:+.2f} d1={out['d1_an']:+.2f}", flush=True)

json.dump(rows, open(os.path.join(HERE, "r_rule.json"), "w"), indent=1)
# a power law r = a (D/r0)^b on the pinned (measured sz2) fits with r > 0.05
for variant in ("untracked", "tiptilt"):
    sel = [w for w in rows if w["variant"] == variant and w["r_meas"] > 0.05]
    X = np.log([w["dr0"] for w in sel]); Y = np.log([w["r_meas"] for w in sel])
    b, a = np.polyfit(X, Y, 1)
    resid = Y - (a + b * X)
    print(f"\n{variant}: ln r = {a:.2f} + {b:.2f} ln(D/r0_gauss) over {len(sel)} points, rms scatter {resid.std():.2f} in ln r")
    # and against the residual phase variance
    key = "s2_full" if variant == "untracked" else "s2_ho"
    X2 = np.log([w[key] for w in sel])
    b2, a2 = np.polyfit(X2, Y, 1); resid2 = Y - (a2 + b2 * X2)
    print(f"{variant}: ln r = {a2:.2f} + {b2:.2f} ln(sigma2_phi) rms scatter {resid2.std():.2f}")


# --- the leave-one-cell-out test -------------------------------------
# the analytic sigma_z^2 of each (tag, D)
sz2_an = {}
for tag, (cols, meta) in tables.items():
    cell = meta["cell"]
    L, cn2 = float(cell["path_length_m"]), float(cell["cn2_m_m23"])
    dz = {float(k): v for k, v in meta["dz_curv_m"].items()}
    for D in (0.10, 0.05):
        f = fd.analytic_feeds(L, cn2, D, dz[D])
        sz2_an[(tag, D)] = float(np.log1p(f["scint"]["sigma2_P"]))

def cell_of(tag):
    return tag.replace("_rapid", "").replace("_standard", "")

for variant, key in (("untracked", "s2_full"), ("tiptilt", "s2_ho")):
    sel = [w for w in rows if w["variant"] == variant]
    print(f"\n=== {variant}: rule r = a * x^b fitted on the OTHER cells (x = D/r0_gauss, and x = sigma2_phi), r > 0.05 points only")
    print(f"{'cell':<26s} {'D':>4s} {'D/r0':>5s} {'s2phi':>6s} {'r fit':>7s} | {'r(D/r0)':>8s} {'d5':>6s} {'d1':>6s} | {'r(s2)':>8s} {'d5':>6s} {'d1':>6s} | emp p5/p1")
    hold = {"dr0": 0, "s2": 0, "n": 0}
    for w in sel:
        train = [v for v in sel if cell_of(v["tag"]) != cell_of(w["tag"]) and v["r_meas"] > 0.05]
        out = [w["tag"], f"{w['D']:.2f}", f"{w['dr0']:5.2f}", f"{w[key]:6.2f}", f"{w['r_meas']:7.2f}"]
        sz2 = sz2_an[(w["tag"], w["D"])]
        for xkey, name in (("dr0", "dr0"), (key, "s2")):
            X = np.log([v[xkey] for v in train]); Y = np.log([v["r_meas"] for v in train])
            b, a = np.polyfit(X, Y, 1)
            r = float(np.exp(a + b * np.log(w[xkey])))
            f = fd.lognormal_rician_fades(r, sz2)
            d5, d1 = f[0.05] - w["p5"], f[0.01] - w["p1"]
            out += [f"{r:8.2f}", f"{d5:+6.2f}", f"{d1:+6.2f}"]
            if abs(d5) <= 0.5 and abs(d1) <= 1.0:
                hold[name] += 1
        hold["n"] += 1
        out += [f"{w['p5']:.1f}/{w['p1']:.1f}"]
        print(f"{out[0]:<26s} {out[1]:>4s} {out[2]:>5s} {out[3]:>6s} {out[4]:>7s} | {out[5]} {out[6]} {out[7]} | {out[8]} {out[9]} {out[10]} | {out[11]}")
    print(f"holds (|d5|<=0.5, |d1|<=1.0): D/r0 rule {hold['dr0']}/{hold['n']}, sigma2_phi rule {hold['s2']}/{hold['n']}")

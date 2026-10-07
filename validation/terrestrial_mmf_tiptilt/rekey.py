"""Re-key the stored campaigns after the MMF `optimal_focus` removal (2026-10-07).

The campaigns of this study ran with `MMF(optimal_focus=True)`, which resolved
to f = F_MMF. The scenario now gives `focal_length_m=F_MMF` explicitly, so
every trial is the SAME (the detector sets only the in-run scalar mmf_eta, and
that f is unchanged), but the scenario text and the fingerprint moved. This
script makes each campaign again in a scratch root with the new scenario,
takes its fingerprint and scenario text, and writes them into the stored
manifest. A manifest that holds no `optimal_focus=True` is left alone.

    python -m validation.terrestrial_mmf_tiptilt.rekey
"""

import glob
import json
import os
import shutil
import tempfile
import warnings

from validation.terrestrial_mmf_tiptilt import run, smoke

ROOT = os.path.join(run.HERE, "campaigns")


def twin(name, scratch):
    """Make the campaign of a stored root name again under scratch."""
    run.HERE = smoke.HERE = scratch           # make_campaign joins HERE
    os.makedirs(os.path.join(scratch, "campaigns"), exist_ok=True)
    if name.startswith("smoke_"):
        cn2, _, n = name[len("smoke_cn2"):].partition("_n")
        return smoke.smoke_campaign(float(cn2), 10, "cupy",
                                    int(n) if n else None)
    cn2, rest = name[len("cn2"):].split("_", 1)
    n = int(rest.rsplit("_n", 1)[1]) if "_n" in rest else None
    return run.make_campaign(float(cn2), "standard", 50, "cupy", n)


def main():
    here = run.HERE
    for man_path in sorted(glob.glob(os.path.join(ROOT, "*", "manifest.json"))):
        name = os.path.basename(os.path.dirname(man_path))
        with open(man_path) as f:
            man = json.load(f)
        if "optimal_focus=True" not in man["scenario"]:
            print(f"{name}: no optimal_focus=True, left alone")
            continue
        if not glob.glob(os.path.join(os.path.dirname(man_path), "block_*.npz")):
            print(f"{name}: no stored trial (a dry run), left alone")
            continue
        scratch = tempfile.mkdtemp()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                new = twin(name, scratch)
            with open(os.path.join(new.root_dir, "manifest.json")) as f:
                new_man = json.load(f)
        finally:
            run.HERE = smoke.HERE = here
            shutil.rmtree(scratch, ignore_errors=True)
        for k in ("grid", "plan", "patch_radius_m", "block_size", "seed"):
            assert new_man[k] == man[k], (name, k)
        man["fingerprint"] = new_man["fingerprint"]
        man["scenario"] = new_man["scenario"]
        with open(man_path, "w") as f:
            json.dump(man, f, indent=1)
        print(f"{name}: re-keyed")


if __name__ == "__main__":
    main()

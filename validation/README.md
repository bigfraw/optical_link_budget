# olb validation scripts

This folder holds the owner's cross-check and validation scripts. They are not
curated user examples. Each one checks one model against another model, or
against a known result.

These scripts can be specific, they can overlap, and they can be rough. Read
[../examples/](../examples/) first if you want the curated set.

Run each script from the repository root as a module:

    python -m validation.uplink_divergence

## The scripts

| File | Purpose |
| --- | --- |
| [uplink_divergence.py](uplink_divergence.py) | A trade study. It widens the uplink transmit beam on purpose, then it finds the divergence with the best 99% margin. |
| [terrestrial_coupling_jitter.py](terrestrial_coupling_jitter.py) | It splits the terrestrial single-mode-fibre coupling loss into three pointing mechanisms, then it sweeps each one. |
| [terrestrial_mmf_na.py](terrestrial_mmf_na.py) | It shows the numerical-aperture angular gate of a terrestrial multimode-fibre link, then it sweeps the focal length. |
| [mmf_coupling_validation.py](mmf_coupling_validation.py) | It plots the multimode-fibre coupled power against the incident angle, for the correct encircled-energy model and for the old, wrong Gaussian roll-off. It saves a PNG. |

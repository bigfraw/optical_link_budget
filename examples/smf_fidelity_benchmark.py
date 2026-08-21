'''
Benchmark the SMF coupling models: reciprocity Strehl proxy vs FAST fidelity-1.

The downlink single-mode-fibre coupling has two turbulence models in olb:

- "reciprocity": a Strehl proxy. eta = eta_max * on-axis Strehl, where the Strehl
  comes from the Dios coupled-flux of the back-projected fibre mode. It captures
  the tip-tilt / angle-of-arrival fade but is NOT a true modal overlap.
- "fast": fidelity-1. The true LP01 Gaussian-mode overlap under turbulence, from
  FAST (the fast-aosim package). This is the reference.

This script runs a no-AO fibre downlink over a set of elevations and prints the
mean and 99% coupling loss from each model, so you can see where the proxy agrees
with the reference and where it does not. The reference needs fast-aosim; install
it with `pip install fast-aosim`.

Run from the repo root:
    python -m examples.smf_fidelity_benchmark
'''

import warnings
warnings.filterwarnings("ignore", category=UserWarning)

import numpy as np

from olb import (Scenario, Channel, Site, CircularOrbit, Terminal, Transmitter,
                 SMF)
from olb.models.coupling import rx_coupling_term


def main():
    wavelength_m = 1550e-9

    # A no-AO fibre ground receiver. The satellite transmits the downlink.
    ground = Terminal(aperture_m=0.7, obscuration_ratio=0.3, wavelength_m=wavelength_m,
                      detector=SMF(sensitivity_dbm=-110.0))
    space = Terminal(aperture_m=0.08, wavelength_m=wavelength_m,
                     transmitter=Transmitter(waist_m=0.05, power_dbm=30.0))
    scenario = Scenario(ground=ground, space=space, direction="downlink",
                        channel=Channel(site=Site(cn2_ground=1.7e-14),
                                        altitude_m=1500e3))

    rng = np.random.default_rng(0)
    print(f"{'elev':>5} | {'reciprocity mean':>16} {'99%':>7} | "
          f"{'FAST mean':>10} {'99%':>7}")
    print("-" * 56)
    for elevation_deg in (30.0, 45.0, 60.0, 90.0):
        geom = CircularOrbit(scenario.channel.altitude_m,
                             elevation_deg=elevation_deg)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            recip = rx_coupling_term(scenario, geom, n_samples=2000,
                                     smf_fidelity="reciprocity")
            fast = rx_coupling_term(scenario, geom, n_samples=2000,
                                    smf_fidelity="fast")
        # The reciprocity term is Monte-Carlo-only; use samples for the 99%.
        recip_99 = float(np.percentile(recip.sample_db(8000, rng), 99))
        print(f"{elevation_deg:5.0f} | {recip.mean_db:16.2f} {recip_99:7.2f} | "
              f"{fast.mean_db:10.2f} {fast.quantile_db(0.99):7.2f}")

    print("\nThe reciprocity Strehl proxy tracks the FAST modal overlap on the "
          "MEAN, but its deep-fade tail (99%) is far heavier: the Dios on-axis "
          "intensity saturates where the FAST aperture-integrated overlap stays "
          "bounded. So the proxy is a fair mean estimate but pessimistic on the "
          "tail. FAST is the reference (with subharmonics, so the tilt is "
          "captured); the proxy is the cheap, wavefront-free estimate. Neither "
          "models point-ahead here.")


if __name__ == "__main__":
    main()

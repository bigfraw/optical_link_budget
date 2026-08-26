# olb documentation

This directory documents `optical_link_budget` (olb) in its current form. The
package builds optical (laser) link budgets with atmospheric propagation, fade
statistics, and Monte Carlo. It models uplink, downlink, retroreflected, and
horizontal terrestrial links.

All documentation uses ASD-STE100 Simplified Technical English. See
[CONVENTIONS.md](../CONVENTIONS.md).

## Contents

1. [Architecture](architecture.md) — the data model, the one-way dependency, and
   the control flow from a scenario to a Monte Carlo result.
2. [Getting started](getting-started.md) — install, the Quickstart, and a
   walkthrough of the four link families.
3. [API — terminals, scenarios, and geometry](api-terminal-scenario.md) — the
   pure-data classes that describe the hardware and the link case.
4. [API — Terms, Budgets, and link entry points](api-budget.md) — the result
   objects, the assumption checks, and the per-link budget functions.
5. [Physics reference](physics.md) — the models and the turbulence kernels, with
   a source DOI for each equation.
6. [Examples](examples.md) — a guide to each runnable script in `examples/`. It
   also points to the Andrews foundation suite in `examples/andrews/` and to the
   wave-optics suite in `examples/waveoptics/`.
7. [Andrews cross-check](andrews-crosscheck.md) — the running comparison of the
   olb equations against Andrews and Phillips, 2nd ed. (2005), with a status per
   entry and the adjudicated conflicts. This is the internal physics record.

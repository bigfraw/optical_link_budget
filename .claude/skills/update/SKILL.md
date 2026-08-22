---
name: update
description: >-
  Reconcile this repo's documentation with the code. Review the changes made
  (working tree plus recent commits) and bring every stale doc back in sync:
  the files in docs/, the main README.md, and the CLAUDE.md roadmap and
  architecture notes. Use this skill whenever the user says "/update", "update
  the docs", "sync the docs", "are the docs current", or after any change to
  the public API, the physics, the module layout, the examples, or the roadmap
  status — even when the user does not name a specific doc file.
---

# Update the documentation

Bring the documentation back in sync with the code as it is now. Do not invent
behavior. Document what the code does, verified against the source.

## 1. Find the delta

Work out what changed since the docs were last correct. Read the real state:

- `git status --short` and `git diff` — the uncommitted working tree.
- `git log --oneline -15` — the recent commits.
- Compare the code against what each doc claims. The doc is stale, not the code.

State the delta in one short list before you edit anything.

## 2. Map each change to its doc

Each code area has one home doc. Update the home doc when the area changes.

| Code area that changed | Doc(s) to update |
|---|---|
| Module layout, one-way dependency, data flow | `docs/architecture.md`, `CLAUDE.md` architecture section |
| Public API: a signature, a default, a new/removed export in `olb/__init__.py` | `docs/api-terminal-scenario.md`, `docs/api-budget.md`, `docs/getting-started.md` |
| A Terminal part, a Scenario, a Channel, a geometry class | `docs/api-terminal-scenario.md` |
| A Term, a Budget method, a `*_budget` entry point | `docs/api-budget.md` |
| An equation, a DOI, a turbulence kernel, a regime limit | `docs/physics.md` |
| A new, changed, or removed script in `examples/` | `docs/examples.md` |
| A link type or a fidelity tier that moved from planned to done | `README.md` roadmap (the mermaid trees), `CLAUDE.md` "Next task" |
| Install, dependency, or the Quickstart flow | `README.md`, `docs/getting-started.md` |

`docs/README.md` is the index. Update it only when a doc file is added or removed.

## 3. Verify, do not trust the prose

The docs were built by reading the source, so check the same way:

- Read the changed source. Copy signatures, keyword defaults, and field values
  from the code, not from memory. The api docs list exact defaults; keep them
  exact.
- For an equation, keep the DOI the code cites. Do not add a DOI the code does
  not give. If the code has no DOI, write "(source cited in <file>)".
- Run the affected module self-check if a value is in doubt:
  `python -m olb.<module>`.

## 4. Keep the conventions

- Write all prose in ASD-STE100 Simplified Technical English. See
  `CONVENTIONS.md`. Short sentences, active voice, present tense, verb-first
  instructions, keep the articles, no -ing verb form, no figurative language.
- Loss is positive dB. Gain is negative dB.
- Keep the existing flags. A flagged limitation (for example the collimated-beam
  Gaussian Fried note in `docs/physics.md`) stays until the code removes the
  limit. Do not delete a flag because the surrounding text moved.

## 5. Delegate the heavy edits (optional)

For a large rewrite of one doc, delegate it to an Opus subagent, one doc per
agent, the way the docs were first written. Give each agent the STE-100 block,
the exact source files to read, and the one doc file to edit. Independent docs
have no shared state, so run the agents in parallel. For a small edit, do it
inline.

## 6. Report

End with a short report:

- What was stale, and what you changed, one line each.
- What needs a human decision, not a silent edit: a roadmap status change (done
  vs planned), a new deliberate deferral, or a claim the code and the docs now
  disagree on. Ask before you flip a roadmap leaf from planned to done.

Do not commit. Leave the staged/unstaged state for the user, unless the user
asks you to commit.

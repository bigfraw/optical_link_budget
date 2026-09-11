"""SessionStart hook: a short repo briefing.

Prints the branch, the last commit, the dirty files, and a warning when the
branch is not ``main`` or when stray untracked sim outputs (``.npz`` / ``.log``)
are sitting in the tree. Claude Code adds this stdout to the session context, so
the agent starts each session aware of the repo state (this repo juggles several
feature branches and its remote desktop can leave sim junk behind).
"""

import subprocess
import sys


def _git(*args: str) -> str:
    try:
        out = subprocess.run(
            ["git", *args], capture_output=True, text=True, timeout=15
        )
        return out.stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return ""


def main() -> int:
    branch = _git("rev-parse", "--abbrev-ref", "HEAD")
    if not branch:
        return 0  # not a git repo / git unavailable; stay silent

    last = _git("log", "-1", "--oneline")
    porcelain = _git("status", "--porcelain")
    lines = [ln for ln in porcelain.splitlines() if ln.strip()]
    dirty = len(lines)
    junk = [
        ln[3:]
        for ln in lines
        if ln.startswith("??") and ln.rstrip().endswith((".npz", ".log"))
    ]

    print("=== olb session briefing ===")
    print(f"branch : {branch}")
    print(f"last   : {last}")
    print(f"dirty  : {dirty} changed/untracked file(s)")
    if branch != "main":
        print(f"NOTE   : on feature branch '{branch}', not main.")
    if junk:
        shown = ", ".join(junk[:5]) + (" ..." if len(junk) > 5 else "")
        print(f"WARNING: {len(junk)} untracked sim output(s) present: {shown}")
        print("         (do NOT force-add these; they are meant to stay ignored)")
    print("============================")
    return 0


if __name__ == "__main__":
    sys.exit(main())

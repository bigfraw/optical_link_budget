"""PostToolUse hook: run a module self-check after Claude Code edits it.

Claude Code calls this after an Edit/Write. It reads the hook payload as JSON on
stdin, finds the edited file, maps it to its olb module, and runs the module's
``if __name__ == '__main__':`` self-check with the olb conda interpreter.

Behaviour:
  - Acts only on ``*.py`` files inside the ``olb`` package.
  - Skips package ``__init__.py`` and any file with no ``__main__`` block.
  - Skips the SKIP_PREFIXES subtrees (heavy Monte-Carlo self-checks).
  - Exit 0 (silent) when the self-check passes or the file is skipped.
  - Exit 2 with the captured output when the self-check fails or times out;
    Claude Code shows that output back to the agent.
"""

import json
import subprocess
import sys
from pathlib import Path

# The olb conda interpreter (not on PATH in tool shells).
PYTHON = r"C:\Users\alexf\anaconda3\envs\olb\python.exe"

# Module prefixes whose self-checks run long sims; skip them.
SKIP_PREFIXES = ("olb.waveoptics",)

# Kill a self-check that runs longer than this (seconds).
TIMEOUT_S = 90


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0  # nothing to act on

    file_path = (payload.get("tool_input") or {}).get("file_path")
    if not file_path:
        return 0

    path = Path(file_path)

    # A README is not code, so there is no self-check to run -- but announce it
    # loudly, because this repo keeps its READMEs in sync with the code.
    if path.stem.upper() == "README":
        banner = (
            "\n"
            "!!!==============================================================!!!\n"
            "!!!  README EDITED  --  no self-check (not a code module)         !!!\n"
            f"!!!  file: {path.name:<49}!!!\n"
            "!!!  Keep the docs in sync with the code (the /update skill).    !!!\n"
            "!!!==============================================================!!!\n"
        )
        print(banner, file=sys.stderr)
        return 2

    if path.suffix != ".py" or path.name == "__init__.py":
        return 0

    # Find the repo root: the first parent that holds the ``olb`` package.
    root = None
    for parent in path.parents:
        if (parent / "olb" / "__init__.py").exists():
            root = parent
            break
    if root is None:
        return 0

    try:
        rel = path.resolve().relative_to(root.resolve())
    except ValueError:
        return 0
    if rel.parts[0] != "olb":
        return 0

    module = ".".join(rel.with_suffix("").parts)
    if any(module == p or module.startswith(p + ".") for p in SKIP_PREFIXES):
        print(f"[self-check] skipped {module} (heavy subtree)")
        return 0

    # Only run modules that actually have a self-check.
    if "__main__" not in path.read_text(encoding="utf-8", errors="ignore"):
        return 0

    try:
        proc = subprocess.run(
            [PYTHON, "-m", module],
            cwd=str(root),
            capture_output=True,
            text=True,
            timeout=TIMEOUT_S,
        )
    except subprocess.TimeoutExpired:
        print(f"[self-check] {module} timed out after {TIMEOUT_S}s", file=sys.stderr)
        return 2

    if proc.returncode == 0:
        print(f"[self-check] {module} OK")
        return 0

    print(f"[self-check] {module} FAILED (exit {proc.returncode})", file=sys.stderr)
    tail = (proc.stdout or "") + (proc.stderr or "")
    print(tail[-3000:], file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())

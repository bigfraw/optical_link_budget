"""PreToolUse hook: enforce this repo's standing Bash rules before a command runs.

Each rule below is a house rule the owner has stated explicitly. The hook inspects
a Bash command and blocks (exit 2) a violation with a clear fix; Claude Code shows
the message back to the agent. It never blocks a command that breaks no rule.

Rules:
  1. No ``git add -f/--force`` of generated outputs (.log/.json/.npz/.csv) or a
     broad target (. / -A / *) that would sweep them in.
  2. No ``&&`` / ``||`` inside an ``ssh desktop "..."`` command: the remote shell
     is Windows PowerShell 5.1, where those are a parse error. Chain with ``;``.
  3. No destructive ``git clean -f/-d/-x`` (the sim .npz data is untracked, so a
     clean would delete it). A ``git clean -n`` dry run is allowed.
  4. No bare ``python`` / ``python3`` / ``pip``: they resolve to the Windows Store
     shim in tool shells. Use the env's absolute interpreter.
  5. No bare ``gh``: the GitHub CLI is not on PATH in tool shells. Use its
     absolute path.
  6. No ``scp -r`` of campaign / .npz data: hundreds of tiny files are latency
     bound. Zip on the source, scp the one archive, expand on the far side.
"""

import json
import re
import sys

OLB_PYTHON = r"C:\Users\alexf\anaconda3\envs\olb\python.exe"
GH_EXE = r"C:\Program Files\GitHub CLI\gh.exe"

FORBIDDEN_EXT = (".log", ".json", ".npz", ".csv")
BROAD_TARGETS = {".", "./", "-A", "--all", "*", ":/"}
BARE_PY = {"python", "python3", "pip", "pip3"}


def _is_force(token: str) -> bool:
    if token == "--force":
        return True
    return bool(re.fullmatch(r"-[a-zA-Z]*f[a-zA-Z]*", token))


def _first_command(segment: str) -> str:
    """The executable a segment runs, skipping VAR=val prefixes and wrappers."""
    for token in segment.split():
        if "=" in token and not token.startswith("-"):
            continue  # VAR=value assignment prefix
        if token in ("time", "sudo", "env", "nohup", "exec", "command"):
            continue
        return token.strip("\"'")
    return ""


def _git_add_violation(segment: str) -> str | None:
    tokens = segment.split()
    if not any(_is_force(t) for t in tokens):
        return None
    low = segment.lower()
    for ext in FORBIDDEN_EXT:
        if ext in low:
            return (
                f"a forced add referencing {ext} files.\n"
                "House rule: never force-add generated .log/.json/.npz/.csv "
                "outputs (they are gitignored on purpose). Stage scripts, docs "
                "and figures only."
            )
    if any(t in BROAD_TARGETS for t in tokens):
        return (
            "a forced add of a broad target that could sweep in outputs.\n"
            "House rule: never force-add generated outputs. Add specific "
            "tracked files without -f."
        )
    return None


def _git_clean_violation(segment: str) -> str | None:
    if not re.search(r"\bgit\s+clean\b", segment):
        return None
    tokens = segment.split()
    if any(t in ("-n", "--dry-run") for t in tokens):
        return None  # dry run is safe
    destructive = any(
        t == "--force" or bool(re.fullmatch(r"-[a-zA-Z]*[fdx][a-zA-Z]*", t))
        for t in tokens
    )
    if destructive:
        return (
            "a destructive 'git clean'.\n"
            "House rule: never git clean this tree -- the sim .npz results under "
            "validation/**/campaigns/ are UNTRACKED and a clean would delete "
            "them. Remove specific files explicitly if that is truly intended, "
            "or use 'git clean -n' to preview."
        )
    return None


def _scp_violation(segment: str) -> str | None:
    if "scp" not in segment:
        return None
    tokens = segment.split()
    if not any(re.fullmatch(r"-[a-zA-Z]*r[a-zA-Z]*", t) for t in tokens):
        return None
    low = segment.lower()
    if "campaign" in low or ".npz" in low:
        return (
            "an 'scp -r' of campaign / .npz data.\n"
            "House rule: a campaign is hundreds of tiny per-block .npz files; "
            "scp -r is per-file-latency bound. Zip on the source "
            "(Compress-Archive), scp the single archive, then expand on the far "
            "side."
        )
    return None


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0

    if payload.get("tool_name") not in (None, "Bash"):
        return 0
    command = (payload.get("tool_input") or {}).get("command", "")
    if not command.strip():
        return 0

    # Whole-command rule: ssh desktop with a PowerShell-illegal separator.
    if re.search(r"\bssh\s+desktop\b", command) and (
        "&&" in command or "||" in command
    ):
        print(
            "BLOCKED: '&&'/'||' inside an ssh desktop command.\n"
            "The remote shell on bigfraw is Windows PowerShell 5.1, where "
            "'&&' is a parse error. Chain steps with ';' instead.",
            file=sys.stderr,
        )
        return 2

    # Per-segment rules.
    for segment in re.split(r"[;&|\n]+", command):
        seg = segment.strip()
        if not seg:
            continue

        if re.search(r"\bgit\s+add\b", seg):
            reason = _git_add_violation(seg)
            if reason:
                print(f"BLOCKED: {reason}", file=sys.stderr)
                return 2

        reason = _git_clean_violation(seg)
        if reason:
            print(f"BLOCKED: {reason}", file=sys.stderr)
            return 2

        reason = _scp_violation(seg)
        if reason:
            print(f"BLOCKED: {reason}", file=sys.stderr)
            return 2

        cmd = _first_command(seg)
        if cmd in BARE_PY:
            print(
                f"BLOCKED: bare '{cmd}'.\n"
                "It resolves to the Windows Store shim in tool shells and errors. "
                f'Use the env interpreter: "{OLB_PYTHON}" '
                f'(for pip: "{OLB_PYTHON}" -m pip ...).',
                file=sys.stderr,
            )
            return 2
        if cmd == "gh":
            print(
                "BLOCKED: bare 'gh'.\n"
                "The GitHub CLI is not on PATH in tool shells. Use its absolute "
                f'path: "{GH_EXE}"',
                file=sys.stderr,
            )
            return 2

    return 0


if __name__ == "__main__":
    sys.exit(main())

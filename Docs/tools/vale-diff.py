#!/usr/bin/env python3
"""Report Vale alerts only on lines this changeset touched, and say which
rules actually bind.

Changed *files* is too coarse: a one-line edit to a long page surfaces all of
that page's pre-existing alerts, which is what makes people ignore a gate.
This filters to changed lines, so adoption costs nothing retroactively.

"Binding" means at or above --gate-level. Anything below it is advisory and
never affects the exit code, so the report separates the two rather than
leaving the reader to infer it from severities.

Usage:
    vale-diff.py <base-ref> [--config DIR] [--gate-level error] [--strict]
"""
import argparse
import collections
import json
import re
import subprocess
import sys
from pathlib import Path

HUNK = re.compile(r"^@@ -\d+(?:,\d+)? \+(\d+)(?:,(\d+))? @@")
RANK = {"suggestion": 0, "warning": 1, "error": 2}


def changed_lines(base: str, patterns: list[str]) -> dict[str, set[int]]:
    """Map path -> line numbers added or modified relative to base."""
    out = subprocess.run(
        ["git", "diff", "-U0", f"{base}...HEAD", "--", *patterns],
        capture_output=True, text=True, check=True,
    ).stdout
    result: dict[str, set[int]] = collections.defaultdict(set)
    path = None
    for line in out.split("\n"):
        if line.startswith("+++ b/"):
            path = line[6:]
        elif line.startswith("@@") and path:
            m = HUNK.match(line)
            if m:
                start, count = int(m.group(1)), int(m.group(2) or 1)
                result[path].update(range(start, start + count))
    return {p: s for p, s in result.items() if s}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("base")
    ap.add_argument("--config", default=".", help="directory holding .vale.ini")
    ap.add_argument("--vale", default="vale")
    ap.add_argument("--gate-level", default="error", choices=list(RANK),
                    help="severity at or above which an alert binds")
    ap.add_argument("--strict", action="store_true",
                    help="exit 1 if any binding alert lands on a changed line")
    ap.add_argument("--glob", nargs="*",
                    default=["*.rst", "*.md", "*.H", "*.h", "*.cpp"])
    args = ap.parse_args()
    gate = RANK[args.gate_level]

    touched = changed_lines(args.base, args.glob)
    if not touched:
        print("no linted files changed")
        return 0

    # Vale runs with cwd=config so it finds .vale.ini and StylesPath, so the
    # paths handed to it must be absolute -- relative ones silently lint
    # nothing and the gate reports a false pass.
    root = subprocess.run(["git", "rev-parse", "--show-toplevel"],
                          capture_output=True, text=True, check=True).stdout.strip()
    abs_paths = [str(Path(root) / p) for p in touched]
    if missing := [p for p in abs_paths if not Path(p).is_file()]:
        print(f"error: {len(missing)} changed path(s) not on disk, "
              f"e.g. {missing[0]}", file=sys.stderr)
        return 2

    proc = subprocess.run(
        [args.vale, "--no-exit", "--minAlertLevel=suggestion",
         "--output=JSON", *abs_paths],
        cwd=args.config, capture_output=True, text=True,
    )
    try:
        data = json.loads(proc.stdout or "{}")
    except json.JSONDecodeError:
        print(proc.stderr.strip() or "vale produced no parsable output",
              file=sys.stderr)
        return 2

    new, pre = [], 0
    for path, alerts in data.items():
        key = next((k for k in touched if path.endswith(k)), None)
        if key is None:
            continue
        for a in alerts:
            if a.get("Line") in touched[key]:
                new.append((key, a))
            else:
                pre += 1

    for key, a in sorted(new, key=lambda t: (t[0], t[1].get("Line", 0))):
        sev = a.get("Severity", "?")
        mark = "*" if RANK.get(sev, 0) >= gate else " "
        print(f"{mark} {key}:{a.get('Line')}:{a.get('Span', [0])[0]}: "
              f"[{sev:10}] {a.get('Check'):32} {a.get('Message')}")

    binding = collections.Counter(
        a.get("Check") for _k, a in new
        if RANK.get(a.get("Severity", "suggestion"), 0) >= gate)
    advisory = collections.Counter(
        a.get("Check") for _k, a in new
        if RANK.get(a.get("Severity", "suggestion"), 0) < gate)

    if binding:
        print(f"\nBINDING (>= {args.gate_level}) -- these fail the gate:")
        for chk, n in binding.most_common():
            print(f"  {n:5}  {chk}")
    else:
        print(f"\nBINDING (>= {args.gate_level}): none")
    if advisory:
        print("ADVISORY (below gate, informational only):")
        for chk, n in advisory.most_common(10):
            print(f"  {n:5}  {chk}")

    lines = sum(len(v) for v in touched.values())
    print(f"\n{len(new)} alert(s) on {lines} changed line(s) across "
          f"{len(touched)} file(s): {sum(binding.values())} binding, "
          f"{sum(advisory.values())} advisory; {pre} pre-existing ignored")
    return 1 if (args.strict and binding) else 0


if __name__ == "__main__":
    sys.exit(main())

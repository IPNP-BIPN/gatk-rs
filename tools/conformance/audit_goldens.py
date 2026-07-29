#!/usr/bin/env python3
"""Check that every committed golden belongs to a declared conformance suite.

A golden that no suite regenerates has never been compared to anything but itself. Before the
manifest existed, seventeen of them were in that state: they were produced once, committed, and
read only by the Rust tests, so the Rust tests were asserting that the port reproduces a file of
unknown provenance. That is a strictly weaker claim than the one the README makes, and nothing in
the repository said which goldens were in which category.

This makes the condition mechanical. Only committed files count (`git ls-files`), so a golden
being written by an in-flight slice does not fail the audit before its suite lands.
"""

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import compare as comparator  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
DATA = "crates/*/tests/data"


def committed_goldens():
    out = subprocess.run(
        ["git", "ls-files", "crates/*/tests/data/*.txt.gz"],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=True,
    )
    return sorted(line for line in out.stdout.splitlines() if line)


def main():
    manifest = comparator.load_manifest()
    # A `golden-pending` case declares `golden: null` on purpose: its dump has never been produced
    # by CI, and committing one from a developer machine is what decision 0008 is about. It is
    # counted, not silently skipped.
    declared = {
        case["golden"]
        for suite in manifest["suites"]
        for case in suite["cases"]
        if case.get("golden")
    }
    pending = [
        (suite["id"], case["dump"])
        for suite in manifest["suites"]
        for case in suite["cases"]
        if not case.get("golden")
    ]
    committed = committed_goldens()

    undeclared = [g for g in committed if g not in declared]
    missing = sorted(g for g in declared if not (REPO / g).exists())

    for golden in undeclared:
        print(f"FAIL undeclared golden: {golden}")
        print("     add a suite to tools/conformance/manifest.json that regenerates it")
    for golden in missing:
        print(f"FAIL declared golden does not exist: {golden}")

    by_status = {}
    for suite in manifest["suites"]:
        by_status[suite["status"]] = by_status.get(suite["status"], 0) + len(suite["cases"])
    summary = " ".join(f"{status}={count}" for status, count in sorted(by_status.items()))
    print(f"goldens committed={len(committed)} declared={len(declared)} ({summary})")

    return 1 if (undeclared or missing) else 0


if __name__ == "__main__":
    sys.exit(main())

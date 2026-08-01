#!/usr/bin/env python3
"""Measure how well the Rust sources explain themselves, and stop finished files from regressing.

The programme claims byte-identical output, and a reader checking that claim has to compare the
Java with the Rust. Not all of them read Rust. `docs/COMMENTING.md` is the standard that follows
from that; this script is the part of it that can be automated.

It reports a ratio per file: comment lines over code lines. It does NOT fail a file for being low,
because a generated table of constants should not be commented line by line and a mechanical
accessor needs nothing. What it fails is regression: a file listed in `commented.txt` was brought
to the standard deliberately, and its recorded floor may not drop.

Usage:
    comment_density.py            report every file, sorted by ratio
    comment_density.py --check    fail if any listed file fell below its recorded floor
    comment_density.py --record   rewrite the floors from the current tree (review the diff)
"""

from __future__ import annotations

import argparse
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
CRATES = ROOT / "crates"
LISTED = pathlib.Path(__file__).resolve().parent / "commented.txt"


def measure(path: pathlib.Path) -> tuple[int, int]:
    """Return (comment lines, code lines) for one Rust source.

    A comment line is one whose first non-whitespace characters are `//` (which covers `///` and
    `//!`). A code line is any other non-blank line. Block comments are not counted, because this
    codebase does not use them and pretending otherwise would invite a way to game the number.

    String literals containing `//` would be miscounted, which is accepted: the ratio is a
    reporting aid, not a proof, and the only file where it decides anything is one a human has
    already read.
    """
    comments = 0
    code = 0
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("//"):
            comments += 1
        else:
            code += 1
    return comments, code


def sources() -> list[pathlib.Path]:
    """Every Rust source under `crates/`, excluding build output, sorted for a stable report."""
    return sorted(
        p
        for p in CRATES.rglob("*.rs")
        if "target" not in p.parts
    )


def ratio(comments: int, code: int) -> float:
    """Comments per line of code. A file with no code at all counts as fully commented."""
    return comments / code if code else float("inf")


def load_floors() -> dict[str, float]:
    """The recorded floor for each file that has been brought to the standard.

    Format is one `<path> <ratio>` pair per line, `#` starting a comment. Paths are relative to the
    repository root so the file reads as a checklist.
    """
    floors: dict[str, float] = {}
    if not LISTED.exists():
        return floors
    for raw in LISTED.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        name, _, value = line.partition(" ")
        floors[name.strip()] = float(value.strip())
    return floors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="fail on regression")
    parser.add_argument("--record", action="store_true", help="rewrite the floors")
    args = parser.parse_args()

    measured: dict[str, tuple[int, int, float]] = {}
    for path in sources():
        comments, code = measure(path)
        name = str(path.relative_to(ROOT))
        measured[name] = (comments, code, ratio(comments, code))

    floors = load_floors()

    if args.record:
        lines = [
            "# Files brought to the standard in docs/COMMENTING.md, with the comment-to-code",
            "# ratio each had when it was done. tools/audit/comment_density.py --check fails if",
            "# one of them drops below its recorded floor, so the work is a ratchet: nothing",
            "# already explained can quietly come undone.",
            "#",
            "# Regenerate with: python3 tools/audit/comment_density.py --record",
            "",
        ]
        for name in sorted(floors):
            if name not in measured:
                print(f"listed file no longer exists: {name}", file=sys.stderr)
                return 1
            lines.append(f"{name} {measured[name][2]:.3f}")
        LISTED.write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(f"recorded {len(floors)} floors")
        return 0

    if args.check:
        failed = False
        for name, floor in sorted(floors.items()):
            if name not in measured:
                print(f"listed file no longer exists: {name}")
                failed = True
                continue
            comments, code, current = measured[name]
            # A hair of tolerance, so that adding a line of code to a well-commented file does not
            # fail the build. Deleting the comments will still fail it.
            if current < floor - 0.02:
                print(
                    f"{name}: {current:.3f} comments per line, recorded floor {floor:.3f} "
                    f"({comments} comment lines, {code} code lines)"
                )
                failed = True
        if failed:
            print()
            print("see docs/COMMENTING.md; a file on the list may not lose its explanations")
            return 1
        print(f"{len(floors)} files at or above their recorded comment density")
        return 0

    total_comments = sum(c for c, _, _ in measured.values())
    total_code = sum(k for _, k, _ in measured.values())
    print(f"{len(measured)} files, {total_comments} comment lines, {total_code} code lines")
    print(f"overall {ratio(total_comments, total_code):.3f} comments per line of code")
    print(f"{len(floors)} of {len(measured)} files brought to the standard")
    print()
    print(f"{'ratio':>7}  {'cmt':>6}  {'code':>6}  file")
    for name, (comments, code, current) in sorted(
        measured.items(), key=lambda item: item[1][2]
    ):
        mark = "*" if name in floors else " "
        print(f"{current:7.3f}  {comments:6d}  {code:6d} {mark}{name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

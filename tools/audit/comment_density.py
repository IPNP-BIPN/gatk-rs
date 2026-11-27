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
    comment_density.py --record   rewrite the floors from the current tree (review the diff);
                                  a floor may rise, but lowering one needs --allow-lower
    comment_density.py --record --add FILE...
                                  also put FILE... on the list at their current density
    comment_density.py --self-test
                                  check the record refusal on a scratch tree
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys
import tempfile

# The repository root. Overridable with --root so the tests can point the script at a throwaway
# tree; every real invocation leaves it alone.
DEFAULT_ROOT = pathlib.Path(__file__).resolve().parents[2]


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


def sources(root: pathlib.Path) -> list[pathlib.Path]:
    """Every Rust source under `crates/`, excluding build output, sorted for a stable report."""
    return sorted(p for p in (root / "crates").rglob("*.rs") if "target" not in p.parts)


def ratio(comments: int, code: int) -> float:
    """Comments per line of code. A file with no code at all counts as fully commented."""
    return comments / code if code else float("inf")


def load_floors(listed: pathlib.Path) -> dict[str, float]:
    """The recorded floor for each file that has been brought to the standard.

    Format is one `<path> <ratio>` pair per line, `#` starting a comment. Paths are relative to the
    repository root so the file reads as a checklist.
    """
    floors: dict[str, float] = {}
    if not listed.exists():
        return floors
    for raw in listed.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        name, _, value = line.partition(" ")
        floors[name.strip()] = float(value.strip())
    return floors


def self_test() -> int:
    """Check, on a scratch tree, that --record raises a floor freely, refuses to lower one, lowers
    it when --allow-lower says so, and adds a file named by --add."""

    def run(root: pathlib.Path, *args: str) -> int:
        command = [sys.executable, __file__, "--root", str(root), *args]
        return subprocess.run(command, capture_output=True, text=True).returncode

    failures = []
    with tempfile.TemporaryDirectory(prefix="comment-density-self-test-") as tmp:
        root = pathlib.Path(tmp)
        source = root / "crates" / "demo" / "src" / "demo.rs"
        source.parent.mkdir(parents=True)
        listed = root / "tools" / "audit" / "commented.txt"
        listed.parent.mkdir(parents=True)
        name = "crates/demo/src/demo.rs"

        source.write_text("// a comment\n" * 10 + "let x = 1;\n" * 10, encoding="utf-8")
        listed.write_text(f"{name} 0.500\n", encoding="utf-8")
        if run(root, "--record") != 0 or "1.000" not in listed.read_text(encoding="utf-8"):
            failures.append("--record did not raise a floor")

        source.write_text("let x = 1;\n" * 10, encoding="utf-8")
        if run(root, "--record") != 1 or "1.000" not in listed.read_text(encoding="utf-8"):
            failures.append("--record lowered a floor without --allow-lower")
        if run(root, "--record", "--allow-lower") != 0 or (
            "0.000" not in listed.read_text(encoding="utf-8")
        ):
            failures.append("--record --allow-lower did not lower a floor")

        other = source.parent / "other.rs"
        other.write_text("// a comment\nlet x = 1;\n", encoding="utf-8")
        added = "crates/demo/src/other.rs 1.000"
        if run(root, "--record", "--add", "crates/demo/src/other.rs") != 0 or (
            added not in listed.read_text(encoding="utf-8")
        ):
            failures.append("--record --add did not list a file")
    for failure in failures:
        print(f"self-test: {failure}")
    if failures:
        return 1
    print("self-test: the record refusal holds")
    return 0


def main() -> int:
    """Parse the command line and run the report, the check or the record."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="fail on regression")
    parser.add_argument("--record", action="store_true", help="rewrite the floors")
    parser.add_argument(
        "--add",
        nargs="+",
        default=[],
        metavar="FILE",
        help="with --record, put these files on the list at their current density",
    )
    parser.add_argument(
        "--allow-lower",
        action="store_true",
        help="with --record, permit a floor to drop; without it, a drop is refused",
    )
    parser.add_argument("--self-test", action="store_true", help="test the record refusal")
    parser.add_argument(
        "--root", default=str(DEFAULT_ROOT), help="repository root; only the tests pass this"
    )
    args = parser.parse_args()
    if args.self_test:
        return self_test()

    root = pathlib.Path(args.root).resolve()
    listed = root / "tools" / "audit" / "commented.txt"

    measured: dict[str, tuple[int, int, float]] = {}
    for path in sources(root):
        comments, code = measure(path)
        name = str(path.relative_to(root))
        measured[name] = (comments, code, ratio(comments, code))

    floors = load_floors(listed)

    if args.record:
        for name in args.add:
            if name not in measured:
                print(f"not a Rust source under crates/: {name}", file=sys.stderr)
                return 1
            floors.setdefault(name, 0.0)
        # A floor may rise freely: a file gained explanations. A floor may only fall on request,
        # because the usual reason for a fall is that someone deleted comments and did not notice.
        lowered = [
            (name, floors[name], measured[name][2])
            for name in sorted(floors)
            if name in measured and measured[name][2] < floors[name] - 0.02
        ]
        if lowered and not args.allow_lower:
            for name, floor, current in lowered:
                print(f"{name}: recording would lower the floor {floor:.3f} to {current:.3f}")
            print()
            print("if that is deliberate, pass --allow-lower and say why in the commit message")
            return 1

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
        listed.write_text("\n".join(lines) + "\n", encoding="utf-8")
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

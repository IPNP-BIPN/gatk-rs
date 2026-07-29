#!/usr/bin/env python3
"""The one comparator for every conformance suite.

Before this file existed, each suite carried its own copy of the same thirty lines of parsing
and diffing inside `.github/workflows/ci.yml`, and the copies had already drifted: some skipped
comment lines, some did not; the metrics group applied a five-tool prefix list, the RnaSeq group
a two-tool one. A comparison that quietly ignores part of a file is worthless, so the rules that
say what is ignored now live in `manifest.json`, declared per suite, and this file is the only
thing that applies them.

The dump format, produced by the `*Dump.java` harnesses:

    <kind>\\t<case>\\t<payload>

one row per line, where the payload joins the record lines of a file with a *literal* backslash-n
(two characters), not a newline. That is why `SEGMENT_SEP` below is `'\\\\n'` and not `'\\n'`: a
comparison that split on real newlines would see one segment and silently canonicalize nothing.

Canonicalization rules are a closed set. Adding one is a reviewable event, because canonicalizing
is how a bit-identity claim quietly weakens.
"""

import gzip
import json
import sys
from pathlib import Path

# The payload separator emitted by the Java harnesses: a literal backslash followed by 'n'.
SEGMENT_SEP = "\\n"


# --------------------------------------------------------------------------------------
# Canonicalization rules. Each takes the payload and its rule spec, and returns the payload
# that will be compared. Every rule must be declared in the manifest with a `why`.
# --------------------------------------------------------------------------------------


def _strip_line_prefixes(payload, spec):
    """Drop whole segments starting with any declared prefix.

    Used for metrics files, whose `# <Tool> <command line>` and `# Started on: <timestamp>`
    headers carry JVM temp paths and a clock reading. The prefix list is per suite and holds
    *every* tool the harness runs, not just the current one: a multi-tool harness emits more
    than one header, and stripping only the current tool's left the others' temp paths in
    place. That bug made CI report 14 mismatches that were not divergences at all.
    """
    prefixes = tuple(spec["prefixes"])
    return SEGMENT_SEP.join(
        seg for seg in payload.split(SEGMENT_SEP) if not seg.startswith(prefixes)
    )


def _strip_pg(payload, spec):
    """Drop `@PG` segments: the provenance record whose `CL:` is the command line."""
    return SEGMENT_SEP.join(
        seg for seg in payload.split(SEGMENT_SEP) if not seg.startswith("@PG")
    )


def _strip_ur(payload, spec):
    """Drop `UR:` tab-fields: the reference's `file:` URI, which is path-dependent."""
    out = []
    for seg in payload.split(SEGMENT_SEP):
        out.append("\t".join(f for f in seg.split("\t") if not f.startswith("UR:")))
    return SEGMENT_SEP.join(out)


def _strip_banner(payload, spec):
    """Drop the metrics banner: `## htsjdk...` and `# ...` segments.

    This is the stripping `CmpCorpus.java` already applies on the Java side; declaring it here
    keeps the two sides describing the same operation.
    """
    return SEGMENT_SEP.join(
        seg
        for seg in payload.split(SEGMENT_SEP)
        if not (seg.startswith("## htsjdk") or seg.startswith("# "))
    )


RULES = {
    "strip_line_prefixes": _strip_line_prefixes,
    "strip_pg": _strip_pg,
    "strip_ur": _strip_ur,
    "strip_banner": _strip_banner,
}


def canonicalize(kind, payload, rules):
    """Apply the suite's declared rules, in order, to one row's payload."""
    for spec in rules:
        name = spec["rule"]
        if name not in RULES:
            raise SystemExit(f"unknown canonicalization rule: {name}")
        # A rule may be restricted to certain row kinds, e.g. `strip_ur` applies to the `dict`
        # row and not to the inputs beside it.
        kinds = spec.get("kinds")
        if kinds is not None and kind not in kinds:
            continue
        payload = RULES[name](payload, spec)
    return payload


# --------------------------------------------------------------------------------------
# Readers
# --------------------------------------------------------------------------------------


def _open(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_keyed(path, compare_spec):
    """Read a dump into {(kind, case): canonical payload}."""
    rules = compare_spec.get("rules", [])
    skip_kinds = set(compare_spec.get("skip_kinds", []))
    skip_comments = compare_spec.get("skip_comment_lines", True)
    rows = {}
    with _open(path) as fh:
        for line in fh:
            if skip_comments and line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t", 2)
            if len(parts) < 3 or parts[0] in skip_kinds:
                continue
            kind, case, payload = parts
            rows[(kind, case)] = canonicalize(kind, payload, rules)
    return rows


def parse_lines(path, compare_spec):
    """Read a dump as a line sequence.

    Suites whose corpus has several rows of one kind per case (the shards of a split, the
    inputs of a merge) cannot use a (kind, case) key: it would collapse those rows and compare
    only the last one.
    """
    skip_comments = compare_spec.get("skip_comment_lines", True)
    with _open(path) as fh:
        return [
            line.rstrip("\n")
            for line in fh
            if line.strip() and not (skip_comments and line.startswith("#"))
        ]


# --------------------------------------------------------------------------------------
# Comparison
# --------------------------------------------------------------------------------------


def compare_case(real_path, golden_path, compare_spec):
    """Compare one regenerated dump against its committed golden.

    Returns (ok, compared_count, message lines).
    """
    mode = compare_spec.get("mode", "keyed")
    out = []

    if mode == "lines":
        real = parse_lines(real_path, compare_spec)
        committed = parse_lines(golden_path, compare_spec)
        if real != committed:
            out.append(f"lines differ: real={len(real)} committed={len(committed)}")
            for i, (r, c) in enumerate(zip(real, committed)):
                if r != c:
                    out.append(f"  first diff at line {i}")
                    out.append(f"    real     ={r[:200]}")
                    out.append(f"    committed={c[:200]}")
                    break
            return False, len(real), out
        return True, len(real), out

    if mode != "keyed":
        raise SystemExit(f"unknown compare mode: {mode}")

    real = parse_keyed(real_path, compare_spec)
    committed = parse_keyed(golden_path, compare_spec)
    if set(real) != set(committed):
        out.append(f"row sets differ: {sorted(set(real) ^ set(committed))}")
        return False, len(real), out
    bad = [k for k in real if real[k] != committed[k]]
    for k in bad[:5]:
        out.append(f"  {k}")
        out.append(f"    real     ={real[k][:200]}")
        out.append(f"    committed={committed[k][:200]}")
    return not bad, len(real), out


def load_manifest(path=None):
    path = Path(path or Path(__file__).with_name("manifest.json"))
    with open(path) as fh:
        return json.load(fh)


def suite_by_id(manifest, suite_id):
    for suite in manifest["suites"]:
        if suite["id"] == suite_id:
            return suite
    raise SystemExit(f"no suite {suite_id!r} in the manifest")


def main(argv):
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest")
    ap.add_argument("--suite", required=True)
    ap.add_argument(
        "--real",
        action="append",
        required=True,
        metavar="DUMP=PATH",
        help="the regenerated dump for one case, e.g. QualityYieldDump=/tmp/qy.txt",
    )
    args = ap.parse_args(argv)

    manifest = load_manifest(args.manifest)
    suite = suite_by_id(manifest, args.suite)
    reals = dict(pair.split("=", 1) for pair in args.real)

    failed, total = 0, 0
    for case in suite["cases"]:
        dump = case["dump"]
        if dump not in reals:
            print(f"FAIL {suite['id']}/{dump}: no regenerated dump supplied")
            failed += 1
            continue
        ok, compared, messages = compare_case(reals[dump], case["golden"], suite["compare"])
        total += compared
        status = "ok  " if ok else "FAIL"
        print(f"{status} {suite['id']}/{dump}: compared={compared}")
        for line in messages:
            print(line)
        failed += 0 if ok else 1

    print(f"suite={suite['id']} cases={len(suite['cases'])} compared={total} failed={failed}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

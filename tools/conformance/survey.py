#!/usr/bin/env python3
"""Run the goldens we already have against another oracle, and classify what moved.

Written for the move to a newer reference (GATK 4.7.0.0, htsjdk 5.0.0, Picard 3.5.0). Moving the
target does not adjust 216 oracle-backed suites, it re-opens them, and the only honest way to size
that is to measure it: run every suite against the candidate image and sort the differences into
the ones that are the reference saying its own version and the ones that are the reference
behaving differently.

    python3 tools/conformance/survey.py --oracle-image gatk-rs-oracle:4.7.0.0 \
        --version-was 4.6.2.0 --version-now 4.7.0.0

Each suite lands in one of three buckets:

  * **identical**, the golden matches under the new oracle with no help at all;
  * **version stamp only**, it matches once the stamp is taken out of both sides: the version
    string, and the arguments the new version adds to every command line it prints
    (`NEW_ARGUMENTS`), in the dump's text and in the BAM headers it carries as base64;
  * **behaviour**, something else moved, and the first differing row is printed.

The second bucket is a *heuristic and is named as one*: a genuine difference that happened to
contain the literal new version string would be swallowed by it. It is here to make the third
bucket small enough to read, not to decide anything. Nothing produced by this script may be
committed as a golden: a golden comes from the pinned container on a real x86-64 runner, and this
runs wherever it is invoked.
"""

import argparse
import base64
import gzip
import json
import re
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import compare as comparator  # noqa: E402
import run_suite  # noqa: E402

REPO = Path(__file__).resolve().parents[2]


# What the new version stamps into a command line: 4.7.0.0 adds one global argument, and both
# `@PG CL:` and `##GATKCommandLine` print every argument with its value.
NEW_ARGUMENTS = ["--output-cram-version 3.1"]
BASE64 = re.compile(r"^[A-Za-z0-9+/=]{40,}$")


def unstamped_bytes(raw, was, now):
    """A file the dump carries as base64, decompressed, with the stamp taken out of it.

    In a BAM the stamp sits in the header text, whose length is a field of its own and whose
    every later byte, and every index offset, it shifts; so a BAM keeps its text and records, and
    an index is reduced to a marker, because a BAM that matches makes its index match too.
    """
    if raw[:2] == b"\x1f\x8b":
        try:
            raw = gzip.decompress(raw)
        except (OSError, EOFError):
            pass
    if raw[:4] in (b"BAI\x01", b"TBI\x01"):
        return b"INDEX"
    raw = raw.replace(now.encode(), was.encode())
    head, text, rest = b"", raw, b""
    if raw[:4] == b"BAM\x01":
        n = struct.unpack("<i", raw[4:8])[0]
        head, text, rest = b"BAM\x01", raw[8 : 8 + n], raw[8 + n :]
    for argument in NEW_ARGUMENTS:
        text = text.replace(b" " + argument.encode(), b"").replace(argument.encode() + b" ", b"")
    return head + text + rest


def unstamped(text, was, now):
    """A dump's text with the version and the new arguments put back to what 4.6.2.0 printed."""
    lines = []
    for line in text.replace(now, was).splitlines():
        fields = line.split("\t")
        for i, field in enumerate(fields):
            if BASE64.match(field):
                try:
                    fields[i] = "unstamped:" + unstamped_bytes(base64.b64decode(field), was, now).hex()
                except ValueError:
                    pass
        line = "\t".join(fields)
        for argument in NEW_ARGUMENTS:
            line = line.replace(" " + argument, "").replace(argument + " ", "")
        lines.append(line)
    return "\n".join(lines) + "\n"


def rewrite(real, golden, was, now):
    """Both sides with the stamp taken out, so only what the reference DOES is left to differ."""
    def read(path):
        return (gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)).read()

    out_real, out_golden = Path(str(real) + ".unstamped"), Path(str(real) + ".golden.unstamped")
    out_real.write_text(unstamped(read(real), was, now))
    out_golden.write_text(unstamped(read(golden), was, now))
    return out_real, out_golden


def survey_suite(manifest, suite, workdir, was, now, fixtures, keep=None):
    """Run one suite under the overridden oracle and say which bucket it falls in."""
    props = suite.get("java_props", manifest.get("default_java_props", []))
    if suite.get("oracle"):
        # It names an image of its own, which the override does not reach.
        return "error", [f"runs in its own image ({manifest[suite['oracle']]['image']}), not surveyed"]
    if not suite.get("needs_fixtures"):
        fixtures = None
    elif fixtures is None:
        return "error", ["could not build the fixtures"]
    verdict, notes = "identical", []
    for case in suite["cases"]:
        dump = case["dump"]
        real = Path(workdir) / f"{suite['id']}.{dump}.txt"
        rc = run_suite.docker_run(
            manifest["oracle"], suite["harness"], dump, props, real, fixtures,
            suite.get("python_image"), suite,
        )
        if rc != 0:
            return "error", [f"{dump}: the dump did not run under this oracle"]
        golden = REPO / case["golden"]
        ok, _, messages = comparator.compare_case(real, golden, suite["compare"])
        if ok:
            continue
        ok, _, messages = comparator.compare_case(*rewrite(real, golden, was, now), suite["compare"])
        if ok:
            verdict = "version" if verdict != "behaviour" else verdict
            notes.append(f"{dump}: the version stamp, and nothing else")
            continue
        verdict = "behaviour"
        notes.append(f"{dump}: " + "; ".join(messages[:3]))
        if keep:
            Path(keep).mkdir(parents=True, exist_ok=True)
            (Path(keep) / real.name).write_text(real.read_text())
    return verdict, notes


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest")
    ap.add_argument("--oracle-image", required=True)
    ap.add_argument("--version-was", required=True)
    ap.add_argument("--version-now", required=True)
    ap.add_argument("--suites", help="space separated ids; default is every oracle-backed suite")
    ap.add_argument("--shard", help="k/n: survey only every n-th suite, starting at the k-th (1-based)")
    ap.add_argument("--json", help="also write the buckets and notes to this file")
    ap.add_argument("--keep", help="copy the regenerated dump of every suite that moved here")
    args = ap.parse_args(argv)

    manifest = comparator.load_manifest(args.manifest)
    # The fixtures are inputs, so they come from the manifest's own image: htsjdk 5.0.0 deflates
    # with libdeflate by default, and fixtures it wrote would move the inputs along with the tool.
    pinned = dict(manifest)
    print(
        f"!! surveying against {args.oracle_image}, not {manifest['oracle']['image']}.\n"
        f"!! Nothing produced here may be committed as a golden.",
        flush=True,
    )
    manifest["oracle"] = dict(manifest["oracle"], image=args.oracle_image)

    ids = (
        args.suites.split()
        if args.suites
        else [s["id"] for s in manifest["suites"] if s["status"] == "oracle-backed"]
    )
    if args.shard:
        k, n = (int(x) for x in args.shard.split("/"))
        ids = ids[k - 1 :: n]
    buckets = {"identical": [], "version": [], "behaviour": [], "error": []}
    all_notes = {}
    with tempfile.TemporaryDirectory() as workdir:
        fixtures = None
        if any(comparator.suite_by_id(manifest, i).get("needs_fixtures") for i in ids):
            fixtures = Path(workdir) / "fixtures"
            if run_suite.build_fixtures(pinned, fixtures) != 0:
                fixtures = None
        for suite_id in ids:
            suite = comparator.suite_by_id(manifest, suite_id)
            verdict, notes = survey_suite(
                manifest, suite, workdir, args.version_was, args.version_now, fixtures, args.keep
            )
            buckets[verdict].append(suite_id)
            all_notes[suite_id] = notes
            print(f"{verdict:10} {suite_id}", flush=True)
            for note in notes:
                print(f"           {note}", flush=True)

    print()
    for name, title in [
        ("identical", "identical under the new oracle"),
        ("version", "the version stamp only"),
        ("behaviour", "something else moved"),
        ("error", "did not run"),
    ]:
        print(f"{len(buckets[name]):4}  {title}")
        if name in ("behaviour", "error") and buckets[name]:
            for suite_id in buckets[name]:
                print(f"        {suite_id}")
    if args.json:
        Path(args.json).write_text(json.dumps({"buckets": buckets, "notes": all_notes}, indent=1))
    return 1 if buckets["behaviour"] or buckets["error"] else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

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
  * **version stamp only**, it matches once the new version string is rewritten to the old one,
    which is the `@PG` line's `VN:` and nothing else;
  * **behaviour**, something else moved, and the first differing row is printed.

The second bucket is a *heuristic and is named as one*: a genuine difference that happened to
contain the literal new version string would be swallowed by it. It is here to make the third
bucket small enough to read, not to decide anything. Nothing produced by this script may be
committed as a golden: a golden comes from the pinned container on a real x86-64 runner, and this
runs wherever it is invoked.
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import compare as comparator  # noqa: E402
import run_suite  # noqa: E402

REPO = Path(__file__).resolve().parents[2]


def rewrite(path, was, now):
    """A copy of the regenerated dump with the new version string put back to the old one."""
    out = Path(str(path) + ".versioned")
    out.write_text(path.read_text().replace(now, was))
    return out


def survey_suite(manifest, suite, workdir, was, now):
    """Run one suite under the overridden oracle and say which bucket it falls in."""
    props = suite.get("java_props", manifest.get("default_java_props", []))
    fixtures = None
    if suite.get("needs_fixtures"):
        fixtures = Path(workdir) / f"fixtures-{suite['id']}"
        if run_suite.build_fixtures(manifest, fixtures) != 0:
            return "error", ["could not build the fixtures"]
    verdict, notes = "identical", []
    for case in suite["cases"]:
        dump = case["dump"]
        real = Path(workdir) / f"{suite['id']}.{dump}.txt"
        if run_suite.docker_run(manifest, suite["harness"], dump, props, real, fixtures) != 0:
            return "error", [f"{dump}: the dump did not run under this oracle"]
        ok, _, messages = comparator.compare_case(real, case["golden"], suite["compare"])
        if ok:
            continue
        ok, _, messages = comparator.compare_case(
            rewrite(real, was, now), case["golden"], suite["compare"]
        )
        if ok:
            verdict = "version" if verdict != "behaviour" else verdict
            notes.append(f"{dump}: the version stamp, and nothing else")
            continue
        verdict = "behaviour"
        notes.append(f"{dump}: " + "; ".join(messages[:3]))
    return verdict, notes


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest")
    ap.add_argument("--oracle-image", required=True)
    ap.add_argument("--version-was", required=True)
    ap.add_argument("--version-now", required=True)
    ap.add_argument("--suites", help="space separated ids; default is every oracle-backed suite")
    args = ap.parse_args(argv)

    manifest = comparator.load_manifest(args.manifest)
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
    buckets = {"identical": [], "version": [], "behaviour": [], "error": []}
    with tempfile.TemporaryDirectory() as workdir:
        for suite_id in ids:
            suite = comparator.suite_by_id(manifest, suite_id)
            verdict, notes = survey_suite(manifest, suite, workdir, args.version_was, args.version_now)
            buckets[verdict].append(suite_id)
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
    return 1 if buckets["behaviour"] or buckets["error"] else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

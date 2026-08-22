#!/usr/bin/env python3
"""Run one conformance suite (or one probe) against the oracle image.

The CI job and a local run go through this same script, so a suite that passes locally and fails
in CI is a real difference in the environment rather than a difference in how it was invoked.

    python3 tools/conformance/run_suite.py --suites metrics
    python3 tools/conformance/run_suite.py --suites "metrics rnaseq snvq"
    python3 tools/conformance/run_suite.py --probe rnaseq-overlap-order
    python3 tools/conformance/run_suite.py --list

`--oracle-image` runs a suite against another image than the manifest's. It exists for one job:
surveying a candidate reference version before anything is migrated to it, by running the goldens
we already have against the newer oracle and reading off which ones move. A suite that fails under
an overridden image is not a regression, it is a measurement, so the override is printed on every
run and refuses to be silent.

The oracle image must exist; build it with

    docker build --platform linux/amd64 -t picard-rs-oracle:3.4.0 tools/oracle

Goldens are only valid when produced on real x86-64 (docs/decisions/0004 and the README's
bit-identity contract), so a local run on Apple Silicon is a smoke test, not a source of goldens.
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import compare as comparator  # noqa: E402

REPO = Path(__file__).resolve().parents[2]


def container_command(cls, props):
    """The in-container script: compile the harness against the pinned jar, then run it.

    htsjdk chatters on stderr, so a successful run's stderr is dropped. It is kept aside rather
    than discarded outright and replayed when the run fails: an exception thrown out of a dump's
    main used to leave nothing behind at all, which reads exactly like a dump that produced no
    rows, and the two need different fixes.

    The whole harness directory is copied and `-sourcepath .` is passed, not just the one class.
    Copying a single file is what the first version did, and it meant three of the six read-filter
    dumps could not compile in CI at all: `CountingFilterDump`, `ReadCoordinateDump` and
    `ReadClipperDump` share `ReadFilterDump.corpus`, which was not there to compile against. Their
    goldens had been produced by an ad-hoc command instead, so nothing re-derived them, which is
    the exact failure picard-rs decision 0008 is about.
    """
    prop_str = (" ".join(props) + " ") if props else ""
    return (
        f'cp /harness/*.java . && javac -cp "$ORACLE_CP" -sourcepath . -d . {cls}.java '
        f'&& {{ java {prop_str}-cp ".:$ORACLE_CP" {cls} 2>/tmp/dump-stderr '
        f'|| {{ cat /tmp/dump-stderr >&2; exit 1; }}; }}'
    )


def build_fixtures(manifest, into):
    """Materialize the shared fixture corpus with tools/coverage/MakeFixtures.java.

    A harness that needs real input files gets them from the same program the covering arrays use,
    so a rejection case and an array row read the same bytes.
    """
    oracle = manifest["oracle"]
    into.mkdir(parents=True, exist_ok=True)
    command = (
        'cp /harness/MakeFixtures.java . && javac -cp "$ORACLE_CP" -d . MakeFixtures.java '
        '&& java -Dsamjdk.try_use_intel_deflater=false -cp ".:$ORACLE_CP" MakeFixtures /out'
    )
    return subprocess.run(
        [
            "docker", "run", "--rm", "--platform", oracle["platform"],
            "-v", f"{REPO}/tools/coverage:/harness:ro",
            "-v", f"{into}:/out",
            "-w", "/work", oracle["image"], command,
        ],
        capture_output=True,
        text=True,
    ).returncode


def docker_run(manifest, harness, cls, props, stdout, fixtures=None):
    oracle = manifest["oracle"]
    mounts = ["-v", f"{REPO}/{harness}:/harness:ro"]
    if fixtures is not None:
        mounts += ["-v", f"{fixtures}:/work/fixtures:ro"]
    cmd = [
        "docker",
        "run",
        "--rm",
        "--platform",
        oracle["platform"],
        *mounts,
        "-w",
        "/work",
        oracle["image"],
        container_command(cls, props),
    ]
    print("+ " + " ".join(cmd[:-1]) + f" '{cmd[-1][:60]}...'", flush=True)
    with open(stdout, "w") as fh:
        return subprocess.run(cmd, stdout=fh).returncode


PENDING_DIR = REPO / "tools" / "conformance" / "pending"


def run_pending(manifest, suite, workdir):
    """A suite with no golden yet: run it, check the shape, and leave the dump for CI to publish.

    The alternative was to generate the golden here and commit it. That is exactly what produced
    the sixteen goldens of decision 0008, whose provenance turned out to be a laptop rather than
    the pinned container, so it is refused: this prints what the reference did, asserts only the
    row count the suite declares, and says plainly that nothing was compared.
    """
    props = suite.get("java_props", manifest.get("default_java_props", []))
    PENDING_DIR.mkdir(parents=True, exist_ok=True)
    fixtures = None
    if suite.get("needs_fixtures"):
        fixtures = Path(workdir) / "fixtures"
        if build_fixtures(manifest, fixtures) != 0:
            print(f"FAIL {suite['id']}: could not build the fixtures")
            return 1
    failed = 0
    for case in suite["cases"]:
        dump = case["dump"]
        out = PENDING_DIR / f"{suite['id']}.{dump}.txt"
        rc = docker_run(manifest, suite["harness"], dump, props, out, fixtures)
        rows = [l for l in open(out) if l.strip() and not l.startswith("#")]
        print(f"--- {suite['id']}/{dump}: {len(rows)} rows, nothing compared (no golden yet)")
        for line in rows:
            print("   ", line.rstrip("\n")[:200])
        expected = suite.get("expect_rows")
        if rc != 0 or (expected is not None and len(rows) != expected):
            print(f"FAIL {suite['id']}/{dump}: exit {rc}, {len(rows)} rows, expected {expected}")
            failed += 1
            continue
        # A row count is not evidence: the first run of this suite produced four rows that all
        # said "Cannot read non-existent file", because the fixtures were not mounted, and the
        # count alone called that a pass. The behaviours the suite exists for are named, and each
        # one must appear.
        body = "".join(rows)
        for phrase in suite.get("expect_contains", []):
            if phrase not in body:
                print(f"FAIL {suite['id']}/{dump}: expected a row containing {phrase!r}")
                failed += 1
    print(
        f"suite={suite['id']} status=golden-pending cases={len(suite['cases'])} failed={failed}; "
        f"dumps in {PENDING_DIR} are the candidate goldens, valid only from a real x86-64 run"
    )
    return failed


def run_suite(manifest, suite, workdir):
    if suite["status"] == "golden-pending":
        return run_pending(manifest, suite, workdir)
    props = suite.get("java_props", manifest.get("default_java_props", []))
    failed = 0
    reals = {}
    fixtures = None
    if suite.get("needs_fixtures"):
        fixtures = Path(workdir) / "fixtures"
        if build_fixtures(manifest, fixtures) != 0:
            print(f"FAIL {suite['id']}: could not build the fixtures")
            return 1
    # The regenerated dumps are kept, not thrown away with the temporary directory. A dump that
    # grows a case needs its golden refreshed, and the refresh may only come from the pinned
    # container on real x86-64; without this the only way to update a golden was to have the image
    # locally, which is how goldens of unknown provenance get committed in the first place.
    PENDING_DIR.mkdir(parents=True, exist_ok=True)
    for case in suite["cases"]:
        dump = case["dump"]
        out = PENDING_DIR / f"{suite['id']}.{dump}.txt"
        rc = docker_run(manifest, suite["harness"], dump, props, out, fixtures)
        lines = sum(1 for _ in open(out))
        print(f"regenerated {dump}: {lines} lines (docker exit {rc})")
        if rc != 0 or lines == 0:
            print(f"FAIL {suite['id']}/{dump}: the oracle produced no dump")
            failed += 1
            continue
        reals[dump] = out

    total = 0
    for case in suite["cases"]:
        dump = case["dump"]
        if dump not in reals:
            continue
        ok, compared, messages = comparator.compare_case(
            reals[dump], REPO / case["golden"], suite["compare"]
        )
        total += compared
        print(f"{'ok  ' if ok else 'FAIL'} {suite['id']}/{dump}: compared={compared}")
        for line in messages:
            print(line)
        failed += 0 if ok else 1

    print(
        f"suite={suite['id']} status={suite['status']} cases={len(suite['cases'])} "
        f"compared={total} failed={failed}"
    )
    return failed


def run_probe(manifest, probe, workdir):
    out = Path(workdir) / f"probe.{probe['id']}.txt"
    props = probe.get("java_props", manifest.get("default_java_props", []))
    docker_run(manifest, probe["harness"], probe["class"], props, out)
    text = open(out).read()
    print(text)
    if probe["expect"] not in text:
        print(f"FAIL probe {probe['id']}: expected {probe['expect']!r}")
        print(probe["on_failure"])
        return 1
    print(f"ok   probe {probe['id']}: {probe['expect']}")
    return 0


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest")
    ap.add_argument(
        "--suites",
        help="one or more suite ids, space separated. Several suites share a CI job because each "
        "job pays the oracle image restore once.",
    )
    ap.add_argument("--probe")
    ap.add_argument("--list", action="store_true")
    ap.add_argument(
        "--oracle-image",
        help="run against this image instead of the manifest's. For surveying a candidate "
        "reference version: a failure here is a measurement, not a regression.",
    )
    args = ap.parse_args(argv)

    manifest = comparator.load_manifest(args.manifest)
    if args.oracle_image:
        print(
            f"!! oracle overridden: {manifest['oracle']['image']} -> {args.oracle_image}.\n"
            f"!! Differences below are measurements of that image, not regressions, and nothing\n"
            f"!! produced under an overridden oracle may be committed as a golden.",
            flush=True,
        )
        manifest["oracle"] = dict(manifest["oracle"], image=args.oracle_image)

    if args.list:
        for suite in manifest["suites"]:
            print(f"{suite['id']:28} {suite['status']:14} {len(suite['cases'])} case(s)")
        for probe in manifest.get("probes", []):
            print(f"{probe['id']:28} {'probe':14} {probe['class']}")
        return 0

    with tempfile.TemporaryDirectory() as workdir:
        if args.suites:
            ids = args.suites.split()
            # Every suite runs even after one fails: the run exists to say which suites diverge,
            # not that at least one does.
            failed = sum(
                run_suite(manifest, comparator.suite_by_id(manifest, suite_id), workdir)
                for suite_id in ids
            )
            print(f"suites={len(ids)} failing={failed}")
            return 1 if failed else 0
        if args.probe:
            for probe in manifest.get("probes", []):
                if probe["id"] == args.probe:
                    return run_probe(manifest, probe, workdir)
            raise SystemExit(f"no probe {args.probe!r} in the manifest")

    ap.error("pass --suites, --probe or --list")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

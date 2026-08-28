#!/usr/bin/env python3
"""Run, locally, every CI gate that does not need the runner's silicon.

Written because a PR went red on `expect_rows`: a number that had been counted by hand from
`wc -l` rather than taken from the harness, which the suite runner would have caught in seconds.
The point of this script is that "I ran the tests" and "I ran what CI runs" stop being different
claims. It mirrors .github/workflows/ci.yml in the workflow's own order and stops at the first
failure.

What it deliberately does NOT do: produce goldens. The oracle jobs run the pinned container and a
golden may only come from a real x86-64 runner; `--suites` here is a smoke test of the harness and
of the manifest's own expectations, which is exactly the part that has broken.

Usage:
    python3 tools/preflight.py                     # every gate, no suite runs
    python3 tools/preflight.py --suites a,b        # and run these suites through the harness
    python3 tools/preflight.py --fast              # skip the release test build
"""

import argparse
import subprocess
import sys
import time

REPO = subprocess.run(["git", "rev-parse", "--show-toplevel"], capture_output=True, text=True,
                      check=True).stdout.strip()


def step(name, command, shell=False):
    """One CI step, printed the way the workflow names it."""
    started = time.time()
    print(f"\n=== {name}", flush=True)
    result = subprocess.run(command, cwd=REPO, shell=shell)
    elapsed = time.time() - started
    if result.returncode != 0:
        print(f"\nFAILED after {elapsed:.0f}s: {name}", file=sys.stderr)
        sys.exit(result.returncode)
    print(f"--- ok ({elapsed:.0f}s)", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--suites", default="",
                        help="comma-separated suite ids to run through the oracle harness")
    parser.add_argument("--fast", action="store_true",
                        help="skip `cargo test --workspace --release`, which CI runs in release")
    arguments = parser.parse_args()

    # guard: repository invariants
    step("LICENSE and README state the same licence",
         'grep -q "Apache License" LICENSE '
         '&& grep -q "Apache License 2.0, matching GATK" README.md '
         '&& grep -qi "not the official GATK" README.md', shell=True)
    step("The workflow matches the conformance manifest",
         ["python3", "tools/conformance/generate_ci.py", "--check"])
    step("Every committed golden belongs to a declared suite",
         ["python3", "tools/conformance/audit_goldens.py"])

    # test: the x86-64 job, minus the parts that need the runner
    step("Every ported symbol comes from a licence-compatible source",
         ["python3", "tools/audit/provenance.py", "crates"])
    step("A file already explained does not lose its explanations",
         ["python3", "tools/audit/comment_density.py", "--check"])
    step("cargo fmt", ["cargo", "fmt", "--all", "--", "--check"])
    step("cargo clippy",
         ["cargo", "clippy", "--workspace", "--all-targets", "--", "-D", "warnings"])
    if not arguments.fast:
        # CI tests in release. A debug-only run has passed while release failed before.
        step("cargo test --release", ["cargo", "test", "--workspace", "--release"])
    step("The ports below are pinned by revision, and every revision resolves",
         'if grep -E \'^(htsjdk|picard|jmath)[a-z-]* = \\{ git\' Cargo.toml | grep -qv \'rev = "\'; '
         'then echo "a dependency is not pinned to a revision"; exit 1; fi; '
         # Checking the form of the pin is not enough: one pin named the wrong repository's SHA for
         # months because no crate depended on it, so cargo never resolved it. Cargo.lock is cargo's
         # own record of having resolved a revision.
         'for revision in $(grep -oE \'rev = "[0-9a-f]+"\' Cargo.toml | grep -oE \'[0-9a-f]+\' | sort -u); do '
         'if ! grep -q "$revision" Cargo.lock; then '
         'echo "revision $revision is pinned in Cargo.toml and resolved nowhere in Cargo.lock"; exit 1; fi; done',
         shell=True)

    # dashboard
    step("docs/STATUS.md is generated, not written",
         ["python3", "tools/dashboard/generate.py", "--check"])

    # The declarations are generated from the golden and cross-checked against the inventory on the
    # way through, so this step is two gates in one: the module is what the generator writes, and
    # the two readings of the reference still agree.
    step("the tool declarations are generated, not written",
         ["python3", "tools/declarations/generate.py", "--check"])

    # The suites, last, because they are the slowest and need the oracle image.
    suites = [s for s in arguments.suites.split(",") if s]
    if suites:
        step(f"oracle harness: {','.join(suites)}",
             ["python3", "tools/conformance/run_suite.py", "--suites", ",".join(suites)])

    print("\npreflight: every gate green")


if __name__ == "__main__":
    main()

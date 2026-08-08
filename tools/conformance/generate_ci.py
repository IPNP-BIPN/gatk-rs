#!/usr/bin/env python3
"""Generate the oracle jobs of `.github/workflows/ci.yml` from `manifest.json`.

The workflow used to hold one hand-written pair of steps per suite: a `docker run` that
regenerated the golden and a thirty-line Python comparator inlined in YAML. At 39 suites that was
1556 lines, the comparators had already drifted from each other, and each new tool made the file
longer than the port it was testing. At 311 tools it would have been unreadable.

Now the suites are data. This script turns them into:

  * `oracle-image`, which builds the digest-pinned image once and caches it, and
  * `oracle`, a matrix with one entry per suite, so suites run in parallel and a broken suite no
    longer hides the ones after it in a single sequential job.

Usage:
    python3 tools/conformance/generate_ci.py            # write .github/workflows/ci.yml
    python3 tools/conformance/generate_ci.py --check    # fail if the committed file is stale
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import compare as comparator  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
TEMPLATE = HERE / "ci.template.yml"
WORKFLOW = REPO / ".github" / "workflows" / "ci.yml"
MARKER = "# {{ORACLE_JOBS}}"

CACHE_KEY = "gatk-rs-oracle-${{ hashFiles('tools/oracle/**') }}"
IMAGE_TAR = "/tmp/oracle-image.tar"


def wrap(text, indent, width=96):
    """Wrap a `why` into YAML comment lines, so the reason stays next to the job that needs it."""
    words = text.split()
    lines, current = [], indent + "#"
    for word in words:
        if len(current) + 1 + len(word) > width:
            lines.append(current)
            current = indent + "#"
        current += " " + word
    lines.append(current)
    return "\n".join(lines)


def oracle_jobs(manifest):
    suites = manifest["suites"]
    probes = manifest.get("probes", [])
    oracle = manifest["oracle"]
    out = []

    out.append(
        f"""  oracle-image:
    name: Build the pinned oracle image
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Record the host CPU
        run: grep -m1 'model name' /proc/cpuinfo; grep -m1 'vendor_id' /proc/cpuinfo

      - id: cache
        uses: actions/cache@v4
        with:
          path: {IMAGE_TAR}
          key: {CACHE_KEY}

      - name: Build the oracle image
        # The probe runs during the build, so a degraded environment cannot produce an image, let
        # alone a golden. The image is built once and shared with every suite below, because
        # rebuilding it per suite would multiply a five-minute build by the number of suites.
        if: steps.cache.outputs.cache-hit != 'true'
        run: |
          docker build --platform {oracle['platform']} -t {oracle['image']} {oracle['context']}
          docker save {oracle['image']} -o {IMAGE_TAR}
"""
    )

    per_job = manifest.get("ci", {}).get("suites_per_job", 5)
    rows = []
    # Grouped by status, so a green run never lets an `unchecked` suite read as the stronger
    # claim, and chunked, because every job pays the oracle image restore once: one job per suite
    # would download the image forty-odd times per run for no extra information.
    for status in ("oracle-backed", "unchecked", "golden-pending"):
        group = [s for s in suites if s["status"] == status]
        chunks = [group[i : i + per_job] for i in range(0, len(group), per_job)]
        for n, chunk in enumerate(chunks, 1):
            ids = " ".join(s["id"] for s in chunk)
            # The artefact is named after the suites in the job, not the job's position in the
            # matrix. Adding a suite renumbers the matrix, and a golden that took the previous
            # slice's name is a golden published under the wrong claim.
            rows.append(
                f"          - group: {status}\n"
                f"            index: \"{n}/{len(chunks)}\"\n"
                f"            slug: \"{'-'.join(s['id'] for s in chunk)}\"\n"
                f"            suites: \"{ids}\""
            )
    matrix = "\n".join(rows)

    out.append(
        f"""  oracle:
    name: "oracle ${{{{ matrix.group }}}} ${{{{ matrix.index }}}}"
    needs: oracle-image
    runs-on: ubuntu-latest
    strategy:
      # One suite's divergence must not cancel the others: the point of the run is to learn which
      # suites diverge, not to learn that at least one does.
      fail-fast: false
      matrix:
        include:
{matrix}
    steps:
      - uses: actions/checkout@v4

      - uses: actions/cache/restore@v4
        with:
          path: {IMAGE_TAR}
          key: {CACHE_KEY}
          fail-on-cache-miss: true

      - run: docker load -i {IMAGE_TAR}

      - name: Regenerate and compare
        # Every rule that canonicalizes anything away is declared in tools/conformance/manifest.json
        # with the reason it is there; tools/conformance/compare.py is the only code that applies
        # them. Run the same thing locally with:
        #   python3 tools/conformance/run_suite.py --suites "${{{{ matrix.suites }}}}"
        run: python3 tools/conformance/run_suite.py --suites "${{{{ matrix.suites }}}}"

      - name: Publish what the container produced
        # Two cases, one artefact. A `golden-pending` suite has no golden to compare against and
        # one may not be committed from a developer machine (decision 0008). An `oracle-backed`
        # suite whose dump has grown a case needs its golden refreshed, and that refresh has to
        # come from here for the same reason. Both land in tools/conformance/pending/, and
        # committing the file is what makes the claim.
        if: always()
        uses: actions/upload-artifact@v4
        with:
          name: candidate-goldens-${{{{ matrix.slug }}}}
          path: tools/conformance/pending/
          if-no-files-found: ignore
"""
    )

    if probes:
        probe_rows = "\n".join(f"          - probe: {p['id']}" for p in probes)
        why_lines = "\n".join(wrap(p["why"], "        ") for p in probes)
        out.append(
            f"""  oracle-probes:
    name: "probe: ${{{{ matrix.probe }}}}"
    needs: oracle-image
    runs-on: ubuntu-latest
    strategy:
      fail-fast: false
      matrix:
        include:
{probe_rows}
    steps:
      - uses: actions/checkout@v4

      - uses: actions/cache/restore@v4
        with:
          path: {IMAGE_TAR}
          key: {CACHE_KEY}
          fail-on-cache-miss: true

      - run: docker load -i {IMAGE_TAR}

      - name: Run the probe
{why_lines}
        run: python3 tools/conformance/run_suite.py --probe ${{{{ matrix.probe }}}}
"""
        )

    return "\n".join(out)


def render(manifest):
    template = TEMPLATE.read_text()
    if MARKER not in template:
        raise SystemExit(f"{TEMPLATE} has lost its {MARKER} marker")
    banner = (
        "# The jobs below this line are generated from tools/conformance/manifest.json by\n"
        "# tools/conformance/generate_ci.py. Edit the manifest, not this file: the `guard` job\n"
        "# fails if the two disagree.\n"
    )
    return template.replace(MARKER, banner + oracle_jobs(manifest))


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="fail if the committed workflow is stale")
    args = ap.parse_args(argv)

    manifest = comparator.load_manifest()
    rendered = render(manifest)

    if args.check:
        current = WORKFLOW.read_text() if WORKFLOW.exists() else ""
        if current != rendered:
            print("ci.yml is stale: regenerate with tools/conformance/generate_ci.py")
            import difflib

            diff = difflib.unified_diff(
                current.splitlines(), rendered.splitlines(), "committed", "generated", lineterm=""
            )
            for line in list(diff)[:60]:
                print(line)
            return 1
        print(
            f"ci.yml matches the manifest: {len(manifest['suites'])} suites, "
            f"{sum(len(s['cases']) for s in manifest['suites'])} cases, "
            f"{len(manifest.get('probes', []))} probe(s)"
        )
        return 0

    WORKFLOW.write_text(rendered)
    print(f"wrote {WORKFLOW} ({rendered.count(chr(10)) + 1} lines)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

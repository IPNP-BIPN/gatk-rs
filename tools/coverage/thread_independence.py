#!/usr/bin/env python3
"""The port's Spark tools at one thread and at N, which must be the same bytes (#72).

A parallel path that is only concordant is not a claim this programme makes. This takes rows of each
Spark tool's covering array and runs the port twice on each, with `RAYON_NUM_THREADS=1` and
`RAYON_NUM_THREADS=N`, which governs bwa-mem4's pool and, through `bwa_threads`, the `-t` it is
given. What is compared is the exit code, stdout and the sha256 of every file written.

A tool with a parallel path runs every row of its array; the others, which spawn no thread to vary,
run a sample, so that a thread added later cannot arrive unobserved.

    python3 tools/coverage/thread_independence.py --port target/release/gatk
    python3 tools/coverage/thread_independence.py --port target/release/gatk --host --fixtures DIR

By default both runs are in the container, at the array's paths, with the binary CI builds. `--host`
runs the binary here with `/work` rewritten, for a laptop; a row that needs ISA-L refuses on both
sides there and counts as a refusal, so only the x86-64 run is the measurement.
"""

import argparse
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import run_array  # noqa: E402

# The tools that reach `bwa_threads`, directly or through `bwa_align`.
PARALLEL = {
    "BwaSpark",
    "BwaAndMarkDuplicatesPipelineSpark",
    "PathSeqBwaSpark",
    "PathSeqFilterSpark",
    "PathSeqPipelineSpark",
    "ReadsPipelineSpark",
}


# The arrays' BWA fixtures hold a few reads, too few to fill a batch per thread: this holds thousands.
LARGER = (
    "--input /work/fixtures/svbp_reads.bam --reference /work/fixtures/svbp_sv.fasta "
    "--bwa-mem-index-image /work/fixtures/svbp_sv.fasta.img --output /work/out/larger.bam "
    "--single-end-alignment true --use-jdk-deflater true --spark-master local[1]"
)


def spark_arrays():
    """Every array of a tool whose arguments name a Spark master, which is what a Spark tool has."""
    found = []
    for path in sorted(run_array.ARRAYS.glob("*.t2.json")):
        array = json.loads(path.read_text())
        rows = array.get("array", [])
        if path.name.split(".")[0].endswith("Spark") and rows:
            found.append((array["tool"], array))
    return found


def sample(rows, count):
    """`count` rows spread over the array, so a sample is not only its first few assignments."""
    if len(rows) <= count:
        return rows
    return [rows[i * len(rows) // count] for i in range(count)]


def digest(path):
    """A file's digest, the wall clock a metrics file opens with taken out."""
    raw = path.read_bytes()
    if raw.startswith(b"## htsjdk.samtools.metrics"):
        raw = run_array.without_start_time(raw.decode("utf-8")).encode()
    return hashlib.sha256(raw).hexdigest()


def tree(directory):
    """Every file under a directory, by relative name, as its digest."""
    return {
        path.relative_to(directory).as_posix(): digest(path)
        for path in sorted(directory.rglob("*"))
        if path.is_file()
    }


def in_container(binary, tool, cli, threads, work, fixtures):
    """One run at the array's own paths. `timeout` is inside, so a hung run frees its container."""
    out = work / "out"
    command = (
        "rm -rf /work/out/* && mkdir -p /work/tmp /work/tmp2 /work/out && "
        f"timeout 300 /work/port-binary/{binary.name} {tool} {cli}"
    )
    result = subprocess.run(
        [
            "docker", "run", "--rm", "--platform", run_array.PLATFORM,
            "-e", f"RAYON_NUM_THREADS={threads}",
            "-v", f"{fixtures}:/work/fixtures:ro",
            "-v", f"{out}:/work/out",
            "-v", f"{binary.parent}:/work/port-binary:ro",
            "-w", "/work", run_array.IMAGE, command,
        ],
        capture_output=True,
        text=True,
    )
    return result.returncode, result.stdout, tree(out)


def on_host(binary, tool, cli, threads, work, fixtures):
    """One run here, with the container's `/work` under a scratch directory."""
    root = work / "host"
    shutil.rmtree(root, ignore_errors=True)
    for name in ("tmp", "tmp2", "out"):
        (root / name).mkdir(parents=True)
    (root / "fixtures").symlink_to(fixtures)
    words = [word.replace("/work/", f"{root}/") for word in shlex.split(cli)]
    result = subprocess.run(
        ["timeout", "300", str(binary), tool, *words],
        capture_output=True,
        text=True,
        cwd=root,
        env={**os.environ, "RAYON_NUM_THREADS": str(threads)},
    )
    return result.returncode, result.stdout.replace(str(root), "/work"), tree(root / "out")


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", required=True)
    parser.add_argument("--threads", type=int, default=4, help="the N compared with one thread")
    parser.add_argument("--rows", type=int, default=4, help="rows of a tool with no parallel path")
    parser.add_argument("--fixtures", help="a corpus already built; else one is built here")
    parser.add_argument("--host", action="store_true")
    parser.add_argument("--tool", action="append", help="only these tools")
    options = parser.parse_args(argv)
    binary = Path(options.port).resolve()
    run = on_host if options.host else in_container

    differing = []
    total = accepted = 0
    with tempfile.TemporaryDirectory() as scratch:
        work = Path(scratch)
        (work / "out").mkdir()
        if options.fixtures:
            fixtures = Path(options.fixtures).resolve()
        else:
            fixtures = run_array.build_fixtures(work)
        for tool, array in spark_arrays():
            if options.tool and tool not in options.tool:
                continue
            rows = array["array"] if tool in PARALLEL else sample(array["array"], options.rows)
            held = array.get("excluded", [])
            positional = run_array.positional_values(tool)
            tagged = run_array.tagged_arguments(tool)
            lists = run_array.list_arguments(tool)
            same = 0
            ran = 0
            for row in rows:
                args = run_array.row_arguments(row, held)
                cli = " ".join(run_array.as_cli(args, positional, tagged, lists))
                one = run(binary, tool, cli, 1, work, fixtures)
                many = run(binary, tool, cli, options.threads, work, fixtures)
                ran += 1
                total += 1
                accepted += 1 if one[0] == 0 else 0
                if one == many:
                    same += 1
                else:
                    differing.append((tool, row["row"], one, many))
            print(f"{tool}: {same}/{ran} equal at 1 and {options.threads} threads", flush=True)

        if not options.tool or "BwaSpark" in options.tool:
            one = run(binary, "BwaSpark", LARGER, 1, work, fixtures)
            many = run(binary, "BwaSpark", LARGER, options.threads, work, fixtures)
            total += 1
            accepted += 1 if one[0] == 0 else 0
            if one != many:
                differing.append(("BwaSpark", "larger", one, many))
            print(f"BwaSpark, larger input: exit {one[0]}, equal {one == many}", flush=True)

        # What a container wrote belongs to root, and only a container can take it away again.
        if not (options.host and options.fixtures):
            subprocess.run(
                [
                    "docker", "run", "--rm", "--platform", run_array.PLATFORM,
                    "-v", f"{work}:/scratch", run_array.IMAGE, "rm -rf /scratch/*",
                ],
                capture_output=True,
            )

    print(f"tools rows={total} accepted={accepted} differing={len(differing)}")
    for tool, row, one, many in differing:
        changed = sorted(
            name for name in set(one[2]) | set(many[2]) if one[2].get(name) != many[2].get(name)
        )
        print(f"DIFFERS {tool} row {row}: exit {one[0]} vs {many[0]}, files {changed}")
    return 1 if differing else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

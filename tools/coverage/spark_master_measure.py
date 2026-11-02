#!/usr/bin/env python3
"""What the REFERENCE's bytes do when the master, the split size or the reducers change (#72).

Each sampled tool runs in the oracle on a BAM large enough to split, once per setting, and each run
is compared with the same tool at `local[1]`. A second `local[1]` run is the control, and a second
`local[4]` run says whether `local[4]` agrees with itself.

A BAM is compared twice: as its raw bytes, and as its decompressed text and records with the
command line taken out of the `@PG`, because the setting itself is on that line. An output file
is otherwise compared by run_array's rendering, which already takes out what no two runs share.

    python3 tools/coverage/spark_master_measure.py --out measured.tsv
"""

import argparse
import concurrent.futures
import gzip
import hashlib
import re
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import run_array  # noqa: E402

BAM = "/work/fixtures/svbp_reads.bam"
REFERENCE = "/work/fixtures/svbp_sv.fasta"
SMALL_SPLIT = "65536"

# Tool, then its arguments apart from the settings measured.
TOOLS = {
    "CountReadsSpark": ["--input", BAM],
    "CountBasesSpark": ["--input", BAM],
    "FlagStatSpark": ["--input", BAM],
    "PrintReadsSpark": ["--input", BAM, "--output", "/work/out/out.bam"],
    "SortSamSpark": ["--input", BAM, "--output", "/work/out/out.bam", "--sort-order", "queryname"],
    "MarkDuplicatesSpark": ["--input", BAM, "--output", "/work/out/out.bam"],
    "MeanQualityByCycleSpark": ["--input", BAM, "--output", "/work/out/out.txt"],
    "QualityScoreDistributionSpark": ["--input", BAM, "--output", "/work/out/out.txt"],
    "CollectInsertSizeMetricsSpark": [
        "--input", BAM, "--output", "/work/out/out.txt",
        "--histogram-output", "/work/out/hist.pdf", "--produce-plot", "false",
    ],
    "CollectQualityYieldMetricsSpark": ["--input", BAM, "--output", "/work/out/out.txt"],
    "PileupSpark": ["--input", BAM, "--output", "/work/out/pileup", "--reference", REFERENCE],
    "FindBadGenomicKmersSpark": ["--reference", REFERENCE, "--output", "/work/out/out.txt"],
}

SETTINGS = {
    "local[1]-again": ["--spark-master", "local[1]"],
    "local[4]": ["--spark-master", "local[4]"],
    "local[4]-again": ["--spark-master", "local[4]"],
    "split-64k": ["--spark-master", "local[1]", "--bam-partition-size", SMALL_SPLIT],
    "split-64k-local[4]": ["--spark-master", "local[4]", "--bam-partition-size", SMALL_SPLIT],
    "reducers-3": ["--spark-master", "local[1]", "--num-reducers", "3"],
}


def bgzf_payload(raw):
    """A BAM's bytes with its compression taken off and its `@PG` command lines blanked."""
    data = b""
    position = 0
    while position < len(raw):
        decoder = gzip.zlib.decompressobj(31)
        data += decoder.decompress(raw[position:])
        position = len(raw) - len(decoder.unused_data)
        if not decoder.unused_data:
            break
    if not data.startswith(b"BAM\x01"):
        return data
    length = struct.unpack("<i", data[4:8])[0]
    text = re.sub(rb"\tCL:[^\t\n]*", b"", data[8 : 8 + length])
    return text + data[8 + length :]


def render(out_dir):
    """Name to `(raw digest, comparable digest)`, for every file a run left."""
    rendered = {}
    for path in sorted(out_dir.rglob("*")):
        if not path.is_file() or path.name.startswith("."):
            continue
        raw = path.read_bytes()
        name = path.relative_to(out_dir).as_posix()
        if path.suffix == ".bam":
            comparable = bgzf_payload(raw)
        else:
            try:
                comparable = run_array.without_start_time(raw.decode("utf-8")).encode()
            except UnicodeDecodeError:
                comparable = raw
        rendered[name] = (
            hashlib.sha256(raw).hexdigest(),
            hashlib.sha256(comparable).hexdigest(),
        )
    return rendered


def run(tool, arguments, fixtures, out_dir):
    """One run in the container, bounded from inside. `(exit, stdout, files)`."""
    out_dir.mkdir(parents=True)
    quoted = " ".join(run_array.shlex.quote(word) for word in arguments)
    command = (
        "mkdir -p /work/tmp /work/out && timeout 600 "
        f'java {run_array.JAVA_OPENS} -cp "$ORACLE_CP" org.broadinstitute.hellbender.Main '
        f"{tool} {quoted}; code=$?; chmod -R a+rwX /work/out; exit $code"
    )
    result = subprocess.run(
        [
            "docker", "run", "--rm", "--platform", run_array.PLATFORM,
            "-v", f"{fixtures}:/work/fixtures:ro",
            "-v", f"{out_dir}:/work/out",
            "-w", "/work", run_array.IMAGE, command,
        ],
        capture_output=True,
        text=True,
    )
    return result.returncode, result.stdout, render(out_dir), run_array.first_error(result.stderr)


def verdict(base, other):
    if base[0] != other[0]:
        return f"exit {base[0]} vs {other[0]}: {other[3][:100]}"
    if base[0] != 0:
        return "both refused" if base[3] == other[3] else f"refused differently: {other[3][:100]}"
    if base[1] != other[1]:
        return "stdout differs"
    if set(base[2]) != set(other[2]):
        return f"files differ: {sorted(set(base[2]) ^ set(other[2]))[:4]}"
    if all(base[2][name][0] == other[2][name][0] for name in base[2]):
        return "same bytes"
    changed = [name for name in base[2] if base[2][name][1] != other[2][name][1]]
    if not changed:
        return "same content, other bytes"
    return f"content differs: {changed[:4]}"


def measure(tool, arguments, fixtures, scratch):
    base = run(tool, arguments + ["--spark-master", "local[1]"], fixtures, scratch / tool / "base")
    lines = [f"{tool}\tlocal[1]\texit {base[0]}, {len(base[2])} files {base[3][:80]}"]
    for label, setting in SETTINGS.items():
        other = run(tool, arguments + setting, fixtures, scratch / tool / label.replace("[", "_"))
        lines.append(f"{tool}\t{label}\t{verdict(base, other)}")
    print("\n".join(lines), flush=True)
    return lines


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True)
    parser.add_argument("--fixtures", required=True, help="a corpus MakeFixtures built")
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--tool", action="append")
    options = parser.parse_args(argv)
    fixtures = Path(options.fixtures).resolve()
    tools = {t: a for t, a in TOOLS.items() if not options.tool or t in options.tool}
    lines = []
    with tempfile.TemporaryDirectory(dir=Path(options.out).resolve().parent) as scratch:
        with concurrent.futures.ThreadPoolExecutor(options.jobs) as pool:
            for result in pool.map(
                lambda item: measure(item[0], item[1], fixtures, Path(scratch)), tools.items()
            ):
                lines += result
    Path(options.out).write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

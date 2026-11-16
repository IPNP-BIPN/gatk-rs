#!/usr/bin/env python3
"""The most a device inflater could save a read-bound tool: the tool on a BAM, then on the same
BAM with stored blocks, whose inflate is a copy. Milestone GPU, GPU.2 (docs/gpu/targets.md, §3).

    python3 tools/speed/bgzf_inflate_bound.py --port target/release/gatk --bam in.bam --work DIR

From one BAM it writes `deflated.bam` (its records `--copies` times behind one header, level 5 as
htsjdk writes) and `stored.bam` (the same bytes, level 0), checks the tool counts the same reads on
both, and prints the best wall and CPU time of `--runs` runs each, at one rayon thread and at all.
The difference at all threads bounds what any inflater, on any device, can take off the tool.
"""

import argparse
import gzip
import os
import resource
import struct
import subprocess
import time
import zlib

BLOCK = 65280


def bgzf(path, data, level):
    """`data` as BGZF blocks at one deflate level, then the empty block BGZF ends with."""
    with open(path, "wb") as out:
        for chunk in [data[i : i + BLOCK] for i in range(0, len(data), BLOCK)] + [b""]:
            compressor = zlib.compressobj(level, zlib.DEFLATED, -15)
            body = compressor.compress(chunk) + compressor.flush()
            header = struct.pack("<BBBBIBBHBBHH", 31, 139, 8, 4, 0, 0, 255, 6, 66, 67, 2, len(body) + 25)
            out.write(header + body + struct.pack("<II", zlib.crc32(chunk), len(chunk)))


def header_length(raw):
    assert raw[:4] == b"BAM\1", "not a BAM"
    at = 8 + struct.unpack_from("<i", raw, 4)[0]
    references = struct.unpack_from("<i", raw, at)[0]
    at += 4
    for _ in range(references):
        at += 4 + struct.unpack_from("<i", raw, at)[0] + 4
    return at


def timed(command, threads):
    environment = dict(os.environ, RAYON_NUM_THREADS=str(threads))
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.perf_counter()
    out = subprocess.run(command, env=environment, capture_output=True, text=True, check=True).stdout
    wall = time.perf_counter() - start
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    return wall, after.ru_utime - before.ru_utime, out.strip().splitlines()[-1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--work", required=True)
    parser.add_argument("--tool", default="CountReads")
    parser.add_argument("--copies", type=int, default=20)
    parser.add_argument("--runs", type=int, default=3)
    arguments = parser.parse_args()

    raw = gzip.decompress(open(arguments.bam, "rb").read())
    split = header_length(raw)
    data = raw[:split] + raw[split:] * arguments.copies
    files = {name: os.path.join(arguments.work, f"{name}.bam") for name in ("deflated", "stored")}
    bgzf(files["deflated"], data, 5)
    bgzf(files["stored"], data, 0)

    cores = os.cpu_count()
    print(f"{len(data):,} bytes inflated; {arguments.tool}; best of {arguments.runs}")
    print("| file | size | threads | wall | CPU | answer |")
    print("|---|---:|---:|---:|---:|---|")
    answers = set()
    for name, path in files.items():
        for threads in (1, cores):
            runs = [timed([arguments.port, arguments.tool, "-I", path], threads) for _ in range(arguments.runs)]
            wall = min(r[0] for r in runs)
            cpu = min(r[1] for r in runs)
            answers.add(runs[0][2])
            print(f"| {name} | {os.path.getsize(path) / 1e6:.0f} MB | {threads} | {wall:.2f} s | {cpu:.2f} s | {runs[0][2]} |")
    assert len(answers) == 1, f"the two files did not give the same answer: {answers}"


if __name__ == "__main__":
    main()

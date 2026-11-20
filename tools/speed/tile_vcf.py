#!/usr/bin/env python3
"""Tile one VCF window across a contig, to time a VCF tool on real records at scale.

    python3 tools/speed/tile_vcf.py in.vcf WINDOW_START WINDOW_LENGTH COPIES out.vcf

Each copy shifts the window's records by WINDOW_LENGTH, starting at position 1, so the output stays
sorted on the same contig; REF is not re-read from the reference. A `.stats` beside the input
(Mutect2's) is written beside the output with `callable` multiplied by COPIES.
"""
import sys
from pathlib import Path

source, start, length, copies, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
lines = Path(source).read_text().splitlines(keepends=True)
header = [line for line in lines if line.startswith("#")]
records = [line.split("\t", 2) for line in lines if not line.startswith("#")]
with open(out, "w") as sink:
    sink.writelines(header)
    for copy in range(copies):
        for contig, position, rest in records:
            sink.write(f"{contig}\t{int(position) - start + 1 + copy * length}\t{rest}")
stats = Path(source + ".stats")
if stats.exists():
    table = [line.split("\t") for line in stats.read_text().splitlines()]
    scaled = [f"{k}\t{float(v) * copies}" if k == "callable" else f"{k}\t{v}" for k, v in table]
    Path(out + ".stats").write_text("\n".join(scaled) + "\n")

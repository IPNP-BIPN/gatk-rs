#!/usr/bin/env python3
"""A deterministic SAM text input for timing CountReads: one 10 Mb contig, single-end 150 bp reads
copied from it with a few substitutions, coordinate sorted, one read group. Seeded, so every runner
times the same bytes. Usage: make_input.py <out-dir> <reads>"""
import random, sys
from pathlib import Path

out, reads = Path(sys.argv[1]), int(sys.argv[2])
rng = random.Random(20261120)
length, read_len = 10_000_000, 150
ref = bytes(rng.choice(b"ACGT") for _ in range(length))
out.mkdir(parents=True, exist_ok=True)
with open(out / "ref.fasta", "wb") as fa:
    fa.write(b">chr1\n")
    for i in range(0, length, 60):
        fa.write(ref[i : i + 60] + b"\n")
starts = sorted(rng.randrange(1, length - read_len) for _ in range(reads))
quals = bytes(range(ord("#"), ord("J")))
with open(out / "reads.sam", "w") as sam:
    sam.write(f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{length}\n@RG\tID:0\tSM:s\tPL:ILLUMINA\n")
    for n, start in enumerate(starts):
        bases = bytearray(ref[start - 1 : start - 1 + read_len])
        for _ in range(rng.randrange(3)):
            bases[rng.randrange(read_len)] = rng.choice(b"ACGT")
        q = bytes(rng.choice(quals) for _ in range(read_len)).decode()
        flag = 16 if n & 1 else 0
        sam.write(f"r{n}\t{flag}\tchr1\t{start}\t60\t{read_len}M\t*\t0\t0\t{bases.decode()}\t{q}\tRG:Z:0\n")

#!/usr/bin/env python3
"""A deterministic SAM text input for timing CountReads, shaped like a 30x WGS BAM: one 10 Mb
contig, 150 bp mate pairs copied from it with a few substitutions, the tags a BWA + MarkDuplicates
record carries (MC, PG, MQ, AS, XS, MD, NM, RG) with MD/NM computed from the substitutions,
coordinate sorted. Seeded, so every runner times the same bytes. Usage: make_input.py <out-dir> <reads>"""
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

# Pairs: the first mate forward, the second reverse, an insert of about 400.
records = []
for n in range(reads // 2):
    start = rng.randrange(1, length - 1000)
    mate = start + rng.randrange(150, 450)
    tlen = mate + read_len - start
    records.append((start, n, 99, mate, tlen))
    records.append((mate, n, 147, start, -tlen))
records.sort()
quals = bytes(range(ord("#"), ord("J")))
group = "NA12878_TTGCCTAG-ACCACTTA_HCLHLDSXX_L003"
with open(out / "reads.sam", "w") as sam:
    sam.write(f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{length}\n")
    sam.write(f"@RG\tID:{group}\tSM:NA12878\tPL:ILLUMINA\tLB:lib1\n@PG\tID:MarkDuplicates\tPN:MarkDuplicates\n")
    for start, n, flag, mate, tlen in records:
        bases = bytearray(ref[start - 1 : start - 1 + read_len])
        subs = sorted({rng.randrange(read_len) for _ in range(rng.randrange(3))})
        for at in subs:
            bases[at] = rng.choice([b for b in b"ACGT" if b != bases[at]])
        md, run = "", 0
        for i in range(read_len):
            if i in subs:
                md += f"{run}{chr(ref[start - 1 + i])}"
                run = 0
            else:
                run += 1
        md += str(run)
        q = bytes(rng.choice(quals) for _ in range(read_len)).decode()
        name = f"HCLHLDSXX:3:{1101 + n % 50}:{n % 32000}:{n // 32000}"
        score = read_len - 5 * len(subs)
        sam.write(
            f"{name}\t{flag}\tchr1\t{start}\t60\t{read_len}M\t=\t{mate}\t{tlen}\t{bases.decode()}\t{q}\t"
            f"MC:Z:{read_len}M\tPG:Z:MarkDuplicates\tMQ:i:60\tAS:i:{score}\tXS:i:{rng.randrange(20, 60)}\t"
            f"MD:Z:{md}\tNM:i:{len(subs)}\tRG:Z:{group}\n"
        )

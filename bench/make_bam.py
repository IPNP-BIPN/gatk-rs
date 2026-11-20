"""A coordinate-sorted SAM of random 150 bp pairs on one 60 Mb contig, on stdout. Seeded."""
import random, sys
rng = random.Random(20261120)
n, length, contig = int(sys.argv[1]), 150, 60_000_000
print("@HD\tVN:1.6\tSO:coordinate"); print(f"@SQ\tSN:chrS\tLN:{contig}"); print("@RG\tID:rg\tSM:s")
step = (contig - 1000) / n
q = "".join(chr(33 + q) for q in range(2, 41))
for i in range(n):
    pos = 1 + int(i * step)
    seq = "".join(rng.choices("ACGT", k=length))
    qual = "".join(rng.choices(q, k=length))
    print(f"r{i}\t0\tchrS\t{pos}\t60\t{length}M\t*\t0\t0\t{seq}\t{qual}\tRG:Z:rg")

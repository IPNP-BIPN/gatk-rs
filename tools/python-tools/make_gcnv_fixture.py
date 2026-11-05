#!/usr/bin/env python3
"""Read counts for the gCNV measurement: a few samples over two autosomes and an X.

Seeded, and free of numpy, so the host and the image write the same bytes. Sample 0 carries a
duplication on chr1 and sample 1 a deletion on chr2, so that the calls have something to say.

Usage: make_gcnv_fixture.py <directory> [samples] [intervals per contig]
"""
import math
import random
import sys


CONTIGS = [("chr1", 2000000), ("chr2", 1500000), ("chrX", 1200000)]


def poisson(rng, mean):
    # Knuth's method, in logs so that a large mean does not underflow.
    if mean < 30:
        limit, k, p = math.exp(-mean), 0, 1.0
        while True:
            p *= rng.random()
            if p <= limit:
                return k
            k += 1
    return max(0, int(round(rng.gauss(mean, math.sqrt(mean)))))


def main():
    directory = sys.argv[1]
    samples = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    per_contig = int(sys.argv[3]) if len(sys.argv) > 3 else 60
    rng = random.Random(20261102)
    dictionary = "@HD\tVN:1.6\n" + "".join(
        "@SQ\tSN:%s\tLN:%d\n" % contig for contig in CONTIGS
    )
    intervals = [
        (contig, 1000 + i * 10000, 1000 + i * 10000 + 999)
        for contig, _ in CONTIGS
        for i in range(per_contig)
    ]
    with open(directory + "/intervals.interval_list", "w") as f:
        f.write(dictionary)
        for contig, start, end in intervals:
            f.write("%s\t%d\t%d\t+\t.\n" % (contig, start, end))
    # An interval's own bias, the same in every sample.
    bias = [math.exp(rng.gauss(0.0, 0.2)) for _ in intervals]
    for sample in range(samples):
        depth = 100.0 * math.exp(rng.gauss(0.0, 0.1))
        female = sample % 2 == 0
        with open("%s/sample%d.counts.tsv" % (directory, sample), "w") as f:
            f.write(dictionary)
            f.write("@RG\tID:GATKCopyNumber\tSM:sample%d\n" % sample)
            f.write("CONTIG\tSTART\tEND\tCOUNT\n")
            for index, (contig, start, end) in enumerate(intervals):
                copies = 2
                if contig == "chrX" and not female:
                    copies = 1
                position = index - (intervals.index((contig, 1000, 1999)))
                if sample == 0 and contig == "chr1" and 10 <= position < 20:
                    copies = 3
                if sample == 1 and contig == "chr2" and 30 <= position < 40:
                    copies = 1
                f.write(
                    "%s\t%d\t%d\t%d\n"
                    % (contig, start, end, poisson(rng, depth * bias[index] * copies / 2.0))
                )
    with open(directory + "/priors.tsv", "w") as f:
        f.write("CONTIG_NAME\tPLOIDY_PRIOR_0\tPLOIDY_PRIOR_1\tPLOIDY_PRIOR_2\tPLOIDY_PRIOR_3\n")
        f.write("chr1\t0.01\t0.01\t0.97\t0.01\n")
        f.write("chr2\t0.01\t0.01\t0.97\t0.01\n")
        f.write("chrX\t0.01\t0.49\t0.49\t0.01\n")


if __name__ == "__main__":
    main()

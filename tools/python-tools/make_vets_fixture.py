#!/usr/bin/env python3
"""The sites the isolation-forest measurement and the suite run on.

A seeded stream of SNPs and indels whose annotations are drawn from two separated clusters, with a
known-good resource over the first cluster and a few values missing, so that the imputer has
something to do. Nothing here depends on numpy, so the host and the image write the same bytes.

Usage: make_vets_fixture.py <directory> [sites]
"""
import random
import sys


HEADER = (
    "##fileformat=VCFv4.2\n"
    '##FILTER=<ID=LowQual,Description="Low quality">\n'
    '##INFO=<ID=QD,Number=1,Type=Float,Description="QD">\n'
    '##INFO=<ID=MQ,Number=1,Type=Float,Description="MQ">\n'
    '##INFO=<ID=FS,Number=1,Type=Float,Description="FS">\n'
    '##INFO=<ID=SOR,Number=1,Type=Float,Description="SOR">\n'
    "##contig=<ID=chr1,length=2000000>\n"
    "##contig=<ID=chr2,length=1500000>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def main():
    directory = sys.argv[1]
    sites = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
    rng = random.Random(20261102)
    calls = [HEADER]
    resource = [HEADER]
    position = 100
    for index in range(sites):
        contig = "chr1" if index < sites * 3 // 4 else "chr2"
        if index == sites * 3 // 4:
            position = 100
        position += rng.randint(5, 400)
        good = rng.random() < 0.6
        indel = rng.random() < 0.25
        ref, alt = rng.choice([("A", "C"), ("C", "T"), ("G", "A"), ("T", "G")])
        if indel:
            ref, alt = rng.choice([(ref + alt, ref), (ref, ref + alt)])
        qd = rng.gauss(20.0 if good else 6.0, 4.0)
        mq = rng.gauss(60.0 if good else 45.0, 3.0)
        fs = abs(rng.gauss(2.0 if good else 20.0, 3.0))
        sor = abs(rng.gauss(1.0 if good else 3.0, 0.5))
        info = ["QD=%.3f" % qd, "MQ=%.2f" % mq, "FS=%.3f" % fs, "SOR=%.3f" % sor]
        if rng.random() < 0.02:
            info.pop(2)
        record = "%s\t%d\t.\t%s\t%s\t50\t.\t%s\n" % (contig, position, ref, alt, ";".join(info))
        calls.append(record)
        if good and rng.random() < 0.8:
            resource.append("%s\t%d\t.\t%s\t%s\t50\t.\t.\n" % (contig, position, ref, alt))
    with open(directory + "/calls.vcf", "w") as f:
        f.writelines(calls)
    with open(directory + "/resource.vcf", "w") as f:
        f.writelines(resource)


if __name__ == "__main__":
    main()

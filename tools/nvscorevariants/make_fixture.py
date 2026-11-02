#!/usr/bin/env python3
"""The small corpus NVScoreVariants is compared on: a FASTA, a VCF, a BAM and its index.

Runs inside `broadinstitute/gatk:4.6.2.0`, which is what has pysam, so the reference and the port
read the same bytes and nothing here has to be committed. Everything comes from one linear
congruential generator, so the same bytes come out every time.

The cases are chosen for what `encoders.py` decides on rather than for realism: windows clipped by
the start and the end of a contig, IUPAC codes and lower case in the reference, reads with
insertions, deletions and clips, every flag `_fill_rows` reads, a pile deeper than the 128 reads the
reservoir keeps, an `artificial` read group, annotations absent, integer, floating and missing, and
every variant type `determine_polymorphic_type` can answer.

    make_fixture.py <output directory>
"""

import os
import sys

import pysam


class Rand:
    """A 48-bit linear congruential generator, so no library's generator is in the fixture."""

    def __init__(self, seed):
        self.state = (seed ^ 0x5DEECE66D) & ((1 << 48) - 1)

    def next(self, bound):
        self.state = (self.state * 0x5DEECE66D + 0xB) & ((1 << 48) - 1)
        return (self.state >> 17) % bound

    def choice(self, items):
        return items[self.next(len(items))]


def contig_sequence(rand, length, ambiguous):
    bases = [rand.choice("ACGT") for _ in range(length)]
    if ambiguous:
        for at in range(length // 2, length // 2 + 40, 3):
            bases[at] = rand.choice("NKMRYSWBVHDXU-")
        for at in range(length // 4, length // 4 + 60):
            bases[at] = bases[at].lower()
    return "".join(bases)


def write_fasta(path, contigs):
    with open(path, "w") as out:
        for name, sequence in contigs.items():
            out.write(f">{name} a fixture contig\n")
            for at in range(0, len(sequence), 60):
                out.write(sequence[at : at + 60] + "\n")


def cigar_of(rand, length):
    kind = rand.next(12)
    if kind == 0:
        return f"{length // 2}M3I{length - length // 2 - 3}M"
    if kind == 1:
        return f"{length // 3}M2D{length - length // 3}M"
    if kind == 2:
        return f"7S{length - 7}M"
    if kind == 3:
        return f"{length - 5}M5S"
    if kind == 4:
        return f"4H4S{length - 4}M"
    if kind == 5:
        return f"{length // 4}M1I{length // 4}M4D{length - length // 4 - length // 4 - 1}M"
    if kind == 6:
        return f"{length // 2}M40N{length - length // 2}M"
    if kind == 7:
        return f"3S{length // 2 - 3}M6I{length - length // 2 - 6}M"
    return f"{length}M"


def read_line(rand, name, contig, start, sequence, group, length):
    """One SAM line starting at the 0-based `start` of `sequence`."""
    cigar = cigar_of(rand, length)
    flag = 0
    if rand.next(2):
        flag |= 16
    if rand.next(10):
        flag |= 1
        if rand.next(2):
            flag |= 32
        flag |= rand.choice([64, 128, 0])
    if not rand.next(25):
        flag |= 1024
    if not rand.next(40):
        flag |= 256
    seq = []
    for at in range(length):
        base = sequence[(start + at) % len(sequence)].upper()
        if base not in "ACGT":
            base = "N"
        if not rand.next(30):
            base = rand.choice("ACGTN")
        seq.append(base)
    qual = "".join(chr(33 + rand.next(61)) for _ in range(length))
    mapq = rand.choice([0, 1, 20, 37, 59, 60, 60, 60])
    return (
        f"{name}\t{flag}\t{contig}\t{start + 1}\t{mapq}\t{cigar}\t*\t0\t0\t{''.join(seq)}\t{qual}"
        f"\tRG:Z:{group}\n"
    )


def write_bam(path, contigs, rand):
    sam = path + ".sam"
    with open(sam, "w") as out:
        out.write("@HD\tVN:1.6\tSO:coordinate\n")
        for name, sequence in contigs.items():
            out.write(f"@SQ\tSN:{name}\tLN:{len(sequence)}\n")
        out.write("@RG\tID:rg1\tSM:S1\n@RG\tID:ArTiFiCiAl_haplotypes\tSM:S1\n")
        number = 0
        for name, sequence in contigs.items():
            for start in range(0, len(sequence) - 40, 23):
                depth = 1 + rand.next(4)
                for _ in range(depth):
                    number += 1
                    group = "ArTiFiCiAl_haplotypes" if not rand.next(30) else "rg1"
                    length = rand.choice([36, 76, 100, 101, 150, 151])
                    if start + length >= len(sequence):
                        length = 36
                    out.write(
                        read_line(rand, f"r{number}", name, start, sequence, group, length)
                    )
        # A pile deeper than the reservoir, on a stretch every read of which covers one site.
        hot = "c1"
        for offset in range(180):
            number += 1
            start = 1500 + offset % 40
            out.write(read_line(rand, f"h{number}", hot, start, contigs[hot], "rg1", 100))
    unsorted_bam = path + ".unsorted.bam"
    pysam.view("-b", "-o", unsorted_bam, sam, catch_stdout=False)
    pysam.sort("-o", path, unsorted_bam, catch_stdout=False)
    pysam.index(path)
    os.remove(sam)
    os.remove(unsorted_bam)


HEADER = """##fileformat=VCFv4.2
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">
##FORMAT=<ID=VF,Number=1,Type=Float,Description="A floating point per-sample value">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##INFO=<ID=FS,Number=1,Type=Float,Description="Fisher strand">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS mapping quality">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Mapping quality rank sum">
##INFO=<ID=QD,Number=1,Type=Float,Description="Quality by depth">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Read position rank sum">
##INFO=<ID=SOR,Number=1,Type=Float,Description="Strand odds ratio">
##INFO=<ID=NOTE,Number=1,Type=String,Description="Free text">
##INFO=<ID=SMALL,Number=1,Type=Float,Description="A value htslib prints through its exponent path">
##INFO=<ID=LONG,Number=1,Type=Float,Description="A value with more digits than its fast parser">
##contig=<ID=c1,length={c1}>
##contig=<ID=c2,length={c2}>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
"""


def annotations(rand):
    """Every annotation, each sometimes absent, at the precision each comes in from a caller."""
    parts = []
    if rand.next(5):
        parts.append(f"AC={rand.next(4)}")
        parts.append(f"AF={rand.next(1000) / 1000:.3f}")
        parts.append(f"AN={2 + rand.next(4)}")
    if not rand.next(3):
        parts.append("DB")
    if rand.next(8):
        parts.append(f"DP={rand.next(400)}")
    if rand.next(8):
        parts.append(f"FS={rand.next(60000) / 1000:.3f}")
    if rand.next(8):
        parts.append(f"MQ={rand.next(70000) / 1000:.2f}")
    if rand.next(4):
        parts.append(f"MQRankSum={(rand.next(10000) - 5000) / 1000:.3f}")
    if rand.next(8):
        parts.append(f"QD={rand.next(40000) / 1000:.2f}")
    if rand.next(4):
        parts.append(f"ReadPosRankSum={(rand.next(8000) - 4000) / 1000:.3f}")
    if rand.next(8):
        parts.append(f"SOR={rand.next(10000) / 1000:.3f}")
    if not rand.next(6):
        parts.append("NOTE=free_text")
    if not rand.next(4):
        parts.append(f"SMALL={rand.next(900) + 100}e-{rand.next(12)}")
    if not rand.next(4):
        parts.append(f"LONG={rand.next(1 << 30)}.{rand.next(1 << 30)}{rand.next(1 << 30)}")
    return ";".join(parts) if parts else "."


def genotypes(rand, alts):
    if not rand.next(10):
        return "GT", "./.", "0/0"
    alleles = 1 + alts
    ad = lambda: ",".join(str(rand.next(50)) for _ in range(alleles))
    pl = lambda: ",".join(str(rand.next(900)) for _ in range(alleles * (alleles + 1) // 2))
    vf = lambda: f"{rand.next(100000) / 977:.5f}"
    return (
        "GT:AD:DP:GQ:PL:VF",
        f"0/1:{ad()}:{rand.next(90)}:{rand.next(99)}:{pl()}:{vf()}",
        f"0|{alts}:{ad()}:{rand.next(90)}:{rand.next(99)}:{pl()}:{vf()}",
    )


def write_vcf(path, contigs, rand):
    """Site list: (contig, 1-based position, alternate alleles against the reference there)."""
    sites = []

    def add(contig, pos, shape):
        sites.append((contig, pos, shape))

    for pos in [1, 2, 3, 10, 40, 63, 64, 65, 66, 100]:
        add("c1", pos, "snp")
    add("c1", 5, "deletion")
    add("c1", 70, "insertion")
    for pos in range(180, 1400, 41):
        add("c1", pos, rand.choice(["snp", "snp", "snp", "deletion", "insertion", "mnp", "mixed",
                                    "multi", "star", "none", "multi_indel"]))
    for pos in range(1505, 1625, 17):
        add("c1", pos, "snp" if pos % 2 else "insertion")
    add("c1", len(contigs["c1"]) - 70, "snp")
    add("c1", len(contigs["c1"]) - 64, "snp")
    add("c1", len(contigs["c1"]) - 63, "deletion")
    add("c1", len(contigs["c1"]) - 10, "snp")
    add("c1", len(contigs["c1"]), "snp")
    for pos in range(90, len(contigs["c2"]) - 100, 67):
        add("c2", pos, rand.choice(["snp", "snp", "insertion", "deletion"]))

    with open(path, "w") as out:
        out.write(HEADER.format(c1=len(contigs["c1"]), c2=len(contigs["c2"])))
        for contig, pos, shape in sorted(sites, key=lambda s: (s[0], s[1])):
            sequence = contigs[contig]
            ref = sequence[pos - 1].upper()
            other = "A" if ref != "A" else "C"
            third = "G" if ref != "G" else "T"
            if shape == "snp":
                ref_field, alt = ref, other
            elif shape == "deletion":
                ref_field, alt = sequence[pos - 1 : pos + 3].upper(), ref
            elif shape == "insertion":
                ref_field, alt = ref, ref + "TTG"
            elif shape == "mnp":
                ref_field = sequence[pos - 1 : pos + 1].upper()
                alt = other + third
            elif shape == "mixed":
                ref_field, alt = ref, f"{other},{ref}TT"
            elif shape == "multi":
                ref_field, alt = ref, f"{other},{third}"
            elif shape == "multi_indel":
                ref_field, alt = sequence[pos - 1 : pos + 2].upper(), f"{ref},{ref}A"
            elif shape == "star":
                ref_field, alt = ref, f"{other},*"
            else:
                ref_field, alt = ref, "."
            alts = 0 if alt == "." else alt.count(",") + 1
            quality = rand.choice([".", "100", f"{rand.next(1000000) / 97:.2f}", "1234567.89", "0.00005"])
            status = rand.choice(["PASS", ".", "LowQual", "PASS"])
            identifier = rand.choice([".", ".", f"rs{pos}"])
            layout, one, two = genotypes(rand, max(alts, 1))
            if alts == 0:
                layout, one, two = "GT", "0/0", "./."
            out.write(
                f"{contig}\t{pos}\t{identifier}\t{ref_field}\t{alt}\t{quality}\t{status}\t"
                f"{annotations(rand)}\t{layout}\t{one}\t{two}\n"
            )


def write_symbolic(path, variants):
    """The same sites with one symbolic allele partway through, which ends the reference early."""
    with open(variants) as lines, open(path, "w") as out:
        record = 0
        for line in lines:
            if not line.startswith("#"):
                record += 1
                if record == 12:
                    fields = line.split("\t")
                    fields[4] = "<NON_REF>"
                    line = "\t".join(fields)
            out.write(line)


def main():
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    rand = Rand(47382911)
    contigs = {
        "c1": contig_sequence(rand, 2600, False),
        "c2": contig_sequence(rand, 1200, True),
    }
    write_fasta(f"{out}/reference.fasta", contigs)
    write_bam(f"{out}/reads.bam", contigs, rand)
    write_vcf(f"{out}/variants.vcf", contigs, rand)
    write_symbolic(f"{out}/symbolic.vcf", f"{out}/variants.vcf")


if __name__ == "__main__":
    main()

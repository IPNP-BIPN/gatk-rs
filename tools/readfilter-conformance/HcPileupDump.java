/*
 * The callers with pileup detection, end to end (G3.2g, #1336).
 *
 * --pileup-detection counts, at every locus of an active region's span, the substitutions (and
 * with --pileup-detection-enable-indel-pileup-calling the insertions and deletions) its pileup
 * holds. An event with a large enough share of the depth is "good": HaplotypeCaller's trimmer
 * sees it, and haplotypes carrying it are built from the five best assembled ones and the reads'
 * k-mers pick which are kept. An event whose reads are mostly "bad" for assembly is "bad": the
 * trimmer does not see it. A read is bad by its pairing, its supplementary alignments, its share
 * of mismatches (the MZ attribute, from NM or from the reference) and its template length.
 *
 * Three read sets: HcEndToEndDump's (unpaired, two samples), HcDragenEndToEndDump's (short tandem
 * repeats, for indels), and a paired set made here, over one contig: a SNP at 400 in a third of
 * the reads, a two-base insertion after 600 in half, a SNP at 800 carried only by reads with many
 * mismatches, some pairs improper, some reads supplementary, some with NM tags, template lengths
 * spread. A k-mer size longer than the reads stops the assembler finding anything, so the pileup
 * is the only source of the alternate haplotypes. Runs in ONE JVM, in this order.
 *
 * Output:
 *
 *     fasta\t<set>:<contig>=<the contig's bases>
 *     sam\t<set>=<the SAM text the BAM was written from, escaped>
 *     out\t<label>=<the output's ##INFO and ##FORMAT lines and its body from #CHROM on, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: HcPileupDump
 */

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.RampedHaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class HcPileupDump {

    static final int PAIRED_LENGTH = 1000;
    static final int READ = 80;

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("hc-pileup-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# HcPileupDump: the callers with pileup detection");

        // HcEndToEndDump's contigs and reads, from the same seed.
        final Random random = new Random(20261010L);
        final String[] references = new String[HcEndToEndDump.CONTIGS.length];
        for (int c = 0; c < HcEndToEndDump.CONTIGS.length; c++) {
            final StringBuilder ref = new StringBuilder();
            for (int i = 0; i < HcEndToEndDump.LENGTHS[c]; i++) {
                ref.append("ACGT".charAt(random.nextInt(4)));
            }
            references[c] = ref.toString();
            System.out.printf("fasta\te2e:%s=%s%n", HcEndToEndDump.CONTIGS[c], references[c]);
        }
        final Path e2eDir = dir.resolve("e2e");
        Files.createDirectories(e2eDir);
        final Path e2eFasta = HcEndToEndDump.writeReference(e2eDir, references);
        final Path e2eBam = HcLinkedGraphDump.bam(e2eDir, "e2e", HcEndToEndDump.sam(references, random));

        // HcDragenEndToEndDump's contig and reads.
        final Random strRandom = new Random(20261013L);
        final StringBuilder ref = new StringBuilder();
        for (int i = 0; i < HcDragenEndToEndDump.LENGTH; i++) {
            ref.append("ACGT".charAt(strRandom.nextInt(4)));
        }
        HcDragenEndToEndDump.splice(ref, 401, "AAAAAAAAAA");
        HcDragenEndToEndDump.splice(ref, 701, "CACACACACACACA");
        HcDragenEndToEndDump.splice(ref, 1001, "TTTTTTTT");
        ref.setCharAt(399, 'G');
        ref.setCharAt(410, 'G');
        ref.setCharAt(699, 'T');
        ref.setCharAt(714, 'T');
        ref.setCharAt(999, 'G');
        ref.setCharAt(1008, 'G');
        final String strReference = ref.toString();
        System.out.printf("fasta\tstr:%s=%s%n", HcDragenEndToEndDump.CONTIG, strReference);
        final Path strDir = dir.resolve("str");
        Files.createDirectories(strDir);
        final Path strFasta = HcDragenEndToEndDump.writeReference(strDir, strReference);
        final Path strBam = HcLinkedGraphDump.bam(strDir, "str", HcDragenEndToEndDump.sam(strReference));

        // The paired set.
        final Random pairedRandom = new Random(20261014L);
        final StringBuilder paired = new StringBuilder();
        for (int i = 0; i < PAIRED_LENGTH; i++) {
            paired.append("ACGT".charAt(pairedRandom.nextInt(4)));
        }
        final String pairedReference = paired.toString();
        System.out.printf("fasta\tpaired:chr1=%s%n", pairedReference);
        final Path pairedDir = dir.resolve("paired");
        Files.createDirectories(pairedDir);
        final Path pairedFasta = writeReference(pairedDir, pairedReference);
        final Path pairedBam = HcLinkedGraphDump.bam(pairedDir, "paired", pairedSam(pairedReference, pairedRandom));

        final String on = "--pileup-detection";
        final String indels = "--pileup-detection-enable-indel-pileup-calling";
        final List<String> noAssembly = List.of("--kmer-size", "101", "--dont-increase-kmer-sizes-for-cycles");
        final List<String> badness = List.of("--pileup-detection-bad-read-tolerance", "0.4",
                "--pileup-detection-filter-assembly-alt-bad-read-tolerance", "0.4");

        HcAllAnnotationsDump.run(dir, "hc-pileup", HaplotypeCaller::new, e2eBam, e2eFasta, List.of(on));
        HcAllAnnotationsDump.run(dir, "hc-pileup-indels", HaplotypeCaller::new, e2eBam, e2eFasta, List.of(on, indels));
        HcAllAnnotationsDump.run(dir, "hc-pileup-no-assembly", HaplotypeCaller::new, e2eBam, e2eFasta,
                concat(List.of(on, indels), noAssembly));
        HcAllAnnotationsDump.run(dir, "hc-pileup-threshold", HaplotypeCaller::new, e2eBam, e2eFasta,
                concat(List.of(on, indels, "--pileup-detection-active-region-phred-threshold", "3.0"), noAssembly));
        HcAllAnnotationsDump.run(dir, "hc-pileup-unpaired-bad", HaplotypeCaller::new, e2eBam, e2eFasta,
                concat(List.of(on, indels), concat(badness, noAssembly)));
        HcAllAnnotationsDump.run(dir, "hc-pileup-gvcf", HaplotypeCaller::new, e2eBam, e2eFasta,
                concat(List.of(on, indels, "-ERC", "GVCF", "--sample-name", "s1", "-L", "chr1"), noAssembly));
        HcAllAnnotationsDump.run(dir, "ramped-pileup", RampedHaplotypeCaller::new, e2eBam, e2eFasta,
                concat(List.of(on, indels), noAssembly));
        HcAllAnnotationsDump.run(dir, "mutect2-pileup", Mutect2::new, e2eBam, e2eFasta, List.of(on, indels));
        HcAllAnnotationsDump.run(dir, "mutect2-pileup-no-assembly", Mutect2::new, e2eBam, e2eFasta,
                concat(List.of(on, indels), noAssembly));
        HcAllAnnotationsDump.run(dir, "str-hc-pileup", HaplotypeCaller::new, strBam, strFasta,
                concat(List.of(on, indels), noAssembly));
        HcAllAnnotationsDump.run(dir, "str-mutect2-pileup", Mutect2::new, strBam, strFasta,
                concat(List.of(on, indels), noAssembly));

        HcAllAnnotationsDump.run(dir, "paired-hc", HaplotypeCaller::new, pairedBam, pairedFasta, List.of());
        HcAllAnnotationsDump.run(dir, "paired-hc-pileup", HaplotypeCaller::new, pairedBam, pairedFasta,
                List.of(on, indels));
        HcAllAnnotationsDump.run(dir, "paired-hc-pileup-no-assembly", HaplotypeCaller::new, pairedBam, pairedFasta,
                concat(List.of(on, indels), noAssembly));
        HcAllAnnotationsDump.run(dir, "paired-hc-badness", HaplotypeCaller::new, pairedBam, pairedFasta,
                concat(List.of(on, indels), badness));
        HcAllAnnotationsDump.run(dir, "paired-hc-badness-no-assembly", HaplotypeCaller::new, pairedBam, pairedFasta,
                concat(List.of(on, indels), concat(badness, noAssembly)));
        HcAllAnnotationsDump.run(dir, "paired-hc-badness-knobs", HaplotypeCaller::new, pairedBam, pairedFasta,
                concat(List.of(on, indels,
                        "--pileup-detection-proper-pair-read-badness", "false",
                        "--pileup-detection-chimeric-read-badness", "false",
                        "--pileup-detection-edit-distance-read-badness-threshold", "0.03",
                        "--pileup-detection-edit-distance-read-badness-for-assembly-filtering-threshold", "0.05",
                        "--pileup-detection-template-mean-badness-threshold", "300",
                        "--pileup-detection-template-std-badness-threshold", "40"), badness));
        HcAllAnnotationsDump.run(dir, "paired-hc-thresholds", HaplotypeCaller::new, pairedBam, pairedFasta,
                concat(List.of(on, indels,
                        "--pileup-detection-snp-alt-threshold", "0.4",
                        "--pileup-detection-indel-alt-threshold", "0.3",
                        "--pileup-detection-absolute-alt-depth", "30",
                        "--pileup-detection-snp-basequality-filter", "25",
                        "--pileup-detection-snp-adjacent-to-assembled-indel-range", "250",
                        "--num-artificial-haplotypes-to-add-per-allele", "1",
                        "--artifical-haplotype-filtering-kmer-size", "25"), noAssembly));
        HcAllAnnotationsDump.run(dir, "paired-hc-dragen-like", HaplotypeCaller::new, pairedBam, pairedFasta,
                concat(List.of(on, indels, "--pileup-detection-active-region-phred-threshold", "3.0"), badness));
        HcAllAnnotationsDump.run(dir, "paired-mutect2-pileup", Mutect2::new, pairedBam, pairedFasta,
                concat(List.of(on, indels), noAssembly));
        HcAllAnnotationsDump.run(dir, "paired-mutect2-badness", Mutect2::new, pairedBam, pairedFasta,
                concat(List.of(on, indels), badness));
    }

    static Path writeReference(final Path dir, final String reference) throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        final StringBuilder bases = new StringBuilder(">chr1\n");
        for (int i = 0; i < reference.length(); i += 60) {
            bases.append(reference, i, Math.min(i + 60, reference.length())).append("\n");
        }
        Files.writeString(fasta, bases.toString());
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        final htsjdk.samtools.SAMFileHeader header = new htsjdk.samtools.SAMFileHeader();
        header.setSequenceDictionary(new htsjdk.samtools.SAMSequenceDictionary(
                List.of(new htsjdk.samtools.SAMSequenceRecord("chr1", reference.length()))));
        try (final java.io.Writer writer = Files.newBufferedWriter(dir.resolve("reference.dict"))) {
            new htsjdk.samtools.SAMTextHeaderCodec().encode(writer, header);
        }
        return fasta;
    }

    static List<String> concat(final List<String> a, final List<String> b) {
        final List<String> out = new ArrayList<>(a);
        out.addAll(b);
        return out;
    }

    /** One fragment's two reads, as SAM rows sorted later. */
    static String pairedSam(final String reference, final Random random) {
        final List<long[]> keys = new ArrayList<>();
        final List<String> rows = new ArrayList<>();
        for (int f = 0; f < 140; f++) {
            final int start = 1 + random.nextInt(PAIRED_LENGTH - 480);
            final int insert = 200 + random.nextInt(200);
            final int mateStart = start + insert - READ;
            final boolean proper = f % 9 != 3;
            final boolean snp400 = f % 3 == 0;
            final boolean insertion600 = f % 2 == 0;
            final boolean noisy = f % 7 == 5;
            for (int end = 0; end < 2; end++) {
                final int pos = end == 0 ? start : mateStart;
                final StringBuilder bases = new StringBuilder(reference.substring(pos - 1, pos - 1 + READ));
                String cigar = READ + "M";
                int mismatches = 0;
                if (snp400 && pos <= 400 && 400 < pos + READ) {
                    final int i = 400 - pos;
                    bases.setCharAt(i, bases.charAt(i) == 'A' ? 'C' : 'A');
                    mismatches++;
                }
                if (noisy) {
                    // Every eleventh base changed, the SNP at 800 among them.
                    for (int i = 0; i < READ; i += 11) {
                        bases.setCharAt(i, bases.charAt(i) == 'G' ? 'T' : 'G');
                        mismatches++;
                    }
                    if (pos <= 800 && 800 < pos + READ) {
                        final int i = 800 - pos;
                        if (i % 11 != 0) {
                            bases.setCharAt(i, bases.charAt(i) == 'A' ? 'C' : 'A');
                            mismatches++;
                        }
                    }
                }
                if (insertion600 && pos + 10 <= 600 && 600 < pos + READ - 10) {
                    final int i = 600 - pos + 1;
                    bases.insert(i, "TG");
                    bases.setLength(READ);
                    cigar = i + "M2I" + (READ - i - 2) + "M";
                    mismatches += 2;
                }
                final StringBuilder quals = new StringBuilder();
                for (int i = 0; i < READ; i++) {
                    quals.append((char) (33 + (i % 13 == 7 ? 8 : 20 + random.nextInt(20))));
                }
                int flag = 1 | (proper ? 2 : 0) | (end == 0 ? 0x40 | 0x20 : 0x80 | 0x10);
                final int tlen = end == 0 ? insert : -insert;
                final StringBuilder tags = new StringBuilder("\tRG:Z:rg1");
                if (f % 10 == 4) {
                    tags.append("\tSA:Z:chr1,").append(pos).append(",+,40M40S,30,0;");
                }
                if (f % 4 == 1) {
                    tags.append("\tNM:i:").append(mismatches);
                }
                rows.add("f" + f + "\t" + flag + "\tchr1\t" + pos + "\t60\t" + cigar + "\t=\t"
                        + (end == 0 ? mateStart : start) + "\t" + tlen + "\t" + bases + "\t" + quals + tags);
                keys.add(new long[]{pos, rows.size() - 1});
            }
        }
        keys.sort((a, b) -> a[0] != b[0] ? Long.compare(a[0], b[0]) : Long.compare(a[1], b[1]));
        final StringBuilder out = new StringBuilder();
        out.append("@HD\tVN:1.6\tSO:coordinate\n");
        out.append("@SQ\tSN:chr1\tLN:").append(PAIRED_LENGTH).append("\n");
        out.append("@RG\tID:rg1\tSM:p1\tPL:ILLUMINA\n");
        for (final long[] key : keys) {
            out.append(rows.get((int) key[1])).append("\n");
        }
        return out.toString();
    }
}

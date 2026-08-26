/*
 * GeneExpressionEvaluation's counts, taken from the reference.
 *
 * RNA-seq fragments counted against gff3 features. What a read contributes is a WEIGHT rather than
 * a count, and the two ways of computing it disagree about what a read that half misses a gene is
 * worth.
 *
 * Nine behaviours this is built to catch.
 *
 *   - PROPORTIONAL COUNTS THE BASES A READ DOES NOT COVER: summedUnNormalizedWeights takes an
 *     extra `1.0 - totalCoveredBases/basesOnReference` before normalising, so a read half of which
 *     is intergenic gives its gene HALF a count, where EQUAL gives it a whole one;
 *   - EQUAL SPLITS BY FEATURE COUNT AND IGNORES OVERLAP LENGTH, so a read that clips one gene by a
 *     single base and covers another entirely gives each of them a half;
 *   - A GOOD PAIR IS COUNTED ONCE, FROM READ ONE: apply returns unless the read is first of pair or
 *     not in a good pair, and the pair's intervals are the union of both mates' alignment blocks;
 *   - A PAIR WITHOUT ITS MATE QUALITY IS A GATKException rather than a skipped read, because
 *     inGoodPair reads the MQ tag without checking it exists;
 *   - --multi-map-method EQUAL SILENTLY DROPS THE MAPPING QUALITY FILTER TO ZERO in
 *     onTraversalStart, so it changes which reads are counted and not only how they are weighted;
 *   - NH DECIDES MULTI-MAPPING AND ITS ABSENCE MEANS ONE, so a read with NH=3 is dropped entirely
 *     under IGNORE and contributes a third under EQUAL;
 *   - AN UNSTRANDED FEATURE EMITS ONE ROW, NOT TWO, and every read over it is called SENSE because
 *     isSense rewrites Strand.NONE as POSITIVE;
 *   - THE COUNT COLUMN IS NAMED AFTER THE SAMPLE, and the writer's column list is built from it,
 *     so the header of the table depends on the read group;
 *   - AND UNSPLICED MODE REPLACES THE ALIGNMENT BLOCKS WITH ONE INTERVAL, which for a good pair is
 *     the fragment length from the leftmost start, so a spliced read's intron is counted as
 *     covered.
 *
 * Output:
 *
 *     gff\t<the whole gff3 file, escaped>
 *     feature\t<label>\tcontig=<c>\tstart=<n>\tend=<n>\tstrand=<+,-,.>\ttype=<t>\
 *         \tid=<id>\tname=<name>
 *     overlap\t<grouping label>\t<contig>\t<start>\t<end>
 *     read\t<name>\tcontig=<c>\tstart=<n>\tcigar=<c>\tflags=<n>\tmq=<n>\tnh=<n>\
 *         \tmate-start=<n>\tmate-cigar=<c>\tmate-mq=<n>\tfragment=<n>
 *     counts\t<label>=<the whole tsv, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: GeneExpressionEvaluationDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.rnaseq.GeneExpressionEvaluation;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GeneExpressionEvaluationDump {

    /**
     * Two overlapping genes on the forward strand, one on the reverse, and one unstranded.
     *
     * geneA covers 100-299 in two exons with an intron between them, so a spliced read can cross
     * it and an unspliced interval cannot tell the difference. geneB starts inside geneA's second
     * exon, which is what gives a read two features to be split between.
     */
    static final String GFF = String.join("\n",
            "##gff-version 3",
            "##sequence-region chr1 1 10000",
            "chr1\ttest\tgene\t100\t299\t.\t+\t.\tID=gene_a;Name=geneA",
            "chr1\ttest\texon\t100\t149\t.\t+\t.\tID=exon_a1;Parent=gene_a",
            "chr1\ttest\texon\t250\t299\t.\t+\t.\tID=exon_a2;Parent=gene_a",
            "chr1\ttest\tgene\t260\t400\t.\t+\t.\tID=gene_b;Name=geneB",
            "chr1\ttest\texon\t260\t400\t.\t+\t.\tID=exon_b1;Parent=gene_b",
            "chr1\ttest\tgene\t1000\t1199\t.\t-\t.\tID=gene_c;Name=geneC",
            "chr1\ttest\texon\t1000\t1199\t.\t-\t.\tID=exon_c1;Parent=gene_c",
            "chr1\ttest\tgene\t2000\t2199\t.\t.\t.\tID=gene_d;Name=geneD",
            "chr1\ttest\texon\t2000\t2199\t.\t.\t.\tID=exon_d1;Parent=gene_d",
            "");

    record Read(String name, int start, String cigar, boolean reverse, Integer nh,
                Integer mateStart, String mateCigar, Integer mateMq, boolean firstOfPair,
                boolean paired, boolean properPair, boolean mateReverse, int fragmentLength) { }

    /**
     * A read that counts on its own, which still has to be PAIRED.
     *
     * inGoodPair calls mateIsUnmapped() before it asks anything else, and that throws on an
     * unpaired read, so a single-end BAM cannot be processed at all. What stands in for a single
     * read here is a paired read that is not PROPERLY paired: inGoodPair then answers false at the
     * second test and the read is counted alone.
     */
    static Read single(final String name, final int start, final String cigar,
                       final boolean reverse, final Integer nh) {
        return new Read(name, start, cigar, reverse, nh, 5000, "50M", 60, true, true, false, false,
                0);
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final Path dir = Path.of("gene-expression-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# GeneExpressionEvaluationDump: fragments counted against gff3 features");

        final Path gff = dir.resolve("features.gff3");
        Files.writeString(gff, GFF, StandardCharsets.UTF_8);
        new IndexFeatureFile().instanceMain(new String[] {"-I", gff.toString()});
        System.out.printf("gff\t%s%n", ReferenceQueryDump.escape(GFF));

        // Coordinate order, because the writer is told the input is presorted. The walker sees
        // them in this order too, which is the order the counts accumulate in.
        final List<Read> reads = sorted(List.of(
                // Wholly inside geneA's first exon: one feature, one whole count.
                single("inside-a", 100, "50M", false, null),
                // Half over geneA's second exon and half intergenic, which is where PROPORTIONAL
                // and EQUAL disagree.
                single("half-out", 280, "40M", false, null),
                // Over both geneA and geneB, clipping geneB by ten bases.
                single("two-genes", 250, "20M", false, null),
                // Spliced across geneA's intron, so its blocks miss the intron entirely.
                single("spliced", 130, "20M100N20M", false, null),
                // On the reverse strand over the forward geneA, which is antisense.
                single("antisense", 100, "50M", true, null),
                // Over the reverse geneC on the forward strand, which is antisense there too.
                single("over-c", 1000, "50M", false, null),
                // Over the unstranded geneD, which has no antisense row at all.
                single("over-d", 2000, "50M", false, null),
                // Three alignments, so IGNORE drops it and EQUAL gives a third.
                single("multi-map", 100, "50M", false, 3),
                // A proper pair over geneA, counted once from read one.
                new Read("pair", 100, "40M", false, null, 200, "40M", 60, true, true, true, true, 140),
                new Read("pair", 200, "40M", true, null, 100, "40M", 60, false, true, true, false, -140)));

        final Path bam = dir.resolve("reads.bam");
        writeBam(bam, reads, "sm1");
        describe(reads);

        run(dir, "default", bam, gff, List.of());
        run(dir, "equal-overlap", bam, gff, List.of("--multi-overlap-method", "EQUAL"));
        run(dir, "equal-multimap", bam, gff, List.of("--multi-map-method", "EQUAL"));
        run(dir, "unspliced", bam, gff, List.of("--unspliced", "true"));
        run(dir, "forward-forward", bam, gff, List.of("--read-strands", "FORWARD_FORWARD"));
        run(dir, "reverse-forward", bam, gff, List.of("--read-strands", "REVERSE_FORWARD"));
        run(dir, "label-id", bam, gff, List.of("--feature-label-key", "ID"));
        // Group by exon rather than by gene, so every exon is its own row and the overlap type
        // finds nothing beneath it.
        run(dir, "group-exon", bam, gff, List.of("--grouping-type", "exon"));

        // A pair whose mate quality was never written, which is a GATKException rather than a
        // skipped read.
        final Path noMq = dir.resolve("no-mq.bam");
        writeBam(noMq, sorted(List.of(
                new Read("pair", 100, "40M", false, null, 200, "40M", null, true, true, true, true, 140),
                new Read("pair", 200, "40M", true, null, 100, "40M", null, false, true, true, false, -140))),
                "sm1");
        run(dir, "no-mate-quality", noMq, gff, List.of());

        // A single-end read, which the tool cannot process at all: inGoodPair asks a read about
        // its mate before it asks whether it has one.
        final Path unpaired = dir.resolve("unpaired.bam");
        writeBam(unpaired, List.of(new Read("unpaired", 100, "50M", false, null, null, null, null,
                false, false, false, false, 0)), "sm1");
        run(dir, "unpaired", unpaired, gff, List.of());
    }

    /** Coordinate order, stable, which is what the writer and the walker both require. */
    static List<Read> sorted(final List<Read> reads) {
        final List<Read> copy = new ArrayList<>(reads);
        copy.sort(java.util.Comparator.comparingInt(Read::start));
        return copy;
    }

    static void describe(final List<Read> reads) {
        for (final Read read : reads) {
            System.out.printf(
                    "read\t%s\tcontig=chr1\tstart=%d\tcigar=%s\treverse=%s\tnh=%s"
                            + "\tpaired=%s\tproper=%s\tfirst=%s\tmate-start=%s"
                            + "\tmate-cigar=%s\tmate-mq=%s"
                            + "\tmate-reverse=%s\tfragment=%d%n",
                    read.name(), read.start(), read.cigar(), read.reverse(),
                    read.nh() == null ? "none" : read.nh(),
                    read.paired(), read.properPair(), read.firstOfPair(),
                    read.mateStart() == null ? "none" : read.mateStart(),
                    read.mateCigar() == null ? "none" : read.mateCigar(),
                    read.mateMq() == null ? "none" : read.mateMq(),
                    read.mateReverse(), read.fragmentLength());
        }
    }

    static void run(final Path dir, final String label, final Path bam, final Path gff,
                    final List<String> extra) throws Exception {
        final Path out = dir.resolve("counts-" + label + ".tsv");
        final List<String> argv = new ArrayList<>(List.of(
                "-I", bam.toString(), "-G", gff.toString(), "-O", out.toString()));
        argv.addAll(extra);
        try {
            new GeneExpressionEvaluation().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            return;
        }
        if (Files.exists(out)) {
            System.out.printf("counts\t%s=%s%n", label,
                    ReferenceQueryDump.escape(masked(Files.readString(out), dir)));
        }
    }

    static void writeBam(final Path file, final List<Read> reads, final String sample) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 10000))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample(sample);
        group.setPlatform("ILLUMINA");
        header.addReadGroup(group);
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .makeBAMWriter(header, true, file.toFile())) {
            for (final Read read : reads) {
                writer.addAlignment(record(header, read));
            }
        }
    }

    static SAMRecord record(final SAMFileHeader header, final Read read) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(read.name());
        record.setReferenceName("chr1");
        record.setAlignmentStart(read.start());
        record.setCigarString(read.cigar());
        record.setMappingQuality(60);
        record.setReadNegativeStrandFlag(read.reverse());
        final int length = length(read.cigar());
        final byte[] bases = new byte[length];
        Arrays.fill(bases, (byte) 'A');
        record.setReadBases(bases);
        final byte[] quals = new byte[length];
        Arrays.fill(quals, (byte) 30);
        record.setBaseQualities(quals);
        record.setAttribute("RG", "rg1");
        if (read.nh() != null) {
            record.setAttribute("NH", read.nh());
        }
        if (read.paired()) {
            record.setReadPairedFlag(true);
            record.setProperPairFlag(read.properPair());
            record.setFirstOfPairFlag(read.firstOfPair());
            record.setSecondOfPairFlag(!read.firstOfPair());
            record.setMateReferenceName("chr1");
            record.setMateAlignmentStart(read.mateStart());
            record.setMateNegativeStrandFlag(read.mateReverse());
            record.setInferredInsertSize(read.fragmentLength());
            record.setAttribute("MC", read.mateCigar());
            if (read.mateMq() != null) {
                record.setAttribute("MQ", read.mateMq());
            }
        }
        return record;
    }

    /** The read length a cigar implies, which is every operator that consumes read bases. */
    static int length(final String cigar) {
        return htsjdk.samtools.TextCigarCodec.decode(cigar).getReadLength();
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}

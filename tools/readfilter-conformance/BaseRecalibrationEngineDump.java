/*
 * BaseRecalibrationEngine, taken from the reference.
 *
 * The counting pass of BQSR: for every read, mark which bases disagree with the reference, weigh
 * each disagreement by the BAQ array, skip the known sites, and increment the datums the covariates
 * point at. It is the other half of the cycle ApplyBQSR closes.
 *
 * Nine behaviours this is built to catch.
 *
 *   - BAQ IS OFF BY DEFAULT. `recalArgs.enableBAQ` is false, so the default run uses flatBAQArray,
 *     an array of the constant 64, and the whole hidden Markov model is skipped. It is also skipped
 *     when the read has NO errors at all, whatever the flag says, "for efficiency reasons";
 *   - THE READ IS TRANSFORMED BEFORE IT IS COUNTED: the cigar is consolidated, default qualities are
 *     filled in, original qualities are restored, and then the adaptor and the SOFT CLIPS are hard
 *     clipped away. So the read the covariates see is not the read in the file;
 *   - AN INSERTION IS MARKED AT A DIFFERENT BASE ON EACH STRAND. Forward marks `readPos - 1` before
 *     advancing, reverse marks `readPos` after: the base on the far side of the insertion. A
 *     deletion is the same, `readPos - 1` forward and `readPos` reverse;
 *   - AND THE MARK IS CLAMPED AWAY AT THE ENDS rather than wrapping: `1D3M` and `3M1D` mark nothing;
 *   - THE FRACTIONAL ERROR ARRAY SPREADS EACH ERROR OVER A BLOCK, and the block starts ONE BASE
 *     BEFORE the first uncertain one: `Math.max(0, blockStartIndex - 1)`. With a flat BAQ array
 *     there are no blocks at all and the fractions are the marks themselves;
 *   - A KNOWN SITE IS SKIPPED BY READ OFFSET, NOT BY REFERENCE POSITION, and the conversion goes
 *     through getReadIndexForReferenceCoordinate, whose deletion case steps BACK one base;
 *   - AND SO IS EVERY BASE BELOW PRESERVE_QSCORES_LESS_THAN AND EVERY NON-REGULAR BASE;
 *   - THE READ GROUP TABLE IS NOT COUNTED, IT IS COLLAPSED. `finalizeData` marginalises the quality
 *     score table over the reported quality, which is the only place RecalDatum.combine runs in
 *     BQSR, so the read group's reported quality is an ESTIMATE and not an average;
 *   - AND THE TABLE IS ROUNDED TO MATCH WHAT A FILE WOULD HOLD, through
 *     `Math.round((in + Math.ulp(in)) * 10^n) / 10^n`. The ulp is added before the rounding, which
 *     is not the same as rounding twice, and it exists so that a table kept in memory equals one
 *     written and read back.
 *
 * Output:
 *
 *     const\t<name>\t<value>
 *     events\t<label>\t<nErrors>\t<isSNP>\t<isInsertion>\t<isDeletion>
 *     fractional\t<label>\t<comma separated bits>
 *     skip\t<label>\t<comma separated booleans>
 *     round\t<value>\t<places>\t<bits>\t<decimal>
 *     datum\t<label>\t<table>\t<keys>\t<observations>\t<errors>\t<reported quality bits>
 *     reads\t<label>\t<numReadsProcessed>
 *     error\t<what>\t<exception>\t<message>
 *
 * Usage: BaseRecalibrationEngineDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.recalibration.BaseRecalibrationEngine;
import org.broadinstitute.hellbender.utils.recalibration.RecalDatum;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationTables;

import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class BaseRecalibrationEngineDump {

    /** The reference every read here is placed against, written out so the port reads the same. */
    static final String REFERENCE = "ACGTACGTACGTTTTTGGGGCCCCAAAAACGTACGTACGTGATTACAGGC";

    public static void main(final String[] args) throws Exception {
        System.out.println("# BaseRecalibrationEngineDump: BQSR's counting pass");

        final RecalibrationArgumentCollection rac = new RecalibrationArgumentCollection();
        System.out.printf("const\tenableBAQ\t%b%n", rac.enableBAQ);
        System.out.printf("const\tPRESERVE_QSCORES_LESS_THAN\t%d%n", rac.PRESERVE_QSCORES_LESS_THAN);
        System.out.printf("const\tdefaultBaseQualities\t%d%n", rac.defaultBaseQualities);
        System.out.printf("const\tuseOriginalBaseQualities\t%b%n", rac.useOriginalBaseQualities);
        System.out.printf("const\treference\t%s%n", REFERENCE);

        final Path dir = Path.of("baserecal-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);
        final Path fasta = writeReference(dir);

        events(fasta);
        fractional();
        rounding();
        engine(fasta);
    }

    /** calculateIsSNPOrIndel, over cigars that put the mark in different places. */
    static void events(final Path fasta) throws Exception {
        final Method method = BaseRecalibrationEngine.class.getDeclaredMethod(
                "calculateIsSNPOrIndel", GATKRead.class, ReferenceDataSource.class, int[].class,
                int[].class, int[].class);
        method.setAccessible(true);

        try (final ReferenceDataSource reference = new ReferenceFileSource(fasta)) {
            final String[][] cases = {
                    // Every base matches, so there is nothing to mark.
                    {"exact", "10M", "1", "ACGTACGTAC", "0"},
                    // One mismatch, forward.
                    {"one-mismatch", "10M", "1", "ACGTTCGTAC", "0"},
                    // The same read on the reverse strand, to show the SNP mark does not move.
                    {"one-mismatch-reverse", "10M", "1", "ACGTTCGTAC", "16"},
                    // A deletion in the middle, forward and reverse: the mark is on the base before
                    // on one strand and the base after on the other.
                    {"deletion", "4M2D6M", "1", "ACGTACGTAC", "0"},
                    {"deletion-reverse", "4M2D6M", "1", "ACGTACGTAC", "16"},
                    // An insertion, both strands.
                    {"insertion", "4M2I4M", "1", "ACGTACGTAC", "0"},
                    {"insertion-reverse", "4M2I4M", "1", "ACGTACGTAC", "16"},
                    // At the very ends, where the mark is clamped away rather than wrapping.
                    {"leading-deletion", "1D9M", "1", "ACGTACGTAC", "0"},
                    {"trailing-deletion", "9M1D", "1", "ACGTACGTAC", "0"},
                    {"leading-insertion", "1I9M", "1", "ACGTACGTAC", "0"},
                    {"trailing-insertion", "9M1I", "1", "ACGTACGTAC", "0"},
                    // A soft clip, which advances the read but not the reference.
                    {"soft-clipped", "3S7M", "1", "ACGTACGTAC", "0"},
                    // A skipped region, which advances the reference but not the read.
                    {"skipped", "4M2N6M", "1", "ACGTACGTAC", "0"},
                    // An N in the read, which is not equal to any reference base.
                    {"n-in-read", "10M", "1", "ACGTNCGTAC", "0"},
                    // Placed where the reference is a homopolymer, so several bases disagree.
                    {"many-mismatches", "10M", "13", "ACGTACGTAC", "0"},
            };
            for (final String[] one : cases) {
                final SAMRecord record = read(one[0], one[1], Integer.parseInt(one[2]), one[3],
                        Integer.parseInt(one[4]));
                final GATKRead gatkRead = new SAMRecordToGATKReadAdapter(record);
                final int[] snp = new int[gatkRead.getLength()];
                final int[] ins = new int[snp.length];
                final int[] del = new int[snp.length];
                try {
                    final int errors = (int) method.invoke(null, gatkRead, reference, snp, ins, del);
                    System.out.printf("events\t%s\t%d\t%s\t%s\t%s%n", one[0], errors, ints(snp),
                            ints(ins), ints(del));
                } catch (final Exception e) {
                    final Throwable cause = e.getCause() == null ? e : e.getCause();
                    System.out.printf("error\tevents@%s\t%s\t%s%n", one[0],
                            cause.getClass().getSimpleName(), cause.getMessage());
                }
            }
        }
    }

    /** calculateFractionalErrorArray, whose blocks start one base early. */
    static void fractional() {
        final int[][] errors = {
                {0, 0, 0, 0, 0, 0},
                {0, 1, 0, 0, 0, 0},
                {1, 0, 0, 0, 0, 1},
                {1, 1, 1, 1, 1, 1},
        };
        // 64 is NO_BAQ_UNCERTAINTY, so a block is a run of anything else.
        final byte[][] baqs = {
                {64, 64, 64, 64, 64, 64},
                {64, 60, 60, 64, 64, 64},
                {60, 60, 60, 60, 60, 60},
                {64, 64, 64, 64, 64, 60},
                {60, 64, 64, 64, 64, 64},
        };
        for (int e = 0; e < errors.length; e++) {
            for (int b = 0; b < baqs.length; b++) {
                final double[] out =
                        BaseRecalibrationEngine.calculateFractionalErrorArray(errors[e], baqs[b]);
                final StringBuilder text = new StringBuilder();
                for (final double value : out) {
                    if (text.length() != 0) {
                        text.append(',');
                    }
                    text.append(Long.toHexString(Double.doubleToRawLongBits(value)));
                }
                System.out.printf("fractional\te%d-b%d\t%s%n", e, b, text);
            }
        }
        // Mismatched lengths, which the function refuses.
        try {
            BaseRecalibrationEngine.calculateFractionalErrorArray(new int[3], new byte[4]);
            System.out.println("error\tfractional-length-mismatch\tnone\t-");
        } catch (final Exception e) {
            System.out.printf("error\tfractional-length-mismatch\t%s\t%s%n",
                    e.getClass().getSimpleName(), e.getMessage());
        }
    }

    /** roundToNDecimalPlaces, whose ulp is added BEFORE the rounding. */
    static void rounding() {
        final double[] values = {
                0.0, 1.0, 0.005, 0.015, 0.025, 0.125, 1.005, 2.675, 30.0, 30.00005,
                1.0 / 3.0, 2.0 / 3.0, 1e-9, 1234.56789,
        };
        for (final double value : values) {
            for (final int places : new int[] {2, 4}) {
                final double rounded = MathUtils.roundToNDecimalPlaces(value, places);
                System.out.printf("round\t%s\t%d\t%s\t%s%n", value, places,
                        Long.toHexString(Double.doubleToRawLongBits(rounded)), rounded);
            }
        }
        try {
            MathUtils.roundToNDecimalPlaces(1.0, 0);
            System.out.println("error\tround-zero-places\tnone\t-");
        } catch (final Exception e) {
            System.out.printf("error\tround-zero-places\t%s\t%s%n", e.getClass().getSimpleName(),
                    e.getMessage());
        }
    }

    /** The whole engine over a small corpus, with and without known sites and BAQ. */
    static void engine(final Path fasta) throws Exception {
        final SAMFileHeader header = header();
        final List<SAMRecord> corpus = new ArrayList<>();
        corpus.add(read("r0", "10M", 1, "ACGTACGTAC", 0));
        corpus.add(read("r1", "10M", 1, "ACGTTCGTAC", 0));
        corpus.add(read("r2", "4M2D6M", 5, "ACGTACGTAC", 0));
        corpus.add(read("r3", "4M2I4M", 9, "ACGTACGTAC", 16));
        corpus.add(read("r4", "10M", 13, "ACGTACGTAC", 0));
        for (final SAMRecord record : corpus) {
            record.setHeader(header);
        }

        run(fasta, header, corpus, "plain", List.of(), false);
        // One known site over the middle of the first read, which the skip array removes.
        run(fasta, header, corpus, "known-site",
                List.of(new SimpleInterval("chr1", 3, 5)), false);
        // A known site covering everything.
        run(fasta, header, corpus, "known-everything",
                List.of(new SimpleInterval("chr1", 1, 50)), false);
        // With BAQ on, which is not the default.
        run(fasta, header, corpus, "baq-enabled", List.of(), true);
    }

    static void run(final Path fasta, final SAMFileHeader header, final List<SAMRecord> corpus,
                    final String label, final List<? extends Locatable> knownSites,
                    final boolean enableBaq) throws Exception {
        final RecalibrationArgumentCollection rac = new RecalibrationArgumentCollection();
        rac.enableBAQ = enableBaq;
        final BaseRecalibrationEngine engine =
                new BaseRecalibrationEngine(rac, header);
        engine.logCovariatesUsed();

        try (final ReferenceDataSource reference = new ReferenceFileSource(fasta)) {
            for (final SAMRecord record : corpus) {
                final SAMRecord copy = record.deepCopy();
                copy.setHeader(header);
                engine.processRead(new SAMRecordToGATKReadAdapter(copy), reference, knownSites);
            }
        }
        System.out.printf("reads\t%s\t%d%n", label, engine.getNumReadsProcessed());
        engine.finalizeData();

        final RecalibrationTables tables = engine.getFinalRecalibrationTables();
        for (int index = 0; index < tables.numTables(); index++) {
            for (final NestedIntegerArray.Leaf<RecalDatum> leaf : tables.getTable(index).getAllLeaves()) {
                System.out.printf("datum\t%s\t%d\t%s\t%d\t%.2f\t%s%n", label, index,
                        ints(leaf.keys), leaf.value.getNumObservations(),
                        leaf.value.getNumMismatches(),
                        Long.toHexString(Double.doubleToRawLongBits(leaf.value.getReportedQuality())));
            }
        }
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(
                List.of(new SAMSequenceRecord("chr1", REFERENCE.length()))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("s1");
        group.setPlatform("ILLUMINA");
        group.setPlatformUnit("unit-rg1");
        header.addReadGroup(group);
        return header;
    }

    static SAMRecord read(final String name, final String cigar, final int start,
                          final String bases, final int flags) {
        final SAMRecord record = new SAMRecord(header());
        record.setReadName(name);
        record.setFlags(flags);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setReadBases(bases.getBytes(StandardCharsets.UTF_8));
        final byte[] quals = new byte[bases.length()];
        for (int i = 0; i < quals.length; i++) {
            // A gradient, so the skip array's quality test has something on both sides of six.
            quals[i] = (byte) (2 + i * 4);
        }
        record.setBaseQualities(quals);
        record.setMappingQuality(60);
        record.setAttribute("RG", "rg1");
        return record;
    }

    static Path writeReference(final Path dir) throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        Files.writeString(fasta, ">chr1\n" + REFERENCE + "\n", StandardCharsets.UTF_8);
        FastaSequenceIndexCreator.create(fasta, true);
        // The sequence dictionary, which ReferenceFileSource wants beside the FASTA.
        final Path dict = dir.resolve("reference.dict");
        Files.writeString(dict, "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:" + REFERENCE.length()
                + "\tM5:0\tUR:file:" + fasta + "\n", StandardCharsets.UTF_8);
        return fasta;
    }

    static void emptyDirectory(final Path dir) throws Exception {
        if (!Files.isDirectory(dir)) {
            return;
        }
        try (final var entries = Files.list(dir)) {
            for (final Path entry : entries.toList()) {
                Files.delete(entry);
            }
        }
    }

    static String ints(final int[] values) {
        final StringBuilder out = new StringBuilder();
        for (final int value : values) {
            if (out.length() != 0) {
                out.append(',');
            }
            out.append(value);
        }
        return out.toString();
    }
}

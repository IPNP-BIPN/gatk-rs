/*
 * CountBases, CountReads and FlagStat, taken from the reference.
 *
 * The first three tools whose output is a NUMBER rather than a BAM. All three are ReadWalkers whose
 * `apply` is one line and whose whole answer is what `onTraversalSuccess` prints, so what they
 * measure is the traversal and the formatting rather than any transform.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THEIR DEFAULT FILTER IS THE ENGINE'S, WellformedReadFilter, so a malformed read is not counted
 *     by any of the three. That is the same pattern PrintReads takes and not the one BaseRecalibrator
 *     or ReadAnonymizer take;
 *   - CountBases COUNTS read.getLength(), WHICH IS THE NUMBER OF BASES AND NOT THE SPAN, so a read
 *     with a deletion counts fewer bases than it covers and one with an insertion counts more;
 *   - FlagStat's PERCENTAGES ARE COMPUTED IN float AND FORMATTED WITH DecimalFormat("#0.00"), which
 *     rounds HALF_EVEN where String.format("%.2f") rounds HALF_UP. Two different roundings, and this
 *     is the one that reaches the output;
 *   - AND WITH NO READS AT ALL THE RATIO IS 0f/0f, so the percentage is NaN and DecimalFormat writes
 *     whatever its symbols say NaN is;
 *   - read2 IS TESTED BEFORE read1 AND THE TWO ARE `else if`, so a read carrying both 0x40 and 0x80
 *     counts as read2 only;
 *   - `singletons` AND `with_itself_and_mate_mapped` ARE BOTH INSIDE the paired branch, so an
 *     unpaired read contributes to neither however it is mapped;
 *   - AND `isUnmapped` IS THE THREE-PART TEST, not the 0x4 flag, so a read with the flag clear, no
 *     reference and a zero start counts as unmapped in `mapped` and in `singletons` alike.
 *
 * Output:
 *
 *     fixture\t<label>\t<the input BAM, base64>
 *     fixtureindex\t<label>\t<the index, base64>
 *     count\t<tool>\t<label>\t<the printed text, escaped>
 *     error\t<tool>\t<label>\t<exception>\t<message>
 *
 * Usage: CountingWalkersDump
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
import org.broadinstitute.hellbender.tools.CountBases;
import org.broadinstitute.hellbender.tools.CountReads;
import org.broadinstitute.hellbender.tools.FlagStat;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CountingWalkersDump {

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final Path dir = Path.of("countingwalkers-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# CountingWalkersDump: CountBases, CountReads and FlagStat");

        // Every flag combination the three tools look at, plus reads whose length is not their span.
        final Path full = dir.resolve("full.bam");
        buildFull(full.toFile());
        emitFixture(dir, full, "full");

        // No reads at all, which is where the percentages divide by zero.
        final Path empty = dir.resolve("empty.bam");
        buildEmpty(empty.toFile());
        emitFixture(dir, empty, "empty");

        // Only reads the default filter drops, so the traversal keeps nothing.
        final Path malformed = dir.resolve("malformed.bam");
        buildMalformed(malformed.toFile());
        emitFixture(dir, malformed, "malformed");

        for (final String label : new String[] {"full", "empty", "malformed"}) {
            final Path bam = dir.resolve(label + ".bam");
            run(dir, bam, "CountBases", label, new String[] {});
            run(dir, bam, "CountReads", label, new String[] {});
            run(dir, bam, "FlagStat", label, new String[] {});
            // One interval that covers only part of the corpus, so the traversal keeps a subset
            // and the percentages are computed over a different denominator.
            run(dir, bam, "CountReads", label + "-interval", new String[] {"-L", "chr1:1-12"});
            run(dir, bam, "FlagStat", label + "-interval", new String[] {"-L", "chr1:1-12"});
            run(dir, bam, "CountBases", label + "-interval", new String[] {"-L", "chr1:1-12"});
            // With the filter disabled, so the malformed reads are counted.
            run(dir, bam, "CountReads", label + "-unfiltered",
                    new String[] {"--disable-read-filter", "WellformedReadFilter"});
            run(dir, bam, "FlagStat", label + "-unfiltered",
                    new String[] {"--disable-read-filter", "WellformedReadFilter"});
        }
    }

    static void emitFixture(final Path dir, final Path bam, final String label) throws Exception {
        System.out.printf("fixture\t%s\t%s%n", label, RecordTransformDump.base64(bam));
        final Path index = dir.resolve(label + ".bai");
        System.out.printf("fixtureindex\t%s\t%s%n", label,
                Files.exists(index) ? RecordTransformDump.base64(index) : "absent");
    }

    /** Every branch of FlagStat's counter, and two reads whose length is not their span. */
    static void buildFull(final File file) {
        final SAMFileHeader header = header();
        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, file)) {
            // Plain, unpaired, mapped: contributes to mapped and to nothing paired.
            writer.addAlignment(read(header, "plain", 1, "10M", 0, 0, 100, 60));
            // A deletion: ten bases, twelve reference positions.
            writer.addAlignment(read(header, "deletion", 3, "4M2D6M", 0, 0, 100, 60));
            // An insertion: ten bases, eight reference positions.
            writer.addAlignment(read(header, "insertion", 5, "4M2I4M", 0, 0, 100, 60));
            // Paired, first of pair, properly paired, mate on the same contig.
            writer.addAlignment(read(header, "first", 7, "10M", 0x1 | 0x2 | 0x40, 0, 40, 60));
            // Paired, second of pair.
            writer.addAlignment(read(header, "second", 9, "10M", 0x1 | 0x2 | 0x80, 0, 40, 60));
            // Both 0x40 and 0x80, which the `else if` counts as read2 only.
            writer.addAlignment(read(header, "both-of-pair", 11, "10M", 0x1 | 0x40 | 0x80, 0, 40, 60));
            // A singleton: paired, mapped, mate unmapped.
            writer.addAlignment(read(header, "singleton", 13, "10M", 0x1 | 0x8, 0, 40, 60));
            // Mate on a different contig, mapping quality above five.
            writer.addAlignment(read(header, "other-chr", 15, "10M", 0x1, 1, 40, 60));
            // The same, below five, so only one of the two counters moves.
            writer.addAlignment(read(header, "other-chr-lowmq", 17, "10M", 0x1, 1, 40, 3));
            // Vendor failure and duplicate, which are counted and still traversed.
            writer.addAlignment(read(header, "vendor-fail", 19, "10M", 0x200, 0, 100, 60));
            writer.addAlignment(read(header, "duplicate", 21, "10M", 0x400, 0, 100, 60));
            // Unmapped by the flag, which is not counted in `mapped`.
            writer.addAlignment(read(header, "unmapped", 23, "10M", 0x4, 0, 100, 60));
        }
    }

    static void buildEmpty(final File file) {
        final SAMFileHeader header = header();
        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, file)) {
            // Nothing at all.
            final int keepTheWriterHonest = 0;
            if (keepTheWriterHonest != 0) {
                writer.addAlignment(read(header, "never", 1, "10M", 0, 0, 100, 60));
            }
        }
    }

    /** Only reads the engine's default filter drops: a cigar that does not match the length. */
    static void buildMalformed(final File file) {
        final SAMFileHeader header = header();
        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, file)) {
            writer.addAlignment(read(header, "short-cigar", 1, "5M", 0, 0, 100, 60));
            writer.addAlignment(read(header, "long-cigar", 3, "15M", 0, 0, 100, 60));
        }
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 200),
                new SAMSequenceRecord("chr2", 200))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("s1");
        header.addReadGroup(group);
        return header;
    }

    static SAMRecord read(final SAMFileHeader header, final String name, final int start,
                          final String cigar, final int flags, final int mateRef,
                          final int mateStart, final int mapq) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setFlags(flags);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setReadBases("ACGTACGTAC".getBytes(StandardCharsets.UTF_8));
        record.setBaseQualities(new byte[] {30, 30, 30, 30, 30, 30, 30, 30, 30, 30});
        record.setMappingQuality(mapq);
        record.setMateReferenceIndex(mateRef);
        record.setMateAlignmentStart(mateStart);
        record.setAttribute("RG", "rg1");
        return record;
    }

    /** One run of one tool, with the text it printed to its optional output file. */
    static void run(final Path dir, final Path input, final String tool, final String label,
                    final String[] extra) throws Exception {
        final Path output = dir.resolve(tool + "." + label.replace(':', '_') + ".txt");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-I", input.toString(), "-O", output.toString(),
                "--use-jdk-inflater", "true"));
        argv.addAll(Arrays.asList(extra));

        try {
            switch (tool) {
                case "CountBases" -> new CountBases().instanceMain(argv.toArray(new String[0]));
                case "CountReads" -> new CountReads().instanceMain(argv.toArray(new String[0]));
                case "FlagStat" -> new FlagStat().instanceMain(argv.toArray(new String[0]));
                default -> throw new IllegalStateException("no tool " + tool);
            }
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s\t%s\t%s%n", tool, label, e.getClass().getSimpleName(),
                    e.getMessage());
            return;
        }
        System.out.printf("count\t%s\t%s\t%s%n", tool, label,
                ReferenceQueryDump.escape(Files.readString(output)));
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
}

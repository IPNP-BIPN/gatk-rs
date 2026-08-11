/*
 * PostProcessReadsForRSEM, taken from the reference.
 *
 * The twelfth whole tool of the record-transform archetype, the third that is not a walker, and the
 * first that groups the traversal into query-name runs of its own accord.
 *
 * Seven behaviours this is built to catch, and two of them are crashes.
 *
 *   - A FOURTH getDefaultReadFilters PATTERN. PrintReads takes GATKTool's default,
 *     UnmarkDuplicates replaces it with ALLOW_ALL_READS, PrintDistantMates extends it with four
 *     more; this one replaces the whole list with a SINGLE filter that is not Wellformed:
 *     `Collections.singletonList(ReadFilterLibrary.NOT_SUPPLEMENTARY_ALIGNMENT)`. So a
 *     supplementary alignment never reaches the tool, and the `supplementaryAlignments` list
 *     `ReadPair` maintains for it is dead code on this path;
 *   - A PAIR WITH NO FIRST-OF-PAIR THROWS NullPointerException. `passesRSEMFilter` opens with
 *     `if (read1 == null || read2 == null) { logger.warn("..." + read1.getName()); return false; }`
 *     and dereferences `read1` inside the branch that exists because it may be null. A file whose
 *     query-name group holds only a second-of-pair therefore kills the run rather than skipping the
 *     pair. Measured in its own fixture, because it takes the whole traversal down with it;
 *   - SO DOES A ONE-SIDED SET OF SECONDARY ALIGNMENTS. `groupSecondaryReads` collects into
 *     `Collectors.groupingBy(GATKRead::isFirstOfPair)`, which produces no `false` key at all when
 *     every secondary is a first-of-pair, and then calls `read2Reads.size()`. The guard that reads
 *     `read1Reads.size() != read2Reads.size()` is written as if both were present;
 *   - THE OUTPUT ORDER IS PRIMARY PAIR, THEN EACH SECONDARY PAIR, each as first-of-pair then
 *     second-of-pair, which is not the order the query-name group arrived in. A secondary is
 *     matched to its mate by contig, `getStart() == mate.getMateStart()` and
 *     `getMateStart() == mate.getStart()`, not by anything in the record that says so;
 *   - THREE REASONS A PAIR IS DROPPED, and all three drop BOTH reads: either read unmapped, the two
 *     on different contigs, or a cigar that is not exactly one `M`. `100=` is a single element and
 *     is not `M`, so it is refused: the test is on the operator, not on the length or the count
 *     alone;
 *   - A PRIMARY PAIR THAT FAILS TAKES ITS SECONDARY ALIGNMENTS WITH IT, because the secondary loop
 *     is inside the primary's `if`. A secondary pair that fails on its own drops only itself;
 *   - AND THE INPUT MUST BE QUERY-NAME SORTED, checked in `onTraversalStart` against the header
 *     rather than against the reads, so a mis-labelled file is trusted and a correctly labelled one
 *     is not re-verified.
 *
 * Neither input nor output is opened with an index: a queryname-sorted BAM cannot have one.
 *
 * Output:
 *
 *     deflater\t<class>
 *     filters\t<class>
 *     fixture\t<label>\t<base64 bam>
 *     header\t<label>\t<escaped SAM header>
 *     commandline\t<label>\t<@PG command line>
 *     output\t<label>\t<base64 bam>
 *     index\t<label>\t<base64 bai or absent>
 *     reads\t<label>\t<name>\t<flags>\t<contig>\t<start>\t<cigar>
 *     error\t<label>\t<class>:<message>
 *
 * Usage: PostProcessReadsForRSEMDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.walkers.qc.PostProcessReadsForRSEM;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;
import java.util.function.Consumer;

public class PostProcessReadsForRSEMDump {

    static final int PAIRED = 0x1;
    static final int PROPER = 0x2;
    static final int UNMAPPED = 0x4;
    static final int FIRST = 0x40;
    static final int SECOND = 0x80;
    static final int SECONDARY = 0x100;
    static final int SUPPLEMENTARY = 0x800;

    public static void main(final String[] args) throws Exception {
        // The factory is static and the first writer wins. This dump calls no Picard entry point,
        // so nothing should replace it; the pin makes that a fact rather than a hope.
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // Relative on purpose: the string handed to -I and -O is the string recorded inside the
        // output BAM's own @PG, so an absolute temporary path would make every output byte
        // unstable and canonicalization cannot reach inside base64.
        final Path dir = Path.of("rsem-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# PostProcessReadsForRSEMDump: PostProcessReadsForRSEM");
        System.out.printf("deflater\t%s%n",
                BlockCompressedOutputStream.getDefaultDeflaterFactory().getClass().getName());
        // The whole default list, which is one filter and is not Wellformed.
        for (final ReadFilter filter : new PostProcessReadsForRSEM().getDefaultReadFilters()) {
            System.out.printf("filters\t%s%n", filter.getClass().getSimpleName());
        }

        final Path plain = dir.resolve("plain.bam");
        buildFixture(plain.toFile(), SAMFileHeader.SortOrder.queryname, writer -> {
            // p1: a clean pair with two secondary alignments that find each other.
            pair(writer, "p1", 100, 300);
            secondaryPair(writer, "p1", 500, 700);
            // p2: the mate is unmapped.
            final SAMRecord unmappedMate = read(writer, "p2", 0, PAIRED | SECOND | UNMAPPED,
                    "chr1", 0, "100M", "chr1", 900);
            writer.addAlignment(read(writer, "p2", 900, PAIRED | FIRST, "chr1", 900, "100M",
                    "chr1", 0));
            writer.addAlignment(unmappedMate);
            // p3: the two halves are on different contigs.
            writer.addAlignment(read(writer, "p3", 1100, PAIRED | FIRST, "chr1", 1100, "100M",
                    "chr2", 100));
            writer.addAlignment(read(writer, "p3", 100, PAIRED | SECOND, "chr2", 100, "100M",
                    "chr1", 1100));
            // p4: two cigar elements.
            writer.addAlignment(read(writer, "p4", 1300, PAIRED | FIRST, "chr1", 1300, "50M50S",
                    "chr1", 1500));
            writer.addAlignment(read(writer, "p4", 1500, PAIRED | SECOND, "chr1", 1500, "100M",
                    "chr1", 1300));
            // p5: one cigar element, and it is not M. `=` is a match and still refused.
            writer.addAlignment(read(writer, "p5", 1700, PAIRED | FIRST, "chr1", 1700, "100=",
                    "chr1", 1900));
            writer.addAlignment(read(writer, "p5", 1900, PAIRED | SECOND, "chr1", 1900, "100M",
                    "chr1", 1700));
            // p6: a clean pair with a supplementary alignment the filter removes before the tool
            // ever sees it.
            pair(writer, "p6", 2100, 2300);
            writer.addAlignment(read(writer, "p6", 2500, PAIRED | FIRST | SUPPLEMENTARY, "chr1",
                    2500, "100M", "chr1", 2300));
            // p7: a first-of-pair with no second. read1 is not null, so the warn does not throw and
            // the pair is simply dropped.
            writer.addAlignment(read(writer, "p7", 2700, PAIRED | FIRST, "chr1", 2700, "100M",
                    "chr1", 2900));
            // p8: a clean pair whose secondary pair is itself chimeric, so only the secondary is
            // dropped and the primary survives.
            pair(writer, "p8", 3100, 3300);
            writer.addAlignment(read(writer, "p8", 3500, PAIRED | FIRST | SECONDARY, "chr1", 3500,
                    "100M", "chr2", 500));
            writer.addAlignment(read(writer, "p8", 500, PAIRED | SECOND | SECONDARY, "chr2", 500,
                    "100M", "chr1", 3500));
        });
        fixture(plain, "plain");

        // Only a second-of-pair: `read1.getName()` inside the null guard.
        final Path noFirst = dir.resolve("no_first.bam");
        buildFixture(noFirst.toFile(), SAMFileHeader.SortOrder.queryname, writer ->
                writer.addAlignment(read(writer, "q1", 100, PAIRED | SECOND, "chr1", 100, "100M",
                        "chr1", 300)));
        fixture(noFirst, "no_first");

        // Every secondary alignment is a first-of-pair: `groupingBy` produces no `false` key.
        final Path oneSided = dir.resolve("one_sided.bam");
        buildFixture(oneSided.toFile(), SAMFileHeader.SortOrder.queryname, writer -> {
            pair(writer, "r1", 100, 300);
            writer.addAlignment(read(writer, "r1", 500, PAIRED | FIRST | SECONDARY, "chr1", 500,
                    "100M", "chr1", 700));
        });
        fixture(oneSided, "one_sided");

        // The mirror image: every secondary is a second-of-pair, so `groupedByRead1.get(true)` is
        // the null one. Java evaluates `read1Reads.size()` first, so the message names the other
        // list, and a port that assumed one message for both shapes would be inventing one.
        final Path otherSided = dir.resolve("other_sided.bam");
        buildFixture(otherSided.toFile(), SAMFileHeader.SortOrder.queryname, writer -> {
            pair(writer, "t1", 100, 300);
            writer.addAlignment(read(writer, "t1", 700, PAIRED | SECOND | SECONDARY, "chr1", 700,
                    "100M", "chr1", 500));
        });
        fixture(otherSided, "other_sided");

        // The one sort order the tool checks, and it checks the header rather than the reads.
        final Path coordinate = dir.resolve("coordinate.bam");
        buildFixture(coordinate.toFile(), SAMFileHeader.SortOrder.coordinate, writer ->
                pair(writer, "s1", 100, 300));
        fixture(coordinate, "coordinate");

        run(dir, plain, "plain");
        run(dir, noFirst, "nofirst");
        run(dir, oneSided, "onesided");
        run(dir, otherSided, "othersided");
        run(dir, coordinate, "coordinate");

        reads(dir, "plain");
        readsOf(plain, "in:plain");
    }

    /** A clean primary pair on chr1. */
    static void pair(final SAMFileWriter writer, final String name, final int first,
                     final int second) {
        writer.addAlignment(read(writer, name, first, PAIRED | PROPER | FIRST, "chr1", first,
                "100M", "chr1", second));
        writer.addAlignment(read(writer, name, second, PAIRED | PROPER | SECOND, "chr1", second,
                "100M", "chr1", first));
    }

    /** A secondary pair that will find each other through their mate starts. */
    static void secondaryPair(final SAMFileWriter writer, final String name, final int first,
                              final int second) {
        writer.addAlignment(read(writer, name, first, PAIRED | FIRST | SECONDARY, "chr1", first,
                "100M", "chr1", second));
        writer.addAlignment(read(writer, name, second, PAIRED | SECOND | SECONDARY, "chr1", second,
                "100M", "chr1", first));
    }

    static SAMRecord read(final SAMFileWriter writer, final String name, final int ignored,
                          final int flags, final String contig, final int start,
                          final String cigar, final String mateContig, final int mateStart) {
        final SAMFileHeader header = writer.getFileHeader();
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setFlags(flags);
        if ((flags & UNMAPPED) == 0) {
            record.setReferenceName(contig);
            record.setAlignmentStart(start);
            record.setCigarString(cigar);
            record.setMappingQuality(60);
        } else {
            // An unmapped read placed on its mate, which is how STAR leaves them.
            record.setReferenceName(mateContig);
            record.setAlignmentStart(mateStart);
            record.setMappingQuality(0);
        }
        record.setMateReferenceName(mateContig);
        record.setMateAlignmentStart(mateStart);
        record.setReadBases("ACGTACGTAC".repeat(10).getBytes());
        final byte[] quals = new byte[100];
        Arrays.fill(quals, (byte) 30);
        record.setBaseQualities(quals);
        record.setAttribute("RG", "rg1");
        return record;
    }

    static void buildFixture(final File file, final SAMFileHeader.SortOrder order,
                             final Consumer<SAMFileWriter> body) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 5000), new SAMSequenceRecord("chr2", 5000))));
        header.setSortOrder(order);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("s1");
        header.addReadGroup(group);
        final SAMProgramRecord existing = new SAMProgramRecord("upstream");
        existing.setProgramVersion("1.0");
        header.addProgramRecord(existing);

        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().makeBAMWriter(header, true, file)) {
            body.accept(writer);
        }
    }

    static void fixture(final Path bam, final String label) throws Exception {
        System.out.printf("fixture\t%s\t%s%n", label, base64(bam));
    }

    static void reads(final Path dir, final String label) throws Exception {
        readsOf(dir.resolve("PostProcessReadsForRSEM." + label + ".bam"), label);
    }

    static void readsOf(final Path bam, final String label) throws Exception {
        try (final SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT).open(bam.toFile())) {
            for (final SAMRecord record : reader) {
                System.out.printf("reads\t%s\t%s\t%d\t%s\t%d\t%s%n", label, record.getReadName(),
                        record.getFlags(), record.getReferenceName(), record.getAlignmentStart(),
                        record.getCigarString());
            }
        }
    }

    static void run(final Path dir, final Path input, final String label) throws Exception {
        final Path output = dir.resolve("PostProcessReadsForRSEM." + label + ".bam");
        // --use-jdk-deflater is the knob that decides which bytes come out, for the same reason
        // PrintReadsDump names it: the GKL deflater's output is not yet reproduced.
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-I", input.toString(), "-O", output.toString(),
                "--use-jdk-deflater", "true", "--use-jdk-inflater", "true"));

        try {
            new PostProcessReadsForRSEM().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            // A crash is the observable behaviour, so it is dumped rather than swallowed.
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    String.valueOf(e.getMessage()).replace('\n', ' '));
            return;
        }

        String commandLine = "";
        try (final SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .open(output.toFile())) {
            final SAMFileHeader header = reader.getFileHeader();
            for (final SAMProgramRecord record : header.getProgramRecords()) {
                if (record.getCommandLine() != null) {
                    commandLine = record.getCommandLine();
                }
            }
            System.out.printf("header\t%s\t%s%n", label,
                    ReferenceQueryDump.escape(header.getSAMString()));
        }
        System.out.printf("commandline\t%s\t%s%n", label, commandLine);
        System.out.printf("output\t%s\t%s%n", label, base64(output));

        final Path index = dir.resolve(output.getFileName().toString().replace(".bam", ".bai"));
        System.out.printf("index\t%s\t%s%n", label,
                Files.exists(index) ? base64(index) : "absent");
    }

    static void emptyDirectory(final Path dir) throws Exception {
        if (!Files.isDirectory(dir)) {
            return;
        }
        try (final var entries = Files.list(dir)) {
            for (final Path entry : entries.toList()) {
                Files.deleteIfExists(entry);
            }
        }
    }

    static String base64(final Path path) throws Exception {
        return Base64.getEncoder().encodeToString(Files.readAllBytes(path));
    }
}

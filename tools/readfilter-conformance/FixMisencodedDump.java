/*
 * FixMisencodedBaseQualityReads, taken from the reference.
 *
 * The fourth whole tool of the record-transform archetype, and the first whose transform can
 * refuse the file it was given. Its `apply` is one line; everything interesting is in what the
 * transform does to a read it cannot fix and in which reads reach it at all.
 *
 * Four behaviours this is built to catch, and none of them is the subtraction:
 *
 *   - THE TOOL DOES NOT OVERRIDE ITS READ FILTERS, so it takes GATKTool's default of
 *     WellformedReadFilter, where UnmarkDuplicates and RevertBaseQualityScores replace the whole
 *     list with ALLOW_ALL_READS. Three tools, one archetype, two different default traversals, and
 *     this is the one that keeps the default;
 *   - A QUALITY BELOW 31 ABORTS THE WHOLE RUN, with a message that says the read "was correctly
 *     encoded". Not skipped, not clamped: the transform throws, so the output is whatever had been
 *     flushed. A port that clamped would produce a healthy-looking file the reference never wrote;
 *   - THE CHECK IS PER BASE AND AFTER THE SUBTRACTION, so the reads before the offending one are
 *     already transformed and the offending read is partly transformed when it throws. What
 *     survives in the output is therefore a function of the writer's buffering rather than of the
 *     transform;
 *   - A READ WITH NO QUALITIES AT ALL PASSES THROUGH. The loop does not run, the empty array is
 *     set back, and nothing refuses: `*` in a SAM is not a quality of zero.
 *
 * The output BAMs travel in the golden in full, base64, indexes included, as the other tools' do,
 * and the deflater is pinned and recorded for the same reason: the factory is static, and whoever
 * touches it first wins for the life of the JVM.
 *
 * Output:
 *
 *     deflater\t<class>
 *     fixture\t<label>\t<base64 bam>
 *     fixtureindex\t<label>\t<base64 bai>
 *     header\t<tool>\t<label>\t<escaped SAM header>
 *     commandline\t<tool>\t<label>\t<@PG command line>
 *     output\t<tool>\t<label>\t<base64 bam>
 *     index\t<tool>\t<label>\t<base64 bai or absent>
 *     error\t<tool>\t<label>\t<class>:<message>
 *
 * Usage: FixMisencodedDump
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

import org.broadinstitute.hellbender.tools.FixMisencodedBaseQualityReads;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;

public class FixMisencodedDump {

    public static void main(final String[] args) throws Exception {
        // Before the fixture is written: the factory is static and first writer wins.
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // Relative on purpose: the string handed to -I and -O is the string recorded inside the
        // output BAM's own @PG, so an absolute temporary path would make every output byte
        // unstable and canonicalization cannot reach inside base64.
        final Path dir = Path.of("fixmisencoded-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# FixMisencodedDump: FixMisencodedBaseQualityReads");
        System.out.printf("deflater\t%s%n",
                BlockCompressedOutputStream.getDefaultDeflaterFactory().getClass().getName());

        // Every quality is at or above 31, so the subtraction succeeds on every read.
        final Path high = dir.resolve("high.bam");
        buildFixture(high.toFile(), Quality.ALL_HIGH);
        fixture(dir, high, "high");
        // One read carries a quality below 31, which is what the transform refuses.
        final Path low = dir.resolve("low.bam");
        buildFixture(low.toFile(), Quality.LOW_ON_ONE);
        fixture(dir, low, "low");
        // One read carries no qualities at all, which the loop never enters.
        final Path none = dir.resolve("no-quals.bam");
        buildFixture(none.toFile(), Quality.MISSING_ON_ONE);
        fixture(dir, none, "no-quals");

        // The success, the interval, the run that asks for no index and the one that adds no @PG.
        fix(dir, high, "all", new String[] {});
        fix(dir, high, "chr1", new String[] {"-L", "chr1"});
        fix(dir, high, "chr1:100-160", new String[] {"-L", "chr1:100-160"});
        fix(dir, high, "noindex", new String[] {"--create-output-bam-index", "false"});
        fix(dir, high, "nopg", new String[] {"--add-output-sam-program-record", "false"});
        // This tool keeps GATKTool's default filter, so naming ALLOW_ALL_READS changes the
        // traversal where it would not for the other two tools of this archetype.
        fix(dir, high, "allowall", new String[] {"--read-filter", "AllowAllReadsReadFilter"});
        // The refusal, and the read the loop never enters.
        fix(dir, low, "low-quality", new String[] {});
        fix(dir, none, "no-quals", new String[] {});
        // An interval that excludes the offending read: the refusal is a property of what the
        // traversal reaches, not of the file.
        fix(dir, low, "low-excluded", new String[] {"-L", "chr1:100-130"});
    }

    /**
     * A fixture and the index written beside it.
     *
     * The index travels too because the port's reader needs one to open the file at all, and a
     * test that built its own would be inventing part of the input rather than reading it.
     */
    static void fixture(final Path dir, final Path bam, final String label) throws Exception {
        System.out.printf("fixture\t%s\t%s%n", label, base64(bam));
        final Path index = dir.resolve(bam.getFileName().toString().replace(".bam", ".bai"));
        System.out.printf("fixtureindex\t%s\t%s%n", label, base64(index));
    }

    /** What the fixture's qualities look like. */
    enum Quality {
        ALL_HIGH,
        LOW_ON_ONE,
        MISSING_ON_ONE
    }

    /**
     * A small coordinate-sorted BAM, with a duplicate-flagged read and an unmapped tail.
     *
     * Built here rather than shared, because what varies per case is the qualities and that is the
     * point of every case.
     */
    static void buildFixture(final File file, final Quality quality) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 1000),
                new SAMSequenceRecord("chr2", 1000))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("s1");
        header.addReadGroup(group);
        final SAMProgramRecord existing = new SAMProgramRecord("upstream");
        existing.setProgramVersion("1.0");
        header.addProgramRecord(existing);

        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, file)) {
            for (int i = 0; i < 6; i++) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName("r" + i);
                record.setReferenceName(i < 4 ? "chr1" : "chr2");
                record.setAlignmentStart(100 + i * 20);
                record.setCigarString("10M");
                final byte[] bases = new byte[10];
                Arrays.fill(bases, (byte) 'A');
                record.setReadBases(bases);
                // At or above 31 everywhere, so the subtraction succeeds, except where a case
                // asks otherwise. 31 itself is the boundary: it becomes 0, which is not negative.
                final byte[] quals = new byte[10];
                Arrays.fill(quals, (byte) (31 + i));
                if (quality == Quality.LOW_ON_ONE && i == 2) {
                    // One base below the threshold in the middle of the read, so the transform
                    // has already rewritten the ones before it when it throws.
                    quals[5] = 30;
                }
                if (quality == Quality.MISSING_ON_ONE && i == 2) {
                    record.setBaseQualities(SAMRecord.NULL_QUALS);
                } else {
                    record.setBaseQualities(quals);
                }
                record.setMappingQuality(60);
                record.setAttribute("RG", "rg1");
                // Two of the six are flagged, so unmarking is visible and so is leaving the rest
                // alone.
                record.setDuplicateReadFlag(i == 1 || i == 3);
                writer.addAlignment(record);
            }
        }
    }

    static void fix(final Path dir, final Path input, final String label, final String[] extra)
            throws Exception {
        run("FixMisencodedBaseQualityReads", dir, input, label, extra, argv -> {
            new FixMisencodedBaseQualityReads().instanceMain(argv);
            return null;
        });
    }

    interface Invocation {
        Void run(String[] argv) throws Exception;
    }

    static void run(final String tool, final Path dir, final Path input, final String label,
                    final String[] extra, final Invocation invocation) throws Exception {
        final Path output = dir.resolve(tool + "." + label.replace(':', '_') + ".bam");
        // --use-jdk-deflater is the knob that decides which bytes come out, for the same reason
        // PrintReadsDump names it: the GKL deflater's output is not yet reproduced.
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-I", input.toString(), "-O", output.toString(),
                "--use-jdk-deflater", "true", "--use-jdk-inflater", "true"));
        argv.addAll(Arrays.asList(extra));

        try {
            invocation.run(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            // The refusal is the observable behaviour, so it is dumped rather than swallowed.
            System.out.printf("error\t%s\t%s\t%s:%s%n", tool, label, e.getClass().getName(),
                    e.getMessage());
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
            System.out.printf("header\t%s\t%s\t%s%n", tool, label,
                    ReferenceQueryDump.escape(header.getSAMString()));
        }
        System.out.printf("commandline\t%s\t%s\t%s%n", tool, label, commandLine);
        System.out.printf("output\t%s\t%s\t%s%n", tool, label, base64(output));

        final Path index = dir.resolve(output.getFileName().toString().replace(".bam", ".bai"));
        System.out.printf("index\t%s\t%s\t%s%n", tool, label,
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

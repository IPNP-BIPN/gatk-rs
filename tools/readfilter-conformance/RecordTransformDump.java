/*
 * UnmarkDuplicates and RevertBaseQualityScores, taken from the reference.
 *
 * The second and third whole tools, and the point is the calibration gate rather than the tools:
 * G2 asks what a member of the largest archetype costs once the first one has paid for the engine.
 * One dump covers both, which is itself part of the answer.
 *
 * Three behaviours this is built to catch, and none of them is the transform:
 *
 *   - BOTH TOOLS REPLACE THE DEFAULT READ FILTERS with ALLOW_ALL_READS, where PrintReads takes
 *     GATKTool's default of WellformedReadFilter. So on a file with a malformed read the three
 *     tools emit different sets of reads, and it is the filter rather than the apply that decides
 *     it. The fixture carries such a read for exactly this reason;
 *   - REVERTBASEQUALITYSCORES ABORTS THE WHOLE RUN on a read with no OQ. Not skips, not passes
 *     through: it throws a UserException, so the output is whatever had been flushed. A port that
 *     passed the read through would produce a larger and healthier-looking file than the
 *     reference;
 *   - AN EMPTY OQ IS THE SAME AS AN ABSENT ONE, because getOriginalBaseQualities returns null for
 *     both, even though fastqToPhred("") is happy and returns an empty array. That conflation is
 *     the reference's and is measured rather than inferred.
 *
 * The output BAMs travel in the golden in full, base64, indexes included, as PrintReadsDump's do.
 * The deflater is pinned and recorded for the same reason: the factory is static, and whoever
 * touches it first wins for the life of the JVM.
 *
 * Output:
 *
 *     deflater\t<class>
 *     fixture\t<base64 bam>
 *     header\t<tool>\t<label>\t<escaped SAM header>
 *     commandline\t<tool>\t<label>\t<the @PG CL the tool recorded>
 *     output\t<tool>\t<label>\t<base64 bam>
 *     index\t<tool>\t<label>\t<base64 bai or `absent`>
 *     error\t<tool>\t<label>\t<exception class>:<message>
 *
 * Usage: RecordTransformDump
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

import org.broadinstitute.hellbender.tools.walkers.RevertBaseQualityScores;
import org.broadinstitute.hellbender.tools.walkers.UnmarkDuplicates;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;

public class RecordTransformDump {

    public static void main(final String[] args) throws Exception {
        // Before the fixture is written: the factory is static and first writer wins.
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // Relative on purpose: the string handed to -I and -O is the string recorded inside the
        // output BAM's own @PG, so an absolute temporary path would make every output byte
        // unstable and canonicalization cannot reach inside base64.
        final Path dir = Path.of("recordtransform-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# RecordTransformDump: UnmarkDuplicates and RevertBaseQualityScores");
        System.out.printf("deflater\t%s%n",
                BlockCompressedOutputStream.getDefaultDeflaterFactory().getClass().getName());

        // Every read carries OQ, so RevertBaseQualityScores can succeed.
        final Path full = dir.resolve("full.bam");
        buildFixture(full.toFile(), Oq.ALL);
        System.out.printf("fixture\tfull\t%s%n", base64(full));
        // One read has no OQ at all, and one has an empty OQ. Both abort the revert, and the
        // duplicate tool does not care.
        final Path partial = dir.resolve("partial.bam");
        buildFixture(partial.toFile(), Oq.MISSING_ON_ONE);
        System.out.printf("fixture\tpartial\t%s%n", base64(partial));
        final Path empty = dir.resolve("empty-oq.bam");
        buildFixture(empty.toFile(), Oq.EMPTY_ON_ONE);
        System.out.printf("fixture\tempty-oq\t%s%n", base64(empty));

        // UnmarkDuplicates: the transform, the interval, the filter that is not the default, and
        // the run that asks for no index.
        unmark(dir, full, "all", new String[] {});
        unmark(dir, full, "chr1", new String[] {"-L", "chr1"});
        unmark(dir, full, "chr1:100-160", new String[] {"-L", "chr1:100-160"});
        unmark(dir, full, "wellformed",
                new String[] {"--read-filter", "WellformedReadFilter"});
        unmark(dir, full, "noindex", new String[] {"--create-output-bam-index", "false"});
        unmark(dir, full, "nopg",
                new String[] {"--add-output-sam-program-record", "false"});
        // The tool that does not read OQ does not mind a file missing it.
        unmark(dir, partial, "partial-input", new String[] {});

        // RevertBaseQualityScores: the success, and the two ways to abort.
        revert(dir, full, "all", new String[] {});
        revert(dir, full, "chr1:100-160", new String[] {"-L", "chr1:100-160"});
        revert(dir, full, "noindex", new String[] {"--create-output-bam-index", "false"});
        revert(dir, partial, "missing-oq", new String[] {});
        revert(dir, empty, "empty-oq", new String[] {});
    }

    /** Which reads of the fixture carry an OQ tag. */
    enum Oq {
        ALL,
        MISSING_ON_ONE,
        EMPTY_ON_ONE
    }

    /**
     * A small coordinate-sorted BAM, with a duplicate-flagged read and an unmapped tail.
     *
     * Built here rather than shared with ReadWalkerDump because these tools need OQ, and a fixture
     * that varies per case is the point of three of the five cases.
     */
    static void buildFixture(final File file, final Oq oq) {
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
                final byte[] quals = new byte[10];
                Arrays.fill(quals, (byte) (20 + i));
                record.setBaseQualities(quals);
                record.setMappingQuality(60);
                record.setAttribute("RG", "rg1");
                // Two of the six are flagged, so unmarking is visible and so is leaving the rest
                // alone.
                record.setDuplicateReadFlag(i == 1 || i == 3);
                // The original qualities differ from the current ones, so a revert is observable
                // rather than a no-op.
                final String oqString = "!!!!!!!!!!".substring(0, 10);
                if (oq == Oq.ALL || i != 2) {
                    record.setAttribute("OQ", oqString);
                } else if (oq == Oq.EMPTY_ON_ONE) {
                    record.setAttribute("OQ", "");
                }
                writer.addAlignment(record);
            }
        }
    }

    static void unmark(final Path dir, final Path input, final String label, final String[] extra)
            throws Exception {
        run("UnmarkDuplicates", dir, input, label, extra, argv -> {
            new UnmarkDuplicates().instanceMain(argv);
            return null;
        });
    }

    static void revert(final Path dir, final Path input, final String label, final String[] extra)
            throws Exception {
        run("RevertBaseQualityScores", dir, input, label, extra, argv -> {
            new RevertBaseQualityScores().instanceMain(argv);
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

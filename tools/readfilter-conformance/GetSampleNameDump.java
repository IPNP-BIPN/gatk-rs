/*
 * GetSampleName's output, taken from the reference.
 *
 * The first tool of the reporting-walker archetype measured here, and the smallest: its `traverse`
 * is empty and everything happens in `onTraversalStart`, which reads the header's read groups and
 * writes their sample names to a file.
 *
 * Seven behaviours this is built to catch.
 *
 *   - IT WRITES EVERY DISTINCT SAMPLE, ONE PER LINE, AND NO TRAILING NEWLINE. `Collectors.joining("\n")`
 *     puts a separator BETWEEN names and nothing after the last, so a one-sample file is a single
 *     line with no line ending at all;
 *   - `distinct()` KEEPS THE FIRST OCCURRENCE AND ITS ORDER, so two read groups naming the same
 *     sample give one line, and two samples come out in READ GROUP ORDER rather than sorted;
 *   - A HEADER WITH NO READ GROUPS TAKES THE SECOND REFUSAL, NOT THE FIRST. `getReadGroups()`
 *     returns an EMPTY LIST rather than null, so the guard about "no header or no read groups" is
 *     never reached and the message is "The given bam input has no sample names.";
 *   - A READ GROUP WITH NO SM WRITES THE FOUR LETTERS `null`. `getSample()` returns null,
 *     `distinct()` keeps it, and `Collectors.joining` stringifies it: the file is not empty and the
 *     tool does not refuse;
 *   - --use-url-encoding HAS NO DEFAULT IN THE SOURCE AND IS STILL OPTIONAL: a run that leaves it
 *     out finishes with the field's own `false`;
 *   - URL ENCODING IS java.net.URLEncoder's, so a SPACE BECOMES `+` AND NOT `%20`, and the
 *     characters it leaves alone are its own set;
 *   - AND THE TOOL READS ONLY THE HEADER: `traverse` is empty and `requiresReads` is true, so an
 *     input with no records at all still produces a sample name.
 *
 * Output:
 *
 *     fixture\t<label>\t<the input BAM, base64>
 *     fixtureindex\t<label>\t<the index, base64>
 *     sample\t<label>\t<the file's whole content, escaped>
 *     bytes\t<label>\t<its length in bytes>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: GetSampleNameDump
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
import org.broadinstitute.hellbender.tools.GetSampleName;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GetSampleNameDump {

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final Path dir = Path.of("getsamplename-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# GetSampleNameDump: GetSampleName's output, from the reference");

        // One read group, one sample: the shortest possible output.
        fixture(dir, "single", new String[][] {{"rg1", "sample1"}}, true);
        // Two read groups naming the same sample, which `distinct` collapses.
        fixture(dir, "repeated", new String[][] {{"rg1", "s1"}, {"rg2", "s1"}}, true);
        // Two samples, given in an order that is not alphabetical.
        fixture(dir, "two", new String[][] {{"rg1", "zebra"}, {"rg2", "alpha"}}, true);
        // A sample name with a space and characters URL encoding treats differently.
        fixture(dir, "special", new String[][] {{"rg1", "a sample/with+odd chars & more"}}, true);
        // No read groups at all.
        fixture(dir, "no-read-groups", new String[][] {}, true);
        // A read group with no SM tag.
        fixture(dir, "no-sample", new String[][] {{"rg1", null}}, true);
        // A file with a header and no records at all.
        fixture(dir, "no-reads", new String[][] {{"rg1", "sample1"}}, false);

        for (final String label : new String[] {
                "single", "repeated", "two", "special", "no-read-groups", "no-sample", "no-reads"}) {
            run(dir, label, label, new String[] {"--use-url-encoding", "false"});
        }
        // The same fixtures with encoding on, where the special characters change.
        run(dir, "single-encoded", "single", new String[] {"--use-url-encoding", "true"});
        run(dir, "special-encoded", "special", new String[] {"--use-url-encoding", "true"});
        run(dir, "two-encoded", "two", new String[] {"--use-url-encoding", "true"});
        // And with the argument left out, which finishes: the field's own false is the default.
        run(dir, "no-encoding-argument", "single", new String[] {});
    }

    static void fixture(final Path dir, final String label, final String[][] groups,
                        final boolean withReads) throws Exception {
        final Path bam = dir.resolve(label + ".bam");
        final SAMFileHeader header = header(groups);
        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, bam.toFile())) {
            if (withReads && groups.length > 0) {
                writer.addAlignment(read(header, "r1", 10, groups[0][0]));
            }
        }
        System.out.printf("fixture\t%s\t%s%n", label, RecordTransformDump.base64(bam));
        final Path index = dir.resolve(label + ".bai");
        System.out.printf("fixtureindex\t%s\t%s%n", label,
                Files.exists(index) ? RecordTransformDump.base64(index) : "absent");
    }

    static SAMFileHeader header(final String[][] groups) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 200))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        for (final String[] group : groups) {
            final SAMReadGroupRecord record = new SAMReadGroupRecord(group[0]);
            if (group[1] != null) {
                record.setSample(group[1]);
            }
            header.addReadGroup(record);
        }
        return header;
    }

    static SAMRecord read(final SAMFileHeader header, final String name, final int start,
                          final String readGroup) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setFlags(0);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString("10M");
        record.setReadBases("ACGTACGTAC".getBytes(StandardCharsets.UTF_8));
        record.setBaseQualities(new byte[] {35, 35, 35, 35, 35, 35, 35, 35, 35, 35});
        record.setMappingQuality(60);
        record.setAttribute("RG", readGroup);
        return record;
    }

    /** One run of the tool, with the file it wrote and that file's length. */
    static void run(final Path dir, final String label, final String fixture, final String[] extra)
            throws Exception {
        final Path output = dir.resolve("GetSampleName." + label + ".txt");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-I", dir.resolve(fixture + ".bam").toString(), "-O", output.toString(),
                "--use-jdk-inflater", "true"));
        argv.addAll(Arrays.asList(extra));

        try {
            new GetSampleName().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(), e.getMessage());
            return;
        }

        final byte[] bytes = Files.readAllBytes(output);
        System.out.printf("sample\t%s\t%s%n", label,
                ReferenceQueryDump.escape(new String(bytes, StandardCharsets.UTF_8)));
        // The length is dumped separately because the absence of a trailing newline is the point
        // and an escaped string ending in nothing looks the same as one ending in a newline.
        System.out.printf("bytes\t%s\t%d%n", label, bytes.length);
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
}

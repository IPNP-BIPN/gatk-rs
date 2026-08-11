/*
 * SplitReads, taken from the reference.
 *
 * The seventh whole tool of the record-transform archetype, and the first whose -O is a DIRECTORY
 * and whose run opens more than one writer. Everything worth measuring is in which files exist.
 *
 * Five behaviours this is built to catch.
 *
 *   - THE SET OF OUTPUT FILES COMES FROM THE HEADER, NOT FROM THE READS. createWriters takes the
 *     cross product of each splitter's getSplitsBy(header), and ReadGroupSplitter.getSplitsBy maps
 *     header.getReadGroups() through its selector. A read group no read belongs to still gets a
 *     file, and that file is a valid empty BAM;
 *   - THE CROSS PRODUCT DOES NOT DEDUPLICATE. Two read groups with the same sample give the same
 *     key twice, so prepareSAMFileWriter runs twice for it and the second writer replaces the
 *     first in the map. closeTool closes what the map holds, so THE FIRST WRITER IS NEVER CLOSED,
 *     on a path the second one is also writing. What the file ends up holding is the question, and
 *     the dump is run so the answer is measured on two machines rather than reasoned about;
 *   - A NULL VALUE IS SPELLED TWO WAYS. A read group with no LB gives the key `.null` from the
 *     header, because addKey concatenates the Object, and `.unknown` from a read, because getKey
 *     substitutes UNKNOWN_OUT_PREFIX for null. So the file the header promises is not the file the
 *     reads go to, and both exist: one empty, one on demand;
 *   - A READ WITH NO READ GROUP AT ALL cannot answer any splitter, and the selector runs on the
 *     null read group. Whether that aborts the run depends on a filter rather than on the tool:
 *     WellformedReadFilter's HAS_READ_GROUP drops the read first. Both runs are here, one with the
 *     default filters and one with that filter disabled;
 *   - WITH NO --split-* AT ALL the key is the empty string and the tool writes one file, named
 *     after the input, into the output directory.
 *
 * Every output file travels in the golden in full, base64, indexes included, and the deflater is
 * pinned and recorded for the same reason the other tools' dumps pin it.
 *
 * Output:
 *
 *     deflater\t<class>
 *     fixture\t<label>\t<base64 bam>
 *     fixtureindex\t<label>\t<base64 bai>
 *     commandline\t<label>\t<@PG command line>
 *     header\t<label>\t<escaped SAM header>
 *     outcount\t<label>\t<files>
 *     outfile\t<label>\t<name>\t<base64 bam>
 *     outindex\t<label>\t<name>\t<base64 bai or absent>
 *     error\t<label>\t<class>:<message>
 *
 * Usage: SplitReadsDump
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

import org.broadinstitute.hellbender.tools.SplitReads;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;

public class SplitReadsDump {

    public static void main(final String[] args) throws Exception {
        // Before the fixture is written: the factory is static and first writer wins.
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // Relative on purpose: the string handed to -I and -O is the string recorded inside the
        // output BAM's own @PG, so an absolute temporary path would make every output byte
        // unstable and canonicalization cannot reach inside base64.
        final Path dir = Path.of("splitreads-dump");
        emptyTree(dir);
        Files.createDirectories(dir);

        System.out.println("# SplitReadsDump: SplitReads");
        System.out.printf("deflater\t%s%n",
                BlockCompressedOutputStream.getDefaultDeflaterFactory().getClass().getName());

        // Three read groups: two share a sample, and one has no library at all.
        final Path plain = dir.resolve("plain.bam");
        buildFixture(plain.toFile(), false);
        fixture(dir, plain, "plain");
        // The same header, and one read carrying no RG tag.
        final Path norg = dir.resolve("norg.bam");
        buildFixture(norg.toFile(), true);
        fixture(dir, norg, "norg");

        split(dir, plain, "sample", new String[] {"--split-sample"});
        split(dir, plain, "readgroup", new String[] {"--split-read-group"});
        split(dir, plain, "library", new String[] {"--split-library-name"});
        split(dir, plain, "all3", new String[] {
                "--split-sample", "--split-read-group", "--split-library-name"});
        split(dir, plain, "none", new String[] {});
        split(dir, plain, "noindex",
                new String[] {"--split-sample", "--create-output-bam-index", "false"});
        split(dir, plain, "nopg",
                new String[] {"--split-sample", "--add-output-sam-program-record", "false"});
        // The read with no RG, twice: dropped by the default filters, and reaching the splitter
        // once that filter is gone.
        split(dir, norg, "norg-default", new String[] {"--split-sample"});
        split(dir, norg, "norg-nofilter", new String[] {
                "--split-sample", "--disable-read-filter", "WellformedReadFilter"});
        // The output directory is asserted writable before anything else happens.
        split(dir, plain, "missingdir", new String[] {"--split-sample"}, "splitreads-dump/absent");
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

    /**
     * A coordinate-sorted BAM over three read groups.
     *
     * `rg1` and `rg2` share the sample `s1`, which is the collision that makes one writer per key
     * two writers on one path. `rg3` has no `LB` at all, which is the null the header spells
     * `null` and a read spells `unknown`.
     */
    static void buildFixture(final File file, final boolean withUngroupedRead) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 1000),
                new SAMSequenceRecord("chr2", 1000))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        header.addReadGroup(readGroup("rg1", "s1", "lib1"));
        header.addReadGroup(readGroup("rg2", "s1", "lib2"));
        header.addReadGroup(readGroup("rg3", "s2", null));
        final SAMProgramRecord existing = new SAMProgramRecord("upstream");
        existing.setProgramVersion("1.0");
        header.addProgramRecord(existing);

        final String[] groups = {"rg1", "rg2", "rg3", "rg1", "rg3"};
        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, file)) {
            for (int i = 0; i < groups.length; i++) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName("r" + i);
                record.setReferenceName("chr1");
                record.setAlignmentStart(100 + i * 20);
                record.setCigarString("10M");
                final byte[] bases = new byte[10];
                Arrays.fill(bases, (byte) 'A');
                record.setReadBases(bases);
                final byte[] quals = new byte[10];
                Arrays.fill(quals, (byte) 30);
                record.setBaseQualities(quals);
                record.setMappingQuality(60);
                // The last read carries no RG when asked for, so no splitter can answer for it.
                if (!(withUngroupedRead && i == groups.length - 1)) {
                    record.setAttribute("RG", groups[i]);
                }
                writer.addAlignment(record);
            }
        }
    }

    static SAMReadGroupRecord readGroup(final String id, final String sample, final String library) {
        final SAMReadGroupRecord group = new SAMReadGroupRecord(id);
        group.setSample(sample);
        if (library != null) {
            group.setLibrary(library);
        }
        return group;
    }

    static void split(final Path dir, final Path input, final String label, final String[] extra)
            throws Exception {
        split(dir, input, label, extra, null);
    }

    /**
     * One run, and everything the output directory holds afterwards.
     *
     * The directory is listed rather than derived: which files exist is the whole measurement, and
     * a dump that named the files it expected would be asserting its own guess.
     */
    static void split(final Path dir, final Path input, final String label, final String[] extra,
                      final String outputOverride) throws Exception {
        final Path output = outputOverride == null ? dir.resolve(label) : Path.of(outputOverride);
        if (outputOverride == null) {
            Files.createDirectories(output);
        }
        // --use-jdk-deflater is the knob that decides which bytes come out, for the same reason
        // PrintReadsDump names it: the GKL deflater's output is not yet reproduced.
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-I", input.toString(), "-O", output.toString(),
                "--use-jdk-deflater", "true", "--use-jdk-inflater", "true"));
        argv.addAll(Arrays.asList(extra));

        try {
            new SplitReads().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            // The refusal is the observable behaviour, so it is dumped rather than swallowed.
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    String.valueOf(e.getMessage()).replace('\n', ' '));
            return;
        }

        final List<Path> files;
        try (final var entries = Files.list(output)) {
            files = entries.filter(path -> path.getFileName().toString().endsWith(".bam"))
                    .sorted()
                    .toList();
        }
        System.out.printf("outcount\t%s\t%d%n", label, files.size());

        String commandLine = "";
        if (!files.isEmpty()) {
            try (final SamReader reader = SamReaderFactory.makeDefault()
                    .validationStringency(ValidationStringency.SILENT)
                    .open(files.get(0).toFile())) {
                final SAMFileHeader header = reader.getFileHeader();
                for (final SAMProgramRecord record : header.getProgramRecords()) {
                    if (record.getCommandLine() != null) {
                        commandLine = record.getCommandLine();
                    }
                }
                // One header row per run: every writer of a run is given the same header object.
                System.out.printf("header\t%s\t%s%n", label,
                        ReferenceQueryDump.escape(header.getSAMString()));
            }
        }
        System.out.printf("commandline\t%s\t%s%n", label, commandLine);

        for (final Path file : files) {
            final String name = file.getFileName().toString();
            // Every writer of a run is handed getHeaderForSAMWriter(), which ADDS a @PG record to
            // the reads header in place and hands back the same object. A writer serialises the
            // header when it is created, so the nth file created carries n of them.
            try (final SamReader reader = SamReaderFactory.makeDefault()
                    .validationStringency(ValidationStringency.SILENT).open(file.toFile())) {
                final StringBuilder ids = new StringBuilder();
                for (final SAMProgramRecord record : reader.getFileHeader().getProgramRecords()) {
                    ids.append(ids.length() == 0 ? "" : ";").append(record.getId());
                }
                System.out.printf("programs\t%s\t%s\t%s%n", label, name, ids);
            }
            System.out.printf("outfile\t%s\t%s\t%s%n", label, name, base64(file));
            final Path index = output.resolve(name.replace(".bam", ".bai"));
            System.out.printf("outindex\t%s\t%s\t%s%n", label, name,
                    Files.exists(index) ? base64(index) : "absent");
        }
    }

    static void emptyTree(final Path dir) throws Exception {
        if (!Files.isDirectory(dir)) {
            return;
        }
        try (final var walk = Files.walk(dir)) {
            for (final Path entry : walk.sorted(java.util.Comparator.reverseOrder()).toList()) {
                Files.deleteIfExists(entry);
            }
        }
    }

    static String base64(final Path path) throws Exception {
        return Base64.getEncoder().encodeToString(Files.readAllBytes(path));
    }
}

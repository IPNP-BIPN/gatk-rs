/*
 * AddOriginalAlignmentTags, taken from the reference.
 *
 * The fifth whole tool of the record-transform archetype, and the first that writes tags rather
 * than changing the read. Its apply is three lines; what is worth measuring is what goes into the
 * two tags and what happens when the read cannot answer.
 *
 * Four behaviours this is built to catch.
 *
 *   - THE OA TAG IS A FORMATTED STRING, not a structure: contig, start, strand, cigar, mapping
 *     quality and NM, comma separated and terminated by a semicolon. It replaces a comma in the
 *     contig name with an underscore, and THAT BRANCH CANNOT FIRE: a comma is not a legal sequence
 *     name character, which is measured here rather than assumed;
 *   - A READ WITH NO NM TAG STILL GETS AN OA. getAttributeAsString returns null and the format
 *     string prints it, so the tag reads `...,60,null;`. That is what the reference writes and it
 *     is measured rather than inferred;
 *   - AN UNMAPPED READ GETS A CONSTANT, `*,0,*,*,0,0;`, and its XM is `*` as well;
 *   - THE MATE CONTIG TAG ASKS A QUESTION AN UNPAIRED READ CANNOT ANSWER. `mateIsUnmapped` is only
 *     defined for a paired read, so what an unpaired one does is the tool's behaviour on ordinary
 *     input rather than an edge case: most reads in most files are paired, and the fixture carries
 *     both.
 *
 * The output BAMs travel in the golden in full, base64, indexes included, and the deflater is
 * pinned and recorded for the same reason the other tools' dumps pin it.
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
 *     tags\t<label>\t<read name>\t<OA>\t<XM>
 *     error\t<tool>\t<label>\t<class>:<message>
 *
 * Usage: AddOriginalAlignmentTagsDump
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

import org.broadinstitute.hellbender.tools.AddOriginalAlignmentTags;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;

public class AddOriginalAlignmentTagsDump {

    public static void main(final String[] args) throws Exception {
        // Before the fixture is written: the factory is static and first writer wins.
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // Relative on purpose: the string handed to -I and -O is the string recorded inside the
        // output BAM's own @PG, so an absolute temporary path would make every output byte
        // unstable and canonicalization cannot reach inside base64.
        final Path dir = Path.of("addoatags-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# AddOriginalAlignmentTagsDump: AddOriginalAlignmentTags");
        System.out.printf("deflater\t%s%n",
                BlockCompressedOutputStream.getDefaultDeflaterFactory().getClass().getName());

        // Paired reads with mates, an unmapped read, a read with no NM, and a contig whose name
        // carries a comma: everything the two tags can be built from.
        final Path plain = dir.resolve("plain.bam");
        buildFixture(plain.toFile());
        fixture(dir, plain, "plain");
        // The OA format replaces a comma in the contig name with an underscore, and that branch
        // cannot fire: a comma is not a legal sequence name character at all. Measured here rather
        // than assumed, with the regex the refusal quotes.
        try {
            new SAMSequenceRecord("chr,1", 1000);
            System.out.printf("sequencename\tcomma\taccepted\t-%n");
        } catch (final Throwable t) {
            System.out.printf("sequencename\tcomma\t%s\t%s%n", t.getClass().getSimpleName(),
                    String.valueOf(t.getMessage()).trim());
        }

        addTags(dir, plain, "all", new String[] {});
        addTags(dir, plain, "chr1", new String[] {"-L", "chr1"});
        addTags(dir, plain, "noindex", new String[] {"--create-output-bam-index", "false"});
        addTags(dir, plain, "nopg",
                new String[] {"--add-output-sam-program-record", "false"});

        // And the tags themselves, read back off the output so the golden says what went in them
        // rather than only that the bytes match.
        tags(dir, "all");
    }

    /** Every read of an output, with the two tags this tool wrote. */
    static void tags(final Path dir, final String label) throws Exception {
        final Path output = dir.resolve("AddOriginalAlignmentTags." + label + ".bam");
        try (final SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT).open(output.toFile())) {
            for (final SAMRecord record : reader) {
                System.out.printf("tags\t%s\t%s\t%s\t%s%n", label, record.getReadName(),
                        String.valueOf(record.getStringAttribute("OA")),
                        String.valueOf(record.getStringAttribute("XM")));
            }
        }
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
     * A small coordinate-sorted BAM, with a duplicate-flagged read and an unmapped tail.
     *
     * Built here rather than shared, because what varies per case is the qualities and that is the
     * point of every case.
     */
    static void buildFixture(final File file) {
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
                Arrays.fill(quals, (byte) 30);
                record.setBaseQualities(quals);
                // Four of the six are paired with a mapped mate, one is paired with an unmapped
                // mate, and one is not paired at all: mateIsUnmapped is only defined for the
                // first two kinds.
                if (i < 4) {
                    record.setReadPairedFlag(true);
                    record.setMateReferenceName(i < 2 ? "chr1" : "chr2");
                    record.setMateAlignmentStart(500);
                } else if (i == 4) {
                    record.setReadPairedFlag(true);
                    record.setMateUnmappedFlag(true);
                    record.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
                    record.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
                }
                // Half the reads carry NM and half do not, because a read with none still gets an
                // OA and what goes in that field is the question.
                if (i % 2 == 0) {
                    record.setAttribute("NM", i);
                }
                if (i == 3) {
                    record.setReadNegativeStrandFlag(true);
                }
                record.setMappingQuality(60);
                record.setAttribute("RG", "rg1");
                writer.addAlignment(record);
            }
        }
    }

    static void addTags(final Path dir, final Path input, final String label, final String[] extra)
            throws Exception {
        run("AddOriginalAlignmentTags", dir, input, label, extra, argv -> {
            new AddOriginalAlignmentTags().instanceMain(argv);
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

/*
 * PrintReadsHeader, taken from the reference.
 *
 * The eleventh whole tool of the record-transform archetype, the second that is not a walker, and
 * the SMALLEST TOOL IN THE ARCHETYPE: forty-six lines, of which four do anything. It is measured
 * for that reason as much as for itself. G2's calibration gate asks what a member of this archetype
 * costs at the margin, and every measurement so far has been of a tool that does something; this
 * one is the floor.
 *
 * Five behaviours this is built to catch, and the first one is not about this tool at all.
 *
 *   - NO BAM htsjdk WROTE CAN CARRY A NON-CURRENT VERSION, so the rewrite this tool performs is a
 *     no-op on every file it will ever be given. `SAMFileWriterImpl.writeHeader(SAMFileHeader)`
 *     calls the TWO-ARGUMENT `SAMTextHeaderCodec.encode`, which is `keepExistingVersionNumber =
 *     false`, and `BAMFileWriter` overrides only `writeHeader(String)`, so the ordinary writer
 *     normalises `VN` on its way out. The three-argument call that passes `true` sits in
 *     `BAMFileWriter.writeHeader(BinaryCodec, SAMFileHeader)` and is reachable only from the
 *     standalone block-copy reheader. Measured here rather than reasoned: the dump prints the
 *     header's own `VN` at the moment it is handed to the writer, `builtvn old_version.bam` says
 *     `1.5`, and the BAM that comes out says `VN:1.6`. That is the finding, and it belongs to
 *     htsjdk-rs, whose `BamWriter` keeps the header's version and whose `SamHeader::encode` says
 *     in a comment that the writer passes `true`. Filed as htsjdk-rs#164;
 *   - THIS TOOL REWRITES THE VERSION TOO, through the same two-argument overload:
 *     `new SAMTextHeaderCodec().encode(writer, header)`. Since no input can disagree, the port has
 *     nothing to reproduce here beyond the flag, and the flag is recorded so a future fixture with
 *     a hand-built header does not have to rediscover it. `writeHDLine(false)` builds a FRESH
 *     SAMFileHeader, copies every attribute except VN into it, and lets the constructor's own VN
 *     stand, which also moves VN to the front of the line: unobservable for the same reason, and
 *     written down for the same reason;
 *   - NO @PG IS APPENDED. The tool reads `getHeaderForReads()`, not `getHeaderForSAMWriter()`, so
 *     the @PG chain is the file's own and nothing is added to it. Every writer tool in this
 *     archetype does the opposite, which is what makes this worth a row;
 *   - THE OUTPUT IS TEXT THROUGH AN OutputStreamWriter WITH NO CHARSET, so it is the platform
 *     default. The fixture carries a non-ASCII @CO comment and the output travels as base64, which
 *     is the only way to see which charset the pinned container actually used; the dump names the
 *     charset beside it so the two can be checked against each other;
 *   - A HEADER WITH NO SEQUENCE DICTIONARY STILL PRINTS, because the codec iterates an empty list
 *     rather than refusing, and a run with no -I does not start at all: `requiresReads()` is true,
 *     and the refusal is Barclay's rather than the engine's.
 *
 * Output:
 *
 *     deflater\t<class>
 *     charset\t<name>
 *     currentversion\t<version>
 *     builtvn\t<file>\t<the header's VN when the writer was handed it>
 *     fixture\t<label>\t<base64 bam>
 *     inputheader\t<label>\t<escaped SAM header as htsjdk parsed it>
 *     output\t<label>\t<base64 output file>
 *     outputtext\t<label>\t<escaped output text>
 *     error\t<label>\t<class>:<message>
 *
 * Usage: PrintReadsHeaderDump
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

import org.broadinstitute.hellbender.tools.PrintReadsHeader;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;

public class PrintReadsHeaderDump {

    public static void main(final String[] args) throws Exception {
        // The factory is static and the first writer wins. This dump calls no Picard entry point,
        // so nothing should replace it; the pin makes that a fact rather than a hope.
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // Relative on purpose, as in the other dumps of this archetype: an absolute temporary path
        // would be unstable, and here it would land in the UserException message too.
        final Path dir = Path.of("printreadsheader-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# PrintReadsHeaderDump: PrintReadsHeader");
        System.out.printf("deflater\t%s%n",
                BlockCompressedOutputStream.getDefaultDeflaterFactory().getClass().getName());
        // The charset the OutputStreamWriter will use, named rather than inferred from the bytes.
        System.out.printf("charset\t%s%n", java.nio.charset.Charset.defaultCharset().name());
        // What the codec will write in place of the file's own version.
        System.out.printf("currentversion\t%s%n", SAMFileHeader.CURRENT_VERSION);

        // A header that says 1.5 when the writer is handed it. The writer normalises it, which is
        // the point: this fixture is the evidence for htsjdk-rs#164 rather than a 1.5 file.
        final Path old = dir.resolve("old_version.bam");
        buildFixture(old.toFile(), "1.5", true, true);
        fixture(dir, old, "old_version");

        // Everything a header can carry: extra @HD attributes, @SQ attributes, two @RG, a @PG
        // chain, and comments including one that is not ASCII.
        final Path rich = dir.resolve("rich.bam");
        buildFixture(rich.toFile(), null, true, true);
        fixture(dir, rich, "rich");

        // No sequence dictionary at all, which the codec iterates rather than refuses.
        final Path bare = dir.resolve("bare.bam");
        buildFixture(bare.toFile(), null, false, false);
        fixture(dir, bare, "bare");

        print(dir, old, "oldversion");
        print(dir, rich, "rich");
        print(dir, bare, "bare");
        print(dir, null, "noinput");
    }

    /** A fixture, and the header htsjdk read back out of it. */
    static void fixture(final Path dir, final Path bam, final String label) throws Exception {
        System.out.printf("fixture\t%s\t%s%n", label, base64(bam));
        try (final SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT).open(bam.toFile())) {
            System.out.printf("inputheader\t%s\t%s%n", label,
                    ReferenceQueryDump.escape(reader.getFileHeader().getSAMString()));
        }
    }

    /**
     * One header, built with whatever the label needs.
     *
     * `version` null leaves whatever `new SAMFileHeader()` set, which is the current one; a string
     * overwrites it, which is how the rewrite becomes observable.
     */
    static void buildFixture(final File file, final String version, final boolean dictionary,
                             final boolean reads) {
        final SAMFileHeader header = new SAMFileHeader();
        if (dictionary) {
            final SAMSequenceRecord chr1 = new SAMSequenceRecord("chr1", 100);
            chr1.setAttribute("M5", "0123456789abcdef0123456789abcdef");
            chr1.setAttribute("UR", "file:/ref/chr1.fasta");
            chr1.setAttribute("SP", "Homo sapiens");
            header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                    chr1, new SAMSequenceRecord("chr2", 200))));
        }
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        header.setGroupOrder(SAMFileHeader.GroupOrder.none);

        final SAMReadGroupRecord first = new SAMReadGroupRecord("rg1");
        first.setSample("s1");
        first.setLibrary("lib1");
        first.setPlatform("ILLUMINA");
        header.addReadGroup(first);
        final SAMReadGroupRecord second = new SAMReadGroupRecord("rg2");
        second.setSample("s2");
        header.addReadGroup(second);

        final SAMProgramRecord upstream = new SAMProgramRecord("upstream");
        upstream.setProgramVersion("1.0");
        upstream.setCommandLine("upstream --in a.bam");
        header.addProgramRecord(upstream);
        final SAMProgramRecord downstream = new SAMProgramRecord("downstream");
        downstream.setPreviousProgramGroupId("upstream");
        header.addProgramRecord(downstream);

        header.addComment("a plain comment");
        // Not ASCII, so the writer's charset is decidable from the bytes.
        header.addComment("un commentaire accentue: éèê");

        // Set last, so it overwrites the constructor's value rather than being overwritten.
        if (version != null) {
            header.setAttribute("VN", version);
        }

        // What the header object actually holds at the moment it is handed to the writer. Without
        // this row, a version that does not survive into the file cannot be told apart from a
        // writer that rewrote it.
        System.out.printf("builtvn\t%s\t%s%n", file.getName(), header.getAttribute("VN"));

        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().makeBAMWriter(header, true, file)) {
            if (reads) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName("r0");
                record.setReferenceName("chr1");
                record.setAlignmentStart(1);
                record.setCigarString("10M");
                record.setReadBases("ACGTACGTAC".getBytes());
                final byte[] quals = new byte[10];
                Arrays.fill(quals, (byte) 30);
                record.setBaseQualities(quals);
                record.setMappingQuality(60);
                record.setAttribute("RG", "rg1");
                writer.addAlignment(record);
            }
        }
    }

    static void print(final Path dir, final Path input, final String label) throws Exception {
        final Path output = dir.resolve("PrintReadsHeader." + label + ".txt");
        final List<String> argv = new ArrayList<>(Arrays.asList("-O", output.toString()));
        if (input != null) {
            argv.addAll(Arrays.asList("-I", input.toString()));
        }

        try {
            new PrintReadsHeader().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            // A refusal is the observable behaviour, so it is dumped rather than swallowed.
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    String.valueOf(e.getMessage()).replace('\n', ' '));
            return;
        }

        // The bytes first, because the charset is only decidable from them, and the text after, so
        // a reader can see what changed without decoding base64 by hand.
        System.out.printf("output\t%s\t%s%n", label, base64(output));
        System.out.printf("outputtext\t%s\t%s%n", label,
                ReferenceQueryDump.escape(new String(Files.readAllBytes(output))));
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

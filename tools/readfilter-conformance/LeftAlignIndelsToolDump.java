/*
 * LeftAlignIndels, taken from the reference.
 *
 * The eighth whole tool of the record-transform archetype, and the FIRST WHOSE requiresReference()
 * IS TRUE. Its apply is eight lines around one call, and the call is measured on its own by
 * LeftAlignIndelsDump; what is left for this dump is what the tool does with the answer.
 *
 * Four behaviours this is built to catch.
 *
 *   - THE REFERENCE A READ WALKER SEES IS THE READ'S OWN SPAN, no padding, so left-alignment is
 *     bounded by the read's own window rather than by the contig. The window handed to each read is
 *     dumped beside the read so a port cannot quietly widen it;
 *   - A LEFT-ALIGNED DELETION CAN MOVE THE READ. CigarBuilder drops a deletion that ends up leading
 *     and reports the reference bases it removed; the tool then moves the read right by that many.
 *     This is the second tool of the archetype that changes a read's position, after
 *     PrintDistantMates, and it does it for an unrelated reason;
 *   - TWO CLASSES OF READ ARE PASSED THROUGH UNTOUCHED, before the call rather than by it: an
 *     unmapped read, and a read whose cigar has one element or none. The fixture carries both;
 *   - A RUN WITH NO -R IS REFUSED, because requiresReference() is true, and the refusal is the
 *     engine's rather than the tool's.
 *
 * The output BAM travels in the golden in full, base64, index included, and the deflater is pinned
 * and recorded for the same reason the other tools' dumps pin it.
 *
 * Output:
 *
 *     deflater\t<class>
 *     fasta\t<escaped FASTA text>
 *     fai\t<escaped .fai text>
 *     fixture\t<label>\t<base64 bam>
 *     fixtureindex\t<label>\t<base64 bai>
 *     refwindow\t<read name>\t<contig>\t<start>\t<end>\t<bases>
 *     header\t<label>\t<escaped SAM header>
 *     commandline\t<label>\t<@PG command line>
 *     output\t<label>\t<base64 bam>
 *     index\t<label>\t<base64 bai or absent>
 *     reads\t<label>\t<name>\t<flags>\t<contig>\t<start>\t<cigar>
 *     error\t<label>\t<class>:<message>
 *
 * Usage: LeftAlignIndelsToolDump
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
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.tools.LeftAlignIndels;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;

public class LeftAlignIndelsToolDump {

    /**
     * Sixty bases in three runs and one two-base repeat, so an indel has somewhere to go.
     *
     * 1-10 A, 11-20 T, 21-30 the AC repeat, 31-40 G, 41-50 C, 51-60 T.
     */
    static final String FASTA =
            ">chr1 repeats to left-align into\n"
            + "AAAAAAAAAATTTTTTTTTT\n"
            + "ACACACACACGGGGGGGGGG\n"
            + "CCCCCCCCCCTTTTTTTTTT\n";

    public static void main(final String[] args) throws Exception {
        // Before the fixture is written: the factory is static and first writer wins.
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // Relative on purpose: the string handed to -I and -O is the string recorded inside the
        // output BAM's own @PG, so an absolute temporary path would make every output byte
        // unstable and canonicalization cannot reach inside base64.
        final Path dir = Path.of("leftalign-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# LeftAlignIndelsToolDump: LeftAlignIndels");
        System.out.printf("deflater\t%s%n",
                BlockCompressedOutputStream.getDefaultDeflaterFactory().getClass().getName());

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        final Path fai = dir.resolve("ref.fasta.fai");
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});
        // AND AGAIN, because that call replaced it. Picard's CommandLineProgram installs its own
        // deflater factory into the same static this dump set at the top, so every BAM written
        // after it comes out with the GKL deflater's bytes rather than the JDK's, whatever
        // --use-jdk-deflater says on the GATK command line. Measured: without this line the output
        // BAMs are four to eleven bytes longer and diverge from the JDK deflater's stream three
        // bytes in.
        System.out.printf("deflaterafterdict\t%s%n",
                BlockCompressedOutputStream.getDefaultDeflaterFactory().getClass().getName());
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());
        System.out.printf("deflaterrestored\t%s%n",
                BlockCompressedOutputStream.getDefaultDeflaterFactory().getClass().getName());
        System.out.printf("fasta\t%s%n",
                ReferenceQueryDump.escape(new String(Files.readAllBytes(fasta))));
        System.out.printf("fai\t%s%n",
                ReferenceQueryDump.escape(new String(Files.readAllBytes(fai))));

        final Path plain = dir.resolve("plain.bam");
        buildFixture(plain.toFile());
        fixture(dir, plain, "plain");

        // The window each read's apply is given, which is the read's own span and nothing more.
        referenceWindows(fasta, plain);

        leftAlign(dir, plain, fasta, "all", new String[] {});
        leftAlign(dir, plain, fasta, "chr1head", new String[] {"-L", "chr1:1-20"});
        leftAlign(dir, plain, fasta, "noindex",
                new String[] {"--create-output-bam-index", "false"});
        leftAlign(dir, plain, fasta, "nopg",
                new String[] {"--add-output-sam-program-record", "false"});
        // requiresReference() is true, so this one never starts.
        leftAlign(dir, plain, null, "noreference", new String[] {});

        reads(dir, "all");
        reads(dir, "chr1head");
    }

    /**
     * The reference bases a ReadWalker hands each read's apply: `new ReferenceContext(reference,
     * new SimpleInterval(read))`, which is the read's own span with no padding.
     *
     * Dumped because it is the bound on how far an indel can move, and a port that queried the
     * contig instead would left-align further than the reference does.
     */
    static void referenceWindows(final Path fasta, final Path bam) throws Exception {
        try (final ReferenceDataSource reference = new ReferenceFileSource(fasta);
             final SamReader reader = SamReaderFactory.makeDefault()
                     .validationStringency(ValidationStringency.SILENT).open(bam.toFile())) {
            for (final SAMRecord record : reader) {
                if (record.getReadUnmappedFlag()) {
                    System.out.printf("refwindow\t%s\t*\t0\t0\t-%n", record.getReadName());
                    continue;
                }
                final SimpleInterval span = new SimpleInterval(record.getReferenceName(),
                        record.getAlignmentStart(), record.getAlignmentEnd());
                System.out.printf("refwindow\t%s\t%s\t%d\t%d\t%s%n", record.getReadName(),
                        span.getContig(), span.getStart(), span.getEnd(),
                        new String(reference.queryAndPrefetch(span).getBases()));
            }
        }
    }

    /** Every read of an output, with what the tool did to its cigar and its position. */
    static void reads(final Path dir, final String label) throws Exception {
        final Path output = dir.resolve("LeftAlignIndels." + label + ".bam");
        try (final SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT).open(output.toFile())) {
            for (final SAMRecord record : reader) {
                System.out.printf("reads\t%s\t%s\t%d\t%s\t%d\t%s%n", label, record.getReadName(),
                        record.getFlags(), record.getReferenceName(), record.getAlignmentStart(),
                        record.getCigarString());
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
     * Six reads over the repeats: one whose deletion walks off the front, one whose insertion does
     * not, one that cannot move at all, one with two indels, and the two kinds the tool passes
     * through before it calls anything.
     */
    static void buildFixture(final File file) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 60))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("s1");
        header.addReadGroup(group);
        final SAMProgramRecord existing = new SAMProgramRecord("upstream");
        existing.setProgramVersion("1.0");
        header.addProgramRecord(existing);

        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, file)) {
            // A deleted A inside the A homopolymer: it left-aligns to the front of the read's own
            // window, is dropped, and the read moves right.
            writer.addAlignment(read(header, "r0", 6, "4M1D5M", "AAAATTTTT"));
            // Two indels close together, which the right-to-left walk may merge.
            writer.addAlignment(read(header, "r1", 6, "3M1D2M1D4M", "AAAATTTTT"));
            // A cigar with one element: passed through before anything is called.
            writer.addAlignment(read(header, "r2", 11, "10M", "TTTTTTTTTT"));
            // A deleted A whose left neighbour is a T: it cannot move.
            writer.addAlignment(read(header, "r3", 17, "4M1D5M", "TTTTCACAC"));
            // An inserted AC inside the AC repeat: the insertion moves but the read does not.
            writer.addAlignment(read(header, "r4", 21, "4M2I4M", "ACACACACAC"));
            // Unmapped, which is the other pass-through. Placed on chr1 so it sorts with the rest.
            final SAMRecord unmapped = read(header, "r5", 41, "10M", "CCCCCCCCCC");
            unmapped.setReadUnmappedFlag(true);
            unmapped.setMappingQuality(0);
            writer.addAlignment(unmapped);
        }
    }

    static SAMRecord read(final SAMFileHeader header, final String name, final int start,
                          final String cigar, final String bases) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setReadBases(bases.getBytes());
        final byte[] quals = new byte[bases.length()];
        Arrays.fill(quals, (byte) 30);
        record.setBaseQualities(quals);
        record.setMappingQuality(60);
        record.setAttribute("RG", "rg1");
        return record;
    }

    static void leftAlign(final Path dir, final Path input, final Path fasta, final String label,
                          final String[] extra) throws Exception {
        final Path output = dir.resolve("LeftAlignIndels." + label + ".bam");
        // --use-jdk-deflater is the knob that decides which bytes come out, for the same reason
        // PrintReadsDump names it: the GKL deflater's output is not yet reproduced.
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-I", input.toString(), "-O", output.toString(),
                "--use-jdk-deflater", "true", "--use-jdk-inflater", "true"));
        if (fasta != null) {
            argv.addAll(Arrays.asList("-R", fasta.toString()));
        }
        argv.addAll(Arrays.asList(extra));

        try {
            new LeftAlignIndels().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            // The refusal is the observable behaviour, so it is dumped rather than swallowed.
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

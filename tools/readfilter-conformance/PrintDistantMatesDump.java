/*
 * PrintDistantMates, taken from the reference.
 *
 * The sixth whole tool of the record-transform archetype, and the first that moves a read rather
 * than editing it in place: the read is rewritten into its mate's position, unmapped, with its old
 * alignment kept in an OA tag.
 *
 * Four behaviours this is built to catch.
 *
 *   - THE DEFAULT READ FILTERS ARE EXTENDED, NOT TAKEN OR REPLACED. getDefaultReadFilters starts
 *     from super's and adds PAIRED, PRIMARY_LINE, NOT_DUPLICATE and MATE_DISTANT, which is a third
 *     pattern inside one archetype. The two lists are dumped side by side rather than described;
 *   - THE OA TAG IS WRITTEN WITH A DIFFERENT CONVENTION FROM AddOriginalAlignmentTags. Both build
 *     `contig,start,strand,cigar,mapq,NM;`, and a missing NM prints as EMPTY here where the other
 *     tool prints the four characters `null`. Both tools are run over the same fixture so the two
 *     spellings of the same missing field land in the same golden, on the same read;
 *   - THE WRITER IS OPENED WITH preSorted = false. Every other tool of this archetype passes true.
 *     The transform moves each read to its mate, so the traversal order is not the output order,
 *     and whether htsjdk re-sorts is the question. The fixture is built so the two orders differ:
 *     the reads leave the traversal chr2:600, chr1:2500, chr2:150;
 *   - undoDistantMateAlterations CLAIMS TO BE THE INVERSE. It is measured as a round trip against
 *     the original record's own SAM text, field by field, including the tag block's order, and on
 *     an OA it cannot parse.
 *
 * The output BAMs travel in the golden in full, base64, indexes included, and the deflater is
 * pinned and recorded for the same reason the other tools' dumps pin it.
 *
 * Output:
 *
 *     deflater\t<class>
 *     filters\t<tool>\t<index>\t<filter class>
 *     fixture\t<label>\t<base64 bam>
 *     fixtureindex\t<label>\t<base64 bai>
 *     header\t<tool>\t<label>\t<escaped SAM header>
 *     commandline\t<tool>\t<label>\t<@PG command line>
 *     output\t<tool>\t<label>\t<base64 bam>
 *     index\t<tool>\t<label>\t<base64 bai or absent>
 *     reads\t<tool>\t<label>\t<name>\t<flags>\t<contig>\t<start>\t<cigar>\t<mapq>\t<OA>\t<DM>\t<NM>
 *     roundtrip\t<name>\t<same|differs>\t<original SAM>\t<recovered SAM>
 *     isdistant\t<name>\t<before>\t<after>\t<undone>
 *     undo\t<label>\t<result>
 *     error\t<tool>\t<label>\t<class>:<message>
 *
 * Usage: PrintDistantMatesDump
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
import org.broadinstitute.hellbender.tools.AddOriginalAlignmentTags;
import org.broadinstitute.hellbender.tools.PrintDistantMates;
import org.broadinstitute.hellbender.tools.PrintReads;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;

public class PrintDistantMatesDump {

    public static void main(final String[] args) throws Exception {
        // Before the fixture is written: the factory is static and first writer wins.
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // Relative on purpose: the string handed to -I and -O is the string recorded inside the
        // output BAM's own @PG, so an absolute temporary path would make every output byte
        // unstable and canonicalization cannot reach inside base64.
        final Path dir = Path.of("distantmates-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# PrintDistantMatesDump: PrintDistantMates");
        System.out.printf("deflater\t%s%n",
                BlockCompressedOutputStream.getDefaultDeflaterFactory().getClass().getName());

        // The third pattern of the archetype, dumped rather than described: PrintReads takes the
        // default list, this tool starts from the same list and adds four filters to it.
        filters("PrintReads", new PrintReads().getDefaultReadFilters());
        filters("PrintDistantMates", new PrintDistantMates().getDefaultReadFilters());

        final Path plain = dir.resolve("plain.bam");
        buildFixture(plain.toFile());
        fixture(dir, plain, "plain");

        distantMates(dir, plain, "all", new String[] {});
        // r5 leaves the traversal last and lands first: whether the output is in traversal order
        // or in coordinate order is what preSorted = false decides.
        distantMates(dir, plain, "chr1head", new String[] {"-L", "chr1:1-250"});
        distantMates(dir, plain, "noindex",
                new String[] {"--create-output-bam-index", "false"});
        distantMates(dir, plain, "nopg",
                new String[] {"--add-output-sam-program-record", "false"});
        // The one filter of the four that takes an argument. At 5 the near-mate read survives, and
        // it is the read with no NM, so the empty field lands in an output BAM rather than only in
        // a static call.
        distantMates(dir, plain, "distant5",
                new String[] {"--mate-too-distant-length", "5"});

        // The same tag, written by the other tool, over the same reads. -L chr1 keeps the unpaired
        // read out of it: AddOriginalAlignmentTags aborts the run on one.
        addOaTags(dir, plain, "chr1", new String[] {"-L", "chr1"});

        for (final String label : new String[] {"all", "chr1head", "noindex", "nopg", "distant5"}) {
            reads("PrintDistantMates", dir, label);
        }
        reads("AddOriginalAlignmentTags", dir, "chr1");

        roundTrip(plain);
        undoRefusals();
    }

    /** A tool's default read filter list, in order. */
    static void filters(final String tool, final List<ReadFilter> filters) {
        for (int i = 0; i < filters.size(); i++) {
            System.out.printf("filters\t%s\t%d\t%s%n", tool, i,
                    filters.get(i).getClass().getSimpleName());
        }
    }

    /**
     * Every read of an output, in the order it is written, with what the transform did to it.
     *
     * The order is the measurement: preSorted = false hands the records to a sorting collection,
     * and the fixture is built so that traversal order and coordinate order disagree.
     */
    static void reads(final String tool, final Path dir, final String label) throws Exception {
        final Path output = dir.resolve(tool + "." + label + ".bam");
        if (!Files.exists(output)) {
            System.out.printf("reads\t%s\t%s\tabsent%n", tool, label);
            return;
        }
        try (final SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT).open(output.toFile())) {
            for (final SAMRecord record : reader) {
                System.out.printf("reads\t%s\t%s\t%s\t%d\t%s\t%d\t%s\t%d\t%s\t%s\t%s%n",
                        tool, label, record.getReadName(), record.getFlags(),
                        record.getReferenceName(), record.getAlignmentStart(),
                        record.getCigarString(), record.getMappingQuality(),
                        String.valueOf(record.getStringAttribute("OA")),
                        String.valueOf(record.getStringAttribute("DM")),
                        String.valueOf(record.getAttribute("NM")));
            }
        }
    }

    /**
     * doDistantMateAlterations and undoDistantMateAlterations over every read of the fixture.
     *
     * Called as statics rather than through a run, because the question is whether the pair is an
     * inverse and a run only shows the forward half. The comparison is the record's whole SAM
     * text, which includes the order of the tag block: a tag cleared and set again does not
     * necessarily come back where it was.
     */
    static void roundTrip(final Path bam) throws Exception {
        try (final SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT).open(bam.toFile())) {
            final SAMFileHeader header = reader.getFileHeader();
            for (final SAMRecord record : reader) {
                final GATKRead read = new SAMRecordToGATKReadAdapter(record.deepCopy());
                final String before = record.getSAMString().trim();
                // Only what the tool's own filters would let through is altered: the transform
                // reads the mate contig, which an unpaired read does not have.
                if (!record.getReadPairedFlag() || record.getMateReferenceIndex() < 0) {
                    System.out.printf("roundtrip\t%s\tunaltered\t%s\t-%n",
                            record.getReadName(), before);
                    continue;
                }
                final GATKRead altered = PrintDistantMates.doDistantMateAlterations(read);
                final GATKRead undone = PrintDistantMates.undoDistantMateAlterations(altered);
                final String after = undone.convertToSAMRecord(header).getSAMString().trim();
                System.out.printf("roundtrip\t%s\t%s\t%s\t%s%n", record.getReadName(),
                        before.equals(after) ? "same" : "differs", before, after);
                System.out.printf("isdistant\t%s\t%s\t%s\t%s%n", record.getReadName(),
                        PrintDistantMates.isDistantMate(read),
                        PrintDistantMates.isDistantMate(altered),
                        PrintDistantMates.isDistantMate(undone));
            }
        }
    }

    /** undoDistantMateAlterations on what it is not given a valid OA for. */
    static void undoRefusals() {
        // No OA at all: the read comes back, and it comes back as the same object rather than a
        // copy, which is observable because a caller may then mutate its input.
        final SAMFileHeader header = fixtureHeader();
        final SAMRecord bare = new SAMRecord(header);
        bare.setReadName("bare");
        bare.setReferenceName("chr1");
        bare.setAlignmentStart(100);
        bare.setCigarString("10M");
        final GATKRead bareRead = new SAMRecordToGATKReadAdapter(bare);
        final GATKRead returned = PrintDistantMates.undoDistantMateAlterations(bareRead);
        System.out.printf("undo\tnooa\t%s%n", returned == bareRead ? "same object" : "copy");

        for (final String oa : new String[] {"garbage", "chr1,100,+,10M,60", "chr1,x,+,10M,60,3;",
                "chr9,100,+,10M,60,3;"}) {
            final SAMRecord record = new SAMRecord(header);
            record.setReadName("oa");
            record.setReferenceName("chr1");
            record.setAlignmentStart(100);
            record.setCigarString("10M");
            record.setAttribute("OA", oa);
            try {
                final GATKRead undone = PrintDistantMates.undoDistantMateAlterations(
                        new SAMRecordToGATKReadAdapter(record));
                System.out.printf("undo\t%s\t%s%n", oa,
                        undone.convertToSAMRecord(header).getSAMString().trim());
            } catch (final Exception | AssertionError e) {
                System.out.printf("undo\t%s\t%s:%s%n", oa, e.getClass().getName(),
                        String.valueOf(e.getMessage()).replace('\n', ' '));
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

    /** The fixture's header, alone, for the reads built outside a file. */
    static SAMFileHeader fixtureHeader() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 3000),
                new SAMSequenceRecord("chr2", 3000))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("s1");
        header.addReadGroup(group);
        final SAMProgramRecord existing = new SAMProgramRecord("upstream");
        existing.setProgramVersion("1.0");
        header.addProgramRecord(existing);
        return header;
    }

    /**
     * A coordinate-sorted BAM holding one read for each of the four added filters to reject and
     * three for them to keep.
     *
     * The three that survive move to chr2:600, chr1:2500 and chr2:150 in that traversal order,
     * which is deliberately not coordinate order: that is what the preSorted = false writer is
     * measured against.
     */
    static void buildFixture(final File file) {
        final SAMFileHeader header = fixtureHeader();
        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, file)) {
            // r0: passes everything, mate on the other contig, and carries an NM.
            writer.addAlignment(read(header, "r0", "chr1", 100, "chr2", 600, 3, 0));
            // r1: mate ten bases away, so MATE_DISTANT rejects it at the default threshold and
            // keeps it at 5. No NM, which is the field whose spelling this tool disagrees about.
            writer.addAlignment(read(header, "r1", "chr1", 120, "chr1", 130, null, 0));
            // r2: duplicate, rejected by NOT_DUPLICATE.
            writer.addAlignment(read(header, "r2", "chr1", 160, "chr2", 200, null, 0x400));
            // r3: secondary, rejected by PRIMARY_LINE.
            writer.addAlignment(read(header, "r3", "chr1", 180, "chr2", 300, 1, 0x100));
            // r4: distant on the SAME contig, 2300 apart, and no NM.
            writer.addAlignment(read(header, "r4", "chr1", 200, "chr1", 2500, null, 0));
            // r5: reverse strand, NM = 0, and it lands before r0 once moved.
            writer.addAlignment(read(header, "r5", "chr1", 300, "chr2", 150, 0, 0x10));
            // r6: supplementary, rejected by PRIMARY_LINE.
            writer.addAlignment(read(header, "r6", "chr1", 400, "chr2", 100, null, 0x800));
            // r7: not paired at all, rejected by PAIRED. Last, and on chr2, so -L chr1 keeps it
            // away from AddOriginalAlignmentTags, which aborts a run on one.
            final SAMRecord unpaired = read(header, "r7", "chr2", 2000, null, 0, 7, 0);
            writer.addAlignment(unpaired);
        }
    }

    /** One 10M read, paired unless the mate contig is null. */
    static SAMRecord read(final SAMFileHeader header, final String name, final String contig,
                          final int start, final String mateContig, final int mateStart,
                          final Integer editDistance, final int extraFlags) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName(contig);
        record.setAlignmentStart(start);
        record.setCigarString("10M");
        final byte[] bases = new byte[10];
        Arrays.fill(bases, (byte) 'A');
        record.setReadBases(bases);
        final byte[] quals = new byte[10];
        Arrays.fill(quals, (byte) 30);
        record.setBaseQualities(quals);
        record.setMappingQuality(60);
        if (mateContig != null) {
            record.setReadPairedFlag(true);
            record.setMateReferenceName(mateContig);
            record.setMateAlignmentStart(mateStart);
        }
        if (editDistance != null) {
            record.setAttribute("NM", editDistance);
        }
        record.setAttribute("RG", "rg1");
        record.setFlags(record.getFlags() | extraFlags);
        return record;
    }

    static void distantMates(final Path dir, final Path input, final String label,
                             final String[] extra) throws Exception {
        run("PrintDistantMates", dir, input, label, extra, argv -> {
            new PrintDistantMates().instanceMain(argv);
            return null;
        });
    }

    static void addOaTags(final Path dir, final Path input, final String label,
                          final String[] extra) throws Exception {
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

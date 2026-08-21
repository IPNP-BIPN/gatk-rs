/*
 * BuildBamIndex, taken from the reference.
 *
 * The .bai of a coordinate-sorted BAM. The index itself is htsjdk's BAMIndexer, already ported and
 * proven byte-identical; what this dump is for is everything AROUND it, which is where the tool's
 * own behaviour lives.
 *
 * Six behaviours this is built to catch.
 *
 *   - THE DEFAULT OUTPUT IS RELATIVE TO THE WORKING DIRECTORY, NOT TO THE INPUT. The path is built
 *     from `inputPath.getFileName()`, so a BAM in a subdirectory with no OUTPUT given writes its
 *     index into the directory the process was started in;
 *   - AND THE EXTENSION IS REPLACED ONLY FOR A NAME ENDING .bam: anything else gets `.bai`
 *     APPENDED, so `reads.bam.copy` becomes `reads.bam.copy.bai`;
 *   - A SAM INPUT IS REFUSED BY THE READER'S TYPE, not by the extension, with a SAMException whose
 *     message is the tool's own;
 *   - A BAM SORTED BY QUERYNAME IS REFUSED, and so is one whose header says `unsorted`, by the
 *     same message: the check is on the header's SO field alone and reads nothing;
 *   - (the check reads the header's SO field and nothing else, but a file whose header lies cannot
 *     be built here: `SAMFileWriterImpl.assertPresorted` refuses to write records out of order
 *     under a coordinate header, so that claim stays unmeasured rather than guessed);
 *   - READS WITH NO ALIGNMENT START ARE COUNTED, not dropped: the unmapped reads at the end of a
 *     coordinate-sorted file become the index's no-coordinate count;
 *   - AND AN EMPTY BAM STILL PRODUCES AN INDEX, one bin-less reference per sequence in the header.
 *
 * Output:
 *
 *     fixture\t<label>\t<the input BAM, base64>
 *     index\t<label>\t<the .bai written, base64>
 *     wrote\t<label>\t<the path the index landed on, relative to the working directory>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: BuildBamIndexDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import picard.sam.BuildBamIndex;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class BuildBamIndexDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("build-bam-index-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# BuildBamIndexDump: the .bai of a coordinate-sorted BAM");

        // A coordinate-sorted BAM with reads on two contigs and two unmapped reads at the end.
        final Path sorted = dir.resolve("sorted.bam");
        buildBam(sorted, SAMFileHeader.SortOrder.coordinate, true);
        fixture("sorted", sorted);

        // The same reads with the header claiming queryname, and with it claiming unsorted.
        final Path queryname = dir.resolve("queryname.bam");
        buildBam(queryname, SAMFileHeader.SortOrder.queryname, false);
        fixture("queryname", queryname);
        final Path unsorted = dir.resolve("unsorted.bam");
        buildBam(unsorted, SAMFileHeader.SortOrder.unsorted, false);
        fixture("unsorted", unsorted);

        // A BAM with no records at all.
        final Path empty = dir.resolve("empty.bam");
        buildEmptyBam(empty);
        fixture("empty", empty);

        // A sam file, which the reader's type refuses.
        final Path sam = dir.resolve("plain.sam");
        buildSam(sam);
        fixture("plain-sam", sam);

        // The same sorted BAM under a name that does not end in .bam.
        final Path oddName = dir.resolve("sorted.bam.copy");
        Files.copy(sorted, oddName);

        run("sorted", sorted, dir.resolve("sorted-out.bai"));
        run("empty", empty, dir.resolve("empty-out.bai"));
        run("queryname", queryname, dir.resolve("queryname-out.bai"));
        run("unsorted", unsorted, dir.resolve("unsorted-out.bai"));
        run("plain-sam", sam, dir.resolve("sam-out.bai"));
        // No OUTPUT at all, which lands beside the process and not beside the input.
        run("default-output", sorted, null);
        // And a name that does not end in .bam, whose default output gains .bai rather than
        // replacing anything.
        run("default-output-odd-name", oddName, null);
    }

    static void fixture(final String label, final Path bam) throws Exception {
        System.out.printf("fixture\t%s\t%s%n", label, RecordTransformDump.base64(bam));
    }

    /** One run, with the index it wrote and where it landed. */
    static void run(final String label, final Path input, final Path output) throws Exception {
        final List<String> argv = new ArrayList<>(List.of("I=" + input));
        if (output != null) {
            argv.add("O=" + output);
        }
        try {
            final Object code = new BuildBamIndex().instanceMain(argv.toArray(new String[0]));
            if (!Integer.valueOf(0).equals(code)) {
                System.out.printf("exit\t%s\t%s%n", label, code);
                return;
            }
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        final Path landed = output != null ? output : defaultOutput(input);
        System.out.printf("wrote\t%s\t%s%n", label, landed);
        System.out.printf("index\t%s\t%s%n", label, RecordTransformDump.base64(landed));
    }

    /** Where the tool puts an index when OUTPUT is not given, which is the working directory. */
    static Path defaultOutput(final Path input) {
        final String name = input.getFileName().toString();
        if (name.endsWith(".bam")) {
            return Path.of(name.substring(0, name.lastIndexOf('.')) + ".bai");
        }
        return Path.of(name + ".bai");
    }

    /** A BAM whose header claims the given order, optionally with unmapped reads at the end. */
    static void buildBam(final Path file, final SAMFileHeader.SortOrder order,
                         final boolean withUnmapped) {
        final SAMFileHeader header = header(order);
        final int[] starts = {100, 300, 500};
        try (SAMFileWriter writer =
                     new SAMFileWriterFactory().makeBAMWriter(header, true, file.toFile())) {
            for (int i = 0; i < starts.length; i++) {
                writer.addAlignment(read(header, "r" + i, "chr1", starts[i]));
            }
            writer.addAlignment(read(header, "s0", "chr2", 200));
            if (withUnmapped) {
                for (int i = 0; i < 2; i++) {
                    final SAMRecord unmapped = new SAMRecord(header);
                    unmapped.setReadName("u" + i);
                    unmapped.setReadUnmappedFlag(true);
                    unmapped.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
                    unmapped.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
                    unmapped.setReadBases(bases());
                    unmapped.setBaseQualities(qualities());
                    writer.addAlignment(unmapped);
                }
            }
        }
    }

    static void buildEmptyBam(final Path file) {
        try (SAMFileWriter writer = new SAMFileWriterFactory()
                .makeBAMWriter(header(SAMFileHeader.SortOrder.coordinate), true, file.toFile())) {
            // Nothing to add: the header alone is the file.
            assert writer != null;
        }
    }

    static void buildSam(final Path file) {
        final SAMFileHeader header = header(SAMFileHeader.SortOrder.coordinate);
        try (SAMFileWriter writer =
                     new SAMFileWriterFactory().makeSAMWriter(header, true, file.toFile())) {
            writer.addAlignment(read(header, "r0", "chr1", 100));
        }
    }

    static SAMFileHeader header(final SAMFileHeader.SortOrder order) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 1000),
                new SAMSequenceRecord("chr2", 1000))));
        header.setSortOrder(order);
        return header;
    }

    static SAMRecord read(final SAMFileHeader header, final String name, final String contig,
                          final int start) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName(contig);
        record.setAlignmentStart(start);
        record.setCigarString("10M");
        record.setReadBases(bases());
        record.setBaseQualities(qualities());
        record.setMappingQuality(60);
        return record;
    }

    static byte[] bases() {
        final byte[] bases = new byte[10];
        Arrays.fill(bases, (byte) 'A');
        return bases;
    }

    static byte[] qualities() {
        final byte[] qualities = new byte[10];
        Arrays.fill(qualities, (byte) 30);
        return qualities;
    }
}

/*
 * CreateHadoopBamSplittingIndex, taken from the reference.
 *
 * The .sbi of a BAM: a list of the virtual offsets a splitting reader may start at, one every
 * `granularity` records. The index writer itself is htsjdk's; what this dump is for is the tool
 * around it, which is where the naming, the granularity and the refusals live.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE DEFAULT OUTPUT IS THE INPUT'S WHOLE NAME PLUS `.sbi`, `reads.bam` becoming
 *     `reads.bam.sbi` rather than `reads.sbi`, which is the opposite of what BuildBamIndex does
 *     with the same kind of argument;
 *   - THE GRANULARITY IS RECORDS AND NOT BYTES, so a granularity of two over five records leaves
 *     three offsets and a granularity above the record count leaves one;
 *   - THE LAST FIELD IS THE OFFSET THE NEXT RECORD WOULD HAVE BEEN WRITTEN AT, taken from the
 *     last record's chunk end and not from the file's length, so it points inside the last BGZF
 *     block and not past it;
 *   - AN EMPTY BAM STILL WRITES AN INDEX, whose next-start field falls back to the FILE LENGTH
 *     because there is no last record to ask;
 *   - --create-bai WRITES A SECOND FILE BESIDE THE FIRST, whose name is the index's with the
 *     extension replaced, so an index named `reads.bam.sbi` puts its companion at `reads.bam.bai`;
 *   - AND IT MAKES THE TOOL READ THE RECORDS rather than the blocks, which is why it is the only
 *     path that refuses a file that is not coordinate sorted;
 *   - WITHOUT IT A QUERYNAME-SORTED BAM IS INDEXED WITHOUT COMPLAINT, the offsets being the same
 *     ones;
 *   - A GRANULARITY OF ZERO OR LESS IS REFUSED before anything is opened;
 *   - A FILE THAT IS NOT A BAM IS REFUSED BY ITS EXTENSION, the message naming the extension it
 *     found rather than the file;
 *   - AND THE MD5 AND THE UUID THE HEADER HAS ROOM FOR ARE WRITTEN AS ZEROES, the writer having
 *     nowhere to get either from.
 *
 * Output:
 *
 *     fixture\t<label>\t<the input BAM, base64>
 *     wrote\t<label>\t<the path the index landed on, relative to the working directory>
 *     index\t<label>\t<the .sbi written, base64>
 *     fields\t<label>\t<magic>,<file length>,<md5>,<uuid>,<records>,<granularity>,<count>,<offsets...>
 *     bai\t<label>\t<the .bai written, base64, or `absent`>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CreateHadoopBamSplittingIndexDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.tools.spark.CreateHadoopBamSplittingIndex;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CreateHadoopBamSplittingIndexDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("splitting-index-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# CreateHadoopBamSplittingIndexDump: the .sbi of a BAM, and the tool "
                + "around htsjdk's writer");

        // Five records on two contigs, coordinate sorted.
        final Path sorted = dir.resolve("sorted.bam");
        buildBam(sorted, SAMFileHeader.SortOrder.coordinate);
        fixture("sorted", sorted);

        // The same records under a queryname header, which only the .bai path refuses.
        final Path queryname = dir.resolve("queryname.bam");
        buildBam(queryname, SAMFileHeader.SortOrder.queryname);
        fixture("queryname", queryname);

        // A BAM with no records at all, whose next-start field has no record to come from.
        final Path empty = dir.resolve("empty.bam");
        buildEmptyBam(empty);
        fixture("empty", empty);

        // A file that is not a BAM, refused by its extension alone.
        final Path sam = dir.resolve("plain.sam");
        buildSam(sam);

        // The granularity, which counts records.
        run("granularity-default", sorted, dir.resolve("default.sbi"), null, false);
        run("granularity-two", sorted, dir.resolve("two.sbi"), 2L, false);
        run("granularity-one", sorted, dir.resolve("one.sbi"), 1L, false);
        run("granularity-above-the-count", sorted, dir.resolve("many.sbi"), 100L, false);

        // The default output name, which appends rather than replaces.
        run("default-output", sorted, null, 2L, false);

        // The .bai companion, and the sort order only it asks about.
        run("with-bai", sorted, dir.resolve("with-bai.sbi"), 2L, true);
        run("queryname-without-bai", queryname, dir.resolve("queryname.sbi"), 2L, false);
        run("queryname-with-bai", queryname, dir.resolve("queryname-bai.sbi"), 2L, true);

        // The empty file, on both paths.
        run("empty", empty, dir.resolve("empty.sbi"), 2L, false);
        run("empty-with-bai", empty, dir.resolve("empty-bai.sbi"), 2L, true);

        // The two refusals.
        run("granularity-zero", sorted, dir.resolve("zero.sbi"), 0L, false);
        run("granularity-negative", sorted, dir.resolve("negative.sbi"), -1L, false);
        run("not-a-bam", sam, dir.resolve("sam.sbi"), 2L, false);
    }

    static void fixture(final String label, final Path bam) throws Exception {
        System.out.printf("fixture\t%s\t%s%n", label, RecordTransformDump.base64(bam));
    }

    /** One run, with the index it wrote, the fields that index holds, and its companion. */
    static void run(final String label, final Path input, final Path output, final Long granularity,
                    final boolean createBai) throws Exception {
        final List<String> argv = new ArrayList<>(List.of("--input", input.toString()));
        if (output != null) {
            argv.add("--output");
            argv.add(output.toString());
        }
        if (granularity != null) {
            argv.add("--splitting-index-granularity");
            argv.add(Long.toString(granularity));
        }
        if (createBai) {
            argv.add("--create-bai");
            argv.add("true");
        }
        try {
            new CreateHadoopBamSplittingIndex().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(cause.getMessage())));
            return;
        }
        final Path landed = output != null ? output : Path.of(input + ".sbi");
        System.out.printf("wrote\t%s\t%s%n", label, landed.getFileName());
        System.out.printf("index\t%s\t%s%n", label, RecordTransformDump.base64(landed));
        System.out.printf("fields\t%s\t%s%n", label, fields(landed));
        final Path bai = Path.of(landed.toString().replaceAll("\\.sbi$", ".bai"));
        System.out.printf("bai\t%s\t%s%n", label,
                Files.exists(bai) ? RecordTransformDump.base64(bai) : "absent");
    }

    /**
     * The .sbi's own fields, read back.
     *
     * The file is the magic `SBI\1`, then a header of the indexed file's LENGTH, an md5 and a
     * uuid of sixteen bytes each, the total record count and the granularity, then the number of
     * offsets and the offsets themselves, all little-endian. The md5 and the uuid are written as
     * zeroes: the writer has nowhere to get either from, and the reader does not check them.
     */
    static String fields(final Path index) throws Exception {
        final ByteBuffer buffer = ByteBuffer.wrap(Files.readAllBytes(index))
                .order(ByteOrder.LITTLE_ENDIAN);
        final byte[] magic = new byte[4];
        buffer.get(magic);
        final long fileLength = buffer.getLong();
        final byte[] md5 = new byte[16];
        buffer.get(md5);
        final byte[] uuid = new byte[16];
        buffer.get(uuid);
        final long records = buffer.getLong();
        final long granularity = buffer.getLong();
        final long count = buffer.getLong();
        final List<String> parts = new ArrayList<>(List.of(
                new String(magic, StandardCharsets.ISO_8859_1).replace("\1", "\\1"),
                Long.toString(fileLength),
                hex(md5),
                hex(uuid),
                Long.toString(records),
                Long.toString(granularity),
                Long.toString(count)));
        for (long i = 0; i < count; i++) {
            parts.add(Long.toString(buffer.getLong()));
        }
        return String.join(",", parts);
    }

    static String hex(final byte[] bytes) {
        final StringBuilder text = new StringBuilder();
        for (final byte value : bytes) {
            text.append(String.format("%02x", value));
        }
        return text.toString();
    }

    static void buildBam(final Path file, final SAMFileHeader.SortOrder order) {
        final SAMFileHeader header = header(order);
        try (SAMFileWriter writer =
                     new SAMFileWriterFactory().makeBAMWriter(header, true, file.toFile())) {
            for (final int start : new int[]{100, 300, 500}) {
                writer.addAlignment(read(header, "r" + start, "chr1", start));
            }
            writer.addAlignment(read(header, "s0", "chr2", 200));
            writer.addAlignment(read(header, "s1", "chr2", 400));
        }
    }

    static void buildEmptyBam(final Path file) {
        try (SAMFileWriter writer = new SAMFileWriterFactory()
                .makeBAMWriter(header(SAMFileHeader.SortOrder.coordinate), true, file.toFile())) {
            // The header alone is the file.
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
        final byte[] bases = new byte[10];
        Arrays.fill(bases, (byte) 'A');
        record.setReadBases(bases);
        final byte[] qualities = new byte[10];
        Arrays.fill(qualities, (byte) 30);
        record.setBaseQualities(qualities);
        record.setMappingQuality(60);
        return record;
    }
}

/*
 * GatherBamFiles's output, taken from the reference.
 *
 * Shards of a scattered run concatenated. The fast path copies gzip blocks and never looks at a
 * record, which is where every surprise in this dump comes from.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE PATH IS CHOSEN BY WHAT THE FILES ARE, not by what they are called: every input must be a
 *     BAM for the block copying to run, and one sam among them sends the whole run through the
 *     record-by-record gather instead;
 *   - AND THE TWO PATHS DO NOT PRODUCE THE SAME BYTES for the same records, the block copy keeping
 *     each input's own compression and the normal gather recompressing everything;
 *   - THE HEADER IS THE FIRST FILE'S AND THE OTHERS ARE DROPPED UNREAD, so a second shard carrying
 *     a read group the first does not declare is concatenated anyway and the output references a
 *     read group its header never mentions;
 *   - NOTHING CHECKS THE ORDER: shards concatenated out of coordinate order produce a file whose
 *     header says `coordinate` and whose records are not, and the tool says nothing;
 *   - AN EMPTY SHARD CONTRIBUTES NOTHING AND IS NOT AN ERROR, its header and its terminator both
 *     being dropped;
 *   - A SINGLE INPUT IS A COPY, which is not a no-op: the terminator is rewritten;
 *   - CREATE_INDEX WRITES A .bai BESIDE THE OUTPUT, and CREATE_MD5_FILE a .md5 holding the digest
 *     of the output as hex with no name and no newline;
 *   - AND AN INPUT THAT IS A LIST OF PATHS IS UNROLLED, one file per line, by `IOUtil.unrollFiles`
 *     keyed on the extension.
 *
 * Output:
 *
 *     deflater\t<class>
 *     fixture\t<label>\t<the input BAM, base64>
 *     output\t<label>\t<the gathered BAM, base64>
 *     sam\t<label>=<the gathered BAM as text, escaped>
 *     index\t<label>\t<the .bai, base64>
 *     md5\t<label>\t<the .md5's contents>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: GatherBamFilesDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import picard.sam.GatherBamFiles;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

public class GatherBamFilesDump {

    public static void main(final String[] args) throws Exception {
        // The factory is static and decides every output byte, so it is pinned before anything is
        // written and recorded beside the goldens.
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final Path dir = Path.of("gather-bam-files-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# GatherBamFilesDump: shards of a scattered run concatenated");
        System.out.printf("deflater\t%s%n",
                BlockCompressedOutputStream.getDefaultDeflaterFactory().getClass().getName());

        // Two shards in coordinate order, sharing a header.
        final Path first = dir.resolve("first.bam");
        buildBam(first, "rg1", new int[] {100, 200}, "chr1");
        fixture("first", first);
        final Path second = dir.resolve("second.bam");
        buildBam(second, "rg1", new int[] {300, 400}, "chr1");
        fixture("second", second);
        // A shard with no records at all.
        final Path empty = dir.resolve("empty.bam");
        buildBam(empty, "rg1", new int[] {}, "chr1");
        fixture("empty", empty);
        // A shard whose header declares a read group the first shard does not.
        final Path other = dir.resolve("other-rg.bam");
        buildBam(other, "rg2", new int[] {500}, "chr1");
        fixture("other-rg", other);
        // A shard whose records come before the first shard's.
        final Path earlier = dir.resolve("earlier.bam");
        buildBam(earlier, "rg1", new int[] {10, 20}, "chr1");
        fixture("earlier", earlier);
        // The same records as `second`, written as a sam file.
        final Path sam = dir.resolve("second.sam");
        buildSam(sam, "rg1", new int[] {300, 400}, "chr1");
        fixture("second-sam", sam);

        // A text file listing two of the BAMs, one per line.
        final Path list = dir.resolve("shards.bam.list");
        Files.writeString(list, first + "\n" + second + "\n", StandardCharsets.UTF_8);

        run(dir, "two-shards", List.of(first, second));
        run(dir, "single", List.of(first));
        run(dir, "with-empty", List.of(first, empty, second));
        // The second file's read group is never read, so the output declares only the first's.
        run(dir, "other-read-group", List.of(first, other));
        // Out of order, which nothing checks.
        run(dir, "out-of-order", List.of(second, earlier));
        // A sam among the inputs, which sends the whole run down the other path.
        run(dir, "with-sam", List.of(first, sam));
        // The index and the digest.
        run(dir, "indexed", List.of(first, second), "CREATE_INDEX=true");
        run(dir, "md5", List.of(first, second), "CREATE_MD5_FILE=true");
        // A list file, unrolled into the two BAMs it names.
        run(dir, "unrolled", List.of(list));
        // And a file that is not there.
        run(dir, "absent", List.of(dir.resolve("absent.bam")));
    }

    static void fixture(final String label, final Path file) throws Exception {
        System.out.printf("fixture\t%s\t%s%n", label, RecordTransformDump.base64(file));
    }

    static void run(final Path dir, final String label, final List<Path> inputs,
                    final String... extra) throws Exception {
        final Path out = dir.resolve("gathered-" + label + ".bam");
        final List<String> argv = new ArrayList<>();
        for (final Path input : inputs) {
            argv.add("I=" + input);
        }
        argv.add("O=" + out);
        argv.add("USE_JDK_DEFLATER=true");
        argv.add("USE_JDK_INFLATER=true");
        argv.addAll(Arrays.asList(extra));
        try {
            final Object code = new GatherBamFiles().instanceMain(argv.toArray(new String[0]));
            if (!Integer.valueOf(0).equals(code)) {
                System.out.printf("exit\t%s\t%s%n", label, code);
                return;
            }
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            return;
        }
        System.out.printf("output\t%s\t%s%n", label, RecordTransformDump.base64(out));
        System.out.printf("sam\t%s=%s%n", label, ReferenceQueryDump.escape(asText(out)));
        final Path index = dir.resolve("gathered-" + label + ".bai");
        if (Files.exists(index)) {
            System.out.printf("index\t%s\t%s%n", label, RecordTransformDump.base64(index));
        }
        final Path digest = dir.resolve("gathered-" + label + ".bam.md5");
        if (Files.exists(digest)) {
            System.out.printf("md5\t%s\t%s%n", label, Files.readString(digest));
        }
    }

    /** The whole file as text, header included, so a divergence reads as a line. */
    static String asText(final Path bam) {
        final StringBuilder text = new StringBuilder();
        try (SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .open(new File(bam.toString()))) {
            text.append(reader.getFileHeader().getSAMString());
            for (final SAMRecord record : reader) {
                text.append(PrintReadsDump.samLine(record));
            }
        } catch (final Exception e) {
            text.append("error: ").append(e);
        }
        return text.toString();
    }

    static SAMFileHeader header(final String readGroup) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 1000),
                new SAMSequenceRecord("chr2", 1000))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord(readGroup);
        group.setSample("s1");
        group.setLibrary("lib1");
        group.setPlatform("illumina");
        group.setPlatformUnit("unit1");
        header.addReadGroup(group);
        return header;
    }

    static void buildBam(final Path file, final String readGroup, final int[] starts,
                         final String contig) {
        final SAMFileHeader header = header(readGroup);
        try (SAMFileWriter writer =
                     new SAMFileWriterFactory().makeBAMWriter(header, true, file.toFile())) {
            for (final SAMRecord record : records(header, readGroup, starts, contig)) {
                writer.addAlignment(record);
            }
        }
    }

    static void buildSam(final Path file, final String readGroup, final int[] starts,
                         final String contig) {
        final SAMFileHeader header = header(readGroup);
        try (SAMFileWriter writer =
                     new SAMFileWriterFactory().makeSAMWriter(header, true, file.toFile())) {
            for (final SAMRecord record : records(header, readGroup, starts, contig)) {
                writer.addAlignment(record);
            }
        }
    }

    static List<SAMRecord> records(final SAMFileHeader header, final String readGroup,
                                   final int[] starts, final String contig) {
        final List<SAMRecord> records = new ArrayList<>();
        for (final int start : starts) {
            final SAMRecord record = new SAMRecord(header);
            record.setReadName("r" + start);
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
            record.setAttribute("RG", readGroup);
            records.add(record);
        }
        records.sort(Comparator.comparingInt(SAMRecord::getAlignmentStart));
        return records;
    }

    /** The dump's own directory, whose absolute path reaches the refusals. */
    static String masked(final String text, final Path dir) {
        return text.replace(dir.toAbsolutePath().toString(), "<dir>");
    }
}

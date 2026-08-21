/*
 * PrintFileDiagnostics' reports, taken from the reference.
 *
 * An index or a CRAM printed as text. Which analyzer runs is decided by the file's NAME, and the
 * BAI branch is the one with a format of its own: `TextualBAMIndexWriter`, which no other tool
 * reaches.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE ANALYZER IS CHOSEN BY THE EXTENSION and nothing else, and anything unrecognised is a
 *     RuntimeException naming the raw argument rather than the resolved path;
 *   - THE BAI BRANCH READS AN EXISTING .bai AND REWRITES IT AS TEXT, so its input is an index and
 *     never a BAM;
 *   - EVERY REFERENCE IS PRINTED, including one with no bins at all, which gets its own two-line
 *     shape AND ITS OWN SPACING: `n_bin=0` and `n_intv=0` with NO SPACE after the `=`, where every
 *     other line writes `n_bin= 4` with one;
 *   - THE BIN COUNT INCLUDES THE METADATA BIN, `n_bin` being the real bins plus one whenever
 *     metadata is present, though that bin is printed separately and out of numeric order;
 *   - EVERY BIN IS PRINTED WITH ITS RANGE SUMMARY, which is the ladder `GenomicIndexUtil` builds,
 *     and every chunk with its offsets IN HEXADECIMAL, lower case and unpadded;
 *   - THE METADATA BIN IS ALWAYS 37450 AND ALWAYS CLAIMS TWO CHUNKS, one of which is the aligned
 *     and unaligned counts rather than a file offset;
 *   - THE LINEAR INDEX PRINTS ONLY ITS NON-ZERO ENTRIES, so the numbers in the left column skip;
 *   - THE NO-COORDINATE COUNT IS THE LAST LINE, and it is printed even when it is zero;
 *   - AND THE CRAI BRANCH PRINTS A HEADER LINE AND ONE LINE PER ENTRY, which is
 *     `CRAIEntry.toString`, a different shape from the file it read.
 *
 * Output:
 *
 *     fixture\t<label>\t<the input file, base64>
 *     report\t<label>=<the whole report, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: PrintFileDiagnosticsDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import org.broadinstitute.hellbender.tools.PrintFileDiagnostics;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;

public class PrintFileDiagnosticsDump {

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final Path dir = Path.of("print-file-diagnostics-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# PrintFileDiagnosticsDump: an index printed as text");

        // A BAM over two contigs, one of which holds no reads at all, with unmapped reads at the
        // end so the no-coordinate count is not zero.
        final Path bam = dir.resolve("reads.bam");
        buildBam(bam);
        final Path bai = dir.resolve("reads.bai");
        System.out.printf("fixture\tbai\t%s%n", RecordTransformDump.base64(bai));

        // A .crai written by hand, which is a gzipped table of six numbers per line.
        final Path crai = dir.resolve("reads.cram.crai");
        Files.write(crai, gzipped("0\t100\t50\t1000\t200\t300\n0\t200\t50\t1300\t200\t300\n"
                + "1\t10\t20\t1600\t200\t300\n"));
        System.out.printf("fixture\tcrai\t%s%n", RecordTransformDump.base64(crai));

        // And a file whose extension no analyzer claims.
        final Path other = dir.resolve("notes.txt");
        Files.writeString(other, "not an index\n", StandardCharsets.UTF_8);

        run(dir, "bai", bai);
        run(dir, "crai", crai);
        run(dir, "unsupported", other);
    }

    static byte[] gzipped(final String text) throws Exception {
        final java.io.ByteArrayOutputStream bytes = new java.io.ByteArrayOutputStream();
        try (java.util.zip.GZIPOutputStream out = new java.util.zip.GZIPOutputStream(bytes)) {
            out.write(text.getBytes(StandardCharsets.UTF_8));
        }
        return bytes.toByteArray();
    }

    static void run(final Path dir, final String label, final Path input) throws Exception {
        final Path out = dir.resolve("report-" + label + ".txt");
        try {
            new PrintFileDiagnostics().instanceMain(new String[] {
                    "-I", input.toString(), "-O", out.toString()});
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            return;
        }
        if (Files.exists(out)) {
            System.out.printf("report\t%s=%s%n", label,
                    ReferenceQueryDump.escape(masked(Files.readString(out), dir)));
        }
    }

    /** A coordinate-sorted BAM with an index beside it, over two contigs. */
    static void buildBam(final Path file) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 100000),
                // A contig no read is on, so its index content is empty.
                new SAMSequenceRecord("chr2", 100000),
                new SAMSequenceRecord("chr3", 100000))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        try (SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .makeBAMWriter(header, true, file.toFile())) {
            // Reads far apart, so the linear index has holes and the bins differ.
            for (final int start : new int[] {100, 20000, 60000}) {
                writer.addAlignment(read(header, "chr1", start));
            }
            writer.addAlignment(read(header, "chr3", 500));
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

    static SAMRecord read(final SAMFileHeader header, final String contig, final int start) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(contig + "-" + start);
        record.setReferenceName(contig);
        record.setAlignmentStart(start);
        record.setCigarString("100M");
        record.setReadBases(bases());
        record.setBaseQualities(qualities());
        record.setMappingQuality(60);
        return record;
    }

    static byte[] bases() {
        final byte[] bases = new byte[100];
        Arrays.fill(bases, (byte) 'A');
        return bases;
    }

    static byte[] qualities() {
        final byte[] qualities = new byte[100];
        Arrays.fill(qualities, (byte) 30);
        return qualities;
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}

/*
 * CompareBaseQualities' output, taken from the reference.
 *
 * The third tool of the reporting-walker archetype, and the first here that reads TWO BAMs at once.
 * It walks both in lockstep, counts every pair of base qualities into a 94x94 matrix, and prints
 * that matrix twice: once as it stands and once through the static quantization mapping.
 *
 * Ten behaviours this is built to catch.
 *
 *   - IT IS NOT A WALKER AT ALL. `PicardCommandLineProgram.doWork`, two positional arguments and no
 *     read filters, so a duplicate and a vendor failure are counted like any other read. Its reader
 *     is STRICT rather than the engine's SILENT, so a record that contradicts itself is a
 *     SAMFormatException where a walker would have carried on;
 *   - BUT A SECONDARY OR SUPPLEMENTARY READ IS SKIPPED, by htsjdk's
 *     SecondaryOrSupplementarySkippingIterator, in BOTH files independently, so the two files line
 *     up on their primary reads and not on their record counts;
 *   - THE READS MUST HAVE THE SAME NAMES IN THE SAME ORDER, and a mismatch names both:
 *     "files do not have the same exact order of reads:A vs B";
 *   - A FILE THAT RUNS OUT FIRST IS "files do not have the same exact number of reads", a different
 *     message from the one above;
 *   - TWO READS OF DIFFERENT QUALITY LENGTHS ARE REFUSED BY CompareMatrix.add, with the two lengths
 *     in the message and no mention of the read's name;
 *   - THE SUMMARY COLLAPSES THE MATRIX ONTO ITS DIAGONALS. When every count is on the main
 *     diagonal it prints one sentence, "all N quality scores are the same"; otherwise it prints a
 *     table of `diff`, `count` and a percentage FORMATTED %.4f;
 *   - THE FULL MATRIX PRINTS ONLY NON-ZERO ENTRIES, in row-major order, with `diff` = QRead1 -
 *     QRead2, and the header row is always printed even when nothing follows it;
 *   - THE BINNED MATRIX IS THE SAME COUNTS THROUGH --static-quantized-quals, and with no
 *     quantization asked for the mapping is the identity, so the two halves of the output agree;
 *   - --round-down-quantized WITHOUT --static-quantized-quals IS A CommandLineException.BadArgumentValue
 *     before anything is read;
 *   - AND --throw-on-diff TURNS A DIFFERENCE INTO A UserException, where without it the tool
 *     finishes and the difference is only in the file.
 *
 * Output:
 *
 *     fixture\t<label>\t<the input BAM, base64>
 *     report\t<label>\t<the whole output file, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CompareBaseQualitiesDump
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
import org.broadinstitute.hellbender.tools.validation.CompareBaseQualities;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CompareBaseQualitiesDump {

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final Path dir = Path.of("comparequals-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# CompareBaseQualitiesDump: CompareBaseQualities' output, from the reference");

        // Two files whose qualities agree everywhere.
        fixture(dir, "same-a", new int[][] {{30, 30, 30, 30}, {20, 20, 20, 20}}, new String[] {"r1", "r2"},
                new int[] {0, 0});
        fixture(dir, "same-b", new int[][] {{30, 30, 30, 30}, {20, 20, 20, 20}}, new String[] {"r1", "r2"},
                new int[] {0, 0});
        // The same reads with some qualities moved, in both directions.
        fixture(dir, "shifted", new int[][] {{31, 30, 29, 30}, {20, 22, 20, 20}}, new String[] {"r1", "r2"},
                new int[] {0, 0});
        // A file whose second read is named differently.
        fixture(dir, "renamed", new int[][] {{30, 30, 30, 30}, {20, 20, 20, 20}}, new String[] {"r1", "other"},
                new int[] {0, 0});
        // A file with one read fewer.
        fixture(dir, "shorter", new int[][] {{30, 30, 30, 30}}, new String[] {"r1"}, new int[] {0});
        // A file whose second read has three qualities where the others have four.
        fixture(dir, "ragged", new int[][] {{30, 30, 30, 30}, {20, 20, 20}}, new String[] {"r1", "r2"},
                new int[] {0, 0});
        // A file carrying a secondary and a supplementary read between the primary ones.
        fixture(dir, "with-secondary",
                new int[][] {{30, 30, 30, 30}, {40, 40, 40, 40}, {45, 45, 45, 45}, {20, 20, 20, 20}},
                new String[] {"r1", "sec", "supp", "r2"},
                new int[] {0, 0x100, 0x800, 0});
        // Flags no filter would keep, which this tool has none of: a duplicate and a vendor
        // failure, both still mapped so the reader's own validation has nothing to say.
        fixture(dir, "flagged",
                new int[][] {{30, 30, 30, 30}, {20, 20, 20, 20}},
                new String[] {"r1", "r2"},
                new int[] {0x400, 0x200});
        // An unmapped read that kept its mapping quality, which is a validation error rather than
        // a filtered read: this tool reads under STRICT stringency, not the engine's SILENT.
        fixture(dir, "unmapped-with-mapq",
                new int[][] {{30, 30, 30, 30}, {20, 20, 20, 20}},
                new String[] {"r1", "r2"},
                new int[] {0, 0x4});

        run(dir, "identical", "same-a", "same-b", new String[] {});
        run(dir, "shifted", "same-a", "shifted", new String[] {});
        // The same difference the other way round, where the diffs change sign.
        run(dir, "shifted-reversed", "shifted", "same-a", new String[] {});
        run(dir, "renamed", "same-a", "renamed", new String[] {});
        run(dir, "shorter", "same-a", "shorter", new String[] {});
        run(dir, "ragged", "same-a", "ragged", new String[] {});
        // Secondary and supplementary reads are skipped in the file that has them.
        run(dir, "with-secondary", "same-a", "with-secondary", new String[] {});
        // Duplicates and vendor failures are counted like anything else: no filters at all.
        run(dir, "flagged", "same-a", "flagged", new String[] {});
        // And the reader is STRICT, so an inconsistent record is a SAMFormatException.
        run(dir, "strict-validation", "same-a", "unmapped-with-mapq", new String[] {});
        // Static quantization, which only changes the second half of the output.
        run(dir, "quantized", "same-a", "shifted",
                new String[] {"--static-quantized-quals", "10", "--static-quantized-quals", "20",
                        "--static-quantized-quals", "30"});
        run(dir, "quantized-round-down", "same-a", "shifted",
                new String[] {"--static-quantized-quals", "10", "--static-quantized-quals", "20",
                        "--static-quantized-quals", "30", "--round-down-quantized", "true"});
        // Rounding down with nothing to round to.
        run(dir, "round-down-alone", "same-a", "shifted",
                new String[] {"--round-down-quantized", "true"});
        // And the argument that turns a difference into a refusal.
        run(dir, "throw-on-diff", "same-a", "shifted", new String[] {"--throw-on-diff", "true"});
        run(dir, "throw-on-same", "same-a", "same-b", new String[] {"--throw-on-diff", "true"});
    }

    static void fixture(final Path dir, final String label, final int[][] qualities,
                        final String[] names, final int[] flags) throws Exception {
        final Path bam = dir.resolve(label + ".bam");
        final SAMFileHeader header = header();
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .makeBAMWriter(header, true, bam.toFile())) {
            for (int i = 0; i < qualities.length; i++) {
                writer.addAlignment(read(header, names[i], 10 + i * 10, qualities[i], flags[i]));
            }
        }
        System.out.printf("fixture\t%s\t%s%n", label, RecordTransformDump.base64(bam));
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 200))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("s1");
        header.addReadGroup(group);
        return header;
    }

    static SAMRecord read(final SAMFileHeader header, final String name, final int start,
                          final int[] qualities, final int flags) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setFlags(flags);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(qualities.length + "M");
        record.setReadBases("ACGTACGTAC".substring(0, qualities.length)
                .getBytes(StandardCharsets.UTF_8));
        final byte[] quals = new byte[qualities.length];
        for (int i = 0; i < qualities.length; i++) {
            quals[i] = (byte) qualities[i];
        }
        record.setBaseQualities(quals);
        record.setMappingQuality(60);
        record.setAttribute("RG", "rg1");
        return record;
    }

    /** One run over two files, with the whole report it wrote. */
    static void run(final Path dir, final String label, final String first, final String second,
                    final String[] extra) throws Exception {
        final Path output = dir.resolve("CompareBaseQualities." + label + ".txt");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                dir.resolve(first + ".bam").toString(), dir.resolve(second + ".bam").toString(),
                "-O", output.toString(), "--use-jdk-inflater", "true"));
        argv.addAll(Arrays.asList(extra));

        final Object result;
        try {
            result = new CompareBaseQualities().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(), e.getMessage());
            return;
        }
        // The return value is the tool's answer as an exit code: zero when the two agree.
        System.out.printf("result\t%s\t%s%n", label, result);
        if (Files.exists(output)) {
            System.out.printf("report\t%s\t%s%n", label,
                    ReferenceQueryDump.escape(Files.readString(output)));
        } else {
            System.out.printf("report\t%s\tabsent%n", label);
        }
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

/*
 * DumpTabixIndex' output, taken from the reference.
 *
 * The fifth tool of the reporting-walker archetype and the smallest of them all: it reads a `.tbi`
 * file, which is a gzipped little-endian structure, and prints it as text.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE MAGIC IS `TBI\1`, and the fourth byte is compared against the NUMBER 1 rather than a
 *     character, so a file whose fourth byte is the letter `1` is refused;
 *   - EVERY INTEGER IS LITTLE-ENDIAN and read a byte at a time, and a long is two ints of which the
 *     low one is masked to unsigned and the high one is NOT: a chunk offset above 2^63 comes out
 *     negative;
 *   - THE BIN NUMBER IS TURNED BACK INTO A RANGE by a ladder of six cases, whose boundaries are
 *     1, 8, 72, 584, 4680 and everything above, and whose units change from M to K at 585;
 *   - A BIN ABOVE 37448 IS A PSEUDOBIN, whose two chunks are read as four longs and printed as a
 *     summary line rather than as chunks;
 *   - AND THAT SUMMARY LINE PRINTS THE WRONG FIELD. `end=` is built from `chunkStart & 0xffff`
 *     where every other place uses the matching value, so the low half of the end offset is the low
 *     half of the START offset. It is in the output, so a port has to reproduce it;
 *   - THE LINEAR INDEX IS PRINTED IN 16K STEPS, counted in a variable that starts at zero and rises
 *     by 16 per entry, printed as `<n>K`;
 *   - WHATEVER FOLLOWS THE LAST CONTIG IS READ AS ONE LONG, "N unplaced reads.", and anything after
 *     that is "Unexpected data follows index.";
 *   - AND THE CONTIG NAMES ARE NUL-TERMINATED AND THEIR TOTAL LENGTH IS CHECKED, so a names block
 *     whose declared length disagrees is refused.
 *
 * Output:
 *
 *     source\t<label>\t<the file that was indexed, escaped>
 *     tbi\t<label>\t<the .tbi, base64>
 *     dump\t<label>\t<the printed text, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: DumpTabixIndexDump
 */

import org.broadinstitute.hellbender.tools.DumpTabixIndex;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.io.PushbackInputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPInputStream;

public class DumpTabixIndexDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("dumptabix-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# DumpTabixIndexDump: DumpTabixIndex' output, from the reference");

        // A small VCF, indexed as tabix: one contig, a handful of records.
        final Path vcf = dir.resolve("small.vcf");
        Files.writeString(vcf, vcfText(new int[] {100, 200, 300, 70_000, 200_000}),
                StandardCharsets.UTF_8);
        index(dir, vcf, "small");

        // One whose records are far enough apart to fill several linear-index entries and to reach
        // a second level of the bin ladder.
        final Path spread = dir.resolve("spread.vcf");
        Files.writeString(spread, vcfText(new int[] {1, 1_000_000, 20_000_000, 100_000_000}),
                StandardCharsets.UTF_8);
        index(dir, spread, "spread");

        // A bed file, whose tabix format number and columns are not a VCF's.
        final Path bed = dir.resolve("regions.bed");
        Files.writeString(bed, "chr1\t99\t120\nchr1\t199\t220\nchr1\t69999\t70020\n",
                StandardCharsets.UTF_8);
        index(dir, bed, "regions");

        for (final String label : new String[] {"small", "spread", "regions"}) {
            run(dir, label);
        }

        // A file whose magic is wrong, which is the tool's first refusal.
        final Path broken = dir.resolve("broken.tbi");
        try (final var out = new java.util.zip.GZIPOutputStream(Files.newOutputStream(broken))) {
            out.write(new byte[] {'T', 'B', 'I', '1', 0, 0, 0, 0});
        }
        runFile(dir, "wrong-magic", broken);

        // And one that is not gzipped at all.
        final Path plain = dir.resolve("plain.tbi");
        Files.write(plain, new byte[] {'T', 'B', 'I', 1, 0, 0, 0, 0});
        runFile(dir, "not-gzipped", plain);
    }

    /** A minimal VCF with one record per position given. */
    static String vcfText(final int[] positions) {
        final StringBuilder text = new StringBuilder();
        text.append("##fileformat=VCFv4.2\n");
        text.append("##contig=<ID=chr1,length=250000000>\n");
        text.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for (final int position : positions) {
            text.append("chr1\t").append(position).append("\t.\tA\tC\t.\t.\t.\n");
        }
        return text.toString();
    }

    /** Block-compress a file and index it as tabix, which is what makes a `.tbi`. */
    static void index(final Path dir, final Path file, final String label) throws Exception {
        final Path compressed = dir.resolve(file.getFileName() + ".gz");
        try (final var in = Files.newInputStream(file);
             final var out = new htsjdk.samtools.util.BlockCompressedOutputStream(
                     compressed.toFile())) {
            in.transferTo(out);
        }
        new org.broadinstitute.hellbender.tools.IndexFeatureFile()
                .instanceMain(new String[] {"-I", compressed.toString()});

        System.out.printf("source\t%s\t%s%n", label,
                ReferenceQueryDump.escape(Files.readString(file)));
        System.out.printf("tbi\t%s\t%s%n", label,
                RecordTransformDump.base64(dir.resolve(compressed.getFileName() + ".tbi")));
    }

    static void run(final Path dir, final String label) throws Exception {
        final Path tbi = Files.list(dir)
                .filter(path -> path.getFileName().toString().startsWith(label)
                        && path.getFileName().toString().endsWith(".tbi"))
                .findFirst()
                .orElseThrow(() -> new IllegalStateException("no index for " + label));
        runFile(dir, label, tbi);
    }

    /** The tool's own dumping routine, called directly so the text is the whole observable. */
    static void runFile(final Path dir, final String label, final Path tbi) {
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        try (final PrintStream writer = new PrintStream(bytes, true, StandardCharsets.UTF_8);
             final PushbackInputStream is = new PushbackInputStream(
                     new GZIPInputStream(Files.newInputStream(tbi)))) {
            DumpTabixIndex.dumpTabixIndex(is, writer);
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(), e.getMessage());
            return;
        }
        System.out.printf("dump\t%s\t%s%n", label,
                ReferenceQueryDump.escape(bytes.toString(StandardCharsets.UTF_8)));
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

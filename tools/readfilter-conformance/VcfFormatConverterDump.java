/*
 * VcfFormatConverter, taken from the reference.
 *
 * A file rewritten in the format its name asks for: vcf, vcf.gz or bcf, with the index and the
 * dictionary deciding what is refused.
 *
 * Seven behaviours this is built to catch.
 *
 *   - REQUIRE_INDEX DEFAULTS TO TRUE, so a plain vcf with no index beside it is refused before a
 *     record is read, by TRIBBLE and not by the tool: the message is the reader's;
 *   - CREATE_INDEX ALSO DEFAULTS TO TRUE, set in the constructor rather than by the parser, and it
 *     is the reason a file with no ##contig lines is refused: the dictionary is null and the
 *     PicardException names the indexing, not the file;
 *   - AND THE TWO ARE INDEPENDENT: REQUIRE_INDEX=false with CREATE_INDEX=true still refuses a file
 *     with no contigs, and CREATE_INDEX=false accepts it;
 *   - THE HEADER IS COPIED THROUGH `new VCFHeader(header)`, which keeps the sample list and the
 *     lines, and the writer emits them in ITS OWN ORDER, so a file whose header is out of order
 *     comes back sorted whatever the format;
 *   - A CONVERSION IS A REWRITE, NOT A COPY: every record goes through the decoder and the encoder,
 *     so an input's spacing, its missing fields and its unphased separators come out the writer's
 *     way;
 *   - THE OUTPUT FORMAT IS THE EXTENSION'S, `.vcf.gz` being block compressed and `.bcf` binary,
 *     and nothing in the arguments can override it;
 *   - AND A ROUND TRIP THROUGH BCF IS THE IDENTITY, measured and not assumed: vcf -> bcf -> vcf
 *     returns the same bytes as vcf -> vcf, header included, for a file of this shape. The binary
 *     format carries the header as text and hands it back unchanged.
 *
 * The BCF output is recorded as a digest and a length rather than as bytes: a BCF codec is a brick
 * of its own, and this golden is meant to pin the text paths a port can reproduce today while
 * leaving the binary ones measured for later.
 *
 * Output:
 *
 *     input\t<label>=<the whole input vcf, escaped>
 *     converted\t<label>=<the whole output vcf, escaped>
 *     sha256\t<label>=<digest>\t<length in bytes>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: VcfFormatConverterDump
 */

import htsjdk.samtools.util.BlockCompressedInputStream;
import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import picard.vcf.VcfFormatConverter;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class VcfFormatConverterDump {

    /** A header whose lines are deliberately out of the writer's order. */
    static final String HEADER =
            "##fileformat=VCFv4.2\n"
            + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
            + "##FILTER=<ID=q10,Description=\"Quality below 10\">\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            + "##ALT=<ID=DEL,Description=\"Deletion\">\n"
            + "##contig=<ID=chr1,length=1000>\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tone\ttwo\n";

    static final String RECORDS =
            "chr1\t100\t.\tA\tC\t30\t.\tDP=10\tGT\t0/1\t0/0\n"
            // A phased genotype and a no-call, which the writer prints its own way.
            + "chr1\t200\trs1\tA\tC\t.\tq10\t.\tGT\t0|1\t./.\n"
            // A record with two alternates and a symbolic one.
            + "chr1\t300\t.\tA\tC,G\t50\tPASS\tDP=20\tGT\t1/2\t0/0\n"
            + "chr1\t400\t.\tA\t<DEL>\t.\t.\t.\tGT\t0/1\t0/1\n";

    static final String HEADER_NO_CONTIGS = HEADER.replace("##contig=<ID=chr1,length=1000>\n", "");

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("vcf-format-converter-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# VcfFormatConverterDump: a vcf rewritten in the format its name asks for");

        final Path indexed = write(dir, "indexed", HEADER + RECORDS, true);
        final Path bare = write(dir, "bare", HEADER + RECORDS, false);
        final Path noContigs = write(dir, "no-contigs", HEADER_NO_CONTIGS + RECORDS, false);
        final Path empty = write(dir, "empty", HEADER, true);

        // The index is there and REQUIRE_INDEX is left at its default.
        run(dir, "vcf-to-vcf", indexed, ".vcf");
        // The same file with no index beside it, which the reader refuses.
        run(dir, "no-index", bare, ".vcf");
        // And the same run with the requirement dropped.
        run(dir, "no-index-allowed", bare, ".vcf", "REQUIRE_INDEX=false");
        // No contig lines, which the indexing refuses even with REQUIRE_INDEX off.
        run(dir, "no-contigs", noContigs, ".vcf", "REQUIRE_INDEX=false");
        // The same file with the index turned off, which is accepted.
        run(dir, "no-contigs-no-index", noContigs, ".vcf",
                "REQUIRE_INDEX=false", "CREATE_INDEX=false");
        // A file with no records.
        run(dir, "empty", empty, ".vcf");
        // Block compressed output, recorded as the text it decompresses to.
        run(dir, "to-gz", indexed, ".vcf.gz");
        // Binary output, recorded as a digest and a length.
        final Path bcf = run(dir, "to-bcf", indexed, ".bcf");
        // And the round trip back, which is where the header can move.
        if (bcf != null) {
            run(dir, "bcf-to-vcf", bcf, ".vcf", "REQUIRE_INDEX=false");
        }
    }

    /** A vcf written by hand, indexed when the run needs an index. */
    static Path write(final Path dir, final String label, final String text, final boolean index)
            throws Exception {
        final Path file = dir.resolve(label + ".vcf");
        Files.writeString(file, text, StandardCharsets.UTF_8);
        if (index) {
            new IndexFeatureFile().instanceMain(new String[] {"-I", file.toString()});
        }
        System.out.printf("input\t%s=%s%n", label, ReferenceQueryDump.escape(text));
        return file;
    }

    /** One run, returning the output when it was written. */
    static Path run(final Path dir, final String label, final Path input, final String extension,
                    final String... extra) throws Exception {
        final Path out = dir.resolve("converted-" + label + extension);
        final List<String> argv = new ArrayList<>(Arrays.asList("I=" + input, "O=" + out));
        argv.addAll(Arrays.asList(extra));
        try {
            final Object code = new VcfFormatConverter().instanceMain(argv.toArray(new String[0]));
            if (!Integer.valueOf(0).equals(code)) {
                System.out.printf("exit\t%s=%s%n", label, code);
                return null;
            }
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            return null;
        }

        final byte[] bytes = Files.readAllBytes(out);
        if (extension.equals(".bcf")) {
            System.out.printf("sha256\t%s=%s\t%d%n", label, digest(bytes), bytes.length);
            return out;
        }
        final String text = extension.endsWith(".gz") ? decompress(out) : Files.readString(out);
        System.out.printf("converted\t%s=%s%n", label, ReferenceQueryDump.escape(text));
        return out;
    }

    /** The dump's own directory, whose absolute path reaches the reader's messages. */
    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }

    static String decompress(final Path file) throws Exception {
        try (BlockCompressedInputStream in = new BlockCompressedInputStream(file.toFile())) {
            return new String(in.readAllBytes(), StandardCharsets.UTF_8);
        }
    }

    static String digest(final byte[] bytes) throws Exception {
        final StringBuilder text = new StringBuilder();
        for (final byte b : MessageDigest.getInstance("SHA-256").digest(bytes)) {
            text.append(String.format("%02x", b));
        }
        return text.toString();
    }
}

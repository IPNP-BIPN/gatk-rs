/*
 * IndexFeatureFile's indexes, taken from the reference.
 *
 * The tool that every other dump in this harness already calls to make its fixtures, measured for
 * itself at last. Which index it builds is decided by the file's NAME, three ways, and the bytes
 * differ accordingly for one and the same set of records.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - THE INDEX TYPE IS CHOSEN BY THE EXTENSION AND NOTHING ELSE: a block compressed file gets a
 *     TABIX index, a name ending `.g.vcf` gets a LINEAR index with a bin width of 128000, and
 *     anything else gets a DYNAMIC index built FOR_SEEK_TIME;
 *   - SO THE SAME RECORDS INDEX TO DIFFERENT BYTES under two names, `reads.vcf` and `reads.g.vcf`
 *     differing in more than their headers;
 *   - THE DEFAULT OUTPUT SITS BESIDE THE INPUT, `Tribble.indexPath` appending `.idx` to the whole
 *     name, so `reads.vcf` indexes to `reads.vcf.idx` and not to `reads.idx`;
 *   - AND FOR A BLOCK COMPRESSED FILE IT IS `.tbi` APPENDED THE SAME WAY;
 *   - AN EXPLICIT OUTPUT IS TAKEN AS GIVEN FOR A PLAIN FILE, whatever it is called, so an index
 *     can be written to a name ending `.tbi` that is not a tabix index at all;
 *   - BUT A BLOCK COMPRESSED INPUT REFUSES AN OUTPUT THAT DOES NOT END `.tbi`, and the message
 *     quotes the INPUT rather than the output;
 *   - THE TOOL RETURNS THE INDEX PATH as its result, which is what a caller sees;
 *   - THE INDEX EMBEDS THE SOURCE FILE'S URI, so the same records indexed from two directories
 *     differ byte for byte and an index is not portable between them;
 *   - AND IT EMBEDS THE SOURCE FILE'S MODIFICATION TIME, in milliseconds, which is why indexing
 *     the same bytes twice produces two different files. The dump ZEROES those eight bytes before
 *     it prints them, the only masking in this golden, because a golden that carried them could
 *     never be re-derived;
 *   - A FILE IN NO SUPPORTED FORMAT IS REFUSED BY THE CODEC LOOKUP as a NoSuitableCodecs, before
 *     any index is built;
 *   - AND A VCF WHOSE RECORDS ARE OUT OF ORDER IS A CouldNotIndexFile wrapping the Tribble
 *     complaint, which names the two positions.
 *
 * Output:
 *
 *     input\t<label>\t<the whole input file, base64>
 *     index\t<label>\t<the index written, base64>
 *     returned\t<label>\t<the path the tool returned>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: IndexFeatureFileDump
 */

import htsjdk.samtools.util.BlockCompressedOutputStream;
import org.broadinstitute.hellbender.tools.IndexFeatureFile;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class IndexFeatureFileDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            + "##contig=<ID=chr1,length=100000>\n"
            + "##contig=<ID=chr2,length=100000>\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tone\n";

    /** Enough records, spread widely enough, that the bin width has something to decide. */
    static String records(final boolean sorted) {
        final StringBuilder text = new StringBuilder();
        for (int i = 1; i <= 40; i++) {
            text.append("chr1\t").append(i * 1000).append("\t.\tA\tC\t.\t.\t.\tGT\t0/1\n");
        }
        for (int i = 1; i <= 10; i++) {
            text.append("chr2\t").append(i * 5000).append("\t.\tA\tC\t.\t.\t.\tGT\t0/1\n");
        }
        if (!sorted) {
            text.append("chr1\t10\t.\tA\tC\t.\t.\t.\tGT\t0/1\n");
        }
        return text.toString();
    }

    static final String BED =
            "chr1\t100\t200\tone\n"
            + "chr1\t5000\t5100\ttwo\n"
            + "chr2\t300\t400\tthree\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("index-feature-file-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# IndexFeatureFileDump: the index of a feature file");

        final String vcf = HEADER + records(true);
        final Path plain = write(dir, "reads.vcf", vcf);
        input("plain", plain);
        // The same bytes under a name the tool reads as a GVCF.
        final Path gvcf = write(dir, "reads.g.vcf", vcf);
        input("gvcf", gvcf);
        final Path bed = write(dir, "regions.bed", BED);
        input("bed", bed);
        // The same records block compressed.
        final Path compressed = writeCompressed(dir, "reads.vcf.gz", vcf);
        input("compressed", compressed);
        // One whose last record is before the first.
        final Path unsorted = write(dir, "unsorted.vcf", HEADER + records(false));
        input("unsorted", unsorted);
        // A file in no format the codecs know.
        final Path unknown = write(dir, "notes.txt", "this is not a feature file\n");
        input("unknown", unknown);

        run(dir, "plain", plain, null);
        run(dir, "gvcf", gvcf, null);
        run(dir, "bed", bed, null);
        run(dir, "compressed", compressed, null);
        // An explicit output, which a plain file takes as given even when it looks like a tabix.
        run(dir, "explicit", plain, dir.resolve("elsewhere.idx"));
        run(dir, "explicit-tbi-name", plain, dir.resolve("misleading.tbi"));
        // And which a block compressed input refuses unless it ends .tbi.
        run(dir, "compressed-wrong-extension", compressed, dir.resolve("wrong.idx"));
        run(dir, "compressed-explicit", compressed, dir.resolve("right.tbi"));
        run(dir, "unsorted", unsorted, null);
        run(dir, "unknown", unknown, null);
        run(dir, "absent", dir.resolve("absent.vcf"), null);
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path file = dir.resolve(name);
        Files.writeString(file, text, StandardCharsets.UTF_8);
        return file;
    }

    static Path writeCompressed(final Path dir, final String name, final String text)
            throws Exception {
        final Path file = dir.resolve(name);
        try (BlockCompressedOutputStream out = new BlockCompressedOutputStream(file.toFile())) {
            out.write(text.getBytes(StandardCharsets.UTF_8));
        }
        return file;
    }

    static void input(final String label, final Path file) throws Exception {
        System.out.printf("input\t%s\t%s%n", label, RecordTransformDump.base64(file));
    }

    static void run(final Path dir, final String label, final Path input, final Path output)
            throws Exception {
        final List<String> argv = new ArrayList<>(Arrays.asList("-I", input.toString()));
        if (output != null) {
            argv.addAll(Arrays.asList("-O", output.toString()));
        }
        final Object result;
        try {
            result = new IndexFeatureFile().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            return;
        }
        System.out.printf("returned\t%s\t%s%n", label,
                masked(String.valueOf(result), dir));
        final Path written = Path.of(String.valueOf(result));
        if (Files.exists(written)) {
            System.out.printf("index\t%s\t%s%n", label,
                    java.util.Base64.getEncoder().encodeToString(withoutTimestamp(Files.readAllBytes(written))));
        }
    }

    /**
     * A tribble index with its embedded modification time zeroed.
     *
     * `AbstractIndex.write` puts the source file's `lastModified` in the header, so the same
     * records indexed twice never produce the same file. The eight bytes sit after the magic, the
     * type, the version, the NUL-terminated path and the file size. A tabix index carries no such
     * field and is left alone.
     */
    static byte[] withoutTimestamp(final byte[] index) {
        // 'TIDX' little-endian, which is what a tribble index starts with.
        if (index.length < 12 || index[0] != 'T' || index[1] != 'I' || index[2] != 'D'
                || index[3] != 'X') {
            return index;
        }
        int position = 12;
        while (position < index.length && index[position] != 0) {
            position++;
        }
        position += 1 + 8;
        for (int i = 0; i < 8 && position + i < index.length; i++) {
            index[position + i] = 0;
        }
        return index;
    }

    /** The dump's own directory, whose absolute path reaches every message and every result. */
    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}

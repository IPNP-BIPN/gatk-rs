/*
 * PrintBGZFBlockInformation's report, taken from the reference.
 *
 * A bgzf file walked block by block, with its offsets, its sizes and its terminator checked. The
 * tool parses the blocks itself rather than through BlockCompressedInputStream, so what it reports
 * is the framing and not the content.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE REPORT NAMES THE FILE AND NOT THE PATH, `getFileName()` being what reaches the first
 *     line while every refusal quotes the whole path;
 *   - A PREMATURE TERMINATOR IS REPORTED TWICE, once inline when the block AFTER it is read and
 *     once in a summary at the end, and the summary joins the numbers with a comma and no space;
 *   - AND THE NUMBER REPORTED IS THE TERMINATOR'S OWN, `blockNumber - 1`, which is the block
 *     before the one being printed;
 *   - THE INLINE MESSAGE COMES BEFORE THE BLOCK THAT TRIGGERED IT, so the report reads as though
 *     the error belonged to the following block;
 *   - A FILE WITH NO TERMINATOR AT ALL EARNS THE MISSING BANNER, and so does an empty file, the
 *     check being `previousBlockInfo == null || uncompressedSize != 0`;
 *   - A FILE THAT IS ONLY A TERMINATOR IS ACCEPTED, block #1 with an uncompressed size of 0;
 *   - A REGULAR GZIP FILE IS REFUSED AT STARTUP by `IOUtil.isBlockCompressed`, with a message that
 *     names both possibilities, and a file that is not compressed at all takes the same path;
 *   - AND A TRUNCATED BLOCK IS AN IOException WRAPPED IN A UserException, so the report written so
 *     far is kept and the run still fails.
 *
 * Output:
 *
 *     file\t<label>\t<the input file, base64>
 *     report\t<label>=<the whole report, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: PrintBgzfBlockInformationDump
 */

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.BlockCompressedStreamConstants;
import org.broadinstitute.hellbender.tools.PrintBGZFBlockInformation;

import java.io.ByteArrayOutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.zip.GZIPOutputStream;

public class PrintBgzfBlockInformationDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("print-bgzf-block-information-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# PrintBgzfBlockInformationDump: a bgzf file walked block by block");

        // Two blocks of real data and a terminator, which is what a well formed file looks like.
        final byte[] whole = bgzf("first block payload\n", "second block payload\n");
        final Path plain = write(dir, "whole.gz", whole);
        fixture("whole", plain);

        // The same file with its terminator dropped.
        final byte[] noTerminator = withoutTerminator(whole);
        final Path bare = write(dir, "no-terminator.gz", noTerminator);
        fixture("no-terminator", bare);

        // Two well formed files glued together, so the first file's terminator sits in the middle.
        final ByteArrayOutputStream glued = new ByteArrayOutputStream();
        glued.write(whole);
        glued.write(bgzf("third block payload\n"));
        final Path premature = write(dir, "premature.gz", glued.toByteArray());
        fixture("premature", premature);

        // Three of them, so the summary has two numbers to join.
        final ByteArrayOutputStream twice = new ByteArrayOutputStream();
        twice.write(whole);
        twice.write(bgzf("third block payload\n"));
        twice.write(bgzf("fourth block payload\n"));
        final Path prematureTwice = write(dir, "premature-twice.gz", twice.toByteArray());
        fixture("premature-twice", prematureTwice);

        // A file that is nothing but a terminator block.
        final Path terminatorOnly = write(dir, "terminator-only.gz",
                BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK);
        fixture("terminator-only", terminatorOnly);

        // A regular gzip file, and a file that is not compressed at all.
        final Path gzip = write(dir, "regular.gz", gzipped("not a bgzf file at all\n"));
        fixture("regular-gzip", gzip);
        final Path text = write(dir, "plain.txt", "not compressed at all\n".getBytes(StandardCharsets.UTF_8));
        fixture("plain-text", text);

        // A file whose last block is cut in half.
        final byte[] cut = new byte[whole.length - 40];
        System.arraycopy(whole, 0, cut, 0, cut.length);
        final Path truncated = write(dir, "truncated.gz", cut);
        fixture("truncated", truncated);

        run(dir, "whole", plain);
        run(dir, "no-terminator", bare);
        run(dir, "premature", premature);
        run(dir, "premature-twice", prematureTwice);
        run(dir, "terminator-only", terminatorOnly);
        run(dir, "regular-gzip", gzip);
        run(dir, "plain-text", text);
        run(dir, "truncated", truncated);
        run(dir, "absent", dir.resolve("absent.gz"));
    }

    static Path write(final Path dir, final String name, final byte[] bytes) throws Exception {
        final Path file = dir.resolve(name);
        Files.write(file, bytes);
        return file;
    }

    static void fixture(final String label, final Path file) throws Exception {
        System.out.printf("file\t%s\t%s%n", label, RecordTransformDump.base64(file));
    }

    /** One payload per block, each flushed so the block boundaries are where they are asked for. */
    static byte[] bgzf(final String... payloads) throws Exception {
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        try (BlockCompressedOutputStream out = new BlockCompressedOutputStream(bytes, (Path) null)) {
            for (final String payload : payloads) {
                out.write(payload.getBytes(StandardCharsets.UTF_8));
                out.flush();
            }
        }
        return bytes.toByteArray();
    }

    static byte[] withoutTerminator(final byte[] file) {
        final int length = file.length - BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK.length;
        final byte[] shorter = new byte[length];
        System.arraycopy(file, 0, shorter, 0, length);
        return shorter;
    }

    static byte[] gzipped(final String text) throws Exception {
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        try (GZIPOutputStream out = new GZIPOutputStream(bytes)) {
            out.write(text.getBytes(StandardCharsets.UTF_8));
        }
        return bytes.toByteArray();
    }

    static void run(final Path dir, final String label, final Path input) throws Exception {
        final Path report = dir.resolve("report-" + label + ".txt");
        try {
            new PrintBGZFBlockInformation().instanceMain(new String[] {
                    "--bgzf-file", input.toString(), "-O", report.toString()});
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            // The report may exist anyway: the parsing failure comes after some blocks were
            // printed.
            if (Files.exists(report)) {
                System.out.printf("report\t%s=%s%n", label,
                        ReferenceQueryDump.escape(masked(Files.readString(report), dir)));
            }
            return;
        }
        System.out.printf("report\t%s=%s%n", label,
                ReferenceQueryDump.escape(masked(Files.readString(report), dir)));
    }

    /** The dump's own directory, whose absolute path reaches every refusal. */
    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}

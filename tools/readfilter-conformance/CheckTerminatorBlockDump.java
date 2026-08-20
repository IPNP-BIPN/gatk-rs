/*
 * CheckTerminatorBlock, and the termination check under it, taken from the reference.
 *
 * The tool is four lines over `BlockCompressedInputStream.checkTermination`, whose three answers
 * are what every BAM writer's "did this file finish" question comes down to.
 *
 * Seven behaviours this is built to catch.
 *
 *   - A FILE SHORTER THAN THE EMPTY BLOCK IS DEFECTIVE without being read at all, the length test
 *     coming before any seek;
 *   - THE TERMINATOR IS COMPARED AS 28 LITERAL BYTES, so a file whose last block is an empty gzip
 *     block written by any other deflater is NOT a terminator;
 *   - A FILE WITH NO TERMINATOR IS SEARCHED BACKWARDS from its end for a block preamble, and the
 *     search window is the whole file when the file is smaller than the maximum block size;
 *   - THE PREAMBLE MATCH IS THE FIRST ONE FOUND WALKING BACKWARDS, and the answer is decided on
 *     that one alone: a healthy block earlier in the file cannot rescue a truncated one after it;
 *   - THE SIZE FIELD IS `BSIZE`, one LESS than the block's length, so the test is
 *     `remaining == totalBlockSizeMinusOne + 1`;
 *   - A FILE THAT IS NOT GZIP AT ALL IS DEFECTIVE rather than an exception, the backwards search
 *     simply finding no preamble;
 *   - AND THE TOOL'S EXIT CODE IS 100 FOR DEFECTIVE and 0 for both other answers, so a file with a
 *     healthy last block but no terminator passes.
 *
 * Output:
 *
 *     fixture\t<label>\t<the file, base64>
 *     termination\t<label>=<HAS_TERMINATOR_BLOCK|HAS_HEALTHY_LAST_BLOCK|DEFECTIVE>
 *     exit\t<label>=<the tool's return code>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CheckTerminatorBlockDump
 */

import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.BlockCompressedStreamConstants;
import picard.sam.CheckTerminatorBlock;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

public class CheckTerminatorBlockDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("check-terminator-block-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# CheckTerminatorBlockDump: whether a block-compressed file finished");

        // A complete bgzipped file: two blocks of payload and the terminator.
        final byte[] complete = bgzip("The quick brown fox jumps over the lazy dog.\n".repeat(40));
        probe(dir, "complete", complete);
        // The same file with its terminator cut off, whose last block is still healthy.
        probe(dir, "no-terminator", Arrays.copyOf(complete,
                complete.length - BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK.length));
        // The same again with one more byte gone, so the last block is short.
        probe(dir, "truncated-block", Arrays.copyOf(complete,
                complete.length - BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK.length - 1));
        // Only the terminator.
        probe(dir, "terminator-only", BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK.clone());
        // One byte short of the terminator's length, which is the length test.
        probe(dir, "too-short", Arrays.copyOf(BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK,
                BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK.length - 1));
        // Nothing at all.
        probe(dir, "empty", new byte[0]);
        // Not gzip, and long enough to be searched.
        probe(dir, "not-gzip", "this is not a block compressed file at all, not one byte of it\n"
                .getBytes());
        // A complete file with a byte flipped inside its last block's payload, which the check
        // does not decompress and therefore does not notice.
        final byte[] flipped = complete.clone();
        flipped[flipped.length - BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK.length - 5] ^= 0xFF;
        probe(dir, "corrupt-payload", flipped);
        // A terminator whose last byte is wrong, so the 28-byte comparison fails and the backwards
        // search takes over.
        final byte[] wrongTerminator = complete.clone();
        wrongTerminator[wrongTerminator.length - 1] ^= 0x01;
        probe(dir, "wrong-terminator", wrongTerminator);
    }

    /** htsjdk's own bgzip writer, so the bytes are the ones a BAM would carry. */
    static byte[] bgzip(final String text) throws Exception {
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        try (BlockCompressedOutputStream out = new BlockCompressedOutputStream(bytes, (File) null)) {
            out.write(text.getBytes());
        }
        return bytes.toByteArray();
    }

    static void probe(final Path dir, final String label, final byte[] contents) throws Exception {
        final Path file = dir.resolve(label + ".bam");
        Files.write(file, contents);
        System.out.printf("fixture\t%s\t%s%n", label, RecordTransformDump.base64(file));
        try {
            final BlockCompressedInputStream.FileTermination termination =
                    BlockCompressedInputStream.checkTermination(file.toFile());
            System.out.printf("termination\t%s=%s%n", label, termination.name());
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
        try {
            final Object code = new CheckTerminatorBlock()
                    .instanceMain(new String[] {"I=" + file});
            System.out.printf("exit\t%s=%s%n", label, code);
        } catch (final Exception e) {
            System.out.printf("error\t%s-tool\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }
}

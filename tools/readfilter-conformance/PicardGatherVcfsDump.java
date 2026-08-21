/*
 * Picard's GatherVcfs, taken from the reference.
 *
 * Concatenating shards of one scattered run. GATK's GatherVcfsCloud is already ported and measured
 * under `gather-vcfs`; this is the OTHER tool of nearly the same name, and the two disagree about
 * how they refuse, about which order checks they run and about what a comment line is called. The
 * dump's class is named PicardGatherVcfsDump because GatherVcfsDump already holds the GATK tool,
 * and on a case-insensitive filesystem two files of one name are one file.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - A REFUSAL IS AN EXIT CODE, NOT AN EXCEPTION. `doWork` wraps everything after the dictionary
 *     check in `try { ... } catch (RuntimeException e)`, logs, DELETES THE OUTPUT and returns 1, so
 *     a caller sees a status and a log line where GATK's tool throws;
 *   - AND AN AssertionError IS NOT A RuntimeException: the dictionary mismatch, which
 *     `assertSameDictionary` raises as an AssertionError, walks straight out of doWork while the
 *     sample mismatch beside it becomes exit 1;
 *   - CREATE_INDEX DEFAULTS TO TRUE, set in the constructor and not by the parser, and an input
 *     with no ##contig lines is then a PicardException RAISED BEFORE THE TRY, so that one is an
 *     exception too;
 *   - WITH CREATE_INDEX=false AND NO CONTIG LINES the dictionary is null and the ordering check
 *     dereferences it, so the same input fails as a NullPointerException turned into exit 1;
 *   - THERE ARE TWO ORDER CHECKS AND THEY COMPARE DIFFERENT THINGS. The first compares the FIRST
 *     record of each file with the first record of the previous file; the second, inside the
 *     gathering, compares the next file's first record with the LAST RECORD WRITTEN. A pair of
 *     shards that overlap passes the first and fails the second, with a different message;
 *   - AN EMPTY SHARD IS SKIPPED BY BOTH, `lastContext` never moving, so a shard with no records
 *     between two ordered ones changes nothing;
 *   - REORDER_INPUT_BY_FIRST_VARIANT SORTS BY THE FIRST RECORD and puts every EMPTY FILE LAST,
 *     the comparator answering 1 whenever the left file is empty;
 *   - ONLY THE FIRST COMMENT SURVIVES. They are added to the first file's header as
 *     `GatherVcfs.comment` lines, a different key from the MergeVcfs.comment the neighbouring tool
 *     writes, and VCFHeader keys an unstructured line BY ITS KEY ALONE: the second CO= is dropped
 *     without a word, so two comments in and one comment out;
 *   - THE HEADER WRITTEN IS THE FIRST FILE'S, whatever the others declare;
 *   - THE MODE IS CHOSEN BY FILE EXTENSION AND NEEDS EVERY INPUT AND THE OUTPUT BLOCK COMPRESSED,
 *     so .vcf.gz inputs gathered into a .vcf take the conventional path;
 *   - AND THE BLOCK COPYING PATH REWRITES ONLY THE BLOCK THAT ENDS EACH LATER FILE'S HEADER,
 *     copying every other block byte for byte and appending one terminator block of its own.
 *
 * Output:
 *
 *     input\t<label>=<the whole input vcf, escaped>
 *     gathered\t<label>=<the whole output vcf, escaped>
 *     blocks\t<label>=<the output's bgzf block sizes, comma separated>
 *     sha256\t<label>=<the output file's digest>
 *     exit\t<label>=<the exit code>\t<the ERROR log lines, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: PicardGatherVcfsDump
 */

import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.Log;
import picard.vcf.GatherVcfs;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class PicardGatherVcfsDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            + "##contig=<ID=chr1,length=1000>\n"
            + "##contig=<ID=chr2,length=1000>\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tone\ttwo\n";

    /** The same header with the two contig lines gone. */
    static final String HEADER_NO_CONTIGS =
            HEADER.replace("##contig=<ID=chr1,length=1000>\n", "")
                    .replace("##contig=<ID=chr2,length=1000>\n", "");

    /** The same header with one contig a different length, which is a different dictionary. */
    static final String HEADER_OTHER_DICTIONARY =
            HEADER.replace("##contig=<ID=chr2,length=1000>", "##contig=<ID=chr2,length=999>");

    /** The same header with one sample instead of two. */
    static final String HEADER_ONE_SAMPLE =
            HEADER.replace("\tone\ttwo\n", "\tone\n");

    static String record(final String contig, final int position, final String genotypes) {
        return contig + "\t" + position + "\t.\tA\tC\t.\t.\t.\tGT\t" + genotypes + "\n";
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("picard-gather-vcfs-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# PicardGatherVcfsDump: shards of a scattered run concatenated");

        // Three shards in genomic order, sharing a dictionary and a sample list.
        final String first = HEADER + record("chr1", 100, "0/1\t0/0") + record("chr1", 200, "0/1\t0/0");
        final String second = HEADER + record("chr1", 300, "0/0\t1/1") + record("chr2", 50, "0/1\t0/1");
        final String third = HEADER + record("chr2", 100, "1/1\t0/0");
        // One with no records at all.
        final String empty = HEADER;
        // One whose first record is after the first shard's first record but before its last, which
        // passes the first order check and fails the second.
        final String overlapping = HEADER + record("chr1", 150, "0/1\t0/0");

        final Path a = write(dir, "a", first);
        final Path b = write(dir, "b", second);
        final Path c = write(dir, "c", third);
        final Path none = write(dir, "none", empty);
        final Path over = write(dir, "over", overlapping);
        final Path oneSample = write(dir, "one-sample",
                HEADER_ONE_SAMPLE + record("chr2", 300, "0/1"));
        final Path otherDictionary = write(dir, "other-dictionary",
                HEADER_OTHER_DICTIONARY + record("chr2", 400, "0/1\t0/0"));
        final Path bareA = write(dir, "bare-a", HEADER_NO_CONTIGS + record("chr1", 100, "0/1\t0/0"));
        final Path bareB = write(dir, "bare-b", HEADER_NO_CONTIGS + record("chr1", 300, "0/0\t1/1"));

        // The same three shards block compressed, for the two mode choices.
        final Path gzA = writeCompressed(dir, "a", first);
        final Path gzB = writeCompressed(dir, "b", second);
        final Path gzC = writeCompressed(dir, "c", third);

        run(dir, "ordered", ".vcf", List.of(a, b, c));
        run(dir, "single", ".vcf", List.of(a));
        // Two comments, which become GatherVcfs.comment lines in the first file's header.
        run(dir, "comments", ".vcf", List.of(a, b), "CO=one comment", "CO=another");
        // An empty shard between two ordered ones, which neither check notices.
        run(dir, "empty-shard", ".vcf", List.of(a, none, b));
        // The shards out of order, refused by the first check.
        run(dir, "unordered", ".vcf", List.of(c, a, b));
        // The same three reordered by the tool, with the empty one which sorts last.
        run(dir, "reordered", ".vcf", List.of(c, none, a, b), "RI=true");
        // Shards whose first records are in order and whose records overlap, refused by the second
        // check with a different message.
        run(dir, "overlapping", ".vcf", List.of(a, over), "CREATE_INDEX=false");
        // A shard with one sample where the others have two.
        run(dir, "different-samples", ".vcf", List.of(a, oneSample));
        // A shard whose dictionary differs, which is an AssertionError and not an exit code.
        run(dir, "different-dictionary", ".vcf", List.of(a, otherDictionary));
        // No contig lines, which CREATE_INDEX refuses before anything is opened.
        run(dir, "no-contigs", ".vcf", List.of(bareA, bareB));
        // The same inputs without the index, which reach the null dictionary instead.
        run(dir, "no-contigs-no-index", ".vcf", List.of(bareA, bareB), "CREATE_INDEX=false");
        // Block compressed inputs gathered into a plain vcf, which is the conventional path.
        run(dir, "gz-in-plain-out", ".vcf", List.of(gzA, gzB, gzC));
        // And block compressed on both sides, which is the block copying path.
        run(dir, "gz-in-gz-out", ".vcf.gz", List.of(gzA, gzB, gzC));
    }

    /** A plain vcf, printed as one of the golden's inputs. */
    static Path write(final Path dir, final String label, final String text) throws Exception {
        final Path file = dir.resolve(label + ".vcf");
        Files.writeString(file, text, StandardCharsets.UTF_8);
        System.out.printf("input\t%s=%s%n", label, ReferenceQueryDump.escape(text));
        return file;
    }

    /** The same text as a bgzf file, whose bytes are the reference's own writer's. */
    static Path writeCompressed(final Path dir, final String label, final String text)
            throws Exception {
        final Path file = dir.resolve(label + ".vcf.gz");
        try (BlockCompressedOutputStream out = new BlockCompressedOutputStream(file.toFile())) {
            out.write(text.getBytes(StandardCharsets.UTF_8));
        }
        return file;
    }

    static void run(final Path dir, final String label, final String extension,
                    final List<Path> inputs, final String... extra) throws Exception {
        final Path out = dir.resolve("gathered-" + label + extension);
        final List<String> argv = new ArrayList<>();
        for (final Path input : inputs) {
            argv.add("I=" + input);
        }
        argv.add("O=" + out);
        argv.addAll(Arrays.asList(extra));

        // htsjdk's Log holds System.err in a static field taken at class load, so redirecting
        // System.err does nothing: the stream has to be replaced through Log itself.
        final PrintStream realErr = Log.getGlobalPrintStream();
        final ByteArrayOutputStream captured = new ByteArrayOutputStream();
        Object code;
        try {
            Log.setGlobalPrintStream(new PrintStream(captured, true, StandardCharsets.UTF_8));
            code = new GatherVcfs().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Log.setGlobalPrintStream(realErr);
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            return;
        } finally {
            Log.setGlobalPrintStream(realErr);
        }

        if (!Integer.valueOf(0).equals(code)) {
            System.out.printf("exit\t%s=%s\t%s%n", label, code,
                    ReferenceQueryDump.escape(errorLines(captured.toString(StandardCharsets.UTF_8), dir)));
            System.out.printf("output-exists\t%s=%s%n", label, Files.exists(out));
            return;
        }

        if (extension.endsWith(".gz")) {
            final byte[] bytes = Files.readAllBytes(out);
            System.out.printf("gathered\t%s=%s%n", label,
                    ReferenceQueryDump.escape(decompress(out)));
            System.out.printf("blocks\t%s=%s%n", label, blockSizes(bytes));
            System.out.printf("sha256\t%s=%s%n", label, digest(bytes));
            return;
        }
        System.out.printf("gathered\t%s=%s%n", label,
                ReferenceQueryDump.escape(Files.readString(out)));
    }

    /** The ERROR lines of the captured log, with the timestamps and the directory masked. */
    static String errorLines(final String log, final Path dir) {
        final StringBuilder text = new StringBuilder();
        for (final String line : log.split("\n")) {
            if (!line.startsWith("ERROR")) {
                continue;
            }
            text.append(masked(line.replaceAll("\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}", "MASKED"), dir))
                    .append("\n");
        }
        return text.toString();
    }

    /** The dump's own directory, whose absolute path reaches several messages. */
    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }

    /** Everything the output decompresses to, which is what the conventional path would write. */
    static String decompress(final Path file) throws Exception {
        try (BlockCompressedInputStream in = new BlockCompressedInputStream(file.toFile())) {
            return new String(in.readAllBytes(), StandardCharsets.UTF_8);
        }
    }

    /** The BSIZE of every bgzf block, which says which blocks were copied and which rewritten. */
    static String blockSizes(final byte[] bytes) {
        final StringBuilder sizes = new StringBuilder();
        int position = 0;
        while (position + 18 <= bytes.length) {
            final int size = ((bytes[position + 16] & 0xff) | ((bytes[position + 17] & 0xff) << 8)) + 1;
            if (sizes.length() > 0) {
                sizes.append(",");
            }
            sizes.append(size);
            position += size;
        }
        return sizes.toString();
    }

    static String digest(final byte[] bytes) throws Exception {
        final StringBuilder text = new StringBuilder();
        for (final byte b : MessageDigest.getInstance("SHA-256").digest(bytes)) {
            text.append(String.format("%02x", b));
        }
        return text.toString();
    }
}

/*
 * BwaMemIndexImageCreator's image, taken from the reference.
 *
 * The tool hands a FASTA to BWA's own index builder through JNI and writes the result. The image
 * is NOT byte-comparable and this dump is what establishes that: building the same reference twice
 * in one process gives two files of the same length that differ in a handful of bytes, and those
 * bytes are in-process POINTERS. `\x70\x1a\x79\xf8\xff\x7f\x00\x00` is an address, not data.
 *
 * So the golden holds what is stable and says so: the size, where the file landed, whether two
 * builds of one reference agree, and the refusal. The bytes are not a claim this repository can
 * make, and a golden that held them would fail on the next run under a different address layout.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE DEFAULT OUTPUT IS THE INPUT'S WHOLE NAME PLUS `.img`, so `default-output.fasta` becomes
 *     `default-output.fasta.img` and not `default-output.img`;
 *   - THE IMAGE IS NOT DETERMINISTIC: two builds of one reference differ in a few bytes;
 *   - THE SIZE IS, so it is what the golden compares;
 *   - THE SIZE GROWS WITH THE SEQUENCE: five repeats of eight bases give 1333 bytes and twenty
 *     give 1551;
 *   - THE BASES ARE UPPER-CASED BEFORE INDEXING, so a reference in lower case gives an image of
 *     the same size;
 *   - THE CONTIG NAME IS IN THE IMAGE, so renaming `chr1` to `other` moves the size by one byte;
 *   - A RUN OF Ns ADDS TO THE SIZE without adding to the indexed bases, and a reference of one
 *     base still produces an image;
 *   - AND A MISSING FASTA IS A CouldNotReadReferenceException from the native side, whose message
 *     names the file and the reason.
 *
 * Output:
 *
 *     fasta\t<name>\t<that reference, escaped>
 *     wrote\t<case>\t<the file name the image landed on>
 *     size\t<case>\t<its size in bytes>
 *     stable\t<case>\t<whether two builds of this reference are byte-identical>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: BwaMemIndexImageCreatorDump
 */

import org.broadinstitute.hellbender.tools.BwaMemIndexImageCreator;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class BwaMemIndexImageCreatorDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("bwa-index-image-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# BwaMemIndexImageCreatorDump: the image BWA's own builder returns");

        final String plain = fasta(List.of("chr1"), List.of(repeat("ACGTTGCA", 5)));
        run(dir, "plain", plain, false);
        // The same reference again, to show the bytes do not move between runs.
        run(dir, "plain-again", plain, false);
        // The same bases in lower case.
        run(dir, "lower-case", fasta(List.of("chr1"), List.of(repeat("acgttgca", 5))), false);
        // A longer sequence, and a second contig.
        run(dir, "longer", fasta(List.of("chr1"), List.of(repeat("ACGTTGCA", 20))), false);
        run(dir, "two-contigs", fasta(
                List.of("chr1", "chr2"),
                List.of(repeat("ACGTTGCA", 5), repeat("TTAACCGG", 5))), false);
        // The same sequence under another contig name.
        run(dir, "renamed-contig", fasta(List.of("other"), List.of(repeat("ACGTTGCA", 5))), false);
        // Ns, which BWA holds apart from the bases it indexes.
        run(dir, "with-ns", fasta(List.of("chr1"),
                List.of(repeat("ACGTTGCA", 2) + "NNNNNNNN" + repeat("ACGTTGCA", 2))), false);
        // One base.
        run(dir, "one-base", fasta(List.of("chr1"), List.of("A")), false);
        // The default output name.
        run(dir, "default-output", plain, true);
        // A reference that is not there.
        missing(dir);
    }

    static String repeat(final String unit, final int times) {
        return unit.repeat(times);
    }

    static String fasta(final List<String> names, final List<String> sequences) {
        final StringBuilder text = new StringBuilder();
        for (int i = 0; i < names.size(); i++) {
            text.append('>').append(names.get(i)).append('\n');
            final String sequence = sequences.get(i);
            for (int at = 0; at < sequence.length(); at += 60) {
                text.append(sequence, at, Math.min(at + 60, sequence.length())).append('\n');
            }
        }
        return text.toString();
    }

    static void run(final Path dir, final String name, final String fasta,
                    final boolean defaultOutput) throws Exception {
        final Path reference = dir.resolve(name + ".fasta");
        Files.writeString(reference, fasta, StandardCharsets.UTF_8);
        System.out.printf("fasta\t%s\t%s%n", name, ReferenceQueryDump.escape(fasta));
        final Path image = dir.resolve(name + ".img");
        final List<String> argv = new ArrayList<>(List.of("-I", reference.toString()));
        if (!defaultOutput) {
            argv.add("-O");
            argv.add(image.toString());
        }
        try {
            new BwaMemIndexImageCreator().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null && cause.getCause() != cause) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", name, cause.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(cause.getMessage())
                            .replace(dir.toString(), "<dir>")));
            return;
        }
        final Path landed = defaultOutput ? Path.of(reference + ".img") : image;
        System.out.printf("wrote\t%s\t%s%n", name, landed.getFileName());
        System.out.printf("size\t%s\t%d%n", name, Files.size(landed));

        // The same reference again, into another file: the two differ in their pointer fields, so
        // the answer is `false` and the bytes stay out of the golden.
        final Path second = dir.resolve(name + "-second.img");
        new BwaMemIndexImageCreator().instanceMain(new String[]{
            "-I", reference.toString(), "-O", second.toString()});
        final byte[] first = Files.readAllBytes(landed);
        final byte[] again = Files.readAllBytes(second);
        System.out.printf("stable\t%s\t%s%n", name, java.util.Arrays.equals(first, again));
    }

    /** A reference that is not there, which fails on the native side. */
    static void missing(final Path dir) {
        try {
            new BwaMemIndexImageCreator().instanceMain(new String[]{
                "-I", dir.resolve("no-such.fasta").toString(),
                "-O", dir.resolve("no-such.img").toString()});
            System.out.printf("error\tmissing-fasta\tno failure%n");
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null && cause.getCause() != cause) {
                cause = cause.getCause();
            }
            System.out.printf("error\tmissing-fasta\t%s:%s%n", cause.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(cause.getMessage())
                            .replace(dir.toString(), "<dir>")));
        }
    }
}

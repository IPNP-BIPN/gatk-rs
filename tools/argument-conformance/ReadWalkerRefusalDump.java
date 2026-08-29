/*
 * What a read walker does with a file that is not a BAM, taken from the reference.
 *
 * `CountReads`' covering array reproduces four of its eighteen rows, and eleven of the fifteen it
 * does not are one shape: the row hands a read walker a file the corpus carries for another tool.
 * The reference refuses those, and the refusals are neither the port's nor each other's.
 *
 * Six behaviours this is built to catch.
 *
 *   - A PLAIN VCF IS NOT REFUSED BY ITS BYTES BUT BY ITS HEADER: htsjdk reads it as a SAM stream,
 *     finds no sequence dictionary in it, and the failure is an `IllegalArgumentException` from
 *     deep inside rather than a `UserException`, so the STATUS is three and not two;
 *   - A BLOCK-COMPRESSED VCF IS THE SAME REFUSAL, because the reader decompresses first;
 *   - A BED IS THE SAME AGAIN, which is what says the refusal is about the dictionary and not
 *     about the format;
 *   - AN EMPTY FILE IS A DIFFERENT ONE;
 *   - A DIRECTORY IS REFUSED BEFORE ANYTHING IS READ;
 *   - AND A PATH THAT DOES NOT EXIST IS REFUSED BY NAME, which is the one the port already
 *     reproduces.
 *
 * The status matters as much as the message: `Main.handleNonUserException` returns three, and a
 * port that answers two is claiming the reference blamed the user when it did not.
 *
 * Output:
 *
 *     status\t<case>\t<the exit status>
 *     class\t<case>\t<the exception class the run threw, or none>
 *     message\t<case>\t<its message, with the directory replaced>
 *
 * Usage: ReadWalkerRefusalDump
 */

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

public class ReadWalkerRefusalDump {

    static final Path FIXTURES = Paths.get("fixtures");

    public static void main(final String[] args) throws Exception {
        System.out.println("# ReadWalkerRefusalDump: what a read walker does with what is not a BAM");
        final Path dir = Files.createTempDirectory("walkerrefusal");

        run("a-bam", FIXTURES.resolve("reads.bam").toString());
        run("a-plain-vcf", FIXTURES.resolve("reads.vcf").toString());
        run("a-block-compressed-vcf", FIXTURES.resolve("reads.vcf.gz").toString());
        run("a-bed", FIXTURES.resolve("regions.bed").toString());

        final Path empty = dir.resolve("empty.bam");
        Files.write(empty, new byte[0]);
        run("an-empty-file", empty.toString());

        final Path text = dir.resolve("notes.txt");
        Files.writeString(text, "not a bam\n", StandardCharsets.UTF_8);
        run("a-text-file", text.toString());

        final Path directory = Files.createDirectory(dir.resolve("adirectory"));
        run("a-directory", directory.toString());

        run("a-path-that-does-not-exist", dir.resolve("nowhere.bam").toString());

        // The same inputs with an INTERVAL, which is the shape the covering array's rows have: an
        // interval needs a sequence dictionary, and the dictionary is asked for before the records
        // are read, so the refusal is a different one from the same file.
        run("a-plain-vcf-with-an-interval", FIXTURES.resolve("reads.vcf").toString(),
                "-L", "chr1:1-1000");
        run("an-empty-file-with-an-interval", empty.toString(), "-L", "chr1:1-1000");
        run("a-bam-with-an-interval", FIXTURES.resolve("reads.bam").toString(),
                "-L", "chr1:1-1000");
    }

    static void run(final String name, final String input, final String... extra) {
        final PrintStream original = System.out;
        String exception = "none";
        String message = "";
        boolean user = false;
        Object returned = null;
        try {
            System.setOut(new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
            final List<String> argv = new ArrayList<>(List.of("CountReads", "--input", input));
            argv.addAll(List.of(extra));
            returned = new org.broadinstitute.hellbender.Main()
                    .instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError thrown) {
            exception = thrown.getClass().getName();
            message = String.valueOf(thrown.getMessage());
            user = thrown instanceof org.broadinstitute.hellbender.exceptions.UserException;
        } finally {
            System.setOut(original);
        }
        // The STATUS is not measured here: `main-entry` already measured that a UserException is
        // two and anything else is three, and this says which of the two a refusal is. What is
        // measured is the class, the wording, and the value a run that did NOT refuse returned.
        System.out.printf("class\t%s\t%s%n", name, exception);
        System.out.printf("kind\t%s\t%s%n", name, exception.equals("none") ? "returned"
                : user ? "user" : "other");
        if (!exception.equals("none")) {
            System.out.printf("message\t%s\t%s%n", name, canonical(message));
        } else {
            System.out.printf("returned\t%s\t%s%n", name, String.valueOf(returned));
        }
    }

    /** A message with the run's own directory taken out of it, and its newlines escaped. */
    static String canonical(final String message) {
        return message
                .replaceAll("/tmp/walkerrefusal[0-9]+", "<dir>")
                .replace("\\", "\\\\")
                .replace("\t", "\\t")
                .replace("\n", "\\n");
    }
}

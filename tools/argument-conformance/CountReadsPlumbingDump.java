/*
 * What `CountReads` writes and returns, taken from the reference.
 *
 * The count itself is measured in `count-reads-and-bases`, which hands the ported function a
 * corpus and compares the number. What is NOT measured there is the layer between that number and
 * a command line: what the tool RETURNS, which `handleResult` prints, and what it WRITES when
 * `-O` names a file. A port with the number and without those two runs no command line.
 *
 * Six behaviours this is built to catch.
 *
 *   - THE TOOL RETURNS THE COUNT, so `Tool returned:` is followed by a number and not by a path;
 *   - THE FILE HOLDS THE NUMBER AND NOTHING ELSE, with no trailing newline: the reference writes
 *     it with `print`, and a port that used `println` would write a file one byte longer;
 *   - `-O` DOES NOT SUPPRESS THE RETURN: the number is written AND returned;
 *   - `-L` COUNTS WHAT THE INTERVAL HOLDS, which is the walker's traversal rather than the tool's
 *     own filter;
 *   - `--read-filter` ADDS TO THE DEFAULTS rather than replacing them, and
 *     `--disable-tool-default-read-filters` is what replaces them;
 *   - AND A MISSING INPUT IS A REFUSAL BEFORE ANY OF IT, whose class and wording are the
 *     reference's.
 *
 * The input is the shared fixture corpus, so a covering array row and a case here read the same
 * bytes: eight reads on one contig, seven hundred bases apart, one of them flagged a duplicate.
 *
 * Output:
 *
 *     returned\t<case>\t<the object the tool returned, as String.valueOf writes it>
 *     file\t<case>\t<the output file's bytes, base64, or absent where none was written>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: CountReadsPlumbingDump
 */

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;

public class CountReadsPlumbingDump {

    static final Path FIXTURES = Paths.get("fixtures");

    public static void main(final String[] args) throws Exception {
        System.out.println("# CountReadsPlumbingDump: what the tool returns and what it writes");

        run("the-whole-file", true);
        run("no-output-file", false);
        run("an-interval", true, "-L", "chr1:1-1000");
        run("an-interval-with-nothing-in-it", true, "-L", "chr1:50000-60000");
        run("a-second-filter", true, "--read-filter", "NotDuplicateReadFilter");
        run("the-defaults-disabled", true, "--disable-tool-default-read-filters");
        run("a-mapping-quality-filter", true,
                "--read-filter", "MappingQualityReadFilter", "--minimum-mapping-quality", "70");
        run("an-absent-input", true, "--input", "/work/nowhere.bam");
        run("no-input-at-all", true);
    }

    static void run(final String name, final boolean withOutput, final String... extra)
            throws Exception {
        final Path dir = Files.createTempDirectory("countreads");
        final Path out = dir.resolve("count.txt");
        final List<String> argv = new ArrayList<>();
        final List<String> tail = Arrays.asList(extra);
        // `an-absent-input` names its own input and `no-input-at-all` names none, so the fixture is
        // added only where neither did.
        if (!tail.contains("--input") && !name.equals("no-input-at-all")) {
            argv.addAll(List.of("--input", FIXTURES.resolve("reads.bam").toString()));
        }
        if (withOutput) {
            argv.addAll(List.of("--output", out.toString()));
        }
        argv.addAll(tail);

        // The tool's own logging goes to a sink: it carries a clock, and what is measured is the
        // value and the file rather than the noise around them.
        final PrintStream original = System.out;
        String returned = null;
        String error = null;
        try {
            System.setOut(new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
            // `instanceMain` returns the tool's own result object, which is what `handleResult`
            // prints under `Tool returned:` and what a port has to produce.
            returned = String.valueOf(new org.broadinstitute.hellbender.Main()
                    .instanceMain(prepend(argv)));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null && cause.getCause() != cause) {
                cause = cause.getCause();
            }
            error = cause.getClass().getName() + ":"
                    + String.valueOf(cause.getMessage()).replace(dir.toString(), "<dir>");
        } finally {
            System.setOut(original);
        }

        if (returned != null) {
            emit("returned", name, returned);
        }
        if (error != null) {
            emit("error", name, error);
        }
        if (withOutput && Files.exists(out)) {
            emit("file", name, Base64.getEncoder().encodeToString(Files.readAllBytes(out)));
        }
    }

    /** `Main`'s own argv: the tool's name first. */
    static String[] prepend(final List<String> argv) {
        final List<String> whole = new ArrayList<>(List.of("CountReads"));
        whole.addAll(argv);
        return whole.toArray(new String[0]);
    }

    static void emit(final String kind, final String name, final String payload) {
        System.out.printf("%s\t%s\t%s%n", kind, name,
                payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"));
    }
}

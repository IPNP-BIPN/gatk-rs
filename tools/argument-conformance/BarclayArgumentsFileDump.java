/*
 * Barclay's `--arguments_file`, taken from the reference.
 *
 * The only argument that changes the command line rather than a field. When
 * `SpecialArgumentsCollection.ARGUMENTS_FILE` is present, `parseArguments` reads the named files,
 * splits every non-comment non-blank line on whitespace, and **calls itself again** on the result.
 * Four things about that are observable and none is obvious:
 *
 *   - THE FILE'S ARGUMENTS COME FIRST. `newArgs.addAll(Arrays.asList(args))` appends the original
 *     command line to the expansion, not the other way round. So a collection ends up with the
 *     file's values before the user's, and a scalar given in both is "cannot be specified more
 *     than once" rather than "the command line wins";
 *   - RECURSION IS BOUNDED BY A SET, NOT A DEPTH. `argumentsFilesLoadedAlready` records every file
 *     named in a pass, including ones that were skipped, so a file that includes itself is read
 *     once and a file named twice on one command line is read once;
 *   - THE SECOND PASS STILL SEES `--arguments_file`. It is not removed from the command line. The
 *     recursion terminates because the expansion is empty the second time, not because the
 *     argument is gone;
 *   - THE TAG SURROGATES ARE RESET BETWEEN PASSES, because the first pass built keys for tagged
 *     options that the second pass will build again. Without the reset the second pass would see
 *     its own keys as duplicates.
 *
 * Splitting is `StringUtils.split(line)`, which splits on runs of whitespace, so alignment inside
 * a file is free and an argument value containing a space cannot be written there at all.
 *
 * Output:
 *
 *     case\t<label>\t<argv, space separated>
 *     result\t<label>\tok|E:<exception class>:<message>
 *     field\t<label>\t<long name>\t<value>
 *
 * Usage: BarclayArgumentsFileDump
 */

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.SpecialArgumentsCollection;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class BarclayArgumentsFileDump {

    public static final class Args {
        /**
         * The collection that carries `--arguments_file`. The parser does not add it: a tool that
         * wants argument files has to declare it, which GATK's CommandLineProgram does.
         */
        @ArgumentCollection
        public SpecialArgumentsCollection special = new SpecialArgumentsCollection();

        @Argument(fullName = "plain-scalar", optional = true, doc = "a scalar")
        public String plainScalar;

        @Argument(fullName = "plain-collection", optional = true, doc = "a collection")
        public List<String> plainCollection = new ArrayList<>();

        @Argument(fullName = "flag", optional = true, doc = "a flag")
        public boolean flag = false;
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("fixtures");
        if (Files.isDirectory(dir)) {
            try (final java.util.stream.Stream<Path> entries = Files.list(dir)) {
                for (final Path entry : entries.collect(java.util.stream.Collectors.toList())) {
                    Files.delete(entry);
                }
            }
        }
        Files.createDirectories(dir);

        write(dir, "scalar.txt", "--plain-scalar fromfile\n");
        write(dir, "collection.txt", "--plain-collection a --plain-collection b\n");
        // Comments, blank lines, and runs of whitespace that the split collapses.
        write(dir, "messy.txt", "# a comment\n\n   --plain-collection    one\t\t--plain-collection two   \n\n");
        write(dir, "flag.txt", "--flag\n");
        // One file that names another.
        write(dir, "outer.txt", "--arguments_file fixtures/inner.txt\n--plain-collection outer\n");
        write(dir, "inner.txt", "--plain-collection inner\n");
        // A file that names itself.
        write(dir, "self.txt", "--arguments_file fixtures/self.txt\n--plain-collection self\n");
        // Two files that name each other.
        write(dir, "ping.txt", "--arguments_file fixtures/pong.txt\n--plain-collection ping\n");
        write(dir, "pong.txt", "--arguments_file fixtures/ping.txt\n--plain-collection pong\n");

        System.out.println("# BarclayArgumentsFileDump: --arguments_file");

        run("no-file", "--plain-scalar", "direct");
        run("scalar-from-file", "--arguments_file", "fixtures/scalar.txt");
        run("collection-from-file", "--arguments_file", "fixtures/collection.txt");
        run("messy-file", "--arguments_file", "fixtures/messy.txt");
        run("flag-from-file", "--arguments_file", "fixtures/flag.txt");

        // The ordering question: the file's values come before the command line's.
        run("file-then-command-line",
                "--arguments_file", "fixtures/collection.txt", "--plain-collection", "cli");
        run("command-line-then-file",
                "--plain-collection", "cli", "--arguments_file", "fixtures/collection.txt");
        // A scalar in both is a duplicate, not an override.
        run("scalar-in-both",
                "--arguments_file", "fixtures/scalar.txt", "--plain-scalar", "cli");

        // Recursion and its bound.
        run("nested-file", "--arguments_file", "fixtures/outer.txt");
        run("self-referencing-file", "--arguments_file", "fixtures/self.txt");
        run("mutually-referencing-files", "--arguments_file", "fixtures/ping.txt");
        run("same-file-twice",
                "--arguments_file", "fixtures/collection.txt",
                "--arguments_file", "fixtures/collection.txt");
        run("two-files",
                "--arguments_file", "fixtures/scalar.txt",
                "--arguments_file", "fixtures/collection.txt");

        run("missing-file", "--arguments_file", "fixtures/nope.txt");
    }

    static void write(final Path dir, final String name, final String body) throws Exception {
        Files.write(dir.resolve(name), body.getBytes());
    }

    static void run(final String label, final String... argv) {
        final Args target = new Args();
        System.out.printf("case\t%s\t%s%n", label, String.join(" ", argv));

        String result;
        try {
            final PrintStream sink = new PrintStream(new ByteArrayOutputStream());
            final boolean parsed = new CommandLineArgumentParser(target).parseArguments(sink, argv);
            result = parsed ? "ok" : "not-parsed";
        } catch (final Exception | AssertionError e) {
            result = "E:" + e.getClass().getName() + ":"
                    + String.valueOf(e.getMessage()).replace("\n", "\\n");
        }
        System.out.printf("result\t%s\t%s%n", label, result);

        if (result.equals("ok")) {
            System.out.printf("field\t%s\t%s\t%s%n", label, "plain-scalar", target.plainScalar);
            System.out.printf("field\t%s\t%s\t%s%n", label, "plain-collection", target.plainCollection);
            System.out.printf("field\t%s\t%s\t%s%n", label, "flag", target.flag);
            System.out.printf("field\t%s\t%s\t%s%n", label, "arguments_file",
                    target.special.ARGUMENTS_FILE);
        }
    }
}

/*
 * `handleNonUserException`, which is the handler `handleUserException` is not.
 *
 * `main-entry` measured `handleUserException` and left the other one alone, so a port with one
 * banner prints it for both: a failure the reference reports as
 *
 *     java.lang.IllegalArgumentException: Dictionary cannot have size zero
 *
 * comes back from the port as
 *
 *     A USER ERROR has occurred: Dictionary cannot have size zero
 *
 * with the right message, the right status and the wrong wrapper. That is one of the two shapes of
 * divergence a CountReads covering-array row still shows (gatk-rs#1020).
 *
 * Four behaviours this is built to catch.
 *
 *   - THE HANDLER IS printStackTrace AND NOTHING ELSE: no banner of asterisks, no
 *     `A USER ERROR has occurred:` prefix, and no notice about a system property. Whatever else a
 *     port prints on this path is a line the reference never writes;
 *   - THE FIRST LINE IS Throwable.toString(), which is the class's BINARY name and, where there is
 *     one, `: ` and the message. A nested class therefore prints with a `$` and not with the dot
 *     its source spells it with;
 *   - A THROWABLE WITH NO MESSAGE IS THE CLASS ALONE, with no trailing colon at all;
 *   - AND THE Error OVERLOAD IS THE SAME HANDLER, so an OutOfMemoryError prints exactly like an
 *     exception and differs only in the status `mainEntry` exits with.
 *
 * The frames after the first line name this harness and its line numbers, so committing them would
 * make the golden fail on any edit above them. What is emitted instead is the first line, which is
 * the whole of what a row-by-row comparison reads, and whether everything after it is a frame,
 * which is the claim that nothing else was printed.
 *
 * Output:
 *
 *     non-user\t<case>\t<the first line, escaped>
 *     shape\t<case>\tframes-only=<true|false>
 *
 * Usage: MainNonUserDump
 */

import org.broadinstitute.hellbender.Main;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;

public class MainNonUserDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** A nested class, so the golden says which of the two names `Throwable.toString` uses. */
    static class Nested extends RuntimeException {
        Nested(final String message) {
            super(message);
        }
    }

    /** The handler is protected, and a subclass would not be `Main`, so it is reached by reflection. */
    static void invoke(final Class<?> type, final Throwable thrown) throws Exception {
        final var method = Main.class.getDeclaredMethod("handleNonUserException", type);
        method.setAccessible(true);
        method.invoke(new Main(), thrown);
    }

    /** The clock log4j puts in front of a line, which is the only thing here that moves. */
    static String masked(final String text) {
        final StringBuilder out = new StringBuilder();
        for (final String line : text.split("\n", -1)) {
            out.append(line.replaceFirst("^\\d\\d:\\d\\d:\\d\\d\\.\\d\\d\\d ", "")).append('\n');
        }
        return out.substring(0, Math.max(0, out.length() - 1));
    }

    static void nonUser(final String name, final Throwable thrown) throws Exception {
        final ByteArrayOutputStream errBytes = new ByteArrayOutputStream();
        final ByteArrayOutputStream outBytes = new ByteArrayOutputStream();
        final PrintStream realErr = System.err;
        final PrintStream realOut = System.out;
        try {
            System.setErr(new PrintStream(errBytes, true, StandardCharsets.UTF_8));
            // Stdout is captured as well, because "the handler writes to stderr" is part of what
            // is claimed: a port that printed the same line to the other stream would agree with
            // a golden that only looked at one of them.
            System.setOut(new PrintStream(outBytes, true, StandardCharsets.UTF_8));
            invoke(thrown instanceof Error ? Error.class : Exception.class, thrown);
        } finally {
            System.setErr(realErr);
            System.setOut(realOut);
        }
        final String[] lines = masked(errBytes.toString(StandardCharsets.UTF_8)).split("\n", -1);
        boolean framesOnly = lines.length > 1;
        for (int i = 1; i < lines.length; i++) {
            if (!lines[i].isEmpty() && !lines[i].startsWith("\tat ")) {
                framesOnly = false;
            }
        }
        emit("non-user", name, lines[0]);
        emit(
                "shape",
                name,
                "frames-only=" + framesOnly + ",stdout=" + outBytes.size());
    }

    public static void main(final String[] args) throws Exception {
        // The row that started this: an interval against a stream with no sequence dictionary.
        nonUser("illegal-argument", new IllegalArgumentException("Dictionary cannot have size zero"));
        // The other refusal a read walker reaches, whose class is htsjdk's rather than the JDK's.
        nonUser("sam-format", new htsjdk.samtools.SAMFormatException(
                "Error parsing text SAM file. Not enough fields; Line 1"));
        // No message: the class alone, and no trailing colon.
        nonUser("no-message", new IllegalStateException());
        // An empty message, which is not the same thing: a colon and nothing after it.
        nonUser("empty-message", new IllegalStateException(""));
        // A nested class, whose binary name carries a `$` where its source name carries a dot.
        nonUser("nested-class", new Nested("thrown from a class inside another"));
        // A GATK exception that is itself nested, so the `$` is measured on the reference's own
        // classes and not only on this harness's.
        nonUser("gatk-nested", new org.broadinstitute.hellbender.exceptions.GATKException
                .ShouldNeverReachHereException("a branch that was supposed to be unreachable"));
        // A message that spans lines, so the golden says the first line is the class and the head
        // of the message rather than the whole of it.
        nonUser("multi-line-message", new IllegalArgumentException("first line\nsecond line"));
        // The Error overload: the same handler, and the status is what differs.
        nonUser("out-of-memory", new OutOfMemoryError("Java heap space"));

        System.out.print(buf);
    }
}

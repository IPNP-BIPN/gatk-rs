/*
 * Main's entry: which stream each path writes to, what it returns, and what exit status an
 * exception carries. Taken from the reference.
 *
 * `mainEntry` cannot be called here: it ends in System.exit, which would take the dump's own JVM
 * with it. What is called is `instanceMain`, which is the whole of the work, plus the exit
 * constants read off the class and the two handlers `mainEntry` would have called. The mapping
 * from an exception to a status is therefore measured as the pair it is, rather than by watching
 * a process die.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - NO ARGUMENTS IS NOT AN ERROR: the usage goes to STDOUT and the run returns null, which is
 *     the one path that reaches `runCommandLineProgram` with no program;
 *   - `-h` AND `--help` ARE THE SAME PATH, and only as the FIRST argument: a tool name followed
 *     by `-h` resolves the tool instead;
 *   - AN UNKNOWN NAME PRINTS THE SAME USAGE TO STDERR and then throws, so the refusal is the
 *     usage and the message both;
 *   - `--version` IS SCANNED OVER EVERY ARGUMENT and not only the first, so a valid tool name
 *     followed by it prints the version and never runs;
 *   - AND ITS FIRST LINE GOES TO System.out WHATEVER STREAM IT IS HANDED, the two version lines
 *     under it going to the stream: printing it to stderr splits it in two;
 *   - THE FIVE EXIT STATUSES ARE 1, 2, 3, 4 AND 137, a command-line exception being 1 and a user
 *     exception 2, so a name that does not resolve exits differently from a name that does but is
 *     given nothing to work with;
 *   - A DEPRECATED TOOL IS A USER EXCEPTION LIKE ANY OTHER unknown name, carrying the notice
 *     rather than a suggestion;
 *   - A GATK TOOL MISSING ITS REQUIRED ARGUMENTS IS A CommandLineException, which is status 1 and
 *     not 2;
 *   - A PICARD TOOL IS WRAPPED, and its non-zero return arrives as a PicardNonZeroExitException
 *     carrying the code, which is status 4;
 *   - `handleResult` PRINTS NOTHING FOR NULL and `Tool returned:` and a line break for anything
 *     else, which is why a tool that returns nothing is silent;
 *   - AND `handleUserException` DECORATES THE MESSAGE with a banner of A characters and tells the
 *     reader which property prints the stack trace.
 *
 * Output:
 *
 *     status\t<name>\t<the exit constant>
 *     out\t<case>\t<what reached stdout, escaped>
 *     err\t<case>\t<what reached stderr, escaped>
 *     shape\t<case>\t<stream>=<line count>,<first line>
 *     result\t<case>\t<the returned object, or `null`>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: MainEntryDump
 */

import org.broadinstitute.hellbender.Main;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.util.List;

public class MainEntryDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /**
     * Main's handlers and its version printer are protected, and the exit statuses they pair with
     * are private, so both are reached by reflection.
     *
     * A subclass would be the obvious way in, and it is the wrong one: `printVersionInfo` asks
     * `this.getClass()` for its manifest, so a subclass loaded from anywhere but the reference's
     * own jar reports no version at all. What is measured has to be Main itself.
     */
    static void invoke(final String name, final Class<?> type, final Object argument)
            throws Exception {
        final var method = Main.class.getDeclaredMethod(name, type);
        method.setAccessible(true);
        method.invoke(new Main(), argument);
    }

    static int constant(final String name) throws Exception {
        final Field field = Main.class.getDeclaredField(name);
        field.setAccessible(true);
        return field.getInt(null);
    }

    /** The clock log4j puts in front of a line, which is the only thing here that moves. */
    static String masked(final String text) {
        final StringBuilder out = new StringBuilder();
        for (final String line : text.split("\n", -1)) {
            out.append(line.replaceFirst("^\\d\\d:\\d\\d:\\d\\d\\.\\d\\d\\d ", "")).append('\n');
        }
        return out.substring(0, Math.max(0, out.length() - 1));
    }

    /** The usage is three hundred tools long, so its shape is what the golden holds. */
    static String shape(final String stream, final String text) {
        if (text.isEmpty()) {
            return stream + "=0,";
        }
        final String[] lines = text.split("\n", -1);
        return stream + "=" + lines.length + "," + lines[0];
    }

    static void run(final String name, final List<String> argv) {
        final ByteArrayOutputStream outBytes = new ByteArrayOutputStream();
        final ByteArrayOutputStream errBytes = new ByteArrayOutputStream();
        final PrintStream realOut = System.out;
        final PrintStream realErr = System.err;
        Object result = null;
        String failure = null;
        try {
            System.setOut(new PrintStream(outBytes, true, StandardCharsets.UTF_8));
            System.setErr(new PrintStream(errBytes, true, StandardCharsets.UTF_8));
            try {
                result = new Main().instanceMain(argv.toArray(new String[0]));
            } catch (final Exception e) {
                Throwable cause = e;
                while (cause.getCause() != null && cause.getCause() != cause) {
                    cause = cause.getCause();
                }
                failure = e.getClass().getName() + ":" + String.valueOf(cause.getMessage());
                // A Picard tool's non-zero code is the exception's own field and not its message,
                // which is why the message is null on that path.
                if (e instanceof org.broadinstitute.hellbender.exceptions.PicardNonZeroExitException) {
                    failure += " code="
                            + ((org.broadinstitute.hellbender.exceptions.PicardNonZeroExitException) e)
                                    .getToolReturnCode();
                }
            }
        } finally {
            System.out.flush();
            System.err.flush();
            System.setOut(realOut);
            System.setErr(realErr);
        }
        final String out = masked(outBytes.toString(StandardCharsets.UTF_8));
        final String err = masked(errBytes.toString(StandardCharsets.UTF_8));
        emit("shape", name, shape("out", out) + " " + shape("err", err));
        if (failure != null) {
            emit("error", name, failure);
        } else {
            emit("result", name, result == null ? "null" : result.toString());
        }
    }

    public static void main(final String[] args) throws Exception {
        for (final String name : List.of("COMMANDLINE_EXCEPTION_EXIT_VALUE", "USER_EXCEPTION_EXIT_VALUE",
                "PICARD_TOOL_EXCEPTION", "ANY_OTHER_EXCEPTION_EXIT_VALUE", "OUT_OF_MEMORY_EXIT_VALUE")) {
            emit("status", name, Integer.toString(constant(name)));
        }

        // The three ways of asking for the usage, and the one way of not asking for it.
        run("no-arguments", List.of());
        run("dash-h", List.of("-h"));
        run("long-help", List.of("--help"));
        run("help-after-a-tool", List.of("CountReads", "-h"));

        // The version, which is scanned over every argument.
        run("version", List.of("--version"));
        run("version-short", List.of("-version"));
        run("version-after-a-tool", List.of("CountReads", "--version"));

        // A name that resolves to nothing, and one that resolves to a tool that went out.
        run("unknown-name", List.of("NoSuchToolAtAll"));
        run("deprecated-name", List.of("IndelRealigner"));

        // A tool given nothing to work with, on each side of the dispatch.
        run("gatk-tool-no-arguments", List.of("CountReads"));
        run("picard-tool-no-arguments", List.of("MarkDuplicates"));

        // The two handlers mainEntry would have called, on their own.
        final ByteArrayOutputStream outBytes = new ByteArrayOutputStream();
        final ByteArrayOutputStream errBytes = new ByteArrayOutputStream();
        final PrintStream realOut = System.out;
        final PrintStream realErr = System.err;
        try {
            System.setOut(new PrintStream(outBytes, true, StandardCharsets.UTF_8));
            System.setErr(new PrintStream(errBytes, true, StandardCharsets.UTF_8));
            invoke("handleResult", Object.class, null);
            invoke("handleResult", Object.class, 0);
            invoke("handleResult", Object.class, "a string");
            invoke("handleUserException", Exception.class,
                    new org.broadinstitute.hellbender.exceptions.UserException(
                            "the message a refusal carries"));
        } finally {
            System.setOut(realOut);
            System.setErr(realErr);
        }
        emit("out", "handlers", masked(outBytes.toString(StandardCharsets.UTF_8)));
        emit("err", "handlers", masked(errBytes.toString(StandardCharsets.UTF_8)));

        // The version printed to a stream that is not System.out, which splits it in two.
        final ByteArrayOutputStream versionOut = new ByteArrayOutputStream();
        final ByteArrayOutputStream versionElsewhere = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(versionOut, true, StandardCharsets.UTF_8));
            invoke("printVersionInfo", PrintStream.class,
                    new PrintStream(versionElsewhere, true, StandardCharsets.UTF_8));
        } finally {
            System.setOut(realOut);
        }
        emit("out", "version-to-stdout", masked(versionOut.toString(StandardCharsets.UTF_8)));
        emit("out", "version-to-elsewhere", masked(versionElsewhere.toString(StandardCharsets.UTF_8)));

        System.out.print(buf);
    }
}

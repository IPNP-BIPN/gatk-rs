/*
 * The main usage listing: `gatk` with no arguments, in full.
 *
 * `main-entry` measured its SHAPE -- three hundred and seventy-three lines on stdout, and the first
 * of them -- because the other three hundred and seventy-two are tool names and one-line summaries
 * that no golden carried. This carries them, which is what C.1 has been waiting for (gatk-rs#818).
 *
 * Six behaviours this is built to catch.
 *
 *   - THE LISTING IS GROUPED, and the groups are the tools' program groups, each with its own
 *     heading and its own summary line;
 *   - THE GROUPS AND THE TOOLS INSIDE THEM ARE SORTED, and neither ordering is the one a reader
 *     would guess from the class names;
 *   - EVERY LINE IS COLOURED, in escapes that are bytes of the output like any other;
 *   - A TOOL'S SUMMARY IS TRUNCATED to the width the renderer chose, so the golden holds what was
 *     PRINTED and not what the annotation says;
 *   - A BETA OR EXPERIMENTAL TOOL CARRIES A MARKER, which is part of the line rather than a
 *     property of it;
 *   - AND THE LISTING ENDS WITH A NOTE, whose text is as much of the output as the tool names.
 *
 * The listing goes to System.out for `gatk` with no arguments and to System.err for a name that
 * does not resolve, and it is the same text either way. Both are recorded, so a port that only
 * built one of them cannot pass.
 *
 * The DECLARATIONS travel with the rendering, and that is what makes the comparison a test of the
 * layout rather than a copy of it. Every tool's `@CommandLineProgramProperties` gives its group,
 * its one-line summary and whether Barclay considers it beta or experimental, all of them BEFORE
 * the renderer truncated or padded anything; every group gives its name and its own description.
 * A port that rendered from those and got the same 373 lines has reproduced the layout. A port
 * handed the finished lines would have reproduced nothing.
 *
 * Output:
 *
 *     line\t<stream>\t<index>\t<the line, escaped>
 *     count\t<stream>\t<how many lines>
 *     tool\t<display name>\t<group class>\t<beta|experimental|released>\t<simple name>\t<one-line summary>
 *     group\t<class>\t<name>\t<description>
 *
 * Usage: MainUsageDump
 */

import org.broadinstitute.hellbender.Main;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.List;

public class MainUsageDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String stream, final String name, final String payload) {
        buf.append(kind).append('\t').append(stream).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** The clock log4j puts in front of a line, which is the only thing here that moves. */
    static String masked(final String text) {
        final StringBuilder out = new StringBuilder();
        for (final String line : text.split("\n", -1)) {
            out.append(line.replaceFirst("^\\d\\d:\\d\\d:\\d\\d\\.\\d\\d\\d ", "")).append('\n');
        }
        return out.substring(0, Math.max(0, out.length() - 1));
    }

    /** One run of `instanceMain`, with both streams captured. */
    static void run(final String label, final List<String> argv) {
        final ByteArrayOutputStream outBytes = new ByteArrayOutputStream();
        final ByteArrayOutputStream errBytes = new ByteArrayOutputStream();
        final PrintStream realOut = System.out;
        final PrintStream realErr = System.err;
        try {
            System.setOut(new PrintStream(outBytes, true, StandardCharsets.UTF_8));
            System.setErr(new PrintStream(errBytes, true, StandardCharsets.UTF_8));
            try {
                new Main().instanceMain(argv.toArray(new String[0]));
            } catch (final Exception ignored) {
                // A name that does not resolve prints the listing and THEN throws. The listing is
                // what is being measured, so the exception is the end of the run and not a failure
                // of it; `main-entry` already measures the exception itself.
            }
        } finally {
            System.setOut(realOut);
            System.setErr(realErr);
        }
        record(label + "-out", masked(outBytes.toString(StandardCharsets.UTF_8)));
        record(label + "-err", masked(errBytes.toString(StandardCharsets.UTF_8)));
    }

    static void record(final String stream, final String text) {
        if (text.isEmpty()) {
            emit("count", stream, "lines", "0");
            return;
        }
        final String[] lines = text.split("\n", -1);
        emit("count", stream, "lines", Integer.toString(lines.length));
        for (int i = 0; i < lines.length; i++) {
            emit("line", stream, Integer.toString(i), lines[i]);
        }
    }

    /** Every tool the listing names, with the annotation the renderer read it from.
     *
     * The discovery mirrors `extractCommandLineProgram`'s, which is not one hierarchy but two:
     * Picard's `CommandLineProgram` and GATK's own, over the same packages, minus the two wrapper
     * classes and anything that cannot be instantiated or asked to be omitted.
     */
    static void declarations() throws Exception {
        final Main main = new Main();
        final java.lang.reflect.Method packages = Main.class.getDeclaredMethod("getPackageList");
        packages.setAccessible(true);
        @SuppressWarnings("unchecked")
        final List<String> packageList = (List<String>) packages.invoke(main);

        final org.broadinstitute.barclay.argparser.ClassFinder finder =
                new org.broadinstitute.barclay.argparser.ClassFinder();
        for (final String pkg : packageList) {
            finder.find(pkg, picard.cmdline.CommandLineProgram.class);
            finder.find(pkg, org.broadinstitute.hellbender.cmdline.CommandLineProgram.class);
        }

        final java.util.TreeMap<String, String[]> tools = new java.util.TreeMap<>();
        final java.util.TreeMap<String, String[]> groups = new java.util.TreeMap<>();
        for (final Class<?> type : finder.getClasses()) {
            if (type.equals(org.broadinstitute.hellbender.cmdline.PicardCommandLineProgramExecutor.class)
                    || type.getSimpleName().equals("CommandLineArgumentValidator")) {
                continue;
            }
            if (!org.broadinstitute.hellbender.utils.ClassUtils.canMakeInstances(type)) {
                continue;
            }
            final var properties = type.getAnnotation(
                    org.broadinstitute.barclay.argparser.CommandLineProgramProperties.class);
            if (properties == null || properties.omitFromCommandLine()) {
                continue;
            }
            final Class<?> group = properties.programGroup();
            final var instance =
                    (org.broadinstitute.barclay.argparser.CommandLineProgramGroup)
                            group.getDeclaredConstructor().newInstance();
            groups.put(group.getName(), new String[] {instance.getName(), instance.getDescription()});
            // `toolDisplayName`: a Picard tool carries a suffix, and it is part of the NAME rather
            // than a decoration on the line, which is why it is also what the tools are sorted by.
            final String display =
                    picard.cmdline.CommandLineProgram.class.isAssignableFrom(type)
                            ? type.getSimpleName() + " (Picard)"
                            : type.getSimpleName();
            tools.put(display,
                    new String[] {group.getName(), maturityOf(type), properties.oneLineSummary(),
                            type.getSimpleName()});
        }
        for (final var entry : groups.entrySet()) {
            emit("group", entry.getKey(), entry.getValue()[0], entry.getValue()[1]);
        }
        for (final var entry : tools.entrySet()) {
            // display name, group, maturity, simple name, one-line summary. The simple name is
            // here because the line's PADDING is decided by its length and not the display name's.
            emit("tool", entry.getKey(),
                    entry.getValue()[0] + "\t" + entry.getValue()[1] + "\t" + entry.getValue()[3],
                    entry.getValue()[2]);
        }
    }

    /** `BetaFeature` and `ExperimentalFeature`, which are what the marker in a line comes from. */
    static String maturityOf(final Class<?> type) {
        if (type.getAnnotation(org.broadinstitute.barclay.argparser.ExperimentalFeature.class) != null) {
            return "experimental";
        }
        if (type.getAnnotation(org.broadinstitute.barclay.argparser.BetaFeature.class) != null) {
            return "beta";
        }
        return "released";
    }

    public static void main(final String[] args) throws Exception {
        declarations();
        // The listing on stdout, which is what no arguments answers with.
        run("no-arguments", List.of());
        // And on stderr, which is what a name that does not resolve answers with. The same text on
        // the other stream, and a port that built only one of them would pass a suite that only
        // looked at one.
        run("unknown-name", List.of("NoSuchToolAtAll"));

        System.out.print(buf);
    }
}

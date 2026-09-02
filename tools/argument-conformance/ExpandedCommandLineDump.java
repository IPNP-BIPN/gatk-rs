/*
 * `CommandLineArgumentParser.getCommandLine()`: the string a tool records in every file it writes.
 *
 * A BAM's @PG carries it in CL and a VCF's ##GATKCommandLine carries it in CommandLine, so it is a
 * byte of the output for every record-transform and variant-transform tool there is. A port that
 * could not build it could not write any of their files (gatk-rs#69).
 *
 * Five behaviours this is built to catch.
 *
 *   - IT OPENS WITH THE CLASS'S SIMPLE NAME, which is the TOOL's and not the command line's: it is
 *     `PrintReads`, never `gatk PrintReads` and never the display name a Picard tool carries;
 *   - THE ARGUMENTS COME IN TWO GROUPS: those the user SET, in the parser's own declaration order,
 *     and then those that were not set but have a non-null default, in that same order. Not the
 *     order they were typed in;
 *   - EVERY ONE IS PRINTED AS `--longName value`, whatever the user typed: a short alias is
 *     expanded and a value given with `=` is separated;
 *   - A COLLECTION BECOMES ONE `--name value` PER ELEMENT rather than one with a list after it;
 *   - AND A DEFAULT OF NULL IS OMITTED while a default of anything else is printed, which is why
 *     the line is long even for a command line that set two arguments.
 *
 * Output:
 *
 *     line\t<tool>\t<case>\t<the expanded command line>
 *     error\t<tool>\t<case>\t<exception class>: <message>
 *
 * Usage: ExpandedCommandLineDump
 */

import org.broadinstitute.hellbender.cmdline.CommandLineProgram;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

public class ExpandedCommandLineDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String tool, final String label, final String payload) {
        buf.append(kind).append('\t').append(tool).append('\t').append(label).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static CommandLineProgram instance(final String tool) {
        return switch (tool) {
            case "CountReads" -> new org.broadinstitute.hellbender.tools.CountReads();
            case "CountVariants" -> new org.broadinstitute.hellbender.tools.walkers.CountVariants();
            case "IndexFeatureFile" -> new org.broadinstitute.hellbender.tools.IndexFeatureFile();
            case "CreateHadoopBamSplittingIndex" ->
                    new org.broadinstitute.hellbender.tools.spark.CreateHadoopBamSplittingIndex();
            default -> new org.broadinstitute.hellbender.tools.PrintReads();
        };
    }

    /** One command line, parsed and then asked what it expanded to. */
    static void run(final String tool, final String label, final String... argv) {
        try {
            final CommandLineProgram target = instance(tool);
            final PrintStream sink = new PrintStream(new ByteArrayOutputStream());
            if (!target.getCommandLineParser().parseArguments(sink, argv)) {
                emit("error", tool, label, "not-parsed");
                return;
            }
            emit("line", tool, label, target.getCommandLineParser().getCommandLine());
        } catch (final Exception | AssertionError e) {
            emit("error", tool, label, e.getClass().getName() + ": " + e.getMessage());
        }
    }

    public static void main(final String[] args) {
        // The shortest line there is: one required argument and every default behind it.
        run("PrintReads", "input-and-output", "-I", "/dev/null", "-O", "/dev/null");
        // The same two arguments under their long names, which changes nothing: the line is always
        // the long form.
        run("PrintReads", "long-names", "--input", "/dev/null", "--output", "/dev/null");
        // And with `=`, which the line separates again.
        run("PrintReads", "equals-form", "--input=/dev/null", "--output=/dev/null");
        // A collection with two values, which becomes two `--input` pairs rather than one.
        run("PrintReads", "two-inputs", "-I", "/dev/null", "-I", "/dev/null", "-O", "/dev/null");
        // An argument set to the value it already had by default, which moves it from the second
        // group to the FIRST and therefore changes where it appears in the line.
        run("PrintReads", "default-set-explicitly",
                "-I", "/dev/null", "-O", "/dev/null", "--create-output-bam-index", "true");
        // A boolean flag with no value, which the line prints with one.
        run("PrintReads", "flag-without-a-value",
                "-I", "/dev/null", "-O", "/dev/null", "--add-output-sam-program-record");
        // An interval, which is a collection of a different type.
        run("PrintReads", "with-an-interval",
                "-I", "/dev/null", "-O", "/dev/null", "-L", "chr1:1-100");

        // The other four tools, so the opening name and the defaults behind it are measured for
        // each rather than assumed to follow one pattern.
        run("CountReads", "one-input", "-I", "/dev/null");
        run("CountVariants", "one-variant", "-V", "/dev/null");
        run("IndexFeatureFile", "one-input", "-I", "/dev/null");
        run("CreateHadoopBamSplittingIndex", "one-input", "-I", "/dev/null");

        System.out.print(buf);
    }
}

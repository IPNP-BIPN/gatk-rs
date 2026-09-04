/*
 * The usage text Barclay renders for a tool, taken from the reference.
 *
 * `gatk <Tool> -h` prints a rendering of the same annotations the argument definitions come from:
 * an ordering, a grouping, a wrapping and a set of colour codes, all of them byte-comparable.
 * What is measured is that rendering, for tools of three shapes.
 *
 * Ten behaviours this is built to catch.
 *
 *   - THE TEXT CARRIES NO ANSI ESCAPES, unlike the top-level usage `Main` prints, which is
 *     coloured: `CommandLineArgumentParser.usage` renders plain text and the colour lives above
 *     it;
 *   - THE ARGUMENTS ARE GROUPED into `Required Arguments`, `Optional Arguments` and
 *     `Advanced Arguments`, with a title line apiece and a blank line between entries;
 *   - THE ORDER INSIDE A GROUP IS ALPHABETICAL BY LONG NAME, which is not the order the parser
 *     reports its definitions in;
 *   - AN ARGUMENT'S LINE CARRIES ITS SHORT NAME, ITS LONG NAME, ITS TYPE AND ITS DEFAULT, and the
 *     documentation is wrapped under it;
 *   - THE WRAPPING IS AT A FIXED WIDTH and breaks on spaces, so a long doc string becomes several
 *     lines with the same indent;
 *   - A COLLECTION ARGUMENT SAYS SO in the same place a scalar says its default, with `This
 *     argument may be specified 0 or more times`;
 *   - AN ENUM ARGUMENT LISTS ITS MEMBERS under its own line, as `Possible values: {...}`;
 *   - THE PLUGIN ARGUMENTS ARE RENDERED LIKE ANY OTHER, so a read walker's usage carries the four
 *     read-filter ones and a tool that is no walker carries none;
 *   - THE VERSION IS PRINTED IN THE HEADER, which is why the golden is pinned to 4.6.2.0;
 *   - AND `-h` AND `--help` RENDER THE SAME TEXT, while a tool given neither renders none.
 *
 * Output:
 *
 *     usage\t<tool>\t<the whole text, escaped>
 *     lines\t<tool>\t<how many lines it is>
 *     same\t<tool>\t<whether -h and --help agree>
 *
 * Usage: UsageTextDump
 */

import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;

public class UsageTextDump {

    public static void main(final String[] args) {
        usage("CountReads", new org.broadinstitute.hellbender.tools.CountReads());
        usage("IndexFeatureFile", new org.broadinstitute.hellbender.tools.IndexFeatureFile());
        usage("GatherVcfsCloud", new org.broadinstitute.hellbender.tools.GatherVcfsCloud());
        // The two the port can run, whose usage the dispatcher answers `-h` with. The other two
        // of that shape, CheckTerminatorBlock and BuildBamIndex, are Picard's and are measured
        // there.
        // Undocumented, like `PrintBGZFBlockInformation` below it: the inventory carries a
        // summary only for a DOCUMENTED tool, so the usage golden is where this one's comes from
        // and the declarations generator refuses without it.
        usage("CompareIntervalLists",
                new org.broadinstitute.hellbender.tools.CompareIntervalLists());
        usage("PrintBGZFBlockInformation",
                new org.broadinstitute.hellbender.tools.PrintBGZFBlockInformation());
        usage("CreateHadoopBamSplittingIndex",
                new org.broadinstitute.hellbender.tools.spark.CreateHadoopBamSplittingIndex());
    }

    /** The text the tool's own parser renders, and the two spellings that ask for it. */
    static void usage(final String tool, final org.broadinstitute.hellbender.cmdline
            .CommandLineProgram target) {
        final CommandLineArgumentParser parser =
                (CommandLineArgumentParser) target.getCommandLineParser();
        final String text = parser.usage(true, true);
        System.out.printf("usage\t%s\t%s%n", tool, escape(text));
        System.out.printf("lines\t%s\t%d%n", tool, text.split("\n", -1).length);

        // `-h` and `--help` reach the same rendering, which is what the two runs below show.
        final String dash = render(target, "-h");
        final String longSpelling = render(target, "--help");
        System.out.printf("same\t%s\t%s%n", tool, dash.equals(longSpelling));
    }

    /** What one spelling of the help argument prints. */
    static String render(final org.broadinstitute.hellbender.cmdline.CommandLineProgram target,
                         final String spelling) {
        final ByteArrayOutputStream sink = new ByteArrayOutputStream();
        try {
            target.getCommandLineParser().parseArguments(
                    new PrintStream(sink, true, StandardCharsets.UTF_8),
                    new String[]{spelling});
        } catch (final Exception e) {
            return "E:" + e.getClass().getName();
        }
        return sink.toString(StandardCharsets.UTF_8);
    }

    static String escape(final String text) {
        return text.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
    }
}

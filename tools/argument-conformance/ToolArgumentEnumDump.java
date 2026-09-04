/*
 * The enum-valued arguments of the ported tools, and the constants they accept, from the reference.
 *
 * The declarations golden says an argument's type is `IntervalSetRule`. What it cannot say is what
 * an `IntervalSetRule` is, and a parser needs that twice: to convert a value at all, and to write
 * the message a bad one produces, which lists every constant in declaration order.
 *
 * Six behaviours this is built to catch.
 *
 *   - THE CONSTANTS ARE IN DECLARATION ORDER and not in any sorted one, which is what the message
 *     prints and what `values()` returns;
 *   - THE CONVERSION IS `Enum.valueOf` AND IS THEREFORE CASE SENSITIVE, so `union` is not `UNION`
 *     and the refusal names both the value and the type;
 *   - THE MESSAGE LISTS EVERY CONSTANT, which is why the list is measured rather than the count;
 *   - AN ENUM THAT IMPLEMENTS `ClpEnum` DOCUMENTS ITS CONSTANTS, and that documentation is part of
 *     the usage text rather than of the refusal, so the two are measured apart;
 *   - THE SAME TYPE APPEARS UNDER MORE THAN ONE TOOL and is one type, so the table is by type and
 *     the arguments point into it;
 *   - AND A DEFAULT IS ONE OF THE CONSTANTS, which is what makes an unset enum argument optional.
 *
 * Output:
 *
 *     enum\t<type>\t<constants, comma separated, in declaration order>
 *     clp\t<type>\t<constant>=<the documentation ClpEnum gives it, escaped>
 *     arg\t<tool>\t<long name>\t<type>|<default>
 *     parse\t<tool>\t<case>\tok|E:<exception class>:<message>
 *
 * Usage: ToolArgumentEnumDump
 */

import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.barclay.argparser.NamedArgumentDefinition;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class ToolArgumentEnumDump {

    /** One entry per enum type, in the order the tools first mention it. */
    static final Map<String, Class<?>> types = new LinkedHashMap<>();
    static final StringBuilder args = new StringBuilder();

    public static void main(final String[] args) {
        declarations("CountReads", new org.broadinstitute.hellbender.tools.CountReads());
        declarations("CountVariants",
                new org.broadinstitute.hellbender.tools.walkers.CountVariants());
        declarations("PrintReads", new org.broadinstitute.hellbender.tools.PrintReads());
        declarations("ApplyBQSR",
                new org.broadinstitute.hellbender.tools.walkers.bqsr.ApplyBQSR());
        declarations("SelectVariants",
                new org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants());
        declarations("IndexFeatureFile",
                new org.broadinstitute.hellbender.tools.IndexFeatureFile());
        declarations("GatherVcfsCloud",
                new org.broadinstitute.hellbender.tools.GatherVcfsCloud());
        // The tools that joined the declarations after this dump was written. Three of them point
        // at enums the seven above already name; `SplitIntervals` brings one of its own, and an
        // argument whose type is missing from this table is DROPPED by the port's parser -- it
        // answered `subdivision-mode is not a recognized option` for a name the reference accepts.
        declarations("CountBases", new org.broadinstitute.hellbender.tools.CountBases());
        declarations("FlagStat", new org.broadinstitute.hellbender.tools.FlagStat());
        declarations("CountBasesInReference",
                new org.broadinstitute.hellbender.tools.walkers.fasta.CountBasesInReference());
        declarations("SplitIntervals",
                new org.broadinstitute.hellbender.tools.walkers.SplitIntervals());
        // A LOCUS walker, whose traversal is the pileup at each position rather than a record or
        // a base, and a second interval utility, which shares `SplitIntervals`' plumbing. Both
        // dumps gain them together: an enum missing from the enums table is an argument the
        // port's parser drops whole, which cost `SplitIntervals` a whole CI round.
        declarations("Pileup",
                new org.broadinstitute.hellbender.tools.walkers.qc.Pileup());
        declarations("PreprocessIntervals",
                new org.broadinstitute.hellbender.tools.copynumber.PreprocessIntervals());
        // A read walker that WRITES a BAM, which is `PrintReads`' plumbing with a filter in front
        // of it, and a `GATKTool` that opens the reads only to read their header. Both dumps
        // together, as ever.
        declarations("PrintDistantMates",
                new org.broadinstitute.hellbender.tools.PrintDistantMates());
        declarations("GetSampleName",
                new org.broadinstitute.hellbender.tools.GetSampleName());
        // A `CommandLineProgram` that is no GATKTool at all -- fifteen arguments, two interval
        // lists in and a verdict out -- and a read walker that rewrites qualities and writes a BAM.
        declarations("CompareIntervalLists",
                new org.broadinstitute.hellbender.tools.CompareIntervalLists());
        declarations("FixMisencodedBaseQualityReads",
                new org.broadinstitute.hellbender.tools.FixMisencodedBaseQualityReads());
        // A second LOCUS walker, which compares the engine's pileup to a samtools one, and a
        // third interval utility, which annotates each interval with its GC content.
        declarations("CheckPileup",
                new org.broadinstitute.hellbender.tools.walkers.qc.CheckPileup());
        declarations("AnnotateIntervals",
                new org.broadinstitute.hellbender.tools.copynumber.AnnotateIntervals());

        // A second REFERENCE walker, which writes a FASTA rather than a number, and a third LOCUS
        // walker, which counts the reads over each interval. `CollectReadCounts` brings an enum of
        // its own -- `--format`, whose two constants are the TSV and the HDF5 it writes -- and an
        // enum missing from the enums table is an argument the port's parser drops whole, so the
        // two dumps gain the pair together as ever.
        declarations("FastaReferenceMaker",
                new org.broadinstitute.hellbender.tools.walkers.fasta.FastaReferenceMaker());
        declarations("CollectReadCounts",
                new org.broadinstitute.hellbender.tools.copynumber.CollectReadCounts());

        // The table first, then the arguments that point into it.
        final List<String> names = new ArrayList<>(types.keySet());
        java.util.Collections.sort(names);
        for (final String name : names) {
            final Class<?> type = types.get(name);
            final List<String> constants = new ArrayList<>();
            for (final Object constant : type.getEnumConstants()) {
                constants.add(((Enum<?>) constant).name());
            }
            System.out.printf("enum\t%s\t%s%n", name, String.join(",", constants));
            // A `ClpEnum` documents each of its constants, and that documentation is the usage
            // text's rather than the refusal's.
            if (CommandLineParser.ClpEnum.class.isAssignableFrom(type)) {
                for (final Object constant : type.getEnumConstants()) {
                    System.out.printf("clp\t%s\t%s=%s%n", name, ((Enum<?>) constant).name(),
                            escape(((CommandLineParser.ClpEnum) constant).getHelpDoc()));
                }
            }
        }
        System.out.print(ToolArgumentEnumDump.args);

        // What a value the enum does not carry costs, which is where the constants are printed.
        parse("CountReads", "a-constant", new String[]{"-I", "/dev/null", "-isr", "UNION"});
        parse("CountReads", "the-other-constant",
                new String[]{"-I", "/dev/null", "-isr", "INTERSECTION"});
        parse("CountReads", "lower-case", new String[]{"-I", "/dev/null", "-isr", "union"});
        parse("CountReads", "not-a-constant", new String[]{"-I", "/dev/null", "-isr", "NEITHER"});
        parse("CountReads", "an-empty-value", new String[]{"-I", "/dev/null", "-isr", ""});
        // A second type, so the message is not one type's own shape.
        parse("CountReads", "a-stringency",
                new String[]{"-I", "/dev/null", "-VS", "LENIENT"});
        parse("CountReads", "not-a-stringency",
                new String[]{"-I", "/dev/null", "-VS", "PERMISSIVE"});
    }

    /** Every enum-valued argument of one tool, and the type it points at. */
    static void declarations(final String tool, final Object target) {
        final List<NamedArgumentDefinition> definitions =
                ((CommandLineArgumentParser) ((org.broadinstitute.hellbender.cmdline
                        .CommandLineProgram) target).getCommandLineParser())
                        .getNamedArgumentDefinitions();
        for (final NamedArgumentDefinition definition : definitions) {
            final Class<?> type = definition.getUnderlyingFieldClass();
            if (!type.isEnum()) {
                continue;
            }
            types.putIfAbsent(type.getSimpleName(), type);
            args.append(String.format("arg\t%s\t%s\t%s|%s%n", tool, definition.getLongName(),
                    type.getSimpleName(),
                    escape(String.valueOf(definition.getDefaultValueAsString()))));
        }
    }

    static String escape(final String text) {
        return text.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
    }

    /** What the parser makes of one command line, before the tool runs at all. */
    static void parse(final String tool, final String label, final String[] argv) {
        final Object target = new org.broadinstitute.hellbender.tools.CountReads();
        String result;
        try {
            final PrintStream sink = new PrintStream(new ByteArrayOutputStream());
            result = ((org.broadinstitute.hellbender.cmdline.CommandLineProgram) target)
                    .getCommandLineParser().parseArguments(sink, argv) ? "ok" : "not-parsed";
        } catch (final Exception | AssertionError e) {
            result = "E:" + e.getClass().getName() + ":" + e.getMessage();
        }
        System.out.printf("parse\t%s\t%s\t%s%n", tool, label, escape(result));
    }
}

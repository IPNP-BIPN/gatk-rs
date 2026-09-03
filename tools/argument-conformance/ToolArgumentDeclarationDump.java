/*
 * The argument declarations of three ported tools, taken from the reference.
 *
 * The Barclay parser is measured already, argument by argument and mechanism by mechanism. What is
 * not measured is what any ONE tool declares: the flattened namespace its own fields and its
 * inherited collections produce, and therefore which command lines it accepts. That is the layer
 * Milestone C's per-tool declarations have to reproduce, and it cannot be read off the parser.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE NAMESPACE IS FLAT AND INHERITED: a read walker declares a handful of arguments of its
 *     own and ends up with dozens, the rest coming from the collections its superclasses hold;
 *   - AND WHICH PARSER IS ASKED DECIDES HOW MANY THERE ARE: a parser built straight from the
 *     instance sees 38 where the one the tool hands out sees 70, the gap being the plugin
 *     descriptors' four read-filter arguments and the standard collections the tool adds when it
 *     builds its own. Every one of the 38 is in the 70, so the instance list is a subset and not
 *     a different reading;
 *   - THE ORDER IS THE PARSER'S OWN, subclass first, and it is what the usage prints in;
 *   - A REQUIRED ARGUMENT IS REQUIRED BY ITS DECLARATION AND NOT BY THE TOOL, and the two
 *     archetypes do not mirror each other: a read walker REQUIRES `--input` and has no `--variant`
 *     at all, while a variant walker requires `--variant` and takes `--input` as an OPTIONAL
 *     argument, so `-V` on the first is "not a recognized option" and `-I` on the second parses;
 *   - THE SHORT NAME IS PART OF THE DECLARATION, and several arguments have none;
 *   - A COLLECTION ARGUMENT IS DECLARED AS ONE, which is what lets it be repeated;
 *   - THE DEFAULT IS THE FIELD'S OWN VALUE at construction time, so it is read off an instance and
 *     not off the annotation;
 *   - AN ARGUMENT NAMED TWICE IS AN ERROR ONLY IF IT IS A SCALAR: `--input` is declared as a
 *     collection, so naming it twice parses, while `--output` is a scalar and naming it twice is
 *     refused;
 *   - AN UNKNOWN ARGUMENT IS REFUSED BY THE PARSER AND NOT BY THE TOOL, with a message that names
 *     it;
 *   - AND A MISSING REQUIRED ARGUMENT IS REFUSED THE SAME WAY, before the tool runs at all.
 *
 * Five more, added when the declarations had to carry enough to BUILD a parser rather than only to
 * count one.
 *
 *   - THE TYPE IS THE UNDERLYING FIELD'S, which for a collection is its ELEMENT class: `--input`
 *     is a `List<GATKPath>` and reports `GATKPath`, so the conversion a value goes through is the
 *     element's and the collection is only how many of them there may be;
 *   - PRIMITIVE IS A SEPARATE QUESTION FROM THE CLASS, the class being boxed either way, and it is
 *     the one the null check asks;
 *   - HIDDEN, ADVANCED AND COMMON ARE THREE DIFFERENT FLAGS, and they are what decides which
 *     section of the usage an argument is printed in, or whether it is printed at all;
 *   - A BOUNDED RANGE IS FOUR NULLABLE DOUBLES and not a pair, the recommended range being
 *     declared beside the hard one rather than instead of it;
 *   - AND THE DOCUMENTATION IS PART OF THE DECLARATION: it is the annotation's own string, which
 *     is what the usage text wraps, so it belongs to the argument and not to the renderer.
 *
 * Output:
 *
 *     count\t<tool>\t<how many named arguments it declares>
 *     def\t<tool>\t<index>\t<longName>|<aliases>|<required>|<collection>|<default>|<type>|
 *         <primitive>|<flag>|<hidden>|<advanced>|<common>|<minElements>|<maxElements>|
 *         <minValue>|<maxValue>|<minRecommended>|<maxRecommended>|<mutex>|<plugin>
 *     doc\t<tool>\t<index>\t<the documentation string, escaped>
 *     parse\t<tool>\t<case>\tok|E:<exception class>:<message>
 *
 * Usage: ToolArgumentDeclarationDump
 */

import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.NamedArgumentDefinition;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class ToolArgumentDeclarationDump {

    public static void main(final String[] args) {
        declarations("CountReads",
                new org.broadinstitute.hellbender.tools.CountReads());
        declarations("CountVariants",
                new org.broadinstitute.hellbender.tools.walkers.CountVariants());
        declarations("PrintReads",
                new org.broadinstitute.hellbender.tools.PrintReads());
        // Four more archetypes, so the list is not three walkers of two kinds: a read walker that
        // writes a recalibrated file, a variant walker with a large argument surface of its own, a
        // tool that is no walker at all, and one that takes a list of files rather than one.
        declarations("ApplyBQSR",
                new org.broadinstitute.hellbender.tools.walkers.bqsr.ApplyBQSR());
        declarations("SelectVariants",
                new org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants());
        declarations("IndexFeatureFile",
                new org.broadinstitute.hellbender.tools.IndexFeatureFile());
        // `GatherVcfs` is Picard's tool of that name, which GATK dispatches to; the GATK one is
        // GatherVcfsCloud, and the two declare different arguments. Naming it wrongly here made
        // the generator's cross-check report a `--COMMENT` the parser did not declare, which is
        // Picard's argument and not this tool's.
        declarations("GatherVcfsCloud",
                new org.broadinstitute.hellbender.tools.GatherVcfsCloud());
        // Two more tools that are no walkers, chosen because the port can RUN them: each takes a
        // file and answers, so a declaration for one of them is a command line the binary can be
        // handed end to end rather than only parsed.
        //
        // `CheckTerminatorBlock` and `BuildBamIndex` would have been two more and are not here:
        // they are PICARD's tools, which GATK dispatches to, so their declarations come from
        // Picard's own parser and belong in picard-rs's measurement rather than this one.
        // The two read walkers that share `CountReads`' plumbing exactly: reads in, a number or a
        // report out. They are here because an archetype is the unit of this milestone -- a tool
        // of one archetype shares its argument shape AND its file plumbing with the others, so the
        // second and third cost a fraction of the first.
        declarations("CountBases",
                new org.broadinstitute.hellbender.tools.CountBases());
        declarations("FlagStat",
                new org.broadinstitute.hellbender.tools.FlagStat());
        // Two archetypes neither the read walkers nor the variant walkers reach: a REFERENCE
        // walker, whose traversal is the FASTA rather than a file of records, and an interval
        // utility, which is a `GATKTool` with no traversal at all and writes a DIRECTORY of files.
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
        declarations("PrintBGZFBlockInformation",
                new org.broadinstitute.hellbender.tools.PrintBGZFBlockInformation());
        declarations("CreateHadoopBamSplittingIndex",
                new org.broadinstitute.hellbender.tools.spark.CreateHadoopBamSplittingIndex());

        // A read walker: its own input, and the arguments it inherits.
        parse("CountReads", "no-arguments", new String[]{});
        parse("CountReads", "input-only", new String[]{"-I", "/dev/null"});
        parse("CountReads", "long-input", new String[]{"--input", "/dev/null"});
        parse("CountReads", "input-twice",
                new String[]{"-I", "/dev/null", "-I", "/dev/null"});
        parse("CountReads", "unknown-argument", new String[]{"--no-such-argument", "1"});
        parse("CountReads", "an-interval", new String[]{"-I", "/dev/null", "-L", "chr1"});
        parse("CountReads", "two-intervals",
                new String[]{"-I", "/dev/null", "-L", "chr1", "-L", "chr2"});
        parse("CountReads", "a-variant-argument",
                new String[]{"-I", "/dev/null", "-V", "/dev/null"});

        // A variant walker, which has no `--input` at all.
        parse("CountVariants", "variant-only", new String[]{"-V", "/dev/null"});
        parse("CountVariants", "an-input", new String[]{"-V", "/dev/null", "-I", "/dev/null"});
        parse("CountVariants", "no-arguments", new String[]{});

        // The output argument a print tool requires on top of its input.
        parse("PrintReads", "input-only", new String[]{"-I", "/dev/null"});
        // A tool that is no walker takes a positional-looking input and nothing else required.
        parse("IndexFeatureFile", "no-arguments", new String[]{});
        parse("IndexFeatureFile", "input-only", new String[]{"-I", "/dev/null"});
        parse("IndexFeatureFile", "an-interval", new String[]{"-I", "/dev/null", "-L", "chr1"});

        parse("PrintReads", "input-and-output",
                new String[]{"-I", "/dev/null", "-O", "/dev/null"});
        // The output is a scalar where the input is a collection, so naming it twice is refused.
        parse("PrintReads", "output-twice", new String[]{
            "-I", "/dev/null", "-O", "/dev/null", "-O", "/dev/null"});
    }

    /**
     * Every named argument the tool declares, in the parser's own order, and how many a parser
     * built straight from the instance would have seen instead.
     *
     * The two are not the same list. A parser constructed from the instance knows nothing about
     * the plugin descriptors or the standard argument collections the tool adds when it builds its
     * own, so it sees 38 arguments where the tool's parser sees 70. The four read-filter arguments
     * are the visible half of that gap; the rest are the common and advanced ones the usage text
     * does not print either.
     */
    static void declarations(final String tool, final Object target) {
        final List<NamedArgumentDefinition> instance =
                new CommandLineArgumentParser(target).getNamedArgumentDefinitions();
        final List<NamedArgumentDefinition> definitions =
                ((CommandLineArgumentParser) ((org.broadinstitute.hellbender.cmdline
                        .CommandLineProgram) target).getCommandLineParser())
                        .getNamedArgumentDefinitions();
        System.out.printf("count\t%s\tinstance=%d tool=%d%n", tool, instance.size(),
                definitions.size());
        final List<String> seen = new ArrayList<>();
        for (final NamedArgumentDefinition definition : instance) {
            seen.add(definition.getLongName());
        }
        for (int i = 0; i < definitions.size(); i++) {
            final NamedArgumentDefinition definition = definitions.get(i);
            final List<String> shorts = new ArrayList<>(definition.getArgumentAliases());
            // The class the parser converts a value to is the UNDERLYING field's, which for a
            // collection is its element class and not the collection's own.
            final String type = definition.getUnderlyingFieldClass().getSimpleName();
            final boolean primitive = definition.getUnderlyingField().getType().isPrimitive();
            // NOT sorted. The message a mutex violation prints joins this list in the order it
            // holds, and that order is neither alphabetical nor the annotation's: `quantize-quals`
            // declares no mutex at all, and Barclay fills its list in reverse as it walks the
            // declarations, so it reads `static-quantized-quals round-down-quantized` -- the order
            // those two are declared in. Sorting it here made the generated declarations hold the
            // alphabetical order and the port print a message the reference never writes.
            final List<String> mutex = new ArrayList<>(definition.getMutexTargetList());
            final Object plugin = definition.getDescriptorForControllingPlugin();
            System.out.printf(
                    "def\t%s\t%d\t%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%d|%d|%s|%s|%s|%s|%s|%s%n",
                    tool, i,
                    definition.getLongName(),
                    String.join(",", shorts),
                    definition.isOptional() ? "optional" : "required",
                    definition.isCollection() ? "collection" : "scalar",
                    escape(String.valueOf(definition.getDefaultValueAsString())),
                    type,
                    primitive ? "primitive" : "boxed",
                    definition.isFlag() ? "flag" : "valued",
                    definition.isHidden() ? "hidden" : "printed",
                    definition.isAdvanced() ? "advanced" : "plain",
                    definition.isCommon() ? "common" : "own",
                    definition.getMinElements(),
                    definition.getMaxElements(),
                    String.valueOf(definition.getMinValue()),
                    String.valueOf(definition.getMaxValue()),
                    String.valueOf(definition.getMinRecommendedValue()),
                    String.valueOf(definition.getMaxRecommendedValue()),
                    mutex.isEmpty() ? "none" : String.join(",", mutex),
                    plugin == null ? "none" : plugin.getClass().getSimpleName());
            // The documentation is a line of its own: it is prose, and prose carries the pipe the
            // line above uses as its separator.
            System.out.printf("doc\t%s\t%d\t%s%n", tool, i,
                    escape(String.valueOf(definition.getDocString())));
            if (!seen.contains(definition.getLongName())) {
                System.out.printf("only-on-the-tool\t%s\t%s%n", tool,
                        definition.getLongName());
            }
        }
    }

    /** The dump's own escaping, since this harness has no shared helper. */
    static String escape(final String text) {
        return text.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
    }

    /** What the parser makes of one command line, before the tool runs at all. */
    static void parse(final String tool, final String label, final String[] argv) {
        final Object target = switch (tool) {
            case "CountReads" -> new org.broadinstitute.hellbender.tools.CountReads();
            case "CountVariants" -> new org.broadinstitute.hellbender.tools.walkers.CountVariants();
            case "IndexFeatureFile" -> new org.broadinstitute.hellbender.tools.IndexFeatureFile();
            default -> new org.broadinstitute.hellbender.tools.PrintReads();
        };
        String result;
        try {
            final PrintStream sink = new PrintStream(new ByteArrayOutputStream());
            result = ((org.broadinstitute.hellbender.cmdline.CommandLineProgram) target)
                    .getCommandLineParser().parseArguments(sink, argv) ? "ok" : "not-parsed";
        } catch (final Exception | AssertionError e) {
            result = "E:" + e.getClass().getName() + ":" + e.getMessage();
        }
        System.out.printf("parse\t%s\t%s\t%s%n", tool, label,
                escape(result));
    }
}

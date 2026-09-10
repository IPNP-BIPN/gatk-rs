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

        // A third LOCUS walker, which is the first tool here whose traversal drives a FEATURE
        // source -- its `-V` is the allele frequencies it summarises the pileup at -- and a read
        // walker that clears one flag bit and writes the reads back, which is `PrintReads`'
        // plumbing with a mutation in the middle.
        declarations("GetPileupSummaries",
                new org.broadinstitute.hellbender.tools.walkers.contamination.GetPileupSummaries());
        declarations("UnmarkDuplicates",
                new org.broadinstitute.hellbender.tools.walkers.UnmarkDuplicates());

        // A `CommandLineProgram` that is no GATKTool -- it reads the table `GetPileupSummaries`
        // writes, so the two together are the first CHAIN of ported tools a command line can run
        // end to end -- and a `GATKTool` with no traversal at all, which opens the reads only to
        // print their header.
        declarations("CalculateContamination",
                new org.broadinstitute.hellbender.tools.walkers.contamination.CalculateContamination());
        declarations("PrintReadsHeader",
                new org.broadinstitute.hellbender.tools.PrintReadsHeader());
        // A third LOCUS walker, which writes TWO files -- a BED of runs and a summary of counts --
        // and a second REFERENCE utility, which writes a shifted FASTA with its dictionary and two
        // interval lists beside it. Neither shape has been declared here before: every tool
        // measured so far writes one file or a directory of them.
        declarations("CallableLoci",
                new org.broadinstitute.hellbender.tools.walkers.coverage.CallableLoci());
        declarations("ShiftFasta",
                new org.broadinstitute.hellbender.tools.walkers.fasta.ShiftFasta());
        declarations("RevertBaseQualityScores",
                new org.broadinstitute.hellbender.tools.walkers.RevertBaseQualityScores());
        declarations("AddOriginalAlignmentTags",
                new org.broadinstitute.hellbender.tools.AddOriginalAlignmentTags());
        // The next two. `LeftAlignIndels` is a read walker that needs a REFERENCE, which is the
        // first required argument in this dump that is not the reads, and `DumpTabixIndex` is a
        // tool that is no walker at all and whose whole namespace is fourteen arguments.
        declarations("LeftAlignIndels",
                new org.broadinstitute.hellbender.tools.LeftAlignIndels());
        declarations("DumpTabixIndex",
                new org.broadinstitute.hellbender.tools.DumpTabixIndex());
        // The next two. `ReadAnonymizer` is a read walker that needs a reference AND rewrites the
        // bases it reads, and `PrintFileDiagnostics` declares fifteen arguments and is no walker.
        declarations("ReadAnonymizer",
                new org.broadinstitute.hellbender.tools.walkers.ReadAnonymizer());
        declarations("PrintFileDiagnostics",
                new org.broadinstitute.hellbender.tools.PrintFileDiagnostics());
        // The next two, and both write MORE than one output: `SplitReads` writes one file per
        // key of up to three splitters, and `ClipReads` writes a statistics file beside its BAM.
        declarations("SplitReads",
                new org.broadinstitute.hellbender.tools.SplitReads());
        declarations("ClipReads",
                new org.broadinstitute.hellbender.tools.ClipReads());
        // The next two, and both need a REFERENCE: `SplitNCigarReads` splits a read at every `N`
        // and `MethylationTypeCaller` writes a VCF rather than reads.
        declarations("SplitNCigarReads",
                new org.broadinstitute.hellbender.tools.walkers.rnaseq.SplitNCigarReads());
        declarations("MethylationTypeCaller",
                new org.broadinstitute.hellbender.tools.walkers.MethylationTypeCaller());
        // The next two. `BaseRecalibrator` writes a GATKReport rather than reads and takes a
        // FeatureInput of known sites, and `GtfToBed` reads an annotation and writes a BED.
        declarations("BaseRecalibrator",
                new org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator());
        declarations("GtfToBed",
                new org.broadinstitute.hellbender.tools.walkers.conversion.GtfToBed());
        // Two tools of the record-transform archetype that are no WALKERS: both extend `GATKTool`
        // and override `traverse()`, so a second reads source is opened by hand and the engine's
        // filter, transformer and interval machinery never runs. `TransferReadTags` walks two files
        // in lockstep and copies tags from the unmapped one, and `PostProcessReadsForRSEM` reorders
        // a query-name-sorted file into the pairs RSEM will read.
        declarations("TransferReadTags",
                new org.broadinstitute.hellbender.tools.walkers.qc.TransferReadTags());
        declarations("PostProcessReadsForRSEM",
                new org.broadinstitute.hellbender.tools.walkers.qc.PostProcessReadsForRSEM());
        // A `VariantWalker` that writes a TABLE rather than a VCF, and a tool that is no GATK tool
        // at all: `CompareBaseQualities` extends `PicardCommandLineProgram`, so its namespace is
        // Picard's argument set and not the engine's, which nothing declared here has been.
        declarations("VariantsToTable",
                new org.broadinstitute.hellbender.tools.walkers.variantutils.VariantsToTable());
        declarations("CompareBaseQualities",
                new org.broadinstitute.hellbender.tools.validation.CompareBaseQualities());
        // The first two tools here that WRITE a VCF, which is a writer this port has not used from a
        // runner before. `RemoveNearbyIndels` buffers one indel at a time and drops any pair closer
        // than a spacing; `UpdateVCFSequenceDictionary` replaces the header's dictionary and passes
        // every record through, and its `--source-dictionary` is the first argument here that takes
        // a dictionary from any of four file kinds.
        declarations("RemoveNearbyIndels",
                new org.broadinstitute.hellbender.tools.walkers.validation.RemoveNearbyIndels());
        declarations("UpdateVCFSequenceDictionary",
                new org.broadinstitute.hellbender.tools.walkers.variantutils.UpdateVCFSequenceDictionary());
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

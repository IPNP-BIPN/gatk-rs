/*
 * The multi-pass walkers, taken from the reference.
 *
 * ReadWalker, LocusWalker, IntervalWalker, VariantWalker and AssemblyRegionWalker each make one
 * traversal. Four classes in the engine package make several, and what they do between passes is
 * not the same in any two of them:
 *
 *   - MultiplePassVariantWalker builds ONE CountingVariantFilter and ONE CountingReadFilter before
 *     the loop and reuses both for every pass, so the filtered counts accumulate across passes: a
 *     three-pass run over a file with two filtered records reports six, not two;
 *   - MultiplePassReadWalker builds its filter in traverse() and then builds a NEW one at the top
 *     of every pass after the first, and resets the reads data source with it, so each pass after
 *     the first counts only its own drops. The two classes are one directory apart and disagree;
 *   - afterNthPass is called after EVERY pass including the last, so a two-pass run calls it twice;
 *   - TwoPassVariantWalker routes n==0 to firstPassApply and n==1 to secondPassApply, and routes
 *     afterNthPass(0) to afterFirstPass and afterNthPass(1) to NOTHING AT ALL. Its guard is
 *     `else if (n > 1) throw`, so there is no afterSecondPass and no error either;
 *   - the two base classes disagree about apply(): MultiplePassVariantWalker makes it an empty
 *     final method, while MultiplePassReadWalker makes it a final method that throws GATKException.
 *
 * numberOfPasses() is consulted once per loop iteration by `n < numberOfPasses()`, so zero passes
 * is a legal traversal that visits nothing and calls afterNthPass never.
 *
 * Output:
 *
 *     event\t<label>\t<n>\t<event>
 *     filter\t<label>\t<n>\t<filtered count of the nth filter this run built>
 *     filters-built\t<label>\t<how many times makeXFilter was called>
 *     summary\t<label>\t<ok|E:class>
 *
 * Usage: MultiPassWalkerDump
 */

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.MultiplePassReadWalker;
import org.broadinstitute.hellbender.engine.MultiplePassVariantWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.TwoPassVariantWalker;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.CountingVariantFilter;
import org.broadinstitute.hellbender.engine.filters.VariantFilterLibrary;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class MultiPassWalkerDump {

    /** The callbacks of the current traversal, in the order they happened. */
    static final List<String> EVENTS = new ArrayList<>();

    /** Every filter object the current traversal built, in the order it built them. */
    static final List<CountingVariantFilter> VARIANT_FILTERS = new ArrayList<>();
    static final List<CountingReadFilter> READ_FILTERS = new ArrayList<>();

    /** How many passes the configurable probes should make. */
    static int PASSES = 2;

    static String where(final VariantContext variant) {
        return variant.getContig() + ":" + variant.getStart();
    }

    /**
     * A two-pass walker.
     *
     * The variant filter is overridden to one that actually drops something, which the base class
     * invites ("Subclasses can extend to provide own filters"): the default allows every variant,
     * and a counter that never counts cannot show whether it was reused.
     */
    @CommandLineProgramProperties(
            summary = "Records the callback order of a TwoPassVariantWalker",
            oneLineSummary = "TwoPassVariantWalker probe",
            programGroup = ReadDataManipulationProgramGroup.class)
    public static final class TwoPassProbe extends TwoPassVariantWalker {
        @Override
        protected CountingVariantFilter makeVariantFilter() {
            final CountingVariantFilter filter =
                    new CountingVariantFilter(VariantFilterLibrary.PASSES_FILTERS);
            VARIANT_FILTERS.add(filter);
            return filter;
        }

        @Override
        protected void firstPassApply(final VariantContext variant, final ReadsContext reads,
                                      final ReferenceContext reference, final FeatureContext features) {
            EVENTS.add("firstPassApply " + where(variant));
        }

        @Override
        protected void afterFirstPass() {
            EVENTS.add("afterFirstPass");
        }

        @Override
        protected void secondPassApply(final VariantContext variant, final ReadsContext reads,
                                       final ReferenceContext reference, final FeatureContext features) {
            EVENTS.add("secondPassApply " + where(variant));
        }

        @Override
        public Object onTraversalSuccess() {
            EVENTS.add("onTraversalSuccess");
            return super.onTraversalSuccess();
        }
    }

    /** An n-pass variant walker, where n is whatever PASSES says. */
    @CommandLineProgramProperties(
            summary = "Records the callback order of a MultiplePassVariantWalker",
            oneLineSummary = "MultiplePassVariantWalker probe",
            programGroup = ReadDataManipulationProgramGroup.class)
    public static final class NPassProbe extends MultiplePassVariantWalker {
        @Override
        protected int numberOfPasses() {
            return PASSES;
        }

        @Override
        protected CountingVariantFilter makeVariantFilter() {
            final CountingVariantFilter filter =
                    new CountingVariantFilter(VariantFilterLibrary.PASSES_FILTERS);
            VARIANT_FILTERS.add(filter);
            return filter;
        }

        @Override
        protected void nthPassApply(final VariantContext variant, final ReadsContext reads,
                                    final ReferenceContext reference, final FeatureContext features,
                                    final int n) {
            EVENTS.add("nthPassApply " + n + " " + where(variant));
        }

        @Override
        protected void afterNthPass(final int n) {
            EVENTS.add("afterNthPass " + n);
        }
    }

    /** An n-pass read walker, which calls forEachRead once per requested pass. */
    @CommandLineProgramProperties(
            summary = "Records the callback order of a MultiplePassReadWalker",
            oneLineSummary = "MultiplePassReadWalker probe",
            programGroup = ReadDataManipulationProgramGroup.class)
    public static final class MultiReadProbe extends MultiplePassReadWalker {
        @Override
        public CountingReadFilter makeReadFilter() {
            final CountingReadFilter filter = super.makeReadFilter();
            READ_FILTERS.add(filter);
            return filter;
        }

        @Override
        public void traverseReads() {
            for (int pass = 0; pass < PASSES; pass++) {
                final int n = pass;
                forEachRead((read, reference, features) ->
                        EVENTS.add("read " + n + " " + read.getName()));
                EVENTS.add("endOfPass " + n);
            }
        }
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("multipass-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        // The variant fixture is VariantWalkerDump's, so the two suites describe the same records.
        // Two of its six are filtered: the LowQual one and the one with a "." FILTER column, which
        // is what makes PASSES_FILTERS count anything at all.
        final Path vcf = dir.resolve("variants.vcf");
        Files.write(vcf, VariantWalkerDump.VCF.getBytes());
        new org.broadinstitute.hellbender.tools.IndexFeatureFile()
                .instanceMain(new String[] {"-I", vcf.toString()});

        // And the read fixture is ReadWalkerDump's: eleven records, of which WellformedReadFilter
        // drops two.
        final Path bam = dir.resolve("reads.bam");
        ReadWalkerDump.buildFixture(bam.toFile());

        System.out.println("# MultiPassWalkerDump: the multi-pass walkers");

        PASSES = 2;
        variants("two-pass", new TwoPassProbe(), vcf);
        variants("npass-2", new NPassProbe(), vcf);
        PASSES = 3;
        variants("npass-3", new NPassProbe(), vcf);
        PASSES = 1;
        variants("npass-1", new NPassProbe(), vcf);
        // `for (n = 0; n < numberOfPasses(); n++)` with zero passes: a traversal that visits
        // nothing, calls afterNthPass never, and still logs both filter summaries.
        PASSES = 0;
        variants("npass-0", new NPassProbe(), vcf);

        // The same three counts over reads, where the filter is rebuilt per pass rather than reused.
        PASSES = 2;
        reads("read-two-pass", bam);
        PASSES = 3;
        reads("read-three-pass", bam);
        PASSES = 1;
        reads("read-one-pass", bam);
        PASSES = 0;
        reads("read-zero-pass", bam);

        // A bounded traversal, to show the data source is reset with its bounds rather than
        // without them: pass 2 sees the same reads as pass 1, not the whole file.
        PASSES = 2;
        reads("read-two-pass-interval", bam, "-L", "chr1:100-200");
    }

    static void variants(final String label, final Object probe, final Path vcf,
                         final String... extra) {
        EVENTS.clear();
        VARIANT_FILTERS.clear();
        READ_FILTERS.clear();
        final List<String> argv = new ArrayList<>(Arrays.asList("-V", vcf.toString()));
        argv.addAll(Arrays.asList(extra));

        String summary;
        try {
            if (probe instanceof TwoPassProbe) {
                ((TwoPassProbe) probe).instanceMain(argv.toArray(new String[0]));
            } else {
                ((NPassProbe) probe).instanceMain(argv.toArray(new String[0]));
            }
            summary = "ok";
        } catch (final Exception | AssertionError e) {
            summary = "E:" + e.getClass().getName();
        }
        report(label, summary, true);
    }

    static void reads(final String label, final Path bam, final String... extra) {
        EVENTS.clear();
        VARIANT_FILTERS.clear();
        READ_FILTERS.clear();
        final List<String> argv = new ArrayList<>(Arrays.asList("-I", bam.toString()));
        argv.addAll(Arrays.asList(extra));

        String summary;
        try {
            new MultiReadProbe().instanceMain(argv.toArray(new String[0]));
            summary = "ok";
        } catch (final Exception | AssertionError e) {
            summary = "E:" + e.getClass().getName();
        }
        report(label, summary, false);
    }

    static void report(final String label, final String summary, final boolean variants) {
        for (int i = 0; i < EVENTS.size(); i++) {
            System.out.printf("event\t%s\t%d\t%s%n", label, i, EVENTS.get(i));
        }
        if (variants) {
            for (int i = 0; i < VARIANT_FILTERS.size(); i++) {
                System.out.printf("filter\t%s\t%d\t%d%n", label, i,
                        VARIANT_FILTERS.get(i).getFilteredCount());
            }
            System.out.printf("filters-built\t%s\t%d%n", label, VARIANT_FILTERS.size());
        } else {
            for (int i = 0; i < READ_FILTERS.size(); i++) {
                System.out.printf("filter\t%s\t%d\t%d%n", label, i,
                        READ_FILTERS.get(i).getFilteredCount());
            }
            System.out.printf("filters-built\t%s\t%d%n", label, READ_FILTERS.size());
        }
        System.out.printf("summary\t%s\t%s%n", label, summary);
    }
}

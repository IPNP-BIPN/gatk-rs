/*
 * `StrandArtifactFilter`'s M step and the Brent optimiser under it, taken from the reference.
 *
 * The other half of the last of the ten filters the `filter-mutect-calls` golden needs.
 * `learnParameters` re-estimates the strand-artifact prior and the beta shape between passes, and
 * the shape comes out of `OptimizationUtils.max`, which is commons-math's `BrentOptimizer`.
 *
 * Six behaviours this is built to catch.
 *
 *   - `GoalType.MAXIMIZE` IS IMPLEMENTED BY NEGATING INSIDE THE SEARCH, not by negating the
 *     objective: `fx = -fx` after every evaluation, and the pair reports `isMinim ? fx : -fx`. The
 *     returned VALUE is therefore the objective's, while every comparison in between is on its
 *     negation;
 *   - THE OPTIMISER RETURNS THE BEST POINT EVER SEEN, not the last one. `best` is folded over the
 *     previous and current pairs at every iteration, and the initial guess is the first candidate;
 *   - THE TOLERANCES ARE CHECKED AT CONSTRUCTION: a relative tolerance below `2 * ulp(1)` and an
 *     absolute tolerance of zero or less are refusals, not clamps;
 *   - THE EVALUATION BUDGET IS A REFUSAL TOO, and the filter's call site allows a hundred;
 *   - `potentialArtifacts` IS `getArtifactProbability() > 0.1` WHILE `totalNonArtifacts` IS OVER
 *     ALL the accumulated E steps, so the two sums are over different sets;
 *   - AND `eSteps.clear()` HAPPENS AT THE END OF `learnParameters` as well as in
 *     `clearAccumulatedData`, so the list is cleared twice per pass.
 *
 * Output:
 *
 *     opt\t<label>\t<min>,<max>,<guess>,<rel>,<abs>,<maxEval>=<point>,<value>
 *     learned\t<label>\t<strandArtifactPrior>,<alphaStrand>,<betaStrand>
 *     accumulated\t<label>\t<number of E steps>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: StrandArtifactMStepDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.ErrorProbabilities;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2Filter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.StrandArtifactFilter;
import org.broadinstitute.hellbender.utils.OptimizationUtils;

import java.io.File;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.function.DoubleUnaryOperator;

public class StrandArtifactMStepDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele SNV = Allele.create("C", false);

    public static void main(final String[] args) throws Exception {
        System.out.println("# StrandArtifactMStepDump: the M step and the Brent optimiser under it");

        // The filter's own call: a hundred evaluations, both tolerances 0.01, over [0.01, 100]
        // starting at 1.0.
        optimise("quadratic-at-three", x -> -(x - 3) * (x - 3), 0.01, 100, 1.0, 0.01, 0.01, 100);
        optimise("quadratic-at-fifty", x -> -(x - 50) * (x - 50), 0.01, 100, 1.0, 0.01, 0.01, 100);
        // The maximum sits on a bound.
        optimise("maximum-at-the-lower-bound", x -> -x, 0.01, 100, 1.0, 0.01, 0.01, 100);
        optimise("maximum-at-the-upper-bound", x -> x, 0.01, 100, 1.0, 0.01, 0.01, 100);
        // Nothing to find at all.
        optimise("flat", x -> 1.0, 0.01, 100, 1.0, 0.01, 0.01, 100);
        // Two maxima in the interval: which one is found depends on the trajectory.
        optimise("two-maxima", Math::sin, 0.01, 20, 1.0, 0.01, 0.01, 100);
        // The guess on each bound.
        optimise("guess-at-the-minimum", x -> -(x - 3) * (x - 3), 0.01, 100, 0.01, 0.01, 0.01, 100);
        optimise("guess-at-the-maximum", x -> -(x - 3) * (x - 3), 0.01, 100, 100, 0.01, 0.01, 100);
        // A tighter tolerance walks further.
        optimise("tight-tolerance", x -> -(x - 3) * (x - 3), 0.01, 100, 1.0, 1e-10, 1e-10, 1000);
        // The interval given the other way round, which the optimiser sorts.
        optimise("reversed-interval", x -> -(x - 3) * (x - 3), 100, 0.01, 1.0, 0.01, 0.01, 100);
        // The refusals.
        optimise("too-few-evaluations", x -> -(x - 3) * (x - 3), 0.01, 100, 1.0, 0.01, 0.01, 3);
        optimise("relative-tolerance-too-small", x -> -(x - 3) * (x - 3), 0.01, 100, 1.0, 1e-17, 0.01, 100);
        optimise("absolute-tolerance-zero", x -> -(x - 3) * (x - 3), 0.01, 100, 1.0, 0.01, 0.0, 100);
        optimise("guess-outside-the-interval", x -> -(x - 3) * (x - 3), 0.01, 100, 200, 0.01, 0.01, 100);
        // An objective that is not a number.
        optimise("nan-objective", x -> Double.NaN, 0.01, 100, 1.0, 0.01, 0.01, 100);

        // The M step itself, over sets of accumulated E steps.
        learn("no-data");
        learn("one-strong-artifact", "50,50|20,0");
        learn("one-weak-site", "50,50|10,10");
        learn("strong-and-weak", "50,50|20,0", "50,50|10,10");
        learn("two-strong-artifacts", "50,50|20,0", "50,50|0,30");
        learn("every-site-weak", "50,50|10,10", "50,50|9,11");
        learn("deep-artifact", "2000,2000|400,0");
        // Learning twice: the second pass has nothing left to learn from.
        learnTwice("learned-twice", "50,50|20,0");
    }

    static void optimise(final String label, final DoubleUnaryOperator objective, final double min,
                         final double max, final double guess, final double relative,
                         final double absolute, final int evaluations) {
        try {
            final var pair = OptimizationUtils.max(objective, min, max, guess, relative, absolute,
                    evaluations);
            System.out.printf("opt\t%s\t%s,%s,%s,%s,%s,%d=%s,%s%n", label, Double.toString(min),
                    Double.toString(max), Double.toString(guess), Double.toString(relative),
                    Double.toString(absolute), evaluations, Double.toString(pair.getPoint()),
                    Double.toString(pair.getValue()));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    static Mutect2FilteringEngine engine() {
        final Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        lines.add(new VCFHeaderLine("normal_sample", "N1"));
        final VCFHeader header = new VCFHeader(lines, List.of("T1", "N1"));
        return new Mutect2FilteringEngine(new M2FiltersArgumentCollection(), header,
                new File("no-such-stats-file.tsv"));
    }

    static VariantContext record(final String table) {
        return new VariantContextBuilder("dump", "chr1", 100, 100, List.of(REF, SNV))
                .attribute("AS_SB_TABLE", table).make();
    }

    /** Accumulate the given records, learn once, and print what was learned. */
    static void learn(final String label, final String... tables) throws Exception {
        final StrandArtifactFilter filter = new StrandArtifactFilter();
        accumulate(filter, tables);
        System.out.printf("accumulated\t%s\t%d%n", label, accumulatedCount(filter));
        invokeLearn(filter);
        System.out.printf("learned\t%s\t%s%n", label, learned(filter));
    }

    /** The same, learned twice: the second pass starts from an empty list. */
    static void learnTwice(final String label, final String... tables) throws Exception {
        final StrandArtifactFilter filter = new StrandArtifactFilter();
        accumulate(filter, tables);
        invokeLearn(filter);
        System.out.printf("learned\t%s-first\t%s%n", label, learned(filter));
        System.out.printf("accumulated\t%s-after\t%d%n", label, accumulatedCount(filter));
        invokeLearn(filter);
        System.out.printf("learned\t%s-second\t%s%n", label, learned(filter));
    }

    static void accumulate(final StrandArtifactFilter filter, final String... tables) throws Exception {
        final Method method = Mutect2Filter.class.getDeclaredMethod("accumulateDataForLearning",
                VariantContext.class, ErrorProbabilities.class, Mutect2FilteringEngine.class);
        method.setAccessible(true);
        for (final String table : tables) {
            final VariantContext vc = record(table);
            final Mutect2FilteringEngine engine = engine();
            method.invoke(filter, vc, new ErrorProbabilities(List.of(filter), vc, engine, null), engine);
        }
    }

    static void invokeLearn(final StrandArtifactFilter filter) throws Exception {
        final Method method =
                Mutect2Filter.class.getDeclaredMethod("learnParametersAndClearAccumulatedData");
        method.setAccessible(true);
        method.invoke(filter);
    }

    static int accumulatedCount(final StrandArtifactFilter filter) throws Exception {
        final Field field = StrandArtifactFilter.class.getDeclaredField("eSteps");
        field.setAccessible(true);
        return ((List<?>) field.get(filter)).size();
    }

    static String learned(final StrandArtifactFilter filter) throws Exception {
        return Double.toString(readDouble(filter, "strandArtifactPrior")) + ","
                + Double.toString(readDouble(filter, "alphaStrand")) + ","
                + Double.toString(readDouble(filter, "betaStrand"));
    }

    static double readDouble(final StrandArtifactFilter filter, final String name) throws Exception {
        final Field field = StrandArtifactFilter.class.getDeclaredField(name);
        field.setAccessible(true);
        return (double) field.get(filter);
    }
}

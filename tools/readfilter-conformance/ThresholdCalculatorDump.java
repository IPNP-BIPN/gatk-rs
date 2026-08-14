/*
 * ThresholdCalculator, taken from the reference.
 *
 * The piece of `FilterMutectCalls` that decides where the error probability is cut: a list of
 * per-variant posteriors goes in, and one threshold comes out, by one of three strategies. Seven
 * behaviours this is built to catch.
 *
 *   - THE LIST IS SORTED IN PLACE, so the caller's own list comes back reordered. Both strategies
 *     call `Collections.sort(posteriors)` on the list they were handed;
 *   - THE OPTIMAL F SCORE KEEPS THE LAST TIE, the comparison being `F >= optimalFScore` rather than
 *     `>`, so a run of equally good cut points resolves to the largest of them;
 *   - AND ITS ANSWER IS THREE-WAY: no index at all is `0`, the last index is `1`, and anything else
 *     is the posterior AT that index rather than between it and the next;
 *   - THE FALSE DISCOVERY RATE WALKS UNTIL THE RATE IS EXCEEDED and then steps BACK one, so the
 *     threshold is the last posterior that kept the cumulative rate acceptable, `0.0` if that was
 *     the very first one, and `1.0` if the rate was never exceeded at all;
 *   - CONSTANT LEAVES THE INITIAL THRESHOLD ALONE, whatever was accumulated;
 *   - RELEARNING CLEARS THE ACCUMULATED PROBABILITIES, so relearning twice computes the second
 *     threshold from an empty list — which the two strategies answer differently, `1.0` for the
 *     false discovery rate and `0.0` for the F score;
 *   - AND A NEGATIVE BETA OR RATE IS REFUSED by `ParamUtils`, with the message the caller passed.
 *
 * Output:
 *
 *     threshold\t<label>\t<the threshold, as Double.toString>
 *     sorted\t<label>\t<the caller's list after the call>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: ThresholdCalculatorDump
 */

import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.ThresholdCalculator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ThresholdCalculatorDump {

    public static void main(final String[] args) {
        System.out.println("# ThresholdCalculatorDump: one list of posteriors, one cut point");

        // A spread of posteriors, deliberately out of order.
        final List<Double> spread = Arrays.asList(0.9, 0.1, 0.5, 0.02, 0.3);
        // Every posterior the same, which is what makes the tie-breaking visible.
        final List<Double> tied = Arrays.asList(0.2, 0.2, 0.2, 0.2);
        // Everything nearly certain to be an error.
        final List<Double> hopeless = Arrays.asList(0.99, 0.98, 0.97);
        // Everything nearly certain to be real.
        final List<Double> clean = Arrays.asList(0.001, 0.002, 0.003);
        final List<Double> single = Arrays.asList(0.4);
        final List<Double> empty = new ArrayList<>();

        for (final ThresholdCalculator.Strategy strategy : ThresholdCalculator.Strategy.values()) {
            run(strategy + "-spread", strategy, spread, 0.05, 1.0);
            run(strategy + "-tied", strategy, tied, 0.05, 1.0);
            run(strategy + "-hopeless", strategy, hopeless, 0.05, 1.0);
            run(strategy + "-clean", strategy, clean, 0.05, 1.0);
            run(strategy + "-single", strategy, single, 0.05, 1.0);
            run(strategy + "-empty", strategy, empty, 0.05, 1.0);
        }

        // The same list under a rate loose enough to accept everything, and one tight enough to
        // reject the first posterior it sees.
        run("FALSE_DISCOVERY_RATE-loose", ThresholdCalculator.Strategy.FALSE_DISCOVERY_RATE,
                spread, 1.0, 1.0);
        run("FALSE_DISCOVERY_RATE-tight", ThresholdCalculator.Strategy.FALSE_DISCOVERY_RATE,
                spread, 0.001, 1.0);
        // Beta weighs recall against precision: zero is precision alone.
        run("OPTIMAL_F_SCORE-beta-zero", ThresholdCalculator.Strategy.OPTIMAL_F_SCORE,
                spread, 0.05, 0.0);
        run("OPTIMAL_F_SCORE-beta-ten", ThresholdCalculator.Strategy.OPTIMAL_F_SCORE,
                spread, 0.05, 10.0);

        // The list is sorted in place: the same list object, before and after.
        final List<Double> mutated = new ArrayList<>(Arrays.asList(0.9, 0.1, 0.5));
        System.out.printf("sorted\tbefore\t%s%n", mutated);
        ThresholdCalculator.calculateThresholdBasedOnFalseDiscoveryRate(mutated, 0.05);
        System.out.printf("sorted\tafter\t%s%n", mutated);

        // Relearning twice: the second time there is nothing left to learn from.
        twice("FALSE_DISCOVERY_RATE", ThresholdCalculator.Strategy.FALSE_DISCOVERY_RATE, spread);
        twice("OPTIMAL_F_SCORE", ThresholdCalculator.Strategy.OPTIMAL_F_SCORE, spread);
        twice("CONSTANT", ThresholdCalculator.Strategy.CONSTANT, spread);

        // The two refusals.
        refusal("negative-beta", ThresholdCalculator.Strategy.OPTIMAL_F_SCORE, spread, 0.05, -1.0);
        refusal("negative-rate", ThresholdCalculator.Strategy.FALSE_DISCOVERY_RATE, spread, -0.5, 1.0);
    }

    static void run(final String label, final ThresholdCalculator.Strategy strategy,
                    final List<Double> posteriors, final double maxFalseDiscoveryRate,
                    final double beta) {
        final ThresholdCalculator calculator =
                new ThresholdCalculator(strategy, 0.123, maxFalseDiscoveryRate, beta);
        // A copy, so that one run's in-place sort is not the next run's input.
        calculator.addCombinedErrorProbabilites(new ArrayList<>(posteriors));
        calculator.relearnThresholdAndClearAcumulatedProbabilities();
        System.out.printf("threshold\t%s\t%s%n", label, Double.toString(calculator.getThreshold()));
    }

    static void twice(final String label, final ThresholdCalculator.Strategy strategy,
                      final List<Double> posteriors) {
        final ThresholdCalculator calculator = new ThresholdCalculator(strategy, 0.123, 0.05, 1.0);
        calculator.addCombinedErrorProbabilites(new ArrayList<>(posteriors));
        calculator.relearnThresholdAndClearAcumulatedProbabilities();
        System.out.printf("threshold\t%s-first\t%s%n", label,
                Double.toString(calculator.getThreshold()));
        calculator.relearnThresholdAndClearAcumulatedProbabilities();
        System.out.printf("threshold\t%s-second\t%s%n", label,
                Double.toString(calculator.getThreshold()));
    }

    static void refusal(final String label, final ThresholdCalculator.Strategy strategy,
                        final List<Double> posteriors, final double maxFalseDiscoveryRate,
                        final double beta) {
        try {
            run(label, strategy, posteriors, maxFalseDiscoveryRate, beta);
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }
}

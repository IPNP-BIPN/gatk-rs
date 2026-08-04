/*
 * NaturalLogUtils, taken from the reference.
 *
 * The natural-log arithmetic every somatic likelihood goes through, and the layer G1.9 stands on.
 * It is dumped as RAW BIT PATTERNS, never as decimal: the whole question here is whether a port
 * built on a permissively-licensed `exp` lands on the same double, and decimal rendering discards
 * exactly the last-place difference that question is about.
 *
 * Three things this suite is designed to separate:
 *
 *   - THE EXACT PATH. logSumExp starts its accumulator at 1.0, because the maximum's own term is
 *     skipped in the loop and folded in as that 1. Then `sum != 1.0` skips the log entirely. So an
 *     array with one non-infinite entry, or one maximum with everything else at -Infinity, returns
 *     maxValue untouched: no exp, no log, exact by construction. The `one-*` and `*-neginf` cases
 *     are those, and a port that differs on them is wrong about the algorithm rather than about
 *     the exponential;
 *   - THE BOUNDED PATH. Everything else calls exp once per element. htsjdk-rs decision 0025 puts
 *     the port's exp within 1 ulp of Math.exp, so these rows measure how far that propagates
 *     through a sum and a log;
 *   - THE REFUSAL, which is on the ACCUMULATOR and not on the inputs. The exception fires after
 *     the loop, on `sum`, so a NaN input reaches it and a -Infinity input does not, though the
 *     message names the inputs either way.
 *
 * Output:
 *
 *     lse\t<label>\t<result bits>\t<inputs, bits, comma separated>
 *     lse\t<label>\tE:<class>:<message>\t<inputs>
 *     norm\t<label>\t<result bits, comma separated>
 *     post\t<label>\t<result bits, comma separated>
 *     l1me\t<label>\t<input bits>\t<result bits>
 *
 * Usage: NaturalLogUtilsDump
 */

import org.broadinstitute.hellbender.utils.NaturalLogUtils;

import java.util.Arrays;
import java.util.stream.Collectors;

public class NaturalLogUtilsDump {

    public static void main(final String[] args) {
        System.out.println("# NaturalLogUtilsDump: logSumExp, normalize, posteriors, log1mexp");

        // The exact path: no exp and no log is called on any of these.
        lse("one-zero", 0.0);
        lse("one-negative", -3.5);
        lse("one-large", 700.0);
        lse("max-and-neginf", -1.5, Double.NEGATIVE_INFINITY);
        lse("neginf-and-max", Double.NEGATIVE_INFINITY, -1.5);
        lse("all-neginf", Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
        lse("one-neginf", Double.NEGATIVE_INFINITY);

        // Ties: maxElementIndex takes the FIRST maximum, and which index is skipped decides which
        // term is folded in as the 1 rather than exponentiated. Same values, same answer only if
        // the port skips the same one.
        lse("tie-first", 2.0, 2.0);
        lse("tie-three", 2.0, 2.0, 2.0);
        lse("tie-with-smaller", 2.0, 1.0, 2.0);

        // The bounded path.
        lse("two-close", 0.0, -0.5);
        lse("two-far", 0.0, -50.0);
        lse("two-very-far", 0.0, -745.0);
        lse("three-equal-ish", -1.0, -1.0000001, -0.9999999);
        lse("phred-scale", -2.302585092994046, -4.605170185988091, -6.907755278982137);
        lse("many", -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1.0);
        lse("negatives-large", -700.0, -701.0, -702.0);
        lse("mixed-with-neginf", -1.0, Double.NEGATIVE_INFINITY, -2.0, Double.NEGATIVE_INFINITY);
        // Ordering: the same multiset, summed in two orders. Addition is not associative, so these
        // may differ from each other in the reference too.
        lse("order-ascending", -3.0, -2.0, -1.0);
        lse("order-descending", -1.0, -2.0, -3.0);

        // The refusal, which is on the accumulator.
        lse("nan-input", 0.0, Double.NaN);
        lse("posinf-input", 0.0, Double.POSITIVE_INFINITY);

        // normalizeFromLogToLinearSpace: every element goes through exp, so no exact path.
        norm("norm-uniform", -1.0, -1.0, -1.0);
        norm("norm-skewed", 0.0, -10.0, -20.0);
        norm("norm-with-neginf", 0.0, Double.NEGATIVE_INFINITY, -1.0);
        norm("norm-single", -4.0);
        norm("norm-phred", -2.302585092994046, -4.605170185988091);

        // posteriors: ebeAdd then normalize, which is the per-read call of the Dirichlet fixed
        // point and therefore the exact shape G1.9.2 will run in a loop.
        post("post-flat", new double[] {-1.0, -1.0}, new double[] {-2.0, -3.0});
        post("post-three", new double[] {-0.5, -1.0, -2.0}, new double[] {-1.0, -1.0, -1.0});
        post("post-with-neginf", new double[] {0.0, 0.0}, new double[] {0.0, Double.NEGATIVE_INFINITY});
        post("post-dirichlet-like",
                new double[] {-1.0986122886681098, -1.0986122886681098, -1.0986122886681098},
                new double[] {-0.6931471805599453, -2.3025850929940455, -6.907755278982137});

        // log1mexp, whose branch is the function.
        for (final double a : new double[] {
                -0.01, -0.1, -0.5, -0.6931471805599453, -0.7, -1.0, -5.0, -50.0, -745.0, 0.0, 1.0}) {
            l1me(a);
        }
    }

    static String bits(final double value) {
        return Long.toHexString(Double.doubleToRawLongBits(value));
    }

    static String bitList(final double[] values) {
        return Arrays.stream(values).mapToObj(NaturalLogUtilsDump::bits).collect(Collectors.joining(","));
    }

    static void lse(final String label, final double... values) {
        String result;
        try {
            result = bits(NaturalLogUtils.logSumExp(values));
        } catch (final Exception | AssertionError e) {
            result = "E:" + e.getClass().getName() + ":" + e.getMessage();
        }
        System.out.printf("lse\t%s\t%s\t%s%n", label, result, bitList(values));
    }

    static void norm(final String label, final double... values) {
        String result;
        try {
            result = bitList(NaturalLogUtils.normalizeFromLogToLinearSpace(values.clone()));
        } catch (final Exception | AssertionError e) {
            result = "E:" + e.getClass().getName() + ":" + e.getMessage();
        }
        System.out.printf("norm\t%s\t%s\t%s%n", label, result, bitList(values));
    }

    static void post(final String label, final double[] priors, final double[] likelihoods) {
        String result;
        try {
            result = bitList(NaturalLogUtils.posteriors(priors.clone(), likelihoods.clone()));
        } catch (final Exception | AssertionError e) {
            result = "E:" + e.getClass().getName() + ":" + e.getMessage();
        }
        System.out.printf("post\t%s\t%s\t%s|%s%n", label, result, bitList(priors), bitList(likelihoods));
    }

    static void l1me(final double a) {
        System.out.printf("l1me\t%s\t%s%n", bits(a), bits(NaturalLogUtils.log1mexp(a)));
    }
}

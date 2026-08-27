/*
 * OptimizationUtils' univariate maximiser, taken from the reference.
 *
 * Apache Commons' Brent optimiser as GATK calls it, which is what fits the beta shape of a
 * somatic panel of normals. It is golden-section search with parabolic interpolation, and where it
 * stops is decided by two tolerances rather than by the function.
 *
 * Seven behaviours this is built to catch.
 *
 *   - IT RETURNS THE BEST POINT IT EVALUATED, not an exactly converged one, so the answer carries
 *     the tolerances in its digits and a symmetric function comes back a hair off centre: a
 *     printed `-0.00000000000000` is a small NEGATIVE number rounded, not a negative zero;
 *   - A FUNCTION THAT IS EXACTLY A PARABOLA IS SOLVED IN ONE STEP by the interpolation, so the
 *     tolerances do not reach it at all: every setting returns the vertex exactly. The tolerances
 *     are only visible on a function the interpolation cannot fit;
 *   - THE SEARCH IS LOCAL TO ITS INTERVAL RATHER THAN TO ITS GUESS: a function with two maxima
 *     answers with the SAME one from a guess near either of them, because the first golden step is
 *     long enough to leave whichever peak the guess was near, and it is bracketing the interval
 *     around one peak that separates them;
 *   - A GUESS AT AN INTERVAL END is still accepted, and the search runs one-sided from it;
 *   - MAXIMISING IS MINIMISING THE NEGATIVE, which the optimiser does internally, so the value
 *     that comes back is the function's own;
 *   - THE EVALUATION BUDGET IS A REFUSAL, not a silent stop: exceeding it throws;
 *   - AND A GUESS OUTSIDE THE INTERVAL IS REFUSED before anything is evaluated.
 *
 * Output:
 *
 *     max\t<label>=<point>,<value>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: BrentOptimizerDump
 */

import org.broadinstitute.hellbender.utils.OptimizationUtils;

import java.util.function.DoubleUnaryOperator;

public class BrentOptimizerDump {

    static void max(final String label, final DoubleUnaryOperator function, final double min,
                    final double max, final double guess, final double relativeTolerance,
                    final double absoluteTolerance, final int maxEvaluations) {
        try {
            final org.apache.commons.math3.optim.univariate.UnivariatePointValuePair result =
                    OptimizationUtils.max(function, min, max, guess, relativeTolerance,
                            absoluteTolerance, maxEvaluations);
            System.out.printf("max\t%s=%.14f,%.14f%n", label, result.getPoint(),
                    result.getValue());
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(cause.getMessage())));
        }
    }

    public static void main(final String[] args) {
        System.out.println("# BrentOptimizerDump: where the univariate maximiser stops");

        // A parabola with its peak at 3, which every setting should find.
        final DoubleUnaryOperator parabola = x -> -(x - 3.0) * (x - 3.0) + 10.0;
        max("parabola-default", parabola, 0.0, 10.0, 1.0, 0.001, 0.001, 1000);
        // The same peak under tighter tolerances, which stops somewhere else.
        max("parabola-tight", parabola, 0.0, 10.0, 1.0, 1e-12, 1e-12, 1000);
        max("parabola-loose", parabola, 0.0, 10.0, 1.0, 0.1, 0.1, 1000);
        // The same function from the other side of the peak.
        max("parabola-guess-high", parabola, 0.0, 10.0, 9.0, 0.001, 0.001, 1000);
        // A guess sitting exactly on the peak.
        max("parabola-guess-exact", parabola, 0.0, 10.0, 3.0, 0.001, 0.001, 1000);
        // A guess at an end of the interval.
        max("parabola-guess-at-min", parabola, 0.0, 10.0, 0.0, 0.001, 0.001, 1000);
        max("parabola-guess-at-max", parabola, 0.0, 10.0, 10.0, 0.001, 0.001, 1000);

        // A peak the interpolation cannot fit, so the tolerances decide where it stops. It is
        // piecewise LINEAR rather than logarithmic on purpose: Math.log and Math.exp are not
        // transcribable between a JVM and a libm under decision 0014, and a golden built on one
        // would be pinning the platform's transcendentals rather than the optimiser.
        final DoubleUnaryOperator kinked = s -> -Math.abs(s - 7.0) - s / 50.0;
        max("kinked-settings", kinked, 0.01, 100.0, 1.0, 0.01, 0.1, 100);
        max("kinked-tight", kinked, 0.01, 100.0, 1.0, 1e-12, 1e-12, 1000);
        max("kinked-loose", kinked, 0.01, 100.0, 1.0, 0.5, 0.5, 1000);

        // Two maxima, where the guess decides which one is found. A quartic rather than a sine,
        // for the same reason as above: it peaks at 4 and at 16, both at zero.
        final DoubleUnaryOperator twoPeaks =
                x -> -(x - 4.0) * (x - 4.0) * (x - 16.0) * (x - 16.0);
        max("two-peaks-low-guess", twoPeaks, 0.0, 20.0, 1.0, 0.001, 0.001, 1000);
        max("two-peaks-high-guess", twoPeaks, 0.0, 20.0, 14.0, 0.001, 0.001, 1000);
        // The same function bracketed around each peak in turn, which is what actually separates
        // them: over the whole interval the first golden step is long enough to leave whichever
        // peak the guess was near.
        max("two-peaks-low-interval", twoPeaks, 0.0, 10.0, 1.0, 0.001, 0.001, 1000);
        max("two-peaks-high-interval", twoPeaks, 10.0, 20.0, 14.0, 0.001, 0.001, 1000);

        // A symmetric function, whose answer is symmetric too.
        max("symmetric", x -> -x * x, -5.0, 5.0, 2.0, 0.001, 0.001, 1000);

        // A budget too small to converge.
        max("too-few-evaluations", parabola, 0.0, 10.0, 1.0, 1e-14, 1e-14, 5);
        // A guess outside the interval.
        max("guess-below-min", parabola, 0.0, 10.0, -1.0, 0.001, 0.001, 1000);
        max("guess-above-max", parabola, 0.0, 10.0, 11.0, 0.001, 0.001, 1000);
        // An interval whose ends are the wrong way round.
        max("inverted-interval", parabola, 10.0, 0.0, 5.0, 0.001, 0.001, 1000);
    }
}

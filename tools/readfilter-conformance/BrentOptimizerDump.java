/*
 * OptimizationUtils' univariate maximiser, taken from the reference.
 *
 * Apache Commons' Brent optimiser as GATK calls it, which is what fits the beta shape of a
 * somatic panel of normals. It is golden-section search with parabolic interpolation, and where it
 * stops is decided by two tolerances rather than by the function.
 *
 * Seven behaviours this is built to catch.
 *
 *   - IT RETURNS THE POINT IT LAST EVALUATED, not an exactly converged one, so the answer carries
 *     the tolerances in its digits;
 *   - A FUNCTION THAT IS EXACTLY A PARABOLA IS SOLVED IN ONE STEP by the interpolation, so the
 *     tolerances do not reach it at all: every setting returns the vertex exactly. The tolerances
 *     are only visible on a function the interpolation cannot fit;
 *   - THE INITIAL GUESS MOVES THE ANSWER for a function with more than one maximum, because the
 *     search is local;
 *   - A GUESS AT AN INTERVAL END is still accepted, and the search runs one-sided from it;
 *   - MAXIMISING IS MINIMISING THE NEGATIVE, so a symmetric function gives a symmetric answer;
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

        // The settings the panel of normals uses, on a function shaped like its likelihood. This
        // one the interpolation cannot fit, so the tolerances decide where it stops.
        final DoubleUnaryOperator panelLike = s -> -Math.abs(Math.log(s) - Math.log(7.0)) - s / 50.0;
        max("panel-settings", panelLike, 0.01, 100.0, 1.0, 0.01, 0.1, 100);
        max("panel-tight", panelLike, 0.01, 100.0, 1.0, 1e-12, 1e-12, 1000);
        max("panel-loose", panelLike, 0.01, 100.0, 1.0, 0.5, 0.5, 1000);

        // Two maxima, where the guess decides which one is found.
        final DoubleUnaryOperator twoPeaks = x -> Math.sin(x);
        max("two-peaks-low-guess", twoPeaks, 0.0, 20.0, 1.0, 0.001, 0.001, 1000);
        max("two-peaks-high-guess", twoPeaks, 0.0, 20.0, 14.0, 0.001, 0.001, 1000);

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

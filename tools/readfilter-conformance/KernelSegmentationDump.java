/*
 * `SingularValueDecomposition` on the matrices `KernelSegmenter` builds, and the changepoints it
 * turns them into, taken from the reference.
 *
 * `CalculateContamination` is blocked on this decomposition and nothing else. What has to be decided
 * before a line is ported is whether byte identity is required THROUGH the decomposition or only
 * through what the segmenter does with it: the changepoints are indices, and a decomposition that
 * differs in the last bits may still give the same ones. This measures both ends so the decision has
 * numbers under it.
 *
 * Four things it is built to catch.
 *
 *   - THE SUBSAMPLE IS SEEDED, so the matrix is deterministic: `new Random(1216)` through
 *     `RandomGeneratorFactory`, and `rng.nextInt(data.size())` per point. A port that subsampled
 *     differently would decompose a different matrix and the comparison would say nothing;
 *   - THE SINGULAR VALUES ARE USED AS `1 / (sqrt(s) + EPSILON)`, so a small one is amplified: two
 *     decompositions agreeing to the last bit in the large values can disagree visibly in the
 *     reduced observation matrix through the small ones;
 *   - `U`'S SIGN CONVENTION REACHES THE OUTPUT. The reduced matrix multiplies by `U` column by
 *     column, so flipping a column's sign flips that column of the observations. Whether the
 *     changepoints survive it is the question, and the golden holds `U` entry by entry so a port can
 *     be judged on the matrix as well as on the indices;
 *   - AND THE CHANGEPOINTS ARE INDICES, so they can agree while the arithmetic under them does not.
 *
 * The dimension is kept small on purpose: a six-point approximation gives thirty-six `U` entries
 * rather than ten thousand, and the question is which bits differ rather than how many there are.
 *
 * Output:
 *
 *     singular\t<label>\t<index>=<bits>,<decimal>
 *     u\t<label>\t<row>,<column>=<bits>
 *     changepoints\t<label>\t<comma-separated indices>
 *
 * Usage: KernelSegmentationDump
 */

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenter;

import java.util.ArrayList;
import java.util.List;
import java.util.function.BiFunction;

public class KernelSegmentationDump {

    /** The kernel `ContaminationSegmenter` uses: a Gaussian of the difference. */
    static BiFunction<Double, Double, Double> gaussian(final double variance) {
        return (x, y) -> Math.exp(-(x - y) * (x - y) / (2.0 * variance));
    }

    public static void main(final String[] args) {
        System.out.println("# KernelSegmentationDump: the decomposition CalculateContamination waits on");

        // A step function with two changes, which is what a segmenter is for.
        final List<Double> twoSteps = new ArrayList<>();
        for (int i = 0; i < 90; i++) {
            twoSteps.add(i < 30 ? 0.1 : (i < 60 ? 0.5 : 0.2));
        }
        // A flat series, which has no changepoint to find.
        final List<Double> flat = new ArrayList<>();
        for (int i = 0; i < 90; i++) {
            flat.add(0.3);
        }
        // A ramp, where every point is a little different from the last.
        final List<Double> ramp = new ArrayList<>();
        for (int i = 0; i < 90; i++) {
            ramp.add(i / 90.0);
        }

        decompose("two-steps", twoSteps);
        decompose("flat", flat);
        decompose("ramp", ramp);

        changepoints("two-steps", twoSteps);
        changepoints("flat", flat);
        changepoints("ramp", ramp);
        // The same series with more changepoints allowed and a tighter penalty.
        changepoints("two-steps-lenient", twoSteps, 10, 0.0, 0.0);
        changepoints("two-steps-strict", twoSteps, 10, 10.0, 10.0);
    }

    /** The kernel matrix the segmenter builds, decomposed the way it decomposes it. */
    static void decompose(final String label, final List<Double> data) {
        final RealMatrix matrix = subKernelMatrix(data, 6, 0.01);
        final SingularValueDecomposition svd = new SingularValueDecomposition(matrix);
        final double[] values = svd.getSingularValues();
        for (int i = 0; i < values.length; i++) {
            System.out.printf("singular\t%s\t%d=%016x,%s%n", label, i,
                    Double.doubleToRawLongBits(values[i]), Double.toString(values[i]));
        }
        final RealMatrix u = svd.getU();
        for (int i = 0; i < u.getRowDimension(); i++) {
            for (int j = 0; j < u.getColumnDimension(); j++) {
                System.out.printf("u\t%s\t%d,%d=%016x%n", label, i, j,
                        Double.doubleToRawLongBits(u.getEntry(i, j)));
            }
        }
    }

    /**
     * The subsampled kernel matrix, built exactly as `calculateReducedObservationMatrix` builds it:
     * the same seed, the same generator, the same order of calls.
     */
    static RealMatrix subKernelMatrix(final List<Double> data, final int dimension,
                                      final double variance) {
        final org.apache.commons.math3.random.RandomGenerator rng =
                org.apache.commons.math3.random.RandomGeneratorFactory
                        .createRandomGenerator(new java.util.Random(1216));
        final int numSubsample = Math.min(dimension, data.size());
        final List<Double> subsample = new ArrayList<>();
        if (numSubsample == data.size()) {
            subsample.addAll(data);
        } else {
            for (int i = 0; i < numSubsample; i++) {
                subsample.add(data.get(rng.nextInt(data.size())));
            }
        }
        final BiFunction<Double, Double, Double> kernel = gaussian(variance);
        final RealMatrix matrix = new Array2DRowRealMatrix(numSubsample, numSubsample);
        for (int i = 0; i < numSubsample; i++) {
            for (int j = 0; j < i; j++) {
                final double value = kernel.apply(subsample.get(i), subsample.get(j));
                matrix.setEntry(i, j, value);
                matrix.setEntry(j, i, value);
            }
            matrix.setEntry(i, i, kernel.apply(subsample.get(i), subsample.get(i)));
        }
        return matrix;
    }

    static void changepoints(final String label, final List<Double> data) {
        changepoints(label, data, 10, 1.0, 1.0);
    }

    static void changepoints(final String label, final List<Double> data, final int maximum,
                             final double linear, final double logLinear) {
        final List<Integer> found = new KernelSegmenter<>(data).findChangepoints(maximum,
                gaussian(0.01), 6, List.of(8, 16), linear, logLinear,
                KernelSegmenter.ChangepointSortOrder.INDEX);
        final StringBuilder text = new StringBuilder();
        for (final int index : found) {
            if (text.length() > 0) {
                text.append(',');
            }
            text.append(index);
        }
        System.out.printf("changepoints\t%s\t%s%n", label,
                text.length() == 0 ? "(none)" : text.toString());
    }
}

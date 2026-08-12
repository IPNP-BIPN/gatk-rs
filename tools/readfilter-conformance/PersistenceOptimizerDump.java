/*
 * PersistenceOptimizer's local minima, taken from the reference.
 *
 * The watershed under the kernel segmenter: given one-dimensional data it returns every local
 * minimum, ordered by topological persistence, which is what CalculateContamination's segmenter
 * uses to pick changepoint candidates.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE ORDER IS BY DECREASING PERSISTENCE and the GLOBAL MINIMUM IS PREPENDED, so the first
 *     index is not part of the sorted run and the first persistence is the whole range, maximum
 *     minus minimum, rather than any pair's;
 *   - THE INDICES ARE SORTED BY Comparator.comparingDouble, WHICH IS Double.compare, so -0.0 SORTS
 *     BELOW 0.0 and a NaN SORTS ABOVE EVERYTHING. A port comparing with `<` would order neither the
 *     same way, and the watershed starts from whichever point that ordering calls lowest;
 *   - THE SORT IS STABLE, so equal values keep their index order and a plateau is entered from the
 *     left;
 *   - A PLATEAU'S LEFTMOST POINT IS THE LOCAL MINIMUM when the plateau follows a maximum or opens
 *     the data, which is the documented consequence of that stability;
 *   - A LOCAL MAXIMUM TAKES THE COLOUR OF THE POINT ON ITS LEFT, read AFTER the two components have
 *     been merged, so the merge order decides the colour;
 *   - MERGING KEEPS THE COMPONENT WITH THE LOWER MINIMUM, and on a tie the LOWER COLOUR, which is
 *     the earlier component and not the earlier index;
 *   - A SINGLE POINT HAS NO PAIRS AT ALL, so the answer is one index and a persistence of zero;
 *   - AND EMPTY DATA IS REFUSED, with a message from the argument check and not from the algorithm.
 *
 * Output:
 *
 *     data\t<label>\t<the values, comma separated>
 *     minima\t<label>\t<the indices, comma separated>
 *     persistence\t<label>\t<the persistences, comma separated>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: PersistenceOptimizerDump
 */

import org.broadinstitute.hellbender.tools.copynumber.utils.optimization.PersistenceOptimizer;

import java.util.ArrayList;
import java.util.List;

public class PersistenceOptimizerDump {

    public static void main(final String[] args) {
        System.out.println("# PersistenceOptimizerDump: local minima by topological persistence, from the reference");

        run("single", new double[] {1.0});
        run("two-equal", new double[] {1.0, 1.0});
        run("increasing", new double[] {0.0, 1.0, 2.0, 3.0, 4.0});
        run("decreasing", new double[] {4.0, 3.0, 2.0, 1.0, 0.0});
        // One deep valley and one shallow one: the deep one is first, the global minimum before it.
        run("two-valleys", new double[] {5.0, 1.0, 4.0, 3.0, 4.5, 0.0, 6.0});
        // A plateau opening the data, and one after a maximum.
        run("plateau-at-start", new double[] {2.0, 2.0, 2.0, 5.0, 1.0});
        run("plateau-after-maximum", new double[] {1.0, 5.0, 3.0, 3.0, 3.0, 6.0});
        run("constant", new double[] {2.0, 2.0, 2.0, 2.0});
        // Zeroes of both signs, which Double.compare separates and `<` does not.
        run("signed-zeroes", new double[] {0.0, -0.0, 0.0, -0.0});
        run("signed-zeroes-with-hill", new double[] {-0.0, 1.0, 0.0, 1.0, -0.0});
        // A NaN, which sorts above everything, and the infinities.
        run("with-nan", new double[] {1.0, Double.NaN, 0.0, 2.0, 0.5});
        run("with-infinities", new double[] {Double.POSITIVE_INFINITY, 1.0, Double.NEGATIVE_INFINITY, 1.0, 2.0});
        // Something long enough for several merges, computed rather than written out.
        final double[] wiggly = new double[41];
        for (int i = 0; i < wiggly.length; i++) {
            wiggly[i] = Math.sin(i / 3.0) + 0.25 * Math.sin(i * 1.7) + i / 40.0;
        }
        run("wiggly", wiggly);
        // Ties between whole valleys, so the merge's tie-break decides.
        run("twin-valleys", new double[] {3.0, 0.0, 3.0, 0.0, 3.0});
        run("twin-valleys-uneven", new double[] {3.0, 0.0, 2.0, 0.0, 4.0});

        // And the refusal.
        runRefused("empty", new double[] {});
    }

    static void run(final String label, final double[] data) {
        final List<String> values = new ArrayList<>();
        for (final double value : data) {
            values.add(String.valueOf(value));
        }
        System.out.printf("data\t%s\t%s%n", label, String.join(",", values));

        final PersistenceOptimizer optimizer = new PersistenceOptimizer(data);
        final List<String> indices = new ArrayList<>();
        for (final int index : optimizer.getMinimaIndices()) {
            indices.add(String.valueOf(index));
        }
        System.out.printf("minima\t%s\t%s%n", label, String.join(",", indices));

        final List<String> persistences = new ArrayList<>();
        for (final double persistence : optimizer.getPersistences()) {
            persistences.add(String.valueOf(persistence));
        }
        System.out.printf("persistence\t%s\t%s%n", label, String.join(",", persistences));
    }

    static void runRefused(final String label, final double[] data) {
        try {
            new PersistenceOptimizer(data);
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(), e.getMessage());
            return;
        }
        System.out.printf("error\t%s\tnone%n", label);
    }
}

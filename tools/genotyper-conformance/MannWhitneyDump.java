/*
 * MannWhitneyU, taken from the reference.
 *
 * The rank-sum test behind BaseQRankSum, MQRankSum, ReadPosRankSum and ClippingRankSum. Four of
 * its decisions are not what "a Mann-Whitney U test" suggests:
 *
 *   - the RANKS ARE FLOATS. The values are doubles, the ranks and the sums of ranks are single
 *     precision, and the averaged rank of a tie band is a float division;
 *   - two different tests, chosen by size: with either series at 10 or more it is the normal
 *     approximation, otherwise every permutation of the group labels is enumerated. So Z comes
 *     from a cumulative probability in one regime and from an INVERSE cumulative probability in
 *     the other, and those two are not inverses of each other in commons-math3;
 *   - transformTies reports ZERO ties when every value is tied, because the sigma formula breaks
 *     down there, and the continuity correction is then dropped as well, so two mechanisms
 *     conspire to produce exactly p = 0.5;
 *   - the permutation p-value is half the observed bin plus everything more extreme, not the
 *     cumulative distribution, "which gives a p-value of 1 in the most extreme case".
 *
 * Output:
 *
 *     mwu\t<label>\t<side>\t<u bits>\t<z bits>\t<p bits>\t<median shift bits>
 *     z\t<u>\t<n1>\t<n2>\t<nties>\t<side>\t<z bits>
 *
 * The ranks themselves are not emitted: Rank is a private nested class, so nothing outside the
 * package can name it. They are measured through U, which is a float sum of them.
 *
 * Usage: MannWhitneyDump
 */

import org.broadinstitute.hellbender.utils.MannWhitneyU;

import java.util.Arrays;

public class MannWhitneyDump {

    public static void main(final String[] args) {
        System.out.println("# MannWhitneyDump: the rank-sum test under the RankSum annotations");

        final MannWhitneyU mwu = new MannWhitneyU();

        // Ordinary shapes, both below and above the size at which the test changes.
        emit(mwu, "tiny-separated", new double[] {1, 2, 3}, new double[] {4, 5, 6});
        emit(mwu, "tiny-overlapping", new double[] {1, 3, 5}, new double[] {2, 4, 6});
        emit(mwu, "tiny-reversed", new double[] {4, 5, 6}, new double[] {1, 2, 3});
        emit(mwu, "one-each", new double[] {1}, new double[] {2});
        emit(mwu, "one-vs-many", new double[] {1}, new double[] {2, 3, 4, 5});
        emit(mwu, "empty-first", new double[] {}, new double[] {1, 2});
        emit(mwu, "empty-second", new double[] {1, 2}, new double[] {});

        // The boundary of the two regimes: 9 and 9 is exact, 10 and 9 is normal.
        emit(mwu, "nine-and-nine", ramp(1, 9), ramp(2, 9));
        emit(mwu, "ten-and-nine", ramp(1, 10), ramp(2, 9));
        emit(mwu, "nine-and-ten", ramp(1, 9), ramp(2, 10));
        emit(mwu, "ten-and-ten", ramp(1, 10), ramp(2, 10));

        // Ties, including the case where every single value is tied.
        emit(mwu, "all-tied", constant(5, 12), constant(5, 12));
        emit(mwu, "all-tied-small", constant(5, 4), constant(5, 4));
        emit(mwu, "half-tied", new double[] {1, 1, 2, 2, 3}, new double[] {1, 1, 2, 2, 3});
        emit(mwu, "tied-across", ramp(1, 12), ramp(1, 12));
        emit(mwu, "one-tie-band", new double[] {1, 2, 2, 2, 3}, new double[] {4, 5, 6, 7, 8});

        // Large enough that a float rank sum stops being exact.
        emit(mwu, "large-ramp", ramp(1, 300), ramp(301, 300));
        emit(mwu, "large-interleaved", evens(600), odds(600));
        emit(mwu, "large-tied", constant(30, 200), constant(30, 200));

        // Quality-like values, which is what the annotations actually pass.
        emit(mwu, "qualities", new double[] {30, 30, 31, 32, 30, 29, 30, 30, 30, 31, 30, 30},
                new double[] {30, 28, 30, 27, 30, 30, 26, 30, 30, 30, 30, 25});
        emit(mwu, "mapping-qualities", constant(60, 15),
                new double[] {60, 60, 60, 59, 60, 60, 60, 60, 57, 60, 60, 60, 60, 60, 60});

        // Negative and non-integral values, which the annotations do not produce and the class
        // accepts.
        emit(mwu, "negative", new double[] {-3, -2, -1}, new double[] {-6, -5, -4});
        emit(mwu, "fractional", new double[] {0.5, 1.5, 2.5}, new double[] {1.0, 2.0, 3.0});

        // calculateZ on its own, so the continuity correction and the tie term are visible
        // without the rest of the test around them.
        for (final double u : new double[] {0, 1, 10, 50, 100, 4.5}) {
            for (final int[] sizes : new int[][] {{10, 10}, {12, 8}, {100, 3}}) {
                for (final double nties : new double[] {0, 6, 120}) {
                    for (final MannWhitneyU.TestType side : MannWhitneyU.TestType.values()) {
                        System.out.printf("z\t%d\t%d\t%d\t%d\t%s\t%d%n",
                                Double.doubleToRawLongBits(u), sizes[0], sizes[1],
                                Double.doubleToRawLongBits(nties), side,
                                Double.doubleToRawLongBits(
                                        mwu.calculateZ(u, sizes[0], sizes[1], nties, side)));
                    }
                }
            }
        }
    }

    static void emit(final MannWhitneyU mwu, final String label, final double[] series1,
                     final double[] series2) {
        for (final MannWhitneyU.TestType side : MannWhitneyU.TestType.values()) {
            // test() sorts its inputs in place, so each call gets its own copy.
            final MannWhitneyU.Result result =
                    mwu.test(series1.clone(), series2.clone(), side);
            System.out.printf("mwu\t%s\t%s\t%d\t%d\t%d\t%d%n", label, side,
                    Double.doubleToRawLongBits(result.getU()),
                    Double.doubleToRawLongBits(result.getZ()),
                    Double.doubleToRawLongBits(result.getP()),
                    Double.doubleToRawLongBits(result.getMedianShift()));
        }
    }

    static double[] ramp(final int from, final int count) {
        final double[] values = new double[count];
        for (int i = 0; i < count; i++) {
            values[i] = from + i;
        }
        return values;
    }

    static double[] constant(final double value, final int count) {
        final double[] values = new double[count];
        Arrays.fill(values, value);
        return values;
    }

    static double[] evens(final int count) {
        final double[] values = new double[count];
        for (int i = 0; i < count; i++) {
            values[i] = 2 * i;
        }
        return values;
    }

    static double[] odds(final int count) {
        final double[] values = new double[count];
        for (int i = 0; i < count; i++) {
            values[i] = 2 * i + 1;
        }
        return values;
    }
}

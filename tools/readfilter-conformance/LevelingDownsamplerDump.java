/*
 * What a LevelingDownsampler keeps, and the permutation underneath it, taken from the reference.
 *
 * LevelingDownsampler is what --max-depth-per-sample reaches when several stacks have to be cut
 * down to one total. It removes items by asking MathUtils.sampleIndicesWithoutReplacement, which
 * is RandomDataGenerator.nextPermutation over the static Well19937c: a different generator from
 * the java.util.Random that ReservoirDownsampler uses, with a different algorithm and a different
 * stream position.
 *
 * Both layers are dumped here because both decide the answer.
 *
 * nextPermutation(n, k):
 *
 *   - it shuffles ALL n entries and returns the first k, so it costs n - 1 draws whatever k is.
 *     Sampling 2 of 1000 spends 999 draws, and every later consumer of the stream moves with them.
 *     The probe reports the generator's next value after each call, which is what measures that;
 *   - the shuffle runs downwards, from n - 1 to 0, with the target drawn from [0, i], so the
 *     bounds shrink. Upwards would be a valid Fisher-Yates and a different permutation;
 *   - the last step takes no draw: at i == 0 the target is 0 without a call. That is why the cost
 *     is n - 1 and not n;
 *   - the returned indices are the shuffled head, so they are not sorted. The caller only uses
 *     them as membership, so sorting them would keep the same items: the raw list is dumped so a
 *     port that sorted them still fails.
 *
 * LevelingDownsampler:
 *
 *   - the plan is arithmetic. A round-robin walk over the sizes decrements one at a time, so how
 *     many go from each stack is fixed before any draw happens, and only which ones is random;
 *   - the walk stops on consecutive refusals rather than on a scan, so a stack that is at the
 *     minimum and then is not is reconsidered;
 *   - a stack that keeps everything takes no draw at all, so stack sizes change the stream
 *     position for the stacks after them;
 *   - LinkedList and ArrayList take different code paths and must keep the same items. Both are
 *     dumped, so a divergence between them is a failure rather than an assumption;
 *   - with minElementsPerStack 0 a stack can be planned down to zero items, and nextPermutation
 *     then throws NotStrictlyPositiveException rather than returning nothing.
 *
 * Utils.resetRandomGenerator() is called before each case, because both static generators are one
 * stream each and a case inheriting the previous one's position would measure the order of the
 * cases rather than the code under test.
 *
 * Output:
 *
 *     perm\t<n>\t<k>\t<comma-separated indices>\t<next int from the Well19937c>
 *     permerror\t<n>\t<k>\t<class>
 *     level\t<label>\t<list kind>\t<groups, stacks separated by |, items by comma>
 *     levelstats\t<label>\t<list kind>\t<size>\t<discarded>\t<next int from the Well19937c>
 *     levelerror\t<label>\t<class>
 *
 * Usage: LevelingDownsamplerDump
 */

import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.downsampling.LevelingDownsampler;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.StringJoiner;

public class LevelingDownsamplerDump {

    public static void main(final String[] args) {
        System.out.println("# LevelingDownsamplerDump: leveling, and the permutation underneath it");

        // The permutation on its own, at sizes either side of the interesting cases: k == n (a
        // full shuffle), k == 1 (still n - 1 draws), and n == 1 (no draw at all).
        for (final int[] nk : new int[][] {
                {1, 1}, {2, 1}, {2, 2}, {3, 1}, {3, 2}, {3, 3},
                {5, 1}, {5, 3}, {5, 5}, {10, 1}, {10, 5}, {10, 10},
                {16, 4}, {17, 4}, {50, 7}, {100, 99}, {1000, 2}}) {
            permutation(nk[0], nk[1]);
        }
        // The two refusals.
        permutationError(3, 4);
        permutationError(3, 0);
        permutationError(3, -1);

        // Leveling. Sizes are given per stack; the target is what the sum may not exceed.
        level("under-target", new int[] {3, 3, 3}, 20, 1);
        level("exactly-target", new int[] {3, 3, 3}, 9, 1);
        level("one-over", new int[] {3, 3, 3}, 8, 1);
        level("even-cut", new int[] {10, 10, 10}, 15, 1);
        // Uneven stacks, where the round-robin plan and the minimum interact.
        level("uneven", new int[] {1, 5, 20}, 10, 1);
        level("floor-blocks", new int[] {1, 1, 20}, 5, 1);
        // A minimum high enough that the walk gives up before reaching the target.
        level("minimum-blocks", new int[] {4, 4, 4}, 3, 3);
        // A single stack, and an empty one among others.
        level("one-stack", new int[] {25}, 4, 1);
        level("empty-among-others", new int[] {0, 6, 6}, 5, 1);
        // No stacks at all, which never enters the plan.
        level("no-stacks", new int[] {}, 5, 1);
        // Target zero with minimum one: the stacks can only fall to one item each.
        level("target-zero-min-one", new int[] {4, 4}, 0, 1);
        // Target zero with minimum zero, which plans a stack down to nothing and throws.
        level("target-zero-min-zero", new int[] {4, 4}, 0, 0);
    }

    static void permutation(final int n, final int k) {
        Utils.resetRandomGenerator();
        final int[] indices = MathUtils.sampleIndicesWithoutReplacement(n, k);
        final StringJoiner values = new StringJoiner(",");
        for (final int index : indices) {
            values.add(Integer.toString(index));
        }
        // Where the shared Well19937c ended up, which is what shows the shuffle cost n - 1 draws
        // rather than k or n.
        System.out.printf("perm\t%d\t%d\t%s\t%d%n",
                n, k, values, Utils.getRandomDataGenerator().getRandomGenerator().nextInt());
    }

    static void permutationError(final int n, final int k) {
        Utils.resetRandomGenerator();
        try {
            MathUtils.sampleIndicesWithoutReplacement(n, k);
            System.out.printf("permerror\t%d\t%d\t%s%n", n, k, "none");
        } catch (final Exception e) {
            System.out.printf("permerror\t%d\t%d\t%s%n", n, k, e.getClass().getName());
        }
    }

    /** Both list kinds, because they take different removal paths and must agree. */
    static void level(final String label, final int[] sizes, final long target, final int minimum) {
        levelWith(label, "linked", sizes, target, minimum, true);
        levelWith(label, "array", sizes, target, minimum, false);
    }

    static void levelWith(final String label, final String kind, final int[] sizes,
                          final long target, final int minimum, final boolean linked) {
        Utils.resetRandomGenerator();

        final LevelingDownsampler<List<String>, String> downsampler =
                new LevelingDownsampler<>(target, minimum);

        int next = 0;
        for (final int size : sizes) {
            final List<String> stack = linked ? new LinkedList<>() : new ArrayList<>();
            for (int i = 0; i < size; i++) {
                stack.add(String.format("s%02d", next++));
            }
            downsampler.submit(stack);
        }

        try {
            downsampler.signalEndOfInput();
        } catch (final Exception e) {
            System.out.printf("levelerror\t%s\t%s\t%s%n", label, kind, e.getClass().getName());
            return;
        }

        final int discarded = downsampler.getNumberOfDiscardedItems();
        final int size = downsampler.size();
        final List<List<String>> groups = downsampler.consumeFinalizedItems();
        final StringJoiner stacks = new StringJoiner("|");
        for (final List<String> group : groups) {
            stacks.add(String.join(",", group));
        }
        System.out.printf("level\t%s\t%s\t%s%n", label, kind, stacks);
        System.out.printf("levelstats\t%s\t%s\t%d\t%d\t%d%n", label, kind, size, discarded,
                Utils.getRandomDataGenerator().getRandomGenerator().nextInt());
    }
}

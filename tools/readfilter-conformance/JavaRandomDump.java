/*
 * The exact sequences java.util.Random produces, taken from the reference.
 *
 * GATK's Utils.getRandomGenerator() is one static `new Random(47382911L)`, seeded at class
 * initialisation, and every draw in a run comes from that single stream. A consumer that takes one
 * value too few or too many does not merely shuffle its own result differently: it shifts the
 * stream for every later consumer as well. So the sequence has to be exact before anything is
 * allowed to draw from it, and this dump exists to establish that before the downsampler is
 * ported.
 *
 * java.util.Random is the case where a port from the *specification* is legitimate: its Javadoc
 * states the algorithm and the constants, and requires every method to produce that exact
 * sequence. That is the opposite of HashMap's iteration order, which is documented as unspecified;
 * see docs/an-unspecified-order-that-reaches-the-output.md.
 *
 * Two methods are probed for their shape rather than their values:
 *
 *   - nextInt(bound) takes a different path for a power-of-two bound, so the same seed yields a
 *     different value at bound 16 than the surrounding bounds would suggest;
 *   - the general path rejects biased values and redraws, so the number of draws it consumes
 *     depends on the values themselves. The sequences below are long enough to contain rejections.
 *
 * Output:
 *
 *     seq\t<seed>\t<method>\t<comma-separated values>
 *     interleave\t<seed>\t<comma-separated values, methods mixed>
 *
 * Usage: JavaRandomDump
 */

import java.util.Random;
import java.util.StringJoiner;

public class JavaRandomDump {

    /** The seed GATK pins in Utils, plus edges: zero, negative, and one that scrambles to zero. */
    static final long[] SEEDS = {47382911L, 0L, 1L, -1L, Long.MAX_VALUE, Long.MIN_VALUE,
            0x5DEECE66DL};

    static final int COUNT = 24;

    public static void main(final String[] args) {
        System.out.println("# JavaRandomDump: the exact sequences java.util.Random produces");

        for (final long seed : SEEDS) {
            ints(seed);
            longs(seed);
            booleans(seed);
            doubles(seed);
            floats(seed);
            // Bounds either side of a power of two, and one large enough that rejections happen.
            for (final int bound : new int[] {2, 3, 7, 15, 16, 17, 100, 1000, Integer.MAX_VALUE}) {
                bounded(seed, bound);
            }
            interleaved(seed);
        }
    }

    static void ints(final long seed) {
        final Random random = new Random(seed);
        final StringJoiner values = new StringJoiner(",");
        for (int i = 0; i < COUNT; i++) {
            values.add(Integer.toString(random.nextInt()));
        }
        System.out.printf("seq\t%d\tnextInt\t%s%n", seed, values);
    }

    static void longs(final long seed) {
        final Random random = new Random(seed);
        final StringJoiner values = new StringJoiner(",");
        for (int i = 0; i < COUNT; i++) {
            values.add(Long.toString(random.nextLong()));
        }
        System.out.printf("seq\t%d\tnextLong\t%s%n", seed, values);
    }

    static void booleans(final long seed) {
        final Random random = new Random(seed);
        final StringBuilder values = new StringBuilder();
        for (int i = 0; i < COUNT; i++) {
            values.append(random.nextBoolean() ? '1' : '0');
        }
        System.out.printf("seq\t%d\tnextBoolean\t%s%n", seed, values);
    }

    static void doubles(final long seed) {
        final Random random = new Random(seed);
        final StringJoiner values = new StringJoiner(",");
        for (int i = 0; i < COUNT; i++) {
            // The raw bits, not the printed decimal: a double compared as text would compare
            // Double.toString rather than the value.
            values.add(Long.toString(Double.doubleToRawLongBits(random.nextDouble())));
        }
        System.out.printf("seq\t%d\tnextDouble\t%s%n", seed, values);
    }

    static void floats(final long seed) {
        final Random random = new Random(seed);
        final StringJoiner values = new StringJoiner(",");
        for (int i = 0; i < COUNT; i++) {
            values.add(Integer.toString(Float.floatToRawIntBits(random.nextFloat())));
        }
        System.out.printf("seq\t%d\tnextFloat\t%s%n", seed, values);
    }

    static void bounded(final long seed, final int bound) {
        final Random random = new Random(seed);
        final StringJoiner values = new StringJoiner(",");
        for (int i = 0; i < COUNT; i++) {
            values.add(Integer.toString(random.nextInt(bound)));
        }
        System.out.printf("seq\t%d\tnextInt(%d)\t%s%n", seed, bound, values);
    }

    /**
     * One generator, methods mixed. This is what a real consumer does, and it catches a port whose
     * individual methods are right but whose draw *counts* are not: nextDouble takes two draws and
     * nextLong takes two, so a port taking one would agree here only by accident.
     */
    static void interleaved(final long seed) {
        final Random random = new Random(seed);
        final StringJoiner values = new StringJoiner(",");
        for (int i = 0; i < COUNT; i++) {
            switch (i % 5) {
                case 0 -> values.add("i" + random.nextInt());
                case 1 -> values.add("d" + Double.doubleToRawLongBits(random.nextDouble()));
                case 2 -> values.add("b" + (random.nextBoolean() ? 1 : 0));
                case 3 -> values.add("l" + random.nextLong());
                default -> values.add("n" + random.nextInt(37));
            }
        }
        System.out.printf("interleave\t%d\t%s%n", seed, values);
    }
}

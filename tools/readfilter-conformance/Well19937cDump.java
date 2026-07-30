/*
 * The exact sequences a Well19937c produces, taken from the reference.
 *
 * GATK carries two static generators, not one:
 *
 *     private static final Random randomGenerator = new Random(GATK_RANDOM_SEED);
 *     private static final RandomDataGenerator randomDataGenerator =
 *             new RandomDataGenerator(new Well19937c(GATK_RANDOM_SEED));
 *
 * They are different algorithms with different streams, and Utils.resetRandomGenerator() resets
 * both. JavaRandomDump established the first; this establishes the second, before anything is
 * allowed to draw from it. Every draw a consumer takes moves the position for every later
 * consumer, so the sequence has to be exact and so does the number of draws each method costs.
 *
 * Four shapes are probed because a port gets them wrong while looking right:
 *
 *   - the two constructors disagree. new Well19937c(long) splits the seed into TWO ints, high word
 *     first, so it seeds a different pool from new Well19937c((int) seed) even when the long fits
 *     in an int: the seed array's length decides how much of the 624-word pool the scrambler
 *     fills. GATK uses the long one. Both are dumped, and they must differ;
 *   - the derived methods are not java.util.Random's. nextDouble is
 *     ((next(26) << 26) | next(26)) * 2^-52 where java.util.Random uses 27 bits and 2^-53;
 *     nextFloat takes 23 bits, not 24; nextLong ors the halves where java.util.Random adds them
 *     signed;
 *   - nextInt(bound) IS the same shape, because commons-math copied it from Apache Harmony. The
 *     power-of-two path and the rejection loop are dumped either side of 16 to catch that;
 *   - nextBytes takes one extra full draw after its four-at-a-time loop, even when the length is a
 *     multiple of four and nothing is left to fill. Lengths 0 through 9 are dumped, which is what
 *     makes that visible.
 *
 * The raw bits of doubles and floats are printed rather than the decimals, so the comparison is of
 * the value and not of Double.toString.
 *
 * Output:
 *
 *     seq\t<constructor>\t<seed>\t<method>\t<comma-separated values>
 *     bytes\t<constructor>\t<seed>\t<length>\t<hex>
 *     interleave\t<constructor>\t<seed>\t<comma-separated values, methods mixed>
 *
 * Usage: Well19937cDump
 */

import org.apache.commons.math3.random.Well19937c;

import java.util.StringJoiner;
import java.util.function.Supplier;

public class Well19937cDump {

    /** The seed GATK pins, plus edges: zero, one, negative, and the extremes of long. */
    static final long[] LONG_SEEDS = {47382911L, 0L, 1L, -1L, Long.MAX_VALUE, Long.MIN_VALUE};

    /** The same values as int seeds, which must NOT produce the same streams. */
    static final int[] INT_SEEDS = {47382911, 0, 1, -1, Integer.MIN_VALUE, Integer.MAX_VALUE};

    /** A seed array longer than the two words the long constructor makes, to reach the copy path. */
    static final int[][] ARRAY_SEEDS = {{1, 2, 3}, {0, 47382911}, {-1, -1, -1, -1}};

    static final int COUNT = 24;

    static final int[] BOUNDS = {2, 3, 7, 15, 16, 17, 100, 1000, Integer.MAX_VALUE};

    public static void main(final String[] args) {
        System.out.println("# Well19937cDump: the exact sequences commons-math3's Well19937c produces");

        for (final long seed : LONG_SEEDS) {
            all("long", Long.toString(seed), () -> new Well19937c(seed));
        }
        for (final int seed : INT_SEEDS) {
            all("int", Integer.toString(seed), () -> new Well19937c(seed));
        }
        for (final int[] seed : ARRAY_SEEDS) {
            final StringJoiner label = new StringJoiner(",");
            for (final int word : seed) {
                label.add(Integer.toString(word));
            }
            all("ints", label.toString(), () -> new Well19937c(seed));
        }
    }

    /** Every method, from a freshly seeded generator each time. */
    static void all(final String ctor, final String seed, final Supplier<Well19937c> make) {
        ints(ctor, seed, make.get());
        longs(ctor, seed, make.get());
        booleans(ctor, seed, make.get());
        doubles(ctor, seed, make.get());
        floats(ctor, seed, make.get());
        for (final int bound : BOUNDS) {
            bounded(ctor, seed, bound, make.get());
        }
        // nextLong(bound), whose rejection loop consumes two draws per attempt.
        for (final long bound : new long[] {7L, 1L << 40, Long.MAX_VALUE}) {
            boundedLong(ctor, seed, bound, make.get());
        }
        for (int length = 0; length <= 9; length++) {
            bytes(ctor, seed, length, make.get());
        }
        interleaved(ctor, seed, make.get());
    }

    static void ints(final String ctor, final String seed, final Well19937c random) {
        final StringJoiner values = new StringJoiner(",");
        for (int i = 0; i < COUNT; i++) {
            values.add(Integer.toString(random.nextInt()));
        }
        System.out.printf("seq\t%s\t%s\tnextInt\t%s%n", ctor, seed, values);
    }

    static void longs(final String ctor, final String seed, final Well19937c random) {
        final StringJoiner values = new StringJoiner(",");
        for (int i = 0; i < COUNT; i++) {
            values.add(Long.toString(random.nextLong()));
        }
        System.out.printf("seq\t%s\t%s\tnextLong\t%s%n", ctor, seed, values);
    }

    static void booleans(final String ctor, final String seed, final Well19937c random) {
        final StringBuilder values = new StringBuilder();
        for (int i = 0; i < COUNT; i++) {
            values.append(random.nextBoolean() ? '1' : '0');
        }
        System.out.printf("seq\t%s\t%s\tnextBoolean\t%s%n", ctor, seed, values);
    }

    static void doubles(final String ctor, final String seed, final Well19937c random) {
        final StringJoiner values = new StringJoiner(",");
        for (int i = 0; i < COUNT; i++) {
            values.add(Long.toString(Double.doubleToRawLongBits(random.nextDouble())));
        }
        System.out.printf("seq\t%s\t%s\tnextDouble\t%s%n", ctor, seed, values);
    }

    static void floats(final String ctor, final String seed, final Well19937c random) {
        final StringJoiner values = new StringJoiner(",");
        for (int i = 0; i < COUNT; i++) {
            values.add(Integer.toString(Float.floatToRawIntBits(random.nextFloat())));
        }
        System.out.printf("seq\t%s\t%s\tnextFloat\t%s%n", ctor, seed, values);
    }

    static void bounded(final String ctor, final String seed, final int bound,
                        final Well19937c random) {
        final StringJoiner values = new StringJoiner(",");
        for (int i = 0; i < COUNT; i++) {
            values.add(Integer.toString(random.nextInt(bound)));
        }
        System.out.printf("seq\t%s\t%s\tnextInt(%d)\t%s%n", ctor, seed, bound, values);
    }

    static void boundedLong(final String ctor, final String seed, final long bound,
                            final Well19937c random) {
        final StringJoiner values = new StringJoiner(",");
        for (int i = 0; i < COUNT; i++) {
            values.add(Long.toString(random.nextLong(bound)));
        }
        System.out.printf("seq\t%s\t%s\tnextLong(%d)\t%s%n", ctor, seed, bound, values);
    }

    /**
     * nextBytes at one length, then one more int from the same generator. The trailing int is the
     * point: it reports where the stream ended up, so a port that took the wrong number of draws
     * fails even at a length where the bytes themselves happen to agree.
     */
    static void bytes(final String ctor, final String seed, final int length,
                      final Well19937c random) {
        final byte[] buffer = new byte[length];
        random.nextBytes(buffer);
        final StringBuilder hex = new StringBuilder();
        for (final byte b : buffer) {
            hex.append(String.format("%02x", b));
        }
        System.out.printf("bytes\t%s\t%s\t%d\t%s\t%d%n",
                ctor, seed, length, hex, random.nextInt());
    }

    /**
     * One generator, methods mixed. This catches a port whose individual methods are right but
     * whose draw counts are not: nextDouble takes two draws and nextLong takes two, so a port
     * taking one would agree with every single-method row and fail here.
     */
    static void interleaved(final String ctor, final String seed, final Well19937c random) {
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
        System.out.printf("interleave\t%s\t%s\t%s%n", ctor, seed, values);
    }
}

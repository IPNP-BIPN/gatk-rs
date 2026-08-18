/*
 * `Double.toString`, over a corpus wide enough to say WHICH values the port disagrees on.
 *
 * `gatk_engine::tsv_table::java_double_to_string` defers to Rust's shortest-round-trip formatter,
 * and the pre-JDK19 `FloatingDecimal` does not always emit the shortest digits. One value was known
 * (`4.9E-324`, from the `string-format` golden); this establishes the rest rather than assuming
 * there are none.
 *
 * The corpus is deterministic on both sides: sixty named values, then a thousand bit patterns from
 * splitmix64 seeded at 1, with the non-finite ones skipped. Nothing here is random.
 *
 * Output:
 *
 *     tostring\t<the bits, as sixteen hex digits>\t<Double.toString>
 *
 * Usage: DoubleToStringDump
 */

public class DoubleToStringDump {

    public static void main(final String[] args) {
        System.out.println("# DoubleToStringDump: Double.toString over a deterministic corpus");

        // The named values: the ends of the range, the powers of two, and the neighbours of one.
        show(0.0);
        show(-0.0);
        show(Double.MIN_VALUE);
        show(-Double.MIN_VALUE);
        show(Double.MIN_NORMAL);
        show(Double.MAX_VALUE);
        show(1.0);
        show(-1.0);
        show(Math.nextUp(1.0));
        show(Math.nextDown(1.0));
        show(0.1);
        show(0.2);
        show(0.3);
        show(1.0 / 3.0);
        show(2.0 / 3.0);
        show(Math.PI);
        show(Math.E);
        show(1e-3);
        show(1e7);
        show(Math.nextDown(1e7));
        show(1e21);
        show(1e22);
        show(1e23);
        show(2.675);
        show(1.005);
        show(1.2345);
        show(0.0625);
        // Powers of two across the whole exponent range.
        for (int exponent = -1074; exponent <= 1023; exponent += 67) {
            show(Math.scalb(1.0, exponent));
        }

        // A thousand bit patterns, from splitmix64 seeded at one.
        long state = 1L;
        int emitted = 0;
        while (emitted < 1000) {
            state += 0x9E3779B97F4A7C15L;
            long z = state;
            z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
            z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
            z = z ^ (z >>> 31);
            final double value = Double.longBitsToDouble(z);
            if (Double.isNaN(value) || Double.isInfinite(value)) {
                continue;
            }
            show(value);
            emitted++;
        }
    }

    static void show(final double value) {
        System.out.printf("tostring\t%016x\t%s%n", Double.doubleToRawLongBits(value),
                Double.toString(value));
    }
}

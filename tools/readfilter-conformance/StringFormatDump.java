/*
 * `String.format("%.Nf", x)`, taken from the reference.
 *
 * Not the exact decimal expansion of the double. `java.util.Formatter` takes the digits
 * `Double.toString` would produce -- the shortest that round-trip -- and pads or rounds THOSE to the
 * requested scale. Four consequences this is built to catch.
 *
 *   - A VALUE WHOSE SHORTEST FORM ENDS IN 5 ROUNDS UP EVEN WHEN ITS EXACT EXPANSION DOES NOT.
 *     `2.675` is `2.67499999999999982...` exactly, and Java prints `2.68` where a formatter working
 *     from the expansion prints `2.67`. C's printf prints `2.67`;
 *   - A LARGE VALUE IS PADDED WITH ZEROS, NOT EXPANDED. `1e300` is a one and three hundred zeros,
 *     not the 1000000000000000052504760255204420248704... the double actually is;
 *   - THE ROUNDING IS HALF_UP ON THOSE DIGITS, not half-even, so `0.0625` at three places is
 *     `0.063` where Rust's own formatter answers `0.062`;
 *   - AND THE THREE NON-FINITE SPELLINGS HAVE NO SIGN ON NaN.
 *
 * Output:
 *
 *     format\t<label>\t<places>=<formatted>
 *     tostring\t<label>\t<Double.toString>
 *
 * Usage: StringFormatDump
 */

import java.util.List;

public class StringFormatDump {

    public static void main(final String[] args) {
        System.out.println("# StringFormatDump: String.format(\"%.Nf\") over the shortest digits");

        // The ties, where the shortest digits and the exact expansion disagree.
        value("two-point-six-seven-five", 2.675, List.of(0, 1, 2, 3, 6));
        value("one-sixteenth", 0.0625, List.of(2, 3, 4, 6));
        value("eight-thousandths", 0.008, List.of(1, 2, 3));
        value("one-point-zero-zero-five", 1.005, List.of(2, 3));
        // Magnitudes, where the padding shows.
        value("one-e-twenty-one", 1e21, List.of(0, 3));
        value("one-e-three-hundred", 1e300, List.of(3));
        value("max-value", Double.MAX_VALUE, List.of(3));
        value("min-normal", Double.MIN_NORMAL, List.of(3, 320));
        value("min-value", Double.MIN_VALUE, List.of(3, 330));
        // The ordinary cases two ported tools reach.
        value("two-sevenths-percent", 100.0 * 2.0 / 7.0, List.of(2));
        value("ten-sixty-fifths-percent", 100.0 * 10.0 / 65.0, List.of(2));
        value("all-of-it", 100.0 * 7.0 / 7.0, List.of(2));
        // Carries that run the whole way up.
        value("nine-nine-nine", 9.999, List.of(2));
        value("zero-nine-nine-nine", 0.999, List.of(2));
        // Zero and its sign.
        value("zero", 0.0, List.of(0, 2));
        value("negative-zero", -0.0, List.of(2));
        value("negative", -2.675, List.of(2));
        // The non-finite spellings.
        value("nan", Double.NaN, List.of(2));
        value("infinity", Double.POSITIVE_INFINITY, List.of(2));
        value("negative-infinity", Double.NEGATIVE_INFINITY, List.of(2));
    }

    static void value(final String label, final double value, final List<Integer> places) {
        System.out.printf("tostring\t%s\t%s%n", label, Double.toString(value));
        for (final int scale : places) {
            System.out.printf("format\t%s\t%d=%s%n", label, scale,
                    String.format("%." + scale + "f", value));
        }
    }
}

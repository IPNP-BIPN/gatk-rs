/*
 * java.text.DecimalFormat, for the two patterns AllelePseudoDepth emits through.
 *
 * The annotation writes strings, not numbers, so the formatter is part of the output and not an
 * afterthought. It is also the piece with the least documentation behind it: the Javadoc specifies
 * HALF_EVEN and the pattern language, and says nothing about the three facts below, each of which
 * changes a printed digit.
 *
 *   - IT ROUNDS THE SHORTEST DECIMAL FORM, NOT THE VALUE. Formatting 0.1 to forty places gives
 *     "0.1", not the double's true 0.1000000000000000055511151231257827.
 *   - THE TIE IS BROKEN BY THE TRUE VALUE. 0.155 rounds DOWN and 0.165 rounds UP, and no rule
 *     reading only the digit string separates them: the doubles sit on opposite sides of the
 *     halfway point. Only when the shortest form IS the value, as for 0.125, is it a real tie and
 *     does HALF_EVEN look at the neighbour.
 *   - THE TWO PATTERNS DISAGREE WITH EACH OTHER AT ONE PLACE. When the rounding position falls
 *     before the first digit, "#.##" follows the value and "#.####" rounds to even. 0.005 and 5e-5
 *     are the same shape one decade apart and go opposite ways. That is where DecimalFormat's
 *     internal fast path stops applying, so the pattern decides the rule.
 *
 * Values are dumped as BIT PATTERNS with their two formatted strings. A decimal literal in a
 * golden would be re-parsed by the reader and could name a different double than the one measured,
 * which for a suite about the last digit is the whole question.
 *
 * Output:
 *
 *     f\t<bits>\t<"#.##">\t<"#.####">
 *
 * Usage: DecimalFormatDump
 */

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

public class DecimalFormatDump {

    static final DecimalFormat TWO = new DecimalFormat("#.##");
    static final DecimalFormat FOUR = new DecimalFormat("#.####");

    public static void main(final String[] args) {
        System.out.println("# DecimalFormatDump: #.## and #.####, the two patterns AllelePseudoDepth uses");

        final List<Double> corpus = new ArrayList<>();

        // Zero, its sign, and the non-finite symbols.
        add(corpus, 0.0, -0.0, Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY);

        // Apparent ties whose doubles fall on opposite sides of the halfway point. These are the
        // rows a digit-string rule gets wrong, and it gets about half of them wrong.
        add(corpus, 0.145, 0.155, 0.165, 0.175, 0.185, 0.195, 0.115, 0.135);
        add(corpus, 1.005, 1.015, 1.025, 1.035, 1.045, 2.675, 8.835, 0.0625);
        add(corpus, 0.12345, 0.12355, 0.12365, 0.99995, 0.99985, 0.00125, 0.00135);

        // Real ties: exactly representable, so HALF_EVEN goes to the even neighbour and the
        // neighbour's parity is what decides.
        add(corpus, 0.125, 0.375, 0.625, 0.875, 1.125, 2.125, 3.375, 0.25, 0.75, 1.5, 2.5, 3.5);
        add(corpus, 0.015625, 0.046875, 0.078125, 0.0078125, 0.00390625);

        // The underflow boundary, where the two patterns part company. Both neighbours of each are
        // included, because the rule is about the exact double and not about the decade.
        for (int k = 1; k <= 8; k++) {
            final double at = Double.parseDouble("5e-" + k);
            add(corpus, at, Math.nextUp(at), Math.nextDown(at));
            add(corpus, Double.parseDouble("1.5e-" + k), Double.parseDouble("2.5e-" + k));
        }

        // The shapes a pseudo-depth and a pseudo-fraction actually take: counts near integers and
        // fractions in [0, 1] that need every one of the four digits.
        add(corpus, 1.0, 2.0, 5.0, 10.0, 25.0, 100.0, 1000.0, 0.5, 0.3333333333333333);
        add(corpus, 0.6666666666666666, 0.14285714285714285, 0.9999, 0.99999, 0.0001, 0.00005);
        add(corpus, 12345.6789, 1234567.891, 3.141592653589793, 2.718281828459045);
        for (int i = 1; i < 40; i++) {
            add(corpus, i / 3.0, i / 7.0, i / 1024.0, i * 0.005, i * 1e-5);
        }

        // Negatives, including ones that round away to nothing and keep their sign.
        add(corpus, -0.5, -0.155, -0.125, -1e-9, -5e-5, -12345.6789, -0.0001);

        // Past the domain the port claims: sixteen significant digits, and above 2^53. Both are
        // Java's digit generation rather than its rounding, and both are recorded so the limit is
        // measured rather than asserted.
        add(corpus, 6.985838094673373e14, 5.936134122025243e14, -2.0097114521676883e17);
        add(corpus, 9.007199254740992e15, 1e16, 1e21, 1e23, Double.MAX_VALUE, Double.MIN_VALUE);

        for (final double value : corpus) {
            System.out.printf("f\t%s\t%s\t%s%n",
                    Long.toHexString(Double.doubleToRawLongBits(value)),
                    TWO.format(value), FOUR.format(value));
        }
    }

    static void add(final List<Double> corpus, final double... values) {
        for (final double value : values) {
            corpus.add(value);
        }
    }
}

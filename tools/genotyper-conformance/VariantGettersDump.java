/*
 * How an INFO attribute becomes a double[], taken from the reference.
 *
 * Everything that reads a numeric annotation off a VariantContext goes through
 * VariantContextGetters.attributeValueToDoubleArray, and the function is held together by an
 * accident:
 *
 *     } else if (value.getClass().isAssignableFrom(Iterable.class)) {
 *
 * The test is BACKWARDS. isAssignableFrom asks whether the argument's type can be assigned to the
 * receiver's, so this asks whether an Iterable fits in an ArrayList variable, which is false. Every
 * list therefore falls through to the branch labelled "as a last resort", which renders the value
 * with String.valueOf and takes the result apart:
 *
 *     Stream.of(String.valueOf(value).trim().replaceAll("\\[|\\]", "").split(","))
 *
 * String.valueOf(List.of(1.0, 2.0)) is "[1.0, 2.0]", so the second field is " 2.0" WITH A LEADING
 * SPACE, and it parses only because Double.parseDouble trims. A port that split on ", ", or used a
 * parser refusing leading whitespace, agrees on a one-element list and diverges on every longer one.
 *
 * Two more things the signature hides: a field equal to "." becomes the caller's missingValue,
 * which is -1 for the overload every annotator uses, so a missing TLOD arrives as a number nobody
 * can recognise as missing; and a field that is neither becomes a GATKException naming the
 * annotation rather than a NumberFormatException.
 *
 * The doubles travel as raw bits, because getTumorLogOdds multiplies them by Math.log(10).
 *
 * Output:
 *
 *     array\t<label>\t<comma-separated raw bits|null|E:class:message>
 *     tlod\t<label>\t<comma-separated raw bits|null|E:class:message>
 *     maxindex\t<label>\t<MathUtils.maxElementIndex>
 *
 * Usage: VariantGettersDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

public class VariantGettersDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele ALT1 = Allele.create("C", false);
    static final Allele ALT2 = Allele.create("G", false);

    public static void main(final String[] args) {
        System.out.println("# VariantGettersDump: an INFO attribute as a double[]");

        // Absent: the null default of the two-argument overload.
        array("absent", null);

        // Strings, which is what a parsed VCF holds.
        array("one-string", "1.5");
        array("two-strings", "1.5,2.5");
        // With the space the list rendering would have produced, to show the parse survives it.
        array("string-with-spaces", "1.5, 2.5");
        array("string-with-brackets", "[1.5, 2.5]");
        array("negative", "-1.5,-2.5");
        array("scientific", "1.5e3,1.5E-3");
        // The missing value, alone and among others.
        array("all-missing", ".");
        array("one-missing", "1.5,.,2.5");
        // Not a double.
        array("not-a-double", "abc");
        array("empty-string", "");
        array("trailing-comma", "1.5,2.5,");
        // The spellings Double.parseDouble takes and Rust's parser does not, and the reverse.
        array("type-suffix", "1.5f,2.5d");
        array("hexadecimal", "0x1p3");
        array("spelled-infinity", "Infinity,-Infinity");
        array("lower-case-inf", "inf");

        // A List attribute, which is the branch the backwards test makes unreachable.
        array("list-of-doubles", List.of(1.5, 2.5));
        array("list-of-strings", List.of("1.5", "2.5"));
        array("list-of-one", List.of(1.5));
        array("list-with-missing", List.of("1.5", ".", "2.5"));
        // A boxed number, which renders through its own toString.
        array("boxed-double", 1.5);
        array("boxed-integer", 3);

        // getTumorLogOdds, which is the getter followed by a log10-to-log conversion in place.
        tlod("absent", null);
        tlod("one", "1.5");
        tlod("two", "1.5,2.5");
        tlod("missing", ".");
        tlod("list", List.of(1.5, 2.5));

        // maxElementIndex, which is where getTumorLogOdds's answer is used.
        maxIndex("one", new double[] {1.0});
        maxIndex("ascending", new double[] {1.0, 2.0, 3.0});
        maxIndex("descending", new double[] {3.0, 2.0, 1.0});
        maxIndex("tie", new double[] {2.0, 2.0, 1.0});
        maxIndex("tie-at-the-end", new double[] {1.0, 2.0, 2.0});
        maxIndex("negatives", new double[] {-3.0, -1.0, -2.0});
        maxIndex("with-nan", new double[] {Double.NaN, 1.0});
        maxIndex("nan-last", new double[] {1.0, Double.NaN});
        maxIndex("all-nan", new double[] {Double.NaN, Double.NaN});
        maxIndex("with-infinity", new double[] {1.0, Double.POSITIVE_INFINITY, 2.0});
    }

    static VariantContext contextWith(final String key, final Object value) {
        final VariantContextBuilder builder = new VariantContextBuilder()
                .chr("chr1").start(100).stop(100).alleles(List.of(REF, ALT1, ALT2));
        if (value != null) {
            final Map<String, Object> attributes = new LinkedHashMap<>();
            attributes.put(key, value);
            builder.attributes(attributes);
        }
        return builder.make();
    }

    static void array(final String label, final Object value) {
        final VariantContext vc = contextWith("X", value);
        try {
            final double[] result = VariantContextGetters.getAttributeAsDoubleArray(vc, "X");
            System.out.printf("array\t%s\t%s%n", label, bits(result));
        } catch (final Exception | AssertionError e) {
            System.out.printf("array\t%s\tE:%s:%s%n", label, e.getClass().getName(),
                    e.getMessage() == null ? "" : e.getMessage().replace('\n', ' '));
        }
    }

    static void tlod(final String label, final Object value) {
        final VariantContext vc = contextWith(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, value);
        try {
            System.out.printf("tlod\t%s\t%s%n", label,
                    bits(Mutect2FilteringEngine.getTumorLogOdds(vc)));
        } catch (final Exception | AssertionError e) {
            System.out.printf("tlod\t%s\tE:%s:%s%n", label, e.getClass().getName(),
                    e.getMessage() == null ? "" : e.getMessage().replace('\n', ' '));
        }
    }

    static void maxIndex(final String label, final double[] values) {
        try {
            System.out.printf("maxindex\t%s\t%d%n", label, MathUtils.maxElementIndex(values));
        } catch (final Exception | AssertionError e) {
            System.out.printf("maxindex\t%s\tE:%s:%s%n", label, e.getClass().getName(),
                    e.getMessage() == null ? "" : e.getMessage().replace('\n', ' '));
        }
    }

    static String bits(final double[] values) {
        if (values == null) {
            return "null";
        }
        final StringJoiner joiner = new StringJoiner(",");
        for (final double value : values) {
            joiner.add(Long.toString(Double.doubleToRawLongBits(value)));
        }
        return joiner.toString();
    }

    static List<Object> unused() {
        return Arrays.asList();
    }
}

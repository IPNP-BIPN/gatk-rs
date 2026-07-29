/*
 * What JexlExpressionReadTagValueFilter decides, taken from the reference.
 *
 * This is the 56th read filter and the only one that runs an expression language. Its answers are
 * not really the filter's: they are commons-jexl 2.1.1's, reached through one detail of the filter
 * that changes everything downstream. The context is
 *
 *     public Object get(final String name) { return read.getAttributeAsString(name); }
 *
 * so every attribute arrives as a **String**, whatever it is in the BAM. Which branch of
 * JexlArithmetic.compare then runs is decided by the *literal's* type rather than the tag's:
 *
 *   - `NM > 3`: 3 is an Integer, isNumberable is true, both sides go through toLong, and
 *     toLong of a non-numeric string throws rather than answering false;
 *   - `NM > 3.5`: 3.5 is a **Float**, because ASTNumberLiteral.setReal defaults to Float, so
 *     isFloatingPoint is true and both sides go through toDouble instead;
 *   - `NM > '3'`: neither side is a number, so it is a lexicographic String comparison, and
 *     "30" > "3" is true for a different reason than 30 > 3 is.
 *
 * The filter's own test is `!v.equals(Boolean.TRUE)`, so a null result is a NullPointerException
 * rather than a false, and a non-boolean result is a quiet false.
 *
 * The matrix is expressions x reads, because a coercion's outcome depends on both.
 *
 * Output:
 *
 *     eval\t<expression>\t<read>\t<ok:value|E:class>
 *     filter\t<expression>\t<read>\t<true|false|E:class>
 *     attr\t<read>\t<tag>=<the string the context hands over>
 *
 * Usage: JexlFilterDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.VariantContextUtils;
import org.apache.commons.jexl2.Expression;
import org.apache.commons.jexl2.JexlContext;
import org.broadinstitute.hellbender.engine.filters.JexlExpressionReadTagValueFilter;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class JexlFilterDump {

    static final int CONTIG_LENGTH = 200;

    /** The tags every probed read carries, in the order they are reported. */
    static final String[] TAGS = {"NM", "AS", "RA", "FL", "EM", "NG", "BI"};

    /**
     * The reads, by label. Each is a map of tag to value, and the value's Java type is chosen so
     * that getAttributeAsString has something different to render each time.
     */
    static Map<String, Map<String, Object>> reads() {
        final Map<String, Map<String, Object>> reads = new LinkedHashMap<>();
        // An ordinary read: small integers, a string that looks like a number.
        reads.put("plain", Map.of("NM", 3, "AS", 100, "RA", "30", "FL", 2.5f, "EM", "", "NG", 0));
        // The one that separates integral from lexicographic comparison: "30" > "3" but 30 < 300.
        reads.put("big", Map.of("NM", 300, "AS", 5, "RA", "300", "FL", 12.0f, "EM", "", "NG", 0));
        // Non-numeric tag values, which make every numeric coercion throw.
        reads.put("text", Map.of("NM", "abc", "AS", "x", "RA", "hello", "FL", "y", "EM", "",
                "NG", "0"));
        // A tag whose value has a decimal point, so isFloatingPointNumber is true for a *String*.
        reads.put("decimal", Map.of("NM", "3.0", "AS", "2.5", "RA", "1e3", "FL", 1.5f, "EM", "",
                "NG", "0"));
        // A read missing NM entirely: the context returns null, which every operator refuses.
        reads.put("missing", Map.of("AS", 7, "RA", "7", "FL", 1.0f, "EM", "", "NG", 0));
        // Boolean-looking strings, since toBoolean of a String is `"true".equals(s)`.
        reads.put("boolish", Map.of("NM", "true", "AS", "false", "RA", "TRUE", "FL", 1.0f,
                "EM", "", "NG", 0));
        return reads;
    }

    /** The expressions probed, chosen for the coercion each one forces. */
    static final String[] EXPRESSIONS = {
        // Integral versus floating versus lexicographic, over the same tag.
        "NM > 3",
        "NM > 3.5",
        "NM > '3'",
        "NM == 3",
        "NM == '3'",
        "NM != 3",
        "NM < 3",
        "NM <= 3",
        "NM >= 3",
        // Long and explicit-double literals, which take different branches again.
        "NM > 3L",
        "NM > 3.0d",
        // Arithmetic, where only `+` falls back to concatenation.
        "NM + 1 == 4",
        "NM - 1 == 2",
        "NM * 2 == 6",
        "NM / 3 == 1",
        "NM % 2 == 1",
        // A String tag with a decimal point sends even a whole-number literal down the double path.
        "RA + 1 == 31",
        "RA > 3",
        // Booleans and logic.
        "NM > 1 && AS > 1",
        "NM > 1000 || AS > 1",
        "!(NM > 1000)",
        "NM",
        "EM",
        // The empty and null cases, where the filter's own equals(Boolean.TRUE) is what throws.
        "empty(EM)",
        "empty(NM)",
        "NG == 0",
        "NG == null",
        // A tag no read carries at all.
        "ZZ == 1",
        "ZZ == null",
        // Numbers that only differ once coerced.
        "AS > NM",
        "AS == NM",
    };

    public static void main(final String[] args) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(
                List.of(new SAMSequenceRecord("chr1", CONTIG_LENGTH))));

        System.out.println("# JexlFilterDump: what the expression filter decides");

        final Map<String, Map<String, Object>> reads = reads();

        // What the context actually hands the engine, per read. Without this the port would be
        // comparing its own rendering of a tag rather than the reference's.
        for (final Map.Entry<String, Map<String, Object>> entry : reads.entrySet()) {
            final GATKRead read = makeRead(header, entry.getKey(), entry.getValue());
            for (final String tag : TAGS) {
                final String value = read.getAttributeAsString(tag);
                System.out.printf("attr\t%s\t%s=%s%n", entry.getKey(), tag,
                        value == null ? "null" : value);
            }
        }

        for (final String expression : EXPRESSIONS) {
            for (final Map.Entry<String, Map<String, Object>> entry : reads.entrySet()) {
                final GATKRead read = makeRead(header, entry.getKey(), entry.getValue());
                probeEvaluation(expression, entry.getKey(), read);
                probeFilter(expression, entry.getKey(), read);
            }
        }
    }

    /** The engine's own answer, before the filter interprets it. */
    static void probeEvaluation(final String expression, final String label, final GATKRead read) {
        String outcome;
        try {
            final Expression compiled = VariantContextUtils.engine.get().createExpression(expression);
            final Object value = compiled.evaluate(new ReadContext(read));
            outcome = value == null
                    ? "ok:null"
                    : "ok:" + value.getClass().getSimpleName() + ":" + value;
        } catch (final Exception e) {
            outcome = "E:" + e.getClass().getName();
        }
        System.out.printf("eval\t%s\t%s\t%s%n", expression, label, outcome);
    }

    /** The filter's answer, which is `!v.equals(Boolean.TRUE)` over the same value. */
    static void probeFilter(final String expression, final String label, final GATKRead read) {
        String outcome;
        try {
            final JexlExpressionReadTagValueFilter filter =
                    new JexlExpressionReadTagValueFilter(Collections.singletonList(expression));
            outcome = Boolean.toString(filter.test(read));
        } catch (final Exception e) {
            outcome = "E:" + e.getClass().getName();
        }
        System.out.printf("filter\t%s\t%s\t%s%n", expression, label, outcome);
    }

    /** The filter's own GATKReadJexlContext, which is private, rebuilt with the same three methods. */
    static final class ReadContext implements JexlContext {
        private final GATKRead read;

        ReadContext(final GATKRead read) {
            this.read = read;
        }

        @Override
        public Object get(final String name) {
            return read.getAttributeAsString(name);
        }

        @Override
        public void set(final String name, final Object value) {
            throw new IllegalArgumentException("setting attributes is not allowed");
        }

        @Override
        public boolean has(final String name) {
            return read.hasAttribute(name);
        }
    }

    static GATKRead makeRead(final SAMFileHeader header, final String label,
                             final Map<String, Object> attributes) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName("read-" + label);
        record.setReferenceName("chr1");
        record.setAlignmentStart(101);
        record.setCigar(TextCigarCodec.decode("10M"));
        record.setReadBases("ACGTACGTAC".getBytes());
        final byte[] quals = new byte[10];
        java.util.Arrays.fill(quals, (byte) 30);
        record.setBaseQualities(quals);
        record.setMappingQuality(60);
        for (final Map.Entry<String, Object> entry : attributes.entrySet()) {
            record.setAttribute(entry.getKey(), entry.getValue());
        }
        return new SAMRecordToGATKReadAdapter(record);
    }
}

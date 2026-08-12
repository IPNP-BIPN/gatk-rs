/*
 * QualQuantizer and QuantizationInfo, taken from the reference.
 *
 * The last table a recalibration report carries and the first one ApplyBQSR reads: a map from every
 * original quality score to the quantized one that replaces it. BaseRecalibrator computes it from
 * the empirical quality histogram and writes it out; ApplyBQSR reads it back and applies it, so the
 * map is the interface between the two tools.
 *
 * Eight behaviours this is built to catch.
 *
 *   - EVERY LEAF INTERVAL IS CREATED WITH A FIXED QUAL, so its error rate is qualToErrorProb of its
 *     own quality score and NOT the (nErrors+1)/(nObservations+1) the class also implements. That
 *     second formula is only ever reached by a MERGED interval, and a reader of getErrorRate would
 *     assume the opposite;
 *   - THE ERROR COUNT OF A LEAF IS floor(nObservations * errorProb) CAST THROUGH AN int and stored
 *     in a long, so a histogram bin large enough to make that product exceed Integer.MAX_VALUE
 *     SATURATES there. Measured: two bins of three billion observations at qualities 0 and 1 merge
 *     into an interval of 6000000000 observations and 4294967294 errors, which is twice
 *     Integer.MAX_VALUE and not the six billion the counts imply;
 *   - THE MERGE SEARCH KEEPS THE FIRST MINIMUM, because the comparison is strictly less-than over
 *     an iteration in qStart order, so ties go to the leftmost pair;
 *   - THE PENALTY OF A MERGED INTERVAL IS A SUM OVER ITS LEAVES, recomputed against the merged
 *     interval's own error rate, and every leaf at or below minInterestingQual contributes ZERO,
 *     which is what makes the low qualities free to merge;
 *   - A GLOBAL ERROR RATE OF ZERO IS A PENALTY OF ZERO, checked before the leaves are walked, so an
 *     empty histogram merges in iteration order rather than by cost;
 *   - errorProbToQual ROUNDS AND THEN CLAMPS TO [1, 93], so an error rate of one comes back as 1
 *     rather than 0 and an error rate of zero comes back as 93 through a saturating cast of
 *     Math.round(Infinity). A QUANTIZED QUALITY OF ZERO IS STILL REACHABLE, because a leaf that was
 *     never merged returns its FIXED qual directly and never goes through that function: a
 *     five-bin histogram at sixteen levels quantizes to 0,1,2,3,4;
 *   - nLevels = 0 IS A NullPointerException, not an empty map: the merge loop runs until one
 *     interval is left and then runs once more with nothing to pair it with;
 *   - AND QuantizationInfo COUNTS ITS OWN LEVELS BY WALKING THE MAP FOR CHANGES, so a map that
 *     returns to a value it already used counts that as another level.
 *
 * Output:
 *
 *     const\t<name>\t<value>
 *     errorprob\t<rate>\t<qual>
 *     map\t<label>\t<nLevels>\t<minInterestingQual>\t<comma separated quantized quals>
 *     interval\t<label>\t<name>\t<nObservations>\t<nErrors>\t<level>\t<fixedQual>\t<qual>\t<errorRate bits>\t<penalty bits>
 *     levels\t<label>\t<quantizationLevels>
 *     noquant\t<label>\t<comma separated quantized quals>
 *     table\t<label>\t<n>\t<the report line, with spaces shown as underscores>
 *     error\t<what>\t<exception>\t<message>
 *
 * Usage: QualQuantizerDump
 */

import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.recalibration.QualQuantizer;
import org.broadinstitute.hellbender.utils.recalibration.QuantizationInfo;
import org.broadinstitute.hellbender.utils.report.GATKReport;
import org.broadinstitute.hellbender.utils.report.GATKReportTable;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class QualQuantizerDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# QualQuantizerDump: QualQuantizer and QuantizationInfo");

        System.out.printf("const\tMIN_USABLE_Q_SCORE\t%d%n", QualityUtils.MIN_USABLE_Q_SCORE);
        System.out.printf("const\tMAX_SAM_QUAL_SCORE\t%d%n", QualityUtils.MAX_SAM_QUAL_SCORE);

        errorProbabilities();

        for (final Map.Entry<String, List<Long>> entry : histograms().entrySet()) {
            for (final int levels : new int[] {1, 2, 4, 16, 93, 94, 1000}) {
                quantize(entry.getKey() + "@" + levels, entry.getValue(), levels,
                        QualityUtils.MIN_USABLE_Q_SCORE);
            }
            // A different minimum interesting quality changes which merges are free.
            quantize(entry.getKey() + "@16,min0", entry.getValue(), 16, 0);
            quantize(entry.getKey() + "@16,min40", entry.getValue(), 16, 40);
        }

        quantizationInfo();
        errors();
    }

    /**
     * errorProbToQual, whose rounding and clamping decide every value in the quantization map.
     */
    static void errorProbabilities() {
        final double[] rates = {
                0.0, 1.0, 0.5, 0.1, 0.01, 0.001, 1.0e-9, 1.0e-10, 1.0e-30,
                // Just either side of the half that Math.round decides.
                0.0316227766016838, 0.03162277660168379,
                // The smallest positive double, whose phred score is far above the clamp.
                Double.MIN_VALUE,
        };
        for (final double rate : rates) {
            System.out.printf("errorprob\t%s\t%d%n", rate, QualityUtils.errorProbToQual(rate));
        }
        // And the arguments it refuses, which are not probabilities.
        for (final double rate : new double[] {-0.1, 1.1, Double.NaN}) {
            try {
                System.out.printf("errorprob\t%s\t%d%n", rate, QualityUtils.errorProbToQual(rate));
            } catch (final Exception e) {
                System.out.printf("error\terrorProbToQual@%s\t%s\t%s%n", rate,
                        e.getClass().getSimpleName(), e.getMessage());
            }
        }
    }

    /** The histograms, each shaped to make a different merge sequence. */
    static Map<String, List<Long>> histograms() {
        final Map<String, List<Long>> out = new LinkedHashMap<>();

        // Every bin equal: the penalty is decided entirely by the error rates, not the counts.
        out.put("flat", constant(94, 1000L));

        // Nothing at all: every error rate is that of the fixed qual, but every count is zero, so
        // the penalty is zero everywhere and the merges happen in iteration order.
        out.put("empty", constant(94, 0L));

        // What Illumina data looks like: a pile of Q2 and a peak at Q30.
        final List<Long> illumina = constant(94, 0L);
        illumina.set(2, 5_000_000L);
        for (int q = 25; q <= 35; q++) {
            illumina.set(q, 1_000_000L * (11 - Math.abs(q - 30)));
        }
        out.put("illumina", illumina);

        // One bin only, so every merge is with an empty neighbour.
        final List<Long> single = constant(94, 0L);
        single.set(40, 1234L);
        out.put("single", single);

        // Five bins: fewer than most level counts, so most runs do no merging at all.
        out.put("short", constant(5, 100L));

        // A count large enough that nObservations * errorProb overflows the int the error count is
        // cast through, at the low qualities where the error probability is near one.
        final List<Long> huge = constant(94, 0L);
        huge.set(0, 3_000_000_000L);
        huge.set(1, 3_000_000_000L);
        out.put("overflowing", huge);

        return out;
    }

    static List<Long> constant(final int size, final long value) {
        final List<Long> out = new ArrayList<>(size);
        for (int i = 0; i < size; i++) {
            out.add(value);
        }
        return out;
    }

    /** One quantization: its map, and every interval of the forest it ended with. */
    static void quantize(final String label, final List<Long> histogram, final int levels,
                         final int minInterestingQual) throws Exception {
        final QualQuantizer quantizer;
        try {
            quantizer = new QualQuantizer(histogram, levels, minInterestingQual);
        } catch (final Exception e) {
            System.out.printf("error\tquantize@%s\t%s\t%s%n", label, e.getClass().getSimpleName(),
                    e.getMessage());
            return;
        }
        System.out.printf("map\t%s\t%d\t%d\t%s%n", label, levels, minInterestingQual,
                join(quantizer.getOriginalToQuantizedMap()));

        final Field field = QualQuantizer.class.getDeclaredField("quantizedIntervals");
        field.setAccessible(true);
        final Collection<?> intervals = (Collection<?>) field.get(quantizer);
        for (final Object interval : intervals) {
            System.out.printf("interval\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s%n", label,
                    call(interval, "getName"),
                    longField(interval, "nObservations"),
                    longField(interval, "nErrors"),
                    intField(interval, "level"),
                    intField(interval, "fixedQual"),
                    (byte) call(interval, "getQual"),
                    bits((double) call(interval, "getErrorRate")),
                    bits((double) call(interval, "getPenalty")));
        }
    }

    /** QuantizationInfo's own arithmetic: the level count, noQuantization, and the report table. */
    static void quantizationInfo() throws Exception {
        for (final Map.Entry<String, List<Long>> entry : histograms().entrySet()) {
            final QualQuantizer quantizer =
                    new QualQuantizer(pad(entry.getValue()), 16, QualityUtils.MIN_USABLE_Q_SCORE);
            final QuantizationInfo info = new QuantizationInfo(
                    quantizer.getOriginalToQuantizedMap(), pad(entry.getValue()));
            System.out.printf("levels\t%s\t%d%n", entry.getKey(), info.getQuantizationLevels());

            // The report table, which is what a recalibration report's third table is.
            final GATKReport report = new GATKReport();
            final GATKReportTable table = info.generateReportTable();
            report.addTable(table);
            final String text = render(report);
            final String[] lines = text.split("\n", -1);
            for (int i = 0; i < lines.length; i++) {
                System.out.printf("table\t%s\t%d\t%s%n", entry.getKey(), i,
                        lines[i].replace(' ', '_'));
            }

            // noQuantization overwrites the map with the identity, for the first 93 entries only.
            info.noQuantization();
            System.out.printf("noquant\t%s\t%s%n", entry.getKey(), join(info.getQuantizedQuals()));
            System.out.printf("levels\t%s-after-noquant\t%d%n", entry.getKey(),
                    info.getQuantizationLevels());
        }

        // A map that returns to a value it already used, to show the level count is a count of
        // CHANGES and not of distinct values.
        final List<Byte> repeating = new ArrayList<>();
        for (final byte q : new byte[] {2, 2, 10, 10, 2, 2, 30}) {
            repeating.add(q);
        }
        final QuantizationInfo repeated = new QuantizationInfo(repeating, constant(7, 1L));
        System.out.printf("levels\trepeating\t%d%n", repeated.getQuantizationLevels());
    }

    /** A histogram at the full quality range, which is what QuantizationInfo always holds. */
    static List<Long> pad(final List<Long> histogram) {
        final List<Long> out = constant(QualityUtils.MAX_SAM_QUAL_SCORE + 1, 0L);
        for (int i = 0; i < histogram.size() && i < out.size(); i++) {
            out.set(i, histogram.get(i));
        }
        return out;
    }

    /** Every argument the quantizer refuses, and the ones it does not refuse but cannot survive. */
    static void errors() throws Exception {
        attempt("negative-counts", () -> {
            final List<Long> negative = constant(10, 1L);
            negative.set(3, -1L);
            return new QualQuantizer(negative, 4, 6);
        });
        attempt("negative-levels", () -> new QualQuantizer(constant(10, 1L), -1, 6));
        attempt("negative-min-interesting", () -> new QualQuantizer(constant(10, 1L), 4, -1));
        // Not refused, and not survivable: the merge loop runs out of pairs.
        attempt("zero-levels", () -> new QualQuantizer(constant(10, 1L), 0, 6));
        // An empty histogram: the minimum of an empty collection.
        attempt("empty-histogram", () -> new QualQuantizer(constant(0, 0L), 4, 6));
        // One bin and one level, which needs no merging at all.
        attempt("one-bin-one-level", () -> new QualQuantizer(constant(1, 5L), 1, 6));
        // The protected constructor, which leaves every field null.
        attempt("protected-constructor", () -> {
            final Constructor<QualQuantizer> constructor =
                    QualQuantizer.class.getDeclaredConstructor(int.class);
            constructor.setAccessible(true);
            final QualQuantizer quantizer = constructor.newInstance(6);
            return quantizer.getOriginalToQuantizedMap();
        });
    }

    interface Attempt {
        Object run() throws Exception;
    }

    static void attempt(final String what, final Attempt attempt) {
        try {
            final Object result = attempt.run();
            System.out.printf("error\t%s\tnone\t%s%n", what,
                    result instanceof QualQuantizer
                            ? join(((QualQuantizer) result).getOriginalToQuantizedMap())
                            : String.valueOf(result));
        } catch (final Exception e) {
            final Throwable cause =
                    e instanceof java.lang.reflect.InvocationTargetException ? e.getCause() : e;
            System.out.printf("error\t%s\t%s\t%s%n", what, cause.getClass().getSimpleName(),
                    cause.getMessage());
        }
    }

    static String join(final List<Byte> values) {
        if (values == null) {
            return "null";
        }
        final StringBuilder out = new StringBuilder();
        for (final byte value : values) {
            if (out.length() != 0) {
                out.append(',');
            }
            out.append(value);
        }
        return out.toString();
    }

    static Object call(final Object target, final String name) throws Exception {
        final Method method = target.getClass().getDeclaredMethod(name);
        method.setAccessible(true);
        return method.invoke(target);
    }

    static long longField(final Object target, final String name) throws Exception {
        final Field field = target.getClass().getDeclaredField(name);
        field.setAccessible(true);
        return field.getLong(target);
    }

    static int intField(final Object target, final String name) throws Exception {
        final Field field = target.getClass().getDeclaredField(name);
        field.setAccessible(true);
        return field.getInt(target);
    }

    static String bits(final double value) {
        return Long.toHexString(Double.doubleToRawLongBits(value));
    }

    static String render(final GATKReport report) {
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        try (final PrintStream out = new PrintStream(bytes, true, StandardCharsets.UTF_8)) {
            report.print(out);
        }
        return bytes.toString(StandardCharsets.UTF_8);
    }
}

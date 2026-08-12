/*
 * RecalDatum and EventType, taken from the reference.
 *
 * The counter every BQSR table cell is, measured before either tool that fills one. A recalibration
 * table is an array of these, and every number ApplyBQSR reads back out of a report comes from one,
 * so the arithmetic is settled here or it is settled twice.
 *
 * Nine behaviours this is built to catch, and every one of them changes a number.
 *
 *   - THE MISMATCH COUNT IS STORED MULTIPLIED BY 100000 and divided back on the way out. The field
 *     comment in the reference says the multiplier that gives sort-insensitive results is 10000.0
 *     and the constant beside it is 100000.0, so the comment is stale and the constant is the
 *     behaviour. A value that is not representable survives the round trip only up to that scaling,
 *     which is why the dump sets a mismatch count and reads it straight back;
 *   - THE EMPIRICAL QUALITY IS AN int, not a double, since February 2025. Every getter returns it
 *     widened, so `%.2f` on it always prints two zeros after the point and no other digits;
 *   - IT IS CACHED AND THE CACHE IS INVALIDATED BY EVERY SETTER, so getEmpiricalQuality(prior)
 *     computed once with one prior keeps answering that value for a different prior until a setter
 *     resets it to -1. The dump asks with two priors in both orders;
 *   - THE SMOOTHING IS ONE ERROR AND TWO OBSERVATIONS, and it is applied twice with different
 *     arithmetic: getEmpiricalErrorRate adds them to doubles, calcEmpiricalQuality adds them to
 *     `(long)(mismatches + 0.5)`, a TRUNCATING CAST after a half, not a rounding;
 *   - A NON-INTEGER QUALITY DIFFERENCE TRUNCATES TOWARD ZERO. getLogPrior computes
 *     `Math.min(Math.abs((int)(q - prior)), 40)`, so the cast happens BEFORE the absolute value and
 *     -0.5 becomes 0 rather than 1;
 *   - THE PRIOR CACHE HAS 41 ENTRIES AND THE SEARCH HAS 61 BINS. The posterior is maximised over
 *     quality scores 0..60 while the Gaussian prior is tabulated only to 40, so every difference
 *     past 40 shares one prior value;
 *   - ZERO OBSERVATIONS SHORT-CIRCUIT to a log likelihood of exactly 0.0, which is a probability of
 *     one, and an infinite or NaN likelihood becomes -Double.MAX_VALUE rather than -Infinity;
 *   - MORE THAN Integer.MAX_VALUE-1 OBSERVATIONS ARE RESCALED, errors included, by a Math.round of
 *     the scaled count;
 *   - AND combine() RECOMPUTES THE REPORTED QUALITY FROM THE EXPECTED ERRORS OF BOTH SIDES rather
 *     than averaging, which for two empty datums is -10*log10(0/0) and therefore NaN, with no
 *     validation to stop it because the field is assigned directly and not through the setter.
 *
 * Every double travels as its raw bits as well as its decimal, because the port has to reproduce
 * the bits and a decimal rendering cannot show the last of them.
 *
 * Output:
 *
 *     const\t<name>\t<value>
 *     logprior\t<i>\t<bits>\t<decimal>
 *     getlogprior\t<quality>\t<prior>\t<bits>\t<decimal>
 *     loglik\t<quality>\t<nObservations>\t<nErrors>\t<bits>\t<decimal>
 *     bayes\t<nObservations>\t<nErrors>\t<prior>\t<result>
 *     datum\t<label>\t<field>\t<bits>\t<decimal>
 *     text\t<label>\t<toString>\t<stringForCSV>
 *     cache\t<label>\t<step>\t<empiricalQuality>
 *     combine\t<label>\t<field>\t<bits>\t<decimal>
 *     narrow\t<what>\t<byte>
 *     error\t<what>\t<exception>\t<message>
 *     event\t<ordinal>\t<toString>\t<prettyPrint>
 *
 * Usage: RecalDatumDump
 */

import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.recalibration.EventType;
import org.broadinstitute.hellbender.utils.recalibration.RecalDatum;

import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

public class RecalDatumDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# RecalDatumDump: RecalDatum and EventType");

        constants();
        priorCache();
        logPriors();
        logLikelihoods();
        bayes();
        datums();
        caching();
        combining();
        narrowing();
        errors();
        events();
    }

    /** The numbers the class is built out of, so the port cannot quietly pick different ones. */
    static void constants() throws Exception {
        System.out.printf("const\tMAX_RECALIBRATED_Q_SCORE\t%d%n", RecalDatum.MAX_RECALIBRATED_Q_SCORE);
        System.out.printf("const\tMAX_GATK_USABLE_Q_SCORE\t%d%n", RecalDatum.MAX_GATK_USABLE_Q_SCORE);
        System.out.printf("const\tMAX_REASONABLE_Q_SCORE\t%d%n", QualityUtils.MAX_REASONABLE_Q_SCORE);
        // Private, and the multiplier is the one number that decides whether two ports agree.
        System.out.printf("const\tMULTIPLIER\t%s%n", privateStatic("MULTIPLIER"));
        System.out.printf("const\tSMOOTHING_CONSTANT\t%s%n", privateStatic("SMOOTHING_CONSTANT"));
        System.out.printf("const\tUNINITIALIZED_EMPIRICAL_QUALITY\t%s%n",
                privateStatic("UNINITIALIZED_EMPIRICAL_QUALITY"));
        System.out.printf("const\tMAX_NUMBER_OF_OBSERVATIONS\t%s%n",
                privateStatic("MAX_NUMBER_OF_OBSERVATIONS"));
    }

    static String privateStatic(final String name) throws Exception {
        final Field field = RecalDatum.class.getDeclaredField(name);
        field.setAccessible(true);
        return String.valueOf(field.get(null));
    }

    /**
     * The 41 entries of logPriorCache, which is a Gaussian log density with mean 0 and sigma 0.5
     * over commons-math's own FastMath rather than java.lang.Math.
     */
    static void priorCache() throws Exception {
        final Field field = RecalDatum.class.getDeclaredField("logPriorCache");
        field.setAccessible(true);
        final double[] cache = (double[]) field.get(null);
        System.out.printf("const\tlogPriorCache.length\t%d%n", cache.length);
        for (int i = 0; i < cache.length; i++) {
            emit("logprior", String.valueOf(i), cache[i]);
        }
    }

    /**
     * getLogPrior, whose cast to int runs BEFORE the absolute value and whose clamp is at 40 while
     * the search that calls it runs to 60.
     */
    static void logPriors() throws Exception {
        final Method method = RecalDatum.class.getDeclaredMethod("getLogPrior", double.class, double.class);
        method.setAccessible(true);
        final double[][] pairs = {
                {0.0, 0.0},         // no difference at all
                {30.0, 30.0},       // the same, away from zero
                {31.0, 30.0},       // one above
                {29.0, 30.0},       // one below, the same cache entry
                {30.5, 30.0},       // half a point: (int) 0.5 is 0, so this is the same as no difference
                {29.5, 30.0},       // (int) -0.5 is 0 too, which is the cast running before the abs
                {30.9, 30.0},       // and so is nine tenths
                {0.0, 40.0},        // exactly at the clamp
                {0.0, 41.0},        // past it, and therefore the same value as at it
                {0.0, 60.0},        // the top of the search range, still the clamped prior
                {60.0, 0.0},        // the same difference the other way round
                {10.0, 30.7},       // a non-integer prior: (int) -20.7 is -20, not -21
        };
        for (final double[] pair : pairs) {
            final double value = (double) method.invoke(null, pair[0], pair[1]);
            System.out.printf("getlogprior\t%s\t%s\t%s\t%s%n", show(pair[0]), show(pair[1]),
                    Long.toHexString(Double.doubleToRawLongBits(value)), value);
        }
    }

    /**
     * getLogBinomialLikelihood, with its three escapes: no observations, the rescaling above
     * Integer.MAX_VALUE-1, and the substitution of -Double.MAX_VALUE for a non-finite result.
     */
    static void logLikelihoods() throws Exception {
        final Method method = RecalDatum.class.getDeclaredMethod(
                "getLogBinomialLikelihood", double.class, long.class, long.class);
        method.setAccessible(true);
        final Object[][] cases = {
                {0.0, 0L, 0L},               // no observations: exactly 0.0, whatever the quality
                {30.0, 0L, 0L},
                {30.0, 100L, 1L},            // the ordinary middle branch
                {30.0, 100L, 0L},            // x == 0, and p = 0.001 is under the 0.1 test
                {2.0, 100L, 0L},             // x == 0 with p = 0.63, which is over it
                {30.0, 100L, 100L},          // x == n with q near one, over the 0.1 test
                {0.0, 100L, 100L},           // x == n with q = 0, which is under it
                {0.0, 100L, 1L},             // quality zero: p = 1, so one error has probability 0
                {60.0, 1000000L, 1L},
                {30.0, 100L, 200L},          // more errors than observations
                {30.0, 3000000000L, 3000L},  // over Integer.MAX_VALUE-1: both counts are rescaled
                {30.0, 2147483646L, 5L},     // exactly at the limit, so no rescaling
                {30.0, 2147483647L, 5L},     // one past it, so rescaling by a fraction just under 1
        };
        for (final Object[] c : cases) {
            final double value = (double) method.invoke(null, c[0], c[1], c[2]);
            System.out.printf("loglik\t%s\t%d\t%d\t%s\t%s%n", show((Double) c[0]), (Long) c[1], (Long) c[2],
                    Long.toHexString(Double.doubleToRawLongBits(value)), value);
        }
    }

    /** The argmax itself, over the 61 bins, including the ties the first-wins rule decides. */
    static void bayes() {
        final long[][] counts = {
                {0, 0}, {1, 0}, {1, 1}, {10, 1}, {100, 1}, {100, 0}, {1000, 10}, {1000, 0},
                {100, 100}, {10000, 3}, {2, 1}, {3000000000L, 3000},
        };
        final double[] priors = {0.0, 10.0, 20.0, 30.0, 45.0, 60.0};
        for (final long[] count : counts) {
            for (final double prior : priors) {
                System.out.printf("bayes\t%d\t%d\t%s\t%d%n", count[0], count[1], show(prior),
                        RecalDatum.bayesianEstimateOfEmpiricalQuality(count[0], count[1], prior));
            }
        }
    }

    /** One datum per shape, with every reader on it. */
    static void datums() {
        emitDatum("plain", new RecalDatum(1000L, 10.0, (byte) 30));
        emitDatum("empty", new RecalDatum(0L, 0.0, (byte) 30));
        emitDatum("perfect", new RecalDatum(1000L, 0.0, (byte) 30));
        emitDatum("all-errors", new RecalDatum(1000L, 1000.0, (byte) 30));
        emitDatum("qual-zero", new RecalDatum(1000L, 10.0, (byte) 0));
        emitDatum("qual-max", new RecalDatum(1000L, 10.0, (byte) 93));
        // A mismatch count that is not a tenth: the multiplier decides what comes back out.
        emitDatum("fractional", new RecalDatum(1000L, 0.1, (byte) 30));
        emitDatum("tiny-fraction", new RecalDatum(1000L, 1.0e-7, (byte) 30));
        // Half a mismatch, which is the value the truncating cast in calcEmpiricalQuality turns on.
        emitDatum("half", new RecalDatum(1000L, 0.5, (byte) 30));
        emitDatum("just-under-half", new RecalDatum(1000L, 0.49999, (byte) 30));
        emitDatum("huge", new RecalDatum(3000000000L, 3000.0, (byte) 30));
        // The empirical quality this one wants is above 93 and is capped there.
        emitDatum("capped", new RecalDatum(100000000L, 0.0, (byte) 93));

        // The copy constructor copies the RAW mismatch field, cache included, without rescaling.
        final RecalDatum source = new RecalDatum(1000L, 10.0, (byte) 30);
        source.getEmpiricalQuality();
        emitDatum("copy-of-computed", new RecalDatum(source));

        // A setter takes the value unscaled and stores it scaled, like the constructor.
        final RecalDatum set = new RecalDatum(1000L, 10.0, (byte) 30);
        set.setNumMismatches(0.1);
        set.setNumObservations(7L);
        set.setReportedQuality(20.5);
        emitDatum("after-setters", set);

        // setEmpiricalQuality writes the cache directly, so the counts no longer explain it.
        final RecalDatum forced = new RecalDatum(1000L, 10.0, (byte) 30);
        forced.setEmpiricalQuality(7);
        emitDatum("forced-empirical", forced);

        // increment(boolean) is one observation and one or no error.
        final RecalDatum incremented = new RecalDatum(0L, 0.0, (byte) 30);
        incremented.increment(true);
        incremented.increment(false);
        incremented.increment(false);
        incremented.incrementNumObservations(10L);
        incremented.incrementNumMismatches(0.25);
        emitDatum("incremented", incremented);
    }

    static void emitDatum(final String label, final RecalDatum datum) {
        emit("datum", label + "\tnumObservations", (double) datum.getNumObservations());
        emit("datum", label + "\tnumMismatches", datum.getNumMismatches());
        emit("datum", label + "\treportedQuality", datum.getReportedQuality());
        emit("datum", label + "\tempiricalErrorRate", datum.getEmpiricalErrorRate());
        System.out.printf("datum\t%s\treportedQualityAsByte\t-\t%d%n", label, datum.getReportedQualityAsByte());
        // toString before any getEmpiricalQuality call, because it calls one itself.
        System.out.printf("text\t%s\t%s\t%s%n", label, datum.toString(), datum.stringForCSV());
        emit("datum", label + "\tempiricalQuality", datum.getEmpiricalQuality());
        System.out.printf("datum\t%s\tempiricalQualityAsByte\t-\t%d%n", label, datum.getEmpiricalQualityAsByte());
    }

    /**
     * The cache and its invalidation: one prior computes a value, a second prior gets the first
     * one back, and a setter makes the second prior take effect.
     */
    static void caching() {
        final RecalDatum datum = new RecalDatum(1000L, 10.0, (byte) 30);
        System.out.printf("cache\tprior-then-prior\tfirst-with-10\t%s%n", datum.getEmpiricalQuality(10.0));
        System.out.printf("cache\tprior-then-prior\tthen-with-45\t%s%n", datum.getEmpiricalQuality(45.0));
        datum.setNumObservations(1000L);  // the same value, but a setter all the same
        System.out.printf("cache\tprior-then-prior\tafter-setter-with-45\t%s%n", datum.getEmpiricalQuality(45.0));

        // The other order, to show neither prior is privileged.
        final RecalDatum other = new RecalDatum(1000L, 10.0, (byte) 30);
        System.out.printf("cache\treversed\tfirst-with-45\t%s%n", other.getEmpiricalQuality(45.0));
        System.out.printf("cache\treversed\tthen-with-10\t%s%n", other.getEmpiricalQuality(10.0));

        // The no-argument getter uses the reported quality as the prior, and it caches the same way.
        final RecalDatum implicitPrior = new RecalDatum(1000L, 10.0, (byte) 30);
        System.out.printf("cache\timplicit\tno-argument\t%s%n", implicitPrior.getEmpiricalQuality());
        System.out.printf("cache\timplicit\tthen-with-0\t%s%n", implicitPrior.getEmpiricalQuality(0.0));

        // setEmpiricalQuality wins over any prior until a count changes.
        final RecalDatum forced = new RecalDatum(1000L, 10.0, (byte) 30);
        forced.setEmpiricalQuality(3);
        System.out.printf("cache\tforced\twith-45\t%s%n", forced.getEmpiricalQuality(45.0));
        forced.incrementNumObservations(0L);  // adds nothing, invalidates anyway
        System.out.printf("cache\tforced\tafter-zero-increment\t%s%n", forced.getEmpiricalQuality(45.0));
    }

    /** combine(), which recomputes the reported quality from both sides' expected errors. */
    static void combining() {
        final RecalDatum a = new RecalDatum(1000L, 10.0, (byte) 30);
        a.combine(new RecalDatum(1000L, 10.0, (byte) 20));
        emitCombine("different-qualities", a);

        final RecalDatum same = new RecalDatum(1000L, 10.0, (byte) 30);
        same.combine(new RecalDatum(1000L, 10.0, (byte) 30));
        emitCombine("same-quality", same);

        // Two empty datums: the expected errors are zero and so are the observations, so the
        // reported quality is -10*log10(0/0), and nothing validates it.
        final RecalDatum empty = new RecalDatum(0L, 0.0, (byte) 30);
        empty.combine(new RecalDatum(0L, 0.0, (byte) 30));
        emitCombine("both-empty", empty);

        // One empty side, which contributes no expected errors and no observations.
        final RecalDatum oneEmpty = new RecalDatum(1000L, 10.0, (byte) 30);
        oneEmpty.combine(new RecalDatum(0L, 0.0, (byte) 30));
        emitCombine("one-empty", oneEmpty);

        // Quality zero on one side: every base is expected to be wrong, so the combined reported
        // quality is dragged most of the way down.
        final RecalDatum withZero = new RecalDatum(1000L, 10.0, (byte) 40);
        withZero.combine(new RecalDatum(1000L, 10.0, (byte) 0));
        emitCombine("with-quality-zero", withZero);

        // A combine after the empirical quality was computed, to show the cache is dropped.
        final RecalDatum computed = new RecalDatum(1000L, 10.0, (byte) 30);
        computed.getEmpiricalQuality();
        computed.combine(new RecalDatum(1000L, 500.0, (byte) 30));
        emitCombine("after-computing", computed);

        // A second combine on the NaN reported quality the first one produced. calcExpectedErrors
        // asks QualityUtils.qualToErrorProb(double), whose guard is `qual >= 0.0`, and that is
        // false for NaN, so the reachable end of the both-empty case is an exception rather than
        // another NaN. This is why the port cannot simply drop the guard.
        error("combine-onto-nan-reported", () -> {
            final RecalDatum nan = new RecalDatum(0L, 0.0, (byte) 30);
            nan.combine(new RecalDatum(0L, 0.0, (byte) 30));
            nan.combine(new RecalDatum(1000L, 10.0, (byte) 30));
            return null;
        });
    }

    /** The two narrowing casts, which are the only place a quality can come back negative. */
    static void narrowing() {
        final RecalDatum big = new RecalDatum(1L, 0.0, (byte) 30);
        big.setReportedQuality(200.0);
        System.out.printf("narrow\treported-200\t%d%n", big.getReportedQualityAsByte());
        big.setReportedQuality(127.4);
        System.out.printf("narrow\treported-127.4\t%d%n", big.getReportedQualityAsByte());
        big.setReportedQuality(127.5);
        System.out.printf("narrow\treported-127.5\t%d%n", big.getReportedQualityAsByte());
        big.setReportedQuality(0.5);
        System.out.printf("narrow\treported-0.5\t%d%n", big.getReportedQualityAsByte());
        big.setReportedQuality(1.5);
        System.out.printf("narrow\treported-1.5\t%d%n", big.getReportedQualityAsByte());

        final RecalDatum forced = new RecalDatum(1L, 0.0, (byte) 30);
        forced.setEmpiricalQuality(200);
        System.out.printf("narrow\tempirical-200\t%d%n", forced.getEmpiricalQualityAsByte());
        forced.setEmpiricalQuality(93);
        System.out.printf("narrow\tempirical-93\t%d%n", forced.getEmpiricalQualityAsByte());
    }

    static void emitCombine(final String label, final RecalDatum datum) {
        emit("combine", label + "\tnumObservations", (double) datum.getNumObservations());
        emit("combine", label + "\tnumMismatches", datum.getNumMismatches());
        emit("combine", label + "\treportedQuality", datum.getReportedQuality());
        emit("combine", label + "\tempiricalQuality", datum.getEmpiricalQuality());
    }

    /** Every argument the class refuses, and the words it refuses it with. */
    static void errors() throws Exception {
        error("constructor-negative-observations", () -> new RecalDatum(-1L, 0.0, (byte) 30));
        error("constructor-negative-mismatches", () -> new RecalDatum(1L, -1.0, (byte) 30));
        error("constructor-negative-quality", () -> new RecalDatum(1L, 0.0, (byte) -1));
        error("set-negative-observations", () -> {
            new RecalDatum(1L, 0.0, (byte) 30).setNumObservations(-1L);
            return null;
        });
        error("set-negative-mismatches", () -> {
            new RecalDatum(1L, 0.0, (byte) 30).setNumMismatches(-1.0);
            return null;
        });
        error("set-negative-reported", () -> {
            new RecalDatum(1L, 0.0, (byte) 30).setReportedQuality(-1.0);
            return null;
        });
        error("set-infinite-reported", () -> {
            new RecalDatum(1L, 0.0, (byte) 30).setReportedQuality(Double.POSITIVE_INFINITY);
            return null;
        });
        error("set-nan-reported", () -> {
            new RecalDatum(1L, 0.0, (byte) 30).setReportedQuality(Double.NaN);
            return null;
        });
        error("set-negative-empirical", () -> {
            new RecalDatum(1L, 0.0, (byte) 30).setEmpiricalQuality(-1);
            return null;
        });
        // A negative increment is not checked at all, so the counts can go negative.
        final RecalDatum negative = new RecalDatum(1L, 0.0, (byte) 30);
        negative.increment(-5L, -5.0);
        System.out.printf("error\tincrement-negative-is-allowed\t-\t%d,%s%n",
                negative.getNumObservations(), negative.getNumMismatches());
        // And so is a reported quality of 127 as a byte, which is above MAX_RECALIBRATED_Q_SCORE.
        System.out.printf("error\tquality-above-max-is-allowed\t-\t%s%n",
                new RecalDatum(1L, 0.0, (byte) 127).getReportedQuality());
        error("event-from-unknown", () -> EventType.eventFrom("X"));
        error("event-from-index-past-end", () -> EventType.eventFrom(3));
        error("event-from-negative-index", () -> EventType.eventFrom(-1));
    }

    interface Attempt {
        Object run() throws Exception;
    }

    static void error(final String what, final Attempt attempt) {
        try {
            attempt.run();
            System.out.printf("error\t%s\tnone\t-%n", what);
        } catch (final InvocationTargetException e) {
            final Throwable cause = e.getCause();
            System.out.printf("error\t%s\t%s\t%s%n", what, cause.getClass().getSimpleName(),
                    cause.getMessage());
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s\t%s%n", what, e.getClass().getSimpleName(), e.getMessage());
        }
    }

    /** The three event types, which are the three tables of a recalibration report. */
    static void events() {
        for (int i = 0; i < 3; i++) {
            final EventType event = EventType.eventFrom(i);
            System.out.printf("event\t%d\t%s\t%s\t%s\t%s%n", i, event.name(), event.toString(),
                    event.prettyPrint(), EventType.eventFrom(event.toString()).name());
        }
    }

    static void emit(final String kind, final String label, final double value) {
        System.out.printf("%s\t%s\t%s\t%s%n", kind, label,
                Long.toHexString(Double.doubleToRawLongBits(value)), value);
    }

    /** A double as Java writes it, so a decimal in a label matches a decimal in a row. */
    static String show(final double value) {
        return String.valueOf(value);
    }
}

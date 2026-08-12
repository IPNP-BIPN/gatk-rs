/*
 * BQSRReadTransformer, taken from the reference.
 *
 * The thing ApplyBQSR is: a read in, the same read out with every base quality replaced. Everything
 * measured before this (RecalDatum, the covariates, the tables, the quantizer, the report reader) is
 * an input to this one function.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE RECALIBRATION OF A READ DEPENDS ON WHICH READS CAME BEFORE IT. The hierarchical estimate
 *     calls RecalDatum.getEmpiricalQuality(prior), which computes on the FIRST call and returns the
 *     cached value on every call after that, whatever prior is asked for. The datums live in the
 *     table and are shared across every read, so the first read to reach a datum fixes its empirical
 *     quality for the whole run. The reference's own TODO says "the prior is ignored if the
 *     empirical quality for the datum is already cached". The dump runs the same two reads in both
 *     orders over freshly built tables and prints both results;
 *   - THE ESTIMATE IS y_3 + y_4 - y_2, not a sum of corrections: each special covariate contributes
 *     its own empirical quality MINUS the two-covariate posterior, and a null datum contributes
 *     nothing at all, which is not the same as contributing zero delta from a prior;
 *   - A COVARIATE KEY OF -1 IS SKIPPED, so the first bases of a read, which have no context, are
 *     recalibrated from fewer covariates than the rest;
 *   - A QUALITY BELOW preserveQLessThan IS LEFT ALONE ENTIRELY, not even quantized;
 *   - THE ROUNDING IS fastRound, `(int)(x + 0.5)`, and not Math.round, so it rounds the double below
 *     a half twice;
 *   - --allow-missing-read-group COVERS A DIFFERENT CASE THAN ITS NAME SUGGESTS. It acts only when
 *     the COVARIATE does not know the read group, which happens when the recalibration report was
 *     written from a different set of read groups than the BAM has; the read is then QUANTIZED BUT
 *     NOT RECALIBRATED. A read whose group the covariate knows but whose TABLE holds no datum is a
 *     GATKException whatever the flag says, because the flag is tested before that lookup;
 *   - constructStaticQuantizedMapping SORTS ITS ARGUMENT LIST IN PLACE, preserves every quality
 *     below MIN_USABLE_Q_SCORE as itself, and rounds in PROBABILITY space rather than in Phred
 *     space, so the midpoint between two static quals is not their arithmetic mean;
 *   - AND WITH STATIC QUANTIZATION THE QUALITY IS QUANTIZED TWICE, first through the dynamic map and
 *     then through the static one, which the reference's own TODO calls out.
 *
 * Output:
 *
 *     static\t<label>\t<comma separated mapping>
 *     sorted\t<label>\t<the argument list after the call>
 *     estimate\t<label>\t<bits>\t<decimal>
 *     order\t<label>\t<read>\t<comma separated recalibrated qualities>
 *     apply\t<label>\t<read>\t<comma separated recalibrated qualities>
 *     oqtag\t<label>\t<read>\t<OQ>
 *     round\t<value>\t<bounded qual>
 *     error\t<what>\t<exception>\t<message>
 *
 * Usage: BqsrTransformerDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.tools.ApplyBQSRArgumentCollection;
import org.broadinstitute.hellbender.transformers.BQSRReadTransformer;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.recalibration.QuantizationInfo;
import org.broadinstitute.hellbender.utils.recalibration.RecalDatum;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;

import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class BqsrTransformerDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# BqsrTransformerDump: BQSRReadTransformer");

        staticMappings();
        rounding();
        estimates();
        applying();
    }

    /** constructStaticQuantizedMapping, which rounds in probability space and sorts in place. */
    static void staticMappings() {
        final List<List<Integer>> cases = List.of(
                List.of(),
                List.of(30),
                List.of(10, 20, 30, 40),
                // Deliberately out of order, to show the call sorts the caller's own list.
                List.of(40, 10, 30, 20),
                List.of(6, 93),
                List.of(2, 3)
        );
        for (final List<Integer> quals : cases) {
            for (final boolean roundDown : new boolean[] {false, true}) {
                final List<Integer> mutable = new ArrayList<>(quals);
                final byte[] mapping =
                        BQSRReadTransformer.constructStaticQuantizedMapping(mutable, roundDown);
                final String label = quals + "@" + (roundDown ? "down" : "nearest");
                System.out.printf("static\t%s\t%s%n", label, join(mapping));
                System.out.printf("sorted\t%s\t%s%n", label, mutable);
            }
        }
    }

    /** getBoundedIntegerQual: fastRound and then the [1, 93] clamp. */
    static void rounding() {
        final double[] values = {
                0.0, 0.4, 0.5, 1.4, 1.5, 2.5, -1.0, -0.5, 92.5, 93.4, 93.5, 200.0,
                // The double just below a half, which fastRound rounds twice and Math.round does not.
                0.49999999999999994,
        };
        for (final double value : values) {
            System.out.printf("round\t%s\t%d\t%d%n", value, MathUtils.fastRound(value),
                    QualityUtils.boundQual(MathUtils.fastRound(value), RecalDatum.MAX_RECALIBRATED_Q_SCORE));
        }
    }

    /** hierarchicalBayesianQualityEstimate over every combination of present and absent datums. */
    static void estimates() {
        final double prior = 25.0;
        System.out.printf("estimate\tall-null\t%s%n",
                bits(BQSRReadTransformer.hierarchicalBayesianQualityEstimate(prior, null, null, null, null)));
        System.out.printf("estimate\tread-group-only\t%s%n",
                bits(BQSRReadTransformer.hierarchicalBayesianQualityEstimate(prior,
                        datum(10000, 100.0, 30), null, null, null)));
        System.out.printf("estimate\ttwo-covariates\t%s%n",
                bits(BQSRReadTransformer.hierarchicalBayesianQualityEstimate(prior,
                        datum(10000, 100.0, 30), datum(5000, 20.0, 30), null, null)));
        System.out.printf("estimate\tone-special\t%s%n",
                bits(BQSRReadTransformer.hierarchicalBayesianQualityEstimate(prior,
                        datum(10000, 100.0, 30), datum(5000, 20.0, 30), datum(1000, 1.0, 30), null)));
        System.out.printf("estimate\ttwo-specials\t%s%n",
                bits(BQSRReadTransformer.hierarchicalBayesianQualityEstimate(prior,
                        datum(10000, 100.0, 30), datum(5000, 20.0, 30), datum(1000, 1.0, 30),
                        datum(2000, 200.0, 30))));
        // A special covariate with no read group or quality datum: the prior flows straight through
        // as the posterior, and the delta is taken against it.
        System.out.printf("estimate\tspecial-only\t%s%n",
                bits(BQSRReadTransformer.hierarchicalBayesianQualityEstimate(prior, null, null,
                        datum(1000, 1.0, 30), null)));

        // The caching: the SAME datum asked twice with different priors answers the first one both
        // times, so the order the estimates are computed in changes the answer.
        final RecalDatum shared = datum(1000, 1.0, 30);
        System.out.printf("estimate\tcached-first-25\t%s%n",
                bits(BQSRReadTransformer.hierarchicalBayesianQualityEstimate(25.0, null, shared)));
        System.out.printf("estimate\tcached-then-45\t%s%n",
                bits(BQSRReadTransformer.hierarchicalBayesianQualityEstimate(45.0, null, shared)));
        final RecalDatum fresh = datum(1000, 1.0, 30);
        System.out.printf("estimate\tfresh-45\t%s%n",
                bits(BQSRReadTransformer.hierarchicalBayesianQualityEstimate(45.0, null, fresh)));
    }

    /** The transformer itself, over a small corpus and a hand-built recalibration table. */
    static void applying() throws Exception {
        final SAMFileHeader header = ReadFilterDump.header();

        // Two reads of the same shape in one read group, differing only in their qualities. Running
        // them in both orders shows whether the datum cache carries between reads.
        final SAMRecord first = read(header, "first", "ACGTACGTAC", new byte[] {30, 30, 30, 30, 30, 30, 30, 30, 30, 30});
        final SAMRecord second = read(header, "second", "ACGTACGTAC", new byte[] {20, 20, 20, 20, 20, 20, 20, 20, 20, 20});
        // A read with qualities below the preserve threshold, which are left untouched.
        final SAMRecord low = read(header, "low", "ACGTACGTAC", new byte[] {2, 2, 5, 5, 6, 6, 30, 30, 2, 40});
        // A read whose group is not in the recalibration table.
        final SAMRecord unknown = read(header, "unknown", "ACGTACGTAC", new byte[] {30, 30, 30, 30, 30, 30, 30, 30, 30, 30});
        unknown.setAttribute("RG", "rg2");

        // The header and the four reads travel, because the read group identifiers are read out of
        // the header and the transformer's answer depends on them.
        ReadFilterDump.printCorpus(header, List.of(first, second, low, unknown));

        apply("in-order", header, List.of(first, second, low), defaults(), false);
        apply("reversed", header, List.of(second, first, low), defaults(), false);

        final ApplyBQSRArgumentCollection quantized = defaults();
        quantized.quantizationLevels = 4;
        apply("quantized-4", header, List.of(first, second), quantized, false);

        final ApplyBQSRArgumentCollection noQuantization = defaults();
        noQuantization.quantizationLevels = 0;
        apply("no-quantization", header, List.of(first, second), noQuantization, false);

        final ApplyBQSRArgumentCollection statics = defaults();
        statics.staticQuantizationQuals = new ArrayList<>(List.of(10, 20, 30, 40));
        apply("static-quals", header, List.of(first, second), statics, false);

        final ApplyBQSRArgumentCollection prior = defaults();
        prior.globalQScorePrior = 20.0;
        apply("global-prior", header, List.of(first, second), prior, false);

        final ApplyBQSRArgumentCollection preserve = defaults();
        preserve.PRESERVE_QSCORES_LESS_THAN = 31;
        apply("preserve-31", header, List.of(first, second, low), preserve, false);

        final ApplyBQSRArgumentCollection original = defaults();
        original.emitOriginalQuals = true;
        apply("emit-original", header, List.of(first), original, true);

        // A read group the TABLE has no datum for, which is not what --allow-missing-read-group
        // covers: the covariate still knows the group, so the key is not the missing code and the
        // flag never gets a chance to act.
        apply("no-datum-for-read-group", header, List.of(unknown), defaults(), false);
        final ApplyBQSRArgumentCollection allow = defaults();
        allow.allowMissingReadGroups = true;
        apply("no-datum-for-read-group-allowed", header, List.of(unknown), allow, false);

        // A read group the COVARIATE does not know, which is the case the flag is for. The
        // covariates come from the recalibration report and the header comes from the BAM, so a
        // report written from one read group and applied to a two-group BAM lands here.
        applyWithOneReadGroup("covariate-missing-read-group", header, List.of(first, unknown),
                defaults());
        final ApplyBQSRArgumentCollection allowCovariate = defaults();
        allowCovariate.allowMissingReadGroups = true;
        applyWithOneReadGroup("covariate-missing-read-group-allowed", header,
                List.of(first, unknown), allowCovariate);
    }

    /**
     * The same run, but with covariates built from ONE read group while the header holds three.
     * This is what applying a recalibration report to a BAM it was not written from looks like.
     */
    static void applyWithOneReadGroup(final String label, final SAMFileHeader header,
                                      final List<SAMRecord> reads,
                                      final ApplyBQSRArgumentCollection args) throws Exception {
        final RecalibrationArgumentCollection rac = new RecalibrationArgumentCollection();
        final StandardCovariateList covariates =
                new StandardCovariateList(rac, List.of("unit-rg1"));
        final RecalibrationTables tables = buildTables(covariates);
        final QuantizationInfo quantization = new QuantizationInfo(tables, rac.QUANTIZING_LEVELS);

        final Constructor<BQSRReadTransformer> constructor =
                BQSRReadTransformer.class.getDeclaredConstructor(SAMFileHeader.class,
                        RecalibrationTables.class, QuantizationInfo.class,
                        StandardCovariateList.class, ApplyBQSRArgumentCollection.class);
        constructor.setAccessible(true);
        final BQSRReadTransformer transformer =
                constructor.newInstance(header, tables, quantization, covariates, args);

        for (final SAMRecord record : reads) {
            final GATKRead read = new SAMRecordToGATKReadAdapter(record.deepCopy());
            try {
                final GATKRead out = transformer.apply(read);
                System.out.printf("apply\t%s\t%s\t%s%n", label, read.getName(),
                        join(out.getBaseQualities()));
            } catch (final Exception e) {
                System.out.printf("apply\t%s\t%s\tE:%s:%s%n", label, read.getName(),
                        e.getClass().getSimpleName(), e.getMessage());
            }
        }
    }

    static ApplyBQSRArgumentCollection defaults() {
        return new ApplyBQSRArgumentCollection();
    }

    /** One run of the transformer over a list of reads, with a table built fresh for each run. */
    static void apply(final String label, final SAMFileHeader header, final List<SAMRecord> reads,
                      final ApplyBQSRArgumentCollection args, final boolean printTags)
            throws Exception {
        final RecalibrationArgumentCollection rac = new RecalibrationArgumentCollection();
        final StandardCovariateList covariates = new StandardCovariateList(rac, header);
        final RecalibrationTables tables = buildTables(covariates);
        final QuantizationInfo quantization = new QuantizationInfo(tables, rac.QUANTIZING_LEVELS);

        final Constructor<BQSRReadTransformer> constructor =
                BQSRReadTransformer.class.getDeclaredConstructor(SAMFileHeader.class,
                        RecalibrationTables.class, QuantizationInfo.class,
                        StandardCovariateList.class, ApplyBQSRArgumentCollection.class);
        constructor.setAccessible(true);
        final BQSRReadTransformer transformer =
                constructor.newInstance(header, tables, quantization, covariates, args);

        for (final SAMRecord record : reads) {
            final GATKRead read = new SAMRecordToGATKReadAdapter(record.deepCopy());
            try {
                final GATKRead out = transformer.apply(read);
                System.out.printf("apply\t%s\t%s\t%s%n", label, read.getName(),
                        join(out.getBaseQualities()));
                if (printTags) {
                    // Not "tag": the shared corpus printer already owns that row kind, and a reader of
                    // the golden would take this for one of its rows.
                    System.out.printf("oqtag\t%s\t%s\t%s%n", label, read.getName(),
                            String.valueOf(out.getAttributeAsString("OQ")));
                }
            } catch (final Exception e) {
                System.out.printf("apply\t%s\t%s\tE:%s:%s%n", label, read.getName(),
                        e.getClass().getSimpleName(), e.getMessage());
            }
        }
    }

    /**
     * A recalibration table with one datum in every place the transformer looks, so no lookup is
     * null and the whole estimate is exercised.
     */
    static RecalibrationTables buildTables(final StandardCovariateList covariates) {
        final RecalibrationTables tables = new RecalibrationTables(covariates);
        // Read group 0, base substitution.
        tables.getReadGroupTable().put(datum(100000, 1000.0, 30), 0, 0);
        // Every reported quality this corpus carries.
        for (final int quality : new int[] {2, 5, 6, 20, 30, 40}) {
            tables.getQualityScoreTable().put(datum(10000, 50.0, quality), 0, quality, 0);
            // And every context and cycle key those reads produce, which is more than they use.
            for (int key = 0; key < 260; key++) {
                tables.getTable(2).put(datum(1000, 5.0 + (key % 7), quality), 0, quality, key, 0);
            }
            for (int key = 0; key < 24; key++) {
                tables.getTable(3).put(datum(1000, 3.0 + (key % 5), quality), 0, quality, key, 0);
            }
        }
        return tables;
    }

    static RecalDatum datum(final long observations, final double mismatches, final int quality) {
        return new RecalDatum(observations, mismatches, (byte) quality);
    }

    static SAMRecord read(final SAMFileHeader header, final String name, final String bases,
                          final byte[] qualities) {
        final SAMRecord record =
                ReadFilterDump.read(header, name, 0, 0, 100, 60, "10M", 0, 200, 100, true);
        record.setReadString(bases);
        record.setBaseQualities(qualities);
        return record;
    }

    static String bits(final double value) {
        return Long.toHexString(Double.doubleToRawLongBits(value)) + "\t" + value;
    }

    static String join(final byte[] values) {
        final StringBuilder out = new StringBuilder();
        for (final byte value : values) {
            if (out.length() != 0) {
                out.append(',');
            }
            out.append(value);
        }
        return out.toString();
    }
}

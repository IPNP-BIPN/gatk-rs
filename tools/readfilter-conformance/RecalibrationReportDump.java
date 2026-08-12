/*
 * RecalibrationReport, taken from the reference.
 *
 * A recalibration table read off disk. It is what ApplyBQSR is given: five named tables in a
 * GATKReport, turned back into a RecalibrationTables, a StandardCovariateList and a QuantizationInfo.
 * The reader, the covariates, the tables and the quantizer are all measured already; this is the
 * assembly, and the assembly is where the keys are chosen.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE READ GROUP KEYS COME FROM THE Arguments AND RecalTable0 TABLES, not from any BAM. The
 *     covariate list is built from a SORTED set of the read groups the report names, so the key a
 *     read group gets depends on the alphabet and not on the file's order;
 *   - AND getReadGroups IS CALLED BEFORE THE TABLES ARE PARSED, so a read group named in RecalTable1
 *     or RecalTable2 but NOT in RecalTable0 gets the missing key -1 and its datum is written at a
 *     negative index, which is an ArrayIndexOutOfBoundsException;
 *   - EVERY DATUM IS BUILT WITH A REPORTED QUALITY OF 1 AND THEN CORRECTED. The constructor takes
 *     `(byte) 1` and `setReportedQuality` overwrites it, which matters because the constructor
 *     refuses a negative quality and the setter refuses NaN and infinity: a report carrying a
 *     negative EstimatedQReported is refused by the setter and not by the constructor;
 *   - THE READ GROUP TABLE READS ITS REPORTED QUALITY FROM EstimatedQReported AND THE OTHERS FROM
 *     QualityScore, which is a different column of a different type: a Double against a Long;
 *   - THE EMPIRICAL QUALITY COLUMN IS IGNORED ON THE WAY IN. Every datum is left uncomputed so it
 *     can be recomputed against a different prior, so a report whose EmpiricalQuality column is
 *     nonsense parses without complaint;
 *   - THE Arguments TABLE IS PARSED BY NAME with a chain of string comparisons, and the value `null`
 *     is turned into a real null first, so `binary_tag_name` of `null` is absent rather than the
 *     four characters;
 *   - NON-STANDARD COVARIATES ARE REFUSED, and so are the two solid_* arguments with any value but
 *     their one allowed one, each with its own message;
 *   - AND THE QUANTIZATION TABLE IS READ BY ROW INDEX AND NOT BY ITS QualityScore COLUMN, so a table
 *     whose rows are out of order is read as though they were in order.
 *
 * Output:
 *
 *     report\t<label>\t<numReadGroups>\t<quantizationLevels>\t<isEmpty>
 *     readgroups\t<label>\t<comma separated, in the order the covariate numbered them>
 *     argument\t<label>\t<name>\t<value>
 *     datum\t<label>\t<table>\t<keys>\t<observations>\t<errors>\t<reported quality bits>
 *     quantized\t<label>\t<comma separated map>
 *     counts\t<label>\t<comma separated counts>
 *     roundtrip\t<label>\t<true|false>
 *     error\t<what>\t<exception>\t<message>
 *
 * Usage: RecalibrationReportDump
 */

import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;
import org.broadinstitute.hellbender.utils.recalibration.EventType;
import org.broadinstitute.hellbender.utils.recalibration.QuantizationInfo;
import org.broadinstitute.hellbender.utils.recalibration.RecalUtils;
import org.broadinstitute.hellbender.utils.recalibration.covariates.ContextCovariate;
import org.broadinstitute.hellbender.utils.recalibration.covariates.CycleCovariate;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.utils.recalibration.RecalDatum;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationReport;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.utils.report.GATKReport;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class RecalibrationReportDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# RecalibrationReportDump: RecalibrationReport");

        for (final Map.Entry<String, String> entry : reports().entrySet()) {
            read(entry.getKey(), entry.getValue());
        }
        errors();
    }

    /**
     * The reports to parse, each shaped for a different decision in the assembly.
     *
     * They are WRITTEN by the reference and then read back by it, because the reader cuts every
     * data line at the header's word starts and a hand-aligned file is cut in the wrong places. The
     * gatk-report-reader suite measures that separately; here the point is the assembly.
     */
    static Map<String, String> reports() {
        final Map<String, String> out = new LinkedHashMap<>();

        // Two read groups, given to the covariate list in the file's order and NOT in alphabetical
        // order, so the keys the report hands back are visibly the sorted ones.
        out.put("two-groups", write(List.of("zebra", "alpha"), true, true));
        // One read group, one datum of every event type.
        out.put("all-events", write(List.of("alpha"), true, true));
        // Nothing beyond the read group table, which is what isEmpty is about.
        out.put("read-group-only", write(List.of("alpha"), false, false));
        return out;
    }

    /**
     * A recalibration report the reference wrote, from tables built here.
     */
    static String write(final List<String> readGroups, final boolean withQualityScores,
                        final boolean withCovariates) {
        final RecalibrationArgumentCollection rac = new RecalibrationArgumentCollection();
        final StandardCovariateList covariates = new StandardCovariateList(rac, readGroups);
        final RecalibrationTables tables = new RecalibrationTables(covariates);

        for (int group = 0; group < readGroups.size(); group++) {
            for (final EventType event : EventType.values()) {
                final RecalDatum datum = new RecalDatum(1000L + group, 10.0 + group, (byte) 1);
                datum.setReportedQuality(30.0 + group + event.ordinal());
                tables.getReadGroupTable().put(datum, group, event.ordinal());
            }
            if (withQualityScores) {
                for (final int quality : new int[] {20, 30}) {
                    for (final EventType event : EventType.values()) {
                        final RecalDatum datum =
                                new RecalDatum(100L + quality, 1.0 + quality / 10.0, (byte) 1);
                        datum.setReportedQuality(quality);
                        tables.getQualityScoreTable().put(datum, group, quality, event.ordinal());
                    }
                }
            }
            if (withCovariates) {
                // One context key and one cycle key, at one reported quality each.
                final RecalDatum context = new RecalDatum(50L, 0.5, (byte) 1);
                context.setReportedQuality(30.0);
                tables.getTable(2).put(context, group, 30, ContextCovariate.keyFromContext("AC"), 0);
                final RecalDatum cycle = new RecalDatum(60L, 0.6, (byte) 1);
                cycle.setReportedQuality(20.0);
                tables.getTable(3).put(cycle, group, 20, CycleCovariate.keyFromCycle(-3, 500), 0);
            }
        }

        final QuantizationInfo quantization = new QuantizationInfo(tables, rac.QUANTIZING_LEVELS);
        final GATKReport report = RecalUtils.createRecalibrationGATKReport(
                rac.generateReportTable(covariates.covariateNames()), quantization, tables,
                covariates);
        return render(report);
    }

    static String render(final GATKReport report) {
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        try (final PrintStream out = new PrintStream(bytes, true, StandardCharsets.UTF_8)) {
            report.print(out);
        }
        return bytes.toString(StandardCharsets.UTF_8);
    }

    /** One report: the covariates it built, every datum it parsed, and its quantization map. */
    static void read(final String label, final String text) {
        final RecalibrationReport report;
        try {
            report = new RecalibrationReport(
                    new ByteArrayInputStream(text.getBytes(StandardCharsets.UTF_8)));
        } catch (final Exception e) {
            System.out.printf("error\tread@%s\t%s\t%s%n", label, e.getClass().getSimpleName(),
                    e.getMessage());
            return;
        }

        final RecalibrationTables tables = report.getRecalibrationTables();
        System.out.printf("report\t%s\t%d\t%d\t%b%n", label,
                report.getCovariates().getReadGroupCovariate().maximumKeyValue() + 1,
                report.getQuantizationInfo().getQuantizationLevels(), report.isEmpty());

        // The order the covariate numbered them, which is the sorted order and not the file's.
        final StringBuilder groups = new StringBuilder();
        for (int key = 0; key <= report.getCovariates().getReadGroupCovariate().maximumKeyValue(); key++) {
            if (groups.length() != 0) {
                groups.append(',');
            }
            groups.append(report.getCovariates().getReadGroupCovariate().formatKey(key));
        }
        System.out.printf("readgroups\t%s\t%s%n", label, groups);

        // The arguments the report round-tripped into the collection.
        final RecalibrationArgumentCollection rac = report.getRAC();
        System.out.printf("argument\t%s\tmismatches_context_size\t%d%n", label, rac.MISMATCHES_CONTEXT_SIZE);
        System.out.printf("argument\t%s\tindels_context_size\t%d%n", label, rac.INDELS_CONTEXT_SIZE);
        System.out.printf("argument\t%s\tmaximum_cycle_value\t%d%n", label, rac.MAXIMUM_CYCLE_VALUE);
        System.out.printf("argument\t%s\tlow_quality_tail\t%d%n", label, rac.LOW_QUAL_TAIL);
        System.out.printf("argument\t%s\tquantizing_levels\t%d%n", label, rac.QUANTIZING_LEVELS);
        System.out.printf("argument\t%s\tmismatches_default_quality\t%d%n", label, rac.MISMATCHES_DEFAULT_QUALITY);
        System.out.printf("argument\t%s\tbinary_tag_name\t%s%n", label, String.valueOf(rac.BINARY_TAG_NAME));
        System.out.printf("argument\t%s\tdefault_platform\t%s%n", label, String.valueOf(rac.DEFAULT_PLATFORM));

        for (int index = 0; index < tables.numTables(); index++) {
            for (final NestedIntegerArray.Leaf<RecalDatum> leaf : tables.getTable(index).getAllLeaves()) {
                System.out.printf("datum\t%s\t%d\t%s\t%d\t%.2f\t%s%n", label, index,
                        join(leaf.keys), leaf.value.getNumObservations(),
                        leaf.value.getNumMismatches(),
                        Long.toHexString(Double.doubleToRawLongBits(leaf.value.getReportedQuality())));
            }
        }

        final StringBuilder map = new StringBuilder();
        for (final byte quantized : report.getQuantizationInfo().getQuantizedQuals()) {
            if (map.length() != 0) {
                map.append(',');
            }
            map.append(quantized);
        }
        System.out.printf("quantized\t%s\t%s%n", label, map);

        // Writing the report back out and reading it again, which is what gathering does.
        try {
            final GATKReport written = report.createGATKReport();
            final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
            try (final PrintStream out = new PrintStream(bytes, true, StandardCharsets.UTF_8)) {
                written.print(out);
            }
            final String second = bytes.toString(StandardCharsets.UTF_8);
            System.out.printf("roundtrip\t%s\t%b%n", label, second.equals(text));
            System.out.printf("rewritten\t%s\t%d%n", label, second.split("\n", -1).length);
        } catch (final Exception e) {
            System.out.printf("error\twrite@%s\t%s\t%s%n", label, e.getClass().getSimpleName(),
                    e.getMessage());
        }
    }

    /** Every report the assembly refuses, made by editing one the reference wrote. */
    static void errors() {
        final String base = write(List.of("alpha"), true, true);

        // A read group named in RecalTable1 but not in RecalTable0: getReadGroups reads only
        // RecalTable0, so this one gets the missing key -1 and is written at a negative index.
        attempt("group-missing-from-read-group-table", base.replaceFirst("(?m)^alpha( +20 +M)", "ghost$1"));

        // A covariate name the list does not know, which is a null dereference.
        attempt("unknown-covariate-name", base.replace("Context", "Nonesuch"));

        // An event type that is not M, I or D.
        attempt("unknown-event-type", base.replaceFirst("(?m)^(alpha +)M", "$1X"));

        // A report with no RecalTable0 at all.
        attempt("no-read-group-table",
                "#:GATKReport.v1.1:1\n"
                        + "#:GATKTable:1:1:%s:;\n"
                        + "#:GATKTable:Quantized:\n"
                        + "QualityScore\n"
                        + "0\n"
                        + "\n");
    }

    static void attempt(final String what, final String text) {
        try {
            final RecalibrationReport report = new RecalibrationReport(
                    new ByteArrayInputStream(text.getBytes(StandardCharsets.UTF_8)));
            System.out.printf("error\t%s\tnone\t%b%n", what, report.isEmpty());
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s\t%s%n", what, e.getClass().getSimpleName(),
                    e.getMessage());
        }
    }

    static String join(final int[] values) {
        final StringBuilder out = new StringBuilder();
        for (final int value : values) {
            if (out.length() != 0) {
                out.append(',');
            }
            out.append(value);
        }
        return out.toString();
    }
}

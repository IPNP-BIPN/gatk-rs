/*
 * AnalyzeCovariates' intermediate csv, taken from the reference.
 *
 * The tool draws plots, and the plots are R's. What the tool itself produces is the csv the R
 * script reads, and that csv is the whole observable of the port: ten columns per row, one row per
 * key of a table built by folding every covariate table into one.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE QUALITY SCORE TABLE IS FILED UNDER AN INDEX THAT IS NOT A COVARIATE'S. The code writes
 *     `newCovs[1] = covariates.size()`, one past the last real index, and calls it a HACK in a
 *     comment; reading it back, `covariateIndex == covariates.size()` is what names the row
 *     "QualityScore". A port that used the covariate's own index would collide with a real one;
 *   - AND THE OPTIONAL COVARIATES DROP THE QUALITY SCORE FROM THEIR KEY, `covs[2] = leaf.keys[2]`
 *     skipping `keys[1]`, so every row of a context or cycle table is SUMMED OVER THE REPORTED
 *     QUALITY. Two data at the same context and different qualities become one row whose
 *     observations are the sum;
 *   - THE HEADER IS PRINTED BY THE PRINTER'S CONSTRUCTOR, so it appears exactly once even though
 *     `writeCsv` takes a `printHeader` flag, which is always false where it is called from;
 *   - THE MODES COME OUT IN THE MAP'S ORDER, BQSR then Before then After, whatever order the
 *     arguments were given in;
 *   - A DATUM IS PRINTED BY String.format WITH FOUR SEPARATE FORMATS: `%d,%.2f,%.2f` from
 *     toString, then the reported quality and the difference at `%.2f` each, so a value that ends
 *     in a 5 is rounded HALF_UP and not to even;
 *   - THE EVENT TYPE IS SPELT OUT, "Base Substitution", "Base Insertion", "Base Deletion", and not
 *     as the single letter the report uses;
 *   - TWO REPORTS WITH DIFFERENT ARGUMENTS ARE REFUSED with a message that names each difference,
 *     joined by "// ", and the key is capitalised. BUT THE CHECK IS PARTLY DEAD: four of its
 *     fifteen comparisons pass the SAME CONSTANT on both sides, `DO_NOT_USE_STANDARD_COVARIATES`
 *     against itself and three more, so they can never fire; and `indels_context_size` is not
 *     compared at all, so two reports built with different indel contexts are combined without a
 *     word;
 *   - AND THE TOOL REFUSES BEFORE IT READS ANYTHING when no report is given or no output is asked
 *     for, which are two different messages naming the short argument names.
 *
 * Output:
 *
 *     report\t<label>\t<the recalibration report, escaped>
 *     csv\t<label>\t<the whole csv, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: AnalyzeCovariatesDump
 */

import org.broadinstitute.hellbender.tools.walkers.bqsr.AnalyzeCovariates;
import org.broadinstitute.hellbender.utils.recalibration.EventType;
import org.broadinstitute.hellbender.utils.recalibration.QuantizationInfo;
import org.broadinstitute.hellbender.utils.recalibration.RecalDatum;
import org.broadinstitute.hellbender.utils.recalibration.RecalUtils;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.utils.recalibration.covariates.ContextCovariate;
import org.broadinstitute.hellbender.utils.recalibration.covariates.CycleCovariate;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.utils.report.GATKReport;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class AnalyzeCovariatesDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("analyzecovariates-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# AnalyzeCovariatesDump: the intermediate csv, from the reference");

        // A report with two data at the same context and different reported qualities, which is
        // what the dropped key aggregates, and one read group.
        final Path plain = write(dir, "plain", List.of("alpha"), 1);
        // Two read groups, so the read group column changes and the keys are the sorted ones.
        final Path groups = write(dir, "two-groups", List.of("zebra", "alpha"), 1);
        // The same shape with different counts, which is what a "before" and an "after" differ by.
        final Path after = write(dir, "after", List.of("alpha"), 7);

        run(dir, "one-report", "-bqsr", plain.toString());
        run(dir, "two-groups", "-bqsr", groups.toString());
        // Before and after together, and then given in the other order, to show the map's order
        // rather than the command line's.
        run(dir, "before-after", "-before", plain.toString(), "-after", after.toString());
        run(dir, "after-before", "-after", after.toString(), "-before", plain.toString());
        // All three roles at once.
        run(dir, "all-three", "-bqsr", plain.toString(), "-before", plain.toString(),
                "-after", after.toString());

        // A report whose arguments differ from the other's, which is the consistency refusal.
        final Path different = write(dir, "different-arguments", List.of("alpha"), 1,
                rac -> rac.MISMATCHES_CONTEXT_SIZE = 3);
        run(dir, "inconsistent", "-before", plain.toString(), "-after", different.toString());

        // And one the check never looks at: the indel context size is missing from the fifteen
        // comparisons, so the two reports are combined without a word.
        final Path unchecked = write(dir, "unchecked-argument", List.of("alpha"), 1,
                rac -> rac.INDELS_CONTEXT_SIZE = 2);
        run(dir, "unchecked-argument", "-before", plain.toString(), "-after", unchecked.toString());

        // The two refusals that happen before anything is read.
        runRefused(dir, "no-report", "-csv", dir.resolve("unused.csv").toString());
        runRefused(dir, "no-output", "-bqsr", plain.toString());
        runRefused(dir, "missing-input", "-bqsr", dir.resolve("absent.table").toString(),
                "-csv", dir.resolve("unused.csv").toString());
    }

    /** A recalibration report the reference wrote, from tables built here. */
    static Path write(final Path dir, final String label, final List<String> readGroups,
                      final long scale) throws Exception {
        return write(dir, label, readGroups, scale, rac -> { });
    }

    static Path write(final Path dir, final String label, final List<String> readGroups,
                      final long scale,
                      final java.util.function.Consumer<RecalibrationArgumentCollection> change)
            throws Exception {
        final RecalibrationArgumentCollection rac = new RecalibrationArgumentCollection();
        change.accept(rac);
        final StandardCovariateList covariates = new StandardCovariateList(rac, readGroups);
        final RecalibrationTables tables = new RecalibrationTables(covariates);

        for (int group = 0; group < readGroups.size(); group++) {
            for (final EventType event : EventType.values()) {
                final RecalDatum datum = new RecalDatum(1000L * scale + group, 10.0 + group, (byte) 1);
                datum.setReportedQuality(30.0 + group + event.ordinal());
                tables.getReadGroupTable().put(datum, group, event.ordinal());
            }
            for (final int quality : new int[] {20, 30}) {
                for (final EventType event : EventType.values()) {
                    final RecalDatum datum =
                            new RecalDatum(100L * scale + quality, 1.0 + quality / 10.0, (byte) 1);
                    datum.setReportedQuality(quality);
                    tables.getQualityScoreTable().put(datum, group, quality, event.ordinal());
                }
            }
            // The SAME context at two different reported qualities, which the csv folds into one
            // row, and one cycle key beside it.
            for (final int quality : new int[] {20, 30}) {
                final RecalDatum context = new RecalDatum(50L * scale + quality, 0.5, (byte) 1);
                context.setReportedQuality(quality);
                tables.getTable(2).put(context, group, quality, ContextCovariate.keyFromContext("AC"), 0);
            }
            final RecalDatum cycle = new RecalDatum(60L * scale, 0.6, (byte) 1);
            cycle.setReportedQuality(20.0);
            tables.getTable(3).put(cycle, group, 20, CycleCovariate.keyFromCycle(-3, 500), 0);
        }

        final QuantizationInfo quantization = new QuantizationInfo(tables, rac.QUANTIZING_LEVELS);
        final GATKReport report = RecalUtils.createRecalibrationGATKReport(
                rac.generateReportTable(covariates.covariateNames()), quantization, tables,
                covariates);
        final Path file = dir.resolve(label + ".table");
        Files.writeString(file, render(report), StandardCharsets.UTF_8);
        System.out.printf("report\t%s\t%s%n", label,
                ReferenceQueryDump.escape(Files.readString(file)));
        return file;
    }

    /** Run the tool for its csv, which is the only output measured here. */
    static void run(final Path dir, final String label, final String... arguments) {
        final Path csv = dir.resolve(label + ".csv");
        final List<String> all = new ArrayList<>(List.of(arguments));
        all.add("-csv");
        all.add(csv.toString());
        try {
            new AnalyzeCovariates().instanceMain(all.toArray(new String[0]));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        try {
            System.out.printf("csv\t%s\t%s%n", label,
                    ReferenceQueryDump.escape(Files.readString(csv)));
        } catch (final Exception e) {
            System.out.printf("error\t%s-read\t%s:%s%n", label, e.getClass().getName(),
                    String.valueOf(e.getMessage()));
        }
    }

    /** Run the tool with arguments it refuses, taking the message exactly as it comes. */
    static void runRefused(final Path dir, final String label, final String... arguments) {
        try {
            new AnalyzeCovariates().instanceMain(arguments);
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("error\t%s\tnone%n", label);
    }

    static String render(final GATKReport report) {
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        try (final PrintStream stream = new PrintStream(bytes, true, StandardCharsets.UTF_8)) {
            report.print(stream);
        }
        return bytes.toString(StandardCharsets.UTF_8);
    }

    static void emptyDirectory(final Path dir) throws Exception {
        if (!Files.isDirectory(dir)) {
            return;
        }
        try (final var entries = Files.list(dir)) {
            for (final Path entry : entries.toList()) {
                Files.deleteIfExists(entry);
            }
        }
    }
}

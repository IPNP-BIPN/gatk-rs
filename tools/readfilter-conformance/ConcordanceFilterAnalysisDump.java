/*
 * Concordance's filter-analysis table, taken from the reference.
 *
 * The second output of the tool: one row per FILTER declared in the EVAL header, counting how many
 * filtered eval records that filter accounts for, split by whether truth had anything at the locus
 * and by whether the filter was the only one on the record.
 *
 * Four behaviours this is built to catch.
 *
 *   - THE COUNTING HAPPENS WITHOUT THE FLAG THAT ASKS FOR IT. The guard is
 *
 *         if (filterAnalysis != null && concordanceState == FILTERED_TRUE_NEGATIVE || concordanceState == FILTERED_FALSE_NEGATIVE)
 *
 *     which parses as `(filterAnalysis != null && FTN) || FFN`, so a FILTERED_FALSE_NEGATIVE walks
 *     into the map on every run and a FILTERED_TRUE_NEGATIVE only on a run that asked for the table.
 *     The counters are invisible when no table is written, so the difference shows as a CRASH: the
 *     lookup is `filterAnalysisRecords::get` keyed by the eval header's own FILTER lines, and a
 *     filter the header does not declare is a null that is then incremented. The same undeclared
 *     filter is fatal on an eval record at a truth locus with no `--filter-analysis` anywhere on the
 *     command line, and harmless on an eval record standing alone;
 *   - A RECORD WITH TWO FILTERS INCREMENTS NEITHER UNIQUE COLUMN, `unique` being
 *     `filters.size() == 1`, a property of the record, computed once and handed to every filter on
 *     it;
 *   - A DECLARED FILTER NOTHING CARRIES STILL GETS A ROW, the map being filled from the header
 *     before the traversal starts;
 *   - AND THE ROW ORDER IS A HashMap'S. The records come out of `filterAnalysisRecords.values()`, so
 *     the order is the bucket order of the filter names' `String.hashCode`, not the header's:
 *     `weak, shallow, noisy, unused` in the header come out `shallow, unused, noisy, weak`.
 *
 * Output:
 *
 *     input\t<label>\t<the whole vcf, escaped>
 *     summary\t<label>\t<the whole summary table, escaped>
 *     table\t<label>\t<the whole filter-analysis table, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: ConcordanceFilterAnalysisDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class ConcordanceFilterAnalysisDump {

    /*
     * Four filters, declared in an order the HashMap does not keep, and one of them carried by no
     * record at all.
     */
    static final String HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##FILTER=<ID=weak,Description=\"the first declared\">\n"
                    + "##FILTER=<ID=shallow,Description=\"the second\">\n"
                    + "##FILTER=<ID=noisy,Description=\"the third\">\n"
                    + "##FILTER=<ID=unused,Description=\"declared and carried by nothing\">\n"
                    + "##contig=<ID=chr1,length=1000>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("concordance-filter-analysis-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# ConcordanceFilterAnalysisDump: what each filter cost, and what counts it");

        // Two filtered records at a truth locus and two standing alone, one of each carrying two
        // filters, plus one plain agreement so the summary is not empty.
        final Path truth = writeVcf(dir, "truth", HEADER,
                "chr1\t100\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t200\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t300\t.\tA\tC\t50\tPASS\t.\tGT\t0/1");
        final Path eval = writeVcf(dir, "eval", HEADER,
                "chr1\t100\t.\tA\tC\t50\tweak\t.\tGT\t0/1",
                "chr1\t200\t.\tA\tC\t50\tweak;shallow\t.\tGT\t0/1",
                "chr1\t300\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t400\t.\tA\tC\t50\tweak\t.\tGT\t0/1",
                "chr1\t500\t.\tA\tC\t50\tshallow;noisy\t.\tGT\t0/1");

        // A filter the header never declares, on a record with nothing at its locus: a filtered true
        // negative, which the guard only lets through when the table was asked for.
        final Path ghostTruth = writeVcf(dir, "ghost-truth", HEADER,
                "chr1\t100\t.\tA\tC\t50\tPASS\t.\tGT\t0/1");
        final Path ghostEval = writeVcf(dir, "ghost-eval", HEADER,
                "chr1\t400\t.\tA\tC\t50\tghost\t.\tGT\t0/1");

        // The same undeclared filter on a record at a truth locus: a filtered false negative, which
        // the guard lets through whatever the command line said.
        final Path ghostLocusTruth = writeVcf(dir, "ghost-locus-truth", HEADER,
                "chr1\t100\t.\tA\tC\t50\tPASS\t.\tGT\t0/1");
        final Path ghostLocusEval = writeVcf(dir, "ghost-locus-eval", HEADER,
                "chr1\t100\t.\tA\tC\t50\tghost\t.\tGT\t0/1");

        run(dir, "baseline", truth, eval, true);
        run(dir, "ftn-undeclared-no-flag", ghostTruth, ghostEval, false);
        run(dir, "ftn-undeclared-with-flag", ghostTruth, ghostEval, true);
        run(dir, "ffn-undeclared-no-flag", ghostLocusTruth, ghostLocusEval, false);
    }

    static Path writeVcf(final Path dir, final String label, final String header,
                         final String... records) throws Exception {
        final StringBuilder text = new StringBuilder(header);
        for (final String record : records) {
            text.append(record).append("\n");
        }
        final Path file = dir.resolve(label + ".vcf");
        Files.writeString(file, text.toString(), StandardCharsets.UTF_8);
        new IndexFeatureFile().instanceMain(new String[] {"-I", file.toString()});
        System.out.printf("input\t%s\t%s%n", label, ReferenceQueryDump.escape(text.toString()));
        return file;
    }

    static void run(final Path dir, final String label, final Path truth, final Path eval,
                    final boolean withFilterAnalysis) {
        final Path summary = dir.resolve(label + ".summary");
        final Path analysis = dir.resolve(label + ".filters");
        final List<String> all = new ArrayList<>(List.of(
                "--truth", truth.toString(),
                "--evaluation", eval.toString(),
                "--summary", summary.toString()));
        if (withFilterAnalysis) {
            all.add("--filter-analysis");
            all.add(analysis.toString());
        }
        try {
            new Concordance().instanceMain(all.toArray(new String[0]));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        print("summary", label, summary);
        if (withFilterAnalysis) {
            print("table", label, analysis);
        }
    }

    static void print(final String kind, final String label, final Path file) {
        try {
            System.out.printf("%s\t%s\t%s%n", kind, label,
                    ReferenceQueryDump.escape(Files.readString(file, StandardCharsets.UTF_8)));
        } catch (final Exception e) {
            System.out.printf("error\t%s-%s\t%s:%s%n", label, kind, e.getClass().getName(),
                    String.valueOf(e.getMessage()));
        }
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

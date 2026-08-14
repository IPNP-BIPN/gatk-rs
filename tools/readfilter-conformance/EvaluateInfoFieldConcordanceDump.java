/*
 * EvaluateInfoFieldConcordance's summary table, taken from the reference.
 *
 * The first tool written on the concordance iterator: it walks a truth VCF and an eval VCF, and for
 * every true positive takes the difference between one INFO field of each, reporting a mean and a
 * standard deviation per variant type.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE STANDARD DEVIATION CAN BE NaN BY ARITHMETIC ALONE. The variance is
 *     `(sumSq - sum * sum / n) / n`, the cancelling form, so a set of equal deltas can produce a
 *     small negative variance whose square root is NaN. Nothing guards it;
 *   - AN EMPTY BUCKET IS TWO NaN COLUMNS, since n = 0 makes both the mean and the variance 0/0: a
 *     run with no indel among its true positives still writes an INDEL row, and that row is NaN;
 *   - THE MEAN IS OF ABSOLUTE DIFFERENCES BUT COMPUTED AS sqrt(delta * delta), which is not the
 *     same double as Math.abs(delta) for every input;
 *   - AND THE VARIANCE IS OF THOSE ABSOLUTE VALUES while the sum of squares is of the signed ones,
 *     which is the same number only because squaring drops the sign;
 *   - A RECORD WHOSE KEY IS ABSENT IS COUNTED BUT CONTRIBUTES NOTHING: the counter is incremented
 *     in `apply` and the delta only inside `infoDifference`, which returns early when either side
 *     lacks the key, so the mean is divided by a count larger than the number of deltas;
 *   - ONLY TRUE POSITIVES ARE LOOKED AT: false positives, false negatives and both filtered states
 *     fall through the switch untouched, and this walker's own filters make the filtered ones
 *     unreachable anyway;
 *   - AND A MISSING KEY IN EITHER HEADER IS A REFUSAL BEFORE ANY RECORD, worded with the file it
 *     was not found in.
 *
 * Output:
 *
 *     input\t<label>\t<the whole vcf, escaped>
 *     table\t<label>\t<the whole summary table, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: EvaluateInfoFieldConcordanceDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.validation.EvaluateInfoFieldConcordance;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class EvaluateInfoFieldConcordanceDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##INFO=<ID=SCORE,Number=1,Type=Float,Description=\"the key both sides use\">\n"
                    + "##INFO=<ID=OTHER,Number=1,Type=Float,Description=\"a second key\">\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##FILTER=<ID=weak,Description=\"Was already there\">\n"
                    + "##contig=<ID=chr1,length=1000>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n";

    /** The same header without SCORE, for the two refusals. */
    static final String HEADER_WITHOUT_SCORE =
            HEADER.replace("##INFO=<ID=SCORE,Number=1,Type=Float,Description=\"the key both sides use\">\n", "");

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("evaluateinfofieldconcordance-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# EvaluateInfoFieldConcordanceDump: two rows, and the arithmetic behind them");

        // Four true positives with deltas of 1.0, 1.0, 1.0 and 1.0: equal deltas are what make the
        // cancelling variance interesting, since the answer should be exactly zero.
        final Path truth = writeVcf(dir, "truth", HEADER,
                "chr1\t100\t.\tA\tC\t50\tPASS\tSCORE=1.0\tGT\t0/1",
                "chr1\t200\t.\tA\tC\t50\tPASS\tSCORE=2.0\tGT\t0/1",
                "chr1\t300\t.\tA\tC\t50\tPASS\tSCORE=3.0\tGT\t0/1",
                "chr1\t400\t.\tA\tC\t50\tPASS\tSCORE=4.0\tGT\t0/1",
                // An indel true positive, so the INDEL row is not empty in the first run.
                "chr1\t500\t.\tACC\tA\t50\tPASS\tSCORE=1.0\tGT\t0/1",
                // A record whose key is absent on the truth side.
                "chr1\t600\t.\tA\tC\t50\tPASS\tOTHER=9.0\tGT\t0/1",
                // A truth record with no eval at all.
                "chr1\t700\t.\tA\tC\t50\tPASS\tSCORE=1.0\tGT\t0/1");
        final Path eval = writeVcf(dir, "eval", HEADER,
                "chr1\t100\t.\tA\tC\t50\tPASS\tSCORE=2.0\tGT\t0/1",
                "chr1\t200\t.\tA\tC\t50\tPASS\tSCORE=3.0\tGT\t0/1",
                "chr1\t300\t.\tA\tC\t50\tPASS\tSCORE=4.0\tGT\t0/1",
                "chr1\t400\t.\tA\tC\t50\tPASS\tSCORE=5.0\tGT\t0/1",
                "chr1\t500\t.\tACC\tA\t50\tPASS\tSCORE=2.5\tGT\t0/1",
                "chr1\t600\t.\tA\tC\t50\tPASS\tSCORE=1.0\tGT\t0/1",
                // An eval record with no truth at all, and a filtered one.
                "chr1\t800\t.\tA\tC\t50\tPASS\tSCORE=1.0\tGT\t0/1",
                "chr1\t900\t.\tA\tC\t50\tweak\tSCORE=1.0\tGT\t0/1");

        // Deltas that do not cancel, so the standard deviation is a real number.
        final Path spread = writeVcf(dir, "spread", HEADER,
                "chr1\t100\t.\tA\tC\t50\tPASS\tSCORE=1.0\tGT\t0/1",
                "chr1\t200\t.\tA\tC\t50\tPASS\tSCORE=1.0\tGT\t0/1",
                "chr1\t300\t.\tA\tC\t50\tPASS\tSCORE=1.0\tGT\t0/1");
        final Path spreadEval = writeVcf(dir, "spread-eval", HEADER,
                "chr1\t100\t.\tA\tC\t50\tPASS\tSCORE=1.5\tGT\t0/1",
                "chr1\t200\t.\tA\tC\t50\tPASS\tSCORE=0.25\tGT\t0/1",
                "chr1\t300\t.\tA\tC\t50\tPASS\tSCORE=101.0\tGT\t0/1");

        // Nothing agrees, so both rows are empty and therefore NaN.
        final Path nothingAgrees = writeVcf(dir, "nothing-agrees", HEADER,
                "chr1\t100\t.\tA\tC\t50\tPASS\tSCORE=1.0\tGT\t0/1");
        final Path nothingAgreesEval = writeVcf(dir, "nothing-agrees-eval", HEADER,
                "chr1\t100\t.\tA\tG\t50\tPASS\tSCORE=2.0\tGT\t0/1");

        // A header that never declares the key.
        final Path noKey = writeVcf(dir, "no-key", HEADER_WITHOUT_SCORE,
                "chr1\t100\t.\tA\tC\t50\tPASS\tOTHER=1.0\tGT\t0/1");

        run(dir, "baseline", truth, eval, "SCORE", "SCORE", "baseline.table");
        run(dir, "spread", spread, spreadEval, "SCORE", "SCORE", "spread.table");
        run(dir, "nothing-agrees", nothingAgrees, nothingAgreesEval, "SCORE", "SCORE",
                "nothing-agrees.table");
        // The two keys need not be the same one.
        run(dir, "different-keys", truth, eval, "SCORE", "OTHER", "different-keys.table");
        run(dir, "missing-eval-key", truth, noKey, "SCORE", "SCORE", "missing-eval-key.table");
        run(dir, "missing-truth-key", noKey, eval, "SCORE", "SCORE", "missing-truth-key.table");
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
                    final String evalKey, final String truthKey, final String output) {
        final Path file = dir.resolve(output);
        final List<String> all = new ArrayList<>(List.of(
                "--truth", truth.toString(),
                "--evaluation", eval.toString(),
                "--summary", file.toString(),
                "--eval-info-key", evalKey,
                "--truth-info-key", truthKey));
        try {
            new EvaluateInfoFieldConcordance().instanceMain(all.toArray(new String[0]));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        try {
            System.out.printf("table\t%s\t%s%n", label,
                    ReferenceQueryDump.escape(Files.readString(file, StandardCharsets.UTF_8)));
        } catch (final Exception e) {
            System.out.printf("error\t%s-read\t%s:%s%n", label, e.getClass().getName(),
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

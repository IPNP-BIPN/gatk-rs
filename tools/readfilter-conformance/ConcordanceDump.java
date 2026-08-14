/*
 * Concordance's summary table, taken from the reference.
 *
 * The tool the concordance iterator was written for: a truth VCF and an eval VCF walked in lockstep,
 * every step labelled with one of five states, and the states tallied into a six-column table with
 * one row for SNPs and one for everything else. Only `--summary` is asked for here; the
 * filter-analysis table and the three annotated VCFs are measured separately.
 *
 * Six behaviours this is built to catch.
 *
 *   - AN EMPTY CALLSET REPORTS ZERO RATHER THAN NaN. Both rates are long divisions, so 0/0 is NaN,
 *     but they are printed through `MathUtils.roundToNDecimalPlaces(x, 3)`, which is
 *     `Math.round((x + Math.ulp(x)) * 1000) / 1000`, and `Math.round` of NaN is 0. The same 0/0 that
 *     `EvaluateInfoFieldConcordance` writes as `NaN` is written here as `0.0`;
 *   - A FILTERED EVAL RECORD ALONE LEAVES NO TRACE. FILTERED_TRUE_NEGATIVE is counted into the enum
 *     map and never read again: the summary's FN column is FALSE_NEGATIVE + FILTERED_FALSE_NEGATIVE
 *     and its FP column is FALSE_POSITIVE alone, so an unmatched filtered eval record moves neither
 *     rate;
 *   - EVERYTHING THAT IS NOT A SNP IS AN INDEL. The stratification is a single `isSNP()` on
 *     `getTruthIfPresentElseEval()`, so an MNP and a symbolic record both land in the INDEL row;
 *   - THE TRUTH SIDE DROPS SYMBOLIC AND SV RECORDS AND THE EVAL SIDE KEEPS EVERYTHING, the truth
 *     filter being `!isFiltered() && !isSymbolicOrSV()` against a default eval filter of `vc -> true`,
 *     so the same symbolic record in both files is a false positive rather than a true positive;
 *   - AGREEMENT NEEDS THE SAME NUMBER OF ALTERNATES BUT ONLY TRUTH'S FIRST. A truth `A/C` against an
 *     eval `A/C,G` disagrees on the count, while a truth `A/C,G` against an eval `A/G,C` agrees,
 *     because `eval.getAlternateAlleles().contains(truth.getAlternateAllele(0))` asks for one allele
 *     and the size test asks for nothing else;
 *   - AND THE ROUNDING IS HALF-UP AFTER AN ULP IS ADDED, so 2/3 prints `0.667`.
 *
 * Output:
 *
 *     input\t<label>\t<the whole vcf, escaped>
 *     table\t<label>\t<the whole summary table, escaped>
 *
 * Usage: ConcordanceDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

public class ConcordanceDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##ALT=<ID=DEL,Description=\"a deletion, symbolically\">\n"
                    + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"the end a symbolic allele needs\">\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##FILTER=<ID=weak,Description=\"Was already there\">\n"
                    + "##contig=<ID=chr1,length=1000>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("concordance-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# ConcordanceDump: six columns, and what each state does to them");

        // Two true positives, one false negative and one false positive on each row, so both rates
        // are 2/3 on the SNP row and 1/2 on the INDEL row.
        final Path truth = writeVcf(dir, "truth", HEADER,
                "chr1\t100\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t200\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t300\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t400\t.\tACC\tA\t50\tPASS\t.\tGT\t0/1",
                "chr1\t500\t.\tACC\tA\t50\tPASS\t.\tGT\t0/1");
        final Path eval = writeVcf(dir, "eval", HEADER,
                "chr1\t100\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t200\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t400\t.\tACC\tA\t50\tPASS\t.\tGT\t0/1",
                "chr1\t600\t.\tA\tG\t50\tPASS\t.\tGT\t0/1",
                "chr1\t700\t.\tACC\tA\t50\tPASS\t.\tGT\t0/1");

        // A filtered eval record at a truth locus and a filtered eval record alone: the first is a
        // filtered false negative and reaches the FN column, the second is a filtered true negative
        // and reaches nothing, so the precision stays 1.0 with an unmatched eval record in the file.
        final Path filteredTruth = writeVcf(dir, "filtered-truth", HEADER,
                "chr1\t100\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t200\t.\tA\tC\t50\tPASS\t.\tGT\t0/1");
        final Path filteredEval = writeVcf(dir, "filtered-eval", HEADER,
                "chr1\t100\t.\tA\tC\t50\tweak\t.\tGT\t0/1",
                "chr1\t200\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t300\t.\tA\tC\t50\tweak\t.\tGT\t0/1");

        // Nothing on either side: four zeroes and four 0/0 divisions.
        final Path emptyTruth = writeVcf(dir, "empty-truth", HEADER);
        final Path emptyEval = writeVcf(dir, "empty-eval", HEADER);

        // An MNP, a symbolic record and a multi-allelic SNP. The symbolic truth record is dropped by
        // the truth filter and the symbolic eval record is not, so it comes out a false positive on
        // the row that is not the SNP row.
        final Path typesTruth = writeVcf(dir, "types-truth", HEADER,
                "chr1\t100\t.\tAT\tGC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t200\t.\tA\t<DEL>\t50\tPASS\tEND=250\tGT\t0/1",
                "chr1\t300\t.\tA\tC,G\t50\tPASS\t.\tGT\t1/2");
        final Path typesEval = writeVcf(dir, "types-eval", HEADER,
                "chr1\t100\t.\tAT\tGC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t200\t.\tA\t<DEL>\t50\tPASS\tEND=250\tGT\t0/1",
                "chr1\t300\t.\tA\tG,C\t50\tPASS\t.\tGT\t1/2");

        // The allele rule: a count that differs, an order that does not matter, and two indels whose
        // reference alleles differ. The first disagreement is what puts the two sides out of step.
        final Path allelesTruth = writeVcf(dir, "alleles-truth", HEADER,
                "chr1\t100\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t200\t.\tA\tC,G\t50\tPASS\t.\tGT\t1/2",
                "chr1\t300\t.\tAT\tA\t50\tPASS\t.\tGT\t0/1");
        final Path allelesEval = writeVcf(dir, "alleles-eval", HEADER,
                "chr1\t100\t.\tA\tC,G\t50\tPASS\t.\tGT\t1/2",
                "chr1\t200\t.\tA\tG,C\t50\tPASS\t.\tGT\t1/2",
                "chr1\t300\t.\tACC\tA\t50\tPASS\t.\tGT\t0/1");

        run(dir, "baseline", truth, eval, "baseline.table");
        run(dir, "filtered", filteredTruth, filteredEval, "filtered.table");
        run(dir, "empty", emptyTruth, emptyEval, "empty.table");
        run(dir, "types", typesTruth, typesEval, "types.table");
        run(dir, "alleles", allelesTruth, allelesEval, "alleles.table");
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
                    final String output) {
        final Path file = dir.resolve(output);
        final List<String> all = List.of(
                "--truth", truth.toString(),
                "--evaluation", eval.toString(),
                "--summary", file.toString());
        try {
            new Concordance().instanceMain(all.toArray(new String[0]));
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

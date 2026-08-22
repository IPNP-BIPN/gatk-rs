/*
 * MergeMutect2CallsWithMC3's output, taken from the reference.
 *
 * Mutect2 calls merged into the MC3 pan-cancer call set, one concordance state at a time. The
 * walking is `AbstractConcordanceWalker`'s and is already measured; what is measured here is what
 * each of the five states writes, which is five different records.
 *
 * Twelve behaviours this is built to catch.
 *
 *   - EACH STATE WRITES A DIFFERENT RECORD, and one writes nothing at all: a FILTERED_TRUE_NEGATIVE
 *     is skipped outright, so an eval-only filtered call disappears without a trace;
 *   - A FALSE POSITIVE IS REBUILT FROM SCRATCH, `new VariantContextBuilder(source, contig, start,
 *     end, alleles)`, so every INFO field, the ID, the QUAL and the FILTER column of the M2 record
 *     are DROPPED and only CENTERS=M2 survives;
 *   - A TRUE POSITIVE AND A FILTERED FALSE NEGATIVE KEEP THE MC3 RECORD WHOLE and only add to it,
 *     so the two sides are treated asymmetrically: MC3's annotations are authoritative and M2's are
 *     discarded;
 *   - CENTERS IS APPENDED TO, not replaced, and an absent CENTERS becomes a list of one, but ONLY
 *     for the three states that add it: a FALSE NEGATIVE is emitted unchanged except for its
 *     genotype and never learns that M2 looked at it;
 *   - THE GENOTYPE'S PLOIDY IS THE NUMBER OF ALLELES AT THE SITE, `new GenotypeBuilder(sample,
 *     variant.getAlleles())` being handed every allele rather than a called pair, so a
 *     multiallelic false positive comes out `0/1/2`;
 *   - AN M2 GENOTYPE WITHOUT AD LEAVES THE OUTPUT GENOTYPE WITHOUT ONE, `getAD()` answering null
 *     and `GenotypeBuilder.AD(null)` setting nothing, rather than throwing or writing zeroes;
 *   - THE GENOTYPE'S ALLELES COME FROM `getTruthIfPresentElseEval()`, so a false positive's
 *     genotype is built from the M2 alleles and everything else from MC3's;
 *   - THE ALLELE DEPTHS COME FROM M2 WHEN M2 IS THERE, and from MC3's NREF and NALT INFO fields
 *     when it is not, each defaulting to ZERO when absent;
 *   - THE OUTPUT CARRIES ONE SAMPLE, named by the EVAL header's tumor_sample line, whatever the two
 *     inputs' own sample columns say;
 *   - AND AN EVAL HEADER WITHOUT THAT LINE IS A NullPointerException rather than a message;
 *   - THE HEADER IS THE TRUTH'S metadata plus the standard GT and AD format lines, the tool's own
 *     lines and M2_FILTERS, so the eval header's own INFO lines never reach the output;
 *   - AND M2_FILTERS IS WRITTEN ONLY FOR A FILTERED FALSE NEGATIVE, as the list of M2's filters.
 *
 * Output:
 *
 *     input\t<label>=<the whole vcf, escaped>
 *     merged\t<label>=<the whole output vcf, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: MergeMutect2WithMC3Dump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.validation.MergeMutect2CallsWithMC3;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class MergeMutect2WithMC3Dump {

    static final String CONTIG = "##contig=<ID=chr1,length=10000>\n";

    /** The MC3 side: NREF and NALT carry the depths, CENTERS the calling centres. */
    static final String TRUTH_HEADER =
            "##fileformat=VCFv4.2\n"
            + "##FILTER=<ID=weak,Description=\"Weak\">\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            + "##INFO=<ID=CENTERS,Number=.,Type=String,Description=\"Calling centers\">\n"
            + "##INFO=<ID=NREF,Number=1,Type=Integer,Description=\"Reference count\">\n"
            + "##INFO=<ID=NALT,Number=1,Type=Integer,Description=\"Alternate count\">\n"
            + "##INFO=<ID=MC3ONLY,Number=1,Type=String,Description=\"An MC3 annotation\">\n"
            + CONTIG
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmc3sample\n";

    /** The M2 side, whose header names the tumour sample. */
    static final String EVAL_HEADER =
            "##fileformat=VCFv4.2\n"
            + "##FILTER=<ID=weak_evidence,Description=\"Weak evidence\">\n"
            + "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            + "##INFO=<ID=M2ONLY,Number=1,Type=String,Description=\"An M2 annotation\">\n"
            + "##tumor_sample=tumour\n"
            + CONTIG
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumour\n";

    /** The same, with no tumor_sample line at all. */
    static final String EVAL_HEADER_NO_SAMPLE = EVAL_HEADER.replace("##tumor_sample=tumour\n", "");

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("merge-mutect2-mc3-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# MergeMutect2WithMC3Dump: M2 calls merged into the MC3 call set");

        // One record per concordance state, plus the variations that matter within each.
        //
        //   100  both, M2 unfiltered            -> true positive
        //   110  MC3 only                       -> false negative
        //   120  M2 only, unfiltered            -> false positive
        //   130  M2 only, filtered              -> filtered true negative, written nowhere
        //   140  both, M2 filtered              -> filtered false negative, M2_FILTERS written
        //   150  both, MC3 already has CENTERS  -> the list is appended to
        //   160  both, MC3 without NREF or NALT -> the depths come from M2 anyway
        //   170  MC3 only, without NREF or NALT  -> both default to zero
        final String truth = TRUTH_HEADER
                + "chr1\t100\t.\tA\tC\t50\tPASS\tNREF=80;NALT=20;MC3ONLY=kept\tGT\t0/1\n"
                + "chr1\t110\t.\tA\tC\t50\tPASS\tNREF=70;NALT=30\tGT\t0/1\n"
                + "chr1\t140\t.\tA\tC\t50\tPASS\tNREF=60;NALT=40\tGT\t0/1\n"
                + "chr1\t150\t.\tA\tC\t50\tPASS\tCENTERS=broad,wustl;NREF=50;NALT=50\tGT\t0/1\n"
                + "chr1\t160\t.\tA\tC\t50\tPASS\t.\tGT\t0/1\n"
                + "chr1\t170\t.\tA\tC\t50\tPASS\t.\tGT\t0/1\n";
        final String eval = EVAL_HEADER
                + "chr1\t100\t.\tA\tC\t50\tPASS\tM2ONLY=dropped\tGT:AD\t0/1:11,12\n"
                + "chr1\t120\trs1\tA\tC\t99\tPASS\tM2ONLY=dropped\tGT:AD\t0/1:13,14\n"
                + "chr1\t130\t.\tA\tC\t50\tweak_evidence\t.\tGT:AD\t0/1:15,16\n"
                + "chr1\t140\t.\tA\tC\t50\tweak_evidence\t.\tGT:AD\t0/1:17,18\n"
                + "chr1\t150\t.\tA\tC\t50\tPASS\t.\tGT:AD\t0/1:19,20\n"
                + "chr1\t160\t.\tA\tC\t50\tPASS\t.\tGT:AD\t0/1:21,22\n";
        run(dir, "every-state", truth, eval);

        // An M2 record whose genotype has no AD at all, against a truth record that does.
        run(dir, "eval-without-ad",
                TRUTH_HEADER + "chr1\t100\t.\tA\tC\t50\tPASS\tNREF=80;NALT=20\tGT\t0/1\n",
                EVAL_HEADER + "chr1\t100\t.\tA\tC\t50\tPASS\t.\tGT\t0/1\n");

        // A multiallelic M2 record containing the MC3 alternate, which is concordant, and one that
        // does not contain it, which is not.
        run(dir, "multiallelic",
                TRUTH_HEADER + "chr1\t100\t.\tA\tC\t50\tPASS\tNREF=80;NALT=20\tGT\t0/1\n"
                + "chr1\t200\t.\tA\tC\t50\tPASS\tNREF=70;NALT=30\tGT\t0/1\n",
                EVAL_HEADER + "chr1\t100\t.\tA\tC,G\t50\tPASS\t.\tGT:AD\t0/1:10,20,30\n"
                + "chr1\t200\t.\tA\tG,T\t50\tPASS\t.\tGT:AD\t0/1:10,20,30\n");

        // And an eval header that never names the tumour sample.
        run(dir, "no-tumor-sample",
                TRUTH_HEADER + "chr1\t100\t.\tA\tC\t50\tPASS\tNREF=80;NALT=20\tGT\t0/1\n",
                EVAL_HEADER_NO_SAMPLE + "chr1\t100\t.\tA\tC\t50\tPASS\t.\tGT:AD\t0/1:11,12\n");
    }

    static void run(final Path dir, final String label, final String truth, final String eval)
            throws Exception {
        final Path truthFile = write(dir, label + "-truth", truth);
        final Path evalFile = write(dir, label + "-eval", eval);
        final Path out = dir.resolve(label + "-merged.vcf");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "--truth", truthFile.toString(),
                "--evaluation", evalFile.toString(),
                "-O", out.toString()));
        try {
            new MergeMutect2CallsWithMC3().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            return;
        }
        if (Files.exists(out)) {
            System.out.printf("merged\t%s=%s%n", label,
                    ReferenceQueryDump.escape(masked(Files.readString(out), dir)));
        }
    }

    static Path write(final Path dir, final String name, final String vcf) throws Exception {
        final Path file = dir.resolve(name + ".vcf");
        Files.writeString(file, vcf, StandardCharsets.UTF_8);
        new IndexFeatureFile().instanceMain(new String[] {"-I", file.toString()});
        System.out.printf("input\t%s=%s%n", name, ReferenceQueryDump.escape(vcf));
        return file;
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>")
                .replaceAll("##GATKCommandLine=<[^\n]*>\n", "")
                .replaceAll("##source=[^\n]*\n", "");
    }
}

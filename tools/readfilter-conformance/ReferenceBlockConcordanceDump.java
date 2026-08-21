/*
 * ReferenceBlockConcordance's three histograms, taken from the reference.
 *
 * Two GVCFs walked side by side by AbstractConcordanceWalker, which is already ported, and the
 * reference blocks they hold turned into three metrics files.
 *
 * Nine behaviours this is built to catch.
 *
 *   - ONLY HOM-REF SITES ARE WALKED. Both filters are `isHomRef`, which reads GENOTYPE 0 and
 *     nothing else, so a variant site is dropped from both inputs before the walk begins;
 *   - AND `isHomRef` READS GENOTYPE 0 OF A FILE THAT MAY HAVE MANY, so a multi-sample record is
 *     filtered on its FIRST sample and only then refused by `extractLengthAndGQ`;
 *   - A FILTERED BLOCK IS STILL WALKED: the tool's filters test the genotype alone, so a block
 *     carrying a FILTER lands in the histogram like any other;
 *   - THE BLOCK HISTOGRAMS ARE KEYED BY A JAVA PAIR'S toString, which is `length,GQ` with a comma
 *     and NO SPACE and no brackets, and the metrics file then SORTS THOSE KEYS AS STRINGS: the
 *     golden's truth histogram runs `1,80`, `100,20`, `50,40`, `50,60`, which is neither the
 *     lengths' order nor the file's;
 *   - THE LENGTH IS `getLengthOnReference()`, so a block's END attribute decides it and a block of
 *     one base is length 1;
 *   - THE CONCORDANCE HISTOGRAM IS INCREMENTED BY THE OVERLAP'S LENGTH, not by one, and it fires
 *     ONLY while both sides currently hold a block;
 *   - AND THE CURRENT BLOCKS ARE CLEARED BY AN OVERLAP TEST AGAINST THE STEP ITSELF, so a step
 *     whose truth block ended before it drops that side and no pair is counted;
 *   - A SITE PRESENT ON ONE SIDE ONLY STILL COUNTS TOWARDS THAT SIDE'S BLOCK HISTOGRAM;
 *   - AND THE THREE FILES ARE METRICS FILES WITH NO METRIC ROWS, so each carries a `## HISTOGRAM`
 *     section alone, under a header whose lines are masked here because they carry the run's date.
 *
 * Output:
 *
 *     truth\t<label>=<the truth gvcf, escaped>
 *     eval\t<label>=<the eval gvcf, escaped>
 *     histogram\t<label>\t<which>=<the metrics file, escaped, dates masked>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: ReferenceBlockConcordanceDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.validation.ReferenceBlockConcordance;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

public class ReferenceBlockConcordanceDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
            + "##ALT=<ID=NON_REF,Description=\"Non-reference allele\">\n"
            + "##FILTER=<ID=q10,Description=\"Quality below 10\">\n"
            + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n"
            + "##contig=<ID=chr1,length=100000>\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n";

    /** A reference block: hom-ref with an END and a GQ. */
    static String block(final int start, final int end, final int gq, final String filter) {
        return "chr1\t" + start + "\t.\tA\t<NON_REF>\t.\t" + filter + "\tEND=" + end
                + "\tGT:GQ\t0/0:" + gq + "\n";
    }

    /** A variant site, which both filters drop. */
    static String variant(final int start, final int gq) {
        return "chr1\t" + start + "\t.\tA\tC,<NON_REF>\t50\t.\t.\tGT:GQ\t0/1:" + gq + "\n";
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("reference-block-concordance-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# ReferenceBlockConcordanceDump: the reference blocks of two GVCFs");

        // Truth: three blocks and a variant site the filter drops.
        final String truth = HEADER
                + block(100, 199, 20, ".")
                + variant(200, 99)
                + block(300, 349, 40, ".")
                // A filtered block, which the genotype filter still lets through.
                + block(400, 449, 60, "q10")
                // A block of one base.
                + block(500, 500, 80, ".");
        // Eval: blocks that overlap the truth's in different ways.
        final String eval = HEADER
                // Overlaps the first truth block for 50 bases, with a different GQ.
                + block(150, 249, 30, ".")
                // Sits entirely inside the second, same length as the overlap.
                + block(310, 329, 40, ".")
                // Starts where the truth's filtered block ends, so they share one base.
                + block(449, 460, 70, ".")
                // And one the truth has nothing near, which still counts for the eval histogram.
                + block(700, 799, 90, ".");
        run(dir, "blocks", truth, eval);

        // Two files that hold the same blocks, so every overlap is whole.
        run(dir, "identical", truth, truth);

        // A truth file with no blocks at all.
        run(dir, "empty-truth", HEADER, eval);

        // A multi-sample record, filtered on genotype 0 and then refused by the length extraction.
        final String twoSamples =
                HEADER.replace("\tsample\n", "\tsample\tother\n")
                + "chr1\t100\t.\tA\t<NON_REF>\t.\t.\tEND=199\tGT:GQ\t0/0:20\t0/0:25\n";
        run(dir, "two-samples", twoSamples, twoSamples);
    }

    static void run(final Path dir, final String label, final String truth, final String eval)
            throws Exception {
        final Path truthFile = write(dir, label + "-truth.vcf", truth);
        final Path evalFile = write(dir, label + "-eval.vcf", eval);
        System.out.printf("truth\t%s=%s%n", label, ReferenceQueryDump.escape(truth));
        System.out.printf("eval\t%s=%s%n", label, ReferenceQueryDump.escape(eval));

        final Path truthHistogram = dir.resolve(label + "-truth-blocks.txt");
        final Path evalHistogram = dir.resolve(label + "-eval-blocks.txt");
        final Path concordance = dir.resolve(label + "-concordance.txt");
        try {
            new ReferenceBlockConcordance().instanceMain(new String[] {
                    "--truth", truthFile.toString(),
                    "--eval", evalFile.toString(),
                    "--truth-block-histogram", truthHistogram.toString(),
                    "--eval-block-histogram", evalHistogram.toString(),
                    "--confidence-concordance-histogram", concordance.toString()});
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            return;
        }
        print(label, "truth-blocks", truthHistogram, dir);
        print(label, "eval-blocks", evalHistogram, dir);
        print(label, "concordance", concordance, dir);
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path file = dir.resolve(name);
        Files.writeString(file, text, StandardCharsets.UTF_8);
        new IndexFeatureFile().instanceMain(new String[] {"-I", file.toString()});
        return file;
    }

    /** One metrics file, with the header lines that carry the run's own identity masked. */
    static void print(final String label, final String which, final Path file, final Path dir)
            throws Exception {
        if (!Files.exists(file)) {
            return;
        }
        final StringBuilder text = new StringBuilder();
        for (final String line : Files.readAllLines(file, StandardCharsets.UTF_8)) {
            if (line.startsWith("# ")) {
                // The command line and the start-up date, neither of which is a result.
                text.append("# MASKED\n");
                continue;
            }
            text.append(masked(line, dir)).append("\n");
        }
        System.out.printf("histogram\t%s\t%s=%s%n", label, which,
                ReferenceQueryDump.escape(text.toString()));
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}

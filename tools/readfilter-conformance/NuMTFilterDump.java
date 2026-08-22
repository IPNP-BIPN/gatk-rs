/*
 * NuMTFilterTool's output, taken from the reference.
 *
 * A mitochondrial call filtered when its alternate depth is low enough to be a nuclear insertion of
 * mitochondrial DNA rather than a real heteroplasmy. One hundred and five lines, and most of them
 * are not what the name suggests.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE TOOL DOES NOTHING AT ITS DEFAULTS. `medianAutosomalCoverage` defaults to ZERO and the
 *     cutoff is only computed when BOTH it and `maxNuMTAutosomalCopies` are above zero, so the
 *     cutoff stays 0, `max(AD) < 0` is never true, and the output is the input with one more
 *     header line;
 *   - THE CUTOFF IS A POISSON QUANTILE, `PoissonDistribution(coverage * copies / 2)
 *     .inverseCumulativeProbability(0.99)`, so it is an integer read off a distribution rather than
 *     a depth the user gave;
 *   - THE FILTER IS APPLIED WHEN NO ALTERNATE ESCAPES IT, `!appliedFilter.contains(FALSE)`, and an
 *     EMPTY list contains no FALSE either, so a record with no alternate allele at all is filtered;
 *   - THE AS_FilterStatus ATTRIBUTE IS WRITTEN WHENEVER ANY alternate is filtered, so a
 *     multiallelic record can carry the attribute without carrying the site filter;
 *   - AND WRITING IT THROWS WHEN THE RECORD HAS NO AS_FilterStatus TO MERGE INTO:
 *     `getMergedASFilterString` validates that the decoded list is the same length as the list of
 *     alternates, and an absent attribute decodes to an EMPTY list, so an ordinary VCF that has
 *     never been through Mutect2's filtering makes the tool throw;
 *   - THE DEPTH COMPARED IS THE MAXIMUM ACROSS SAMPLES for that allele, not the sum;
 *   - A GENOTYPE WITHOUT AD IS SKIPPED ENTIRELY by the `Genotype::hasAD` precondition, so an allele
 *     nobody reports leaves an empty list whose maximum is taken as ZERO and is therefore filtered;
 *   - AN EXISTING ALLELE FILTER IS UNIONED THROUGH A LinkedHashSet, so order is insertion order and
 *     a repeat is dropped, and the placeholder SITE is REPLACED rather than added to;
 *   - AND THE HEADER GAINS THE FILTER LINE whether or not anything is filtered.
 *
 * Output:
 *
 *     cutoff\t<coverage>\t<copies>\t<the integer cutoff>
 *     input\t<label>=<the whole input vcf, escaped>
 *     filtered\t<label>=<the whole output vcf, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: NuMTFilterDump
 */

import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.NuMTFilterTool;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class NuMTFilterDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
            + "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            + "##INFO=<ID=AS_FilterStatus,Number=A,Type=String,Description=\"Filter status for each allele\">\n"
            + "##contig=<ID=chrM,length=16569>\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tone\ttwo\n";

    /** Records that already carry AS_FilterStatus, which is the only shape the tool can merge into. */
    static final String WITH_STATUS =
            // A deep alternate, above any cutoff these runs use.
            "chrM\t100\t.\tA\tC\t50\t.\tAS_FilterStatus=SITE\tGT:AD\t0/1:100,100\t0/0:200,0\n"
            // A shallow one, below them.
            + "chrM\t200\t.\tA\tC\t50\t.\tAS_FilterStatus=SITE\tGT:AD\t0/1:200,5\t0/0:200,0\n"
            // Two alternates, one deep and one shallow, so the attribute is written and the site
            // filter is not.
            + "chrM\t300\t.\tA\tC,G\t50\t.\tAS_FilterStatus=SITE|SITE\tGT:AD\t0/1:100,100,4\t0/0:200,0,0\n"
            // The maximum across samples rather than the sum: two samples of fifty each, so the
            // maximum is under the cutoff of seventy-nine while the sum is over it.
            + "chrM\t400\t.\tA\tC\t50\t.\tAS_FilterStatus=SITE\tGT:AD\t0/1:100,50\t0/1:100,50\n"
            // No AD anywhere, so the allele's list is empty and its maximum is taken as zero.
            + "chrM\t500\t.\tA\tC\t50\t.\tAS_FilterStatus=SITE\tGT\t0/1\t0/0\n"
            // An allele filter already set, which is unioned rather than replaced.
            + "chrM\t600\t.\tA\tC\t50\t.\tAS_FilterStatus=weak_evidence\tGT:AD\t0/1:200,5\t0/0:200,0\n"
            // The same filter already set, which the set drops.
            + "chrM\t700\t.\tA\tC\t50\t.\tAS_FilterStatus=possible_numt\tGT:AD\t0/1:200,5\t0/0:200,0\n";

    /** The same records without the attribute, which is what an ordinary VCF looks like. */
    static final String WITHOUT_STATUS =
            "chrM\t100\t.\tA\tC\t50\t.\t.\tGT:AD\t0/1:100,100\t0/0:200,0\n"
            + "chrM\t200\t.\tA\tC\t50\t.\t.\tGT:AD\t0/1:200,5\t0/0:200,0\n";

    /** A record whose ALT is absent, so the list of alternates is empty. */
    static final String NO_ALT =
            "chrM\t100\t.\tA\t.\t50\t.\tAS_FilterStatus=SITE\tGT:AD\t0/0:200\t0/0:200\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("numt-filter-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# NuMTFilterDump: mitochondrial calls filtered by a Poisson depth cutoff");

        // The cutoff on its own, which is the only arithmetic in the tool.
        for (final double coverage : new double[] {1, 5, 10, 30, 100, 1000, 0.5}) {
            for (final double copies : new double[] {1, 4, 0.5}) {
                System.out.printf("cutoff\t%s\t%s\t%d%n", coverage, copies,
                        cutoff(copies, coverage));
            }
        }

        // The defaults, where the coverage is zero and nothing at all happens.
        run(dir, "defaults", HEADER + WITH_STATUS);
        // A coverage that puts the cutoff between the shallow and the deep alternates.
        run(dir, "coverage-30", HEADER + WITH_STATUS, "--autosomal-coverage", "30");
        // Copies of zero, which keeps the cutoff at zero however deep the coverage.
        run(dir, "copies-0", HEADER + WITH_STATUS,
                "--autosomal-coverage", "30", "--max-numt-autosomal-copies", "0");
        // A cutoff so high that every alternate falls under it.
        run(dir, "coverage-1000", HEADER + WITH_STATUS, "--autosomal-coverage", "1000");
        // The same records without AS_FilterStatus: harmless while nothing is filtered.
        run(dir, "no-status-defaults", HEADER + WITHOUT_STATUS);
        // And the same, once something is.
        run(dir, "no-status-filtered", HEADER + WITHOUT_STATUS, "--autosomal-coverage", "30");
        // A record with no alternate allele, whose list of decisions is empty.
        run(dir, "no-alt", HEADER + NO_ALT, "--autosomal-coverage", "30");
    }

    /**
     * `getMaxAltDepthCutoff`, which is package-private and `@VisibleForTesting`, so it is reached
     * reflectively rather than by moving this dump into the tool's package.
     */
    static int cutoff(final double copies, final double coverage) throws Exception {
        final java.lang.reflect.Method method = NuMTFilterTool.class
                .getDeclaredMethod("getMaxAltDepthCutoff", double.class, double.class);
        method.setAccessible(true);
        return (Integer) method.invoke(null, copies, coverage);
    }

    static void run(final Path dir, final String label, final String vcf, final String... extra)
            throws Exception {
        final Path in = dir.resolve(label + ".vcf");
        Files.writeString(in, vcf, StandardCharsets.UTF_8);
        System.out.printf("input\t%s=%s%n", label, ReferenceQueryDump.escape(vcf));

        final Path out = dir.resolve(label + "-filtered.vcf");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-V", in.toString(), "-O", out.toString()));
        argv.addAll(Arrays.asList(extra));
        try {
            new NuMTFilterTool().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            return;
        }
        if (Files.exists(out)) {
            System.out.printf("filtered\t%s=%s%n", label,
                    ReferenceQueryDump.escape(masked(Files.readString(out), dir)));
        }
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>")
                .replaceAll("##GATKCommandLine=<[^\n]*>\n", "");
    }
}

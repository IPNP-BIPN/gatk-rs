/*
 * What a VariantWalker hands to apply(), taken from the reference.
 *
 * The traversal every variant-based tool inherits, run through the real command line rather than
 * reconstructed, so the defaults are measured. Two of its behaviours are not what "restrict the
 * traversal to -L" suggests, and both are visible only when the intervals are not the sorted,
 * non-overlapping list the class asks for:
 *
 *   - the de-duplication remembers ONE interval, not all of them. FeatureIntervalIterator's
 *     featureIsNovel tests the feature against previousInterval only, so a variant covered by
 *     intervals 1 and 3 but not by 2 is handed to apply twice, while the same variant covered by
 *     two adjacent intervals arrives once. Nothing enforces the precondition: breaking it produces
 *     duplicates rather than an error;
 *   - an empty interval list is not an empty traversal. setIntervalsForTraversal maps both null
 *     and an empty list to null, which means no restriction at all.
 *
 * The record fields are dumped alongside, because the walker is also where the decoded record
 * becomes what a tool sees: the interval it builds is SimpleInterval(variant), whose end is the
 * variant's END attribute when it has one and start + ref.length() - 1 otherwise.
 *
 * Output:
 *
 *     apply\t<label>\t<n>\t<contig>:<start>-<stop>|<id>|<alleles>|<filters>
 *     summary\t<label>\t<ok|E:class>
 *     count\t<label>\t<number of apply calls>
 *
 * Usage: VariantWalkerDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;

public class VariantWalkerDump {

    /** Every apply() call of the current traversal, filled by the probe walker. */
    static final List<String> APPLIED = new ArrayList<>();

    /**
     * The fixture, sorted and single-sample.
     *
     * Two records are deliberately adjacent (200 and 201) so an interval boundary can fall between
     * them, and one carries an END so its stop is not derived from the reference allele.
     */
    static final String VCF = String.join("\n",
            "##fileformat=VCFv4.2",
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">",
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "##FILTER=<ID=LowQual,Description=\"Low quality\">",
            "##contig=<ID=chr1,length=1000>",
            "##contig=<ID=chr2,length=1000>",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1",
            "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT\t0/1",
            "chr1\t150\trs1\tACGT\tA\t50\tLowQual\tDP=20\tGT\t1/1",
            "chr1\t200\t.\tG\tC\t50\t.\tDP=30\tGT\t0/1",
            "chr1\t201\t.\tT\tG\t50\tPASS\tDP=40\tGT\t0/0",
            "chr1\t300\t.\tA\t<DEL>\t50\tPASS\tEND=400\tGT\t0/1",
            "chr2\t100\t.\tC\tG\t50\tPASS\tDP=50\tGT\t0/1",
            "");

    @CommandLineProgramProperties(
            summary = "Records what a VariantWalker hands to apply()",
            oneLineSummary = "VariantWalker traversal probe",
            programGroup = ReadDataManipulationProgramGroup.class)
    public static class ProbeWalker extends VariantWalker {
        @Override
        public void apply(final VariantContext variant, final ReadsContext reads,
                          final ReferenceContext reference, final FeatureContext features) {
            final StringJoiner alleles = new StringJoiner(",");
            for (final Allele allele : variant.getAlleles()) {
                alleles.add(allele.getDisplayString());
            }
            final String filters;
            if (!variant.filtersWereApplied()) {
                filters = "unfiltered";
            } else if (variant.getFilters().isEmpty()) {
                filters = "PASS";
            } else {
                filters = String.join(",", new java.util.TreeSet<>(variant.getFilters()));
            }
            APPLIED.add(String.format("%s:%d-%d|%s|%s|%s",
                    variant.getContig(), variant.getStart(), variant.getEnd(), variant.getID(),
                    alleles, filters));
        }
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("variantwalker-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path vcf = dir.resolve("variants.vcf");
        Files.write(vcf, VCF.getBytes());
        // A traversal restricted by -L needs an index, which is what makes the per-interval query
        // path run at all.
        new org.broadinstitute.hellbender.tools.IndexFeatureFile()
                .instanceMain(new String[] {"-I", vcf.toString()});

        System.out.println("# VariantWalkerDump: what a VariantWalker hands to apply()");

        traverse("all", vcf);
        traverse("one-interval", vcf, "-L", "chr1:100-160");
        // A boundary between two adjacent records.
        traverse("boundary", vcf, "-L", "chr1:150-200");
        // Two intervals, sorted and non-overlapping, which is the documented precondition.
        traverse("two-sorted", vcf, "-L", "chr1:100-160", "-L", "chr1:200-210");
        // Two intervals that overlap, where the one-interval memory suppresses the repeat.
        traverse("two-overlapping", vcf, "-L", "chr1:100-200", "-L", "chr1:150-250");
        // Three intervals where the first and third cover the same record and the second does not,
        // which is the case the one-interval memory cannot suppress.
        traverse("gap-then-repeat", vcf,
                "-L", "chr1:100-160", "-L", "chr1:500-600", "-L", "chr1:100-160");
        // Intervals out of order, which the class asks the caller not to do.
        traverse("unsorted", vcf, "-L", "chr1:200-210", "-L", "chr1:100-160");
        // The END record, queried by an interval that touches only its tail.
        traverse("end-tail", vcf, "-L", "chr1:390-410");
        // A second contig, and an interval on a contig with no records.
        traverse("chr2", vcf, "-L", "chr2");
        traverse("empty-contig-interval", vcf, "-L", "chr2:500-600");
        // Excluding an interval rather than including one.
        traverse("exclude", vcf, "-XL", "chr1");
        // A record-level filter argument, to confirm the traversal itself applies none by default.
        traverse("select-passing", vcf, "-L", "chr1");
    }

    static void traverse(final String label, final Path vcf, final String... extra) {
        APPLIED.clear();
        final List<String> argv = new ArrayList<>(Arrays.asList("-V", vcf.toString()));
        argv.addAll(Arrays.asList(extra));

        String summary;
        try {
            new ProbeWalker().instanceMain(argv.toArray(new String[0]));
            summary = "ok";
        } catch (final Exception | AssertionError e) {
            summary = "E:" + e.getClass().getName();
        }
        for (int i = 0; i < APPLIED.size(); i++) {
            System.out.printf("apply\t%s\t%d\t%s%n", label, i, APPLIED.get(i));
        }
        System.out.printf("summary\t%s\t%s%n", label, summary);
        System.out.printf("count\t%s\t%d%n", label, APPLIED.size());
    }
}

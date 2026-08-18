/*
 * `StrandArtifactFilter`'s E step, taken from the reference.
 *
 * The last of the ten filters the `filter-mutect-calls` golden needs, split in two: this is
 * `calculateArtifactProbabilities` and the `strandArtifactProbability` under it, at the initial
 * parameters. The M step needs a Brent optimizer and is measured separately.
 *
 * Six behaviours this is built to catch.
 *
 *   - THE STRAND COUNTS COME OUT OF A STRING. `AS_SB_TABLE` is split on `|` and then on `,`, with
 *     brackets stripped and each field trimmed before `Integer.parseInt`. A non-integer field is a
 *     `NumberFormatException` out of a filter, an empty annotation is an empty list, and A TABLE
 *     WITH ONE ENTRY IS DROPPED: `sbs.size() <= 1` guards the whole filter, not one allele;
 *   - TWO BRANCHES ANSWER A HARD ZERO RATHER THAN A PROBABILITY. An alternate with no reads on
 *     either strand, and an indel longer than `LONGEST_STRAND_ARTIFACT_INDEL_SIZE = 4`, both build
 *     an `EStep` with both responsibilities zero. The filter still reports for that allele;
 *   - THE SEQUENCING-ERROR PRIOR TAKES THREE STEPS IN THE INDEL SIZE: 1000 at zero, 5000 below
 *     three, 50000 from three to four, and the size is `Math.abs(ref.length() - alt.length())`, so
 *     an insertion and a deletion of the same length take the same branch;
 *   - THE TOTALS ARE OVER EVERY ALLELE INCLUDING THE REFERENCE while the responsibilities are per
 *     alternate;
 *   - THE NORMALISATION IS IN LOG10 AND THE LIKELIHOODS ARE NATURAL, converted by multiplying by
 *     `MathUtils.LOG10_E` rather than by changing the base of the logarithm;
 *   - AND `binomialCoefficientLog` APPEARS THREE TIMES IN ONE EXPRESSION in the no-artifact branch,
 *     which is where a deep site reaches the sum-of-logs route.
 *
 * Output:
 *
 *     default\t<name>\t<value>
 *     name\tstrand-artifact\t<filterName>,<errorType>,<annotation>,<required annotations>
 *     estep\t<label>\t<forward>,<reverse>,<fwdCount>,<revCount>,<fwdAlt>,<revAlt>
 *     prob\t<label>\t<one error probability per alternate allele>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: StrandArtifactEStepDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.StrandArtifactFilter;

import java.io.File;
import java.lang.reflect.Method;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

public class StrandArtifactEStepDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele SNV = Allele.create("C", false);
    static final Allele SNV2 = Allele.create("G", false);
    static final Allele SYMBOLIC = Allele.create("<NON_REF>", false);

    public static void main(final String[] args) throws Exception {
        System.out.println("# StrandArtifactEStepDump: the strand-artifact filter's E step");

        // The initial parameters, which are private fields.
        System.out.println("default\tinitialAlphaStrand\t1.0");
        System.out.println("default\tinitialBetaStrand\t20.0");
        System.out.println("default\tinitialStrandArtifactPrior\t0.001");

        final StrandArtifactFilter filter = new StrandArtifactFilter();
        System.out.printf("name\tstrand-artifact\t%s,%s,%s,%s%n", filter.filterName(), filter.errorType(),
                filter.phredScaledPosteriorAnnotationName().orElse("none"), "AS_SB_TABLE");

        // A biallelic SNV whose alternate is balanced across the strands.
        run("balanced", biallelic("50,50|10,10"));
        // Forward only, and reverse only: the two artifact hypotheses in turn.
        run("forward-only", biallelic("50,50|20,0"));
        run("reverse-only", biallelic("50,50|0,20"));
        // One read on one strand, which is where the prior dominates.
        run("one-forward-read", biallelic("50,50|1,0"));
        // No alternate reads at all: the hard zero.
        run("no-alt-reads", biallelic("50,50|0,0"));
        // Depth, which is what moves the beta binomial.
        run("shallow", biallelic("5,5|3,0"));
        run("deep", biallelic("2000,2000|400,0"));

        // The indel-size branches of the sequencing-error prior.
        run("snv", indel("A", "C", "50,50|20,0"));
        run("one-base-deletion", indel("AA", "A", "50,50|20,0"));
        run("two-base-insertion", indel("A", "AAA", "50,50|20,0"));
        run("three-base-insertion", indel("A", "AAAA", "50,50|20,0"));
        run("four-base-deletion", indel("AAAAA", "A", "50,50|20,0"));
        // One past the longest strand-artifact indel: the other hard zero.
        run("five-base-deletion", indel("AAAAAA", "A", "50,50|20,0"));

        // Two alternates, each with its own strand counts, and totals over all three entries.
        run("two-alternates", new VariantContextBuilder("dump", "chr1", 100, 100,
                List.of(REF, SNV, SNV2))
                .attribute("AS_SB_TABLE", "50,50|20,0|5,5").make());

        // A symbolic alternate, whose data the filter removes.
        run("symbolic-alternate", new VariantContextBuilder("dump", "chr1", 100, 100,
                List.of(REF, SNV, SYMBOLIC))
                .attribute("AS_SB_TABLE", "50,50|20,0|1,1").make());

        // The table's own shapes.
        run("one-entry-table", biallelic("50,50"));
        run("empty-table", biallelic(""));
        run("no-table", new VariantContextBuilder("dump", "chr1", 100, 100, List.of(REF, SNV)).make());
        run("bracketed-table", biallelic("[50,50|20,0]"));
        run("spaced-table", biallelic("50, 50 | 20, 0"));
        run("non-integer-table", biallelic("50,50|twenty,0"));
        run("empty-field", biallelic("50,50|,0"));

        // `strandArtifactProbability` directly, which is package-private, over the prior and the
        // indel size the record cannot vary independently.
        direct("direct-default-prior", 0.001, 50, 50, 20, 0, 0);
        direct("direct-high-prior", 0.5, 50, 50, 20, 0, 0);
        direct("direct-zero-prior", 0.0, 50, 50, 20, 0, 0);
        direct("direct-prior-of-one", 1.0, 50, 50, 20, 0, 0);
        direct("direct-long-indel", 0.001, 50, 50, 20, 0, 3);
        direct("direct-no-reads", 0.001, 0, 0, 0, 0, 0);
    }

    static VariantContext biallelic(final String table) {
        final VariantContextBuilder builder =
                new VariantContextBuilder("dump", "chr1", 100, 100, List.of(REF, SNV));
        return builder.attribute("AS_SB_TABLE", table).make();
    }

    static VariantContext indel(final String ref, final String alt, final String table) {
        final Allele reference = Allele.create(ref, true);
        final Allele alternate = Allele.create(alt, false);
        return new VariantContextBuilder("dump", "chr1", 100, 100 + ref.length() - 1,
                List.of(reference, alternate))
                .attribute("AS_SB_TABLE", table).make();
    }

    static Mutect2FilteringEngine engine() {
        final Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        lines.add(new VCFHeaderLine("normal_sample", "N1"));
        final VCFHeader header = new VCFHeader(lines, List.of("T1", "N1"));
        return new Mutect2FilteringEngine(new M2FiltersArgumentCollection(), header,
                new File("no-such-stats-file.tsv"));
    }

    /** Both the E steps and the probabilities they become, from one fresh filter. */
    static void run(final String label, final VariantContext vc) {
        final StrandArtifactFilter filter = new StrandArtifactFilter();
        try {
            for (final StrandArtifactFilter.EStep step : filter.calculateArtifactProbabilities(vc, engine())) {
                System.out.printf("estep\t%s\t%s%n", label, render(step));
            }
            System.out.printf("prob\t%s\t%s%n", label, filter.errorProbabilities(vc, engine(), null));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    /** `strandArtifactProbability`, which is package-private and `@VisibleForTesting`. */
    static void direct(final String label, final double prior, final int forwardCount,
                       final int reverseCount, final int forwardAltCount, final int reverseAltCount,
                       final int indelSize) {
        try {
            final Method method = StrandArtifactFilter.class.getDeclaredMethod(
                    "strandArtifactProbability", double.class, int.class, int.class, int.class,
                    int.class, int.class);
            method.setAccessible(true);
            final StrandArtifactFilter.EStep step = (StrandArtifactFilter.EStep) method.invoke(
                    new StrandArtifactFilter(), prior, forwardCount, reverseCount, forwardAltCount,
                    reverseAltCount, indelSize);
            System.out.printf("estep\t%s\t%s%n", label, render(step));
        } catch (final Exception e) {
            final Throwable cause = e.getCause() == null ? e : e.getCause();
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(cause.getMessage())));
        }
    }

    static String render(final StrandArtifactFilter.EStep step) {
        return Double.toString(step.getForwardArtifactResponsibility()) + ","
                + Double.toString(step.getReverseArtifactResponsibility()) + ","
                + step.getForwardCount() + "," + step.getReverseCount() + ","
                + step.getForwardAltCount() + "," + step.getReverseAltCount();
    }
}

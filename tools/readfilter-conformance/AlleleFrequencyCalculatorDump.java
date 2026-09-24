/*
 * AlleleFrequencyCalculator and AFCalculationResult, taken from the reference.
 *
 * What GenotypeGVCFs asks before it calls anything: how likely the site is to hold no variant at
 * all, which is its QUAL, how likely each alternate is to be absent, which decides the alleles it
 * keeps, and the integer allele counts that become MLEAC. The answer comes out of an iterative
 * Dirichlet fit, so every number is printed as its raw bits.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE FIRST ITERATION USES A FLAT PRIOR, one over the allele count for every allele, and the
 *     prior pseudocounts only enter from the second, so a single strong het is not drowned by the
 *     reference pseudocount;
 *   - THE PSEUDOCOUNTS DEPEND ON THE ALLELE'S LENGTH, not on its kind: an alternate as long as the
 *     reference is a SNP and gets the SNP pseudocount, any other length the indel one;
 *   - A BIALLELIC SITE SHORT-CIRCUITS the per-allele sum, and its one allele-absent probability is
 *     the no-variant probability itself;
 *   - A SPANNING DELETION COUNTS AS NO VARIANT: genotypes made only of the reference and `*` are
 *     summed into the no-variant probability, capped at zero;
 *   - A HOM-REF WITH NO LIKELIHOODS BUT A GQ, and only if diploid, is given likelihoods invented
 *     from the GQ (0, GQ, 20 GQ, every het at GQ), and any other genotype without likelihoods is
 *     skipped;
 *   - THE MLE ALLELE COUNTS ARE THE ROUNDED EFFECTIVE COUNTS of the last iteration, not a maximum;
 *   - AND `passesThreshold` ADDS AN EPSILON of 1e-10 before comparing with the phred threshold.
 *
 * Output:
 *
 *     calc\t<label>\t<log10 P(no variant)>\t<log10 P(variant)>\t<MLE alt counts>\t<per alternate: log10 P(absent)>\t<per alternate: passes 10, 30, 50>
 *     error\t<label>\t<exception class>:<message>
 *
 * Doubles are printed as the sixteen hexadecimal digits of their raw bits.
 *
 * Usage: AlleleFrequencyCalculatorDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AFCalculationResult;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AlleleFrequencyCalculator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class AlleleFrequencyCalculatorDump {

    static final Allele A = Allele.create("A", true);
    static final Allele C = Allele.create("C", false);
    static final Allele G = Allele.create("G", false);
    static final Allele T = Allele.create("T", false);
    static final Allele AT = Allele.create("AT", false);
    static final Allele DEL_REF = Allele.create("ACG", true);
    static final Allele DEL_ALT = Allele.create("A", false);
    static final Allele SPAN = Allele.SPAN_DEL;

    public static void main(final String[] args) {
        System.out.println("# AlleleFrequencyCalculatorDump: the Dirichlet fit behind QUAL, from the reference");

        final AlleleFrequencyCalculator standard =
                AlleleFrequencyCalculator.makeCalculator(new GenotypeCalculationArgumentCollection());

        // One diploid sample, biallelic: a het, a hom-var, the reference best, and a marginal call.
        run(standard, "het", List.of(A, C), pl("s1", A, C, 300, 0, 500));
        run(standard, "hom-var", List.of(A, C), pl("s1", C, C, 900, 60, 0));
        run(standard, "reference-best", List.of(A, C), pl("s1", A, A, 0, 30, 400));
        run(standard, "marginal", List.of(A, C), pl("s1", A, C, 8, 0, 300));
        run(standard, "flat", List.of(A, C), pl("s1", A, C, 0, 0, 0));

        // A cohort: a het, a hom-ref, a hom-var and a second het.
        run(standard, "cohort", List.of(A, C),
                pl("s1", A, C, 250, 0, 480),
                pl("s2", A, A, 0, 45, 700),
                pl("s3", C, C, 820, 55, 0),
                pl("s4", A, C, 120, 0, 350));
        // Many hom-refs, which is what a joint call of a rare het looks like.
        final List<Genotype> rare = new ArrayList<>();
        rare.add(pl("s0", A, C, 200, 0, 400));
        for (int i = 1; i < 20; i++) {
            rare.add(pl("s" + i, A, A, 0, 30 + i, 450));
        }
        run(standard, "rare-het", List.of(A, C), rare.toArray(new Genotype[0]));
        run(standard, "all-hom-ref", List.of(A, C),
                pl("s1", A, A, 0, 60, 900), pl("s2", A, A, 0, 42, 630), pl("s3", A, A, 0, 21, 315));

        // Three and five alleles, a SNP beside an indel, which take different pseudocounts.
        run(standard, "snp-and-snp", List.of(A, C, G), pl("s1", C, G, 900, 400, 500, 450, 0, 800));
        run(standard, "snp-and-indel", List.of(A, C, AT),
                pl("s1", A, C, 300, 0, 500, 280, 450, 900),
                pl("s2", A, AT, 260, 350, 700, 0, 380, 640));
        run(standard, "indel-only", List.of(DEL_REF, DEL_ALT), pl("s1", DEL_REF, DEL_ALT, 400, 0, 600));
        run(standard, "five-alleles", List.of(A, C, G, T, AT),
                pl("s1", C, G, 900, 400, 500, 450, 0, 800, 600, 700, 650, 900, 620, 710, 660, 910, 990),
                pl("s2", A, T, 300, 350, 600, 380, 610, 640, 0, 400, 420, 500, 330, 440, 460, 520, 800));

        // Other ploidies: haploid, triploid, tetraploid, and a haploid beside a diploid.
        run(standard, "haploid", List.of(A, C), pl("s1", new Allele[] {C}, 200, 0));
        run(standard, "triploid", List.of(A, C), pl("s1", new Allele[] {A, A, C}, 150, 0, 90, 400));
        run(standard, "tetraploid", List.of(A, C), pl("s1", new Allele[] {A, A, C, C}, 400, 60, 0, 70, 500));
        run(standard, "mixed-ploidy", List.of(A, C),
                pl("s1", new Allele[] {C}, 200, 0),
                pl("s2", A, C, 300, 0, 500));

        // A spanning deletion, alone and beside a SNP.
        run(standard, "span-del-only", List.of(A, SPAN), pl("s1", A, SPAN, 300, 0, 500));
        run(standard, "span-del-and-snp", List.of(A, SPAN, C),
                pl("s1", A, SPAN, 300, 0, 500, 320, 520, 900),
                pl("s2", A, C, 280, 400, 700, 0, 350, 600));

        // A hom-ref with only a GQ, which is given likelihoods, beside one that has them.
        run(standard, "gq-hom-ref", List.of(A, C),
                pl("s1", A, C, 300, 0, 500),
                gq("s2", 2, 40));
        run(standard, "gq-hom-ref-triallelic", List.of(A, C, G),
                pl("s1", C, G, 900, 400, 500, 450, 0, 800),
                gq("s2", 2, 25));
        // A haploid hom-ref with only a GQ is not usable and is skipped.
        run(standard, "gq-haploid-skipped", List.of(A, C),
                pl("s1", A, C, 300, 0, 500),
                gq("s2", 1, 40));
        // A no-call with likelihoods is used; one with neither is skipped.
        run(standard, "no-call-with-pl", List.of(A, C),
                pl("s1", A, C, 300, 0, 500),
                new GenotypeBuilder("s2", List.of(Allele.NO_CALL, Allele.NO_CALL)).PL(new int[] {0, 10, 100}).make());
        run(standard, "no-call-without-pl", List.of(A, C),
                pl("s1", A, C, 300, 0, 500),
                new GenotypeBuilder("s2", List.of(Allele.NO_CALL, Allele.NO_CALL)).make());

        // A calculator with other priors: rarer SNPs and indels, a wider spread.
        final GenotypeCalculationArgumentCollection other = new GenotypeCalculationArgumentCollection();
        other.snpHeterozygosity = 0.01;
        other.indelHeterozygosity = 0.0005;
        other.heterozygosityStandardDeviation = 0.05;
        final AlleleFrequencyCalculator wide = AlleleFrequencyCalculator.makeCalculator(other);
        run(wide, "wide-prior-het", List.of(A, C), pl("s1", A, C, 300, 0, 500));
        run(wide, "wide-prior-snp-and-indel", List.of(A, C, AT),
                pl("s1", A, C, 300, 0, 500, 280, 450, 900),
                pl("s2", A, AT, 260, 350, 700, 0, 380, 640));

        // The refusals.
        run(standard, "only-reference", List.of(A), new GenotypeBuilder("s1", List.of(A, A)).PL(new int[] {0}).make());
        run(standard, "no-likelihoods", List.of(A, C), gq("s1", 2, 40));
        run(standard, "inconsistent-pl", List.of(A, C, G), pl("s1", A, C, 300, 0, 500));
    }

    static Genotype pl(final String sample, final Allele first, final Allele second, final int... pls) {
        return pl(sample, new Allele[] {first, second}, pls);
    }

    static Genotype pl(final String sample, final Allele[] alleles, final int... pls) {
        return new GenotypeBuilder(sample, Arrays.asList(alleles)).PL(pls).make();
    }

    /** A hom-ref of the given ploidy carrying a GQ and no likelihoods. */
    static Genotype gq(final String sample, final int ploidy, final int gq) {
        final List<Allele> alleles = new ArrayList<>();
        for (int i = 0; i < ploidy; i++) {
            alleles.add(A);
        }
        return new GenotypeBuilder(sample, alleles).GQ(gq).make();
    }

    static void run(final AlleleFrequencyCalculator calculator, final String label, final List<Allele> alleles,
                    final Genotype... genotypes) {
        final Allele reference = alleles.get(0);
        final VariantContext vc = new VariantContextBuilder("dump", "chr1", 100, 100 + reference.length() - 1, alleles)
                .genotypes(genotypes).make();
        try {
            final AFCalculationResult result = calculator.calculate(vc);
            final List<Allele> alternates = alleles.subList(1, alleles.size());
            final String counts = Arrays.stream(result.getAlleleCountsOfMLE())
                    .mapToObj(Integer::toString).collect(Collectors.joining(","));
            final String absent = alternates.stream()
                    .map(a -> bits(result.getLog10PosteriorOfAlleleAbsent(a))).collect(Collectors.joining(","));
            final String passes = alternates.stream()
                    .map(a -> (result.passesThreshold(a, 10) ? "1" : "0")
                            + (result.passesThreshold(a, 30) ? "1" : "0")
                            + (result.passesThreshold(a, 50) ? "1" : "0"))
                    .collect(Collectors.joining(","));
            System.out.printf("calc\t%s\t%s\t%s\t%s\t%s\t%s%n", label,
                    bits(result.log10ProbOnlyRefAlleleExists()), bits(result.log10ProbVariantPresent()),
                    counts, absent, passes);
        } catch (final RuntimeException e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(), e.getMessage());
        }
    }

    static String bits(final double value) {
        return String.format("%016x", Double.doubleToRawLongBits(value));
    }
}

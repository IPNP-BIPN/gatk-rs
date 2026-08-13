/*
 * GenotypeIndexCalculator and GenotypeAlleleCounts, taken from the reference.
 *
 * The combinatorics under every PL array: which genotype each index of a likelihood vector means,
 * and where a genotype lands when the allele list changes. It is what AlleleSubsettingUtils needs
 * before it can subset a single PL, and it is the piece between LeftAlignAndTrimVariants and a
 * genotyped vcf.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE INDEX IS A SUM OF BINOMIALS, not a position in a nested loop, and the canonical order is
 *     the one `GenotypeAlleleCounts.first(ploidy).next()` walks: for two alleles diploid it is
 *     0/0, 0/1, 1/1, and for three it continues 0/2, 1/2, 2/2;
 *   - THE NUMBER OF GENOTYPES IS `C(alleles + ploidy - 1, ploidy)`, so a triploid site with four
 *     alleles has twenty and not sixty-four;
 *   - A GENOTYPE IS A MULTISET, so `alleleCountsToIndex` takes pairs of allele and count and the
 *     ORDER OF THE PAIRS DOES NOT MATTER, while a repeated allele does;
 *   - AN ODD-LENGTH COUNT ARRAY IS REFUSED, which is the only argument check;
 *   - PLOIDY ZERO HAS EXACTLY ONE GENOTYPE, the empty one, and its index is 0;
 *   - THE ITERATION IS DESTRUCTIVE IN THE REFERENCE: `GenotypeAlleleCounts.next()` returns a new
 *     object but `increase()` mutates, and the iterable hands out the SAME object each time, so a
 *     caller that keeps one keeps the last;
 *   - `distinctAlleleCount` COUNTS ALLELES PRESENT, not ploidy, so 1/1 has one and 0/1 has two;
 *   - AND `alleleCountFor` IS ZERO FOR AN ALLELE THAT IS NOT THERE, which is what makes the
 *     subsetting loop work at all.
 *
 * Output:
 *
 *     count\t<ploidy>\t<alleles>\t<number of genotypes>
 *     order\t<ploidy>\t<alleles>\t<index>\t<the genotype as allele:count pairs>
 *     index\t<label>\t<the allele counts given>\t<the index>
 *     subset\t<label>\t<ploidy>\t<old alleles>\t<kept alleles>\t<the pl indices>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: GenotypeIndexDump
 */

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.walkers.genotyper.AlleleSubsettingUtils;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAlleleCounts;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeIndexCalculator;

import java.util.ArrayList;
import java.util.List;

public class GenotypeIndexDump {

    public static void main(final String[] args) {
        System.out.println("# GenotypeIndexDump: the genotype index combinatorics, from the reference");

        // How many genotypes, over the shapes a real vcf reaches and a few beyond.
        for (final int ploidy : new int[] {0, 1, 2, 3, 4}) {
            for (final int alleles : new int[] {1, 2, 3, 4, 6}) {
                System.out.printf("count\t%d\t%d\t%d%n", ploidy, alleles,
                        GenotypeIndexCalculator.genotypeCount(ploidy, alleles));
            }
        }

        // The canonical order, which is what a PL array is indexed by.
        order(2, 2);
        order(2, 3);
        order(2, 4);
        order(3, 3);
        order(1, 4);
        order(0, 3);

        // The index of a genotype given as allele/count pairs.
        index("diploid-ref", 0, 2);
        index("diploid-het", 0, 1, 1, 1);
        index("diploid-hom-alt", 1, 2);
        index("diploid-het-non-ref", 1, 1, 2, 1);
        // The same genotype with its pairs the other way round.
        index("diploid-het-non-ref-reversed", 2, 1, 1, 1);
        index("triploid-mixed", 0, 1, 1, 1, 2, 1);
        index("ploidy-zero");
        // A count of zero for an allele, which is a pair the caller may pass.
        index("zero-count", 0, 2, 1, 0);
        // And the one argument check.
        indexRefused("odd-length", 0, 2, 1);

        // The subsetting itself: which old index each new index takes its likelihood from.
        subset("keep-first-alt", 2, 3, new int[] {0, 1});
        subset("keep-second-alt", 2, 3, new int[] {0, 2});
        subset("keep-ref-only", 2, 3, new int[] {0});
        subset("keep-both-swapped", 2, 3, new int[] {0, 2, 1});
        subset("four-to-two", 2, 4, new int[] {0, 3});
        subset("triploid-keep-one", 3, 3, new int[] {0, 1});
        subset("keep-everything", 2, 3, new int[] {0, 1, 2});
    }

    /** Every genotype of one shape, in the order the reference walks them. */
    static void order(final int ploidy, final int alleles) {
        int index = 0;
        for (final GenotypeAlleleCounts counts : GenotypeAlleleCounts.iterable(ploidy, alleles)) {
            final List<String> pairs = new ArrayList<>();
            for (int position = 0; position < counts.distinctAlleleCount(); position++) {
                pairs.add(counts.alleleIndexAt(position) + ":" + counts.alleleCountAt(position));
            }
            System.out.printf("order\t%d\t%d\t%d\t%s%n", ploidy, alleles, index,
                    String.join(",", pairs));
            index++;
        }
    }

    static void index(final String label, final int... alleleCounts) {
        final List<String> given = new ArrayList<>();
        for (final int value : alleleCounts) {
            given.add(String.valueOf(value));
        }
        System.out.printf("index\t%s\t%s\t%d%n", label, String.join(",", given),
                GenotypeIndexCalculator.alleleCountsToIndex(alleleCounts));
    }

    static void indexRefused(final String label, final int... alleleCounts) {
        final List<String> given = new ArrayList<>();
        for (final int value : alleleCounts) {
            given.add(String.valueOf(value));
        }
        try {
            GenotypeIndexCalculator.alleleCountsToIndex(alleleCounts);
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("error\t%s\tnone%n", label);
    }

    /** `subsettedPLIndices`, over allele lists built from single bases. */
    static void subset(final String label, final int ploidy, final int alleleCount,
                       final int[] keep) {
        final List<Allele> original = alleles(alleleCount);
        final List<Allele> kept = new ArrayList<>();
        for (final int index : keep) {
            kept.add(original.get(index));
        }
        final int[] indices = AlleleSubsettingUtils.subsettedPLIndices(ploidy, original, kept);
        final List<String> out = new ArrayList<>();
        for (final int value : indices) {
            out.add(String.valueOf(value));
        }
        final List<String> keptNames = new ArrayList<>();
        for (final int index : keep) {
            keptNames.add(String.valueOf(index));
        }
        System.out.printf("subset\t%s\t%d\t%d\t%s\t%s%n", label, ploidy, alleleCount,
                String.join(",", keptNames), String.join(",", out));
    }

    /** `alleleCount` alleles, the first being the reference. */
    static List<Allele> alleles(final int alleleCount) {
        final String bases = "ACGTNM";
        final List<Allele> out = new ArrayList<>();
        for (int index = 0; index < alleleCount; index++) {
            out.add(Allele.create(bases.substring(index, index + 1), index == 0));
        }
        return out;
    }
}

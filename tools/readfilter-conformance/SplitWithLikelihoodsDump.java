/*
 * splitVariantContextToBiallelics over genotypes that carry likelihoods, taken from the reference.
 *
 * The split-biallelics suite measured the splitter over records whose genotypes had no PLs, and the
 * subset-alleles suite measured what the subsetting does to one genotype in isolation. This is the
 * junction: the two together, which is what a real genotyped vcf goes through.
 *
 * Seven behaviours this is built to catch.
 *
 *   - EVERY OUTPUT RECORD CARRIES ITS OWN SUBSET of the PLs, so one record with two alternates
 *     leaves two records whose likelihood vectors are different permutations of the same six
 *     numbers, each rescaled on its own;
 *   - AND THE CALL DOES NOT FOLLOW THE PLs. Under BEST_MATCH_TO_ORIGINAL a sample whose alternate
 *     this record dropped is called HOM-REF, even where its own subset likelihoods make the
 *     heterozygote the most likely genotype by 10 phred: `A/G` with `G` dropped comes out `A/A`
 *     beside PLs of `10,0,30`. The call and the likelihoods disagree in the record that is
 *     written;
 *   - THE AC AND AF ARE RECOMPUTED FROM THE SUBSET CALLS, so the same sample contributes to one
 *     record's count and not the other's;
 *   - ONE HET-NON-REF GENOTYPE STILL EMPTIES EVERYTHING, likelihoods included, which is the rule
 *     from the split meeting the rule from the subsetting;
 *   - THE GQ IS RECOMPUTED PER RECORD, so a sample confident about one alternate is not confident
 *     about the other;
 *   - THE AD IS SPLIT BY ALLELE, and the depth of the dropped alternate is simply gone rather than
 *     folded into the reference;
 *   - AND THE TRIMMING STILL HAPPENS AFTER ALL OF IT, so a record can move while its genotypes are
 *     being rewritten.
 *
 * Output:
 *
 *     in\t<label>\t<start>-<end>\t<alleles>\t<genotypes>
 *     out\t<label>\t<n>\t<start>-<end>\t<alleles>\t<attributes>\t<genotypes>
 *     error\t<label>\t<exception class>:<message>
 *
 * A genotype is printed as `sample=alleles|PL|GQ|AD|DP`, samples separated by `;`.
 *
 * Usage: SplitWithLikelihoodsDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class SplitWithLikelihoodsDump {

    /** One sample's call: two allele indices, PLs, GQ, AD and DP, any of the last four absent. */
    record Spec(String sample, int first, int second, int[] pl, Integer gq, int[] ad, Integer dp) {}

    public static void main(final String[] args) {
        System.out.println("# SplitWithLikelihoodsDump: splitting a genotyped record, from the reference");

        // Three alleles, two samples, each het for a different alternate.
        run("two-samples", 100, List.of("A", "C", "G"), List.of(
                new Spec("s0", 0, 1, new int[] {50, 0, 60, 40, 30, 70}, 50, new int[] {10, 12, 8}, 30),
                new Spec("s1", 0, 2, new int[] {45, 35, 65, 0, 25, 55}, 25, new int[] {9, 4, 11}, 24)));

        // A sample that is hom for the second alternate, so one record loses its call entirely.
        run("hom-alt", 100, List.of("A", "C", "G"), List.of(
                new Spec("s0", 2, 2, new int[] {70, 60, 80, 30, 20, 0}, 20, new int[] {2, 3, 15}, 20)));

        // A het-non-ref, which empties every sample of every record.
        run("het-non-ref", 100, List.of("A", "C", "G"), List.of(
                new Spec("s0", 1, 2, new int[] {60, 30, 40, 20, 0, 50}, 20, new int[] {1, 8, 9}, 18),
                new Spec("s1", 0, 1, new int[] {50, 0, 60, 40, 30, 70}, 50, new int[] {10, 12, 8}, 30)));

        // Indels of different lengths, so the trimming moves the two records differently.
        run("indels", 100, List.of("ACGT", "A", "ACGTT"), List.of(
                new Spec("s0", 0, 1, new int[] {40, 0, 50, 30, 20, 60}, 20, new int[] {6, 7, 5}, 18)));

        // A sample with no likelihoods beside one that has them.
        run("mixed-samples", 100, List.of("A", "C", "G"), List.of(
                new Spec("s0", 0, 1, new int[] {50, 0, 60, 40, 30, 70}, 50, new int[] {10, 12, 8}, 30),
                new Spec("s1", 0, 2, null, 35, null, 15)));

        // Four alleles, so a record drops two alternates rather than one.
        run("four-alleles", 100, List.of("A", "C", "G", "T"), List.of(
                new Spec("s0", 0, 3, new int[] {60, 50, 70, 40, 30, 80, 0, 20, 45, 55}, 20,
                        new int[] {5, 4, 3, 12}, 24)));
    }

    static void run(final String label, final int start, final List<String> bases,
                    final List<Spec> specs) {
        final List<Allele> alleles = new ArrayList<>();
        for (int index = 0; index < bases.size(); index++) {
            alleles.add(Allele.create(bases.get(index), index == 0));
        }
        final List<Genotype> genotypes = new ArrayList<>();
        for (final Spec spec : specs) {
            final GenotypeBuilder builder = new GenotypeBuilder(spec.sample(),
                    List.of(alleles.get(spec.first()), alleles.get(spec.second())));
            if (spec.pl() != null) {
                builder.PL(spec.pl());
            }
            if (spec.gq() != null) {
                builder.GQ(spec.gq());
            }
            if (spec.ad() != null) {
                builder.AD(spec.ad());
            }
            if (spec.dp() != null) {
                builder.DP(spec.dp());
            }
            genotypes.add(builder.make());
        }
        run(label, new VariantContextBuilder("x", "chr1", start,
                start + bases.get(0).length() - 1, alleles).genotypes(genotypes).make());
    }

    static void run(final String label, final VariantContext variant) {
        System.out.printf("in\t%s\t%s\t%s\t%s%n", label, place(variant), alleles(variant),
                genotypes(variant));
        final List<VariantContext> results;
        try {
            results = GATKVariantContextUtils.splitVariantContextToBiallelics(variant, false,
                    GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, false);
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        for (int index = 0; index < results.size(); index++) {
            final VariantContext result = results.get(index);
            System.out.printf("out\t%s\t%d\t%s\t%s\t%s\t%s%n", label, index, place(result),
                    alleles(result), attributes(result), genotypes(result));
        }
    }

    static String place(final VariantContext variant) {
        return variant.getStart() + "-" + variant.getEnd();
    }

    static String alleles(final VariantContext variant) {
        final List<String> out = new ArrayList<>();
        for (final Allele allele : variant.getAlleles()) {
            out.add(allele.getDisplayString() + (allele.isReference() ? "(ref)" : ""));
        }
        return String.join(",", out);
    }

    static String attributes(final VariantContext variant) {
        final Map<String, Object> sorted = new LinkedHashMap<>();
        variant.getAttributes().keySet().stream().sorted()
                .forEach(key -> sorted.put(key, variant.getAttribute(key)));
        final List<String> out = new ArrayList<>();
        for (final Map.Entry<String, Object> entry : sorted.entrySet()) {
            out.add(entry.getKey() + "=" + String.valueOf(entry.getValue()));
        }
        return String.join(";", out);
    }

    static String genotypes(final VariantContext variant) {
        final List<String> out = new ArrayList<>();
        for (final Genotype genotype : variant.getGenotypes()) {
            final List<String> called = new ArrayList<>();
            for (final Allele allele : genotype.getAlleles()) {
                called.add(allele.getDisplayString());
            }
            out.add(genotype.getSampleName() + "=" + String.join("/", called)
                    + "|" + (genotype.hasPL() ? ints(genotype.getPL()) : "")
                    + "|" + (genotype.hasGQ() ? String.valueOf(genotype.getGQ()) : "")
                    + "|" + (genotype.hasAD() ? ints(genotype.getAD()) : "")
                    + "|" + (genotype.hasDP() ? String.valueOf(genotype.getDP()) : ""));
        }
        return String.join(";", out);
    }

    static String ints(final int[] values) {
        final List<String> out = new ArrayList<>();
        for (final int value : values) {
            out.add(String.valueOf(value));
        }
        return String.join(",", out);
    }
}

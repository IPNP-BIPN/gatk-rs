/*
 * GATKVariantContextUtils.splitVariantContextToBiallelics, taken from the reference.
 *
 * What --split-multiallelics does in LeftAlignAndTrimVariants: one record with N alternates becomes
 * N records with one each, every one of them right trimmed.
 *
 * Nine behaviours this is built to catch.
 *
 *   - A NON-VARIANT RECORD BECOMES AN EMPTY LIST, not a list of one, so a record with no alternate
 *     at all disappears rather than passing through;
 *   - AND A BIALLELIC ONE IS RETURNED AS ITSELF, the same object, untrimmed, so the trimming that
 *     every split record gets is not applied to a record that did not need splitting;
 *   - ONE HET-NON-REF GENOTYPE ANYWHERE CHANGES THE RULE FOR EVERY RECORD, `hasHetNonRef` swapping
 *     the assignment method to SET_TO_NO_CALL_NO_ANNOTATIONS for all of them, so a single 1/2 call
 *     in one sample empties the calls of every sample and every output record;
 *   - THE ATTRIBUTE FILTER IS A NEGATED DISJUNCTION, `!(AC || AF || AN) || method ==
 *     SET_TO_NO_CALL_NO_ANNOTATIONS`, so AC, AF and AN survive a normal split and are removed
 *     along with everything else once that method is in force;
 *   - AND THE SURVIVING AC, AF AND AN ARE THEN RECOMPUTED FROM THE GENOTYPES, so their values
 *     change rather than being carried over, and A RECORD WITH NO GENOTYPES LOSES THEM ALTOGETHER:
 *     the three that the filter kept are dropped by the recomputation that has nothing to count;
 *   - EVERY OUTPUT IS RIGHT TRIMMED and only left trimmed when asked, which is trimAlleles with
 *     trimReverse always true;
 *   - SO TWO ALTERNATES OF DIFFERENT LENGTHS COME OUT AT DIFFERENT POSITIONS, the trim being
 *     computed for each pair on its own;
 *   - THE ORDER IS THE ALTERNATE ORDER of the input, whatever the trimming does to the positions;
 *   - AND THE GENOTYPES ARE SUBSET, not merely remapped: with no likelihoods to go on,
 *     BEST_MATCH_TO_ORIGINAL turns an allele that is not in the pair into the REFERENCE, so `A/G`
 *     becomes `A/A` in the record that kept `C`. Nothing here carries PLs, so the reindexing that
 *     AlleleSubsettingUtils does to them is not measured by this suite.
 *
 * Output:
 *
 *     in\t<label>\t<start>-<end>\t<alleles>\t<attributes>\t<genotypes>
 *     out\t<label>\t<n>\t<start>-<end>\t<alleles>\t<attributes>\t<genotypes>
 *     count\t<label>\t<how many records came back>
 *     same\t<label>\t<true when the one record back is the very same object>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: SplitBiallelicsDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class SplitBiallelicsDump {

    public static void main(final String[] args) {
        System.out.println("# SplitBiallelicsDump: splitVariantContextToBiallelics, from the reference");

        // The two shapes that never split.
        run("no-alternate", record(100, List.of("A")), false);
        run("biallelic", record(100, List.of("AC", "A")), false);

        // Two alternates of the same length, then of different lengths, so the trim differs.
        run("two-snps", record(100, List.of("A", "C", "G")), false);
        run("two-deletions", record(100, List.of("ACGT", "ACG", "A")), false);
        run("deletion-and-insertion", record(100, List.of("AC", "A", "ACGG")), false);
        // Three alternates, so the order of the output can be seen.
        run("three-alternates", record(100, List.of("AC", "A", "ACG", "G")), false);

        // Left trimming asked for, which can move a record forward.
        run("left-trim", record(100, List.of("AACG", "AAC", "AACGT")), true);
        run("no-left-trim", record(100, List.of("AACG", "AAC", "AACGT")), false);

        // Attributes: AC, AF and AN survive a normal split and are recomputed; anything else goes.
        run("attributes", withAttributes(record(100, List.of("A", "C", "G")),
                Map.of("AC", List.of(1, 2), "AF", List.of(0.25, 0.5), "AN", 4, "DP", 30)), false);

        // Genotypes without likelihoods, which is what a hand written vcf carries.
        run("called-genotypes", withGenotypes(record(100, List.of("A", "C", "G")),
                List.of("0/1", "1/2", "2/2")), false);
        run("no-het-non-ref", withGenotypes(record(100, List.of("A", "C", "G")),
                List.of("0/1", "0/2", "0/0")), false);
        // And the same records with attributes, so the two rules interact.
        run("het-non-ref-with-attributes",
                withAttributes(withGenotypes(record(100, List.of("A", "C", "G")),
                        List.of("0/1", "1/2")), Map.of("AC", List.of(1, 2), "AN", 4, "DP", 30)),
                false);
    }

    /** The first string is the reference allele, the rest alternates. */
    static VariantContext record(final int start, final List<String> bases) {
        final List<Allele> alleles = new ArrayList<>();
        for (int index = 0; index < bases.size(); index++) {
            alleles.add(Allele.create(bases.get(index), index == 0));
        }
        return new VariantContextBuilder("x", "chr1", start,
                start + alleles.get(0).length() - 1, alleles).make();
    }

    static VariantContext withAttributes(final VariantContext variant,
                                         final Map<String, Object> attributes) {
        final VariantContextBuilder builder = new VariantContextBuilder(variant);
        for (final Map.Entry<String, Object> entry : attributes.entrySet()) {
            builder.attribute(entry.getKey(), entry.getValue());
        }
        return builder.make();
    }

    /** Genotypes given as `0/1`, indices into the record's allele list. */
    static VariantContext withGenotypes(final VariantContext variant, final List<String> calls) {
        final List<Genotype> genotypes = new ArrayList<>();
        for (int sample = 0; sample < calls.size(); sample++) {
            final List<Allele> called = new ArrayList<>();
            for (final String index : calls.get(sample).split("/")) {
                called.add(variant.getAlleles().get(Integer.parseInt(index)));
            }
            genotypes.add(new GenotypeBuilder("s" + sample, called).make());
        }
        return new VariantContextBuilder(variant).genotypes(genotypes).make();
    }

    static void run(final String label, final VariantContext variant, final boolean trimLeft) {
        System.out.printf("in\t%s\t%s\t%s\t%s\t%s%n", label, place(variant), alleles(variant),
                attributes(variant), genotypes(variant));
        final List<VariantContext> results;
        try {
            results = GATKVariantContextUtils.splitVariantContextToBiallelics(variant, trimLeft,
                    GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, false);
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("count\t%s\t%d%n", label, results.size());
        for (int index = 0; index < results.size(); index++) {
            final VariantContext result = results.get(index);
            System.out.printf("out\t%s\t%d\t%s\t%s\t%s\t%s%n", label, index, place(result),
                    alleles(result), attributes(result), genotypes(result));
        }
        System.out.printf("same\t%s\t%s%n", label,
                results.size() == 1 && results.get(0) == variant);
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

    /** The attributes in a fixed order, so the row does not depend on a map's iteration. */
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
            out.add(genotype.getSampleName() + "=" + String.join("/", called));
        }
        return String.join(";", out);
    }
}

/*
 * GATKVariantContextUtils.trimAlleles, taken from the reference.
 *
 * The other half of what LeftAlignAndTrimVariants does, and the piece the multiallelic splitter
 * calls on every biallelic it produces: cut the bases every allele shares, from the front, from the
 * back, or from both.
 *
 * Nine behaviours this is built to catch.
 *
 *   - ONE ALLELE OF LENGTH ONE STOPS EVERYTHING. `anyMatch(a -> a.length() == 1 && !a.equals(
 *     Allele.SPAN_DEL))` returns the input untouched, so a record carrying a snp beside a deletion
 *     is never trimmed at all, however much its other alleles share;
 *   - AND THE SPANNING DELETION IS THE ONE EXCEPTION to that test, being length one and excluded by
 *     name;
 *   - SYMBOLIC ALLELES AND `*` ARE NOT FED TO THE COMPARISON but are kept in the output, so a
 *     record can be trimmed on the strength of the alleles that are left;
 *   - THE FORWARD TRIM IS A NEGATIVE SHIFT, `startTrim = -shifts.getLeft()`, which is normalizeAlleles
 *     being asked to move alleles that start at index 0 and can only move backwards;
 *   - AN ALLELE TRIMMED TO NOTHING GETS ONE BASE BACK, at the END when nothing was trimmed from the
 *     front and at the START when something was, which is the only reason a vcf record always keeps
 *     an anchor base;
 *   - THE INNER trimAlleles TAKES AN INCLUSIVE INDEX, so the caller passes `startBasesToClip - 1`
 *     and -1 means "trim nothing from the front";
 *   - AND IT RETURNS THE INPUT ITSELF when that index is -1 and the reverse trim is 0;
 *   - THE START MOVES BY `fwdTrimEnd + 1` AND THE END IS RECOMPUTED FROM THE REFERENCE ALLELE's new
 *     length, not by subtracting the reverse trim;
 *   - THE ALLELE MAP IS A LinkedHashMap here, where leftAlignAndTrim uses a HashMap, so the output
 *     order is the input order and the genotypes are remapped through it.
 *
 * Output:
 *
 *     in\t<label>\t<start>-<end>\t<alleles>\t<genotypes>
 *     out\t<label>\t<start>-<end>\t<alleles>\t<genotypes>
 *     same\t<label>\t<true when the very same object came back>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: TrimAllelesDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.ArrayList;
import java.util.List;

public class TrimAllelesDump {

    public static void main(final String[] args) {
        System.out.println("# TrimAllelesDump: trimAlleles, from the reference");

        // Nothing shared beyond the anchor base, which is what a well formed record looks like.
        run("anchored-deletion", record(100, "AC", "A"), true, true);
        // A shared suffix, which only the reverse trim removes.
        run("shared-suffix", record(100, "ACGT", "AGT"), true, true);
        run("shared-suffix-forward-only", record(100, "ACGT", "AGT"), true, false);
        run("shared-suffix-reverse-only", record(100, "ACGT", "AGT"), false, true);
        // A shared prefix of more than one base, which only the forward trim removes.
        run("shared-prefix", record(100, "AACG", "AACGT"), true, true);
        run("shared-prefix-forward-only", record(100, "AACG", "AACGT"), true, false);
        run("shared-prefix-reverse-only", record(100, "AACG", "AACGT"), false, true);
        // Shared at both ends at once.
        run("shared-both-ends", record(100, "AACGTT", "AATT"), true, true);

        // An allele of length one anywhere stops the whole thing, however much the others share.
        run("one-base-allele", record(100, List.of("ACGT", "AGT", "A")), true, true);
        run("snp-beside-indel", record(100, List.of("AC", "A", "GC")), true, true);

        // The spanning deletion is length one and is excluded by name, so trimming still runs.
        run("spanning-deletion", record(100, List.of("ACGT", "AGT", "*")), true, true);
        // A symbolic allele is not compared either, and is kept.
        run("symbolic", record(100, List.of("ACGT", "AGT", "<DEL>")), true, true);

        // Alleles that trim to nothing: an insertion whose alleles share everything but the
        // inserted bases, in both directions.
        run("empty-after-reverse", record(100, "AAA", "AA"), true, true);
        run("empty-after-forward", record(100, "TAAA", "TAA"), true, true);

        // Genotypes, which are remapped through the map the alleles came out of.
        run("with-genotypes", withGenotypes(record(100, "ACGT", "AGT")), true, true);
        // And a genotype holding the spanning deletion, which maps to itself.
        run("genotypes-with-star",
                withStarGenotype(record(100, List.of("ACGT", "AGT", "*"))), true, true);

        // Neither direction asked for, which is the inner function's own early return.
        run("no-trim-requested", record(100, "ACGT", "AGT"), false, false);
    }

    static VariantContext record(final int start, final String reference, final String alternate) {
        return record(start, List.of(reference, alternate));
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

    static VariantContext withGenotypes(final VariantContext variant) {
        final List<Genotype> genotypes = List.of(
                new GenotypeBuilder("s1",
                        List.of(variant.getReference(), variant.getAlternateAllele(0))).make(),
                new GenotypeBuilder("s2",
                        List.of(variant.getAlternateAllele(0), variant.getAlternateAllele(0))).make());
        return new VariantContextBuilder(variant).genotypes(genotypes).make();
    }

    static VariantContext withStarGenotype(final VariantContext variant) {
        final List<Genotype> genotypes = List.of(new GenotypeBuilder("s1",
                List.of(variant.getReference(), Allele.SPAN_DEL)).make());
        return new VariantContextBuilder(variant).genotypes(genotypes).make();
    }

    static void run(final String label, final VariantContext variant, final boolean trimForward,
                    final boolean trimReverse) {
        System.out.printf("in\t%s\t%s\t%s\t%s%n", label, place(variant), alleles(variant),
                genotypes(variant));
        final VariantContext result;
        try {
            result = GATKVariantContextUtils.trimAlleles(variant, trimForward, trimReverse);
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("out\t%s\t%s\t%s\t%s%n", label, place(result), alleles(result),
                genotypes(result));
        System.out.printf("same\t%s\t%s%n", label, result == variant);
    }

    static String place(final VariantContext variant) {
        return variant.getStart() + "-" + variant.getEnd();
    }

    /** The alleles in the order the record carries them, the reference marked with a `*`. */
    static String alleles(final VariantContext variant) {
        final List<String> out = new ArrayList<>();
        for (final Allele allele : variant.getAlleles()) {
            out.add(allele.getDisplayString() + (allele.isReference() ? "(ref)" : ""));
        }
        return String.join(",", out);
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

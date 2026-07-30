/*
 * The two axes of every likelihood matrix, taken from the reference.
 *
 * AlleleLikelihoods is defined over an AlleleList and a SampleList, so anything surprising about
 * those two is inherited by every likelihood GATK computes. Three things are surprising.
 *
 *   - an indexed list is a SET. IndexedSet drops a duplicate silently and the survivor keeps the
 *     index of its FIRST occurrence, so a list built from [A, C, A] has two entries and nothing
 *     says so. A caller that built its list from a variant context with a repeated allele gets a
 *     matrix with fewer rows than it asked for;
 *   - membership is Allele.equals, which is bases AND the reference flag, so the reference A and
 *     a non-reference A are two different entries in the same list, and indexOfReference finds
 *     only the first entry whose flag is set;
 *   - a permutation is a directional SUBSET map. It may drop and reorder, never invent. It refuses
 *     a longer target and a target holding an unknown allele, with the same message both ways, and
 *     it reports isPartial / isNonPermuted / isKept, which callers use to skip work. permutation()
 *     short circuits on equality before constructing anything, which is the only path where a
 *     separately built but identical target comes back non-permuted without a scan.
 *
 * Output:
 *
 *     list\t<label>\t<size>\t<comma-separated alleles as A* or A>\t<indexOfReference>
 *     index\t<label>\t<query>\t<indexOfAllele>\t<containsAllele>
 *     samples\t<label>\t<size>\t<comma-separated names>
 *     perm\t<label>\t<isPartial>\t<isNonPermuted>\t<fromSize>\t<toSize>\t<fromIndex[]>\t<toIndex[]>\t<kept[]>
 *     perm\t<label>\tE:<class>:<message>
 *
 * Usage: AlleleListDump
 */

import htsjdk.variant.variantcontext.Allele;

import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.AlleleListPermutation;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;

import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;

public class AlleleListDump {

    static final Allele REF = Allele.create("A", true);
    /** The same bases, not flagged reference: a different key to equals and to hashCode. */
    static final Allele REF_BASES_AS_ALT = Allele.create("A", false);
    static final Allele ALT1 = Allele.create("C", false);
    static final Allele ALT2 = Allele.create("G", false);
    static final Allele ALT3 = Allele.create("T", false);
    static final Allele SECOND_REF = Allele.create("C", true);
    static final Allele NO_CALL = Allele.NO_CALL;
    static final Allele SPAN_DEL = Allele.SPAN_DEL;

    public static void main(final String[] args) {
        System.out.println("# AlleleListDump: the allele and sample axes, and the permutation");

        list("empty");
        list("one-ref", REF);
        list("ref-and-alt", REF, ALT1);
        list("alt-first", ALT1, REF);
        list("no-reference", ALT1, ALT2);
        // A duplicate, which is dropped and keeps its first index.
        list("duplicate-adjacent", REF, REF, ALT1);
        list("duplicate-separated", REF, ALT1, REF);
        list("all-duplicates", ALT1, ALT1, ALT1);
        // Same bases, different reference flag: two entries.
        list("ref-flag-pair", REF, REF_BASES_AS_ALT);
        // Two alleles both flagged reference: indexOfReference reports only the first.
        list("two-references", REF, SECOND_REF);
        // The special alleles, which are ordinary members here.
        list("with-no-call", REF, NO_CALL, ALT1);
        list("with-span-del", REF, SPAN_DEL, ALT1);

        // indexOfAllele on a list holding both flavours of A.
        final AlleleList<Allele> pair = new IndexedAlleleList<>(REF, REF_BASES_AS_ALT, ALT1);
        index("ref-flag-pair", pair, "A*", REF);
        index("ref-flag-pair", pair, "A", REF_BASES_AS_ALT);
        index("ref-flag-pair", pair, "C", ALT1);
        index("ref-flag-pair", pair, "G", ALT2);
        index("ref-flag-pair", pair, "no-call", NO_CALL);

        samples("empty");
        samples("one", "s1");
        samples("two", "s1", "s2");
        // The same set semantics as the allele list.
        samples("duplicate", "s1", "s2", "s1");
        samples("out-of-order", "b", "A", "a");

        // Permutations, from a four-allele original.
        final AlleleList<Allele> original = new IndexedAlleleList<>(REF, ALT1, ALT2, ALT3);
        perm("identity", original, new IndexedAlleleList<>(REF, ALT1, ALT2, ALT3));
        perm("reordered", original, new IndexedAlleleList<>(ALT3, ALT2, ALT1, REF));
        perm("swap-two", original, new IndexedAlleleList<>(REF, ALT2, ALT1, ALT3));
        perm("drop-last", original, new IndexedAlleleList<>(REF, ALT1, ALT2));
        perm("drop-first", original, new IndexedAlleleList<>(ALT1, ALT2, ALT3));
        perm("keep-one", original, new IndexedAlleleList<>(ALT2));
        perm("keep-none", original, new IndexedAlleleList<>());
        // A subset that happens to be in order but shorter: partial, and NOT non-permuted even
        // though every kept allele is where it started.
        perm("prefix", original, new IndexedAlleleList<>(REF, ALT1));
        // Refusals.
        perm("longer-target", original,
                new IndexedAlleleList<>(REF, ALT1, ALT2, ALT3, Allele.create("AA", false)));
        perm("unknown-allele", original, new IndexedAlleleList<>(REF, Allele.create("AA", false)));
        // The reference flag again: the target's A is not the original's A*.
        perm("wrong-ref-flag", original, new IndexedAlleleList<>(REF_BASES_AS_ALT, ALT1));
        // A target with a duplicate, which the target's own constructor collapses first, so this
        // is a permutation of size two and not a refusal.
        perm("duplicate-in-target", original, new IndexedAlleleList<>(ALT1, ALT1, ALT2));
        // The empty original, whose only legal target is empty.
        perm("empty-to-empty", new IndexedAlleleList<>(), new IndexedAlleleList<>());
        perm("empty-to-one", new IndexedAlleleList<>(), new IndexedAlleleList<>(REF));
    }

    static void list(final String label, final Allele... alleles) {
        final AlleleList<Allele> list = new IndexedAlleleList<>(alleles);
        final StringJoiner joiner = new StringJoiner(",");
        for (int i = 0; i < list.numberOfAlleles(); i++) {
            joiner.add(show(list.getAllele(i)));
        }
        System.out.printf("list\t%s\t%d\t%s\t%d%n", label, list.numberOfAlleles(), joiner,
                list.indexOfReference());
    }

    static void index(final String label, final AlleleList<Allele> list, final String query,
                      final Allele allele) {
        System.out.printf("index\t%s\t%s\t%d\t%b%n", label, query, list.indexOfAllele(allele),
                list.containsAllele(allele));
    }

    static void samples(final String label, final String... names) {
        final IndexedSampleList list = new IndexedSampleList(Arrays.asList(names));
        final StringJoiner joiner = new StringJoiner(",");
        for (int i = 0; i < list.numberOfSamples(); i++) {
            joiner.add(list.getSample(i));
        }
        System.out.printf("samples\t%s\t%d\t%s%n", label, list.numberOfSamples(), joiner);
    }

    static void perm(final String label, final AlleleList<Allele> from,
                     final AlleleList<Allele> to) {
        try {
            final AlleleListPermutation<Allele> permutation = from.permutation(to);
            final StringJoiner fromIndices = new StringJoiner(",");
            for (int i = 0; i < permutation.toSize(); i++) {
                fromIndices.add(Integer.toString(permutation.fromIndex(i)));
            }
            final StringJoiner toIndices = new StringJoiner(",");
            final StringJoiner kept = new StringJoiner(",");
            for (int i = 0; i < permutation.fromSize(); i++) {
                toIndices.add(Integer.toString(permutation.toIndex(i)));
                kept.add(Boolean.toString(permutation.isKept(i)));
            }
            System.out.printf("perm\t%s\t%b\t%b\t%d\t%d\t%s\t%s\t%s%n", label,
                    permutation.isPartial(), permutation.isNonPermuted(), permutation.fromSize(),
                    permutation.toSize(), fromIndices, toIndices, kept);
        } catch (final Exception | AssertionError e) {
            System.out.printf("perm\t%s\tE:%s:%s%n", label, e.getClass().getName(),
                    e.getMessage() == null ? "" : e.getMessage().replace('\n', ' '));
        }
    }

    /** `A*` for a reference allele, its bases otherwise, which is htsjdk's own rendering. */
    static String show(final Allele allele) {
        return allele.getDisplayString() + (allele.isReference() ? "*" : "");
    }

    static List<Allele> unused() {
        return List.of();
    }
}

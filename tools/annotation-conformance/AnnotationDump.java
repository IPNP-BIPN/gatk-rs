/*
 * Three INFO annotations, taken from the reference through their own interface.
 *
 * ChromosomeCounts, SampleList and RawGtCount are the members of the counting family that need no
 * AlleleLikelihoods, so they can be asked directly. What is measured is not only the numbers but
 * the SHAPE of what an annotation returns, which is where a port goes wrong without noticing:
 *
 *   - the map is Map<String, Object>, so the Java CLASS of each value is part of the answer.
 *     ChromosomeCounts puts an Integer for one alternate allele and an ArrayList for two or more,
 *     and both render the same way;
 *   - an empty map is not zeroes. It means the keys are ABSENT from the record, and every one of
 *     these annotations has at least one guard reaching that branch;
 *   - RawGtCount.annotate returns NULL rather than an empty map, which is a third thing again. Its
 *     real work is combineRawData, whose output drops one of the three counts it just summed and
 *     writes a literal "." in its place, so combining a value with itself is not that value
 *     doubled;
 *   - SampleList's guard is isMonomorphicInSamples, which is not "every genotype is hom-ref": a
 *     site with alternate alleles and no genotypes at all is not monomorphic, and a site where
 *     everybody is a no-call is.
 *
 * Output:
 *
 *     anno\t<annotation>\t<label>\t<key>=<value>[java.lang.X];...    (empty for an empty map, null for null)
 *     keys\t<annotation>\t<comma-separated getKeyNames()>
 *     combine\t<label>\t<combined value|E:class:message>
 *
 * Usage: AnnotationDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import org.broadinstitute.hellbender.tools.walkers.annotator.ChromosomeCounts;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.RawGtCount;
import org.broadinstitute.hellbender.tools.walkers.annotator.SampleList;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotationData;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

public class AnnotationDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele REF_BASES_AS_ALT = Allele.create("A", false);
    static final Allele ALT1 = Allele.create("C", false);
    static final Allele ALT2 = Allele.create("G", false);
    static final Allele NO_CALL = Allele.NO_CALL;

    public static void main(final String[] args) {
        System.out.println("# AnnotationDump: ChromosomeCounts, SampleList and RawGtCount");

        final InfoFieldAnnotation counts = new ChromosomeCounts();
        final InfoFieldAnnotation samples = new SampleList();
        final InfoFieldAnnotation rawGt = new RawGtCount();

        keys("ChromosomeCounts", counts);
        keys("SampleList", samples);
        keys("RawGtCount", rawGt);

        for (final Object[] fixture : fixtures()) {
            final String label = (String) fixture[0];
            final VariantContext vc = (VariantContext) fixture[1];
            annotate("ChromosomeCounts", label, counts, vc);
            annotate("SampleList", label, samples, vc);
            annotate("RawGtCount", label, rawGt, vc);
        }

        // combineRawData, which is the only path in RawGtCount that produces anything.
        combine("single", "[1, 2, 3]");
        combine("single-no-brackets", "1,2,3");
        combine("single-no-spaces", "[1,2,3]");
        combine("doubled", "[1, 2, 3]", "[1, 2, 3]");
        combine("three-ways", "[1, 2, 3]", "[10, 20, 30]", "[100, 200, 300]");
        combine("zeroes", "[0, 0, 0]");
        // The output of a combine, fed back in: the "." is not an integer, so this fails, which is
        // what makes the operation non-associative in practice and not only in principle.
        combine("round-trip", ".,2,3");
        // Whitespace: trim takes the ends, and the split takes a comma FOLLOWED by spaces.
        combine("leading-space", "  [1, 2, 3]  ");
        combine("space-before-comma", "1 , 2, 3");
        combine("many-spaces", "1,    2,    3");
        combine("tab-separated", "1,\t2,\t3");
        // Brackets anywhere, since replaceAll is global.
        combine("brackets-inside", "[1, 2], [3]");
        // Arities that are not three.
        combine("two-values", "[1, 2]");
        combine("four-values", "[1, 2, 3, 4]");
        combine("empty", "");
        combine("trailing-comma", "1,2,3,");
        // Not integers.
        combine("negative", "[-1, -2, -3]");
        combine("plus-sign", "[+1, 2, 3]");
        combine("not-a-number", "[a, 2, 3]");
        combine("overflow", "[2147483648, 2, 3]");
        combine("max-int", "[2147483647, 2147483647, 2147483647]");
        // Two maximum ints summed, which is a silent overflow rather than an error.
        combine("overflowing-sum", "[2147483647, 2147483647, 2147483647]",
                "[1, 1, 1]");
    }

    static Object[][] fixtures() {
        return new Object[][] {
            {"no-genotypes", build(List.of(REF, ALT1))},
            {"ref-only-no-genotypes", build(List.of(REF))},
            {"one-het", build(List.of(REF, ALT1), gt("s1", REF, ALT1))},
            {"one-hom-var", build(List.of(REF, ALT1), gt("s1", ALT1, ALT1))},
            {"all-hom-ref", build(List.of(REF, ALT1), gt("s1", REF, REF), gt("s2", REF, REF))},
            {"two-alts", build(List.of(REF, ALT1, ALT2),
                    gt("s1", REF, ALT1), gt("s2", ALT1, ALT2))},
            {"ref-only-site", build(List.of(REF), gt("s1", REF, REF))},
            {"all-no-call", build(List.of(REF, ALT1), gt("s1", NO_CALL, NO_CALL))},
            {"half-no-call", build(List.of(REF, ALT1), gt("s1", ALT1, NO_CALL))},
            // A partially called genotype is MIXED, and therefore called, so SampleList lists it.
            {"mixed-and-hom-ref", build(List.of(REF, ALT1),
                    gt("s1", REF, NO_CALL), gt("s2", ALT1, REF))},
            // Every non-hom-ref genotype filtered: the called counts ignore it, so the site is
            // monomorphic and SampleList says nothing even though a genotype carries an alt.
            {"filtered-het", build(List.of(REF, ALT1),
                    filtered(gt("s1", REF, ALT1), "LowGQ"), gt("s2", REF, REF))},
            // A genotype that is HET only because the reference bases are also present unflagged.
            {"het-by-ref-flag", build(List.of(REF, ALT1),
                    gt("s1", REF, REF_BASES_AS_ALT), gt("s2", REF, REF))},
            // Sample names out of storage order, to show the iteration order of SampleList.
            {"name-order", build(List.of(REF, ALT1),
                    gt("b", REF, ALT1), gt("A", REF, ALT1), gt("a", REF, ALT1),
                    gt("B", REF, REF))},
            {"three-hets", build(List.of(REF, ALT1),
                    gt("m0", REF, ALT1), gt("m1", REF, ALT1), gt("m2", REF, ALT1))},
        };
    }

    static Genotype gt(final String sample, final Allele... alleles) {
        return new GenotypeBuilder(sample, Arrays.asList(alleles)).make();
    }

    static Genotype filtered(final Genotype genotype, final String filter) {
        return new GenotypeBuilder(genotype).filters(filter).make();
    }

    static VariantContext build(final List<Allele> alleles, final Genotype... genotypes) {
        final VariantContextBuilder builder = new VariantContextBuilder()
                .chr("chr1").start(100).stop(100).alleles(alleles);
        if (genotypes.length > 0) {
            builder.genotypes(Arrays.asList(genotypes));
        }
        return builder.make();
    }

    static void keys(final String name, final InfoFieldAnnotation annotation) {
        System.out.printf("keys\t%s\t%s%n", name, String.join(",", annotation.getKeyNames()));
    }

    /** The reference context is null, which every annotation here accepts. */
    static void annotate(final String name, final String label,
                         final InfoFieldAnnotation annotation, final VariantContext vc) {
        try {
            final Map<String, Object> result = annotation.annotate(null, vc, null);
            if (result == null) {
                System.out.printf("anno\t%s\t%s\tnull%n", name, label);
                return;
            }
            final StringJoiner joiner = new StringJoiner(";");
            for (final Map.Entry<String, Object> entry : result.entrySet()) {
                final Object value = entry.getValue();
                joiner.add(String.format("%s=%s[%s]", entry.getKey(), render(value),
                        value == null ? "null" : value.getClass().getName()));
            }
            System.out.printf("anno\t%s\t%s\t%s%n", name, label, joiner);
        } catch (final Exception | AssertionError e) {
            System.out.printf("anno\t%s\t%s\tE:%s:%s%n", name, label, e.getClass().getName(),
                    oneLine(e.getMessage()));
        }
    }

    /** A Double travels as raw bits, since AF is a division. */
    static String render(final Object value) {
        if (value instanceof Double) {
            return Long.toString(Double.doubleToRawLongBits((Double) value));
        }
        if (value instanceof List) {
            final StringJoiner joiner = new StringJoiner(",", "(", ")");
            for (final Object element : (List<?>) value) {
                joiner.add(render(element));
            }
            return joiner.toString();
        }
        return String.valueOf(value);
    }

    static void combine(final String label, final String... rawValues) {
        try {
            final List<ReducibleAnnotationData<?>> data = new ArrayList<>();
            for (final String raw : rawValues) {
                data.add(new ReducibleAnnotationData<>(raw));
            }
            final Map<String, Object> combined =
                    new RawGtCount().combineRawData(List.of(REF, ALT1), data);
            System.out.printf("combine\t%s\t%s%n", label,
                    combined.get("RAW_GT_COUNT"));
        } catch (final Exception | AssertionError e) {
            System.out.printf("combine\t%s\tE:%s:%s%n", label, e.getClass().getName(),
                    oneLine(e.getMessage()));
        }
    }

    static String oneLine(final String message) {
        return message == null ? "" : message.replace('\n', ' ').replace('\t', ' ');
    }
}

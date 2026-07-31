/*
 * PossibleDeNovo and TransmittedSingleton, and the MendelianViolation under them, taken from the
 * reference.
 *
 *   - TransmittedSingleton reads the CHILD's depth three times and names the three results
 *     childIsHighDepth, momIsHighDepth and dadIsHighDepth. The parents' depths are never looked at,
 *     so the documented caveat about all three samples being deep is not what the code does;
 *   - Genotype.getDP() is -1 when absent, not 0, and PossibleDeNovo's depth threshold defaults to
 *     0, so a trio with no DP fails the high-confidence branch and falls into the low one;
 *   - isViolation indexes the child's alleles at 0 and 1 and asks each parent to CONTAIN one, so a
 *     half-called parent transmits whatever it does carry;
 *   - the de novo allele-frequency cutoff is max(4, nSamples/1000), so the flat four wins for every
 *     cohort below four thousand samples;
 *   - TransmittedSingleton's call-rate gate counts genotypes whose GQ is STRICTLY greater than 20,
 *     where every other threshold in the pair is greater-or-equal.
 *
 * Output:
 *
 *     type\t<label>\t<GenotypeType>
 *     violation\t<label>\t<true|false>
 *     denovo\t<label>\t<key>=<value>[<class>];...
 *     singleton\t<label>\t<key>=<value>[<class>];...
 *
 * Usage: PedigreeAnnotationDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import org.broadinstitute.hellbender.tools.walkers.annotator.PossibleDeNovo;
import org.broadinstitute.hellbender.tools.walkers.annotator.TransmittedSingleton;
import org.broadinstitute.hellbender.utils.samples.MendelianViolation;
import org.broadinstitute.hellbender.utils.samples.Sample;
import org.broadinstitute.hellbender.utils.samples.Sex;
import org.broadinstitute.hellbender.utils.samples.Trio;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

public class PedigreeAnnotationDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele ALT = Allele.create("C", false);
    static final int START = 105;

    public static void main(final String[] args) {
        System.out.println("# PedigreeAnnotationDump: hiConfDeNovo, loConfDeNovo and the singletons");

        for (final String label : new String[] {
                "hom-ref", "het", "hom-var", "no-call", "half-called", "empty"}) {
            genotypeType(label);
        }

        for (final String label : new String[] {
                "de-novo-het", "de-novo-hom-var", "inherited-het", "all-hom-ref",
                "mother-no-call", "father-no-call", "child-no-call", "mother-half-called",
                "both-parents-no-call", "hom-var-parent"}) {
            violation(label);
        }

        for (final String label : new String[] {
                "high-confidence", "low-confidence-gq", "no-depth", "shallow-child",
                "low-parent-gq", "inherited", "no-gq", "multiallelic", "common-allele",
                "two-trios"}) {
            denovo(label);
        }

        for (final String label : new String[] {
                "transmitted", "non-transmitted", "shallow-parents", "shallow-child",
                "low-call-rate", "ac-three", "child-hom-var", "both-parents-het"}) {
            singleton(label);
        }
    }

    static Genotype genotype(final String name, final List<Allele> alleles,
                             final Integer gq, final Integer dp) {
        final GenotypeBuilder gb = new GenotypeBuilder(name, alleles);
        if (gq != null) {
            gb.GQ(gq);
        }
        if (dp != null) {
            gb.DP(dp);
        }
        return gb.make();
    }

    static List<Allele> allelesFor(final String kind) {
        switch (kind) {
            case "hom-ref": return List.of(REF, REF);
            case "het": return List.of(REF, ALT);
            case "hom-var": return List.of(ALT, ALT);
            case "no-call": return List.of(Allele.NO_CALL, Allele.NO_CALL);
            case "half-called": return List.of(REF, Allele.NO_CALL);
            case "empty": return List.of();
            default: throw new IllegalArgumentException(kind);
        }
    }

    static void genotypeType(final String label) {
        final Genotype g = genotype("s", allelesFor(label), 30, 20);
        System.out.printf("type\t%s\t%s%n", label, g.getType());
    }

    /** mother, father, child, in that order. */
    static String[] trioKinds(final String label) {
        switch (label) {
            case "de-novo-het": return new String[] {"hom-ref", "hom-ref", "het"};
            case "de-novo-hom-var": return new String[] {"hom-ref", "hom-ref", "hom-var"};
            case "inherited-het": return new String[] {"het", "hom-ref", "het"};
            case "all-hom-ref": return new String[] {"hom-ref", "hom-ref", "hom-ref"};
            case "mother-no-call": return new String[] {"no-call", "hom-ref", "hom-var"};
            case "father-no-call": return new String[] {"hom-ref", "no-call", "hom-var"};
            case "child-no-call": return new String[] {"hom-ref", "hom-ref", "no-call"};
            case "mother-half-called": return new String[] {"half-called", "hom-ref", "het"};
            case "both-parents-no-call": return new String[] {"no-call", "no-call", "het"};
            case "hom-var-parent": return new String[] {"hom-var", "hom-ref", "hom-ref"};
            default: throw new IllegalArgumentException(label);
        }
    }

    static void violation(final String label) {
        final String[] kinds = trioKinds(label);
        final Genotype mom = genotype("mom", allelesFor(kinds[0]), 50, 30);
        final Genotype dad = genotype("dad", allelesFor(kinds[1]), 50, 30);
        final Genotype kid = genotype("kid", allelesFor(kinds[2]), 50, 30);
        System.out.printf("violation\t%s\t%b%n", label,
                MendelianViolation.isViolation(mom, dad, kid));
    }

    static Set<Trio> trios(final int count) {
        final Set<Trio> trios = new LinkedHashSet<>();
        for (int i = 0; i < count; i++) {
            final String suffix = i == 0 ? "" : Integer.toString(i);
            final Sample mother = new Sample("mom" + suffix, "fam" + i, null, null, Sex.FEMALE);
            final Sample father = new Sample("dad" + suffix, "fam" + i, null, null, Sex.MALE);
            final Sample child = new Sample("kid" + suffix, "fam" + i,
                    "dad" + suffix, "mom" + suffix, Sex.MALE);
            trios.add(new Trio(mother, father, child));
        }
        return trios;
    }

    static VariantContext denovoContext(final String label) {
        final List<Genotype> genotypes = new ArrayList<>();
        final List<Allele> alleles = "multiallelic".equals(label)
                ? List.of(REF, ALT, Allele.create("G", false)) : List.of(REF, ALT);
        switch (label) {
            case "high-confidence":
                genotypes.add(genotype("mom", allelesFor("hom-ref"), 50, 30));
                genotypes.add(genotype("dad", allelesFor("hom-ref"), 50, 30));
                genotypes.add(genotype("kid", allelesFor("het"), 50, 30));
                break;
            case "low-confidence-gq":
                genotypes.add(genotype("mom", allelesFor("hom-ref"), 5, 30));
                genotypes.add(genotype("dad", allelesFor("hom-ref"), 5, 30));
                genotypes.add(genotype("kid", allelesFor("het"), 15, 30));
                break;
            case "no-depth":
                genotypes.add(genotype("mom", allelesFor("hom-ref"), 50, null));
                genotypes.add(genotype("dad", allelesFor("hom-ref"), 50, null));
                genotypes.add(genotype("kid", allelesFor("het"), 50, null));
                break;
            case "shallow-child":
                genotypes.add(genotype("mom", allelesFor("hom-ref"), 50, 30));
                genotypes.add(genotype("dad", allelesFor("hom-ref"), 50, 30));
                genotypes.add(genotype("kid", allelesFor("het"), 50, 0));
                break;
            case "low-parent-gq":
                genotypes.add(genotype("mom", allelesFor("hom-ref"), 15, 30));
                genotypes.add(genotype("dad", allelesFor("hom-ref"), 15, 30));
                genotypes.add(genotype("kid", allelesFor("het"), 50, 30));
                break;
            case "inherited":
                genotypes.add(genotype("mom", allelesFor("het"), 50, 30));
                genotypes.add(genotype("dad", allelesFor("hom-ref"), 50, 30));
                genotypes.add(genotype("kid", allelesFor("het"), 50, 30));
                break;
            case "no-gq":
                genotypes.add(genotype("mom", allelesFor("hom-ref"), null, 30));
                genotypes.add(genotype("dad", allelesFor("hom-ref"), null, 30));
                genotypes.add(genotype("kid", allelesFor("het"), null, 30));
                break;
            case "multiallelic":
                genotypes.add(genotype("mom", allelesFor("hom-ref"), 50, 30));
                genotypes.add(genotype("dad", allelesFor("hom-ref"), 50, 30));
                genotypes.add(genotype("kid", allelesFor("het"), 50, 30));
                break;
            case "common-allele":
                genotypes.add(genotype("mom", allelesFor("hom-ref"), 50, 30));
                genotypes.add(genotype("dad", allelesFor("hom-ref"), 50, 30));
                genotypes.add(genotype("kid", allelesFor("het"), 50, 30));
                // Four more het samples, so the alternate is carried five times and clears the
                // flat cutoff of four.
                for (int i = 0; i < 4; i++) {
                    genotypes.add(genotype("x" + i, allelesFor("het"), 50, 30));
                }
                break;
            case "two-trios":
                genotypes.add(genotype("mom", allelesFor("hom-ref"), 50, 30));
                genotypes.add(genotype("dad", allelesFor("hom-ref"), 50, 30));
                genotypes.add(genotype("kid", allelesFor("het"), 50, 30));
                genotypes.add(genotype("mom1", allelesFor("hom-ref"), 50, 30));
                genotypes.add(genotype("dad1", allelesFor("hom-ref"), 50, 30));
                genotypes.add(genotype("kid1", allelesFor("het"), 15, 30));
                break;
            default: throw new IllegalArgumentException(label);
        }
        return new VariantContextBuilder().chr("chr1").start(START).stop(START)
                .alleles(alleles).genotypes(genotypes).make();
    }

    static void denovo(final String label) {
        final int trioCount = "two-trios".equals(label) ? 2 : 1;
        final PossibleDeNovo annotation =
                new PossibleDeNovo(trios(trioCount), PossibleDeNovo.DEFAULT_MIN_GENOTYPE_QUALITY_P);
        try {
            emitMap("denovo", label, annotation.annotate(null, denovoContext(label), null));
        } catch (final Exception | AssertionError e) {
            System.out.printf("denovo\t%s\tE:%s%n", label, e.getClass().getName());
        }
    }

    static VariantContext singletonContext(final String label) {
        final List<Genotype> genotypes = new ArrayList<>();
        int alleleCount = 2;
        switch (label) {
            case "transmitted":
                genotypes.add(genotype("mom", allelesFor("het"), 50, 30));
                genotypes.add(genotype("dad", allelesFor("hom-ref"), 50, 30));
                genotypes.add(genotype("kid", allelesFor("het"), 50, 30));
                break;
            case "non-transmitted":
                genotypes.add(genotype("mom", allelesFor("het"), 50, 30));
                genotypes.add(genotype("dad", allelesFor("hom-ref"), 50, 30));
                genotypes.add(genotype("kid", allelesFor("hom-ref"), 50, 30));
                alleleCount = 1;
                break;
            case "shallow-parents":
                // The parents are shallow and the child is deep: the reference passes this,
                // because all three depth tests read the child.
                genotypes.add(genotype("mom", allelesFor("het"), 50, 1));
                genotypes.add(genotype("dad", allelesFor("hom-ref"), 50, 1));
                genotypes.add(genotype("kid", allelesFor("het"), 50, 30));
                break;
            case "shallow-child":
                genotypes.add(genotype("mom", allelesFor("het"), 50, 30));
                genotypes.add(genotype("dad", allelesFor("hom-ref"), 50, 30));
                genotypes.add(genotype("kid", allelesFor("het"), 50, 1));
                break;
            case "low-call-rate":
                genotypes.add(genotype("mom", allelesFor("het"), 50, 30));
                genotypes.add(genotype("dad", allelesFor("hom-ref"), 50, 30));
                genotypes.add(genotype("kid", allelesFor("het"), 50, 30));
                for (int i = 0; i < 7; i++) {
                    genotypes.add(genotype("x" + i, allelesFor("hom-ref"), 5, 30));
                }
                break;
            case "ac-three":
                genotypes.add(genotype("mom", allelesFor("het"), 50, 30));
                genotypes.add(genotype("dad", allelesFor("hom-ref"), 50, 30));
                genotypes.add(genotype("kid", allelesFor("het"), 50, 30));
                alleleCount = 3;
                break;
            case "child-hom-var":
                genotypes.add(genotype("mom", allelesFor("het"), 50, 30));
                genotypes.add(genotype("dad", allelesFor("hom-ref"), 50, 30));
                genotypes.add(genotype("kid", allelesFor("hom-var"), 50, 30));
                break;
            case "both-parents-het":
                genotypes.add(genotype("mom", allelesFor("het"), 50, 30));
                genotypes.add(genotype("dad", allelesFor("het"), 50, 30));
                genotypes.add(genotype("kid", allelesFor("het"), 50, 30));
                break;
            default: throw new IllegalArgumentException(label);
        }
        return new VariantContextBuilder().chr("chr1").start(START).stop(START)
                .alleles(List.of(REF, ALT)).genotypes(genotypes)
                .attribute("AC", alleleCount).make();
    }

    static void singleton(final String label) {
        final TransmittedSingleton annotation = new TransmittedSingleton(trios(1));
        try {
            emitMap("singleton", label, annotation.annotate(null, singletonContext(label), null));
        } catch (final Exception | AssertionError e) {
            System.out.printf("singleton\t%s\tE:%s%n", label, e.getClass().getName());
        }
    }

    static void emitMap(final String kind, final String label, final Map<String, Object> result) {
        final StringJoiner joiner = new StringJoiner(";");
        if (result != null) {
            for (final Map.Entry<String, Object> entry : result.entrySet()) {
                joiner.add(String.format("%s=%s[%s]", entry.getKey(), entry.getValue(),
                        entry.getValue().getClass().getName()));
            }
        }
        System.out.printf("%s\t%s\t%s%n", kind, label, joiner);
    }
}

/*
 * GnarlyGenotyper's finalized genotypes, taken from the reference.
 *
 * The tool reads what GenomicsDB combined and finishes the job: it decides whether a site is worth
 * calling at all, recomputes QD and MQ from the raw annotations, and rewrites the genotypes. The
 * ENGINE IS ASKED DIRECTLY here rather than through the tool, which is what unblocks the
 * measurement: a hand-written combined GVCF never gets past the tool's own reader, but
 * `finalizeGenotype` is public and takes a variant context.
 *
 * Twelve behaviours this is built to catch.
 *
 *   - THE QUALITY FLOOR IS NOT THE CONFIDENCE ARGUMENT: it is that argument LESS ten times the
 *     logarithm of the site prior, so the SNP floor is exactly 60 and the indel floor about
 *     69.03, neither of them the 30 the argument names;
 *   - A SITE UNDER ITS FLOOR IS DROPPED ENTIRELY, the engine returning null, unless
 *     --keep-all-sites asks for it back;
 *   - A SITE KEPT THAT WAY COMES BACK FILTERED `LowQual` WITH AN ADJUSTED ALLELE COUNT OF ZERO
 *     and none of the other annotations: no QD, no AC, no QUAL, and its `<NON_REF>` still in the
 *     alternates;
 *   - WHICH FLOOR APPLIES TURNS ON WHETHER ANY ALTERNATE IS THE REFERENCE'S OWN LENGTH, so a site
 *     with one SNP and one indel is judged as a SNP;
 *   - A SPANNING DELETION DOES NOT COUNT AS A SNP for that test;
 *   - QUALapprox IS READ FROM THE PLAIN KEY WHEN THERE IS ONE and summed from the allele-specific
 *     list otherwise, and a site with neither scores zero and is always under the floor;
 *   - A CALLED SITE LOSES ITS `<NON_REF>` and gains AC, AF, AN, QD, FS, SOR, ExcessHet and MQ;
 *   - THE SITE'S QUAL IS `QUALapprox / 10 - log10(prior)`, which for a QUALapprox of 900 is 93
 *     and is written out as 870 after the Phred scaling: it is not QUALapprox;
 *   - QD IS QUALapprox OVER THE VARIANT DEPTH, so 900 over 30 is 30;
 *   - MQ IS FINALIZED FROM THE RAW SUM AND ITS DEPTH, and a raw key whose depth is zero drops MQ
 *     from the output exactly as having no raw key at all does;
 *   - --strip-allele-specific-annotations LEAVES THE AS_ KEYS PRESENT BUT NULL rather than
 *     removing them, so the output carries `AS_QD=null` either way;
 *   - AND AN ALLELE-SPECIFIC LIST WHOSE LENGTH DOES NOT MATCH THE ALLELES IS AN
 *     IllegalStateException, which is how a site that keeps its annotations and its `<NON_REF>`
 *     ends.
 *
 * Output:
 *
 *     in\t<case>\t<the input variant context, as its fields>
 *     out\t<case>\t<the finalized context, or `null`>
 *     db\t<case>\t<the annotation database context>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: GnarlyGenotyperDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.tools.walkers.gnarlyGenotyper.GnarlyGenotyperEngine;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class GnarlyGenotyperDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** One variant context, rendered as the fields the engine reads and writes. */
    static String render(final VariantContext vc) {
        if (vc == null) {
            return "null";
        }
        final StringBuilder text = new StringBuilder();
        text.append(vc.getContig()).append('\t').append(vc.getStart()).append('\t');
        text.append(vc.getReference().getDisplayString()).append('\t');
        final List<String> alts = new ArrayList<>();
        for (final Allele allele : vc.getAlternateAlleles()) {
            alts.add(allele.getDisplayString());
        }
        text.append(String.join(",", alts)).append('\t');
        text.append(vc.isFiltered() ? String.join(";", vc.getFilters()) : ".").append('\t');
        text.append(vc.hasLog10PError() ? String.format("%.4f", vc.getPhredScaledQual()) : ".")
                .append('\t');
        final List<String> info = new ArrayList<>();
        for (final Map.Entry<String, Object> entry
                : new java.util.TreeMap<>(vc.getAttributes()).entrySet()) {
            info.add(entry.getKey() + "=" + entry.getValue());
        }
        text.append(String.join(";", info)).append('\t');
        final List<String> genotypes = new ArrayList<>();
        for (final Genotype genotype : vc.getGenotypes()) {
            genotypes.add(genotype.getSampleName() + ":" + genotype.getGenotypeString()
                    + ":" + (genotype.hasPL() ? java.util.Arrays.toString(genotype.getPL()) : ".")
                    + ":" + (genotype.hasDP() ? genotype.getDP() : "."));
        }
        text.append(String.join(" ", genotypes));
        return text.toString();
    }

    /** A site with the attributes and alternates given. */
    static VariantContext site(final String reference, final List<String> alternates,
                               final Map<String, Object> attributes, final int[] pls,
                               final int depth) {
        final List<Allele> alleles = new ArrayList<>();
        alleles.add(Allele.create(reference, true));
        for (final String alternate : alternates) {
            alleles.add(Allele.create(alternate, false));
        }
        final GenotypeBuilder genotype = new GenotypeBuilder("sample1");
        genotype.alleles(List.of(alleles.get(0), alleles.get(0)));
        if (pls != null) {
            genotype.PL(pls);
        }
        if (depth >= 0) {
            genotype.DP(depth);
        }
        return new VariantContextBuilder("fixture", "chr1", 1000,
                1000 + reference.length() - 1L, alleles)
                .attributes(attributes)
                .genotypes(genotype.make())
                .make();
    }

    static Map<String, Object> attributes(final Object... pairs) {
        final Map<String, Object> map = new LinkedHashMap<>();
        for (int i = 0; i < pairs.length; i += 2) {
            map.put((String) pairs[i], pairs[i + 1]);
        }
        return map;
    }

    static void run(final String name, final VariantContext input, final boolean keepAllSites,
                    final int maxAlternates, final boolean stripAlleleSpecific) {
        emit("in", name, render(input));
        final GnarlyGenotyperEngine engine =
                new GnarlyGenotyperEngine(keepAllSites, maxAlternates, stripAlleleSpecific);
        final VariantContextBuilder database = new VariantContextBuilder(input);
        try {
            final VariantContext out = engine.finalizeGenotype(input, database);
            emit("out", name, render(out));
            emit("db", name, render(database.make()));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            emit("error", name, cause.getClass().getName() + ":" + cause.getMessage());
        }
    }

    public static void main(final String[] args) {
        System.out.println("# GnarlyGenotyperDump: the engine's own rules, asked directly");

        final int[] pls = {900, 0, 1200};

        // A SNP well above its floor.
        run("snp-called",
                site("A", List.of("G", "<NON_REF>"),
                        attributes("QUALapprox", 900, "VarDP", 30, "RAW_MQandDP", "108000,30"),
                        pls, 30),
                false, 6, false);
        // The same site with the quality just under the SNP floor, dropped and kept.
        run("snp-under-floor",
                site("A", List.of("G", "<NON_REF>"),
                        attributes("QUALapprox", 20, "VarDP", 30, "RAW_MQandDP", "108000,30"),
                        pls, 30),
                false, 6, false);
        run("snp-under-floor-kept",
                site("A", List.of("G", "<NON_REF>"),
                        attributes("QUALapprox", 20, "VarDP", 30, "RAW_MQandDP", "108000,30"),
                        pls, 30),
                true, 6, false);
        // A quality between the two floors, which a SNP passes and an indel does not.
        run("between-the-floors-snp",
                site("A", List.of("G", "<NON_REF>"),
                        attributes("QUALapprox", 40, "VarDP", 30, "RAW_MQandDP", "108000,30"),
                        pls, 30),
                true, 6, false);
        run("between-the-floors-indel",
                site("A", List.of("AT", "<NON_REF>"),
                        attributes("QUALapprox", 40, "VarDP", 30, "RAW_MQandDP", "108000,30"),
                        pls, 30),
                true, 6, false);
        // A site with one SNP and one indel, which is judged as a SNP.
        run("mixed-site",
                site("A", List.of("G", "AT", "<NON_REF>"),
                        attributes("QUALapprox", 40, "VarDP", 30, "RAW_MQandDP", "108000,30"),
                        new int[] {900, 0, 1200, 1500, 1800, 2100}, 30),
                true, 6, false);
        // A spanning deletion, which does not count as a SNP however long it is.
        run("spanning-deletion",
                site("A", List.of("*", "<NON_REF>"),
                        attributes("QUALapprox", 40, "VarDP", 30, "RAW_MQandDP", "108000,30"),
                        pls, 30),
                true, 6, false);

        // The allele-specific quality list, summed.
        run("allele-specific-qual",
                site("A", List.of("G", "<NON_REF>"),
                        attributes("AS_QUALapprox", "|500|400", "VarDP", 30,
                                "RAW_MQandDP", "108000,30"),
                        pls, 30),
                false, 6, false);
        // Neither key at all, which scores zero and is always under the floor.
        run("no-qual-key",
                site("A", List.of("G", "<NON_REF>"),
                        attributes("VarDP", 30, "RAW_MQandDP", "108000,30"), pls, 30),
                true, 6, false);

        // The allele-specific annotations stripped.
        run("stripped-annotations",
                site("A", List.of("G", "<NON_REF>"),
                        attributes("QUALapprox", 900, "AS_QUALapprox", "|500|400",
                                "AS_SB_TABLE", "10,10|5,5", "VarDP", 30,
                                "RAW_MQandDP", "108000,30"),
                        pls, 30),
                false, 6, true);
        run("kept-annotations",
                site("A", List.of("G", "<NON_REF>"),
                        attributes("QUALapprox", 900, "AS_QUALapprox", "|500|400",
                                "AS_SB_TABLE", "10,10|5,5", "VarDP", 30,
                                "RAW_MQandDP", "108000,30"),
                        pls, 30),
                false, 6, false);

        // More alternates than the cap allows.
        run("too-many-alternates",
                site("A", List.of("G", "C", "T", "<NON_REF>"),
                        attributes("QUALapprox", 900, "VarDP", 30, "RAW_MQandDP", "108000,30"),
                        new int[] {900, 0, 1200, 1500, 1800, 2100, 2400, 2700, 3000, 3300}, 30),
                false, 1, false);

        // A site with the raw MQ key and no depth beside it.
        run("raw-mq-without-depth",
                site("A", List.of("G", "<NON_REF>"),
                        attributes("QUALapprox", 900, "VarDP", 30, "RAW_MQandDP", "108000,0"),
                        pls, 30),
                false, 6, false);
        // And one with no raw MQ at all.
        run("no-raw-mq",
                site("A", List.of("G", "<NON_REF>"),
                        attributes("QUALapprox", 900, "VarDP", 30), pls, 30),
                false, 6, false);

        System.out.print(buf);
    }
}

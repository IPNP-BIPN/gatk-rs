/*
 * Mutect2FilteringEngine as constructed, taken from the reference.
 *
 * What the engine knows before any record reaches it: which samples are normal, what the threshold
 * starts at, and what a somatic clustering model with nothing learned yet believes. Five behaviours
 * this is built to catch.
 *
 *   - THE NORMAL SAMPLES ARE A HEADER LINE, `##normal_sample=`, and EVERY OTHER SAMPLE IS A TUMOUR
 *     SAMPLE. `isTumor` is `!isNormal`, so a sample the header never mentions — including one that
 *     is not in the VCF at all — is treated as tumour;
 *   - THE KEY IS MATCHED EXACTLY, so a header line whose key differs in case names no normal sample;
 *   - THE THRESHOLD STARTS AT THE ARGUMENT COLLECTION'S DEFAULT and not at anything learned;
 *   - A CLUSTERING MODEL WITH NO DATA STILL HAS PRIORS, one for a somatic variant and one for a
 *     variant against an artifact, and they are what every filter's posterior is weighed against on
 *     the first pass;
 *   - AND THE MISSING STATS FILE IS NOT AN ERROR: the engine reads it only `if (exists)`, so a run
 *     without one starts from an empty list rather than refusing.
 *
 * Output:
 *
 *     sample\t<label>\t<normal|tumour>
 *     value\t<label>\t<a double>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: MutectEngineConstructionDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine;

import java.io.File;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

public class MutectEngineConstructionDump {

    public static void main(final String[] args) {
        System.out.println("# MutectEngineConstructionDump: what the engine knows before any record");

        final Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        lines.add(new VCFHeaderLine("normal_sample", "N1"));
        lines.add(new VCFHeaderLine("normal_sample", "N2"));
        // A key that differs only in case, which the engine compares exactly.
        lines.add(new VCFHeaderLine("Normal_Sample", "N3"));
        final VCFHeader header = new VCFHeader(lines, List.of("T1", "N1", "N2", "N3"));

        final Mutect2FilteringEngine engine =
                new Mutect2FilteringEngine(new M2FiltersArgumentCollection(), header,
                        new File("no-such-stats-file.tsv"));

        for (final String sample : new String[] {"T1", "N1", "N2", "N3", "never-mentioned"}) {
            final Genotype genotype = new GenotypeBuilder(sample, List.of(Allele.REF_A)).make();
            System.out.printf("sample\t%s\t%s%n", sample,
                    engine.isNormal(genotype) ? "normal" : "tumour");
        }

        System.out.printf("value\tinitial-threshold\t%s%n",
                Double.toString(engine.getThreshold()));
        System.out.printf("value\tlog-prior-variant-versus-artifact\t%s%n",
                Double.toString(engine.getSomaticClusteringModel().getLogPriorOfVariantVersusArtifact()));

        final VariantContext snp = new VariantContextBuilder("dump", "chr1", 100, 100,
                List.of(Allele.REF_A, Allele.ALT_C)).make();
        System.out.printf("value\tlog-somatic-prior-snp\t%s%n",
                Double.toString(engine.getLogSomaticPrior(snp, 0)));

        final VariantContext indel = new VariantContextBuilder("dump", "chr1", 100, 102,
                List.of(Allele.create("ATT", true), Allele.create("A", false))).make();
        System.out.printf("value\tlog-somatic-prior-indel\t%s%n",
                Double.toString(engine.getLogSomaticPrior(indel, 0)));

        // The two posteriors the engine offers on top of its own priors.
        for (final double odds : new double[] {-10.0, 0.0, 10.0}) {
            System.out.printf("value\tposterior-normal-artifact-%s\t%s%n", Double.toString(odds),
                    Double.toString(engine.posteriorProbabilityOfNormalArtifact(odds)));
            System.out.printf("value\tposterior-error-snp-%s\t%s%n", Double.toString(odds),
                    Double.toString(engine.posteriorProbabilityOfError(snp, odds, 0)));
        }
    }
}

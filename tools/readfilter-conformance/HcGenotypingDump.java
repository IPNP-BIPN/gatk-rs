/*
 * HaplotypeCallerGenotypingEngine.assignGenotypeLikelihoods without annotations or phasing: the
 * calls HaplotypeCaller makes from a region's haplotypes and read likelihoods.
 *
 * At each event start in the window: the events, `replaceSpanDels`, the merged context, the allele
 * mapper, the read-by-allele matrix (`marginalize`), the reads overlapping the context widened by
 * `--informative-read-overlap-margin` (`retainEvidence`), each sample's genotype likelihoods from
 * that matrix (`IndependentSampleGenotypesModel`, `GenotypeLikelihoodCalculator`) stored as PLs,
 * `calculateGenotypes` (QUAL, alternates kept, LowQual, MLEAC/MLEAF, GT/GQ/PL), and the allele
 * reverse trim when alleles were dropped. The annotation engine has no annotations and physical
 * phasing is off; both are measured by their own suites.
 *
 * Output:
 *
 *     case\t<label>\t<samples,...>\t<ploidy>\t<window start>\t<window end>\t<standard confidence>\t<max alt alleles>\t<spanning genotyping>
 *     haplotype\t<label>\t<index>\t<ref>\t<bases>\t<cigar>\t<score bits>
 *     read\t<label>\t<sample>\t<name>\t<start>\t<length>
 *     lk\t<label>\t<sample index>\t<haplotype index>\t<values,...>
 *     call\t<label>\t<VCF line>
 *     calls\t<label>\t<count>
 *     error\t<label>\t<exception class>: <message>
 *
 * The reference is REF at chr1:1000; reads are 40-base 40M records whose position only decides
 * which sites they overlap.
 *
 * Usage: HcGenotypingDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.CalledHaplotypes;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class HcGenotypingDump {

    static final String REF = "GATTACACGTTGCAAGTCCGATAGCTTAGGCATCGATCGGATCCATGCAAGTCAGTAAAAAGTCCATGAC";
    static final int START = 1000;

    record Hap(String bases, String cigar, double score) {
    }

    record Read(String sample, int start, int support, int strength) {
    }

    public static void main(final String[] args) {
        System.out.println("# HcGenotypingDump: HaplotypeCallerGenotypingEngine.assignGenotypeLikelihoods");
        final int n = REF.length();
        final Hap ref = new Hap(REF, n + "M", -1.0);
        final Hap snp20 = new Hap(snp(REF, 20), n + "M", -2.0);
        final Hap snp20b = new Hap(snp(snp(REF, 20), 20), n + "M", -3.0);
        final Hap snp20and40 = new Hap(snp(snp(REF, 20), 40), n + "M", -2.5);
        final Hap del20 = new Hap(REF.substring(0, 21) + REF.substring(24), "21M3D" + (n - 24) + "M", -2.2);
        final Hap del18 = new Hap(REF.substring(0, 19) + REF.substring(23), "19M4D" + (n - 23) + "M", -2.4);
        final Hap ins30 = new Hap(REF.substring(0, 31) + "GG" + REF.substring(31), "31M2I" + (n - 31) + "M", -2.6);
        final Hap homDel = new Hap(REF.substring(0, 56) + REF.substring(58), "56M2D" + (n - 58) + "M", -2.1);
        final Hap homDel1 = new Hap(REF.substring(0, 56) + REF.substring(57), "56M1D" + (n - 57) + "M", -2.3);

        final List<String> one = List.of("sA");
        final List<String> two = List.of("sA", "sB");
        run("het-snp", one, 2, 1000, 1069, 30.0, 6, true, List.of(ref, snp20), reads("sA", new int[]{0, 1, 0, 1, 0, 1, 0, 1, 1, 0}, 8));
        run("hom-snp", one, 2, 1000, 1069, 30.0, 6, true, List.of(ref, snp20), reads("sA", new int[]{1, 1, 1, 1, 1, 1, 1, 1}, 8));
        run("hom-ref", one, 2, 1000, 1069, 30.0, 6, true, List.of(ref, snp20), reads("sA", new int[]{0, 0, 0, 0, 0, 0, 0, 0}, 8));
        run("low-qual", one, 2, 1000, 1069, 30.0, 6, true, List.of(ref, snp20), reads("sA", new int[]{0, 0, 0, 1, 0, 0}, 1));
        run("emit-threshold", one, 2, 1000, 1069, 10.0, 6, true, List.of(ref, snp20), reads("sA", new int[]{0, 0, 0, 1, 0, 0}, 1));
        run("triallelic", one, 2, 1000, 1069, 30.0, 6, true, List.of(ref, snp20, snp20b), reads("sA", new int[]{1, 2, 1, 2, 1, 2, 1, 2}, 6));
        run("triallelic-ploidy-3", one, 3, 1000, 1069, 30.0, 6, true, List.of(ref, snp20, snp20b), reads("sA", new int[]{0, 1, 2, 0, 1, 2, 0, 1, 2}, 6));
        run("two-sites", one, 2, 1000, 1069, 30.0, 6, true, List.of(ref, snp20, snp20and40), reads("sA", new int[]{0, 2, 0, 2, 0, 2, 1, 2}, 6));
        run("deletion", one, 2, 1000, 1069, 30.0, 6, true, List.of(ref, del20), reads("sA", new int[]{0, 1, 0, 1, 0, 1, 0, 1}, 6));
        run("snp-and-deletion", one, 2, 1000, 1069, 30.0, 6, true, List.of(ref, snp20, del20), reads("sA", new int[]{1, 2, 1, 2, 1, 2, 1, 2}, 6));
        run("spanning-deletion", one, 2, 1000, 1069, 30.0, 6, true, List.of(ref, snp20, del18), reads("sA", new int[]{1, 2, 1, 2, 1, 2, 1, 2}, 6));
        run("spanning-deletion-off", one, 2, 1000, 1069, 30.0, 6, false, List.of(ref, snp20, del18), reads("sA", new int[]{1, 2, 1, 2, 1, 2, 1, 2}, 6));
        run("insertion", one, 2, 1000, 1069, 30.0, 6, true, List.of(ref, ins30), reads("sA", new int[]{0, 1, 1, 1, 0, 1}, 6));
        run("homopolymer-deletions-trim", one, 2, 1000, 1069, 30.0, 6, true, List.of(ref, homDel, homDel1), reads("sA", new int[]{0, 2, 0, 2, 0, 2, 0, 2}, 6));
        run("max-alt-1", one, 2, 1000, 1069, 30.0, 1, true, List.of(ref, snp20, snp20b), reads("sA", new int[]{1, 2, 1, 1, 1, 2, 1, 1}, 6));
        run("two-samples", two, 2, 1000, 1069, 30.0, 6, true, List.of(ref, snp20, del20),
                concat(reads("sA", new int[]{0, 1, 0, 1, 0, 1}, 6), reads("sB", new int[]{2, 2, 2, 2, 2, 0}, 6)));
        run("window", one, 2, 1030, 1069, 30.0, 6, true, List.of(ref, snp20, ins30), reads("sA", new int[]{1, 2, 1, 2, 1, 2}, 6));
        run("distant-reads", one, 2, 1000, 1069, 30.0, 6, true, List.of(ref, snp20), farReads("sA"));
    }

    static List<Read> concat(final List<Read> a, final List<Read> b) {
        final List<Read> all = new ArrayList<>(a);
        all.addAll(b);
        return all;
    }

    static String snp(final String s, final int at) {
        final char c = s.charAt(at);
        return s.substring(0, at) + "ACGT".charAt(("ACGT".indexOf(c) + 1) % 4) + s.substring(at + 1);
    }

    /** One read per entry, best supporting that haplotype by `strength` (log10) over the rest. */
    static List<Read> reads(final String sample, final int[] support, final int strength) {
        final List<Read> out = new ArrayList<>();
        for (int i = 0; i < support.length; i++) {
            out.add(new Read(sample, 995 + (i * 5) % 30, support[i], strength));
        }
        return out;
    }

    /** Reads that end before the first event and start after the last. */
    static List<Read> farReads(final String sample) {
        return List.of(new Read(sample, 950, 1, 6), new Read(sample, 951, 1, 6), new Read(sample, 1060, 1, 6), new Read(sample, 1001, 1, 6));
    }

    static String bits(final double d) {
        return Long.toHexString(Double.doubleToRawLongBits(d));
    }

    static void run(final String label, final List<String> samples, final int ploidy, final int windowStart, final int windowEnd,
                    final double confidence, final int maxAlt, final boolean spanning, final List<Hap> haps, final List<Read> reads) {
        System.out.println("case\t" + label + "\t" + String.join(",", samples) + "\t" + ploidy + "\t" + windowStart + "\t" + windowEnd
                + "\t" + bits(confidence) + "\t" + maxAlt + "\t" + spanning);
        try {
            final SAMFileHeader header = new SAMFileHeader();
            header.setSequenceDictionary(new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", 5000))));
            for (final String s : samples) {
                final SAMReadGroupRecord g = new SAMReadGroupRecord("rg" + s);
                g.setSample(s);
                header.addReadGroup(g);
            }
            final List<Haplotype> haplotypes = new ArrayList<>();
            for (int i = 0; i < haps.size(); i++) {
                final Hap h = haps.get(i);
                final boolean isRef = h.bases().equals(REF);
                System.out.println("haplotype\t" + label + "\t" + i + "\t" + isRef + "\t" + h.bases() + "\t" + h.cigar() + "\t" + bits(h.score()));
                final Haplotype hap = new Haplotype(h.bases().getBytes(StandardCharsets.US_ASCII), isRef);
                hap.setCigar(TextCigarCodec.decode(h.cigar()));
                hap.setAlignmentStartHapwrtRef(0);
                hap.setGenomeLocation(new SimpleInterval("chr1", START, START + REF.length() - 1));
                hap.setScore(h.score());
                haplotypes.add(hap);
            }
            final Map<String, List<GATKRead>> evidence = new LinkedHashMap<>();
            final Map<String, List<Read>> specs = new LinkedHashMap<>();
            for (final String s : samples) {
                evidence.put(s, new ArrayList<>());
                specs.put(s, new ArrayList<>());
            }
            int count = 0;
            for (final Read r : reads) {
                final String name = r.sample() + "r" + count++;
                System.out.println("read\t" + label + "\t" + r.sample() + "\t" + name + "\t" + r.start() + "\t40");
                final GATKRead read = ArtificialReadUtils.createArtificialRead(header, name, 0, r.start(), 40);
                read.setReadGroup("rg" + r.sample());
                evidence.get(r.sample()).add(read);
                specs.get(r.sample()).add(r);
            }
            final AlleleLikelihoods<GATKRead, Haplotype> likelihoods = new AlleleLikelihoods<>(new IndexedSampleList(samples),
                    new IndexedAlleleList<>(haplotypes), evidence);
            for (int s = 0; s < samples.size(); s++) {
                final LikelihoodMatrix<GATKRead, Haplotype> m = likelihoods.sampleMatrix(s);
                final List<Read> sampleSpecs = specs.get(samples.get(s));
                for (int a = 0; a < haplotypes.size(); a++) {
                    final StringBuilder b = new StringBuilder();
                    for (int r = 0; r < sampleSpecs.size(); r++) {
                        final Read spec = sampleSpecs.get(r);
                        final double v = spec.support() == a ? -0.5 - 0.125 * (r % 3) : -0.5 - spec.strength() - 0.25 * ((a + r) % 4);
                        m.set(a, r, v);
                        b.append(r == 0 ? "" : ",").append(bits(v));
                    }
                    System.out.println("lk\t" + label + "\t" + s + "\t" + a + "\t" + b);
                }
            }
            final HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
            hcArgs.standardArgs.genotypeArgs.samplePloidy = ploidy;
            hcArgs.standardArgs.genotypeArgs.standardConfidenceForCalling = confidence;
            hcArgs.standardArgs.genotypeArgs.maxAlternateAlleles = maxAlt;
            hcArgs.disableSpanningEventGenotyping = !spanning;
            final HaplotypeCallerGenotypingEngine engine = new HaplotypeCallerGenotypingEngine(hcArgs, new IndexedSampleList(samples), false, false);
            engine.setAnnotationEngine(new VariantAnnotatorEngine(Collections.emptyList(), null, Collections.emptyList(), false, false));
            final CalledHaplotypes called = engine.assignGenotypeLikelihoods(haplotypes, likelihoods, new LinkedHashMap<>(),
                    REF.getBytes(StandardCharsets.US_ASCII), new SimpleInterval("chr1", START, START + REF.length() - 1),
                    new SimpleInterval("chr1", windowStart, windowEnd), new FeatureContext(), Collections.emptyList(), false, 0, header,
                    false, new HashSet<>(), null);
            final VCFEncoder encoder = new VCFEncoder(new VCFHeader(Collections.emptySet(), samples), true, false);
            for (final VariantContext vc : called.getCalls()) {
                System.out.println("call\t" + label + "\t" + encoder.encode(vc));
            }
            System.out.println("calls\t" + label + "\t" + called.getCalls().size());
        } catch (final Exception e) {
            System.out.println("error\t" + label + "\t" + e.getClass().getSimpleName() + ": " + e.getMessage());
        }
    }
}

/*
 * HaplotypeCaller's physical phasing: AssemblyBasedCallerUtils.phaseCalls over the calls of
 * assignGenotypeLikelihoods (phasing on, no annotations).
 *
 * A call with exactly one site-specific alternate is tied to the called haplotypes whose event map
 * holds that alternate at its start; two calls on the same haplotypes (or one on every haplotype
 * that carries a called variant) are phased 0|1 with each other, two on complementary haplotypes
 * 0|1 against 1|0, and a call whose partner is already in another group empties the whole
 * mapping. Each group's genotypes become phased, a het one reversed when its phased alternate
 * index holds no site-specific alternate, with PID (start_ref_alt of the group's first call),
 * PGT (0|1 or 1|0) and PS (that call's start).
 *
 * Output: as HcGenotypingDump (case, haplotype, read, lk, call, calls, error rows).
 *
 * Usage: HcPhasingDump
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

public class HcPhasingDump {

    static final String REF = "GATTACACGTTGCAAGTCCGATAGCTTAGGCATCGATCGGATCCATGCAAGTCAGTAAAAAGTCCATGAC";
    static final int START = 1000;

    record Hap(String bases, String cigar, double score) {
    }

    record Read(String sample, int start, int support, int strength) {
    }

    public static void main(final String[] args) {
        System.out.println("# HcPhasingDump: HaplotypeCallerGenotypingEngine.assignGenotypeLikelihoods");
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
        final Hap snp40 = new Hap(snp(REF, 40), n + "M", -2.7);
        final Hap snp20and40and60 = new Hap(snp(snp(snp(REF, 20), 40), 60), n + "M", -2.8);
        final Hap snp60 = new Hap(snp(REF, 60), n + "M", -2.9);
        final Hap del20and40 = new Hap(snp(REF.substring(0, 21) + REF.substring(24), 37), "21M3D" + (n - 24) + "M", -2.9);
        // In cis: both SNPs on one haplotype, the other the reference.
        run("cis", one, 2, 1000, 1069, 30.0, 6, true, List.of(ref, snp20and40), reads("sA", new int[]{0, 1, 0, 1, 0, 1, 0, 1}, 8));
        // In trans: each SNP on its own haplotype.
        run("trans", one, 2, 1000, 1069, 30.0, 6, true, List.of(snp20, snp40), reads("sA", new int[]{0, 1, 0, 1, 0, 1, 0, 1}, 8));
        // Trans against the reference too: three haplotypes, the two SNPs apart.
        run("trans-with-ref", one, 2, 1000, 1069, 30.0, 6, true, List.of(ref, snp20, snp40), reads("sA", new int[]{1, 2, 1, 2, 1, 2, 1, 2}, 8));
        // A hom SNP with a het one: the hom call is on every called haplotype.
        run("hom-and-het", one, 2, 1000, 1069, 30.0, 6, true, List.of(snp20, snp20and40), reads("sA", new int[]{0, 1, 0, 1, 0, 1, 0, 1}, 8));
        // Three variants on one haplotype.
        run("three-cis", one, 2, 1000, 1069, 30.0, 6, true, List.of(ref, snp20and40and60), reads("sA", new int[]{0, 1, 0, 1, 0, 1, 0, 1}, 8));
        // Two cis, one trans.
        run("cis-and-trans", one, 2, 1000, 1069, 30.0, 6, true, List.of(snp20and40, snp60), reads("sA", new int[]{0, 1, 0, 1, 0, 1, 0, 1}, 8));
        // A triallelic site is never phased.
        run("triallelic-not-phased", one, 2, 1000, 1069, 30.0, 6, true, List.of(snp20, snp20b, snp20and40), reads("sA", new int[]{0, 1, 2, 0, 1, 2, 0, 1}, 8));
        // A deletion and a SNP in cis.
        run("deletion-cis", one, 2, 1000, 1069, 30.0, 6, true, List.of(ref, del20and40), reads("sA", new int[]{0, 1, 0, 1, 0, 1, 0, 1}, 8));
        // Two samples sharing the haplotypes.
        run("two-samples", two, 2, 1000, 1069, 30.0, 6, true, List.of(ref, snp20and40, snp40),
                concat(reads("sA", new int[]{0, 1, 0, 1, 0, 1}, 8), reads("sB", new int[]{0, 2, 0, 2, 0, 2}, 8)));
        // Unrelated haplotypes: neither the same nor complementary.
        run("unrelated", one, 2, 1000, 1069, 30.0, 6, true, List.of(ref, snp20, snp20and40, snp40), reads("sA", new int[]{1, 2, 3, 1, 2, 3, 0, 1}, 8));
        // A single call: nothing to phase.
        run("single", one, 2, 1000, 1069, 30.0, 6, true, List.of(ref, snp20), reads("sA", new int[]{0, 1, 0, 1, 0, 1}, 8));
        // Ploidy 3.
        run("ploidy-3", one, 3, 1000, 1069, 30.0, 6, true, List.of(ref, snp20and40), reads("sA", new int[]{0, 1, 1, 0, 1, 1, 0, 1, 1}, 8));
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
            final HaplotypeCallerGenotypingEngine engine = new HaplotypeCallerGenotypingEngine(hcArgs, new IndexedSampleList(samples), true, false);
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

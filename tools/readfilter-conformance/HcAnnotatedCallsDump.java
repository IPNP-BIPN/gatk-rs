/*
 * HaplotypeCaller's calls with its default annotations (StandardAnnotation and
 * StandardHCAnnotation): HaplotypeCallerGenotypingEngine.assignGenotypeLikelihoods with phasing
 * off, the annotations computed by `makeAnnotatedCall` from the read-by-allele matrix the call was
 * genotyped with, before the reverse trim.
 *
 * Reads are cut from the haplotypes and aligned to the reference through the haplotype's own
 * CIGAR, with varied base qualities, mapping qualities and strands, because the rank-sum, strand
 * and mapping-quality annotations read them. The random generator is reset before each case, as a
 * tool run resets it.
 *
 * Output:
 *
 *     case\t<label>\t<samples,...>\t<ploidy>
 *     haplotype\t<label>\t<index>\t<ref>\t<bases>\t<cigar>\t<score bits>
 *     read\t<label>\t<sample>\t<name>\t<start>\t<cigar>\t<bases>\t<quals phred+33>\t<mapq>\t<reverse>
 *     lk\t<label>\t<sample index>\t<haplotype index>\t<values bits,...>
 *     call\t<label>\t<VCF line>
 *     calls\t<label>\t<count>
 *     error\t<label>\t<exception class>: <message>
 *
 * The reference is REF at chr1:1000, the window chr1:1000-1069.
 *
 * Usage: HcAnnotatedCallsDump
 */

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.CalledHaplotypes;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
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

public class HcAnnotatedCallsDump {

    static final String REF = "GATTACACGTTGCAAGTCCGATAGCTTAGGCATCGATCGGATCCATGCAAGTCAGTAAAAAGTCCATGACTTGCA";
    static final int START = 1000;
    static final String[] ANNOTATIONS = {"BaseQualityRankSumTest", "ChromosomeCounts", "Coverage", "DepthPerAlleleBySample",
            "DepthPerSampleHC", "ExcessHet", "FisherStrand", "InbreedingCoeff", "MappingQualityRankSumTest", "QualByDepth",
            "RMSMappingQuality", "ReadPosRankSumTest", "StrandOddsRatio"};

    record Hap(String bases, String cigar, double score) {
    }

    /** A read: its sample, the haplotype it is cut from, where on the reference it starts. */
    record ReadSpec(String sample, int hap, int refStart, int length, int baseQual, int mapq, boolean reverse) {
    }

    public static void main(final String[] args) throws Exception {
        System.out.println("# HcAnnotatedCallsDump: assignGenotypeLikelihoods with HaplotypeCaller's default annotations");
        final int n = REF.length();
        final Hap ref = new Hap(REF, n + "M", -1.0);
        final Hap snp20 = new Hap(snp(REF, 20), n + "M", -2.0);
        final Hap snp20b = new Hap(snp(snp(REF, 20), 20), n + "M", -3.0);
        final Hap snp20and40 = new Hap(snp(snp(REF, 20), 40), n + "M", -2.5);
        final Hap del20 = new Hap(REF.substring(0, 21) + REF.substring(24), "21M3D" + (n - 24) + "M", -2.2);
        final Hap ins30 = new Hap(REF.substring(0, 31) + "GG" + REF.substring(31), "31M2I" + (n - 31) + "M", -2.6);
        final Hap homDel = new Hap(REF.substring(0, 56) + REF.substring(58), "56M2D" + (n - 58) + "M", -2.1);
        final Hap homDel1 = new Hap(REF.substring(0, 56) + REF.substring(57), "56M1D" + (n - 57) + "M", -2.3);

        run("het-snp", List.of("sA"), 2, List.of(ref, snp20), tile("sA", new int[]{0, 1}, 1000, 3));
        run("hom-snp", List.of("sA"), 2, List.of(ref, snp20), tile("sA", new int[]{1}, 1000, 3));
        run("strand-biased", List.of("sA"), 2, List.of(ref, snp20), biased("sA"));
        run("triallelic", List.of("sA"), 2, List.of(ref, snp20, snp20b), tile("sA", new int[]{1, 2}, 1000, 3));
        run("two-sites", List.of("sA"), 2, List.of(ref, snp20, snp20and40), tile("sA", new int[]{0, 2, 1}, 1000, 3));
        run("deletion", List.of("sA"), 2, List.of(ref, del20), tile("sA", new int[]{0, 1}, 1000, 3));
        run("insertion", List.of("sA"), 2, List.of(ref, ins30), tile("sA", new int[]{0, 1}, 1005, 3));
        run("homopolymer-trim", List.of("sA"), 2, List.of(ref, homDel, homDel1), tile("sA", new int[]{0, 2}, 1020, 2));
        run("two-samples", List.of("sA", "sB"), 2, List.of(ref, snp20, del20),
                concat(tile("sA", new int[]{0, 1}, 1000, 3), tile("sB", new int[]{2, 2, 0}, 1001, 4)));
        run("three-samples-hwe", List.of("sA", "sB", "sC"), 2, List.of(ref, snp20),
                concat(concat(tile("sA", new int[]{1}, 1000, 4), tile("sB", new int[]{1}, 1001, 4)), tile("sC", new int[]{0}, 1002, 4)));
        run("deep", List.of("sA"), 2, List.of(ref, snp20), tile("sA", new int[]{1, 1, 1, 0}, 995, 1));
        run("ploidy-1", List.of("sA"), 1, List.of(ref, snp20), tile("sA", new int[]{1}, 1000, 3));
    }

    static <T> List<T> concat(final List<T> a, final List<T> b) {
        final List<T> all = new ArrayList<>(a);
        all.addAll(b);
        return all;
    }

    static String snp(final String s, final int at) {
        final char c = s.charAt(at);
        return s.substring(0, at) + "ACGT".charAt(("ACGT".indexOf(c) + 1) % 4) + s.substring(at + 1);
    }

    /** Reads every `step` bases from `from` to 1030, cycling through `haps`, with varied qualities and strands. */
    static List<ReadSpec> tile(final String sample, final int[] haps, final int from, final int step) {
        final List<ReadSpec> out = new ArrayList<>();
        int i = 0;
        for (int start = from; start <= 1030; start += step) {
            out.add(new ReadSpec(sample, haps[i % haps.length], start, 36, 20 + (i * 7) % 21, 40 + (i * 13) % 21, i % 3 == 0));
            i++;
        }
        return out;
    }

    /** The alternate on the forward strand only, the reference on both. */
    static List<ReadSpec> biased(final String sample) {
        final List<ReadSpec> out = new ArrayList<>();
        int i = 0;
        for (int start = 1000; start <= 1015; start += 1) {
            final int hap = i % 2;
            out.add(new ReadSpec(sample, hap, start, 36, 30, 60, hap == 0 && i % 4 == 0));
            i++;
        }
        return out;
    }

    /**
     * The read of `length` haplotype bases starting at the haplotype base aligned to `refStart`,
     * with the CIGAR the haplotype's own alignment gives it.
     */
    static Object[] cut(final Hap hap, final int refStart, final int length) {
        final List<Integer> refOf = new ArrayList<>(); // the reference position of each haplotype base, -1 if inserted
        int r = START;
        for (final CigarElement e : TextCigarCodec.decode(hap.cigar()).getCigarElements()) {
            for (int k = 0; k < e.getLength(); k++) {
                if (e.getOperator() == CigarOperator.M) {
                    refOf.add(r++);
                } else if (e.getOperator() == CigarOperator.I) {
                    refOf.add(-1);
                } else if (e.getOperator() == CigarOperator.D) {
                    r++;
                }
            }
        }
        // The first haplotype base aligned at or after `refStart`: a start inside a deletion moves on.
        int first = 0;
        while (first < refOf.size() && (refOf.get(first) == -1 || refOf.get(first) < refStart)) {
            first++;
        }
        final int last = Math.min(first + length, refOf.size());
        final List<CigarElement> elements = new ArrayList<>();
        int prevRef = -2;
        for (int k = first; k < last; k++) {
            final int pos = refOf.get(k);
            if (pos == -1) {
                elements.add(new CigarElement(1, CigarOperator.I));
            } else {
                if (prevRef >= 0 && pos > prevRef + 1) {
                    elements.add(new CigarElement(pos - prevRef - 1, CigarOperator.D));
                }
                elements.add(new CigarElement(1, CigarOperator.M));
                prevRef = pos;
            }
        }
        final Cigar cigar = new org.broadinstitute.hellbender.utils.read.CigarBuilder().addAll(elements).make();
        return new Object[]{hap.bases().substring(first, last), cigar.toString(), refOf.get(first)};
    }

    static String bits(final double d) {
        return Long.toHexString(Double.doubleToRawLongBits(d));
    }

    static void run(final String label, final List<String> samples, final int ploidy, final List<Hap> haps, final List<ReadSpec> specs) throws Exception {
        System.out.println("case\t" + label + "\t" + String.join(",", samples) + "\t" + ploidy);
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
            final Map<String, List<ReadSpec>> bySample = new LinkedHashMap<>();
            for (final String s : samples) {
                evidence.put(s, new ArrayList<>());
                bySample.put(s, new ArrayList<>());
            }
            int count = 0;
            for (final ReadSpec spec : specs) {
                final Object[] cut = cut(haps.get(spec.hap()), spec.refStart(), spec.length());
                final String bases = (String) cut[0];
                final String cigar = (String) cut[1];
                final int readStart = (Integer) cut[2];
                final byte[] q = new byte[bases.length()];
                final StringBuilder qs = new StringBuilder();
                for (int k = 0; k < q.length; k++) {
                    q[k] = (byte) (spec.baseQual() + (k % 5));
                    qs.append((char) (q[k] + 33));
                }
                final String name = spec.sample() + "r" + count++;
                final GATKRead read = ArtificialReadUtils.createArtificialRead(header, name, 0, readStart,
                        bases.getBytes(StandardCharsets.US_ASCII), q, cigar);
                read.setReadGroup("rg" + spec.sample());
                read.setMappingQuality(spec.mapq());
                read.setIsReverseStrand(spec.reverse());
                evidence.get(spec.sample()).add(read);
                bySample.get(spec.sample()).add(spec);
                System.out.println("read\t" + label + "\t" + spec.sample() + "\t" + name + "\t" + readStart + "\t" + cigar + "\t" + bases
                        + "\t" + qs + "\t" + spec.mapq() + "\t" + spec.reverse());
            }
            final AlleleLikelihoods<GATKRead, Haplotype> likelihoods = new AlleleLikelihoods<>(new IndexedSampleList(samples),
                    new IndexedAlleleList<>(haplotypes), evidence);
            for (int s = 0; s < samples.size(); s++) {
                final LikelihoodMatrix<GATKRead, Haplotype> m = likelihoods.sampleMatrix(s);
                final List<ReadSpec> sampleSpecs = bySample.get(samples.get(s));
                for (int a = 0; a < haplotypes.size(); a++) {
                    final StringBuilder b = new StringBuilder();
                    for (int r = 0; r < sampleSpecs.size(); r++) {
                        final double v = sampleSpecs.get(r).hap() == a ? -0.5 - 0.125 * (r % 3) : -6.5 - 0.25 * ((a + r) % 4);
                        m.set(a, r, v);
                        b.append(r == 0 ? "" : ",").append(bits(v));
                    }
                    System.out.println("lk\t" + label + "\t" + s + "\t" + a + "\t" + b);
                }
            }
            final List<Annotation> annotations = new ArrayList<>();
            for (final String name : ANNOTATIONS) {
                annotations.add((Annotation) Class.forName("org.broadinstitute.hellbender.tools.walkers.annotator." + name)
                        .getDeclaredConstructor().newInstance());
            }
            final HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
            hcArgs.standardArgs.genotypeArgs.samplePloidy = ploidy;
            final HaplotypeCallerGenotypingEngine engine = new HaplotypeCallerGenotypingEngine(hcArgs, new IndexedSampleList(samples), false, false);
            engine.setAnnotationEngine(new VariantAnnotatorEngine(annotations, null, Collections.emptyList(), false, false));
            Utils.resetRandomGenerator();
            final CalledHaplotypes called = engine.assignGenotypeLikelihoods(haplotypes, likelihoods, new LinkedHashMap<>(),
                    REF.getBytes(StandardCharsets.US_ASCII), new SimpleInterval("chr1", START, START + REF.length() - 1),
                    new SimpleInterval("chr1", 1000, 1069), new FeatureContext(), Collections.emptyList(), false, 0, header,
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

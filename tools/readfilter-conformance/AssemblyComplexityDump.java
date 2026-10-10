/*
 * AssemblyComplexity: how a variant's haplotypes relate to the germline ones.
 *
 * `annotate(vc, haplotypeLikelihoods, germlineMode)` counts, for each haplotype, the reads that
 * best support it (ties broken), groups haplotypes by their events away from the site (HEC: each
 * group's support, descending), takes the most supported haplotype as germline (and the second if
 * it has at least half the first's support; the reference in germline mode), and per alternate
 * allele gives the edit distance (events in one but not the other, away from the site) from the
 * most supported haplotype carrying it to the closest germline one (HAPCOMP), and the share of the
 * carriers' support that haplotype holds (HAPDOM; 1 / haplotypes when they hold none). An event
 * map's allele matches a variant context allele with the context's longer reference trimmed.
 *
 * Output:
 *
 *     case\t<label>\t<germline mode>\t<vc start>\t<vc ref>\t<vc alts,...>
 *     haplotype\t<label>\t<index>\t<ref>\t<bases>\t<cigar>
 *     lk\t<label>\t<sample>\t<allele index>\t<values,...>      (evidence count = values)
 *     hec\t<label>\t<counts,...>
 *     hapcomp\t<label>\t<distances,...>
 *     hapdom\t<label>\t<bits:value,...>
 *     error\t<label>\t<exception class>: <message>
 *
 * The reference is REF at chr1:1000; event maps use maxMnpDistance 0, as HaplotypeCaller does.
 *
 * Usage: AssemblyComplexityDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.lang3.tuple.Triple;
import org.broadinstitute.hellbender.tools.walkers.annotator.AssemblyComplexity;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class AssemblyComplexityDump {

    static final String REF = "GATTACACGTTGCAAGTCCGATAGCTTAGGCATCGATCGGATCCATGCAAGTCAGT";
    static final int START = 1000;
    static final SAMFileHeader HEADER = ArtificialReadUtils.createArtificialSamHeader();

    record Hap(String bases, String cigar, boolean ref) {
    }

    public static void main(final String[] args) {
        System.out.println("# AssemblyComplexityDump: AssemblyComplexity.annotate");
        final Hap ref = new Hap(REF, REF.length() + "M", true);
        final Hap snp20 = new Hap(snp(REF, 20), REF.length() + "M", false);
        final Hap snp20b = new Hap(snp(snp(REF, 20), 20), REF.length() + "M", false);
        final Hap snp20and40 = new Hap(snp(snp(REF, 20), 40), REF.length() + "M", false);
        final Hap snp40 = new Hap(snp(REF, 40), REF.length() + "M", false);
        final Hap snp10and20 = new Hap(snp(snp(REF, 10), 20), REF.length() + "M", false);
        final Hap snp10and20and45 = new Hap(snp(snp(snp(REF, 10), 20), 45), REF.length() + "M", false);
        final Hap del20 = new Hap(REF.substring(0, 21) + REF.substring(24), "21M3D" + (REF.length() - 24) + "M", false);
        final Hap del20long = new Hap(REF.substring(0, 21) + REF.substring(26), "21M5D" + (REF.length() - 26) + "M", false);

        // A SNP at 1020 (offset 20), carried by one haplotype, the reference the germline.
        run("one-snp", false, 1020, List.of(alt(20)), List.of(ref, snp20), support(2, 0, 0, 0, 0, 0, 0, 1, 1));
        // The variant haplotype also carries a second event: HAPCOMP 1.
        run("snp-with-passenger", false, 1020, List.of(alt(20)), List.of(ref, snp20and40), support(2, 0, 0, 0, 0, 0, 0, 1, 1));
        // The variant haplotype is the most supported: it is germline itself, unless in reference mode.
        run("variant-dominant", false, 1020, List.of(alt(20)), List.of(ref, snp20and40), support(2, 1, 1, 1, 1, 1, 1, 0, 0));
        run("variant-dominant-germline-mode", true, 1020, List.of(alt(20)), List.of(ref, snp20and40), support(2, 1, 1, 1, 1, 1, 1, 0, 0));
        // A germline het elsewhere: the second haplotype is germline too, and the closest.
        run("germline-het", false, 1020, List.of(alt(20)), List.of(ref, snp40, snp20and40), support(3, 0, 0, 0, 0, 1, 1, 1, 2, 2));
        run("germline-het-weak-second", false, 1020, List.of(alt(20)), List.of(ref, snp40, snp20and40), support(3, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2));
        run("germline-het-exactly-half", false, 1020, List.of(alt(20)), List.of(ref, snp40, snp20and40), support(3, 0, 0, 0, 0, 0, 1, 1, 1, 2));
        run("germline-het-odd-half", false, 1020, List.of(alt(20)), List.of(ref, snp40, snp20and40), support(3, 0, 0, 0, 0, 0, 1, 1, 2));
        // Two variant haplotypes carrying the site, several events apart.
        run("far-variant", false, 1020, List.of(alt(20)), List.of(ref, snp10and20and45), support(2, 0, 0, 0, 0, 0, 1, 1));
        // Ties between haplotypes: the first index wins.
        run("ties", false, 1020, List.of(alt(20)), List.of(ref, snp20, snp10and20),
                lk(new int[][]{{0, 0, -5, -5}, {-5, -1, 0, 0}, {-5, -1, 0, -2}}));
        // Two carriers of the allele: dominance is the larger's share.
        run("two-carriers", false, 1020, List.of(alt(20)), List.of(ref, snp20, snp10and20, snp10and20and45),
                support(4, 0, 0, 0, 0, 0, 0, 1, 1, 2, 3));
        run("two-carriers-tied", false, 1020, List.of(alt(20)), List.of(ref, snp10and20, snp20, snp10and20and45),
                support(4, 0, 0, 0, 0, 0, 0, 1, 2, 3));
        // A carrier with no read support at all: 1 / haplotypes.
        run("unsupported-carrier", false, 1020, List.of(alt(20)), List.of(ref, snp20),
                lk(new int[][]{{0, 0, 0}, {-5, -5, -5}}));
        // Two alternate alleles at the site.
        run("biallelic-two-alts", false, 1020, List.of(alt(20), alt2(20)), List.of(ref, snp20, snp20b),
                lk(new int[][]{{0, 0, -5, -5, -5}, {-5, -5, 0, 0, -5}, {-5, -5, -5, -5, 0}}));
        // A deletion; the context's reference is longer than the event map's.
        run("deletion", false, 1020, List.of(Allele.create(REF.substring(20, 21))), List.of(ref, del20),
                lk(new int[][]{{0, 0, -5, -5}, {-5, -5, 0, 0}}), REF.substring(20, 24));
        run("deletion-merged-context", false, 1020, List.of(Allele.create(REF.substring(20, 21) + REF.substring(24, 26)), Allele.create(REF.substring(20, 21))),
                List.of(ref, del20, del20long), lk(new int[][]{{0, 0, -5, -5, -5}, {-5, -5, 0, 0, -5}, {-5, -5, -5, -5, 0}}), REF.substring(20, 26));
        // A spanning deletion allele and a symbolic one: zero.
        run("star-allele", false, 1020, List.of(alt(20), Allele.SPAN_DEL), List.of(ref, snp20),
                lk(new int[][]{{0, 0, -5}, {-5, -5, 0}}));
        run("non-ref-allele", false, 1020, List.of(alt(20), Allele.NON_REF_ALLELE), List.of(ref, snp20),
                lk(new int[][]{{0, 0, -5}, {-5, -5, 0}}));
        // An allele no haplotype carries: findFirst() on nothing.
        run("allele-not-carried", false, 1020, List.of(alt2(20)), List.of(ref, snp20),
                lk(new int[][]{{0, 0, -5}, {-5, -5, 0}}));
        // Two samples.
        run("two-samples", false, 1020, List.of(alt(20)), List.of(ref, snp20, snp20and40),
                lk(new int[][]{{0, 0, -5}, {-5, -5, 0}, {-5, -5, -5}}), lk(new int[][]{{0, -5}, {-5, -5}, {-5, 0}}));
    }

    static String snp(final String s, final int at) {
        final char c = s.charAt(at);
        return s.substring(0, at) + "ACGT".charAt(("ACGT".indexOf(c) + 1) % 4) + s.substring(at + 1);
    }

    static Allele alt(final int at) {
        return Allele.create(snp(REF, at).substring(at, at + 1));
    }

    static Allele alt2(final int at) {
        return Allele.create(snp(snp(REF, at), at).substring(at, at + 1));
    }

    /** One read per entry, best supporting that haplotype index (0) against the others (-5). */
    static int[][] support(final int haplotypes, final int... best) {
        final int[][] v = new int[haplotypes][best.length];
        for (int a = 0; a < haplotypes; a++) {
            for (int r = 0; r < best.length; r++) {
                v[a][r] = best[r] == a ? 0 : -5;
            }
        }
        return v;
    }

    static int[][] lk(final int[][] values) {
        return values;
    }

    static void run(final String label, final boolean germline, final int vcStart, final List<Allele> alts, final List<Hap> haps,
                    final int[][]... perSample) {
        run(label, germline, vcStart, alts, haps, perSample, REF.substring(vcStart - START, vcStart - START + 1));
    }

    static void run(final String label, final boolean germline, final int vcStart, final List<Allele> alts, final List<Hap> haps,
                    final int[][] values, final String refAllele) {
        run(label, germline, vcStart, alts, haps, new int[][][]{values}, refAllele);
    }

    static void run(final String label, final boolean germline, final int vcStart, final List<Allele> alts, final List<Hap> haps,
                    final int[][][] perSample, final String refAllele) {
        System.out.println("case\t" + label + "\t" + germline + "\t" + vcStart + "\t" + refAllele + "\t"
                + alts.stream().map(Allele::getDisplayString).collect(Collectors.joining(",")));
        try {
            final List<Haplotype> haplotypes = new ArrayList<>();
            for (int i = 0; i < haps.size(); i++) {
                final Hap h = haps.get(i);
                System.out.println("haplotype\t" + label + "\t" + i + "\t" + h.ref() + "\t" + h.bases() + "\t" + h.cigar());
                final Haplotype hap = new Haplotype(h.bases().getBytes(StandardCharsets.US_ASCII), h.ref());
                hap.setCigar(TextCigarCodec.decode(h.cigar()));
                hap.setAlignmentStartHapwrtRef(0);
                hap.setEventMap(EventMap.fromHaplotype(hap, REF.getBytes(StandardCharsets.US_ASCII),
                        new SimpleInterval("chr1", START, START + REF.length() - 1), 0));
                haplotypes.add(hap);
            }
            final List<String> samples = new ArrayList<>();
            final Map<String, List<GATKRead>> evidence = new LinkedHashMap<>();
            for (int s = 0; s < perSample.length; s++) {
                final String sample = "s" + (s + 1);
                samples.add(sample);
                final List<GATKRead> reads = new ArrayList<>();
                for (int r = 0; r < perSample[s][0].length; r++) {
                    reads.add(ArtificialReadUtils.createArtificialRead(HEADER, sample + "r" + r, 0, 1, 10));
                }
                evidence.put(sample, reads);
            }
            final AlleleLikelihoods<GATKRead, Haplotype> likelihoods = new AlleleLikelihoods<>(new IndexedSampleList(samples),
                    new IndexedAlleleList<>(haplotypes), evidence);
            for (int s = 0; s < perSample.length; s++) {
                final LikelihoodMatrix<GATKRead, Haplotype> m = likelihoods.sampleMatrix(s);
                for (int a = 0; a < perSample[s].length; a++) {
                    final StringBuilder b = new StringBuilder();
                    for (int r = 0; r < perSample[s][a].length; r++) {
                        m.set(a, r, perSample[s][a][r]);
                        b.append(r == 0 ? "" : ",").append(perSample[s][a][r]);
                    }
                    System.out.println("lk\t" + label + "\t" + s + "\t" + a + "\t" + b);
                }
            }
            final List<Allele> alleles = new ArrayList<>();
            alleles.add(Allele.create(refAllele, true));
            alleles.addAll(alts);
            final VariantContext vc = new VariantContextBuilder("dump", "chr1", vcStart, vcStart + refAllele.length() - 1, alleles).make();
            final Triple<int[], int[], double[]> result = AssemblyComplexity.annotate(vc, likelihoods, germline);
            System.out.println("hec\t" + label + "\t" + Arrays.stream(result.getLeft()).mapToObj(String::valueOf).collect(Collectors.joining(",")));
            System.out.println("hapcomp\t" + label + "\t" + Arrays.stream(result.getMiddle()).mapToObj(String::valueOf).collect(Collectors.joining(",")));
            System.out.println("hapdom\t" + label + "\t" + Arrays.stream(result.getRight())
                    .mapToObj(d -> Long.toHexString(Double.doubleToRawLongBits(d)) + ":" + d).collect(Collectors.joining(",")));
        } catch (final Exception e) {
            System.out.println("error\t" + label + "\t" + e.getClass().getSimpleName() + ": " + e.getMessage());
        }
    }
}

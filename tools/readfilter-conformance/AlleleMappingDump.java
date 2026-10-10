/*
 * From haplotypes to the alleles HaplotypeCaller genotypes at each event start.
 *
 * For every start in `EventMap.getEventStartPositions` inside the active window, as
 * `HaplotypeCallerGenotypingEngine.assignGenotypeLikelihoods` does:
 * `AssemblyBasedCallerUtils.getVariantsFromActiveHaplotypes` (the distinct events overlapping the
 * start, spanning ones included unless spanning-event genotyping is disabled, each a biallelic
 * context named after its haplotype), `replaceSpanDels` (an event starting earlier becomes ref/*),
 * `makeMergedVariantContext` (`GATKVariantContextUtils.simpleMerge`: the longest reference, every
 * alternate extended to it), `createAlleleMapper` (which haplotypes support each merged allele,
 * spanning events going to * or to the reference), and `AlleleLikelihoods.marginalize` (each
 * allele's likelihood the maximum over its haplotypes).
 *
 * Output:
 *
 *     case\t<label>\t<max mnp distance>\t<spanning genotyping>\t<window start>\t<window end>
 *     haplotype\t<label>\t<index>\t<ref>\t<bases>\t<cigar>
 *     lk\t<label>\t<haplotype index>\t<values,...>          (one sample, values as given)
 *     starts\t<label>\t<positions,...>
 *     event\t<label>\t<loc>\t<source>\t<start>\t<end>\t<ref>\t<alt>   (after replaceSpanDels)
 *     merged\t<label>\t<loc>\t<source>\t<start>\t<end>\t<alleles,...>
 *     mapper\t<label>\t<loc>\t<allele>\t<haplotype indices,...>
 *     marginal\t<label>\t<loc>\t<allele index>\t<allele>\t<bits,...>
 *     error\t<label>\t<loc>\t<exception class>: <message>
 *
 * Alleles print as `Allele.toString()` (`*` after a reference allele). The reference is REF at
 * chr1:1000.
 *
 * Usage: AlleleMappingDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingEngine;
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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.stream.Collectors;

public class AlleleMappingDump {

    static final String REF = "GATTACACGTTGCAAGTCCGATAGCTTAGGCATCGATCGGATCCATGCAAGTCAGTAAAAAGTC";
    static final int START = 1000;
    static final SAMFileHeader HEADER = ArtificialReadUtils.createArtificialSamHeader();

    public static void main(final String[] args) throws Exception {
        System.out.println("# AlleleMappingDump: getVariantsFromActiveHaplotypes, makeMergedVariantContext, createAlleleMapper, marginalize");
        final int n = REF.length();
        final String[] ref = {REF, n + "M"};
        final String[] snp20 = {snp(REF, 20), n + "M"};
        final String[] snp20b = {snp(snp(REF, 20), 20), n + "M"};
        final String[] snp30 = {snp(REF, 30), n + "M"};
        final String[] snp20and30 = {snp(snp(REF, 20), 30), n + "M"};
        final String[] del20 = {REF.substring(0, 21) + REF.substring(24), "21M3D" + (n - 24) + "M"};
        final String[] del20long = {REF.substring(0, 21) + REF.substring(26), "21M5D" + (n - 26) + "M"};
        final String[] ins20 = {REF.substring(0, 21) + "TT" + REF.substring(21), "21M2I" + (n - 21) + "M"};
        final String[] del18 = {REF.substring(0, 19) + REF.substring(23), "19M4D" + (n - 23) + "M"};
        final String[] mnp20 = {snp(snp(REF, 20), 21), n + "M"};
        final String[] snpDel20 = {snp(REF, 20).substring(0, 21) + REF.substring(24), "21M3D" + (n - 24) + "M"};
        final String[] homDel = {REF.substring(0, 56) + REF.substring(58), "56M2D" + (n - 58) + "M"};
        final String[] homDel1 = {REF.substring(0, 56) + REF.substring(57), "56M1D" + (n - 57) + "M"};

        run("snp", 0, true, 1000, 1063, List.of(ref, snp20));
        run("triallelic-snp", 0, true, 1000, 1063, List.of(ref, snp20, snp20b));
        run("two-sites", 0, true, 1000, 1063, List.of(ref, snp20, snp30, snp20and30));
        run("deletion-and-snp", 0, true, 1000, 1063, List.of(ref, snp20, del20));
        run("two-deletions", 0, true, 1000, 1063, List.of(ref, del20, del20long));
        run("insertion-and-deletion", 0, true, 1000, 1063, List.of(ref, ins20, del20));
        run("spanning-deletion", 0, true, 1000, 1063, List.of(ref, snp20, del18));
        run("spanning-deletion-off", 0, false, 1000, 1063, List.of(ref, snp20, del18));
        run("mnp-split", 0, true, 1000, 1063, List.of(ref, mnp20, snp20));
        run("mnp-merged", 1, true, 1000, 1063, List.of(ref, mnp20, snp20));
        run("snp-plus-deletion-compound", 0, true, 1000, 1063, List.of(ref, snpDel20, del20, snp20));
        run("homopolymer-deletions", 0, true, 1000, 1063, List.of(ref, homDel, homDel1));
        run("window-excludes", 0, true, 1025, 1063, List.of(ref, snp20, snp30, del18));
        run("no-reference-haplotype", 0, true, 1000, 1063, List.of(snp20, snp30));
    }

    static String snp(final String s, final int at) {
        final char c = s.charAt(at);
        return s.substring(0, at) + "ACGT".charAt(("ACGT".indexOf(c) + 1) % 4) + s.substring(at + 1);
    }

    static String bits(final double d) {
        return Long.toHexString(Double.doubleToRawLongBits(d));
    }

    static void run(final String label, final int maxMnp, final boolean spanning, final int windowStart, final int windowEnd,
                    final List<String[]> haps) throws Exception {
        System.out.println("case\t" + label + "\t" + maxMnp + "\t" + spanning + "\t" + windowStart + "\t" + windowEnd);
        final List<Haplotype> haplotypes = new ArrayList<>();
        for (int i = 0; i < haps.size(); i++) {
            final boolean isRef = haps.get(i)[0].equals(REF);
            System.out.println("haplotype\t" + label + "\t" + i + "\t" + isRef + "\t" + haps.get(i)[0] + "\t" + haps.get(i)[1]);
            final Haplotype h = new Haplotype(haps.get(i)[0].getBytes(StandardCharsets.US_ASCII), isRef);
            h.setCigar(TextCigarCodec.decode(haps.get(i)[1]));
            h.setAlignmentStartHapwrtRef(0);
            h.setGenomeLocation(new SimpleInterval("chr1", START, START + REF.length() - 1));
            haplotypes.add(h);
        }
        final SimpleInterval refLoc = new SimpleInterval("chr1", START, START + REF.length() - 1);
        final byte[] refBytes = REF.getBytes(StandardCharsets.US_ASCII);
        EventMap.buildEventMapsForHaplotypes(haplotypes, refBytes, refLoc, false, maxMnp);

        // One sample, five reads, values in [-6, 0) varying by haplotype and read.
        final List<GATKRead> reads = new ArrayList<>();
        for (int r = 0; r < 5; r++) {
            reads.add(ArtificialReadUtils.createArtificialRead(HEADER, "r" + r, 0, 1, 10));
        }
        final AlleleLikelihoods<GATKRead, Haplotype> likelihoods = new AlleleLikelihoods<>(new IndexedSampleList(List.of("s1")),
                new IndexedAlleleList<>(haplotypes), Map.of("s1", reads));
        final LikelihoodMatrix<GATKRead, Haplotype> m = likelihoods.sampleMatrix(0);
        for (int a = 0; a < haplotypes.size(); a++) {
            final StringBuilder b = new StringBuilder();
            for (int r = 0; r < reads.size(); r++) {
                final double v = -0.25 - ((a * 7 + r * 3) % 11) * 0.5;
                m.set(a, r, v);
                b.append(r == 0 ? "" : ",").append(v);
            }
            System.out.println("lk\t" + label + "\t" + a + "\t" + b);
        }

        final SortedSet<Integer> starts = EventMap.getEventStartPositions(haplotypes);
        System.out.println("starts\t" + label + "\t" + starts.stream().map(String::valueOf).collect(Collectors.joining(",")));
        for (final int loc : starts) {
            if (loc < windowStart || loc > windowEnd) {
                continue;
            }
            try {
                final List<VariantContext> events = AssemblyBasedCallerUtils.getVariantsFromActiveHaplotypes(loc, haplotypes, spanning);
                final List<VariantContext> replaced = HaplotypeCallerGenotypingEngine.replaceSpanDels(events,
                        Allele.create(refBytes[loc - START], true), loc);
                for (final VariantContext vc : replaced) {
                    System.out.println("event\t" + label + "\t" + loc + "\t" + vc.getSource() + "\t" + vc.getStart() + "\t" + vc.getEnd()
                            + "\t" + vc.getReference() + "\t" + vc.getAlternateAllele(0));
                }
                final VariantContext merged = AssemblyBasedCallerUtils.makeMergedVariantContext(replaced);
                if (merged == null) {
                    System.out.println("merged\t" + label + "\t" + loc + "\tnull");
                    continue;
                }
                System.out.println("merged\t" + label + "\t" + loc + "\t" + merged.getSource() + "\t" + merged.getStart() + "\t" + merged.getEnd()
                        + "\t" + merged.getAlleles().stream().map(Allele::toString).collect(Collectors.joining(",")));
                final Map<Allele, List<Haplotype>> mapper = AssemblyBasedCallerUtils.createAlleleMapper(merged, loc, haplotypes, spanning);
                for (final Map.Entry<Allele, List<Haplotype>> e : mapper.entrySet()) {
                    System.out.println("mapper\t" + label + "\t" + loc + "\t" + e.getKey() + "\t"
                            + e.getValue().stream().map(h -> String.valueOf(haplotypes.indexOf(h))).collect(Collectors.joining(",")));
                }
                final AlleleLikelihoods<GATKRead, Allele> marginal = likelihoods.marginalize(mapper);
                for (int a = 0; a < marginal.numberOfAlleles(); a++) {
                    final StringBuilder b = new StringBuilder();
                    for (int r = 0; r < reads.size(); r++) {
                        b.append(r == 0 ? "" : ",").append(bits(marginal.sampleMatrix(0).get(a, r)));
                    }
                    System.out.println("marginal\t" + label + "\t" + loc + "\t" + a + "\t" + marginal.getAllele(a) + "\t" + b);
                }
            } catch (final Exception e) {
                System.out.println("error\t" + label + "\t" + loc + "\t" + e.getClass().getSimpleName() + ": " + e.getMessage());
            }
        }
    }
}

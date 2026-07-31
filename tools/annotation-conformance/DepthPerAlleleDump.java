/*
 * DepthPerAlleleBySample, AlleleFraction and DepthPerSampleHC, taken from the reference.
 *
 * AD, AF and the HaplotypeCaller's DP, all counted off a MARGINALISED likelihood matrix.
 *
 *   - the allele order the marginalisation uses is a HashMap's. Collectors.toMap builds one and
 *     marginalize takes its key set as the new allele array, so the new matrix's allele order
 *     follows Allele.hashCode. That order is observable, because searchBestAllele breaks a tie by
 *     keeping the first index;
 *   - marginalize takes the MAXIMUM likelihood of the old alleles a new one stands for, not their
 *     sum, so collapsing two alleles a read supports equally leaves the value unchanged;
 *   - AD is keyed on vc.getAlleles() and counted over the matrix, so an allele the matrix holds
 *     and the variant does not is counted into a bucket nobody reads;
 *   - AF drops the first entry, so it is one shorter than AD and the two do not line up;
 *   - DP counts INFORMATIVE reads only, so it can be lower than the INFO-level DP over the same
 *     site.
 *
 * Output:
 *
 *     ad\t<label>\t<counts, comma-separated>
 *     af\t<label>\t<fractions, bits, comma-separated>
 *     dp\t<label>\t<count>
 *     order\t<alleles in, comma-separated>\t<alleles out>
 *     marginal\t<label>\t<sample>\t<allele>\t<likelihood bits, comma-separated>
 *
 * Usage: DepthPerAlleleDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import org.broadinstitute.hellbender.tools.walkers.annotator.AlleleFraction;
import org.broadinstitute.hellbender.tools.walkers.annotator.DepthPerAlleleBySample;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class DepthPerAlleleDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele ALT = Allele.create("C", false);
    static final Allele ALT2 = Allele.create("G", false);
    static final Allele ALT3 = Allele.create("T", false);
    static final SAMFileHeader HEADER = makeHeader();
    static final int VARIANT_START = 105;

    public static void main(final String[] args) {
        System.out.println("# DepthPerAlleleDump: AD, AF, DP and the marginalisation under them");

        // The HashMap order, on every allele set the rest of the dump uses.
        order(List.of(REF, ALT));
        order(List.of(REF, ALT, ALT2));
        order(List.of(REF, ALT, ALT2, ALT3));
        order(List.of(Allele.create("AT", true), Allele.create("A", false)));
        order(List.of(Allele.create("ACGTACGTAC", true), Allele.create("A", false),
                Allele.create("ACGTACGTACGT", false)));

        emit("two-and-two", biallelic(2, 2, false));
        emit("all-ref", biallelic(4, 0, false));
        emit("all-alt", biallelic(0, 4, false));
        emit("empty", biallelic(0, 0, false));
        emit("uninformative", biallelic(2, 2, true));
        emit("triallelic", triallelic());
        emit("matrix-has-extra-allele", extraAllele());
    }

    static SAMFileHeader makeHeader() {
        final SAMSequenceDictionary dictionary =
                new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", 1000)));
        final SAMFileHeader header = new SAMFileHeader(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return header;
    }

    static void order(final List<Allele> alleles) {
        final Map<Allele, List<Allele>> map = new LinkedHashSet<>(alleles).stream()
                .collect(Collectors.toMap(a -> a, Arrays::asList));
        final StringBuilder in = new StringBuilder();
        for (final Allele a : alleles) {
            if (in.length() > 0) { in.append(','); }
            in.append(a.getDisplayString()).append(a.isReference() ? "*" : "");
        }
        final StringBuilder out = new StringBuilder();
        for (final Allele a : map.keySet()) {
            if (out.length() > 0) { out.append(','); }
            out.append(a.getDisplayString()).append(a.isReference() ? "*" : "");
        }
        System.out.printf("order\t%s\t%s%n", in, out);
    }

    static GATKRead read(final String name, final int quality) {
        final SAMRecord record = new SAMRecord(HEADER);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(100);
        record.setCigarString("20M");
        final byte[] bases = new byte[20];
        Arrays.fill(bases, (byte) 'A');
        record.setReadBases(bases);
        final byte[] qualities = new byte[20];
        Arrays.fill(qualities, (byte) quality);
        record.setBaseQualities(qualities);
        record.setMappingQuality(60);
        return new SAMRecordToGATKReadAdapter(record);
    }

    static AlleleLikelihoods<GATKRead, Allele> biallelic(final int refReads, final int altReads,
                                                          final boolean uninformative) {
        final List<GATKRead> reads = new ArrayList<>();
        for (int i = 0; i < refReads; i++) { reads.add(read("r" + i, 30)); }
        for (int i = 0; i < altReads; i++) { reads.add(read("a" + i, 30)); }
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", reads);
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(List.of("s1")), new IndexedAlleleList<>(REF, ALT), bySample);
        final LikelihoodMatrix<GATKRead, Allele> m = likelihoods.sampleMatrix(0);
        for (int e = 0; e < reads.size(); e++) {
            final boolean isRef = e < refReads;
            if (uninformative) {
                // Within the informative threshold of each other, so nothing is counted.
                m.set(0, e, -1.0);
                m.set(1, e, -1.1);
            } else {
                m.set(0, e, isRef ? -1 : -10);
                m.set(1, e, isRef ? -10 : -1);
            }
        }
        return likelihoods;
    }

    static AlleleLikelihoods<GATKRead, Allele> triallelic() {
        final List<GATKRead> reads = List.of(read("r0", 30), read("a0", 30), read("b0", 30),
                read("t0", 30));
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", new ArrayList<>(reads));
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(List.of("s1")),
                new IndexedAlleleList<>(List.of(REF, ALT, ALT2)), bySample);
        final LikelihoodMatrix<GATKRead, Allele> m = likelihoods.sampleMatrix(0);
        for (int a = 0; a < 3; a++) {
            for (int e = 0; e < reads.size(); e++) {
                // The last read ties between two alleles, which is where the HashMap order shows.
                m.set(a, e, e == 3 ? (a == 1 || a == 2 ? -1 : -10) : (a == e ? -1 : -10));
            }
        }
        return likelihoods;
    }

    static AlleleLikelihoods<GATKRead, Allele> extraAllele() {
        final List<GATKRead> reads = List.of(read("r0", 30), read("a0", 30), read("x0", 30));
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", new ArrayList<>(reads));
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(List.of("s1")),
                new IndexedAlleleList<>(List.of(REF, ALT, ALT3)), bySample);
        final LikelihoodMatrix<GATKRead, Allele> m = likelihoods.sampleMatrix(0);
        for (int a = 0; a < 3; a++) {
            for (int e = 0; e < reads.size(); e++) {
                m.set(a, e, a == e ? -1 : -10);
            }
        }
        return likelihoods;
    }

    static void emit(final String label, final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        final List<Allele> alleles = likelihoods.alleles().size() > 2
                && likelihoods.alleles().contains(ALT2)
                ? List.of(REF, ALT, ALT2) : List.of(REF, ALT);
        final VariantContext vc = new VariantContextBuilder().chr("chr1").start(VARIANT_START)
                .stop(VARIANT_START).alleles(alleles)
                .genotypes(List.of(new GenotypeBuilder("s1", List.of(REF, ALT)).make())).make();
        final Genotype g = vc.getGenotype(0);

        try {
            final int[] ad = DepthPerAlleleBySample.annotateWithLikelihoods(vc, g,
                    new LinkedHashSet<>(vc.getAlleles()), likelihoods);
            final StringBuilder counts = new StringBuilder();
            for (final int count : ad) {
                if (counts.length() > 0) { counts.append(','); }
                counts.append(count);
            }
            System.out.printf("ad\t%s\t%s%n", label, counts);
        } catch (final Exception | AssertionError e) {
            System.out.printf("ad\t%s\tE:%s%n", label, e.getClass().getName());
        }

        try {
            final GenotypeBuilder gb = new GenotypeBuilder("s1", List.of(REF, ALT));
            new AlleleFraction().annotate(null, vc, g, gb, likelihoods);
            final Object value = gb.make().getAnyAttribute("AF");
            final StringBuilder fractions = new StringBuilder();
            if (value instanceof double[]) {
                for (final double f : (double[]) value) {
                    if (fractions.length() > 0) { fractions.append(','); }
                    fractions.append(Double.doubleToRawLongBits(f));
                }
            }
            System.out.printf("af\t%s\t%s%n", label, fractions);
        } catch (final Exception | AssertionError e) {
            System.out.printf("af\t%s\tE:%s%n", label, e.getClass().getName());
        }

        // DP as DepthPerSampleHC computes it: informative reads only.
        final long informative = likelihoods.bestAllelesBreakingTies("s1").stream()
                .filter(ba -> ba.isInformative()).count();
        System.out.printf("dp\t%s\t%d%n", label, informative);

        // The marginalised matrix itself, so the maximum rule is measured rather than inferred.
        final Map<Allele, List<Allele>> subset = new LinkedHashSet<>(vc.getAlleles()).stream()
                .collect(Collectors.toMap(a -> a, Arrays::asList));
        final AlleleLikelihoods<GATKRead, Allele> marginal = likelihoods.marginalize(subset);
        for (int a = 0; a < marginal.numberOfAlleles(); a++) {
            final StringBuilder values = new StringBuilder();
            final LikelihoodMatrix<GATKRead, Allele> m = marginal.sampleMatrix(0);
            for (int e = 0; e < m.evidenceCount(); e++) {
                if (values.length() > 0) { values.append(','); }
                values.append(Double.doubleToRawLongBits(m.get(a, e)));
            }
            System.out.printf("marginal\t%s\ts1\t%s\t%s%n", label,
                    marginal.getAllele(a).getDisplayString()
                            + (marginal.getAllele(a).isReference() ? "*" : ""), values);
        }
    }
}

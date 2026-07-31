/*
 * QualByDepth, GenotypeSummaries and LikelihoodRankSumTest, taken from the reference.
 *
 *   - QD above 35 is RANDOMISED: the reference replaces it with 30 plus a Gaussian jitter, so the
 *     value written to a VCF is a draw from a seeded generator. This dump records where the
 *     boundary is and what the raw ratio was, so the port can refuse the branch and still be
 *     measured up to it;
 *   - the depth QD divides by is not DP. Only het and hom-var genotypes count, the whole AD total
 *     is added for each, and a separate AD-restricted tally collects the same totals only where
 *     the alternate part exceeds one. If that tally ends non-zero it REPLACES the depth, so one
 *     sample with two alternate reads can discard every other sample's depth;
 *   - GQ_MEAN and GQ_STDDEV are strings, and the deviation is written only with more than one GQ;
 *   - NCC counts no-call ALLELES across genotypes, not no-call samples;
 *   - LikelihoodRankSumTest is the one member of the rank-sum family whose value comes from the
 *     matrix rather than from the read.
 *
 * Output:
 *
 *     qd\t<label>\t<value or E:class>
 *     depth\t<label>\t<depth>
 *     summaries\t<label>\t<key>=<value>[<class>];...
 *     lrs\t<label>\t<key>=<value>[<class>]
 *
 * Usage: SiteStatisticsDump
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

import org.broadinstitute.hellbender.tools.walkers.annotator.GenotypeSummaries;
import org.broadinstitute.hellbender.tools.walkers.annotator.LikelihoodRankSumTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.QualByDepth;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

public class SiteStatisticsDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele ALT = Allele.create("C", false);
    static final SAMFileHeader HEADER = makeHeader();
    static final int START = 105;

    public static void main(final String[] args) {
        System.out.println("# SiteStatisticsDump: QD, NCC/GQ_MEAN/GQ_STDDEV and LikelihoodRankSum");

        // QD, across the threshold and around the depth rules.
        qd("ordinary", -25.0, new int[][] {{5, 5}}, null, null);
        qd("just-below-threshold", -34.9, new int[][] {{5, 5}}, null, null);
        qd("at-threshold", -35.0, new int[][] {{5, 5}}, null, null);
        qd("above-threshold", -100.0, new int[][] {{5, 5}}, null, null);
        qd("no-qual", null, new int[][] {{5, 5}}, null, null);
        qd("raw-qual-approx", null, new int[][] {{5, 5}}, 300, null);
        qd("hom-ref-only", -25.0, new int[][] {{10, 0}}, null, "0/0");
        qd("one-alt-read", -20.0, new int[][] {{9, 1}}, null, null);
        qd("two-alt-reads", -20.0, new int[][] {{8, 2}}, null, null);
        qd("zero-ad", -25.0, new int[][] {{0, 0}}, null, null);
        qd("no-ad-with-dp", -25.0, null, null, null);
        qd("two-samples", -20.0, new int[][] {{9, 1}, {8, 2}}, null, null);

        // GenotypeSummaries.
        summaries("no-genotypes", new int[0], 0);
        summaries("one-gq", new int[] {50}, 0);
        summaries("two-gqs", new int[] {50, 70}, 0);
        summaries("three-gqs", new int[] {10, 50, 99}, 0);
        summaries("no-gq", new int[0], 2);
        summaries("with-no-calls", new int[] {50, 70}, 2);

        // LikelihoodRankSum, whose values are the likelihoods themselves.
        lrs("separated");
        lrs("overlapping");
        lrs("ref-only");
    }

    static SAMFileHeader makeHeader() {
        final SAMSequenceDictionary dictionary =
                new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", 1000)));
        final SAMFileHeader header = new SAMFileHeader(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return header;
    }

    static void qd(final String label, final Double log10PError, final int[][] ads,
                   final Integer rawQualApprox, final String genotypeString) {
        final VariantContextBuilder builder = new VariantContextBuilder().chr("chr1").start(START)
                .stop(START).alleles(List.of(REF, ALT));
        if (log10PError != null) {
            builder.log10PError(log10PError);
        }
        if (rawQualApprox != null) {
            builder.attribute("QUALapprox", rawQualApprox);
        }
        final List<Genotype> genotypes = new ArrayList<>();
        final int sampleCount = ads == null ? 1 : ads.length;
        for (int s = 0; s < sampleCount; s++) {
            final List<Allele> called = "0/0".equals(genotypeString)
                    ? List.of(REF, REF) : List.of(REF, ALT);
            final GenotypeBuilder gb = new GenotypeBuilder("s" + s, called);
            if (ads != null) {
                gb.AD(ads[s]);
            } else {
                gb.DP(17);
            }
            genotypes.add(gb.make());
        }
        final VariantContext vc = builder.genotypes(genotypes).make();

        System.out.printf("depth\t%s\t%d%n", label,
                QualByDepth.getDepth(vc.getGenotypes(), null));
        try {
            final Map<String, Object> result = new QualByDepth().annotate(null, vc, null);
            final StringJoiner joiner = new StringJoiner(";");
            for (final Map.Entry<String, Object> entry : result.entrySet()) {
                joiner.add(String.format("%s=%s", entry.getKey(), entry.getValue()));
            }
            System.out.printf("qd\t%s\t%s%n", label, joiner);
        } catch (final Exception | AssertionError e) {
            System.out.printf("qd\t%s\tE:%s%n", label, e.getClass().getName());
        }
    }

    static void summaries(final String label, final int[] gqs, final int noCallSamples) {
        final List<Genotype> genotypes = new ArrayList<>();
        for (int i = 0; i < gqs.length; i++) {
            genotypes.add(new GenotypeBuilder("s" + i, List.of(REF, ALT)).GQ(gqs[i]).make());
        }
        for (int i = 0; i < noCallSamples; i++) {
            genotypes.add(new GenotypeBuilder("n" + i,
                    List.of(Allele.NO_CALL, Allele.NO_CALL)).make());
        }
        final VariantContext vc = new VariantContextBuilder().chr("chr1").start(START).stop(START)
                .alleles(List.of(REF, ALT)).genotypes(genotypes).make();
        final Map<String, Object> result = new GenotypeSummaries().annotate(null, vc, null);
        final StringJoiner joiner = new StringJoiner(";");
        for (final Map.Entry<String, Object> entry : result.entrySet()) {
            joiner.add(String.format("%s=%s[%s]", entry.getKey(), entry.getValue(),
                    entry.getValue().getClass().getName()));
        }
        System.out.printf("summaries\t%s\t%s%n", label, joiner);
    }

    static void lrs(final String label) {
        final List<GATKRead> reads = new ArrayList<>();
        for (int i = 0; i < 12; i++) {
            reads.add(read("r" + i));
        }
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", reads);
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(List.of("s1")), new IndexedAlleleList<>(REF, ALT), bySample);
        final LikelihoodMatrix<GATKRead, Allele> m = likelihoods.sampleMatrix(0);
        for (int e = 0; e < reads.size(); e++) {
            final boolean isRef = "ref-only".equals(label) || e < 6;
            final double strong = "overlapping".equals(label) ? -1 - (e * 0.1) : -1;
            final double weak = "overlapping".equals(label) ? -5 - (e * 0.1) : -10;
            m.set(0, e, isRef ? strong : weak);
            m.set(1, e, isRef ? weak : strong);
        }
        final VariantContext vc = new VariantContextBuilder().chr("chr1").start(START).stop(START)
                .alleles(List.of(REF, ALT))
                .genotypes(List.of(new GenotypeBuilder("s1", List.of(REF, ALT)).make())).make();
        final Map<String, Object> result = new LikelihoodRankSumTest().annotate(null, vc, likelihoods);
        final StringJoiner joiner = new StringJoiner(";");
        for (final Map.Entry<String, Object> entry : result.entrySet()) {
            joiner.add(String.format("%s=%s[%s]", entry.getKey(), entry.getValue(),
                    entry.getValue().getClass().getName()));
        }
        System.out.printf("lrs\t%s\t%s%n", label, joiner);
    }

    static GATKRead read(final String name) {
        final SAMRecord record = new SAMRecord(HEADER);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(100);
        record.setCigarString("20M");
        final byte[] bases = new byte[20];
        Arrays.fill(bases, (byte) 'A');
        record.setReadBases(bases);
        final byte[] qualities = new byte[20];
        Arrays.fill(qualities, (byte) 30);
        record.setBaseQualities(qualities);
        record.setMappingQuality(60);
        return new SAMRecordToGATKReadAdapter(record);
    }
}

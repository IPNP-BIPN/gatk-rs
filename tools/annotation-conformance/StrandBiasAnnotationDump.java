/*
 * FisherStrand, StrandOddsRatio and StrandBiasBySample, taken from the reference.
 *
 * FS, SOR and SB: the same 2x2 contingency table of forward and reverse reads on the reference and
 * the alternate, read three ways. Four things decide what a consumer sees:
 *
 *   - the table is built PER SAMPLE and only then added up, and a sample under the threshold
 *     contributes nothing rather than contributing to a pooled total that would pass. FS uses
 *     MIN_COUNT = 2 and SOR uses 0, so they do not see the same table;
 *   - the GENOTYPE FIELD WINS: if any genotype carries SB, both are computed from those arrays and
 *     the likelihood matrix is never consulted, which makes SB load-bearing rather than
 *     diagnostic;
 *   - FS normalises above 400 reads by scaling to about 200 and TRUNCATING to int, so a deep site
 *     is scored from a table that no longer sums to its coverage. SOR does not normalise;
 *   - FS is phred-scaled off a floor of 1e-320 "to prevent INFINITYs", so the largest value a site
 *     can report is about 3200 whatever the evidence.
 *
 * Output:
 *
 *     anno\t<annotation>\t<label>\t<key>=<value>[<class>];...
 *     table\t<label>\t<minCount>\t<a,b,c,d>
 *     fisher\t<a,b,c,d>\t<p-value bits>
 *     sor\t<a,b,c,d>\t<sor bits>
 *
 * Usage: StrandBiasAnnotationDump
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

import org.broadinstitute.hellbender.tools.walkers.annotator.FisherStrand;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.StrandBiasBySample;
import org.broadinstitute.hellbender.tools.walkers.annotator.StrandBiasTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.StrandOddsRatio;
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

public class StrandBiasAnnotationDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele ALT = Allele.create("C", false);
    static final SAMFileHeader HEADER = makeHeader();
    static final int VARIANT_START = 105;

    public static void main(final String[] args) {
        System.out.println("# StrandBiasAnnotationDump: FS, SOR and SB");

        emit("null-likelihoods", site(null), null);
        emit("monomorphic", monomorphic(), balanced(5, 5, 5, 5));
        emit("balanced", site(null), balanced(5, 5, 5, 5));
        emit("skewed", site(null), balanced(10, 0, 0, 10));
        emit("one-read", site(null), balanced(1, 0, 0, 0));
        emit("two-reads", site(null), balanced(1, 1, 0, 0));
        emit("three-reads", site(null), balanced(2, 1, 0, 0));
        emit("deep-balanced", site(null), balanced(200, 200, 200, 200));
        emit("deep-skewed", site(null), balanced(400, 1, 1, 400));
        emit("ref-only", site(null), balanced(5, 5, 0, 0));
        emit("alt-only", site(null), balanced(0, 0, 5, 5));
        emit("empty-matrix", site(null), balanced(0, 0, 0, 0));

        // The genotype field wins: the matrix says one thing and SB says another.
        emit("sb-field-string", site("1,2,3,4"), balanced(50, 50, 50, 50));
        emit("sb-field-list", siteWithList(), balanced(50, 50, 50, 50));
        emit("sb-field-below-threshold", site("0,1,0,0"), balanced(50, 50, 50, 50));

        // The table on its own, at both thresholds, so the per-sample rule is visible.
        for (final int minCount : new int[] {0, 2}) {
            table("balanced", minCount, balanced(5, 5, 5, 5));
            table("one-read", minCount, balanced(1, 0, 0, 0));
            table("two-reads", minCount, balanced(1, 1, 0, 0));
            table("three-reads", minCount, balanced(2, 1, 0, 0));
            table("two-samples-each-small", minCount, twoSamples());
        }

        // The two statistics on tables directly, including the ones the normalisation reshapes.
        for (final int[] cells : new int[][] {
                {0, 0, 0, 0}, {1, 0, 0, 0}, {5, 5, 5, 5}, {10, 0, 0, 10}, {0, 10, 10, 0},
                {100, 100, 100, 100}, {200, 200, 200, 200}, {201, 200, 200, 200},
                {1000, 1000, 1000, 1000}, {1000, 1, 1, 1000}, {3, 1, 1, 3}, {50, 30, 20, 40},
                {1, 1, 1, 1}, {2, 0, 0, 2}, {7, 0, 0, 0}, {0, 0, 7, 0}}) {
            final int[][] t = {{cells[0], cells[1]}, {cells[2], cells[3]}};
            System.out.printf("fisher\t%d,%d,%d,%d\t%d%n", cells[0], cells[1], cells[2], cells[3],
                    Double.doubleToRawLongBits(FisherStrand.pValueForContingencyTable(t)));
            System.out.printf("sor\t%d,%d,%d,%d\t%d%n", cells[0], cells[1], cells[2], cells[3],
                    Double.doubleToRawLongBits(StrandOddsRatio.calculateSOR(t)));
        }
    }

    static SAMFileHeader makeHeader() {
        final SAMSequenceDictionary dictionary =
                new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", 1000)));
        final SAMFileHeader header = new SAMFileHeader(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return header;
    }

    static VariantContext site(final String sbField) {
        final GenotypeBuilder builder = new GenotypeBuilder("s1", List.of(REF, ALT));
        if (sbField != null) {
            builder.attribute("SB", sbField);
        }
        return new VariantContextBuilder().chr("chr1").start(VARIANT_START).stop(VARIANT_START)
                .alleles(List.of(REF, ALT)).genotypes(List.of(builder.make())).make();
    }

    static VariantContext siteWithList() {
        final Genotype genotype = new GenotypeBuilder("s1", List.of(REF, ALT))
                .attribute("SB", new ArrayList<>(List.of(4, 3, 2, 1))).make();
        return new VariantContextBuilder().chr("chr1").start(VARIANT_START).stop(VARIANT_START)
                .alleles(List.of(REF, ALT)).genotypes(List.of(genotype)).make();
    }

    static VariantContext monomorphic() {
        return new VariantContextBuilder().chr("chr1").start(VARIANT_START).stop(VARIANT_START)
                .alleles(List.of(REF)).make();
    }

    static GATKRead read(final String name, final boolean reverse) {
        final SAMRecord record = new SAMRecord(HEADER);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(100);
        record.setCigarString("10M");
        final byte[] bases = new byte[10];
        Arrays.fill(bases, (byte) 'A');
        record.setReadBases(bases);
        final byte[] qualities = new byte[10];
        Arrays.fill(qualities, (byte) 30);
        record.setBaseQualities(qualities);
        record.setMappingQuality(60);
        record.setReadNegativeStrandFlag(reverse);
        return new SAMRecordToGATKReadAdapter(record);
    }

    /** A matrix with the four counts: ref-forward, ref-reverse, alt-forward, alt-reverse. */
    static AlleleLikelihoods<GATKRead, Allele> balanced(final int refFwd, final int refRev,
                                                        final int altFwd, final int altRev) {
        final List<GATKRead> reads = new ArrayList<>();
        final List<Boolean> isRef = new ArrayList<>();
        for (int i = 0; i < refFwd; i++) { reads.add(read("rf" + i, false)); isRef.add(true); }
        for (int i = 0; i < refRev; i++) { reads.add(read("rr" + i, true)); isRef.add(true); }
        for (int i = 0; i < altFwd; i++) { reads.add(read("af" + i, false)); isRef.add(false); }
        for (int i = 0; i < altRev; i++) { reads.add(read("ar" + i, true)); isRef.add(false); }
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", reads);
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(List.of("s1")), new IndexedAlleleList<>(REF, ALT), bySample);
        final LikelihoodMatrix<GATKRead, Allele> matrix = likelihoods.sampleMatrix(0);
        for (int e = 0; e < reads.size(); e++) {
            matrix.set(0, e, isRef.get(e) ? -1 : -10);
            matrix.set(1, e, isRef.get(e) ? -10 : -1);
        }
        return likelihoods;
    }

    /** Two samples with one read each: below the FS threshold separately, above it pooled. */
    static AlleleLikelihoods<GATKRead, Allele> twoSamples() {
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", new ArrayList<>(List.of(read("a", false))));
        bySample.put("s2", new ArrayList<>(List.of(read("b", true))));
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(List.of("s1", "s2")), new IndexedAlleleList<>(REF, ALT),
                bySample);
        for (int s = 0; s < 2; s++) {
            final LikelihoodMatrix<GATKRead, Allele> matrix = likelihoods.sampleMatrix(s);
            matrix.set(0, 0, -1);
            matrix.set(1, 0, -10);
        }
        return likelihoods;
    }

    static void table(final String label, final int minCount,
                      final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        final int[][] t = StrandBiasTest.getContingencyTable(likelihoods, site(null), minCount);
        System.out.printf("table\t%s\t%d\t%d,%d,%d,%d%n", label, minCount,
                t[0][0], t[0][1], t[1][0], t[1][1]);
    }

    static void emit(final String label, final VariantContext vc,
                     final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        one("FisherStrand", label, new FisherStrand(), vc, likelihoods);
        one("StrandOddsRatio", label, new StrandOddsRatio(), vc, likelihoods);
        // SB is a genotype annotation, so it is asked through its own signature.
        try {
            final GenotypeBuilder gb = new GenotypeBuilder("s1", List.of(REF, ALT));
            final Genotype genotype = vc.getGenotypes().isEmpty()
                    ? new GenotypeBuilder("s1", List.of(REF, ALT)).make()
                    : vc.getGenotype(0);
            new StrandBiasBySample().annotate(null, vc, genotype, gb, likelihoods);
            final Object value = gb.make().getAnyAttribute("SB");
            System.out.printf("anno\tStrandBiasBySample\t%s\t%s%n", label,
                    value == null ? "" : String.format("SB=%s[%s]", value,
                            value.getClass().getName()));
        } catch (final Exception | AssertionError e) {
            System.out.printf("anno\tStrandBiasBySample\t%s\tE:%s%n", label,
                    e.getClass().getName());
        }
    }

    static void one(final String name, final String label, final InfoFieldAnnotation annotation,
                    final VariantContext vc,
                    final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        try {
            final Map<String, Object> result = annotation.annotate(null, vc, likelihoods);
            final StringJoiner joiner = new StringJoiner(";");
            for (final Map.Entry<String, Object> entry : result.entrySet()) {
                final Object value = entry.getValue();
                joiner.add(String.format("%s=%s[%s]", entry.getKey(), value,
                        value == null ? "null" : value.getClass().getName()));
            }
            System.out.printf("anno\t%s\t%s\t%s%n", name, label, joiner);
        } catch (final Exception | AssertionError e) {
            System.out.printf("anno\t%s\t%s\tE:%s:%s%n", name, label, e.getClass().getName(),
                    e.getMessage() == null ? "" : e.getMessage().replace('\n', ' '));
        }
    }
}

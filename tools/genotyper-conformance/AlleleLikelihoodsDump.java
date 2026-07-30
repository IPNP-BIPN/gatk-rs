/*
 * The read side of AlleleLikelihoods, taken from the reference.
 *
 * Every likelihood-reading annotation goes through searchBestAllele, and four of its decisions are
 * invisible in its signature.
 *
 *   - bestAlleleIndex and secondBestIndex BOTH start at 0. Nothing means "no second best", so a
 *     one-allele matrix ends with the two equal, which the tail then turns into a second-best
 *     likelihood of negative infinity;
 *   - the comparison is strictly greater, so among equal likelihoods the EARLIEST allele wins and
 *     the allele order is observable in a result that looks order-independent;
 *   - isInformative() compares the confidence against LOG_10_INFORMATIVE_THRESHOLD whatever base
 *     the matrix is in, while the tie-breaking pass inside searchBestAllele converts its threshold
 *     with the base. After switchToNaturalLog() the two ask their question at different thresholds;
 *   - confidence is guarded by `likelihood == secondBestLikelihood ? 0 : ...`, whose load-bearing
 *     case is the pair of negative infinities, where the subtraction would be NaN.
 *
 * The likelihoods travel as raw bits, because they are read back out of the matrix after a
 * conversion and a decimal rendering would hide a divergence in the last place.
 *
 * Output:
 *
 *     counts\t<label>\t<numberOfSamples>\t<numberOfAlleles>\t<evidenceCount>\t<per-sample counts>
 *     best\t<label>\t<sample>:<evidence>\t<best allele>\t<second best>\t<likelihood bits>\t<second bits>\t<confidence bits>\t<informative>
 *     natural\t<label>\t<threshold bits before>\t<threshold bits after>\t<value bits after>
 *
 * Usage: AlleleLikelihoodsDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;

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

public class AlleleLikelihoodsDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele ALT1 = Allele.create("C", false);
    static final Allele ALT2 = Allele.create("G", false);

    static final SAMFileHeader HEADER = makeHeader();

    public static void main(final String[] args) {
        System.out.println("# AlleleLikelihoodsDump: the matrix, the counts and the best-allele search");

        // One sample, one allele: the search ends with best and second best equal.
        emit("one-allele", List.of(REF), Map.of("s1", 1), new double[][][] {{{-1.0}}});

        // The ordinary shape: one sample, two alleles, three pieces of evidence.
        emit("two-alleles", List.of(REF, ALT1), Map.of("s1", 3),
                new double[][][] {{{-1.0, -3.0, -0.5}, {-2.0, -0.5, -0.5}}});

        // A tie, which goes to the earlier allele because the comparison is strict.
        emit("tie", List.of(REF, ALT1, ALT2), Map.of("s1", 1),
                new double[][][] {{{-1.0}, {-1.0}, {-1.0}}});

        // A tie between the two alternates only, with the reference worse.
        emit("tie-between-alts", List.of(REF, ALT1, ALT2), Map.of("s1", 1),
                new double[][][] {{{-5.0}, {-1.0}, {-1.0}}});

        // Differences just above and just below the informative threshold of 0.2.
        emit("just-informative", List.of(REF, ALT1), Map.of("s1", 1),
                new double[][][] {{{-1.0}, {-1.2000001}}});
        emit("just-uninformative", List.of(REF, ALT1), Map.of("s1", 1),
                new double[][][] {{{-1.0}, {-1.1999999}}});
        emit("exactly-threshold", List.of(REF, ALT1), Map.of("s1", 1),
                new double[][][] {{{-1.0}, {-1.2}}});

        // Infinities, where the confidence guard earns its keep.
        emit("all-infinite", List.of(REF, ALT1), Map.of("s1", 1),
                new double[][][] {{{Double.NEGATIVE_INFINITY}, {Double.NEGATIVE_INFINITY}}});
        emit("one-infinite", List.of(REF, ALT1), Map.of("s1", 1),
                new double[][][] {{{-1.0}, {Double.NEGATIVE_INFINITY}}});
        // A NaN, which loses every comparison, so the first allele stays best.
        emit("with-nan", List.of(REF, ALT1), Map.of("s1", 1),
                new double[][][] {{{Double.NaN}, {-1.0}}});

        // Two samples, to show the traversal order and the total evidence count.
        final Map<String, Integer> twoSamples = new LinkedHashMap<>();
        twoSamples.put("s1", 2);
        twoSamples.put("s2", 1);
        emit("two-samples", List.of(REF, ALT1), twoSamples,
                new double[][][] {{{-1.0, -2.0}, {-2.0, -1.0}}, {{-1.0}, {-3.0}}});

        // A sample with no evidence at all, which contributes nothing and is not an error.
        final Map<String, Integer> emptySample = new LinkedHashMap<>();
        emptySample.put("s1", 0);
        emit("empty-sample", List.of(REF, ALT1), emptySample,
                new double[][][] {{{}, {}}});

        // No alleles: the search returns the missing-index answer.
        emit("no-alleles", List.of(), Map.of("s1", 1), new double[][][] {{}});

        // The natural-log switch, and the two thresholds that stop agreeing after it.
        natural("switch", List.of(REF, ALT1), new double[][][] {{{-1.0}, {-1.1}}});
    }

    static SAMFileHeader makeHeader() {
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(
                List.of(new SAMSequenceRecord("chr1", 1000)));
        final SAMFileHeader header = new SAMFileHeader(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return header;
    }

    /** One read, named so the dump can identify it, mapped somewhere harmless. */
    static GATKRead read(final String name) {
        final SAMRecord record = new SAMRecord(HEADER);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(100);
        record.setCigarString("10M");
        record.setReadBases("ACGTACGTAC".getBytes());
        record.setBaseQualities(new byte[] {30, 30, 30, 30, 30, 30, 30, 30, 30, 30});
        return new SAMRecordToGATKReadAdapter(record);
    }

    static AlleleLikelihoods<GATKRead, Allele> build(final List<Allele> alleles,
                                                     final Map<String, Integer> evidenceCounts,
                                                     final double[][][] values) {
        final List<String> sampleNames = new ArrayList<>(evidenceCounts.keySet());
        final Map<String, List<GATKRead>> evidenceBySample = new LinkedHashMap<>();
        for (final String sample : sampleNames) {
            final List<GATKRead> reads = new ArrayList<>();
            for (int i = 0; i < evidenceCounts.get(sample); i++) {
                reads.add(read(sample + "-r" + i));
            }
            evidenceBySample.put(sample, reads);
        }

        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(sampleNames),
                new IndexedAlleleList<>(alleles.toArray(new Allele[0])),
                evidenceBySample);

        for (int s = 0; s < sampleNames.size(); s++) {
            final LikelihoodMatrix<GATKRead, Allele> matrix = likelihoods.sampleMatrix(s);
            for (int a = 0; a < alleles.size(); a++) {
                for (int e = 0; e < values[s][a].length; e++) {
                    matrix.set(a, e, values[s][a][e]);
                }
            }
        }
        return likelihoods;
    }

    static void emit(final String label, final List<Allele> alleles,
                     final Map<String, Integer> evidenceCounts, final double[][][] values) {
        try {
            final AlleleLikelihoods<GATKRead, Allele> likelihoods =
                    build(alleles, evidenceCounts, values);

            final StringJoiner perSample = new StringJoiner(",");
            for (int s = 0; s < likelihoods.numberOfSamples(); s++) {
                perSample.add(Integer.toString(likelihoods.sampleEvidenceCount(s)));
            }
            System.out.printf("counts\t%s\t%d\t%d\t%d\t%s%n", label, likelihoods.numberOfSamples(),
                    likelihoods.numberOfAlleles(), likelihoods.evidenceCount(), perSample);

            for (final AlleleLikelihoods<GATKRead, Allele>.BestAllele best
                    : likelihoods.bestAllelesBreakingTies()) {
                System.out.printf("best\t%s\t%s:%s\t%s\t%s\t%d\t%d\t%d\t%b%n", label, best.sample,
                        best.evidence.getName(), show(best.allele), show(best.second_best_allele),
                        Double.doubleToRawLongBits(best.likelihood),
                        Double.doubleToRawLongBits(best.secondBestLikelihood),
                        Double.doubleToRawLongBits(best.confidence), best.isInformative());
            }
        } catch (final Exception | AssertionError e) {
            System.out.printf("counts\t%s\tE:%s:%s%n", label, e.getClass().getName(),
                    e.getMessage() == null ? "" : e.getMessage().replace('\n', ' '));
        }
    }

    /** The natural-log switch, with the value and both thresholds around it. */
    static void natural(final String label, final List<Allele> alleles, final double[][][] values) {
        final AlleleLikelihoods<GATKRead, Allele> likelihoods =
                build(alleles, Map.of("s1", values[0][0].length), values);

        final AlleleLikelihoods<GATKRead, Allele>.BestAllele before =
                likelihoods.bestAllelesBreakingTies().iterator().next();
        likelihoods.switchToNaturalLog();
        final AlleleLikelihoods<GATKRead, Allele>.BestAllele after =
                likelihoods.bestAllelesBreakingTies().iterator().next();

        System.out.printf("natural\t%s\t%d\t%d\t%b\t%b\t%b%n", label,
                Double.doubleToRawLongBits(before.confidence),
                Double.doubleToRawLongBits(after.confidence),
                before.isInformative(), after.isInformative(), likelihoods.isNaturalLog());

        // Switching twice is refused.
        try {
            likelihoods.switchToNaturalLog();
            System.out.printf("natural\t%s-twice\tno-refusal%n", label);
        } catch (final Exception | AssertionError e) {
            System.out.printf("natural\t%s-twice\tE:%s:%s%n", label, e.getClass().getName(),
                    e.getMessage() == null ? "" : e.getMessage().replace('\n', ' '));
        }
    }

    static String show(final Allele allele) {
        if (allele == null) {
            return "null";
        }
        return allele.getDisplayString() + (allele.isReference() ? "*" : "");
    }

    static List<Allele> unused() {
        return Arrays.asList();
    }
}

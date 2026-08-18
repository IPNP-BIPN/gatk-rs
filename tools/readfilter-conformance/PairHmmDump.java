/*
 * The three `PairHMM` implementations over the same read and haplotype pairs, taken from the
 * reference.
 *
 * `HaplotypeCaller` picks `FASTEST_AVAILABLE`, which resolves to `VectorLoglessPairHMM` -- Intel GKL,
 * AVX, through JNI -- whenever the native library loads, and it loads on the x86-64 runner the
 * goldens come from. So a golden taken today would pin the vectorised results while the port targets
 * the readable Java one, and nobody has checked whether the two produce the same bytes.
 *
 * What this decides:
 *
 *   - IF THE IMPLEMENTATIONS AGREE BIT FOR BIT, the choice is free and the port targets the readable
 *     one;
 *   - IF THEY DISAGREE, the oracle contract has to name one, every dump has to force it with
 *     `--pair-hmm-implementation`, and the divergence is a measured defect like any other.
 *
 * Each likelihood is printed as its raw bits as well as its decimal form, because a comparison of
 * decimal renderings would hide exactly the difference this is looking for.
 *
 * WHICH IMPLEMENTATIONS RAN IS PART OF THE ANSWER. The native library loads on the runner and may
 * not load elsewhere, so the dump prints what it could build and what it could not, and a `loaded`
 * row of `no` is a result rather than a failure.
 *
 * Output:
 *
 *     loaded\t<implementation>\t<yes|no:reason>
 *     likelihood\t<pair>\t<implementation>=<bits>,<decimal>
 *     agree\t<pair>\t<yes|no>
 *
 * Usage: PairHmmDump
 */

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pairhmm.Log10PairHMM;
import org.broadinstitute.hellbender.utils.pairhmm.LoglessPairHMM;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMMInputScoreImputation;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMMInputScoreImputator;
import org.broadinstitute.hellbender.utils.pairhmm.VectorLoglessPairHMM;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class PairHmmDump {

    /** One read against one haplotype, with the four quality tracks the model reads. */
    record Pair(String label, String haplotype, String read, byte baseQual, byte insertionQual,
                byte deletionQual, byte gapContinuation) {
    }

    public static void main(final String[] args) throws Exception {
        System.out.println("# PairHmmDump: LoglessPairHMM, Log10PairHMM and VectorLoglessPairHMM");

        final List<Pair> pairs = List.of(
                // The read is the haplotype: the alignment is exact.
                new Pair("identical", "ACGTACGTACGT", "ACGTACGTACGT", (byte) 30, (byte) 45,
                        (byte) 45, (byte) 10),
                // One substitution in the middle.
                new Pair("one-mismatch", "ACGTACGTACGT", "ACGTAAGTACGT", (byte) 30, (byte) 45,
                        (byte) 45, (byte) 10),
                // The read is shorter than the haplotype.
                new Pair("short-read", "ACGTACGTACGTACGT", "ACGTACGT", (byte) 30, (byte) 45,
                        (byte) 45, (byte) 10),
                // A deletion: the read skips four bases of the haplotype.
                new Pair("deletion", "ACGTACGTACGT", "ACGTACGT", (byte) 30, (byte) 45, (byte) 45,
                        (byte) 10),
                // An insertion: the read carries four bases the haplotype does not.
                new Pair("insertion", "ACGTACGT", "ACGTTTTTACGT", (byte) 30, (byte) 45, (byte) 45,
                        (byte) 10),
                // Base qualities at the bottom and the top of the range.
                new Pair("low-base-quality", "ACGTACGTACGT", "ACGTAAGTACGT", (byte) 2, (byte) 45,
                        (byte) 45, (byte) 10),
                new Pair("high-base-quality", "ACGTACGTACGT", "ACGTAAGTACGT", (byte) 60, (byte) 45,
                        (byte) 45, (byte) 10),
                // The base-quality threshold the model squashes below.
                new Pair("at-the-quality-threshold", "ACGTACGTACGT", "ACGTAAGTACGT", (byte) 18,
                        (byte) 45, (byte) 45, (byte) 10),
                new Pair("below-the-quality-threshold", "ACGTACGTACGT", "ACGTAAGTACGT", (byte) 17,
                        (byte) 45, (byte) 45, (byte) 10),
                // Cheap gaps, which move the indel hypotheses.
                new Pair("cheap-gaps", "ACGTACGTACGT", "ACGTACGT", (byte) 30, (byte) 10, (byte) 10,
                        (byte) 5),
                // A long haplotype, where the vectorised kernel has something to vectorise.
                new Pair("long", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", (byte) 30,
                        (byte) 45, (byte) 45, (byte) 10),
                // A homopolymer, where an indel has many equally good placements.
                new Pair("homopolymer", "AAAAAAAAAAAA", "AAAAAAAAAA", (byte) 30, (byte) 45,
                        (byte) 45, (byte) 10));

        // What could be built. The native library loads on the runner and may not load elsewhere.
        final Map<String, PairHMM> implementations = new LinkedHashMap<>();
        implementations.put("LoglessPairHMM", new LoglessPairHMM());
        implementations.put("Log10PairHMM", new Log10PairHMM(true));
        for (final Map.Entry<String, PairHMM> entry : implementations.entrySet()) {
            System.out.printf("loaded\t%s\tyes%n", entry.getKey());
        }
        for (final VectorLoglessPairHMM.Implementation flavour :
                VectorLoglessPairHMM.Implementation.values()) {
            final String name = "VectorLoglessPairHMM." + flavour.name();
            try {
                implementations.put(name,
                        new VectorLoglessPairHMM(flavour, new PairHMMNativeArguments()));
                System.out.printf("loaded\t%s\tyes%n", name);
            } catch (final Throwable failure) {
                System.out.printf("loaded\t%s\tno:%s%n", name,
                        ReferenceQueryDump.escape(failure.getClass().getSimpleName()));
            }
        }

        for (final Pair pair : pairs) {
            final List<String> renderings = new ArrayList<>();
            for (final Map.Entry<String, PairHMM> entry : implementations.entrySet()) {
                final double likelihood = likelihood(entry.getValue(), pair);
                renderings.add(Double.toString(likelihood));
                System.out.printf("likelihood\t%s\t%s=%016x,%s%n", pair.label(), entry.getKey(),
                        Double.doubleToRawLongBits(likelihood), Double.toString(likelihood));
            }
            final boolean agree = renderings.stream().distinct().count() == 1;
            System.out.printf("agree\t%s\t%s%n", pair.label(), agree ? "yes" : "no");
        }

        for (final PairHMM hmm : implementations.values()) {
            hmm.close();
        }
    }

    static double likelihood(final PairHMM hmm, final Pair pair) {
        final Haplotype haplotype = new Haplotype(pair.haplotype().getBytes(), true);
        final byte[] bases = pair.read().getBytes();
        final byte[] baseQuals = filled(bases.length, pair.baseQual());
        final byte[] insertionQuals = filled(bases.length, pair.insertionQual());
        final byte[] deletionQuals = filled(bases.length, pair.deletionQual());
        final byte[] gcp = filled(bases.length, pair.gapContinuation());

        final GATKRead read =
                ArtificialReadUtils.createArtificialRead(bases, baseQuals, bases.length + "M");
        ReadUtils.setInsertionBaseQualities(read, insertionQuals);
        ReadUtils.setDeletionBaseQualities(read, deletionQuals);

        final PairHMMInputScoreImputator imputator = ignored -> new PairHMMInputScoreImputation() {
            @Override
            public byte[] delOpenPenalties() {
                return deletionQuals;
            }

            @Override
            public byte[] insOpenPenalties() {
                return insertionQuals;
            }

            @Override
            public byte[] gapContinuationPenalties() {
                return gcp;
            }
        };

        // `VectorLoglessPairHMM` overrides the four-argument form to set up its JNI data and works
        // out the maximum lengths itself; the scalar ones inherit a version that only forwards the
        // two zeros, which their own validation then refuses. So each is initialised the way its own
        // class expects.
        if (hmm instanceof VectorLoglessPairHMM) {
            hmm.initialize(Collections.singletonList(haplotype), null, 0, 0);
        } else {
            hmm.initialize(bases.length, pair.haplotype().length());
        }
        hmm.computeLog10Likelihoods(matrix(Collections.singletonList(haplotype)),
                Collections.singletonList(read), imputator);
        return hmm.getLogLikelihoodArray()[0];
    }

    static byte[] filled(final int length, final byte value) {
        final byte[] values = new byte[length];
        java.util.Arrays.fill(values, value);
        return values;
    }

    /** The one method `computeLog10Likelihoods` needs, and nothing else. */
    static LikelihoodMatrix<GATKRead, Haplotype> matrix(final List<Haplotype> haplotypes) {
        return new LikelihoodMatrix<GATKRead, Haplotype>() {
            @Override
            public List<GATKRead> evidence() {
                throw new UnsupportedOperationException();
            }

            @Override
            public List<Haplotype> alleles() {
                return haplotypes;
            }

            @Override
            public void set(final int alleleIndex, final int evidenceIndex, final double value) {
                // The likelihoods are read from `getLogLikelihoodArray`, not from here.
            }

            @Override
            public double get(final int alleleIndex, final int evidenceIndex) {
                throw new UnsupportedOperationException();
            }

            @Override
            public int indexOfAllele(final Allele allele) {
                throw new UnsupportedOperationException();
            }

            @Override
            public int indexOfEvidence(final GATKRead evidence) {
                throw new UnsupportedOperationException();
            }

            @Override
            public int numberOfAlleles() {
                return haplotypes.size();
            }

            @Override
            public int evidenceCount() {
                throw new UnsupportedOperationException();
            }

            @Override
            public Haplotype getAllele(final int alleleIndex) {
                throw new UnsupportedOperationException();
            }

            @Override
            public GATKRead getEvidence(final int evidenceIndex) {
                throw new UnsupportedOperationException();
            }

            @Override
            public void copyAlleleLikelihoods(final int alleleIndex, final double[] dest,
                                              final int offset) {
                throw new UnsupportedOperationException();
            }
        };
    }
}

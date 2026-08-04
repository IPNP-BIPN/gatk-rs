/*
 * SomaticLikelihoodsEngine.alleleFractionsPosterior and the Dirichlet under it, from the reference.
 *
 * The variational fixed point AllelePseudoDepth ends in, and where G1.9's open question lives.
 * NaturalLogUtils established that one ulp is what enters here; this measures what the iteration
 * does with it.
 *
 * THE ITERATION COUNT IS DUMPED ALONGSIDE THE RESULT, and that is the point of the harness rather
 * than a convenience. Convergence is a threshold test:
 *
 *     distance1(old, new) / sum(new) < 0.001
 *
 * so a difference far too small to see in the values can still put one side of that comparison on
 * the other side of the threshold, after which the two runs do a different amount of work. A
 * divergence has to be attributable to amplification or to a different iteration count, and
 * without the count the two are indistinguishable.
 *
 * The count is not exposed by the engine, so it is recovered by re-running the loop's body here
 * with the same calls in the same order: getEffectiveCounts, ebeAdd, then the same convergence
 * test. The result is checked against the engine's own answer on every case, and a row where they
 * disagree would mean this harness has drifted from the code it is measuring.
 *
 * Values are RAW BIT PATTERNS. Decimal rendering discards exactly the difference being measured.
 *
 * Output:
 *
 *     post\t<label>\t<iterations>\t<result bits, comma separated>
 *     agree\t<label>\t<true if the engine's own answer matches this harness's replay>
 *     weights\t<label>\t<effectiveLogMultinomialWeights bits, comma separated>
 *
 * Usage: SomaticLikelihoodsDump
 */

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.MathArrays;

import org.broadinstitute.hellbender.utils.Dirichlet;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.NaturalLogUtils;
import org.broadinstitute.hellbender.tools.walkers.mutect.SomaticLikelihoodsEngine;

import java.util.Arrays;
import java.util.stream.Collectors;

public class SomaticLikelihoodsDump {

    public static void main(final String[] args) {
        System.out.println("# SomaticLikelihoodsDump: the Dirichlet fixed point, with its iteration count");

        // The flat prior, which is what AllelePseudoDepth composes by default.
        final double[] flat2 = {1.0, 1.0};
        final double[] flat3 = {1.0, 1.0, 1.0};

        // One read, strongly on the first allele: converges immediately.
        run("one-read-clean", new double[][] {{-0.001}, {-10.0}}, flat2, null);
        // One read that says nothing: the two alleles stay level.
        run("one-read-flat", new double[][] {{-1.0}, {-1.0}}, flat2, null);
        // Ten reads all on the first allele.
        run("ten-reads-clean", column(10, new double[] {-0.001, -10.0}), flat2, null);
        // Ten reads split five and five, which is the slowest shape to settle.
        run("ten-reads-split", split(5, 5), flat2, null);
        // A near-tie, which is where the convergence threshold is most likely to be decided by a
        // last bit.
        run("near-tie", split(50, 49), flat2, null);
        // Three alleles.
        run("three-alleles", new double[][] {
                {-0.1, -0.2, -5.0}, {-3.0, -0.3, -4.0}, {-6.0, -7.0, -0.05}}, flat3, null);
        // A prior that is not flat, so the addition inside the loop actually does something.
        run("skewed-prior", split(5, 5), new double[] {0.5, 3.0}, null);
        run("large-prior", split(5, 5), new double[] {100.0, 100.0}, null);
        // Weighted, the other branch AllelePseudoDepth can take.
        run("weighted-uniform", split(5, 5), flat2, ones(10));
        run("weighted-decaying", split(5, 5), flat2, decay(10));
        // Likelihoods far apart, so posteriors are saturated and the fixed point barely moves.
        run("saturated", column(4, new double[] {-1e-9, -700.0}), flat2, null);
        // Many reads, to give amplification the most room this harness offers.
        run("fifty-reads", split(25, 25), flat2, null);

        // The Dirichlet weights on their own: digamma, no exp, so these carry no 1-ulp bound and
        // must be bit-identical.
        weights("weights-flat2", flat2);
        weights("weights-flat3", flat3);
        weights("weights-skewed", new double[] {0.5, 3.0});
        weights("weights-large", new double[] {100.0, 100.0});
        weights("weights-tiny", new double[] {1e-3, 1e-3});
    }

    /** `reads` copies of the same per-allele likelihood column. */
    static double[][] column(final int reads, final double[] perAllele) {
        final double[][] matrix = new double[perAllele.length][reads];
        for (int a = 0; a < perAllele.length; a++) {
            Arrays.fill(matrix[a], perAllele[a]);
        }
        return matrix;
    }

    /** `first` reads favouring allele 0, then `second` favouring allele 1. */
    static double[][] split(final int first, final int second) {
        final double[][] matrix = new double[2][first + second];
        for (int r = 0; r < first + second; r++) {
            final boolean firstAllele = r < first;
            matrix[0][r] = firstAllele ? -0.01 : -4.0;
            matrix[1][r] = firstAllele ? -4.0 : -0.01;
        }
        return matrix;
    }

    static double[] ones(final int n) {
        final double[] w = new double[n];
        Arrays.fill(w, 1.0);
        return w;
    }

    /** Weights that decay with read index, which is the shape AllelePseudoDepth computes. */
    static double[] decay(final int n) {
        final double[] w = new double[n];
        for (int i = 0; i < n; i++) {
            w[i] = 1.0 / (1.0 + i);
        }
        return w;
    }

    static String bits(final double value) {
        return Long.toHexString(Double.doubleToRawLongBits(value));
    }

    static String bitList(final double[] values) {
        return Arrays.stream(values).mapToObj(SomaticLikelihoodsDump::bits)
                .collect(Collectors.joining(","));
    }

    static void run(final String label, final double[][] logLikelihoods, final double[] prior,
                    final double[] weights) {
        final RealMatrix matrix = new Array2DRowRealMatrix(logLikelihoods);

        // The engine's own answer.
        final double[] engine = SomaticLikelihoodsEngine.alleleFractionsPosterior(matrix, prior, weights);

        // The same loop, replayed here so the iteration count can be observed. Same calls, same
        // order: if this ever stops agreeing with the engine, the harness has drifted.
        double[] posterior = new double[prior.length];
        Arrays.fill(posterior, 1.0);
        int iterations = 0;
        boolean converged = false;
        while (!converged) {
            final double[] counts = effectiveCounts(matrix, posterior, weights);
            final double[] next = MathArrays.ebeAdd(counts, prior);
            iterations++;
            converged = MathArrays.distance1(posterior, next) / MathUtils.sum(next) < 0.001;
            posterior = next;
        }

        System.out.printf("post\t%s\t%d\t%s%n", label, iterations, bitList(posterior));
        System.out.printf("agree\t%s\t%s%n", label, Arrays.equals(engine, posterior));
    }

    /** `SomaticLikelihoodsEngine.getEffectiveCounts`, which is package-visible, rewritten here. */
    static double[] effectiveCounts(final RealMatrix logLikelihoods, final double[] dirichletPrior,
                                    final double[] weights) {
        final double[] effectiveLogWeights = new Dirichlet(dirichletPrior).effectiveLogMultinomialWeights();
        return MathUtils.sumArrayFunction(0, logLikelihoods.getColumnDimension(),
                read -> {
                    final double[] unweighted =
                            NaturalLogUtils.posteriors(effectiveLogWeights, logLikelihoods.getColumn(read));
                    return weights == null
                            ? unweighted
                            : MathUtils.applyToArrayInPlace(unweighted, d -> d * weights[read]);
                });
    }

    static void weights(final String label, final double[] alpha) {
        System.out.printf("weights\t%s\t%s%n", label,
                bitList(new Dirichlet(alpha).effectiveLogMultinomialWeights()));
    }
}

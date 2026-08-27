/*
 * LearnReadOrientationModel's mixture model, taken from the reference.
 *
 * Which of twelve states a site is in: eight orientation artefacts, one per alternate base per
 * read orientation, and four real states. The tool fits the prior over those states by EM. What is
 * decidable without the counts file is the two functions the fit is built from: the flat prior it
 * starts from, and the responsibilities one site gets given a prior.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE FLAT PRIOR IS NOT FLAT OVER TWELVE STATES: the two artefact states whose alternate base
 *     IS the reference base are given zero, and the remaining ten share the mass, so the value
 *     depends on the reference base;
 *   - A REF-TO-REF ARTEFACT IS IMPOSSIBLE, and is given zero responsibility however the counts
 *     fall;
 *   - AN ARTEFACT STATE WHOSE BASE IS NOT THE OBSERVED ALTERNATE IS ALSO IMPOSSIBLE, so only two
 *     of the eight artefact states can ever be non-zero for one site;
 *   - THE TWO SURVIVING ARTEFACT STATES ARE F1R2 AND F2R1 OF THE OBSERVED BASE, and which of them
 *     takes the mass is decided by the F1R2 count alone;
 *   - `givenNotHomRef` ZEROES HOM_REF AFTER the posteriors are computed, so the remaining states
 *     are renormalised rather than rescaled;
 *   - THE RESPONSIBILITIES ARE NORMALISED FROM LOG SPACE, so they sum to one whatever the prior;
 *   - A STATE WITH A ZERO PRIOR STAYS AT ZERO, which is what makes the flat prior's zeros stick
 *     through every iteration;
 *   - DEPTH AND ALT DEPTH MOVE THE ANSWER SEPARATELY: the same alt fraction at a greater depth is
 *     more certain;
 *   - AND THE STATES HAVE A FIXED ORDER, the eight artefacts first and the four real states after,
 *     which is the order every prior array is indexed in.
 *
 * Output:
 *
 *     order\tstates=<the twelve state names, comma separated>
 *     flat\t<base>=<the twelve prior entries>
 *     resp\t<label>=<the twelve responsibilities>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: LearnReadOrientationModelDump
 */

import htsjdk.samtools.util.Histogram;
import org.broadinstitute.hellbender.tools.walkers.readorientation.ArtifactState;
import org.broadinstitute.hellbender.tools.walkers.readorientation.LearnReadOrientationModelEngine;
import org.broadinstitute.hellbender.utils.Nucleotide;

import java.util.ArrayList;
import java.util.List;

public class LearnReadOrientationModelDump {

    /** Every double printed the same way, so a change of one part in a million shows. */
    static String format(final double[] values) {
        final List<String> parts = new ArrayList<>();
        for (final double value : values) {
            parts.add(String.format("%.10f", value));
        }
        return String.join(",", parts);
    }

    static void flat(final Nucleotide base) {
        System.out.printf("flat\t%s=%s%n", base,
                format(LearnReadOrientationModelEngine.getFlatPrior(base)));
    }

    static void responsibilities(final String label, final Nucleotide reference,
                                 final Nucleotide alternate, final int altDepth,
                                 final int f1r2AltCount, final int depth, final double[] prior,
                                 final boolean givenNotHomRef) {
        try {
            System.out.printf("resp\t%s=%s%n", label,
                    format(LearnReadOrientationModelEngine.computeResponsibilities(reference,
                            alternate, altDepth, f1r2AltCount, depth, prior, givenNotHomRef)));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(cause.getMessage())));
        }
    }

    public static void main(final String[] args) {
        System.out.println("# LearnReadOrientationModelDump: which of twelve states a site is in");

        final List<String> names = new ArrayList<>();
        for (final ArtifactState state : ArtifactState.values()) {
            names.add(state.name());
        }
        System.out.printf("order\tstates=%s%n", String.join(",", names));

        // The flat prior, which is not flat: the two states whose alternate IS the reference get
        // zero and the other ten share the mass.
        for (final Nucleotide base : new Nucleotide[] {Nucleotide.A, Nucleotide.C, Nucleotide.G,
                Nucleotide.T}) {
            flat(base);
        }

        final double[] flatA = LearnReadOrientationModelEngine.getFlatPrior(Nucleotide.A);
        final double[] flatC = LearnReadOrientationModelEngine.getFlatPrior(Nucleotide.C);

        // A site with no alternate read at all: hom ref should take almost everything.
        responsibilities("no-alt", Nucleotide.A, Nucleotide.C, 0, 0, 50, flatA, false);
        // The same site with hom ref forbidden, which renormalises over what is left.
        responsibilities("no-alt-not-hom-ref", Nucleotide.A, Nucleotide.C, 0, 0, 50, flatA, true);
        // A low alt fraction entirely on one orientation, which is what an artefact looks like.
        responsibilities("artifact-f1r2", Nucleotide.A, Nucleotide.C, 5, 5, 50, flatA, false);
        // The same counts on the other orientation.
        responsibilities("artifact-f2r1", Nucleotide.A, Nucleotide.C, 5, 0, 50, flatA, false);
        // The same alt fraction split evenly, which is not an artefact.
        responsibilities("balanced", Nucleotide.A, Nucleotide.C, 6, 3, 50, flatA, false);
        // A half alt fraction, which is a germline het.
        responsibilities("het", Nucleotide.A, Nucleotide.C, 25, 12, 50, flatA, false);
        // Every read alternate, which is hom var.
        responsibilities("hom-var", Nucleotide.A, Nucleotide.C, 50, 25, 50, flatA, false);
        // The same alt fraction at ten times the depth, which is more certain.
        responsibilities("artifact-deep", Nucleotide.A, Nucleotide.C, 50, 50, 500, flatA, false);
        // A different observed alternate, which moves the mass to a different pair of states.
        responsibilities("artifact-g", Nucleotide.A, Nucleotide.G, 5, 5, 50, flatA, false);
        // The reference base as the alternate, which is the ref-to-ref case.
        responsibilities("alt-is-ref", Nucleotide.A, Nucleotide.A, 5, 5, 50, flatA, false);
        // A different reference base, which zeroes a different pair of states.
        responsibilities("ref-c", Nucleotide.C, Nucleotide.A, 5, 5, 50, flatC, false);
        // A prior that already rules out both surviving artefact states.
        final double[] noArtifacts = flatA.clone();
        noArtifacts[ArtifactState.F1R2_C.ordinal()] = 0;
        noArtifacts[ArtifactState.F2R1_C.ordinal()] = 0;
        responsibilities("prior-without-artifacts", Nucleotide.A, Nucleotide.C, 5, 5, 50,
                noArtifacts, false);
        // A prior concentrated entirely on one artefact state.
        final double[] onlyF1R2 = new double[ArtifactState.values().length];
        onlyF1R2[ArtifactState.F1R2_C.ordinal()] = 1.0;
        responsibilities("prior-only-f1r2", Nucleotide.A, Nucleotide.C, 5, 5, 50, onlyF1R2, false);
        // And the same prior with the observation pointing the other way.
        responsibilities("prior-only-f1r2-wrong-way", Nucleotide.A, Nucleotide.C, 5, 0, 50,
                onlyF1R2, false);

        // A histogram whose label is not a canonical kmer, which the engine refuses.
        try {
            final Histogram<Integer> histogram = new Histogram<>("depth", "ACGT");
            new LearnReadOrientationModelEngine(histogram, List.of(), List.of(), 1e-4, 20, 100,
                    null);
            System.out.println("error\tbad-context\tNONE:no exception");
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\tbad-context\t%s:%s%n", cause.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(cause.getMessage())));
        }
        // A canonical label with NO alt sites, which passes both validations and then builds a
        // matrix with zero rows.
        try {
            final Histogram<Integer> histogram = new Histogram<>("depth", "TCA");
            new LearnReadOrientationModelEngine(histogram, List.of(), List.of(), 1e-4, 20, 100,
                    null);
            System.out.println("error\tempty-design-matrix\tNONE:no exception");
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\tempty-design-matrix\t%s:%s%n", cause.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(cause.getMessage())));
        }
        // A label of the right length whose middle base is not canonical, which IS refused by the
        // validation the empty matrix never reaches.
        try {
            final Histogram<Integer> histogram = new Histogram<>("depth", "AGT");
            new LearnReadOrientationModelEngine(histogram, List.of(), List.of(), 1e-4, 20, 100,
                    null);
            System.out.println("error\tnon-canonical\tNONE:no exception");
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\tnon-canonical\t%s:%s%n", cause.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(cause.getMessage())));
        }
    }
}

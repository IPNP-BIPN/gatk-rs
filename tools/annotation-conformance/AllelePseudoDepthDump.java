/*
 * AllelePseudoDepth, taken from the reference.
 *
 * The annotation G1.9 was opened for: DD and DF, an allele depth and an allele fraction read off
 * the Dirichlet posterior and emitted as STRINGS. That is what makes this suite possible at all.
 * The chain underneath calls Math.exp, which cannot be reproduced bit for bit (htsjdk-rs #71), but
 * G1.9.1 measured the port within 1 ulp of it and G1.9.2 measured the fixed point at zero ulp, and
 * the output then goes through DecimalFormat at two and four decimals. So the suite compares the
 * FORMATTED STRINGS, which is the only thing the reference actually writes.
 *
 * The golden carries the inputs as well as the outputs. Every matrix is dumped as raw bit
 * patterns, after any conversion the setup performed, so the Rust side feeds exactly the doubles
 * the reference saw rather than a decimal re-parse of them.
 *
 * Four behaviours this is built to catch, none of them the arithmetic:
 *
 *   - THE MEMO IS WRITTEN THROUGH. composePriorPseudoCounts hands out the array it stores, not a
 *     copy. On the empty-evidence branch the posteriors ARE that array, so the final
 *     `posteriors[i] -= prior[i]` zeroes the cache, and every later genotype with the same allele
 *     count gets a prior of zeros. Every case is therefore annotated TWICE with the same
 *     annotation object, and both answers are dumped;
 *   - THE LOG10 BRANCH INDEXES THE EVIDENCE BY ALLELE. The visitor receives (row, column) as
 *     (allele, read) and looks the mapping quality up at evidence().get(row). More alleles than
 *     reads throws; fewer floors each allele's row with an unrelated read's quality;
 *   - THE DECAY RATE CHANGES WHICH INTRINSICS RUN. Math.pow(10, ...) always runs; the second
 *     Math.pow only when weightDecay != 1.0; neither when weightDecay == 0, which returns null
 *     weights before any of it;
 *   - THE GUARDS EMIT NOTHING. A null likelihoods object and a single allele both return without
 *     touching the builder, and an absent key is not an empty one.
 *
 * Output:
 *
 *     case\t<label>\t<naturalLog>\t<prior>\t<keepPrior>\t<weightDecay>\t<emitted allele indices>
 *              \t<mapping qualities>\t<matrix bits, alleles separated by ';', reads by ','>
 *     out\t<label>\t<call>\t<DD>\t<DF>
 *     err\t<label>\t<call>\t<exception class>:<message>
 *
 * `-` in a DD or DF field means the key was absent, which is not the same as empty.
 *
 * Usage: AllelePseudoDepthDump
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

import org.broadinstitute.hellbender.tools.walkers.annotator.AllelePseudoDepth;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class AllelePseudoDepthDump {

    static final SAMFileHeader HEADER = header();
    static final String SAMPLE = "s1";
    static final int VARIANT_START = 100;

    static final Allele REF = Allele.create("A", true);
    static final Allele ALT1 = Allele.create("C", false);
    static final Allele ALT2 = Allele.create("G", false);

    public static void main(final String[] args) {
        System.out.println("# AllelePseudoDepthDump: DD and DF as the reference formats them");

        // The default settings, on a clean two-allele site. The baseline everything else moves off.
        run("default-two-alleles", setup(2, 10, true, "alternating"), List.of(REF, ALT1), 1.0, false, 1.0);
        run("default-three-alleles", setup(3, 10, true, "alternating"), List.of(REF, ALT1, ALT2),
                1.0, false, 1.0);
        // A site where the reads say nothing, so the posterior barely moves off the prior.
        run("uninformative", setup(2, 10, true, "flat"), List.of(REF, ALT1), 1.0, false, 1.0);

        // EMPTY EVIDENCE. The branch that aliases the memo, and the reason every case runs twice.
        run("empty-evidence", setup(2, 0, true, "alternating"), List.of(REF, ALT1), 1.0, false, 1.0);
        run("empty-evidence-keep-prior", setup(2, 0, true, "alternating"), List.of(REF, ALT1),
                1.0, true, 1.0);
        run("empty-evidence-three", setup(3, 0, true, "alternating"), List.of(REF, ALT1, ALT2),
                1.0, false, 1.0);

        // THE PRIOR. Not one, and kept or subtracted.
        run("keep-prior", setup(2, 9, true, "graded"), List.of(REF, ALT1), 1.0, true, 1.0);
        run("small-prior", setup(2, 9, true, "graded"), List.of(REF, ALT1), 0.1, false, 1.0);
        run("large-prior", setup(2, 9, true, "graded"), List.of(REF, ALT1), 100.0, false, 1.0);
        run("zero-prior", setup(2, 9, true, "graded"), List.of(REF, ALT1), 0.0, false, 1.0);

        // THE DECAY RATE, at the three settings that reach three different sets of intrinsics.
        run("decay-zero", setup(2, 9, true, "graded"), List.of(REF, ALT1), 1.0, false, 0.0);
        run("decay-quadratic", setup(2, 9, true, "graded"), List.of(REF, ALT1), 1.0, false, 2.0);
        run("decay-half", setup(2, 9, true, "graded"), List.of(REF, ALT1), 1.0, false, 0.5);
        run("decay-tiny", setup(2, 9, true, "graded"), List.of(REF, ALT1), 1.0, false, 0.001);

        // THE LOG10 BRANCH, where the mapping quality is looked up by allele index.
        run("log10-two-by-two", setup(2, 2, false, "graded"), List.of(REF, ALT1), 1.0, false, 1.0);
        run("log10-two-by-ten", setup(2, 10, false, "graded"), List.of(REF, ALT1), 1.0, false, 1.0);
        // Three alleles and two reads: evidence().get(2) is out of bounds.
        run("log10-three-by-two", setup(3, 2, false, "graded"), List.of(REF, ALT1, ALT2),
                1.0, false, 1.0);

        // SUBSETTING. The matrix holds three alleles and the variant emits two, so the annotation
        // wraps the matrix rather than using it directly.
        run("subsetted-drop-last", setup(3, 9, true, "graded"), List.of(REF, ALT1), 1.0, false, 1.0);
        run("subsetted-drop-middle", setup(3, 9, true, "graded"), List.of(REF, ALT2), 1.0, false, 1.0);

        // THE GUARDS, which put no key in the genotype at all.
        run("single-allele", setup(2, 10, true, "alternating"), List.of(REF), 1.0, false, 1.0);
        runNullLikelihoods("null-likelihoods");

        // A single read, where the fixed point converges immediately, and many reads, where it
        // does not.
        run("one-read", setup(2, 1, true, "graded"), List.of(REF, ALT1), 1.0, false, 1.0);
        run("fifty-reads", setup(2, 50, true, "graded"), List.of(REF, ALT1), 1.0, false, 1.0);
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(
                List.of(new SAMSequenceRecord("chr1", 100000))));
        return header;
    }

    /**
     * A read whose mapping quality varies with its index, so the log10 branch's allele-indexed
     * lookup produces visibly different floors rather than the same one twice.
     */
    static GATKRead read(final int index) {
        final SAMRecord record = new SAMRecord(HEADER);
        record.setReadName("r" + index);
        record.setReferenceName("chr1");
        record.setAlignmentStart(VARIANT_START + index);
        record.setCigarString("20M");
        final byte[] bases = new byte[20];
        Arrays.fill(bases, (byte) 'A');
        record.setReadBases(bases);
        final byte[] qualities = new byte[20];
        Arrays.fill(qualities, (byte) 30);
        record.setBaseQualities(qualities);
        record.setMappingQuality(20 + 10 * index);
        return new SAMRecordToGATKReadAdapter(record);
    }

    /**
     * A likelihood matrix with `alleles` rows and `reads` columns.
     *
     * The values are in log10 as the engine's own matrices are, and switchToNaturalLog scales them
     * when the natural-log branch is wanted.
     *
     * The shape matters more than it looks. `alternating` gives every read the same margin, and a
     * site like that is SYMMETRIC: the weights are all equal, so every setting of weightDecay
     * produces the same answer and the cases meant to separate the three intrinsic paths separate
     * nothing. `graded` gives each read a different margin, which is what makes a weight a weight.
     */
    static AlleleLikelihoods<GATKRead, Allele> setup(final int alleles, final int reads,
                                                     final boolean naturalLog, final String shape) {
        final List<Allele> alleleList = List.of(REF, ALT1, ALT2).subList(0, alleles);
        final List<GATKRead> readList = new ArrayList<>();
        for (int i = 0; i < reads; i++) {
            readList.add(read(i));
        }
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put(SAMPLE, readList);
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(List.of(SAMPLE)), new IndexedAlleleList<>(alleleList), bySample);

        final LikelihoodMatrix<GATKRead, Allele> matrix = likelihoods.sampleMatrix(0);
        for (int a = 0; a < alleles; a++) {
            for (int e = 0; e < reads; e++) {
                // Reads alternate which allele they favour, so a site is genuinely mixed rather
                // than unanimous, and the fixed point has somewhere to go.
                final double value;
                switch (shape) {
                    case "flat":
                        value = -1.0;
                        break;
                    case "graded":
                        // Each read's margin differs, so the weights differ and weightDecay has
                        // something to act on.
                        value = (e % alleles == a) ? -0.01 : -(0.2 + 0.35 * e + a);
                        break;
                    default:
                        value = (e % alleles == a) ? -0.01 : -4.0 - a;
                        break;
                }
                matrix.set(a, e, value);
            }
        }
        if (naturalLog) {
            likelihoods.switchToNaturalLog();
        }
        return likelihoods;
    }

    static VariantContext site(final List<Allele> alleles) {
        return new VariantContextBuilder().chr("chr1").start(VARIANT_START).stop(VARIANT_START)
                .alleles(alleles).make();
    }

    static void run(final String label, final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                    final List<Allele> emitted, final double prior, final boolean keepPrior,
                    final double weightDecay) {
        describe(label, likelihoods, emitted, prior, keepPrior, weightDecay);

        final AllelePseudoDepth annotation = new AllelePseudoDepth();
        annotation.prior = prior;
        annotation.keepPriorInCount = keepPrior;
        annotation.weightDecay = weightDecay;

        // TWICE, with the same annotation object. The second call is where a memo that was written
        // through shows up, and it is invisible to a suite that only ever calls once.
        for (int call = 1; call <= 2; call++) {
            emit(label, call, annotation, likelihoods, emitted);
        }
    }

    static void runNullLikelihoods(final String label) {
        System.out.printf("case\t%s\t-\t1.0\tfalse\t1.0\t-\t-\t-%n", label);
        final AllelePseudoDepth annotation = new AllelePseudoDepth();
        for (int call = 1; call <= 2; call++) {
            emit(label, call, annotation, null, List.of(REF, ALT1));
        }
    }

    static void emit(final String label, final int call, final AllelePseudoDepth annotation,
                     final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                     final List<Allele> emitted) {
        final GenotypeBuilder builder = new GenotypeBuilder(SAMPLE, List.of(REF, REF));
        try {
            annotation.annotate(null, site(emitted), new GenotypeBuilder(SAMPLE,
                    List.of(REF, REF)).make(), builder, likelihoods);
        } catch (final Exception | AssertionError e) {
            System.out.printf("err\t%s\t%d\t%s:%s%n", label, call, e.getClass().getName(),
                    e.getMessage());
            return;
        }
        final Genotype genotype = builder.make();
        System.out.printf("out\t%s\t%d\t%s\t%s%n", label, call,
                value(genotype, GATKVCFConstants.PSEUDO_DEPTH_KEY),
                value(genotype, GATKVCFConstants.PSEUDO_FRACTION_KEY));
    }

    /** `-` when the key is absent, which is not the same as an empty string. */
    static String value(final Genotype genotype, final String key) {
        final Object value = genotype.getExtendedAttribute(key);
        return value == null ? "-" : value.toString();
    }

    static void describe(final String label, final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                         final List<Allele> emitted, final double prior, final boolean keepPrior,
                         final double weightDecay) {
        final LikelihoodMatrix<GATKRead, Allele> matrix = likelihoods.sampleMatrix(0);
        // The emitted alleles by their row in the matrix, which is what SubsettedLikelihoodMatrix
        // renumbers and what the port needs to reproduce the same subset.
        final String indices = emitted.stream()
                .map(a -> Integer.toString(likelihoods.indexOfAllele(a)))
                .collect(Collectors.joining(","));
        final String qualities = matrix.evidence().stream()
                .map(r -> Integer.toString(r.getMappingQuality()))
                .collect(Collectors.joining(","));
        final StringBuilder rows = new StringBuilder();
        for (int a = 0; a < likelihoods.numberOfAlleles(); a++) {
            if (a > 0) {
                rows.append(';');
            }
            for (int e = 0; e < matrix.evidenceCount(); e++) {
                if (e > 0) {
                    rows.append(',');
                }
                rows.append(Long.toHexString(Double.doubleToRawLongBits(matrix.get(a, e))));
            }
        }
        System.out.printf("case\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%n", label,
                likelihoods.isNaturalLog(), Double.toString(prior), keepPrior,
                Double.toString(weightDecay), indices, qualities.isEmpty() ? "-" : qualities,
                rows.length() == 0 ? "-" : rows.toString());
    }
}

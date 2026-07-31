/*
 * ExcessHet, InbreedingCoeff, RMSMappingQuality, MappingQualityZero and TandemRepeat, taken from
 * the reference.
 *
 *   - the same genotype counts feed ExcessHet (rounded) and InbreedingCoeff (not rounded), so the
 *     two annotations on one site can disagree about how many hets there are. Both flags are
 *     dumped for every cohort;
 *   - the rounded branch guards against counting a GQ of zero twice by comparing two normalised
 *     likelihoods for exact equality, which is reachable from PLs of [0,0,X] and [X,0,0];
 *   - a hom-ref genotype with no PLs is given an invented distribution from its GQ, three
 *     different ways depending on the rounding flag and on whether the GQ is zero;
 *   - ExcessHet saturates at 160.0000 below a p-value of 10e-60, which is 1e-59;
 *   - MQ drops reads whose mapping quality is 255 from BOTH the sum and the count, so an all-255
 *     matrix divides by zero and the annotation is the four characters NaN;
 *   - MQ0 returns 0 for an empty matrix where every other annotation returns nothing;
 *   - findRepeatedSubstring cannot see a partial trailing repeat, because Arrays.copyOfRange pads
 *     with zero bytes past the end, and returns 1 with a zero-byte unit on an empty input.
 *
 * Output:
 *
 *     counts\t<label>\t<rounded|raw>\t<refs>\t<hets>\t<homs>
 *     eh\t<label>\t<value or E:class>
 *     ic\t<label>\t<value or E:class>
 *     mq\t<label>\t<key>=<value>[<class>];...
 *     rawmq\t<label>\t<key>=<value>[<class>];...
 *     mq0\t<label>\t<key>=<value>[<class>];...
 *     repunit\t<bases>\t<length>
 *     reps\t<unit>\t<test>\t<leading>\t<count or E:class>
 *     str\t<label>\t<key>=<value>[<class>];...
 *
 * Usage: HeterozygosityAndMqDump
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

import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.tools.walkers.annotator.ExcessHet;
import org.broadinstitute.hellbender.tools.walkers.annotator.InbreedingCoeff;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.MappingQualityZero;
import org.broadinstitute.hellbender.tools.walkers.annotator.RMSMappingQuality;
import org.broadinstitute.hellbender.tools.walkers.annotator.TandemRepeat;
import org.broadinstitute.hellbender.utils.GenotypeCounts;
import org.broadinstitute.hellbender.utils.GenotypeUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

public class HeterozygosityAndMqDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele ALT = Allele.create("C", false);
    static final Allele ALT2 = Allele.create("G", false);
    static final SAMFileHeader HEADER = makeHeader();
    static final int START = 105;

    /** PLs for a hom-ref, a het and a hom-var call at a biallelic site. */
    static final int[] HOM_REF = {0, 60, 600};
    static final int[] HET = {60, 0, 60};
    static final int[] HOM_VAR = {600, 60, 0};
    /** [0,0,X]: the two likelihoods normalise to exactly one half each. */
    static final int[] GQ_ZERO_REF_HET = {0, 0, 60};
    /** [X,0,0]: het and hom-var tie. */
    static final int[] GQ_ZERO_HET_VAR = {60, 0, 0};
    /** [0,0,0]: all three tie, which is why the variant count is tracked separately. */
    static final int[] FLAT = {0, 0, 0};

    public static void main(final String[] args) {
        System.out.println("# HeterozygosityAndMqDump: ExcessHet, InbreedingCoeff, MQ, MQ0, STR");

        // Cohorts in and out of Hardy-Weinberg equilibrium.
        cohort("equilibrium", repeat(HOM_REF, 25), repeat(HET, 50), repeat(HOM_VAR, 25));
        cohort("all-het", repeat(HET, 20));
        cohort("all-hom-ref", repeat(HOM_REF, 20));
        cohort("all-hom-var", repeat(HOM_VAR, 20));
        cohort("excess-het-small", repeat(HOM_REF, 2), repeat(HET, 8));
        cohort("excess-het-large", repeat(HOM_REF, 5), repeat(HET, 40), repeat(HOM_VAR, 5));
        cohort("saturating", repeat(HET, 200));
        cohort("nine-samples", repeat(HET, 5), repeat(HOM_REF, 4));
        cohort("ten-samples", repeat(HET, 5), repeat(HOM_REF, 5));
        cohort("gq-zero-ref-het", repeat(GQ_ZERO_REF_HET, 10));
        cohort("gq-zero-het-var", repeat(GQ_ZERO_HET_VAR, 10));
        cohort("flat", repeat(FLAT, 10));
        cohort("mixed-ties", repeat(GQ_ZERO_REF_HET, 5), repeat(GQ_ZERO_HET_VAR, 5));

        // Hom-ref without likelihoods, which is given a distribution invented from its GQ.
        gqOnly("gq-only-zero", 0, 12);
        gqOnly("gq-only-thirty", 30, 12);
        gqOnly("gq-only-ninetynine", 99, 12);

        // Multiallelic, where the best PL decides which pair of columns is used.
        multiallelic("multiallelic-het", 12);
        multiallelic("multiallelic-no-ref", 12);

        // Empty and monomorphic.
        emptyCohort("no-genotypes");
        monomorphic("monomorphic");

        // MQ and MQ0.
        mq("ordinary", new int[] {60, 60, 60});
        mq("with-zeroes", new int[] {60, 0, 60});
        mq("all-zero", new int[] {0, 0, 0});
        mq("all-unavailable", new int[] {255, 255});
        mq("mixed-unavailable", new int[] {60, 255, 30});
        mq("single-read", new int[] {37});
        mq("empty", new int[0]);
        mqNullLikelihoods("null-likelihoods");

        // The raw tuple and its finalisation.
        rawRoundTrip("raw-ordinary", 60 * 60 * 3, 3);
        rawRoundTrip("raw-zero-depth", 0, 0);
        rawRoundTrip("raw-one", 3600, 1);

        // The repeat arithmetic, on its own.
        repeatUnit("");
        repeatUnit("A");
        repeatUnit("AT");
        repeatUnit("ATAT");
        repeatUnit("ACTACT");
        repeatUnit("ACTACTAC");
        repeatUnit("CCCCCCCC");
        repeatUnit("ACTACA");

        repetitions("AT", "ATATG", true);
        repetitions("AT", "GATAT", true);
        repetitions("AT", "GATAT", false);
        repetitions("A", "ATATG", true);
        repetitions("CCC", "CCCCCCCC", true);
        repetitions("CCC", "CCCCCCCC", false);
        repetitions("AT", "", true);

        // TandemRepeat, over a real reference window.
        str("deletion-of-one-unit", "GATCCACCACCAGTCGA", 100, 102, "TCCA", "T");
        str("insertion-of-one-unit", "GATCCACCACCAGTCGA", 100, 102, "T", "TCCA");
        str("not-a-repeat", "GATCCACCACCAGTCGA", 100, 102, "TC", "T");
        str("snp-is-not-an-indel", "GATCCACCACCAGTCGA", 100, 102, "T", "G");
        str("homopolymer", "GAAAAAAAAAAAAGTCGA", 100, 102, "AA", "A");
        str("multiallelic-indel", "GATCCACCACCAGTCGA", 100, 102, "TCCA", "T,TCCACCA");
    }

    static SAMFileHeader makeHeader() {
        final SAMSequenceDictionary dictionary =
                new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", 1000)));
        final SAMFileHeader header = new SAMFileHeader(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return header;
    }

    static List<int[]> repeat(final int[] pls, final int count) {
        final List<int[]> out = new ArrayList<>();
        for (int i = 0; i < count; i++) {
            out.add(pls);
        }
        return out;
    }

    @SafeVarargs
    static void cohort(final String label, final List<int[]>... groups) {
        final List<Genotype> genotypes = new ArrayList<>();
        int index = 0;
        for (final List<int[]> group : groups) {
            for (final int[] pls : group) {
                genotypes.add(new GenotypeBuilder("s" + index++, callFor(pls)).PL(pls).make());
            }
        }
        emitCohort(label, new VariantContextBuilder().chr("chr1").start(START).stop(START)
                .alleles(List.of(REF, ALT)).genotypes(genotypes).make());
    }

    /** The called alleles the PLs imply, chosen so the genotype is called and diploid. */
    static List<Allele> callFor(final int[] pls) {
        int best = 0;
        for (int i = 1; i < pls.length; i++) {
            if (pls[i] < pls[best]) {
                best = i;
            }
        }
        if (best == 0) {
            return List.of(REF, REF);
        }
        return best == 1 ? List.of(REF, ALT) : List.of(ALT, ALT);
    }

    static void gqOnly(final String label, final int gq, final int count) {
        final List<Genotype> genotypes = new ArrayList<>();
        for (int i = 0; i < count; i++) {
            genotypes.add(new GenotypeBuilder("s" + i, List.of(REF, REF)).GQ(gq).make());
        }
        emitCohort(label, new VariantContextBuilder().chr("chr1").start(START).stop(START)
                .alleles(List.of(REF, ALT)).genotypes(genotypes).make());
    }

    static void multiallelic(final String label, final int count) {
        final boolean noRef = label.endsWith("no-ref");
        // Six PL entries for three alleles: AA, AB, BB, AC, BC, CC.
        final int[] pls = noRef ? new int[] {600, 300, 60, 300, 0, 60}
                : new int[] {60, 0, 60, 60, 60, 600};
        final List<Genotype> genotypes = new ArrayList<>();
        for (int i = 0; i < count; i++) {
            genotypes.add(new GenotypeBuilder("s" + i,
                    noRef ? List.of(ALT, ALT2) : List.of(REF, ALT)).PL(pls).make());
        }
        emitCohort(label, new VariantContextBuilder().chr("chr1").start(START).stop(START)
                .alleles(List.of(REF, ALT, ALT2)).genotypes(genotypes).make());
    }

    static void emptyCohort(final String label) {
        emitCohort(label, new VariantContextBuilder().chr("chr1").start(START).stop(START)
                .alleles(List.of(REF, ALT)).make());
    }

    static void monomorphic(final String label) {
        final List<Genotype> genotypes = new ArrayList<>();
        for (int i = 0; i < 12; i++) {
            genotypes.add(new GenotypeBuilder("s" + i, List.of(REF, REF)).PL(HOM_REF).make());
        }
        emitCohort(label, new VariantContextBuilder().chr("chr1").start(START).stop(START)
                .alleles(List.of(REF)).genotypes(genotypes).make());
    }

    static void emitCohort(final String label, final VariantContext vc) {
        counts(label, vc, true);
        counts(label, vc, false);
        try {
            emitMap("eh", label, new ExcessHet().annotate(null, vc, null));
        } catch (final Exception | AssertionError e) {
            System.out.printf("eh\t%s\tE:%s%n", label, e.getClass().getName());
        }
        try {
            emitMap("ic", label, new InbreedingCoeff().annotate(null, vc, null));
        } catch (final Exception | AssertionError e) {
            System.out.printf("ic\t%s\tE:%s%n", label, e.getClass().getName());
        }
    }

    static void counts(final String label, final VariantContext vc, final boolean rounded) {
        try {
            final GenotypeCounts counts =
                    GenotypeUtils.computeDiploidGenotypeCounts(vc, vc.getGenotypes(), rounded);
            // The raw bits, because these are the numbers the two annotations divide and index by
            // and a decimal rendering would hide the last ulp.
            System.out.printf("counts\t%s\t%s\t%d\t%d\t%d%n", label, rounded ? "rounded" : "raw",
                    Double.doubleToRawLongBits(counts.getRefs()),
                    Double.doubleToRawLongBits(counts.getHets()),
                    Double.doubleToRawLongBits(counts.getHoms()));
        } catch (final Exception | AssertionError e) {
            System.out.printf("counts\t%s\t%s\tE:%s%n", label, rounded ? "rounded" : "raw",
                    e.getClass().getName());
        }
    }

    static void mq(final String label, final int[] mappingQualities) {
        final List<GATKRead> reads = new ArrayList<>();
        for (int i = 0; i < mappingQualities.length; i++) {
            reads.add(read("r" + i, mappingQualities[i]));
        }
        emitMq(label, matrix(reads));
    }

    static void mqNullLikelihoods(final String label) {
        emitMq(label, null);
    }

    static void emitMq(final String label, final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        final VariantContext vc = new VariantContextBuilder().chr("chr1").start(START).stop(START)
                .alleles(List.of(REF, ALT)).make();
        one("mq", label, new RMSMappingQuality(), vc, likelihoods);
        one("mq0", label, new MappingQualityZero(), vc, likelihoods);
        try {
            emitMap("rawmq", label, new RMSMappingQuality().annotateRawData(null, vc, likelihoods));
        } catch (final Exception | AssertionError e) {
            System.out.printf("rawmq\t%s\tE:%s%n", label, e.getClass().getName());
        }
    }

    static void rawRoundTrip(final String label, final long squareSum, final long depth) {
        final String raw = String.format("%d,%d", squareSum, depth);
        final VariantContext vc = new VariantContextBuilder().chr("chr1").start(START).stop(START)
                .alleles(List.of(REF, ALT)).attribute("RAW_MQandDP", raw).make();
        try {
            emitMap("finalized", label, new RMSMappingQuality().finalizeRawData(vc, vc));
        } catch (final Exception | AssertionError e) {
            System.out.printf("finalized\t%s\tE:%s%n", label, e.getClass().getName());
        }
    }

    static void repeatUnit(final String bases) {
        System.out.printf("repunit\t%s\t%d%n", bases,
                GATKVariantContextUtils.findRepeatedSubstring(bases.getBytes()));
    }

    static void repetitions(final String unit, final String test, final boolean leading) {
        try {
            System.out.printf("reps\t%s\t%s\t%b\t%d%n", unit, test, leading,
                    GATKVariantContextUtils.findNumberOfRepetitions(
                            unit.getBytes(), test.getBytes(), leading));
        } catch (final Exception | AssertionError e) {
            System.out.printf("reps\t%s\t%s\t%b\tE:%s%n", unit, test, leading,
                    e.getClass().getName());
        }
    }

    static void str(final String label, final String bases, final int windowStart,
                    final int variantStart, final String refAllele, final String altAlleles) {
        final SimpleInterval window = new SimpleInterval("chr1", windowStart,
                windowStart + bases.length() - 1);
        final org.broadinstitute.hellbender.utils.reference.ReferenceBases sequence =
                new org.broadinstitute.hellbender.utils.reference.ReferenceBases(
                        bases.getBytes(), window);
        final ReferenceContext context = new ReferenceContext(
                new ReferenceMemorySource(sequence, HEADER.getSequenceDictionary()), window);
        final List<Allele> alleles = new ArrayList<>();
        alleles.add(Allele.create(refAllele, true));
        for (final String alt : altAlleles.split(",")) {
            alleles.add(Allele.create(alt, false));
        }
        final VariantContext vc = new VariantContextBuilder().chr("chr1").start(variantStart)
                .stop(variantStart + refAllele.length() - 1).alleles(alleles).make();
        try {
            emitMap("str", label, new TandemRepeat().annotate(context, vc, null));
        } catch (final Exception | AssertionError e) {
            System.out.printf("str\t%s\tE:%s%n", label, e.getClass().getName());
        }
    }

    static void one(final String kind, final String label, final InfoFieldAnnotation annotation,
                    final VariantContext vc,
                    final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        try {
            emitMap(kind, label, annotation.annotate(null, vc, likelihoods));
        } catch (final Exception | AssertionError e) {
            System.out.printf("%s\t%s\tE:%s%n", kind, label, e.getClass().getName());
        }
    }

    static void emitMap(final String kind, final String label, final Map<String, Object> result) {
        final StringJoiner joiner = new StringJoiner(";");
        if (result != null) {
            for (final Map.Entry<String, Object> entry : result.entrySet()) {
                joiner.add(String.format("%s=%s[%s]", entry.getKey(), entry.getValue(),
                        entry.getValue().getClass().getName()));
            }
        }
        System.out.printf("%s\t%s\t%s%n", kind, label, joiner);
    }

    static AlleleLikelihoods<GATKRead, Allele> matrix(final List<GATKRead> reads) {
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", reads);
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(List.of("s1")), new IndexedAlleleList<>(REF, ALT), bySample);
        final LikelihoodMatrix<GATKRead, Allele> m = likelihoods.sampleMatrix(0);
        for (int e = 0; e < reads.size(); e++) {
            m.set(0, e, e % 2 == 0 ? -1 : -10);
            m.set(1, e, e % 2 == 0 ? -10 : -1);
        }
        return likelihoods;
    }

    static GATKRead read(final String name, final int mappingQuality) {
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
        record.setMappingQuality(mappingQuality);
        return new SAMRecordToGATKReadAdapter(record);
    }
}

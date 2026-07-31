/*
 * AS_QualByDepth, AS_RMSMappingQuality and AS_InbreedingCoeff, taken from the reference.
 *
 *   - AS_QD has no raw data of its own: annotateRawData returns an empty map and annotate is
 *     finalizeRawData(vc, vc) with the same context passed twice;
 *   - AS_QD adds the REFERENCE depth to every alternate's denominator, "to match biallelic case",
 *     so a triallelic site's denominators overlap;
 *   - AS_QD inherits fixTooHighQD once per alternate, so one entry of a comma-separated field can
 *     be deterministic and the next drawn from a random generator;
 *   - AS_MQ divides a sum of squares taken from the likelihoods by an AD summed from the genotypes,
 *     and nothing checks the two describe the same reads;
 *   - AS_MQ's getADcounts counts every genotype, while AS_QD's getAlleleDepths counts only het and
 *     hom-var genotypes whose alternate depth exceeds one. Two definitions of the same depth;
 *   - AS_InbreedingCoeff derives its allele frequency from the allele counts when biallelic and
 *     from getCalledChrCount otherwise, and the two branches do not agree.
 *
 * Output:
 *
 *     asqd\t<label>\t<key>=<value>[<class>];...
 *     asqdraw\t<label>\t<key>=<value>[<class>];...
 *     asdepths\t<label>\t<comma-separated or null>
 *     asmq\t<label>\t<key>=<value>[<class>];...
 *     asmqraw\t<label>\t<key>=<value>[<class>];...
 *     asmqfinal\t<label>\t<key>=<value>[<class>];...
 *     asic\t<label>\t<key>=<value>[<class>];...
 *     hetcounts\t<label>\t<sampleCount>\t<allele>=<bits>;...\t<allele>=<bits>;...
 *
 * Usage: AlleleSpecificSiteStatisticsDump
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

import org.broadinstitute.hellbender.tools.walkers.annotator.HeterozygosityCalculator;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_InbreedingCoeff;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_QualByDepth;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_RMSMappingQuality;
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

public class AlleleSpecificSiteStatisticsDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele ALT = Allele.create("C", false);
    static final Allele ALT2 = Allele.create("G", false);
    static final SAMFileHeader HEADER = makeHeader();
    static final int START = 105;

    public static void main(final String[] args) {
        System.out.println("# AlleleSpecificSiteStatisticsDump: AS_QD, AS_MQ and AS_InbreedingCoeff");

        for (final String label : new String[] {
                "biallelic", "triallelic", "as-qual-approx", "as-vardp", "no-qual", "no-genotypes",
                "no-ad", "one-alt-read", "high-qd", "mixed-qd", "empty-slot", "hom-ref-only"}) {
            asqd(label);
        }

        for (final String label : new String[] {
                "biallelic", "triallelic", "ref-only", "alt-only", "all-unavailable",
                "mixed-unavailable", "no-ad", "null-likelihoods"}) {
            asmq(label);
        }

        for (final String raw : new String[] {
                "3600.00|1600.00", "3600.00|", "|1600.00", "3600.00|1600.00|400.00", "0.00|0.00"}) {
            asmqFinal(raw);
        }

        for (final String label : new String[] {
                "ten-het", "ten-hom-ref", "nine-samples", "twenty-mixed", "triallelic",
                "gq-only", "no-genotypes"}) {
            asic(label);
            hetCounts(label);
        }
    }

    static SAMFileHeader makeHeader() {
        final SAMSequenceDictionary dictionary =
                new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", 1000)));
        final SAMFileHeader header = new SAMFileHeader(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return header;
    }

    // ---- AS_QualByDepth ----------------------------------------------------------------------

    static VariantContext qdContext(final String label) {
        final boolean triallelic = "triallelic".equals(label) || "mixed-qd".equals(label)
                || "empty-slot".equals(label);
        final List<Allele> alleles = triallelic ? List.of(REF, ALT, ALT2) : List.of(REF, ALT);
        final VariantContextBuilder builder = new VariantContextBuilder().chr("chr1")
                .start(START).stop(START).alleles(alleles);
        switch (label) {
            case "as-qual-approx": builder.attribute("AS_QUALapprox", "0|300"); break;
            case "empty-slot": builder.attribute("AS_QUALapprox", "0|300|"); break;
            case "as-vardp":
                builder.attribute("AS_QUAL", "300").attribute("AS_VarDP", "10|8");
                break;
            case "no-qual": break;
            case "high-qd": builder.attribute("AS_QUAL", "5000"); break;
            case "mixed-qd": builder.attribute("AS_QUAL", "70,5000"); break;
            case "triallelic": builder.attribute("AS_QUAL", "70,140"); break;
            default: builder.attribute("AS_QUAL", "300"); break;
        }
        if ("no-genotypes".equals(label)) {
            return builder.make();
        }
        final List<Genotype> genotypes = new ArrayList<>();
        if ("hom-ref-only".equals(label)) {
            genotypes.add(new GenotypeBuilder("s0", List.of(REF, REF))
                    .AD(triallelic ? new int[] {10, 0, 0} : new int[] {10, 0}).make());
        } else if ("no-ad".equals(label)) {
            genotypes.add(new GenotypeBuilder("s0", List.of(REF, ALT)).DP(17).make());
        } else if ("one-alt-read".equals(label)) {
            genotypes.add(new GenotypeBuilder("s0", List.of(REF, ALT))
                    .AD(new int[] {10, 1}).make());
        } else {
            genotypes.add(new GenotypeBuilder("s0", List.of(REF, ALT))
                    .AD(triallelic ? new int[] {10, 4, 6} : new int[] {10, 8}).make());
        }
        return builder.genotypes(genotypes).make();
    }

    static void asqd(final String label) {
        final VariantContext vc = qdContext(label);
        final AS_QualByDepth annotation = new AS_QualByDepth();
        try {
            emitMap("asqd", label, annotation.annotate(null, vc, null));
        } catch (final Exception | AssertionError e) {
            System.out.printf("asqd\t%s\tE:%s%n", label, e.getClass().getName());
        }
        try {
            emitMap("asqdraw", label, annotation.annotateRawData(null, vc, null));
        } catch (final Exception | AssertionError e) {
            System.out.printf("asqdraw\t%s\tE:%s%n", label, e.getClass().getName());
        }
        try {
            final List<Integer> depths = AS_QualByDepth.getAlleleDepths(vc.getGenotypes());
            final StringJoiner joiner = new StringJoiner(",");
            if (depths == null) {
                System.out.printf("asdepths\t%s\tnull%n", label);
            } else {
                for (final Integer depth : depths) {
                    joiner.add(depth.toString());
                }
                System.out.printf("asdepths\t%s\t%s%n", label, joiner);
            }
        } catch (final Exception | AssertionError e) {
            System.out.printf("asdepths\t%s\tE:%s%n", label, e.getClass().getName());
        }
    }

    // ---- AS_RMSMappingQuality ----------------------------------------------------------------

    static int[][] mqComposition(final String label) {
        switch (label) {
            case "biallelic": return new int[][] {{60, 60, 60}, {30, 30}};
            case "triallelic": return new int[][] {{60, 60}, {30, 30}, {20, 20}};
            case "ref-only": return new int[][] {{60, 60, 60}, {}};
            case "alt-only": return new int[][] {{}, {30, 30}};
            case "all-unavailable": return new int[][] {{255, 255}, {255, 255}};
            case "mixed-unavailable": return new int[][] {{60, 255}, {30, 255}};
            case "no-ad": return new int[][] {{60, 60}, {30, 30}};
            case "null-likelihoods": return null;
            default: throw new IllegalArgumentException(label);
        }
    }

    static VariantContext mqContext(final String label) {
        final List<Allele> alleles = "triallelic".equals(label)
                ? List.of(REF, ALT, ALT2) : List.of(REF, ALT);
        final GenotypeBuilder gb = new GenotypeBuilder("s1", List.of(REF, ALT));
        if (!"no-ad".equals(label)) {
            gb.AD("triallelic".equals(label) ? new int[] {2, 2, 2} : new int[] {3, 2});
        }
        return new VariantContextBuilder().chr("chr1").start(START).stop(START)
                .alleles(alleles).genotypes(List.of(gb.make())).make();
    }

    static AlleleLikelihoods<GATKRead, Allele> mqLikelihoods(final String label) {
        final int[][] composition = mqComposition(label);
        if (composition == null) {
            return null;
        }
        final List<Allele> alleles = "triallelic".equals(label)
                ? List.of(REF, ALT, ALT2) : List.of(REF, ALT);
        final List<GATKRead> reads = new ArrayList<>();
        final List<Integer> best = new ArrayList<>();
        for (int a = 0; a < composition.length; a++) {
            for (int i = 0; i < composition[a].length; i++) {
                reads.add(read("a" + a + "i" + i, composition[a][i]));
                best.add(a);
            }
        }
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", reads);
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(List.of("s1")), new IndexedAlleleList<>(alleles), bySample);
        final LikelihoodMatrix<GATKRead, Allele> m = likelihoods.sampleMatrix(0);
        for (int e = 0; e < reads.size(); e++) {
            for (int a = 0; a < alleles.size(); a++) {
                m.set(a, e, a == best.get(e) ? -1 : -10);
            }
        }
        return likelihoods;
    }

    static void asmq(final String label) {
        final VariantContext vc = mqContext(label);
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = mqLikelihoods(label);
        final AS_RMSMappingQuality annotation = new AS_RMSMappingQuality();
        try {
            emitMap("asmq", label, annotation.annotate(null, vc, likelihoods));
        } catch (final Exception | AssertionError e) {
            System.out.printf("asmq\t%s\tE:%s%n", label, e.getClass().getName());
        }
        try {
            emitMap("asmqraw", label, annotation.annotateRawData(null, vc, likelihoods));
        } catch (final Exception | AssertionError e) {
            System.out.printf("asmqraw\t%s\tE:%s%n", label, e.getClass().getName());
        }
    }

    static void asmqFinal(final String raw) {
        final int slots = raw.split("\\|", -1).length;
        final List<Allele> alleles = slots > 2 ? List.of(REF, ALT, ALT2) : List.of(REF, ALT);
        final GenotypeBuilder gb = new GenotypeBuilder("s1", List.of(REF, ALT));
        gb.AD(slots > 2 ? new int[] {2, 2, 2} : new int[] {3, 2});
        final VariantContext vc = new VariantContextBuilder().chr("chr1").start(START).stop(START)
                .alleles(alleles).genotypes(List.of(gb.make()))
                .attribute("AS_RAW_MQ", raw).make();
        try {
            emitMap("asmqfinal", raw, new AS_RMSMappingQuality().finalizeRawData(vc, vc));
        } catch (final Exception | AssertionError e) {
            System.out.printf("asmqfinal\t%s\tE:%s%n", raw, e.getClass().getName());
        }
    }

    // ---- AS_InbreedingCoeff ------------------------------------------------------------------

    static final int[] HOM_REF = {0, 60, 600};
    static final int[] HET = {60, 0, 60};
    static final int[] HOM_VAR = {600, 60, 0};

    static VariantContext icContext(final String label) {
        final boolean triallelic = "triallelic".equals(label);
        final List<Allele> alleles = triallelic ? List.of(REF, ALT, ALT2) : List.of(REF, ALT);
        final List<Genotype> genotypes = new ArrayList<>();
        switch (label) {
            case "ten-het":
                for (int i = 0; i < 10; i++) {
                    genotypes.add(new GenotypeBuilder("s" + i, List.of(REF, ALT)).PL(HET).make());
                }
                break;
            case "ten-hom-ref":
                for (int i = 0; i < 10; i++) {
                    genotypes.add(new GenotypeBuilder("s" + i, List.of(REF, REF))
                            .PL(HOM_REF).make());
                }
                break;
            case "nine-samples":
                for (int i = 0; i < 9; i++) {
                    genotypes.add(new GenotypeBuilder("s" + i, List.of(REF, ALT)).PL(HET).make());
                }
                break;
            case "twenty-mixed":
                for (int i = 0; i < 5; i++) {
                    genotypes.add(new GenotypeBuilder("r" + i, List.of(REF, REF))
                            .PL(HOM_REF).make());
                }
                for (int i = 0; i < 10; i++) {
                    genotypes.add(new GenotypeBuilder("h" + i, List.of(REF, ALT)).PL(HET).make());
                }
                for (int i = 0; i < 5; i++) {
                    genotypes.add(new GenotypeBuilder("v" + i, List.of(ALT, ALT))
                            .PL(HOM_VAR).make());
                }
                break;
            case "triallelic":
                for (int i = 0; i < 12; i++) {
                    genotypes.add(new GenotypeBuilder("s" + i, List.of(REF, ALT))
                            .PL(new int[] {60, 0, 60, 60, 60, 600}).make());
                }
                break;
            case "gq-only":
                for (int i = 0; i < 12; i++) {
                    genotypes.add(new GenotypeBuilder("s" + i, List.of(REF, ALT)).GQ(30).make());
                }
                break;
            case "no-genotypes":
                break;
            default: throw new IllegalArgumentException(label);
        }
        final VariantContextBuilder builder = new VariantContextBuilder().chr("chr1")
                .start(START).stop(START).alleles(alleles);
        return genotypes.isEmpty() ? builder.make() : builder.genotypes(genotypes).make();
    }

    static void asic(final String label) {
        final VariantContext vc = icContext(label);
        try {
            emitMap("asic", label, new AS_InbreedingCoeff().annotate(null, vc, null));
        } catch (final Exception | AssertionError e) {
            System.out.printf("asic\t%s\tE:%s%n", label, e.getClass().getName());
        }
    }

    static void hetCounts(final String label) {
        final VariantContext vc = icContext(label);
        try {
            final HeterozygosityCalculator calculator = new HeterozygosityCalculator(vc);
            final StringJoiner hets = new StringJoiner(";");
            for (final Allele a : vc.getAlternateAlleles()) {
                hets.add(a.getDisplayString() + "="
                        + Double.doubleToRawLongBits(calculator.getHetCount(a)));
            }
            final StringJoiner counts = new StringJoiner(";");
            for (final Allele a : vc.getAlleles()) {
                counts.add(a.getDisplayString() + "="
                        + Double.doubleToRawLongBits(calculator.getAlleleCount(a)));
            }
            System.out.printf("hetcounts\t%s\t%d\t%s\t%s%n", label, calculator.getSampleCount(),
                    hets, counts);
        } catch (final Exception | AssertionError e) {
            System.out.printf("hetcounts\t%s\tE:%s%n", label, e.getClass().getName());
        }
    }

    // ---- shared -------------------------------------------------------------------------------

    static void emitMap(final String kind, final String label, final Map<String, Object> result) {
        final StringJoiner joiner = new StringJoiner(";");
        if (result != null) {
            final List<String> keys = new ArrayList<>(result.keySet());
            keys.sort(String::compareTo);
            for (final String key : keys) {
                final Object value = result.get(key);
                joiner.add(String.format("%s=%s[%s]", key, value,
                        value == null ? "null" : value.getClass().getName()));
            }
        }
        System.out.printf("%s\t%s\t%s%n", kind, label, joiner);
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

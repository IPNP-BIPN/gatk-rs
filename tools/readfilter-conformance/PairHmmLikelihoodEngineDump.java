/*
 * PairHMMLikelihoodCalculationEngine, taken from the reference: every read against every
 * haplotype, the way HaplotypeCaller and Mutect2 score their reads.
 *
 * `computeReadLikelihoods` prepares each read (soft clips hard-clipped unless soft-clipped bases
 * are kept, the PCR indel error model capping the insertion and deletion qualities by the tandem
 * repeat around each base, base qualities capped by the mapping quality and squashed to Q6 under
 * the threshold, indel qualities squashed to Q6 under Q6), runs the Java LOGLESS_CACHING PairHMM
 * with a constant gap continuation penalty, normalizes each read's likelihoods to the best allele
 * (`normalizeLikelihoods`), and removes the reads no haplotype explains
 * (`filterPoorlyModeledEvidence`, a fixed or a dynamic threshold).
 *
 * Output:
 *
 *     pcr\t<rate factor>\t<repeat length>\t<adjusted qual>
 *     repeat\t<sequence>\t<offset>\t<unit>\t<length>
 *     case\t<label>\t<parameters>
 *     processed\t<label>\t<sample>\t<read>\t<bases>\t<quals>\t<ins quals>\t<del quals>
 *     evidence\t<label>\t<sample>\t<index>\t<read>\t<HMM quals length>
 *     lk\t<label>\t<sample>\t<allele>\t<evidence>\t<raw bits>\t<value>
 *     filtered\t<label>\t<sample>\t<read>
 *     error\t<label>\t<exception class>: <message>
 *
 * Likelihoods are printed as their raw IEEE bits: the normalization cap and the filtering
 * threshold are compared against them.
 *
 * Usage: PairHmmLikelihoodEngineDump
 */

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.PairHMMLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class PairHmmLikelihoodEngineDump {

    /** A read: quality strings are phred+33; a null indel string leaves the BI/BD tag off. */
    record Read(String sample, String name, int start, String bases, String quals, String cigar, int mapq,
                String ins, String del) {
    }

    record Params(byte gcp, double mismapping, PairHMMLikelihoodCalculationEngine.PCRErrorModel pcr, byte bqThreshold,
                  boolean dynamic, double scale, double errorRate, boolean symmetric, boolean disableCap,
                  boolean keepSoftClips) {
        String describe() {
            return String.format("gcp=%d\tmismapping=%s\tpcr=%s\tbq=%d\tdynamic=%b\tscale=%s\terror=%s\tsymmetric=%b\tdisablecap=%b\tkeepsoftclips=%b",
                    gcp, mismapping, pcr, bqThreshold, dynamic, scale, errorRate, symmetric, disableCap, keepSoftClips);
        }
    }

    record Case(String label, Params params, List<String> haplotypes, List<String> samples, List<Read> reads) {
    }

    static final String REF = "ACGTTGCATGTCAAAAAAAGATGCACACACAGAGCTCAGTCTAGGCTTAC";

    static final Params DEFAULT = new Params((byte) 10, -4.5, PairHMMLikelihoodCalculationEngine.PCRErrorModel.CONSERVATIVE,
            PairHMM.BASE_QUALITY_SCORE_THRESHOLD, false, 1.0, 0.02, true, false, false);

    static String q(final int length, final char c) {
        return String.valueOf(c).repeat(length);
    }

    public static void main(final String[] args) throws Exception {
        System.out.println("# PairHmmLikelihoodEngineDump: PairHMMLikelihoodCalculationEngine.computeReadLikelihoods");
        for (final PairHMMLikelihoodCalculationEngine.PCRErrorModel model : PairHMMLikelihoodCalculationEngine.PCRErrorModel.values()) {
            if (!model.hasRateFactor()) {
                continue;
            }
            final Method adjusted = PairHMMLikelihoodCalculationEngine.class
                    .getDeclaredMethod("getErrorModelAdjustedQual", int.class, double.class);
            adjusted.setAccessible(true);
            for (int i = 0; i <= ReadLikelihoodCalculationEngine.MAX_REPEAT_LENGTH; i++) {
                System.out.printf("pcr\t%s\t%d\t%d%n", model.getRateFactor(), i, (byte) adjusted.invoke(null, i, model.getRateFactor()));
            }
        }
        final Method repeats = ReadLikelihoodCalculationEngine.class.getDeclaredMethod("findTandemRepeatUnits", byte[].class, int.class);
        repeats.setAccessible(true);
        for (final String s : List.of(REF, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "ACACACACACACAC", "TTCTTCCCCAGT", "ACGT", "A",
                "AACAACAACAACGTGTGT", "GATATATCGCGCGCGCGCGA")) {
            final byte[] bytes = s.getBytes(StandardCharsets.US_ASCII);
            for (int offset = 0; offset < bytes.length; offset++) {
                @SuppressWarnings("unchecked") final org.apache.commons.lang3.tuple.Pair<byte[], Integer> p =
                        (org.apache.commons.lang3.tuple.Pair<byte[], Integer>) repeats.invoke(null, bytes, offset);
                System.out.printf("repeat\t%s\t%d\t%s\t%d%n", s, offset, new String(p.getLeft(), StandardCharsets.US_ASCII), p.getRight());
            }
        }

        final String snp = REF.substring(0, 25) + 'G' + REF.substring(26);
        final String homDel = REF.substring(0, 12) + REF.substring(13);
        final String caIns = REF.substring(0, 23) + "CA" + REF.substring(23);
        final List<String> haps = List.of(REF, snp, homDel, caIns);
        final List<Read> reads = List.of(
                new Read("s1", "exact", 6, REF.substring(5, 45), q(40, '?'), "40M", 60, null, null),
                new Read("s1", "snp", 6, snp.substring(5, 45), q(40, '?'), "40M", 60, null, null),
                new Read("s1", "lowqual", 6, REF.substring(5, 45), "?????#####((((((5555)))))+++++????????&&", "40M", 60, null, null),
                new Read("s1", "mapq20", 6, snp.substring(5, 45), q(40, 'I'), "40M", 20, null, null),
                new Read("s1", "softclip", 6, "GGGGG" + REF.substring(10, 40) + "TTTTT", q(40, '?'), "5S30M5S", 60, null, null),
                new Read("s1", "indelquals", 6, homDel.substring(5, 40), q(35, '?'), "35M", 60,
                        "####%%%%''''))))++++----////1111333", "DDDD::::0000&&&&!!!!######$$$$%%%&&"),
                new Read("s2", "garbage", 6, "TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGG", q(40, 'I'), "40M", 60, null, null),
                new Read("s2", "homdel", 6, homDel.substring(5, 45), q(40, '?'), "40M", 60, null, null),
                new Read("s2", "ins", 18, caIns.substring(17, 47), q(30, '?'), "30M", 60, null, null),
                new Read("s2", "short", 11, REF.substring(10, 20), q(10, '5'), "10M", 60, null, null),
                new Read("s2", "mapq0", 6, REF.substring(5, 45), q(40, '?'), "40M", 0, null, null),
                new Read("s2", "leftclip", 1, "CCCCCCC" + REF.substring(0, 25), q(32, '?'), "7S25M", 60, null, null));
        final List<String> both = List.of("s1", "s2");

        final List<Case> cases = new ArrayList<>(List.of(
                new Case("default", DEFAULT, haps, both, reads),
                new Case("pcr-none", with(DEFAULT, "pcr", PairHMMLikelihoodCalculationEngine.PCRErrorModel.NONE), haps, both, reads),
                new Case("pcr-hostile", with(DEFAULT, "pcr", PairHMMLikelihoodCalculationEngine.PCRErrorModel.HOSTILE), haps, both, reads),
                new Case("pcr-aggressive", with(DEFAULT, "pcr", PairHMMLikelihoodCalculationEngine.PCRErrorModel.AGGRESSIVE), haps, both, reads),
                new Case("dynamic-1", with(DEFAULT, "dynamic", true), haps, both, reads),
                new Case("dynamic-0.5", with(with(DEFAULT, "dynamic", true), "scale", 0.5), haps, both, reads),
                new Case("error-0.05", with(DEFAULT, "errorRate", 0.05), haps, both, reads),
                new Case("asymmetric", with(DEFAULT, "symmetric", false), haps, both, reads),
                new Case("disable-cap", with(DEFAULT, "disableCap", true), haps, both, reads),
                new Case("keep-softclips", with(DEFAULT, "keepSoftClips", true), haps, both, reads),
                new Case("no-mismapping-cap", with(DEFAULT, "mismapping", Double.NEGATIVE_INFINITY), haps, both, reads),
                new Case("mismapping-1", with(DEFAULT, "mismapping", -1.0), haps, both, reads),
                new Case("gcp-5", with(DEFAULT, "gcp", (byte) 5), haps, both, reads),
                new Case("bq-6", with(DEFAULT, "bqThreshold", (byte) 6), haps, both, reads),
                new Case("bq-25", with(DEFAULT, "bqThreshold", (byte) 25), haps, both, reads),
                new Case("one-haplotype", DEFAULT, List.of(REF), both, reads),
                new Case("ref-second", DEFAULT, List.of(snp, REF), both, reads),
                new Case("empty-sample", DEFAULT, haps, List.of("s1", "s3"), reads),
                new Case("negative-gcp", with(DEFAULT, "gcp", (byte) -1), haps, both, reads),
                new Case("positive-mismapping", with(DEFAULT, "mismapping", 0.5), haps, both, reads),
                new Case("bq-5", with(DEFAULT, "bqThreshold", (byte) 5), haps, both, reads)));
        for (final Case c : cases) {
            run(c);
        }
    }

    static Params with(final Params p, final String field, final Object value) {
        return new Params(
                field.equals("gcp") ? (byte) value : p.gcp(),
                field.equals("mismapping") ? (double) value : p.mismapping(),
                field.equals("pcr") ? (PairHMMLikelihoodCalculationEngine.PCRErrorModel) value : p.pcr(),
                field.equals("bqThreshold") ? (byte) value : p.bqThreshold(),
                field.equals("dynamic") ? (boolean) value : p.dynamic(),
                field.equals("scale") ? (double) value : p.scale(),
                field.equals("errorRate") ? (double) value : p.errorRate(),
                field.equals("symmetric") ? (boolean) value : p.symmetric(),
                field.equals("disableCap") ? (boolean) value : p.disableCap(),
                field.equals("keepSoftClips") ? (boolean) value : p.keepSoftClips());
    }

    static byte[] quals(final String s) {
        final byte[] q = s.getBytes(StandardCharsets.US_ASCII);
        for (int i = 0; i < q.length; i++) {
            q[i] -= 33;
        }
        return q;
    }

    static String show(final byte[] q) {
        final StringBuilder sb = new StringBuilder();
        for (final byte b : q) {
            if (sb.length() > 0) {
                sb.append(',');
            }
            sb.append(b);
        }
        return sb.toString();
    }

    @SuppressWarnings("unchecked")
    static void run(final Case c) throws Exception {
        System.out.printf("case\t%s\t%s%n", c.label(), c.params().describe());
        final Params p = c.params();
        final PairHMMLikelihoodCalculationEngine engine;
        try {
            engine = new PairHMMLikelihoodCalculationEngine(p.gcp(), null, null, PairHMM.Implementation.LOGLESS_CACHING,
                    null, p.mismapping(), p.pcr(), p.bqThreshold(), p.dynamic(), p.scale(), p.errorRate(), p.symmetric(),
                    p.disableCap(), p.keepSoftClips());
        } catch (final RuntimeException e) {
            System.out.printf("error\t%s\t%s: %s%n", c.label(), e.getClass().getSimpleName(), e.getMessage());
            return;
        }
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final Map<String, List<GATKRead>> perSample = new LinkedHashMap<>();
        for (final String s : c.samples()) {
            perSample.put(s, new ArrayList<>());
        }
        for (final Read r : c.reads()) {
            if (!perSample.containsKey(r.sample())) {
                continue;
            }
            final GATKRead read = ArtificialReadUtils.createArtificialRead(header, r.name(), 0, r.start(),
                    r.bases().getBytes(StandardCharsets.US_ASCII), quals(r.quals()), r.cigar());
            read.setMappingQuality(r.mapq());
            if (r.ins() != null) {
                ReadUtils.setInsertionBaseQualities(read, quals(r.ins()));
                ReadUtils.setDeletionBaseQualities(read, quals(r.del()));
            }
            perSample.get(r.sample()).add(read);
        }
        final List<Haplotype> haplotypes = new ArrayList<>();
        for (int i = 0; i < c.haplotypes().size(); i++) {
            final String h = c.haplotypes().get(i);
            haplotypes.add(new Haplotype(h.getBytes(StandardCharsets.US_ASCII), h.equals(REF)));
        }

        final Method modify = PairHMMLikelihoodCalculationEngine.class.getDeclaredMethod("modifyReadQualities", List.class);
        modify.setAccessible(true);
        for (final Map.Entry<String, List<GATKRead>> e : perSample.entrySet()) {
            for (final GATKRead processed : (List<GATKRead>) modify.invoke(engine, e.getValue())) {
                System.out.printf("processed\t%s\t%s\t%s\t%s\t%s\t%s\t%s%n", c.label(), e.getKey(), processed.getName(),
                        processed.getBasesString(), show(processed.getBaseQualities()),
                        show(ReadUtils.getBaseInsertionQualities(processed)), show(ReadUtils.getBaseDeletionQualities(processed)));
            }
        }

        final AlleleLikelihoods<GATKRead, Haplotype> result = engine.computeReadLikelihoods(haplotypes, header,
                new IndexedSampleList(c.samples()), perSample, true);
        for (int s = 0; s < result.numberOfSamples(); s++) {
            final String sample = result.getSample(s);
            final LikelihoodMatrix<GATKRead, Haplotype> matrix = result.sampleMatrix(s);
            for (int r = 0; r < matrix.evidenceCount(); r++) {
                final GATKRead read = matrix.evidence().get(r);
                final byte[] hmm = read.getTransientAttribute(PairHMMLikelihoodCalculationEngine.HMM_BASE_QUALITIES_TAG) == null
                        ? null : (byte[]) read.getTransientAttribute(PairHMMLikelihoodCalculationEngine.HMM_BASE_QUALITIES_TAG);
                System.out.printf("evidence\t%s\t%s\t%d\t%s\t%s%n", c.label(), sample, r, read.getName(),
                        hmm == null ? "-" : String.valueOf(hmm.length));
            }
            for (int a = 0; a < matrix.numberOfAlleles(); a++) {
                for (int r = 0; r < matrix.evidenceCount(); r++) {
                    final double v = matrix.get(a, r);
                    System.out.printf("lk\t%s\t%s\t%d\t%d\t%016x\t%s%n", c.label(), sample, a, r, Double.doubleToRawLongBits(v), v);
                }
            }
        }
        final Field filtered = AlleleLikelihoods.class.getDeclaredField("filteredEvidenceBySampleIndex");
        filtered.setAccessible(true);
        final List<List<GATKRead>> removed = (List<List<GATKRead>>) filtered.get(result);
        for (int s = 0; s < removed.size(); s++) {
            for (final GATKRead read : removed.get(s)) {
                System.out.printf("filtered\t%s\t%s\t%s%n", c.label(), result.getSample(s), read.getName());
            }
        }
        engine.close();
    }
}

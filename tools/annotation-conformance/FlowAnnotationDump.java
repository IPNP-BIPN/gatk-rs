/*
 * The eight flow-based annotations, taken from the reference.
 *
 * VARIANT_TYPE, INDEL_CLASSIFY, INDEL_LENGTH, HMER_INDEL_LENGTH, HMER_INDEL_NUC, LEFT/RIGHT_MOTIF,
 * GC_CONTENT and CYCLESKIP_STATUS describe a variant in FLOW SPACE rather than in base space.
 *
 *   - a flow key is a run-length encoding against a repeating flow order (TGCA by default), and
 *     N matches EVERY flow base, so an ambiguous base is absorbed into whatever run it sits in;
 *   - an hmer indel is one whose reference and alternate flow keys differ in exactly one flow;
 *     the length reported is the larger of the two, and the nucleotide is the order's base there;
 *   - VARIANT_TYPE asks three questions in order, and one ordinary indel beside an hmer one makes
 *     the whole site non-h-indel;
 *   - a motif that would run off the contig is not truncated: the annotation is DROPPED;
 *   - the left motif SHIFTS for an indel, dropping its first base and appending the reference's;
 *   - GC_CONTENT is a float, not a double.
 *
 * Output:
 *
 *     flow\t<label>\t<key>=<value>;...
 *     key\t<bases>\t<flow order>\t<key, comma-separated or E>
 *
 * Usage: FlowAnnotationDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.flow.CycleSkipStatus;
import org.broadinstitute.hellbender.tools.walkers.annotator.flow.GcContent;
import org.broadinstitute.hellbender.tools.walkers.annotator.flow.HmerIndelLength;
import org.broadinstitute.hellbender.tools.walkers.annotator.flow.HmerIndelNuc;
import org.broadinstitute.hellbender.tools.walkers.annotator.flow.HmerMotifs;
import org.broadinstitute.hellbender.tools.walkers.annotator.flow.IndelClassify;
import org.broadinstitute.hellbender.tools.walkers.annotator.flow.IndelLength;
import org.broadinstitute.hellbender.tools.walkers.annotator.flow.VariantType;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.FlowBasedKeyCodec;

import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

public class FlowAnnotationDump {

    static final SAMFileHeader HEADER = makeHeader();
    /** A window with a homopolymer run in it, so the hmer paths are reachable. */
    static final String BASES = "ACGTACGTACAAAAAGGGGTTTTCCCCACGTACGTACGT";
    static final int WINDOW_START = 90;

    public static void main(final String[] args) {
        System.out.println("# FlowAnnotationDump: the eight flow-based annotations");

        // The key codec on its own, including the two cases the port's tests name.
        for (final String bases : new String[] {
                "TAA", "TNA", "ACGT", "TGCA", "AAAA", "", "T", "TTTTGGGG", "NNNN", "ACGTN"}) {
            key(bases, "TGCA");
            key(bases, "ACGT");
        }
        key("X", "TGCA");

        // Sites: a SNP, an insertion, a deletion, an hmer indel, a mixed site, a spanning
        // deletion, and one at each edge of the window.
        site("snp", 100, "A", List.of("C"));
        site("insertion", 100, "A", List.of("AC"));
        site("deletion", 100, "AC", List.of("A"));
        site("hmer-insertion", 99, "A", List.of("AA"));
        site("hmer-deletion", 99, "AA", List.of("A"));
        site("mixed", 100, "A", List.of("C", "AC"));
        site("spanning-deletion", 100, "A", List.of("C", "*"));
        site("non-ref", 100, "A", List.of("C", "<NON_REF>"));
        site("mnp", 100, "AC", List.of("GT"));
        site("near-window-start", 91, "A", List.of("C"));
        site("near-window-end", 127, "A", List.of("C"));
        site("long-insertion", 105, "A", List.of("ACGTACGT"));
    }

    static SAMFileHeader makeHeader() {
        final SAMSequenceDictionary dictionary =
                new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", 1000)));
        return new SAMFileHeader(dictionary);
    }

    static void key(final String bases, final String flowOrder) {
        try {
            final int[] result = FlowBasedKeyCodec.baseArrayToKey(bases.getBytes(), flowOrder);
            final StringJoiner joiner = new StringJoiner(",");
            for (final int value : result) {
                joiner.add(Integer.toString(value));
            }
            System.out.printf("key\t%s\t%s\t%s%n", bases, flowOrder, joiner);
        } catch (final Exception | AssertionError e) {
            System.out.printf("key\t%s\t%s\tE%n", bases, flowOrder);
        }
    }

    static void site(final String label, final int start, final String ref,
                     final List<String> alts) {
        final VariantContextBuilder builder = new VariantContextBuilder().chr("chr1").start(start)
                .stop(start + ref.length() - 1);
        final java.util.List<Allele> alleles = new java.util.ArrayList<>();
        alleles.add(Allele.create(ref, true));
        for (final String alt : alts) {
            alleles.add(alt.equals("*") ? Allele.SPAN_DEL
                    : alt.equals("<NON_REF>") ? Allele.NON_REF_ALLELE : Allele.create(alt, false));
        }
        final VariantContext vc = builder.alleles(alleles).make();

        final SimpleInterval window =
                new SimpleInterval("chr1", WINDOW_START, WINDOW_START + BASES.length() - 1);
        final org.broadinstitute.hellbender.utils.reference.ReferenceBases bases =
                new org.broadinstitute.hellbender.utils.reference.ReferenceBases(
                        BASES.getBytes(), window);
        final ReferenceContext context = new ReferenceContext(
                new ReferenceMemorySource(bases, HEADER.getSequenceDictionary()), window);

        one("VariantType", label, new VariantType(), context, vc);
        one("IndelClassify", label, new IndelClassify(), context, vc);
        one("IndelLength", label, new IndelLength(), context, vc);
        one("HmerIndelLength", label, new HmerIndelLength(), context, vc);
        one("HmerIndelNuc", label, new HmerIndelNuc(), context, vc);
        one("HmerMotifs", label, new HmerMotifs(), context, vc);
        one("GcContent", label, new GcContent(), context, vc);
        one("CycleSkipStatus", label, new CycleSkipStatus(), context, vc);
    }

    static void one(final String name, final String label, final InfoFieldAnnotation annotation,
                    final ReferenceContext ref, final VariantContext vc) {
        try {
            final Map<String, Object> result = annotation.annotate(ref, vc, null);
            final StringJoiner joiner = new StringJoiner(";");
            for (final Map.Entry<String, Object> entry : result.entrySet()) {
                joiner.add(String.format("%s=%s", entry.getKey(), entry.getValue()));
            }
            System.out.printf("flow\t%s\t%s\t%s%n", name, label, joiner);
        } catch (final Exception | AssertionError e) {
            System.out.printf("flow\t%s\t%s\tE:%s%n", name, label, e.getClass().getName());
        }
    }
}

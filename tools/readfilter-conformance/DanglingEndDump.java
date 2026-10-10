/*
 * Dangling-end recovery on the read threading graph, taken from the reference.
 *
 * `ReadThreadingAssembler.getAssemblyResult` prunes low-weight chains, then (dangling branches are
 * recovered by default) calls `recoverDanglingTails` and `recoverDanglingHeads`, then
 * `removePathsNotConnectedToRef`. Recovery walks from each dangling end back to the reference,
 * aligns the branch against the reference with Smith-Waterman (the Java aligner and
 * `STANDARD_NGS`, as HaplotypeCaller passes them), and, when the CIGAR is simple enough, adds one
 * edge from the branch back into the reference. A dangling head whose matching bases run into its
 * source k-mer is first extended with new vertices that are not in the k-mer map.
 *
 * Two code paths for heads: `numDanglingMatchingPrefixBases` of -1 (the default) takes
 * `mergeDanglingHeadLegacy`, any value from 0 takes `mergeDanglingHead`, and the same value is the
 * minimum suffix a tail must match.
 *
 * Order: vertices and edges are printed in the order the graph iterates them, as in
 * ChainPrunerDump; a vertex's index is its position in the CURRENT vertex set.
 *
 * Output, per case and per stage (pruned, tails, heads, connected):
 *
 *     case\t<label>\tk=<k>\tprunefactor=<f>\tminbranch=<n>\trecoverall=<b>\tminmatching=<m>
 *     vertex\t<label>\t<stage>\t<index>\t<sequence>
 *     edge\t<label>\t<stage>\t<source>\t<target>\t<multiplicity>\t<pruning multiplicity>\t<ref>
 *     summary\t<label>\t<stage>\tvertices=<n>\tedges=<n>\tcycles=<b>\trefsource=<i>\trefsink=<i>
 *     error\t<label>\t<stage>\t<exception class>
 *
 * Usage: DanglingEndDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.LowWeightChainPruner;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.AbstractReadThreadingGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingGraph;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignmentConstants;

import java.lang.reflect.Constructor;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;

public class DanglingEndDump {

    record Read(String sample, String bases, int copies) {
    }

    record Case(String label, int k, int pruneFactor, int minBranch, boolean recoverAll, int minMatching,
                List<Read> reads) {
    }

    static final String REF = "ATGCGTACCTGAGTCAAGCTTGGACTCTAGAGCAATCGGTACCATTGACGTAGCTAGGCA";

    /** A base that is not the reference's. */
    static char alt(final char c) {
        return switch (c) {
            case 'A' -> 'C';
            case 'C' -> 'G';
            case 'G' -> 'T';
            default -> 'A';
        };
    }

    static String mutate(final String s, final int at) {
        return s.substring(0, at) + alt(s.charAt(at)) + s.substring(at + 1);
    }

    static Case tails(final String label, final List<Read> reads) {
        return new Case(label, 10, 2, 4, false, -1, reads);
    }

    static List<Read> one(final String bases) {
        return List.of(new Read("s1", bases, 3));
    }

    public static void main(final String[] args) throws Exception {
        System.out.println("# DanglingEndDump: recoverDanglingTails, recoverDanglingHeads, removePathsNotConnectedToRef");
        final String r = REF;
        // Tails: the read leaves the reference at 40 and ends before a k-mer rejoins it.
        final String tailSnp4 = r.substring(0, 40) + alt(r.charAt(40)) + r.substring(41, 45);
        final String tailSnp1 = r.substring(0, 40) + alt(r.charAt(40)) + r.substring(41, 42);
        final String tailSnp8 = r.substring(0, 40) + alt(r.charAt(40)) + r.substring(41, 49);
        final String tailIns = r.substring(0, 40) + "GA" + r.substring(40, 46);
        final String tailDel = r.substring(0, 40) + r.substring(43, 49);
        final String tailMnp = r.substring(0, 40) + alt(r.charAt(40)) + alt(r.charAt(41)) + r.substring(42, 47);
        final String tailJunk = r.substring(0, 40) + "TTTTTTTT";
        final String tailOther = r.substring(0, 40) + alt(r.charAt(40)) + "GGGG";
        final String tailNearStart = r.substring(0, 12) + alt(r.charAt(12)) + r.substring(13, 17);
        final String tailAtSource = r.substring(0, 10) + alt(r.charAt(10)) + r.substring(11, 16);
        // Heads: the read starts on the reference, mutated at offset p, and runs to the end.
        final String head2 = mutate(r.substring(20), 2);
        final String head4 = mutate(r.substring(20), 4);
        final String head6 = mutate(r.substring(20), 6);
        final String head9 = mutate(r.substring(20), 9);
        final String headIns = r.substring(20, 26) + "GA" + r.substring(26);
        final String headDel = r.substring(20, 26) + r.substring(29);
        final String headJunk = "TTTTTTTT" + r.substring(28);
        final String headLong = mutate(r.substring(5), 3);

        final List<Case> cases = List.of(
                new Case("ref-only", 10, 2, 4, false, -1, List.of()),
                tails("tail-snp-4", one(tailSnp4)),
                tails("tail-snp-1", one(tailSnp1)),
                new Case("tail-snp-1-minbranch-0", 10, 2, 0, false, -1, one(tailSnp1)),
                tails("tail-snp-8", one(tailSnp8)),
                tails("tail-insertion", one(tailIns)),
                tails("tail-deletion", one(tailDel)),
                tails("tail-mnp", one(tailMnp)),
                tails("tail-junk", one(tailJunk)),
                new Case("tail-junk-recover-all", 10, 2, 4, true, -1, one(tailJunk)),
                tails("tail-near-start", one(tailNearStart)),
                tails("tail-lca-is-ref-source", one(tailAtSource)),
                new Case("tail-snp-4-minmatching-3", 10, 2, 4, false, 3, one(tailSnp4)),
                new Case("tail-snp-4-minmatching-5", 10, 2, 4, false, 5, one(tailSnp4)),
                new Case("tail-deletion-minmatching-2", 10, 2, 4, false, 2, one(tailDel)),
                // A fork in the branch: given up on unless every branch is recovered.
                tails("tail-fork", List.of(new Read("s1", tailSnp4, 3), new Read("s1", tailOther, 3))),
                new Case("tail-fork-recover-all", 10, 2, 4, true, -1,
                        List.of(new Read("s1", tailSnp4, 3), new Read("s1", tailOther, 3))),
                // The last edges of the branch are under the prune factor, so the walk restarts.
                tails("tail-light-end", List.of(new Read("s1", tailSnp4, 3), new Read("s1", tailSnp8, 1))),
                new Case("tail-light-end-factor-0", 10, 0, 4, false, -1,
                        List.of(new Read("s1", tailSnp4, 3), new Read("s1", tailSnp8, 1))),
                tails("head-2", one(head2)),
                tails("head-4", one(head4)),
                tails("head-6", one(head6)),
                tails("head-9", one(head9)),
                tails("head-insertion", one(headIns)),
                tails("head-deletion", one(headDel)),
                tails("head-junk", one(headJunk)),
                tails("head-long", one(headLong)),
                new Case("head-2-minbranch-0", 10, 2, 0, false, -1, one(head2)),
                new Case("head-4-minmatching-3", 10, 2, 4, false, 3, one(head4)),
                new Case("head-6-minmatching-3", 10, 2, 4, false, 3, one(head6)),
                new Case("head-9-minmatching-3", 10, 2, 4, false, 3, one(head9)),
                new Case("head-9-minmatching-0", 10, 2, 4, false, 0, one(head9)),
                new Case("head-insertion-minmatching-3", 10, 2, 4, false, 3, one(headIns)),
                new Case("head-deletion-minmatching-3", 10, 2, 4, false, 3, one(headDel)),
                new Case("head-junk-recover-all", 10, 2, 4, true, -1, one(headJunk)),
                new Case("head-6-recover-all", 10, 2, 4, true, -1, one(head6)),
                // Both ends at once, and from two samples.
                tails("head-and-tail", List.of(new Read("s1", head6, 3), new Read("s1", tailSnp4, 3))),
                new Case("head-and-tail-minmatching-3", 10, 2, 4, false, 3,
                        List.of(new Read("s1", head6, 3), new Read("s1", tailSnp4, 3))),
                tails("two-samples", List.of(new Read("s1", head9, 2), new Read("s2", tailIns, 2))),
                new Case("k5-tail-snp-4", 5, 2, 4, false, -1, one(tailSnp4)),
                new Case("k5-head-6", 5, 2, 4, false, -1, one(head6)),
                new Case("negative-prune-factor", 10, -1, 4, false, -1, one(tailSnp4)),
                new Case("negative-min-branch", 10, 2, -1, false, -1, one(tailSnp4)));

        for (final Case c : cases) {
            run(c);
        }
    }

    static void run(final Case c) throws Exception {
        System.out.printf("case\t%s\tk=%d\tprunefactor=%d\tminbranch=%d\trecoverall=%b\tminmatching=%d%n", c.label(),
                c.k(), c.pruneFactor(), c.minBranch(), c.recoverAll(), c.minMatching());
        final Constructor<ReadThreadingGraph> ctor = ReadThreadingGraph.class.getDeclaredConstructor(
                int.class, boolean.class, byte.class, int.class, int.class);
        ctor.setAccessible(true);
        final ReadThreadingGraph graph = ctor.newInstance(c.k(), false, (byte) 10, 1, c.minMatching());
        graph.addSequence("ref", REF.getBytes(StandardCharsets.US_ASCII), true);
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final Method addRead = AbstractReadThreadingGraph.class
                .getDeclaredMethod("addRead", GATKRead.class, SAMFileHeader.class);
        addRead.setAccessible(true);
        int n = 0;
        for (final Read r : c.reads()) {
            if (header.getReadGroup(r.sample()) == null) {
                final SAMReadGroupRecord group = new SAMReadGroupRecord(r.sample());
                group.setSample(r.sample());
                header.addReadGroup(group);
            }
            final byte[] bases = r.bases().getBytes(StandardCharsets.US_ASCII);
            final byte[] quals = new byte[bases.length];
            java.util.Arrays.fill(quals, (byte) 30);
            for (int copy = 0; copy < r.copies(); copy++) {
                final GATKRead read = ArtificialReadUtils.createArtificialRead(bases, quals, bases.length + "M");
                read.setName("read" + n++);
                read.setReadGroup(r.sample());
                addRead.invoke(graph, read, header);
            }
        }
        graph.buildGraphIfNecessary();
        new LowWeightChainPruner<MultiDeBruijnVertex, MultiSampleEdge>(Math.max(c.pruneFactor(), 0))
                .pruneLowWeightChains(graph);
        print(c.label(), "pruned", graph);

        final SmithWatermanAligner aligner = SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.JAVA);
        try {
            graph.recoverDanglingTails(c.pruneFactor(), c.minBranch(), c.recoverAll(), aligner,
                    SmithWatermanAlignmentConstants.STANDARD_NGS);
            print(c.label(), "tails", graph);
            graph.recoverDanglingHeads(c.pruneFactor(), c.minBranch(), c.recoverAll(), aligner,
                    SmithWatermanAlignmentConstants.STANDARD_NGS);
            print(c.label(), "heads", graph);
            graph.removePathsNotConnectedToRef();
            print(c.label(), "connected", graph);
        } catch (final RuntimeException e) {
            // Messages carry vertices' identity hash codes; the class is what is compared.
            System.out.printf("error\t%s\t%s%n", c.label(), e.getClass().getSimpleName());
        }
    }

    static void print(final String label, final String stage, final ReadThreadingGraph graph) {
        final Map<MultiDeBruijnVertex, Integer> index = new IdentityHashMap<>();
        for (final MultiDeBruijnVertex v : graph.vertexSet()) {
            index.put(v, index.size());
        }
        for (final MultiDeBruijnVertex v : graph.vertexSet()) {
            System.out.printf("vertex\t%s\t%s\t%d\t%s%n", label, stage, index.get(v), v.getSequenceString());
        }
        for (final MultiSampleEdge e : graph.edgeSet()) {
            System.out.printf("edge\t%s\t%s\t%d\t%d\t%d\t%d\t%b%n", label, stage, index.get(graph.getEdgeSource(e)),
                    index.get(graph.getEdgeTarget(e)), e.getMultiplicity(), e.getPruningMultiplicity(), e.isRef());
        }
        final MultiDeBruijnVertex source = graph.getReferenceSourceVertex();
        final MultiDeBruijnVertex sink = graph.getReferenceSinkVertex();
        System.out.printf("summary\t%s\t%s\tvertices=%d\tedges=%d\tcycles=%b\trefsource=%s\trefsink=%s%n",
                label, stage, graph.vertexSet().size(), graph.edgeSet().size(), graph.hasCycles(),
                source == null ? "-" : String.valueOf(index.get(source)),
                sink == null ? "-" : String.valueOf(index.get(sink)));
    }
}

/*
 * The sequence graph and its simplification, taken from the reference.
 *
 * After pruning, `getAssemblyResult` removes the read threading graph's paths not connected to the
 * reference and turns it into a `SeqGraph`: one vertex per k-mer vertex, a source keeping its whole
 * k-mer and every other vertex only its last base. `cleanupSeqGraph` then zips linear chains,
 * removes orphans and vertices not connected to the reference in either direction, simplifies
 * (`MergeDiamonds`, `MergeTails`, `SplitCommonSuffices`, `MergeCommonSuffices`, zipping, until
 * nothing changes), removes paths not connected to the reference, simplifies again, and adds an
 * empty vertex when everything collapsed into one.
 *
 * Part one runs that pipeline on graphs threaded from reads and prints the graph after every
 * stage. Part two builds sequence graphs by hand, the way the reference's own unit tests do, to
 * reach each transform where reads alone would not: merged tails, split common suffixes, merged
 * shared sequences, diamonds whose middles are wholly prefix and suffix.
 *
 * ORDER. Sequence vertices compare by identity, JGraphT 1.1.0 keeps vertices, edges and each
 * vertex's edge lists in insertion order, and every transform walks those orders, so vertices and
 * edges are printed as the graph iterates them. A vertex's index is its position in the current
 * vertex set.
 *
 * Output:
 *
 *     case\t<label>
 *     vertex\t<label>\t<stage>\t<index>\t<sequence>
 *     edge\t<label>\t<stage>\t<source>\t<target>\t<multiplicity>\t<ref>
 *     summary\t<label>\t<stage>\tvertices=<n>\tedges=<n>\trefsource=<i>\trefsink=<i>
 *     status\t<label>\t<JUST_ASSEMBLED_REFERENCE|ASSEMBLED_SOME_VARIATION>
 *     result\t<label>\t<stage>\t<true|false>                     (what a transform answered)
 *     error\t<label>\t<stage>\t<exception class>: <message>
 *
 * Usage: SeqGraphDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.LowWeightChainPruner;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqVertex;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.AbstractReadThreadingGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingGraph;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.lang.reflect.Constructor;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;

public class SeqGraphDump {

    record Read(String bases, int copies) {
    }

    record Case(String label, int k, String ref, List<Read> reads) {
    }

    static final String REF = "ACGTTGCATGTCGCATGATGCATGAGAGCTCAGTCTAGGC";

    public static void main(final String[] args) throws Exception {
        System.out.println("# SeqGraphDump: toSequenceGraph, cleanupSeqGraph and SeqGraph.simplifyGraph");
        final String snp = REF.substring(0, 20) + 'T' + REF.substring(21);
        final String snp2 = REF.substring(0, 12) + 'A' + REF.substring(13);
        final String snpB = REF.substring(0, 20) + 'C' + REF.substring(21);
        final String ins = REF.substring(0, 20) + "GG" + REF.substring(20);
        final String longIns = REF.substring(0, 20) + "TTACCGGT" + REF.substring(20);
        final String del = REF.substring(0, 18) + REF.substring(21);
        final String mnp = REF.substring(0, 19) + "TT" + REF.substring(21);
        final String before = "TTCCAAGGTT" + REF.substring(0, 20);
        final List<Case> cases = List.of(
                new Case("ref-only", 10, REF, List.of()),
                new Case("snp", 10, REF, List.of(new Read(REF, 3), new Read(snp, 3))),
                new Case("two-snps-apart", 10, REF, List.of(new Read(snp, 3), new Read(snp2, 3))),
                new Case("triallelic-snp", 10, REF, List.of(new Read(snp, 3), new Read(snpB, 3))),
                new Case("insertion", 10, REF, List.of(new Read(ins, 3))),
                new Case("long-insertion", 10, REF, List.of(new Read(longIns, 3))),
                new Case("deletion", 10, REF, List.of(new Read(del, 3))),
                new Case("snp-and-deletion", 10, REF, List.of(new Read(snp, 3), new Read(del, 3))),
                new Case("mnp", 10, REF, List.of(new Read(mnp, 3))),
                new Case("prefix-before-reference", 10, REF, List.of(new Read(before, 3))),
                new Case("matching-reads-only", 10, REF, List.of(new Read(REF, 5))),
                new Case("k5-snp", 5, REF, List.of(new Read(snp, 3))),
                new Case("k25-snp", 25, REF, List.of(new Read(snp, 3))));
        for (final Case c : cases) {
            fromReads(c);
        }
        byHand();
    }

    static void fromReads(final Case c) throws Exception {
        System.out.printf("case\t%s%n", c.label());
        final Constructor<ReadThreadingGraph> ctor = ReadThreadingGraph.class.getDeclaredConstructor(
                int.class, boolean.class, byte.class, int.class, int.class);
        ctor.setAccessible(true);
        final ReadThreadingGraph rtg = ctor.newInstance(c.k(), false, (byte) 10, 1, -1);
        rtg.addSequence("ref", c.ref().getBytes(StandardCharsets.US_ASCII), true);
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final SAMReadGroupRecord group = new SAMReadGroupRecord("s1");
        group.setSample("s1");
        header.addReadGroup(group);
        final Method addRead = AbstractReadThreadingGraph.class
                .getDeclaredMethod("addRead", GATKRead.class, SAMFileHeader.class);
        addRead.setAccessible(true);
        int n = 0;
        for (final Read r : c.reads()) {
            final byte[] bases = r.bases().getBytes(StandardCharsets.US_ASCII);
            final byte[] quals = new byte[bases.length];
            java.util.Arrays.fill(quals, (byte) 30);
            for (int copy = 0; copy < r.copies(); copy++) {
                final GATKRead read = ArtificialReadUtils.createArtificialRead(bases, quals, bases.length + "M");
                read.setName("read" + n++);
                read.setReadGroup("s1");
                addRead.invoke(rtg, read, header);
            }
        }
        rtg.buildGraphIfNecessary();
        new LowWeightChainPruner<MultiDeBruijnVertex, MultiSampleEdge>(2).pruneLowWeightChains(rtg);
        rtg.removePathsNotConnectedToRef();

        final SeqGraph seq = rtg.toSequenceGraph();
        print(c.label(), "seq", seq);
        seq.cleanNonRefPaths();
        print(c.label(), "cleaned", seq);
        cleanup(c.label(), seq);
    }

    /** `ReadThreadingAssembler.cleanupSeqGraph`, stage by stage. */
    static void cleanup(final String label, final SeqGraph seq) {
        try {
            seq.zipLinearChains();
            print(label, "zipped", seq);
            seq.removeSingletonOrphanVertices();
            seq.removeVerticesNotConnectedToRefRegardlessOfEdgeDirection();
            print(label, "pruned", seq);
            seq.simplifyGraph();
            print(label, "merged", seq);
            if (seq.getReferenceSourceVertex() == null || seq.getReferenceSinkVertex() == null) {
                System.out.printf("status\t%s\tJUST_ASSEMBLED_REFERENCE%n", label);
                return;
            }
            seq.removePathsNotConnectedToRef();
            seq.simplifyGraph();
            if (seq.vertexSet().size() == 1) {
                final SeqVertex complete = seq.vertexSet().iterator().next();
                final SeqVertex dummy = new SeqVertex("");
                seq.addVertex(dummy);
                seq.addEdge(complete, dummy, new BaseEdge(true, 0));
            }
            print(label, "final", seq);
            System.out.printf("status\t%s\tASSEMBLED_SOME_VARIATION%n", label);
        } catch (final RuntimeException e) {
            System.out.printf("error\t%s\tcleanup\t%s: %s%n", label, e.getClass().getSimpleName(),
                    String.valueOf(e.getMessage()).replaceAll("_id_-?\\d+_", "_id_"));
        }
    }

    // ---------------------------------------------------------------------------------------------
    // Part two: graphs built by hand.

    static SeqVertex v(final SeqGraph g, final String seq) {
        final SeqVertex vertex = new SeqVertex(seq);
        g.addVertex(vertex);
        return vertex;
    }

    static void e(final SeqGraph g, final SeqVertex a, final SeqVertex b, final boolean ref, final int mult) {
        g.addEdge(a, b, new BaseEdge(ref, mult));
    }

    static void byHand() {
        final String tail = "GGGGCCCCAAAATTTT"; // sixteen bases, past MergeTails' minimum of ten

        // Three sinks sharing a long suffix under one top: MergeTails.
        SeqGraph g = new SeqGraph(10);
        SeqVertex top = v(g, "ACGTACGTAC");
        SeqVertex t1 = v(g, "A" + tail);
        SeqVertex t2 = v(g, "CC" + tail);
        SeqVertex t3 = v(g, "TTT" + tail);
        e(g, top, t1, true, 5);
        e(g, top, t2, false, 2);
        e(g, top, t3, false, 1);
        simplify("tails", g);

        // Sinks sharing only a short suffix: left alone.
        g = new SeqGraph(10);
        top = v(g, "ACGTACGTAC");
        t1 = v(g, "AGG");
        t2 = v(g, "CGG");
        e(g, top, t1, true, 5);
        e(g, top, t2, false, 2);
        simplify("short-tails", g);

        // A bottom whose incoming vertices share a suffix: SplitCommonSuffices.
        g = new SeqGraph(10);
        SeqVertex s1 = v(g, "ACGTTT");
        SeqVertex s2 = v(g, "CGATTT");
        SeqVertex m1 = v(g, "AAACCC");
        SeqVertex m2 = v(g, "GGACCC");
        SeqVertex bottom = v(g, "TTGCA");
        e(g, s1, m1, true, 4);
        e(g, s2, m2, false, 3);
        e(g, m1, bottom, true, 4);
        e(g, m2, bottom, false, 3);
        simplify("common-suffix", g);

        // Incoming vertices with the same sequence: SharedSequenceMerger.
        g = new SeqGraph(10);
        SeqVertex a = v(g, "ACGT");
        SeqVertex b = v(g, "TTGG");
        SeqVertex p1 = v(g, "CCA");
        SeqVertex p2 = v(g, "CCA");
        bottom = v(g, "GAGA");
        e(g, a, p1, true, 2);
        e(g, b, p2, false, 2);
        e(g, p1, bottom, true, 2);
        e(g, p2, bottom, false, 2);
        simplify("shared-sequence", g);

        // Two same-sequence predecessors fed by ONE source: the second edge to the merged vertex is
        // refused by the graph, so its multiplicity is lost.
        g = new SeqGraph(10);
        a = v(g, "ACGT");
        p1 = v(g, "CCA");
        p2 = v(g, "CCA");
        bottom = v(g, "GAGA");
        e(g, a, p1, true, 2);
        e(g, a, p2, false, 3);
        e(g, p1, bottom, true, 2);
        e(g, p2, bottom, false, 3);
        simplify("shared-sequence-one-source", g);

        // A diamond whose middles are wholly the shared prefix and suffix: a direct prefix->suffix edge.
        g = new SeqGraph(10);
        top = v(g, "ACGTA");
        SeqVertex x1 = v(g, "CCGG");
        SeqVertex x2 = v(g, "CCAGG");
        SeqVertex x3 = v(g, "CCTTGG");
        bottom = v(g, "TACGT");
        e(g, top, x1, true, 3);
        e(g, top, x2, false, 2);
        e(g, top, x3, false, 1);
        e(g, x1, bottom, true, 3);
        e(g, x2, bottom, false, 2);
        e(g, x3, bottom, false, 1);
        simplify("diamond-prefix-suffix", g);

        // A diamond with nothing shared: left as it is.
        g = new SeqGraph(10);
        top = v(g, "ACGTA");
        x1 = v(g, "C");
        x2 = v(g, "G");
        bottom = v(g, "TACGT");
        e(g, top, x1, true, 3);
        e(g, top, x2, false, 2);
        e(g, x1, bottom, true, 3);
        e(g, x2, bottom, false, 2);
        simplify("diamond-nothing-shared", g);

        // A chain crossing from reference to non-reference: zipping stops at the boundary.
        g = new SeqGraph(10);
        a = v(g, "AAA");
        b = v(g, "CCC");
        SeqVertex c = v(g, "GGG");
        SeqVertex d = v(g, "TTT");
        e(g, a, b, true, 1);
        e(g, b, c, true, 1);
        e(g, c, d, false, 1);
        simplify("zip-ref-boundary", g);

        // A self-loop, which no transform may zip through.
        g = new SeqGraph(10);
        a = v(g, "AAA");
        b = v(g, "CCC");
        c = v(g, "GGG");
        e(g, a, b, true, 1);
        e(g, b, b, false, 1);
        e(g, b, c, true, 1);
        simplify("self-loop", g);
    }

    static void simplify(final String label, final SeqGraph g) {
        System.out.printf("case\t%s%n", label);
        print(label, "built", g);
        try {
            System.out.printf("result\t%s\tzip\t%b%n", label, g.zipLinearChains());
            print(label, "zipped", g);
            g.simplifyGraph();
            print(label, "simplified", g);
        } catch (final RuntimeException e) {
            System.out.printf("error\t%s\tsimplify\t%s: %s%n", label, e.getClass().getSimpleName(),
                    String.valueOf(e.getMessage()).replaceAll("_id_-?\\d+_", "_id_"));
        }
    }

    static void print(final String label, final String stage, final SeqGraph g) {
        final Map<SeqVertex, Integer> index = new IdentityHashMap<>();
        for (final SeqVertex v : g.vertexSet()) {
            index.put(v, index.size());
            System.out.printf("vertex\t%s\t%s\t%d\t%s%n", label, stage, index.get(v), v.getSequenceString());
        }
        for (final BaseEdge e : g.edgeSet()) {
            System.out.printf("edge\t%s\t%s\t%d\t%d\t%d\t%b%n", label, stage, index.get(g.getEdgeSource(e)),
                    index.get(g.getEdgeTarget(e)), e.getMultiplicity(), e.isRef());
        }
        final SeqVertex source = g.getReferenceSourceVertex();
        final SeqVertex sink = g.getReferenceSinkVertex();
        System.out.printf("summary\t%s\t%s\tvertices=%d\tedges=%d\trefsource=%s\trefsink=%s%n", label, stage,
                g.vertexSet().size(), g.edgeSet().size(),
                source == null ? "-" : String.valueOf(index.get(source)),
                sink == null ? "-" : String.valueOf(index.get(sink)));
    }
}

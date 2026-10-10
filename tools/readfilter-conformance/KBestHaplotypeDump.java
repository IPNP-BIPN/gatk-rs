/*
 * The k best haplotypes of a sequence graph, taken from the reference.
 *
 * `GraphBasedKBestHaplotypeFinder` is Dijkstra's k-shortest-paths over the sequence graph: paths
 * leave the sources, each step scores `log10(edge multiplicity) - log10(total outgoing
 * multiplicity)`, a `java.util.PriorityQueue` orders them by score (highest first) and then by
 * bases (`BASES_COMPARATOR` reversed), a vertex is extended at most `k` times, and a path that
 * reaches a sink is a haplotype. A graph with a cycle is first cut by
 * `removeCyclesAndVerticesThatDontLeadToSinks`, whose depth-first walk never forgets a visited
 * vertex, so it cuts cross edges as well as back edges.
 *
 * Scores are printed as raw bits: a zero multiplicity makes them infinite, and two infinities
 * make NaN, which `Double.compare` puts first once reversed.
 *
 * Output:
 *
 *     case\t<label>
 *     vertex\t<label>\t<index>\t<sequence>                      (the graph searched)
 *     edge\t<label>\t<source>\t<target>\t<multiplicity>\t<ref>
 *     haplotype\t<label>\t<mode>\t<k>\t<rank>\t<bases>\t<score bits>,<score>\t<ref>\t<vertex,...>
 *     error\t<label>\t<mode>\t<k>\t<exception class>: <message>
 *
 * `mode` is `refpath` (the reference source and sink, as the assembler asks) or `all` (every
 * source and sink, the default constructor).
 *
 * Usage: KBestHaplotypeDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.GraphBasedKBestHaplotypeFinder;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.KBestHaplotype;
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
import java.util.ArrayList;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;

public class KBestHaplotypeDump {

    record Read(String bases, int copies) {
    }

    static final String REF = "ACGTTGCATGTCGCATGATGCATGAGAGCTCAGTCTAGGC";
    static final int[] KS = {1, 2, 3, 128};

    public static void main(final String[] args) throws Exception {
        System.out.println("# KBestHaplotypeDump: GraphBasedKBestHaplotypeFinder");
        final String snp = REF.substring(0, 20) + 'T' + REF.substring(21);
        final String snp2 = REF.substring(0, 12) + 'A' + REF.substring(13);
        final String snpB = REF.substring(0, 20) + 'C' + REF.substring(21);
        final String ins = REF.substring(0, 20) + "GG" + REF.substring(20);
        final String del = REF.substring(0, 18) + REF.substring(21);
        fromReads("ref-only", 10, List.of());
        fromReads("snp", 10, List.of(new Read(REF, 3), new Read(snp, 3)));
        fromReads("snp-minor", 10, List.of(new Read(REF, 9), new Read(snp, 2)));
        fromReads("two-snps-apart", 10, List.of(new Read(snp, 3), new Read(snp2, 4)));
        fromReads("triallelic-snp", 10, List.of(new Read(snp, 3), new Read(snpB, 3)));
        fromReads("insertion", 10, List.of(new Read(ins, 3)));
        fromReads("snp-and-deletion", 10, List.of(new Read(snp, 3), new Read(del, 5)));
        fromReads("k5-snp", 5, List.of(new Read(snp, 3)));
        byHand();
    }

    static void fromReads(final String label, final int k, final List<Read> reads) throws Exception {
        final Constructor<ReadThreadingGraph> ctor = ReadThreadingGraph.class.getDeclaredConstructor(
                int.class, boolean.class, byte.class, int.class, int.class);
        ctor.setAccessible(true);
        final ReadThreadingGraph rtg = ctor.newInstance(k, false, (byte) 10, 1, -1);
        rtg.addSequence("ref", REF.getBytes(StandardCharsets.US_ASCII), true);
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final SAMReadGroupRecord group = new SAMReadGroupRecord("s1");
        group.setSample("s1");
        header.addReadGroup(group);
        final Method addRead = AbstractReadThreadingGraph.class
                .getDeclaredMethod("addRead", GATKRead.class, SAMFileHeader.class);
        addRead.setAccessible(true);
        int n = 0;
        for (final Read r : reads) {
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
        seq.cleanNonRefPaths();
        seq.zipLinearChains();
        seq.removeSingletonOrphanVertices();
        seq.removeVerticesNotConnectedToRefRegardlessOfEdgeDirection();
        seq.simplifyGraph();
        seq.removePathsNotConnectedToRef();
        seq.simplifyGraph();
        if (seq.vertexSet().size() == 1) {
            final SeqVertex complete = seq.vertexSet().iterator().next();
            final SeqVertex dummy = new SeqVertex("");
            seq.addVertex(dummy);
            seq.addEdge(complete, dummy, new BaseEdge(true, 0));
        }
        search(label, seq);
    }

    static SeqVertex v(final SeqGraph g, final String seq) {
        final SeqVertex vertex = new SeqVertex(seq);
        g.addVertex(vertex);
        return vertex;
    }

    static void e(final SeqGraph g, final SeqVertex a, final SeqVertex b, final boolean ref, final int mult) {
        g.addEdge(a, b, new BaseEdge(ref, mult));
    }

    static void byHand() {
        // Two equally weighted branches: the tie is broken by the bases, reversed.
        SeqGraph g = new SeqGraph(10);
        SeqVertex top = v(g, "ACGT");
        SeqVertex x = v(g, "A");
        SeqVertex y = v(g, "C");
        SeqVertex bottom = v(g, "TTTT");
        e(g, top, x, true, 5);
        e(g, top, y, false, 5);
        e(g, x, bottom, true, 5);
        e(g, y, bottom, false, 5);
        search("tie-on-score", g);

        // A zero-multiplicity branch: a score of minus infinity.
        g = new SeqGraph(10);
        top = v(g, "ACGT");
        x = v(g, "A");
        y = v(g, "C");
        bottom = v(g, "TTTT");
        e(g, top, x, true, 5);
        e(g, top, y, false, 0);
        e(g, x, bottom, true, 5);
        e(g, y, bottom, false, 0);
        search("zero-multiplicity", g);

        // Every outgoing multiplicity zero: log10(0) - log10(0) is NaN.
        g = new SeqGraph(10);
        top = v(g, "ACGT");
        x = v(g, "A");
        y = v(g, "C");
        bottom = v(g, "TTTT");
        e(g, top, x, true, 0);
        e(g, top, y, false, 0);
        e(g, x, bottom, true, 0);
        e(g, y, bottom, false, 0);
        search("all-zero", g);

        // Three nested bubbles: more haplotypes than the smaller k values keep.
        g = new SeqGraph(10);
        final SeqVertex a = v(g, "AAAA");
        final SeqVertex b1 = v(g, "C");
        final SeqVertex b2 = v(g, "G");
        final SeqVertex m = v(g, "TTTT");
        final SeqVertex c1 = v(g, "A");
        final SeqVertex c2 = v(g, "C");
        final SeqVertex c3 = v(g, "G");
        final SeqVertex z = v(g, "GGGG");
        e(g, a, b1, true, 6);
        e(g, a, b2, false, 3);
        e(g, b1, m, true, 6);
        e(g, b2, m, false, 3);
        e(g, m, c1, true, 4);
        e(g, m, c2, false, 3);
        e(g, m, c3, false, 2);
        e(g, c1, z, true, 4);
        e(g, c2, z, false, 3);
        e(g, c3, z, false, 2);
        search("six-haplotypes", g);

        // A cycle, cut before the search; the cut also drops a cross edge.
        g = new SeqGraph(10);
        final SeqVertex s = v(g, "AAAA");
        final SeqVertex p = v(g, "CC");
        final SeqVertex q = v(g, "GG");
        final SeqVertex t = v(g, "TTTT");
        e(g, s, p, true, 4);
        e(g, p, q, true, 4);
        e(g, q, p, false, 1);
        e(g, s, q, false, 2);
        e(g, q, t, true, 4);
        search("cycle", g);

        // Two sources and two sinks for the default constructor.
        g = new SeqGraph(10);
        final SeqVertex s1 = v(g, "AAAA");
        final SeqVertex s2 = v(g, "CCCC");
        final SeqVertex mid = v(g, "GG");
        final SeqVertex t1 = v(g, "TTTT");
        final SeqVertex t2 = v(g, "TATA");
        e(g, s1, mid, true, 4);
        e(g, s2, mid, false, 2);
        e(g, mid, t1, true, 4);
        e(g, mid, t2, false, 2);
        search("two-sources-two-sinks", g);
    }

    static void search(final String label, final SeqGraph g) {
        System.out.printf("case\t%s%n", label);
        final Map<SeqVertex, Integer> index = new IdentityHashMap<>();
        for (final SeqVertex v : g.vertexSet()) {
            index.put(v, index.size());
            System.out.printf("vertex\t%s\t%d\t%s%n", label, index.get(v), v.getSequenceString());
        }
        for (final BaseEdge e : g.edgeSet()) {
            System.out.printf("edge\t%s\t%d\t%d\t%d\t%b%n", label, index.get(g.getEdgeSource(e)),
                    index.get(g.getEdgeTarget(e)), e.getMultiplicity(), e.isRef());
        }
        for (final String mode : new String[] {"refpath", "all"}) {
            for (final int k : KS) {
                try {
                    final GraphBasedKBestHaplotypeFinder<SeqVertex, BaseEdge> finder = mode.equals("refpath")
                            ? new GraphBasedKBestHaplotypeFinder<>(g, g.getReferenceSourceVertex(), g.getReferenceSinkVertex())
                            : new GraphBasedKBestHaplotypeFinder<>(g);
                    final List<KBestHaplotype<SeqVertex, BaseEdge>> found = finder.findBestHaplotypes(k);
                    int rank = 0;
                    for (final KBestHaplotype<SeqVertex, BaseEdge> h : found) {
                        final List<String> path = new ArrayList<>();
                        for (final SeqVertex v : h.getVertices()) {
                            path.add(String.valueOf(index.get(v)));
                        }
                        System.out.printf("haplotype\t%s\t%s\t%d\t%d\t%s\t%016x,%s\t%b\t%s%n", label, mode, k, rank++,
                                new String(h.getBases(), StandardCharsets.US_ASCII),
                                Double.doubleToRawLongBits(h.score()), Double.toString(h.score()), h.isReference(),
                                String.join(",", path));
                    }
                } catch (final RuntimeException e) {
                    System.out.printf("error\t%s\t%s\t%d\t%s: %s%n", label, mode, k, e.getClass().getSimpleName(),
                            String.valueOf(e.getMessage()).replaceAll("_id_-?\\d+_", "_id_"));
                }
            }
        }
    }
}

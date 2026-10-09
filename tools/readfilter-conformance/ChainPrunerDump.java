/*
 * Low-weight chain pruning on the read threading graph, and the removal of paths not connected to
 * the reference, taken from the reference.
 *
 * After `buildGraphIfNecessary`, `ReadThreadingAssembler.createGraph` prunes before it checks for
 * cycles: `LowWeightChainPruner` (HaplotypeCaller's default) splits the graph into maximal linear
 * chains, removes every chain whose edges are all non-reference with a pruning multiplicity under
 * the prune factor, then `removeSingletonOrphanVertices`. `getAssemblyResult` later calls
 * `removePathsNotConnectedToRef`. Both are here, each printed with the graph it leaves.
 *
 * `ReadThreadingGraph`'s orphan removal is its own: unlike `BaseGraph`'s it does not spare an
 * isolated reference source, and removing a vertex removes its k-mer from the k-mer map, so the
 * map is printed after each stage too.
 *
 * Order: vertices and edges are printed in the order the graph iterates them, as in
 * ReadThreadingGraphDump, and a vertex's index is its position in the CURRENT vertex set, so the
 * indices of one stage are not those of the last. Chains are printed in `findAllChains` order.
 *
 * Output, per case and per stage (built, pruned, connected):
 *
 *     case\t<label>\tk=<k>\tminbq=<q>\tpruning=<n>\tprunefactor=<f>
 *     chain\t<label>\t<vertex index,...>\t<multiplicities>\t<pruning multiplicities>\t<removed>
 *     vertex\t<label>\t<stage>\t<index>\t<sequence>\t<additional info>
 *     edge\t<label>\t<stage>\t<source>\t<target>\t<multiplicity>\t<pruning multiplicity>\t<ref>
 *     kmer\t<label>\t<stage>\t<kmer>\t<vertex index>
 *     summary\t<label>\t<stage>\tvertices=<n>\tedges=<n>\tcycles=<b>\trefsource=<i>\trefsink=<i>
 *     error\t<label>\t<stage>\t<exception class>: <message>
 *
 * Usage: ChainPrunerDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.Kmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.ChainPruner;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.LowWeightChainPruner;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.Path;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.AbstractReadThreadingGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingGraph;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;

public class ChainPrunerDump {

    record Read(String sample, String bases, int copies) {
    }

    record Case(String label, int k, int pruningSamples, int pruneFactor, String ref, List<Read> reads) {
    }

    static final String REF = "ACGTTGCATGTCGCATGATGCATGAGAGCTCAGTCTAGGC";

    public static void main(final String[] args) throws Exception {
        System.out.println("# ChainPrunerDump: LowWeightChainPruner and removePathsNotConnectedToRef");
        final String snp = REF.substring(0, 20) + 'T' + REF.substring(21);
        final String snp2 = REF.substring(0, 12) + 'A' + REF.substring(13);
        final String tail = REF.substring(0, 32) + "TTTTTTTT";
        final String head = "GGGGGGGG" + REF.substring(8);
        final String before = "TTCCAAGGTT" + REF.substring(0, 20);
        final List<Case> cases = List.of(
                // A weight-one bubble under a prune factor of two: removed, the reference intact.
                new Case("snp-once-factor-2", 10, 1, 2, REF, List.of(new Read("s1", snp, 1))),
                // The same bubble seen twice survives.
                new Case("snp-twice-factor-2", 10, 1, 2, REF, List.of(new Read("s1", snp, 2))),
                // A prune factor of zero removes nothing; one removes nothing either (no count is 0).
                new Case("snp-once-factor-0", 10, 1, 0, REF, List.of(new Read("s1", snp, 1))),
                new Case("snp-once-factor-1", 10, 1, 1, REF, List.of(new Read("s1", snp, 1))),
                // Two bubbles of different weights under a factor of three.
                new Case("two-bubbles-factor-3", 10, 1, 3, REF,
                        List.of(new Read("s1", snp, 3), new Read("s1", snp2, 2))),
                // A dangling tail and a dangling head of weight one.
                new Case("dangling-tail", 10, 1, 2, REF, List.of(new Read("s1", tail, 1))),
                new Case("dangling-head", 10, 1, 2, REF, List.of(new Read("s1", head, 1))),
                // Heavy enough to survive pruning, then cut by removePathsNotConnectedToRef.
                new Case("heavy-dangling-tail", 10, 1, 2, REF, List.of(new Read("s1", tail, 4))),
                new Case("heavy-prefix-before-reference", 10, 1, 2, REF, List.of(new Read("s1", before, 3))),
                // Pruning multiplicity across samples: each sample saw the bubble once.
                new Case("two-samples-two-pruning-samples", 10, 2, 2, REF,
                        List.of(new Read("s1", snp, 1), new Read("s2", snp, 1))),
                new Case("two-samples-one-pruning-sample", 10, 1, 2, REF,
                        List.of(new Read("s1", snp, 1), new Read("s2", snp, 1))),
                // A chain mixing weights: kept if any edge reaches the factor.
                new Case("mixed-weights", 10, 1, 2, REF,
                        List.of(new Read("s1", snp, 1), new Read("s1", snp.substring(0, 28), 1))),
                // The reference alone, and a read that closes a cycle.
                new Case("ref-only", 10, 1, 2, REF, List.of()),
                new Case("rearranged-read", 5, 1, 2, REF,
                        List.of(new Read("s1", REF.substring(25) + REF.substring(0, 15), 1))),
                new Case("rearranged-read-heavy", 5, 1, 2, REF,
                        List.of(new Read("s1", REF.substring(25) + REF.substring(0, 15), 3))),
                // k = 5 duplicates the reference's repeated k-mers.
                new Case("k5-snp-once", 5, 1, 2, REF, List.of(new Read("s1", snp, 1))));

        for (final Case c : cases) {
            run(c);
        }
    }

    static void run(final Case c) throws Exception {
        System.out.printf("case\t%s\tk=%d\tminbq=10\tpruning=%d\tprunefactor=%d%n", c.label(), c.k(),
                c.pruningSamples(), c.pruneFactor());
        final Constructor<ReadThreadingGraph> ctor = ReadThreadingGraph.class.getDeclaredConstructor(
                int.class, boolean.class, byte.class, int.class, int.class);
        ctor.setAccessible(true);
        final ReadThreadingGraph graph = ctor.newInstance(c.k(), false, (byte) 10, c.pruningSamples(), -1);
        graph.addSequence("ref", c.ref().getBytes(StandardCharsets.US_ASCII), true);
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
        print(c.label(), "built", graph);

        final LowWeightChainPruner<MultiDeBruijnVertex, MultiSampleEdge> pruner =
                new LowWeightChainPruner<>(c.pruneFactor());
        chains(c.label(), graph, pruner, c.pruneFactor());
        pruner.pruneLowWeightChains(graph);
        print(c.label(), "pruned", graph);

        try {
            graph.removePathsNotConnectedToRef();
            print(c.label(), "connected", graph);
        } catch (final RuntimeException e) {
            // The message lists vertices by `toString`, which carries an identity hash code; that
            // part is not reproducible, so it is masked.
            System.out.printf("error\t%s\tconnected\t%s: %s%n", c.label(), e.getClass().getSimpleName(),
                    String.valueOf(e.getMessage()).replaceAll("_id_-?\\d+_", "_id_"));
        }
    }

    @SuppressWarnings("unchecked")
    static void chains(final String label, final ReadThreadingGraph graph,
                       final LowWeightChainPruner<MultiDeBruijnVertex, MultiSampleEdge> pruner, final int factor)
            throws Exception {
        final Map<MultiDeBruijnVertex, Integer> index = index(graph);
        final Method find = ChainPruner.class.getDeclaredMethod("findAllChains", BaseGraph.class);
        find.setAccessible(true);
        final List<Path<MultiDeBruijnVertex, MultiSampleEdge>> chains =
                (List<Path<MultiDeBruijnVertex, MultiSampleEdge>>) find.invoke(pruner, graph);
        for (final Path<MultiDeBruijnVertex, MultiSampleEdge> chain : chains) {
            final List<String> vertices = new ArrayList<>();
            for (final MultiDeBruijnVertex v : chain.getVertices()) {
                vertices.add(String.valueOf(index.get(v)));
            }
            final List<String> multiplicities = new ArrayList<>();
            final List<String> pruning = new ArrayList<>();
            boolean removed = true;
            for (final MultiSampleEdge e : chain.getEdges()) {
                multiplicities.add(String.valueOf(e.getMultiplicity()));
                pruning.add(String.valueOf(e.getPruningMultiplicity()));
                removed &= e.getPruningMultiplicity() < factor && !e.isRef();
            }
            System.out.printf("chain\t%s\t%s\t%s\t%s\t%b%n", label, String.join(",", vertices),
                    String.join(",", multiplicities), String.join(",", pruning), removed);
        }
    }

    static Map<MultiDeBruijnVertex, Integer> index(final ReadThreadingGraph graph) {
        final Map<MultiDeBruijnVertex, Integer> index = new IdentityHashMap<>();
        for (final MultiDeBruijnVertex v : graph.vertexSet()) {
            index.put(v, index.size());
        }
        return index;
    }

    @SuppressWarnings("unchecked")
    static void print(final String label, final String stage, final ReadThreadingGraph graph) throws Exception {
        final Map<MultiDeBruijnVertex, Integer> index = index(graph);
        for (final MultiDeBruijnVertex v : graph.vertexSet()) {
            System.out.printf("vertex\t%s\t%s\t%d\t%s\t%s%n", label, stage, index.get(v), v.getSequenceString(),
                    v.getAdditionalInfo());
        }
        for (final MultiSampleEdge e : graph.edgeSet()) {
            System.out.printf("edge\t%s\t%s\t%d\t%d\t%d\t%d\t%b%n", label, stage, index.get(graph.getEdgeSource(e)),
                    index.get(graph.getEdgeTarget(e)), e.getMultiplicity(), e.getPruningMultiplicity(), e.isRef());
        }
        final Field kmers = AbstractReadThreadingGraph.class.getDeclaredField("kmerToVertexMap");
        kmers.setAccessible(true);
        for (final Map.Entry<Kmer, MultiDeBruijnVertex> e : ((Map<Kmer, MultiDeBruijnVertex>) kmers.get(graph)).entrySet()) {
            System.out.printf("kmer\t%s\t%s\t%s\t%s%n", label, stage,
                    new String(e.getKey().bases(), StandardCharsets.US_ASCII),
                    index.containsKey(e.getValue()) ? String.valueOf(index.get(e.getValue())) : "removed");
        }
        final MultiDeBruijnVertex source = graph.getReferenceSourceVertex();
        final MultiDeBruijnVertex sink = graph.getReferenceSinkVertex();
        System.out.printf("summary\t%s\t%s\tvertices=%d\tedges=%d\tcycles=%b\trefsource=%s\trefsink=%s%n",
                label, stage, graph.vertexSet().size(), graph.edgeSet().size(), graph.hasCycles(),
                source == null ? "-" : String.valueOf(index.get(source)),
                sink == null ? "-" : String.valueOf(index.get(sink)));
    }
}

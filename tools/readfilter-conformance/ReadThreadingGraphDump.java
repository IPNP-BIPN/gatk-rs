/*
 * The read threading graph, as `ReadThreadingGraph.buildGraphIfNecessary` leaves it, taken from the
 * reference.
 *
 * This is the first stage of the assembler under `HaplotypeCaller` and `Mutect2`: the reference
 * haplotype is threaded k-mer by k-mer, then each read's runs of usable bases, sample by sample.
 * A k-mer that repeats within any one sequence is non-unique and never a merge point; every other
 * k-mer is one vertex, and an edge counts the sequences that walked it. Pruning, dangling-end
 * recovery and the conversion to a sequence graph come after, and are not here.
 *
 * WHAT ORDER MEANS. The graph is a JGraphT `DefaultDirectedGraph` (jgrapht-core 1.1.0): vertices
 * and edges compare by identity, the vertex and edge sets keep insertion order, and so do each
 * vertex's incoming and outgoing edge lists. The threading loop takes the FIRST outgoing edge
 * whose target's suffix matches, so that order decides the graph, and the dump prints vertices and
 * edges in the order the graph iterates them rather than sorted.
 *
 * The constructor that sets the minimum base quality and the number of pruning samples,
 * `addRead`, `getNonUniqueKmers` and the per-sample `addSequence` are package-private or
 * protected, as the assembler's own tests find them, so they are reached by reflection.
 *
 * Output, per case:
 *
 *     case\t<label>\tk=<k>\tminbq=<q>\tpruning=<n>\texisting=<b>\tbranches=<b>
 *     vertex\t<label>\t<index>\t<sequence>\t<additional info>
 *     edge\t<label>\t<source index>\t<target index>\t<multiplicity>\t<pruning multiplicity>\t<ref>
 *     refpath\t<label>\t<index,index,...>
 *     kmer\t<label>\t<kmer>\t<vertex index>          (kmerToVertexMap, in its order)
 *     nonunique\t<label>\t<kmer,kmer,...>            (sorted: a HashSet, only ever asked contains)
 *     summary\t<label>\tvertices=<n>\tedges=<n>\tcycles=<b>\tlowquality=<b>\trefsource=<i>\trefsink=<i>
 *     error\t<label>\t<exception class>: <message>
 *
 * Usage: ReadThreadingGraphDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.Kmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.AbstractReadThreadingGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingGraph;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

public class ReadThreadingGraphDump {

    /** A read: its sample, its bases, and its qualities (`null` for thirty everywhere). */
    record Read(String sample, String bases, byte[] quals) {
    }

    record Case(String label, int k, byte minBaseQuality, int pruningSamples,
                boolean onlyAtExisting, boolean throughBranches, String ref, List<Read> reads) {
    }

    // Forty bases with no repeated 10-mer and three repeated 5-mers (CATGA, GCATG, TGCAT), so at
    // k = 5 the reference itself carries non-unique k-mers, each copy its own vertex.
    static final String REF = "ACGTTGCATGTCGCATGATGCATGAGAGCTCAGTCTAGGC";

    public static void main(final String[] args) throws Exception {
        System.out.println("# ReadThreadingGraphDump: ReadThreadingGraph after buildGraphIfNecessary");
        final String snp = REF.substring(0, 20) + 'T' + REF.substring(21);
        final String ins = REF.substring(0, 20) + "GG" + REF.substring(20);
        final String del = REF.substring(0, 18) + REF.substring(21);
        final String repeat = "ACGTACGTACGTACGTTGCATGTCGCATG";
        final List<Case> cases = List.of(
                // The reference alone: one chain, every edge a reference edge of multiplicity one.
                new Case("ref-only", 5, (byte) 10, 1, false, false, REF, List.of()),
                // Three reads identical to the reference add to its edges' counts.
                new Case("ref-and-matching-reads", 5, (byte) 10, 1, false, false, REF,
                        List.of(new Read("s1", REF, null), new Read("s1", REF.substring(5, 35), null),
                                new Read("s1", REF.substring(10), null))),
                // A substitution in the middle: a bubble that leaves and rejoins the reference.
                new Case("snp-bubble", 5, (byte) 10, 1, false, false, REF,
                        List.of(new Read("s1", snp, null), new Read("s1", snp, null))),
                new Case("insertion-bubble", 5, (byte) 10, 1, false, false, REF,
                        List.of(new Read("s1", ins, null))),
                new Case("deletion-bubble", 5, (byte) 10, 1, false, false, REF,
                        List.of(new Read("s1", del, null))),
                // An N splits the read into two runs; each run long enough is threaded on its own.
                new Case("read-split-at-n", 5, (byte) 10, 1, false, false, REF,
                        List.of(new Read("s1", REF.substring(0, 15) + 'N' + REF.substring(16), null))),
                // Two bases under the minimum quality split it the same way, and a run shorter
                // than k is dropped.
                new Case("read-split-at-low-quality", 5, (byte) 10, 1, false, false, REF,
                        List.of(new Read("s1", REF, quals(REF.length(), 30, new int[] {3, 20, 21})))),
                // The same read at a minimum quality of zero is not split at all.
                new Case("min-quality-zero", 5, (byte) 0, 1, false, false, REF,
                        List.of(new Read("s1", REF, quals(REF.length(), 30, new int[] {3, 20, 21})))),
                // A reference with a repeated 4-mer: each copy is its own vertex.
                new Case("repeat-in-reference", 4, (byte) 10, 1, false, false, repeat,
                        List.of(new Read("s1", repeat.substring(4), null))),
                // Most of the k-mers repeat: the graph is low quality.
                new Case("low-complexity", 3, (byte) 10, 1, false, false, "AAAAAAAAAACAAAAAAAAAA",
                        List.of(new Read("s1", "AAAAAAAACAAAAAAA", null))),
                // A read with no unique k-mer to start from is not threaded at all.
                new Case("read-without-a-start", 3, (byte) 10, 1, false, false, REF,
                        List.of(new Read("s1", "ATATATATATAT", null))),
                // A read that starts off the reference creates its own source, then merges.
                new Case("read-before-reference", 5, (byte) 10, 1, false, false, REF.substring(10),
                        List.of(new Read("s1", REF.substring(0, 25), null))),
                // The reference rearranged: the read goes from the end of the second half back to
                // the start of the first, which closes a cycle.
                new Case("rearranged-read", 5, (byte) 10, 1, false, false, REF,
                        List.of(new Read("s1", REF.substring(25) + REF.substring(0, 15), null))),
                // Starting only at existing vertices: the off-reference prefix is not threaded.
                new Case("start-only-at-existing", 5, (byte) 10, 1, true, false, REF.substring(10),
                        List.of(new Read("s1", REF.substring(0, 25), null))),
                // Two samples, two pruning samples: the pruning multiplicity is the per-sample
                // count kept in the priority queue, not the total.
                new Case("two-samples", 5, (byte) 10, 2, false, false, REF,
                        List.of(new Read("s1", snp, null), new Read("s1", snp, null),
                                new Read("s2", snp, null), new Read("s2", REF, null))),
                new Case("two-samples-one-pruning-sample", 5, (byte) 10, 1, false, false, REF,
                        List.of(new Read("s1", snp, null), new Read("s1", snp, null),
                                new Read("s2", snp, null), new Read("s2", REF, null))),
                // The backwards count increase stops at a branch unless told to go through it.
                new Case("counts-through-branches", 5, (byte) 10, 1, false, true, REF,
                        List.of(new Read("s1", snp, null), new Read("s1", snp.substring(12), null))),
                new Case("counts-stop-at-branches", 5, (byte) 10, 1, false, false, REF,
                        List.of(new Read("s1", snp, null), new Read("s1", snp.substring(12), null))),
                // A larger k, as the assembler's first attempt uses.
                new Case("k10-snp", 10, (byte) 10, 1, false, false, REF,
                        List.of(new Read("s1", snp, null), new Read("s1", REF, null))),
                new Case("k25-snp", 25, (byte) 10, 1, false, false, REF,
                        List.of(new Read("s1", snp, null))));

        for (final Case c : cases) {
            run(c);
        }
    }

    static byte[] quals(final int length, final int value, final int[] low) {
        final byte[] q = new byte[length];
        java.util.Arrays.fill(q, (byte) value);
        for (final int i : low) {
            q[i] = 2;
        }
        return q;
    }

    static void run(final Case c) throws Exception {
        System.out.printf("case\t%s\tk=%d\tminbq=%d\tpruning=%d\texisting=%b\tbranches=%b%n", c.label(),
                c.k(), c.minBaseQuality(), c.pruningSamples(), c.onlyAtExisting(), c.throughBranches());
        try {
            final Constructor<ReadThreadingGraph> ctor = ReadThreadingGraph.class.getDeclaredConstructor(
                    int.class, boolean.class, byte.class, int.class, int.class);
            ctor.setAccessible(true);
            final ReadThreadingGraph graph =
                    ctor.newInstance(c.k(), false, c.minBaseQuality(), c.pruningSamples(), -1);
            graph.setThreadingStartOnlyAtExistingVertex(c.onlyAtExisting());
            if (c.throughBranches()) {
                final Method branches = AbstractReadThreadingGraph.class
                        .getDeclaredMethod("setIncreaseCountsThroughBranches", boolean.class);
                branches.setAccessible(true);
                branches.invoke(graph, true);
            }

            // As `ReadThreadingAssembler.createGraph`: the reference first, then the reads.
            graph.addSequence("ref", c.ref().getBytes(StandardCharsets.US_ASCII), true);

            final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
            for (final Read r : c.reads()) {
                if (header.getReadGroup(r.sample()) == null) {
                    final SAMReadGroupRecord group = new SAMReadGroupRecord(r.sample());
                    group.setSample(r.sample());
                    header.addReadGroup(group);
                }
            }
            final Method addRead = AbstractReadThreadingGraph.class
                    .getDeclaredMethod("addRead", GATKRead.class, SAMFileHeader.class);
            addRead.setAccessible(true);
            int n = 0;
            for (final Read r : c.reads()) {
                final byte[] bases = r.bases().getBytes(StandardCharsets.US_ASCII);
                final byte[] quals = r.quals() != null ? r.quals() : quals(bases.length, 30, new int[0]);
                final GATKRead read =
                        ArtificialReadUtils.createArtificialRead(bases, quals, bases.length + "M");
                read.setName("read" + n++);
                read.setReadGroup(r.sample());
                addRead.invoke(graph, read, header);
            }

            graph.buildGraphIfNecessary();
            print(c.label(), graph);
        } catch (final InvocationTargetException e) {
            final Throwable cause = e.getCause();
            System.out.printf("error\t%s\t%s: %s%n", c.label(), cause.getClass().getSimpleName(), cause.getMessage());
        } catch (final RuntimeException e) {
            System.out.printf("error\t%s\t%s: %s%n", c.label(), e.getClass().getSimpleName(), e.getMessage());
        }
    }

    @SuppressWarnings("unchecked")
    static void print(final String label, final ReadThreadingGraph graph) throws Exception {
        final Map<MultiDeBruijnVertex, Integer> index = new IdentityHashMap<>();
        for (final MultiDeBruijnVertex v : graph.vertexSet()) {
            index.put(v, index.size());
            System.out.printf("vertex\t%s\t%d\t%s\t%s%n", label, index.get(v), v.getSequenceString(),
                    v.getAdditionalInfo());
        }
        for (final MultiSampleEdge e : graph.edgeSet()) {
            System.out.printf("edge\t%s\t%d\t%d\t%d\t%d\t%b%n", label, index.get(graph.getEdgeSource(e)),
                    index.get(graph.getEdgeTarget(e)), e.getMultiplicity(), e.getPruningMultiplicity(), e.isRef());
        }

        final Field refPath = AbstractReadThreadingGraph.class.getDeclaredField("referencePath");
        refPath.setAccessible(true);
        final List<MultiDeBruijnVertex> path = (List<MultiDeBruijnVertex>) refPath.get(graph);
        final List<String> steps = new ArrayList<>();
        if (path != null) {
            for (final MultiDeBruijnVertex v : path) {
                steps.add(String.valueOf(index.get(v)));
            }
        }
        System.out.printf("refpath\t%s\t%s%n", label, String.join(",", steps));

        final Field kmers = AbstractReadThreadingGraph.class.getDeclaredField("kmerToVertexMap");
        kmers.setAccessible(true);
        for (final Map.Entry<Kmer, MultiDeBruijnVertex> e : ((Map<Kmer, MultiDeBruijnVertex>) kmers.get(graph)).entrySet()) {
            System.out.printf("kmer\t%s\t%s\t%d%n", label, new String(e.getKey().bases(), StandardCharsets.US_ASCII),
                    index.get(e.getValue()));
        }

        final Method nonUniques = ReadThreadingGraph.class.getDeclaredMethod("getNonUniqueKmers");
        nonUniques.setAccessible(true);
        final Set<String> sorted = new TreeSet<>();
        for (final Kmer k : (Set<Kmer>) nonUniques.invoke(graph)) {
            sorted.add(new String(k.bases(), StandardCharsets.US_ASCII));
        }
        System.out.printf("nonunique\t%s\t%s%n", label, String.join(",", sorted));

        final MultiDeBruijnVertex source = graph.getReferenceSourceVertex();
        final MultiDeBruijnVertex sink = graph.getReferenceSinkVertex();
        System.out.printf("summary\t%s\tvertices=%d\tedges=%d\tcycles=%b\tlowquality=%b\trefsource=%s\trefsink=%s%n",
                label, graph.vertexSet().size(), graph.edgeSet().size(), graph.hasCycles(), graph.isLowQualityGraph(),
                source == null ? "-" : String.valueOf(index.get(source)),
                sink == null ? "-" : String.valueOf(index.get(sink)));
    }
}

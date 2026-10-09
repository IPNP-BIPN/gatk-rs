/*
 * Adaptive chain pruning on the read threading graph, and the log likelihood ratio it rests on,
 * taken from the reference.
 *
 * `AdaptiveChainPruner` is Mutect2's default pruner (HaplotypeCaller's with --adaptive-pruning).
 * Instead of a fixed prune factor it scores each chain's log odds of being real against the other
 * edges at its two ends, with `Mutect2Engine.logLikelihoodRatio`, grows a subgraph of good chains
 * from seeds through a priority queue, and removes the rest. It does that twice: once at the
 * initial error rate, to estimate the graph's own error rate, and once at that estimate.
 *
 * Doubles are printed as their raw bits as well as their text, because the log odds are compared
 * against thresholds and an ulp is enough to move a chain across one.
 *
 * Output:
 *
 *     llr\t<ref count>\t<alt count>\t<error probability>\t<bits>,<value>|exception:<class>
 *     case\t<label>\tk=<k>\tinitial=<p>\tthreshold=<t>\tseeding=<t>\tmaxvariants=<n>
 *     chain\t<label>\t<pass>\t<vertex index,...>\t<multiplicities>\t<left bits>\t<right bits>\t<error chain>
 *     errorrate\t<label>\t<bits>,<value>
 *     vertex / edge / summary \t<label>\tpruned\t...      (as ChainPrunerDump)
 *     error\t<label>\t<exception class>: <message>
 *
 * Usage: AdaptiveChainPrunerDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadThreadingAssemblerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.AdaptiveChainPruner;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.ChainPruner;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.Path;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.AbstractReadThreadingGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingGraph;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collection;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class AdaptiveChainPrunerDump {

    record Read(String sample, String bases, int copies) {
    }

    record Case(String label, int k, double initial, double threshold, double seeding, int maxVariants,
                String ref, List<Read> reads) {
    }

    static final String REF = "ACGTTGCATGTCGCATGATGCATGAGAGCTCAGTCTAGGC";
    static final double LOD = ReadThreadingAssemblerArgumentCollection.DEFAULT_PRUNING_LOG_ODDS_THRESHOLD;
    static final double SEED = ReadThreadingAssemblerArgumentCollection.DEFAULT_PRUNING_SEEDING_LOG_ODDS_THRESHOLD;

    public static void main(final String[] args) throws Exception {
        System.out.println("# AdaptiveChainPrunerDump: Mutect2Engine.logLikelihoodRatio and AdaptiveChainPruner");
        for (final int ref : new int[] {0, 1, 2, 5, 10, 30, 100}) {
            for (final int alt : new int[] {1, 2, 3, 10, 40}) {
                for (final double p : new double[] {0.001, 0.01, 0.05, 0.2, 0.5, 0.0, 1.0, Double.NaN, 1.5}) {
                    String shown;
                    try {
                        shown = render(Mutect2Engine.logLikelihoodRatio(ref, alt, p));
                    } catch (final RuntimeException e) {
                        shown = "exception:" + e.getClass().getSimpleName();
                    }
                    System.out.printf("llr\t%d\t%d\t%s\t%s%n", ref, alt, p, shown);
                }
            }
        }

        final String snp = REF.substring(0, 20) + 'T' + REF.substring(21);
        final String snp2 = REF.substring(0, 12) + 'A' + REF.substring(13);
        final String tail = REF.substring(0, 32) + "TTTTTTTT";
        final List<Case> cases = List.of(
                // Mutect2's defaults: a real bubble at a fifth of the depth survives.
                new Case("real-bubble", 10, 0.001, LOD, SEED, 100, REF,
                        List.of(new Read("s1", REF, 8), new Read("s1", snp, 2))),
                // One read in twenty: an error by the graph's own rate.
                new Case("error-bubble", 10, 0.001, LOD, SEED, 100, REF,
                        List.of(new Read("s1", REF, 20), new Read("s1", snp, 1))),
                new Case("two-bubbles", 10, 0.001, LOD, SEED, 100, REF,
                        List.of(new Read("s1", REF, 10), new Read("s1", snp, 5), new Read("s1", snp2, 1))),
                new Case("dangling-tail", 10, 0.001, LOD, SEED, 100, REF,
                        List.of(new Read("s1", REF, 10), new Read("s1", tail, 1))),
                new Case("heavy-dangling-tail", 10, 0.001, LOD, SEED, 100, REF,
                        List.of(new Read("s1", REF, 4), new Read("s1", tail, 4))),
                // No variant may stay unpruned.
                new Case("max-variants-zero", 10, 0.001, LOD, SEED, 0, REF,
                        List.of(new Read("s1", REF, 8), new Read("s1", snp, 4), new Read("s1", snp2, 4))),
                new Case("max-variants-one", 10, 0.001, LOD, SEED, 1, REF,
                        List.of(new Read("s1", REF, 8), new Read("s1", snp, 4), new Read("s1", snp2, 4))),
                // A threshold no chain reaches, and one every chain does.
                new Case("high-threshold", 10, 0.001, 50.0, 60.0, 100, REF,
                        List.of(new Read("s1", REF, 8), new Read("s1", snp, 3))),
                new Case("zero-threshold", 10, 0.001, 0.0, 0.0, 100, REF,
                        List.of(new Read("s1", REF, 8), new Read("s1", snp, 1))),
                new Case("initial-error-rate-high", 10, 0.2, LOD, SEED, 100, REF,
                        List.of(new Read("s1", REF, 8), new Read("s1", snp, 2))),
                new Case("ref-only", 10, 0.001, LOD, SEED, 100, REF, List.of()),
                new Case("k5-real-bubble", 5, 0.001, LOD, SEED, 100, REF,
                        List.of(new Read("s1", REF, 6), new Read("s1", snp, 3))));
        for (final Case c : cases) {
            run(c);
        }
    }

    static String render(final double value) {
        return String.format("%016x,%s", Double.doubleToRawLongBits(value), Double.toString(value));
    }

    static void run(final Case c) throws Exception {
        System.out.printf("case\t%s\tk=%d\tinitial=%s\tthreshold=%s\tseeding=%s\tmaxvariants=%d%n", c.label(), c.k(),
                c.initial(), c.threshold(), c.seeding(), c.maxVariants());
        final Constructor<ReadThreadingGraph> ctor = ReadThreadingGraph.class.getDeclaredConstructor(
                int.class, boolean.class, byte.class, int.class, int.class);
        ctor.setAccessible(true);
        final ReadThreadingGraph graph = ctor.newInstance(c.k(), false, (byte) 10, 1, -1);
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

        final AdaptiveChainPruner<MultiDeBruijnVertex, MultiSampleEdge> pruner =
                new AdaptiveChainPruner<>(c.initial(), c.threshold(), c.seeding(), c.maxVariants());
        try {
            report(c, graph, pruner);
            pruner.pruneLowWeightChains(graph);
            print(c.label(), graph);
        } catch (final InvocationTargetException e) {
            System.out.printf("error\t%s\t%s: %s%n", c.label(), e.getCause().getClass().getSimpleName(),
                    e.getCause().getMessage());
        } catch (final RuntimeException e) {
            System.out.printf("error\t%s\t%s: %s%n", c.label(), e.getClass().getSimpleName(), e.getMessage());
        }
    }

    /** The two passes `chainsToRemove` makes, chain by chain, before the pruner removes anything. */
    @SuppressWarnings("unchecked")
    static void report(final Case c, final ReadThreadingGraph graph,
                       final AdaptiveChainPruner<MultiDeBruijnVertex, MultiSampleEdge> pruner) throws Exception {
        final Method find = ChainPruner.class.getDeclaredMethod("findAllChains", BaseGraph.class);
        find.setAccessible(true);
        final List<Path<MultiDeBruijnVertex, MultiSampleEdge>> chains =
                (List<Path<MultiDeBruijnVertex, MultiSampleEdge>>) find.invoke(pruner, graph);
        if (chains.isEmpty()) {
            return;
        }
        final Method odds = AdaptiveChainPruner.class.getDeclaredMethod("chainLogOdds", Path.class, BaseGraph.class,
                double.class);
        odds.setAccessible(true);
        final Method likely = AdaptiveChainPruner.class.getDeclaredMethod("likelyErrorChains", List.class,
                BaseGraph.class, double.class);
        likely.setAccessible(true);

        final Collection<Path<MultiDeBruijnVertex, MultiSampleEdge>> first =
                (Collection<Path<MultiDeBruijnVertex, MultiSampleEdge>>) likely.invoke(pruner, chains, graph, c.initial());
        pass(c.label(), "initial", graph, chains, odds, pruner, c.initial(), first);

        final int errorCount = first.stream().mapToInt(p -> p.getLastEdge().getMultiplicity()).sum();
        final int totalBases = chains.stream()
                .mapToInt(p -> p.getEdges().stream().mapToInt(MultiSampleEdge::getMultiplicity).sum()).sum();
        final double errorRate = (double) errorCount / totalBases;
        System.out.printf("errorrate\t%s\t%s%n", c.label(), render(errorRate));

        final Collection<Path<MultiDeBruijnVertex, MultiSampleEdge>> second =
                (Collection<Path<MultiDeBruijnVertex, MultiSampleEdge>>) likely.invoke(pruner, chains, graph, errorRate);
        pass(c.label(), "estimated", graph, chains, odds, pruner, errorRate, second);
    }

    @SuppressWarnings("unchecked")
    static void pass(final String label, final String name, final ReadThreadingGraph graph,
                     final List<Path<MultiDeBruijnVertex, MultiSampleEdge>> chains, final Method odds,
                     final AdaptiveChainPruner<MultiDeBruijnVertex, MultiSampleEdge> pruner, final double rate,
                     final Collection<Path<MultiDeBruijnVertex, MultiSampleEdge>> errors) throws Exception {
        final Map<MultiDeBruijnVertex, Integer> index = index(graph);
        final Set<Path<MultiDeBruijnVertex, MultiSampleEdge>> errorSet =
                java.util.Collections.newSetFromMap(new IdentityHashMap<>());
        errorSet.addAll(errors);
        for (final Path<MultiDeBruijnVertex, MultiSampleEdge> chain : chains) {
            final List<String> vertices = new ArrayList<>();
            for (final MultiDeBruijnVertex v : chain.getVertices()) {
                vertices.add(String.valueOf(index.get(v)));
            }
            final List<String> multiplicities = new ArrayList<>();
            for (final MultiSampleEdge e : chain.getEdges()) {
                multiplicities.add(String.valueOf(e.getMultiplicity()));
            }
            final Pair<Double, Double> lr = (Pair<Double, Double>) odds.invoke(pruner, chain, graph, rate);
            System.out.printf("chain\t%s\t%s\t%s\t%s\t%s\t%s\t%b%n", label, name, String.join(",", vertices),
                    String.join(",", multiplicities), render(lr.getLeft()), render(lr.getRight()),
                    errorSet.contains(chain));
        }
    }

    static Map<MultiDeBruijnVertex, Integer> index(final ReadThreadingGraph graph) {
        final Map<MultiDeBruijnVertex, Integer> index = new IdentityHashMap<>();
        for (final MultiDeBruijnVertex v : graph.vertexSet()) {
            index.put(v, index.size());
        }
        return index;
    }

    static void print(final String label, final ReadThreadingGraph graph) {
        final Map<MultiDeBruijnVertex, Integer> index = index(graph);
        for (final MultiDeBruijnVertex v : graph.vertexSet()) {
            System.out.printf("vertex\t%s\tpruned\t%d\t%s%n", label, index.get(v), v.getSequenceString());
        }
        for (final MultiSampleEdge e : graph.edgeSet()) {
            System.out.printf("edge\t%s\tpruned\t%d\t%d\t%d\t%b%n", label, index.get(graph.getEdgeSource(e)),
                    index.get(graph.getEdgeTarget(e)), e.getMultiplicity(), e.isRef());
        }
        System.out.printf("summary\t%s\tpruned\tvertices=%d\tedges=%d%n", label, graph.vertexSet().size(),
                graph.edgeSet().size());
    }
}

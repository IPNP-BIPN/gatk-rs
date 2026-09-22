/*
 * The file PathSeqBuildReferenceTaxonomy writes, byte for byte, and the map orders that decide it.
 *
 * `PSTaxonomyDatabase` is a Kryo stream whose nodes, children and contig-to-taxon entries are
 * written in their maps' own iteration order. `kryo-stream` pins the encoding over a tree built by
 * hand; what it cannot pin is the order a map reaches after the tool has built it: sized at
 * construction, copied, trimmed three times by `retainNodes`, and iterated again. That order is
 * unspecified by `HashMap`'s contract, so it is measured here the way
 * docs/an-unspecified-order-that-reaches-the-output.md measures a HashSet's, and the port's table
 * is a hypothesis this golden either confirms or refutes (IPNP-BIPN/gatk-rs#1181).
 *
 * Two halves.
 *
 * PROBES build one map with the constructor the tool uses, insert, remove, and print the order it
 * iterates in. Each row carries its own inputs, so a port replays it without knowing the case:
 *
 *   - `new HashMap<>()` over integer keys, below and across both growth points, over keys that
 *     share low bits, and over keys whose HIGH bits decide the bucket once folded in;
 *   - `new HashMap<>(n)` for the sizes `retainNodes` and `PSPathogenReferenceTaxonProperties` ask
 *     for, including two buckets and one;
 *   - `new HashSet<>(collection)`, which is `PSTreeNode.copy`, followed by removals;
 *   - a map grown large and then emptied down to three keys, which is what the properties map
 *     is after `removeUnusedTaxIds`;
 *   - and the same shapes over contig names, which are the output map's keys.
 *
 * RUNS execute the tool on a taxonomy and a reference and print the database file as hex. The
 * small fixture is the one `pathseq-build-reference-taxonomy` reads back; the wide one is built to
 * push every table the file passes through past a growth point: twenty species under one genus, a
 * taxon holding three contigs, an output map of more than forty entries, a branch with no contigs
 * that is trimmed away, and a subtree hanging from a parent nodes.dmp never declares, which is
 * unreachable and removed.
 *
 * Output:
 *
 *     probe\t<label>\t<int|string>\t<constructor>\t<inserted, comma-separated>\t<removed>=<order>
 *     fixture\t<name>=<the whole file, escaped>
 *     contig\t<reference>\t<index>=<name as the dictionary holds it>\t<length>
 *     len\t<label>=<the file's length in bytes>
 *     stream\t<label>=<the file, in lower-case hex>
 *     error\t<label>\t<exception class>:<message>
 *
 * where <constructor> is `default`, `sized:<n>` or `copy`, the last meaning `new HashSet<>(list)`
 * over the inserted keys in that order.
 *
 * Usage: PathSeqTaxonomyKryoDump
 */

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.apache.commons.compress.archivers.tar.TarArchiveOutputStream;
import org.broadinstitute.hellbender.tools.spark.pathseq.PathSeqBuildReferenceTaxonomy;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.zip.GZIPOutputStream;

public class PathSeqTaxonomyKryoDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# PathSeqTaxonomyKryoDump: the taxonomy database's bytes and the orders behind them");
        probes();
        runs();
    }

    // ---------------------------------------------------------------------------------------
    // Probes.

    static void probes() {
        // The seven taxa of the small fixture, in the order names.dmp lists them.
        intProbe("fixture-taxa", "default", List.of(1, 2, 10239, 562, 11234, 9606, 40674), List.of());
        // Twelve keys stay in sixteen buckets; the thirteenth doubles, the twenty-fifth again.
        intProbe("twelve", "default", spread(12, 13, 5), List.of());
        intProbe("thirteen", "default", spread(13, 13, 5), List.of());
        intProbe("twenty-five", "default", spread(25, 13, 5), List.of());
        // Keys sixty-four apart share a bucket until the table passes sixty-four buckets.
        intProbe("shared-low-bits", "default", List.of(4096, 4160, 4224, 4288, 32, 96, 160, 7), List.of());
        // Above 65535 the high half is folded into the bucket.
        intProbe("high-bits", "default",
                List.of(131567, 2759, 1224, 1236, 91347, 65536, 65537, 196608, 2147483647, 1117, 543),
                List.of());
        intProbe("sized-one", "sized:1", List.of(5, 3, 1, 7), List.of());
        intProbe("sized-two", "sized:2", List.of(4, 8, 12, 16, 1), List.of());
        intProbe("sized-three", "sized:3", List.of(4, 8, 12, 16, 1, 2), List.of());
        intProbe("sized-seventeen", "sized:17", List.of(48, 16, 32, 1, 33, 17), List.of());
        intProbe("sized-forty", "sized:40", spread(40, 29, 3), List.of());
        // `PSTreeNode.copy` and then `removeChild`.
        final List<Integer> twenty = spread(20, 1009, 562);
        intProbe("copy-twenty", "copy", twenty, List.of());
        intProbe("copy-twenty-trimmed", "copy", twenty, twenty.subList(3, 18));
        intProbe("copy-three", "copy", List.of(40, 24, 8), List.of());
        // Grown to thirty, left with three: the table does not shrink.
        final List<Integer> thirty = spread(30, 7, 1);
        intProbe("grown-then-trimmed", "default", thirty, thirty.subList(3, 30));
        // Contig names, into two buckets and into the default.
        final List<String> six = IntStream.range(0, 6).mapToObj(i -> "ref|NC_S" + i + ".1|")
                .collect(Collectors.toList());
        stringProbe("names-sized-two", "sized:" + SVUtils.hashMapCapacity(1), six, List.of());
        final List<String> thirtyNames = IntStream.range(0, 30).mapToObj(i -> "taxid|" + (1000 + 37 * i) + "|")
                .collect(Collectors.toList());
        stringProbe("names-default", "default", thirtyNames, List.of());
        stringProbe("names-trimmed", "default", thirtyNames, thirtyNames.subList(2, 28));
    }

    /** `count` keys `step` apart from `start`, inserted in a scrambled but fixed order. */
    static List<Integer> spread(final int count, final int step, final int start) {
        final List<Integer> keys = new ArrayList<>();
        for (int i = 0; i < count; i++) {
            keys.add(start + step * ((i * 7) % count));
        }
        return keys;
    }

    static void intProbe(final String label, final String constructor, final List<Integer> inserted,
                         final List<Integer> removed) {
        probe(label, "int", constructor, inserted, removed);
    }

    static void stringProbe(final String label, final String constructor, final List<String> inserted,
                            final List<String> removed) {
        probe(label, "string", constructor, inserted, removed);
    }

    static <K> void probe(final String label, final String kind, final String constructor,
                          final List<K> inserted, final List<K> removed) {
        final Collection<K> keys;
        if (constructor.equals("copy")) {
            keys = new HashSet<>(inserted);
        } else {
            final Map<K, Boolean> map = constructor.equals("default")
                    ? new HashMap<>()
                    : new HashMap<>(Integer.parseInt(constructor.substring("sized:".length())));
            for (final K key : inserted) {
                map.put(key, Boolean.TRUE);
            }
            keys = map.keySet();
        }
        keys.removeAll(removed);
        System.out.printf("probe\t%s\t%s\t%s\t%s\t%s=%s%n", label, kind, constructor, joined(inserted),
                joined(removed), joined(keys));
    }

    static String joined(final Collection<?> values) {
        return values.stream().map(String::valueOf).collect(Collectors.joining(","));
    }

    // ---------------------------------------------------------------------------------------
    // Runs.

    static void runs() throws Exception {
        final Path dir = Path.of("pathseq-taxonomy-kryo-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        // The small fixture, exactly as pathseq-build-reference-taxonomy builds it.
        final Path smallTaxdump = dir.resolve("small-taxdump.tar.gz");
        writeTaxdump(smallTaxdump, PathSeqBuildReferenceTaxonomyDump.NAMES,
                PathSeqBuildReferenceTaxonomyDump.NODES);
        fixture("small-names.dmp", PathSeqBuildReferenceTaxonomyDump.NAMES);
        fixture("small-nodes.dmp", PathSeqBuildReferenceTaxonomyDump.NODES);
        final Path smallFasta = reference(dir, "small", List.of(
                PathSeqBuildReferenceTaxonomyDump.contig("ref|NC_VIRUS.1| a virus", 300),
                PathSeqBuildReferenceTaxonomyDump.contig("ref|NC_BACT.1| a bacterium", 1000),
                PathSeqBuildReferenceTaxonomyDump.contig("ref|NC_SHORT.1| a short bacterium", 100),
                PathSeqBuildReferenceTaxonomyDump.contig("taxid|562| named by its taxon", 800),
                PathSeqBuildReferenceTaxonomyDump.contig("ACC_PLAIN.1 named by its first word", 900),
                PathSeqBuildReferenceTaxonomyDump.contig("gi|9|ref|NC_BOTH.1|taxid|11234|", 700)));
        final String smallRefseq = String.join("\n",
                "11234\tsomething\tNC_VIRUS.1\tmore\tcolumns\tignored",
                "562\tsomething\tNC_BACT.1\tmore\tcolumns\tignored",
                "562\tsomething\tNC_SHORT.1\tmore\tcolumns\tignored") + "\n";
        final String smallGenbank = "a\tACC_PLAIN.1\tc\td\te\tf\t9606\th\n";
        final Path smallRefseqPath = catalog(dir, "small-refseq.catalog", smallRefseq);
        final Path smallGenbankPath = catalog(dir, "small-genbank.catalog", smallGenbank);
        run(dir, "small-both-catalogs", smallFasta, smallTaxdump, smallRefseqPath, smallGenbankPath, 0);
        run(dir, "small-min-length-500", smallFasta, smallTaxdump, smallRefseqPath, smallGenbankPath, 500);
        run(dir, "small-refseq-only", smallFasta, smallTaxdump, smallRefseqPath, null, 0);

        // The wide fixture.
        final StringBuilder names = new StringBuilder();
        final StringBuilder nodes = new StringBuilder();
        final List<PathSeqBuildReferenceTaxonomyDump.Contig> contigs = new ArrayList<>();
        final StringBuilder refseq = new StringBuilder();
        final StringBuilder genbank = new StringBuilder();

        taxon(names, nodes, 1, "root", 1, "no rank");
        taxon(names, nodes, 2, "Bacteria", 1, "superkingdom");
        taxon(names, nodes, 10239, "Viruses", 1, "superkingdom");
        taxon(names, nodes, 2157, "Archaea", 1, "superkingdom");
        // A branch no contig sits on, which the trim removes.
        taxon(names, nodes, 28890, "Euryarchaeota", 2157, "phylum");
        taxon(names, nodes, 2172, "Methanobrevibacter", 28890, "genus");

        // Twenty species under one genus: the genus's children cross sixteen buckets.
        taxon(names, nodes, 561, "Escherichia", 2, "genus");
        for (int k = 0; k < 20; k++) {
            final int id = 1000 + 37 * k * k + 11 * k;
            taxon(names, nodes, id, "Escherichia species " + k, 561, "species");
            contigs.add(PathSeqBuildReferenceTaxonomyDump.contig("taxid|" + id + "| species " + k, 300 + 50 * k));
            if (k % 3 == 0) {
                // Two more by accession, so this taxon's own map outgrows its two buckets.
                contigs.add(PathSeqBuildReferenceTaxonomyDump.contig("ref|NC_E" + k + "A.1| first", 450 + k));
                contigs.add(PathSeqBuildReferenceTaxonomyDump.contig("ref|NC_E" + k + "B.1| second", 520 + k));
                refseq.append(id).append("\tx\tNC_E").append(k).append("A.1\n");
                refseq.append(id).append("\tx\tNC_E").append(k).append("B.1\n");
            }
        }
        // Keys sixty-four apart, which share low bits in the node map.
        taxon(names, nodes, 1279, "Staphylococcus", 2, "genus");
        for (int k = 0; k < 5; k++) {
            final int id = 8192 + 64 * k;
            taxon(names, nodes, id, "Staphylococcus species " + k, 1279, "species");
            contigs.add(PathSeqBuildReferenceTaxonomyDump.contig("GB_S" + k + ".1 by GenBank", 700));
            genbank.append("a\tGB_S").append(k).append(".1\tc\td\te\tf\t").append(id).append("\n");
        }
        // Thirteen short viruses, which the length filter never drops.
        taxon(names, nodes, 11157, "Paramyxoviridae", 10239, "family");
        for (int k = 0; k < 13; k++) {
            final int id = 11234 + 97 * k;
            taxon(names, nodes, id, "Virus species " + k, 11157, "species");
            contigs.add(PathSeqBuildReferenceTaxonomyDump.contig("ref|NC_V" + k + ".1| virus", 100 + k));
            refseq.append(id).append("\tx\tNC_V").append(k).append(".1\n");
        }
        // A large id, whose high half decides its bucket.
        taxon(names, nodes, 2697049, "Severe acute respiratory syndrome coronavirus 2", 10239, "species");
        contigs.add(PathSeqBuildReferenceTaxonomyDump.contig("taxid|2697049| sars", 29903));
        // A subtree under a parent nodes.dmp never declares: unreachable, and removed with its contig.
        taxon(names, nodes, 777777, "Orphan", 888888, "species");
        contigs.add(PathSeqBuildReferenceTaxonomyDump.contig("taxid|777777| orphan", 900));

        final Path wideTaxdump = dir.resolve("wide-taxdump.tar.gz");
        writeTaxdump(wideTaxdump, names.toString(), nodes.toString());
        fixture("wide-names.dmp", names.toString());
        fixture("wide-nodes.dmp", nodes.toString());
        final Path wideFasta = reference(dir, "wide", contigs);
        final Path wideRefseq = catalog(dir, "wide-refseq.catalog", refseq.toString());
        final Path wideGenbank = catalog(dir, "wide-genbank.catalog", genbank.toString());
        run(dir, "wide-both-catalogs", wideFasta, wideTaxdump, wideRefseq, wideGenbank, 0);
        run(dir, "wide-min-length-500", wideFasta, wideTaxdump, wideRefseq, wideGenbank, 500);
        run(dir, "wide-refseq-only", wideFasta, wideTaxdump, wideRefseq, null, 0);
    }

    static void taxon(final StringBuilder names, final StringBuilder nodes, final int id,
                      final String name, final int parent, final String rank) {
        names.append(id).append("\t|\t").append(name).append("\t|\t\t|\tscientific name\t|\n");
        nodes.append(id).append("\t|\t").append(parent).append("\t|\t").append(rank).append("\t|\n");
    }

    static void fixture(final String name, final String text) {
        System.out.printf("fixture\t%s=%s%n", name, ReferenceQueryDump.escape(text));
    }

    static Path catalog(final Path dir, final String name, final String text) throws Exception {
        fixture(name, text);
        return PathSeqBuildReferenceTaxonomyDump.writeGz(dir, name + ".gz", text);
    }

    /** The reference, and its dictionary's names as the tool will read them. */
    static Path reference(final Path dir, final String label,
                          final List<PathSeqBuildReferenceTaxonomyDump.Contig> contigs) throws Exception {
        final Path fasta = PathSeqBuildReferenceTaxonomyDump.writeReference(dir, label, contigs);
        int index = 0;
        for (final SAMSequenceRecord record : ReferenceSequenceFileFactory
                .getReferenceSequenceFile(fasta).getSequenceDictionary().getSequences()) {
            System.out.printf("contig\t%s\t%d=%s\t%d%n", label, index++, record.getSequenceName(),
                    record.getSequenceLength());
        }
        return fasta;
    }

    static void writeTaxdump(final Path path, final String names, final String nodes) throws Exception {
        try (final TarArchiveOutputStream tar =
                     new TarArchiveOutputStream(new GZIPOutputStream(Files.newOutputStream(path)))) {
            PathSeqBuildReferenceTaxonomyDump.addEntry(tar, "names.dmp", names);
            PathSeqBuildReferenceTaxonomyDump.addEntry(tar, "nodes.dmp", nodes);
        }
    }

    static void run(final Path dir, final String label, final Path fasta, final Path taxdump,
                    final Path refseq, final Path genbank, final int minLength) {
        final Path out = dir.resolve(label + ".db");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-R", fasta.toString(),
                "-O", out.toString(),
                "--tax-dump", taxdump.toString(),
                "--min-non-virus-contig-length", Integer.toString(minLength)));
        if (refseq != null) {
            argv.addAll(Arrays.asList("--refseq-catalog", refseq.toString()));
        }
        if (genbank != null) {
            argv.addAll(Arrays.asList("--genbank-catalog", genbank.toString()));
        }
        try {
            new PathSeqBuildReferenceTaxonomy().instanceMain(argv.toArray(new String[0]));
            final byte[] bytes = Files.readAllBytes(out);
            System.out.printf("len\t%s=%d%n", label, bytes.length);
            final StringBuilder hex = new StringBuilder(bytes.length * 2);
            for (final byte value : bytes) {
                hex.append(String.format("%02x", value));
            }
            System.out.printf("stream\t%s=%s%n", label, hex);
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage()).replace(dir.toString(), "<dir>")));
        }
    }
}

/*
 * The bytes Kryo writes for the two PathSeq builders' outputs, from the reference.
 *
 * `PathSeqBuildKmers` and `PathSeqBuildReferenceTaxonomy` write their whole answer through
 * `PSKmerUtils.writeKryoObject`, which is `new Kryo()` and `kryo.writeObject(output, obj)` and
 * nothing else. A runner that cannot produce those bytes cannot be compared against the reference at
 * all (IPNP-BIPN/gatk-rs#1181), so the stream is measured here the way an unspecified iteration
 * order is measured: as an observable of the pinned container rather than as a specification.
 *
 * GATK pins Kryo at `strictly [4,5)`, so this is Kryo 4's encoding and not Kryo 5's.
 *
 * The cases climb a ladder: `Output`'s own primitives first, then each container on its own, then
 * the two objects a tool actually writes. A port that matches the first rung and not the second has
 * a defect in the container rather than in the encoding, which a whole-file comparison cannot say.
 *
 * Nine behaviours this is built to catch.
 *
 *   - `Output.writeInt` IS FOUR BYTES BIG-ENDIAN and never a varint, which is what the serializers
 *     here call, while `writeInt(value, true)` is the varint the containers do not use;
 *   - `Output.writeLong` IS EIGHT BYTES BIG-ENDIAN by the same rule;
 *   - `Output.writeString` IS KRYO'S OWN and not Java's: a length written as a varint, and an
 *     all-ASCII string written as its bytes with the TOP BIT SET on the last one;
 *   - A NULL STRING IS ONE BYTE, which is how the length distinguishes it from an empty one;
 *   - `kryo.writeObject` WRITES NO CLASS, because the reader is told the class, so a nested object
 *     contributes its serializer's bytes and nothing more;
 *   - `LongHopscotchSet` WRITES ITS WHOLE TABLE, capacity included, so the bytes depend on the
 *     legal size above the requested capacity and on FNV-1a placing each value;
 *   - `LargeLongHopscotchSet` IS A COUNT AND THEN ITS PARTITIONS, each a `LongHopscotchSet`;
 *   - `PSKmerSet` IS `kmerSize`, THE MASK AS A LONG, AND THE SET, in that order;
 *   - AND `PSTaxonomyDatabase` WRITES ITS MAP IN ITS OWN ITERATION ORDER, which is a `HashMap`'s:
 *     the keys are printed here in the order they were written, so the order is data rather than an
 *     assumption. `PSTree` holds its nodes in a `HashMap<Integer, PSTreeNode>` and has the same
 *     property over integer keys, and its order is reported beside the tree's bytes.
 *
 * Output:
 *
 *     note\t<label>=<what was constructed>
 *     len\t<label>=<the stream's length in bytes>
 *     stream\t<label>=<the stream, in lower-case hex>
 *     key\t<label>\t<position>=<the map key written there>
 *     node\t<label>\t<position>=<the tree node id written there>
 *     error\t<label>=<exception class>:<message>
 *
 * Usage: KryoStreamDump
 */

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Output;

import org.broadinstitute.hellbender.tools.spark.pathseq.PSKmerSet;
import org.broadinstitute.hellbender.tools.spark.pathseq.PSTaxonomyDatabase;
import org.broadinstitute.hellbender.tools.spark.pathseq.PSTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerShort;
import org.broadinstitute.hellbender.tools.spark.utils.LargeLongHopscotchSet;
import org.broadinstitute.hellbender.tools.spark.utils.LongHopscotchSet;

import java.util.LinkedHashMap;
import java.util.Map;

public class KryoStreamDump {

    /** The k-mers every set case holds, chosen to be positive and to spread over the table. */
    static final long[] KMERS = {1L, 2L, 3L, 17L, 1024L, 65535L, 1048577L, 123456789L};

    public static void main(final String[] args) {
        System.out.println("# KryoStreamDump: the bytes Kryo 4 writes for the PathSeq outputs");

        primitives();
        containers();
        objects();
    }

    /**
     * `Output`'s own encodings, which every serializer here is built out of.
     *
     * The values are chosen so that a port cannot pass by accident: zero and minus one distinguish
     * big-endian from little-endian, the maxima pin the width, and the strings cross the ASCII
     * boundary and the one-byte length.
     */
    static void primitives() {
        stream("int-zero", "Output.writeInt(0)", output -> output.writeInt(0));
        stream("int-one", "Output.writeInt(1)", output -> output.writeInt(1));
        stream("int-minus-one", "Output.writeInt(-1)", output -> output.writeInt(-1));
        stream("int-max", "Output.writeInt(Integer.MAX_VALUE)",
                output -> output.writeInt(Integer.MAX_VALUE));
        stream("int-min", "Output.writeInt(Integer.MIN_VALUE)",
                output -> output.writeInt(Integer.MIN_VALUE));
        // The varint form, which the containers here never ask for and a port must not confuse
        // with the fixed one.
        stream("varint-300", "Output.writeInt(300, true)", output -> output.writeInt(300, true));
        stream("varint-minus-300", "Output.writeInt(-300, true)",
                output -> output.writeInt(-300, true));

        stream("long-zero", "Output.writeLong(0)", output -> output.writeLong(0L));
        stream("long-minus-one", "Output.writeLong(-1)", output -> output.writeLong(-1L));
        stream("long-max", "Output.writeLong(Long.MAX_VALUE)",
                output -> output.writeLong(Long.MAX_VALUE));
        stream("long-kmer-mask", "Output.writeLong(1048575)", output -> output.writeLong(1048575L));

        stream("string-null", "Output.writeString(null)", output -> output.writeString(null));
        stream("string-empty", "Output.writeString(\"\")", output -> output.writeString(""));
        stream("string-one-char", "Output.writeString(\"A\")", output -> output.writeString("A"));
        stream("string-ascii", "Output.writeString(\"1\")", output -> output.writeString("1"));
        stream("string-number", "Output.writeString(String.valueOf(131567))",
                output -> output.writeString(String.valueOf(131567)));
        stream("string-long-ascii", "Output.writeString(a 40-character name)",
                output -> output.writeString("Homo sapiens neanderthalensis and more xx"));
        // Not ASCII, which is the branch the fast path does not take.
        stream("string-utf8", "Output.writeString(\"Crème\")", output -> output.writeString("Crème"));
    }

    /** Each container on its own, so a divergence names the container rather than the file. */
    static void containers() {
        final Kryo kryo = new Kryo();

        final LongHopscotchSet small = new LongHopscotchSet(8);
        for (final long value : KMERS) {
            small.add(value);
        }
        stream("hopscotch-eight", "LongHopscotchSet(8) holding " + KMERS.length + " k-mers",
                output -> kryo.writeObject(output, small));

        // The same values in a set asked for a larger capacity: the table is bigger and the bytes
        // are longer, which is what says the capacity is written rather than the size.
        final LongHopscotchSet large = new LongHopscotchSet(64);
        for (final long value : KMERS) {
            large.add(value);
        }
        stream("hopscotch-sixty-four", "LongHopscotchSet(64) holding the same k-mers",
                output -> kryo.writeObject(output, large));

        final LongHopscotchSet empty = new LongHopscotchSet(8);
        stream("hopscotch-empty", "LongHopscotchSet(8) holding nothing",
                output -> kryo.writeObject(output, empty));

        final LargeLongHopscotchSet partitioned = new LargeLongHopscotchSet(KMERS.length);
        for (final long value : KMERS) {
            partitioned.add(value);
        }
        stream("large-hopscotch", "LargeLongHopscotchSet(" + KMERS.length + ") holding them",
                output -> kryo.writeObject(output, partitioned));

        final PSTree tree = new PSTree(1);
        tree.addNode(2, "Bacteria", 1, 0L, "superkingdom");
        tree.addNode(3, "Escherichia coli", 2, 4641652L, "species");
        stream("tree-three-nodes", "PSTree(1) with Bacteria and E. coli under it",
                output -> kryo.writeObject(output, tree));
        reportNodeOrder("tree-three-nodes", tree);
    }

    /** The two objects a tool writes, whole, which is what a runner has to produce. */
    static void objects() {
        final Kryo kryo = new Kryo();

        final LargeLongHopscotchSet set = new LargeLongHopscotchSet(KMERS.length);
        for (final long value : KMERS) {
            set.add(value);
        }
        final PSKmerSet kmerSet = new PSKmerSet(set, 31, new SVKmerShort(1048575L));
        stream("kmer-set", "PSKmerSet(k=31, mask=1048575) over the same k-mers",
                output -> kryo.writeObject(output, kmerSet));

        final PSTree tree = new PSTree(1);
        tree.addNode(2, "Bacteria", 1, 0L, "superkingdom");
        tree.addNode(3, "Escherichia coli", 2, 4641652L, "species");

        // A `LinkedHashMap` here, so that what the golden pins is Kryo's encoding and not a
        // HashMap's bucket order. The tool builds a `HashMap`, and the order THAT produces is the
        // second half of the problem: the case below writes one and reports the order it took.
        final Map<String, Integer> ordered = new LinkedHashMap<>();
        ordered.put("NC_000913.3", 3);
        ordered.put("NC_002695.2", 3);
        ordered.put("NZ_CP009273.1", 2);
        final PSTaxonomyDatabase orderedDatabase = new PSTaxonomyDatabase(tree, ordered);
        stream("taxonomy-linked-map", "PSTaxonomyDatabase over a LinkedHashMap of three accessions",
                output -> kryo.writeObject(output, orderedDatabase));
        reportOrder("taxonomy-linked-map", ordered);

        final Map<String, Integer> hashed = new java.util.HashMap<>();
        hashed.put("NC_000913.3", 3);
        hashed.put("NC_002695.2", 3);
        hashed.put("NZ_CP009273.1", 2);
        final PSTaxonomyDatabase hashedDatabase = new PSTaxonomyDatabase(tree, hashed);
        stream("taxonomy-hash-map", "PSTaxonomyDatabase over a HashMap of the same three",
                output -> kryo.writeObject(output, hashedDatabase));
        reportOrder("taxonomy-hash-map", hashed);
    }

    /** The order the tree's own map iterates, which is the order its nodes reached the bytes. */
    static void reportNodeOrder(final String label, final PSTree tree) {
        int position = 0;
        for (final int id : tree.getNodeIDs()) {
            System.out.printf("node\t%s\t%d=%d%n", label, position++, id);
        }
    }

    /** The order the map handed the serializer, which is the order its keys reached the bytes. */
    static void reportOrder(final String label, final Map<String, Integer> map) {
        int position = 0;
        for (final String key : map.keySet()) {
            System.out.printf("key\t%s\t%d=%s%n", label, position++, key);
        }
    }

    /** What one write leaves in the stream, as hex, with its length beside it. */
    static void stream(final String label, final String note, final Written written) {
        System.out.printf("note\t%s=%s%n", label, note);
        final Output output = new Output(1024, -1);
        try {
            written.write(output);
            output.flush();
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s=%s:%s%n", label, e.getClass().getName(),
                    String.valueOf(e.getMessage()).replace("\n", "\\n").replace("\t", "\\t"));
            return;
        }
        final byte[] bytes = output.toBytes();
        System.out.printf("len\t%s=%d%n", label, bytes.length);
        final StringBuilder hex = new StringBuilder(bytes.length * 2);
        for (final byte value : bytes) {
            hex.append(String.format("%02x", value));
        }
        System.out.printf("stream\t%s=%s%n", label, hex);
    }

    /** One write against an `Output`, so each case is a line rather than a block. */
    interface Written {
        void write(Output output) throws Exception;
    }
}

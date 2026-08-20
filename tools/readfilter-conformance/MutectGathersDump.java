/*
 * The three gathers a scattered Mutect run ends with, taken from the reference.
 *
 * `MergeMutectStats`, `GatherPileupSummaries` and `GatherNormalArtifactData` each take the shards
 * of one scattered run and write one file. Their record formats are already ported; what is
 * measured here is what each does with several files, including the empty ones and the ones out of
 * order.
 *
 * Eight behaviours this is built to catch.
 *
 *   - `MergeMutectStats` SUMS EVERY STATISTIC IT KNOWS AND REFUSES EVERY OTHER, its aggregation
 *     map holding `callable` alone, so a shard carrying any other statistic ends the run;
 *   - THE SUM IS OVER DOUBLES AND WRITTEN AS ONE, so three integer counts come back with a decimal
 *     point;
 *   - THE STATS FILES ARE A `Set`, so the same file given twice is read once;
 *   - `GatherPileupSummaries` SORTS ITS FILES BY THEIR FIRST RECORD against the sequence
 *     dictionary, not by the order they were given, and does not sort within a file;
 *   - IT DROPS A FILE WITH NO RECORDS before sorting, so an empty shard is not a first record to
 *     compare;
 *   - AND IT REFUSES A FILE WHOSE SAMPLE DIFFERS, the writer taking the sample name from the first
 *     file it saw;
 *   - `GatherNormalArtifactData` CONCATENATES IN THE ORDER GIVEN, with no sort and no empty check;
 *   - AND EACH TOOL WRITES ITS OWN HEADER, so the output of a gather is a file of the same shape as
 *     its inputs rather than a concatenation of them.
 *
 * Output:
 *
 *     stats\t<label>=<the whole stats file, escaped>
 *     pileup\t<label>=<the whole pileup-summaries file, escaped>
 *     artifact\t<label>=<the whole normal-artifact file, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: MutectGathersDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.walkers.contamination.GatherPileupSummaries;
import org.broadinstitute.hellbender.tools.walkers.mutect.MergeMutectStats;
import org.broadinstitute.hellbender.tools.walkers.mutect.GatherNormalArtifactData;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class MutectGathersDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("mutect-gathers-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, PreprocessIntervalsDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        final Path dict = dir.resolve("ref.dict");
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dict});

        System.out.println("# MutectGathersDump: the three gathers a scattered Mutect run ends with");

        // MergeMutectStats over three shards, one of them zero.
        stats("three-shards", dir,
                "statistic\tvalue\ncallable\t100.0\n",
                "statistic\tvalue\ncallable\t250.0\n",
                "statistic\tvalue\ncallable\t0.0\n");
        // Integer-looking counts, which come back with a decimal point.
        stats("integers", dir, "statistic\tvalue\ncallable\t1\n", "statistic\tvalue\ncallable\t2\n");
        // One shard only.
        stats("one-shard", dir, "statistic\tvalue\ncallable\t7.5\n");
        // A shard with no rows at all.
        stats("empty-shard", dir, "statistic\tvalue\ncallable\t7.5\n", "statistic\tvalue\n");
        // A statistic the aggregation map does not know.
        stats("unknown-statistic", dir, "statistic\tvalue\ncallable\t1.0\nother\t2.0\n");

        // GatherPileupSummaries over shards given out of order, one of them empty.
        pileup("out-of-order", dir, dict,
                "#<METADATA>SAMPLE=sample\n"
                + "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n"
                + "chr2\t10\t10\t5\t0\t0.5\n",
                "#<METADATA>SAMPLE=sample\n"
                + "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n"
                + "chr1\t10\t20\t2\t0\t0.1\n"
                + "chr1\t50\t18\t4\t1\t0.2\n");
        // An empty shard, which is dropped before the sort.
        pileup("empty-shard", dir, dict,
                "#<METADATA>SAMPLE=sample\n"
                + "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n",
                "#<METADATA>SAMPLE=sample\n"
                + "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n"
                + "chr1\t10\t20\t2\t0\t0.1\n");
        // Two samples, which the writer refuses.
        pileup("two-samples", dir, dict,
                "#<METADATA>SAMPLE=first\n"
                + "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n"
                + "chr1\t10\t20\t2\t0\t0.1\n",
                "#<METADATA>SAMPLE=second\n"
                + "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n"
                + "chr1\t50\t18\t4\t1\t0.2\n");

        // GatherNormalArtifactData, which concatenates in the order given.
        artifact("two-shards", dir,
                "normal_alt\tnormal_dp\ttumor_alt\ttumor_dp\tdownsampling\ttype\n"
                + "1\t20\t5\t30\t1.0\tSNV\n",
                "normal_alt\tnormal_dp\ttumor_alt\ttumor_dp\tdownsampling\ttype\n"
                + "0\t25\t7\t35\t1.0\tINDEL\n");
        // The same two the other way round, to show the order is the caller's.
        artifact("reversed", dir,
                "normal_alt\tnormal_dp\ttumor_alt\ttumor_dp\tdownsampling\ttype\n"
                + "0\t25\t7\t35\t1.0\tINDEL\n",
                "normal_alt\tnormal_dp\ttumor_alt\ttumor_dp\tdownsampling\ttype\n"
                + "1\t20\t5\t30\t1.0\tSNV\n");
        // An empty shard, which is written through rather than dropped.
        artifact("empty-shard", dir,
                "normal_alt\tnormal_dp\ttumor_alt\ttumor_dp\tdownsampling\ttype\n",
                "normal_alt\tnormal_dp\ttumor_alt\ttumor_dp\tdownsampling\ttype\n"
                + "1\t20\t5\t30\t1.0\tSNV\n");
    }

    static void stats(final String label, final Path dir, final String... shards) throws Exception {
        final List<String> argv = new ArrayList<>();
        for (int i = 0; i < shards.length; i++) {
            final Path shard = dir.resolve("stats-" + label + "-" + i + ".tsv");
            Files.write(shard, shards[i].getBytes());
            argv.addAll(Arrays.asList("--stats", shard.toString()));
        }
        final Path out = dir.resolve("stats-" + label + ".tsv");
        argv.addAll(Arrays.asList("-O", out.toString()));
        emit("stats", label, out, () -> new MergeMutectStats().instanceMain(argv.toArray(new String[0])));
    }

    static void pileup(final String label, final Path dir, final Path dict, final String... shards)
            throws Exception {
        final List<String> argv = new ArrayList<>();
        for (int i = 0; i < shards.length; i++) {
            final Path shard = dir.resolve("pileup-" + label + "-" + i + ".tsv");
            Files.write(shard, shards[i].getBytes());
            argv.addAll(Arrays.asList("-I", shard.toString()));
        }
        final Path out = dir.resolve("pileup-" + label + ".tsv");
        argv.addAll(Arrays.asList("--sequence-dictionary", dict.toString(), "-O", out.toString()));
        emit("pileup", label, out,
                () -> new GatherPileupSummaries().instanceMain(argv.toArray(new String[0])));
    }

    static void artifact(final String label, final Path dir, final String... shards)
            throws Exception {
        final List<String> argv = new ArrayList<>();
        for (int i = 0; i < shards.length; i++) {
            final Path shard = dir.resolve("artifact-" + label + "-" + i + ".tsv");
            Files.write(shard, shards[i].getBytes());
            argv.addAll(Arrays.asList("-I", shard.toString()));
        }
        final Path out = dir.resolve("artifact-" + label + ".tsv");
        argv.addAll(Arrays.asList("-O", out.toString()));
        emit("artifact", label, out,
                () -> new GatherNormalArtifactData().instanceMain(argv.toArray(new String[0])));
    }

    interface Run {
        void go();
    }

    static void emit(final String kind, final String label, final Path out, final Run run)
            throws Exception {
        try {
            run.go();
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s-%s\t%s:%s%n", kind, label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("%s\t%s=%s%n", kind, label,
                ReferenceQueryDump.escape(new String(Files.readAllBytes(out))));
    }
}

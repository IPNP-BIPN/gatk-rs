/*
 * The seventeen Spark tools this port runs, end to end through Main, against `--spark-master local[1]`.
 *
 * The covering arrays compare the port against the reference row by row on every CI run; what
 * they cannot do is hold a golden, so nothing re-derives a fixed answer for these tools. This dump
 * is that answer: every run goes through Main, the door a user takes, so GATKConfig's level-2
 * BGZF and its other injected properties are in force, and every path is relative, so the command
 * line PrintReadsSpark records in its @PG line is the same in any directory.
 *
 * What is measured:
 *
 *   - THE THREE COUNTS ARE PRINTED, NOT RETURNED: `runTool` println's the number and `doWork`
 *     returns null, so there is no `Tool returned:` line, and `-O` gets the same text;
 *   - COUNTREADSSPARK DEFAULTS TO ALLOW_ALL_READS, the other two to the wellformed filter, so the
 *     read with no read group is counted by one and not the others;
 *   - WITHOUT AN .sbi, DISQ GUESSES THE FIRST RECORD, scanning from the header's own block for
 *     ten records that chain. A mapped read with no cigar is never a start, so the fixture's chr1
 *     reads, which all lie within ten records of one, are never read: the counts are 3, not 11;
 *   - COUNTVARIANTSSPARK WITH -L QUERIES THROUGH THE INDEX, so the same VCF without its `.idx` is
 *     refused while a whole-file count is not;
 *   - PRINTREADSSPARK WRITES DISQ'S BAM, header and records in separate blocks, with `.bai` and
 *     `.sbi` beside a single file and a directory of one part and `_SUCCESS` for sharded output;
 *     through the JDK deflater, so the bytes replay on a host without ISA-L; GKL's level 2 is the
 *     covering array's to compare;
 *   - BASERECALIBRATORSPARK IS THE WALKER'S ENGINE ON ONE PARTITION, over the reads Disq's guess
 *     leaves; APPLYBQSRSPARK APPLIES ITS TABLE through PrintReadsSpark's writer, and
 *     BQSRPIPELINESPARK CHAINS THE TWO without writing the table;
 *   - AND PILEUPSPARK SAVES A TEXT FILE DIRECTORY, `part-00000` and `_SUCCESS`, with Hadoop's own
 *     checksum files beside them.
 *
 * Output:
 *
 *     input\t<name>\t<base64 of a fixture file>
 *     argv\t<label>\t<the command line, tab-separated, escaped>
 *     stdout\t<label>\t<what the tool printed, escaped>
 *     file\t<label>\t<path under the output>\t<base64>
 *     error\t<label>\t<exception class>:<message, escaped, cwd masked>
 *
 * Usage: SparkWalkersDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.Main;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class SparkWalkersDump {

    static final Path DIR = Path.of("sparkwalkers-dump");

    static final String VCF = "##fileformat=VCFv4.2\n"
            + "##contig=<ID=chr1,length=200>\n"
            + "##contig=<ID=chr2,length=200>\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            + "chr1\t10\t.\tA\tC\t50\t.\t.\n"
            + "chr1\t60\t.\tA\tG\t50\t.\t.\n"
            + "chr1\t150\t.\tC\tT\t50\t.\t.\n"
            + "chr2\t20\t.\tG\tA\t50\t.\t.\n";

    /** The same sites with one sample's calls, which PrintVariantsSpark cannot carry through Kryo. */
    static final String READ_NAMES = "@r001\nr002/1\nr003/2\nnobody\n";

    static final String GENOTYPED_VCF = "##fileformat=VCFv4.2\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            + "##contig=<ID=chr1,length=200>\n"
            + "##contig=<ID=chr2,length=200>\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n"
            + "chr1\t10\t.\tA\tC\t50\t.\t.\tGT\t0/1\n"
            + "chr2\t20\t.\tG\tA\t50\t.\t.\tGT\t1/1\n";

    public static void main(final String[] args) throws Exception {
        PrintReadsDump.emptyDirectory(DIR);
        Files.createDirectories(DIR);

        System.out.println("# SparkWalkersDump: twenty-one Spark tools through Main, --spark-master local[1]");

        final Path fasta = DIR.resolve("ref.fasta");
        Files.write(fasta, ReadWalkerDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + DIR.resolve("ref.dict")});
        final Path bam = DIR.resolve("reads.bam");
        ReadWalkerDump.buildFixture(bam.toFile());
        final Path plainVcf = DIR.resolve("plain.vcf");
        Files.write(plainVcf, VCF.getBytes());
        final Path indexedVcf = DIR.resolve("indexed.vcf");
        Files.write(indexedVcf, VCF.getBytes());
        new Main().instanceMain(new String[] {"IndexFeatureFile", "-I", indexedVcf.toString()});
        final Path genotypedVcf = DIR.resolve("genotyped.vcf");
        Files.write(genotypedVcf, GENOTYPED_VCF.getBytes());
        Files.write(DIR.resolve("names.txt"), READ_NAMES.getBytes());
        Files.write(DIR.resolve("revert.map.tsv"), ("READ_GROUP_ID\tOUTPUT\nrg1\t"
                + DIR.resolve("out/revert-sam-map/mapped.bam") + "\n").getBytes());
        Files.write(DIR.resolve("revert.bad-map.tsv"), "ID\tFILE\nrg1\tx.bam\n".getBytes());

        // Pairs of BAMs that agree or disagree in one way each, for CompareDuplicatesSpark.
        final String[][] pairs = {
                {"dup-a.bam", "", "", "", ""},
                // The same reads, with the marks moved and a tag the tool throws away.
                {"dup-marks.bam", "a1=none,b1=dup,b2=none,c1=dup,a2=tag", "", "", ""},
                // One name changed where the group is otherwise whole.
                {"dup-renamed.bam", "", "a2", "", ""},
                // A read moved out of its group: the groups differ in size, and the new one is
                // only in the second BAM, where nothing looks at it.
                {"dup-uneven.bam", "", "", "a2", ""},
                // A group of one moved: the first BAM's key finds nothing in the second.
                {"dup-moved.bam", "", "", "d1", ""},
                // A read with no group, and one with a group the first header does not know.
                {"dup-no-group.bam", "", "", "", "a1=none"},
                {"dup-other-group.bam", "", "", "", "a1=rg3"},
                // Fewer reads altogether.
                {"dup-fewer.bam", "drop=e1", "", "", ""},
        };
        for (final String[] pair : pairs) {
            duplicateFixture(DIR.resolve(pair[0]), pair[1], pair[2], pair[3], pair[4]);
        }

        markingFixture(DIR.resolve("marking.bam"), SAMFileHeader.SortOrder.coordinate, "");
        markingFixture(DIR.resolve("marking-queryname.bam"), SAMFileHeader.SortOrder.queryname, "");
        markingFixture(DIR.resolve("marking-unsorted.bam"), SAMFileHeader.SortOrder.unsorted, "");
        markingFixture(DIR.resolve("marking-one-library.bam"), SAMFileHeader.SortOrder.coordinate,
                "one-library");
        markingFixture(DIR.resolve("marking-clean.bam"), SAMFileHeader.SortOrder.coordinate, "clean");
        markingFixture(DIR.resolve("marking-no-groups.bam"), SAMFileHeader.SortOrder.coordinate,
                "no-groups");
        markingFixture(DIR.resolve("marking-no-read-group.bam"), SAMFileHeader.SortOrder.coordinate,
                "no-read-group");
        markingFixture(DIR.resolve("marking-unknown-read-group.bam"),
                SAMFileHeader.SortOrder.coordinate, "unknown-read-group");
        markingFixture(DIR.resolve("marking-three.bam"), SAMFileHeader.SortOrder.coordinate, "three");
        markingFixture(DIR.resolve("marking-collision.bam"), SAMFileHeader.SortOrder.coordinate,
                "collision");
        markingFixture(DIR.resolve("marking-scattered.bam"), SAMFileHeader.SortOrder.unsorted,
                "scattered");
        insertSizeFixture(DIR.resolve("pairs.bam"), SAMFileHeader.SortOrder.coordinate, false);
        insertSizeFixture(DIR.resolve("pairs-away.bam"), SAMFileHeader.SortOrder.coordinate, true);
        insertSizeFixture(DIR.resolve("pairs-unsorted.bam"), SAMFileHeader.SortOrder.unsorted,
                false);

        // The fixtures travel, so the port runs over these bytes rather than its own rebuild.
        for (final String name : new String[] {"ref.fasta", "ref.fasta.fai", "ref.dict",
                "reads.bam", "reads.bai", "plain.vcf", "genotyped.vcf", "revert.map.tsv",
                "revert.bad-map.tsv", "names.txt", "dup-a.bam", "dup-a.bai", "dup-marks.bam",
                "dup-marks.bai", "dup-renamed.bam", "dup-uneven.bam", "dup-moved.bam",
                "dup-no-group.bam", "dup-other-group.bam", "dup-fewer.bam", "pairs.bam",
                "pairs.bai", "pairs-away.bam", "pairs-unsorted.bam", "marking.bam", "marking.bai",
                "marking-queryname.bam", "marking-unsorted.bam", "marking-one-library.bam",
                "marking-clean.bam", "marking-no-groups.bam", "marking-no-read-group.bam",
                "marking-unknown-read-group.bam", "marking-three.bam", "marking-collision.bam",
                "marking-scattered.bam"}) {
            System.out.printf("input\t%s\t%s%n", name,
                    Base64.getEncoder().encodeToString(Files.readAllBytes(DIR.resolve(name))));
        }

        // COMPAREDUPLICATESSPARK prints what it found; only the second BAM skips the read filter.
        final String compare = "CompareDuplicatesSpark";
        final String dupA = DIR.resolve("dup-a.bam").toString();
        final String[][] comparisons = {
                {"compare-same", "dup-a.bam"},
                {"compare-marks", "dup-marks.bam"},
                {"compare-marks-throw", "dup-marks.bam", "--throw-on-diff", "true"},
                {"compare-same-throw", "dup-a.bam", "--throw-on-diff", "true"},
                {"compare-no-summary", "dup-marks.bam", "--print-summary", "false"},
                {"compare-renamed", "dup-renamed.bam"},
                {"compare-uneven", "dup-uneven.bam"},
                {"compare-moved", "dup-moved.bam"},
                {"compare-moved-reducers", "dup-moved.bam", "--num-reducers", "3"},
                {"compare-no-group", "dup-no-group.bam"},
                {"compare-other-group", "dup-other-group.bam"},
                {"compare-fewer", "dup-fewer.bam"},
                {"compare-filter", "dup-a.bam", "-RF", "MappingQualityReadFilter",
                        "--minimum-mapping-quality", "30"},
                {"compare-intervals", "dup-marks.bam", "-L", "chr1:100-300"},
                {"compare-only-output", "dup-a.bam", "-O", "only.bam"},
        };
        for (final String[] row : comparisons) {
            final List<String> argv = new ArrayList<>(Arrays.asList("-I", dupA,
                    "--input2", DIR.resolve(row[1]).toString()));
            argv.addAll(Arrays.asList(row).subList(2, row.length));
            run(row[0], null, compare, argv.toArray(new String[0]));
        }
        run("compare-output", "first.bam", compare, "-I", dupA,
                "-I2", DIR.resolve("dup-marks.bam").toString(), "--use-jdk-deflater", "true",
                "-O", out("compare-output", "first.bam"),
                "-O2", out("compare-output", "second.bam"));

        // COLLECTINSERTSIZEMETRICSSPARK is GATK's collector, not Picard's: ten filters, FR, TANDEM,
        // RF, a width the first category fixes for the others, and no file when none passes.
        final String insert = "CollectInsertSizeMetricsSpark";
        final String pairsBam = DIR.resolve("pairs.bam").toString();
        final String[][] inserts = {
                {"insert-size"},
                {"insert-size-levels", "--metric-accumulation-level", "SAMPLE",
                        "--metric-accumulation-level", "READ_GROUP",
                        "--metric-accumulation-level", "LIBRARY"},
                {"insert-size-all-levels", "--metric-accumulation-level", "ALL_READS",
                        "--metric-accumulation-level", "SAMPLE",
                        "--metric-accumulation-level", "LIBRARY",
                        "--metric-accumulation-level", "READ_GROUP"},
                {"insert-size-width", "--width", "125"},
                {"insert-size-width-one", "--width", "1"},
                {"insert-size-tolerance", "-TOL", "1.0"},
                {"insert-size-minimum", "-M", "0.2"},
                {"insert-size-all-kept", "-M", "0.0"},
                {"insert-size-interval", "-L", "chr1:100-600"},
                {"insert-size-filter", "-RF", "ReadLengthReadFilter", "--min-read-length", "1",
                        "--max-read-length", "5"},
                {"insert-size-second", "-DF", "FirstOfPairReadFilter", "-RF",
                        "SecondOfPairReadFilter"},
                {"insert-size-no-defaults", "--disable-tool-default-read-filters", "true"},
                {"insert-size-sample", "--metric-accumulation-level", "SAMPLE",
                        "--disable-tool-default-read-filters", "true"},
        };
        for (final String[] row : inserts) {
            final List<String> argv = new ArrayList<>(Arrays.asList("-I", pairsBam, "-H",
                    "unused.pdf",
                    "-O", out(row[0], "insert.metrics")));
            argv.addAll(Arrays.asList(row).subList(1, row.length));
            run(row[0], "insert.metrics", insert, argv.toArray(new String[0]));
        }
        run("insert-size-away", null, insert, "-I", DIR.resolve("pairs-away.bam").toString(),
                "-H", "unused.pdf", "-O", out("insert-size-away", "insert.metrics"));
        run("insert-size-none-pass", "insert.metrics", insert, "-I", pairsBam, "-H", "unused.pdf",
                "-M", "1.0", "-O", out("insert-size-none-pass", "insert.metrics"));
        run("insert-size-unsorted", null, insert, "-I", DIR.resolve("pairs-unsorted.bam").toString(),
                "-H", "unused.pdf", "-O", out("insert-size-unsorted", "insert.metrics"));
        run("insert-size-no-reads", null, insert, "-I", bam.toString(), "-H", "unused.pdf",
                "-O", out("insert-size-no-reads", "insert.metrics"));
        run("insert-size-no-plot-file", null, insert, "-I", pairsBam,
                "-O", out("insert-size-no-plot-file", "insert.metrics"));

        // COLLECTMULTIPLEMETRICSSPARK runs the two collectors over the same reads, in the order
        // asked for, each through its own default filters: the command line's are parsed and then
        // bypassed. -O is a base name, and each collector's failure names its enum constant.
        final String multiple = "CollectMultipleMetricsSpark";
        final String insertOnly = "CollectInsertSizeMetrics";
        final String yieldOnly = "CollectQualityYieldMetrics";
        final String[][] multiples = {
                {"multiple"},
                {"multiple-insert", "--collectors", insertOnly},
                {"multiple-yield", "--collectors", yieldOnly},
                {"multiple-reversed", "--collectors", yieldOnly, "--collectors", insertOnly},
                {"multiple-twice", "--collectors", insertOnly, "--collectors", insertOnly},
                {"multiple-levels", "--metric-accumulation-level", "SAMPLE",
                        "--metric-accumulation-level", "READ_GROUP"},
                {"multiple-interval", "-L", "chr1:100-600"},
                {"multiple-bypassed", "-RF", "ReadLengthReadFilter", "--min-read-length", "1000",
                        "--max-read-length", "2000"},
                {"multiple-no-defaults", "--disable-tool-default-read-filters", "true"},
                {"multiple-not-assumed", "--assume-sorted", "false"},
        };
        for (final String[] row : multiples) {
            final List<String> argv = new ArrayList<>(Arrays.asList("-I", pairsBam,
                    "-O", out(row[0], "metrics")));
            argv.addAll(Arrays.asList(row).subList(1, row.length));
            run(row[0], "metrics", multiple, argv.toArray(new String[0]));
        }
        final String pairsUnsorted = DIR.resolve("pairs-unsorted.bam").toString();
        final String pairsAway = DIR.resolve("pairs-away.bam").toString();
        run("multiple-unsorted", null, multiple, "-I", pairsUnsorted,
                "-O", out("multiple-unsorted", "metrics"));
        run("multiple-unsorted-yield", null, multiple, "-I", pairsUnsorted,
                "--collectors", yieldOnly, "-O", out("multiple-unsorted-yield", "metrics"));
        run("multiple-away", null, multiple, "-I", pairsAway, "-O", out("multiple-away", "metrics"));
        run("multiple-no-reads", null, multiple, "-I", bam.toString(),
                "-O", out("multiple-no-reads", "metrics"));
        run("multiple-no-reads-yield", "metrics", multiple, "-I", bam.toString(),
                "--collectors", yieldOnly, "-O", out("multiple-no-reads-yield", "metrics"));
        run("multiple-bad-collector", null, multiple, "-I", pairsBam, "--collectors", "Nope",
                "-O", out("multiple-bad-collector", "metrics"));
        run("multiple-no-output", null, multiple, "-I", pairsBam);

        // MARKDUPLICATESSPARK answers with the names that are not duplicates, so every read of a
        // name follows it, and it always writes a coordinate-sorted file.
        final String marks = "MarkDuplicatesSpark";
        final String marking = DIR.resolve("marking.bam").toString();
        final String[][] markRows = {
                {"mark"},
                {"mark-all", "--duplicate-tagging-policy", "All"},
                {"mark-optical", "--duplicate-tagging-policy", "OpticalOnly"},
                {"mark-remove-all", "--remove-all-duplicates", "true"},
                {"mark-remove-sequencing", "--remove-sequencing-duplicates", "true"},
                {"mark-length", "-DS", "TOTAL_MAPPED_REFERENCE_LENGTH"},
                {"mark-flow-score", "-DS", "FLOW_SUM_OF_BASE_QUALITIES"},
                {"mark-keep-mates", "--do-not-mark-unmapped-mates", "true"},
                {"mark-no-names", "--read-name-regex", "null", "--duplicate-tagging-policy", "All"},
                {"mark-regex", "--read-name-regex", "[^:]*:[^:]*:([0-9]+):([0-9]+):([0-9]+)",
                        "--duplicate-tagging-policy", "All"},
                {"mark-distance", "--optical-duplicate-pixel-distance", "5",
                        "--duplicate-tagging-policy", "All"},
                {"mark-wide", "--optical-duplicate-pixel-distance", "100000",
                        "--duplicate-tagging-policy", "All"},
                {"mark-interval", "-L", "chr1:100-300"},
                {"mark-reducers", "--num-reducers", "2"},
                {"mark-filter", "-RF", "MappingQualityReadFilter", "--minimum-mapping-quality", "30"},
        };
        for (final String[] row : markRows) {
            final List<String> argv = new ArrayList<>(Arrays.asList("-I", marking,
                    "--use-jdk-deflater", "true",
                    "-O", out(row[0], "marked.bam"), "-M", out(row[0], "marked.metrics")));
            argv.addAll(Arrays.asList(row).subList(1, row.length));
            run(row[0], "marked.bam", marks, argv.toArray(new String[0]));
        }
        final String[][] markInputs = {
                {"mark-queryname", "marking-queryname.bam"},
                {"mark-unsorted", "marking-unsorted.bam"},
                {"mark-one-library", "marking-one-library.bam"},
                {"mark-clean", "marking-clean.bam"},
                {"mark-no-groups", "marking-no-groups.bam"},
                {"mark-no-read-group", "marking-no-read-group.bam"},
                {"mark-unknown-read-group", "marking-unknown-read-group.bam"},
                {"mark-three", "marking-three.bam"},
                {"mark-collision", "marking-collision.bam"},
                {"mark-scattered", "marking-scattered.bam"},
        };
        for (final String[] row : markInputs) {
            run(row[0], "marked.bam", marks, "-I", DIR.resolve(row[1]).toString(),
                    "--use-jdk-deflater", "true",
                    "-O", out(row[0], "marked.bam"), "-M", out(row[0], "marked.metrics"));
        }
        run("mark-no-metrics", "marked.bam", marks, "-I", marking, "--use-jdk-deflater", "true",
                "-O", out("mark-no-metrics", "marked.bam"));
        run("mark-three-no-metrics", null, marks, "-I", DIR.resolve("marking-three.bam").toString(),
                "-O", out("mark-three-no-metrics", "marked.bam"));
        run("mark-treated", "marked.bam", marks, "-I", DIR.resolve("marking-unsorted.bam").toString(),
                "--treat-unsorted-as-querygroup-ordered", "true", "--use-jdk-deflater", "true",
                "-O", out("mark-treated", "marked.bam"), "-M", out("mark-treated", "marked.metrics"));
        run("mark-mutex", null, marks, "-I", marking, "--remove-all-duplicates", "true",
                "--duplicate-tagging-policy", "All", "-O", out("mark-mutex", "marked.bam"));
        run("mark-no-output", null, marks, "-I", marking);
        run("mark-interval-no-index", null, marks, "-I", DIR.resolve("marking-queryname.bam").toString(),
                "-L", "chr1:100-300", "-O", out("mark-interval-no-index", "marked.bam"));

        final String r = fasta.toString();
        final String i = bam.toString();
        run("count-reads", null, "CountReadsSpark", "-I", i);
        run("count-reads-chr1", null, "CountReadsSpark", "-I", i, "-L", "chr1");
        run("count-reads-wellformed", null, "CountReadsSpark", "-I", i,
                "--read-filter", "WellformedReadFilter");
        run("count-reads-output", "count-reads.txt", "CountReadsSpark", "-I", i,
                "-O", out("count-reads-output", "count-reads.txt"));
        run("count-bases", null, "CountBasesSpark", "-I", i);
        run("count-bases-chr2", null, "CountBasesSpark", "-I", i, "-R", r, "-L", "chr2");
        run("count-bases-no-filters", null, "CountBasesSpark", "-I", i,
                "--disable-tool-default-read-filters");
        run("flagstat", null, "FlagStatSpark", "-I", i);
        run("flagstat-output", "flagstat.txt", "FlagStatSpark", "-I", i,
                "-O", out("flagstat-output", "flagstat.txt"));
        run("count-variants", null, "CountVariantsSpark", "-V", plainVcf.toString());
        run("count-variants-chr1", null, "CountVariantsSpark", "-V", indexedVcf.toString(),
                "-R", r, "-L", "chr1");
        run("count-variants-unindexed", null, "CountVariantsSpark", "-V", plainVcf.toString(),
                "-R", r, "-L", "chr1");
        run("print-reads", "print-reads.bam", "PrintReadsSpark", "-I", i, "--use-jdk-deflater", "true",
                "-O", out("print-reads", "print-reads.bam"));
        run("print-reads-chr1", "print-reads-chr1.bam", "PrintReadsSpark", "-I", i, "--use-jdk-deflater", "true", "-L", "chr1",
                "-O", out("print-reads-chr1", "print-reads-chr1.bam"));
        run("print-reads-sharded", "print-reads-sharded.bam", "PrintReadsSpark", "-I", i, "--use-jdk-deflater", "true",
                "--sharded-output", "true", "-O", out("print-reads-sharded", "print-reads-sharded.bam"));
        run("pileup", "pileup", "PileupSpark", "-I", i, "-R", r, "-O", out("pileup", "pileup"));
        run("pileup-chr1", "pileup-chr1", "PileupSpark", "-I", i, "-R", r, "-L", "chr1:100-200",
                "-O", out("pileup-chr1", "pileup-chr1"));
        run("pileup-shuffle", "pileup-shuffle", "PileupSpark", "-I", i, "-R", r, "--shuffle",
                "true", "-O", out("pileup-shuffle", "pileup-shuffle"));
        // BQSR: the table over the indexed VCF's sites, the table applied by the next run, and the
        // two passes chained in one job with the table never written.
        final String sites = indexedVcf.toString();
        run("base-recalibrator", "recal.table", "BaseRecalibratorSpark", "-I", i, "-R", r,
                "--known-sites", sites, "-O", out("base-recalibrator", "recal.table"));
        run("apply-bqsr", "apply-bqsr.bam", "ApplyBQSRSpark", "-I", i, "--use-jdk-deflater",
                "true", "--bqsr-recal-file", out("base-recalibrator", "recal.table"),
                "-O", out("apply-bqsr", "apply-bqsr.bam"));
        run("bqsr-pipeline", "bqsr-pipeline.bam", "BQSRPipelineSpark", "-I", i, "-R", r,
                "--use-jdk-deflater", "true", "--known-sites", sites,
                "-O", out("bqsr-pipeline", "bqsr-pipeline.bam"));
        // SortSamSpark: PrintReadsSpark's writer asked to sort to the header, by coordinate (the
        // default), by name, into shards, and over an interval; an order Spark has no sort for is
        // Barclay's to refuse.
        run("sort-sam", "sort-sam.bam", "SortSamSpark", "-I", i, "--use-jdk-deflater", "true",
                "-O", out("sort-sam", "sort-sam.bam"));
        run("sort-sam-queryname", "sort-sam-queryname.bam", "SortSamSpark", "-I", i,
                "--use-jdk-deflater", "true", "--sort-order", "queryname",
                "-O", out("sort-sam-queryname", "sort-sam-queryname.bam"));
        run("sort-sam-sharded", "sort-sam-sharded.bam", "SortSamSpark", "-I", i,
                "--use-jdk-deflater", "true", "--sharded-output", "true",
                "-O", out("sort-sam-sharded", "sort-sam-sharded.bam"));
        run("sort-sam-chr2", "sort-sam-chr2.bam", "SortSamSpark", "-I", i, "--use-jdk-deflater",
                "true", "-L", "chr2", "-O", out("sort-sam-chr2", "sort-sam-chr2.bam"));
        run("sort-sam-unsorted", "sort-sam-unsorted.bam", "SortSamSpark", "-I", i,
                "--use-jdk-deflater", "true", "--sort-order", "unsorted",
                "-O", out("sort-sam-unsorted", "sort-sam-unsorted.bam"));
        // RevertSamSpark: the reversion and each switch that holds part of it back, the arguments
        // the tool has quirks on (`--remove-default-attributes-to-clear` fails, and the output
        // checks' answers are thrown away), `--sanitize` with its sort, and one file per read group.
        final String revert = "RevertSamSpark";
        run("revert-sam", "revert-sam.bam", revert, "-I", i, "--use-jdk-deflater", "true",
                "-O", out("revert-sam", "revert-sam.bam"));
        run("revert-sam-keep-alignment", "revert-sam-keep-alignment.bam", revert, "-I", i,
                "--use-jdk-deflater", "true", "--keep-alignment-information", "true",
                "--sort-order", "coordinate",
                "-O", out("revert-sam-keep-alignment", "revert-sam-keep-alignment.bam"));
        run("revert-sam-qualities", "revert-sam-qualities.bam", revert, "-I", i,
                "--use-jdk-deflater", "true", "--dont-restore-original-qualities", "true",
                "-O", out("revert-sam-qualities", "revert-sam-qualities.bam"));
        run("revert-sam-duplicates", "revert-sam-duplicates.bam", revert, "-I", i,
                "--use-jdk-deflater", "true", "--remove-duplicate-information", "true",
                "-O", out("revert-sam-duplicates", "revert-sam-duplicates.bam"));
        run("revert-sam-attributes", "revert-sam-attributes.bam", revert, "-I", i,
                "--use-jdk-deflater", "true", "--attributes-to-clear", "NM",
                "--attributes-to-clear", "RG",
                "-O", out("revert-sam-attributes", "revert-sam-attributes.bam"));
        run("revert-sam-defaults", "revert-sam-defaults.bam", revert, "-I", i,
                "--remove-default-attributes-to-clear", "true",
                "-O", out("revert-sam-defaults", "revert-sam-defaults.bam"));
        run("revert-sam-alias", "revert-sam-alias.bam", revert, "-I", i, "--use-jdk-deflater",
                "true", "--sample-alias", "alias", "--library-name", "library",
                "-O", out("revert-sam-alias", "revert-sam-alias.bam"));
        run("revert-sam-sanitize", "revert-sam-sanitize.bam", revert, "-I", i,
                "--use-jdk-deflater", "true", "--sanitize", "true",
                "-O", out("revert-sam-sanitize", "revert-sam-sanitize.bam"));
        run("revert-sam-sanitize-first", "revert-sam-sanitize-first.bam", revert, "-I", i,
                "--use-jdk-deflater", "true", "--sanitize", "true", "--keep-first-duplicate", "true",
                "--num-reducers", "2",
                "-O", out("revert-sam-sanitize-first", "revert-sam-sanitize-first.bam"));
        run("revert-sam-first-alone", "revert-sam-first-alone.bam", revert, "-I", i,
                "--keep-first-duplicate", "true",
                "-O", out("revert-sam-first-alone", "revert-sam-first-alone.bam"));
        run("revert-sam-sharded", "revert-sam-sharded.bam", revert, "-I", i, "--use-jdk-deflater",
                "true", "--sharded-output", "true",
                "-O", out("revert-sam-sharded", "revert-sam-sharded.bam"));
        run("revert-sam-sam", "revert-sam-sam.sam", revert, "-I", i,
                "-O", out("revert-sam-sam", "revert-sam-sam.sam"));
        run("revert-sam-by-group", "by-group", revert, "-I", i, "--use-jdk-deflater", "true",
                "--output-by-readgroup", "true", "-O", out("revert-sam-by-group", "by-group"));
        run("revert-sam-by-group-sam", "by-group", revert, "-I", i,
                "--output-by-readgroup", "true", "--output-by-readgroup-file-format", "sam",
                "-O", out("revert-sam-by-group-sam", "by-group"));
        run("revert-sam-no-output", null, revert, "-I", i, "--output-by-readgroup", "true");
        final String map = DIR.resolve("revert.map.tsv").toString();
        run("revert-sam-map", "mapped.bam", revert, "-I", i, "--use-jdk-deflater", "true",
                "--output-by-readgroup", "true", "--output-map", map);
        out("revert-sam-map", "mapped.bam");
        run("revert-sam-map-alone", null, revert, "-I", i, "--output-map", map);
        run("revert-sam-map-bad", null, revert, "-I", i, "--output-by-readgroup", "true",
                "--output-map", DIR.resolve("revert.bad-map.tsv").toString());
        run("revert-sam-map-missing", null, revert, "-I", i, "--output-by-readgroup", "true",
                "--output-map", DIR.resolve("no-such.tsv").toString());
        // A filter by read name, whose `@`, `/1` and `/2` come off the file's lines, in both senses
        // and with a read filter it parses and never applies.
        final Path names = DIR.resolve("names.txt");
        final String extract = "ExtractOriginalAlignmentRecordsByNameSpark";
        run("extract-by-name", "extract-by-name.bam", extract, "-I", i, "--use-jdk-deflater",
                "true", "--read-name-file", names.toString(),
                "-O", out("extract-by-name", "extract-by-name.bam"));
        run("extract-by-name-invert", "extract-by-name-invert.bam", extract, "-I", i,
                "--use-jdk-deflater", "true", "--read-name-file", names.toString(),
                "--invert-match", "true",
                "-O", out("extract-by-name-invert", "extract-by-name-invert.bam"));
        run("extract-by-name-missing", null, extract, "-I", i,
                "--read-name-file", DIR.resolve("no-such-names.txt").toString(),
                "-O", out("extract-by-name-missing", "extract-by-name-missing.bam"));
        // The metrics the Spark pipelines compute themselves. Their files are compared with the
        // time of the run taken out, and `--chart` needs an R the container does not have.
        run("mean-quality", "mean-quality.metrics", "MeanQualityByCycleSpark", "-I", i,
                "-O", out("mean-quality", "mean-quality.metrics"));
        run("mean-quality-pf", "mean-quality-pf.metrics", "MeanQualityByCycleSpark", "-I", i,
                "--pfReadsOnly", "true", "--alignedReadsOnly", "true",
                "-O", out("mean-quality-pf", "mean-quality-pf.metrics"));
        run("mean-quality-no-output", null, "MeanQualityByCycleSpark", "-I", i);
        run("quality-distribution", "quality-distribution.metrics",
                "QualityScoreDistributionSpark", "-I", i,
                "-O", out("quality-distribution", "quality-distribution.metrics"));
        run("quality-distribution-no-calls", "quality-distribution-no-calls.metrics",
                "QualityScoreDistributionSpark", "-I", i, "--includeNoCalls", "true",
                "--alignedReadsOnly", "true",
                "-O", out("quality-distribution-no-calls", "quality-distribution-no-calls.metrics"));
        run("base-distribution", "base-distribution.metrics",
                "CollectBaseDistributionByCycleSpark", "-I", i,
                "-O", out("base-distribution", "base-distribution.metrics"));
        run("base-distribution-aligned", "base-distribution-aligned.metrics",
                "CollectBaseDistributionByCycleSpark", "-I", i, "--aligned-reads-only", "true",
                "--pf-reads-only", "true",
                "-O", out("base-distribution-aligned", "base-distribution-aligned.metrics"));
        run("quality-yield", "quality-yield.metrics", "CollectQualityYieldMetricsSpark", "-I", i,
                "-O", out("quality-yield", "quality-yield.metrics"));
        run("quality-yield-original", "quality-yield-original.metrics",
                "CollectQualityYieldMetricsSpark", "-I", i, "--use-original-qualities", "true",
                "-O", out("quality-yield-original", "quality-yield-original.metrics"));
        // PrintVariantsSpark: Disq's VCF writer, as text, block-compressed with its tabix index,
        // over an interval of the indexed copy, refused for BCF, and failing inside the job on a
        // record with genotypes.
        run("print-variants", "print-variants.vcf", "PrintVariantsSpark", "-V",
                plainVcf.toString(), "-O", out("print-variants", "print-variants.vcf"));
        run("print-variants-gz", "print-variants-gz.vcf.gz", "PrintVariantsSpark", "--use-jdk-deflater", "true", "-V",
                plainVcf.toString(), "-O", out("print-variants-gz", "print-variants-gz.vcf.gz"));
        run("print-variants-chr1", "print-variants-chr1.vcf", "PrintVariantsSpark", "-V",
                indexedVcf.toString(), "-R", r, "-L", "chr1",
                "-O", out("print-variants-chr1", "print-variants-chr1.vcf"));
        run("print-variants-no-index", "print-variants-no-index.vcf.gz", "PrintVariantsSpark",
                "--use-jdk-deflater", "true", "-V", plainVcf.toString(), "--create-output-variant-index", "false",
                "-O", out("print-variants-no-index", "print-variants-no-index.vcf.gz"));
        run("print-variants-bcf", "print-variants-bcf.bcf", "PrintVariantsSpark", "-V",
                plainVcf.toString(), "-O", out("print-variants-bcf", "print-variants-bcf.bcf"));
        run("print-variants-genotyped", "print-variants-genotyped.vcf", "PrintVariantsSpark",
                "-V", genotypedVcf.toString(),
                "-O", out("print-variants-genotyped", "print-variants-genotyped.vcf"));
    }

    /**
     * Pairs in two read groups of one sample: rg1 has a library and a platform unit, rg2 neither.
     * Mostly FR with one far outlier, then RF and TANDEM pairs, a duplicate and a pair with no
     * group; `away` adds a read on another contig than its mate whose flag still says proper.
     */
    /**
     * Reads for every record type MarkDuplicatesSpark makes: optical and library duplicates among
     * pairs, a set of four for the finder's graph, FR against RF at one key, fragments alone and
     * beside a pair, an unmapped mate, secondaries and supplementaries, a mate missing from the file,
     * names with no location, and a failed vendor check.
     */
    static void markingFixture(final Path bam, final SAMFileHeader.SortOrder order, final String variant) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 100000));
        dictionary.addSequence(new SAMSequenceRecord("chr2", 100000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(order);
        if (variant.equals("scattered")) {
            header.setGroupOrder(SAMFileHeader.GroupOrder.query);
        }
        if (!variant.equals("no-groups")) {
            for (final String[] group : new String[][] {{"rg1", "lib1"}, {"rg2", "lib2"}, {"rg3", null}}) {
                final SAMReadGroupRecord record = new SAMReadGroupRecord(group[0]);
                record.setSample("sample1");
                record.setLibrary(group[1]);
                header.addReadGroup(record);
            }
        }
        final List<SAMRecord> records = new ArrayList<>();
        // name, contig, start, flags, cigar, mate contig, mate start, quality, group, mapping quality
        final Object[][] reads = variant.equals("clean") ? new Object[][] {
                {"c:1:1101:10:10", 0, 100, 0x63, "10M", 0, 200, 30, "rg1", 60},
                {"c:1:1101:10:10", 0, 200, 0x93, "10M", 0, 100, 30, "rg1", 60},
                {"c:1:1101:20:20", 0, 300, 0x63, "10M", 0, 400, 30, "rg1", 60},
                {"c:1:1101:20:20", 0, 400, 0x93, "10M", 0, 300, 30, "rg1", 60},
        } : new Object[][] {
                // Optical and library duplicates of one pair, and the best in the middle.
                {"x:1:1101:1000:1000", 0, 100, 0x63, "10M", 0, 200, 30, "rg1", 60},
                {"x:1:1101:1000:1000", 0, 200, 0x93, "10M", 0, 100, 30, "rg1", 60},
                {"x:1:1101:1010:1010", 0, 100, 0x63, "10M", 0, 200, 40, "rg1", 60},
                {"x:1:1101:1010:1010", 0, 200, 0x93, "10M", 0, 100, 40, "rg1", 60},
                {"x:1:1102:5000:5000", 0, 100, 0x63, "10M", 0, 200, 20, "rg1", 60},
                {"x:1:1102:5000:5000", 0, 200, 0x93, "10M", 0, 100, 20, "rg1", 60},
                // The same place in the other library is no duplicate of it.
                {"y:1:1101:1000:1000", 0, 100, 0x63, "10M", 0, 200, 30, "rg2", 60},
                {"y:1:1101:1000:1000", 0, 200, 0x93, "10M", 0, 100, 30, "rg2", 60},
                // Four and more: the graph, a chain of three and one apart.
                {"g:1:1101:1000:1000", 0, 1000, 0x63, "10M", 0, 1100, 20, "rg1", 60},
                {"g:1:1101:1000:1000", 0, 1100, 0x93, "10M", 0, 1000, 20, "rg1", 60},
                {"g:1:1101:1090:1000", 0, 1000, 0x63, "10M", 0, 1100, 20, "rg1", 60},
                {"g:1:1101:1090:1000", 0, 1100, 0x93, "10M", 0, 1000, 20, "rg1", 60},
                {"g:1:1101:1180:1000", 0, 1000, 0x63, "10M", 0, 1100, 25, "rg1", 60},
                {"g:1:1101:1180:1000", 0, 1100, 0x93, "10M", 0, 1000, 25, "rg1", 60},
                {"g:1:1101:9000:9000", 0, 1000, 0x63, "10M", 0, 1100, 40, "rg1", 60},
                {"g:1:1101:9000:9000", 0, 1100, 0x93, "10M", 0, 1000, 40, "rg1", 60},
                {"g:1:1101:9050:9000", 0, 1000, 0x63, "10M", 0, 1100, 10, "rg1", 60},
                {"g:1:1101:9050:9000", 0, 1100, 0x93, "10M", 0, 1000, 10, "rg1", 60},
                // FR and RF at one key: which read is first of pair flips the optical orientation.
                {"o:1:1101:1000:1000", 0, 2000, 0x63, "10M", 0, 2100, 30, "rg1", 60},
                {"o:1:1101:1000:1000", 0, 2100, 0x93, "10M", 0, 2000, 30, "rg1", 60},
                {"o:1:1101:1020:1000", 0, 2000, 0xA3, "10M", 0, 2100, 20, "rg1", 60},
                {"o:1:1101:1020:1000", 0, 2100, 0x53, "10M", 0, 2000, 20, "rg1", 60},
                {"o:1:1101:1040:1000", 0, 2000, 0x63, "10M", 0, 2100, 10, "rg1", 60},
                {"o:1:1101:1040:1000", 0, 2100, 0x93, "10M", 0, 2000, 10, "rg1", 60},
                {"o:1:1101:1060:1000", 0, 2000, 0xA3, "10M", 0, 2100, 10, "rg1", 60},
                {"o:1:1101:1060:1000", 0, 2100, 0x53, "10M", 0, 2000, 10, "rg1", 60},
                // Both reads on one start, and a pair across contigs.
                {"s:1:1101:1:1", 0, 3000, 0x73, "10M", 0, 3000, 30, "rg1", 60},
                {"s:1:1101:1:1", 0, 3000, 0xB3, "10M", 0, 3000, 30, "rg1", 60},
                {"t:1:1101:1:1", 0, 3100, 0x43, "10M", 1, 3200, 30, "rg1", 60},
                {"t:1:1101:1:1", 1, 3200, 0x83, "10M", 0, 3100, 30, "rg1", 60},
                {"t:1:1101:2:2", 1, 3200, 0x43, "10M", 0, 3100, 40, "rg1", 60},
                {"t:1:1101:2:2", 0, 3100, 0x83, "10M", 1, 3200, 40, "rg1", 60},
                // Fragments on their own, and clipping that moves the unclipped start.
                {"f:1:1101:1:1", 0, 4000, 0, "10M", -1, 0, 30, "rg1", 60},
                {"f:1:1101:2:2", 0, 4000, 0, "10M", -1, 0, 35, "rg1", 60},
                {"f:1:1101:3:3", 0, 4002, 0, "2S8M", -1, 0, 39, "rg1", 60},
                {"f:1:1101:4:4", 0, 4100, 0x10, "10M", -1, 0, 30, "rg1", 60},
                {"f:1:1101:5:5", 0, 4095, 0x10, "5M5S", -1, 0, 30, "rg1", 60},
                // A fragment on a pair's read's place loses to the pair whatever it scores.
                {"h:1:1101:1:1", 0, 5000, 0x63, "10M", 0, 5100, 10, "rg1", 60},
                {"h:1:1101:1:1", 0, 5100, 0x93, "10M", 0, 5000, 10, "rg1", 60},
                {"h:1:1101:2:2", 0, 5000, 0, "10M", -1, 0, 40, "rg1", 60},
                // A mapped read with its unmapped mate, beside a duplicate of it.
                {"u:1:1101:1:1", 0, 6000, 0x49, "10M", 0, 6000, 40, "rg1", 60},
                {"u:1:1101:1:1", 0, 6000, 0x85, "*", 0, 6000, 40, "rg1", 0},
                {"u:1:1101:2:2", 0, 6000, 0x49, "10M", 0, 6000, 20, "rg1", 60},
                {"u:1:1101:2:2", 0, 6000, 0x85, "*", 0, 6000, 20, "rg1", 0},
                // Both unmapped, and unmapped alone.
                {"n:1:1101:1:1", -1, 0, 0x4D, "*", -1, 0, 30, "rg1", 0},
                {"n:1:1101:1:1", -1, 0, 0x8D, "*", -1, 0, 30, "rg1", 0},
                {"n:1:1101:2:2", -1, 0, 0x4, "*", -1, 0, 30, "rg1", 0},
                // A primary that lost keeps its secondary and supplementary company in the verdict.
                {"p:1:1101:1:1", 0, 7000, 0x63, "10M", 0, 7100, 30, "rg1", 60},
                {"p:1:1101:1:1", 0, 7100, 0x93, "10M", 0, 7000, 30, "rg1", 60},
                {"p:1:1101:1:1", 0, 7500, 0x163, "10M", 0, 7100, 30, "rg1", 60},
                {"p:1:1101:2:2", 0, 7000, 0x63, "10M", 0, 7100, 40, "rg1", 60},
                {"p:1:1101:2:2", 0, 7100, 0x93, "10M", 0, 7000, 40, "rg1", 60},
                {"p:1:1101:2:2", 0, 7600, 0x863, "10M", 0, 7100, 40, "rg1", 60},
                // Nothing but a supplementary, and a mate the file does not hold.
                {"q:1:1101:1:1", 0, 7700, 0x800, "10M", -1, 0, 30, "rg1", 60},
                {"l:1:1101:1:1", 0, 7800, 0x63, "10M", 0, 7900, 30, "rg1", 60},
                {"l:1:1101:2:2", 0, 7800, 0x63, "10M", 0, 7900, 30, "rg1", 60},
                // Names with no location, and one where a failed check loses it the score.
                {"plainA", 0, 8000, 0x63, "10M", 0, 8100, 30, "rg1", 60},
                {"plainA", 0, 8100, 0x93, "10M", 0, 8000, 30, "rg1", 60},
                {"plainB", 0, 8000, 0x63, "10M", 0, 8100, 31, "rg1", 60},
                {"plainB", 0, 8100, 0x93, "10M", 0, 8000, 31, "rg1", 60},
                {"v:1:1101:1:1", 0, 8200, 0x263, "10M", 0, 8300, 40, "rg1", 60},
                {"v:1:1101:1:1", 0, 8300, 0x293, "10M", 0, 8200, 40, "rg1", 60},
                {"v:1:1101:2:2", 0, 8200, 0x63, "10M", 0, 8300, 10, "rg1", 60},
                {"v:1:1101:2:2", 0, 8300, 0x93, "10M", 0, 8200, 10, "rg1", 60},
                // The group with no library.
                {"w:1:1101:1:1", 0, 9000, 0, "10M", -1, 0, 30, "rg3", 60},
                {"w:1:1101:2:2", 0, 9000, 0, "10M", -1, 0, 20, "rg3", 60},
                {"w:1:1101:3:3", 0, 9300, 0x400, "20M", -1, 0, 20, "rg3", 60},
        };
        for (final Object[] read : reads) {
            if (variant.equals("one-library") && !"rg1".equals(read[8])) {
                continue;
            }
            final SAMRecord record = new SAMRecord(header);
            record.setReadName((String) read[0]);
            record.setFlags((int) read[3]);
            record.setReferenceIndex((int) read[1]);
            record.setAlignmentStart((int) read[2]);
            record.setCigarString((String) read[4]);
            record.setMateReferenceIndex((int) read[5]);
            record.setMateAlignmentStart((int) read[6]);
            record.setMappingQuality((int) read[9]);
            record.setReadString("ACGTACGTACACGTACGTAC".substring(0, ((String) read[4]).equals("20M") ? 20 : 10));
            final char[] qualities = new char[record.getReadLength()];
            Arrays.fill(qualities, (char) (33 + (int) read[7]));
            record.setBaseQualityString(new String(qualities));
            if (!variant.equals("no-groups")) {
                record.setAttribute("RG", read[8]);
            }
            // A marking the input already carries is taken back.
            if (read[0].equals("x:1:1101:1000:1000")) {
                record.setAttribute("DT", "LB");
            }
            records.add(record);
        }
        final SAMRecord first = records.get(0);
        if (variant.equals("no-read-group")) {
            records.get(4).setAttribute("RG", null);
        }
        if (variant.equals("unknown-read-group")) {
            records.get(4).setAttribute("RG", "elsewhere");
        }
        if (variant.equals("three")) {
            final SAMRecord third = first.deepCopy();
            third.setAlignmentStart(50);
            records.add(third);
        }
        if (variant.equals("collision")) {
            // One name: a read whose mate is mapped and a fragment, both answered as not duplicate.
            final SAMRecord alone = new SAMRecord(header);
            alone.setReadName("l:1:1101:1:1");
            alone.setFlags(0);
            alone.setReferenceIndex(0);
            alone.setAlignmentStart(7850);
            alone.setCigarString("10M");
            alone.setMappingQuality(60);
            alone.setReadString("ACGTACGTAC");
            alone.setBaseQualityString("IIIIIIIIII");
            alone.setAttribute("RG", "rg1");
            records.add(alone);
        }
        final boolean scattered = variant.equals("scattered");
        final boolean presorted = order == SAMFileHeader.SortOrder.unsorted || scattered;
        if (scattered) {
            // Ordered by position under a header that says grouped: a name comes back later.
            records.sort(java.util.Comparator.comparingInt(SAMRecord::getAlignmentStart));
        }
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(order == SAMFileHeader.SortOrder.coordinate)
                .makeBAMWriter(header, presorted, bam.toFile())) {
            records.forEach(writer::addAlignment);
        }
    }

    static void insertSizeFixture(final Path bam, final SAMFileHeader.SortOrder order,
                                  final boolean away) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 100000));
        dictionary.addSequence(new SAMSequenceRecord("chr2", 100000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(order);
        final SAMReadGroupRecord one = new SAMReadGroupRecord("rg1");
        one.setSample("sample1");
        one.setLibrary("lib1");
        one.setPlatformUnit("unit1");
        header.addReadGroup(one);
        final SAMReadGroupRecord two = new SAMReadGroupRecord("rg2");
        two.setSample("sample1");
        header.addReadGroup(two);
        // start, mate start, insert, first's flags beyond paired/first/proper, group.
        final Object[][] pairs = {
                {100, 190, 100, 0x20, "rg1"}, {200, 300, 110, 0x20, "rg1"},
                {300, 390, 100, 0x20, "rg1"}, {400, 510, 120, 0x20, "rg1"},
                {500, 620, 130, 0x20, "rg2"}, {600, 690, 100, 0x20, "rg2"},
                {700, 5690, 5000, 0x20, "rg1"}, {800, 910, 120, 0x420, "rg1"},
                {900, 1010, 120, 0x20, null},
                // The reverse read first: RF, on either side of the forward mate.
                {1100, 1050, -40, 0x10, "rg1"}, {1200, 1140, -50, 0x10, "rg2"},
                // Both forward.
                {1300, 1400, 110, 0, "rg1"},
        };
        final List<SAMRecord> records = new ArrayList<>();
        for (int index = 0; index < pairs.length; index++) {
            final Object[] pair = pairs[index];
            final int flags = (int) pair[3];
            final int start = (int) pair[0];
            final int mate = (int) pair[1];
            final int insert = (int) pair[2];
            for (final boolean first : new boolean[] {true, false}) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName("p" + index);
                final boolean reverse = (flags & 0x10) != 0;
                final boolean mateReverse = (flags & 0x20) != 0;
                record.setFlags(0x1 | 0x2 | (first ? 0x40 : 0x80) | (flags & 0x400)
                        | ((first ? reverse : mateReverse) ? 0x10 : 0)
                        | ((first ? mateReverse : reverse) ? 0x20 : 0));
                record.setReferenceIndex(0);
                record.setAlignmentStart(first ? start : mate);
                record.setMateReferenceIndex(0);
                record.setMateAlignmentStart(first ? mate : start);
                record.setInferredInsertSize(first ? insert : -insert);
                record.setCigarString("10M");
                record.setMappingQuality(60);
                record.setReadString("ACGTACGTAC");
                record.setBaseQualityString("IIIIIIIIII");
                if (pair[4] != null) {
                    record.setAttribute("RG", pair[4]);
                }
                records.add(record);
            }
        }
        if (away) {
            final SAMRecord other = new SAMRecord(header);
            other.setReadName("away");
            other.setFlags(0x1 | 0x2 | 0x40 | 0x20);
            other.setReferenceIndex(1);
            other.setAlignmentStart(50);
            other.setMateReferenceIndex(0);
            other.setMateAlignmentStart(50);
            other.setInferredInsertSize(80);
            other.setCigarString("10M");
            other.setMappingQuality(60);
            other.setReadString("ACGTACGTAC");
            other.setBaseQualityString("IIIIIIIIII");
            other.setAttribute("RG", "rg1");
            records.add(other);
        }
        records.sort(java.util.Comparator.comparingInt(SAMRecord::getReferenceIndex)
                .thenComparingInt(SAMRecord::getAlignmentStart));
        try (final SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(
                order == SAMFileHeader.SortOrder.coordinate)
                .makeBAMWriter(header, true, bam.toFile())) {
            records.forEach(writer::addAlignment);
        }
    }

    /**
     * Twelve reads in groups of the MarkDuplicates key: a1/a2, b1/b2 and c1/c2 share a start, d1 is
     * reverse and clipped, e1 sits on chr2 in a read group with no library, and u1 is unmapped.
     */
    static void duplicateFixture(final Path bam, final String changes, final String renamed,
                                 final String moved, final String groups) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 1000));
        dictionary.addSequence(new SAMSequenceRecord("chr2", 1000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord withLibrary = new SAMReadGroupRecord("rg1");
        withLibrary.setSample("sample1");
        withLibrary.setLibrary("lib1");
        header.addReadGroup(withLibrary);
        final SAMReadGroupRecord noLibrary = new SAMReadGroupRecord("rg2");
        noLibrary.setSample("sample1");
        header.addReadGroup(noLibrary);
        final SAMReadGroupRecord third = new SAMReadGroupRecord("rg3");
        third.setSample("sample1");
        third.setLibrary("lib3");
        if (groups.contains("rg3")) {
            header.addReadGroup(third);
        }
        final List<String> edits = Arrays.asList(changes.split(","));
        final Object[][] reads = {
                {"a1", 0, 100, 0x400, "10M"}, {"a2", 0, 100, 0, "10M"},
                {"b1", 0, 200, 0, "10M"}, {"b2", 0, 200, 0x400, "10M"},
                {"c1", 0, 300, 0, "10M"}, {"c2", 0, 300, 0, "10M"},
                {"d1", 0, 402, 0x10, "2S8M"},
                {"e1", 1, 50, 0, "10M"},
                {"u1", -1, 0, 0x4, "*"},
        };
        try (final SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true)
                .makeBAMWriter(header, false, bam.toFile())) {
            for (final Object[] read : reads) {
                final String name = (String) read[0];
                if (edits.contains("drop=" + name)) {
                    continue;
                }
                final SAMRecord record = new SAMRecord(header);
                record.setReadName(name.equals(renamed) ? "zz" : name);
                int flags = (int) read[3];
                if (edits.contains(name + "=dup")) {
                    flags |= 0x400;
                }
                if (edits.contains(name + "=none")) {
                    flags &= ~0x400;
                }
                record.setFlags(flags);
                if ((int) read[1] >= 0) {
                    record.setReferenceIndex((int) read[1]);
                    record.setAlignmentStart((int) read[2] + (name.equals(moved) ? 20 : 0));
                    record.setCigarString((String) read[4]);
                    record.setMappingQuality(name.equals("c2") ? 10 : 60);
                }
                record.setReadString("ACGTACGTAC");
                record.setBaseQualityString("IIIIIIIIII");
                final String group = name.equals("e1") ? "rg2" : "rg1";
                if (groups.equals(name + "=none")) {
                    // No group at all.
                } else if (groups.startsWith(name + "=")) {
                    record.setAttribute("RG", groups.substring(name.length() + 1));
                } else {
                    record.setAttribute("RG", group);
                }
                if (edits.contains(name + "=tag")) {
                    record.setAttribute("XX", 7);
                }
                writer.addAlignment(record);
            }
        }
    }

    /** Each run writes into a directory of its own, so whatever lands there is its output. */
    static String out(final String label, final String name) throws Exception {
        final Path dir = DIR.resolve("out").resolve(label);
        Files.createDirectories(dir);
        return dir.resolve(name).toString();
    }

    /** One run through Main; `output` says whether it writes, into `out/<label>/`. */
    static void run(final String label, final String output, final String tool,
                    final String... rest) throws Exception {
        final List<String> argv = new ArrayList<>();
        argv.add(tool);
        argv.addAll(Arrays.asList(rest));
        argv.addAll(Arrays.asList("--spark-master", "local[1]"));
        System.out.printf("argv\t%s\t%s%n", label,
                ReferenceQueryDump.escape(String.join("\t", argv)));
        // The flag never restores the JDK factory: GATK installs Intel's when it is false, the
        // setter is static, and every earlier run here left Intel's in place.
        if (argv.contains("--use-jdk-deflater")) {
            htsjdk.samtools.util.BlockCompressedOutputStream.setDefaultDeflaterFactory(
                    new htsjdk.samtools.util.zip.DeflaterFactory());
        }
        final ByteArrayOutputStream outBytes = new ByteArrayOutputStream();
        final PrintStream realOut = System.out;
        String failure = null;
        try {
            System.setOut(new PrintStream(outBytes, true, StandardCharsets.UTF_8));
            try {
                new Main().instanceMain(argv.toArray(new String[0]));
            } catch (final Exception e) {
                Throwable cause = e;
                while (cause.getCause() != null && cause.getCause() != cause) {
                    cause = cause.getCause();
                }
                failure = e.getClass().getName() + ":" + String.valueOf(cause.getMessage());
            }
        } finally {
            System.out.flush();
            System.setOut(realOut);
        }
        System.out.printf("stdout\t%s\t%s%n", label,
                ReferenceQueryDump.escape(masked(outBytes.toString(StandardCharsets.UTF_8))));
        if (failure != null) {
            System.out.printf("error\t%s\t%s%n", label, ReferenceQueryDump.escape(masked(failure)));
            return;
        }
        if (output == null) {
            return;
        }
        // Indexes and Hadoop's hidden checksums included: the directory is the answer.
        final Path outDir = DIR.resolve("out").resolve(label);
        final List<Path> produced;
        try (Stream<Path> walk = Files.walk(outDir)) {
            produced = walk.filter(Files::isRegularFile)
                    .sorted()
                    .collect(Collectors.toList());
        }
        for (final Path path : produced) {
            System.out.printf("file\t%s\t%s\t%s%n", label, outDir.relativize(path),
                    Base64.getEncoder().encodeToString(withoutTheTime(path)));
        }
    }

    /** A metrics file with the one line two runs cannot share taken out. */
    static byte[] withoutTheTime(final Path path) throws Exception {
        final byte[] bytes = Files.readAllBytes(path);
        if (!path.toString().endsWith(".metrics") && !path.toString().endsWith("Metrics.txt")) {
            return bytes;
        }
        return new String(bytes, StandardCharsets.UTF_8)
                .replaceAll("(?m)^# Started on: .*$", "# Started on: <start time>")
                .getBytes(StandardCharsets.UTF_8);
    }

    /** The path, and CompareDuplicatesSpark's summary in the one order two JVMs can share. */
    static String masked(final String text) {
        final String masked = text.replace(Path.of("").toAbsolutePath().toString(), "<cwd>");
        final String rule = "##############################\n";
        final int at = masked.indexOf(rule);
        if (at < 0) {
            return masked;
        }
        final int from = at + rule.length();
        return masked.substring(0, from) + Arrays.stream(masked.substring(from).split("\n"))
                .sorted().map(line -> line + "\n").collect(Collectors.joining());
    }
}

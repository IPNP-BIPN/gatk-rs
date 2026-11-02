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

        System.out.println("# SparkWalkersDump: thirty-seven Spark tools through Main, --spark-master local[1]");

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

        flowMarkingFixture(DIR.resolve("marking-flow.bam"));
        final String[] twoContigs = {"chr1", "chr2"};
        final String[][] twoGroups = {{"rg1", "lib1"}, {"rg4", "lib4"}};
        final String[][] noPrograms = {};
        final SAMFileHeader.SortOrder byName = SAMFileHeader.SortOrder.queryname;
        final SAMFileHeader.SortOrder byPlace = SAMFileHeader.SortOrder.coordinate;
        final SAMFileHeader.SortOrder noOrder = SAMFileHeader.SortOrder.unsorted;
        mergeFixture(DIR.resolve("more-queryname.bam"), byName, false, twoContigs, twoGroups, noPrograms, "from more", "rg4");
        mergeFixture(DIR.resolve("more-copy.bam"), byName, false, twoContigs, twoGroups, noPrograms, "from more", "rg4");
        mergeFixture(DIR.resolve("more-coordinate.bam"), byPlace, false, twoContigs, twoGroups, noPrograms, "from more", "rg4");
        mergeFixture(DIR.resolve("more-unsorted.bam"), noOrder, false, twoContigs, twoGroups, noPrograms, "from more", "rg4");
        mergeFixture(DIR.resolve("more-grouped.bam"), noOrder, true, twoContigs, twoGroups, noPrograms, "from more", "rg4");
        mergeFixture(DIR.resolve("more-scattered.bam"), noOrder, true, twoContigs, twoGroups, noPrograms, "from more", "rg4");
        // One id with two libraries, whose reads keep the id and so take the first file's library.
        mergeFixture(DIR.resolve("more-collision.bam"), byName, false, twoContigs,
                new String[][] {{"rg1", "other"}, {"rg2", "lib2"}}, noPrograms, null, "rg1");
        mergeFixture(DIR.resolve("more-programs-a.bam"), byName, false, twoContigs, twoGroups,
                new String[][] {{"pg1", "first", null}, {"pg2", "second", "pg1"}}, null, "rg4");
        mergeFixture(DIR.resolve("more-programs-b.bam"), byName, false, twoContigs, twoGroups,
                new String[][] {{"pg1", "other", null}, {"pg2", "second", "pg1"}, {"pg3", "third", "pg2"}}, null, "rg4");
        mergeFixture(DIR.resolve("more-programs-dangling.bam"), byName, false, twoContigs, twoGroups,
                new String[][] {{"pg1", "first", "gone"}}, null, "rg4");
        // A contig between the two others, which the merge places in front of chr2.
        mergeFixture(DIR.resolve("more-contig.bam"), byName, false, new String[] {"chr1", "chrM", "chr2"}, twoGroups, noPrograms, null, "rg4");
        mergeFixture(DIR.resolve("more-contig-place.bam"), byPlace, false, new String[] {"chr1", "chrM", "chr2"}, twoGroups, noPrograms, null, "rg4");
        mergeFixture(DIR.resolve("more-swapped.bam"), byName, false, new String[] {"chr2", "chr1"}, twoGroups, noPrograms, null, "rg4");
        // One input still goes through the merger: groups in id order, nothing else of its @HD.
        mergeFixture(DIR.resolve("more-header.bam"), byPlace, false, twoContigs,
                new String[][] {{"rg4", "lib4"}, {"rg1", "lib1"}},
                new String[][] {{"pg9", "last", null}, {"pg1", "first", null}}, "from more", "rg4");
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

        kmerFixtures();
        callerFixtures();
        filterFixtures();
        scoreFixtures();
        bwaFixtures(DIR);

        // The fixtures travel, so the port runs over these bytes rather than its own rebuild.
        svEvidenceBam(DIR.resolve("sv-evidence.bam"), 300, 1L, false);
        svEvidenceBam(DIR.resolve("sv-evidence2.bam"), 450, 2L, false);
        svEvidenceBam(DIR.resolve("sv-evidence-lone.bam"), 300, 3L, true);
        Files.write(DIR.resolve("sv-ignore.txt"), "alt\n".getBytes(StandardCharsets.UTF_8));
        Files.write(DIR.resolve("sv-ignore-bad.txt"), "nosuch\n".getBytes(StandardCharsets.UTF_8));
        svContigReads(DIR);
        for (final String name : new String[] {"caller.fasta", "caller.fasta.fai", "caller.dict",
                "caller.bam", "caller.bai", "caller-names.bam", "caller-dups.bam", "caller-dups.bai", "caller-dups-names.bam",
                "caller-known.vcf", "caller-known.vcf.idx", "caller-known.bed", "caller-no-groups.bam", "kmers.fasta", "kmers.fasta.fai", "kmers.dict",
                "kmers-big.fasta", "kmers-big.fasta.fai", "kmers-big.dict", "high-copy.fasta",
                "high-copy-blank.fasta", "high-copy.fasta.gz", "ref.fasta", "ref.fasta.fai", "ref.dict",
                "reads.bam", "reads.bai", "plain.vcf", "genotyped.vcf", "revert.map.tsv",
                "revert.bad-map.tsv", "names.txt", "dup-a.bam", "dup-a.bai", "dup-marks.bam",
                "dup-marks.bai", "dup-renamed.bam", "dup-uneven.bam", "dup-moved.bam",
                "dup-no-group.bam", "dup-other-group.bam", "dup-fewer.bam", "pairs.bam",
                "pairs.bai", "pairs-away.bam", "pairs-unsorted.bam", "marking.bam", "marking.bai",
                "marking-queryname.bam", "marking-unsorted.bam", "marking-one-library.bam",
                "marking-clean.bam", "marking-no-groups.bam", "marking-no-groups.bai", "marking-no-read-group.bam",
                "marking-unknown-read-group.bam", "marking-three.bam", "marking-collision.bam",
                "marking-scattered.bam", "marking-flow.bam", "more-queryname.bam", "more-copy.bam",
                "more-coordinate.bam", "more-coordinate.bai", "more-unsorted.bam", "more-grouped.bam",
                "more-scattered.bam",
                "more-collision.bam", "more-programs-a.bam", "more-programs-b.bam",
                "more-programs-dangling.bam", "more-contig.bam", "more-contig-place.bam",
                "more-contig-place.bai", "more-swapped.bam", "more-header.bam", "filter.bam",
                "filter-aligned.bam", "filter-empty.bam", "host.hss", "host.bfi", "host.txt", "taxa.db",
                "score-paired.bam", "score-unpaired.bam", "score-no-nm.bam", "score-bad-xa.bam",
                "score-bad-nm.bam", "score-bad-cigar.bam", "score-lonely.bam", "score-not-queryname.bam",
                "sv-evidence.bam", "sv-evidence.bai", "sv-evidence2.bam", "sv-evidence2.bai",
                "sv-evidence-lone.bam", "sv-evidence-lone.bai", "sv-ignore.txt", "sv-ignore-bad.txt",
                "sv_contigs.fasta", "sv_contigs.fasta.fai", "sv_contigs.dict", "sv_contigs.bam", "sv_contigs.bai",
                "sv_contigs_small.bam", "sv_contigs_small.bai", "sv_contigs_adjacent_reverse.bam",
                "sv_contigs_no_primary.bam", "sv_contigs_hard_primary.bam", "sv_contigs_no_sample.bam",
                "sv_contigs_two_samples.bam", "sv_cnv.vcf", "sv_cnv_other_sample.vcf", "sv_cnv_no_dictionary.vcf",
                "sv_cnv_other_dictionary.vcf", "sv_cnv_two_samples.vcf", "sv_truth.vcf", "sv_truth_near_start.vcf",
                "sv_truth_unknown_type.vcf", "sv_truth_unknown_contig.vcf", "sv_assembly.bam", "sv_assembly.bai",
                "sv_assembly_too_many.bam", "sv_assembly_overlap_contained.bam", "sv_assembly_no_breakend.bam",
                "sv_alt_tigs.txt", "sv_quality_missing.bam", "sv_cpx.vcf", "sv_cpx_empty.vcf",
                "sv_cpx_no_dictionary.vcf", "sv_cpx_no_dictionary_one.vcf", "sv_cpx_sparse.vcf",
                "sv_cpx_bad_segment.vcf", "sv_cpx_whole_contig.vcf", "sv_cpx_bad_length.vcf",
                "bwa.fasta", "bwa.fasta.fai", "bwa.dict", "bwa.fasta.img", "bwa_named.img", "bwa_not_image.img",
                "rp_known.vcf", "rp_known.vcf.idx", "rp_reads.bam", "rp_single.bam", "ps_microbe.img", "ps_microbe.dict", "ps_host.img", "ps_taxa.db", "bwa_no_image.fasta", "bwa_no_image.fasta.fai", "bwa_no_image.dict", "bwa_reads.bam",
                "bwa_reads_unsorted.bam", "bwa_reads_no_group.bam", "bwa_reads_dictionary.bam",
                "bwa_reads_other_dictionary.bam", "bwa_empty.bam", "bwa_iupac.bam", "bwa_no_sequence.bam",
                "bwa_odd.bam", "bwa_names.bam", "bwa_no_qualities.bam"}) {
            final byte[] bytes = Files.readAllBytes(DIR.resolve(name));
            System.out.printf("input\t%s\t%s%n", name, Base64.getEncoder().encodeToString(
                    name.endsWith(".img") && !name.equals("bwa_not_image.img") ? withoutAddresses(bytes) : bytes));
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
        final String flows = DIR.resolve("marking-flow.bam").toString();
        final String[][] flowRows = {
                {"mark-flow-plain"},
                {"mark-flowbased", "--flowbased", "true"},
                {"mark-flowbased-wide", "--flowbased", "true", "--flow-end-pos-uncertainty", "3"},
                {"mark-flowbased-narrow", "--flowbased", "true", "--flow-end-pos-uncertainty", "0"},
                {"mark-flow-end", "--single-end-reads-end-position-significant", "true"},
                {"mark-flow-end-wide", "--single-end-reads-end-position-significant", "true",
                        "--flow-end-pos-uncertainty", "2"},
                {"mark-flow-clipped", "--single-end-reads-clipping-is-end", "true"},
                {"mark-flow-sum", "--flow-quality-sum-strategy", "true"},
                {"mark-flow-sum-end", "--flow-quality-sum-strategy", "true",
                        "--single-end-reads-end-position-significant", "true"},
                {"mark-flow-skip", "--flow-skip-start-homopolymers", "1"},
                {"mark-flow-skip-two", "--flow-skip-start-homopolymers", "2"},
                {"mark-flow-skip-clipped", "--flow-skip-start-homopolymers", "2",
                        "--single-end-reads-clipping-is-end", "true"},
                {"mark-flow-q", "--single-end-reads-end-position-significant", "true",
                        "--flow-q-is-known-end", "true"},
                {"mark-flow-uncertainty-alone", "--flow-end-pos-uncertainty", "5"},
        };
        for (final String[] row : flowRows) {
            final List<String> argv = new ArrayList<>(Arrays.asList("-I", flows,
                    "--use-jdk-deflater", "true", "--duplicate-tagging-policy", "All",
                    "-O", out(row[0], "marked.bam"), "-M", out(row[0], "marked.metrics")));
            argv.addAll(Arrays.asList(row).subList(1, row.length));
            run(row[0], "marked.bam", marks, argv.toArray(new String[0]));
        }
        // SEVERAL INPUTS are a union, one partition an input, under the header their merge makes.
        final String allowed = "--allow-multiple-sort-orders-in-input";
        final String treat = "--treat-unsorted-as-querygroup-ordered";
        final String[][] severalRows = {
                {"mark-two", "marking-queryname.bam,more-queryname.bam"},
                {"mark-two-grouped", "more-grouped.bam,more-copy.bam"},
                {"mark-two-scattered", "more-scattered.bam,more-queryname.bam"},
                {"mark-two-same-path", "marking-queryname.bam,marking-queryname.bam"},
                {"mark-two-same-reads", "more-queryname.bam,more-copy.bam"},
                {"mark-two-coordinate", "marking.bam,more-coordinate.bam"},
                {"mark-two-second-bad", "marking-queryname.bam,more-coordinate.bam"},
                {"mark-two-unsorted", "marking-unsorted.bam,more-unsorted.bam"},
                {"mark-two-treated", "marking-unsorted.bam,more-unsorted.bam", treat, "true"},
                {"mark-two-treated-mixed", "marking-queryname.bam,more-unsorted.bam", treat, "true"},
                {"mark-two-allowed", "marking.bam,more-coordinate.bam", allowed, "true"},
                {"mark-two-allowed-mixed", "marking-queryname.bam,more-coordinate.bam", allowed, "true"},
                {"mark-two-allowed-grouped", "marking-queryname.bam,more-queryname.bam", allowed, "true"},
                {"mark-two-same-reads-sorted", "more-queryname.bam,more-copy.bam", allowed, "true"},
                {"mark-two-collision", "marking-queryname.bam,more-collision.bam"},
                {"mark-two-programs", "more-programs-a.bam,more-programs-b.bam"},
                {"mark-two-dangling", "more-queryname.bam,more-programs-dangling.bam"},
                {"mark-two-dictionaries", "marking-queryname.bam,more-contig.bam"},
                {"mark-two-dictionaries-sorted", "marking.bam,more-contig-place.bam", allowed, "true"},
                {"mark-two-dictionaries-swapped", "marking-queryname.bam,more-swapped.bam"},
                {"mark-two-interval", "marking.bam,more-coordinate.bam", allowed, "true", "-L", "chr1:100-300"},
                {"mark-two-interval-contig", "marking.bam,more-contig-place.bam", allowed, "true", "-L", "chrM:1-500"},
                {"mark-two-interval-no-index", "marking.bam,more-unsorted.bam", allowed, "true", "-L", "chr1:100-300"},
                {"mark-two-no-groups", "marking-no-groups.bam,more-coordinate.bam", allowed, "true"},
                {"mark-two-unknown-group", "marking-unknown-read-group.bam,more-coordinate.bam", allowed, "true"},
                {"mark-two-missing", "marking-queryname.bam,nowhere.bam"},
                {"mark-one-merged-header", "more-header.bam"},
        };
        for (final String[] row : severalRows) {
            final List<String> argv = new ArrayList<>();
            for (final String input : row[1].split(",")) {
                argv.addAll(Arrays.asList("-I", DIR.resolve(input).toString()));
            }
            argv.addAll(Arrays.asList("--use-jdk-deflater", "true", "--duplicate-tagging-policy", "All",
                    "-O", out(row[0], "marked.bam"), "-M", out(row[0], "marked.metrics")));
            argv.addAll(Arrays.asList(row).subList(2, row.length));
            run(row[0], "marked.bam", marks, argv.toArray(new String[0]));
        }
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
        // COLLECTALLELICCOUNTSSPARK writes a row per locus of the intervals, reads or none.
        final String allelic = "CollectAllelicCountsSpark";
        final String[][] allelicRows = {
                {"allelic-counts", "-L", "chr1:95-125"},
                {"allelic-counts-quality", "-L", "chr1:95-125", "--minimum-base-quality", "60"},
                {"allelic-counts-two", "-L", "chr1:95-125", "-L", "chr2:95-125"},
                {"allelic-counts-shards", "-L", "chr1:95-125", "--read-shard-size", "7"},
                {"allelic-counts-shuffle", "-L", "chr1:95-125", "--shuffle", "true"},
                {"allelic-counts-no-filters", "-L", "chr1:95-125", "--disable-tool-default-read-filters", "true"},
                {"allelic-counts-mapping", "-L", "chr1:95-125", "--minimum-mapping-quality", "0"},
                {"allelic-counts-away", "-L", "chr2:190-200"},
                {"allelic-counts-no-intervals"},
        };
        for (final String[] row : allelicRows) {
            final List<String> argv = new ArrayList<>(Arrays.asList("-I", i, "-R", r,
                    "-O", out(row[0], "allelic.tsv")));
            argv.addAll(Arrays.asList(row).subList(1, row.length));
            run(row[0], "allelic.tsv", allelic, argv.toArray(new String[0]));
        }

        kmerRows();
        callerRows();
        pipelineRows();
        filterRows();
        scoreRows();
        evidenceRows();
        discoveryRows();
        assemblyRows();
        reinterpreterRows();
        bwaRows();
        pathSeqBwaRows();

        // CALCMETADATASPARK reads everything, filters nothing, and names one partition.
        final String metadata = "CalcMetadataSpark";
        final String[][] metadataRows = {
                {"metadata-pairs", "pairs.bam"},
                {"metadata-pairs-binary", "pairs.bam", "--write-as-binary"},
                {"metadata-pairs-away", "pairs-away.bam"},
                {"metadata-pairs-unsorted", "pairs-unsorted.bam"},
                {"metadata-pairs-interval", "pairs.bam", "-L", "chr1:1-150"},
                {"metadata-marking", "marking.bam"},
                {"metadata-marking-binary", "marking.bam", "--write-as-binary"},
                {"metadata-marking-queryname", "marking-queryname.bam"},
                {"metadata-no-groups", "marking-no-groups.bam"},
                {"metadata-unknown-group", "marking-unknown-read-group.bam"},
                {"metadata-one-library", "marking-one-library.bam"},
                {"metadata-more-programs", "more-programs-a.bam"},
                {"metadata-reads", "reads.bam"},
                {"metadata-dup", "dup-a.bam"},
                {"metadata-away", "marking.bam", "-L", "chr2:190-200"},
                {"metadata-filter-ignored", "pairs.bam", "--read-filter", "ReadLengthReadFilter",
                        "--min-read-length", "1000", "--max-read-length", "2000"},
                {"metadata-two-inputs", "pairs.bam", "-I", DIR.resolve("pairs-away.bam").toString()},
                {"metadata-missing", "nowhere.bam"},
        };
        for (final String[] row : metadataRows) {
            final List<String> argv = new ArrayList<>(Arrays.asList(
                    "-I", DIR.resolve(row[1]).toString(), "-O", out(row[0], "metadata.txt")));
            argv.addAll(Arrays.asList(row).subList(2, row.length));
            run(row[0], "metadata.txt", metadata, argv.toArray(new String[0]));
        }
        run("metadata-no-directory", null, metadata, "-I", DIR.resolve("pairs.bam").toString(),
                "-O", DIR.resolve("nowhere/metadata.txt").toString());
        run("metadata-binary-no-directory", null, metadata, "-I", DIR.resolve("pairs.bam").toString(),
                "-O", DIR.resolve("nowhere/metadata.txt").toString(), "--write-as-binary");
        run("allelic-counts-no-reference", null, allelic, "-I", i, "-L", "chr1:95-125",
                "-O", out("allelic-counts-no-reference", "allelic.tsv"));
        // marking.bam's contigs are longer than the reference's, so these turn validation off.
        // They are the rows with reads under the loci: reads.bam's intervals above come back empty.
        final String[][] markedRows = {
                {"allelic-counts-other-dictionary", "marking.bam"},
                {"allelic-marked-quality", "marking.bam", "--minimum-base-quality", "60"},
                {"allelic-marked-no-quality", "marking.bam", "--minimum-base-quality", "0"},
                {"allelic-marked-two", "marking.bam", "-L", "chr2:95-125"},
                {"allelic-marked-shards", "marking.bam", "--read-shard-size", "7"},
                {"allelic-marked-shuffle", "marking.bam", "--shuffle", "true"},
                {"allelic-marked-shuffle-shards", "marking.bam", "--shuffle", "true", "--read-shard-size", "7"},
                {"allelic-marked-no-filters", "marking.bam", "--disable-tool-default-read-filters", "true"},
                {"allelic-marked-mapping", "marking.bam", "--minimum-mapping-quality", "0"},
                {"allelic-marked-depth", "marking.bam", "--max-depth-per-sample", "1"},
                {"allelic-counts-no-reads", "marking.bam", "--minimum-mapping-quality", "100"},
                {"allelic-counts-no-reads-shuffle", "marking.bam", "--shuffle", "true",
                        "--minimum-mapping-quality", "100"},
                {"allelic-counts-no-sample", "marking-no-groups.bam"},
                {"allelic-counts-no-sample-unfiltered", "marking-no-groups.bam",
                        "--disable-tool-default-read-filters", "true"},
        };
        for (final String[] row : markedRows) {
            final List<String> argv = new ArrayList<>(Arrays.asList("-I", DIR.resolve(row[1]).toString(),
                    "-R", r, "-L", "chr1:95-125", "--disable-sequence-dictionary-validation", "true",
                    "-O", out(row[0], "allelic.tsv")));
            argv.addAll(Arrays.asList(row).subList(2, row.length));
            run(row[0], "allelic.tsv", allelic, argv.toArray(new String[0]));
        }
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
    /** A few reads over the given contigs, with what the header merger has to reconcile. */
    static void mergeFixture(final Path bam, final SAMFileHeader.SortOrder order, final boolean grouped,
                             final String[] contigs, final String[][] groups, final String[][] programs,
                             final String comment, final String group) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        for (final String contig : contigs) {
            dictionary.addSequence(new SAMSequenceRecord(contig, 100000));
        }
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(order);
        if (grouped) {
            header.setGroupOrder(SAMFileHeader.GroupOrder.query);
        }
        for (final String[] entry : groups) {
            final SAMReadGroupRecord record = new SAMReadGroupRecord(entry[0]);
            record.setSample("sample1");
            record.setLibrary(entry[1]);
            header.addReadGroup(record);
        }
        for (final String[] entry : programs) {
            final htsjdk.samtools.SAMProgramRecord record = new htsjdk.samtools.SAMProgramRecord(entry[0]);
            record.setProgramName(entry[1]);
            record.setPreviousProgramGroupId(entry[2]);
            header.addProgramRecord(record);
        }
        if (comment != null) {
            header.addComment(comment);
        }
        final int chr1 = dictionary.getSequenceIndex("chr1");
        final int chr2 = dictionary.getSequenceIndex("chr2");
        final int last = contigs.length > 2 ? 1 : chr2;
        // name, contig, start, flags, mate contig, mate start, quality, group
        final Object[][] reads = {
                // At the place of marking.bam's first pair, in its library when the group is rg1.
                {"m:1:1101:1:1", chr1, 100, 0x63, chr1, 200, 50, group},
                {"m:1:1101:1:1", chr1, 200, 0x93, chr1, 100, 50, group},
                {"m:1:1101:1000:1000", chr1, 100, 0x63, chr1, 200, 15, "rg1"},
                {"m:1:1101:1000:1000", chr1, 200, 0x93, chr1, 100, 15, "rg1"},
                // A name marking.bam also has, with reads of its own.
                {"f:1:1101:1:1", chr1, 4000, 0, -1, 0, 37, "rg1"},
                {"m:1:1101:7:7", last, 300, 0x63, last, 400, 30, group},
                {"m:1:1101:7:7", last, 400, 0x93, last, 300, 30, group},
                {"m:1:1101:8:8", last, 300, 0x63, last, 400, 20, group},
                {"m:1:1101:8:8", last, 400, 0x93, last, 300, 20, group},
                {"m:1:1101:9:9", chr2, 500, 0, -1, 0, 30, group},
        };
        // One name that comes back under a header that says grouped: two fragments, each best at
        // its place, so one name is answered twice and only that name.
        final Object[][] back = {
                {"z:1:1101:1:1", chr1, 8000, 0, -1, 0, 30, group},
                {"z:1:1101:2:2", chr1, 8500, 0, -1, 0, 30, group},
                {"z:1:1101:1:1", chr1, 9000, 0, -1, 0, 30, group},
        };
        final List<Object[]> all = new ArrayList<>(Arrays.asList(reads));
        if (bam.getFileName().toString().contains("scattered")) {
            all.addAll(Arrays.asList(back));
        }
        final List<SAMRecord> records = new ArrayList<>();
        for (final Object[] read : all) {
            final SAMRecord record = new SAMRecord(header);
            record.setReadName((String) read[0]);
            record.setFlags((int) read[3]);
            record.setReferenceIndex((int) read[1]);
            record.setAlignmentStart((int) read[2]);
            record.setCigarString("10M");
            record.setMateReferenceIndex((int) read[4]);
            record.setMateAlignmentStart((int) read[5]);
            record.setMappingQuality(60);
            record.setReadString("ACGTACGTAC");
            final char[] qualities = new char[10];
            Arrays.fill(qualities, (char) (33 + (int) read[6]));
            record.setBaseQualityString(new String(qualities));
            record.setAttribute("RG", read[7]);
            records.add(record);
        }
        final SAMFileWriterFactory factory = new SAMFileWriterFactory()
                .setCreateIndex(order == SAMFileHeader.SortOrder.coordinate);
        try (SAMFileWriter writer = factory.makeBAMWriter(header, false, bam)) {
            records.forEach(writer::addAlignment);
        }
    }

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

    /**
     * Fragments at one start with the ends flow mode tells apart: ends a base or two apart, clips on
     * either side that a clipped location sees, the tm tag's A, Q and Z, homopolymers to skip, and a
     * read in a group with no flow order.
     */
    static void flowMarkingFixture(final Path bam) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 100000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        final SAMReadGroupRecord flow = new SAMReadGroupRecord("flow");
        flow.setSample("sample1");
        flow.setLibrary("lib1");
        flow.setAttribute("FO", "TGCA");
        header.addReadGroup(flow);
        final SAMReadGroupRecord plain = new SAMReadGroupRecord("plain");
        plain.setSample("sample1");
        plain.setLibrary("lib1");
        header.addReadGroup(plain);
        // name, start, flags, cigar, bases, quality, tm, flow tags, group
        final Object[][] reads = {
                {"a:1:1:1:1", 100, 0, "12M", "AACCGGTTAACC", 30, null, true, "flow"},
                {"a:1:1:2:2", 100, 0, "13M", "AACCGGTTAACCG", 31, null, true, "flow"},
                {"a:1:1:3:3", 100, 0, "14M", "AACCGGTTAACCGG", 32, null, false, "flow"},
                {"a:1:1:4:4", 100, 0, "20M", "AACCGGTTAACCGGTTAACC", 20, null, true, "flow"},
                {"a:1:1:5:5", 100, 0, "20M", "AACCGGTTAACCGGTTAACC", 25, "Q", true, "flow"},
                {"a:1:1:6:6", 100, 0, "20M", "AACCGGTTAACCGGTTAACC", 26, "Z", true, "flow"},
                {"b:1:1:1:1", 200, 0, "2S12M", "TTAACCGGTTAACC", 30, null, true, "flow"},
                {"b:1:1:2:2", 202, 0, "12M2S", "AACCGGTTAACCGG", 31, null, true, "flow"},
                {"b:1:1:3:3", 202, 0, "12M", "AACCGGTTAACC", 32, "A", true, "flow"},
                {"b:1:1:4:4", 200, 0, "2S12M", "TTAACCGGTTAACC", 33, "A", true, "flow"},
                {"c:1:1:1:1", 300, 0x10, "12M", "AACCGGTTAACC", 30, null, true, "flow"},
                {"c:1:1:2:2", 301, 0x10, "11M", "AACCGGTTAAC", 31, null, true, "flow"},
                {"c:1:1:3:3", 296, 0x10, "4S12M", "AACCAACCGGTTAACC", 32, null, true, "flow"},
                {"d:1:1:1:1", 400, 0, "12M", "AAAACCGGTTAA", 30, null, true, "flow"},
                {"d:1:1:2:2", 402, 0, "10M", "AACCGGTTAA", 31, null, true, "flow"},
                {"d:1:1:3:3", 404, 0, "8M", "CCGGTTAA", 32, null, true, "flow"},
                {"e:1:1:1:1", 500, 0, "12M", "TTGGCCAATTGG", 30, null, false, "plain"},
                {"e:1:1:2:2", 501, 0, "11M", "TGGCCAATTGG", 31, null, false, "plain"},
        };
        final List<SAMRecord> records = new ArrayList<>();
        for (final Object[] read : reads) {
            final SAMRecord record = new SAMRecord(header);
            record.setReadName((String) read[0]);
            record.setFlags((int) read[2]);
            record.setReferenceIndex(0);
            record.setAlignmentStart((int) read[1]);
            record.setCigarString((String) read[3]);
            record.setMappingQuality(60);
            final String bases = (String) read[4];
            record.setReadString(bases);
            final char[] qualities = new char[bases.length()];
            for (int index = 0; index < qualities.length; index++) {
                qualities[index] = (char) (33 + (int) read[5] + (index % 3));
            }
            record.setBaseQualityString(new String(qualities));
            if (read[6] != null) {
                record.setAttribute("tm", read[6]);
            }
            if ((boolean) read[7]) {
                record.setAttribute("tp", new byte[qualities.length]);
            }
            record.setAttribute("RG", read[8]);
            records.add(record);
        }
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .makeBAMWriter(header, false, bam.toFile())) {
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

    /**
     * BwaSpark: BWA through the JNI, the supplementary lines a paired run filters away, and the
     * refusals of the engine. Its pipeline marks duplicates over the same alignment.
     */
    static void bwaRows() throws Exception {
        final String reference = DIR.resolve("bwa.fasta").toString();
        final String[][] runs = {
                {"bwa-paired", "BwaSpark", "bwa_reads.bam"},
                {"bwa-single", "BwaSpark", "bwa_reads.bam", "--single-end-alignment"},
                {"bwa-named-image", "BwaSpark", "bwa_reads.bam", "--bwa-mem-index-image", DIR.resolve("bwa_named.img").toString()},
                {"bwa-image-missing", "BwaSpark", "bwa_reads.bam", "--bwa-mem-index-image", DIR.resolve("nowhere.img").toString()},
                {"bwa-unsorted", "BwaSpark", "bwa_reads_unsorted.bam"},
                {"bwa-unsorted-single", "BwaSpark", "bwa_reads_unsorted.bam", "--single-end-alignment"},
                {"bwa-no-group", "BwaSpark", "bwa_reads_no_group.bam"},
                {"bwa-dictionary", "BwaSpark", "bwa_reads_dictionary.bam"},
                {"bwa-other-dictionary", "BwaSpark", "bwa_reads_other_dictionary.bam"},
                {"bwa-empty", "BwaSpark", "bwa_empty.bam"},
                {"bwa-empty-single", "BwaSpark", "bwa_empty.bam", "--single-end-alignment"},
                {"bwa-iupac", "BwaSpark", "bwa_iupac.bam"},
                {"bwa-iupac-single", "BwaSpark", "bwa_iupac.bam", "--single-end-alignment"},
                {"bwa-no-sequence", "BwaSpark", "bwa_no_sequence.bam"},
                {"bwa-no-sequence-single", "BwaSpark", "bwa_no_sequence.bam", "--single-end-alignment"},
                {"bwa-odd", "BwaSpark", "bwa_odd.bam"},
                {"bwa-odd-single", "BwaSpark", "bwa_odd.bam", "--single-end-alignment"},
                {"bwa-names", "BwaSpark", "bwa_names.bam"},
                {"bwa-no-qualities", "BwaSpark", "bwa_no_qualities.bam"},
                {"bwa-no-qualities-single", "BwaSpark", "bwa_no_qualities.bam", "--single-end-alignment"},
                {"bwa-filter", "BwaSpark", "bwa_reads.bam", "--read-filter", "MappingQualityReadFilter",
                        "--minimum-mapping-quality", "30"},
                {"bwa-no-filter", "BwaSpark", "bwa_reads.bam", "--disable-tool-default-read-filters", "true"},
                {"bwa-pipeline", "BwaAndMarkDuplicatesPipelineSpark", "bwa_reads.bam"},
                {"bwa-pipeline-tagging", "BwaAndMarkDuplicatesPipelineSpark", "bwa_reads.bam", "--duplicate-tagging-policy", "All"},
                {"bwa-pipeline-remove", "BwaAndMarkDuplicatesPipelineSpark", "bwa_reads.bam", "--remove-all-duplicates", "true"},
                {"bwa-pipeline-length", "BwaAndMarkDuplicatesPipelineSpark", "bwa_reads.bam", "--duplicate-scoring-strategy",
                        "TOTAL_MAPPED_REFERENCE_LENGTH"},
                {"bwa-pipeline-mates", "BwaAndMarkDuplicatesPipelineSpark", "bwa_reads.bam", "--do-not-mark-unmapped-mates", "true"},
                {"bwa-pipeline-no-group", "BwaAndMarkDuplicatesPipelineSpark", "bwa_reads_no_group.bam"},
                {"bwa-pipeline-unsorted", "BwaAndMarkDuplicatesPipelineSpark", "bwa_reads_unsorted.bam"},
                {"bwa-pipeline-empty", "BwaAndMarkDuplicatesPipelineSpark", "bwa_empty.bam"},
        };
        for (final String[] one : runs) {
            final List<String> rest = new ArrayList<>(Arrays.asList(
                    "-I", DIR.resolve(one[2]).toString(), "-R", reference,
                    "--use-jdk-deflater", "true", "-O", out(one[0], "aligned.bam")));
            rest.addAll(Arrays.asList(one).subList(3, one.length));
            run(one[0], "aligned.bam", one[1], rest.toArray(new String[0]));
        }
        run("bwa-default-image-missing", "aligned.bam", "BwaSpark", "-I", DIR.resolve("bwa_reads.bam").toString(),
                "-R", DIR.resolve("bwa_no_image.fasta").toString(), "--use-jdk-deflater", "true",
                "-O", out("bwa-default-image-missing", "aligned.bam"));
        run("bwa-sam", "aligned.sam", "BwaSpark", "-I", DIR.resolve("bwa_reads.bam").toString(), "-R", reference,
                "-O", out("bwa-sam", "aligned.sam"));
        run("bwa-sharded", "aligned.bam", "BwaSpark", "-I", DIR.resolve("bwa_reads.bam").toString(), "-R", reference,
                "--use-jdk-deflater", "true", "--sharded-output", "true", "-O", out("bwa-sharded", "aligned.bam"));
    }
    /**
     * A BWA image without the addresses `bwa_idx2mem` leaves in it, which ASLR changes on every run
     * and nothing reads back: `bwt_t.sa`, then past the BWT and the suffix array `bntseq_t.anns`,
     * `.ambs` and `.fp_pac`, and each `bntann1_t`'s `name` and `anno`. Only the printed row is
     * cleaned; the reference ran on the file as BWA wrote it.
     */
    private static byte[] withoutAddresses(final byte[] image) {
        final byte[] cleaned = image.clone();
        final java.nio.ByteBuffer read = java.nio.ByteBuffer.wrap(image).order(java.nio.ByteOrder.LITTLE_ENDIAN);
        final int bwt = Math.toIntExact(1120 + read.getLong(56) * 4 + read.getLong(1104) * 8);
        final List<Integer> addresses = new ArrayList<>(List.of(1112, bwt + 16, bwt + 32, bwt + 40));
        final int annotations = bwt + 48 + 16 * read.getInt(bwt + 24);
        for (int sequence = 0; sequence < read.getInt(bwt + 8); sequence++) {
            addresses.add(annotations + 40 * sequence + 24);
            addresses.add(annotations + 40 * sequence + 32);
        }
        addresses.forEach(address -> Arrays.fill(cleaned, address, address + 8, (byte) 0));
        return cleaned;
    }

    /**
     * PathSeq with BWA: the aligner alone, the filter's host image, and the whole pipeline, whose
     * batches a small partition size makes several of.
     */
    static void pathSeqBwaRows() throws Exception {
        final String[] microbe = {"--microbe-bwa-image", DIR.resolve("ps_microbe.img").toString(),
                "--microbe-dict", DIR.resolve("ps_microbe.dict").toString()};
        final String reads = DIR.resolve("bwa_reads.bam").toString();
        final String unsorted = DIR.resolve("bwa_reads_unsorted.bam").toString();
        final String[][] aligner = {
                {"psb-paired", "--paired-input", reads},
                {"psb-unpaired", "--unpaired-input", unsorted},
                {"psb-unpaired-queryname", "--unpaired-input", reads},
                {"psb-both", "--paired-input", reads, "--unpaired-input", unsorted},
                {"psb-other-missing", "--paired-input", reads, "--unpaired-input", DIR.resolve("nowhere.bam").toString()},
                {"psb-seed", "--paired-input", reads, "--microbe-min-seed-length", "40"},
                {"psb-seed-short", "--unpaired-input", unsorted, "--microbe-min-seed-length", "11"},
                {"psb-score", "--paired-input", reads, "--bwa-score-threshold", "90"},
                {"psb-score-zero", "--unpaired-input", unsorted, "--bwa-score-threshold", "0"},
                {"psb-hits-zero", "--paired-input", reads, "--max-alternate-hits", "0"},
                {"psb-hits", "--unpaired-input", unsorted, "--max-alternate-hits", "1"},
                {"psb-empty", "--paired-input", DIR.resolve("bwa_empty.bam").toString()},
                {"psb-odd", "--paired-input", DIR.resolve("bwa_odd.bam").toString()},
                {"psb-odd-unpaired", "--unpaired-input", DIR.resolve("bwa_odd.bam").toString()},
                {"psb-no-qualities", "--paired-input", DIR.resolve("bwa_no_qualities.bam").toString()},
                {"psb-iupac", "--paired-input", DIR.resolve("bwa_iupac.bam").toString()},
                {"psb-not-queryname", "--paired-input", unsorted},
                {"psb-aligned", "--paired-input", DIR.resolve("bwa_reads_dictionary.bam").toString()},
                {"psb-input", "-I", reads},
                {"psb-no-inputs"},
                {"psb-all-missing", "--paired-input", DIR.resolve("nowhere.bam").toString()},
        };
        for (final String[] one : aligner) {
            final String name = one[0];
            final List<String> rest = new ArrayList<>(Arrays.asList(microbe));
            rest.addAll(Arrays.asList("--use-jdk-deflater", "true", "--paired-output", out(name, name + "-paired.bam"),
                    "--unpaired-output", out(name, name + "-unpaired.bam")));
            rest.addAll(Arrays.asList(one).subList(1, one.length));
            run(name, name + "-paired.bam", "PathSeqBwaSpark", rest.toArray(new String[0]));
        }
        run("psb-no-output", null, "PathSeqBwaSpark", microbe[0], microbe[1], microbe[2], microbe[3], "--paired-input", reads);
        run("psb-sam-output", null, "PathSeqBwaSpark", microbe[0], microbe[1], microbe[2], microbe[3], "--paired-input", reads,
                "--paired-output", out("psb-sam-output", "p.sam"));
        run("psb-image-missing", null, "PathSeqBwaSpark", microbe[0], DIR.resolve("nowhere.img").toString(), microbe[2], microbe[3],
                "--paired-input", reads, "--paired-output", out("psb-image-missing", "p.bam"));
        run("psb-dict-missing", null, "PathSeqBwaSpark", microbe[0], microbe[1], microbe[2], DIR.resolve("nowhere.dict").toString(),
                "--paired-input", reads, "--paired-output", out("psb-dict-missing", "p.bam"));

        final String[][] filter = {
                {"psf-host"},
                {"psf-host-skip", "--skip-quality-filters", "true"},
                {"psf-host-seed", "--filter-bwa-seed-length", "60"},
                {"psf-host-identity", "--host-min-identity", "95"},
                {"psf-host-no-duplicates", "--filter-duplicates", "false"},
                {"psf-host-missing", "--filter-bwa-image", DIR.resolve("nowhere.img").toString()},
        };
        for (final String[] one : filter) {
            final String name = one[0];
            final List<String> rest = new ArrayList<>(Arrays.asList("-I", reads, "--use-jdk-deflater", "true",
                    "--paired-output", out(name, name + "-paired.bam"), "--unpaired-output", out(name, name + "-unpaired.bam"),
                    "--filter-metrics", out(name, name + "Metrics.txt")));
            if (one.length < 2 || !one[1].equals("--filter-bwa-image")) {
                rest.addAll(Arrays.asList("--filter-bwa-image", DIR.resolve("ps_host.img").toString()));
            }
            rest.addAll(Arrays.asList(one).subList(1, one.length));
            run(name, name + "-paired.bam", "PathSeqFilterSpark", rest.toArray(new String[0]));
        }

        final String db = DIR.resolve("ps_taxa.db").toString();
        final String[][] pipeline = {
                {"psp-default", reads},
                {"psp-skip", reads, "--skip-quality-filters", "true"},
                {"psp-no-duplicates", reads, "--filter-duplicates", "false", "--skip-quality-filters", "true"},
                {"psp-host", reads, "--filter-bwa-image", DIR.resolve("ps_host.img").toString()},
                {"psp-unsorted", unsorted},
                {"psp-batches", reads, "--pipeline-reads-per-partition", "100", "--skip-quality-filters", "true"},
                {"psp-batches-three", reads, "--pipeline-reads-per-partition", "256", "--skip-quality-filters", "true"},
                {"psp-output-partition", reads, "--readsPerPartitionOutput", "100"},
                {"psp-seed", reads, "--microbe-min-seed-length", "40"},
                {"psp-score", reads, "--bwa-score-threshold", "90"},
                {"psp-hits", reads, "--max-alternate-hits", "0"},
                {"psp-identity", reads, "--min-score-identity", "0.99"},
                {"psp-margin", reads, "--identity-margin", "0.3"},
                {"psp-genome-length", reads, "--divide-by-genome-length", "true"},
                {"psp-no-kingdom", reads, "--not-normalized-by-kingdom", "true"},
                {"psp-empty", DIR.resolve("bwa_empty.bam").toString()},
                {"psp-odd", DIR.resolve("bwa_odd.bam").toString()},
                {"psp-reducers", reads, "--num-reducers", "2"},
        };
        for (final String[] one : pipeline) {
            final String name = one[0];
            final List<String> rest = new ArrayList<>(Arrays.asList(microbe));
            rest.addAll(Arrays.asList("-I", one[1], "--taxonomy-file", db, "--use-jdk-deflater", "true",
                    "--scores-output", out(name, name + ".txt"), "--output", out(name, name + ".bam"),
                    "--filter-metrics", out(name, name + "FilterMetrics.txt"),
                    "--score-metrics", out(name, name + "ScoreMetrics.txt"),
                    "--score-warnings", out(name, name + "-warnings.txt")));
            rest.addAll(Arrays.asList(one).subList(2, one.length));
            run(name, name + ".txt", "PathSeqPipelineSpark", rest.toArray(new String[0]));
        }
        final String copy = "ParallelCopyGCSDirectoryIntoHDFSSpark";
        run("copy-local-input", null, copy, "--input-gcs-path", DIR.resolve("bwa.fasta").toString(),
                "--output-hdfs-directory", "hdfs://nowhere:1/out");
        run("copy-local-output", null, copy, "--input-gcs-path", "gs://nowhere/in",
                "--output-hdfs-directory", DIR.resolve("out").toString());
        final String bwaReference = DIR.resolve("bwa.fasta").toString();
        final String platformed = DIR.resolve("rp_reads.bam").toString();
        final String[][] readsPipeline = {
                {"rp-align", platformed},
                {"rp-align-no-platform", reads},
                {"rp-align-single", DIR.resolve("rp_single.bam").toString(), "--single-end-alignment"},
                {"rp-align-dictionary", DIR.resolve("bwa_reads_dictionary.bam").toString()},
                {"rp-align-unsorted", unsorted},
                {"rp-align-named", platformed, "--bwa-mem-index-image", DIR.resolve("bwa_named.img").toString()},
                {"rp-align-image-missing", platformed, "--bwa-mem-index-image", DIR.resolve("nowhere.img").toString()},
                {"rp-align-no-group", DIR.resolve("bwa_reads_no_group.bam").toString()},
                {"rp-align-empty", DIR.resolve("bwa_empty.bam").toString()},
        };
        for (final String[] one : readsPipeline) {
            final String name = one[0];
            final List<String> rest = new ArrayList<>(Arrays.asList("-I", one[1], "-R", bwaReference, "--align", "true",
                    "--known-sites", DIR.resolve("rp_known.vcf").toString(), "--use-jdk-deflater", "true",
                    "-O", out(name, name + ".vcf"), "--output-bam", out(name, name + ".bam")));
            rest.addAll(Arrays.asList(one).subList(2, one.length));
            if (rest.contains("--single-end-alignment")) {
                rest.add("true");
            }
            run(name, name + ".vcf", "ReadsPipelineSpark", rest.toArray(new String[0]));
        }
        final String[][] missing = {
                {"psp-image-missing", DIR.resolve("nowhere.img").toString(), microbe[3], db},
                {"psp-dict-missing", microbe[1], DIR.resolve("nowhere.dict").toString(), db},
                {"psp-db-missing", microbe[1], microbe[3], DIR.resolve("nowhere.db").toString()},
        };
        for (final String[] one : missing) {
            run(one[0], null, "PathSeqPipelineSpark", microbe[0], one[1], microbe[2], one[2], "-I", reads,
                    "--taxonomy-file", one[3], "--scores-output", out(one[0], "s.txt"));
        }
        run("psp-only-scores", "s.txt", "PathSeqPipelineSpark", microbe[0], microbe[1], microbe[2], microbe[3], "-I", reads,
                "--taxonomy-file", db, "--scores-output", out("psp-only-scores", "s.txt"));
    }

    /**
     * BwaSpark and its pipeline: a reference with its image, whose second contig repeats pieces of
     * the first and holds a run of N, and reads as a sequencer would leave them, queryname sorted
     * and unaligned. Some pairs are the same fragment twice, for the duplicate marking, and some
     * are chimeras, junk or too short to place.
     */
    static void bwaFixtures(final Path dir) throws Exception {
        final java.util.Random random = new java.util.Random(77L);
        final java.util.function.IntFunction<String> bases = length -> {
            final StringBuilder made = new StringBuilder(length);
            for (int i = 0; i != length; ++i) {
                made.append("ACGT".charAt(random.nextInt(4)));
            }
            return made.toString();
        };
        final String first = bases.apply(40000);
        final String second = bases.apply(6000) + first.substring(1000, 3500) + "N".repeat(50) + bases.apply(4000)
                + "R" + first.substring(9000, 11000) + bases.apply(5000);
        final java.util.function.BiFunction<String, String, String> contig = (name, sequence) -> {
            final StringBuilder lines = new StringBuilder(">" + name + "\n");
            for (int at = 0; at < sequence.length(); at += 60) {
                lines.append(sequence, at, Math.min(at + 60, sequence.length())).append('\n');
            }
            return lines.toString();
        };
        final Path fasta = dir.resolve("bwa.fasta");
        Files.writeString(fasta, contig.apply("bwa1", first) + contig.apply("bwa2", second), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        Files.writeString(dir.resolve("bwa.dict"), "@HD\tVN:1.6\n@SQ\tSN:bwa1\tLN:" + first.length()
                + "\n@SQ\tSN:bwa2\tLN:" + second.length() + "\n", StandardCharsets.UTF_8);
        new org.broadinstitute.hellbender.Main().instanceMain(new String[] {"BwaMemIndexImageCreator",
                "-I", fasta.toString(), "-O", dir.resolve("bwa.fasta.img").toString()});
        Files.copy(dir.resolve("bwa.fasta.img"), dir.resolve("bwa_named.img"),
                java.nio.file.StandardCopyOption.REPLACE_EXISTING);
        Files.writeString(dir.resolve("bwa_not_image.img"), "not an image\n", StandardCharsets.UTF_8);
        Files.copy(fasta, dir.resolve("bwa_no_image.fasta"), java.nio.file.StandardCopyOption.REPLACE_EXISTING);
        Files.copy(dir.resolve("bwa.fasta.fai"), dir.resolve("bwa_no_image.fasta.fai"),
                java.nio.file.StandardCopyOption.REPLACE_EXISTING);
        Files.copy(dir.resolve("bwa.dict"), dir.resolve("bwa_no_image.dict"),
                java.nio.file.StandardCopyOption.REPLACE_EXISTING);

        final java.util.function.UnaryOperator<String> flip = sequence -> {
            final StringBuilder flipped = new StringBuilder(sequence.length());
            for (int at = sequence.length() - 1; at >= 0; --at) {
                flipped.append("TGCAN".charAt(Math.max("ACGTN".indexOf(sequence.charAt(at)), 0)));
            }
            return flipped.toString();
        };
        final java.util.function.UnaryOperator<String> noisy = sequence -> {
            final StringBuilder changed = new StringBuilder(sequence);
            for (int at = 0; at != changed.length(); ++at) {
                if (random.nextInt(60) == 0) {
                    changed.setCharAt(at, "ACGT".charAt(random.nextInt(4)));
                }
            }
            return changed.toString();
        };
        final List<String[]> pairs = new ArrayList<>();
        String[] last = null;
        for (int i = 0; i != 300; ++i) {
            final String template = random.nextInt(3) == 0 ? second.replace('N', 'A').replace('R', 'A') : first;
            final int length = 200 + random.nextInt(300);
            final int start = random.nextInt(template.length() - length);
            final String fragment = template.substring(start, start + length);
            final String[] pair;
            final int kind = random.nextInt(20);
            if (kind == 0) {
                pair = new String[] {bases.apply(101), bases.apply(101)};
            } else if (kind == 1) {
                pair = new String[] {noisy.apply(fragment.substring(0, 101)), bases.apply(101)};
            } else if (kind == 2) {
                final int other = random.nextInt(first.length() - 60);
                pair = new String[] {fragment.substring(0, 50) + flip.apply(first.substring(other, other + 60)),
                        flip.apply(fragment.substring(length - 101))};
            } else if (kind == 3) {
                pair = new String[] {fragment.substring(0, 15), flip.apply(fragment.substring(length - 18))};
            } else if (kind <= 6 && last != null) {
                pair = last;
            } else {
                pair = new String[] {noisy.apply(fragment.substring(0, 101)),
                        noisy.apply(flip.apply(fragment.substring(length - 101)))};
            }
            pairs.add(pair);
            last = pair;
        }
        bwaReadBam(dir.resolve("bwa_reads.bam"), pairs, SAMFileHeader.SortOrder.queryname, true, null);
        bwaReadBam(dir.resolve("bwa_reads_unsorted.bam"), pairs, SAMFileHeader.SortOrder.unsorted, true, null);
        bwaReadBam(dir.resolve("bwa_reads_no_group.bam"), pairs, SAMFileHeader.SortOrder.queryname, false, null);
        bwaReadBam(dir.resolve("bwa_reads_dictionary.bam"), pairs, SAMFileHeader.SortOrder.queryname, true,
                new String[] {"bwa1", String.valueOf(first.length()), "bwa2", String.valueOf(second.length())});
        bwaReadBam(dir.resolve("bwa_reads_other_dictionary.bam"), pairs, SAMFileHeader.SortOrder.queryname, true,
                new String[] {"bwa1", String.valueOf(first.length())});
        bwaReadBam(dir.resolve("bwa_empty.bam"), new ArrayList<>(), SAMFileHeader.SortOrder.queryname, true, null);
        final List<String[]> iupac = new ArrayList<>();
        for (final String[] pair : pairs.subList(0, 40)) {
            iupac.add(new String[] {pair[0].substring(0, 30) + "RYN" + pair[0].substring(33), pair[1]});
        }
        bwaReadBam(dir.resolve("bwa_iupac.bam"), iupac, SAMFileHeader.SortOrder.queryname, true, null);
        final List<String[]> noSequence = new ArrayList<>(pairs.subList(0, 6));
        noSequence.add(new String[] {"", pairs.get(7)[1]});
        noSequence.addAll(pairs.subList(8, 12));
        bwaReadBam(dir.resolve("bwa_no_sequence.bam"), noSequence, SAMFileHeader.SortOrder.queryname, true, null);

        final SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("s");
        header.addReadGroup(group);
        try (SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(header, true, dir.resolve("bwa_odd.bam"))) {
            for (int i = 0; i != 3; ++i) {
                writer.addAlignment(bwaRead(header, "odd", 77, pairs.get(i)[0], true));
            }
        }
        try (SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(header, true, dir.resolve("bwa_names.bam"))) {
            writer.addAlignment(bwaRead(header, "a", 77, pairs.get(0)[0], true));
            writer.addAlignment(bwaRead(header, "b", 141, pairs.get(0)[1], true));
        }
        try (SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(header, true, dir.resolve("bwa_no_qualities.bam"))) {
            for (int i = 0; i != 10; ++i) {
                writer.addAlignment(bwaRead(header, String.format("q%03d", i), 77, pairs.get(i)[0], false));
                writer.addAlignment(bwaRead(header, String.format("q%03d", i), 141, pairs.get(i)[1], false));
            }
        }

        // PathSeq's microbes are those contigs under names a taxonomy reads, one of them repeated
        // for alternate hits and one no read can reach; its host is the second contig alone.
        final String repeated = first.substring(20000, 24000) + bases.apply(1500);
        final String lonely = bases.apply(1500);
        final Path microbes = dir.resolve("ps_microbe.fasta");
        Files.writeString(microbes, contig.apply("taxid|562|", first) + contig.apply("taxid|1280|", second)
                + contig.apply("taxid|11234|", repeated) + contig.apply("unknown_contig", first.substring(30000, 32000))
                + contig.apply("taxid|300|", lonely), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(microbes, true);
        Files.writeString(dir.resolve("ps_microbe.dict"), "@HD\tVN:1.6\n@SQ\tSN:taxid|562|\tLN:" + first.length()
                + "\n@SQ\tSN:taxid|1280|\tLN:" + second.length() + "\n@SQ\tSN:taxid|11234|\tLN:" + repeated.length()
                + "\n@SQ\tSN:unknown_contig\tLN:2000\n@SQ\tSN:taxid|300|\tLN:" + lonely.length() + "\n",
                StandardCharsets.UTF_8);
        new org.broadinstitute.hellbender.Main().instanceMain(new String[] {"BwaMemIndexImageCreator",
                "-I", microbes.toString(), "-O", dir.resolve("ps_microbe.img").toString()});
        final Path host = dir.resolve("ps_host.fasta");
        Files.writeString(host, contig.apply("host1", second), StandardCharsets.UTF_8);
        new org.broadinstitute.hellbender.Main().instanceMain(new String[] {"BwaMemIndexImageCreator",
                "-I", host.toString(), "-O", dir.resolve("ps_host.img").toString()});
        final String names = String.join("\n", "1\t|\troot\t|\t\t|\tscientific name\t|",
                "2\t|\tBacteria\t|\t\t|\tscientific name\t|", "10239\t|\tViruses\t|\t\t|\tscientific name\t|",
                "562\t|\tEscherichia coli\t|\t\t|\tscientific name\t|", "1280\t|\tStaphylococcus aureus\t|\t\t|\tscientific name\t|",
                "11234\t|\tMeasles morbillivirus\t|\t\t|\tscientific name\t|", "300\t|\tPseudomonas\t|\t\t|\tscientific name\t|") + "\n";
        final String nodes = String.join("\n", "1\t|\t1\t|\tno rank\t|", "2\t|\t1\t|\tsuperkingdom\t|",
                "10239\t|\t1\t|\tsuperkingdom\t|", "562\t|\t2\t|\tspecies\t|", "1280\t|\t2\t|\tspecies\t|",
                "11234\t|\t10239\t|\tspecies\t|", "300\t|\t1\t|\tno rank\t|") + "\n";
        try (final org.apache.commons.compress.archivers.tar.TarArchiveOutputStream tar =
                     new org.apache.commons.compress.archivers.tar.TarArchiveOutputStream(
                             new java.util.zip.GZIPOutputStream(Files.newOutputStream(dir.resolve("ps_taxa.tar.gz"))))) {
            for (final String[] entry : new String[][] {{"names.dmp", names}, {"nodes.dmp", nodes}}) {
                final byte[] bytes = entry[1].getBytes(StandardCharsets.UTF_8);
                final org.apache.commons.compress.archivers.tar.TarArchiveEntry made =
                        new org.apache.commons.compress.archivers.tar.TarArchiveEntry(entry[0]);
                made.setSize(bytes.length);
                made.setModTime(0L);
                tar.putArchiveEntry(made);
                tar.write(bytes);
                tar.closeArchiveEntry();
            }
        }
        Files.writeString(dir.resolve("rp_known.vcf"), "##fileformat=VCFv4.2\n##contig=<ID=bwa1,length=" + first.length()
                + ">\n##contig=<ID=bwa2,length=" + second.length() + ">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                + "bwa1\t5000\t.\t" + first.charAt(4999) + "\t" + (first.charAt(4999) == 'A' ? 'C' : 'A') + "\t.\t.\t.\n",
                StandardCharsets.UTF_8);
        Files.setLastModifiedTime(dir.resolve("rp_known.vcf"), java.nio.file.attribute.FileTime.fromMillis(0L));
        new org.broadinstitute.hellbender.Main().instanceMain(new String[] {"IndexFeatureFile",
                "-I", dir.resolve("rp_known.vcf").toString()});
        try (final htsjdk.samtools.SamReader source = htsjdk.samtools.SamReaderFactory.makeDefault()
                .open(dir.resolve("bwa_reads.bam"))) {
            final SAMFileHeader platform = source.getFileHeader().clone();
            platform.getReadGroups().forEach(line -> line.setPlatform("ILLUMINA"));
            try (SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(platform, true,
                    dir.resolve("rp_reads.bam"))) {
                source.forEach(read -> {
                    read.setHeader(platform);
                    writer.addAlignment(read);
                });
            }
        }
        final SAMFileHeader singles = new SAMFileHeader();
        singles.setSortOrder(SAMFileHeader.SortOrder.queryname);
        final SAMReadGroupRecord singlesGroup = new SAMReadGroupRecord("rg1");
        singlesGroup.setSample("s");
        singlesGroup.setPlatform("ILLUMINA");
        singles.addReadGroup(singlesGroup);
        try (SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(singles, true, dir.resolve("rp_single.bam"))) {
            for (int i = 0; i != 60; ++i) {
                writer.addAlignment(bwaRead(singles, String.format("s%03d", i), 4, pairs.get(i)[0].substring(0, 60), true));
            }
        }
        Files.writeString(dir.resolve("ps_taxa.catalog"), "a\tNOTHING.1\tc\td\te\tf\t562\th\n", StandardCharsets.UTF_8);
        new org.broadinstitute.hellbender.Main().instanceMain(new String[] {"PathSeqBuildReferenceTaxonomy",
                "-R", microbes.toString(), "--tax-dump", dir.resolve("ps_taxa.tar.gz").toString(),
                "-O", dir.resolve("ps_taxa.db").toString(), "--genbank-catalog", dir.resolve("ps_taxa.catalog").toString(),
                "--min-non-virus-contig-length", "0"});
    }

    static void bwaReadBam(final Path bam, final List<String[]> pairs, final SAMFileHeader.SortOrder order,
                           final boolean grouped, final String[] dictionary) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(order);
        if (grouped) {
            final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
            group.setSample("s");
            group.setLibrary("l");
            header.addReadGroup(group);
        }
        if (dictionary != null) {
            for (int at = 0; at < dictionary.length; at += 2) {
                header.addSequence(new SAMSequenceRecord(dictionary[at], Integer.parseInt(dictionary[at + 1])));
            }
        }
        try (SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(header, true, bam)) {
            for (int i = 0; i != pairs.size(); ++i) {
                final String name = String.format("p%04d", i);
                writer.addAlignment(bwaRead(header, name, 77, pairs.get(i)[0], true));
                writer.addAlignment(bwaRead(header, name, 141, pairs.get(i)[1], true));
            }
        }
    }

    static SAMRecord bwaRead(final SAMFileHeader header, final String name, final int flags,
                             final String sequence, final boolean qualities) {
        final SAMRecord read = new SAMRecord(header);
        read.setReadName(name);
        read.setFlags(flags);
        read.setReadString(sequence.isEmpty() ? SAMRecord.NULL_SEQUENCE_STRING : sequence);
        if (qualities && !sequence.isEmpty()) {
            final byte[] scores = new byte[sequence.length()];
            for (int at = 0; at != scores.length; ++at) {
                scores[at] = (byte) (10 + (name.hashCode() + at * 7 & 0x1F));
            }
            read.setBaseQualities(scores);
        }
        if (!header.getReadGroups().isEmpty()) {
            read.setAttribute("RG", "rg1");
        }
        return read;
    }

    /** A sequence of random bases that the same seed always gives. */
    static String randomBases(final java.util.Random random, final int length) {
        final StringBuilder bases = new StringBuilder(length);
        for (int i = 0; i != length; ++i) {
            bases.append("ACGT".charAt(random.nextInt(4)));
        }
        return bases.toString();
    }

    static void kmerReference(final String stem, final String fasta) throws Exception {
        final Path path = DIR.resolve(stem + ".fasta");
        Files.write(path, fasta.getBytes());
        FastaSequenceIndexCreator.create(path, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + path, "O=" + DIR.resolve(stem + ".dict")});
    }

    /** Repeats a few times over and many times over, either case, broken by Ns, and low complexity. */
    static void kmerFixtures() throws Exception {
        final java.util.Random random = new java.util.Random(8675309);
        final String unit = randomBases(random, 300);
        final StringBuilder first = new StringBuilder(">chr1\n");
        for (int copy = 0; copy != 5; ++copy) {
            first.append(copy == 2 ? unit.toLowerCase() : unit).append(randomBases(random, 40));
        }
        first.append("\n>chr2\n").append(unit).append('N').append(unit).append('N').append(unit)
                .append('N').append(unit).append('\n');
        first.append(">chr3\n");
        for (int copy = 0; copy != 40; ++copy) {
            first.append("AC");
        }
        first.append('\n');
        kmerReference("kmers", first.toString());

        // Over 104 records of 9950 bases, so two partitions, and over the maps' first sizes.
        final String big = randomBases(random, 1500);
        final StringBuilder second = new StringBuilder(">big\n");
        for (int copy = 0; copy != 720; ++copy) {
            second.append(big);
        }
        second.append('\n');
        kmerReference("kmers-big", second.toString());

        final String highCopy = ">one\n" + randomBases(random, 60) + "\n" + randomBases(random, 60)
                + "\n>two\r\n" + unit.substring(0, 90) + "\r\n>three\n>four\n"
                + unit.substring(0, 90) + "n" + randomBases(random, 30);
        Files.write(DIR.resolve("high-copy.fasta"), highCopy.getBytes());
        Files.write(DIR.resolve("high-copy-blank.fasta"), ">one\nACGTACGT\n\nACGT\n".getBytes());
        try (final java.util.zip.GZIPOutputStream zipped = new java.util.zip.GZIPOutputStream(
                Files.newOutputStream(DIR.resolve("high-copy.fasta.gz")))) {
            zipped.write(highCopy.getBytes());
        }
    }

    /** Reads over a reference no shard boundary respects: variants beside 5000 and past it. */
    static void callerFixtures() throws Exception {
        final java.util.Random random = new java.util.Random(4242);
        final String bases = randomBases(random, 12000);
        kmerReference("caller", ">chr1\n" + bases + "\n");
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new htsjdk.samtools.SAMSequenceDictionary(
                java.util.List.of(new htsjdk.samtools.SAMSequenceRecord("chr1", 12000))));
        final htsjdk.samtools.SAMReadGroupRecord group = new htsjdk.samtools.SAMReadGroupRecord("rg1");
        group.setSample("caller-sample");
        group.setLibrary("caller-library");
        header.addReadGroup(group);
        final int[] snps = {300, 4985, 5040, 9000};
        final List<SAMRecord> records = new ArrayList<>();
        int serial = 0;
        for (final int[] span : new int[][] {{120, 480}, {4800, 5200}, {8850, 9150}}) {
            for (int start = span[0]; start < span[1]; start += 3) {
                final char[] read = bases.substring(start - 1, start + 99).toCharArray();
                for (final int snp : snps) {
                    if (snp >= start && snp < start + 100 && (start + snp) % 2 == 0) {
                        final int at = snp - start;
                        read[at] = read[at] == 'A' ? 'C' : 'A';
                    }
                }
                final SAMRecord record = new SAMRecord(header);
                record.setReadName("c" + serial++);
                record.setReferenceIndex(0);
                record.setAlignmentStart(start);
                record.setCigarString("100M");
                record.setMappingQuality(60);
                record.setReadString(new String(read));
                final byte[] qualities = new byte[100];
                Arrays.fill(qualities, (byte) 35);
                record.setBaseQualities(qualities);
                record.setAttribute("RG", "rg1");
                record.setReadNegativeStrandFlag(serial % 2 == 0);
                records.add(record);
            }
        }
        for (final String[] one : new String[][] {{"caller.bam", "coordinate"}, {"caller-names.bam", "queryname"}}) {
            final SAMFileHeader own = header.clone();
            own.setSortOrder(SAMFileHeader.SortOrder.valueOf(one[1]));
            final htsjdk.samtools.SAMFileWriterFactory factory = new htsjdk.samtools.SAMFileWriterFactory()
                    .setCreateIndex(one[1].equals("coordinate"));
            try (final htsjdk.samtools.SAMFileWriter writer =
                         factory.makeBAMWriter(own, false, DIR.resolve(one[0]).toFile())) {
                records.forEach(writer::addAlignment);
            }
        }

        // Every twentieth read again under another name and worse qualities, so the pipeline's
        // marking has something to mark, the recalibration something to leave out and, with the
        // second name form, an optical tile to read.
        final List<SAMRecord> withDuplicates = new ArrayList<>();
        int copies = 0;
        for (final SAMRecord record : records) {
            withDuplicates.add(record);
            if ((record.getAlignmentStart() % 20) == 0) {
                final SAMRecord copy = record.deepCopy();
                copy.setReadName("H0:1:1:" + (1000 + copies) + ":1000");
                final byte[] worse = new byte[100];
                Arrays.fill(worse, (byte) 20);
                copy.setBaseQualities(worse);
                withDuplicates.add(copy);
                copies++;
            }
        }
        // The recalibration wants a platform, and `caller-names.bam`'s absence of one is its own row.
        final SAMFileHeader dupHeader = header.clone();
        dupHeader.getReadGroup("rg1").setPlatform("ILLUMINA");
        dupHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        try (final htsjdk.samtools.SAMFileWriter writer = new htsjdk.samtools.SAMFileWriterFactory()
                .setCreateIndex(true).makeBAMWriter(dupHeader, false, DIR.resolve("caller-dups.bam").toFile())) {
            withDuplicates.forEach(writer::addAlignment);
        }
        final SAMFileHeader dupNames = dupHeader.clone();
        dupNames.setSortOrder(SAMFileHeader.SortOrder.queryname);
        try (final htsjdk.samtools.SAMFileWriter writer = new htsjdk.samtools.SAMFileWriterFactory()
                .makeBAMWriter(dupNames, false, DIR.resolve("caller-dups-names.bam").toFile())) {
            withDuplicates.forEach(writer::addAlignment);
        }
        final SAMFileHeader noGroups = header.clone();
        noGroups.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        noGroups.setReadGroups(new ArrayList<>());
        try (final htsjdk.samtools.SAMFileWriter writer = new htsjdk.samtools.SAMFileWriterFactory()
                .makeBAMWriter(noGroups, false, DIR.resolve("caller-no-groups.bam").toFile())) {
            for (final SAMRecord record : records) {
                final SAMRecord copy = record.deepCopy();
                copy.setAttribute("RG", null);
                copy.setHeader(noGroups);
                writer.addAlignment(copy);
            }
        }
        Files.write(DIR.resolve("caller-known.vcf"), ("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=12000>\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                + "chr1\t300\trs1\t" + bases.charAt(299) + "\t" + (bases.charAt(299) == 'A' ? 'C' : 'A') + "\t.\t.\t.\n"
                + "chr1\t5040\trs2\t" + bases.charAt(5039) + "\t" + (bases.charAt(5039) == 'A' ? 'C' : 'A') + "\t.\t.\t.\n")
                .getBytes(java.nio.charset.StandardCharsets.UTF_8));
        // Tribble writes the file's time into the index, which would change the golden on every run.
        Files.setLastModifiedTime(DIR.resolve("caller-known.vcf"), java.nio.file.attribute.FileTime.fromMillis(0L));
        new Main().instanceMain(new String[] {"IndexFeatureFile", "-I", DIR.resolve("caller-known.vcf").toString()});
        Files.write(DIR.resolve("caller-known.bed"), "chr1\t299\t300\n".getBytes(java.nio.charset.StandardCharsets.UTF_8));
    }

    /** READSPIPELINESPARK: the marking, the recalibration, then the caller's shards. */
    static void pipelineRows() throws Exception {
        final String tool = "ReadsPipelineSpark";
        final String ref = DIR.resolve("caller.fasta").toString();
        final String dups = DIR.resolve("caller-dups.bam").toString();
        final String known = DIR.resolve("caller-known.vcf").toString();
        final String[][] runs = {
                {"pipeline-default", "caller.vcf"},
                {"pipeline-gz", "caller.vcf.gz", "--use-jdk-deflater", "true"},
                {"pipeline-bam", "caller.vcf", "--use-jdk-deflater", "true", "--output-bam", DIR.resolve("out/pipeline-bam/out.bam").toString()},
                {"pipeline-bam-remove", "caller.vcf", "--use-jdk-deflater", "true", "--remove-all-duplicates", "--output-bam",
                        DIR.resolve("out/pipeline-bam-remove/out.bam").toString()},
                {"pipeline-gvcf", "caller.g.vcf", "-ERC", "GVCF"},
                {"pipeline-interval", "caller.vcf", "-L", "chr1:8000-9500"},
                {"pipeline-empty-interval", "caller.vcf", "-L", "chr1:2000-3000"},
                {"pipeline-strict", "caller.vcf", "--strict"},
                {"pipeline-small-shards", "caller.vcf", "--read-shard-size", "1000"},
                {"pipeline-tagging", "caller.vcf", "--use-jdk-deflater", "true", "--duplicate-tagging-policy", "All", "--output-bam",
                        DIR.resolve("out/pipeline-tagging/out.bam").toString()},
                {"pipeline-scoring", "caller.vcf", "--duplicate-scoring-strategy", "TOTAL_MAPPED_REFERENCE_LENGTH"},
                {"pipeline-preserve", "caller.vcf", "--preserve-qscores-less-than", "30"},
                {"pipeline-static", "caller.vcf", "--use-jdk-deflater", "true", "--static-quantized-quals", "10", "--static-quantized-quals", "30",
                        "--output-bam", DIR.resolve("out/pipeline-static/out.bam").toString()},
                {"pipeline-original", "caller.vcf", "--use-jdk-deflater", "true", "--emit-original-quals", "--output-bam",
                        DIR.resolve("out/pipeline-original/out.bam").toString()},
                {"pipeline-indels", "caller.vcf", "--indels-context-size", "4", "--mismatches-context-size", "3"},
                {"pipeline-baq-gap", "caller.vcf", "--bqsr-baq-gap-open-penalty", "30.0"},
                {"pipeline-quality", "caller.vcf", "--standard-min-confidence-threshold-for-calling", "2000"},
                {"pipeline-bcf", "caller.bcf"},
                {"pipeline-bamout", "caller.vcf", "-bamout", DIR.resolve("out/pipeline-bamout/out.bam").toString()},
                {"pipeline-dbsnp", "caller.vcf", "--dbsnp", known},
                {"pipeline-comp", "caller.vcf", "--comparison:known", known},
        };
        for (final String[] one : runs) {
            final List<String> rest = new ArrayList<>(Arrays.asList("-I", dups, "-R", ref, "--known-sites", known,
                    "-O", out(one[0], one[1])));
            rest.addAll(Arrays.asList(one).subList(2, one.length));
            run(one[0], one[1], tool, rest.toArray(new String[0]));
        }
        run("pipeline-names", "caller.vcf", tool, "-I", DIR.resolve("caller-names.bam").toString(), "-R", ref,
                "--known-sites", known, "-O", out("pipeline-names", "caller.vcf"));
        run("pipeline-query-sorted", "caller.vcf", tool, "-I", DIR.resolve("caller-dups-names.bam").toString(), "-R", ref,
                "--known-sites", known, "--use-jdk-deflater", "true", "--output-bam",
                DIR.resolve("out/pipeline-query-sorted/out.bam").toString(),
                "-O", out("pipeline-query-sorted", "caller.vcf"));
        run("pipeline-no-groups", "caller.vcf", tool, "-I", DIR.resolve("caller-no-groups.bam").toString(), "-R", ref,
                "--known-sites", known, "-O", out("pipeline-no-groups", "caller.vcf"));
        run("pipeline-bed", "caller.vcf", tool, "-I", dups, "-R", ref,
                "--known-sites", DIR.resolve("caller-known.bed").toString(), "-O", out("pipeline-bed", "caller.vcf"));
        run("recal-no-platform", "recal.table", "BaseRecalibratorSpark", "-I", DIR.resolve("caller.bam").toString(),
                "-R", ref, "--known-sites", known, "-O", out("recal-no-platform", "recal.table"));
        run("pipeline-no-known", null, tool, "-I", dups, "-R", ref, "-O", out("pipeline-no-known", "caller.vcf"));
        run("pipeline-no-reference", null, tool, "-I", dups, "--known-sites", known,
                "-O", out("pipeline-no-reference", "caller.vcf"));
    }

    /** HAPLOTYPECALLERSPARK is HaplotypeCaller's engine, shard by shard. */
    static void callerRows() throws Exception {
        final String tool = "HaplotypeCallerSpark";
        final String ref = DIR.resolve("caller.fasta").toString();
        final String bam = DIR.resolve("caller.bam").toString();
        final String[][] runs = {
                {"caller-default", "caller.vcf"},
                {"caller-gz", "caller.vcf.gz", "--use-jdk-deflater", "true"},
                {"caller-gvcf", "caller.g.vcf", "-ERC", "GVCF"},
                {"caller-bp", "caller.g.vcf", "-ERC", "BP_RESOLUTION"},
                {"caller-small-shards", "caller.vcf", "--read-shard-size", "1000"},
                {"caller-no-padding", "caller.vcf", "--read-shard-padding", "0"},
                {"caller-boundary", "caller.vcf", "--read-shard-size", "5000", "-L", "chr1:4900-5100"},
                {"caller-interval", "caller.vcf", "-L", "chr1:8000-9500"},
                {"caller-empty-interval", "caller.vcf", "-L", "chr1:2000-3000"},
                {"caller-strict", "caller.vcf", "--strict"},
                {"caller-shuffle", "caller.vcf", "--shuffle"},
                {"caller-strict-shuffle", "caller.vcf", "--strict", "--shuffle"},
                {"caller-no-index", "caller.vcf", "--create-output-variant-index", "false"},
                {"caller-downsample", "caller.vcf", "--max-reads-per-alignment-start", "1"},
                {"caller-quality", "caller.vcf", "--standard-min-confidence-threshold-for-calling", "2000"},
                {"caller-bcf", "caller.bcf"},
                {"caller-bamout", "caller.vcf", "-bamout",
                        DIR.resolve("out/caller-bamout/out.bam").toString()},
                {"caller-dbsnp", "caller.vcf", "--dbsnp", DIR.resolve("plain.vcf").toString()},
                {"caller-comp", "caller.vcf", "--comp", DIR.resolve("plain.vcf").toString()},
        };
        for (final String[] one : runs) {
            final List<String> rest = new ArrayList<>(Arrays.asList("-I", bam, "-R", ref, "-O",
                    out(one[0], one[1])));
            rest.addAll(Arrays.asList(one).subList(2, one.length));
            run(one[0], one[1], tool, rest.toArray(new String[0]));
        }
        run("caller-names", "caller.vcf", tool, "-I", DIR.resolve("caller-names.bam").toString(),
                "-R", ref, "-O", out("caller-names", "caller.vcf"));
        run("caller-no-reference", null, tool, "-I", bam, "-O", out("caller-no-reference", "caller.vcf"));
    }

    /** A read of PATHSEQFILTERSPARK's fixtures, with a tag the filter is there to drop. */
    static SAMRecord filterRead(final SAMFileHeader header, final String name, final int flags,
            final String bases, final int quality) {
        final SAMRecord read = new SAMRecord(header);
        read.setReadName(name);
        read.setFlags(flags);
        read.setReadString(bases);
        final byte[] qualities = new byte[bases.length()];
        Arrays.fill(qualities, (byte) quality);
        read.setBaseQualities(qualities);
        read.setAttribute("RG", "rg1");
        read.setAttribute("XX", "dropped");
        return read;
    }

    /** Enough of each thing the filter looks for that each transformer changes what is written. */
    static void filterFixtures() throws Exception {
        final java.util.Random random = new java.util.Random(11);
        final java.util.function.IntFunction<String> sequence = length -> {
            final StringBuilder text = new StringBuilder();
            for (int i = 0; i < length; i++) {
                text.append("ACGT".charAt(random.nextInt(4)));
            }
            return text.toString();
        };
        final SAMFileHeader header = new SAMFileHeader();
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("s");
        header.addReadGroup(group);
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        final List<SAMRecord> reads = new ArrayList<>();
        for (int i = 0; i < 90; i++) {
            final String bases = sequence.apply(70 + i % 11);
            reads.add(filterRead(header, "pair" + i + "/1", 0x1 | 0x4 | 0x8 | 0x40, bases, 30));
            reads.add(filterRead(header, "pair" + i + "/2", 0x1 | 0x4 | 0x8 | 0x80, sequence.apply(75), 30));
        }
        for (int i = 0; i < 30; i++) {
            reads.add(filterRead(header, "lone" + i, 0x4, sequence.apply(80), 30));
        }
        // Each transformer's food: an adapter, repeats, low complexity, quality tails, ambiguity.
        reads.add(filterRead(header, "adapter", 0x4,
                sequence.apply(50) + "AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG", 30));
        reads.add(filterRead(header, "adapter-mismatch", 0x4,
                sequence.apply(50) + "AGATCGGTAGAGCACACGTCTGAACTCCAGTCAC", 30));
        reads.add(filterRead(header, "at-rich", 0x4, "AT".repeat(30) + sequence.apply(30), 30));
        reads.add(filterRead(header, "polya", 0x4, "A".repeat(90), 30));
        reads.add(filterRead(header, "dust-middle", 0x4, sequence.apply(30) + "CAG".repeat(12) + sequence.apply(30), 30));
        reads.add(filterRead(header, "with-n", 0x4, sequence.apply(30) + "N" + sequence.apply(40), 30));
        reads.add(filterRead(header, "many-n", 0x4, sequence.apply(30) + "NNNNN" + sequence.apply(40), 30));
        final SAMRecord tails = filterRead(header, "tails", 0x4, sequence.apply(80), 30);
        final byte[] tailQualities = tails.getBaseQualities();
        Arrays.fill(tailQualities, 0, 6, (byte) 3);
        Arrays.fill(tailQualities, 72, 80, (byte) 2);
        tailQualities[40] = 10;
        tails.setBaseQualities(tailQualities);
        reads.add(tails);
        reads.add(filterRead(header, "all-low", 0x4, sequence.apply(80), 4));
        reads.add(filterRead(header, "short", 0x4, sequence.apply(20), 30));
        // Duplicates: the same bases, the reverse complement, an unpaired one among pairs.
        final String copied = sequence.apply(80);
        reads.add(filterRead(header, "dup-a", 0x4, copied, 30));
        reads.add(filterRead(header, "dup-b", 0x4, copied, 25));
        reads.add(filterRead(header, "dup-rc", 0x4 | 0x10,
                htsjdk.samtools.util.SequenceUtil.reverseComplement(copied), 20));
        final String pairedCopy = sequence.apply(80);
        reads.add(filterRead(header, "dupp/1", 0x1 | 0x4 | 0x8 | 0x40, pairedCopy, 30));
        reads.add(filterRead(header, "dupp/2", 0x1 | 0x4 | 0x8 | 0x80, sequence.apply(75), 30));
        reads.add(filterRead(header, "dupq", 0x4, pairedCopy, 30));
        // What the first pass drops.
        reads.add(filterRead(header, "second", 0x4 | 0x100, sequence.apply(80), 30));
        reads.add(filterRead(header, "supplement", 0x4 | 0x800, sequence.apply(80), 30));
        writeFilterBam(header, reads, "filter.bam");

        // Aligned: the host filter's three answers, a reverse strand turned around, a proper pair.
        final SAMFileHeader aligned = header.clone();
        aligned.addSequence(new htsjdk.samtools.SAMSequenceRecord("chr1", 100000));
        final List<SAMRecord> mapped = new ArrayList<>();
        final String[][] rows = {{"host", "50"}, {"host-weak", "45"}, {"no-nm", null}};
        int at = 100;
        for (final String[] row : rows) {
            for (final int strand : new int[] {0, 0x10}) {
                final SAMRecord read = filterRead(aligned, row[0] + strand, strand, sequence.apply(80), 30);
                read.setReferenceName("chr1");
                read.setAlignmentStart(at += 100);
                read.setMappingQuality(60);
                read.setCigarString("80M");
                if (row[1] != null) {
                    read.setAttribute("NM", 80 - Integer.parseInt(row[1]));
                }
                mapped.add(read);
            }
        }
        final SAMRecord deleted = filterRead(aligned, "deletion", 0, sequence.apply(80), 30);
        deleted.setReferenceName("chr1");
        deleted.setAlignmentStart(at += 100);
        deleted.setMappingQuality(60);
        deleted.setCigarString("40M5D40M");
        deleted.setAttribute("NM", 5);
        mapped.add(deleted);
        final SAMRecord unmapped = filterRead(aligned, "unmapped", 0x4, sequence.apply(80), 30);
        mapped.add(unmapped);
        final SAMRecord proper = filterRead(aligned, "proper/1", 0x1 | 0x2 | 0x20 | 0x40, sequence.apply(80), 30);
        proper.setReferenceName("chr1");
        proper.setAlignmentStart(at += 100);
        proper.setMappingQuality(60);
        proper.setCigarString("80M");
        proper.setMateReferenceName("chr1");
        proper.setMateAlignmentStart(at + 50);
        proper.setAttribute("NM", 70);
        mapped.add(proper);
        writeFilterBam(aligned, mapped, "filter-aligned.bam");

        final SAMFileHeader empty = header.clone();
        writeFilterBam(empty, List.of(), "filter-empty.bam");

        // Host k-mers, from the reference's own builder: the first lone read's and two pairs'.
        final String host = reads.get(180).getReadString() + "\n" + reads.get(0).getReadString() + "\n"
                + reads.get(2).getReadString();
        final StringBuilder fasta = new StringBuilder();
        int contig = 0;
        for (final String line : host.split("\n")) {
            fasta.append(">host").append(contig++).append('\n').append(line).append('\n');
        }
        Files.write(DIR.resolve("host.fasta"), fasta.toString().getBytes(StandardCharsets.UTF_8));
        new Main().instanceMain(new String[] {"PathSeqBuildKmers", "-R", DIR.resolve("host.fasta").toString(),
                "--kmer-size", "31", "--bloom-false-positive-probability", "0",
                "-O", DIR.resolve("host.hss").toString()});
        new Main().instanceMain(new String[] {"PathSeqBuildKmers", "-R", DIR.resolve("host.fasta").toString(),
                "--kmer-size", "31", "--bloom-false-positive-probability", "0.001",
                "-O", DIR.resolve("host.bfi").toString()});
        Files.write(DIR.resolve("host.txt"), new byte[0]);
    }

    static void writeFilterBam(final SAMFileHeader header, final List<SAMRecord> reads, final String name) {
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .makeBAMWriter(header, true, DIR.resolve(name).toFile())) {
            reads.forEach(writer::addAlignment);
        }
    }

    /** One read of PATHSEQSCORESPARK's fixtures: on a contig with an NM, or unmapped. */
    static SAMRecord scoreRead(final SAMFileHeader header, final String name, final int flags, final String contig,
            final int nm, final String... tags) {
        final SAMRecord read = filterRead(header, name, flags | (contig == null ? 0x4 : 0),
                "ACGTTGCAAGGCTTAACCGATCGATGCATCGGATTACGCTAGCTAGGATCCTAGTCCAGTTACGGATCCAATTGGC".substring(0, 70), 30);
        read.setAttribute("XX", null);
        if (contig != null) {
            read.setReferenceName(contig);
            read.setAlignmentStart(100);
            read.setMappingQuality(60);
            read.setCigarString("70M");
            if (nm >= 0) {
                read.setAttribute("NM", nm);
            }
        }
        for (int i = 0; i + 1 < tags.length; i += 2) {
            read.setAttribute(tags[i], tags[i + 1]);
        }
        return read;
    }

    /** A taxonomy built by the reference's own tool, and reads aligned to its contigs. */
    static void scoreFixtures() throws Exception {
        final String[][] contigs = {{"taxid|562|", "800"}, {"taxid|1280|", "900"}, {"taxid|11234|", "700"},
                {"taxid|9606|", "600"}, {"unknown_contig", "500"}, {"taxid|300|", "400"}};
        final StringBuilder fasta = new StringBuilder();
        for (final String[] contig : contigs) {
            fasta.append('>').append(contig[0]).append('\n').append("ACGT".repeat(Integer.parseInt(contig[1]) / 4)).append('\n');
        }
        kmerReference("taxa", fasta.toString());
        final String names = String.join("\n",
                "1\t|\troot\t|\t\t|\tscientific name\t|",
                "2\t|\tBacteria\t|\t\t|\tscientific name\t|",
                "10239\t|\tViruses\t|\t\t|\tscientific name\t|",
                "2759\t|\tEukaryota\t|\t\t|\tscientific name\t|",
                "33208\t|\tMetazoa\t|\t\t|\tscientific name\t|",
                "562\t|\tEscherichia coli\t|\t\t|\tscientific name\t|",
                "1280\t|\tStaphylococcus aureus\t|\t\t|\tscientific name\t|",
                "11234\t|\tMeasles morbillivirus\t|\t\t|\tscientific name\t|",
                "9606\t|\tHomo sapiens\t|\t\t|\tscientific name\t|",
                "300\t|\tPseudomonas\t|\t\t|\tscientific name\t|") + "\n";
        final String nodes = String.join("\n",
                "1\t|\t1\t|\tno rank\t|",
                "2\t|\t1\t|\tsuperkingdom\t|",
                "10239\t|\t1\t|\tsuperkingdom\t|",
                "2759\t|\t1\t|\tsuperkingdom\t|",
                "33208\t|\t2759\t|\tkingdom\t|",
                "562\t|\t2\t|\tspecies\t|",
                "1280\t|\t2\t|\tspecies\t|",
                "11234\t|\t10239\t|\tspecies\t|",
                "9606\t|\t33208\t|\tspecies\t|",
                "300\t|\t1\t|\tno rank\t|") + "\n";
        try (final org.apache.commons.compress.archivers.tar.TarArchiveOutputStream tar =
                     new org.apache.commons.compress.archivers.tar.TarArchiveOutputStream(
                             new java.util.zip.GZIPOutputStream(Files.newOutputStream(DIR.resolve("taxa.tar.gz"))))) {
            for (final String[] entry : new String[][] {{"names.dmp", names}, {"nodes.dmp", nodes}}) {
                final byte[] bytes = entry[1].getBytes(StandardCharsets.UTF_8);
                final org.apache.commons.compress.archivers.tar.TarArchiveEntry header =
                        new org.apache.commons.compress.archivers.tar.TarArchiveEntry(entry[0]);
                header.setSize(bytes.length);
                header.setModTime(0L);
                tar.putArchiveEntry(header);
                tar.write(bytes);
                tar.closeArchiveEntry();
            }
        }
        Files.write(DIR.resolve("taxa.catalog"), "a\tNOTHING.1\tc\td\te\tf\t562\th\n".getBytes(StandardCharsets.UTF_8));
        new Main().instanceMain(new String[] {"PathSeqBuildReferenceTaxonomy", "-R", DIR.resolve("taxa.fasta").toString(),
                "--tax-dump", DIR.resolve("taxa.tar.gz").toString(), "-O", DIR.resolve("taxa.db").toString(),
                "--genbank-catalog", DIR.resolve("taxa.catalog").toString(),
                "--min-non-virus-contig-length", "0"});

        final SAMFileHeader paired = new SAMFileHeader();
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("s");
        paired.addReadGroup(group);
        for (final String[] contig : contigs) {
            paired.addSequence(new htsjdk.samtools.SAMSequenceRecord(contig[0], Integer.parseInt(contig[1])));
        }
        paired.setSortOrder(SAMFileHeader.SortOrder.queryname);
        final String[][] pairs = {
                // Both mates on one taxon, then on two, then disagreeing, then one unmapped.
                {"p0", "taxid|562|", "0", "taxid|562|", "1"},
                {"p1", "taxid|1280|", "0", "taxid|1280|", "2"},
                {"p2", "taxid|562|", "0", "taxid|1280|", "0"},
                {"p3", "taxid|11234|", "1", null, "0"},
                {"p4", "taxid|9606|", "0", "taxid|9606|", "0"},
                {"p5", "unknown_contig", "0", "unknown_contig", "0"},
                {"p6", "taxid|300|", "0", "taxid|300|", "0"},
        };
        final List<SAMRecord> pairReads = new ArrayList<>();
        for (final String[] pair : pairs) {
            pairReads.add(scoreRead(paired, pair[0], 0x1 | 0x40 | (pair[3] == null ? 0x8 : 0), pair[1], Integer.parseInt(pair[2])));
            pairReads.add(scoreRead(paired, pair[0], 0x1 | 0x80 | (pair[3] == null ? 0x4 : 0), pair[3], Integer.parseInt(pair[4])));
        }
        final String xa = "taxid|1280|,+100,70M,1;taxid|11234|,+200,70M,2;";
        pairReads.add(scoreRead(paired, "p7", 0x1 | 0x40, "taxid|562|", 0, "XA", xa));
        pairReads.add(scoreRead(paired, "p7", 0x1 | 0x80, "taxid|562|", 1, "XA", xa));
        writeFilterBam(paired, pairReads, "score-paired.bam");

        final SAMFileHeader unpaired = new SAMFileHeader();
        unpaired.addReadGroup(new SAMReadGroupRecord("rg1"));
        unpaired.addReadGroup(new SAMReadGroupRecord("rg2"));
        unpaired.addSequence(new htsjdk.samtools.SAMSequenceRecord("taxid|1280|", 900));
        unpaired.addSequence(new htsjdk.samtools.SAMSequenceRecord("taxid|562|", 800));
        unpaired.addSequence(new htsjdk.samtools.SAMSequenceRecord("extra_contig", 300));
        unpaired.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        final List<SAMRecord> single = new ArrayList<>();
        single.add(scoreRead(unpaired, "u0", 0, "taxid|562|", 0));
        single.add(scoreRead(unpaired, "u1", 0, "taxid|1280|", 0, "SA", "taxid|562|,50,+,70M,60,3;"));
        single.add(scoreRead(unpaired, "u2", 0, "taxid|1280|", 30));
        single.add(scoreRead(unpaired, "u3", 0, "extra_contig", 0));
        single.add(scoreRead(unpaired, "u4", 0, null, 0));
        single.add(scoreRead(unpaired, "u5", 0x100, "taxid|562|", 0));
        single.add(scoreRead(unpaired, "u6", 0x800, "taxid|562|", 0));
        writeFilterBam(unpaired, single, "score-unpaired.bam");

        writeFilterBam(unpaired, List.of(scoreRead(unpaired, "n", 0, "taxid|562|", -1)), "score-no-nm.bam");
        writeFilterBam(unpaired, List.of(scoreRead(unpaired, "bad", 0, "taxid|562|", 0, "XA", "a,1")), "score-bad-xa.bam");
        writeFilterBam(unpaired, List.of(scoreRead(unpaired, "bad", 0, "taxid|562|", 0, "XA", "a,1,70M,x;")), "score-bad-nm.bam");
        writeFilterBam(unpaired, List.of(scoreRead(unpaired, "bad", 0, "taxid|562|", 0, "XA", "a,1,7Q,1;")), "score-bad-cigar.bam");
        final SAMFileHeader lonely = paired.clone();
        writeFilterBam(lonely, List.of(scoreRead(lonely, "only", 0x1 | 0x40, "taxid|562|", 0)), "score-lonely.bam");
        final SAMFileHeader coordinate = paired.clone();
        coordinate.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        writeFilterBam(coordinate, List.of(), "score-not-queryname.bam");
    }

    /** PATHSEQSCORESPARK: the identity arguments, the normalizations and the outputs. */
    /** EXTRACTSVEVIDENCESPARK's text, the metadata it takes on the way and its high-depth regions. */
    static void evidenceRows() throws Exception {
        final String tool = "ExtractSVEvidenceSpark";
        final String ignore = DIR.resolve("sv-ignore.txt").toString();
        final String[][] runs = {
                {"evidence", "sv-evidence.bam"},
                {"evidence-other-sizes", "sv-evidence2.bam"},
                {"evidence-unpaired", "sv-evidence-lone.bam"},
                {"evidence-ignore-alt", "sv-evidence.bam", "--cross-contigs-to-ignore", ignore},
                {"evidence-ignore-unknown", "sv-evidence.bam", "--cross-contigs-to-ignore",
                        DIR.resolve("sv-ignore-bad.txt").toString()},
                {"evidence-ignore-missing", "sv-evidence.bam", "--cross-contigs-to-ignore",
                        DIR.resolve("nowhere.txt").toString()},
                {"evidence-strict-filter", "sv-evidence.bam", "--min-evidence-mapq", "40",
                        "--min-evidence-match-length", "80"},
                {"evidence-no-overhang", "sv-evidence.bam", "--allowed-short-fragment-overhang", "0",
                        "--max-tracked-fragment-length", "500"},
                {"evidence-low-depth", "sv-evidence2.bam", "--high-depth-coverage-factor", "1",
                        "--high-depth-coverage-peak-factor", "2"},
                {"evidence-interval", "sv-evidence.bam", "-L", "chr2"},
                {"evidence-interval-pile", "sv-evidence.bam", "-L", "chr1:7000-9000"},
                {"evidence-queryname", "marking-queryname.bam"},
                {"evidence-pairs", "pairs.bam"},
                {"evidence-no-reads", "sv-evidence.bam", "-L", "alt:4000-4100"},
        };
        for (final String[] one : runs) {
            final String name = one[0];
            final List<String> rest = new ArrayList<>(Arrays.asList(
                    "-I", DIR.resolve(one[1]).toString(), "-O", out(name, "evidence"),
                    "--aligner-index-image", "image.img", "--kmers-to-ignore", "kmers.txt",
                    "--read-metadata", out(name, "metadata.txt"),
                    "--high-coverage-intervals", out(name, "high.bed")));
            rest.addAll(Arrays.asList(one).subList(2, one.length));
            run(name, "evidence", tool, rest.toArray(new String[0]));
        }
    }

    /** Each signature, each way a contig or an argument is refused, and what the reads are cut to. */
    static void discoveryRows() throws Exception {
        final String tool = "DiscoverVariantsFromContigAlignmentsSAMSpark";
        final String cnv = DIR.resolve("sv_cnv.vcf").toString();
        final String truth = DIR.resolve("sv_truth.vcf").toString();
        final String[][] runs = {
                {"discover", "sv_contigs.bam"},
                {"discover-small", "sv_contigs_small.bam"},
                {"discover-cnv", "sv_contigs.bam", "--cnv-calls", cnv},
                {"discover-truth", "sv_contigs.bam", "--truth-vcf", truth},
                {"discover-truth-near", "sv_contigs.bam", "--truth-vcf", truth, "--truth-interval-padding", "1"},
                {"discover-interval", "sv_contigs.bam", "-L", "chr1:1-30000"},
                {"discover-excluded", "sv_contigs.bam", "-XL", "chr2:13000-20000"},
                {"discover-cut-contig", "sv_contigs.bam", "-L", "chr2"},
                {"discover-no-default-filters", "sv_contigs.bam", "--disable-tool-default-read-filters", "true"},
                {"discover-no-supplementary", "sv_contigs.bam", "--read-filter", "NotSupplementaryAlignmentReadFilter"},
                {"discover-ignored-arguments", "sv_contigs.bam", "--min-align-length", "500", "--min-mq", "-1",
                        "--debug-mode", "true"},
                {"discover-adjacent-reverse", "sv_contigs_adjacent_reverse.bam"},
                {"discover-no-primary", "sv_contigs_no_primary.bam"},
                {"discover-hard-primary", "sv_contigs_hard_primary.bam"},
                {"discover-no-sample", "sv_contigs_no_sample.bam"},
                {"discover-two-samples", "sv_contigs_two_samples.bam"},
                {"discover-bad-length", "sv_contigs.bam", "--min-align-length", "0"},
                {"discover-bad-uncertainty", "sv_contigs.bam", "--assembly-imprecise-evidence-overlap-uncertainty", "0"},
                {"discover-bad-threshold", "sv_contigs.bam", "--imprecise-variant-evidence-threshold", "0"},
                {"discover-bad-padding", "sv_contigs.bam", "--truth-interval-padding", "0"},
                {"discover-bad-deletion", "sv_contigs.bam", "--max-callable-imprecise-deletion-size", "0"},
                {"discover-cnv-missing", "sv_contigs.bam", "--cnv-calls", DIR.resolve("nowhere.vcf").toString()},
                {"discover-cnv-other-sample", "sv_contigs.bam", "--cnv-calls", DIR.resolve("sv_cnv_other_sample.vcf").toString()},
                {"discover-cnv-no-dictionary", "sv_contigs.bam", "--cnv-calls", DIR.resolve("sv_cnv_no_dictionary.vcf").toString()},
                {"discover-cnv-other-dictionary", "sv_contigs.bam", "--cnv-calls", DIR.resolve("sv_cnv_other_dictionary.vcf").toString()},
                {"discover-cnv-two-samples", "sv_contigs.bam", "--cnv-calls", DIR.resolve("sv_cnv_two_samples.vcf").toString()},
                {"discover-truth-missing", "sv_contigs.bam", "--truth-vcf", DIR.resolve("nowhere.vcf").toString()},
                {"discover-truth-near-start", "sv_contigs.bam", "--truth-vcf", DIR.resolve("sv_truth_near_start.vcf").toString()},
                {"discover-truth-unknown-type", "sv_contigs.bam", "--truth-vcf", DIR.resolve("sv_truth_unknown_type.vcf").toString()},
                {"discover-truth-unknown-contig", "sv_contigs.bam", "--truth-vcf", DIR.resolve("sv_truth_unknown_contig.vcf").toString()},
        };
        for (final String[] one : runs) {
            final String name = one[0];
            final List<String> rest = new ArrayList<>(Arrays.asList(
                    "-I", DIR.resolve(one[1]).toString(), "-R", DIR.resolve("sv_contigs.fasta").toString(),
                    "-O", out(name, "discovered"), "--add-output-vcf-command-line", "false"));
            rest.addAll(Arrays.asList(one).subList(2, one.length));
            run(name, "discovered", tool, rest.toArray(new String[0]));
        }
        // The prefix's other two shapes: a directory, and a file it will not overwrite.
        run("discover-directory", "discovered", tool, "-I", DIR.resolve("sv_contigs.bam").toString(),
                "-R", DIR.resolve("sv_contigs.fasta").toString(), "--add-output-vcf-command-line", "false",
                "-O", java.nio.file.Paths.get(out("discover-directory", "discovered")).getParent().toString());
        run("discover-existing-file", "discovered", tool, "-I", DIR.resolve("sv_contigs.bam").toString(),
                "-R", DIR.resolve("sv_contigs.fasta").toString(), "-O", DIR.resolve("sv_contigs.dict").toString());
    }

    /**
     * SvDiscoverFromLocalAssemblyContigAlignmentsSpark: five VCFs under the prefix, the sample's
     * name in each, whatever the calls, and with `--write-sam` a BAM of the contigs each class
     * could not use.
     */
    static void assemblyRows() throws Exception {
        final String tool = "SvDiscoverFromLocalAssemblyContigAlignmentsSpark";
        final String cnv = DIR.resolve("sv_cnv.vcf").toString();
        final String truth = DIR.resolve("sv_truth.vcf").toString();
        final String altTigs = DIR.resolve("sv_alt_tigs.txt").toString();
        final String[][] runs = {
                {"assembly", "sv_assembly.bam"},
                {"assembly-contigs", "sv_contigs.bam"},
                {"assembly-small", "sv_contigs_small.bam"},
                {"assembly-write-sam", "sv_assembly.bam", "--write-sam", "--use-jdk-deflater", "true"},
                {"assembly-alt-tigs", "sv_assembly.bam", "--non-canonical-contig-names-file", altTigs},
                {"assembly-alt-tigs-missing", "sv_assembly.bam", "--non-canonical-contig-names-file",
                        DIR.resolve("nowhere.txt").toString()},
                {"assembly-truth", "sv_assembly.bam", "--truth-vcf", truth},
                {"assembly-truth-near", "sv_assembly.bam", "--truth-vcf", truth, "--truth-interval-padding", "1"},
                {"assembly-truth-missing", "sv_assembly.bam", "--truth-vcf", DIR.resolve("nowhere.vcf").toString()},
                {"assembly-strict", "sv_assembly.bam", "--min-mq", "70", "--min-align-length", "90"},
                {"assembly-loose", "sv_assembly.bam", "--min-mq", "10", "--min-align-length", "20"},
                {"assembly-cnv", "sv_assembly_no_breakend.bam", "--cnv-calls", cnv},
                {"assembly-cnv-breakend", "sv_assembly.bam", "--cnv-calls", cnv},
                {"assembly-cnv-missing", "sv_assembly_no_breakend.bam", "--cnv-calls", DIR.resolve("nowhere.vcf").toString()},
                {"assembly-no-breakend", "sv_assembly_no_breakend.bam"},
                {"assembly-too-many", "sv_assembly_too_many.bam"},
                {"assembly-overlap-contained", "sv_assembly_overlap_contained.bam"},
                {"assembly-interval", "sv_assembly.bam", "-L", "chr1:1-30000"},
                {"assembly-excluded", "sv_assembly.bam", "-XL", "chr2:13000-20000"},
                {"assembly-no-default-filters", "sv_assembly.bam", "--disable-tool-default-read-filters", "true"},
                {"assembly-no-supplementary", "sv_assembly.bam", "--read-filter", "NotSupplementaryAlignmentReadFilter"},
                {"assembly-no-sample", "sv_contigs_no_sample.bam"},
                {"assembly-two-samples", "sv_contigs_two_samples.bam"},
                {"assembly-bad-length", "sv_assembly.bam", "--min-align-length", "0"},
                {"assembly-bad-uncertainty", "sv_assembly.bam", "--assembly-imprecise-evidence-overlap-uncertainty", "0"},
                {"assembly-bad-threshold", "sv_assembly.bam", "--imprecise-variant-evidence-threshold", "0"},
                {"assembly-bad-padding", "sv_assembly.bam", "--truth-interval-padding", "0"},
                {"assembly-bad-deletion", "sv_assembly.bam", "--max-callable-imprecise-deletion-size", "0"},
        };
        for (final String[] one : runs) {
            final String name = one[0];
            final List<String> rest = new ArrayList<>(Arrays.asList(
                    "-I", DIR.resolve(one[1]).toString(), "-R", DIR.resolve("sv_contigs.fasta").toString(),
                    "-O", out(name, "assembled"), "--add-output-vcf-command-line", "false"));
            rest.addAll(Arrays.asList(one).subList(2, one.length));
            run(name, "assembled", tool, rest.toArray(new String[0]));
        }
        // The prefix's other two shapes: a directory, and a file it will not overwrite.
        run("assembly-directory", "assembled", tool, "-I", DIR.resolve("sv_assembly.bam").toString(),
                "-R", DIR.resolve("sv_contigs.fasta").toString(), "--add-output-vcf-command-line", "false",
                "-O", java.nio.file.Paths.get(out("assembly-directory", "assembled")).getParent().toString());
        run("assembly-existing-file", "assembled", tool, "-I", DIR.resolve("sv_assembly.bam").toString(),
                "-R", DIR.resolve("sv_contigs.fasta").toString(), "-O", DIR.resolve("sv_contigs.dict").toString());
        // PrintReadsSpark takes GATKSparkTool's wellformed filter: the read with no qualities goes.
        run("print-reads-no-qualities", "no-qualities.bam", "PrintReadsSpark", "-I",
                DIR.resolve("sv_quality_missing.bam").toString(), "--use-jdk-deflater", "true",
                "-O", out("print-reads-no-qualities", "no-qualities.bam"));
    }

    /**
     * CpxVariantReInterpreterSpark: the assembly step's last two files, from a VCF of complex calls
     * and the contigs that made them, under the VCF's dictionary rather than the reads'.
     */
    static void reinterpreterRows() throws Exception {
        final String tool = "CpxVariantReInterpreterSpark";
        final String altTigs = DIR.resolve("sv_alt_tigs.txt").toString();
        final String[][] runs = {
                {"cpx", "sv_assembly.bam", "sv_cpx.vcf"},
                {"cpx-other-reads", "sv_contigs.bam", "sv_cpx.vcf"},
                {"cpx-no-reads-of-its-contigs", "sv_contigs_small.bam", "sv_cpx.vcf"},
                {"cpx-alt-tigs", "sv_assembly.bam", "sv_cpx.vcf", "--non-canonical-contig-names-file", altTigs},
                {"cpx-alt-tigs-missing", "sv_assembly.bam", "sv_cpx.vcf", "--non-canonical-contig-names-file",
                        DIR.resolve("nowhere.txt").toString()},
                {"cpx-empty", "sv_assembly.bam", "sv_cpx_empty.vcf"},
                {"cpx-sparse", "sv_assembly.bam", "sv_cpx_sparse.vcf"},
                {"cpx-no-dictionary", "sv_assembly.bam", "sv_cpx_no_dictionary.vcf"},
                {"cpx-no-dictionary-one", "sv_assembly.bam", "sv_cpx_no_dictionary_one.vcf"},
                {"cpx-bad-segment", "sv_assembly.bam", "sv_cpx_bad_segment.vcf"},
                {"cpx-whole-contig", "sv_assembly.bam", "sv_cpx_whole_contig.vcf"},
                {"cpx-bad-length", "sv_assembly.bam", "sv_cpx_bad_length.vcf"},
                {"cpx-missing", "sv_assembly.bam", "nowhere.vcf"},
                {"cpx-interval-no-index", "sv_assembly.bam", "sv_cpx.vcf", "-L", "chr1:1-30000"},
                {"cpx-no-default-filters", "sv_assembly.bam", "sv_cpx.vcf", "--disable-tool-default-read-filters", "true"},
                {"cpx-no-supplementary", "sv_assembly.bam", "sv_cpx.vcf", "--read-filter", "NotSupplementaryAlignmentReadFilter"},
                {"cpx-unvalidated-length", "sv_assembly.bam", "sv_cpx.vcf", "--min-align-length", "0"},
                {"cpx-unvalidated-uncertainty", "sv_assembly.bam", "sv_cpx.vcf", "--assembly-imprecise-evidence-overlap-uncertainty", "0"},
                {"cpx-ignored-truth", "sv_assembly.bam", "sv_cpx.vcf", "--truth-vcf", DIR.resolve("nowhere.vcf").toString()},
                {"cpx-ignored-cnv", "sv_assembly.bam", "sv_cpx.vcf", "--cnv-calls", DIR.resolve("nowhere.vcf").toString()},
        };
        for (final String[] one : runs) {
            final String name = one[0];
            final List<String> rest = new ArrayList<>(Arrays.asList(
                    "-I", DIR.resolve(one[1]).toString(), "-R", DIR.resolve("sv_contigs.fasta").toString(),
                    "--cpx-vcf", DIR.resolve(one[2]).toString(),
                    "--prefix-out-vcf", out(name, "reinterpreted"), "--add-output-vcf-command-line", "false"));
            rest.addAll(Arrays.asList(one).subList(3, one.length));
            run(name, "reinterpreted", tool, rest.toArray(new String[0]));
        }
        // A prefix whose directory is not there, and the two required arguments.
        run("cpx-prefix-nowhere", "reinterpreted", tool, "-I", DIR.resolve("sv_assembly.bam").toString(),
                "-R", DIR.resolve("sv_contigs.fasta").toString(), "--cpx-vcf", DIR.resolve("sv_cpx.vcf").toString(),
                "--prefix-out-vcf", DIR.resolve("nowhere").resolve("reinterpreted").toString(),
                "--add-output-vcf-command-line", "false");
        run("cpx-no-vcf", "reinterpreted", tool, "-I", DIR.resolve("sv_assembly.bam").toString(),
                "-R", DIR.resolve("sv_contigs.fasta").toString(), "--prefix-out-vcf", out("cpx-no-vcf", "reinterpreted"));
        run("cpx-no-prefix", "reinterpreted", tool, "-I", DIR.resolve("sv_assembly.bam").toString(),
                "-R", DIR.resolve("sv_contigs.fasta").toString(), "--cpx-vcf", DIR.resolve("sv_cpx.vcf").toString());
    }

    /**
     * What `DiscoverVariantsFromContigAlignmentsSAMSpark` reads: a reference, and assembled contigs
     * whose alignments carry each signature its interpreter tells apart.
     */
    static void svContigReads(final Path dir) throws Exception {
        final int[] lengths = {30000, 20000};
        final String[] names = {"chr1", "chr2"};
        final java.util.Random random = new java.util.Random(11L);
        final StringBuilder fasta = new StringBuilder();
        final StringBuilder index = new StringBuilder();
        final StringBuilder dictionary = new StringBuilder("@HD\tVN:1.6\n");
        final java.util.List<SAMSequenceRecord> sequences = new java.util.ArrayList<>();
        for (int contig = 0; contig < names.length; contig++) {
            final String line = ">" + names[contig] + "\n";
            fasta.append(line);
            index.append(names[contig]).append('\t').append(lengths[contig]).append('\t')
                    .append(fasta.length()).append("\t60\t61\n");
            for (int at = 0; at < lengths[contig]; at++) {
                fasta.append("ACGT".charAt(random.nextInt(4)));
                if (at % 60 == 59 || at == lengths[contig] - 1) {
                    fasta.append('\n');
                }
            }
            dictionary.append("@SQ\tSN:").append(names[contig]).append("\tLN:").append(lengths[contig]).append('\n');
            sequences.add(new SAMSequenceRecord(names[contig], lengths[contig]));
        }
        Files.writeString(dir.resolve("sv_contigs.fasta"), fasta.toString(), StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_contigs.fasta.fai"), index.toString(), StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_contigs.dict"), dictionary.toString(), StandardCharsets.UTF_8);

        // Name, then per alignment: flags, contig, start, cigar, mapping quality.
        final Object[][] contigs = {
            {"del", 0, "chr1", 1000, "100M100S", 60, 2048, "chr1", 2000, "100S100M", 60},
            {"del2", 0, "chr1", 1000, "100M100S", 60, 2048, "chr1", 2000, "100H100M", 60},
            {"delRev", 16, "chr1", 3000, "100M100S", 60, 2064, "chr1", 4000, "100S100M", 60},
            {"delHom", 0, "chr1", 5000, "105M95S", 60, 2048, "chr1", 6000, "100S100M", 60},
            {"delHomTie", 0, "chr1", 5000, "107M93S", 60, 2048, "chr1", 6002, "100S100M", 60},
            {"replace", 0, "chr1", 7000, "100M160S", 60, 2048, "chr1", 8000, "160S100M", 60},
            {"replaceShort", 0, "chr1", 9000, "100M110S", 60, 2048, "chr1", 9500, "110S100M", 60},
            {"ins", 0, "chr1", 10000, "100M160S", 60, 2048, "chr1", 10100, "160S100M", 60},
            {"insRev", 16, "chr1", 10500, "100M160S", 60, 2064, "chr1", 10600, "160S100M", 60},
            {"dup", 0, "chr1", 11000, "100M100S", 60, 2048, "chr1", 11080, "100S100M", 60},
            {"dupRev", 16, "chr1", 11500, "100M100S", 60, 2064, "chr1", 11580, "100S100M", 60},
            {"dupIns", 0, "chr1", 12000, "100M130S", 60, 2048, "chr1", 12080, "130S100M", 60},
            {"dupGap", 0, "chr1", 12500, "50M2I48M100S", 60, 2048, "chr1", 12580, "100S60M1D40M", 60},
            {"contract", 0, "chr1", 13000, "100M80S", 60, 2048, "chr1", 13100, "80S100M", 60},
            {"cpxUp", 0, "chr1", 14000, "100M90S", 60, 2048, "chr1", 14070, "90S100M", 60},
            {"cpxUpRev", 16, "chr1", 14500, "100M90S", 60, 2064, "chr1", 14570, "90S100M", 60},
            {"cpxDown", 0, "chr1", 15000, "100M70S", 60, 2048, "chr1", 15090, "70S100M", 60},
            {"inv55", 0, "chr1", 16000, "100M100S", 60, 2064, "chr1", 17000, "100M100S", 60},
            {"inv33", 16, "chr1", 18000, "100S100M", 60, 2048, "chr1", 19000, "100S100M", 60},
            {"inv55Ins", 0, "chr1", 20000, "100M120S", 60, 2064, "chr1", 21000, "100M120S", 60},
            {"inv33Hom", 16, "chr1", 22000, "95S105M", 60, 2048, "chr1", 23000, "100S100M", 60},
            {"mapped", 0, "chr1", 24000, "100M200S", 60, 2048, "chr2", 500, "100S100M100S", 20, 2048, "chr1", 25000, "200S100M", 60},
            {"threes", 0, "chr1", 26000, "100M200S", 60, 2048, "chr1", 26500, "100S100M100S", 60, 2048, "chr1", 27000, "200S100M", 60},
            {"gapped", 0, "chr2", 1000, "100M60D100M", 60},
            {"gappedIns", 0, "chr2", 2000, "100M60I100M", 60},
            {"gappedSmall", 0, "chr2", 3000, "100M10D100M", 60},
            {"gappedClip", 16, "chr2", 4000, "20S100M60D100M30S", 60},
            {"lowFirst", 0, "chr2", 5000, "100M100S", 10, 2048, "chr2", 6000, "100S100M", 60},
            {"lowEnd", 0, "chr2", 7000, "100M100S", 60, 2048, "chr2", 8000, "100S100M", 10},
            {"tooShort", 0, "chr2", 9000, "30M100S", 60, 2048, "chr2", 9500, "30S100M", 60},
            {"nextShort", 0, "chr2", 9800, "100M30S", 60, 2048, "chr2", 9950, "100S30M", 60},
            {"inter", 0, "chr1", 28000, "100M100S", 60, 2048, "chr2", 10000, "100S100M", 60},
            {"swap", 0, "chr2", 11000, "100M100S", 60, 2048, "chr2", 10500, "100S100M", 60},
            {"contained", 0, "chr2", 12000, "100M100S", 60, 2048, "chr2", 12020, "100S60M40S", 60},
            {"stitched", 0, "chr2", 13000, "100M100S", 60, 2048, "chr2", 13100, "100S100M", 60},
            {"single", 0, "chr2", 14000, "100M", 60},
            {"withSecondary", 0, "chr2", 15000, "100M100S", 60, 2048, "chr2", 16000, "100S100M", 60, 256, "chr2", 17000, "100M100S", 60},
        };
        svContigBam(dir.resolve("sv_contigs.bam"), sequences, "s", contigs, 12L);
        svContigBam(dir.resolve("sv_contigs_small.bam"), sequences, "s", new Object[][] {contigs[0], contigs[17]}, 13L);
        svContigBam(dir.resolve("sv_contigs_adjacent_reverse.bam"), sequences, "s", new Object[][] {
            contigs[0],
            {"adjacentReverse", 16, "chr2", 18000, "100M100S", 60, 2064, "chr2", 18100, "100S100M", 60}}, 14L);
        svContigBam(dir.resolve("sv_contigs_no_primary.bam"), sequences, "s", new Object[][] {
            {"orphan", 2048, "chr1", 1000, "100M100S", 60}}, 15L);
        svContigBam(dir.resolve("sv_contigs_hard_primary.bam"), sequences, "s", new Object[][] {
            {"hard", 0, "chr1", 1000, "100M100H", 60}}, 16L);
        svContigBam(dir.resolve("sv_contigs_no_sample.bam"), sequences, null, new Object[][] {contigs[0]}, 17L);
        svContigBam(dir.resolve("sv_contigs_two_samples.bam"), sequences, "s,t", new Object[][] {contigs[0]}, 18L);
        svAssemblyContigs(dir, sequences);
        svComplexCalls(dir);

        final StringBuilder vcfHead = new StringBuilder("##fileformat=VCFv4.2\n"
                + "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number\">\n"
                + "##FORMAT=<ID=CNQ,Number=1,Type=Float,Description=\"Copy number quality\">\n"
                + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop\">\n"
                + "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">\n");
        for (int contig = 0; contig < names.length; contig++) {
            vcfHead.append("##contig=<ID=").append(names[contig]).append(",length=").append(lengths[contig]).append(">\n");
        }
        final String columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
        final String calls = "chr1\t900\tcnvA\tN\t<DEL>\t.\t.\tEND=2500;SVTYPE=DEL\tGT:CN:CNQ\t0/1:1:40\n"
                + "chr1\t1500\tcnvB\tN\t<DEL>\t.\t.\tEND=4000;SVTYPE=DEL\tGT:CN:CNQ\t0/1:1:30.5\n"
                + "chr2\t1\tcnvC\tN\t<DUP>\t.\t.\tEND=9000;SVTYPE=DUP\tGT:CN:CNQ\t0/1:3:20\n";
        Files.writeString(dir.resolve("sv_cnv.vcf"), vcfHead + columns + "s\n" + calls, StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_cnv_other_sample.vcf"), vcfHead + columns + "t\n" + calls, StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_cnv_no_dictionary.vcf"),
                vcfHead.toString().replaceAll("##contig=[^\n]*\n", "") + columns + "s\n" + calls, StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_cnv_other_dictionary.vcf"),
                vcfHead.toString().replace("length=20000", "length=20001") + columns + "s\n" + calls, StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_cnv_two_samples.vcf"),
                vcfHead + columns + "s\tt\n" + calls.replace("\n", "\t0/1:1:40\n"), StandardCharsets.UTF_8);
        final String truthColumns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        Files.writeString(dir.resolve("sv_truth.vcf"),
                vcfHead + truthColumns
                        + "chr1\t1100\ttruthA\tN\t<DEL>\t.\t.\tEND=2000;SVTYPE=DEL\n"
                        + "chr2\t100\ttruthB\tN\t<INS>\t.\t.\tEND=100;SVTYPE=INS\n"
                        + "chr2\t200\ttruthC\tN\t.\t.\t.\tEND=200\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_truth_near_start.vcf"),
                vcfHead + truthColumns + "chr2\t1\ttruthB\tN\t<INS>\t.\t.\tEND=1;SVTYPE=INS\n", StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_truth_unknown_type.vcf"),
                vcfHead + truthColumns + "chr2\t100\ttruthB\tN\t<INS>\t.\t.\tEND=100;SVTYPE=ODD\n", StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_truth_unknown_contig.vcf"),
                vcfHead + truthColumns + "chr9\t100\ttruthB\tN\t<INS>\t.\t.\tEND=100;SVTYPE=INS\n", StandardCharsets.UTF_8);
    }

    static void svAssemblyContigs(final Path dir, final java.util.List<SAMSequenceRecord> sequences) throws Exception {
        final Object[][] contigs = {
            {"cpxMapped", 0, "chr1", 1000, "100M200S", 60, 2048, "chr2", 500, "100S100M100S", 60, 2048, "chr1", 1500, "200S100M", 60},
            {"cpxInverted", 0, "chr1", 2000, "100M200S", 60, 2064, "chr1", 2300, "100S100M100S", 60, 2048, "chr1", 2600, "200S100M", 60},
            {"cpxDuplicated", 0, "chr1", 3000, "100M200S", 60, 2048, "chr1", 3200, "100S100M100S", 60, 2048, "chr1", 3100, "200S300M", 60},
            {"cpxReverse", 16, "chr1", 5000, "200S100M", 60, 2064, "chr1", 4700, "100S100M100S", 60, 2064, "chr1", 4400, "100M200S", 60},
            {"cpxReverseFlip", 16, "chr1", 5900, "200S100M", 60, 2048, "chr1", 5600, "100S100M100S", 60, 2064, "chr1", 5300, "100M200S", 60},
            {"cpxGapped", 0, "chr1", 6000, "100M320S", 60, 2048, "chr1", 6400, "120S100M200S", 60, 2048, "chr1", 6900, "320S100M", 60},
            {"cpxOverlap", 0, "chr1", 8000, "110M200S", 60, 2048, "chr1", 8300, "100S110M100S", 60, 2048, "chr1", 8700, "200S110M", 60},
            {"cpxFour", 0, "chr1", 9500, "100M400S", 60, 2048, "chr1", 9800, "100S110M290S", 60, 2048, "chr1", 10200, "200S110M190S", 60, 2048, "chr1", 10800, "400S100M", 60},
            {"cpxFourSwitch", 0, "chr1", 11000, "100M400S", 60, 2048, "chr1", 11200, "100S110M290S", 60, 2064, "chr1", 11500, "200S110M190S", 60, 2048, "chr1", 11800, "400S100M", 60},
            {"cpxFourOther", 0, "chr1", 12000, "100M400S", 60, 2048, "chr2", 1100, "100S110M290S", 60, 2048, "chr1", 12300, "200S110M190S", 60, 2048, "chr1", 12800, "400S100M", 60},
            {"cpxClipDeletion", 0, "chr1", 13000, "100M3D10M200S", 60, 2048, "chr1", 13400, "100S10M3I100M100S", 60, 2048, "chr1", 13800, "210S103M", 60},
            {"cpxShared", 0, "chr1", 14000, "100M200S", 60, 2048, "chr2", 3000, "100S100M100S", 60, 2048, "chr1", 14099, "200S100M", 60},
            {"cpxSharedFlip", 0, "chr1", 14500, "100M250S", 60, 2064, "chr2", 3400, "130S100M120S", 60, 2048, "chr1", 14599, "250S100M", 60},
            {"cpxInvertedOne", 0, "chr1", 15000, "100M300S", 60, 2064, "chr1", 15099, "100S200M100S", 60, 2048, "chr1", 15298, "300S100M", 60},
            {"cpxInvertedOneWide", 0, "chr1", 15600, "100M500S", 60, 2064, "chr1", 15699, "200S200M200S", 60, 2048, "chr1", 15898, "500S100M", 60},
            {"cpxInvertedSmall", 0, "chr1", 16200, "100M240S", 60, 2064, "chr1", 16299, "100S20M220S", 60, 2048, "chr1", 16318, "240S100M", 60},
            {"cpxAsIsWide", 0, "chr1", 16600, "100M500S", 60, 2048, "chr1", 16699, "200S200M200S", 60, 2048, "chr1", 16898, "500S100M", 60},
            {"cpxAsIsFront", 0, "chr1", 17200, "100M400S", 60, 2048, "chr1", 17299, "200S200M100S", 60, 2048, "chr1", 17498, "400S100M", 60},
            {"cpxAsIsRear", 0, "chr1", 17800, "100M400S", 60, 2048, "chr1", 17899, "100S200M200S", 60, 2048, "chr1", 18098, "400S100M", 60},
            {"cpxAllAsIs", 0, "chr1", 18400, "100M700S", 60, 2048, "chr1", 18499, "200S150M450S", 60, 2048, "chr1", 18648, "350S150M300S", 60, 2048, "chr1", 18797, "600S100M100S", 60},
            {"cpxSegmentGone", 0, "chr1", 19200, "100M200S", 60, 2048, "chr1", 19400, "100S100M100S", 60, 2048, "chr1", 19700, "200S100M", 60},
            {"cpxLowMiddle", 0, "chr1", 20000, "100M400S", 60, 2048, "chr1", 20400, "100S200M200S", 30, 2048, "chr1", 20900, "400S100M", 60},
            {"cpxSmallGone", 0, "chr1", 21400, "100M300S", 60, 2048, "chr2", 2500, "100S200M100S", 30, 2048, "chr1", 21510, "300S100M", 60},
            {"cpxTiny", 0, "chr1", 21800, "100M200S", 60, 2048, "chr1", 21910, "100S100M100S", 60, 2048, "chr1", 22000, "200S100M", 60},
            {"incompleteContig", 0, "chr1", 22200, "100M200S", 60, 2048, "chr2", 6000, "100S100M100S", 60, 2048, "chr2", 6500, "200S100M", 60},
            {"incompleteStrand", 0, "chr1", 22400, "100M200S", 60, 2048, "chr1", 22600, "100S100M100S", 60, 2064, "chr1", 22800, "200S100M", 60},
            {"incompleteOrder", 0, "chr1", 23200, "100M200S", 60, 2048, "chr1", 23000, "100S100M100S", 60, 2048, "chr1", 22900, "200S100M", 60},
            {"incompleteMiddle", 0, "chr1", 23400, "100M200S", 60, 2048, "chr1", 23380, "100S100M100S", 60, 2048, "chr1", 23700, "200S100M", 60},
            {"uninformative", 0, "chr1", 24000, "200M", 60, 2048, "chr1", 24500, "50S100M50S", 30},
            {"ambiguous", 0, "chr1", 25000, "100M100S", 60, 2048, "chr2", 7000, "100M100S", 60, 2048, "chr1", 25500, "100S100M", 60},
            {"altTig", 0, "chr1", 26000, "100M100S", 60, 2048, "chr1", 26500, "100S100M", 60, 2048, "chr2", 8000, "200M", 60},
            {"manyLow", 0, "chr1", 27000, "100M1000S", 60, 2048, "chr1", 27200, "100S100M900S", 60, 2048, "chr1", 27400, "200S100M800S", 5,
                    2048, "chr1", 27600, "300S100M700S", 60, 2048, "chr1", 27800, "400S100M600S", 5, 2048, "chr1", 28000, "500S100M500S", 60,
                    2048, "chr1", 28200, "600S100M400S", 5, 2048, "chr1", 28400, "700S100M300S", 60, 2048, "chr1", 28600, "800S100M200S", 5,
                    2048, "chr1", 28800, "900S100M100S", 5, 2048, "chr1", 29000, "1000S100M", 60},
            {"gappedDropped", 0, "chr2", 1000, "100M60D100M100S", 60, 2048, "chr2", 1500, "150S50M100S", 60, 2048, "chr2", 1700, "200S100M", 60},
            {"cpxSplit", 0, "chr2", 2000, "100M60D100M200S", 60, 2048, "chr2", 2500, "200S100M100S", 60, 2048, "chr2", 2900, "300S100M", 60},
            {"interForward", 0, "chr1", 29200, "100M100S", 60, 2064, "chr2", 10000, "100M100S", 60},
            {"interForwardIns", 0, "chr1", 29400, "100M130S", 60, 2064, "chr2", 10200, "100M130S", 60},
            {"interReverse", 16, "chr1", 29600, "100S100M", 60, 2048, "chr2", 10400, "100S100M", 60},
            {"interReverseIns", 16, "chr2", 10600, "130S100M", 60, 2048, "chr1", 29800, "130S100M", 60},
            {"interSecond", 0, "chr2", 11000, "100M100S", 60, 2048, "chr1", 28800, "100S100M", 60},
            {"interSecondRev", 16, "chr1", 28400, "100M100S", 60, 2064, "chr2", 11200, "100S100M", 60},
            {"interFirstIns", 0, "chr1", 28200, "100M130S", 60, 2048, "chr2", 11400, "130S100M", 60},
            {"swapIns", 0, "chr2", 12000, "100M130S", 60, 2048, "chr2", 11700, "130S100M", 60},
            {"swapRev", 16, "chr2", 12300, "100M100S", 60, 2064, "chr2", 12600, "100S100M", 60},
            {"fatInsertion", 0, "chr2", 13000, "100M160S", 60, 2048, "chr2", 13120, "160S100M", 60},
            {"fatInsertionRev", 16, "chr2", 13400, "100M160S", 60, 2064, "chr2", 13520, "160S100M", 60},
            {"dupWide", 0, "chr2", 14000, "100M100S", 60, 2048, "chr2", 14040, "100S100M", 60},
            {"dupWideIns", 0, "chr2", 14300, "100M130S", 60, 2048, "chr2", 14340, "130S100M", 60},
            {"dupComplexWide", 0, "chr2", 14600, "100M90S", 60, 2048, "chr2", 14640, "90S100M", 60},
            {"dupComplexDown", 0, "chr2", 15000, "100M50S", 60, 2048, "chr2", 15040, "50S100M", 60},
            {"dupComplexShort", 0, "chr2", 15300, "100M90S", 60, 2048, "chr2", 15380, "90S100M", 60},
            {"weakPair", 0, "chr2", 16000, "100M100S", 15, 2048, "chr2", 16500, "100S100M", 60},
            {"shortPair", 0, "chr2", 16800, "20M100S", 60, 2048, "chr2", 17000, "20S100M", 60},
            {"lengthFiltered", 0, "chr2", 17300, "40M40S", 60, 2048, "chr2", 18000, "40S40M", 60},
            {"qualityFiltered", 0, "chr2", 18300, "100M100S", 25, 2048, "chr2", 19000, "100S100M", 25},
            {"normal", 0, "chr2", 19200, "100M", 60},
        };
        svContigBam(dir.resolve("sv_assembly.bam"), sequences, "s", contigs, 19L, true);
        final Object[][] tooMany = new Object[1][1 + 5 * 31];
        tooMany[0][0] = "tooMany";
        for (int at = 0; at < 31; at++) {
            tooMany[0][1 + 5 * at] = at == 0 ? 0 : 2048;
            tooMany[0][2 + 5 * at] = "chr1";
            tooMany[0][3 + 5 * at] = 100 + 300 * at;
            tooMany[0][4 + 5 * at] = at + "S20M" + (30 - at) + "S";
            tooMany[0][5 + 5 * at] = 60;
        }
        tooMany[0][4] = "20M30S";
        tooMany[0][154] = "30S20M";
        svContigBam(dir.resolve("sv_assembly_too_many.bam"), sequences, "s", tooMany, 20L, true);
        svContigBam(dir.resolve("sv_assembly_overlap_contained.bam"), sequences, "s", new Object[][] {
            {"halfOverlap", 0, "chr1", 1000, "100M400S", 60, 2048, "chr1", 1200, "90S210M200S", 60, 2048, "chr1", 1600, "200S120M180S", 60, 2048, "chr1", 2000, "400S100M", 60}}, 21L, true);
        // No breakend among its calls, or `--cnv-calls` stops on the first one's negative end.
        svContigBam(dir.resolve("sv_assembly_no_breakend.bam"), sequences, "s", new Object[][] {
            {"del", 0, "chr1", 1000, "100M100S", 60, 2048, "chr1", 2000, "100S100M", 60},
            {"ins", 0, "chr1", 10000, "100M160S", 60, 2048, "chr1", 10100, "160S100M", 60},
            {"dup", 0, "chr1", 11000, "100M100S", 60, 2048, "chr1", 11080, "100S100M", 60},
            {"cpxGapped", 0, "chr1", 6000, "100M320S", 60, 2048, "chr1", 6400, "120S100M200S", 60, 2048, "chr1", 6900, "320S100M", 60},
            {"cpxOverlap", 0, "chr1", 8000, "110M200S", 60, 2048, "chr1", 8300, "100S110M100S", 60, 2048, "chr1", 8700, "200S110M", 60},
            {"cpxFour", 0, "chr1", 12000, "100M400S", 60, 2048, "chr1", 12300, "100S110M290S", 60, 2048, "chr1", 12700, "200S110M190S", 60, 2048, "chr1", 13300, "400S100M", 60},
            {"cpxAsIsWide", 0, "chr1", 16600, "100M500S", 60, 2048, "chr1", 16699, "200S200M200S", 60, 2048, "chr1", 16898, "500S100M", 60}}, 22L, true);
        Files.writeString(dir.resolve("sv_alt_tigs.txt"), "chr2\n", StandardCharsets.UTF_8);

        // A read with no qualities, which the wellformed filter drops, beside one it keeps.
        final SAMFileHeader header = new SAMFileHeader();
        sequences.forEach(header::addSequence);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("s");
        header.addReadGroup(group);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        try (SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true)
                .makeBAMWriter(header, true, dir.resolve("sv_quality_missing.bam").toFile())) {
            for (int at = 0; at < 2; at++) {
                final htsjdk.samtools.SAMRecord read = new htsjdk.samtools.SAMRecord(header);
                read.setReadName("read" + at);
                read.setReferenceName("chr1");
                read.setAlignmentStart(100 + at);
                read.setCigarString("4M");
                read.setMappingQuality(60);
                read.setReadBases("ACGT".getBytes(StandardCharsets.US_ASCII));
                read.setBaseQualities(at == 0 ? new byte[] {30, 30, 30, 30} : htsjdk.samtools.SAMRecord.NULL_QUALS);
                read.setAttribute("RG", "rg1");
                writer.addAlignment(read);
            }
        }
    }

    /**
     * What `CpxVariantReInterpreterSpark` reads: `sv_assembly.bam`'s complex calls as the assembly
     * step wrote them, each haplotype reduced to its length, since the length is all it asks of one,
     * and five records of no contig's making: none or one segment, three contigs, one contig twice.
     */
    static void svComplexCalls(final Path dir) throws Exception {
        // Contig, start, end, contig names, segments, arrangement, haplotype length, qualities, longest.
        final Object[][] calls = {
            {"chr1", 1099, 1500, "cpxMapped", "chr1:1099-1500", "-chr2:500-599", 102, "60", 100},
            {"chr1", 2099, 2600, "cpxInverted", "chr1:2099-2300,chr1:2300-2399,chr1:2399-2600", "-2", 102, "60", 100},
            {"chr1", 3099, 3299, "cpxDuplicated", "chr1:3100-3200,chr1:3200-3299", "2,1,2", 301, "60", 100},
            {"chr1", 4499, 5000, "cpxReverse", "chr1:4499-4700,chr1:4700-4799,chr1:4799-5000", "2", 102, "60", 100},
            {"chr1", 5399, 5900, "cpxReverseFlip", "chr1:5399-5600,chr1:5600-5699,chr1:5699-5900", "-2", 102, "60", 100},
            {"chr1", 6099, 6900, "cpxGapped", "chr1:6099-6400,chr1:6400-6499,chr1:6499-6900", "UINS-20,2,UINS-100", 222, "60", 100},
            {"chr1", 8109, 8700, "cpxOverlap", "chr1:8109-8310,chr1:8310-8399,chr1:8399-8700", "2", 92, "60", 110},
            {"chr1", 9599, 10800, "cpxFour", "chr1:9599-9800,chr1:9800-9899,chr1:9899-10200,chr1:10200-10309,chr1:10309-10800", "2,4,UINS-90", 302, "60", 100},
            {"chr1", 11099, 11800, "cpxFourSwitch", "chr1:11099-11200,chr1:11200-11289,chr1:11289-11500,chr1:11500-11609,chr1:11609-11800", "2,-4,UINS-100", 302, "60", 100},
            {"chr1", 12099, 12800, "cpxFourOther", "chr1:12099-12300,chr1:12300-12409,chr1:12409-12800", "-chr2:1100-1199,2,UINS-90", 302, "60", 100},
            {"chr1", 13112, 13800, "cpxClipDeletion", "chr1:13112-13410,chr1:13410-13506,chr1:13506-13800", "UINS-3,2", 102, "60", 103},
            {"chr1", 14099, 14099, "cpxShared", "chr1:14099-14099", "1,chr2:3000-3099,1", 102, "60", 100},
            {"chr1", 14599, 14599, "cpxSharedFlip", "chr1:14599-14599", "1,UINS-20,-chr2:3400-3499,UINS-30,1", 152, "60", 100},
            {"chr1", 15099, 15298, "cpxInvertedOne", "chr1:15099-15298", "-1", 202, "60", 100},
            {"chr1", 15699, 15898, "cpxInvertedOneWide", "chr1:15699-15898", "UINS-100,-1,UINS-100", 402, "60", 100},
            {"chr1", 16299, 16318, "cpxInvertedSmall", "chr1:16299-16318", "UINS-120,-1", 142, "60", 100},
            {"chr1", 16699, 16898, "cpxAsIsWide", "chr1:16699-16898", "UINS-100,1,UINS-100", 402, "60", 100},
            {"chr1", 17299, 17498, "cpxAsIsFront", "chr1:17299-17498", "UINS-100,1", 302, "60", 100},
            {"chr1", 17899, 18098, "cpxAsIsRear", "chr1:17899-18098", "1,UINS-100", 302, "60", 100},
            {"chr1", 18499, 18797, "cpxAllAsIs", "chr1:18499-18648,chr1:18648-18797", "UINS-100,1,2,UINS-100", 502, "60", 100},
            {"chr1", 19299, 19700, "cpxSegmentGone", "chr1:19299-19400,chr1:19400-19499,chr1:19499-19700", "2", 102, "60", 100},
            {"chr1", 20099, 20900, "cpxLowMiddle", "chr1:20099-20400,chr1:20400-20599,chr1:20599-20900", "2,UINS-100", 302, "60", 100},
            {"chr1", 21499, 21510, "cpxSmallGone", "chr1:21499-21510", "-chr2:2500-2699", 202, "60", 100},
            {"chr1", 21899, 22009, "cpxTiny", "chr1:21899-21910,chr1:21910-22000,chr1:22000-22009", "2,3,3", 111, "60", 100},
            {"chr1", 25000, 25000, "zeroSegments", "", "UINS-100", 102, "60", 100},
            {"chr1", 25200, 25200, "zeroSegmentsShort", "", "UINS-10", 12, "60", 10},
            {"chr1", 25400, 25900, "oneSegmentDeleted", "chr1:25400-25900", "UINS-200", 200, "60", 100},
            {"chr1", 26000, 26300, "cpxFour,cpxFourOther,nobody", "chr1:26000-26100,chr1:26100-26300", "2,1", 300, "60,50,40", 120},
            {"chr1", 27099, 29000, "manyLow", "chr1:27099-27200,chr1:27200-27299,chr1:27299-27600,chr1:27600-27699,chr1:27699-28000,chr1:28000-28099,chr1:28099-28400,chr1:28400-28499,chr1:28499-29000", "2,UINS-100,4,UINS-100,6,UINS-100,8,UINS-200", 902, "60", 100},
            {"chr2", 2099, 2900, "cpxSplit", "chr2:2099-2160,chr2:2160-2259,chr2:2259-2500,chr2:2500-2599,chr2:2599-2900", "2,4", 202, "60", 100},
            {"chr2", 4000, 4600, "cpxGapped", "chr2:4000-4300,chr2:4300-4600", "-2,1", 600, "30", 80},
        };
        final StringBuilder lines = new StringBuilder("##fileformat=VCFv4.2\n");
        lines.append("##INFO=<ID=ALT_ARRANGEMENT,Number=.,Type=String,Description=\"For CPX variants only; specifies how reference segments given in SEGMENTS are re-arranged\">\n");
        lines.append("##INFO=<ID=CTG_NAMES,Number=.,Type=String,Description=\"Name of local assembly contigs supporting this variant\">\n");
        lines.append("##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">\n");
        lines.append("##INFO=<ID=MAPPING_QUALITIES,Number=.,Type=Integer,Description=\"Mapping qualities of the contig alignments that support the variant\">\n");
        lines.append("##INFO=<ID=MAX_ALIGN_LENGTH,Number=1,Type=Integer,Description=\"Maximum of the values listed in ALIGN_LENGTHS\">\n");
        lines.append("##INFO=<ID=SEGMENTS,Number=.,Type=String,Description=\"For CPX variants only; segments of reference that are rearranged\">\n");
        lines.append("##INFO=<ID=SEQ_ALT_HAPLOTYPE,Number=A,Type=Character,Description=\"Alt haplotype sequence, one per alt allele\">\n");
        lines.append("##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n");
        lines.append("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
        final String contigs = "##contig=<ID=chr1,length=30000>\n##contig=<ID=chr2,length=20000>\n";
        final String head = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        final StringBuilder records = new StringBuilder();
        for (final Object[] call : calls) {
            records.append(call[0]).append('\t').append(call[1]).append('\t')
                    .append("CPX_").append(call[0]).append(':').append(call[1]).append('-').append(call[2])
                    .append("\tN\t<CPX>\t.\t.\tALT_ARRANGEMENT=").append(call[5])
                    .append(";CTG_NAMES=").append(call[3]).append(";END=").append(call[2])
                    .append(";MAPPING_QUALITIES=").append(call[7]).append(";MAX_ALIGN_LENGTH=").append(call[8]);
            if (!((String) call[4]).isEmpty()) {
                records.append(";SEGMENTS=").append(call[4]);
            }
            records.append(";SEQ_ALT_HAPLOTYPE=").append("A".repeat((Integer) call[6]))
                    .append(";SVTYPE=CPX\n");
        }
        final String record = "chr1\t100\tbad\tN\t<CPX>\t.\t.\t";
        Files.writeString(dir.resolve("sv_cpx.vcf"), lines + contigs + head + records, StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_cpx_empty.vcf"), lines + contigs + head, StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_cpx_no_dictionary.vcf"), lines + head + records, StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_cpx_no_dictionary_one.vcf"),
                lines + head + records.substring(0, records.indexOf("\n") + 1), StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_cpx_sparse.vcf"), lines + contigs + head + record
                + "END=200;SEGMENTS=chr1:100-200;SVTYPE=CPX\n", StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_cpx_bad_segment.vcf"), lines + contigs + head + record
                + "ALT_ARRANGEMENT=1;CTG_NAMES=del;END=200;SEGMENTS=chr1:x-5;SEQ_ALT_HAPLOTYPE=AAA;SVTYPE=CPX\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_cpx_whole_contig.vcf"), lines + contigs + head + record
                + "ALT_ARRANGEMENT=1;CTG_NAMES=del;END=200;SEGMENTS=chr1;SEQ_ALT_HAPLOTYPE=AAA;SVTYPE=CPX\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_cpx_bad_length.vcf"), lines + contigs + head + record
                + "ALT_ARRANGEMENT=1;CTG_NAMES=del;END=200;MAX_ALIGN_LENGTH=many;SEGMENTS=chr1:100-200;SEQ_ALT_HAPLOTYPE=AAA;SVTYPE=CPX\n",
                StandardCharsets.UTF_8);
    }

    static void svContigBam(final Path bam, final java.util.List<SAMSequenceRecord> sequences, final String sample,
                            final Object[][] contigs, final long seed) {
        svContigBam(bam, sequences, sample, contigs, seed, false);
    }

    static void svContigBam(final Path bam, final java.util.List<SAMSequenceRecord> sequences, final String sample,
                            final Object[][] contigs, final long seed, final boolean scoreEachAlignment) {
        final SAMFileHeader header = new SAMFileHeader();
        sequences.forEach(header::addSequence);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        if (sample != null) {
            group.setSample(sample.split(",")[0]);
        }
        header.addReadGroup(group);
        if (sample != null && sample.contains(",")) {
            final SAMReadGroupRecord second = new SAMReadGroupRecord("rg2");
            second.setSample(sample.split(",")[1]);
            header.addReadGroup(second);
        }
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final java.util.Random random = new java.util.Random(seed);
        try (SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, bam.toFile())) {
            final java.util.List<htsjdk.samtools.SAMRecord> reads = new java.util.ArrayList<>();
            for (final Object[] contig : contigs) {
                for (int at = 1; at < contig.length; at += 5) {
                    final htsjdk.samtools.SAMRecord read = new htsjdk.samtools.SAMRecord(header);
                    read.setReadName((String) contig[0]);
                    read.setFlags((Integer) contig[at]);
                    read.setReferenceName((String) contig[at + 1]);
                    read.setAlignmentStart((Integer) contig[at + 2]);
                    read.setCigarString((String) contig[at + 3]);
                    read.setMappingQuality((Integer) contig[at + 4]);
                    final int length = read.getCigar().getReadLength();
                    final byte[] bases = new byte[length];
                    final byte[] qualities = new byte[length];
                    for (int base = 0; base < length; base++) {
                        bases[base] = (byte) "ACGT".charAt(random.nextInt(4));
                        qualities[base] = 30;
                    }
                    read.setReadBases(bases);
                    read.setBaseQualities(qualities);
                    read.setAttribute("RG", "rg1");
                    if (at == 1 || scoreEachAlignment) {
                        read.setAttribute("NM", 2);
                        read.setAttribute("AS", scoreEachAlignment ? read.getCigar().getReferenceLength() : 90);
                    }
                    reads.add(read);
                }
            }
            reads.sort(new htsjdk.samtools.SAMRecordCoordinateComparator());
            reads.forEach(writer::addAlignment);
        }
    }

    static void svEvidenceBam(final Path bam, final int template, final long seed, final boolean lone) {
        final SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("chr1", 30000));
        header.addSequence(new SAMSequenceRecord("chr2", 20000));
        header.addSequence(new SAMSequenceRecord("alt", 5000));
        for (final String[] group : new String[][] {{"rg1", "libA"}, {"rg2", "libB"}, {"rg3", null}}) {
            final SAMReadGroupRecord record = new SAMReadGroupRecord(group[0]);
            record.setSample("s");
            record.setLibrary(group[1]);
            header.addReadGroup(record);
        }
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final java.util.Random random = new java.util.Random(seed);
        final java.util.List<htsjdk.samtools.SAMRecord> reads = new java.util.ArrayList<>();
        final String[] groups = {"rg1", "rg2", "rg3"};
        final java.util.function.BiFunction<String, Integer, htsjdk.samtools.SAMRecord> read = (name, flags) -> {
            final htsjdk.samtools.SAMRecord record = new htsjdk.samtools.SAMRecord(header);
            record.setReadName(name);
            record.setFlags(flags);
            record.setMappingQuality(60);
            final byte[] bases = new byte[100];
            final byte[] quals = new byte[100];
            for (int at = 0; at < 100; at++) {
                bases[at] = (byte) "ACGT".charAt(random.nextInt(4));
                quals[at] = 30;
            }
            record.setReadBases(bases);
            record.setBaseQualities(quals);
            record.setAttribute("RG", groups[random.nextInt(3)]);
            return record;
        };
        // One read at a place, with its mate's place and its template size.
        final java.util.function.Consumer<Object[]> add = row -> {
            final htsjdk.samtools.SAMRecord record = read.apply((String) row[0], (Integer) row[1]);
            record.setReferenceName((String) row[2]);
            record.setAlignmentStart((Integer) row[3]);
            record.setCigarString((String) row[4]);
            if (row[5] != null) {
                record.setMateReferenceName((String) row[5]);
                record.setMateAlignmentStart((Integer) row[6]);
            }
            record.setInferredInsertSize((Integer) row[7]);
            for (int at = 8; at + 1 < row.length; at += 2) {
                record.setAttribute((String) row[at], row[at + 1]);
            }
            reads.add(record);
        };
        // The steady pairs, two reads per 40 bases over both big contigs.
        for (final String contig : new String[] {"chr1", "chr2"}) {
            for (int start = 200; start < 15000; start += 40) {
                final int size = template + random.nextInt(41) - 20;
                final String name = contig + "p" + start;
                add.accept(new Object[] {name, 99, contig, start, "100M", contig, start + size - 100, size, "MC", "100M", "MQ", 60});
                add.accept(new Object[] {name, 147, contig, start + size - 100, "100M", contig, start, -size, "MC", "100M", "MQ", 60});
            }
        }
        // A pile to the peak depth, then a shallower shoulder against it.
        for (int copy = 0; copy < 120; copy++) {
            add.accept(new Object[] {"pile" + copy, 99, "chr1", 8000 + copy % 7, "100M", "chr1", 8200, 300});
        }
        for (int copy = 0; copy < 30; copy++) {
            add.accept(new Object[] {"shoulder" + copy, 99, "chr1", 8100 + copy, "100M", "chr1", 8300, 300});
        }
        // Long templates in one stretch, for the K-S windows.
        for (int start = 16000; start < 17000; start += 10) {
            final int size = template + 1500 + random.nextInt(50);
            add.accept(new Object[] {"long" + start, 99, "chr1", start, "100M", "chr1", start + size - 100, size});
        }
        // Clips: a good SA, two of them, a poor one, one on this very spot, and clips that do not count.
        add.accept(new Object[] {"clipA", 99, "chr1", 1000, "40S60M", "chr1", 1200, 300, "SA", "chr2,5000,+,60M40S,60,0;"});
        add.accept(new Object[] {"clipB", 83, "chr1", 1100, "60M40S", "chr1", 900, -300, "SA", "chr2,5100,-,40S60M,60,0;"});
        add.accept(new Object[] {"clipC", 99, "chr1", 1200, "35S65M", "chr1", 1400, 300, "SA", "chr2,5200,+,60M,60,0;alt,300,+,60M,60,0;"});
        add.accept(new Object[] {"clipD", 99, "chr1", 1300, "35S65M", "chr1", 1500, 300, "SA", "chr2,5300,+,60M,5,0;"});
        add.accept(new Object[] {"clipE", 99, "chr1", 1400, "35S65M", "chr1", 1600, 300, "SA", "chr1,1402,+,60M,60,0;"});
        add.accept(new Object[] {"clipF", 99, "chr1", 1500, "20S80M", "chr1", 1700, 300});
        add.accept(new Object[] {"clipG", 99, "chr1", 1600, "10H35S65M", "chr1", 1800, 300});
        add.accept(new Object[] {"clipH", 99, "chr1", 1700, "65M35S5H", "chr1", 1900, 300, "SA", "alt,400,-,35M,60,0;"});
        add.accept(new Object[] {"clipI", 99, "chr1", 1800, "32S36M32S", "chr1", 2000, 300});
        add.accept(new Object[] {"clipJ", 2147, "chr1", 1900, "40H100M", "chr1", 2100, 300, "SA", "chr1,1000,+,40M60S,60,0;"});
        final htsjdk.samtools.SAMRecord dull = read.apply("clipK", 99);
        dull.setReferenceName("chr1");
        dull.setAlignmentStart(1950);
        dull.setCigarString("40S60M");
        dull.setMateReferenceName("chr1");
        dull.setMateAlignmentStart(2150);
        dull.setInferredInsertSize(300);
        final byte[] poor = dull.getBaseQualities();
        java.util.Arrays.fill(poor, 12, 20, (byte) 5);
        dull.setBaseQualities(poor);
        reads.add(dull);
        // Indels big and small.
        add.accept(new Object[] {"insA", 99, "chr1", 2200, "30M50I20M", "chr1", 2400, 300});
        add.accept(new Object[] {"delA", 99, "chr1", 2300, "50M60D50M", "chr1", 2600, 400});
        add.accept(new Object[] {"delB", 99, "chr1", 2400, "50M20D50M", "chr1", 2700, 400});
        add.accept(new Object[] {"insB", 99, "chr1", 2500, "20M10N30M45I5M", "chr1", 2800, 300});
        // Discordant pairs.
        add.accept(new Object[] {"mateOff", 73, "chr1", 3000, "100M", null, 0, 0});
        add.accept(new Object[] {"mateOff2", 89, "chr1", 3100, "2M1I97M", null, 0, 0, "MD", "A1C97"});
        add.accept(new Object[] {"interA", 97, "chr1", 3200, "100M", "chr2", 6000, 0, "MC", "50M50S", "MQ", 60});
        add.accept(new Object[] {"interB", 81, "chr1", 3300, "100M", "chr2", 6100, 0, "MQ", 3});
        add.accept(new Object[] {"interC", 97, "chr1", 3400, "100M", "alt", 700, 0});
        add.accept(new Object[] {"interD", 97, "chr2", 15500, "100M", "chr1", 3500, 0});
        add.accept(new Object[] {"interE", 97, "alt", 800, "100M", "chr1", 3600, 0});
        add.accept(new Object[] {"sameA", 65, "chr1", 3700, "100M", "chr1", 4000, 400, "MD", "100"});
        add.accept(new Object[] {"sameB", 113, "chr1", 3800, "100M", "chr1", 4100, 400, "MD", "T^AC98"});
        add.accept(new Object[] {"outA", 97, "chr1", 4200, "100M", "chr1", 3900, -300, "MD", "97T1G0"});
        add.accept(new Object[] {"outB", 81, "chr1", 4300, "100M", "chr1", 4600, 400});
        add.accept(new Object[] {"farA", 99, "chr1", 4400, "100M", "chr1", 9000, 4700, "MQ", 60});
        add.accept(new Object[] {"farB", 147, "chr1", 4500, "100M", "chr1", 100, -4500});
        add.accept(new Object[] {"nearA", 99, "chr1", 4600, "100M", "chr1", 4610, 110});
        add.accept(new Object[] {"second", 355, "chr1", 4700, "100M", "chr1", 4900, 300});
        add.accept(new Object[] {"supp", 2195, "chr1", 4800, "100M", "chr1", 5000, 300});
        add.accept(new Object[] {"dup", 1123, "chr1", 4900, "40S60M", "chr1", 5100, 300});
        add.accept(new Object[] {"fail", 611, "chr1", 5000, "40S60M", "chr1", 5200, 300});
        add.accept(new Object[] {"unmapped", 77, "*", 0, "*", null, 0, 0});
        // Reads the filter's thresholds decide.
        final htsjdk.samtools.SAMRecord low = read.apply("lowMapq", 99);
        low.setReferenceName("chr1");
        low.setAlignmentStart(5100);
        low.setCigarString("40S60M");
        low.setMateReferenceName("chr1");
        low.setMateAlignmentStart(5300);
        low.setInferredInsertSize(300);
        low.setMappingQuality(30);
        reads.add(low);
        add.accept(new Object[] {"shortMatch", 99, "chr1", 5200, "50S50M", "chr1", 5400, 300, "SA", "chr2,5400,+,50M50S,60,0;"});
        add.accept(new Object[] {"altClip", 99, "alt", 900, "40S60M", "alt", 1100, 300, "SA", "chr1,5500,+,60M40S,60,0;"});
        if (lone) {
            add.accept(new Object[] {"alone", 0, "chr1", 5500, "100M", null, 0, 0});
        }
        // Mapped or not, a read with no mate place has its unmapped flags right.
        for (final htsjdk.samtools.SAMRecord record : reads) {
            if (record.getReadUnmappedFlag()) {
                record.setMappingQuality(0);
            }
            if (record.getReadPairedFlag() && record.getMateUnmappedFlag()) {
                record.setMateReferenceName(record.getReferenceName());
                record.setMateAlignmentStart(record.getAlignmentStart());
            }
        }
        try (SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, bam.toFile())) {
            reads.forEach(writer::addAlignment);
        }
    }

    static void scoreRows() throws Exception {
        final String tool = "PathSeqScoreSpark";
        final String db = DIR.resolve("taxa.db").toString();
        final String pairedBam = DIR.resolve("score-paired.bam").toString();
        final String lone = DIR.resolve("score-unpaired.bam").toString();
        final String[][] runs = {
                {"score-both", "--paired-input", pairedBam, "--unpaired-input", lone},
                {"score-paired", "--paired-input", pairedBam},
                {"score-unpaired", "--unpaired-input", lone},
                {"score-paired-missing-other", "--paired-input", pairedBam, "--unpaired-input",
                        DIR.resolve("nowhere.bam").toString()},
                {"score-genome-length", "--paired-input", pairedBam, "--unpaired-input", lone,
                        "--divide-by-genome-length", "true"},
                {"score-no-kingdom", "--paired-input", pairedBam, "--unpaired-input", lone,
                        "--not-normalized-by-kingdom", "true"},
                {"score-identity", "--paired-input", pairedBam, "--unpaired-input", lone,
                        "--min-score-identity", "0.99"},
                {"score-identity-low", "--paired-input", pairedBam, "--unpaired-input", lone,
                        "--min-score-identity", "0.1"},
                {"score-margin", "--paired-input", pairedBam, "--unpaired-input", lone,
                        "--identity-margin", "0.3"},
                {"score-margin-zero", "--paired-input", pairedBam, "--unpaired-input", lone,
                        "--identity-margin", "0"},
                {"score-estimate", "--paired-input", pairedBam, "--score-reads-per-partition-estimate", "2"},
        };
        for (final String[] one : runs) {
            final String name = one[0];
            final List<String> rest = new ArrayList<>(Arrays.asList("--taxonomy-file", db, "--use-jdk-deflater", "true",
                    "--scores-output", out(name, name + ".txt"),
                    "--output", out(name, name + ".bam"),
                    "--score-metrics", out(name, name + "Metrics.txt"),
                    "--score-warnings", out(name, name + "-warnings.txt")));
            rest.addAll(Arrays.asList(one).subList(1, one.length));
            run(name, name + ".txt", tool, rest.toArray(new String[0]));
        }
        run("score-only-scores", "s.txt", tool, "--taxonomy-file", db, "--scores-output", out("score-only-scores", "s.txt"),
                "--unpaired-input", lone);
        final String[][] refused = {
                {"score-input", "-I", lone},
                {"score-no-input"},
                {"score-all-missing", "--paired-input", DIR.resolve("nowhere.bam").toString()},
                {"score-not-queryname", "--paired-input", DIR.resolve("score-not-queryname.bam").toString()},
                {"score-lonely", "--paired-input", DIR.resolve("score-lonely.bam").toString()},
                {"score-no-nm", "--unpaired-input", DIR.resolve("score-no-nm.bam").toString()},
                {"score-bad-xa", "--unpaired-input", DIR.resolve("score-bad-xa.bam").toString()},
                {"score-bad-nm", "--unpaired-input", DIR.resolve("score-bad-nm.bam").toString()},
                {"score-bad-cigar", "--unpaired-input", DIR.resolve("score-bad-cigar.bam").toString()},
        };
        for (final String[] one : refused) {
            final List<String> rest = new ArrayList<>(Arrays.asList("--taxonomy-file", db,
                    "--scores-output", out(one[0], "s.txt")));
            rest.addAll(Arrays.asList(one).subList(1, one.length));
            run(one[0], null, tool, rest.toArray(new String[0]));
        }
        run("score-no-database", null, tool, "--taxonomy-file", DIR.resolve("nowhere.db").toString(),
                "--scores-output", out("score-no-database", "s.txt"), "--unpaired-input", lone);
    }

    /** PATHSEQFILTERSPARK: each step's arguments, and the two lists it hands the writer. */
    static void filterRows() throws Exception {
        final String tool = "PathSeqFilterSpark";
        final String bam = DIR.resolve("filter.bam").toString();
        final String aligned = DIR.resolve("filter-aligned.bam").toString();
        final String[][] runs = {
                {"filter-default", bam},
                {"filter-skip", bam, "--skip-quality-filters", "true"},
                {"filter-keep-duplicates", bam, "--filter-duplicates", "false"},
                {"filter-short-table", bam, "--filter-reads-per-partition", "3"},
                {"filter-adapter-strict", bam, "--max-adapter-mismatches", "0"},
                {"filter-adapter-loose", bam, "--max-adapter-mismatches", "4", "--min-adapter-length", "8"},
                {"filter-min-length", bam, "--min-clipped-read-length", "75"},
                {"filter-masked", bam, "--max-masked-bases", "0"},
                {"filter-masked-many", bam, "--max-masked-bases", "40"},
                {"filter-base-quality", bam, "--min-base-quality", "31", "--max-masked-bases", "100"},
                {"filter-trim", bam, "--quality-threshold", "31"},
                {"filter-dust-off", bam, "--dust-t", "1000"},
                {"filter-dust-tight", bam, "--dust-t", "1", "--max-masked-bases", "100"},
                {"filter-dust-window", bam, "--dust-window", "16", "--max-masked-bases", "100"},
                {"filter-dust-quality", bam, "--dust-mask-quality", "9", "--max-masked-bases", "100",
                        "--min-base-quality", "5"},
                {"filter-repartition", bam, "--skip-pre-bwa-repartition", "true"},
                {"filter-hss", bam, "--kmer-file", DIR.resolve("host.hss").toString()},
                {"filter-bfi", bam, "--kmer-file", DIR.resolve("host.bfi").toString()},
                {"filter-kmer-thresh", bam, "--kmer-file", DIR.resolve("host.hss").toString(),
                        "--host-kmer-thresh", "30"},
                {"filter-kmer-extension", bam, "--kmer-file", DIR.resolve("host.txt").toString()},
                {"filter-kmer-missing", bam, "--kmer-file", DIR.resolve("nowhere.hss").toString()},
                {"filter-aligned", aligned, "--is-host-aligned", "true", "--skip-quality-filters", "true"},
                {"filter-aligned-off", aligned, "--skip-quality-filters", "true"},
                {"filter-aligned-identity", aligned, "--is-host-aligned", "true", "--host-min-identity", "46",
                        "--skip-quality-filters", "true"},
                {"filter-aligned-ignored", aligned, "--is-host-aligned", "true",
                        "--ignore-alignment-contigs", "chr1"},
                {"filter-ignored-unaligned", aligned, "--ignore-alignment-contigs", "chr1",
                        "--skip-quality-filters", "true"},
                {"filter-empty", DIR.resolve("filter-empty.bam").toString()},
        };
        for (final String[] one : runs) {
            final String name = one[0];
            final List<String> rest = new ArrayList<>(Arrays.asList("-I", one[1], "--use-jdk-deflater", "true",
                    "--paired-output", out(name, name + "-paired.bam"),
                    "--unpaired-output", out(name, name + "-unpaired.bam"),
                    "--filter-metrics", out(name, name + "Metrics.txt")));
            rest.addAll(Arrays.asList(one).subList(2, one.length));
            run(name, name + "-paired.bam", tool, rest.toArray(new String[0]));
        }
        run("filter-bwa", null, tool, "-I", bam, "--paired-output", out("filter-bwa", "p.bam"),
                "--filter-bwa-image", DIR.resolve("nowhere.img").toString());
        run("filter-no-paired", null, tool, "-I", bam, "--use-jdk-deflater", "true", "--unpaired-output", out("filter-no-paired", "u.bam"));
        run("filter-no-unpaired", null, tool, "-I", bam, "--use-jdk-deflater", "true", "--paired-output", out("filter-no-unpaired", "p.bam"));
    }

    /** FINDBADGENOMICKMERSSPARK writes what it collected, in the order its tables hold it. */
    static void kmerRows() throws Exception {
        final String tool = "FindBadGenomicKmersSpark";
        final String ref = DIR.resolve("kmers.fasta").toString();
        final String[][] runs = {
                {"kmers-default"},
                {"kmers-21", "--k-size", "21"},
                {"kmers-dust", "--k-size", "21", "--kmer-max-dust-score", "3"},
                {"kmers-dust-off", "--k-size", "21", "--kmer-max-dust-score", "1000"},
                {"kmers-odd-small", "--k-size", "5"},
                {"kmers-long", "--k-size", "63"},
                {"kmers-65", "--k-size", "65"},
                {"kmers-high-copy", "--k-size", "21", "--high-copy-fasta",
                        DIR.resolve("high-copy.fasta").toString()},
                {"kmers-high-copy-gz", "--k-size", "21", "--high-copy-fasta",
                        DIR.resolve("high-copy.fasta.gz").toString()},
                {"kmers-high-copy-blank", "--high-copy-fasta",
                        DIR.resolve("high-copy-blank.fasta").toString()},
                {"kmers-high-copy-missing", "--high-copy-fasta",
                        DIR.resolve("nowhere.fasta").toString()},
                {"kmers-high-copy-even", "--k-size", "4000", "--high-copy-fasta",
                        DIR.resolve("high-copy.fasta").toString()},
                {"kmers-even", "--k-size", "20"},
                {"kmers-too-small", "--k-size", "3"},
                {"kmers-zero", "--k-size", "0"},
                {"kmers-negative", "--k-size", "-1"},
                {"kmers-no-dust", "--kmer-max-dust-score", "0"},
        };
        for (final String[] one : runs) {
            final List<String> rest = new ArrayList<>(Arrays.asList("-R", ref, "-O",
                    out(one[0], one[0] + ".txt")));
            rest.addAll(Arrays.asList(one).subList(1, one.length));
            run(one[0], one[0] + ".txt", tool, rest.toArray(new String[0]));
        }
        run("kmers-big", "kmers-big.txt", tool, "-R", DIR.resolve("kmers-big.fasta").toString(),
                "-O", out("kmers-big", "kmers-big.txt"), "--k-size", "31");
        run("kmers-reads", "kmers-reads.txt", tool, "-R", DIR.resolve("ref.fasta").toString(),
                "-I", DIR.resolve("reads.bam").toString(), "-O", out("kmers-reads", "kmers-reads.txt"));
        run("kmers-reads-mismatch", "kmers-reads-mismatch.txt", tool, "-R", ref,
                "-I", DIR.resolve("reads.bam").toString(),
                "-O", out("kmers-reads-mismatch", "kmers-reads-mismatch.txt"));
        run("kmers-no-directory", null, tool, "-R", ref,
                "-O", DIR.resolve("nowhere/kmers.txt").toString());
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

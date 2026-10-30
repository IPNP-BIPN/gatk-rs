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

        System.out.println("# SparkWalkersDump: eighteen Spark tools through Main, --spark-master local[1]");

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

        // The fixtures travel, so the port runs over these bytes rather than its own rebuild.
        for (final String name : new String[] {"ref.fasta", "ref.fasta.fai", "ref.dict",
                "reads.bam", "reads.bai", "plain.vcf", "genotyped.vcf", "revert.map.tsv",
                "revert.bad-map.tsv", "names.txt", "dup-a.bam", "dup-a.bai", "dup-marks.bam",
                "dup-marks.bai", "dup-renamed.bam", "dup-uneven.bam", "dup-moved.bam",
                "dup-no-group.bam", "dup-other-group.bam", "dup-fewer.bam"}) {
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
        if (!path.toString().endsWith(".metrics")) {
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

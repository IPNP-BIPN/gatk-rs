/*
 * The eleven Spark tools this port runs, end to end through Main, against `--spark-master local[1]`.
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

        System.out.println("# SparkWalkersDump: eleven Spark tools through Main, --spark-master local[1]");

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

        // The fixtures travel, so the port runs over these bytes rather than its own rebuild.
        for (final String name : new String[] {"ref.fasta", "ref.fasta.fai", "ref.dict",
                "reads.bam", "reads.bai", "plain.vcf", "genotyped.vcf"}) {
            System.out.printf("input\t%s\t%s%n", name,
                    Base64.getEncoder().encodeToString(Files.readAllBytes(DIR.resolve(name))));
        }

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
                    Base64.getEncoder().encodeToString(Files.readAllBytes(path)));
        }
    }

    static String masked(final String text) {
        return text.replace(Path.of("").toAbsolutePath().toString(), "<cwd>");
    }
}

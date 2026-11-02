/*
 * The eight tools that stop at an interpreter the pinned image does not have, through Main.
 *
 * The image has neither Rscript, nor `python`, nor gcnvkernel or sklearn. What the reference does
 * before it reaches them is all that can be observed there, and this dump is that much:
 *
 *   - THE PLOTS READ AND CHECK EVERYTHING FIRST: `RScriptExecutor` looks for its executable only
 *     when the script runs, so every input refusal comes before the missing Rscript, and the
 *     output directory exists by then;
 *   - THE GCNV TOOLS AND NVSCOREVARIANTS IMPORT THEIR PACKAGE IN `onStartup`, before `doWork`, so
 *     nothing but the parser can refuse a command line ahead of the missing python;
 *   - POSTPROCESSGERMLINECNVCALLS CALLS `super.onStartup()` FIRST, so GATKTool's reference, reads
 *     and `-L` come ahead of it, and `-L` reads the shards' interval lists;
 *   - TRAINVARIANTANNOTATIONSMODEL ONLY IMPORTS AFTER READING THE ANNOTATION NAMES, and JAVA_BGMM
 *     and PYTHON_SCRIPT read every matrix before their first refusal of a missing backend;
 *   - SCOREVARIANTANNOTATIONS IS A WALKER: GATKTool's startup, then `onTraversalStart` and the
 *     VCF's header, then the backend, and under PYTHON_SCRIPT a whole pass and the annotations
 *     file before the first script runs. A pass that extracts nothing never runs one and ends
 *     with a VCF.
 *
 * Output:
 *
 *     input\t<name>\t<base64 of a fixture file>
 *     folder\t<name>
 *     argv\t<label>\t<the command line, tab-separated, escaped>
 *     stdout\t<label>\t<what the tool printed, escaped>
 *     error\t<label>\t<exception class>:<message, escaped, cwd masked>
 *     dir\t<label>\t<a directory under the run's output>
 *     file\t<label>\t<path under the run's output>\t<base64>
 *     hdf5\t<label>\t<path under the run's output>\t<the annotations file's contents, as text>
 *
 * Unlike the Spark dump this one lists what a refused run left behind: what the reference writes
 * before it stops is the half of the answer the interpreter does not touch. An HDF5 file is
 * listed as text because two runs of the library write different bytes. Runs labelled `setup-`
 * make those files, and the replay makes them with the port.
 *
 * Usage: InterpreterToolsDump
 */

import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.LabeledVariantAnnotationsData;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.VariantAnnotationsScorer;

import java.io.ByteArrayOutputStream;
import java.io.File;
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

public class InterpreterToolsDump {

    static final Path DIR = Path.of("interpreter-dump");

    static final String SEQUENCES = "@SQ\tSN:chr1\tLN:2000000\n"
            + "@SQ\tSN:chr2\tLN:1500000\n"
            + "@SQ\tSN:chrShort\tLN:1000\n";
    static final String DICT = "@HD\tVN:1.6\n" + SEQUENCES;
    static final String HEADER = DICT + "@RG\tID:GATKCopyNumber\tSM:sampleA\n";
    static final String COPY_RATIOS = HEADER + "CONTIG\tSTART\tEND\tLOG2_COPY_RATIO\n"
            + "chr1\t1000\t1999\t0.5\n"
            + "chr1\t2000\t2999\t-0.25\n"
            + "chr2\t500\t999\t1.0\n";
    static final String SEGMENT_COLUMNS = "CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tNUM_POINTS_ALLELE_FRACTION\t"
            + "LOG2_COPY_RATIO_POSTERIOR_10\tLOG2_COPY_RATIO_POSTERIOR_50\tLOG2_COPY_RATIO_POSTERIOR_90\t"
            + "MINOR_ALLELE_FRACTION_POSTERIOR_10\tMINOR_ALLELE_FRACTION_POSTERIOR_50\tMINOR_ALLELE_FRACTION_POSTERIOR_90\n";
    static final String POSTERIORS = "\t0.1\t0.2\t0.3\t0.4\t0.45\t0.5\n";
    static final String SEGMENTS = HEADER + SEGMENT_COLUMNS
            + "chr1\t1000\t2999\t2\t2" + POSTERIORS
            + "chr2\t500\t999\t1\t1" + POSTERIORS;
    static final String ALLELIC_COUNTS = HEADER + "CONTIG\tPOSITION\tREF_COUNT\tALT_COUNT\tREF_NUCLEOTIDE\tALT_NUCLEOTIDE\n"
            + "chr1\t1500\t10\t12\tA\tC\n"
            + "chr1\t2500\t8\t9\tG\tT\n"
            + "chr2\t700\t4\t3\tT\tN\n";

    static final String VCF_HEADER = "##fileformat=VCFv4.2\n"
            + "##FILTER=<ID=LowQual,Description=\"Low quality\">\n"
            + "##INFO=<ID=QD,Number=1,Type=Float,Description=\"QD\">\n"
            + "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"MQ\">\n"
            + "##INFO=<ID=NEVER,Number=1,Type=Float,Description=\"Declared and on no record\">\n"
            + "##INFO=<ID=AS_QD,Number=A,Type=Float,Description=\"Allele-specific QD\">\n"
            + "##contig=<ID=chr1,length=2000000>\n"
            + "##contig=<ID=chr2,length=1500000>\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    static final String CALLS = VCF_HEADER
            + "chr1\t100\t.\tA\tC\t50\t.\tQD=1.5;MQ=60.0;AS_QD=1.5\n"
            + "chr1\t200\t.\tA\tG\t50\t.\tQD=2.5;MQ=50.0;AS_QD=2.5\n"
            + "chr1\t300\t.\tC\tT\t50\tLowQual\tQD=3.5;MQ=40.0;AS_QD=3.5\n"
            + "chr1\t400\t.\tCA\tC\t50\t.\tQD=4.5;MQ=30.0;AS_QD=4.5\n"
            + "chr1\t500\t.\tG\tGT\t50\t.\tQD=5.5;MQ=20.0;AS_QD=5.5\n"
            + "chr2\t100\t.\tG\tA,T\t50\t.\tQD=6.5;MQ=10.0;AS_QD=6.5,7.5\n";
    static final String RESOURCE = VCF_HEADER
            + "chr1\t100\t.\tA\tC\t50\t.\t.\n"
            + "chr1\t400\t.\tCA\tC\t50\t.\t.\n"
            + "chr2\t100\t.\tG\tA\t50\t.\t.\n";
    static final String RESOURCE_SNPS = VCF_HEADER
            + "chr1\t100\t.\tA\tC\t50\t.\t.\n";
    static final String RESOURCE_EVERY_SNP = VCF_HEADER
            + "chr1\t100\t.\tA\tC\t50\t.\t.\n"
            + "chr1\t200\t.\tA\tG\t50\t.\t.\n"
            + "chr1\t300\t.\tC\tT\t50\t.\t.\n"
            + "chr2\t100\t.\tG\tA\t50\t.\t.\n";

    public static void main(final String[] args) throws Exception {
        PrintReadsDump.emptyDirectory(DIR);
        Files.createDirectories(DIR);

        System.out.println("# InterpreterToolsDump: eight tools through Main, up to the interpreter");

        input("ref.dict", DICT);
        input("sequences-with-delimiter.dict", "@HD\tVN:1.6\n@SQ\tSN:chrCONTIG_DELIMITER1\tLN:2000000\n");
        input("only-chr1.dict", "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:2000000\n");
        input("short-chr1.dict", "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:2000\n@SQ\tSN:chr2\tLN:1500000\n");
        input("no-sequences.dict", "@HD\tVN:1.6\n");
        input("standardized.tsv", COPY_RATIOS);
        input("denoised.tsv", COPY_RATIOS.replace("0.5", "0.75"));
        input("denoised-other-intervals.tsv", COPY_RATIOS.replace("chr2\t500", "chr2\t600"));
        input("denoised-other-sample.tsv", COPY_RATIOS.replace("sampleA", "sampleB"));
        input("denoised-two-samples.tsv", COPY_RATIOS.replace("@RG\tID:GATKCopyNumber\tSM:sampleA\n",
                "@RG\tID:one\tSM:sampleA\n@RG\tID:two\tSM:sampleB\n"));
        input("denoised-no-sample.tsv", COPY_RATIOS.replace("@RG\tID:GATKCopyNumber\tSM:sampleA\n", ""));
        input("denoised-no-column.tsv", COPY_RATIOS.replace("LOG2_COPY_RATIO", "OTHER"));
        input("denoised-unsorted.tsv", HEADER + "CONTIG\tSTART\tEND\tLOG2_COPY_RATIO\n"
                + "chr1\t2000\t2999\t0.5\nchr1\t1000\t1999\t0.5\n");
        input("denoised-overlapping.tsv", HEADER + "CONTIG\tSTART\tEND\tLOG2_COPY_RATIO\n"
                + "chr1\t1000\t1999\t0.5\nchr1\t1500\t2999\t0.5\n");
        input("denoised-off-dictionary.tsv", HEADER + "CONTIG\tSTART\tEND\tLOG2_COPY_RATIO\n"
                + "chrOther\t1000\t1999\t0.5\n");
        input("denoised-bad-number.tsv", HEADER + "CONTIG\tSTART\tEND\tLOG2_COPY_RATIO\n"
                + "chr1\t1000\t1999\tmany\n");
        input("denoised-bad-start.tsv", HEADER + "CONTIG\tSTART\tEND\tLOG2_COPY_RATIO\n"
                + "chr1\tfirst\t1999\t0.5\n");
        input("denoised-short-line.tsv", HEADER + "CONTIG\tSTART\tEND\tLOG2_COPY_RATIO\n"
                + "chr1\t1000\t1999\n");
        input("empty.tsv", "");
        input("segments.seg", SEGMENTS);
        input("segments-other-sample.seg", SEGMENTS.replace("sampleA", "sampleB"));
        input("segments-wrong-copy-ratio-points.seg", HEADER + SEGMENT_COLUMNS
                + "chr1\t1000\t2999\t3\t2" + POSTERIORS + "chr2\t500\t999\t1\t1" + POSTERIORS);
        input("segments-wrong-allele-points.seg", HEADER + SEGMENT_COLUMNS
                + "chr1\t1000\t2999\t2\t3" + POSTERIORS + "chr2\t500\t999\t1\t1" + POSTERIORS);
        input("segments-extra-contig.seg", HEADER + SEGMENT_COLUMNS
                + "chr1\t1000\t2999\t2\t2" + POSTERIORS + "chr2\t500\t999\t1\t1" + POSTERIORS
                + "chrShort\t1\t900\t0\t0" + POSTERIORS);
        input("allelic.tsv", ALLELIC_COUNTS);
        input("allelic-fewer.tsv", ALLELIC_COUNTS.replace("chr2\t700\t4\t3\tT\tN\n", ""));

        input("hyperparameters.json", "{}\n");
        input("score.py", "# never run: the image has no python\n");
        input("not-hdf5.hdf5", "This is not an HDF5 file.\n");
        input("calls.vcf", CALLS);
        input("resource.vcf", RESOURCE);
        input("resource-unindexed.vcf", RESOURCE);
        input("resource-snps.vcf", RESOURCE_SNPS);
        input("resource-every-snp.vcf", RESOURCE_EVERY_SNP);
        for (final String indexed : List.of("calls.vcf", "resource.vcf", "resource-snps.vcf", "resource-every-snp.vcf")) {
            new Main().instanceMain(new String[] {"IndexFeatureFile", "-I", DIR.resolve(indexed).toString()});
        }

        // The shards PostprocessGermlineCNVCalls reads `-L` against, and the model files
        // ScoreVariantAnnotations only asks can be read.
        final String shardIntervals = DICT + "CONTIG\tSTART\tEND\nchr1\t1000\t1999\nchr1\t2000\t2999\n";
        folder("calls-shard");
        input("calls-shard/interval_list.tsv", shardIntervals);
        folder("model-shard");
        input("model-shard/interval_list.tsv", shardIntervals);
        folder("model-other-dictionary");
        input("model-other-dictionary/interval_list.tsv", shardIntervals.replace("LN:2000000", "LN:2000001"));
        folder("empty-folder");
        input("both.snp.scorer.pkl", "no pickle\n");
        input("both.indel.scorer.pkl", "no pickle\n");
        input("snp-only.snp.scorer.pkl", "no pickle\n");
        input("indel-only.indel.scorer.pkl", "no pickle\n");
        input("ser.snp.scorer.ser", "no serialized scorer\n");
        input("bgmm.indel.bgmmScorer.ser", "a scorer the module cannot read\n");
        input("broken-calibration.snp.scorer.pkl", "no pickle\n");
        input("broken-calibration.indel.scorer.pkl", "no pickle\n");
        input("broken-calibration.snp.calibrationScores.hdf5", "This is not an HDF5 file.\n");

        // Scores written by the library here and by the port's writer in the replay.
        VariantAnnotationsScorer.writeScores(DIR.resolve("calibrated.snp.calibrationScores.hdf5").toFile(), new double[] {1., 2., 3.});
        VariantAnnotationsScorer.writeScores(DIR.resolve("infinite.snp.calibrationScores.hdf5").toFile(), new double[] {1., Double.POSITIVE_INFINITY});
        input("calibrated.snp.scorer.pkl", "no pickle\n");
        input("calibrated.indel.scorer.pkl", "no pickle\n");
        input("infinite.snp.scorer.pkl", "no pickle\n");
        input("infinite.indel.scorer.pkl", "no pickle\n");

        setupRows();
        plotDenoisedRows();
        plotSegmentsRows();
        ploidyRows();
        callerRows();
        postprocessRows();
        nvScoreRows();
        trainRows();
        scoreRows();
    }

    static String at(final String name) {
        return DIR.resolve(name).toString();
    }

    static String out(final String label, final String name) {
        return DIR.resolve("out").resolve(label).resolve(name).toString();
    }

    static void input(final String name, final String content) throws Exception {
        final byte[] bytes = content.getBytes(StandardCharsets.UTF_8);
        Files.write(DIR.resolve(name), bytes);
        System.out.printf("input\t%s\t%s%n", name, Base64.getEncoder().encodeToString(bytes));
    }

    static void folder(final String name) throws Exception {
        Files.createDirectories(DIR.resolve(name));
        System.out.printf("folder\t%s%n", name);
    }

    /** One run through Main, and everything it left under its output directory, refused or not. */
    static void run(final String label, final String tool, final String... rest) throws Exception {
        final List<String> argv = new ArrayList<>();
        argv.add(tool);
        argv.addAll(Arrays.asList(rest));
        System.out.printf("argv\t%s\t%s%n", label,
                ReferenceQueryDump.escape(String.join("\t", argv)));
        final Path outDir = DIR.resolve("out").resolve(label);
        Files.createDirectories(outDir);
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
        }
        final List<Path> produced;
        try (Stream<Path> walk = Files.walk(outDir)) {
            produced = walk.filter(path -> !path.equals(outDir)).sorted().collect(Collectors.toList());
        }
        for (final Path path : produced) {
            final Path relative = outDir.relativize(path);
            if (Files.isDirectory(path)) {
                System.out.printf("dir\t%s\t%s%n", label, relative);
            } else if (path.toString().endsWith(".hdf5")) {
                System.out.printf("hdf5\t%s\t%s\t%s%n", label, relative, annotations(path.toFile()));
            } else if (path.toString().endsWith(".idx")) {
                // A Tribble index holds the output's absolute path and its modification time.
                System.out.printf("index\t%s\t%s%n", label, relative);
            } else {
                System.out.printf("file\t%s\t%s\t%s%n", label, relative,
                        Base64.getEncoder().encodeToString(Files.readAllBytes(path)));
            }
        }
    }

    /** An annotations file as text; a label the file does not hold is a dash. */
    static String annotations(final File file) {
        final StringBuilder text = new StringBuilder();
        text.append("names=").append(LabeledVariantAnnotationsData.readAnnotationNames(file));
        text.append(";annotations=").append(Arrays.deepToString(LabeledVariantAnnotationsData.readAnnotations(file)));
        for (final String label : List.of("snp", "training", "calibration")) {
            String values;
            try {
                values = LabeledVariantAnnotationsData.readLabel(file, label).toString();
            } catch (final RuntimeException e) {
                values = "-";
            }
            text.append(';').append(label).append('=').append(values);
        }
        return text.toString();
    }

    static String masked(final String text) {
        return text.replace(Path.of("").toAbsolutePath().toString(), "<cwd>");
    }

    /** The annotation files the two walkers' rows read. */
    static void setupRows() throws Exception {
        final String[][] setups = {
                {"setup-labelled", "QD MQ", "train,training=true,calibration=true", "resource.vcf"},
                {"setup-other-names", "QD", "train,training=true,calibration=true", "resource.vcf"},
                {"setup-missing-values", "QD NEVER", "train,training=true,calibration=true", "resource.vcf"},
                {"setup-calibration-only", "QD MQ", "train,calibration=true", "resource.vcf"},
                {"setup-training-only", "QD MQ", "train,training=true", "resource.vcf"},
                {"setup-only-snps", "QD MQ", "train,training=true,calibration=true", "resource-snps.vcf"},
                {"setup-every-snp", "QD MQ", "train,training=true,calibration=true", "resource-every-snp.vcf"},
        };
        for (final String[] setup : setups) {
            final List<String> rest = new ArrayList<>(List.of(
                    "-V", at("calls.vcf"),
                    "--resource:" + setup[2], at(setup[3]),
                    "-O", out(setup[0], "extract"),
                    "--maximum-number-of-unlabeled-variants", "100",
                    "--add-output-vcf-command-line", "false",
                    "--do-not-gzip-vcf-output",
                    "--create-output-variant-index", "false"));
            for (final String annotation : setup[1].split(" ")) {
                rest.add("-A");
                rest.add(annotation);
            }
            run(setup[0], "ExtractVariantAnnotations", rest.toArray(new String[0]));
        }
    }

    static String[] plotDenoised(final String label, final String standardized, final String denoised,
                                 final String dictionary, final String... extra) {
        final List<String> rest = new ArrayList<>(List.of(
                "--standardized-copy-ratios", at(standardized),
                "--denoised-copy-ratios", at(denoised),
                "--sequence-dictionary", at(dictionary),
                "--output-prefix", "sample",
                "--output", out(label, "plots")));
        rest.addAll(Arrays.asList(extra));
        return rest.toArray(new String[0]);
    }

    static void denoised(final String label, final String denoised, final String... extra) throws Exception {
        run("plot-denoised-" + label, "PlotDenoisedCopyRatios",
                plotDenoised("plot-denoised-" + label, "standardized.tsv", denoised, "ref.dict", extra));
    }

    static void plotDenoisedRows() throws Exception {
        final String tool = "PlotDenoisedCopyRatios";
        denoised("happy", "denoised.tsv");
        denoised("tuned", "denoised.tsv", "--minimum-contig-length", "0", "--maximum-copy-ratio", "Infinity",
                "--point-size-copy-ratio", "1.0");
        denoised("no-such-file", "nothing.tsv");
        denoised("folder", "empty-folder");
        denoised("other-intervals", "denoised-other-intervals.tsv");
        denoised("other-sample", "denoised-other-sample.tsv");
        denoised("two-samples", "denoised-two-samples.tsv");
        denoised("no-sample", "denoised-no-sample.tsv");
        denoised("no-column", "denoised-no-column.tsv");
        denoised("unsorted", "denoised-unsorted.tsv");
        denoised("overlapping", "denoised-overlapping.tsv");
        denoised("off-dictionary", "denoised-off-dictionary.tsv");
        denoised("bad-number", "denoised-bad-number.tsv");
        denoised("bad-start", "denoised-bad-start.tsv");
        denoised("short-line", "denoised-short-line.tsv");
        denoised("empty-file", "empty.tsv");
        denoised("negative-length", "denoised.tsv", "--minimum-contig-length", "-1");
        denoised("negative-ratio", "denoised.tsv", "--maximum-copy-ratio", "-1");
        denoised("negative-point", "denoised.tsv", "--point-size-copy-ratio", "-0.5");
        denoised("every-contig-short", "denoised.tsv", "--minimum-contig-length", "2000001");
        run("plot-denoised-no-dictionary", tool, plotDenoised("plot-denoised-no-dictionary",
                "standardized.tsv", "denoised.tsv", "nothing.dict"));
        run("plot-denoised-empty-dictionary", tool, plotDenoised("plot-denoised-empty-dictionary",
                "standardized.tsv", "denoised.tsv", "no-sequences.dict"));
        run("plot-denoised-delimiter", tool, plotDenoised("plot-denoised-delimiter",
                "standardized.tsv", "denoised.tsv", "sequences-with-delimiter.dict"));
        run("plot-denoised-subset", tool, plotDenoised("plot-denoised-subset",
                "standardized.tsv", "denoised.tsv", "only-chr1.dict"));
        run("plot-denoised-exceeds", tool, plotDenoised("plot-denoised-exceeds",
                "standardized.tsv", "denoised.tsv", "short-chr1.dict"));
        run("plot-denoised-missing-prefix", tool,
                "--standardized-copy-ratios", at("standardized.tsv"), "--denoised-copy-ratios", at("denoised.tsv"),
                "--sequence-dictionary", at("ref.dict"), "--output", out("plot-denoised-missing-prefix", "plots"));
        run("plot-denoised-empty-prefix", tool,
                "--standardized-copy-ratios", at("standardized.tsv"), "--denoised-copy-ratios", at("denoised.tsv"),
                "--sequence-dictionary", at("ref.dict"), "--output-prefix", "",
                "--output", out("plot-denoised-empty-prefix", "plots"));
        run("plot-denoised-output-under-file", tool,
                "--standardized-copy-ratios", at("standardized.tsv"), "--denoised-copy-ratios", at("denoised.tsv"),
                "--sequence-dictionary", at("ref.dict"), "--output-prefix", "sample",
                "--output", at("ref.dict") + "/plots");
        run("plot-denoised-output-is-file", tool,
                "--standardized-copy-ratios", at("standardized.tsv"), "--denoised-copy-ratios", at("denoised.tsv"),
                "--sequence-dictionary", at("ref.dict"), "--output-prefix", "sample",
                "--output", at("empty.tsv"));
    }

    static void segments(final String label, final String... rest) throws Exception {
        final String name = "plot-segments-" + label;
        final List<String> argv = new ArrayList<>(List.of("--output-prefix", "sample", "--output", out(name, "plots")));
        argv.addAll(Arrays.asList(rest));
        run(name, "PlotModeledSegments", argv.toArray(new String[0]));
    }

    static void plotSegmentsRows() throws Exception {
        final String[] dict = {"--sequence-dictionary", at("ref.dict")};
        final String denoised = at("denoised.tsv");
        final String allelic = at("allelic.tsv");
        final String seg = at("segments.seg");
        segments("happy", "--denoised-copy-ratios", denoised, "--allelic-counts", allelic, "--segments", seg, dict[0], dict[1]);
        segments("ratios-only", "--denoised-copy-ratios", denoised, "--segments", seg, dict[0], dict[1]);
        segments("counts-only", "--allelic-counts", allelic, "--segments", seg, dict[0], dict[1]);
        segments("tuned", "--allelic-counts", allelic, "--segments", seg, dict[0], dict[1],
                "--point-size-allele-fraction", "1.5", "--minimum-contig-length", "0");
        segments("neither", "--segments", seg, dict[0], dict[1]);
        segments("no-segments-file", "--denoised-copy-ratios", denoised, "--segments", at("nothing.seg"), dict[0], dict[1]);
        segments("no-counts-file", "--allelic-counts", at("nothing.tsv"), "--segments", seg, dict[0], dict[1]);
        segments("other-sample", "--denoised-copy-ratios", denoised, "--segments", at("segments-other-sample.seg"), dict[0], dict[1]);
        segments("wrong-ratio-points", "--denoised-copy-ratios", denoised,
                "--segments", at("segments-wrong-copy-ratio-points.seg"), dict[0], dict[1]);
        segments("wrong-allele-points", "--allelic-counts", allelic,
                "--segments", at("segments-wrong-allele-points.seg"), dict[0], dict[1]);
        segments("fewer-counts", "--allelic-counts", at("allelic-fewer.tsv"), "--segments", seg, dict[0], dict[1]);
        segments("zero-points-elsewhere", "--denoised-copy-ratios", denoised, "--allelic-counts", allelic,
                "--segments", at("segments-extra-contig.seg"), dict[0], dict[1]);
        segments("no-dictionary", "--denoised-copy-ratios", denoised, "--segments", seg, dict[0], at("nothing.dict"));
        segments("exceeds", "--denoised-copy-ratios", denoised, "--segments", seg, dict[0], at("short-chr1.dict"));
        segments("negative-fraction", "--denoised-copy-ratios", denoised, "--segments", seg, dict[0], dict[1],
                "--point-size-allele-fraction", "-1");
        segments("missing-segments", "--denoised-copy-ratios", denoised, dict[0], dict[1]);
    }

    static void ploidyRows() throws Exception {
        final String tool = "DetermineGermlineContigPloidy";
        // None of these inputs exists: the import comes first, and the output directory is not made.
        run("ploidy-happy", tool, "--input", at("nothing.tsv"), "--output-prefix", "p", "--output", out("ploidy-happy", "made"));
        run("ploidy-empty-prefix", tool, "--input", at("nothing.tsv"), "--output-prefix", "", "--output", out("ploidy-empty-prefix", "made"));
        run("ploidy-duplicates", tool, "--input", at("nothing.tsv"), "--input", at("nothing.tsv"),
                "--output-prefix", "p", "--output", out("ploidy-duplicates", "made"));
        run("ploidy-missing-input", tool, "--output-prefix", "p", "--output", out("ploidy-missing-input", "made"));
        run("ploidy-missing-prefix", tool, "--input", at("nothing.tsv"), "--output", out("ploidy-missing-prefix", "made"));
        run("ploidy-missing-output", tool, "--input", at("nothing.tsv"), "--output-prefix", "p");
        run("ploidy-bad-number", tool, "--input", at("nothing.tsv"), "--output-prefix", "p", "--output", out("ploidy-bad-number", "made"),
                "--mean-bias-standard-deviation", "many");
        run("ploidy-below-range", tool, "--input", at("nothing.tsv"), "--output-prefix", "p", "--output", out("ploidy-below-range", "made"),
                "--mean-bias-standard-deviation", "-1");
        run("ploidy-unknown-argument", tool, "--input", at("nothing.tsv"), "--output-prefix", "p", "--output", out("ploidy-unknown-argument", "made"),
                "--no-such-argument", "1");
    }

    static void callerRows() throws Exception {
        final String tool = "GermlineCNVCaller";
        final String[] required = {"--input", at("nothing.tsv"), "--contig-ploidy-calls", at("nothing"), "--output-prefix", "p"};
        run("caller-cohort", tool, concat(required, "--run-mode", "COHORT", "--output", out("caller-cohort", "made")));
        run("caller-case", tool, concat(required, "--run-mode", "CASE", "--output", out("caller-case", "made")));
        run("caller-bad-mode", tool, concat(required, "--run-mode", "OTHER", "--output", out("caller-bad-mode", "made")));
        run("caller-missing-mode", tool, concat(required, "--output", out("caller-missing-mode", "made")));
        run("caller-missing-ploidy", tool, "--input", at("nothing.tsv"), "--output-prefix", "p", "--run-mode", "CASE",
                "--output", out("caller-missing-ploidy", "made"));
        run("caller-below-range", tool, concat(required, "--run-mode", "CASE", "--output", out("caller-below-range", "made"),
                "--max-copy-number", "-1"));
    }

    static String[] concat(final String[] first, final String... second) {
        final List<String> all = new ArrayList<>(Arrays.asList(first));
        all.addAll(Arrays.asList(second));
        return all.toArray(new String[0]);
    }

    static String[] postprocess(final String label, final String calls, final String model) {
        return new String[] {
                "--calls-shard-path", at(calls), "--model-shard-path", at(model),
                "--contig-ploidy-calls", at("nothing"), "--sample-index", "0",
                "--output-genotyped-intervals", out(label, "intervals.vcf.gz"),
                "--output-genotyped-segments", out(label, "segments.vcf.gz"),
                "--output-denoised-copy-ratios", out(label, "ratios.tsv"),
        };
    }

    static void postprocessRows() throws Exception {
        final String tool = "PostprocessGermlineCNVCalls";
        run("post-happy", tool, postprocess("post-happy", "nothing", "nothing"));
        run("post-copy-number", tool, concat(postprocess("post-copy-number", "nothing", "nothing"), "--autosomal-ref-copy-number", "1"));
        run("post-negative-index", tool, concat(postprocess("post-negative-index", "nothing", "nothing"), "--sample-index", "-1"));
        run("post-negative-copy-number", tool, concat(postprocess("post-negative-copy-number", "nothing", "nothing"), "--autosomal-ref-copy-number", "-1"));
        run("post-missing-calls", tool, "--model-shard-path", at("nothing"), "--contig-ploidy-calls", at("nothing"), "--sample-index", "0",
                "--output-genotyped-intervals", out("post-missing-calls", "a"), "--output-genotyped-segments", out("post-missing-calls", "b"),
                "--output-denoised-copy-ratios", out("post-missing-calls", "c"));
        // GATKTool's startup comes before the import.
        run("post-no-reference", tool, concat(postprocess("post-no-reference", "nothing", "nothing"), "-R", at("nothing.fasta")));
        run("post-no-reads", tool, concat(postprocess("post-no-reads", "nothing", "nothing"), "-I", at("nothing.bam")));
        run("post-no-master", tool, concat(postprocess("post-no-master", "nothing", "nothing"), "--sequence-dictionary", at("nothing.dict")));
        run("post-filtering-without-intervals", tool, concat(postprocess("post-filtering-without-intervals", "nothing", "nothing"),
                "--variant-output-filtering", "CONTAINED"));
        run("post-shards", tool, postprocess("post-shards", "calls-shard", "model-shard"));
        run("post-shards-undisturbed", tool, concat(postprocess("post-shards-undisturbed", "nothing", "nothing"), "--disable-sequence-dictionary-validation"));
        run("post-shards-other-dictionary", tool, postprocess("post-shards-other-dictionary", "calls-shard", "model-other-dictionary"));
        run("post-intervals-no-shard", tool, concat(postprocess("post-intervals-no-shard", "nothing", "nothing"), "-L", "chr1"));
        run("post-intervals-empty-shard", tool, concat(postprocess("post-intervals-empty-shard", "empty-folder", "empty-folder"), "-L", "chr1"));
        run("post-intervals", tool, concat(postprocess("post-intervals", "calls-shard", "model-shard"), "-L", "chr1"));
        run("post-intervals-unknown-contig", tool, concat(postprocess("post-intervals-unknown-contig", "calls-shard", "model-shard"), "-L", "chrOther"));
        run("post-intervals-other-dictionary", tool, concat(postprocess("post-intervals-other-dictionary", "calls-shard", "model-other-dictionary"), "-L", "chr1"));
        run("post-intervals-with-master", tool, concat(postprocess("post-intervals-with-master", "nothing", "nothing"),
                "--sequence-dictionary", at("ref.dict"), "-L", "chr1"));
    }

    static void nvScoreRows() throws Exception {
        final String tool = "NVScoreVariants";
        final String[] required = {"-V", at("nothing.vcf"), "-R", at("nothing.fasta")};
        run("nv-happy", tool, concat(required, "-O", out("nv-happy", "scored.vcf")));
        // The two refusals `doWork` holds never come.
        run("nv-reads-without-tensor", tool, concat(required, "-O", out("nv-reads-without-tensor", "scored.vcf"), "-I", at("nothing.bam")));
        run("nv-tensor-without-reads", tool, concat(required, "-O", out("nv-tensor-without-reads", "scored.vcf"), "--tensor-type", "read_tensor"));
        run("nv-bad-tensor", tool, concat(required, "-O", out("nv-bad-tensor", "scored.vcf"), "--tensor-type", "other"));
        run("nv-bad-batch", tool, concat(required, "-O", out("nv-bad-batch", "scored.vcf"), "--batch-size", "many"));
        run("nv-missing-output", tool, required);
        run("nv-missing-reference", tool, "-V", at("nothing.vcf"), "-O", out("nv-missing-reference", "scored.vcf"));
    }

    static void train(final String label, final String annotations, final String... rest) throws Exception {
        final String name = "train-" + label;
        run(name, "TrainVariantAnnotationsModel",
                concat(new String[] {"--annotations-hdf5", annotations, "-O", out(name, "model")}, rest));
    }

    static void trainRows() throws Exception {
        final String labelled = out("setup-labelled", "extract.annot.hdf5");
        final String unlabeled = out("setup-labelled", "extract.unlabeled.annot.hdf5");
        final String json = at("hyperparameters.json");
        final String script = at("score.py");
        final String[] byScript = {"--model-backend", "PYTHON_SCRIPT", "--python-script", script, "--hyperparameters-json", json};

        train("forest", labelled);
        // The default backend reads nothing of the file before it imports.
        train("forest-not-hdf5", at("not-hdf5.hdf5"));
        train("forest-no-file", at("nothing.hdf5"));
        train("forest-folder", at("empty-folder"));
        train("forest-script", labelled, "--python-script", script);
        train("forest-json", labelled, "--hyperparameters-json", json);
        train("forest-no-json", labelled, "--hyperparameters-json", at("nothing.json"));
        train("forest-unlabeled", labelled, "--unlabeled-annotations-hdf5", unlabeled);
        train("forest-unlabeled-no-file", labelled, "--unlabeled-annotations-hdf5", at("nothing.hdf5"));
        train("forest-unlabeled-not-hdf5", labelled, "--unlabeled-annotations-hdf5", at("not-hdf5.hdf5"));
        train("forest-other-names", labelled, "--unlabeled-annotations-hdf5", out("setup-other-names", "extract.annot.hdf5"));
        train("forest-first-not-hdf5", at("not-hdf5.hdf5"), "--unlabeled-annotations-hdf5", unlabeled);
        train("bad-backend", labelled, "--model-backend", "OTHER");
        train("bad-mode", labelled, "--mode", "OTHER");

        train("bgmm-script", labelled, "--model-backend", "JAVA_BGMM", "--python-script", script);
        train("bgmm-no-json", labelled, "--model-backend", "JAVA_BGMM");
        train("bgmm-missing-json", labelled, "--model-backend", "JAVA_BGMM", "--hyperparameters-json", at("nothing.json"));
        train("bgmm", labelled, "--model-backend", "JAVA_BGMM", "--hyperparameters-json", json);
        train("bgmm-indel", labelled, "--model-backend", "JAVA_BGMM", "--hyperparameters-json", json, "--mode", "INDEL");
        train("bgmm-no-training-label", out("setup-calibration-only", "extract.annot.hdf5"), "--model-backend", "JAVA_BGMM", "--hyperparameters-json", json);
        train("bgmm-no-indels", out("setup-only-snps", "extract.annot.hdf5"), "--model-backend", "JAVA_BGMM", "--hyperparameters-json", json, "--mode", "INDEL");
        train("bgmm-missing-values", out("setup-missing-values", "extract.annot.hdf5"), "--model-backend", "JAVA_BGMM", "--hyperparameters-json", json);

        train("script-no-json", labelled, "--model-backend", "PYTHON_SCRIPT", "--python-script", script);
        train("script-no-script", labelled, "--model-backend", "PYTHON_SCRIPT", "--hyperparameters-json", json);
        train("script-missing-script", labelled, "--model-backend", "PYTHON_SCRIPT", "--python-script", at("nothing.py"), "--hyperparameters-json", json);
        train("script-missing-json", labelled, "--model-backend", "PYTHON_SCRIPT", "--python-script", script, "--hyperparameters-json", at("nothing.json"));
        train("script", labelled, byScript);
        train("script-snp", labelled, concat(byScript, "--mode", "SNP"));
        train("script-indel", labelled, concat(byScript, "--mode", "INDEL"));
        train("script-not-hdf5", at("not-hdf5.hdf5"), byScript);
        train("script-no-calibration-label", out("setup-training-only", "extract.annot.hdf5"), byScript);
        train("script-no-training-label", out("setup-calibration-only", "extract.annot.hdf5"), byScript);
        train("script-no-indels", out("setup-only-snps", "extract.annot.hdf5"), concat(byScript, "--mode", "INDEL"));
        train("script-missing-values", out("setup-missing-values", "extract.annot.hdf5"), byScript);
        train("script-unlabeled", labelled, concat(byScript, "--unlabeled-annotations-hdf5", unlabeled));
        train("script-no-unlabeled-snps", out("setup-every-snp", "extract.annot.hdf5"),
                concat(byScript, "--unlabeled-annotations-hdf5", out("setup-every-snp", "extract.unlabeled.annot.hdf5")));
        train("script-unlabeled-values", labelled,
                concat(byScript, "--unlabeled-annotations-hdf5", out("setup-missing-values", "extract.annot.hdf5")));
    }

    static void score(final String label, final String prefix, final String... rest) throws Exception {
        final String name = "score-" + label;
        run(name, "ScoreVariantAnnotations", concat(new String[] {
                "-V", at("calls.vcf"), "-A", "QD", "-A", "MQ",
                "--model-prefix", at(prefix), "-O", out(name, "scored"),
                "--add-output-vcf-command-line", "false", "--use-jdk-deflater"}, rest));
    }

    static void scoreRows() throws Exception {
        final String script = at("score.py");
        final String[] byScript = {"--model-backend", "PYTHON_SCRIPT", "--python-script", script};
        final String resource = "--resource:train,training=true,calibration=true";

        // PYTHON_IFOREST, the default: the header, then the import.
        score("forest", "both");
        score("forest-plain", "both", "--do-not-gzip-vcf-output");
        score("forest-plain-no-index", "both", "--do-not-gzip-vcf-output", "--create-output-variant-index", "false");
        score("forest-no-models", "nothing");
        score("forest-script", "both", "--python-script", script);
        score("forest-resource", "both", resource, at("resource.vcf"), "--do-not-gzip-vcf-output");
        score("forest-reserved", "both", "--resource:typed,snp=true", at("resource.vcf"));
        score("forest-keys", "both", "--do-not-gzip-vcf-output", "--score-key", "S", "--calibration-sensitivity-key", "C",
                "--low-score-filter-name", "L", "--snp-key", "null");
        score("forest-null-score-key", "both", "--do-not-gzip-vcf-output", "--score-key", "null");
        score("forest-null-calibration-key", "both", "--do-not-gzip-vcf-output", "--calibration-sensitivity-key", "null");
        score("forest-null-filter", "both", "--do-not-gzip-vcf-output", "--low-score-filter-name", "null");
        score("forest-null-format", "both", "--do-not-gzip-vcf-output", "--double-format", "null");
        score("forest-allele-specific", "both", "--do-not-gzip-vcf-output", "-A", "AS_QD");
        score("forest-undeclared", "both", "-A", "NOT_DECLARED");
        run("score-missing-folder", "ScoreVariantAnnotations", "-V", at("calls.vcf"), "-A", "QD",
                "--model-prefix", at("both"), "-O", out("score-missing-folder", "no/such/folder/scored"),
                "--add-output-vcf-command-line", "false");
        run("score-missing-prefix", "ScoreVariantAnnotations", "-V", at("calls.vcf"), "-A", "QD",
                "-O", out("score-missing-prefix", "scored"));
        run("score-no-annotation", "ScoreVariantAnnotations", "-V", at("calls.vcf"), "--model-prefix", at("both"),
                "-O", out("score-no-annotation", "scored"));
        run("score-no-variants-file", "ScoreVariantAnnotations", "-V", at("nothing.vcf"), "-A", "QD",
                "--model-prefix", at("both"), "-O", out("score-no-variants-file", "scored"));
        score("threshold-high", "both", "--snp-calibration-sensitivity-threshold", "1.5");
        score("bad-backend", "both", "--model-backend", "OTHER");

        // JAVA_BGMM: no scorer is a refusal and a scorer is the reference's own missing module.
        score("bgmm-script", "ser", "--model-backend", "JAVA_BGMM", "--python-script", script);
        score("bgmm-no-models", "nothing", "--model-backend", "JAVA_BGMM");
        score("bgmm-scorer", "ser", "--model-backend", "JAVA_BGMM");
        score("bgmm-real-scorer", "bgmm", "--model-backend", "JAVA_BGMM");

        // PYTHON_SCRIPT: the scorers are only asked whether they can be read.
        score("script-no-script", "both", "--model-backend", "PYTHON_SCRIPT");
        score("script-missing-script", "both", "--model-backend", "PYTHON_SCRIPT", "--python-script", at("nothing.py"));
        score("script-no-models", "nothing", byScript);
        score("script-only-snps", "snp-only", byScript);
        score("script-only-indels", "indel-only", byScript);
        score("script-threshold-snp", "both", concat(byScript, "--snp-calibration-sensitivity-threshold", "0.5"));
        score("script-threshold-indel", "both", concat(byScript, "--indel-calibration-sensitivity-threshold", "0.5"));
        score("script-broken-calibration", "broken-calibration", byScript);
        score("script-infinite-calibration", "infinite", byScript);
        score("script-calibration", "calibrated", concat(byScript, "--snp-calibration-sensitivity-threshold", "0.5", "--do-not-gzip-vcf-output"));
        score("script", "both", concat(byScript, "--do-not-gzip-vcf-output"));
        score("script-gz", "both", byScript);
        score("script-omit-alleles", "both", concat(byScript, "--omit-alleles-in-hdf5", "--do-not-gzip-vcf-output"));
        score("script-resource", "both", concat(byScript, resource, at("resource.vcf"), "--do-not-gzip-vcf-output"));
        score("script-resource-unindexed", "both", concat(byScript, resource, at("resource-unindexed.vcf"), "--do-not-gzip-vcf-output"));
        score("script-only-snps-mode", "snp-only", concat(byScript, "--mode", "SNP", "--do-not-gzip-vcf-output"));
        score("script-allele-specific", "both", concat(byScript, "-A", "AS_QD", "--do-not-gzip-vcf-output"));
        score("script-ignore-filter", "both", concat(byScript, "--ignore-filter", "LowQual", "--do-not-gzip-vcf-output"));
        score("script-ignore-all", "both", concat(byScript, "--ignore-all-filters", "--do-not-gzip-vcf-output"));

        // A pass that extracts nothing: no script runs, and the VCF is the input's records.
        score("script-nothing", "both", concat(byScript, "-L", "chr2:1-50", "--do-not-gzip-vcf-output", "--create-output-variant-index", "false"));
        score("script-nothing-gz", "both", concat(byScript, "-L", "chr2:1-50"));
        score("script-nothing-indexed", "both", concat(byScript, "-L", "chr2:1-50", "--do-not-gzip-vcf-output"));
        score("script-only-filtered", "both", concat(byScript, "-L", "chr1:300-300", "--do-not-gzip-vcf-output", "--create-output-variant-index", "false"));
        score("script-other-type", "both", concat(byScript, "--mode", "INDEL", "-L", "chr1:100-200", "--do-not-gzip-vcf-output", "--create-output-variant-index", "false"));
    }
}

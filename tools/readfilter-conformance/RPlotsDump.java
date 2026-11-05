/*
 * The three tools that hand their work to an R script, with R there to run it.
 *
 * This is the one dump that does not run in gatk-rs-oracle: that image has no Rscript, on purpose,
 * so `tools/conformance/manifest.json` sends it to `oracle_with_r`, the Broad's own image of the
 * same release, where the reference's R and its packages are the ones the release was built with.
 *
 * What it holds is everything around the image, because the PNG bytes are not comparable: cairo's
 * rasterisation, the font metrics and the ggplot geometry all differ, and the rewrite is on ggplot2
 * besides. No row here says anything about a pixel.
 *
 *   - THE SIDE FILES, BYTE FOR BYTE: PlotDenoisedCopyRatios' four MADs, each one number rounded to
 *     three places, where the MAD is `median(abs(diff(copy_ratio)))`;
 *   - A ONE-POINT FILE WRITES NA and a flat one divides nought by nought into NaN, and neither
 *     fails the run;
 *   - WHAT EACH RUN LEAVES IN THE WORKING DIRECTORY, since an R error there writes
 *     `plotting_dump.rda` beside whatever the output directory got;
 *   - THE STATUS OF A SCRIPT RUN ON ITS OWN, whose bare `quit(status=1)` on a missing input the
 *     Java side's own refusal would otherwise always get ahead of;
 *   - THE NUMBERS THE SCRIPTS PLOT, which they keep nowhere: RPlotsTrace.R masks the graphics calls
 *     ahead of the reference's own unmodified scripts and records their arguments. That is where
 *     `maximum_copy_ratio = Infinity` giving each panel its own limit shows, and the running index
 *     of the modeled segments: consecutive segments SHARING A POINT, a zero-point segment skipped
 *     without advancing, the last segment RUNNING OFF THE TABLE into an NA, and a point shared
 *     across a contig boundary taking the later contig's offset;
 *   - THE INNER JOIN OF plotAlleleFrequencyQC.R, on tables a VariantEval report never gives it, which
 *     is why those rows run the script directly;
 *   - AND THE PNG's PIXEL SIZE AND RESOLUTION, read from its header, which the rewrite keeps.
 *
 * Output:
 *
 *     input\t<name>\t<the file, escaped>
 *     argv\t<label>\t<the command line, tab-separated, escaped>
 *     stdout\t<label>\t<what the tool printed, escaped, cwd masked>
 *     error\t<label>\t<exception class>:<the message's first line, escaped, cwd masked>
 *     produced\t<label>\t<path under the run's output, or cwd/ and a path>
 *     side\t<label>\t<path under the run's output>\t<the file, escaped>
 *     png\t<label>\t<path under the run's output>\t<width>x<height>\t<pixels per metre>
 *     metrics\t<label>\t<the metrics table without its comments, escaped>
 *     script\t<label>\t<denoised, modeled or qc>\t<the script's arguments, tab-separated, escaped>
 *     script-status\t<label>\t<the status of that run>
 *     trace\t<label>\t<one row RPlotsTrace.R recorded, escaped>
 *
 * A `trace` row under an `argv` label comes from a run of the scripts on the arguments Java passes,
 * since a trace needs the mask in front and Java's executor cannot put it there.
 *
 * Usage: RPlotsDump
 */

import org.broadinstitute.hellbender.Main;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.InputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class RPlotsDump {

    static final Path DIR = Path.of("r-plots-dump");
    static final Path RESOURCES = DIR.resolve("reference-scripts");
    static final String TRACE = "/harness/RPlotsTrace.R";
    static final String PLOTTING = "/org/broadinstitute/hellbender/tools/copynumber/plotting/";
    static final String VARIANTEVAL = "/org/broadinstitute/hellbender/tools/walkers/varianteval/";

    static final String SEQUENCES = "@SQ\tSN:chr1\tLN:2000000\n"
            + "@SQ\tSN:chr2\tLN:1500000\n"
            + "@SQ\tSN:chrShort\tLN:1000\n";
    static final String DICT = "@HD\tVN:1.6\n" + SEQUENCES;
    static final String HEADER = DICT + "@RG\tID:GATKCopyNumber\tSM:sampleA\n";
    static final String RATIO_COLUMNS = "CONTIG\tSTART\tEND\tLOG2_COPY_RATIO\n";
    static final String COUNT_COLUMNS = "CONTIG\tPOSITION\tREF_COUNT\tALT_COUNT\tREF_NUCLEOTIDE\tALT_NUCLEOTIDE\n";
    static final String SEGMENT_COLUMNS = "CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tNUM_POINTS_ALLELE_FRACTION\t"
            + "LOG2_COPY_RATIO_POSTERIOR_10\tLOG2_COPY_RATIO_POSTERIOR_50\tLOG2_COPY_RATIO_POSTERIOR_90\t"
            + "MINOR_ALLELE_FRACTION_POSTERIOR_10\tMINOR_ALLELE_FRACTION_POSTERIOR_50\tMINOR_ALLELE_FRACTION_POSTERIOR_90\n";

    /** Six points on chr1 and four on chr2, with a midpoint that rounds to even and a large ratio. */
    static final String RATIOS = RATIO_COLUMNS
            + "chr1\t1000\t1999\t0.5\n"
            + "chr1\t2000\t2999\t-0.25\n"
            + "chr1\t3000\t3999\t0.125\n"
            + "chr1\t5001\t5002\t-1.5\n"
            + "chr1\t6000\t6999\t3.25\n"
            + "chr1\t7000\t7999\t0.0625\n"
            + "chr2\t500\t999\t1.0\n"
            + "chr2\t1000\t1499\t-0.75\n"
            + "chr2\t1500\t1999\t0.3333333333333333\n"
            + "chr2\t2000\t2499\t-3.5\n";
    static final String DENOISED_RATIOS = RATIO_COLUMNS
            + "chr1\t1000\t1999\t0.25\n"
            + "chr1\t2000\t2999\t-0.125\n"
            + "chr1\t3000\t3999\t0.0625\n"
            + "chr1\t5001\t5002\t-0.75\n"
            + "chr1\t6000\t6999\t1.625\n"
            + "chr1\t7000\t7999\t0.03125\n"
            + "chr2\t500\t999\t0.5\n"
            + "chr2\t1000\t1499\t-0.375\n"
            + "chr2\t1500\t1999\t0.16666666666666666\n"
            + "chr2\t2000\t2499\t-1.75\n";
    static final String COUNTS = COUNT_COLUMNS
            + "chr1\t1500\t10\t12\tA\tC\n"
            + "chr1\t2500\t8\t9\tG\tT\n"
            + "chr1\t3500\t20\t3\tC\tA\n"
            + "chr1\t5500\t0\t0\tC\tA\n"
            + "chr2\t700\t4\t3\tT\tA\n"
            + "chr2\t1200\t7\t7\tT\tG\n";
    static final String POSTERIORS = "\t-0.1\t0.2\t0.3\t0.1\t0.25\t0.4\n";
    static final String OTHER_POSTERIORS = "\t-1.5\t-1.25\t-1.0\t0.3\t0.35\t0.5\n";

    /** Three segments, so that two boundaries are shared and the last one runs off the table. */
    static final String SEGMENTS = SEGMENT_COLUMNS
            + "chr1\t1000\t3999\t3\t3" + POSTERIORS
            + "chr1\t5001\t7999\t3\t1" + OTHER_POSTERIORS
            + "chr2\t500\t2499\t4\t2" + POSTERIORS;
    /** A segment with no copy-ratio points in the middle, and one with no allele fractions: no both, which Java refuses. */
    static final String SEGMENTS_ZERO_POINTS = SEGMENT_COLUMNS
            + "chr1\t1000\t3999\t3\t3" + POSTERIORS
            + "chr1\t4000\t4999\t0\t1" + OTHER_POSTERIORS
            + "chr1\t5001\t7999\t3\t0" + OTHER_POSTERIORS
            + "chr2\t500\t2499\t4\t2" + POSTERIORS;
    static final String SEGMENTS_NO_POINTS = SEGMENT_COLUMNS
            + "chr1\t1000\t3999\t0\t0" + POSTERIORS;

    static final String QC_COLUMNS = "VariantAFEvaluator  CompFeatureInput  EvalFeatureInput  AlleleFrequency  Filter  avgVarAF\n";

    public static void main(final String[] args) throws Exception {
        PrintReadsDump.emptyDirectory(DIR);
        Files.createDirectories(DIR.resolve("out"));
        Files.createDirectories(RESOURCES);

        System.out.println("# RPlotsDump: PlotDenoisedCopyRatios, PlotModeledSegments and AlleleFrequencyQC with R");

        extract(PLOTTING, "CNVPlottingLibrary.R");
        extract(PLOTTING, "PlotDenoisedCopyRatios.R");
        extract(PLOTTING, "PlotModeledSegments.R");
        extract(VARIANTEVAL, "plotAlleleFrequencyQC.R");

        input("ref.dict", DICT);
        input("only-long.dict", "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:2000000\n");
        input("standardized.tsv", HEADER + RATIOS);
        input("denoised.tsv", HEADER + DENOISED_RATIOS);
        input("one-standardized.tsv", HEADER + RATIO_COLUMNS + "chr1\t1000\t1999\t0.5\n");
        input("one-denoised.tsv", HEADER + RATIO_COLUMNS + "chr1\t1000\t1999\t0.25\n");
        input("flat.tsv", HEADER + RATIO_COLUMNS
                + "chr1\t1000\t1999\t0.5\nchr1\t2000\t2999\t0.5\nchr1\t3000\t3999\t0.5\n");
        input("no-points.tsv", HEADER + RATIO_COLUMNS);
        input("allelic.tsv", HEADER + COUNTS);
        input("segments.seg", HEADER + SEGMENTS);
        input("segments-zero-points.seg", HEADER + SEGMENTS_ZERO_POINTS);
        input("segments-no-points.seg", HEADER + SEGMENTS_NO_POINTS);
        input("segments-none.seg", HEADER + SEGMENT_COLUMNS);

        denoisedRows();
        modeledRows();
        qcRows();
    }

    static void denoisedRows() throws Exception {
        denoised("default", "standardized.tsv", "denoised.tsv", "ref.dict");
        denoised("infinity", "standardized.tsv", "denoised.tsv", "ref.dict", "--maximum-copy-ratio", "Infinity");
        denoised("low-maximum", "standardized.tsv", "denoised.tsv", "ref.dict", "--maximum-copy-ratio", "1.5",
                "--point-size-copy-ratio", "0.5");
        denoised("short-contig-kept", "standardized.tsv", "denoised.tsv", "ref.dict", "--minimum-contig-length", "0");
        denoised("contig-dropped", "standardized.tsv", "denoised.tsv", "only-long.dict");
        denoised("one-point", "one-standardized.tsv", "one-denoised.tsv", "ref.dict");
        denoised("flat", "flat.tsv", "flat.tsv", "ref.dict");
        denoised("no-points", "no-points.tsv", "no-points.tsv", "ref.dict");
        denoised("no-points-infinity", "no-points.tsv", "no-points.tsv", "ref.dict", "--maximum-copy-ratio", "Infinity");
        denoised("no-such-file", "nothing.tsv", "denoised.tsv", "ref.dict");
        script("script-denoised-no-such-file", "denoised", "--sample_name=sampleA",
                "--standardized_copy_ratios_file=" + DIR.resolve("nothing.tsv"),
                "--denoised_copy_ratios_file=" + DIR.resolve("denoised.tsv"),
                "--contig_names=chr1CONTIG_DELIMITERchr2", "--contig_lengths=2000000CONTIG_DELIMITER1500000",
                "--maximum_copy_ratio=4.0", "--point_size_copy_ratio=0.2",
                "--output_dir=" + DIR.resolve("out/script-denoised-no-such-file") + "/", "--output_prefix=sample");
    }

    static void modeledRows() throws Exception {
        final String denoised = "denoised.tsv";
        final String allelic = "allelic.tsv";
        modeled("both", denoised, allelic, "segments.seg", "ref.dict");
        modeled("ratios-only", denoised, null, "segments.seg", "ref.dict");
        modeled("counts-only", null, allelic, "segments.seg", "ref.dict");
        modeled("infinity", denoised, allelic, "segments.seg", "ref.dict", "--maximum-copy-ratio", "Infinity");
        modeled("zero-points", denoised, allelic, "segments-zero-points.seg", "ref.dict");
        modeled("contig-dropped", denoised, allelic, "segments.seg", "only-long.dict");
        modeled("no-segments", "no-points.tsv", null, "segments-none.seg", "ref.dict");
        modeled("all-zero", "no-points.tsv", null, "segments-no-points.seg", "ref.dict");
        modeled("all-zero-infinity", "no-points.tsv", null, "segments-no-points.seg", "ref.dict",
                "--maximum-copy-ratio", "Infinity");
        modeled("neither", null, null, "segments.seg", "ref.dict");
        modeled("no-such-segments", denoised, allelic, "nothing.seg", "ref.dict");
        final List<String> common = List.of(
                "--contig_names=chr1CONTIG_DELIMITERchr2", "--contig_lengths=2000000CONTIG_DELIMITER1500000",
                "--maximum_copy_ratio=4.0", "--point_size_copy_ratio=0.2", "--point_size_allele_fraction=0.4");
        script("script-modeled-no-such-file", "modeled", concat(List.of("--sample_name=sampleA",
                "--denoised_copy_ratios_file=" + DIR.resolve(denoised), "--allelic_counts_file=null",
                "--modeled_segments_file=" + DIR.resolve("nothing.seg")), common,
                List.of("--output_file=" + DIR.resolve("out/script-modeled-no-such-file/sample.modeled.png"))));
        script("script-modeled-no-inputs", "modeled", concat(List.of("--sample_name=sampleA",
                "--denoised_copy_ratios_file=null", "--allelic_counts_file=null",
                "--modeled_segments_file=" + DIR.resolve("segments.seg")), common,
                List.of("--output_file=" + DIR.resolve("out/script-modeled-no-inputs/sample.modeled.png"))));
    }

    static void qcRows() throws Exception {
        final Path fasta = VariantEvalDump.writeReference(DIR);
        input("eval.vcf", AlleleFrequencyQCDump.evalVcf(List.of(
                AlleleFrequencyQCDump.called(1000, "0/1"), AlleleFrequencyQCDump.called(2000, "0/1"),
                AlleleFrequencyQCDump.called(3000, "0/1"), AlleleFrequencyQCDump.called(4000, "1/1")), true));
        input("comp.vcf", AlleleFrequencyQCDump.compVcf(List.of(
                AlleleFrequencyQCDump.site(1000, "0.005"), AlleleFrequencyQCDump.site(2000, "0.05"),
                AlleleFrequencyQCDump.site(3000, "0.25"), AlleleFrequencyQCDump.site(4000, "0.5"))));
        for (final String name : List.of(fasta.getFileName().toString(), fasta.getFileName() + ".fai",
                fasta.getFileName().toString().replaceAll("\\.[^.]*$", "") + ".dict")) {
            System.out.printf("input\t%s\t%s%n", name, ReferenceQueryDump.escape(Files.readString(DIR.resolve(name))));
        }
        for (final String name : List.of("eval.vcf", "comp.vcf")) {
            new Main().instanceMain(new String[]{"IndexFeatureFile", "-I", at(name)});
        }
        // The comparison set twice, once as --comp and once as a second eval track: see
        // AlleleFrequencyQCDump. The report is kept, so that the traced run can read it.
        for (final String[] variant : new String[][]{{"qc", "0.05"}, {"qc-threshold", "1.1"}}) {
            final String label = variant[0];
            final String metrics = out(label, "qc.metrics");
            final String report = out(label, "report.txt");
            run(label, "AlleleFrequencyQC", "-O", metrics, "-R", fasta.toString(),
                    "--eval", at("eval.vcf"), "--comp", at("comp.vcf"),
                    "-eval:thousand_genomes", at("comp.vcf"), "-L", at("comp.vcf"),
                    "-debug-file", report, "-pvalue-threshold", variant[1]);
            final Path traceOut = DIR.resolve("traces").resolve(label);
            Files.createDirectories(traceOut);
            traced(label, "qc", report, traceOut.resolve("qc.metrics").toString(), "NA12878");
        }

        input("qc-same-bins.txt", QC_COLUMNS
                + "VariantAFEvaluator  none  eval  1.0  called  0.25\n"
                + "VariantAFEvaluator  none  eval  3.0  called  0.5\n"
                + "VariantAFEvaluator  none  eval  3.0  filtered  0.75\n"
                + "VariantAFEvaluator  none  thousand_genomes  1.0  called  0.2\n"
                + "VariantAFEvaluator  none  thousand_genomes  3.0  called  0.4\n");
        input("qc-other-bins.txt", QC_COLUMNS
                + "VariantAFEvaluator  none  eval  1.0  called  0.25\n"
                + "VariantAFEvaluator  none  eval  3.0  called  0.5\n"
                + "VariantAFEvaluator  none  thousand_genomes  3.0  called  0.4\n"
                + "VariantAFEvaluator  none  thousand_genomes  7.0  called  0.1\n");
        input("qc-no-called.txt", QC_COLUMNS
                + "VariantAFEvaluator  none  eval  1.0  filtered  0.25\n"
                + "VariantAFEvaluator  none  thousand_genomes  1.0  filtered  0.2\n");
        for (final String table : new String[]{"qc-same-bins", "qc-other-bins", "qc-no-called", "qc-no-such-file"}) {
            final String name = "script-" + table;
            script(name, "qc", at(table + ".txt"), out(name, "qc.metrics"), "NA12878");
        }
    }

    static void denoised(final String label, final String standardized, final String denoised,
                         final String dictionary, final String... extra) throws Exception {
        final String name = "denoised-" + label;
        final List<String> rest = new ArrayList<>(List.of(
                "--standardized-copy-ratios", at(standardized),
                "--denoised-copy-ratios", at(denoised),
                "--sequence-dictionary", at(dictionary),
                "--output-prefix", "sample",
                "--output", out(name, "plots")));
        rest.addAll(Arrays.asList(extra));
        run(name, "PlotDenoisedCopyRatios", rest.toArray(new String[0]));
        if (!Files.exists(DIR.resolve(standardized))) {
            return;
        }
        // The arguments Java passes, in a directory of their own so that the run above stays whole.
        final Path traceOut = DIR.resolve("traces").resolve(name);
        Files.createDirectories(traceOut);
        traced(name, "denoised", "--sample_name=sampleA",
                "--standardized_copy_ratios_file=" + DIR.resolve(standardized),
                "--denoised_copy_ratios_file=" + DIR.resolve(denoised),
                "--contig_names=" + contigs(dictionary, extra), "--contig_lengths=" + lengths(dictionary, extra),
                "--maximum_copy_ratio=" + option(extra, "--maximum-copy-ratio", "4.0"),
                "--point_size_copy_ratio=" + option(extra, "--point-size-copy-ratio", "0.2"),
                "--output_dir=" + traceOut + "/", "--output_prefix=sample");
    }

    static void modeled(final String label, final String denoised, final String allelic, final String segments,
                        final String dictionary, final String... extra) throws Exception {
        final String name = "modeled-" + label;
        final List<String> rest = new ArrayList<>(List.of(
                "--segments", at(segments),
                "--sequence-dictionary", at(dictionary),
                "--output-prefix", "sample",
                "--output", out(name, "plots")));
        if (denoised != null) {
            rest.addAll(List.of("--denoised-copy-ratios", at(denoised)));
        }
        if (allelic != null) {
            rest.addAll(List.of("--allelic-counts", at(allelic)));
        }
        rest.addAll(Arrays.asList(extra));
        run(name, "PlotModeledSegments", rest.toArray(new String[0]));
        if (!Files.exists(DIR.resolve(segments)) || (denoised == null && allelic == null)) {
            return;
        }
        final Path traceOut = DIR.resolve("traces").resolve(name);
        Files.createDirectories(traceOut);
        traced(name, "modeled", "--sample_name=sampleA",
                "--denoised_copy_ratios_file=" + (denoised == null ? "null" : DIR.resolve(denoised)),
                "--allelic_counts_file=" + (allelic == null ? "null" : DIR.resolve(allelic)),
                "--modeled_segments_file=" + DIR.resolve(segments),
                "--contig_names=" + contigs(dictionary, extra), "--contig_lengths=" + lengths(dictionary, extra),
                "--maximum_copy_ratio=" + option(extra, "--maximum-copy-ratio", "4.0"),
                "--point_size_copy_ratio=0.2", "--point_size_allele_fraction=0.4",
                "--output_file=" + traceOut.resolve("sample.modeled.png"));
    }

    /** `PlottingUtils.getContigLengthMap` over this dump's dictionaries: chrShort only below the default. */
    static String contigs(final String dictionary, final String[] extra) {
        return !dictionary.equals("ref.dict") ? "chr1"
                : "chr1CONTIG_DELIMITERchr2" + (option(extra, "--minimum-contig-length", "").equals("0") ? "CONTIG_DELIMITERchrShort" : "");
    }

    static String lengths(final String dictionary, final String[] extra) {
        return !dictionary.equals("ref.dict") ? "2000000"
                : "2000000CONTIG_DELIMITER1500000" + (option(extra, "--minimum-contig-length", "").equals("0") ? "CONTIG_DELIMITER1000" : "");
    }

    static String option(final String[] extra, final String name, final String fallback) {
        final int at = Arrays.asList(extra).indexOf(name);
        return at < 0 ? fallback : extra[at + 1];
    }

    static String[] concat(final List<String> a, final List<String> b, final List<String> c) {
        final List<String> all = new ArrayList<>(a);
        all.addAll(b);
        all.addAll(c);
        return all.toArray(new String[0]);
    }

    static void extract(final String path, final String name) throws Exception {
        try (InputStream in = RPlotsDump.class.getResourceAsStream(path + name)) {
            Files.copy(in, RESOURCES.resolve(name));
        }
    }

    static void input(final String name, final String text) throws Exception {
        Files.writeString(DIR.resolve(name), text);
        System.out.printf("input\t%s\t%s%n", name, ReferenceQueryDump.escape(text));
    }

    static String at(final String name) {
        return DIR.resolve(name).toString();
    }

    static String out(final String label, final String name) throws Exception {
        final Path directory = DIR.resolve("out").resolve(label);
        Files.createDirectories(directory);
        return directory.resolve(name).toString();
    }

    static String masked(final String text) {
        return text.replace(new File("").getAbsolutePath(), "<cwd>");
    }

    /** The message's first line with anything in it: the rest of a script's is its temporary paths. */
    static String firstLine(final String message) {
        if (message == null) {
            return "null";
        }
        return Arrays.stream(message.split("\n")).filter(line -> !line.isBlank()).findFirst().orElse("");
    }

    /** Everything in the working directory a run added: what an R error leaves is `cwd/`. */
    static Set<String> workingDirectory() {
        return Arrays.stream(new File(".").list())
                .filter(name -> !name.endsWith(".java") && !name.endsWith(".class") && !name.equals(DIR.toString()))
                .collect(Collectors.toSet());
    }

    static boolean run(final String label, final String tool, final String... rest) throws Exception {
        final List<String> argv = new ArrayList<>();
        argv.add(tool);
        argv.addAll(Arrays.asList(rest));
        System.out.printf("argv\t%s\t%s%n", label, ReferenceQueryDump.escape(String.join("\t", argv)));
        final Path outDir = DIR.resolve("out").resolve(label);
        Files.createDirectories(outDir);
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
                failure = e.getClass().getName() + ":" + firstLine(cause.getMessage());
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
        listOutput(label, outDir);
        return failure == null;
    }

    static void listOutput(final String label, final Path outDir) throws Exception {
        final List<Path> produced;
        try (Stream<Path> walk = Files.walk(outDir)) {
            produced = walk.filter(Files::isRegularFile).sorted().collect(Collectors.toList());
        }
        for (final Path path : produced) {
            final Path relative = outDir.relativize(path);
            System.out.printf("produced\t%s\t%s%n", label, relative);
            final String name = relative.toString();
            if (name.endsWith(".png")) {
                final byte[] bytes = Files.readAllBytes(path);
                System.out.printf("png\t%s\t%s\t%s%n", label, relative, pngShape(bytes));
            } else if (name.endsWith(".txt") && name.contains("MAD")) {
                System.out.printf("side\t%s\t%s\t%s%n", label, relative,
                        ReferenceQueryDump.escape(Files.readString(path)));
            } else if (name.equals("qc.metrics")) {
                final String table = Files.readAllLines(path).stream()
                        .filter(line -> !line.startsWith("#"))
                        .collect(Collectors.joining("\n"));
                System.out.printf("metrics\t%s\t%s%n", label, ReferenceQueryDump.escape(table));
            }
        }
        for (final String name : workingDirectory()) {
            System.out.printf("produced\t%s\tcwd/%s%n", label, name);
            final File file = new File(name);
            if (file.isDirectory()) {
                PrintReadsDump.emptyDirectory(file.toPath());
            }
            file.delete();
        }
    }

    /** IHDR's width and height, and pHYs's pixels per metre, or a dash where a PNG has none. */
    static String pngShape(final byte[] bytes) {
        final java.nio.ByteBuffer buffer = java.nio.ByteBuffer.wrap(bytes);
        String perMetre = "-";
        int at = 8;
        while (at + 8 <= bytes.length) {
            final int length = buffer.getInt(at);
            final String type = new String(bytes, at + 4, 4, StandardCharsets.US_ASCII);
            if (type.equals("pHYs")) {
                perMetre = Integer.toString(buffer.getInt(at + 8));
            }
            if (type.equals("IDAT")) {
                break;
            }
            at += 12 + length;
        }
        return buffer.getInt(16) + "x" + buffer.getInt(20) + "\t" + perMetre;
    }

    /** Runs the reference's scripts on their own, ahead of them the mask, and returns the status. */
    static int rscript(final String kind, final Path trace, final String... arguments) throws Exception {
        final StringBuilder expression = new StringBuilder("source('" + TRACE + "');");
        if (kind.equals("qc")) {
            expression.append("source('").append(RESOURCES.resolve("plotAlleleFrequencyQC.R")).append("');");
        } else {
            expression.append("source('").append(RESOURCES.resolve("CNVPlottingLibrary.R")).append("');");
            expression.append("source('").append(RESOURCES.resolve(
                    kind.equals("denoised") ? "PlotDenoisedCopyRatios.R" : "PlotModeledSegments.R")).append("');");
        }
        final List<String> command = new ArrayList<>(List.of("Rscript", "-e", expression.toString()));
        command.addAll(Arrays.asList(arguments));
        Files.deleteIfExists(trace);
        Files.createFile(trace);
        final ProcessBuilder builder = new ProcessBuilder(command);
        builder.environment().put("R_PLOTS_TRACE_FILE", trace.toString());
        builder.redirectErrorStream(true);
        builder.redirectOutput(ProcessBuilder.Redirect.DISCARD);
        return builder.start().waitFor();
    }

    static void emitTrace(final String label, final Path trace) throws Exception {
        for (final String row : Files.readAllLines(trace)) {
            System.out.printf("trace\t%s\t%s%n", label, ReferenceQueryDump.escape(row));
        }
        for (final String name : workingDirectory()) {
            final File file = new File(name);
            if (file.isFile()) {
                file.delete();
            }
        }
    }

    /** The trace of an `argv` row. */
    static void traced(final String label, final String kind, final String... arguments) throws Exception {
        final Path trace = DIR.resolve(label + ".trace");
        rscript(kind, trace, arguments);
        emitTrace(label, trace);
    }

    /** A script run on its own: its status, what it left, and its trace. */
    static void script(final String label, final String kind, final String... arguments) throws Exception {
        System.out.printf("script\t%s\t%s\t%s%n", label, kind, ReferenceQueryDump.escape(String.join("\t", arguments)));
        final Path outDir = DIR.resolve("out").resolve(label);
        Files.createDirectories(outDir);
        final Path trace = DIR.resolve(label + ".trace");
        final int status = rscript(kind, trace, arguments);
        System.out.printf("script-status\t%s\t%d%n", label, status);
        for (final String row : Files.readAllLines(trace)) {
            System.out.printf("trace\t%s\t%s%n", label, ReferenceQueryDump.escape(row));
        }
        listOutput(label, outDir);
    }
}

/*
 * VariantRecalibrator's tranches file and its recalibration table, taken from the reference.
 *
 * The training half of VQSR. A Gaussian mixture is fitted to the annotations of the variants a
 * resource marks as training sites, every variant is scored against it, and the scores are then
 * cut into tranches by TRUTH SENSITIVITY: what fraction of the truth sites a cutoff still calls.
 * The tranche arithmetic is the half of the tool that does not depend on the model, and it is
 * what this measures against the scores the model produced.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - A TRANCHE IS THE LARGEST SET OF VARIANTS WHOSE RUNNING SENSITIVITY REACHES THE TARGET,
 *     found by walking the LOD-sorted variants upwards from the worst;
 *   - THE RUNNING SENSITIVITY IS COMPUTED FROM THE TOP DOWN, each entry being one less the
 *     truth sites at or above it over the truth sites in total;
 *   - A TRANCHE'S COUNTS ARE OVER EVERY VARIANT AT OR ABOVE ITS minVQSLod rather than over the
 *     ones the walk had reached, so the tranches nest rather than partition;
 *   - Ti/Tv IS COUNTED OVER SNPs ALONE AND ITS DENOMINATOR IS FLOORED AT ONE, so a tranche with
 *     no transversion in it reports its transition COUNT where a ratio should be;
 *   - THE FILE IS SORTED BY CALLS AT TRUTH SITES AND NOT BY THE TARGET SENSITIVITY, which is
 *     the same order only when the targets were given in increasing order: given 100, 99.9, 99
 *     and 90 the file comes out 90, 99.9, 99, 100;
 *   - AND THE SORT IS STABLE, so two targets that found the same tranche keep the order they
 *     were given in;
 *   - EACH ROW'S FILTER NAME NAMES THE PREVIOUS ROW'S TARGET AS ITS LOWER BOUND, whatever that
 *     target is, so an unsorted target list produces a band that runs backwards:
 *     `VQSRTrancheSNP99.90to99.00`;
 *   - THE FIRST ROW'S LOWER BOUND IS 0.00, there being no previous row;
 *   - A TARGET OF 0 IS REACHABLE rather than a refusal, and produces a tranche that calls no
 *     truth site at all;
 *   - THE VQSLOD TRANCHES ARE CUT ON THE SCORE ITSELF, are written under a DIFFERENT version
 *     number and a different first column, and take their minVQSLod FROM THE REQUEST rather
 *     than from the data, so a threshold no variant reaches still produces a row, an empty one;
 *   - AND THE RECALIBRATION TABLE IS A VCF of `<VQSR>` records carrying VQSLOD, the culprit
 *     annotation and the training-site flags.
 *
 * Output:
 *
 *     vcf\t<label>=<that vcf, escaped>
 *     tranches\t<label>=<the whole tranches file, escaped>
 *     recal\t<label>=<pos, ref, alt, VQSLOD, culprit and flags per record, escaped>
 *     none\t<label>=<what was not written>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: VariantRecalibratorDump
 */

import org.broadinstitute.hellbender.tools.walkers.vqsr.VariantRecalibrator;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class VariantRecalibratorDump {

    static final int CONTIG_LENGTH = 1999980;
    static final int VARIANTS = 300;

    static List<String> header() {
        return new ArrayList<>(List.of(
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=" + CONTIG_LENGTH + ">",
                "##FILTER=<ID=LOW,Description=\"Low\">",
                "##INFO=<ID=QD,Number=1,Type=Float,Description=\"Quality by depth\">",
                "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"Mapping quality\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample"));
    }

    /**
     * The fixture's variants, built from the index alone so the file is the same every time.
     *
     * Two thirds are drawn towards a high-quality cluster and one third towards a low-quality
     * one, with a spread that repeats every eleven variants. Every even variant is a transition
     * and every odd one a transversion, so the Ti/Tv of any subset is readable off the positions.
     */
    static String[] variant(final int i) {
        final int position = 1000 + i * 100;
        final boolean transition = i % 2 == 0;
        final String alternate = transition ? "G" : "C";
        final boolean good = i % 3 != 0;
        // The two annotations get DIFFERENT spreads: drawn from the same one they are perfectly
        // correlated, the covariance matrix is singular, and the negative model finds no data.
        final double qd = (good ? 22.0 : 6.0) + (((i * 37) % 11) - 5.0) * 0.5;
        final double mq = (good ? 59.0 : 41.0) + (((i * 53) % 13) - 6.0) * 0.4;
        // Every seventh variant carries a filter, which the tool drops unless told not to.
        final String filter = i % 7 == 3 ? "LOW" : "PASS";
        return new String[] {
                "chr1", String.valueOf(position), ".", "A", alternate, "100.00", filter,
                String.format("QD=%.3f;MQ=%.3f", qd, mq), "GT", "0/1"};
    }

    static String row(final String[] fields) {
        return String.join("\t", fields);
    }

    static String vcf(final List<String> sites) {
        final List<String> lines = header();
        lines.addAll(sites);
        lines.add("");
        return String.join("\n", lines);
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("variant-recalibrator-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# VariantRecalibratorDump: the tranches a truth sensitivity target "
                + "cuts out of a scored callset");

        final List<String> input = new ArrayList<>();
        final List<String> training = new ArrayList<>();
        final List<String> truth = new ArrayList<>();
        final List<String> known = new ArrayList<>();
        for (int i = 0; i < VARIANTS; i++) {
            final String[] fields = variant(i);
            input.add(row(fields));
            // The three resources overlap deliberately: a variant can be a training site, a truth
            // site, both, or neither, and `known` is the widest of the three.
            final String bare = row(new String[] {fields[0], fields[1], ".", fields[3], fields[4],
                    "100.00", "PASS", ".", "GT", "0/1"});
            if (i % 3 == 1) {
                training.add(bare);
            }
            if (i % 3 == 2) {
                truth.add(bare);
            }
            if (i % 4 != 0) {
                known.add(bare);
            }
        }
        final Path inputPath = index(write(dir, "input.vcf", vcf(input)));
        final Path trainingPath = index(write(dir, "training.vcf", vcf(training)));
        final Path truthPath = index(write(dir, "truth.vcf", vcf(truth)));
        final Path knownPath = index(write(dir, "known.vcf", vcf(known)));
        System.out.printf("vcf\tinput=%s%n", ReferenceQueryDump.escape(vcf(input)));
        System.out.printf("vcf\ttraining=%s%n", ReferenceQueryDump.escape(vcf(training)));
        System.out.printf("vcf\ttruth=%s%n", ReferenceQueryDump.escape(vcf(truth)));
        System.out.printf("vcf\tknown=%s%n", ReferenceQueryDump.escape(vcf(known)));

        // Every run below carries --minimum-bad-variants and --bad-lod-score-cutoff. The negative
        // model is fitted to the WORST variants, and the defaults that decide which are worst are
        // meant for a whole-genome callset: a thousand bad variants and a cutoff of -5. A fixture
        // this size reaches neither, and the run dies with "No data found".
        final List<String> resources = List.of(
                "--resource:tr,known=false,training=true,truth=false,prior=15.0",
                trainingPath.toString(),
                "--resource:tv,known=false,training=false,truth=true,prior=12.0",
                truthPath.toString(),
                "--resource:kn,known=true,training=false,truth=false,prior=2.0",
                knownPath.toString());

        run(dir, "four-tranches", inputPath, resources, List.of(
                "-an", "QD", "-an", "MQ", "-mode", "SNP", "--max-gaussians", "2", "--minimum-bad-variants", "30", "--bad-lod-score-cutoff", "5.0",
                "-tranche", "100.0", "-tranche", "99.9", "-tranche", "99.0", "-tranche", "90.0"));
        // The same targets given out of order, which the file does not keep.
        run(dir, "targets-out-of-order", inputPath, resources, List.of(
                "-an", "QD", "-an", "MQ", "-mode", "SNP", "--max-gaussians", "2", "--minimum-bad-variants", "30", "--bad-lod-score-cutoff", "5.0",
                "-tranche", "90.0", "-tranche", "100.0", "-tranche", "99.0", "-tranche", "99.9"));
        // One target only, whose lower bound is therefore 0.00.
        run(dir, "one-tranche", inputPath, resources, List.of(
                "-an", "QD", "-an", "MQ", "-mode", "SNP", "--max-gaussians", "2", "--minimum-bad-variants", "30", "--bad-lod-score-cutoff", "5.0",
                "-tranche", "99.0"));
        // A target no tranche reaches, first in the list and then last.
        run(dir, "target-too-low-first", inputPath, resources, List.of(
                "-an", "QD", "-an", "MQ", "-mode", "SNP", "--max-gaussians", "2", "--minimum-bad-variants", "30", "--bad-lod-score-cutoff", "5.0",
                "-tranche", "0.0"));
        run(dir, "target-too-low-last", inputPath, resources, List.of(
                "-an", "QD", "-an", "MQ", "-mode", "SNP", "--max-gaussians", "2", "--minimum-bad-variants", "30", "--bad-lod-score-cutoff", "5.0",
                "-tranche", "100.0", "-tranche", "0.0"));
        // The filtered variants brought back in, which changes every count.
        run(dir, "ignore-all-filters", inputPath, resources, List.of(
                "-an", "QD", "-an", "MQ", "-mode", "SNP", "--max-gaussians", "2", "--minimum-bad-variants", "30", "--bad-lod-score-cutoff", "5.0",
                "--ignore-all-filters", "true",
                "-tranche", "100.0", "-tranche", "99.0"));
        run(dir, "ignore-filter", inputPath, resources, List.of(
                "-an", "QD", "-an", "MQ", "-mode", "SNP", "--max-gaussians", "2", "--minimum-bad-variants", "30", "--bad-lod-score-cutoff", "5.0",
                "--ignore-filter", "LOW",
                "-tranche", "100.0", "-tranche", "99.0"));
        // The other tranche file, cut on the score rather than on sensitivity.
        run(dir, "vqslod-tranches", inputPath, resources, List.of(
                "-an", "QD", "-an", "MQ", "-mode", "SNP", "--max-gaussians", "2", "--minimum-bad-variants", "30", "--bad-lod-score-cutoff", "5.0",
                "--output-tranches-for-scatter", "true",
                "--vqslod-tranche", "10.0", "--vqslod-tranche", "0.0",
                "--vqslod-tranche", "-10.0"));
        // A threshold no variant reaches, which still produces a row.
        run(dir, "vqslod-unreachable", inputPath, resources, List.of(
                "-an", "QD", "-an", "MQ", "-mode", "SNP", "--max-gaussians", "2", "--minimum-bad-variants", "30", "--bad-lod-score-cutoff", "5.0",
                "--output-tranches-for-scatter", "true",
                "--vqslod-tranche", "100000.0", "--vqslod-tranche", "0.0"));
        // One annotation rather than two, which moves every score.
        run(dir, "one-annotation", inputPath, resources, List.of(
                "-an", "QD", "-mode", "SNP", "--max-gaussians", "2", "--minimum-bad-variants", "30", "--bad-lod-score-cutoff", "5.0",
                "-tranche", "100.0", "-tranche", "99.0"));
        // No truth resource at all, so no tranche has any sensitivity to reach.
        run(dir, "no-truth", inputPath, List.of(
                "--resource:tr,known=false,training=true,truth=false,prior=15.0",
                trainingPath.toString()), List.of(
                "-an", "QD", "-an", "MQ", "-mode", "SNP", "--max-gaussians", "2", "--minimum-bad-variants", "30", "--bad-lod-score-cutoff", "5.0",
                "-tranche", "100.0", "-tranche", "99.0"));
        // No training resource at all, which is a refusal rather than an empty model.
        run(dir, "no-training", inputPath, List.of(
                "--resource:tv,known=false,training=false,truth=true,prior=12.0",
                truthPath.toString()), List.of(
                "-an", "QD", "-an", "MQ", "-mode", "SNP", "--max-gaussians", "2", "--minimum-bad-variants", "30", "--bad-lod-score-cutoff", "5.0",
                "-tranche", "100.0"));
        // INDEL mode over a callset with no indel in it.
        run(dir, "indel-mode", inputPath, resources, List.of(
                "-an", "QD", "-an", "MQ", "-mode", "INDEL", "--max-gaussians", "2", "--minimum-bad-variants", "30", "--bad-lod-score-cutoff", "5.0",
                "-tranche", "100.0"));
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    /** Every VCF here is queried by interval, so each needs an index beside it. */
    static Path index(final Path path) throws Exception {
        htsjdk.tribble.index.IndexFactory.createLinearIndex(path.toFile(),
                new htsjdk.variant.vcf.VCFCodec()).writeBasedOnFeatureFile(path.toFile());
        return path;
    }

    static void run(final Path dir, final String label, final Path input,
                    final List<String> resources, final List<String> extra) throws Exception {
        final Path recal = dir.resolve("out-" + label + ".recal");
        final Path tranches = dir.resolve("out-" + label + ".tranches");
        final List<String> argv = new ArrayList<>(List.of(
                "-V", input.toString(),
                "-O", recal.toString(),
                "--tranches-file", tranches.toString(),
                "--dont-run-rscript", "true"));
        argv.addAll(resources);
        argv.addAll(extra);
        try {
            new VariantRecalibrator().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(cause.getMessage()), dir)));
            return;
        }
        if (Files.exists(tranches)) {
            System.out.printf("tranches\t%s=%s%n", label,
                    ReferenceQueryDump.escape(masked(Files.readString(tranches), dir)));
        } else {
            System.out.printf("none\t%s=no tranches file%n", label);
        }
        if (!Files.exists(recal)) {
            System.out.printf("none\t%s=no recalibration table%n", label);
            return;
        }
        // The table is a VCF, reported as the fields the tranche arithmetic reads off it.
        final StringBuilder body = new StringBuilder();
        for (final String line : Files.readString(recal).split("\n", -1)) {
            if (line.startsWith("#") || line.isEmpty()) {
                continue;
            }
            final String[] columns = line.split("\t");
            final StringBuilder flags = new StringBuilder();
            String lod = ".";
            String culprit = ".";
            for (final String part : columns[7].split(";")) {
                if (part.startsWith("VQSLOD=")) {
                    lod = part.substring("VQSLOD=".length());
                } else if (part.startsWith("culprit=")) {
                    culprit = part.substring("culprit=".length());
                } else if (part.endsWith("TRAIN_SITE")) {
                    flags.append(flags.length() == 0 ? "" : ",").append(part);
                }
            }
            body.append(String.join("\t", columns[1], columns[3], columns[4], lod, culprit,
                    flags.toString())).append("\n");
        }
        System.out.printf("recal\t%s=%s%n", label, ReferenceQueryDump.escape(body.toString()));
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}

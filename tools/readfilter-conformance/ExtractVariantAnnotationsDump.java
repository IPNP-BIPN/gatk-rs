/*
 * ExtractVariantAnnotations' annotation matrices and its sites-only VCF, taken from the reference.
 *
 * Which variants a scalable-VQSR extraction keeps, what labels it gives them, and what numbers it
 * writes for their annotations. One prefix produces up to THREE files, and none of them is a
 * subset of another: a labelled HDF5 matrix, an unlabelled one, and a sites-only VCF.
 *
 * Thirteen behaviours this is built to catch.
 *
 *   - THE MODE DECIDES WHICH VARIANT TYPES ARE KEPT, and a record of the wrong type is dropped
 *     from every file rather than written with no label;
 *   - A LABEL COMES FROM THE RESOURCE'S TAG, and only from a tag whose value is the string
 *     `true`: `training=false` labels nothing and so extracts nothing;
 *   - `snp` IS RESERVED, because the matrix carries a `snp` label of its own for every run
 *     whatever the resources say, and a resource tagged with it is refused;
 *   - AN UNLABELLED VARIANT IS DROPPED unless --maximum-number-of-unlabeled-variants asks for
 *     some, and then it goes to a SEPARATE FILE, `<prefix>.unlabeled.annot.hdf5`, never into the
 *     labelled one;
 *   - THE RESERVOIR IS NOT IN GENOMIC ORDER and the seed decides both which records it keeps and
 *     what order they land in: seeds 0, 1 and 100 keep the same two here and seed 42 keeps a
 *     different pair, later position first;
 *   - THE ANNOTATION COLUMNS ARE SORTED BY NAME, whatever order they were asked for in;
 *   - AN ANNOTATION A RECORD DOES NOT CARRY IS NaN, and so is an infinite one: the two absences
 *     are not told apart;
 *   - A FILTERED RECORD IS DROPPED, and --ignore-filter names one filter while
 *     --ignore-all-filters takes them all;
 *   - THE MATCHING STRATEGY DECIDES WHAT COUNTS AS THE SAME VARIANT, and THE DEFAULT IS
 *     START_POSITION, which matches on the position alone: a resource record with an entirely
 *     different alternate still labels the input;
 *   - ONLY THE MINIMAL REPRESENTATION RECONCILES A PADDED ALLELE, so the same insertion written
 *     as `C>CAT` and as `CG>CATG` matches under that strategy and under no other;
 *   - ASKING FOR AN ALLELE-SPECIFIC ANNOTATION SWITCHES THE WHOLE RUN to one row per alternate,
 *     and a multiallelic record becomes two rows at the same position;
 *   - THE ALT ALLELE ARRAY IS FLAT, so a multiallelic record contributes both alternates and the
 *     array is longer than the reference one beside it;
 *   - AND --omit-alleles-in-hdf5 REMOVES THE ALLELE DATASETS FROM THE MATRIX ALONE, while a run
 *     that extracts nothing writes no matrix at all and still writes its VCF.
 *
 * Output:
 *
 *     vcf\t<label>=<that vcf, escaped>
 *     out\t<label>=<the sites-only vcf without its header, escaped>
 *     names\t<label>=<the annotation columns, comma-separated>
 *     intervals\t<label>=<one interval per row, comma-separated>
 *     alleles\t<label><path>=<that allele dataset, comma-separated>
 *     rows\t<label>=<one row per record, the annotations rendered by Double.toString>
 *     labels\t<label>=<one line per label, `<label>: <true/false per row>`>
 *     none\t<label>=<what was not written>
 *     error\t<label>\t<exception class>:<message>
 *
 * The unlabelled matrix is reported under `<label>.unlabeled`.
 *
 * Usage: ExtractVariantAnnotationsDump
 */

import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.ExtractVariantAnnotations;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.LabeledVariantAnnotationsData;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class ExtractVariantAnnotationsDump {

    static final int CONTIG_LENGTH = 199980;

    static List<String> header() {
        return new ArrayList<>(List.of(
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=" + CONTIG_LENGTH + ">",
                "##FILTER=<ID=LOW,Description=\"Low\">",
                "##FILTER=<ID=BAD,Description=\"Bad\">",
                "##INFO=<ID=QD,Number=1,Type=Float,Description=\"Quality by depth\">",
                "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"Mapping quality\">",
                "##INFO=<ID=FS,Number=1,Type=Float,Description=\"Strand bias\">",
                "##INFO=<ID=AS_QD,Number=A,Type=Float,Description=\"Quality by depth per alt\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample"));
    }

    static String site(final int position, final String reference, final String alternate,
                       final String filter, final String info) {
        return "chr1\t" + position + "\t.\t" + reference + "\t" + alternate
                + "\t100.00\t" + filter + "\t" + info + "\tGT\t0/1";
    }

    static String vcf(final List<String> sites) {
        final List<String> lines = header();
        lines.addAll(sites);
        lines.add("");
        return String.join("\n", lines);
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("extract-variant-annotations-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# ExtractVariantAnnotationsDump: which variants a scalable-VQSR "
                + "extraction keeps, and what it writes for them");

        final Path fasta = writeReference(dir);

        // The input. Two SNPs and two indels, one record with a filter, one whose QD is missing
        // and one whose MQ is infinite, and one multiallelic site for the allele-specific runs.
        final String input = vcf(List.of(
                site(1000, "A", "C", "PASS", "QD=1.5;MQ=60.0;FS=0.5;AS_QD=1.5"),
                site(2000, "G", "T", "PASS", "QD=2.5;MQ=59.0;FS=1.5;AS_QD=2.5"),
                site(3000, "C", "CAT", "PASS", "QD=3.5;MQ=58.0;FS=2.5;AS_QD=3.5"),
                site(4000, "TGG", "T", "PASS", "QD=4.5;MQ=57.0;FS=3.5;AS_QD=4.5"),
                site(5000, "A", "G", "LOW", "QD=5.5;MQ=56.0;FS=4.5;AS_QD=5.5"),
                // No QD at all, and an MQ the reader turns into a number it cannot hold.
                site(6000, "A", "T", "PASS", "MQ=Infinity;FS=5.5;AS_QD=6.5"),
                site(7000, "A", "C,G", "PASS", "QD=7.5;MQ=55.0;FS=6.5;AS_QD=7.5,8.5")));
        final Path inputPath = index(write(dir, "input.vcf", input));
        System.out.printf("vcf\tinput=%s%n", ReferenceQueryDump.escape(input));

        // The training resource: the first SNP and the first indel, one of them with a DIFFERENT
        // representation of the same event, so the matching strategies can be told apart.
        final String training = vcf(List.of(
                site(1000, "A", "C", "PASS", "."),
                // The SAME event as the input's insertion at 3000, written with a padding base,
                // which only a reconciled minimal representation recognises.
                site(3000, "CG", "CATG", "PASS", "."),
                site(7000, "A", "T", "PASS", ".")));
        final Path trainingPath = index(write(dir, "training.vcf", training));
        System.out.printf("vcf\ttraining=%s%n", ReferenceQueryDump.escape(training));

        // A second resource under a second label, overlapping the first at one site only.
        final String calibration = vcf(List.of(
                site(1000, "A", "C", "PASS", "."),
                site(2000, "G", "T", "PASS", ".")));
        final Path calibrationPath = index(write(dir, "calibration.vcf", calibration));
        System.out.printf("vcf\tcalibration=%s%n", ReferenceQueryDump.escape(calibration));

        final String trainingTag = "--resource:train,training=true";
        final String calibrationTag = "--resource:cal,calibration=true";

        run(dir, "snp", inputPath, fasta, List.of(
                trainingTag, trainingPath.toString(),
                "-A", "QD", "-A", "MQ", "--mode", "SNP"));
        run(dir, "indel", inputPath, fasta, List.of(
                trainingTag, trainingPath.toString(),
                "-A", "QD", "-A", "MQ", "--mode", "INDEL"));
        run(dir, "both-modes", inputPath, fasta, List.of(
                trainingTag, trainingPath.toString(),
                "-A", "QD", "-A", "MQ", "--mode", "SNP", "--mode", "INDEL"));
        // Two labels, and the annotations asked for out of order.
        run(dir, "two-labels", inputPath, fasta, List.of(
                trainingTag, trainingPath.toString(),
                calibrationTag, calibrationPath.toString(),
                "-A", "MQ", "-A", "FS", "-A", "QD", "--mode", "SNP", "--mode", "INDEL"));
        // A tag whose value is not the string `true`, which is not a label.
        run(dir, "false-tag", inputPath, fasta, List.of(
                "--resource:train,training=false", trainingPath.toString(),
                "-A", "QD", "--mode", "SNP"));
        // The reserved label.
        run(dir, "reserved-label", inputPath, fasta, List.of(
                "--resource:train,snp=true", trainingPath.toString(),
                "-A", "QD", "--mode", "SNP"));
        // The unlabelled variants, which are dropped by default and reservoir-sampled when asked
        // for, with the seed deciding which.
        // MQ is asked for here too, because the record with no QD at all is unlabelled and so is
        // the one whose MQ the reader turns into an infinity.
        run(dir, "unlabeled-all", inputPath, fasta, List.of(
                trainingTag, trainingPath.toString(),
                "-A", "QD", "-A", "MQ", "--mode", "SNP", "--mode", "INDEL",
                "--maximum-number-of-unlabeled-variants", "10"));
        // Two of the four unlabelled records, under two seeds.
        run(dir, "unlabeled-two-seed-0", inputPath, fasta, List.of(
                trainingTag, trainingPath.toString(),
                "-A", "QD", "--mode", "SNP", "--mode", "INDEL",
                "--ignore-all-filters", "true",
                "--maximum-number-of-unlabeled-variants", "2",
                "--reservoir-sampling-random-seed", "0"));
        run(dir, "unlabeled-two-seed-1", inputPath, fasta, List.of(
                trainingTag, trainingPath.toString(),
                "-A", "QD", "--mode", "SNP", "--mode", "INDEL",
                "--ignore-all-filters", "true",
                "--maximum-number-of-unlabeled-variants", "2",
                "--reservoir-sampling-random-seed", "1"));
        run(dir, "unlabeled-two-seed-42", inputPath, fasta, List.of(
                trainingTag, trainingPath.toString(),
                "-A", "QD", "--mode", "SNP", "--mode", "INDEL",
                "--ignore-all-filters", "true",
                "--maximum-number-of-unlabeled-variants", "2",
                "--reservoir-sampling-random-seed", "42"));
        run(dir, "unlabeled-two-seed-100", inputPath, fasta, List.of(
                trainingTag, trainingPath.toString(),
                "-A", "QD", "--mode", "SNP", "--mode", "INDEL",
                "--ignore-all-filters", "true",
                "--maximum-number-of-unlabeled-variants", "2",
                "--reservoir-sampling-random-seed", "100"));
        // The filters.
        run(dir, "ignore-filter", inputPath, fasta, List.of(
                trainingTag, trainingPath.toString(),
                "-A", "QD", "--mode", "SNP", "--mode", "INDEL",
                "--maximum-number-of-unlabeled-variants", "10",
                "--ignore-filter", "LOW"));
        run(dir, "ignore-all-filters", inputPath, fasta, List.of(
                trainingTag, trainingPath.toString(),
                "-A", "QD", "--mode", "SNP", "--mode", "INDEL",
                "--maximum-number-of-unlabeled-variants", "10",
                "--ignore-all-filters", "true"));
        // The three matching strategies over the site whose representation differs.
        for (final String strategy : List.of("START_POSITION",
                "START_POSITION_AND_GIVEN_REPRESENTATION",
                "START_POSITION_AND_MINIMAL_REPRESENTATION")) {
            run(dir, "match-" + strategy.toLowerCase(), inputPath, fasta, List.of(
                    trainingTag, trainingPath.toString(),
                    "-A", "QD", "--mode", "SNP", "--mode", "INDEL",
                    "--resource-matching-strategy", strategy));
        }
        // An allele-specific annotation, which switches the whole run to one row per alternate.
        run(dir, "allele-specific", inputPath, fasta, List.of(
                trainingTag, trainingPath.toString(),
                "-A", "AS_QD", "--mode", "SNP", "--mode", "INDEL",
                "--maximum-number-of-unlabeled-variants", "10"));
        // The alleles left out of the HDF5, which is a different file for the same extraction.
        run(dir, "omit-alleles", inputPath, fasta, List.of(
                trainingTag, trainingPath.toString(),
                "-A", "QD", "--mode", "SNP", "--omit-alleles-in-hdf5", "true"));
        // A run that keeps nothing at all: no resource, and no unlabelled variants asked for.
        run(dir, "extracts-nothing", inputPath, fasta, List.of(
                "-A", "QD", "--mode", "SNP", "--mode", "INDEL"));
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    /** Every VCF here is QUERIED by interval, so each needs an index beside it. */
    static Path index(final Path path) throws Exception {
        htsjdk.tribble.index.IndexFactory.createLinearIndex(path.toFile(),
                new htsjdk.variant.vcf.VCFCodec()).writeBasedOnFeatureFile(path.toFile());
        return path;
    }

    static void run(final Path dir, final String label, final Path input, final Path fasta,
                    final List<String> extra) throws Exception {
        final Path prefix = dir.resolve("out-" + label);
        final List<String> argv = new ArrayList<>(List.of(
                "-V", input.toString(),
                "-O", prefix.toString(),
                "-R", fasta.toString(),
                "--do-not-gzip-vcf-output", "true"));
        argv.addAll(extra);
        try {
            new ExtractVariantAnnotations().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(cause.getMessage()), dir)));
            return;
        }

        final Path vcf = dir.resolve("out-" + label + ".vcf");
        if (Files.exists(vcf)) {
            final StringBuilder body = new StringBuilder();
            for (final String line : Files.readString(vcf).split("\n", -1)) {
                if (!line.startsWith("##") && !line.isEmpty()) {
                    body.append(line).append("\n");
                }
            }
            System.out.printf("out\t%s=%s%n", label,
                    ReferenceQueryDump.escape(masked(body.toString(), dir)));
        } else {
            System.out.printf("none\t%s=no vcf%n", label);
        }

        // The labelled matrix, and beside it the reservoir of unlabelled ones, which is a
        // SEPARATE file written only when some were asked for.
        matrix(label, "", new File(prefix + ".annot.hdf5"));
        matrix(label, ".unlabeled", new File(prefix + ".unlabeled.annot.hdf5"));
    }

    /** One HDF5 matrix: its columns, its intervals, its alleles, its rows and its labels. */
    static void matrix(final String label, final String suffix, final File hdf5) {
        final String name = label + suffix;
        if (!hdf5.exists()) {
            // A run that extracted nothing writes no matrix, which is a warning and not an error.
            System.out.printf("none\t%s=no annotations hdf5%n", name);
            return;
        }
        System.out.printf("names\t%s=%s%n", name,
                String.join(",", LabeledVariantAnnotationsData.readAnnotationNames(hdf5)));

        try (final org.broadinstitute.hdf5.HDF5File file =
                     new org.broadinstitute.hdf5.HDF5File(hdf5)) {
            final List<org.broadinstitute.hellbender.utils.SimpleInterval> intervals =
                    org.broadinstitute.hellbender.tools.copynumber.utils.HDF5Utils.readIntervals(
                            file, LabeledVariantAnnotationsData.INTERVALS_PATH);
            final List<String> places = new ArrayList<>();
            for (final org.broadinstitute.hellbender.utils.SimpleInterval interval : intervals) {
                places.add(interval.getContig() + ":" + interval.getStart() + "-"
                        + interval.getEnd());
            }
            System.out.printf("intervals\t%s=%s%n", name, String.join(",", places));

            // The alleles, which --omit-alleles-in-hdf5 leaves out of this file alone.
            for (final String path : List.of(LabeledVariantAnnotationsData.ALLELES_REF_PATH,
                    LabeledVariantAnnotationsData.ALLELES_ALT_PATH)) {
                try {
                    System.out.printf("alleles\t%s%s=%s%n", name, path,
                            String.join(",", file.readStringArray(path)));
                } catch (final Exception e) {
                    System.out.printf("none\t%s%s=no alleles%n", name, path);
                }
            }
        }

        final double[][] rows = LabeledVariantAnnotationsData.readAnnotations(hdf5);
        final StringBuilder cells = new StringBuilder();
        for (final double[] row : rows) {
            final List<String> values = new ArrayList<>();
            for (final double value : row) {
                values.add(Double.toString(value));
            }
            cells.append(String.join("\t", values)).append("\n");
        }
        System.out.printf("rows\t%s=%s%n", name, ReferenceQueryDump.escape(cells.toString()));

        // The labels the matrix carries, `snp` among them: it is written for every run whatever
        // the resources say, which is why the name is reserved.
        final StringBuilder flags = new StringBuilder();
        for (final String key : List.of("snp", "training", "calibration")) {
            final List<Boolean> values;
            try {
                values = LabeledVariantAnnotationsData.readLabel(hdf5, key);
            } catch (final Exception e) {
                continue;
            }
            final List<String> written = new ArrayList<>();
            for (final Boolean value : values) {
                written.add(String.valueOf(value));
            }
            flags.append(key).append(": ").append(String.join(",", written)).append("\n");
        }
        System.out.printf("labels\t%s=%s%n", name, ReferenceQueryDump.escape(flags.toString()));
    }

    static Path writeReference(final Path dir) throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        final StringBuilder bases = new StringBuilder(">chr1\n");
        for (int i = 0; i < CONTIG_LENGTH / 60; i++) {
            bases.append("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n");
        }
        Files.writeString(fasta, bases.toString(), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        final htsjdk.samtools.SAMFileHeader header = new htsjdk.samtools.SAMFileHeader();
        header.setSequenceDictionary(new htsjdk.samtools.SAMSequenceDictionary(List.of(
                new htsjdk.samtools.SAMSequenceRecord("chr1", CONTIG_LENGTH))));
        try (final java.io.Writer writer = Files.newBufferedWriter(dir.resolve("reference.dict"))) {
            new htsjdk.samtools.SAMTextHeaderCodec().encode(writer, header);
        }
        return fasta;
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}

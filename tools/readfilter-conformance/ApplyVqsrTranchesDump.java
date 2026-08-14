/*
 * The tranches file, as ApplyVQSR reads it, taken from the reference.
 *
 * `TruthSensitivityTranche.readTranches` and what `onTraversalStart` makes of the tranches that
 * survive `--truth-sensitivity-filter-level`: a set of FILTER header lines, which is the whole of
 * the tool's output when the input VCF has no record in it.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE TWO OPTIONAL COLUMNS ARE NOT OPTIONAL. `numKnown` is read with `getOptionalInteger(...,
 *     -1)` and the constructor then refuses `numKnown < 0` with a GATKException, so a header that
 *     does not name the column is fatal rather than defaulted. `accessibleTruthSites` and
 *     `callsAtTruthSites` default to -1 and are not checked at all;
 *   - THE HEADER NAMES THE COLUMNS AND THERE MUST BE ELEVEN OF THEM. The bindings are built
 *     header-to-value, so the order is free, but a header of any other length is refused before any
 *     row is read and a row of a different length is refused against the header's;
 *   - THE SAME READER WORDS ITS TWO REFUSALS DIFFERENTLY: the missing-key path builds its
 *     `MalformedFile` with no file at all, so it reads "Unknown file is malformed", while the two
 *     length checks name the file they were given;
 *   - `model` AND `filterName` ARE READ WITH A BARE `bindings.get`, so a header that does not name
 *     the model column reaches `Mode.valueOf(null)` and comes out `NullPointerException: Name is
 *     null`, while `numNovel` is a `long` field parsed with `Integer.valueOf`, so a count the writer
 *     could have produced cannot be read back past 2^31;
 *   - THE FILTER IDs ARE THE filterName COLUMN, not anything synthesized: a tranches file naming its
 *     tranches something else produces FILTER lines with those names;
 *   - EACH TRANCHE'S INTERVAL IS DESCRIBED WITH THE NEXT TRANCHE'S minVQSLod, and the first tranche
 *     gets a SECOND line with a `+` appended to its own name and an open-ended description, so one
 *     tranche can appear twice in the header under two IDs, while THE LAST TRANCHE NEVER BECOMES A
 *     FILTER AT ALL, the loop running to `size - 1`: the tranche whose variants are all kept is the
 *     one with no line;
 *   - THE ORDER IS BY targetTruthSensitivity, then filtered by the level and REVERSED, so the tranche
 *     that is kept whole is the last one and the `+` line belongs to the most specific;
 *   - A LEVEL ABOVE EVERY TRANCHE IS A UserException rather than an empty filter set, and the
 *     mutual exclusion of `--truth-sensitivity-filter-level` and `--lod-score-cutoff` is checked
 *     AFTER the tranches file has been read, so a broken file wins over the mutex;
 *   - AND WITH NO LEVEL AT ALL there is one filter line, `LOW_VQSLOD`, described with the default
 *     cutoff of 0.0.
 *
 * Output:
 *
 *     tranches\t<label>\t<the whole tranches file, escaped>
 *     filter\t<run>\t<one ##FILTER line of the output vcf, escaped>
 *     error\t<run>\t<exception class>:<message>
 *
 * Usage: ApplyVqsrTranchesDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.vqsr.ApplyVQSR;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class ApplyVqsrTranchesDump {

    /** No record: the header the tool builds is the whole measurement. */
    static final String VCF =
            "##fileformat=VCFv4.2\n"
                    + "##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=\"the score\">\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##contig=<ID=chr1,length=1000>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n";

    static final String COLUMNS =
            "targetTruthSensitivity,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,model,accessibleTruthSites,callsAtTruthSites,truthSensitivity";

    static final String PREAMBLE =
            "# Variant quality score tranches file\n# Version number 5\n";

    static final String ROW_90 =
            "90.00,10,5,2.1000,1.9000,3.5000,VQSRTrancheSNP0.00to90.00,SNP,100,90,0.9000";
    static final String ROW_99 =
            "99.00,20,9,2.0000,1.8000,1.5000,VQSRTrancheSNP90.00to99.00,SNP,100,99,0.9900";
    static final String ROW_100 =
            "100.00,30,15,1.9000,1.7000,-0.5000,VQSRTrancheSNP99.00to100.00,SNP,100,100,1.0000";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("apply-vqsr-tranches-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# ApplyVqsrTranchesDump: eleven columns, and the filters they become");

        final Path variants = writeVcf(dir, "variants");
        final Path recal = writeVcf(dir, "recal");

        // Three tranches, deliberately written out of sensitivity order so that the sort is visible.
        final Path tranches = writeTranches(dir, "tranches", PREAMBLE + COLUMNS + "\n"
                + ROW_100 + "\n" + ROW_90 + "\n" + ROW_99 + "\n");

        // The same three under names of their own.
        final Path names = writeTranches(dir, "custom-names", PREAMBLE + COLUMNS + "\n"
                + ROW_90.replace("VQSRTrancheSNP0.00to90.00", "loose")
                + "\n" + ROW_99.replace("VQSRTrancheSNP90.00to99.00", "middling")
                + "\n" + ROW_100.replace("VQSRTrancheSNP99.00to100.00", "tight") + "\n");

        // A required key the header does not name.
        final Path missingRequired = writeTranches(dir, "missing-required",
                PREAMBLE + COLUMNS.replace("minVQSLod", "minVQSLOD") + "\n" + ROW_99 + "\n");

        // An OPTIONAL key the header does not name, which is the one that is not optional.
        final Path missingOptional = writeTranches(dir, "missing-optional",
                PREAMBLE + COLUMNS.replace("numKnown", "numknown") + "\n" + ROW_99 + "\n");

        // A model the enum does not have.
        final Path badModel = writeTranches(dir, "bad-model",
                PREAMBLE + COLUMNS + "\n" + ROW_99.replace(",SNP,", ",GERMLINE,") + "\n");

        // A sensitivity the constructor calls unreasonable.
        final Path unreasonable = writeTranches(dir, "unreasonable",
                PREAMBLE + COLUMNS + "\n" + ROW_99.replace("99.00,20", "150.00,20") + "\n");

        // A value that is not a number at all.
        final Path invalidValue = writeTranches(dir, "invalid-value",
                PREAMBLE + COLUMNS + "\n" + ROW_99.replace(",20,9,", ",20,many,") + "\n");

        // A count `numNovel` holds as a long and reads back with Integer.valueOf.
        final Path bigNovel = writeTranches(dir, "novel-past-int",
                PREAMBLE + COLUMNS + "\n" + ROW_99.replace(",20,9,", ",20,3000000000,") + "\n");

        // A header that does not name the model column, which nothing checks for null.
        final Path missingModel = writeTranches(dir, "missing-model",
                PREAMBLE + COLUMNS.replace(",model,", ",mdl,") + "\n" + ROW_99 + "\n");

        // A row with one field fewer than the header.
        final Path shortRow = writeTranches(dir, "short-row",
                PREAMBLE + COLUMNS + "\n" + ROW_99.substring(0, ROW_99.lastIndexOf(',')) + "\n");

        // A header with one column fewer.
        final Path shortHeader = writeTranches(dir, "short-header",
                PREAMBLE + COLUMNS.substring(0, COLUMNS.lastIndexOf(',')) + "\n" + ROW_99 + "\n");

        run(dir, "level-99", variants, recal, tranches, "99.0", null);
        run(dir, "level-0", variants, recal, tranches, "0.0", null);
        run(dir, "level-above-everything", variants, recal, tranches, "100.1", null);
        run(dir, "custom-names", variants, recal, names, "0.0", null);
        run(dir, "no-level", variants, recal, tranches, null, null);
        run(dir, "both-cutoffs", variants, recal, tranches, "99.0", "0.5");
        // The mutex is checked after the file is read, so this one reports the file's problem.
        run(dir, "both-cutoffs-broken-file", variants, recal, missingRequired, "99.0", "0.5");
        run(dir, "missing-required", variants, recal, missingRequired, "99.0", null);
        run(dir, "missing-optional", variants, recal, missingOptional, "99.0", null);
        run(dir, "bad-model", variants, recal, badModel, "99.0", null);
        run(dir, "missing-model", variants, recal, missingModel, "99.0", null);
        run(dir, "invalid-value", variants, recal, invalidValue, "99.0", null);
        run(dir, "novel-past-int", variants, recal, bigNovel, "99.0", null);
        run(dir, "unreasonable", variants, recal, unreasonable, "99.0", null);
        run(dir, "short-row", variants, recal, shortRow, "99.0", null);
        run(dir, "short-header", variants, recal, shortHeader, "99.0", null);
    }

    static Path writeVcf(final Path dir, final String label) throws Exception {
        final Path file = dir.resolve(label + ".vcf");
        Files.writeString(file, VCF, StandardCharsets.UTF_8);
        new IndexFeatureFile().instanceMain(new String[] {"-I", file.toString()});
        return file;
    }

    static Path writeTranches(final Path dir, final String label, final String text) throws Exception {
        final Path file = dir.resolve(label + ".tranches");
        Files.writeString(file, text, StandardCharsets.UTF_8);
        System.out.printf("tranches\t%s\t%s%n", label, ReferenceQueryDump.escape(text));
        return file;
    }

    static void run(final Path dir, final String label, final Path variants, final Path recal,
                    final Path tranches, final String level, final String lodCutoff) {
        final Path output = dir.resolve(label + ".vcf");
        final List<String> all = new ArrayList<>(List.of(
                "-V", variants.toString(),
                "--recal-file", recal.toString(),
                "--tranches-file", tranches.toString(),
                "-O", output.toString()));
        if (level != null) {
            all.add("--truth-sensitivity-filter-level");
            all.add(level);
        }
        if (lodCutoff != null) {
            all.add("--lod-score-cutoff");
            all.add(lodCutoff);
        }
        try {
            new ApplyVQSR().instanceMain(all.toArray(new String[0]));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        try {
            for (final String line : Files.readAllLines(output, StandardCharsets.UTF_8)) {
                if (line.startsWith("##FILTER=")) {
                    System.out.printf("filter\t%s\t%s%n", label, ReferenceQueryDump.escape(line));
                }
            }
        } catch (final Exception e) {
            System.out.printf("error\t%s-read\t%s:%s%n", label, e.getClass().getName(),
                    String.valueOf(e.getMessage()));
        }
    }

    static void emptyDirectory(final Path dir) throws Exception {
        if (!Files.isDirectory(dir)) {
            return;
        }
        try (final var entries = Files.list(dir)) {
            for (final Path entry : entries.toList()) {
                Files.deleteIfExists(entry);
            }
        }
    }
}

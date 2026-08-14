/*
 * ApplyVQSR's per-record filtering, taken from the reference.
 *
 * `apply` decides whether a record is recalibrated at all, `doSiteSpecificFiltering` finds its recal
 * record and reads the LOD out of it, and `generateFilterString` turns that LOD into a FILTER value.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE TRANCHES ARE WALKED BACKWARDS AND THE LAST ONE MEANS PASS. The loop runs from the end of
 *     the reversed list, takes the first tranche whose `minVQSLod` the record reaches, and answers
 *     `PASS` only when that tranche is the last one; anything below every tranche gets the FIRST
 *     tranche's name with a `+`. The bands are therefore exactly the intervals the header lines
 *     describe, and the boundary belongs to the tranche below it, `lod >= tranche.minVQSLod`;
 *   - THE RECAL RECORD MUST AGREE ON BOTH ENDS, BY TWO DIFFERENT MECHANISMS. The query is
 *     `featureContext.getValues(recal, vc.getStart())`, which keeps only records STARTING at the
 *     input record's start, and `getMatchingRecalVC` then takes the first of those whose END agrees.
 *     A recal record at the same start with another length is skipped, and one with no partner that
 *     ends where the input record does is a refusal;
 *   - A RECORD OF THE WRONG CLASS IS EMITTED UNTOUCHED, with no VQSLOD and no filter, so a run in
 *     SNP mode passes indels through exactly as they came in;
 *   - AND SO IS A RECORD THAT WAS ALREADY FILTERED, unless `--ignore-all-filters` was given or
 *     `--ignore-filter` names every filter it carries. `--exclude-filtered` NEVER DROPS IT: the
 *     exclusion sits inside the branch that recalibrates, so a run that drops every record the tool
 *     filtered still writes out the records that arrived filtered;
 *   - EVERY NEGATIVE VQSLOD IS WRITTEN IN SCIENTIFIC NOTATION. The value is reparsed with
 *     `Double.valueOf` and written back as a double, so it goes through htsjdk's `formatVCFDouble`,
 *     whose branch is on the SIGNED value rather than the magnitude: `-3.0` comes out
 *     `-3.000e+00` while `5.0` comes out `5.00` and `0.0` comes out `0.00`;
 *   - THE CULPRIT IS COPIED WITH getAttribute AND NO DEFAULT, so a recal record that does not carry
 *     one writes `culprit=.`;
 *   - THE TWO TRAINING LABELS ARE COPIED AS `true` when present, whatever they held;
 *   - AND THE THREE REFUSALS ARE WRAPPED: what `apply` throws is a `UserException`, and the walker
 *     rethrows it as a `GATKException` whose message is the locus and the whole `VariantContext`, so
 *     the tool's own wording survives only in the cause. The three are a missing recal record, a
 *     recal record with no LOD, and a LOD that will not parse.
 *
 * Output:
 *
 *     input\t<label>\t<the whole vcf, escaped>
 *     vcfline\t<run>\t<one record line of the output vcf, escaped>
 *     error\t<run>\t<exception class>:<message>
 *     cause\t<run>\t<the wrapped exception's class>:<message>
 *
 * Usage: ApplyVqsrSiteFilteringDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.vqsr.ApplyVQSR;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class ApplyVqsrSiteFilteringDump {

    static final String INPUT_HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##FILTER=<ID=weak,Description=\"Was already there\">\n"
                    + "##contig=<ID=chr1,length=1000>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n";

    /** The recal file VariantRecalibrator writes: a VQSLOD, a culprit and the two training labels. */
    static final String RECAL_HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=\"the score\">\n"
                    + "##INFO=<ID=culprit,Number=1,Type=String,Description=\"the worst annotation\">\n"
                    + "##INFO=<ID=POSITIVE_TRAIN_SITE,Number=0,Type=Flag,Description=\"a positive training site\">\n"
                    + "##INFO=<ID=NEGATIVE_TRAIN_SITE,Number=0,Type=Flag,Description=\"a negative training site\">\n"
                    + "##contig=<ID=chr1,length=1000>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    static final String TRANCHES =
            "# Variant quality score tranches file\n"
                    + "# Version number 5\n"
                    + "targetTruthSensitivity,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,model,accessibleTruthSites,callsAtTruthSites,truthSensitivity\n"
                    + "90.00,10,5,2.1000,1.9000,3.5000,VQSRTrancheSNP0.00to90.00,SNP,100,90,0.9000\n"
                    + "99.00,20,9,2.0000,1.8000,1.5000,VQSRTrancheSNP90.00to99.00,SNP,100,99,0.9900\n"
                    + "100.00,30,15,1.9000,1.7000,-0.5000,VQSRTrancheSNP99.00to100.00,SNP,100,100,1.0000\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("apply-vqsr-site-filtering-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# ApplyVqsrSiteFilteringDump: one lod, one band, one filter");

        // One record per band, a boundary, an indel, an already-filtered record, and a record whose
        // recal partner begins five bases earlier and ends where it does.
        final Path variants = writeVcf(dir, "variants", INPUT_HEADER,
                "chr1\t100\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t200\t.\tA\tC\t50\t.\t.\tGT\t0/1",
                "chr1\t300\t.\tA\tC\t50\t.\t.\tGT\t0/1",
                "chr1\t400\t.\tA\tC\t50\t.\t.\tGT\t0/1",
                "chr1\t500\t.\tACC\tA\t50\t.\t.\tGT\t0/1",
                "chr1\t600\t.\tA\tC\t50\tweak\t.\tGT\t0/1",
                "chr1\t700\t.\tA\tC\t50\t.\t.\tGT\t0/1");
        final Path recal = writeVcf(dir, "recal", RECAL_HEADER,
                "chr1\t100\t.\tA\tC\t.\t.\tVQSLOD=5.0;culprit=QD;POSITIVE_TRAIN_SITE",
                "chr1\t200\t.\tA\tC\t.\t.\tVQSLOD=2.0;culprit=MQ",
                "chr1\t300\t.\tA\tC\t.\t.\tVQSLOD=0.0;culprit=FS;NEGATIVE_TRAIN_SITE",
                // No culprit at all, which is copied with no default.
                "chr1\t400\t.\tA\tC\t.\t.\tVQSLOD=-3.0",
                "chr1\t500\t.\tACC\tA\t.\t.\tVQSLOD=4.0;culprit=QD",
                "chr1\t600\t.\tA\tC\t.\t.\tVQSLOD=4.0;culprit=QD",
                // Exactly on a tranche's minVQSLod.
                "chr1\t700\t.\tA\tC\t.\t.\tVQSLOD=1.5;culprit=SOR");

        // The end test, on top of a query that already fixed the start: two recal records begin
        // where the input record does and only the second ends where it does.
        final Path ends = writeVcf(dir, "ends", INPUT_HEADER,
                "chr1\t800\t.\tA\tC\t50\t.\t.\tGT\t0/1");
        // And one whose only candidate begins where it does and ends five bases later.
        final Path endsMismatch = writeVcf(dir, "ends-mismatch", INPUT_HEADER,
                "chr1\t850\t.\tA\tC\t50\t.\t.\tGT\t0/1");
        final Path endsRecal = writeVcf(dir, "ends-recal", RECAL_HEADER,
                "chr1\t800\t.\tACCCCC\tA\t.\t.\tVQSLOD=9.0;culprit=SKIPPED",
                "chr1\t800\t.\tA\tC\t.\t.\tVQSLOD=5.0;culprit=TAKEN",
                // Nothing at 850 ends where the input record does.
                "chr1\t850\t.\tACCCCC\tA\t.\t.\tVQSLOD=9.0;culprit=NEVER");

        // One record with no recal partner at all.
        final Path orphan = writeVcf(dir, "orphan", INPUT_HEADER,
                "chr1\t900\t.\tA\tC\t50\t.\t.\tGT\t0/1");

        // A recal record carrying no VQSLOD, and one whose VQSLOD will not parse.
        final Path noLod = writeVcf(dir, "no-lod", RECAL_HEADER,
                "chr1\t900\t.\tA\tC\t.\t.\tculprit=QD");
        final Path badLod = writeVcf(dir, "bad-lod", RECAL_HEADER,
                "chr1\t900\t.\tA\tC\t.\t.\tVQSLOD=nonsense;culprit=QD");

        final Path tranches = writeTranches(dir, "tranches", TRANCHES);

        run(dir, "snp-mode", variants, recal, tranches, "0.0", null, List.of());
        run(dir, "indel-mode", variants, recal, tranches, "0.0", null, List.of("-mode", "INDEL"));
        run(dir, "ignore-all-filters", variants, recal, tranches, "0.0", null,
                List.of("--ignore-all-filters"));
        run(dir, "ignore-named-filter", variants, recal, tranches, "0.0", null,
                List.of("--ignore-filter", "weak"));
        run(dir, "exclude-filtered", variants, recal, tranches, "0.0", null,
                List.of("--exclude-filtered"));
        run(dir, "lod-cutoff", variants, recal, tranches, null, "1.0", List.of());
        run(dir, "ends", ends, endsRecal, tranches, "0.0", null, List.of());
        run(dir, "ends-mismatch", endsMismatch, endsRecal, tranches, "0.0", null, List.of());
        run(dir, "no-recal-record", orphan, recal, tranches, "0.0", null, List.of());
        run(dir, "no-lod", orphan, noLod, tranches, "0.0", null, List.of());
        run(dir, "bad-lod", orphan, badLod, tranches, "0.0", null, List.of());
    }

    static Path writeVcf(final Path dir, final String label, final String header,
                         final String... records) throws Exception {
        final StringBuilder text = new StringBuilder(header);
        for (final String record : records) {
            text.append(record).append("\n");
        }
        final Path file = dir.resolve(label + ".vcf");
        Files.writeString(file, text.toString(), StandardCharsets.UTF_8);
        new IndexFeatureFile().instanceMain(new String[] {"-I", file.toString()});
        System.out.printf("input\t%s\t%s%n", label, ReferenceQueryDump.escape(text.toString()));
        return file;
    }

    static Path writeTranches(final Path dir, final String label, final String text) throws Exception {
        final Path file = dir.resolve(label + ".tranches");
        Files.writeString(file, text, StandardCharsets.UTF_8);
        System.out.printf("input\t%s\t%s%n", label, ReferenceQueryDump.escape(text));
        return file;
    }

    static void run(final Path dir, final String label, final Path variants, final Path recal,
                    final Path tranches, final String level, final String lodCutoff,
                    final List<String> extra) {
        // A name of its own: an output that collided with an input would truncate the file the
        // run is reading.
        final Path output = dir.resolve(label + ".out.vcf");
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
        all.addAll(extra);
        try {
            new ApplyVQSR().instanceMain(all.toArray(new String[0]));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            // The walker wraps what `apply` threw, and the wrapper's message is the locus and the
            // record: the tool's own wording is only in the cause.
            for (Throwable cause = e.getCause(); cause != null; cause = cause.getCause()) {
                System.out.printf("cause\t%s\t%s:%s%n", label, cause.getClass().getName(),
                        ReferenceQueryDump.escape(String.valueOf(cause.getMessage())));
            }
            return;
        }
        try {
            for (final String line : Files.readAllLines(output, StandardCharsets.UTF_8)) {
                if (!line.startsWith("#")) {
                    System.out.printf("vcfline\t%s\t%s%n", label, ReferenceQueryDump.escape(line));
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

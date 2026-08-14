/*
 * ApplyVQSR's allele-specific mode, taken from the reference.
 *
 * With `-AS` every alternate allele is scored on its own against a recal record of its own, and the
 * three per-allele annotations are written as comma-separated lists. Seven behaviours this is built
 * to catch.
 *
 *   - THE PER-ALLELE CLASS TEST IS A LENGTH COMPARISON, not a type: SNP mode keeps an allele whose
 *     length equals the reference's and INDEL mode keeps the rest, so a SPANNING DELETION COUNTS AS
 *     A SNP, as the reference's own comment says;
 *   - AN ALLELE OF THE OTHER MODE IS NOT SKIPPED BUT PADDED, with `NA` for its filter and its
 *     culprit and `NaN` for its LOD, so the three lists always have one entry per alternate allele
 *     whatever mode is running;
 *   - A SPANNING DELETION IS PADDED THE SAME WAY WITHOUT BEING LOOKED UP, so it needs no recal
 *     record even though the class test kept it;
 *   - A MIXED SITE IS LEFT UNFILTERED UNTIL BOTH MODES HAVE RUN. `generateFilterStringFromAlleles`
 *     returns the unfiltered value when neither `bothModesWereRun` nor `onlyOneModeNeeded` holds, so
 *     the FILTER column of a mixed record is `.` while its AS_FilterStatus already names a tranche;
 *   - EVERY RECORD IS EVALUATED IN AS MODE, `evaluateThisVariant` being `useASannotations || ...`,
 *     so a pure indel is annotated by a SNP-mode run rather than passed through;
 *   - THERE IS NO SITE-LEVEL VQSLOD OR CULPRIT AT ALL in this mode: `doAlleleSpecificFiltering`
 *     writes only the three lists, each LOD formatted with `%.4f` rather than through the writer's
 *     double format, while the two training labels are still written at the site;
 *   - AND AN ALLELE WITH NO RECAL RECORD IS A REFUSAL OF ITS OWN, worded for `-AS` and quoting the
 *     whole record, which the walker rethrows as a GATKException naming the locus.
 *
 * Output:
 *
 *     input\t<label>\t<the whole vcf, escaped>
 *     vcfline\t<run>\t<one record line of the output vcf, escaped>
 *     error\t<run>\t<exception class>:<message>
 *     cause\t<run>\t<the wrapped exception's class>:<message>
 *
 * Usage: ApplyVqsrAlleleSpecificDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.vqsr.ApplyVQSR;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class ApplyVqsrAlleleSpecificDump {

    static final String INPUT_HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##contig=<ID=chr1,length=1000>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n";

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
        final Path dir = Path.of("apply-vqsr-allele-specific-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# ApplyVqsrAlleleSpecificDump: one list entry per alternate allele");

        // A multi-allelic SNP, a mixed site, a spanning deletion beside a SNP, and a pure indel.
        final Path variants = writeVcf(dir, "variants", INPUT_HEADER,
                "chr1\t100\t.\tA\tC,G\t50\t.\t.\tGT\t1/2",
                "chr1\t200\t.\tA\tC,ACC\t50\t.\t.\tGT\t1/2",
                "chr1\t300\t.\tA\t*,C\t50\t.\t.\tGT\t1/2",
                "chr1\t400\t.\tACC\tA\t50\t.\t.\tGT\t0/1");

        // One recal record per allele, told apart by its first alternate allele.
        final Path recal = writeVcf(dir, "recal", RECAL_HEADER,
                "chr1\t100\t.\tA\tC\t.\t.\tVQSLOD=5.0000;culprit=QD;POSITIVE_TRAIN_SITE",
                "chr1\t100\t.\tA\tG\t.\t.\tVQSLOD=2.0000;culprit=MQ",
                "chr1\t200\t.\tA\tC\t.\t.\tVQSLOD=-3.0000;culprit=FS",
                "chr1\t200\t.\tA\tACC\t.\t.\tVQSLOD=4.0000;culprit=SOR",
                "chr1\t300\t.\tA\tC\t.\t.\tVQSLOD=2.0000;culprit=MQ",
                "chr1\t400\t.\tACC\tA\t.\t.\tVQSLOD=4.0000;culprit=QD");

        // A record whose only alternate allele has no recal record of its own.
        final Path orphan = writeVcf(dir, "orphan", INPUT_HEADER,
                "chr1\t500\t.\tA\tC\t50\t.\t.\tGT\t0/1");
        final Path orphanRecal = writeVcf(dir, "orphan-recal", RECAL_HEADER,
                // The right locus, the wrong allele.
                "chr1\t500\t.\tA\tG\t.\t.\tVQSLOD=5.0000;culprit=QD");

        final Path tranches = writeTranches(dir, "tranches", TRANCHES);

        run(dir, "as-snp-mode", variants, recal, tranches, "SNP");
        run(dir, "as-indel-mode", variants, recal, tranches, "INDEL");
        run(dir, "as-missing-allele", orphan, orphanRecal, tranches, "SNP");
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
                    final Path tranches, final String mode) {
        // A name of its own: an output that collided with an input would truncate the file the run
        // is reading.
        final Path output = dir.resolve(label + ".out.vcf");
        final List<String> all = new ArrayList<>(List.of(
                "-V", variants.toString(),
                "--recal-file", recal.toString(),
                "--tranches-file", tranches.toString(),
                "-O", output.toString(),
                "--truth-sensitivity-filter-level", "0.0",
                "-AS",
                "-mode", mode));
        try {
            new ApplyVQSR().instanceMain(all.toArray(new String[0]));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
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

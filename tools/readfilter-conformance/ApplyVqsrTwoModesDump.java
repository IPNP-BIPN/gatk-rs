/*
 * ApplyVQSR run twice, taken from the reference.
 *
 * The second run reads the first run's own FILTER lines back out of the header to decide whether the
 * other mode has been applied, and then compares the two modes' filters to pick the most lenient.
 * Six behaviours this is built to catch.
 *
 *   - THE PREVIOUS RUN IS DETECTED BY THREE MAGIC NUMBERS. A FILTER line counts as a tranche when
 *     its ID is at least 12 characters, its first 11 are `VQSRTranche` ignoring case, and its 12th
 *     is `S` or `I`; the interval is then `substring(14)` or `substring(16)`, which is the length of
 *     `VQSRTrancheSNP` and of `VQSRTrancheINDEL`;
 *   - CASE-INSENSITIVE FOR THE PREFIX AND CASE-SENSITIVE FOR THE MODE. `equalsIgnoreCase` accepts
 *     `vqsrtranchesnp...`, and `charAt(11) == 'S'` then rejects it, so a lowercase tranche name
 *     passes the first test and is not a tranche;
 *   - A NAME WHOSE INTERVAL WILL NOT PARSE IS A REFUSAL, but only when it splits into exactly two
 *     pieces on `to`: `trancheIntervalIsValid` returns false for any other count and throws for two
 *     that are not numbers;
 *   - THE SITE FILTER OF THE SECOND RUN IS THE MOST LENIENT OF BOTH MODES, and `most lenient` means
 *     the smallest lower limit pulled out of the filter NAME by a regex rather than any number the
 *     tool holds;
 *   - THAT REGEX IS GREEDY AND LOSES A DIGIT. `VQSRTranche\S+(\d+\.\d+)to(\d+\.\d+)` lets `\S+` eat
 *     as much as it can, so `VQSRTrancheINDEL90.00to99.00` has a lower limit of `0.00` rather than
 *     `90.00` while `VQSRTrancheSNP5.00to90.00` keeps its `5.00`. The second pair of runs is built
 *     on exactly that: the site is filtered with the LESS lenient of the two names, the one whose
 *     true interval starts at 90 rather than the one starting at 5;
 *   - AND A MIXED SITE IS FILTERED ONLY ONCE BOTH MODES HAVE RUN, which is what the pair of runs
 *     shows: the same record is `.` after the first and carries a filter after the second.
 *
 * Output:
 *
 *     input\t<label>\t<the whole vcf, escaped>
 *     vcfline\t<run>\t<one record line of the output vcf, escaped>
 *     filter\t<run>\t<one ##FILTER line of the output vcf, escaped>
 *     error\t<run>\t<exception class>:<message>
 *     cause\t<run>\t<the wrapped exception's class>:<message>
 *
 * Usage: ApplyVqsrTwoModesDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.vqsr.ApplyVQSR;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class ApplyVqsrTwoModesDump {

    static final String INPUT_HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##contig=<ID=chr1,length=1000>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n";

    static final String RECAL_HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=\"the score\">\n"
                    + "##INFO=<ID=culprit,Number=1,Type=String,Description=\"the worst annotation\">\n"
                    + "##contig=<ID=chr1,length=1000>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    static String tranches(final String model) {
        return tranches(model, "VQSRTranche" + model + "90.00to99.00");
    }

    /** The same three tranches with a name of one's own for the middle one. */
    static String tranches(final String model, final String middleName) {
        return "# Variant quality score tranches file\n"
                + "# Version number 5\n"
                + "targetTruthSensitivity,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,model,accessibleTruthSites,callsAtTruthSites,truthSensitivity\n"
                + "90.00,10,5,2.1000,1.9000,3.5000,VQSRTranche" + model + "0.00to90.00," + model + ",100,90,0.9000\n"
                + "99.00,20,9,2.0000,1.8000,1.5000," + middleName + "," + model + ",100,99,0.9900\n"
                + "100.00,30,15,1.9000,1.7000,-0.5000,VQSRTranche" + model + "99.00to100.00," + model + ",100,100,1.0000\n";
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("apply-vqsr-two-modes-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# ApplyVqsrTwoModesDump: the second run reads the first run's header");

        // A mixed site, which is the record that needs both modes, and a plain SNP beside it.
        final Path variants = writeVcf(dir, "variants", INPUT_HEADER,
                "chr1\t100\t.\tA\tC,ACC\t50\t.\t.\tGT\t1/2",
                "chr1\t200\t.\tA\tC\t50\t.\t.\tGT\t0/1");
        final Path recal = writeVcf(dir, "recal", RECAL_HEADER,
                // The SNP allele of the mixed site sits in the middle tranche of the SNP run.
                "chr1\t100\t.\tA\tC\t.\t.\tVQSLOD=2.0000;culprit=QD",
                // Its insertion sits in the widest tranche of the INDEL run.
                "chr1\t100\t.\tA\tACC\t.\t.\tVQSLOD=0.0000;culprit=FS",
                "chr1\t200\t.\tA\tC\t.\t.\tVQSLOD=2.0000;culprit=MQ");

        final Path snpTranches = writeTranches(dir, "tranches-snp", tranches("SNP"));
        final Path indelTranches = writeTranches(dir, "tranches-indel", tranches("INDEL"));

        // A header whose tranche name is spelled in lower case, which the prefix test accepts and
        // the mode test does not.
        final Path lowercase = writeVcf(dir, "lowercase-tranche",
                INPUT_HEADER.replace("##contig=",
                        "##FILTER=<ID=vqsrtranchesnp90.00to99.00,Description=\"lower case\">\n##contig="),
                "chr1\t200\t.\tA\tC\t50\t.\t.\tGT\t0/1");

        // A header whose tranche name splits into two pieces that are not numbers.
        final Path malformed = writeVcf(dir, "malformed-tranche",
                INPUT_HEADER.replace("##contig=",
                        "##FILTER=<ID=VQSRTrancheSNPaatobb,Description=\"two pieces, neither a number\">\n##contig="),
                "chr1\t200\t.\tA\tC\t50\t.\t.\tGT\t0/1");

        // The same two runs again, with two tranche names whose lower limits the greedy regex
        // reorders: 5.00 parses as 5.00 and 90.00 parses as 0.00.
        final Path invertedVariants = writeVcf(dir, "variants-inverted", INPUT_HEADER,
                "chr1\t300\t.\tA\tC,ACC\t50\t.\t.\tGT\t1/2");
        final Path invertedRecal = writeVcf(dir, "recal-inverted", RECAL_HEADER,
                "chr1\t300\t.\tA\tC\t.\t.\tVQSLOD=2.0000;culprit=QD",
                "chr1\t300\t.\tA\tACC\t.\t.\tVQSLOD=2.0000;culprit=FS");
        final Path invertedSnpTranches = writeTranches(dir, "tranches-snp-inverted",
                tranches("SNP", "VQSRTrancheSNP5.00to90.00"));
        final Path invertedIndelTranches = writeTranches(dir, "tranches-indel-inverted",
                tranches("INDEL", "VQSRTrancheINDEL90.00to99.00"));

        final Path first = run(dir, "first-snp", variants, recal, snpTranches, "SNP");
        if (first != null) {
            // The second run reads the first run's output, header and all.
            new IndexFeatureFile().instanceMain(new String[] {"-I", first.toString()});
            run(dir, "second-indel", first, recal, indelTranches, "INDEL");
        }
        final Path firstInverted = run(dir, "first-snp-inverted", invertedVariants, invertedRecal,
                invertedSnpTranches, "SNP");
        if (firstInverted != null) {
            new IndexFeatureFile().instanceMain(new String[] {"-I", firstInverted.toString()});
            run(dir, "second-indel-inverted", firstInverted, invertedRecal, invertedIndelTranches,
                    "INDEL");
        }
        run(dir, "lowercase-tranche", lowercase, recal, indelTranches, "INDEL");
        run(dir, "malformed-tranche", malformed, recal, indelTranches, "INDEL");
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

    static Path run(final Path dir, final String label, final Path variants, final Path recal,
                    final Path tranches, final String mode) {
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
            return null;
        }
        try {
            for (final String line : Files.readAllLines(output, StandardCharsets.UTF_8)) {
                if (line.startsWith("##FILTER=")) {
                    System.out.printf("filter\t%s\t%s%n", label, ReferenceQueryDump.escape(line));
                } else if (!line.startsWith("#")) {
                    System.out.printf("vcfline\t%s\t%s%n", label, ReferenceQueryDump.escape(line));
                }
            }
        } catch (final Exception e) {
            System.out.printf("error\t%s-read\t%s:%s%n", label, e.getClass().getName(),
                    String.valueOf(e.getMessage()));
            return null;
        }
        return output;
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

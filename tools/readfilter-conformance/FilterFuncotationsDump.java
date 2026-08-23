/*
 * FilterFuncotations's output, taken from the reference.
 *
 * A Funcotated VCF read back and marked with the clinical-significance rules that match it. The
 * tool is five filters over the funcotations of each transcript, and what comes out is one INFO
 * field naming the filters that matched and one FILTER value when none did.
 *
 * Twelve behaviours this is built to catch.
 *
 *   - A FILTER MATCHES ONLY WHEN EVERY ONE OF ITS RULES DOES, the rules being reduced with
 *     logicalAnd and a filter with no rules at all reducing to false;
 *   - THE CLINSIG VALUE IS A HashSet OF THE MATCHING NAMES JOINED WITH A COMMA, so its order is
 *     not the order the filters were registered in but the set's own, which over AR, CLINVAR, LOF
 *     and LMM happens to read alphabetically, and the two autosomal recessive filters share one
 *     name and so contribute one entry between them;
 *   - A VARIANT THAT MATCHES NOTHING GETS CLINSIG=NONE AND THE FILTER NOT_CLINSIG, and one that
 *     matches anything is passed rather than left with whatever filters it arrived with;
 *   - CLINVAR'S SIGNIFICANCE TEST IS AN EXACT STRING MATCH against three values, so a ClinVar
 *     value carrying anything else, a qualifier or a second term, does not match;
 *   - AN ExAC SUB-POPULATION WITH AN ALLELE NUMBER OF ZERO IS A MAF OF ZERO rather than a
 *     division, so a variant never seen in ExAC PASSES the frequency rule;
 *   - AND SO IS ONE WHOSE COUNTS DO NOT PARSE: the NumberFormatException is caught and turned into
 *     zero, so a malformed allele count makes the variant pass;
 *   - ONLY SUB-POPULATIONS WHOSE ALLELE COUNT KEY IS PRESENT ARE CONSIDERED, and the maximum over
 *     no sub-population at all is zero;
 *   - THE GNOMAD PATH IS NOT THE SAME CODE: a dataset counts as present when it has no FILTER
 *     funcotation or when its FILTER says PASS, and when neither dataset is present the rule fails
 *     whatever the frequencies say;
 *   - AND A GNOMAD FREQUENCY THAT DOES NOT PARSE IS NOT CAUGHT, unlike ExAC's;
 *   - THE LMM FILTER READS ITS FLAG WITH Boolean.valueOf, which is case-insensitive on "true" and
 *     false for everything else, including "yes" and "1";
 *   - THE HOMVAR RULE ASKS THE VARIANT, NOT THE FUNCOTATIONS, for its hom-var count, and fires only
 *     for the two genes the tool knows;
 *   - AND THE HETVAR RULE IS A SECOND PASS: a gene needs MORE THAN ONE het variant before any of
 *     them is compound het, so a lone het variant in the same gene matches nothing.
 *
 * Output:
 *
 *     input\t<label>=<the whole VCF, escaped>
 *     output\t<label>=<the whole VCF, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: FilterFuncotationsDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.funcotator.FilterFuncotations;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class FilterFuncotationsDump {

    /** The funcotation keys the fixture declares, in the order the values are written. */
    static final List<String> KEYS = List.of(
            "Gencode_19_hugoSymbol",
            "Gencode_19_variantClassification",
            "Gencode_19_annotationTranscript",
            "ACMG_recommendation_Disease_Name",
            "ClinVar_VCF_CLNSIG",
            "ACMGLMMLof_LOF_Mechanism",
            "LMMKnown_LMM_FLAGGED",
            "ExAC_AC_AFR",
            "ExAC_AN_AFR",
            "ExAC_AC_NFE",
            "ExAC_AN_NFE",
            "gnomAD_genome_AF_afr",
            "gnomAD_genome_FILTER",
            "gnomAD_exome_AF_nfe",
            "gnomAD_exome_FILTER");

    /** One transcript's funcotations, in key order, missing values written empty. */
    static String funcotation(final String... values) {
        final String[] fields = new String[KEYS.size()];
        Arrays.fill(fields, "");
        for (int index = 0; index < values.length; index += 2) {
            final int position = KEYS.indexOf(values[index]);
            if (position < 0) {
                throw new IllegalArgumentException("no such key: " + values[index]);
            }
            fields[position] = values[index + 1];
        }
        return "[" + String.join("|", fields) + "]";
    }

    record Variant(int position, String genotype, String funcotation) {}

    static Variant variant(final int position, final String genotype, final String funcotation) {
        return new Variant(position, genotype, funcotation);
    }

    static String vcf(final List<Variant> variants) {
        final StringBuilder text = new StringBuilder();
        text.append("##fileformat=VCFv4.2\n");
        text.append("##contig=<ID=chr1,length=100000>\n");
        text.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        text.append("##INFO=<ID=FUNCOTATION,Number=A,Type=String,Description=\"Funcotation fields are: ")
                .append(String.join("|", KEYS)).append("\">\n");
        text.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n");
        for (final Variant variant : variants) {
            text.append("chr1\t").append(variant.position()).append("\t.\tA\tC\t100\t.\tFUNCOTATION=")
                    .append(variant.funcotation()).append("\tGT\t").append(variant.genotype())
                    .append('\n');
        }
        return text.toString();
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("filter-funcotations-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# FilterFuncotationsDump: a Funcotated VCF marked with the rules that match it");

        // Every ExAC-side rule, one variant each.
        run(dir, "exac", "hg19", "exac", List.of(
                // ClinVar: the gene is on the ACMG list, the significance is one of the three
                // matching values and the frequency is under five per cent.
                variant(100, "0/1", funcotation(
                        "ACMG_recommendation_Disease_Name", "Cardiomyopathy",
                        "ClinVar_VCF_CLNSIG", "Pathogenic",
                        "ExAC_AC_AFR", "1", "ExAC_AN_AFR", "1000")),
                // The same, with the third matching value.
                variant(200, "0/1", funcotation(
                        "ACMG_recommendation_Disease_Name", "Cardiomyopathy",
                        "ClinVar_VCF_CLNSIG", "Pathogenic/Likely_pathogenic",
                        "ExAC_AC_AFR", "1", "ExAC_AN_AFR", "1000")),
                // A significance that carries a qualifier, which is not one of the three.
                variant(300, "0/1", funcotation(
                        "ACMG_recommendation_Disease_Name", "Cardiomyopathy",
                        "ClinVar_VCF_CLNSIG", "Pathogenic_low_penetrance",
                        "ExAC_AC_AFR", "1", "ExAC_AN_AFR", "1000")),
                // A frequency above the ClinVar threshold, taken as the maximum over two
                // sub-populations.
                variant(400, "0/1", funcotation(
                        "ACMG_recommendation_Disease_Name", "Cardiomyopathy",
                        "ClinVar_VCF_CLNSIG", "Pathogenic",
                        "ExAC_AC_AFR", "1", "ExAC_AN_AFR", "1000",
                        "ExAC_AC_NFE", "60", "ExAC_AN_NFE", "1000")),
                // An allele number of zero, which is a frequency of zero and not a division.
                variant(500, "0/1", funcotation(
                        "ACMG_recommendation_Disease_Name", "Cardiomyopathy",
                        "ClinVar_VCF_CLNSIG", "Pathogenic",
                        "ExAC_AC_AFR", "5", "ExAC_AN_AFR", "0")),
                // An allele count that does not parse, which is caught and read as zero.
                variant(600, "0/1", funcotation(
                        "ACMG_recommendation_Disease_Name", "Cardiomyopathy",
                        "ClinVar_VCF_CLNSIG", "Pathogenic",
                        "ExAC_AC_AFR", "many", "ExAC_AN_AFR", "1000")),
                // An allele count with no allele number beside it, which is the same zero.
                variant(700, "0/1", funcotation(
                        "ACMG_recommendation_Disease_Name", "Cardiomyopathy",
                        "ClinVar_VCF_CLNSIG", "Pathogenic",
                        "ExAC_AC_AFR", "5")),
                // No ACMG gene at all, so the first ClinVar rule fails.
                variant(800, "0/1", funcotation(
                        "ClinVar_VCF_CLNSIG", "Pathogenic",
                        "ExAC_AC_AFR", "1", "ExAC_AN_AFR", "1000")),
                // Loss of function: a matching classification, the mechanism and a frequency under
                // one per cent.
                variant(900, "0/1", funcotation(
                        "Gencode_19_variantClassification", "FRAME_SHIFT_DEL",
                        "ACMGLMMLof_LOF_Mechanism", "YES",
                        "ExAC_AC_AFR", "1", "ExAC_AN_AFR", "1000")),
                // A classification that is not one of the five.
                variant(1000, "0/1", funcotation(
                        "Gencode_19_variantClassification", "MISSENSE",
                        "ACMGLMMLof_LOF_Mechanism", "YES",
                        "ExAC_AC_AFR", "1", "ExAC_AN_AFR", "1000")),
                // A frequency between the two thresholds, which passes ClinVar and fails LOF.
                variant(1100, "0/1", funcotation(
                        "ACMG_recommendation_Disease_Name", "Cardiomyopathy",
                        "ClinVar_VCF_CLNSIG", "Pathogenic",
                        "Gencode_19_variantClassification", "NONSENSE",
                        "ACMGLMMLof_LOF_Mechanism", "YES",
                        "ExAC_AC_AFR", "30", "ExAC_AN_AFR", "1000")),
                // LMM, whose flag is read with Boolean.valueOf.
                variant(1200, "0/1", funcotation("LMMKnown_LMM_FLAGGED", "true")),
                variant(1300, "0/1", funcotation("LMMKnown_LMM_FLAGGED", "TRUE")),
                variant(1400, "0/1", funcotation("LMMKnown_LMM_FLAGGED", "yes")),
                variant(1500, "0/1", funcotation("LMMKnown_LMM_FLAGGED", "1")),
                // A hom-var in one of the two genes the tool knows, and a het in the same gene,
                // which is alone and therefore not compound.
                variant(1600, "1/1", funcotation("Gencode_19_hugoSymbol", "ATP7B")),
                variant(1700, "0/1", funcotation("Gencode_19_hugoSymbol", "ATP7B")),
                // A hom-var in a gene the tool does not know.
                variant(1800, "1/1", funcotation("Gencode_19_hugoSymbol", "BRCA1")),
                // Two hets in the other known gene, which makes both of them compound.
                variant(1900, "0/1", funcotation("Gencode_19_hugoSymbol", "MUTYH")),
                variant(2000, "0/1", funcotation("Gencode_19_hugoSymbol", "MUTYH")),
                // Everything at once, to show what the joined CLINSIG value looks like.
                variant(2100, "1/1", funcotation(
                        "Gencode_19_hugoSymbol", "ATP7B",
                        "ACMG_recommendation_Disease_Name", "Wilson",
                        "ClinVar_VCF_CLNSIG", "Likely_pathogenic",
                        "Gencode_19_variantClassification", "SPLICE_SITE",
                        "ACMGLMMLof_LOF_Mechanism", "YES",
                        "LMMKnown_LMM_FLAGGED", "true",
                        "ExAC_AC_AFR", "1", "ExAC_AN_AFR", "1000")),
                // Nothing at all.
                variant(2200, "0/1", funcotation("Gencode_19_hugoSymbol", "BRCA1"))));

        // The gnomAD path, which is different code with different rules.
        run(dir, "gnomad", "hg19", "gnomad", List.of(
                // No FILTER funcotation at all, so both datasets count as present.
                variant(100, "0/1", funcotation(
                        "ACMG_recommendation_Disease_Name", "Cardiomyopathy",
                        "ClinVar_VCF_CLNSIG", "Pathogenic",
                        "gnomAD_genome_AF_afr", "0.01")),
                // A FILTER that says PASS, which also counts as present.
                variant(200, "0/1", funcotation(
                        "ACMG_recommendation_Disease_Name", "Cardiomyopathy",
                        "ClinVar_VCF_CLNSIG", "Pathogenic",
                        "gnomAD_genome_AF_afr", "0.01",
                        "gnomAD_genome_FILTER", "PASS")),
                // A frequency over the threshold.
                variant(300, "0/1", funcotation(
                        "ACMG_recommendation_Disease_Name", "Cardiomyopathy",
                        "ClinVar_VCF_CLNSIG", "Pathogenic",
                        "gnomAD_genome_AF_afr", "0.06",
                        "gnomAD_genome_FILTER", "PASS")),
                // Both datasets filtered out, which fails the rule whatever the frequency says.
                variant(400, "0/1", funcotation(
                        "ACMG_recommendation_Disease_Name", "Cardiomyopathy",
                        "ClinVar_VCF_CLNSIG", "Pathogenic",
                        "gnomAD_genome_AF_afr", "0.01",
                        "gnomAD_genome_FILTER", "RF",
                        "gnomAD_exome_FILTER", "AC0")),
                // One filtered out and one not.
                variant(500, "0/1", funcotation(
                        "ACMG_recommendation_Disease_Name", "Cardiomyopathy",
                        "ClinVar_VCF_CLNSIG", "Pathogenic",
                        "gnomAD_genome_AF_afr", "0.01",
                        "gnomAD_genome_FILTER", "RF",
                        "gnomAD_exome_AF_nfe", "0.02",
                        "gnomAD_exome_FILTER", "PASS"))));

        // A gnomAD frequency that does not parse, which nothing catches, so it takes the whole
        // run down and has to stand on its own.
        run(dir, "gnomad-unparseable", "hg19", "gnomad", List.of(
                variant(100, "0/1", funcotation(
                        "ACMG_recommendation_Disease_Name", "Cardiomyopathy",
                        "ClinVar_VCF_CLNSIG", "Pathogenic",
                        "gnomAD_genome_AF_afr", "many",
                        "gnomAD_genome_FILTER", "PASS"))));

        // The other reference versions, whose gencode number is part of every key the tool looks
        // for, so the same funcotations stop matching.
        run(dir, "hg38-keys", "hg38", "exac", List.of(
                variant(100, "0/1", funcotation(
                        "Gencode_19_variantClassification", "NONSENSE",
                        "ACMGLMMLof_LOF_Mechanism", "YES",
                        "ExAC_AC_AFR", "1", "ExAC_AN_AFR", "1000"))));

        // A VCF with no FUNCOTATION header line at all.
        final Path bare = dir.resolve("bare.vcf");
        Files.writeString(bare, "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                + "chr1\t100\t.\tA\tC\t100\t.\t.\n", StandardCharsets.UTF_8);
        execute(dir, "no-funcotation-header", bare, "hg19", "exac");
    }

    static void run(final Path dir, final String label, final String reference,
                    final String dataSource, final List<Variant> variants) throws Exception {
        final Path input = dir.resolve(label + ".vcf");
        Files.writeString(input, vcf(variants), StandardCharsets.UTF_8);
        execute(dir, label, input, reference, dataSource);
    }

    static void execute(final Path dir, final String label, final Path input,
                        final String reference, final String dataSource) throws Exception {
        System.out.printf("input\t%s=%s%n", label,
                ReferenceQueryDump.escape(Files.readString(input)));
        try {
            new IndexFeatureFile().instanceMain(new String[] {"-I", input.toString()});
        } catch (final Exception e) {
            System.out.printf("index\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
        }
        final Path out = dir.resolve(label + ".filtered.vcf");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-V", input.toString(), "-O", out.toString(),
                "--ref-version", reference,
                "--allele-frequency-data-source", dataSource));
        try {
            new FilterFuncotations().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            return;
        }
        if (Files.exists(out)) {
            System.out.printf("output\t%s=%s%n", label,
                    ReferenceQueryDump.escape(masked(Files.readString(out), dir)));
        }
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}

/*
 * Funcotator's annotated VCF, taken from the reference.
 *
 * A VCF annotated from a folder of data sources. The GENCODE source is the one FuncotateSegments'
 * fixture already builds, reused whole: what is measured here is what changes when the thing
 * being annotated is a variant rather than a segment.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - THE ANNOTATIONS LAND IN ONE INFO FIELD, `FUNCOTATION`, whose value is the source's fields
 *     joined by `|` and wrapped in brackets, in an order the header line declares;
 *   - EVERY ALTERNATE GETS ITS OWN BRACKETED FUNCOTATION, joined by a comma, so a two-alternate
 *     site carries two and a reader who splits on the comma alone splits inside neither;
 *   - THE VARIANT CLASSIFICATION IS COMPUTED FROM THE TRANSCRIPT: the fixture's one exon runs
 *     from 1000 to 1200 and its CDS with it, so 1050 is SILENT, 1060 and 1070 are MISSENSE, 1500
 *     is INTRON and 5000 is IGR;
 *   - A VARIANT OUTSIDE EVERY GENE STILL GETS A FUNCOTATION, whose gene name is the literal
 *     `Unknown` rather than an empty field;
 *   - THE CODON CHANGE AND THE PROTEIN CHANGE ARE WRITTEN OUT, so the two alternates of one site
 *     differ in `p.A21T` and `p.A21S`;
 *   - --annotation-override REPLACES A FIELD THE SOURCES PRODUCED, by its prefixed name;
 *   - --annotation-default ADDS ONE with a fixed value;
 *   - --remove-filtered-variants DROPS A FILTERED RECORD rather than annotating it, taking the
 *     file from five records to four;
 *   - THE MAF OUTPUT IS A DIFFERENT FILE ENTIRELY, opening `#version 2.4` and carrying one row
 *     per alternate with no VCF header at all;
 *   - A REFERENCE VERSION WITH NO DIRECTORY IS REFUSED, the same way it is for segments;
 *   - AND THE LOCATABLE XSV SOURCE THAT FuncotateSegments READS STRAIGHT THROUGH HAS TO BE
 *     INDEXED FOR FUNCOTATOR, which queries it by interval: the folder here holds GENCODE alone
 *     for that reason.
 *
 * Output:
 *
 *     vcf\t<label>=<that vcf, escaped>
 *     out\t<label>=<the whole output, escaped>
 *     header\t<label>\t<the `##INFO=<ID=FUNCOTATION` line>
 *     none\t<label>=<what was not written>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: FuncotatorDump
 */

import org.broadinstitute.hellbender.tools.funcotator.Funcotator;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class FuncotatorDump {

    static List<String> header() {
        return new ArrayList<>(List.of(
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=" + FuncotateSegmentsDump.CONTIG_LENGTH + ">",
                "##FILTER=<ID=LOW,Description=\"Low\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample"));
    }

    static String site(final int position, final String reference, final String alternate,
                       final String filter) {
        return "chr1\t" + position + "\t.\t" + reference + "\t" + alternate
                + "\t100.00\t" + filter + "\t.\tGT\t0/1";
    }

    static String vcf(final List<String> sites) {
        final List<String> lines = header();
        lines.addAll(sites);
        lines.add("");
        return String.join("\n", lines);
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("funcotator-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# FuncotatorDump: a VCF annotated from the same folder of data "
                + "sources FuncotateSegments uses");

        final Path fasta = FuncotateSegmentsDump.writeReference(dir);
        // Only the GENCODE source. FuncotateSegments' folder also holds a locatable XSV, which
        // Funcotator QUERIES BY INTERVAL and so wants an index beside: the segment tool reads it
        // straight through and never asks for one.
        final Path good = gencodeOnly(dir, "ds-good", "hg38");
        final Path hg19Only = gencodeOnly(dir, "ds-hg19", "hg19");

        // The fixture's gene sits at chr1:1000-2000 with its one exon at 1000-1200 and its CDS
        // over the same span, so a substitution at 1050 is inside the coding exon and one at
        // 1500 is inside the gene but outside it.
        // The records have to be in position order or the indexer refuses the file before the
        // tool ever sees it.
        final String variants = vcf(List.of(
                site(1050, "C", "A", "PASS"),
                site(1060, "G", "A,T", "PASS"),
                site(1070, "T", "A", "LOW"),
                site(1500, "T", "A", "PASS"),
                site(5000, "A", "G", "PASS")));
        final Path vcfPath = index(FuncotateSegmentsDump.write(dir, "variants.vcf", variants));
        System.out.printf("vcf\tvariants=%s%n", ReferenceQueryDump.escape(variants));

        run(dir, "annotated", vcfPath, fasta, good, "hg38", List.of());
        // The MAF output, which is a different file entirely.
        run(dir, "maf", vcfPath, fasta, good, "hg38",
                List.of("--output-file-format", "MAF"));
        // A default annotation and an override.
        run(dir, "annotation-default", vcfPath, fasta, good, "hg38",
                List.of("--annotation-default", "NOTE:added"));
        run(dir, "annotation-override", vcfPath, fasta, good, "hg38",
                List.of("--annotation-override", "Gencode_1_hugoSymbol:OVERRIDDEN"));
        // The filtered record dropped.
        run(dir, "remove-filtered", vcfPath, fasta, good, "hg38",
                List.of("--remove-filtered-variants", "true"));
        // The transcript selection mode.
        run(dir, "canonical-transcripts", vcfPath, fasta, good, "hg38",
                List.of("--transcript-selection-mode", "CANONICAL"));
        run(dir, "all-transcripts", vcfPath, fasta, good, "hg38",
                List.of("--transcript-selection-mode", "ALL"));
        // A reference version the folder has no directory for.
        run(dir, "wrong-ref-version", vcfPath, fasta, hg19Only, "hg38", List.of());
    }

    /** A data-source root holding the manifest and the GENCODE source alone. */
    static Path gencodeOnly(final Path dir, final String name, final String refVersion)
            throws Exception {
        final Path root = dir.resolve(name);
        Files.createDirectories(root);
        Files.writeString(root.resolve("MANIFEST.txt"),
                "Version: 1.7.hg38.20220101\nSource: test\nAlternate Source: test\n",
                StandardCharsets.UTF_8);
        FuncotateSegmentsDump.writeGencode(root, refVersion);
        return root;
    }

    static Path index(final Path path) throws Exception {
        htsjdk.tribble.index.IndexFactory.createLinearIndex(path.toFile(),
                new htsjdk.variant.vcf.VCFCodec()).writeBasedOnFeatureFile(path.toFile());
        return path;
    }

    static void run(final Path dir, final String label, final Path input, final Path fasta,
                    final Path dataSources, final String refVersion, final List<String> extra)
            throws Exception {
        final boolean maf = extra.contains("MAF");
        final Path out = dir.resolve("out-" + label + (maf ? ".maf" : ".vcf"));
        final List<String> argv = new ArrayList<>(List.of(
                "-V", input.toString(),
                "-O", out.toString(),
                "-R", fasta.toString(),
                "--data-sources-path", dataSources.toString(),
                "--ref-version", refVersion));
        if (!maf) {
            argv.addAll(List.of("--output-file-format", "VCF"));
        }
        argv.addAll(extra);
        try {
            new Funcotator().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(cause.getMessage()), dir)));
            return;
        }
        if (!Files.exists(out)) {
            System.out.printf("none\t%s=no output%n", label);
            return;
        }
        final StringBuilder body = new StringBuilder();
        for (final String line : Files.readString(out).split("\n", -1)) {
            if (line.isEmpty()) {
                continue;
            }
            if (line.startsWith("##INFO=<ID=FUNCOTATION")) {
                // The header line names every field in order, so it is its own row.
                System.out.printf("header\t%s\t%s%n", label, masked(line, dir));
            } else if (!line.startsWith("##")) {
                body.append(line).append("\n");
            }
        }
        System.out.printf("out\t%s=%s%n", label,
                ReferenceQueryDump.escape(masked(body.toString(), dir)));
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}

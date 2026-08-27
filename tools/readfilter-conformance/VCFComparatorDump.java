/*
 * VCFComparator's complaints, taken from the reference.
 *
 * What counts as two VCFs disagreeing. The tool's output IS its exceptions: it walks two files in
 * step and throws on the first difference it is not told to tolerate, so what is measurable is
 * which differences it notices, what it says about each, and which argument silences it.
 *
 * Ten behaviours this is built to catch.
 *
 *   - THE EXPECTED FILE IS IDENTIFIED BY ITS TAG, not by its order, so both inputs must be tagged
 *     and exactly one of them `expected`;
 *   - A VARIANT IN ONE FILE AND NOT THE OTHER IS ONLY A COMPLAINT WHEN A GENOTYPE HAS QUALITY
 *     ZERO: a confidently called variant present on one side alone passes in silence, and
 *     --ignore-gq0 silences the other case too;
 *   - AN EXTRA ALLELE CHANGES `AC` WITH IT, so the ATTRIBUTE complaint is reached first;
 *   - AND WITH `AC` IGNORED THE EXTRA ALLELE IS STILL NOT CHECKED, because the guard in front of
 *     the allele comparison is `actualHasNewAlleles(expected, actual)` with its arguments the
 *     other way round: it asks whether EXPECTED has an allele actual lacks. An allele ADDED to
 *     actual therefore passes in silence, and --allow-extra-alleles is unreachable through it;
 *   - EVERY DIFFERENCE IS REPORTED WITH ITS POSITION, wrapped around the message rather than
 *     inside it;
 *   - THE QUAL COMPARISON IS A TOLERANCE, and its message names both the difference and the
 *     tolerance it exceeded;
 *   - DIFFERENT FILTERS AND UNAPPLIED FILTERS ARE TWO DIFFERENT COMPLAINTS;
 *   - AN EXTRA ALLELE IS A MISMATCH UNLESS IT IS ALLOWED, and allowing it is a different argument
 *     from allowing a new star;
 *   - AN INFO ATTRIBUTE DIFFERENCE IS ONE COMPLAINT FOR THE WHOLE RECORD, whatever attribute it
 *     was, and --ignore-attribute takes a key at a time;
 *   - --positions-only IS MUTUALLY EXCLUSIVE WITH ALMOST EVERY OTHER TOLERANCE, because it already
 *     ignores what they tolerate;
 *   - --warn-on-errors TURNS EVERY COMPLAINT INTO A WARNING and the run then succeeds;
 *   - AND A SINGLE-SAMPLE PAIR IS EXPECTED TO BE GVCFs, so two ordinary single-sample VCFs are
 *     refused for having no <NON_REF>.
 *
 * Output:
 *
 *     vcf\t<label>=<that vcf, escaped>
 *     ok\t<label>=<the run succeeded>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: VCFComparatorDump
 */

import org.broadinstitute.hellbender.tools.walkers.variantutils.VCFComparator;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class VCFComparatorDump {

    static final int CONTIG_LENGTH = 199980;

    static List<String> header() {
        return new ArrayList<>(List.of(
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=" + CONTIG_LENGTH + ">",
                "##FILTER=<ID=LOW,Description=\"Low\">",
                "##FILTER=<ID=OTHER,Description=\"Other\">",
                "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">",
                "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2"));
    }

    static String site(final int position, final String id, final String reference,
                       final String alternate, final String qual, final String filter,
                       final String info, final String g1, final String g2) {
        return "chr1\t" + position + "\t" + id + "\t" + reference + "\t" + alternate + "\t" + qual
                + "\t" + filter + "\t" + info + "\tGT:GQ\t" + g1 + "\t" + g2;
    }

    static String vcf(final List<String> sites) {
        final List<String> lines = header();
        lines.addAll(sites);
        lines.add("");
        return String.join("\n", lines);
    }

    /**
     * The baseline every `actual` is a single change away from.
     *
     * The site at 4000 carries TWO alternates, so an actual can drop one: that is the only
     * direction the allele check is reachable from.
     */
    static List<String> baseline() {
        return new ArrayList<>(List.of(
                site(1000, ".", "A", "C", "100.00", "PASS", "AC=1;DP=20", "0/1:50", "0/0:50"),
                site(2000, ".", "G", "T", "200.00", "PASS", "AC=2;DP=30", "0/1:60", "0/1:60"),
                site(4000, ".", "C", "A,T", "400.00", "PASS", "AC=1,1;DP=50", "0/1:80",
                        "0/2:80")));
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("vcf-comparator-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# VCFComparatorDump: what counts as two VCFs disagreeing");

        final Path fasta = writeReference(dir);
        final String expected = vcf(baseline());
        final Path expectedPath = write(dir, "expected.vcf", expected);
        System.out.printf("vcf\texpected=%s%n", ReferenceQueryDump.escape(expected));

        // Identical.
        final Path same = writeActual(dir, "same", baseline());
        run(dir, "identical", expectedPath, same, List.of());

        // A variant only in actual, whose genotypes are confidently called. The unmatched-variant
        // complaint is guarded on a genotype quality of ZERO, so this one passes.
        final List<String> extraSite = baseline();
        // INSERTED before the site at 4000, not appended: the driving iterator refuses a VCF whose
        // records are not in position order.
        extraSite.add(2, site(3000, ".", "T", "A", "300.00", "PASS", "AC=1;DP=40", "0/1:70",
                "0/0:70"));
        run(dir, "extra-variant", expectedPath, writeActual(dir, "extra", extraSite), List.of());
        // The same variant with a genotype quality of zero, which is what the complaint is looking
        // for.
        final List<String> extraGq0 = baseline();
        extraGq0.add(2, site(3000, ".", "T", "A", "300.00", "PASS", "AC=1;DP=40", "0/1:0",
                "0/0:0"));
        final Path extraGq0Path = writeActual(dir, "extra-gq0", extraGq0);
        run(dir, "extra-variant-gq0", expectedPath, extraGq0Path, List.of());
        run(dir, "extra-variant-gq0-ignored", expectedPath, extraGq0Path,
                List.of("--ignore-gq0", "true"));

        // A different QUAL, with and without the tolerance that allows it.
        final List<String> qual = baseline();
        qual.set(0, site(1000, ".", "A", "C", "105.00", "PASS", "AC=1;DP=20", "0/1:50", "0/0:50"));
        final Path qualPath = writeActual(dir, "qual", qual);
        run(dir, "qual-differs", expectedPath, qualPath, List.of());
        run(dir, "qual-tolerated", expectedPath, qualPath,
                List.of("--qual-change-allowed", "10"));
        run(dir, "qual-ignored", expectedPath, qualPath, List.of("--ignore-quals", "true"));

        // A different filter, and one applied on only one side.
        final List<String> filtered = baseline();
        filtered.set(0, site(1000, ".", "A", "C", "100.00", "LOW", "AC=1;DP=20", "0/1:50",
                "0/0:50"));
        final Path filteredPath = writeActual(dir, "filtered", filtered);
        run(dir, "filters-differ", expectedPath, filteredPath, List.of());
        run(dir, "filters-ignored", expectedPath, filteredPath, List.of("--ignore-filters", "true"));
        final List<String> unfiltered = baseline();
        unfiltered.set(0, site(1000, ".", "A", "C", "100.00", ".", "AC=1;DP=20", "0/1:50",
                "0/0:50"));
        run(dir, "filters-unapplied", expectedPath, writeActual(dir, "unfiltered", unfiltered),
                List.of());

        // An extra allele, and the two arguments that allow one.
        final List<String> extraAllele = baseline();
        extraAllele.set(0, site(1000, ".", "A", "C,G", "100.00", "PASS", "AC=1,0;DP=20", "0/1:50",
                "0/0:50"));
        final Path extraAllelePath = writeActual(dir, "extra-allele", extraAllele);
        // AC is Number=A, so an extra allele changes it too and the ATTRIBUTE complaint is reached
        // first: the allele check only runs once AC is ignored. That ordering is the measurement.
        run(dir, "alleles-differ-hits-ac", expectedPath, extraAllelePath, List.of());
        run(dir, "alleles-differ", expectedPath, extraAllelePath,
                List.of("--ignore-attribute", "AC"));
        run(dir, "alleles-allowed", expectedPath, extraAllelePath,
                List.of("--ignore-attribute", "AC", "--allow-extra-alleles", "true"));
        // An allele MISSING from actual, which is the direction the guard actually reaches. The
        // guard is `actualHasNewAlleles(expected, actual)`, its arguments the other way round, so
        // it asks whether EXPECTED has an allele actual lacks.
        final List<String> missingAllele = baseline();
        missingAllele.set(2, site(4000, ".", "C", "A", "400.00", "PASS", "AC=1;DP=50", "0/1:80",
                "0/0:80"));
        final Path missingAllelePath = writeActual(dir, "missing-allele", missingAllele);
        run(dir, "allele-missing", expectedPath, missingAllelePath,
                List.of("--ignore-attribute", "AC"));
        run(dir, "allele-missing-allowed", expectedPath, missingAllelePath,
                List.of("--ignore-attribute", "AC", "--allow-extra-alleles", "true"));

        // A different INFO attribute, and the key-at-a-time argument that ignores one.
        final List<String> info = baseline();
        info.set(0, site(1000, ".", "A", "C", "100.00", "PASS", "AC=1;DP=25", "0/1:50", "0/0:50"));
        final Path infoPath = writeActual(dir, "info", info);
        run(dir, "info-differs", expectedPath, infoPath, List.of());
        run(dir, "info-ignored-key", expectedPath, infoPath, List.of("--ignore-attribute", "DP"));
        run(dir, "info-ignored-wrong-key", expectedPath, infoPath,
                List.of("--ignore-attribute", "AC"));
        run(dir, "info-ignored-all", expectedPath, infoPath, List.of("--ignore-annotations", "true"));

        // A different dbSNP id.
        final List<String> ids = baseline();
        ids.set(0, site(1000, "rs1", "A", "C", "100.00", "PASS", "AC=1;DP=20", "0/1:50", "0/0:50"));
        run(dir, "ids-differ", expectedPath, writeActual(dir, "ids", ids), List.of());

        // Positions only, which ignores every one of the above.
        run(dir, "positions-only-alleles", expectedPath, extraAllelePath,
                List.of("--positions-only", "true"));
        // And is refused beside the arguments it subsumes.
        run(dir, "positions-only-with-quals", expectedPath, qualPath,
                List.of("--positions-only", "true", "--ignore-quals", "true"));

        // Warn instead of throw.
        run(dir, "warn-on-errors", expectedPath, qualPath, List.of("--warn-on-errors", "true"));

        // One input, and two inputs neither of which is tagged `expected`.
        runRaw(dir, "one-input", List.of("-V:expected", expectedPath.toString(),
                "-R", fasta.toString()));
        runRaw(dir, "no-expected", List.of(
                "-V:first", expectedPath.toString(),
                "-V:second", qualPath.toString(),
                "-R", fasta.toString()));
    }

    static Path writeActual(final Path dir, final String name, final List<String> sites)
            throws Exception {
        final String text = vcf(sites);
        System.out.printf("vcf\t%s=%s%n", name, ReferenceQueryDump.escape(text));
        return write(dir, name + ".vcf", text);
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static void run(final Path dir, final String label, final Path expected, final Path actual,
                    final List<String> extra) throws Exception {
        final List<String> argv = new ArrayList<>(List.of(
                "-V:expected", expected.toString(),
                "-V:actual", actual.toString(),
                "-R", dir.resolve("reference.fasta").toString()));
        argv.addAll(extra);
        runRaw(dir, label, argv);
    }

    /** The tool requires a reference, whatever it compares. */
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

    static void runRaw(final Path dir, final String label, final List<String> argv) {
        try {
            new VCFComparator().instanceMain(argv.toArray(new String[0]));
            System.out.printf("ok\t%s=succeeded%n", label);
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(cause.getMessage()), dir)));
        }
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}

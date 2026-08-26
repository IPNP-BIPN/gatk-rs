/*
 * SVConcordance's eval-versus-truth matching, taken from the reference.
 *
 * Which truth variant an evaluation variant is judged to be, and what that judgement writes on it.
 * The match is not the clustering of SVCluster: it is a single closest truth record per eval
 * record, chosen by a chain of tiebreakers, under a linkage whose one override is asymmetric.
 *
 * Nine behaviours this is built to catch.
 *
 *   - EACH EVAL RECORD TAKES ONE TRUTH RECORD, not a cluster, and the one it takes is the closest
 *     by TOTAL breakend distance;
 *   - THE TIEBREAKER IS THE CLOSEST SINGLE BREAKEND: two truth records at the same total distance
 *     are separated by the smaller of their two ends, so the one that is exact on one side wins
 *     over the one that is off by the same amount on both;
 *   - THE CNV OVERRIDE IS ASYMMETRIC: the comment says CNV/DEL and CNV/DUP matching is not allowed,
 *     and `(aType == CNV || bType != CNV) && aType != bType` only refuses it when the CNV is the
 *     EVAL record. An eval DEL against a truth CNV falls through to the base linkage and matches;
 *   - AN UNMATCHED EVAL RECORD IS STILL WRITTEN, as a false positive with its truth allele counts
 *     set to nothing;
 *   - A MULTIALLELIC CNV IS SCORED ON COPY STATE INSTEAD, per sample as TRUTH_CN_EQUAL and over the
 *     record as CNV_CONCORDANCE, and gets no genotype concordance at all;
 *   - CONCORDANCE IS COMPUTED ON THE COMMON SAMPLES ONLY, so a sample present in one VCF and not
 *     the other is annotated on neither side;
 *   - A TRUTH RECORD THAT CARRIES AC/AF/AN HAS THEM COPIED, and one that does not has them
 *     RECOUNTED FROM THE EVAL RECORD'S ALT ALLELES over the truth genotypes;
 *   - --do-not-sort CHANGES NOTHING FOR AN IN-ORDER INPUT: the output buffer it removes sorts by
 *     position, and the records complete in the order they were added, so the flag only drops the
 *     index. It is measured here as the control that shows that;
 *   - AND THE SEQUENCE DICTIONARY IS REQUIRED, from --sequence-dictionary rather than from either
 *     VCF, because the VCFs' own dictionaries are sometimes out of order.
 *
 * Output:
 *
 *     vcf\teval=<the whole eval vcf, escaped>
 *     vcf\ttruth=<the whole truth vcf, escaped>
 *     out\t<label>=<the whole output vcf without its header, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: SVConcordanceDump
 */

import org.broadinstitute.hellbender.tools.walkers.sv.SVConcordance;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class SVConcordanceDump {

    /** The eval samples. s1 is eval-only, so it is outside the common set. */
    static final List<String> EVAL_SAMPLES = List.of("s1", "s2", "s3");
    /** The truth samples. s4 is truth-only, so it is outside the common set too. */
    static final List<String> TRUTH_SAMPLES = List.of("s2", "s3", "s4");

    static List<String> header(final List<String> samples) {
        return new ArrayList<>(List.of(
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=199980>",
                "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">",
                "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length\">",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">",
                "##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Algorithms\">",
                "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">",
                "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">",
                "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=ECN,Number=1,Type=Integer,Description=\"Expected copy number\">",
                "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number\">",
                "##ALT=<ID=DEL,Description=\"Deletion\">",
                "##ALT=<ID=DUP,Description=\"Duplication\">",
                "##ALT=<ID=CNV,Description=\"Copy number variant\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                        + String.join("\t", samples)));
    }

    /** A record whose genotypes are GT:ECN. */
    static String record(final String id, final int start, final String type, final int end,
                         final String info, final String... genotypes) {
        final StringBuilder attributes = new StringBuilder("SVTYPE=" + type + ";END=" + end
                + ";SVLEN=" + (end - start + 1) + ";ALGORITHMS=depth");
        if (info != null) {
            attributes.append(';').append(info);
        }
        final List<String> withEcn = new ArrayList<>();
        for (final String genotype : genotypes) {
            withEcn.add(genotype + ":2");
        }
        return "chr1\t" + start + "\t" + id + "\tN\t<" + type + ">\t.\t.\t" + attributes
                + "\tGT:ECN\t" + String.join("\t", withEcn);
    }

    /** A multiallelic CNV, whose genotypes carry a copy number instead of an informative GT. */
    static String cnv(final String id, final int start, final int end, final int... copyNumbers) {
        final List<String> genotypes = new ArrayList<>();
        for (final int copyNumber : copyNumbers) {
            genotypes.add("./.:2:" + copyNumber);
        }
        return "chr1\t" + start + "\t" + id + "\tN\t<CNV>\t.\t.\tSVTYPE=CNV;END=" + end
                + ";SVLEN=" + (end - start + 1) + ";ALGORITHMS=depth\tGT:ECN:CN\t"
                + String.join("\t", genotypes);
    }

    static String buildEval() {
        final List<String> lines = header(EVAL_SAMPLES);
        // An exact match, genotype for genotype, over the common samples.
        lines.add(record("exact", 1000, "DEL", 2000, null, "0/1", "0/1", "0/0"));
        // The same locus judged against a truth record whose genotypes disagree on both common
        // samples, which is what moves the concordance off 1.
        lines.add(record("discordant", 20000, "DEL", 21000, null, "0/1", "0/1", "1/1"));
        // Nothing within reach in the truth VCF.
        lines.add(record("nomatch", 40000, "DEL", 41000, null, "0/1", "0/0", "0/0"));
        // Two truth records at the SAME total distance, separated only by the closer breakend.
        lines.add(record("tie", 60000, "DEL", 61000, null, "0/1", "0/1", "0/0"));
        // An eval CNV against a truth DEL: refused, because the CNV is the eval record.
        lines.add(cnv("eval-cnv", 80000, 81000, 2, 1, 2));
        // An eval DEL against a truth CNV: allowed, because the CNV is the truth record. The same
        // pair of types, the other way round.
        lines.add(record("eval-del", 100000, "DEL", 101000, null, "0/1", "0/1", "0/0"));
        // Two CNVs, scored on copy state: s2 agrees, s3 does not.
        lines.add(cnv("cnv-pair", 120000, 121000, 3, 1, 3));
        // A truth record that carries its own allele counts.
        lines.add(record("af-given", 140000, "DEL", 141000, null, "0/1", "0/1", "0/0"));
        // And one that does not, so they are recounted from the truth genotypes.
        lines.add(record("af-missing", 160000, "DEL", 161000, null, "0/1", "0/1", "0/0"));
        lines.add("");
        return String.join("\n", lines);
    }

    static String buildTruth() {
        final List<String> lines = header(TRUTH_SAMPLES);
        lines.add(record("t-exact", 1000, "DEL", 2000, null, "0/1", "0/0", "0/1"));
        lines.add(record("t-discordant", 20000, "DEL", 21000, null, "1/1", "0/0", "0/1"));
        // Exact on the start and off by 200 at the end: total 200, closest breakend 0.
        lines.add(record("t-tie-one", 60000, "DEL", 61200, null, "0/0", "0/0", "0/0"));
        // Off by 100 on BOTH ends: the same total 200, closest breakend 100. Written second
        // because the finder refuses records that are not in position order, which is what the
        // `unsorted` run is for.
        lines.add(record("t-tie-both", 60100, "DEL", 61100, null, "0/1", "0/1", "0/0"));
        lines.add(record("t-eval-cnv", 80000, "DEL", 81000, null, "0/1", "0/1", "0/0"));
        lines.add(cnv("t-eval-del", 100000, 101000, 1, 1, 2));
        // s2 agrees with the eval copy state and s3 does not, which is what puts the record
        // concordance at a half rather than at either end.
        lines.add(cnv("t-cnv-pair", 120000, 121000, 1, 2, 1));
        lines.add(record("t-af-given", 140000, "DEL", 141000, "AC=3;AF=0.5;AN=6",
                "0/1", "0/1", "0/1"));
        lines.add(record("t-af-missing", 160000, "DEL", 161000, null, "1/1", "0/1", "0/0"));
        lines.add("");
        return String.join("\n", lines);
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("sv-concordance-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# SVConcordanceDump: which truth variant an eval variant is judged to "
                + "be");

        final Path dict = SVClusterDump.writeDictionary(dir);
        final String eval = buildEval();
        final String truth = buildTruth();
        final Path evalPath = write(dir, "eval.vcf", eval);
        final Path truthPath = write(dir, "truth.vcf", truth);
        System.out.printf("vcf\teval=%s%n", ReferenceQueryDump.escape(eval));
        System.out.printf("vcf\ttruth=%s%n", ReferenceQueryDump.escape(truth));

        run(dir, "default", evalPath, truthPath, dict, List.of());
        // The same run without the output sort. For an in-order input the records complete in
        // the order they were added, so this is a control: it says the flag is a memory setting
        // and not a behaviour.
        run(dir, "do-not-sort", evalPath, truthPath, dict, List.of("--do-not-sort", "true"));
        // A reciprocal overlap that only ONE pair fails: the tie pair is the only one whose
        // truth record is not an exact interval match, so raising the threshold to 0.999 breaks
        // that match alone and leaves the other eight untouched.
        run(dir, "overlap-high", evalPath, truthPath, dict,
                List.of("--depth-interval-overlap", "0.999"));

        // The dictionary is not taken from either VCF.
        runRaw(dir, "no-dictionary", List.of(
                "--eval", evalPath.toString(),
                "--truth", truthPath.toString(),
                "-O", dir.resolve("out-no-dictionary.vcf").toString()));

        // An eval VCF whose records are not in position order, which the finder refuses.
        final List<String> unsortedLines = header(EVAL_SAMPLES);
        unsortedLines.add(record("late", 5000, "DEL", 6000, null, "0/1", "0/1", "0/0"));
        unsortedLines.add(record("early", 1000, "DEL", 2000, null, "0/1", "0/1", "0/0"));
        unsortedLines.add("");
        final Path unsorted = write(dir, "unsorted.vcf", String.join("\n", unsortedLines));
        run(dir, "unsorted", unsorted, truthPath, dict, List.of());
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static void run(final Path dir, final String label, final Path eval, final Path truth,
                    final Path dict, final List<String> extra) throws Exception {
        final List<String> argv = new ArrayList<>(List.of(
                "--eval", eval.toString(),
                "--truth", truth.toString(),
                "-O", dir.resolve("out-" + label + ".vcf").toString(),
                "--sequence-dictionary", dict.toString()));
        argv.addAll(extra);
        runRaw(dir, label, argv);
    }

    static void runRaw(final Path dir, final String label, final List<String> argv)
            throws Exception {
        try {
            new SVConcordance().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(cause.getMessage()), dir)));
            return;
        }
        final Path out = dir.resolve("out-" + label + ".vcf");
        if (!Files.exists(out)) {
            return;
        }
        final StringBuilder body = new StringBuilder();
        for (final String line : Files.readString(out).split("\n", -1)) {
            if (!line.startsWith("##") && !line.isEmpty()) {
                body.append(line).append("\n");
            }
        }
        System.out.printf("out\t%s=%s%n", label,
                ReferenceQueryDump.escape(masked(body.toString(), dir)));
    }

    /** The temporary directory is not the same twice, so it never reaches the golden. */
    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}

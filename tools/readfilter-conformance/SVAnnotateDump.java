/*
 * SVAnnotate's predicted consequences, taken from the reference.
 *
 * What a structural variant is predicted to do to a gene. The answer is a word per transcript, and
 * which word depends on the SV type, on which features the two breakpoints land in, and on whether
 * the variant SPANS a feature or merely overlaps it.
 *
 * Ten behaviours this is built to catch.
 *
 *   - SPANNING AND OVERLAPPING ARE DIFFERENT QUESTIONS: a duplication that contains the whole
 *     transcript is COPY_GAIN, one that contains a whole coding exon without a breakpoint in one is
 *     INT_EXON_DUP, and one with a breakpoint in a coding exon is PARTIAL_EXON_DUP;
 *   - TWO BREAKPOINTS IN CODING SEQUENCE MAKE IT LOF, and so does one in coding sequence with one
 *     in a UTR, which is the only place the UTR count is read;
 *   - THE TRANSCRIPTION START SITE IS A SPECIAL CASE for deletions, inversions and duplications,
 *     and it is the transcript's END on the minus strand;
 *   - A DELETION AND AN INVERSION SHARE THEIR RULE, the inversion only adding the spanning case;
 *   - A BREAKEND USES THE DELETION RULE AND THEN DOWNGRADES ITS LOF to BREAKEND_EXONIC, because a
 *     low-confidence breakend should not be called loss of function;
 *   - A MULTIALLELIC CNV IS ANNOTATED AS A DUPLICATION AND THEN RECLASSIFIED: six of the
 *     duplication answers become MSV_EXON_OVERLAP and the rest are kept;
 *   - A TRANSLOCATION IS LOF WITHOUT LOOKING AT ANYTHING, because breaking a gene at all is
 *     predicted to break it;
 *   - THE PROMOTER IS INFERRED, not read: it is the window upstream of the TSS, and a window of
 *     zero removes the annotation entirely;
 *   - A VARIANT THAT HITS NOTHING GETS THE NEAREST TSS instead, which is a different field;
 *   - AND THE CONSEQUENCE LISTS ARE SORTED, which is the only reason two runs can be compared at
 *     all.
 *
 * Output:
 *
 *     gtf\tmain=<the protein-coding GTF, escaped>
 *     bed\tmain=<the non-coding BED, escaped>
 *     vcf\tinput=<the whole input vcf, escaped>
 *     out\t<label>=<the whole output vcf without its header, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: SVAnnotateDump
 */

import org.broadinstitute.hellbender.tools.walkers.sv.SVAnnotate;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class SVAnnotateDump {

    /** The attributes every feature of one transcript repeats. */
    static String attributes(final String gene, final String number, final String extra) {
        return "gene_id \"ENSG0000000000" + number + ".1\"; transcript_id \"ENST0000000000"
                + number + ".1\"; gene_type \"protein_coding\"; gene_name \"" + gene
                + "\"; transcript_type \"protein_coding\"; transcript_name \"" + gene
                + "-201\"; tag \"MANE_Select\";" + extra;
    }

    static String feature(final String type, final int start, final int end, final String strand,
                          final String gene, final String number) {
        return "chr1\tensembl_havana\t" + type + "\t" + start + "\t" + end + "\t.\t" + strand
                + "\t" + ("CDS".equals(type) ? "0" : ".") + "\t" + attributes(gene, number, "");
    }

    /**
     * Two genes, one per strand, with the same shape: five exons, the middle two coding.
     *
     * ALPHA runs forward from 10000, so its transcription start site is 10000. BETA runs backward
     * to 21000, so its transcription start site is 21000, which is the only way to tell that the
     * rule reads the strand rather than the smaller coordinate.
     */
    static String buildGtf() {
        final List<String> lines = new ArrayList<>();
        // ALPHA, plus strand.
        lines.add("chr1\tensembl_havana\tgene\t10000\t11000\t.\t+\t.\t"
                + "gene_id \"ENSG00000000001.1\"; gene_type \"protein_coding\"; "
                + "gene_name \"ALPHA\";");
        lines.add(feature("transcript", 10000, 11000, "+", "ALPHA", "1"));
        lines.add(feature("exon", 10000, 10200, "+", "ALPHA", "1"));
        lines.add(feature("exon", 10300, 10400, "+", "ALPHA", "1"));
        lines.add(feature("exon", 10500, 10600, "+", "ALPHA", "1"));
        lines.add(feature("exon", 10700, 10800, "+", "ALPHA", "1"));
        lines.add(feature("exon", 10900, 11000, "+", "ALPHA", "1"));
        lines.add(feature("CDS", 10300, 10400, "+", "ALPHA", "1"));
        lines.add(feature("CDS", 10500, 10600, "+", "ALPHA", "1"));
        lines.add(feature("start_codon", 10300, 10302, "+", "ALPHA", "1"));
        lines.add(feature("stop_codon", 10601, 10603, "+", "ALPHA", "1"));
        lines.add(feature("UTR", 10000, 10200, "+", "ALPHA", "1"));
        lines.add(feature("UTR", 10700, 10800, "+", "ALPHA", "1"));
        lines.add(feature("UTR", 10900, 11000, "+", "ALPHA", "1"));
        // BETA, minus strand.
        lines.add("chr1\tensembl_havana\tgene\t20000\t21000\t.\t-\t.\t"
                + "gene_id \"ENSG00000000002.1\"; gene_type \"protein_coding\"; "
                + "gene_name \"BETA\";");
        lines.add(feature("transcript", 20100, 21000, "-", "BETA", "2"));
        lines.add(feature("exon", 20900, 21000, "-", "BETA", "2"));
        lines.add(feature("exon", 20700, 20800, "-", "BETA", "2"));
        lines.add(feature("exon", 20500, 20600, "-", "BETA", "2"));
        lines.add(feature("exon", 20300, 20400, "-", "BETA", "2"));
        lines.add(feature("exon", 20100, 20200, "-", "BETA", "2"));
        lines.add(feature("CDS", 20500, 20600, "-", "BETA", "2"));
        lines.add(feature("CDS", 20300, 20400, "-", "BETA", "2"));
        lines.add(feature("start_codon", 20598, 20600, "-", "BETA", "2"));
        lines.add(feature("stop_codon", 20297, 20299, "-", "BETA", "2"));
        lines.add(feature("UTR", 20900, 21000, "-", "BETA", "2"));
        lines.add(feature("UTR", 20700, 20800, "-", "BETA", "2"));
        lines.add(feature("UTR", 20100, 20200, "-", "BETA", "2"));
        lines.add("");
        return String.join("\n", lines);
    }

    /**
     * BED is half-open, so a start of 10399 is base 10400.
     *
     * With NO header, whatever `--non-coding-bed`'s own documentation says: see
     * {@link #buildBedWithHeader}.
     */
    static String buildBed() {
        return String.join("\n",
                "chr1\t10399\t10500\tDNase\t.\t+",
                "chr1\t11199\t11600\tEnhancer\t.\t+",
                "chr1\t49999\t50500\tTAD\t.\t+",
                "");
    }

    /**
     * The same file with the header the argument documentation asks for, which is refused: the BED
     * codec reads the header row as a feature and fails to parse `start` as a number.
     */
    static String buildBedWithHeader() {
        return "chrom\tstart\tend\tname\tscore\tstrand\n" + buildBed();
    }

    static String record(final String id, final int start, final String type, final int end,
                         final String extra) {
        final int length = "INS".equals(type) || "BND".equals(type) || "CTX".equals(type)
                ? -1 : end - start + 1;
        final StringBuilder info = new StringBuilder("SVTYPE=" + type + ";END=" + end
                + ";SVLEN=" + length + ";ALGORITHMS=depth");
        if (extra != null) {
            info.append(';').append(extra);
        }
        return "chr1\t" + start + "\t" + id + "\tN\t<" + type + ">\t.\t.\t" + info
                + "\tGT:ECN\t0/1:2";
    }

    /** A breakend on chr1 with a real length, so it can also be read as a DEL or a DUP. */
    static String breakend(final String id, final int start, final int length,
                           final String strands) {
        return "chr1\t" + start + "\t" + id + "\tN\t<BND>\t.\t.\tSVTYPE=BND;END=" + start
                + ";SVLEN=" + length + ";ALGORITHMS=depth;CHR2=chr1;END2=" + (start + length)
                + ";STRANDS=" + strands + "\tGT:ECN\t0/1:2";
    }

    static String buildVcf() {
        final List<String> lines = new ArrayList<>(List.of(
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=199980>",
                "##contig=<ID=chr2,length=199980>",
                "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">",
                "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length\">",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">",
                "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Second contig\">",
                "##INFO=<ID=END2,Number=1,Type=Integer,Description=\"Second position\">",
                "##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Strands\">",
                "##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Algorithms\">",
                "##INFO=<ID=CPX_TYPE,Number=1,Type=String,Description=\"Complex type\">",
                "##INFO=<ID=CPX_INTERVALS,Number=.,Type=String,Description=\"Complex intervals\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=ECN,Number=1,Type=Integer,Description=\"Expected copy number\">",
                "##ALT=<ID=DEL,Description=\"Deletion\">",
                "##ALT=<ID=DUP,Description=\"Duplication\">",
                "##ALT=<ID=CNV,Description=\"Copy number variant\">",
                "##ALT=<ID=INS,Description=\"Insertion\">",
                "##ALT=<ID=INV,Description=\"Inversion\">",
                "##ALT=<ID=BND,Description=\"Breakend\">",
                "##ALT=<ID=CTX,Description=\"Translocation\">",
                "##ALT=<ID=CPX,Description=\"Complex\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1"));
        // Three thousand bases upstream of ALPHA's start: outside the default window and inside a
        // wide one, which is the only thing that makes the window's width observable.
        lines.add(record("far-upstream", 7000, "DEL", 7100, null));
        // The promoter, which is inferred as the window upstream of ALPHA's start.
        lines.add(record("promoter", 9200, "DEL", 9400, null));
        // A deletion over the transcription start site, which is LOF whatever else it hits.
        lines.add(record("del-tss", 9950, "DEL", 10050, null));
        // A deletion over coding sequence.
        lines.add(record("del-cds", 10350, "DEL", 10450, null));
        // A deletion inside an intron, which also overlaps the DNase element.
        lines.add(record("del-intron", 10420, "DEL", 10480, null));
        // A deletion over a UTR and nothing else.
        lines.add(record("del-utr", 10720, "DEL", 10780, null));
        // A duplication containing the whole transcript.
        lines.add(record("dup-span", 9000, "DUP", 12000, null));
        // A duplication over the start site, which does not contain the transcript.
        lines.add(record("dup-tss", 9950, "DUP", 10500, null));
        // A duplication with one breakpoint inside and one past the end.
        lines.add(record("dup-partial", 10500, "DUP", 12000, null));
        // A duplication with one breakpoint in coding sequence.
        lines.add(record("dup-partial-exon", 10350, "DUP", 10450, null));
        // A duplication containing both coding exons whole, with neither breakpoint in one.
        lines.add(record("dup-int-exon", 10250, "DUP", 10650, null));
        // A duplication with a breakpoint in EACH coding exon.
        lines.add(record("dup-lof", 10350, "DUP", 10550, null));
        // A duplication with both breakpoints in a UTR.
        lines.add(record("dup-utr", 10720, "DUP", 10780, null));
        // A duplication inside an intron.
        lines.add(record("dup-intron", 10420, "DUP", 10480, null));
        // A multiallelic CNV over the same coding exon as dup-partial-exon.
        lines.add(record("cnv-cds", 10350, "CNV", 10450, null));
        // And one inside an intron, whose answer is not in the reclassified set.
        lines.add(record("cnv-intron", 10420, "CNV", 10480, null));
        // An inversion containing the transcript, and one over coding sequence.
        lines.add(record("inv-span", 9000, "INV", 12000, "STRANDS=++"));
        lines.add(record("inv-cds", 10350, "INV", 10450, "STRANDS=++"));
        // An insertion inside coding sequence.
        lines.add(record("ins-cds", 10350, "INS", 10350, null));
        // A breakend in coding sequence, and one in an intron. Both carry a POSITIVE SVLEN,
        // because a BND read as a deletion is given the interval pos..pos+SVLEN and the
        // conventional -1 makes that interval empty: see the `bnd-no-length` run.
        lines.add(breakend("bnd-cds", 10350, 4650, "+-"));
        lines.add(breakend("bnd-intron", 10450, 4550, "+-"));
        // A breakend whose strands make it a duplication rather than a deletion.
        lines.add(breakend("bnd-dup", 10350, 300, "-+"));
        // A translocation out of the gene, which is LOF without looking.
        lines.add(record("ctx", 10450, "CTX", 10450, "CHR2=chr2;END2=1000;STRANDS=+-"));
        // A variant that hits no gene at all, which overlaps the TAD element.
        lines.add(record("intergenic", 50100, "DEL", 50300, null));
        // The minus-strand gene's transcription start site is its END.
        lines.add(record("del-tss-minus", 20950, "DEL", 21050, null));
        // And its smaller coordinate is not a start site, so this is only UTR.
        lines.add(record("del-end-minus", 20120, "DEL", 20180, null));
        lines.add("");
        return String.join("\n", lines);
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("sv-annotate-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# SVAnnotateDump: what a structural variant is predicted to do to a "
                + "gene");

        final String gtf = buildGtf();
        final String bed = buildBed();
        final String vcf = buildVcf();
        final Path gtfPath = write(dir, "genes.gtf", gtf);
        final Path bedPath = write(dir, "noncoding.bed", bed);
        final Path input = write(dir, "input.vcf", vcf);
        System.out.printf("gtf\tmain=%s%n", ReferenceQueryDump.escape(gtf));
        System.out.printf("bed\tmain=%s%n", ReferenceQueryDump.escape(bed));
        System.out.printf("vcf\tinput=%s%n", ReferenceQueryDump.escape(vcf));

        run(dir, "default", input, List.of(
                "--protein-coding-gtf", gtfPath.toString(),
                "--non-coding-bed", bedPath.toString()));
        // A promoter window of zero, which the argument accepts (minValue = 0) and which then
        // builds an empty interval and throws.
        run(dir, "zero-promoter", input, List.of(
                "--protein-coding-gtf", gtfPath.toString(),
                "--non-coding-bed", bedPath.toString(),
                "--promoter-window-length", "0"));
        // A wider window, which reaches further upstream.
        run(dir, "wide-promoter", input, List.of(
                "--protein-coding-gtf", gtfPath.toString(),
                "--non-coding-bed", bedPath.toString(),
                "--promoter-window-length", "5000"));
        // Breakends short enough to be read as deletions or duplications instead.
        run(dir, "breakend-as-cnv", input, List.of(
                "--protein-coding-gtf", gtfPath.toString(),
                "--non-coding-bed", bedPath.toString(),
                "--max-breakend-as-cnv-length", "1000000"));
        // The GTF alone, and the BED alone.
        run(dir, "gtf-only", input, List.of("--protein-coding-gtf", gtfPath.toString()));
        run(dir, "bed-only", input, List.of("--non-coding-bed", bedPath.toString()));
        // Neither, which annotates nothing at all rather than refusing.
        run(dir, "neither", input, List.of());
        // The BED with the header its own argument documentation asks for.
        final Path headered = write(dir, "headered.bed", buildBedWithHeader());
        run(dir, "bed-header", input, List.of("--non-coding-bed", headered.toString()));

        // A complex variant with no intervals, and one with no type.
        final Path noIntervals = write(dir, "no-intervals.vcf", complex(null, "delINV"));
        run(dir, "cpx-no-intervals", noIntervals,
                List.of("--protein-coding-gtf", gtfPath.toString()));
        final Path noType = write(dir, "no-type.vcf",
                complex("chr1:10300-10400", null));
        run(dir, "cpx-no-type", noType, List.of("--protein-coding-gtf", gtfPath.toString()));
        // A breakend carrying the conventional SVLEN of -1, which is accepted as being under the
        // maximum and then given an empty interval.
        final Path noLength = write(dir, "no-length.vcf", singleRecord(
                record("bnd-no-length", 10350, "BND", 10350, "CHR2=chr1;END2=15000;STRANDS=+-")));
        run(dir, "bnd-no-length", noLength, List.of(
                "--protein-coding-gtf", gtfPath.toString(),
                "--max-breakend-as-cnv-length", "1000000"));
        // And a translocation with no second contig.
        final Path noChr2 = write(dir, "no-chr2.vcf", singleRecord(
                record("ctx-bare", 10450, "CTX", 10450, "STRANDS=+-")));
        run(dir, "ctx-no-chr2", noChr2, List.of("--protein-coding-gtf", gtfPath.toString()));
    }

    /** One record in a file with the same header as the main input. */
    static String singleRecord(final String line) {
        final String[] all = buildVcf().split("\n");
        final List<String> lines = new ArrayList<>();
        for (final String header : all) {
            if (header.startsWith("#")) {
                lines.add(header);
            }
        }
        lines.add(line);
        lines.add("");
        return String.join("\n", lines);
    }

    static String complex(final String intervals, final String type) {
        final StringBuilder extra = new StringBuilder();
        if (intervals != null) {
            extra.append("CPX_INTERVALS=DEL_").append(intervals);
        }
        if (type != null) {
            if (extra.length() > 0) {
                extra.append(';');
            }
            extra.append("CPX_TYPE=").append(type);
        }
        return singleRecord(record("cpx", 10300, "CPX", 10400,
                extra.length() == 0 ? null : extra.toString()));
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static void run(final Path dir, final String label, final Path input, final List<String> extra)
            throws Exception {
        final Path out = dir.resolve("out-" + label + ".vcf");
        final List<String> argv = new ArrayList<>(List.of(
                "-V", input.toString(),
                "-O", out.toString()));
        argv.addAll(extra);
        try {
            new SVAnnotate().instanceMain(argv.toArray(new String[0]));
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

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}

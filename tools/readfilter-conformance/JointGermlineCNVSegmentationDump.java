/*
 * JointGermlineCNVSegmentation's segments, taken from the reference.
 *
 * How single-sample gCNV calls become one joint call set. Two engines run in sequence: a
 * defragmenter that joins one sample's adjacent segments, and a max-clique cluster that joins
 * different samples' events into one site. What is measurable is which calls survive the entry
 * filter, what ploidy each genotype is given, and which segments end up as one record.
 *
 * Nine behaviours this is built to catch.
 *
 *   - FOUR SEPARATE REASONS DROP A SINGLE-SAMPLE CALL BEFORE ANYTHING ELSE RUNS: a hom-ref
 *     genotype, a no-call with no CN, a QS below --minimum-qs-score, and a "null call", which is a
 *     no-call whose CN is exactly 0;
 *   - THE QS TEST IS STRICTLY LESS THAN, so a call exactly at the threshold survives;
 *   - PLOIDY ON AN AUTOSOME IS AN ARGUMENT and on an allosome is the PEDIGREE'S SEX: chrX is 2 for
 *     a female and 1 for a male, chrY is 0 for a female and 1 for a male, and an UNKNOWN sex is 1
 *     on BOTH, so it is not the female answer on either;
 *   - AN ECN ALREADY ON THE INPUT GENOTYPE DOES NOT REACH THE OUTPUT: the ploidy is derived again
 *     for every sample when the site is written, so a genotype carrying ECN=7 on chrX still comes
 *     out haploid for a male and the site still reports AN=4;
 *   - THE GENOTYPE IS PADDED TO ITS PLOIDY with reference alleles, and a single no-call allele is
 *     a special case that becomes a no-call of the full ploidy instead;
 *   - DEFRAGMENTATION PADS BY A FRACTION OF EACH RECORD'S OWN LENGTH and compares against every
 *     MEMBER of a joined run rather than against the run's span, so the run reaches only as far as
 *     its longest member: at 0 nothing joins, at the default 0.25 both the short pair and the long
 *     pair join across the same hundred-base gap while 90000 stays out, and at 1.0 everything from
 *     20000 to 91000 is one record because the twenty-thousand-base members reach that far alone;
 *   - --min-sample-set-fraction-overlap IS ABOUT THE DEFRAGMENTER, which only ever sees one
 *     sample, so a value of 1.0 changes nothing: it is measured as the control that says so;
 *   - A MULTI-SAMPLE INPUT SKIPS DEFRAGMENTATION ENTIRELY, being assumed pre-clustered;
 *   - AND THE PEDIGREE IS VALIDATED STRICTLY, so a sample in the VCFs and not in the pedigree is
 *     refused before a record is read.
 *
 * Output:
 *
 *     ped\tmain=<the pedigree, escaped>
 *     vcf\t<sample>=<that sample's whole gCNV vcf, escaped>
 *     out\t<label>=<the whole output vcf without its header, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: JointGermlineCNVSegmentationDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.tools.walkers.sv.JointGermlineCNVSegmentation;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class JointGermlineCNVSegmentationDump {

    static final int CONTIG_LENGTH = 199980;

    /** A gCNV segments VCF for one sample. */
    static String vcf(final String sample, final List<String> records) {
        final List<String> lines = new ArrayList<>(List.of(
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=" + CONTIG_LENGTH + ">",
                "##contig=<ID=chrX,length=" + CONTIG_LENGTH + ">",
                "##contig=<ID=chrY,length=" + CONTIG_LENGTH + ">",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number\">",
                "##FORMAT=<ID=NP,Number=1,Type=Integer,Description=\"Number of points\">",
                "##FORMAT=<ID=QA,Number=1,Type=Integer,Description=\"Quality, any\">",
                "##FORMAT=<ID=QS,Number=1,Type=Integer,Description=\"Quality, some\">",
                "##FORMAT=<ID=QSE,Number=1,Type=Integer,Description=\"Quality, end\">",
                "##FORMAT=<ID=QSS,Number=1,Type=Integer,Description=\"Quality, start\">",
                "##FORMAT=<ID=ECN,Number=1,Type=Integer,Description=\"Expected copy number\">",
                "##ALT=<ID=DEL,Description=\"Deletion\">",
                "##ALT=<ID=DUP,Description=\"Duplication\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample));
        lines.addAll(records);
        lines.add("");
        return String.join("\n", lines);
    }

    /** One gCNV segment. `format` and `values` carry whatever FORMAT fields the case needs. */
    static String segment(final String contig, final int start, final int end, final String alt,
                          final String genotype, final String format, final String values) {
        return contig + "\t" + start + "\tseg\tN\t" + alt + "\t.\t.\tEND=" + end + "\tGT:" + format
                + "\t" + genotype + ":" + values;
    }

    /** The ordinary case: a called segment with a copy number and a quality. */
    static String call(final String contig, final int start, final int end, final String alt,
                       final int copyNumber, final int quality) {
        return segment(contig, start, end, alt, "0/1", "CN:QS", copyNumber + ":" + quality);
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("joint-germline-cnv-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# JointGermlineCNVSegmentationDump: how single-sample gCNV calls "
                + "become one joint call set");

        final Path dict = writeDictionary(dir);
        final Path fasta = writeReference(dir);

        // FAMILY, sample, father, mother, sex (1 male, 2 female, 0 unknown), phenotype.
        final String pedigree = String.join("\n",
                "fam\tmale1\t0\t0\t1\t0",
                "fam\tfemale1\t0\t0\t2\t0",
                "fam\tunknown1\t0\t0\t0\t0",
                "");
        final Path ped = write(dir, "family.ped", pedigree);
        System.out.printf("ped\tmain=%s%n", ReferenceQueryDump.escape(pedigree));

        // male1: the four entry filters, the QS boundary, and two short segments to defragment.
        final List<String> male = new ArrayList<>();
        // Dropped: hom-ref.
        male.add(call("chr1", 1000, 2000, "<DEL>", 1, 100).replace("0/1:", "0/0:"));
        // Dropped: a no-call with no CN at all.
        male.add(segment("chr1", 3000, 4000, "<DEL>", "./.", "QS", "100"));
        // Dropped: a null call, which is a no-call whose CN is exactly zero.
        male.add(segment("chr1", 5000, 6000, "<DEL>", "./.", "CN:QS", "0:100"));
        // Dropped: a quality one below the threshold the runs use.
        male.add(call("chr1", 7000, 8000, "<DEL>", 1, 49));
        // Kept: a quality exactly at it, which is what makes the test strictly less than.
        male.add(call("chr1", 9000, 10000, "<DEL>", 1, 50));
        // Two short segments a hundred bases apart, and below them two long ones the same hundred
        // apart, so the two pairs differ only in the length the padding is a fraction OF.
        male.add(call("chr1", 20000, 20500, "<DEL>", 1, 100));
        male.add(call("chr1", 20600, 21100, "<DEL>", 1, 100));
        male.add(call("chr1", 40000, 60000, "<DEL>", 1, 100));
        male.add(call("chr1", 60100, 80100, "<DEL>", 1, 100));
        // A single no-call allele, which is padded to a full no-call rather than to a ref call.
        male.add(segment("chr1", 90000, 91000, "<DEL>", ".", "CN:QS", "1:100"));
        // The allosomes, whose ploidy comes from the pedigree. A male is haploid on both, and a
        // genotype carrying more alleles than its ploidy is refused outright, so these are written
        // with ONE allele. The records are in dictionary order because the driving iterator over
        // the three VCFs requires it.
        male.add(segment("chrX", 1000, 2000, "<DEL>", "1", "CN:QS", "0:100"));
        // A genotype that already carries an ECN, which is read before the contig is.
        male.add(segment("chrX", 5000, 6000, "<DEL>", "0/1", "CN:QS:ECN", "0:100:7"));
        male.add(segment("chrY", 1000, 2000, "<DEL>", "1", "CN:QS", "0:100"));

        // female1 and unknown1: the same allosomal segments, for the sex rule.
        final List<String> female = new ArrayList<>();
        female.add(call("chr1", 20000, 20500, "<DEL>", 1, 100));
        // Diploid on chrX, which is the answer a male does not get.
        female.add(call("chrX", 1000, 2000, "<DEL>", 1, 100));
        final List<String> unknown = new ArrayList<>();
        // An unknown sex is haploid on chrX, so it is NOT given the female answer.
        unknown.add(segment("chrX", 1000, 2000, "<DEL>", "1", "CN:QS", "0:100"));
        unknown.add(segment("chrY", 1000, 2000, "<DEL>", "1", "CN:QS", "0:100"));

        final Path maleVcf = write(dir, "male1.vcf", vcf("male1", male));
        final Path femaleVcf = write(dir, "female1.vcf", vcf("female1", female));
        final Path unknownVcf = write(dir, "unknown1.vcf", vcf("unknown1", unknown));
        System.out.printf("vcf\tmale1=%s%n", ReferenceQueryDump.escape(vcf("male1", male)));
        System.out.printf("vcf\tfemale1=%s%n", ReferenceQueryDump.escape(vcf("female1", female)));
        System.out.printf("vcf\tunknown1=%s%n", ReferenceQueryDump.escape(vcf("unknown1", unknown)));

        final List<String> inputs = List.of(
                "-V", maleVcf.toString(),
                "-V", femaleVcf.toString(),
                "-V", unknownVcf.toString());

        run(dir, "default", inputs, ped, fasta, dict, List.of("--minimum-qs-score", "50"));
        // A threshold above every quality in the input, which drops every call and writes a header
        // with no records rather than refusing.
        run(dir, "high-quality", inputs, ped, fasta, dict,
                List.of("--minimum-qs-score", "101"));

        // ONE input VCF, which is the only way the defragmenter runs at all: more than one sample
        // is assumed to be pre-clustered and defragmentation is skipped entirely. The three runs
        // above are therefore clustering only, and these three are the defragmenter.
        final List<String> alone = List.of("-V", maleVcf.toString());
        run(dir, "single", alone, ped, fasta, dict, List.of("--minimum-qs-score", "50"));
        run(dir, "single-no-padding", alone, ped, fasta, dict,
                List.of("--minimum-qs-score", "50", "--defragmentation-padding-fraction", "0"));
        run(dir, "single-wide-padding", alone, ped, fasta, dict,
                List.of("--minimum-qs-score", "50", "--defragmentation-padding-fraction", "1.0"));
        // A sample overlap of one, which one sample always reaches, so it is the control that says
        // the argument is about DIFFERENT samples.
        run(dir, "single-sample-overlap", alone, ped, fasta, dict,
                List.of("--minimum-qs-score", "50", "--min-sample-set-fraction-overlap", "1.0"));
        // A different autosomal reference copy number, which is the ploidy on chr1.
        run(dir, "haploid-autosomes", inputs, ped, fasta, dict,
                List.of("--minimum-qs-score", "50", "--autosomal-ref-copy-number", "1"));

        // A sample in the VCFs and not in the pedigree.
        final Path shortPed = write(dir, "short.ped", "fam\tmale1\t0\t0\t1\t0\n");
        run(dir, "missing-from-pedigree", inputs, shortPed, fasta, dict,
                List.of("--minimum-qs-score", "50"));
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static void run(final Path dir, final String label, final List<String> inputs, final Path ped,
                    final Path fasta, final Path dict, final List<String> extra) throws Exception {
        final Path out = dir.resolve("out-" + label + ".vcf");
        final List<String> argv = new ArrayList<>(inputs);
        argv.addAll(List.of(
                "-O", out.toString(),
                "-R", fasta.toString(),
                "--sequence-dictionary", dict.toString(),
                "--pedigree", ped.toString()));
        argv.addAll(extra);
        try {
            new JointGermlineCNVSegmentation().instanceMain(argv.toArray(new String[0]));
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

    static Path writeReference(final Path dir) throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        final StringBuilder bases = new StringBuilder();
        for (final String contig : new String[] {"chr1", "chrX", "chrY"}) {
            bases.append(">").append(contig).append("\n");
            for (int i = 0; i < CONTIG_LENGTH / 60; i++) {
                bases.append("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n");
            }
        }
        Files.writeString(fasta, bases.toString(), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        return fasta;
    }

    static Path writeDictionary(final Path dir) throws Exception {
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", CONTIG_LENGTH),
                new SAMSequenceRecord("chrX", CONTIG_LENGTH),
                new SAMSequenceRecord("chrY", CONTIG_LENGTH)));
        final Path path = dir.resolve("reference.dict");
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(dictionary);
        try (final java.io.Writer writer = Files.newBufferedWriter(path)) {
            new htsjdk.samtools.SAMTextHeaderCodec().encode(writer, header);
        }
        return path;
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}

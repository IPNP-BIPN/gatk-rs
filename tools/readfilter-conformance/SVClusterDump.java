/*
 * SVCluster's clusters, taken from the reference.
 *
 * Which structural variants are judged the same event. What decides it is a predicate over a pair,
 * and which of three parameter sets that predicate uses is decided by the records' own ALGORITHMS
 * field rather than by any argument.
 *
 * Ten behaviours this is built to catch.
 *
 *   - THE PARAMETER SET IS CHOSEN BY THE PAIR, NOT BY THE USER: two depth-only records take the
 *     depth parameters, two PESR records the PESR ones, and one of each the MIXED ones, so a single
 *     run applies three different thresholds to different pairs;
 *   - SINGLE LINKAGE AND MAX CLIQUE DISAGREE ON A CHAIN: three records where A clusters with B and
 *     B with C but A not with C become ONE cluster under single linkage and more than one under max
 *     clique;
 *   - THE DEPTH PARAMETERS REQUIRE OVERLAP AND PROXIMITY, the PESR ones require either, which is
 *     what `requiresOverlapAndProximity` decides from the parameter set's own validity rule;
 *   - RECIPROCAL OVERLAP AND SIZE SIMILARITY ARE SEPARATE TESTS, and both must pass: two records
 *     can overlap reciprocally and still fail on size;
 *   - AN INSERTION IS GIVEN AN ASSUMED LENGTH for both tests, because it has none of its own;
 *   - INTERCHROMOSOMAL PAIRS SKIP OVERLAP ENTIRELY and are judged on breakend proximity alone;
 *   - ONLY BND AND INV REQUIRE MATCHING STRANDS, and a record with null strands matches anything;
 *   - A DELETION AND A DUPLICATION DO NOT CLUSTER unless one of them is a CNV, or --enable-cnv is
 *     given;
 *   - SAMPLE OVERLAP IS TESTED LAST because it is the slowest, and it is a fraction of the smaller
 *     carrier set;
 *   - AND THE MEMBER IDS ARE WRITTEN IN THE OUTPUT, which is how the clustering is observable
 *     without reading the collapsed representative.
 *
 * Output:
 *
 *     vcf\tinput=<the whole input vcf, escaped>
 *     out\t<label>=<the whole output vcf without its header, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: SVClusterDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.tools.walkers.sv.SVCluster;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class SVClusterDump {

    /**
     * The measured input.
     *
     * The deletions on chr1 are a chain: d1 and d2 overlap well, d2 and d3 overlap well, d1 and d3
     * do not. That is the pair that separates single linkage from max clique. The rest are one
     * behaviour each.
     */
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
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                // ECN, the expected copy number, is not optional: SVCallRecord refuses a genotype
                // without it.
                "##FORMAT=<ID=ECN,Number=1,Type=Integer,Description=\"Expected copy number\">",
                "##ALT=<ID=DEL,Description=\"Deletion\">",
                "##ALT=<ID=DUP,Description=\"Duplication\">",
                "##ALT=<ID=CNV,Description=\"Copy number variant\">",
                "##ALT=<ID=INS,Description=\"Insertion\">",
                "##ALT=<ID=INV,Description=\"Inversion\">",
                "##ALT=<ID=BND,Description=\"Breakend\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3"));
        // The chain, all depth-only so one parameter set decides it.
        lines.add(del("d1", 1000, 11000, "depth", "0/1", "0/1", "0/0"));
        lines.add(del("d2", 1600, 11600, "depth", "0/1", "0/1", "0/0"));
        lines.add(del("d3", 2200, 12200, "depth", "0/1", "0/1", "0/0"));
        // A deletion whose caller is PESR, so a pair with a depth one takes the MIXED parameters.
        lines.add(del("p1", 20000, 30000, "manta", "0/1", "0/0", "0/0"));
        lines.add(del("p2", 20200, 30200, "manta", "0/1", "0/0", "0/0"));
        lines.add(del("m1", 20400, 30400, "depth", "0/1", "0/0", "0/0"));
        // Reciprocal overlap without size similarity: same start, very different length.
        lines.add(del("z1", 50000, 70000, "depth", "0/1", "0/1", "0/0"));
        lines.add(del("z2", 50000, 51000, "depth", "0/1", "0/1", "0/0"));
        // A duplication next to a deletion, and a CNV between them.
        lines.add(record("u1", 80000, "DUP", 90000, 10001, null, null, null, "depth",
                "0/1", "0/1", "0/0"));
        lines.add(del("e1", 80000, 90000, "depth", "0/1", "0/1", "0/0"));
        lines.add(record("c1", 80000, "CNV", 90000, 10001, null, null, null, "depth",
                "0/1", "0/1", "0/0"));
        // Two insertions at the same locus, which have no length of their own.
        lines.add("chr1\t100000\ti1\tN\t<INS>\t.\t.\tSVTYPE=INS;END=100000;SVLEN=-1;"
                + "ALGORITHMS=manta\tGT:ECN\t0/1:2\t0/1:2\t0/0:2");
        lines.add("chr1\t100050\ti2\tN\t<INS>\t.\t.\tSVTYPE=INS;END=100050;SVLEN=-1;"
                + "ALGORITHMS=manta\tGT:ECN\t0/1:2\t0/1:2\t0/0:2");
        // Two inversions whose strands disagree.
        lines.add(record("v1", 120000, "INV", 130000, 10001, null, null, "++", "manta",
                "0/1", "0/1", "0/0"));
        lines.add(record("v2", 120100, "INV", 130100, 10001, null, null, "--", "manta",
                "0/1", "0/1", "0/0"));
        // Two breakends into chr2, close together.
        lines.add(record("b1", 150000, "BND", 150000, -1, "chr2", 50000, "+-", "manta",
                "0/1", "0/1", "0/0"));
        lines.add(record("b2", 150100, "BND", 150100, -1, "chr2", 50100, "+-", "manta",
                "0/1", "0/1", "0/0"));
        // Two deletions that overlap but share no carrier.
        lines.add(del("s1", 170000, 180000, "depth", "0/1", "0/0", "0/0"));
        lines.add(del("s2", 170100, 180100, "depth", "0/0", "0/0", "0/1"));
        lines.add("");
        return String.join("\n", lines);
    }

    static String del(final String id, final int start, final int end, final String algorithm,
                      final String... genotypes) {
        return record(id, start, "DEL", end, end - start + 1, null, null, null, algorithm,
                genotypes);
    }

    static String record(final String id, final int start, final String type, final int end,
                         final int length, final String contig2, final Integer end2,
                         final String strands, final String algorithm, final String... genotypes) {
        final StringBuilder info = new StringBuilder("SVTYPE=" + type + ";END=" + end
                + ";SVLEN=" + length + ";ALGORITHMS=" + algorithm);
        if (contig2 != null) {
            info.append(";CHR2=").append(contig2).append(";END2=").append(end2);
        }
        if (strands != null) {
            info.append(";STRANDS=").append(strands);
        }
        final List<String> withEcn = new ArrayList<>();
        for (final String genotype : genotypes) {
            withEcn.add(genotype + ":2");
        }
        return "chr1\t" + start + "\t" + id + "\tN\t<" + type + ">\t.\t.\t" + info
                + "\tGT:ECN\t" + String.join("\t", withEcn);
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("sv-cluster-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# SVClusterDump: which structural variants are judged the same event");

        final Path dict = writeDictionary(dir);
        // A real reference is required too: the tool reads the base at each record's position for
        // the collapsed representative's REF allele.
        final Path fasta = writeReference(dir);
        // The ploidy table is required: a header of SAMPLE plus one column per contig, and one
        // row per sample. Without it the run is refused before a record is read.
        final Path ploidy = write(dir, "ploidy.tsv", String.join("\n",
                "SAMPLE\tchr1\tchr2",
                "s1\t2\t2",
                "s2\t2\t2",
                "s3\t2\t2",
                ""));
        System.out.printf("ploidy\ttable=%s%n",
                ReferenceQueryDump.escape(Files.readString(ploidy)));
        final String vcf = buildVcf();
        final Path input = write(dir, "input.vcf", vcf);
        System.out.printf("vcf\tinput=%s%n", ReferenceQueryDump.escape(vcf));

        // The two algorithms over the same records, which is the chain test.
        run(dir, "single-linkage", input, dict, List.of("--algorithm", "SINGLE_LINKAGE"));
        run(dir, "max-clique", input, dict, List.of("--algorithm", "MAX_CLIQUE"));
        // Deletions and duplications allowed to cluster with each other.
        run(dir, "enable-cnv", input, dict, List.of("--algorithm", "SINGLE_LINKAGE",
                "--enable-cnv", "true"));
        // --enable-cnv is only observable under max clique: under single linkage the CNV already
        // chains the deletion to the duplication, so allowing them to link directly changes
        // nothing.
        run(dir, "max-clique-enable-cnv", input, dict, List.of("--algorithm", "MAX_CLIQUE",
                "--enable-cnv", "true"));
        // A depth overlap threshold nothing reaches, and one everything reaches.
        run(dir, "depth-overlap-high", input, dict, List.of("--algorithm", "SINGLE_LINKAGE",
                "--depth-interval-overlap", "0.99"));
        run(dir, "depth-overlap-low", input, dict, List.of("--algorithm", "SINGLE_LINKAGE",
                "--depth-interval-overlap", "0.1", "--depth-size-similarity", "0.1"));
        // Sample overlap, which is the last test the predicate makes.
        run(dir, "sample-overlap", input, dict, List.of("--algorithm", "SINGLE_LINKAGE",
                "--depth-sample-overlap", "0.5"));
        // A breakend window wide enough to reach the pairs proximity alone would not.
        run(dir, "wide-window", input, dict, List.of("--algorithm", "SINGLE_LINKAGE",
                "--pesr-breakend-window", "5000"));
        // The member ids omitted, which is the only thing that hides the clustering.
        run(dir, "omit-members", input, dict, List.of("--algorithm", "SINGLE_LINKAGE",
                "--omit-members", "true"));
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static void run(final Path dir, final String label, final Path input, final Path dict,
                    final List<String> extra) throws Exception {
        final Path out = dir.resolve("out-" + label + ".vcf");
        final List<String> argv = new ArrayList<>(List.of(
                "-V", input.toString(),
                "-O", out.toString(),
                "--sequence-dictionary", dict.toString(),
                "--ploidy-table", dir.resolve("ploidy.tsv").toString(),
                "-R", dir.resolve("reference.fasta").toString()));
        argv.addAll(extra);
        try {
            new SVCluster().instanceMain(argv.toArray(new String[0]));
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

    /**
     * Two contigs of 199980 bases, which is 3333 lines of 60.
     *
     * The length is not round on purpose: it has to be a whole number of FASTA lines, and a
     * dictionary claiming 200000 against a file holding 199980 is refused with `Index length does
     * not match dictionary length for contig: chr1`.
     */
    static Path writeReference(final Path dir) throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        final StringBuilder bases = new StringBuilder();
        for (final String contig : new String[] {"chr1", "chr2"}) {
            bases.append(">").append(contig).append("\n");
            for (int i = 0; i < 199980 / 60; i++) {
                bases.append("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n");
            }
        }
        Files.writeString(fasta, bases.toString(), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        // No CreateSequenceDictionary here: writeDictionary already wrote reference.dict beside
        // the fasta, which is where htsjdk looks for it.
        return fasta;
    }

    static Path writeDictionary(final Path dir) throws Exception {
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 199980),
                new SAMSequenceRecord("chr2", 199980)));
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

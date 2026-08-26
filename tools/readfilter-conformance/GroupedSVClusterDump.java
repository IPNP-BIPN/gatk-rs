/*
 * GroupedSVCluster's clusters, taken from the reference.
 *
 * The stratification engine and the clustering engine bolted together: each record is sorted into
 * exactly one stratum, and each stratum clusters with its own thresholds. What the tool adds over
 * the two engines is the pairing of the two configurations, and what it refuses when they disagree.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE TWO CONFIGURATIONS MUST NAME THE SAME GROUPS, checked twice: once by count and once by
 *     name, with two different messages;
 *   - AN EMPTY STRATIFICATION CONFIGURATION IS REFUSED before either engine is built;
 *   - A RECORD MATCHING TWO STRATA IS REFUSED HERE, with a DIFFERENT message from SVStratify's:
 *     this tool says the groups must be mutually exclusive rather than offering a flag, because
 *     allowing it would proliferate variants;
 *   - AN UNMATCHED RECORD IS NOT CLUSTERED AT ALL: it is written straight out carrying its own id
 *     as its only member and `default` as its stratum;
 *   - EACH STRATUM CLUSTERS WITH ITS OWN THRESHOLDS, so two records that would cluster under one
 *     group's parameters do not under another's;
 *   - THE STRATUM IS WRITTEN INTO EVERY RECORD, matched or not;
 *   - THE OUTPUT IS NOT SORTED and the tool says so by disabling its own index;
 *   - AND THE COLUMN-COUNT MESSAGE PRINTS THE SAME NUMBER TWICE, in a second copy of the same
 *     mistake the stratification parser makes.
 *
 * Output:
 *
 *     vcf\tinput=<the whole input vcf, escaped>
 *     strata\t<label>=<the stratification table, escaped>
 *     cluster\t<label>=<the clustering table, escaped>
 *     out\t<label>=<the whole output vcf without its header, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: GroupedSVClusterDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.tools.walkers.sv.GroupedSVCluster;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class GroupedSVClusterDump {

    /** Two strata over deletions, split by size, plus one over duplications. */
    static final String STRATA = String.join("\n",
            "NAME\tSVTYPE\tMIN_SIZE\tMAX_SIZE\tTRACKS",
            "DEL_small\tDEL\t-1\t5000\t-1",
            "DEL_large\tDEL\t5000\t-1\t-1",
            "DUP_any\tDUP\t-1\t-1\t-1",
            "");

    /** One row per stratum, with thresholds that differ between them on purpose. */
    static final String CLUSTERING = String.join("\n",
            "NAME\tRECIPROCAL_OVERLAP\tSIZE_SIMILARITY\tBREAKEND_WINDOW\tSAMPLE_OVERLAP",
            "DEL_small\t0.5\t0\t500\t0",
            "DEL_large\t0.99\t0\t100\t0",
            "DUP_any\t0.5\t0\t500\t0",
            "");

    static String buildVcf() {
        final List<String> lines = new ArrayList<>(List.of(
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=199980>",
                "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">",
                "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length\">",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">",
                "##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Algorithms\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=ECN,Number=1,Type=Integer,Description=\"Expected copy number\">",
                "##ALT=<ID=DEL,Description=\"Deletion\">",
                "##ALT=<ID=DUP,Description=\"Duplication\">",
                "##ALT=<ID=INS,Description=\"Insertion\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2"));
        // Two small deletions that cluster under DEL_small's half-overlap threshold.
        lines.add(record("small1", 1000, "DEL", 2000, 1001));
        lines.add(record("small2", 1300, "DEL", 2300, 1001));
        // Two large deletions that would cluster under a half-overlap threshold and do NOT under
        // DEL_large's 0.99, which is what shows each stratum using its own parameters.
        lines.add(record("large1", 20000, "DEL", 30000, 10001));
        lines.add(record("large2", 21000, "DEL", 31000, 10001));
        // Two duplications, in their own stratum.
        lines.add(record("dup1", 50000, "DUP", 51000, 1001));
        lines.add(record("dup2", 50300, "DUP", 51300, 1001));
        // An insertion, which no stratum claims.
        lines.add("chr1\t80000\tunmatched\tN\t<INS>\t.\t.\tSVTYPE=INS;END=80000;SVLEN=-1;"
                + "ALGORITHMS=depth\tGT:ECN\t0/1:2\t0/1:2");
        lines.add("");
        return String.join("\n", lines);
    }

    static String record(final String id, final int start, final String type, final int end,
                         final int length) {
        return "chr1\t" + start + "\t" + id + "\tN\t<" + type + ">\t.\t.\tSVTYPE=" + type
                + ";END=" + end + ";SVLEN=" + length + ";ALGORITHMS=depth\tGT:ECN\t0/1:2\t0/1:2";
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("grouped-sv-cluster-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# GroupedSVClusterDump: one stratum per record, one threshold set per stratum");

        final Path dict = writeDictionary(dir);
        final Path fasta = writeReference(dir);
        final Path ploidy = write(dir, "ploidy.tsv", "SAMPLE\tchr1\ns1\t2\ns2\t2\n");
        final String vcf = buildVcf();
        final Path input = write(dir, "input.vcf", vcf);
        System.out.printf("vcf\tinput=%s%n", ReferenceQueryDump.escape(vcf));

        final Path strata = write(dir, "strata.tsv", STRATA);
        final Path clustering = write(dir, "clustering.tsv", CLUSTERING);
        System.out.printf("strata\tmain=%s%n", ReferenceQueryDump.escape(STRATA));
        System.out.printf("cluster\tmain=%s%n", ReferenceQueryDump.escape(CLUSTERING));

        run(dir, "default", input, dict, fasta, ploidy, strata, clustering, List.of());
        run(dir, "max-clique", input, dict, fasta, ploidy, strata, clustering,
                List.of("--algorithm", "MAX_CLIQUE"));

        // A clustering configuration with one row too few, and one whose row is named differently.
        final Path shortTable = write(dir, "short.tsv", String.join("\n",
                "NAME\tRECIPROCAL_OVERLAP\tSIZE_SIMILARITY\tBREAKEND_WINDOW\tSAMPLE_OVERLAP",
                "DEL_small\t0.5\t0\t500\t0",
                "DEL_large\t0.99\t0\t100\t0",
                ""));
        System.out.printf("cluster\tshort=%s%n",
                ReferenceQueryDump.escape(Files.readString(shortTable)));
        run(dir, "too-few-groups", input, dict, fasta, ploidy, strata, shortTable, List.of());

        final Path renamed = write(dir, "renamed.tsv", CLUSTERING.replace("DUP_any", "DUP_other"));
        System.out.printf("cluster\trenamed=%s%n",
                ReferenceQueryDump.escape(Files.readString(renamed)));
        run(dir, "group-not-found", input, dict, fasta, ploidy, strata, renamed, List.of());

        // A stratification configuration with a header and no rows.
        final Path empty = write(dir, "empty-strata.tsv",
                "NAME\tSVTYPE\tMIN_SIZE\tMAX_SIZE\tTRACKS\n");
        System.out.printf("strata\tempty=%s%n",
                ReferenceQueryDump.escape(Files.readString(empty)));
        run(dir, "no-strata", input, dict, fasta, ploidy, empty, clustering, List.of());

        // Two strata a single record matches, which this tool refuses outright.
        final Path overlapping = write(dir, "overlapping.tsv", String.join("\n",
                "NAME\tSVTYPE\tMIN_SIZE\tMAX_SIZE\tTRACKS",
                "DEL_a\tDEL\t-1\t-1\t-1",
                "DEL_b\tDEL\t-1\t-1\t-1",
                ""));
        final Path twoGroups = write(dir, "two-groups.tsv", String.join("\n",
                "NAME\tRECIPROCAL_OVERLAP\tSIZE_SIMILARITY\tBREAKEND_WINDOW\tSAMPLE_OVERLAP",
                "DEL_a\t0.5\t0\t500\t0",
                "DEL_b\t0.5\t0\t500\t0",
                ""));
        System.out.printf("strata\toverlapping=%s%n",
                ReferenceQueryDump.escape(Files.readString(overlapping)));
        run(dir, "multiple-matches", input, dict, fasta, ploidy, overlapping, twoGroups, List.of());

        // A clustering table with an extra column, which is the second copy of the doubled-number
        // message.
        final Path extra = write(dir, "extra.tsv", String.join("\n",
                "NAME\tRECIPROCAL_OVERLAP\tSIZE_SIMILARITY\tBREAKEND_WINDOW\tSAMPLE_OVERLAP\tEXTRA",
                "DEL_small\t0.5\t0\t500\t0\tx",
                ""));
        run(dir, "extra-column", input, dict, fasta, ploidy, strata, extra, List.of());
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static void run(final Path dir, final String label, final Path input, final Path dict,
                    final Path fasta, final Path ploidy, final Path strata, final Path clustering,
                    final List<String> extra) throws Exception {
        final Path out = dir.resolve("out-" + label + ".vcf");
        final List<String> argv = new ArrayList<>(List.of(
                "-V", input.toString(),
                "-O", out.toString(),
                "-R", fasta.toString(),
                "--sequence-dictionary", dict.toString(),
                "--ploidy-table", ploidy.toString(),
                "--stratify-config", strata.toString(),
                "--clustering-config", clustering.toString()));
        argv.addAll(extra);
        try {
            new GroupedSVCluster().instanceMain(argv.toArray(new String[0]));
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

    /** One contig of 199980 bases, which is 3333 lines of 60. */
    static Path writeReference(final Path dir) throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        final StringBuilder bases = new StringBuilder(">chr1\n");
        for (int i = 0; i < 3333; i++) {
            bases.append("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n");
        }
        Files.writeString(fasta, bases.toString(), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        return fasta;
    }

    static Path writeDictionary(final Path dir) throws Exception {
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 199980)));
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

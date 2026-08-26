/*
 * SVStratify's strata, taken from the reference.
 *
 * Structural variants sorted into groups by type, size and reference-track overlap. What a stratum
 * matches is three independent tests, and two of the three ignore some of the thresholds they are
 * handed.
 *
 * Ten behaviours this is built to catch.
 *
 *   - MIN IS INCLUSIVE AND MAX IS EXCLUSIVE, so a stratum of [1000, 5000) takes a 1000 bp deletion
 *     and refuses a 5000 bp one;
 *   - A RECORD WITH NO LENGTH MATCHES ONLY A STRATUM WITH NEITHER BOUND, whatever its type, which
 *     is how BND records reach a stratum at all;
 *   - `-1` IS THE ONLY NEGATIVE THAT MEANS NULL: the null set is {"-1", "", "NULL", "NA"}, so -2
 *     parses as a number and is then refused as a negative bound;
 *   - AN INSERTION IGNORES --stratify-num-breakpoint-overlaps and requires exactly ONE end in a
 *     track, whatever the argument said;
 *   - A BREAKPOINT COUNTS ONCE PER END, NOT ONCE PER TRACK: countAnyTrackOverlap returns 1 as soon
 *     as any track overlaps, so two tracks over one end still count one;
 *   - THE OVERLAP FRACTION IS MEASURED AGAINST THE MERGED UNION OF EVERY NAMED TRACK, and compared
 *     with >=, so a stratum naming two tracks is easier to match than either alone;
 *   - BND AND CTX CANNOT CARRY A SIZE, and a config that gives them one is refused by name;
 *   - THE COLUMN-COUNT MESSAGE PRINTS THE SAME NUMBER TWICE, because both halves of it read
 *     columns.columnCount();
 *   - WITHOUT --split-output EVERY STRATUM WRITES TO THE SAME FILE, so a record matching two
 *     strata under --allow-multiple-matches appears TWICE, once per stratum name;
 *   - AND THE RESERVED NAME `default` IS REFUSED IN THE CONFIG, because it is the group unmatched
 *     records go to.
 *
 * Output:
 *
 *     config\t<label>=<the whole stratification table, escaped>
 *     track\t<name>=<the whole bed, escaped>
 *     vcf\tinput=<the whole input vcf, escaped>
 *     out\t<label>\t<file name>=<the whole output vcf, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: SVStratifyDump
 */

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.tools.walkers.sv.SVStratify;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Stream;

public class SVStratifyDump {

    /**
     * Two tracks over chr1, chosen so that one variant's ends land in different tracks.
     *
     * RM covers the left half of the measured deletions and SD the right, which is what separates
     * a stratum naming one track from a stratum naming both.
     */
    static final String RM_BED = String.join("\n",
            "chr1\t900\t1100",
            "chr1\t4900\t5100",
            "");

    static final String SD_BED = String.join("\n",
            "chr1\t1900\t2100",
            "chr1\t2500\t3500",
            "");

    /** The strata every successful run is measured against. */
    static final String CONFIG = String.join("\n",
            "NAME\tSVTYPE\tMIN_SIZE\tMAX_SIZE\tTRACKS",
            "DEL_small_RM\tDEL\t-1\t5000\tRM",
            "DEL_large_RM\tDEL\t5000\t-1\tRM",
            "DEL_small_both\tDEL\t-1\t5000\tRM,SD",
            "DUP_any\tDUP\t-1\t-1\t-1",
            "INS_RM\tINS\t-1\t-1\tRM",
            "BND_RM\tBND\t-1\t-1\tRM",
            "");

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("sv-stratify-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# SVStratifyDump: structural variants sorted into strata");

        final Path dict = writeDictionary(dir);
        final Path rm = write(dir, "rm.bed", RM_BED);
        final Path sd = write(dir, "sd.bed", SD_BED);
        System.out.printf("track\tRM=%s%n", ReferenceQueryDump.escape(RM_BED));
        System.out.printf("track\tSD=%s%n", ReferenceQueryDump.escape(SD_BED));

        final Path config = write(dir, "strata.tsv", CONFIG);
        System.out.printf("config\tmain=%s%n", ReferenceQueryDump.escape(CONFIG));

        final String vcf = buildVcf();
        final Path input = write(dir, "input.vcf", vcf);
        System.out.printf("vcf\tinput=%s%n", ReferenceQueryDump.escape(vcf));

        final List<String> tracks = List.of(
                "--track-name", "RM", "--track-intervals", rm.toString(),
                "--track-name", "SD", "--track-intervals", sd.toString());

        // The default: one output file, the breakpoint threshold at one and no overlap fraction.
        run(dir, "default", input, dict, config, tracks, List.of(), false);
        // The same with the overlap fraction raised, which is the test the breakpoint count does
        // not make.
        run(dir, "overlap-half", input, dict, config, tracks,
                List.of("--stratify-overlap-fraction", "0.5"), false);
        run(dir, "overlap-full", input, dict, config, tracks,
                List.of("--stratify-overlap-fraction", "1.0"), false);
        // Two breakpoints required, which an insertion still ignores.
        run(dir, "two-breakpoints", input, dict, config, tracks,
                List.of("--stratify-num-breakpoint-overlaps", "2"), false);
        // Multiple matches, allowed and refused.
        run(dir, "multiple-allowed", input, dict, config, tracks,
                List.of("--allow-multiple-matches", "true"), false);
        run(dir, "multiple-refused", input, dict, config, tracks, List.of(), false);
        // One file per stratum, which is the only mode where the stratum's own writer is used.
        // Refused on the first record, which still leaves every file created and headed.
        run(dir, "split", input, dict, config, tracks,
                List.of("--split-output", "true", "--output-prefix", "strat"), true);
        // And the same allowing multiple matches, which is where the records actually land in
        // more than one file.
        run(dir, "split-allowed", input, dict, config, tracks,
                List.of("--split-output", "true", "--output-prefix", "strat",
                        "--allow-multiple-matches", "true"), true);

        // The configuration refusals, each its own table.
        refuse(dir, "reserved-name", input, dict, tracks, String.join("\n",
                "NAME\tSVTYPE\tMIN_SIZE\tMAX_SIZE\tTRACKS",
                "default\tDEL\t-1\t-1\tRM",
                ""));
        refuse(dir, "unknown-track", input, dict, tracks, String.join("\n",
                "NAME\tSVTYPE\tMIN_SIZE\tMAX_SIZE\tTRACKS",
                "DEL_XX\tDEL\t-1\t-1\tXX",
                ""));
        refuse(dir, "negative-min", input, dict, tracks, String.join("\n",
                "NAME\tSVTYPE\tMIN_SIZE\tMAX_SIZE\tTRACKS",
                "DEL_neg\tDEL\t-2\t5000\tRM",
                ""));
        refuse(dir, "min-over-max", input, dict, tracks, String.join("\n",
                "NAME\tSVTYPE\tMIN_SIZE\tMAX_SIZE\tTRACKS",
                "DEL_bad\tDEL\t5000\t1000\tRM",
                ""));
        refuse(dir, "bnd-with-size", input, dict, tracks, String.join("\n",
                "NAME\tSVTYPE\tMIN_SIZE\tMAX_SIZE\tTRACKS",
                "BND_sized\tBND\t100\t5000\tRM",
                ""));
        refuse(dir, "extra-column", input, dict, tracks, String.join("\n",
                "NAME\tSVTYPE\tMIN_SIZE\tMAX_SIZE\tTRACKS\tEXTRA",
                "DEL_x\tDEL\t-1\t-1\tRM\tx",
                ""));
        refuse(dir, "missing-column", input, dict, tracks, String.join("\n",
                "NAME\tSVTYPE\tMIN_SIZE\tMAX_SIZE",
                "DEL_x\tDEL\t-1\t-1",
                ""));

        // A duplicate track name, and a count mismatch between names and files.
        run(dir, "duplicate-track", input, dict, config, List.of(
                "--track-name", "RM", "--track-intervals", rm.toString(),
                "--track-name", "RM", "--track-intervals", sd.toString()), List.of(), false);
        run(dir, "track-count-mismatch", input, dict, config, List.of(
                "--track-name", "RM", "--track-intervals", rm.toString(),
                "--track-name", "SD"), List.of(), false);

        // Both thresholds at zero, which the engine refuses rather than matching everything.
        run(dir, "both-thresholds-zero", input, dict, config, tracks,
                List.of("--stratify-overlap-fraction", "0", "--stratify-num-breakpoint-overlaps",
                        "0"), false);
    }

    /**
     * The measured variants.
     *
     * The deletions differ in where their ends fall and how long they are; the duplication carries
     * no track at all in its stratum; the insertion has no length of its own; and the breakend has
     * neither length nor a second contig on chr1.
     */
    static String buildVcf() {
        final List<String> lines = new ArrayList<>(List.of(
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=100000>",
                "##contig=<ID=chr2,length=100000>",
                "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">",
                "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length\">",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">",
                "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Second contig\">",
                "##INFO=<ID=END2,Number=1,Type=Integer,Description=\"Second position\">",
                "##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Strands\">",
                "##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Algorithms\">",
                "##ALT=<ID=DEL,Description=\"Deletion\">",
                "##ALT=<ID=DUP,Description=\"Duplication\">",
                "##ALT=<ID=INS,Description=\"Insertion\">",
                "##ALT=<ID=BND,Description=\"Breakend\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"));
        // Both ends in RM, 4001 long: small, and its interval covers only a little of any track.
        lines.add(record("del_small_rm", 1000, "DEL", 5000, 4001, null, null));
        // Both ends in RM as well but long enough for the large stratum.
        lines.add(record("del_large_rm", 1000, "DEL", 9000, 8001, null, null));
        // One end in RM and the other in SD, which is the record that matches two strata.
        lines.add(record("del_both", 1000, "DEL", 2000, 1001, null, null));
        // Neither end in any track.
        lines.add(record("del_no_track", 20000, "DEL", 21000, 1001, null, null));
        // A duplication, whose stratum names no track at all.
        lines.add(record("dup_any", 30000, "DUP", 31000, 1001, null, null));
        // An insertion, which has one locus and no length.
        lines.add("chr1\t1000\tins_rm\tN\t<INS>\t.\t.\t"
                + "SVTYPE=INS;END=1000;SVLEN=-1;ALGORITHMS=depth");
        // A breakend from chr1 into chr2, which has no length either.
        lines.add("chr1\t1000\tbnd_rm\tN\t<BND>\t.\t.\t"
                + "SVTYPE=BND;END=1000;CHR2=chr2;END2=5000;STRANDS=+-;ALGORITHMS=depth");
        lines.add("");
        return String.join("\n", lines);
    }

    static String record(final String id, final int start, final String type, final int end,
                         final int length, final String contig2, final Integer end2) {
        // ALGORITHMS is not optional: SVCallRecordUtils.create refuses a record without it.
        final StringBuilder info = new StringBuilder("SVTYPE=" + type + ";END=" + end
                + ";SVLEN=" + length + ";ALGORITHMS=depth");
        if (contig2 != null) {
            info.append(";CHR2=").append(contig2).append(";END2=").append(end2);
        }
        return "chr1\t" + start + "\t" + id + "\tN\t<" + type + ">\t.\t.\t" + info;
    }

    static void run(final Path dir, final String label, final Path input, final Path dict,
                    final Path config, final List<String> tracks, final List<String> extra,
                    final boolean split) throws Exception {
        final Path output = split ? dir.resolve("split-" + label) : dir.resolve(label + ".vcf");
        if (split) {
            Files.createDirectories(output);
        }
        final List<String> argv = new ArrayList<>(List.of(
                "-V", input.toString(),
                "-O", output.toString(),
                "--sequence-dictionary", dict.toString(),
                "--stratify-config", config.toString()));
        argv.addAll(tracks);
        argv.addAll(extra);
        invoke(dir, label, argv);
        report(dir, label, output, split);
    }

    static void refuse(final Path dir, final String label, final Path input, final Path dict,
                       final List<String> tracks, final String table) throws Exception {
        final Path config = write(dir, label + ".tsv", table);
        System.out.printf("config\t%s=%s%n", label, ReferenceQueryDump.escape(table));
        run(dir, label, input, dict, config, tracks, List.of(), false);
    }

    static void invoke(final Path dir, final String label, final List<String> argv) {
        try {
            new SVStratify().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(cause.getMessage()), dir)));
        }
    }

    static void report(final Path dir, final String label, final Path output, final boolean split)
            throws Exception {
        if (!split) {
            if (Files.exists(output)) {
                System.out.printf("out\t%s\t%s=%s%n", label, output.getFileName(),
                        ReferenceQueryDump.escape(masked(body(output), dir)));
            }
            return;
        }
        try (final Stream<Path> entries = Files.list(output)) {
            for (final Path entry : entries.sorted(Comparator.comparing(Path::toString)).toList()) {
                if (entry.toString().endsWith(".vcf.gz")) {
                    System.out.printf("out\t%s\t%s=%s%n", label, entry.getFileName(),
                            ReferenceQueryDump.escape(masked(body(entry), dir)));
                }
            }
        }
    }

    /** A VCF without its header, which is where the stratum annotation lands. */
    static String body(final Path path) throws Exception {
        final String text = path.toString().endsWith(".gz")
                ? new String(new java.util.zip.GZIPInputStream(
                        Files.newInputStream(path)).readAllBytes(), StandardCharsets.UTF_8)
                : Files.readString(path);
        final StringBuilder out = new StringBuilder();
        for (final String line : text.split("\n", -1)) {
            if (!line.startsWith("##") && !line.isEmpty()) {
                out.append(line).append("\n");
            }
        }
        return out.toString();
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static Path writeDictionary(final Path dir) throws Exception {
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 100000),
                new SAMSequenceRecord("chr2", 100000)));
        final Path path = dir.resolve("reference.dict");
        final htsjdk.samtools.SAMFileHeader header = new htsjdk.samtools.SAMFileHeader();
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

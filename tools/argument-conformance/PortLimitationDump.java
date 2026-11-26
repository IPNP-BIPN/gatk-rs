/*
 * What GATK does with the tool and argument combinations the port used to refuse as its own
 * limitation, taken from the reference.
 *
 * Each case is one command line handed to `Main.instanceMain`, run from a directory of its own so
 * the files it leaves are the case's alone. The golden carries the command line itself, so the
 * port's test replays exactly what ran here and a new case is a new line in `main` and nothing
 * else.
 *
 * Output:
 *
 *     cwd\t<the directory the cases ran from>
 *     input\t<input path>\t<its bytes, base64>
 *     directory\t<an input that is a directory>
 *     argv\t<case>\t<argument>\t<argument>...
 *     error\t<case>\t<the exception class>:<its message>
 *     kind\t<case>\t<commandline, user or other>
 *     dir\t<case>\t<a directory the run left, under the case's>
 *     file\t<case>\t<a file the run left, under the case's>\t<its bytes, base64>
 *
 * Usage: PortLimitationDump
 */

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

public class PortLimitationDump {

    static final Path OUT = Paths.get("port-limitations");
    static final Set<String> INPUTS = new LinkedHashSet<>();

    public static void main(final String[] args) throws Exception {
        System.out.println("# PortLimitationDump: the reference's answer to what the port used to refuse");
        // Some messages name a file by its absolute path; the port's test puts its own directory
        // in place of this one.
        System.out.printf("cwd\t%s%n", Paths.get("").toAbsolutePath());

        // AlleleFrequencyQC's chi-squared over fewer than two bins, which no command line is
        // known to reach: commons-math's own refusal for zero and minus one degrees of freedom.
        for (final int df : new int[] {0, -1}) {
            try {
                new org.apache.commons.math3.distribution.ChiSquaredDistribution(df);
                System.out.printf("probe\tchi-squared%d\tnone%n", df);
            } catch (final Exception e) {
                System.out.printf("probe\tchi-squared%d\t%s:%s%n", df, e.getClass().getName(), e.getMessage());
            }
        }

        // A walker whose `requiresReference()` is true declares -R required, so Barclay refuses a
        // run without it before the tool opens anything; the other inputs need not exist.
        run("svcluster-no-reference", "SVCluster", "-V", "absent.vcf", "--ploidy-table",
                "absent.tsv", "-O", "{out}/out.vcf");
        run("jointgermlinecnvsegmentation-no-reference", "JointGermlineCNVSegmentation", "-V",
                "absent.vcf", "--pedigree", "absent.ped", "-O", "{out}/out.vcf");

        // DownsampleByDuplicateSet reads MI lazily: a read's number when it arrives, its strand
        // only when it starts a set or when the set it joined is closed by a later molecule, and
        // never for the reads of the last set after its first. What a malformed tag throws is the
        // JVM's own exception from the place it is first read, and `closeTool` still closes the
        // writer, so the sets applied before the failure are written.
        downsample("no-reads");
        downsample("mi-missing-first", read(1, 100, null));
        downsample("mi-not-a-number", read(1, 100, "MI:Z:x/A"));
        downsample("mi-no-strand-first", read(1, 100, "MI:Z:1"));
        downsample("mi-slash-only", read(1, 100, "MI:Z:/"));
        downsample("mi-integer", read(1, 100, "MI:i:3"));
        downsample("mi-missing-later", read(1, 100, "MI:Z:1/A"), read(1, 110, "MI:Z:1/B"),
                read(1, 120, null));
        downsample("mi-no-strand-in-a-closed-set", read(1, 100, "MI:Z:1/A"), read(1, 110, "MI:Z:1/B"),
                read(1, 120, "MI:Z:2/A"), read(1, 130, "MI:Z:2"), read(1, 140, "MI:Z:3/A"),
                read(1, 150, "MI:Z:3/B"));
        downsample("mi-no-strand-in-the-last-set", read(1, 100, "MI:Z:1/A"), read(1, 110, "MI:Z:1/B"),
                read(1, 120, "MI:Z:2/A"), read(1, 130, "MI:Z:2"));
        downsample("mi-another-contig", read(1, 100, "MI:Z:1/A"), read(2, 100, "MI:Z:1/B"));
        downsample("mi-not-sorted", read(1, 100, "MI:Z:1/A"), read(1, 110, "MI:Z:1/B"),
                read(1, 120, "MI:Z:3/A"), read(1, 130, "MI:Z:3/B"), read(1, 140, "MI:Z:2/A"));

        // GatherPileupSummaries loads --sequence-dictionary in `onStartup`, through a GATKPath
        // made absolute, before any shard is read.
        Files.createDirectories(Paths.get("inputs", "a-directory.dict"));
        INPUTS.add("inputs/a-directory.dict");
        run("gatherpileupsummaries-missing-dictionary", "GatherPileupSummaries",
                "--sequence-dictionary", "inputs/absent.dict", "-I", "absent.table",
                "-O", "{out}/out.table");
        run("gatherpileupsummaries-dictionary-is-a-directory", "GatherPileupSummaries",
                "--sequence-dictionary", "inputs/a-directory.dict", "-I", "absent.table",
                "-O", "{out}/out.table");

        // RevertSamSpark's output map is read on the driver as a TableFeature stream, every group
        // of the header must be in it, and the files are written in the map's HashMap order, so a
        // group only the map names fails on its null header after the groups before it are
        // written. An illegal --attributes-to-clear fails inside the task that writes the first
        // output: Hadoop's "Job aborted." and its empty .parts directory. With --sanitize an
        // earlier task fails, and its message names the executor's host: no stable bytes there.
        final String reverted = bam("revertsamspark", read(1, 100, "XX:Z:a"), read(1, 110, "XX:Z:b"));
        run("revertsamspark-short-map-row", "RevertSamSpark", "-I", reverted,
                "--output-by-readgroup", "--output-map",
                text("revertsamspark-short-map-row.tsv", "READ_GROUP_ID\tOUTPUT\nrg\n"));
        run("revertsamspark-unknown-group", "RevertSamSpark", "-I", reverted,
                "--output-by-readgroup", "--output-map",
                text("revertsamspark-unknown-group.tsv", "READ_GROUP_ID\tOUTPUT\nother\t"
                        + OUT.resolve("revertsamspark-unknown-group").resolve("other.bam") + "\n"));
        run("revertsamspark-extra-group", "RevertSamSpark", "-I", reverted,
                "--output-by-readgroup", "--output-map",
                text("revertsamspark-extra-group.tsv", "READ_GROUP_ID\tOUTPUT\nrg\t"
                        + OUT.resolve("revertsamspark-extra-group").resolve("rg.bam") + "\nother\t"
                        + OUT.resolve("revertsamspark-extra-group").resolve("other.bam") + "\n"));
        run("revertsamspark-bad-attribute", "RevertSamSpark", "-I", reverted, "-O",
                "{out}/out.bam", "--attributes-to-clear", "XXX");
        run("revertsamspark-bad-attribute-by-group", "RevertSamSpark", "-I", reverted,
                "--output-by-readgroup", "-O", "{out}", "--attributes-to-clear", "XXX");

        // ConvertCountsToDepthFile reads its counts through SimpleCountCodec, whose decode turns
        // a NumberFormatException into MalformedFile and lets every other failure through, after
        // the writer has been opened.
        final String countsHeader = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10000\n@RG\tID:rg\tSM:sample\n"
                + "CONTIG\tSTART\tEND\tCOUNT\n";
        counts("counts-not-a-number", countsHeader + "chr1\t1\t100\t5\nchr1\t101\t200\tx\n");
        // Tribble decodes one record ahead, so the record before the refused line is not written.
        counts("counts-not-a-number-third", countsHeader
                + "chr1\t1\t100\t5\nchr1\t101\t200\t6\nchr1\t201\t300\tx\n");
        counts("counts-short-row", countsHeader + "chr1\t1\t100\n");
        counts("counts-start-after-end", countsHeader + "chr1\t200\t100\t5\n");
        counts("counts-start-zero", countsHeader + "chr1\t0\t100\t5\n");

        // DepthOfCoverage's --include-deletions puts the deletions in the pileup, where they count
        // into the depth and into a D column whatever their quality; --ignore-deletion-sites takes
        // them out of the pileup again but leaves the column, at zero.
        final String fasta = reference();
        final String covered = bam("depthofcoverage",
                cigarRead(1, 100, "10M", "IIIIIIIIII"),
                cigarRead(1, 102, "4M2D6M", "IIII#IIIII"),
                cigarRead(1, 104, "2S3M1D5M", "IIIIIIIIII"),
                cigarRead(1, 105, "3M1I6M", "#IIIIIIIII"));
        for (final String[] variant : new String[][] {
                {"default"}, {"include-deletions", "--include-deletions"},
                {"include-deletions-ignore-sites", "--include-deletions", "--ignore-deletion-sites"},
                {"include-deletions-min-quality", "--include-deletions", "--min-base-quality", "20"}}) {
            final List<String> argv = new ArrayList<>(List.of("DepthOfCoverage", "-R", fasta, "-I",
                    covered, "-L", "chr1:98-118", "-O", "{out}/out", "--omit-locus-table",
                    "--omit-interval-statistics", "--omit-per-sample-statistics"));
            argv.addAll(Arrays.asList(variant).subList(1, variant.length));
            run("depthofcoverage-" + variant[0], argv.toArray(new String[0]));
            argv.add("--print-base-counts");
            run("depthofcoverage-" + variant[0] + "-base-counts", argv.toArray(new String[0]));
        }

        for (final String input : INPUTS) {
            final Path path = Paths.get(input);
            if (Files.isRegularFile(path)) {
                System.out.printf("input\t%s\t%s%n", input, base64(path));
            } else if (Files.isDirectory(path)) {
                System.out.printf("directory\t%s%n", input);
            }
        }
    }

    /** A DownsampleByDuplicateSet case over a BAM of the given reads. */
    static void downsample(final String name, final String... reads) throws Exception {
        final String input = bam(name, reads);
        run(name, "DownsampleByDuplicateSet", "-I", input, "--fraction-to-keep", "1.0",
                "-O", "{out}/out.bam");
    }

    /** One mapped read on chr1 or chr2 in the one read group, its MI tag given whole (or none). */
    static String read(final int contig, final int start, final String tag) {
        return String.join("\t", "r" + contig + "_" + start, "0", "chr" + contig,
                String.valueOf(start), "60", "10M", "*", "0", "0", "ACGTACGTAC", "IIIIIIIIII", "RG:Z:rg")
                + (tag == null ? "" : "\t" + tag);
    }

    /** One mapped read on a contig with its own cigar and qualities, its bases cut to fit. */
    static String cigarRead(final int contig, final int start, final String cigar, final String quals) {
        return String.join("\t", "r" + contig + "_" + start, "0", "chr" + contig,
                String.valueOf(start), "60", cigar, "*", "0", "0",
                "ACGTACGTACGT".substring(0, quals.length()), quals, "RG:Z:rg");
    }

    /** The two contigs of the BAMs' header as a FASTA, with its index and dictionary. */
    static String reference() throws Exception {
        final Path dir = Paths.get("inputs");
        Files.createDirectories(dir);
        final StringBuilder fasta = new StringBuilder();
        final StringBuilder fai = new StringBuilder();
        final StringBuilder dict = new StringBuilder("@HD\tVN:1.6\n");
        for (final String contig : new String[] {"chr1", "chr2"}) {
            fasta.append('>').append(contig).append('\n');
            fai.append(String.format("%s\t10000\t%d\t60\t61%n", contig, fasta.length()));
            for (int at = 0; at < 10000; at++) {
                fasta.append("ACGTTGCA".charAt(at % 8));
                if (at % 60 == 59 || at == 9999) {
                    fasta.append('\n');
                }
            }
            dict.append("@SQ\tSN:").append(contig).append("\tLN:10000\n");
        }
        final String[][] files = {{"reference.fasta", fasta.toString()},
                {"reference.fasta.fai", fai.toString()}, {"reference.dict", dict.toString()}};
        for (final String[] file : files) {
            Files.writeString(dir.resolve(file[0]), file[1], StandardCharsets.US_ASCII);
            INPUTS.add(dir.resolve(file[0]).toString());
        }
        return dir.resolve("reference.fasta").toString();
    }

    /**
     * The reads as a coordinate-sorted BAM, written by GATK's own PrintReads so the compression
     * level is GATK's: htsjdk writing first would fix `Defaults` at 5 for the rest of the JVM.
     */
    static String bam(final String name, final String... reads) throws Exception {
        final Path dir = Paths.get("inputs");
        Files.createDirectories(dir);
        final Path sam = dir.resolve(name + ".sam");
        final StringBuilder text = new StringBuilder(
                "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10000\n@SQ\tSN:chr2\tLN:10000\n@RG\tID:rg\tSM:sample\n");
        for (final String read : reads) {
            text.append(read).append('\n');
        }
        Files.writeString(sam, text.toString(), StandardCharsets.UTF_8);
        final String bam = dir.resolve(name + ".bam").toString();
        final PrintStream original = System.out;
        try {
            System.setOut(new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
            new org.broadinstitute.hellbender.Main().instanceMain(
                    new String[] {"PrintReads", "-I", sam.toString(), "-O", bam});
        } finally {
            System.setOut(original);
        }
        INPUTS.add(bam);
        INPUTS.add(dir.resolve(name + ".bai").toString());
        return bam;
    }

    /** A ConvertCountsToDepthFile case over a counts file of the given text. */
    static void counts(final String name, final String content) throws Exception {
        run(name, "ConvertCountsToDepthFile", "--counts-filename",
                text(name + ".counts.tsv", content), "-O", "{out}/out.rd.txt");
    }

    /** A text input, recorded with the others. */
    static String text(final String name, final String content) throws Exception {
        final Path path = Paths.get("inputs", name);
        Files.createDirectories(path.getParent());
        Files.writeString(path, content, StandardCharsets.UTF_8);
        INPUTS.add(path.toString());
        return path.toString();
    }

    static void run(final String name, final String... command) throws Exception {
        final Path dir = OUT.resolve(name);
        Files.createDirectories(dir);
        final List<String> argv = new ArrayList<>();
        for (final String argument : command) {
            final String resolved = argument.replace("{out}", dir.toString());
            if (resolved.startsWith("fixtures/")) {
                INPUTS.add(resolved);
            }
            argv.add(resolved);
        }
        System.out.printf("argv\t%s\t%s%n", name,
                String.join("\t", argv.stream().map(PortLimitationDump::escape).toList()));

        // The tool's own logging carries a clock; what is measured is the files and the refusal.
        final PrintStream original = System.out;
        try {
            System.setOut(new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
            new org.broadinstitute.hellbender.Main().instanceMain(argv.toArray(new String[0]));
        } catch (final Throwable e) {
            System.setOut(original);
            System.out.printf("error\t%s\t%s%n", name, escape(e.getClass().getName() + ":" + e.getMessage()));
            // Which handler `mainEntry` gives it to, which decides the status and the banner.
            System.out.printf("kind\t%s\t%s%n", name,
                    e instanceof org.broadinstitute.barclay.argparser.CommandLineException ? "commandline"
                            : e instanceof org.broadinstitute.hellbender.exceptions.UserException ? "user"
                            : "other");
        } finally {
            System.setOut(original);
        }
        // Everything the run left, directories included, by its path under the case's directory.
        try (final var entries = Files.walk(dir)) {
            for (final Path path : entries.filter(p -> !p.equals(dir)).sorted().toList()) {
                final String relative = dir.relativize(path).toString();
                if (Files.isDirectory(path)) {
                    System.out.printf("dir\t%s\t%s%n", name, relative);
                } else {
                    System.out.printf("file\t%s\t%s\t%s%n", name, relative, base64(path));
                }
            }
        }
    }

    static String base64(final Path path) throws Exception {
        return Base64.getEncoder().encodeToString(Files.readAllBytes(path));
    }

    static String escape(final String text) {
        return text.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
    }
}

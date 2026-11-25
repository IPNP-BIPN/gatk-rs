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
 *     file\t<case>\t<file name>\t<its bytes, base64>
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
        } catch (final Exception | AssertionError e) {
            System.setOut(original);
            System.out.printf("error\t%s\t%s%n", name, escape(e.getClass().getName() + ":" + e.getMessage()));
        } finally {
            System.setOut(original);
        }
        try (final var entries = Files.list(dir)) {
            for (final Path file : entries.sorted().toList()) {
                System.out.printf("file\t%s\t%s\t%s%n", name, file.getFileName(), base64(file));
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

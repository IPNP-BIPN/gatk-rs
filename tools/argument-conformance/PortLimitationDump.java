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
 *     input\t<fixture path>\t<its bytes, base64>
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

        // A walker whose `requiresReference()` is true declares -R required, so Barclay refuses a
        // run without it before the tool opens anything; the other inputs need not exist.
        run("svcluster-no-reference", "SVCluster", "-V", "absent.vcf", "--ploidy-table",
                "absent.tsv", "-O", "{out}/out.vcf");
        run("jointgermlinecnvsegmentation-no-reference", "JointGermlineCNVSegmentation", "-V",
                "absent.vcf", "--pedigree", "absent.ped", "-O", "{out}/out.vcf");

        for (final String input : INPUTS) {
            final Path path = Paths.get(input);
            if (Files.isRegularFile(path)) {
                System.out.printf("input\t%s\t%s%n", input, base64(path));
            }
        }
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

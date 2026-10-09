import java.io.IOException;
import java.io.PrintWriter;
import java.lang.management.ManagementFactory;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Comparator;
import java.util.stream.Stream;

import org.broadinstitute.hellbender.Main;

/**
 * The reference's steady state: one JVM running the same command line {@code n} times.
 *
 * <p>A cold run pays class loading and JIT warm-up every time, which is what a user invoking GATK
 * once per file pays, and it is not what the reference costs once warm. This calls
 * {@code Main.instanceMain}, the door the dumps use, so the loop stays inside one process: the
 * later iterations are the JIT-compiled tool, and their median is the steady-state number.
 *
 * <p>Usage: {@code Steady <iterations> <report> <Tool> <args...>}. Each iteration writes one line,
 * {@code index wall_ns cpu_ns status}, and the loop stops at the first iteration that throws: a
 * tool that keeps static state between runs is reported as having no steady state rather than
 * timed on a second run that does something else.
 */
public final class Steady {
    private Steady() {}

    public static void main(final String[] argv) throws IOException {
        final int iterations = Integer.parseInt(argv[0]);
        final Path report = Paths.get(argv[1]);
        final String[] args = Arrays.copyOfRange(argv, 2, argv.length);
        final com.sun.management.OperatingSystemMXBean os =
                (com.sun.management.OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean();
        try (PrintWriter out = new PrintWriter(Files.newBufferedWriter(report))) {
            for (int i = 0; i < iterations; i++) {
                // The cold runs start from an empty output directory, so every iteration does too:
                // a tool that refuses an existing output would otherwise fail its second run.
                reset(Paths.get("/work/out"));
                Files.createDirectories(Paths.get("/work/tmp"));
                Files.createDirectories(Paths.get("/work/tmp2"));
                final long wall0 = System.nanoTime();
                final long cpu0 = os.getProcessCpuTime();
                String status = "ok";
                try {
                    new Main().instanceMain(args);
                } catch (final Throwable t) {
                    status = t.getClass().getName();
                }
                final long wall = System.nanoTime() - wall0;
                final long cpu = os.getProcessCpuTime() - cpu0;
                out.println(i + "\t" + wall + "\t" + cpu + "\t" + status);
                out.flush();
                if (!status.equals("ok")) {
                    break;
                }
            }
        }
        // A tool can leave a non-daemon thread behind; the loop's answer is already written.
        System.exit(0);
    }

    private static void reset(final Path dir) throws IOException {
        if (Files.exists(dir)) {
            try (Stream<Path> walk = Files.walk(dir)) {
                walk.sorted(Comparator.reverseOrder())
                        .filter(p -> !p.equals(dir))
                        .forEach(p -> p.toFile().delete());
            }
        }
        Files.createDirectories(dir);
    }
}

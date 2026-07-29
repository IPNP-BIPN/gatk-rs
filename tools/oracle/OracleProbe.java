/*
 * Asserts that this container satisfies gatk-rs's oracle contract, and records what it found.
 *
 * The principle, inherited from htsjdk-rs and picard-rs: a degraded environment must FAIL rather
 * than produce a golden that looks exactly like a good one. Every check exists because its failure
 * mode is silent.
 *
 * Two checks here have no counterpart in the picard-rs probe, because they are properties of GATK
 * rather than of Picard:
 *
 *   - the PairHMM implementation. HaplotypeCaller and Mutect2 default to FASTEST_AVAILABLE, which
 *     resolves per host: an AVX-512 machine, an AVX2 machine and a machine with neither run
 *     different code and can produce different bytes. A golden produced under an unrecorded
 *     resolution describes no machine in particular, so the resolution is recorded here and the
 *     vector library's availability is asserted.
 *   - the `gatk` wrapper. The bit-identity claim is defined against `gatk <Tool>`, because that
 *     fixes the argument parser to Barclay rather than Picard's legacy syntax. An image whose
 *     wrapper cannot run would still be able to produce goldens through the jar directly, under a
 *     different parser, and nothing downstream would notice.
 *
 * Exits 2 on violation, so the Docker build itself cannot complete in a bad environment.
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.TreeSet;

public class OracleProbe {

    private static final String EXPECTED_ARCH = "amd64";
    private static final String EXPECTED_JAVA_MAJOR = "17";
    private static final String EXPECTED_GATK_VERSION = "4.6.2.0";

    /**
     * htsjdk-rs decision 0011: FormatUtil reaches NumberFormat.getNumberInstance(), which takes the
     * default locale, and nothing pins it. Under fr-FR a metrics file has commas for decimal
     * points. GATK writes the same metrics files through the same htsjdk, so the locale is part of
     * this contract too.
     */
    private static final String EXPECTED_LOCALE = "en_US";

    private static final String[] REQUIRED_CPU_FLAGS = {"avx", "avx2", "sse4_2"};

    public static void main(final String[] args) throws Exception {
        final List<String> failures = new ArrayList<>();

        final String arch = System.getProperty("os.arch");
        final String javaVersion = System.getProperty("java.version");
        final String javaVendor = System.getProperty("java.vendor");
        final String javaMajor = javaVersion.split("\\.")[0];
        final Locale locale = Locale.getDefault();
        final String decimalSample = new java.text.DecimalFormat("0.0#####").format(1.0 / 3.0);
        final TreeSet<String> cpuFlags = readCpuFlags();

        // The version GATK reports about itself, not the one the build asked for: a jar from the
        // wrong tag would otherwise be invisible.
        String gatkVersion = "unknown";
        try {
            final Class<?> main = Class.forName("org.broadinstitute.hellbender.Main");
            final String implementation = main.getPackage().getImplementationVersion();
            gatkVersion = implementation == null ? "absent" : implementation;
        } catch (final Throwable t) {
            failures.add("GATK classes are not on the classpath: " + t);
        }

        // Picard is bundled inside the GATK jar and is what `gatk MarkDuplicates` runs.
        //
        // Recorded, not asserted: the local jar is shaded, so every package reports the *jar's*
        // implementation version (4.6.2.0), not the bundled library's. This field therefore says
        // "Picard classes are present", and the authority on which Picard is inside remains
        // build.gradle at the pinned tag, which says 3.4.0. Asserting "3.4.0" here would fail on a
        // correct image, and asserting what it actually prints would assert nothing.
        String picardVersion = "unknown";
        try {
            picardVersion = Class.forName("picard.cmdline.CommandLineProgram")
                    .getPackage().getImplementationVersion();
            if (picardVersion == null) picardVersion = "absent";
        } catch (final Throwable t) {
            failures.add("Picard classes are not on the classpath: " + t);
        }

        // GKL degrades silently. The contract pins the JDK deflater, so what matters is that GKL is
        // present and working: its silent absence would mean the pin is doing nothing.
        boolean gklPresent = false;
        try {
            final Class<?> f = Class.forName("com.intel.gkl.compression.IntelDeflaterFactory");
            gklPresent =
                (Boolean) f.getMethod("usingIntelDeflater").invoke(f.getConstructor().newInstance());
        } catch (final Throwable t) {
            failures.add("Intel GKL is not usable: " + t
                    + ". It degrades silently, which is why this is checked.");
        }

        // The vector PairHMM: available or not, and which resolution FASTEST_AVAILABLE would take.
        // Recorded rather than required, because a machine without it is a valid oracle as long as
        // the goldens say so; what is not valid is not knowing.
        String pairHmmVector = "unavailable";
        try {
            final Class<?> lib =
                Class.forName("com.intel.gkl.pairhmm.IntelPairHmm");
            final Object instance = lib.getConstructor().newInstance();
            pairHmmVector = instance.getClass().getName();
        } catch (final Throwable t) {
            pairHmmVector = "unavailable: " + t.getClass().getSimpleName();
        }
        final String pairHmmResolution =
            cpuFlags.contains("avx512f") ? "AVX512 path available"
                : cpuFlags.contains("avx2") ? "AVX2 path available"
                : "scalar only";

        // The entry point the claim is defined against.
        boolean wrapperPresent = new File("/opt/gatk/gatk").canExecute();
        if (!wrapperPresent) {
            failures.add("the `gatk` wrapper is missing or not executable. The bit-identity claim"
                    + " is defined against `gatk <Tool>`, which fixes the parser to Barclay.");
        }

        if (!EXPECTED_ARCH.equals(arch)) {
            failures.add("os.arch is '" + arch + "', expected '" + EXPECTED_ARCH + "'");
        }
        if (!EXPECTED_JAVA_MAJOR.equals(javaMajor)) {
            failures.add("java major is '" + javaMajor + "', expected '" + EXPECTED_JAVA_MAJOR + "'");
        }
        if (!EXPECTED_GATK_VERSION.equals(gatkVersion)) {
            failures.add("GATK reports version '" + gatkVersion + "', expected '"
                    + EXPECTED_GATK_VERSION + "'");
        }
        if (!EXPECTED_LOCALE.equals(locale.toString())) {
            failures.add("default locale is '" + locale + "', expected '" + EXPECTED_LOCALE
                    + "'. Metrics number formatting is locale-dependent; see htsjdk-rs"
                    + " decision 0011.");
        }
        if (!"0.333333".equals(decimalSample)) {
            failures.add("a decimal formats as '" + decimalSample + "', expected '0.333333'.");
        }
        if (!gklPresent) {
            failures.add("usingIntelDeflater is false. The oracle pins the JDK deflater, but a"
                    + " GKL that cannot load means that pin is untested.");
        }
        for (final String flag : REQUIRED_CPU_FLAGS) {
            if (!cpuFlags.isEmpty() && !cpuFlags.contains(flag)) {
                failures.add("CPU flag '" + flag + "' is absent");
            }
        }

        final StringBuilder json = new StringBuilder();
        json.append("{\n");
        json.append("  \"os_arch\": \"").append(arch).append("\",\n");
        json.append("  \"java_version\": \"").append(javaVersion).append("\",\n");
        json.append("  \"java_vendor\": \"").append(javaVendor).append("\",\n");
        json.append("  \"gatk_version\": \"").append(gatkVersion).append("\",\n");
        json.append("  \"picard_classes_present\": ").append(!"unknown".equals(picardVersion))
            .append(",\n");
        json.append("  \"jar_implementation_version\": \"").append(picardVersion).append("\",\n");
        json.append("  \"default_locale\": \"").append(locale).append("\",\n");
        json.append("  \"decimal_sample\": \"").append(decimalSample).append("\",\n");
        json.append("  \"using_intel_deflater\": ").append(gklPresent).append(",\n");
        json.append("  \"gatk_wrapper\": ").append(wrapperPresent).append(",\n");
        json.append("  \"pairhmm_vector\": \"").append(pairHmmVector.replace("\"", "'")).append("\",\n");
        json.append("  \"pairhmm_resolution\": \"").append(pairHmmResolution).append("\",\n");
        json.append("  \"avx\": ").append(cpuFlags.contains("avx")).append(",\n");
        json.append("  \"avx2\": ").append(cpuFlags.contains("avx2")).append(",\n");
        json.append("  \"avx512f\": ").append(cpuFlags.contains("avx512f")).append(",\n");
        json.append("  \"contract_satisfied\": ").append(failures.isEmpty()).append("\n");
        json.append("}");

        if (!failures.isEmpty()) {
            System.err.println("ORACLE CONTRACT VIOLATED. No golden produced here may be trusted.");
            System.err.println(json);
            for (final String f : failures) System.err.println("  - " + f);
            System.exit(2);
        }
        System.out.println(json);
    }

    private static TreeSet<String> readCpuFlags() {
        final TreeSet<String> flags = new TreeSet<>();
        try (BufferedReader r = new BufferedReader(new FileReader("/proc/cpuinfo"))) {
            String line;
            while ((line = r.readLine()) != null) {
                if (line.startsWith("flags")) {
                    for (final String flag : line.split(":", 2)[1].trim().split("\\s+")) {
                        flags.add(flag);
                    }
                    break;
                }
            }
        } catch (final Exception ignored) {
            // No /proc/cpuinfo means the flag checks are skipped rather than failed: the file is
            // absent on non-Linux hosts, and its absence is not evidence of a bad CPU.
        }
        return flags;
    }
}

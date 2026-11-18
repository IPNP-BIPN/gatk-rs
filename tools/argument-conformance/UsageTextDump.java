/*
 * What `gatk <Tool> --help` prints, for every tool `Main` can name.
 *
 * The usage is a rendering of the same annotations the argument definitions come from: an
 * ordering, a grouping, a wrapping, conditional blocks per plugin descriptor, and possible values
 * that a descriptor rather than an enum answers for. All of it is byte-comparable, and none of it
 * can be measured on a handful of tools: a walker's `--disable-read-filter` lists THAT tool's
 * defaults, an annotation descriptor adds its own conditional block, and a Spark tool prints a
 * BETA banner above its usage line. So every tool is measured, GATK's and the Picard tools GATK
 * dispatches to, through `Main.instanceMain` exactly as the command line reaches it.
 *
 * What this is built to catch:
 *
 *   - THE TEXT IS WHAT `Main` WRITES TO STDERR for `<Tool> --help`, not a parser built on the side;
 *   - A PICARD TOOL PRINTS ITS USAGE THROUGH PICARD'S OWN `instanceMain`, which then returns 1;
 *   - `--disable-read-filter`, `--annotation` and their siblings print the DESCRIPTOR'S allowed
 *     values, which differ per tool and are recorded on their own rows;
 *   - THE MUTEX SENTENCE NAMES THE TARGET'S FIELD, recorded where it is not the long name;
 *   - A DEPRECATED ARGUMENT SAYS SO before its documentation, which no declaration records;
 *   - AND `-h` AND `--help` RENDER THE SAME TEXT.
 *
 * Output, tools sorted by name:
 *
 *     usage\t<tool>\t<what --help wrote to stderr, escaped>
 *     lines\t<tool>\t<how many lines it is>
 *     same\t<tool>\t<whether -h wrote the same>
 *     allowed\t<tool>\t<long name>\t<the descriptor's allowed values, sorted, comma-separated>
 *     field\t<tool>\t<long name>\t<the field name, for a mutex target whose field is not its name>
 *     deprecated\t<tool>\t<long name>\t<the deprecation detail the usage prints before the doc>
 *
 * Usage: UsageTextDump
 */

import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLinePluginDescriptor;
import org.broadinstitute.barclay.argparser.NamedArgumentDefinition;
import org.broadinstitute.hellbender.Main;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

public class UsageTextDump {

    public static void main(final String[] args) throws Exception {
        // What GATK's `Main` sets before anything else: Picard tools parse with Barclay.
        System.setProperty("picard.useLegacyParser", "false");
        for (final var tool : tools().entrySet()) {
            final String name = tool.getKey();
            final String text = help(name, "--help");
            System.out.printf("usage\t%s\t%s%n", name, escape(text));
            System.out.printf("lines\t%s\t%d%n", name, text.split("\n", -1).length);
            System.out.printf("same\t%s\t%s%n", name, text.equals(help(name, "-h")));
            descriptors(name, parser(tool.getValue().getDeclaredConstructor().newInstance()));
        }
    }

    /** Every tool the listing names, by the discovery `extractCommandLineProgram` runs. */
    static TreeMap<String, Class<?>> tools() throws Exception {
        final Main main = new Main();
        final java.lang.reflect.Method packages = Main.class.getDeclaredMethod("getPackageList");
        packages.setAccessible(true);
        @SuppressWarnings("unchecked")
        final List<String> packageList = (List<String>) packages.invoke(main);
        final org.broadinstitute.barclay.argparser.ClassFinder finder =
                new org.broadinstitute.barclay.argparser.ClassFinder();
        for (final String pkg : packageList) {
            finder.find(pkg, picard.cmdline.CommandLineProgram.class);
            finder.find(pkg, org.broadinstitute.hellbender.cmdline.CommandLineProgram.class);
        }
        final TreeMap<String, Class<?>> tools = new TreeMap<>();
        for (final Class<?> type : finder.getClasses()) {
            if (type.equals(org.broadinstitute.hellbender.cmdline.PicardCommandLineProgramExecutor.class)
                    || type.getSimpleName().equals("CommandLineArgumentValidator")
                    || !org.broadinstitute.hellbender.utils.ClassUtils.canMakeInstances(type)) {
                continue;
            }
            final var properties = type.getAnnotation(
                    org.broadinstitute.barclay.argparser.CommandLineProgramProperties.class);
            if (properties == null || properties.omitFromCommandLine()) {
                continue;
            }
            tools.put(type.getSimpleName(), type);
        }
        return tools;
    }

    /** What one spelling of the help argument makes `Main` write to stderr. */
    static String help(final String tool, final String spelling) {
        final ByteArrayOutputStream errBytes = new ByteArrayOutputStream();
        final PrintStream realOut = System.out;
        final PrintStream realErr = System.err;
        try {
            System.setOut(new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
            System.setErr(new PrintStream(errBytes, true, StandardCharsets.UTF_8));
            try {
                new Main().instanceMain(new String[]{tool, spelling});
            } catch (final Exception ignored) {
                // A Picard tool's `instanceMain` returns 1 after the usage, which `Main` turns into
                // a PicardNonZeroExitException; the usage is what is measured here.
            }
        } finally {
            System.setOut(realOut);
            System.setErr(realErr);
        }
        return errBytes.toString(StandardCharsets.UTF_8);
    }

    /** The parser a tool hands out, whichever of the two `CommandLineProgram`s it extends. */
    static CommandLineArgumentParser parser(final Object target) {
        return (CommandLineArgumentParser) (target instanceof picard.cmdline.CommandLineProgram
                ? ((picard.cmdline.CommandLineProgram) target).getCommandLineParser()
                : ((org.broadinstitute.hellbender.cmdline.CommandLineProgram) target)
                        .getCommandLineParser());
    }

    /** The descriptors' allowed values, the mutex sentence's field names, the deprecations. */
    static void descriptors(final String tool, final CommandLineArgumentParser parser)
            throws Exception {
        final java.lang.reflect.Field registry = CommandLineArgumentParser.class
                .getDeclaredField("pluginDescriptors");
        registry.setAccessible(true);
        @SuppressWarnings("unchecked")
        final java.util.Map<String, CommandLinePluginDescriptor<?>> descriptors =
                (java.util.Map<String, CommandLinePluginDescriptor<?>>) registry.get(parser);
        final Set<String> targets = new HashSet<>();
        for (final NamedArgumentDefinition definition : parser.getNamedArgumentDefinitions()) {
            targets.addAll(definition.getMutexTargetList());
        }
        for (final NamedArgumentDefinition definition : parser.getNamedArgumentDefinitions()) {
            final String longName = definition.getLongName();
            // `usage` asks every descriptor about every argument and prints what one answers;
            // a descriptor refuses a name that is not its own, which is no answer.
            for (final CommandLinePluginDescriptor<?> descriptor : descriptors.values()) {
                final Set<String> allowed;
                try {
                    allowed = descriptor.getAllowedValuesForDescriptorHelp(longName);
                } catch (final RuntimeException notItsOwn) {
                    continue;
                }
                if (allowed != null) {
                    // Sorted as the usage prints them, which ignores case: `AlleleFraction`
                    // before `AS_BaseQualityRankSumTest`.
                    final Set<String> sorted = new java.util.TreeSet<>(String.CASE_INSENSITIVE_ORDER);
                    sorted.addAll(allowed);
                    System.out.printf("allowed\t%s\t%s\t%s%n", tool, longName,
                            String.join(",", sorted));
                }
            }
            final String field = definition.getUnderlyingField().getName();
            if (targets.contains(longName) && !field.equals(longName)) {
                System.out.printf("field\t%s\t%s\t%s%n", tool, longName, field);
            }
            if (definition.isDeprecated()) {
                System.out.printf("deprecated\t%s\t%s\t%s%n", tool, longName,
                        escape(definition.getDeprecationDetail()));
            }
        }
    }

    static String escape(final String text) {
        return text.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
    }
}

/*
 * Which plugin owns which argument, taken from the reference.
 *
 * `barclay-plugin-descriptors` measured what the trim DOES: a controlled argument whose plugin
 * nobody named is dropped before the required check, so a required argument of an unselected
 * plugin does not fire. What that suite deliberately did not measure is the table the trim runs
 * over, because it is a property of GATK's filter library rather than of the mechanism: WHICH
 * class each of the controlled arguments belongs to, and which of them are required.
 *
 * Without that table a port has the rule and cannot apply it. It cannot be derived either: the
 * declarations golden records that an argument is controlled by `GATKReadFilterPluginDescriptor`
 * and not which read filter contributed it, and the argument names do not name their filters
 * (`--library` belongs to `LibraryReadFilter`, `--platform-filter-name` to `PlatformReadFilter`,
 * `--black-listed-lanes` to `ReadGroupBlackListReadFilter`).
 *
 * Six things this is built to record.
 *
 *   - THE OWNERSHIP TABLE: one row per controlled argument, with the SIMPLE CLASS NAME of the
 *     instance that declared it, its short name, and whether it is required;
 *   - THAT A REQUIRED CONTROLLED ARGUMENT EXISTS AT ALL, which is what makes the trim load-bearing
 *     rather than decorative;
 *   - THAT A TOOL DEFAULT COUNTS AS SELECTED: a descriptor constructed with default instances
 *     allows their dependent arguments with no `--read-filter` on the command line, which is why
 *     a plain walker command line parses;
 *   - AND THAT A DEFAULT IS STILL SELECTED ONCE DISABLED, or is not: `--disable-read-filter` and
 *     `--disable-tool-default-read-filters` are answered by `isDependentArgumentAllowed` after the
 *     parse, not before it;
 *   - THE REQUIRED CHECK ON A SELECTED PLUGIN, which does fire: naming `LibraryReadFilter` without
 *     `--library` is the error the trim does not swallow;
 *   - AND THE THREE-WAY ANSWER on a required argument that belongs to nobody selected: dropped
 *     when absent, an error naming the CLASS when given.
 *
 * Output:
 *
 *     owner\t<plugin simple name>\t<long name>\t<short name or ->\t<required|optional>
 *     descriptor\t<display name>\t<count of controlled arguments>
 *     allowed\t<label>\t<plugin simple name>\t<true|false>
 *     case\t<label>\t<argv, space separated>
 *     result\t<label>\tok|not-parsed|E:<exception class>:<message>
 *     resolved\t<label>\t<the filters the descriptor resolved, in its order>
 *
 * Usage: PluginArgumentOwnershipDump
 */

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.NamedArgumentDefinition;

import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary.MappedReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class PluginArgumentOwnershipDump {

    /** A stand-in for the tool's own arguments, which no descriptor controls. */
    public static final class Args {
        @Argument(fullName = "tool-arg", optional = true, doc = "the tool's own")
        public String toolArg;
    }

    static CommandLineArgumentParser parser(final Args target, final List<ReadFilter> defaults) {
        return new CommandLineArgumentParser(
                target,
                Collections.singletonList(new GATKReadFilterPluginDescriptor(defaults)),
                Collections.emptySet());
    }

    public static void main(final String[] args) {
        System.out.println("# PluginArgumentOwnershipDump: which plugin owns which argument");

        // The table itself, over a descriptor with no defaults: discovery registers every read
        // filter's arguments whatever the tool asked for, so the ownership does not depend on the
        // defaults and this list is the whole of it.
        final CommandLineArgumentParser parser = parser(new Args(), new ArrayList<ReadFilter>());
        final List<String> rows = new ArrayList<>();
        int controlled = 0;
        for (final NamedArgumentDefinition definition : parser.getNamedArgumentDefinitions()) {
            if (!definition.isControlledByPlugin()) {
                continue;
            }
            controlled++;
            final String owner = definition.getContainingObject().getClass().getSimpleName();
            final String shortName = definition.getShortName() == null
                    || definition.getShortName().isEmpty() ? "-" : definition.getShortName();
            rows.add(String.format("owner\t%s\t%s\t%s\t%s", owner, definition.getLongName(),
                    shortName, definition.isOptional() ? "optional" : "required"));
        }
        Collections.sort(rows);
        rows.forEach(System.out::println);
        System.out.printf("descriptor\t%s\t%d%n",
                new GATKReadFilterPluginDescriptor(new ArrayList<ReadFilter>()).getDisplayName(),
                controlled);

        // A tool default counts as selected, and a filter nobody named does not. The question is
        // asked of the descriptor rather than of the parser, which is where the trim asks it.
        allowed("no-defaults", new ArrayList<ReadFilter>());
        allowed("wellformed-and-mapped",
                List.of(new WellformedReadFilter(), ReadFilterLibrary.MAPPED));

        // The trim, over the arguments the table above says are required.
        run("nothing-named", new ArrayList<ReadFilter>());
        run("a-required-argument-of-nobody", new ArrayList<ReadFilter>(), "--library", "lib1");
        run("its-filter-named", new ArrayList<ReadFilter>(),
                "--read-filter", "LibraryReadFilter", "--library", "lib1");
        run("its-filter-named-without-it", new ArrayList<ReadFilter>(),
                "--read-filter", "LibraryReadFilter");
        run("a-default-filter-argument", List.of(new WellformedReadFilter(), ReadFilterLibrary.MAPPED),
                "--ambig-filter-frac", "0.1");
        run("a-default-filters-own-argument",
                List.of(new WellformedReadFilter(), ReadFilterLibrary.MAPPED),
                "--read-filter", "AmbiguousBaseReadFilter", "--ambig-filter-frac", "0.1");
        run("the-defaults-disabled", List.of(new WellformedReadFilter(), ReadFilterLibrary.MAPPED),
                "--disable-tool-default-read-filters", "true");
        run("one-default-disabled", List.of(new WellformedReadFilter(), ReadFilterLibrary.MAPPED),
                "--disable-read-filter", "MappedReadFilter");
        run("a-second-required-argument", new ArrayList<ReadFilter>(),
                "--read-filter", "PlatformReadFilter", "--platform-filter-name", "ILLUMINA");
        run("a-second-required-argument-missing", new ArrayList<ReadFilter>(),
                "--read-filter", "PlatformReadFilter");
    }

    static void allowed(final String label, final List<ReadFilter> defaults) {
        final GATKReadFilterPluginDescriptor descriptor =
                new GATKReadFilterPluginDescriptor(defaults);
        for (final Class<?> plugin : List.of(WellformedReadFilter.class, MappedReadFilter.class,
                org.broadinstitute.hellbender.engine.filters.LibraryReadFilter.class)) {
            System.out.printf("allowed\t%s\t%s\t%s%n", label, plugin.getSimpleName(),
                    descriptor.isDependentArgumentAllowed(plugin));
        }
    }

    static void run(final String label, final List<ReadFilter> defaults, final String... argv) {
        System.out.printf("case\t%s\t%s%n", label, String.join(" ", argv));
        String result;
        GATKReadFilterPluginDescriptor descriptor = new GATKReadFilterPluginDescriptor(defaults);
        try {
            final PrintStream sink = new PrintStream(new ByteArrayOutputStream());
            final CommandLineArgumentParser parser = new CommandLineArgumentParser(
                    new Args(), Collections.singletonList(descriptor), Collections.emptySet());
            result = parser.parseArguments(sink, argv) ? "ok" : "not-parsed";
        } catch (final Exception | AssertionError e) {
            result = "E:" + e.getClass().getName() + ":"
                    + String.valueOf(e.getMessage()).replace("\n", "\\n");
        }
        System.out.printf("result\t%s\t%s%n", label, result);
        if (result.equals("ok")) {
            final List<String> names = new ArrayList<>();
            for (final ReadFilter filter : descriptor.getResolvedInstances()) {
                names.add(filter.getClass().getSimpleName());
            }
            System.out.printf("resolved\t%s\t%s%n", label, String.join(" ", names));
        }
    }
}

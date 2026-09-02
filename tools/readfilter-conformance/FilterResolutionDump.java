/*
 * `GATKReadFilterPluginDescriptor`: which read filters a command line actually ends up with.
 *
 * The runner reads --read-filter and --disable-tool-default-read-filters and ignores the other two
 * of the four the descriptor owns: --disable-read-filter and --inverted-read-filter change nothing
 * in the port and change the filters, and therefore the count, in the reference (gatk-rs#69).
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE ORDER IS DEFAULTS FIRST, then the user's enabled filters in the order given, then the
 *     inverted ones; disabling removes from the DEFAULTS only, before any of that;
 *   - AN ENABLED FILTER THAT IS ALREADY A DEFAULT IS NOT ADDED TWICE, because the list is checked
 *     with contains() -- on filter INSTANCES, which is why the answer depends on their equality
 *     and not on their names;
 *   - AN INVERTED FILTER IS ALWAYS ADDED, defaults or not, and it is added as a NEGATED instance
 *     rather than replacing anything, so inverting a tool default leaves both in the list;
 *   - --disable-tool-default-read-filters DROPS THE DEFAULTS WHOLESALE, and a --disable-read-filter
 *     naming one of them is then redundant rather than an error;
 *   - ENABLING OR DISABLING THE SAME FILTER TWICE IS A REFUSAL, and the two messages differ;
 *   - DISABLING A FILTER THAT DOES NOT EXIST IS A REFUSAL, where disabling one the tool never
 *     enabled is only a warning;
 *   - INVERTING A TOOL DEFAULT IS A REFUSAL OF ITS OWN, unless the defaults were disabled, "so we
 *     do not inadvertently filter all reads from the input";
 *   - AND THE "enabled and inverted" REFUSAL PRINTS THE WRONG VARIABLE: the reference builds
 *     `enabledAndInverted` and then formats `enabledAndDisabled`, so the message lists the empty
 *     set. That is the reference's behaviour and a port has to reproduce it.
 *
 * Output:
 *
 *     filters\t<case>\t<class simple names, space separated, or (empty)>
 *     error\t<case>\t<exception class>: <message>
 *
 * Usage: FilterResolutionDump
 */

import org.broadinstitute.hellbender.cmdline.GATKPlugin.DefaultGATKReadFilterArgumentCollection;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class FilterResolutionDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** The defaults every case starts from: what a read walker takes. */
    static List<ReadFilter> toolDefaults() {
        return new ArrayList<>(Arrays.asList(
                ReadFilterLibrary.MAPPED,
                ReadFilterLibrary.NOT_DUPLICATE,
                ReadFilterLibrary.PRIMARY_LINE));
    }

    static List<String> list(final String... values) {
        return new ArrayList<>(Arrays.asList(values));
    }

    /**
     * One command line's worth of filter arguments, resolved.
     *
     * A negated filter is a DIFFERENT INSTANCE and not a flag on one: `negate()` answers a
     * `ReadFilterNegate`, whose class name is what the list then carries. So the resolved list can
     * hold a filter and its negation at once, and does when the two are asked for.
     */
    static void run(final String name, final List<String> enabled, final List<String> disabled,
                    final List<String> inverted, final boolean disableDefaults) {
        try {
            final DefaultGATKReadFilterArgumentCollection userArgs =
                    new DefaultGATKReadFilterArgumentCollection();
            userArgs.userEnabledReadFilterNames.addAll(enabled);
            userArgs.userDisabledReadFilterNames.addAll(disabled);
            userArgs.userEnabledInvertedReadFilterNames.addAll(inverted);
            userArgs.disableToolDefaultReadFilters = disableDefaults;

            final GATKReadFilterPluginDescriptor descriptor =
                    new GATKReadFilterPluginDescriptor(userArgs, toolDefaults());
            // The descriptor has to see the library before it can resolve a name to an instance.
            // The parser does this while it walks the plugin packages; here it is done the same
            // way, because `allDiscoveredReadFilters` is otherwise empty and every name refuses.
            final org.broadinstitute.barclay.argparser.ClassFinder finder =
                    new org.broadinstitute.barclay.argparser.ClassFinder();
            for (final String pkg : descriptor.getPackageNames()) {
                finder.find(pkg, descriptor.getPluginBaseClass());
            }
            for (final Class<?> type : finder.getConcreteClasses()) {
                if (descriptor.includePluginClass(type)) {
                    descriptor.createInstanceForPlugin(type);
                }
            }
            descriptor.validateAndResolvePlugins();

            final List<String> names = new ArrayList<>();
            for (final ReadFilter filter : descriptor.getResolvedInstances()) {
                final String simple = filter.getClass().getSimpleName();
                names.add(simple.isEmpty() ? "(negated)" : simple);
            }
            emit("filters", name, names.isEmpty() ? "(empty)" : String.join(" ", names));
        } catch (final Exception e) {
            emit("error", name, e.getClass().getName() + ": " + e.getMessage());
        }
    }

    static final List<String> NONE = list();

    public static void main(final String[] args) {
        // The baseline: the tool's own defaults, in the order it declares them.
        run("defaults-only", NONE, NONE, NONE, false);
        // One more filter, appended after the defaults.
        run("one-enabled", list("GoodCigarReadFilter"), NONE, NONE, false);
        // Two more, in the order given rather than alphabetically.
        run("two-enabled", list("GoodCigarReadFilter", "FirstOfPairReadFilter"), NONE, NONE, false);
        // A filter that is already a default, which is not added twice.
        run("enabled-is-a-default", list("MappedReadFilter"), NONE, NONE, false);

        // Disabling a default removes it, and the rest keep their order.
        run("disable-a-default", NONE, list("MappedReadFilter"), NONE, false);
        run("disable-two-defaults", NONE, list("MappedReadFilter", "PrimaryLineReadFilter"), NONE, false);
        // Disabling one the tool never enabled is a WARNING, not a refusal.
        run("disable-a-non-default", NONE, list("GoodCigarReadFilter"), NONE, false);
        // Disabling a name no filter carries is a refusal, and so is ENABLING one, with a
        // different message and a different exception class.
        run("disable-unknown", NONE, list("NoSuchReadFilter"), NONE, false);
        run("enabled-unknown", list("NoSuchReadFilter"), NONE, NONE, false);
        run("inverted-unknown", NONE, NONE, list("NoSuchReadFilter"), false);

        // Inverting: always added, and added NEGATED rather than replacing anything.
        run("invert-a-non-default", NONE, NONE, list("GoodCigarReadFilter"), false);
        run("invert-a-default", NONE, NONE, list("MappedReadFilter"), false);
        // Enabled and inverted at once, which is the refusal whose message names the wrong set.
        run("enabled-and-inverted", list("GoodCigarReadFilter"), NONE, list("GoodCigarReadFilter"), false);

        // The defaults dropped wholesale, with and without a filter of the user's own.
        run("no-defaults", NONE, NONE, NONE, true);
        run("no-defaults-one-enabled", list("GoodCigarReadFilter"), NONE, NONE, true);
        // Disabling a default that is already gone.
        run("no-defaults-and-disable", NONE, list("MappedReadFilter"), NONE, true);

        // The three duplicate refusals.
        run("enabled-twice", list("GoodCigarReadFilter", "GoodCigarReadFilter"), NONE, NONE, false);
        run("disabled-twice", NONE, list("MappedReadFilter", "MappedReadFilter"), NONE, false);
        run("enabled-and-disabled", list("GoodCigarReadFilter"), list("GoodCigarReadFilter"), NONE, false);

        System.out.print(buf);
    }
}

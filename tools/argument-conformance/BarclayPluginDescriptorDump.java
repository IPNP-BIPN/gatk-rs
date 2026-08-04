/*
 * Barclay's plugin descriptors, taken from the reference.
 *
 * This is the mechanism behind `--read-filter` and `--annotation`: arguments that do not exist
 * until another argument names the thing that owns them. A descriptor declares a base class and a
 * set of packages; the parser finds every implementation, instantiates it, and registers its
 * `@Argument` fields as *controlled by* that descriptor. So `--minimum-mapping-quality` is
 * registered whether or not anybody asked for `MappingQualityReadFilter`.
 *
 * What decides whether it is usable is `validatePluginArgumentValues`, which runs before the
 * required check and rewrites the definition list:
 *
 *   - a controlled argument whose predecessor was NOT named and which was NOT given is DROPPED. It
 *     is not merely optional: it leaves the list, so a REQUIRED argument belonging to a plugin
 *     nobody selected does not fire. That is the row worth having;
 *   - a controlled argument whose predecessor was not named and which WAS given is an error naming
 *     the predecessor's SIMPLE CLASS NAME, not the argument that would have selected it;
 *   - the error is built from getShortName() + "/" + getLongName(), so an argument with no short
 *     name reports a leading slash, the same way a tag on a non-taggable argument does.
 *
 * The descriptor here is GATK's own `GATKReadFilterPluginDescriptor` over GATK's own read
 * filters, not a mock: `--read-filter` IS this mechanism, so measuring it with a stand-in would be
 * measuring the stand-in. That also means the discovery half runs for real — `ClassFinder` scans
 * `org.broadinstitute.hellbender.engine.filters` in the pinned jar.
 *
 * Output:
 *
 *     defs\t<label>\t<index>\t<long name>
 *     case\t<label>\t<argv, space separated>
 *     result\t<label>\tok|E:<exception class>:<message>
 *     field\t<label>\t<long name>\t<value>
 *
 * Usage: BarclayPluginDescriptorDump
 */

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.NamedArgumentDefinition;

import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class BarclayPluginDescriptorDump {

    /** The tool's own arguments, which are not controlled by any descriptor. */
    public static final class Args {
        @Argument(fullName = "tool-arg", optional = true, doc = "the tool's own")
        public String toolArg;
    }

    public static void main(final String[] args) {
        System.out.println("# BarclayPluginDescriptorDump: plugin descriptors");

        // Three arguments of three read filters, chosen because their shapes differ: one required,
        // one optional with a short name, one optional without.
        definitions("discovered");

        // No filter named. Every controlled argument is dropped before the required check, so
        // AmbiguousBaseReadFilter's required threshold does not fire.
        run("no-filter");
        run("tool-arg-only", "--tool-arg", "t");

        // A dependent argument without its filter: the error names the CLASS.
        run("dependent-without-filter", "--ambig-filter-frac", "0.1");
        run("dependent-without-filter-mapping-quality", "--minimum-mapping-quality", "30");

        // The filter named, so its argument is allowed.
        run("filter-then-dependent",
                "--read-filter", "AmbiguousBaseReadFilter", "--ambig-filter-frac", "0.1");
        run("filter-without-dependent", "--read-filter", "AmbiguousBaseReadFilter");
        run("mapping-quality-filter",
                "--read-filter", "MappingQualityReadFilter", "--minimum-mapping-quality", "30");

        // Two filters at once.
        run("two-filters",
                "--read-filter", "AmbiguousBaseReadFilter",
                "--read-filter", "MappingQualityReadFilter",
                "--minimum-mapping-quality", "30");
        // One named, the other's argument given.
        run("one-filter-other-dependent",
                "--read-filter", "AmbiguousBaseReadFilter", "--minimum-mapping-quality", "30");

        // A name no filter has.
        run("unknown-filter-name", "--read-filter", "NoSuchReadFilter");
    }

    /** The arguments this dump reports on, and the filter each belongs to. */
    static final String[] REPORTED = {
        "read-filter",
        "ambig-filter-frac",
        "minimum-mapping-quality",
        "tool-arg",
    };

    static CommandLineArgumentParser parser(final Args target) {
        return new CommandLineArgumentParser(
                target,
                Collections.singletonList(
                        new GATKReadFilterPluginDescriptor(new ArrayList<ReadFilter>())),
                Collections.emptySet());
    }

    /**
     * Which of the reported arguments the parser holds, and whether each is controlled.
     *
     * Not the whole list: a real descriptor registers every read filter's arguments at once, and
     * that count is a property of the filter library rather than of this mechanism.
     */
    static void definitions(final String label) {
        try {
            final List<NamedArgumentDefinition> defs = parser(new Args()).getNamedArgumentDefinitions();
            int index = 0;
            for (final String name : REPORTED) {
                for (final NamedArgumentDefinition def : defs) {
                    if (def.getLongName().equals(name)) {
                        System.out.printf("defs\t%s\t%d\t%s\t%s%n", label, index++, name,
                                def.isControlledByPlugin()
                                        ? def.getContainingObject().getClass().getSimpleName()
                                        : "tool");
                    }
                }
            }
        } catch (final Exception | AssertionError e) {
            System.out.printf("defs\t%s\tE:%s:%s%n", label, e.getClass().getName(), e.getMessage());
        }
    }

    static void run(final String label, final String... argv) {
        final Args target = new Args();
        System.out.printf("case\t%s\t%s%n", label, String.join(" ", argv));

        String result;
        List<NamedArgumentDefinition> defs = null;
        try {
            final PrintStream sink = new PrintStream(new ByteArrayOutputStream());
            final CommandLineArgumentParser parser = parser(target);
            defs = parser.getNamedArgumentDefinitions();
            result = parser.parseArguments(sink, argv) ? "ok" : "not-parsed";
        } catch (final Exception | AssertionError e) {
            result = "E:" + e.getClass().getName() + ":"
                    + String.valueOf(e.getMessage()).replace("\n", "\\n");
        }
        System.out.printf("result\t%s\t%s%n", label, result);

        if (result.equals("ok")) {
            for (final String name : REPORTED) {
                for (final NamedArgumentDefinition def : defs) {
                    if (def.getLongName().equals(name)) {
                        System.out.printf("field\t%s\t%s\t%s%n",
                                label, name, String.valueOf(def.getArgumentValue()));
                    }
                }
            }
        }
    }
}

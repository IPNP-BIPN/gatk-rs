/*
 * The read filter catalogue and each tool's defaults, taken from the reference.
 *
 * Two sets a walker's usage prints and a port cannot derive.
 *
 *   - THE CATALOGUE: every read filter `ClassFinder` discovered, which is what `--read-filter`
 *     lists as its possible values. It is the filter LIBRARY, and a filter with no arguments is in
 *     it and is therefore not in the ownership table;
 *   - AND THE TOOL'S OWN DEFAULTS: what `getDefaultReadFilters()` returned, which is what
 *     `--disable-read-filter` lists as its possible values. That set is per tool: `CountReads`
 *     takes `ReadWalker`'s one filter and `ApplyBQSR` takes the same, while a tool that is no
 *     walker has no descriptor at all.
 *
 * The second set is the one the trim is still short of. A default counts as SELECTED, so a default
 * filter's own arguments are accepted with no `--read-filter` on the command line, and a port
 * without the list refuses a command line the reference takes.
 *
 * Both are read off the descriptor the tool's own parser built, through
 * `getAllowedValuesForDescriptorHelp`, which is the same call the usage rendering makes. Sets have
 * no order, so both are sorted here; the usage sorts them too, and `usage-text` is where that
 * order is compared.
 *
 * Output:
 *
 *     catalogue\t<count>\t<every filter name, space separated, sorted>
 *     descriptor\t<tool>\t<the descriptor's display name, or none>
 *     defaults\t<tool>\t<the tool's default filter names, space separated, sorted>
 *     allowed\t<tool>\t<argument long name>\t<its possible values, space separated, sorted>
 *
 * Usage: ReadFilterCatalogueDump
 */

import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

public class ReadFilterCatalogueDump {

    public static void main(final String[] args) {
        System.out.println("# ReadFilterCatalogueDump: the filter library and each tool's defaults");

        // The catalogue, off the first tool that has a descriptor at all: discovery is the
        // descriptor's and does not depend on which tool asked.
        final GATKReadFilterPluginDescriptor first =
                descriptor(new org.broadinstitute.hellbender.tools.CountReads());
        final List<String> catalogue =
                sorted(first.getAllowedValuesForDescriptorHelp("read-filter"));
        System.out.printf("catalogue\t%d\t%s%n", catalogue.size(), String.join(" ", catalogue));

        // The nine tools whose declarations this port carries, walkers and not.
        tool("CountReads", new org.broadinstitute.hellbender.tools.CountReads());
        tool("CountVariants", new org.broadinstitute.hellbender.tools.walkers.CountVariants());
        tool("PrintReads", new org.broadinstitute.hellbender.tools.PrintReads());
        tool("ApplyBQSR", new org.broadinstitute.hellbender.tools.walkers.bqsr.ApplyBQSR());
        tool("SelectVariants",
                new org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants());
        tool("IndexFeatureFile", new org.broadinstitute.hellbender.tools.IndexFeatureFile());
        tool("GatherVcfsCloud", new org.broadinstitute.hellbender.tools.GatherVcfsCloud());
        tool("PrintBGZFBlockInformation",
                new org.broadinstitute.hellbender.tools.PrintBGZFBlockInformation());
        tool("CreateHadoopBamSplittingIndex",
                new org.broadinstitute.hellbender.tools.spark.CreateHadoopBamSplittingIndex());
    }

    /** The read filter descriptor the tool's own parser built, or null where it built none. */
    static GATKReadFilterPluginDescriptor descriptor(final CommandLineProgram target) {
        final CommandLineArgumentParser parser =
                (CommandLineArgumentParser) target.getCommandLineParser();
        return parser.getPluginDescriptor(GATKReadFilterPluginDescriptor.class);
    }

    static void tool(final String name, final CommandLineProgram target) {
        final GATKReadFilterPluginDescriptor plugin = descriptor(target);
        if (plugin == null) {
            System.out.printf("descriptor\t%s\tnone%n", name);
            return;
        }
        System.out.printf("descriptor\t%s\t%s%n", name, plugin.getDisplayName());
        // `--disable-read-filter`'s possible values ARE the tool's defaults, which is how the
        // usage prints them and the only place the list is visible from outside the tool.
        final List<String> defaults =
                sorted(plugin.getAllowedValuesForDescriptorHelp("disable-read-filter"));
        System.out.printf("defaults\t%s\t%s%n", name, String.join(" ", defaults));
        for (final String argument : new String[]{"read-filter", "disable-read-filter"}) {
            final List<String> allowed = sorted(plugin.getAllowedValuesForDescriptorHelp(argument));
            System.out.printf("allowed\t%s\t%s\t%s%n", name, argument, String.join(" ", allowed));
        }
    }

    static List<String> sorted(final Set<String> values) {
        final List<String> list = new ArrayList<>(values == null ? Set.of() : values);
        Collections.sort(list);
        return list;
    }
}

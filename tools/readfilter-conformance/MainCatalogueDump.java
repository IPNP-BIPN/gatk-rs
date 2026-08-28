/*
 * Main's own tool catalogue and the suggestions it makes over it, taken from the reference.
 *
 * MainDispatchDump measured the SEARCH over a catalogue of five named classes. This measures the
 * CATALOGUE: the set of tool names `gatk <Tool>` will resolve, discovered exactly the way Main
 * discovers it, and what the suggestion search answers when a near miss is asked of the whole set
 * rather than of a handful.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE CATALOGUE IS THE TWO PACKAGES `org.broadinstitute.hellbender` AND `picard`, scanned for
 *     both flavours of CommandLineProgram, so the Picard tools GATK dispatches to are in it and
 *     `SortVcf` resolves;
 *   - THE CATALOGUE IS BIGGER THAN THE DOCUMENTED TOOL LIST: the scan finds 331 classes where the
 *     CLI reports 311 tools, so twenty names resolve that no inventory of the documentation names;
 *   - A CLASS WITHOUT THE ANNOTATION IS NOT IN IT, and PicardCommandLineProgramExecutor and
 *     CommandLineArgumentValidator are excluded by name whatever their annotation says;
 *   - THE KEY IS THE SIMPLE CLASS NAME, so two tools of the same simple name in different
 *     packages would be a collision the reference refuses to start with;
 *   - A DEPRECATED TOOL IS NOT IN THE CATALOGUE AT ALL. CNNScoreVariants has been taken out of the
 *     code, so its name falls through to the search, and the deprecation registry answering before
 *     the search is the ONLY reason the notice is ever seen;
 *   - THE SPARK TOOLS ARE IN IT, so `PrintReadsSpark` resolves as `PrintReads` does;
 *   - A NAME ONE CHARACTER OFF A REAL TOOL FINDS IT AND ONLY IT, over three hundred tools and not
 *     only over a handful;
 *   - A PREFIX OF SEVERAL TOOLS FINDS ALL OF THEM at distance zero, and a substring of five
 *     characters or more does the same wherever it sits: `Fingerprint` finds five;
 *   - AND A NAME FAR FROM EVERYTHING FINDS NOTHING even against the whole catalogue.
 *
 * The suggestion names are emitted SORTED and the raw message only where there is nothing to
 * suggest. The reference walks its class set in hash order, which is not something a byte
 * comparison should be asked to reproduce; what is being pinned here is WHICH tools a query
 * finds, not the order the message lists them in.
 *
 * Output:
 *
 *     count\tcatalogue\t<how many tools the scan found>
 *     catalogue\t<n>\t<the names, sorted, comma-separated, one line per hundred>
 *     suggests\t<query>\t<the names it suggests, sorted, comma-separated>
 *     message\t<query>\t<the whole message, where nothing is suggested>
 *     holds\t names\t<name>=<whether the catalogue holds it>, comma-separated>
 *     error\t<query>\t<exception class>:<message>
 *
 * Usage: MainCatalogueDump
 */

import org.broadinstitute.barclay.argparser.ClassFinder;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.utils.ClassUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

public class MainCatalogueDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** Main's own package list, which is protected and therefore reachable from a subclass. */
    static final class Reachable extends Main {
        List<String> packages() {
            return getPackageList();
        }

        List<Class<? extends CommandLineProgram>> classes() {
            return getClassList();
        }
    }

    /**
     * `extractCommandLineProgram`'s scan, transcribed: both flavours of CommandLineProgram in
     * every package, keeping the instantiable classes that carry the annotation and dropping the
     * two the reference names.
     */
    static Set<Class<?>> catalogue(final Reachable main) {
        final ClassFinder finder = new ClassFinder();
        for (final String pkg : main.packages()) {
            finder.find(pkg, picard.cmdline.CommandLineProgram.class);
            finder.find(pkg, CommandLineProgram.class);
        }
        final Set<Class<?>> toCheck = new LinkedHashSet<>(finder.getClasses());
        toCheck.addAll(main.classes());
        final Set<Class<?>> kept = new LinkedHashSet<>();
        for (final Class<?> clazz : toCheck) {
            if (clazz.getSimpleName().equals("PicardCommandLineProgramExecutor")
                    || clazz.getSimpleName().equals("CommandLineArgumentValidator")) {
                continue;
            }
            if (!ClassUtils.canMakeInstances(clazz)) {
                continue;
            }
            final CommandLineProgramProperties property = Main.getProgramProperty(clazz);
            if (property == null) {
                continue;
            }
            kept.add(clazz);
        }
        return kept;
    }

    /** The names a message suggests, sorted, or nothing where it suggests none. */
    static List<String> suggestions(final String message) {
        final int mean = message.indexOf("Did you mean");
        if (mean < 0) {
            return List.of();
        }
        final int newline = message.indexOf('\n', mean);
        final String tail = message.substring(newline + 1);
        final List<String> names = new ArrayList<>();
        for (final String part : tail.split(" {8}")) {
            if (!part.isEmpty()) {
                names.add(part);
            }
        }
        return new ArrayList<>(new TreeSet<>(names));
    }

    static void query(final Main main, final Set<Class<?>> classes, final String name) {
        final String message;
        try {
            message = main.getUnknownCommandMessage(classes, name);
        } catch (final Exception e) {
            emit("error", name, e.getClass().getName() + ":" + String.valueOf(e.getMessage()));
            return;
        }
        final List<String> names = suggestions(message);
        if (names.isEmpty()) {
            emit("message", name, message);
        } else {
            emit("suggests", name, String.join(",", names));
        }
    }

    public static void main(final String[] args) {
        final Reachable main = new Reachable();
        final Set<Class<?>> classes = catalogue(main);

        final List<String> names = new ArrayList<>();
        for (final Class<?> clazz : classes) {
            names.add(clazz.getSimpleName());
        }
        final List<String> sorted = new ArrayList<>(new TreeSet<>(names));
        emit("count", "catalogue", Integer.toString(sorted.size()));
        for (int i = 0; i < sorted.size(); i += 100) {
            emit("catalogue", Integer.toString(i / 100),
                    String.join(",", sorted.subList(i, Math.min(i + 100, sorted.size()))));
        }

        // A name one character off a real tool, over the whole catalogue.
        query(main, classes, "PrintRead");
        query(main, classes, "PrintReadz");
        query(main, classes, "HaplotypeCallr");
        query(main, classes, "MarkDuplicate");
        // A prefix several tools share.
        query(main, classes, "CollectQuality");
        query(main, classes, "PathSeq");
        // A substring of many, five characters or more.
        query(main, classes, "Fingerprint");
        // A name far from everything.
        query(main, classes, "Zzzzzzzzzzzzzzzzzzzz");
        // A deprecated name. It is not in the catalogue at all, so it reaches the search, and the
        // registry answering first is the only reason the notice is seen.
        query(main, classes, "IndelRealigner");
        // A Picard tool, to show the second package is scanned.
        query(main, classes, "SortVcf");

        // The names the two guards drop, and a handful the catalogue must hold.
        final List<String> expected = Arrays.asList(
                "PrintReads", "HaplotypeCaller", "MarkDuplicates", "SortSam", "SortVcfs",
                "CNNScoreVariants", "PrintReadsSpark", "PicardCommandLineProgramExecutor",
                "CommandLineArgumentValidator");
        final List<String> present = new ArrayList<>();
        for (final String name : expected) {
            present.add(name + "=" + sorted.contains(name));
        }
        emit("holds", "names", String.join(",", present));

        System.out.print(buf);
    }
}

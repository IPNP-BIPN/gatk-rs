/*
 * Main's tool resolution and its refusals, taken from the reference.
 *
 * `gatk <Tool> <args>` resolves the first token to a tool and hands the rest to Barclay. What is
 * measured is what happens when the first token names nothing: the deprecation notice, the
 * suggestion search and the message they are folded into. The catalogue is not the reference's
 * own here, which would be its whole class path; it is five named classes, so what is pinned is
 * the SEARCH and not the catalogue.
 *
 * Ten behaviours this is built to catch.
 *
 *   - A DEPRECATED TOOL SHORT-CIRCUITS EVERYTHING, its notice naming the version it went out in
 *     and never a suggestion, and a tool that is not deprecated answers null rather than a
 *     message;
 *   - A NAME THAT IS A PREFIX OF A TOOL SCORES ZERO, whatever its length;
 *   - A SUBSTRING SCORES ZERO ONLY FROM FIVE CHARACTERS: `ountR` scores zero against CountReads
 *     while `ount` falls back on the distance, which still finds it;
 *   - THE DISTANCE IS NOT THE PLAIN LEVENSHTEIN ONE. Dropping a character from the command costs
 *     FOUR and adding one costs ONE, so `PrintReadsxy` is over the floor at eight while
 *     `PrtReads` is well under it at two;
 *   - THE FLOOR IS SEVEN, which is what those two straddle;
 *   - `Did you mean this?` BECOMES `Did you mean one of these?` AT TWO MATCHES, counted over the
 *     BEST distance alone, so a prefix of one tool beats a near neighbour and is suggested by
 *     itself;
 *   - THE SUGGESTIONS ARE PRINTED WITH NO SEPARATOR BETWEEN THEM, eight spaces before each and no
 *     line break after, so two of them run together on one line;
 *   - WHEN EVERY TOOL SCORES ZERO THE SUGGESTION IS SUPPRESSED, the distance being bumped over
 *     the floor rather than every tool being listed. A one-tool catalogue the command is a prefix
 *     of is exactly that case;
 *   - THE MESSAGE ALWAYS OPENS ON `'<name>' is not a valid command.` and a line break, suggestion
 *     or not, so a message with nothing to suggest still ends on one;
 *   - AND A NAME THAT DOES MATCH A TOOL IS A RuntimeException rather than an answer, the search
 *     being reached only once resolution has failed.
 *
 * Output:
 *
 *     tools\t<set>\t<the class names, comma-separated>
 *     message\t<case>\t<the message, escaped>
 *     deprecated\t<name>\t<the notice, escaped>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: MainDispatchDump
 */

import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.cmdline.DeprecatedToolsRegistry;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

public class MainDispatchDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** The classes the search scores against, named by their simple names alone. */
    static Set<Class<?>> tools(final List<Class<?>> classes) {
        return new LinkedHashSet<>(classes);
    }

    // A fixed set of real tool classes, whose simple names are what the search reads.
    static final List<Class<?>> CATALOGUE = List.of(
            org.broadinstitute.hellbender.tools.PrintReads.class,
            org.broadinstitute.hellbender.tools.CountReads.class,
            org.broadinstitute.hellbender.tools.CountBases.class,
            org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants.class,
            org.broadinstitute.hellbender.tools.walkers.filters.VariantFiltration.class);

    static void message(final String name, final String command, final List<Class<?>> classes) {
        try {
            emit("message", name, new Main().getUnknownCommandMessage(tools(classes), command));
        } catch (final Exception e) {
            emit("error", name, e.getClass().getName() + ":" + String.valueOf(e.getMessage()));
        }
    }

    public static void main(final String[] args) {
        final List<String> names = new ArrayList<>();
        for (final Class<?> clazz : CATALOGUE) {
            names.add(clazz.getSimpleName());
        }
        emit("tools", "catalogue", String.join(",", names));

        // Every deprecated tool the registry names.
        for (final String tool : List.of("IndelRealigner", "RealignerTargetCreator", "CNNScoreVariants",
                "CNNVariantTrain", "CNNVariantWriteTensors")) {
            emit("deprecated", tool, String.valueOf(DeprecatedToolsRegistry.getToolDeprecationInfo(tool)));
        }
        // And one that is not.
        emit("deprecated", "PrintReads", String.valueOf(DeprecatedToolsRegistry.getToolDeprecationInfo("PrintReads")));

        // A deprecated name short-circuits the search even though the catalogue holds neighbours.
        message("deprecated-short-circuits", "IndelRealigner", CATALOGUE);

        // A prefix of one tool, of one character and of four.
        message("prefix-one-character", "P", CATALOGUE);
        message("prefix-four-characters", "Prin", CATALOGUE);

        // Four characters that are a substring but not a prefix, against five that are.
        message("substring-four-characters", "ount", CATALOGUE);
        message("substring-five-characters", "ountR", CATALOGUE);

        // One character too many and one too few, which are both under the floor.
        message("one-insertion", "PrintReadss", CATALOGUE);
        message("one-deletion", "PrntReads", CATALOGUE);
        message("one-substitution", "PrintReadz", CATALOGUE);

        // Two characters too many against two too few, which the weights price very differently:
        // dropping a character from the command costs four and adding one costs one.
        message("two-too-many", "PrintReadsxy", CATALOGUE);
        message("two-too-few", "PrtReads", CATALOGUE);

        // A name that resolves: the search is never reached, and asking it anyway throws.
        message("name-matches-a-tool", "PrintReads", CATALOGUE);

        // A name far from everything, so nothing is suggested.
        message("nothing-close", "Zzzzzzzzzzzzzzzzzzzz", CATALOGUE);

        // A name that scores zero against every tool at once, which suppresses the suggestion.
        message("every-tool-scores-zero", "", CATALOGUE);

        // A prefix of one tool that is close to another: only the prefix is suggested.
        message("prefix-beats-a-neighbour", "CountRead", CATALOGUE);
        message("substring-of-two", "ount", List.of(
                org.broadinstitute.hellbender.tools.CountReads.class,
                org.broadinstitute.hellbender.tools.CountBases.class));

        // A one-tool catalogue the command is a prefix of, which is still EVERY tool scoring
        // zero, so the suggestion is suppressed there too.
        message("single-tool-zero", "Print", List.of(
                org.broadinstitute.hellbender.tools.PrintReads.class));

        System.out.print(buf);
    }
}

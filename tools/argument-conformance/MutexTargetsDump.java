/*
 * The two names a mutex target has, taken from the reference.
 *
 * `tool-argument-declarations` measured `getMutexTargetList()`, which is the LONG NAME the target
 * resolves to: `--ambig-filter-bases` reports `ambig-filter-frac`. The usage text prints a
 * different name for the same relation, the FIELD name the annotation was written with:
 * `Cannot be used in conjunction with argument(s) maxAmbiguousBaseFraction`. Both are in the
 * reference and only one of them was measured, which is why a walker's composed usage differs from
 * the `usage-text` golden in exactly four lines.
 *
 * Three things this records.
 *
 *   - THE TWO NAMES SIDE BY SIDE, per argument that has a mutex at all, so the map between them is
 *     a measurement rather than a guess;
 *   - THAT THE TWO DIFFER ONLY SOMETIMES: a field whose name IS its long name reports the same
 *     string twice, which is why the sentence looked right everywhere it had been compared;
 *   - AND THE SENTENCE ITSELF, taken off the rendered usage, so the port has the text it has to
 *     produce and not only the names in it.
 *
 * The arguments are reached through the tools' own parsers rather than through a hand-built list:
 * a plugin's arguments exist only once a descriptor has discovered them, and four of the five
 * mutex arguments here are a read filter's.
 *
 * Output:
 *
 *     order\t<tool>\t<long name>\t<the target list, space-joined, in the order it holds>
 *     mutex\t<tool>\t<long name>\t<target long name>\t<target field name>\t<annotated|resolved>
 *     sentence\t<tool>\t<long name>\t<the rendered sentence, escaped>
 *
 * Usage: MutexTargetsDump
 */

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.NamedArgumentDefinition;

import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class MutexTargetsDump {

    public static void main(final String[] args) {
        System.out.println("# MutexTargetsDump: the long name and the field name of a mutex target");
        tool("CountReads", new org.broadinstitute.hellbender.tools.CountReads());
        tool("IndexFeatureFile", new org.broadinstitute.hellbender.tools.IndexFeatureFile());
        tool("SelectVariants",
                new org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants());
        tool("GatherVcfsCloud", new org.broadinstitute.hellbender.tools.GatherVcfsCloud());
    }

    static void tool(final String name,
                     final org.broadinstitute.hellbender.cmdline.CommandLineProgram target) {
        final CommandLineArgumentParser parser =
                (CommandLineArgumentParser) target.getCommandLineParser();
        final String usage = parser.usage(true, true);
        final List<String> rows = new ArrayList<>();
        final List<String> sentences = new ArrayList<>();
        for (final NamedArgumentDefinition definition : parser.getNamedArgumentDefinitions()) {
            if (definition.getMutexTargetList().isEmpty()) {
                continue;
            }
            final List<String> ownOrder = new ArrayList<>(definition.getMutexTargetList());
            // The list in the order it HOLDS, which is the order the refusal joins it in. It is
            // its own row because the rows below are sorted for stability, and sorting the list
            // itself is what hid this: `quantize-quals` reads `static-quantized-quals
            // round-down-quantized`, the order those two declare their own mutex in, and nothing
            // recovers that from an alphabetical list.
            rows.add(String.format("order\t%s\t%s\t%s", name, definition.getLongName(),
                    String.join(" ", ownOrder)));
            final List<String> targets = new ArrayList<>(ownOrder);
            Collections.sort(targets);
            for (final String targetName : targets) {
                // The name the usage prints is the TARGET DEFINITION'S FIELD name, which is
                // neither the target's long name nor anything the annotation holds: the
                // annotation's own `mutex()` list is long names too, as the rows below show.
                String fieldName = "-";
                for (final NamedArgumentDefinition other : parser.getNamedArgumentDefinitions()) {
                    if (other.getLongName().equals(targetName)) {
                        fieldName = other.getUnderlyingField().getName();
                        break;
                    }
                }
                final Field field = definition.getUnderlyingField();
                final Argument annotation = field.getAnnotation(Argument.class);
                final List<String> annotated = new ArrayList<>(List.of(annotation.mutex()));
                Collections.sort(annotated);
                rows.add(String.format("mutex\t%s\t%s\t%s\t%s\t%s", name,
                        definition.getLongName(), targetName, fieldName,
                        annotated.contains(targetName) ? "annotated" : "resolved"));
            }
            // The sentence as the usage renders it, over a text whose wrapping is collapsed: the
            // line breaks fall inside it, and what the port has to produce is the sentence.
            final String collapsed = usage.replaceAll("\\s+", " ");
            final String marker = "Cannot be used in conjunction with argument(s)";
            int at = collapsed.indexOf(marker);
            while (at >= 0) {
                // The sentence runs to the next argument or the next heading, whichever comes
                // first: it is the last piece of a description, so nothing follows it inside one.
                int end = collapsed.length();
                for (final String stop : new String[]{"--", "Valid only if", "Conditional Arguments"}) {
                    final int found = collapsed.indexOf(stop, at);
                    if (found >= 0 && found < end) {
                        end = found;
                    }
                }
                final String sentence = collapsed.substring(at, end).trim();
                if (!sentences.contains(sentence)) {
                    sentences.add(sentence);
                }
                at = collapsed.indexOf(marker, at + 1);
            }
        }
        Collections.sort(rows);
        rows.forEach(System.out::println);
        Collections.sort(sentences);
        for (final String sentence : sentences) {
            System.out.printf("sentence\t%s\t%s%n", name,
                    sentence.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"));
        }
    }
}

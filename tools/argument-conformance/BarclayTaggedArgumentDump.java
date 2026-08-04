/*
 * Barclay's tagged arguments and its collection-file expansion, taken from the reference.
 *
 * Two mechanisms that sit either side of the value model, and that a covering-array vector reaches
 * as soon as it names a file or a logical name:
 *
 *   - TAGS. `--argument:logical_name,key=value raw_value` is rewritten BEFORE jopt-simple sees the
 *     command line: `TaggedArgumentParser.preprocessTaggedOptions` peels the tag off, stores
 *     Pair(tag, value) in a map under a surrogate key built as `option_string:value`, and hands
 *     the parser `--argument` and that key. So the *value* jopt-simple parses is a synthetic
 *     string, and the real one is retrieved later. Three consequences are visible: the surrogate
 *     key is the uniqueness test, so the same option with the same tag and the same value twice is
 *     "duplicated on the command line" and not two values; a tag on a field whose type does not
 *     implement TaggedArgument is an error raised at value-setting time, naming the argument as
 *     "shortName/fullName" even when the short name is empty; and the tag is populated onto the
 *     value BEFORE the value is stored;
 *   - EXPANSION. A COLLECTION argument whose value ends in `.list` or `.args` is replaced by the
 *     lines of that file, trimmed, with blanks and `#` comments dropped. A scalar argument with
 *     the same value is not expanded at all, and neither is a collection value ending in anything
 *     else. A file whose surviving lines start with `@` produces a warning and is expanded anyway.
 *
 * Output:
 *
 *     case\t<label>\t<argv, space separated>
 *     result\t<label>\tok|E:<exception class>:<message>
 *     field\t<label>\t<field name>\t<value>
 *     tag\t<label>\t<field name>\t<index>\t<tag>\t<attributes>
 *
 * Usage: BarclayTaggedArgumentDump
 */

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.TaggedArgument;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

public class BarclayTaggedArgumentDump {

    /**
     * A value type that accepts tags. Its String constructor is what `constructFromString` calls,
     * and `populateArgumentTags` writes onto the instance afterwards.
     */
    public static final class TaggedPath implements TaggedArgument {
        private final String value;
        private String tag;
        private Map<String, String> attributes = Collections.emptyMap();

        public TaggedPath(final String value) { this.value = value; }

        @Override public void setTag(final String tagName) { this.tag = tagName; }
        @Override public String getTag() { return tag; }
        @Override public void setTagAttributes(final Map<String, String> attributes) { this.attributes = attributes; }
        @Override public Map<String, String> getTagAttributes() { return attributes; }
        @Override public String toString() { return value; }
    }

    public static final class Args {
        @Argument(fullName = "tagged-collection", optional = true, doc = "a collection that accepts tags")
        public List<TaggedPath> taggedCollection = new ArrayList<>();

        @Argument(fullName = "tagged-scalar", shortName = "T", optional = true, doc = "a scalar that accepts tags")
        public TaggedPath taggedScalar;

        @Argument(fullName = "plain-collection", optional = true, doc = "a collection of plain strings")
        public List<String> plainCollection = new ArrayList<>();

        @Argument(fullName = "plain-scalar", optional = true, doc = "a scalar that does not accept tags")
        public String plainScalar;

        @Argument(fullName = "no-expansion", optional = true, suppressFileExpansion = true,
                doc = "a collection that refuses expansion")
        public List<String> noExpansion = new ArrayList<>();
    }

    public static void main(final String[] args) throws Exception {
        // A fixed relative directory, not a temporary one: the paths appear in the golden, and a
        // temporary directory's name changes every run, so the golden would never match itself.
        final Path dir = Path.of("fixtures");
        if (Files.isDirectory(dir)) {
            try (final java.util.stream.Stream<Path> entries = Files.list(dir)) {
                for (final Path entry : entries.collect(java.util.stream.Collectors.toList())) {
                    Files.delete(entry);
                }
            }
        }
        Files.createDirectories(dir);

        // Three files with the same body and three extensions: two expand, one does not.
        final String body = String.join("\n",
                "first",
                "",
                "# a comment",
                "  second  ",
                "third",
                "");
        final Path list = dir.resolve("values.list");
        final Path argsFile = dir.resolve("values.args");
        final Path text = dir.resolve("values.txt");
        Files.write(list, body.getBytes());
        Files.write(argsFile, body.getBytes());
        Files.write(text, body.getBytes());
        // A file whose lines look like a sequence dictionary: expanded, with a warning.
        final Path dictionary = dir.resolve("looks-like-a-dict.list");
        Files.write(dictionary, "@HD\tVN:1.6\nchr1\n".getBytes());
        final Path empty = dir.resolve("empty.list");
        Files.write(empty, "\n# only a comment\n".getBytes());

        System.out.println("# BarclayTaggedArgumentDump: tags and collection-file expansion");

        // Tags on a field whose type accepts them.
        run("tag-name-only", "--tagged-scalar:tumour", "a.bam");
        run("tag-with-attribute", "--tagged-scalar:tumour,kind=wgs", "a.bam");
        run("tag-with-two-attributes", "--tagged-scalar:tumour,kind=wgs,lane=3", "a.bam");
        run("tag-on-short-name", "-T:tumour", "a.bam");
        run("untagged", "--tagged-scalar", "a.bam");
        run("tag-on-collection", "--tagged-collection:one", "a.bam", "--tagged-collection:two", "b.bam");
        // The surrogate key is option string plus value, so this pair collides with itself.
        run("same-option-tag-and-value-twice",
                "--tagged-collection:one", "a.bam", "--tagged-collection:one", "a.bam");
        // Same tag, different value: two distinct keys, so both survive.
        run("same-tag-different-values",
                "--tagged-collection:one", "a.bam", "--tagged-collection:one", "b.bam");

        // Malformed tags.
        run("zero-length-tag", "--tagged-scalar:", "a.bam");
        run("zero-length-argument-name", "--:tumour", "a.bam");
        run("tag-with-no-value", "--tagged-scalar:tumour");
        run("tag-followed-by-option", "--tagged-scalar:tumour", "--plain-scalar", "x");
        run("empty-attribute", "--tagged-scalar:tumour,", "a.bam");
        run("attribute-without-equals", "--tagged-scalar:tumour,kind", "a.bam");
        run("attribute-with-empty-value", "--tagged-scalar:tumour,kind=", "a.bam");
        run("duplicate-attribute-key", "--tagged-scalar:tumour,kind=wgs,kind=wes", "a.bam");
        run("tag-name-containing-equals", "--tagged-scalar:kind=wgs", "a.bam");

        // A tag on a field whose type does not implement TaggedArgument.
        run("tag-on-plain-scalar", "--plain-scalar:tumour", "a.bam");
        run("tag-on-plain-collection", "--plain-collection:tumour", "a.bam");

        // Expansion, which is a collection-only mechanism.
        run("expand-list", "--plain-collection", list.toString());
        run("expand-args", "--plain-collection", argsFile.toString());
        run("no-expand-other-extension", "--plain-collection", text.toString());
        run("no-expand-scalar", "--plain-scalar", list.toString());
        run("expansion-suppressed", "--no-expansion", list.toString());
        run("expand-with-at-sign", "--plain-collection", dictionary.toString());
        run("expand-empty-file", "--plain-collection", empty.toString());
        run("expand-missing-file", "--plain-collection", dir.resolve("nope.list").toString());
        // Expansion and ordinary values in one argument, and the order they end up in.
        run("expand-among-values",
                "--plain-collection", "before",
                "--plain-collection", list.toString(),
                "--plain-collection", "after");
        // A tagged collection whose value is an expansion file: the tag is propagated to every
        // value the file produced.
        run("expand-tagged", "--tagged-collection:one", list.toString());
    }

    static void run(final String label, final String... argv) {
        final Args target = new Args();
        System.out.printf("case\t%s\t%s%n", label, String.join(" ", argv));

        String result;
        try {
            final PrintStream sink = new PrintStream(new ByteArrayOutputStream());
            final boolean parsed = new CommandLineArgumentParser(target).parseArguments(sink, argv);
            result = parsed ? "ok" : "not-parsed";
        } catch (final Exception | AssertionError e) {
            // One message embeds a newline ("Duplicate key %s\n" + USAGE). Escaped rather than
            // printed, so a row stays a row and the golden stays line-comparable.
            result = "E:" + e.getClass().getName() + ":"
                    + String.valueOf(e.getMessage()).replace("\n", "\\n");
        }
        System.out.printf("result\t%s\t%s%n", label, result);

        if (result.equals("ok")) {
            System.out.printf("field\t%s\t%s\t%s%n", label, "tagged-collection", target.taggedCollection);
            System.out.printf("field\t%s\t%s\t%s%n", label, "tagged-scalar", target.taggedScalar);
            System.out.printf("field\t%s\t%s\t%s%n", label, "plain-collection", target.plainCollection);
            System.out.printf("field\t%s\t%s\t%s%n", label, "plain-scalar", target.plainScalar);
            System.out.printf("field\t%s\t%s\t%s%n", label, "no-expansion", target.noExpansion);
            for (int i = 0; i < target.taggedCollection.size(); i++) {
                tag(label, "tagged-collection", i, target.taggedCollection.get(i));
            }
            if (target.taggedScalar != null) {
                tag(label, "tagged-scalar", 0, target.taggedScalar);
            }
        }
    }

    static void tag(final String label, final String name, final int index, final TaggedPath value) {
        final StringJoiner attributes = new StringJoiner(",");
        // The attribute map is a LinkedHashMap, so the order is the command line's.
        new LinkedHashMap<>(value.getTagAttributes()).forEach((k, v) -> attributes.add(k + "=" + v));
        System.out.printf("tag\t%s\t%s\t%d\t%s\t%s%n",
                label, name, index, String.valueOf(value.getTag()), attributes);
    }
}

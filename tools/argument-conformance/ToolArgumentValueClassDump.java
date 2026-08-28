/*
 * The four value classes a tool declares and no golden converts, taken from the reference.
 *
 * `ValueClass` carries the conversions that are measured: an integer, a double, a string, a boolean
 * and an enum. The declarations name four more, and every one of the seven ported tools declares at
 * least one of them, which is why no command line reaches the ported parser yet.
 *
 * Ten behaviours this is built to catch.
 *
 *   - `constructFromString` USES THE STRING CONSTRUCTOR, so what a class accepts is what its
 *     constructor accepts, and a class that accepts anything refuses nothing;
 *   - A `File` ACCEPTS ANY STRING, including one that names nothing and one that is empty, so a
 *     bad path is not a bad value;
 *   - A `GATKPath` ACCEPTS A URI AS WELL AS A PATH, and renders as the string it was given rather
 *     than as a resolved path;
 *   - A `GATKPath` IS TAGGABLE AND A `File` IS NOT, which is the difference the same command line
 *     turns into a refusal on one and a tag on the other;
 *   - THE REFUSAL FOR A TAG ON AN UNTAGGABLE FIELD NAMES THE ARGUMENT as `shortName/fullName`;
 *   - A `FeatureInput` IS TAGGABLE AND ITS TAG IS ITS NAME, which is what a walker then looks the
 *     feature up by, and an untagged one takes a name from the file;
 *   - A `Float` REFUSES WHAT `Float.valueOf` REFUSES, and its message names `Float` where a
 *     double's names `Double`;
 *   - A `Float` ACCEPTS WHAT `Double` ACCEPTS in exponent form, and rounds to the nearer float;
 *   - AND `Float.valueOf` IS NOT A NAIVE PARSE: it takes a trailing `f` or `d`, a hexadecimal
 *     literal, leading and trailing whitespace and a leading plus, and it spells its infinity and
 *     its not-a-number with capitals, refusing `inf` and `nan`. An underscore is refused too.
 *     Every one of those is a spelling a port would otherwise guess at, and Rust's own `parse`
 *     disagrees with Java on six of them;
 *   - THE DEFAULT RENDERING IS `String.valueOf(field)`, so an initialised File renders as its path
 *     and an uninitialised one as `null`;
 *   - AND A COLLECTION OF ANY OF THEM RENDERS AS `null` WHEN EMPTY, which is what makes it
 *     optional.
 *
 * Output:
 *
 *     case\t<label>\t<argv, space separated>
 *     result\t<label>\tok|E:<exception class>:<message>
 *     field\t<label>\t<field name>\t<value>
 *     tag\t<label>\t<field name>\t<tag>\t<attributes>
 *     default\t<field name>\t<the rendering a constructed instance has>
 *
 * Usage: ToolArgumentValueClassDump
 */

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.NamedArgumentDefinition;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKPath;

import htsjdk.variant.variantcontext.VariantContext;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

public class ToolArgumentValueClassDump {

    /** One field per class, plus a collection of each, so both renderings are visible. */
    public static final class Target {
        @Argument(fullName = "path", shortName = "P", doc = "a GATKPath", optional = true)
        public GATKPath path;

        @Argument(fullName = "paths", doc = "a collection of them", optional = true)
        public List<GATKPath> paths = new ArrayList<>();

        @Argument(fullName = "file", shortName = "F", doc = "a File", optional = true)
        public File file;

        @Argument(fullName = "initialised-file", doc = "a File with a value", optional = true)
        public File initialisedFile = new File("already/here");

        @Argument(fullName = "fraction", shortName = "FR", doc = "a Float", optional = true)
        public Float fraction;

        @Argument(fullName = "primitive-fraction", doc = "a float", optional = true)
        public float primitiveFraction = 0.5f;

        @Argument(fullName = "feature", doc = "a FeatureInput", optional = true)
        public FeatureInput<VariantContext> feature;
    }

    static void emit(final String kind, final String... parts) {
        final StringJoiner joined = new StringJoiner("\t");
        joined.add(kind);
        for (final String part : parts) {
            joined.add(escape(part));
        }
        System.out.println(joined);
    }

    static String escape(final String text) {
        return text == null ? "null"
                : text.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
    }

    /** One command line against a fresh target, and what each field holds afterwards. */
    static void run(final String label, final String... argv) {
        final Target target = new Target();
        final CommandLineArgumentParser parser = new CommandLineArgumentParser(target);
        emit("case", label, String.join(" ", argv));
        String result;
        try {
            final PrintStream sink = new PrintStream(new ByteArrayOutputStream());
            result = parser.parseArguments(sink, argv) ? "ok" : "not-parsed";
        } catch (final Exception | AssertionError e) {
            result = "E:" + e.getClass().getName() + ":" + e.getMessage();
        }
        emit("result", label, result);
        if (!result.equals("ok")) {
            return;
        }
        emit("field", label, "path", String.valueOf(target.path));
        emit("field", label, "paths", String.valueOf(target.paths));
        emit("field", label, "file", String.valueOf(target.file));
        emit("field", label, "initialised-file", String.valueOf(target.initialisedFile));
        emit("field", label, "fraction", String.valueOf(target.fraction));
        emit("field", label, "primitive-fraction", String.valueOf(target.primitiveFraction));
        emit("field", label, "feature", String.valueOf(target.feature));
        if (target.path != null) {
            emit("tag", label, "path", String.valueOf(target.path.getTag()),
                    attributes(target.path.getTagAttributes()));
        }
        if (target.feature != null) {
            // A FeatureInput's tag is its NAME, which is what a walker looks the feature up by.
            emit("tag", label, "feature", String.valueOf(target.feature.getTag()),
                    attributes(target.feature.getTagAttributes()));
            emit("field", label, "feature-name", String.valueOf(target.feature.getName()));
            emit("field", label, "feature-path",
                    String.valueOf(target.feature.getFeaturePath()));
        }
    }

    static String attributes(final Map<String, String> attributes) {
        if (attributes == null || attributes.isEmpty()) {
            return "";
        }
        final StringJoiner joined = new StringJoiner(",");
        // The map's own order is what the reference exposes; sorting it here would hide a change.
        attributes.forEach((key, value) -> joined.add(key + "=" + value));
        return joined.toString();
    }

    public static void main(final String[] args) {
        // What a constructed instance renders as, which is what decides optionality.
        final CommandLineArgumentParser parser = new CommandLineArgumentParser(new Target());
        for (final NamedArgumentDefinition definition : parser.getNamedArgumentDefinitions()) {
            emit("default", definition.getLongName(),
                    String.valueOf(definition.getDefaultValueAsString()));
        }

        // A path, a URI, a relative path and one that names nothing.
        run("a-path", "--path", "/tmp/reads.bam");
        run("a-uri", "--path", "gs://bucket/reads.bam");
        run("a-relative-path", "--path", "reads.bam");
        run("a-path-that-is-not-there", "--path", "/no/such/file.bam");
        run("an-empty-path", "--path", "");

        // A file, which accepts the same strings and refuses nothing.
        run("a-file", "--file", "/tmp/reads.bam");
        run("a-file-that-is-not-there", "--file", "/no/such/file.bam");
        run("an-empty-file", "--file", "");

        // The tag, on the taggable class and on the one that is not.
        run("a-tagged-path", "--path:name,key=value", "/tmp/reads.bam");
        run("a-tagged-path-without-attributes", "--path:name", "/tmp/reads.bam");
        run("a-tagged-file", "--file:name", "/tmp/reads.bam");
        run("a-tagged-path-by-short-name", "-P:name", "/tmp/reads.bam");

        // A collection of paths, which is where the empty rendering matters.
        run("no-paths");
        run("one-path", "--paths", "/tmp/a.bam");
        run("two-paths", "--paths", "/tmp/a.bam", "--paths", "/tmp/b.bam");

        // A feature input, tagged and not.
        run("a-feature", "--feature", "/tmp/sites.vcf");
        run("a-tagged-feature", "--feature:known", "/tmp/sites.vcf");

        // The float, on both sides of what `Float.valueOf` accepts.
        run("a-fraction", "--fraction", "0.25");
        run("a-fraction-in-exponent-form", "--fraction", "1e3");
        run("a-fraction-that-is-not-a-number", "--fraction", "abc");
        run("a-fraction-that-is-a-word", "--fraction", "NaN");
        run("a-fraction-out-of-a-floats-range", "--fraction", "1e40");
        // `Float.valueOf` is not `str::parse`: Java's grammar takes a trailing type suffix and a
        // hex literal, and spells its infinity and its not-a-number with capitals. Every one of
        // these is a spelling a port would otherwise have to guess at.
        run("a-fraction-with-a-float-suffix", "--fraction", "1.5f");
        run("a-fraction-with-a-double-suffix", "--fraction", "1.5d");
        run("a-fraction-in-hexadecimal", "--fraction", "0x1p3");
        run("a-fraction-spelled-inf", "--fraction", "inf");
        run("a-fraction-spelled-infinity", "--fraction", "Infinity");
        run("a-fraction-spelled-negative-infinity", "--fraction", "-Infinity");
        run("a-fraction-spelled-nan-in-lower-case", "--fraction", "nan");
        run("a-fraction-with-a-leading-space", "--fraction", " 1.5");
        run("a-fraction-with-a-trailing-space", "--fraction", "1.5 ");
        run("a-fraction-with-an-underscore", "--fraction", "1_000");
        run("a-fraction-with-a-leading-plus", "--fraction", "+1.5");
        run("a-fraction-that-rounds-to-a-float", "--fraction", "0.1000000000000000055511151231257827");

        run("a-primitive-fraction", "--primitive-fraction", "0.25");
        run("a-primitive-fraction-that-is-not-a-number", "--primitive-fraction", "abc");
    }
}

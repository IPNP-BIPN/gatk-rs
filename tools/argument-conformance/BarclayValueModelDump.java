/*
 * Barclay's argument value model, taken from the reference.
 *
 * This is the layer a covering-array vector is interpreted by: `CommandLineArgumentParser` over a
 * set of `@Argument` fields, run at library level, with no GATK tool around it. The unified CLI
 * dispatcher is out of scope; what is in scope is which vectors are accepted, what value each
 * field ends up holding, and which exception the rejected ones produce.
 *
 * Six of the rules are not what the annotation names suggest:
 *
 *   - `optional()` is not what decides whether an argument is optional. `isOptional` is
 *     `annotation.optional() || !defaultValueAsString.equals("null")`, so a field declared
 *     `optional = false` that was initialised to anything is optional anyway. An empty collection
 *     counts as "null" for this test, and a non-empty one does not;
 *   - `"null"` is a value rather than a token, with three different outcomes: on a collection it
 *     clears the collection (and warns, if it is not the first value), on a non-optional argument
 *     it throws, and on a scalar whose FIELD (not boxed class) is primitive it throws a different
 *     exception;
 *   - the bounds check treats a null value as out of range: `isValueOutOfRange` starts with
 *     `value == null ||`, so `--bounded-int null` is an OutOfRangeArgumentValue rather than an
 *     accepted null;
 *   - the RECOMMENDED range is checked with `isValueOutOfRange`, which compares against the HARD
 *     minValue/maxValue. For an argument with a recommended range and no hard range those are
 *     infinities, so the warning can fire only for a null value; for an argument with both, the
 *     hard check has already thrown. The recommended-range warning is nearly unreachable, and the
 *     `recommended-out` case is here to record that nothing happens;
 *   - a scalar refuses a second occurrence with `getHasBeenSet() || originalValues.size() > 1`,
 *     which is a BadArgumentValue and not a "last one wins";
 *   - a collection is cleared once before the first value unless the parser is in
 *     APPEND_TO_COLLECTIONS mode, so its declared default is discarded the moment the user names
 *     the argument at all.
 *
 * And `validateValues` treats a mutex partner having been set as satisfying "required", which is
 * why the required test is `!getHasBeenSet() && providedMutexArguments.isEmpty()`.
 *
 * Output:
 *
 *     case\t<label>\t<argv, space separated>
 *     result\t<label>\tok|E:<exception class>:<message>
 *     field\t<label>\t<field name>\t<value>
 *
 * Usage: BarclayValueModelDump
 */

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class BarclayValueModelDump {

    /** An enum argument, to reach `Enum.valueOf` and the message its failure produces. */
    public enum Mode { FAST, SLOW }

    /**
     * One field per rule under test.
     *
     * `defaultedRequired` is the field that shows `optional()` is not the whole story: it is
     * declared without `optional = true` and is optional anyway, because it was initialised.
     */
    public static final class Args {
        @Argument(fullName = "required-string", shortName = "R", doc = "required, uninitialised")
        public String requiredString;

        @Argument(fullName = "required-collection", doc = "required: an empty collection reads as null")
        public List<String> requiredCollection = new ArrayList<>();

        @Argument(fullName = "optional-string", shortName = "S", optional = true, doc = "optional scalar")
        public String optionalString;

        @Argument(fullName = "defaulted-required", doc = "declared required, initialised, therefore optional")
        public String defaultedRequired = "preset";

        @Argument(fullName = "flag", optional = true, doc = "a boolean, which may appear with no value")
        public boolean flag = false;

        @Argument(fullName = "bounded-int", optional = true, minValue = 1, maxValue = 10,
                doc = "a hard range")
        public Integer boundedInt;

        @Argument(fullName = "primitive-int", optional = true, doc = "a primitive field, so not nullable")
        public int primitiveInt = 7;

        @Argument(fullName = "recommended-int", optional = true, minRecommendedValue = 5,
                maxRecommendedValue = 8, doc = "a recommended range and no hard range")
        public Integer recommendedInt;

        @Argument(fullName = "bounded-double", optional = true, minValue = 0.0,
                doc = "a minimum only, so the message names one bound")
        public Double boundedDouble;

        @Argument(fullName = "collection", optional = true, doc = "an optional collection")
        public List<String> collection = new ArrayList<>(Arrays.asList("declared"));

        @Argument(fullName = "mutex-a", optional = true, mutex = {"mutex-b"}, doc = "mutex with b")
        public String mutexA;

        @Argument(fullName = "mutex-b", optional = true, mutex = {"mutex-a"}, doc = "mutex with a")
        public String mutexB;

        @Argument(fullName = "enum-arg", optional = true, doc = "an enum, built by Enum.valueOf")
        public Mode enumArg;
    }

    /** The two arguments every accepted vector has to carry, unless the case is about them. */
    static final String[] REQUIRED = {"--required-string", "r", "--required-collection", "c"};

    public static void main(final String[] args) {
        System.out.println("# BarclayValueModelDump: Barclay's argument value model");

        // Nothing but the required pair: every other field shows its declared default.
        run("minimal", with());

        // The required arguments, and what "required" means.
        run("nothing-given");
        run("only-required-string", "--required-string", "r");
        run("only-required-collection", "--required-collection", "c");
        // `defaulted-required` is declared without `optional = true` and is never given. It does
        // not fail, because it was initialised.
        run("required-collection-null", "--required-string", "r", "--required-collection", "null");

        // Scalars.
        run("scalar-once", with("--optional-string", "value"));
        run("scalar-twice", with("--optional-string", "a", "--optional-string", "b"));
        run("scalar-equals-syntax", with("--optional-string=value"));
        run("scalar-short-name", with("-S", "value"));
        run("scalar-null-boxed", with("--optional-string", "null"));
        run("scalar-null-primitive", with("--primitive-int", "null"));
        run("scalar-empty-value", with("--optional-string", ""));

        // Collections. The declared default is discarded as soon as the argument is named.
        run("collection-untouched", with());
        run("collection-one", with("--collection", "a"));
        run("collection-two", with("--collection", "a", "--collection", "b"));
        run("collection-null-first", with("--collection", "null"));
        run("collection-null-after-values", with("--collection", "a", "--collection", "null"));
        run("collection-values-after-null", with("--collection", "null", "--collection", "a"));

        // Flags.
        run("flag-bare", with("--flag"));
        run("flag-true", with("--flag", "true"));
        run("flag-false", with("--flag", "false"));
        run("flag-bad-value", with("--flag", "maybe"));

        // The hard range.
        run("bounded-in-range", with("--bounded-int", "5"));
        run("bounded-at-min", with("--bounded-int", "1"));
        run("bounded-at-max", with("--bounded-int", "10"));
        run("bounded-below", with("--bounded-int", "0"));
        run("bounded-above", with("--bounded-int", "11"));
        run("bounded-null", with("--bounded-int", "null"));
        run("bounded-not-a-number", with("--bounded-int", "x"));
        // A minimum with no maximum, on a Double: the message names one bound and formats it as a
        // double rather than as an integer.
        run("double-below-min", with("--bounded-double", "-1"));
        run("double-in-range", with("--bounded-double", "0.5"));

        // The recommended range, which is checked against the hard bounds.
        run("recommended-in-range", with("--recommended-int", "6"));
        run("recommended-far-out", with("--recommended-int", "100"));

        // Enums.
        run("enum-valid", with("--enum-arg", "FAST"));
        run("enum-wrong-case", with("--enum-arg", "fast"));
        run("enum-unknown", with("--enum-arg", "NOPE"));

        // Mutex.
        run("mutex-one", with("--mutex-a", "x"));
        run("mutex-both", with("--mutex-a", "x", "--mutex-b", "y"));

        // Unknown and malformed.
        run("unknown-argument", with("--nope", "1"));
        run("missing-value", with("--optional-string"));
        run("positional-value", with("bare"));
    }

    static String[] with(final String... extra) {
        final List<String> argv = new ArrayList<>(Arrays.asList(REQUIRED));
        argv.addAll(Arrays.asList(extra));
        return argv.toArray(new String[0]);
    }

    static void run(final String label, final String... argv) {
        final Args target = new Args();
        System.out.printf("case\t%s\t%s%n", label, String.join(" ", argv));

        String result;
        try {
            // The parser writes usage and warnings to the stream it is given. It is swallowed
            // here: what it prints is the help layer, which is a separate slice, and letting it
            // reach stdout would put a usage block in the middle of the golden.
            final PrintStream sink = new PrintStream(new ByteArrayOutputStream());
            final boolean parsed = new CommandLineArgumentParser(target).parseArguments(sink, argv);
            result = parsed ? "ok" : "not-parsed";
        } catch (final Exception | AssertionError e) {
            result = "E:" + e.getClass().getName() + ":" + e.getMessage();
        }
        System.out.printf("result\t%s\t%s%n", label, result);

        if (result.equals("ok")) {
            field(label, "required-string", target.requiredString);
            field(label, "required-collection", target.requiredCollection);
            field(label, "optional-string", target.optionalString);
            field(label, "defaulted-required", target.defaultedRequired);
            field(label, "flag", target.flag);
            field(label, "bounded-int", target.boundedInt);
            field(label, "primitive-int", target.primitiveInt);
            field(label, "recommended-int", target.recommendedInt);
            field(label, "bounded-double", target.boundedDouble);
            field(label, "collection", target.collection);
            field(label, "mutex-a", target.mutexA);
            field(label, "mutex-b", target.mutexB);
            field(label, "enum-arg", target.enumArg);
        }
    }

    /** `String.valueOf`, so a null field is the four characters `null` and not an empty column. */
    static void field(final String label, final String name, final Object value) {
        System.out.printf("field\t%s\t%s\t%s%n", label, name, String.valueOf(value));
    }
}

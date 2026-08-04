/*
 * Barclay's @ArgumentCollection flattening, taken from the reference.
 *
 * This is how `-L`, `-XL` and the read-filter arguments actually reach a tool: they are not
 * declared on the tool at all. They live on collection objects the tool holds, and
 * `createArgumentDefinitions` recurses into those objects and adds their `@Argument` fields to
 * ONE flat namespace. Nothing in the resulting command line says which object an argument came
 * from.
 *
 * Three things about that flattening are observable and none is stated anywhere:
 *
 *   - THE ORDER IS SUBCLASS-FIRST. `CommandLineParserUtilities.getAllFields` walks
 *     `clazz.getDeclaredFields()` and then climbs to the superclass, so a subclass's own fields
 *     are registered BEFORE the fields it inherits. That order is the order values are propagated
 *     in and the order `validateArgumentValues` reports a missing required argument in, so which
 *     of two missing arguments a user is told about depends on which class declared it;
 *   - THE RECURSION IS DEPTH-FIRST, AT THE POINT THE FIELD APPEARS. A collection declared between
 *     two `@Argument` fields inserts all of its arguments between them, not after them;
 *   - A DUPLICATE ALIAS IS A CONSTRUCTION FAILURE, not a shadowing rule. Two collections that
 *     happen to declare the same name make the parser unusable before it sees a command line, and
 *     the message names the alias display string rather than the field.
 *
 * Two more refusals happen while the definitions are built rather than while a command line is
 * parsed: an `@ArgumentCollection` field left null, and a field carrying both annotations.
 *
 * Output:
 *
 *     defs\t<label>\t<index>\t<long name>
 *     case\t<label>\t<argv, space separated>
 *     result\t<label>\tok|E:<exception class>:<message>
 *     field\t<label>\t<long name>\t<value>
 *
 * Usage: BarclayArgumentCollectionDump
 */

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.NamedArgumentDefinition;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class BarclayArgumentCollectionDump {

    /** The innermost collection, to reach a recursion two levels deep. */
    public static final class Inner {
        @Argument(fullName = "inner-one", optional = true, doc = "innermost")
        public String innerOne;

        @Argument(fullName = "inner-two", optional = true, doc = "innermost")
        public String innerTwo;
    }

    /** A collection that itself holds a collection. */
    public static final class Middle {
        @Argument(fullName = "middle-before", optional = true, doc = "declared before the nesting")
        public String middleBefore;

        @ArgumentCollection
        public Inner inner = new Inner();

        @Argument(fullName = "middle-after", optional = true, doc = "declared after the nesting")
        public String middleAfter;
    }

    /** A base class, whose fields are registered *after* the subclass's own. */
    public static class Base {
        @Argument(fullName = "base-required", doc = "required, declared on the base class")
        public String baseRequired;

        @Argument(fullName = "base-optional", optional = true, doc = "optional, on the base class")
        public String baseOptional;
    }

    /** The ordinary case: a subclass with its own fields and two nested collections. */
    public static final class Derived extends Base {
        @Argument(fullName = "derived-required", doc = "required, declared on the subclass")
        public String derivedRequired;

        @ArgumentCollection
        public Middle middle = new Middle();

        @Argument(fullName = "derived-last", optional = true, doc = "declared after the nesting")
        public String derivedLast;
    }

    /** Two collections declaring the same argument name. */
    public static final class ClashingA {
        @Argument(fullName = "clash", optional = true, doc = "one of two")
        public String clash;
    }

    public static final class ClashingB {
        @Argument(fullName = "clash", optional = true, doc = "the other of two")
        public String clash;
    }

    public static final class Clashing {
        @ArgumentCollection
        public ClashingA first = new ClashingA();

        @ArgumentCollection
        public ClashingB second = new ClashingB();
    }

    /** An `@ArgumentCollection` field the constructor left null. */
    public static final class Uninitialised {
        @ArgumentCollection
        public Inner inner;
    }

    /** A field carrying both annotations, which is refused. */
    public static final class BothAnnotations {
        @Argument(fullName = "both", optional = true, doc = "annotated twice")
        @ArgumentCollection
        public Inner both = new Inner();
    }

    /** A collection whose own field shadows one the outer object declares. */
    public static final class ShadowingInner {
        @Argument(fullName = "derived-last", optional = true, doc = "the same name the outer uses")
        public String derivedLast;
    }

    public static final class Shadowing {
        @Argument(fullName = "derived-last", optional = true, doc = "declared on the outer object")
        public String derivedLast;

        @ArgumentCollection
        public ShadowingInner inner = new ShadowingInner();
    }

    public static void main(final String[] args) {
        System.out.println("# BarclayArgumentCollectionDump: @ArgumentCollection flattening");

        // The definition order, which is what everything else follows from.
        definitions("derived", new Derived());

        // Every argument of every nesting level, named without any prefix.
        run("all-levels", new Derived(),
                "--base-required", "b", "--derived-required", "d",
                "--base-optional", "bo", "--derived-last", "dl",
                "--middle-before", "mb", "--middle-after", "ma",
                "--inner-one", "i1", "--inner-two", "i2");
        // Nothing but the two required ones, to see which of the two is checked first when both
        // are missing.
        run("nothing-given", new Derived());
        run("only-derived-required", new Derived(), "--derived-required", "d");
        run("only-base-required", new Derived(), "--base-required", "b");
        // An inner argument alone, with the required ones supplied.
        run("inner-only", new Derived(),
                "--base-required", "b", "--derived-required", "d", "--inner-two", "i2");

        // The construction failures, which happen before any command line is seen.
        run("clashing-aliases", new Clashing());
        run("uninitialised-collection", new Uninitialised());
        run("both-annotations", new BothAnnotations());
        run("shadowing", new Shadowing(), "--derived-last", "x");
    }

    /** The flat list of definitions, in the order the parser built them. */
    static void definitions(final String label, final Object target) {
        try {
            final List<NamedArgumentDefinition> defs =
                    new CommandLineArgumentParser(target).getNamedArgumentDefinitions();
            for (int i = 0; i < defs.size(); i++) {
                System.out.printf("defs\t%s\t%d\t%s%n", label, i, defs.get(i).getLongName());
            }
        } catch (final Exception | AssertionError e) {
            System.out.printf("defs\t%s\tE:%s:%s%n", label, e.getClass().getName(), e.getMessage());
        }
    }

    static void run(final String label, final Object target, final String... argv) {
        System.out.printf("case\t%s\t%s%n", label, String.join(" ", argv));

        String result;
        List<NamedArgumentDefinition> defs = null;
        try {
            final PrintStream sink = new PrintStream(new ByteArrayOutputStream());
            final CommandLineArgumentParser parser = new CommandLineArgumentParser(target);
            defs = parser.getNamedArgumentDefinitions();
            result = parser.parseArguments(sink, argv) ? "ok" : "not-parsed";
        } catch (final Exception | AssertionError e) {
            result = "E:" + e.getClass().getName() + ":"
                    + String.valueOf(e.getMessage()).replace("\n", "\\n");
        }
        System.out.printf("result\t%s\t%s%n", label, result);

        if (result.equals("ok")) {
            for (final NamedArgumentDefinition def : defs) {
                System.out.printf("field\t%s\t%s\t%s%n",
                        label, def.getLongName(), String.valueOf(def.getArgumentValue()));
            }
        }
    }

    static List<String> argv(final String... values) {
        return new ArrayList<>(Arrays.asList(values));
    }
}

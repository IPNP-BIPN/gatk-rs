/*
 * NestedIntegerArray and RecalibrationTables, taken from the reference.
 *
 * The container a recalibration table is. Four sparse arrays of RecalDatum indexed by covariate key,
 * two of them "special" because the report writes them differently, and the shapes come from the
 * covariates' own maximumKeyValue. What BaseRecalibrator fills and ApplyBQSR looks up lives here.
 *
 * Eight behaviours this is built to catch, and most of them are about which lookups answer null and
 * which throw.
 *
 *   - THE VARARGS get() DOES NOT BOUNDS-CHECK ITS LAST KEY. It tests `keys[i] >= dimensions[i]` for
 *     the NESTED dimensions only and then indexes the leaf array directly, so an out-of-range last
 *     key is an ArrayIndexOutOfBoundsException where every other out-of-range key is null;
 *   - AND THE SPECIALISED getNKeys() DO. get2Keys, get3Keys and get4Keys test EVERY key including
 *     the last, so get(a, b) and get2Keys(a, b) disagree on exactly the case above. They are
 *     documented as a performance specialisation of the same function and they are not the same
 *     function;
 *   - EXCEPT get1Key, WHICH CHECKS NOTHING AT ALL. Its comment says "bounds check is done in the
 *     caller" and no caller does it;
 *   - A NEGATIVE KEY IS NEVER CHECKED ANYWHERE, because every test is `>=` against the dimension, so
 *     it reaches the array index and throws;
 *   - THE FIRST TWO DIMENSIONS ARE PRE-ALLOCATED and the rest are created on demand in put(), which
 *     is why an empty four-dimensional table already holds a tree of empty arrays and why
 *     getAllLeaves walks it;
 *   - put() BOUNDS-CHECKS ONLY THE NESTED DIMENSIONS TOO, with a message naming the dimension and
 *     its maximum, and refuses a wrong number of keys with a different message;
 *   - THE TABLE SHAPES COME FROM THE COVARIATES: read group is (numReadGroups x 3), quality score is
 *     (numReadGroups x 94 x 3), and each additional covariate is (numReadGroups x 94 x
 *     maximumKeyValue+1 x 3), so the context table is 1012 wide and the cycle table 1002;
 *   - combineTables ADDS THE OTHER TABLE'S DATUMS IN PLACE, calling RecalDatum.combine where a
 *     datum already exists and STORING THE OTHER TABLE'S OBJECT ITSELF where one does not, so the
 *     two tables end up sharing that datum;
 *   - AND safeCombine IS NOT SAFE. It allocates a new set of tables and then combines BOTH
 *     arguments into it, so the first combine moves the left table's datum objects into the new
 *     table and the second combine MUTATES them. Measured: after safeCombine(one, two) the new
 *     table's datum at (0,0) is the same object as one's, and it now holds the sum of both. A port
 *     that copied the datums would leave `one` unchanged and be wrong about what the caller holds
 *     afterwards.
 *
 * Output:
 *
 *     readgroups\t<comma separated read group identifiers>
 *     shape\t<table>\t<comma separated dimensions>
 *     const\t<name>\t<value>
 *     get\t<label>\t<keys>\t<result>
 *     put\t<label>\t<keys>\t<result>
 *     leaf\t<label>\t<keys>\t<datum>
 *     values\t<label>\t<count>
 *     combine\t<label>\t<keys>\t<datum>
 *     shared\t<label>\t<true|false>
 *     error\t<what>\t<exception>\t<message>
 *
 * Usage: RecalibrationTablesDump
 */

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;
import org.broadinstitute.hellbender.utils.recalibration.RecalDatum;
import org.broadinstitute.hellbender.utils.recalibration.RecalUtils;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.utils.recalibration.covariates.Covariate;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;

import java.util.Arrays;

public class RecalibrationTablesDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# RecalibrationTablesDump: NestedIntegerArray and RecalibrationTables");

        shapes();
        lookups();
        insertions();
        traversal();
        combining();
    }

    /** The shapes, which come from the covariates rather than from any constant here. */
    static void shapes() {
        final SAMFileHeader header = ReadFilterDump.header();
        final RecalibrationArgumentCollection rac = new RecalibrationArgumentCollection();
        final StandardCovariateList covariates = new StandardCovariateList(rac, header);
        final RecalibrationTables tables = new RecalibrationTables(covariates);

        // The read groups the shapes are computed from, so a port can build the same covariates
        // without the whole corpus travelling with them.
        System.out.printf("readgroups\t%s%n",
                String.join(",", org.broadinstitute.hellbender.utils.recalibration.covariates
                        .ReadGroupCovariate.getReadGroupIDs(header)));
        System.out.printf("const\tnumTables\t%d%n", tables.numTables());
        System.out.printf("const\teventDimension\t%d%n", 3);
        System.out.printf("const\tnumReadGroups\t%d%n",
                covariates.getReadGroupCovariate().maximumKeyValue() + 1);
        System.out.printf("const\tqualDimension\t%d%n",
                covariates.getQualityScoreCovariate().maximumKeyValue() + 1);

        System.out.printf("shape\treadGroup\t%s%n", join(tables.getReadGroupTable().getDimensions()));
        System.out.printf("shape\tqualityScore\t%s%n",
                join(tables.getQualityScoreTable().getDimensions()));
        int index = 0;
        for (final Covariate covariate : covariates.getAdditionalCovariates()) {
            System.out.printf("shape\t%s\t%s%n", covariate.parseNameForReport(),
                    join(tables.getAdditionalTables().get(index).getDimensions()));
            index++;
        }
        // Which table is which, by the reference-identity tests the class exposes.
        System.out.printf("const\tisReadGroupTable\t%b%n",
                tables.isReadGroupTable(tables.getReadGroupTable()));
        System.out.printf("const\tisQualityScoreTable\t%b%n",
                tables.isQualityScoreTable(tables.getQualityScoreTable()));
        System.out.printf("const\tqualityScoreIsReadGroup\t%b%n",
                tables.isReadGroupTable(tables.getQualityScoreTable()));
        // A freshly made quality score table has the same shape and is NOT the same table.
        System.out.printf("const\tmadeQualityScoreIsQualityScore\t%b%n",
                tables.isQualityScoreTable(tables.makeQualityScoreTable()));
        System.out.printf("const\tisEmpty\t%b%n", tables.isEmpty());
    }

    /**
     * The four ways to read a value, which do not agree.
     */
    static void lookups() {
        final NestedIntegerArray<RecalDatum> table = new NestedIntegerArray<>(3, 4, 5);
        final RecalDatum datum = new RecalDatum(1000L, 10.0, (byte) 30);
        table.put(datum, 1, 2, 3);

        // In range, every way.
        get("in-range-varargs", () -> table.get(1, 2, 3));
        get("in-range-3keys", () -> table.get3Keys(1, 2, 3));
        // Set but read at a neighbouring leaf: null, not an error.
        get("unset-varargs", () -> table.get(1, 2, 4));
        get("unset-3keys", () -> table.get3Keys(1, 2, 4));
        // A branch that was never created: null both ways.
        get("unset-branch-varargs", () -> table.get(0, 0, 0));
        get("unset-branch-3keys", () -> table.get3Keys(0, 0, 0));

        // The divergence: an out-of-range LAST key. The varargs form indexes the leaf array
        // directly and throws; the specialised form checks it and answers null.
        get("last-key-too-big-varargs", () -> table.get(1, 2, 5));
        get("last-key-too-big-3keys", () -> table.get3Keys(1, 2, 5));

        // An out-of-range nested key is null both ways.
        get("first-key-too-big-varargs", () -> table.get(3, 0, 0));
        get("first-key-too-big-3keys", () -> table.get3Keys(3, 0, 0));
        get("second-key-too-big-varargs", () -> table.get(0, 4, 0));
        get("second-key-too-big-3keys", () -> table.get3Keys(0, 4, 0));

        // A negative key is checked by nothing, because every test is `>=`.
        get("negative-first-varargs", () -> table.get(-1, 0, 0));
        get("negative-first-3keys", () -> table.get3Keys(-1, 0, 0));
        get("negative-last-varargs", () -> table.get(1, 2, -1));
        get("negative-last-3keys", () -> table.get3Keys(1, 2, -1));

        // get1Key checks nothing at all, on a one-dimensional table.
        final NestedIntegerArray<RecalDatum> flat = new NestedIntegerArray<>(2);
        flat.put(datum, 0);
        get("flat-1key-in-range", () -> flat.get1Key(0));
        get("flat-1key-unset", () -> flat.get1Key(1));
        get("flat-1key-too-big", () -> flat.get1Key(2));
        get("flat-1key-negative", () -> flat.get1Key(-1));
        get("flat-varargs-too-big", () -> flat.get(2));

        // Two and four dimensions, for the other two specialisations.
        final NestedIntegerArray<RecalDatum> two = new NestedIntegerArray<>(2, 3);
        two.put(datum, 1, 2);
        get("two-2keys", () -> two.get2Keys(1, 2));
        get("two-2keys-last-too-big", () -> two.get2Keys(1, 3));
        get("two-varargs-last-too-big", () -> two.get(1, 3));

        final NestedIntegerArray<RecalDatum> four = new NestedIntegerArray<>(2, 3, 4, 5);
        four.put(datum, 1, 2, 3, 4);
        get("four-4keys", () -> four.get4Keys(1, 2, 3, 4));
        get("four-4keys-last-too-big", () -> four.get4Keys(1, 2, 3, 5));
        get("four-varargs-last-too-big", () -> four.get(1, 2, 3, 5));
        // The third dimension of a four-dimensional table is created on demand, so a key into a
        // branch that put() never touched is null rather than an error.
        get("four-unset-branch", () -> four.get(0, 0, 0, 0));
        get("four-4keys-unset-branch", () -> four.get4Keys(0, 0, 0, 0));
    }

    /** put(), whose checks are the mirror of get()'s and whose messages are exact. */
    static void insertions() {
        final RecalDatum datum = new RecalDatum(1000L, 10.0, (byte) 30);

        put("wrong-key-count-too-few", () -> {
            new NestedIntegerArray<RecalDatum>(3, 4, 5).put(datum, 1, 2);
            return null;
        });
        put("wrong-key-count-too-many", () -> {
            new NestedIntegerArray<RecalDatum>(3, 4, 5).put(datum, 1, 2, 3, 4);
            return null;
        });
        put("first-key-too-big", () -> {
            new NestedIntegerArray<RecalDatum>(3, 4, 5).put(datum, 3, 0, 0);
            return null;
        });
        put("second-key-too-big", () -> {
            new NestedIntegerArray<RecalDatum>(3, 4, 5).put(datum, 0, 4, 0);
            return null;
        });
        // The last key is not checked here either, so this is an index error and not a message.
        put("last-key-too-big", () -> {
            new NestedIntegerArray<RecalDatum>(3, 4, 5).put(datum, 0, 0, 5);
            return null;
        });
        put("negative-first-key", () -> {
            new NestedIntegerArray<RecalDatum>(3, 4, 5).put(datum, -1, 0, 0);
            return null;
        });
        // Overwriting is silent.
        put("overwrite", () -> {
            final NestedIntegerArray<RecalDatum> table = new NestedIntegerArray<>(2, 2);
            table.put(datum, 0, 0);
            table.put(new RecalDatum(7L, 1.0, (byte) 20), 0, 0);
            return table.get(0, 0);
        });
        // Zero dimensions is refused at construction.
        put("no-dimensions", () -> new NestedIntegerArray<RecalDatum>());
        // A zero-length dimension makes an array nothing fits in.
        put("zero-length-dimension", () -> {
            final NestedIntegerArray<RecalDatum> table = new NestedIntegerArray<>(0, 3);
            return table.getAllValues().size();
        });
    }

    /** getAllValues and getAllLeaves, which walk the pre-allocated tree. */
    static void traversal() {
        final NestedIntegerArray<RecalDatum> table = new NestedIntegerArray<>(2, 3, 4, 5);
        System.out.printf("values\tempty-four\t%d%n", table.getAllValues().size());
        System.out.printf("values\tempty-four-leaves\t%d%n", table.getAllLeaves().size());

        table.put(new RecalDatum(10L, 1.0, (byte) 30), 0, 0, 0, 0);
        table.put(new RecalDatum(20L, 2.0, (byte) 30), 1, 2, 3, 4);
        table.put(new RecalDatum(30L, 3.0, (byte) 30), 0, 1, 0, 1);
        System.out.printf("values\tthree-values\t%d%n", table.getAllValues().size());
        for (final NestedIntegerArray.Leaf<RecalDatum> leaf : table.getAllLeaves()) {
            System.out.printf("leaf\tthree-values\t%s\t%s%n", join(leaf.keys), leaf.value);
        }
        // The order getAllValues returns them in, which is the tree walk's order and not the
        // insertion order.
        final StringBuilder order = new StringBuilder();
        for (final RecalDatum value : table.getAllValues()) {
            if (order.length() != 0) {
                order.append(';');
            }
            order.append(value);
        }
        System.out.printf("values\tthree-values-order\t%s%n", order);

        // A one-dimensional table has no nesting to walk.
        final NestedIntegerArray<RecalDatum> flat = new NestedIntegerArray<>(3);
        flat.put(new RecalDatum(5L, 0.0, (byte) 30), 2);
        for (final NestedIntegerArray.Leaf<RecalDatum> leaf : flat.getAllLeaves()) {
            System.out.printf("leaf\tflat\t%s\t%s%n", join(leaf.keys), leaf.value);
        }
    }

    /** combineTables, which shares objects rather than copying them. */
    static void combining() {
        final NestedIntegerArray<RecalDatum> left = new NestedIntegerArray<>(2, 3);
        final NestedIntegerArray<RecalDatum> right = new NestedIntegerArray<>(2, 3);

        // One position in both, so the datums are combined.
        left.put(new RecalDatum(1000L, 10.0, (byte) 30), 0, 0);
        right.put(new RecalDatum(2000L, 20.0, (byte) 20), 0, 0);
        // One position only in the right table, so its datum is MOVED, not copied.
        final RecalDatum onlyRight = new RecalDatum(500L, 5.0, (byte) 25);
        right.put(onlyRight, 1, 1);
        // And one only in the left, which the merge never touches.
        left.put(new RecalDatum(300L, 3.0, (byte) 35), 1, 2);

        RecalUtils.combineTables(left, right);
        for (final NestedIntegerArray.Leaf<RecalDatum> leaf : left.getAllLeaves()) {
            System.out.printf("combine\tmerged\t%s\t%s%n", join(leaf.keys), leaf.value);
        }
        // The two tables now hold the SAME object at that key, so a change to one is a change to
        // both. This is the behaviour a port that copies would lose.
        System.out.printf("shared\tonly-right\t%b%n", left.get(1, 1) == right.get(1, 1));
        System.out.printf("shared\tcombined\t%b%n", left.get(0, 0) == right.get(0, 0));

        // Mismatched shapes are refused, with both shapes in the message.
        try {
            RecalUtils.combineTables(new NestedIntegerArray<RecalDatum>(2, 3),
                    new NestedIntegerArray<RecalDatum>(2, 4));
            System.out.println("error\tcombine-different-shapes\tnone\t-");
        } catch (final Exception e) {
            System.out.printf("error\tcombine-different-shapes\t%s\t%s%n",
                    e.getClass().getSimpleName(), e.getMessage());
        }

        // And the whole-table combine, which checks the table COUNT and not the shapes.
        final SAMFileHeader header = ReadFilterDump.header();
        final RecalibrationArgumentCollection rac = new RecalibrationArgumentCollection();
        final StandardCovariateList covariates = new StandardCovariateList(rac, header);
        final RecalibrationTables one = new RecalibrationTables(covariates);
        final RecalibrationTables two = new RecalibrationTables(covariates);
        one.getReadGroupTable().put(new RecalDatum(100L, 1.0, (byte) 30), 0, 0);
        two.getReadGroupTable().put(new RecalDatum(200L, 2.0, (byte) 30), 0, 0);
        two.getReadGroupTable().put(new RecalDatum(300L, 3.0, (byte) 30), 1, 1);
        final RecalibrationTables combined = RecalibrationTables.safeCombine(one, two);
        for (final NestedIntegerArray.Leaf<RecalDatum> leaf : combined.getReadGroupTable().getAllLeaves()) {
            System.out.printf("combine\tsafeCombine\t%s\t%s%n", join(leaf.keys), leaf.value);
        }
        System.out.printf("const\tsafeCombineIsEmpty\t%b%n", combined.isEmpty());
        // safeCombine reads from both, so the left table is untouched and the right one is not
        // either; the shared datums are the new table's.
        System.out.printf("shared\tsafeCombine-left\t%b%n",
                combined.getReadGroupTable().get(0, 0) == one.getReadGroupTable().get(0, 0));
    }

    interface Attempt {
        Object run() throws Exception;
    }

    static void get(final String label, final Attempt attempt) {
        emit("get", label, attempt);
    }

    static void put(final String label, final Attempt attempt) {
        emit("put", label, attempt);
    }

    static void emit(final String kind, final String label, final Attempt attempt) {
        try {
            System.out.printf("%s\t%s\t%s%n", kind, label, String.valueOf(attempt.run()));
        } catch (final Exception e) {
            System.out.printf("%s\t%s\tE:%s:%s%n", kind, label, e.getClass().getSimpleName(),
                    e.getMessage());
        }
    }

    static String join(final int[] values) {
        return Arrays.toString(values).replace("[", "").replace("]", "").replace(" ", "");
    }
}

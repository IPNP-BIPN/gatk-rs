/*
 * `IntervalArgumentCollection`: the five arguments that turn -L and -XL into the list a walker
 * traverses.
 *
 * The runners read `--intervals` and ignore the other four, so `--interval-padding`,
 * `--interval-exclusion-padding`, `--interval-set-rule` and `--exclude-intervals` change nothing
 * in the port and change the answer in the reference. This measures what they change.
 *
 * Eight behaviours this is built to catch.
 *
 *   - PADDING IS APPLIED PER -L ARGUMENT, before the set operator, and each padded batch is sorted
 *     and merged with ALL whatever the merging rule is;
 *   - THE SET OPERATOR FOLDS ONE ARGUMENT AT A TIME into an accumulator, so INTERSECTION over three
 *     arguments is not the intersection of all three at once;
 *   - AN EMPTY LIST SHORT-CIRCUITS the operator: intersecting with nothing returns the other side
 *     rather than nothing, which is why the FIRST -L is never intersected with anything;
 *   - AN EMPTY INTERSECTION IS A REFUSAL, not an empty traversal;
 *   - -XL WITH NO -L means the whole reference, contig by contig, and not an empty set;
 *   - -XL THAT REMOVES EVERYTHING IS A REFUSAL, with its own message;
 *   - THE EXCLUSION PADDING IS A SEPARATE ARGUMENT from the inclusion padding, and each applies to
 *     its own side only;
 *   - AND PADDING IS CLAMPED TO THE CONTIG, so a padded interval never runs past the dictionary's
 *     length or below one.
 *
 * Output:
 *
 *     intervals\t<case>\t<contig:start-end, space separated, or (empty)>
 *     unmapped\t<case>\t<whether the traversal includes unmapped records>
 *     error\t<case>\t<exception class>: <message>
 *
 * Usage: IntervalArgumentsDump
 */

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.engine.TraversalParameters;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalSetRule;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class IntervalArgumentsDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** The dictionary every case resolves against: three contigs of different lengths. */
    static SAMSequenceDictionary dictionary() {
        return new SAMSequenceDictionary(Arrays.asList(
                new SAMSequenceRecord("chr1", 1000),
                new SAMSequenceRecord("chr2", 500),
                new SAMSequenceRecord("chr10", 200)));
    }

    static void set(final Object target, final String field, final Object value) throws Exception {
        Field declared = null;
        for (Class<?> type = target.getClass(); type != null; type = type.getSuperclass()) {
            try {
                declared = type.getDeclaredField(field);
                break;
            } catch (final NoSuchFieldException ignored) {
                // keep walking up: the strings live on the subclass and the rest on the base.
            }
        }
        declared.setAccessible(true);
        if (value instanceof List) {
            @SuppressWarnings("unchecked")
            final List<String> into = (List<String>) declared.get(target);
            into.addAll((List<String>) value);
        } else {
            declared.set(target, value);
        }
    }

    /** One command line's worth of interval arguments, resolved. */
    static void run(final String name, final List<String> include, final List<String> exclude,
                    final IntervalSetRule setRule, final IntervalMergingRule mergingRule,
                    final int padding, final int exclusionPadding) {
        try {
            final IntervalArgumentCollection collection = new OptionalIntervalArgumentCollection();
            set(collection, "intervalStrings", include);
            set(collection, "excludeIntervalStrings", exclude);
            set(collection, "intervalSetRule", setRule);
            set(collection, "intervalMergingRule", mergingRule);
            set(collection, "intervalPadding", padding);
            set(collection, "intervalExclusionPadding", exclusionPadding);

            final TraversalParameters parameters = collection.getTraversalParameters(dictionary());
            final List<String> rendered = new ArrayList<>();
            for (final SimpleInterval interval : parameters.getIntervalsForTraversal()) {
                rendered.add(interval.getContig() + ":" + interval.getStart() + "-" + interval.getEnd());
            }
            emit("intervals", name, rendered.isEmpty() ? "(empty)" : String.join(" ", rendered));
            emit("unmapped", name, Boolean.toString(parameters.traverseUnmappedReads()));
        } catch (final Exception e) {
            emit("error", name, e.getClass().getName() + ": " + e.getMessage());
        }
    }

    static List<String> list(final String... values) {
        return Arrays.asList(values);
    }

    static final List<String> NONE = list();
    static final IntervalSetRule UNION = IntervalSetRule.UNION;
    static final IntervalSetRule INTERSECTION = IntervalSetRule.INTERSECTION;
    static final IntervalMergingRule ALL = IntervalMergingRule.ALL;
    static final IntervalMergingRule OVERLAPPING = IntervalMergingRule.OVERLAPPING_ONLY;

    public static void main(final String[] args) {
        // The baseline: one interval, nothing else set.
        run("one-interval", list("chr1:100-200"), NONE, UNION, ALL, 0, 0);

        // Padding, and padding clamped by the contig at both ends.
        run("padded", list("chr1:100-200"), NONE, UNION, ALL, 50, 0);
        run("padded-past-the-start", list("chr1:10-20"), NONE, UNION, ALL, 50, 0);
        run("padded-past-the-end", list("chr2:480-490"), NONE, UNION, ALL, 50, 0);
        // Padding that makes two separate intervals touch, and therefore merge.
        run("padding-merges-two", list("chr1:100-200", "chr1:260-300"), NONE, UNION, ALL, 30, 0);

        // The merging rule, on intervals that are adjacent rather than overlapping.
        run("adjacent-all", list("chr1:100-200", "chr1:201-300"), NONE, UNION, ALL, 0, 0);
        run("adjacent-overlapping-only", list("chr1:100-200", "chr1:201-300"), NONE, UNION, OVERLAPPING, 0, 0);

        // The set rule over two arguments, and over three, which is folded one at a time.
        run("union-two", list("chr1:100-200", "chr1:150-300"), NONE, UNION, ALL, 0, 0);
        run("intersection-two", list("chr1:100-200", "chr1:150-300"), NONE, INTERSECTION, ALL, 0, 0);
        run("intersection-three", list("chr1:100-400", "chr1:200-500", "chr1:300-600"),
                NONE, INTERSECTION, ALL, 0, 0);
        // An intersection that is empty, which is a refusal rather than an empty traversal.
        run("intersection-empty", list("chr1:100-200", "chr1:300-400"), NONE, INTERSECTION, ALL, 0, 0);
        // Intersection across contigs, which cannot overlap at all.
        run("intersection-across-contigs", list("chr1:100-200", "chr2:100-200"),
                NONE, INTERSECTION, ALL, 0, 0);

        // Exclusion: a hole in the middle, a bite off one end, and the whole thing.
        run("exclude-middle", list("chr1:100-300"), list("chr1:150-200"), UNION, ALL, 0, 0);
        run("exclude-prefix", list("chr1:100-300"), list("chr1:50-150"), UNION, ALL, 0, 0);
        run("exclude-everything", list("chr1:100-200"), list("chr1:1-1000"), UNION, ALL, 0, 0);
        // Exclusion padding, which is its own argument and applies to -XL only.
        run("exclusion-padded", list("chr1:100-300"), list("chr1:180-190"), UNION, ALL, 0, 20);
        run("both-paddings", list("chr1:100-300"), list("chr1:180-190"), UNION, ALL, 10, 20);

        // -XL with no -L at all, which is the whole reference minus the exclusion.
        run("exclude-only", NONE, list("chr1:1-900"), UNION, ALL, 0, 0);
        run("exclude-only-whole-contig", NONE, list("chr2"), UNION, ALL, 0, 0);

        // `unmapped`, which is a traversal flag rather than an interval, and is refused on -XL.
        run("unmapped-included", list("chr1:100-200", "unmapped"), NONE, UNION, ALL, 0, 0);
        run("unmapped-only", list("unmapped"), NONE, UNION, ALL, 0, 0);
        run("unmapped-excluded", list("chr1:100-200"), list("unmapped"), UNION, ALL, 0, 0);

        // A whole contig by name, and the dictionary's own order, which is not alphabetical.
        run("whole-contig", list("chr2"), NONE, UNION, ALL, 0, 0);
        run("three-contigs-out-of-order", list("chr10:1-10", "chr2:1-10", "chr1:1-10"),
                NONE, UNION, ALL, 0, 0);

        System.out.print(buf);
    }
}

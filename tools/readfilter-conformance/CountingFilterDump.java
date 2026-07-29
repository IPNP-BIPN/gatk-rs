/*
 * CountingReadFilter's counts and its summary text, taken from the reference.
 *
 * This is the first *output text* gatk-rs has to reproduce. Every GATK tool that reads reads ends
 * its run by printing this block, so the summary is not a debugging aid: it is part of the bytes a
 * run produces, down to the trailing space that CountingReadFilter puts after the filter name in
 * one branch and not in the other.
 *
 * It is also the only place where the *order* of a conjunction becomes observable. WellformedReadFilter
 * is eight filters and'ed together; the boolean it returns cannot tell you which of the eight
 * rejected a read, but the counts can, because `and` short-circuits and only the first failing
 * filter increments.
 *
 * Output:
 *
 *     composition\t<id>\t<filter.getName()>
 *     summary\t<id>\t<getSummaryLine(), newlines escaped as \n>
 *
 * The composition row is the reference's own name for the tree, and the port rebuilds the tree by
 * parsing it, so the two sides cannot disagree about what was composed.
 *
 * The corpus is ReadFilterDump's, unchanged, for the same reason: one corpus, judged once.
 *
 * Usage: CountingFilterDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class CountingFilterDump {

    /**
     * The compositions, each exercising one shape of the summary.
     *
     * Only filters that decide for every record in the corpus appear here. A filter that throws
     * would stop the walk part way and leave the counters describing a prefix of the corpus, which
     * is a real behaviour but a different one, and not the one this suite is about.
     */
    static Map<String, CountingReadFilter> compositions(final SAMFileHeader header) {
        final Map<String, CountingReadFilter> map = new LinkedHashMap<>();

        map.put("single", new CountingReadFilter(ReadFilterLibrary.MAPPED));

        // fromList is what the engine itself builds from --read-filter arguments: a left-nested
        // chain of ANDs, which is the shape the simplified summary flattens.
        final List<ReadFilter> chain = Arrays.asList(
                ReadFilterLibrary.MAPPED,
                ReadFilterLibrary.NOT_DUPLICATE,
                ReadFilterLibrary.PRIMARY_LINE,
                ReadFilterLibrary.GOOD_CIGAR);
        map.put("from-list", CountingReadFilter.fromList(chain, header));

        // One filter that rejects nothing, to hold the "0 read(s) filtered by" branch, which
        // formats its count as a literal rather than through Long.toString.
        map.put("nothing-filtered", CountingReadFilter.fromList(
                Arrays.asList(ReadFilterLibrary.ALLOW_ALL_READS, ReadFilterLibrary.MAPPED), header));

        // An OR inside the tree: the flattening bails out on anything that is not an AND, so this
        // takes the nested, indented branch instead.
        final CountingReadFilter nested = new CountingReadFilter(ReadFilterLibrary.MAPPED)
                .and(new CountingReadFilter(ReadFilterLibrary.NOT_DUPLICATE)
                        .or(new CountingReadFilter(ReadFilterLibrary.PRIMARY_LINE)));
        map.put("with-or", nested);

        // A negation, whose name is composed rather than a class name.
        map.put("negated", new CountingReadFilter(ReadFilterLibrary.MAPPED)
                .and(new CountingReadFilter(ReadFilterLibrary.NOT_DUPLICATE).negate()));

        // WellformedReadFilter, the conjunction whose order only the counts can reveal. It is a
        // single ReadFilter to the engine, so wrapping it counts the whole conjunction as one; the
        // composition below spells the same eight filters out in the same order, and the two
        // summaries agreeing on the total is what says the order was ported and not guessed.
        map.put("wellformed-opaque", new CountingReadFilter(new WellformedReadFilter(header)));
        map.put("wellformed-spelled-out", CountingReadFilter.fromList(Arrays.asList(
                ReadFilterLibrary.VALID_ALIGNMENT_START,
                ReadFilterLibrary.VALID_ALIGNMENT_END,
                new org.broadinstitute.hellbender.engine.filters.AlignmentAgreesWithHeaderReadFilter(),
                ReadFilterLibrary.HAS_READ_GROUP,
                ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS,
                ReadFilterLibrary.READLENGTH_EQUALS_CIGARLENGTH,
                ReadFilterLibrary.SEQ_IS_STORED,
                new ReadFilterLibrary.CigarContainsNoNOperator()), header));

        return map;
    }

    public static void main(final String[] args) throws Exception {
        final SAMFileHeader header = ReadFilterDump.header();
        final List<SAMRecord> corpus = ReadFilterDump.corpus(header);

        System.out.println("# CountingFilterDump: counts and summary text, from the reference");
        // The corpus travels here too rather than being read from the other golden: this file has
        // to say which records produced these counts, or the counts describe nothing in particular.
        ReadFilterDump.printCorpus(header, corpus);

        for (final Map.Entry<String, CountingReadFilter> entry : compositions(header).entrySet()) {
            final CountingReadFilter filter = entry.getValue();
            filter.setHeader(header);
            final StringBuilder decisions = new StringBuilder();
            for (final SAMRecord record : corpus) {
                final GATKRead read = new SAMRecordToGATKReadAdapter(record);
                decisions.append(filter.test(read) ? '1' : '0');
            }
            System.out.printf("composition\t%s\t%s%n", entry.getKey(), filter.getName());
            System.out.printf("decisions\t%s\t%s%n", entry.getKey(), decisions);
            System.out.printf("summary\t%s\t%s%n",
                    entry.getKey(), filter.getSummaryLine().replace("\n", "\\n"));
        }
    }
}

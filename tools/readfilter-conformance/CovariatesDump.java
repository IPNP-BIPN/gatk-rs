/*
 * The four BQSR covariates, taken from the reference.
 *
 * A recalibration table is indexed by covariate keys, so what BaseRecalibrator counts and what
 * ApplyBQSR looks up are both decided here. Four covariates, one matrix per read of shape
 * (event type) x (read position) x (covariate), and every cell is an integer key.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE KEY MATRIX IS REUSED ACROSS READS OF THE SAME LENGTH AND IT DOES NOT LEAK.
 *     PerReadCovariateMatrix takes a CovariateKeyCache and, on a hit, hands back THE SAME int[][][]
 *     a previous read filled in, with nothing clearing it in between. That looks like stale data
 *     waiting to happen, and it is not: the dump runs the whole corpus twice, once with one shared
 *     cache as BQSR does and once with a fresh cache per read, and prints both matrices so the two
 *     can be compared cell by cell. They agree everywhere, because every covariate writes a key at
 *     every position of the read and ContextCovariate zeroes the whole array first on the one path
 *     that would otherwise leave a gap, the low-quality clip shortening the read. A port is
 *     therefore free to allocate per read, and this is the measurement that says so;
 *   - THE READ GROUP IS IDENTIFIED BY PU AND NOT BY ID, falling back to ID only when PU is absent,
 *     so two read groups with the same PU would collide and a table keyed by ID would be a different
 *     table;
 *   - A READ WHOSE GROUP IS NOT IN THE HEADER IS A NullPointerException, not the -1 key. The -1 is
 *     documented as the missing-read-group code, and the path that produces it is
 *     keyForReadGroup on an identifier the covariate's own table does not hold; a read whose RG
 *     names a group the HEADER does not declare never gets that far, because
 *     ReadUtils.getSAMReadGroupRecord answers null and getReadGroupIdentifier dereferences it. The
 *     dump carries both: a read like that, and keyFromValue("nonesuch") answering -1;
 *   - THE CONTEXT IS OF THE READ'S OWN BASES, REVERSE-COMPLEMENTED FOR A NEGATIVE-STRAND READ, and
 *     the low-quality tail is overwritten with N first. The first contextSize-1 positions have no
 *     context and get -1;
 *   - AN N IN THE CONTEXT POISONS THE FOLLOWING contextSize-1 POSITIONS, and the recovery loop that
 *     rebuilds the key after the first -1 walks BACKWARDS from the penalty position, which is a
 *     different thing from restarting;
 *   - THE CONTEXT KEY CARRIES ITS LENGTH IN THE LOW FOUR BITS, so a key is not just packed bases and
 *     contextFromKey reads the length back out of it;
 *   - THE CYCLE IS SIGNED AND THE SIGN IS THE LOW BIT of the key, so cycle -3 and cycle 3 are
 *     different keys; a second-of-pair read counts the other way, and a negative-strand read counts
 *     from the far end;
 *   - INDEL CYCLE KEYS ARE -1 WITHIN FOUR BASES OF EITHER END, which is a cushion the substitution
 *     keys do not have;
 *   - AND A READ WITH BASES AND NO QUALITIES IS AN EXCEPTION, not a covariate value. The quality
 *     covariate reads the base quality COUNT and writes nothing, and then ContextCovariate's
 *     low-quality clip indexes an empty quality array by the read's length and throws
 *     ArrayIndexOutOfBoundsException. Nothing in the covariates guards it.
 *
 * Output:
 *
 *     const\t<name>\t<value>
 *     covariate\t<index>\t<simpleName>\t<parseNameForReport>\t<maximumKeyValue>
 *     names\t<covariateNames>
 *     rgids\t<comma separated identifiers from the header>
 *     matrix\t<cache>\t<readIndex>\t<eventType>\t<covIndex>\t<comma separated keys>
 *     clipped\t<readIndex>\t<lowQualTail>\t<stranded clipped bases>
 *     context\t<dna>\t<key>\t<contextFromKey>
 *     cycle\t<cycle>\t<maxCycle>\t<key>\t<cycleFromKey>
 *     format\t<covariate>\t<key>\t<formatted>
 *     fromvalue\t<covariate>\t<value>\t<key>
 *     error\t<what>\t<exception>\t<message>
 *
 * Usage: CovariatesDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.recalibration.EventType;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.covariates.ContextCovariate;
import org.broadinstitute.hellbender.utils.recalibration.covariates.Covariate;
import org.broadinstitute.hellbender.utils.recalibration.covariates.CovariateKeyCache;
import org.broadinstitute.hellbender.utils.recalibration.covariates.CycleCovariate;
import org.broadinstitute.hellbender.utils.recalibration.covariates.PerReadCovariateMatrix;
import org.broadinstitute.hellbender.utils.recalibration.covariates.ReadGroupCovariate;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

public class CovariatesDump {

    static final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();

    public static void main(final String[] args) throws Exception {
        final SAMFileHeader header = ReadFilterDump.header();
        final List<SAMRecord> corpus = corpus(header);

        System.out.println("# CovariatesDump: the four BQSR covariates");
        ReadFilterDump.printCorpus(header, corpus);

        constants(header);

        final StandardCovariateList covariates = new StandardCovariateList(RAC, header);
        describe(covariates);

        // The same reads twice. "shared" is what BQSR does: one CovariateKeyCache for the whole
        // traversal, so reads of equal length are handed the same array. "fresh" gives every read
        // its own, which is the only way to see what the shared run is carrying over.
        matrices(header, corpus, covariates, "shared", new CovariateKeyCache(), true);
        matrices(header, corpus, covariates, "fresh", null, true);
        // And once without indel values, which takes the other half of every duplicated loop.
        matrices(header, corpus, covariates, "fresh-no-indels", null, false);

        clippedBases(header, corpus);
        contexts();
        cycles();
        formatting(covariates);
        errors(header, covariates);
    }

    static void constants(final SAMFileHeader header) {
        System.out.printf("const\tMISSING_READ_GROUP_KEY\t%d%n", ReadGroupCovariate.MISSING_READ_GROUP_KEY);
        System.out.printf("const\tUNKNOWN_OR_ERROR_CONTEXT_CODE\t%d%n", ContextCovariate.UNKNOWN_OR_ERROR_CONTEXT_CODE);
        System.out.printf("const\tCUSHION_FOR_INDELS\t%d%n", CycleCovariate.CUSHION_FOR_INDELS);
        System.out.printf("const\tREAD_GROUP_COVARIATE_DEFAULT_INDEX\t%d%n", StandardCovariateList.READ_GROUP_COVARIATE_DEFAULT_INDEX);
        System.out.printf("const\tBASE_QUALITY_COVARIATE_DEFAULT_INDEX\t%d%n", StandardCovariateList.BASE_QUALITY_COVARIATE_DEFAULT_INDEX);
        System.out.printf("const\tCONTEXT_COVARIATE_DEFAULT_INDEX\t%d%n", StandardCovariateList.CONTEXT_COVARIATE_DEFAULT_INDEX);
        System.out.printf("const\tCYCLE_COVARIATE_DEFAULT_INDEX\t%d%n", StandardCovariateList.CYCLE_COVARIATE_DEFAULT_INDEX);
        System.out.printf("const\tNUM_REQUIRED_COVARITES\t%d%n", StandardCovariateList.NUM_REQUIRED_COVARITES);
        System.out.printf("const\tMISMATCHES_CONTEXT_SIZE\t%d%n", RAC.MISMATCHES_CONTEXT_SIZE);
        System.out.printf("const\tINDELS_CONTEXT_SIZE\t%d%n", RAC.INDELS_CONTEXT_SIZE);
        System.out.printf("const\tMAXIMUM_CYCLE_VALUE\t%d%n", RAC.MAXIMUM_CYCLE_VALUE);
        System.out.printf("const\tLOW_QUAL_TAIL\t%d%n", RAC.LOW_QUAL_TAIL);
        // PU and not ID, which is the whole point of getReadGroupIdentifier.
        System.out.printf("rgids\t%s%n", String.join(",", ReadGroupCovariate.getReadGroupIDs(header)));
        for (final SAMReadGroupRecord rg : header.getReadGroups()) {
            System.out.printf("rgidentifier\t%s\t%s\t%s%n", rg.getId(),
                    String.valueOf(rg.getPlatformUnit()), ReadGroupCovariate.getReadGroupIdentifier(rg));
        }
        // A read group with no PU falls back to its ID.
        final SAMReadGroupRecord noPlatformUnit = new SAMReadGroupRecord("rg-no-pu");
        System.out.printf("rgidentifier\t%s\t%s\t%s%n", noPlatformUnit.getId(),
                String.valueOf(noPlatformUnit.getPlatformUnit()),
                ReadGroupCovariate.getReadGroupIdentifier(noPlatformUnit));
    }

    static void describe(final StandardCovariateList covariates) {
        System.out.printf("names\t%s%n", covariates.covariateNames());
        System.out.printf("const\tsize\t%d%n", covariates.size());
        System.out.printf("const\tnumberOfSpecialCovariates\t%d%n", covariates.numberOfSpecialCovariates());
        for (int i = 0; i < covariates.size(); i++) {
            final Covariate covariate = covariates.get(i);
            System.out.printf("covariate\t%d\t%s\t%s\t%d\t%d%n", i,
                    covariate.getClass().getSimpleName(), covariate.parseNameForReport(),
                    covariate.maximumKeyValue(), covariates.indexByClass(covariate.getClass()));
        }
        // Looked up by the name the report writes, which is the class name minus "Covariate".
        for (final String name : new String[] {"ReadGroup", "QualityScore", "Context", "Cycle", "Nonesuch"}) {
            final Covariate found = covariates.getCovariateByParsedName(name);
            System.out.printf("byname\t%s\t%s%n", name,
                    found == null ? "null" : found.getClass().getSimpleName());
        }
    }

    /**
     * Every key of every read, under one cache policy.
     *
     * @param cache the shared cache, or null for a fresh one per read
     */
    static void matrices(final SAMFileHeader header, final List<SAMRecord> corpus,
                         final StandardCovariateList covariates, final String label,
                         final CovariateKeyCache cache, final boolean recordIndelValues) {
        for (int i = 0; i < corpus.size(); i++) {
            final GATKRead read = new SAMRecordToGATKReadAdapter(corpus.get(i));
            final int readLength = read.getLength();
            if (readLength == 0) {
                System.out.printf("matrix\t%s\t%d\t-\t-\tempty%n", label, i);
                continue;
            }
            final CovariateKeyCache keysCache = cache == null ? new CovariateKeyCache() : cache;
            final PerReadCovariateMatrix matrix =
                    new PerReadCovariateMatrix(readLength, covariates.size(), keysCache);
            try {
                covariates.populatePerReadCovariateMatrix(read, header, matrix, recordIndelValues);
            } catch (final Exception e) {
                System.out.printf("matrix\t%s\t%d\t-\t-\tE:%s:%s%n", label, i,
                        e.getClass().getSimpleName(), e.getMessage());
                continue;
            }
            for (final EventType event : EventType.values()) {
                final int[][] byPosition = matrix.getMatrixForErrorModel(event);
                for (int cov = 0; cov < covariates.size(); cov++) {
                    final StringBuilder keys = new StringBuilder();
                    for (int position = 0; position < readLength; position++) {
                        if (position > 0) {
                            keys.append(',');
                        }
                        keys.append(byPosition[position][cov]);
                    }
                    System.out.printf("matrix\t%s\t%d\t%s\t%d\t%s%n", label, i, event, cov, keys);
                }
            }
        }
    }

    /**
     * The bases the context covariate actually sees: low-quality tail overwritten with N, then
     * reverse-complemented for a negative-strand read.
     */
    static void clippedBases(final SAMFileHeader header, final List<SAMRecord> corpus) throws Exception {
        final Method method = ContextCovariate.class.getDeclaredMethod(
                "getStrandedClippedBytes", GATKRead.class, byte.class);
        method.setAccessible(true);
        for (int i = 0; i < corpus.size(); i++) {
            final GATKRead read = new SAMRecordToGATKReadAdapter(corpus.get(i));
            for (final byte lowQual : new byte[] {0, 2, 20, 30}) {
                try {
                    final byte[] bases = (byte[]) method.invoke(null, read, lowQual);
                    System.out.printf("clipped\t%d\t%d\t%s%n", i, lowQual,
                            new String(bases, StandardCharsets.UTF_8));
                } catch (final InvocationTargetException e) {
                    // A read with bases and no qualities: clipLowQualEnds indexes the quality array
                    // by the read's length and runs off the end of an empty one. The covariate has
                    // no guard, so this is a real end of the code path and not a fixture mistake.
                    final Throwable cause = e.getCause();
                    System.out.printf("clipped\t%d\t%d\tE:%s:%s%n", i, lowQual,
                            cause.getClass().getSimpleName(), cause.getMessage());
                }
            }
        }
    }

    /** The context key's encoding, including the length it carries in its low four bits. */
    static void contexts() {
        final String[] sequences = {
                "A", "C", "G", "T", "AA", "AC", "TT", "ACG", "ACGT", "TTTTTTTTTTTTT",
                // Non-ACGT, which makes the key -1 wherever it appears.
                "N", "AN", "NA", "ACN", "a", "ac",
                // The maximum context the covariate allows, and one base of each at the ends.
                "AAAAAAAAAAAAA", "ACGTACGTACGTA",
        };
        for (final String dna : sequences) {
            final int key = ContextCovariate.keyFromContext(dna);
            String back;
            try {
                back = key < 0 ? "E" : ContextCovariate.contextFromKey(key);
            } catch (final Exception e) {
                back = "E:" + e.getClass().getSimpleName();
            }
            System.out.printf("context\t%s\t%d\t%s%n", dna, key, back);
        }
        // A key whose length nibble says more bases than the key holds, which is what an unchecked
        // integer from a table would look like.
        for (final int key : new int[] {0, 1, 2, 16, 17, 4095, 65535}) {
            System.out.printf("contextfromkey\t%d\t%s%n", key, ContextCovariate.contextFromKey(key));
        }
    }

    /** The cycle key's encoding, whose sign lives in the low bit. */
    static void cycles() {
        final int[] values = {0, 1, -1, 2, -2, 3, 250, -250, 500, -500};
        for (final int cycle : values) {
            final int key = CycleCovariate.keyFromCycle(cycle, RAC.MAXIMUM_CYCLE_VALUE);
            System.out.printf("cycle\t%d\t%d\t%d\t%d%n", cycle, RAC.MAXIMUM_CYCLE_VALUE, key,
                    CycleCovariate.cycleFromKey(key));
        }
        // Past the maximum, which is a UserException rather than a wrapped key.
        for (final int cycle : new int[] {501, -501}) {
            try {
                CycleCovariate.keyFromCycle(cycle, RAC.MAXIMUM_CYCLE_VALUE);
                System.out.printf("cycle\t%d\t%d\tnone\tnone%n", cycle, RAC.MAXIMUM_CYCLE_VALUE);
            } catch (final Exception e) {
                System.out.printf("error\tkeyFromCycle@%d\t%s\t%s%n", cycle,
                        e.getClass().getSimpleName(), e.getMessage());
            }
        }
        // The decoder on keys that no encoder produces.
        for (final int key : new int[] {0, 1, 2, 3, 1000, 1001}) {
            System.out.printf("cyclefromkey\t%d\t%d%n", key, CycleCovariate.cycleFromKey(key));
        }
    }

    /** formatKey and keyFromValue, which are what the recalibration report is written and read by. */
    static void formatting(final StandardCovariateList covariates) {
        final Covariate readGroup = covariates.getReadGroupCovariate();
        final Covariate quality = covariates.getQualityScoreCovariate();
        final List<Covariate> additional = new ArrayList<>();
        covariates.getAdditionalCovariates().forEach(additional::add);
        final Covariate context = additional.get(0);
        final Covariate cycle = additional.get(1);

        for (final int key : new int[] {0, 1, 2}) {
            System.out.printf("format\tReadGroup\t%d\t%s%n", key, readGroup.formatKey(key));
        }
        for (final int key : new int[] {0, 2, 30, 93}) {
            System.out.printf("format\tQualityScore\t%d\t%s%n", key, quality.formatKey(key));
        }
        for (final int key : new int[] {-1, 18, 33, 4095}) {
            System.out.printf("format\tContext\t%d\t%s%n", key, String.valueOf(context.formatKey(key)));
        }
        for (final int key : new int[] {0, 1, 2, 3, 1000, 1001}) {
            System.out.printf("format\tCycle\t%d\t%s%n", key, cycle.formatKey(key));
        }

        for (final String value : new String[] {"unit-rg1", "unit-rg2", "unit-rg3", "nonesuch"}) {
            System.out.printf("fromvalue\tReadGroup\t%s\t%d%n", value, readGroup.keyFromValue(value));
        }
        // A String, a Long and a Byte all reach different branches of the same method.
        System.out.printf("fromvalue\tQualityScore\tstring:30\t%d%n", quality.keyFromValue("30"));
        System.out.printf("fromvalue\tQualityScore\tlong:30\t%d%n", quality.keyFromValue(30L));
        System.out.printf("fromvalue\tQualityScore\tbyte:30\t%d%n", quality.keyFromValue((byte) 30));
        System.out.printf("fromvalue\tQualityScore\tbyte:-1\t%d%n", quality.keyFromValue((byte) -1));
        for (final String value : new String[] {"AC", "ACG", "TT"}) {
            System.out.printf("fromvalue\tContext\t%s\t%d%n", value, context.keyFromValue(value));
        }
        System.out.printf("fromvalue\tCycle\tstring:3\t%d%n", cycle.keyFromValue("3"));
        System.out.printf("fromvalue\tCycle\tstring:-3\t%d%n", cycle.keyFromValue("-3"));
        System.out.printf("fromvalue\tCycle\tinteger:3\t%d%n", cycle.keyFromValue(3));
        System.out.printf("fromvalue\tCycle\tinteger:-3\t%d%n", cycle.keyFromValue(-3));
    }

    /** Every argument the covariates refuse, and the words they refuse it in. */
    static void errors(final SAMFileHeader header, final StandardCovariateList covariates) {
        // A context size above the maximum, and a non-positive one.
        for (final int size : new int[] {14, 0, -1}) {
            final RecalibrationArgumentCollection rac = new RecalibrationArgumentCollection();
            rac.MISMATCHES_CONTEXT_SIZE = size;
            try {
                new ContextCovariate(rac);
                System.out.printf("error\tmismatches-context-size@%d\tnone\t-%n", size);
            } catch (final Exception e) {
                System.out.printf("error\tmismatches-context-size@%d\t%s\t%s%n", size,
                        e.getClass().getSimpleName(), e.getMessage());
            }
        }
        for (final int size : new int[] {14, 0}) {
            final RecalibrationArgumentCollection rac = new RecalibrationArgumentCollection();
            rac.INDELS_CONTEXT_SIZE = size;
            try {
                new ContextCovariate(rac);
                System.out.printf("error\tindels-context-size@%d\tnone\t-%n", size);
            } catch (final Exception e) {
                System.out.printf("error\tindels-context-size@%d\t%s\t%s%n", size,
                        e.getClass().getSimpleName(), e.getMessage());
            }
        }
        // A key the read group does not know, which is an exception on the way out and -1 on the
        // way in. The asymmetry is the point.
        try {
            covariates.getReadGroupCovariate().formatKey(99);
            System.out.println("error\treadgroup-format-unknown\tnone\t-");
        } catch (final Exception e) {
            System.out.printf("error\treadgroup-format-unknown\t%s\t%s%n",
                    e.getClass().getSimpleName(), e.getMessage());
        }
        // A negative context key, which contextFromKey refuses outright.
        try {
            ContextCovariate.contextFromKey(-1);
            System.out.println("error\tcontext-from-negative\tnone\t-");
        } catch (final Exception e) {
            System.out.printf("error\tcontext-from-negative\t%s\t%s%n",
                    e.getClass().getSimpleName(), e.getMessage());
        }
        // A read longer than the maximum cycle, which is where keyFromCycle throws during a run.
        final SAMRecord longRead = ReadFilterDump.read(header, "too_many_cycles", 0, 0, 100, 60,
                "6M", 0, 200, 100, true);
        longRead.setReadString("ACGTAC");
        longRead.setBaseQualityString("IIIIII");
        final RecalibrationArgumentCollection tiny = new RecalibrationArgumentCollection();
        tiny.MAXIMUM_CYCLE_VALUE = 3;
        final StandardCovariateList tinyList = new StandardCovariateList(tiny, header);
        final GATKRead read = new SAMRecordToGATKReadAdapter(longRead);
        final PerReadCovariateMatrix matrix =
                new PerReadCovariateMatrix(read.getLength(), tinyList.size(), new CovariateKeyCache());
        try {
            tinyList.populatePerReadCovariateMatrix(read, header, matrix, true);
            System.out.println("error\tcycle-past-maximum\tnone\t-");
        } catch (final Exception e) {
            System.out.printf("error\tcycle-past-maximum\t%s\t%s%n",
                    e.getClass().getSimpleName(), e.getMessage());
        }
    }

    /**
     * The shared read-filter corpus, plus reads that exist only to separate the covariates.
     */
    static List<SAMRecord> corpus(final SAMFileHeader header) {
        final List<SAMRecord> out = new ArrayList<>(ReadFilterDump.corpus(header));

        // A low-quality tail on each end, which is what the context covariate clips to N.
        final SAMRecord lowQualEnds = ReadFilterDump.read(header, "cov_low_qual_ends", 0, 0, 100, 60,
                "10M", 0, 200, 100, true);
        lowQualEnds.setReadString("ACGTACGTAC");
        lowQualEnds.setBaseQualityString("!!IIIIII!!");
        out.add(lowQualEnds);

        // An N in the middle: the key is -1 there and for the following contextSize-1 positions,
        // and the recovery loop after it walks backwards.
        final SAMRecord withN = ReadFilterDump.read(header, "cov_n_in_middle", 0, 0, 100, 60,
                "10M", 0, 200, 100, true);
        withN.setReadString("ACGTNACGTA");
        withN.setBaseQualityString("IIIIIIIIII");
        out.add(withN);

        // An N in the FIRST context, which is the branch that sets currentNPenalty before the loop.
        final SAMRecord nAtStart = ReadFilterDump.read(header, "cov_n_at_start", 0, 0, 100, 60,
                "10M", 0, 200, 100, true);
        nAtStart.setReadString("NACGTACGTA");
        nAtStart.setBaseQualityString("IIIIIIIIII");
        out.add(nAtStart);

        // The same read on the negative strand, so the context is of the reverse complement and the
        // cycle counts from the far end.
        final SAMRecord reverse = ReadFilterDump.read(header, "cov_reverse", 16, 0, 100, 60,
                "10M", 0, 200, 100, true);
        reverse.setReadString("ACGTNACGTA");
        reverse.setBaseQualityString("IIIIIIIIII");
        out.add(reverse);

        // Second of pair, which flips the sign of every cycle.
        final SAMRecord secondOfPair = ReadFilterDump.read(header, "cov_second_of_pair",
                0x1 | 0x80, 0, 100, 60, "10M", 0, 200, 100, true);
        secondOfPair.setReadString("ACGTACGTAC");
        secondOfPair.setBaseQualityString("IIIIIIIIII");
        out.add(secondOfPair);

        // Second of pair AND negative strand, which is the fourth corner of the cycle sign.
        final SAMRecord secondReverse = ReadFilterDump.read(header, "cov_second_reverse",
                0x1 | 0x80 | 0x10, 0, 100, 60, "10M", 0, 200, 100, true);
        secondReverse.setReadString("ACGTACGTAC");
        secondReverse.setBaseQualityString("IIIIIIIIII");
        out.add(secondReverse);

        // Nine bases: shorter than twice the indel cushion plus one, so EVERY indel cycle key is -1.
        final SAMRecord shortRead = ReadFilterDump.read(header, "cov_short", 0, 0, 100, 60,
                "9M", 0, 200, 100, true);
        shortRead.setReadString("ACGTACGTA");
        shortRead.setBaseQualityString("IIIIIIIII");
        out.add(shortRead);

        // Two bases: shorter than the indel context size, so getReadContextAtEachPosition returns
        // early with fewer keys than the read has positions.
        final SAMRecord twoBases = ReadFilterDump.read(header, "cov_two_bases", 0, 0, 100, 60,
                "2M", 0, 200, 100, true);
        twoBases.setReadString("AC");
        twoBases.setBaseQualityString("II");
        out.add(twoBases);

        // Every base below the low-quality tail, so the clip leaves nothing and the context
        // covariate zeroes the whole matrix.
        final SAMRecord allLowQual = ReadFilterDump.read(header, "cov_all_low_qual", 0, 0, 100, 60,
                "10M", 0, 200, 100, true);
        allLowQual.setReadString("ACGTACGTAC");
        allLowQual.setBaseQualityString("!!!!!!!!!!");
        out.add(allLowQual);

        // A read in a read group that is not in the header, which is the -1 key.
        final SAMRecord unknownGroup = ReadFilterDump.read(header, "cov_unknown_group", 0, 0, 100,
                60, "10M", 0, 200, 100, true);
        unknownGroup.setReadString("ACGTACGTAC");
        unknownGroup.setBaseQualityString("IIIIIIIIII");
        unknownGroup.setAttribute("RG", "rg-not-in-header");
        out.add(unknownGroup);

        return out;
    }
}

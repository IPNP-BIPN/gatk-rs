/*
 * OverhangFixingManager, taken from the reference.
 *
 * The second piece of SplitNCigarReads, and the one that decides what the output reads look like.
 * The tool splits a read at every N and hands the family here; this class holds the families in a
 * priority queue, remembers where the splices were, and soft-clips the overhang a read leaves on
 * the far side of a splice when the bases there disagree with the reference.
 *
 * Nine behaviours this is built to catch.
 *
 *   - overhangingBasesMismatch REFUSES THREE KINDS OF SPAN BEFORE COUNTING ANYTHING: a span below
 *     1, a span above --max-bases-in-overhang, and a span STRICTLY LONGER than half the read's
 *     non-clipped length. `readLength / 2` is integer division and the comparison is `>`, so on ten
 *     bases a span of five is still tested and six is refused, while on NINE bases five is already
 *     refused: the same span, two different answers, one base of read length apart;
 *   - AND IT HAS TWO WAYS TO SAY YES. More than --max-mismatches-in-overhang mismatches returns
 *     early, and a span where at least `(span+1)/2` bases mismatch returns true at the end even
 *     when the first rule was never tripped. One mismatch out of two is a mismatch by the second
 *     rule while the tolerance is one;
 *   - isLeftOverhang AND isRightOverhang ARE THREE STRICT AND NON-STRICT COMPARISONS EACH, so a
 *     read that starts exactly at the splice start is not a left overhang and one that starts one
 *     base later is;
 *   - addSplicePosition RETURNS null FOR A SPLICE IT HAS ALREADY SEEN, so the reference bases are
 *     fetched once per splice and the waiting reads are run against it once;
 *   - A SPLICE ON A NEW CONTIG CLEARS EVERY SPLICE SEEN SO FAR, and the check is against the FIRST
 *     splice in the sorted set rather than the last one added;
 *   - THE QUEUE IS FLUSHED HALFWAY ON PRESSURE AND COMPLETELY ON A NEW CONTIG: `maxRecordsInMemory
 *     / 2` against 0;
 *   - THE FAMILY IS REPAIRED ON THE WAY OUT, not on the way in: NM, MD and NH are cleared from
 *     every read and the SA tags are written, and only when writing is active;
 *   - activateWriting A SECOND TIME IS A GATKException, and the first call flushes the queue and
 *     clears the splices;
 *   - AND THE MATE KEY IS `name@0-or-1@start`, where the digit is `!isFirstOfPair` when the key is
 *     stored and `isFirstOfPair` when it is looked up: the two calls are deliberately opposite,
 *     because a read stores the key its MATE will search for. A first-of-pair read that stored a
 *     key is found by its second-of-pair mate and not by another first-of-pair read, and the start
 *     in the key is the OLD start, so a mate pointing at the clipped read's new position misses.
 *
 * Output:
 *
 *     reference\t<the reference bases of chr1>
 *     mismatch\t<label>\t<true|false>
 *     overhang\t<label>\t<left|right|neither>
 *     key\t<label>\t<the mate key>
 *     splice\t<label>\t<new|null>\t<the splices held afterwards>
 *     queue\t<label>\t<reads waiting>
 *     written\t<label>\t<one written read, as SAM>
 *     mate\t<label>\t<edited|untouched>\t<the read, as SAM>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: OverhangFixingManagerDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.broadinstitute.hellbender.tools.walkers.rnaseq.OverhangFixingManager;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.GATKReadWriter;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class OverhangFixingManagerDump {

    /** 120 bases, non-repeating enough that a mismatch is a mismatch and not an alignment. */
    static final String CHR1 =
            "ACGTACGTACGTTTTTGGGGCCCCAAAAACGTACGTACGTGATTACAGGCTCTAGCATCGATCGATCGATTAGCTAGCTAGCTAACCGGTTACGTAGGCTTACCGGATCGATCGATCGAT";
    static final String CHR2 =
            "TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGG";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("overhang-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# OverhangFixingManagerDump: the overhang clipper SplitNCigarReads writes through");
        System.out.printf("reference\t%s%n", CHR1);

        final Path fasta = writeReference(dir);
        final SAMFileHeader header = header();
        final GenomeLocParser parser = new GenomeLocParser(header.getSequenceDictionary());

        predicates(header, parser, fasta);
        keys();
        lifecycle(header, parser, fasta);
        clipping(header, parser, fasta);
        mates(header, parser, fasta);
    }

    /** The two pure predicates, reached by reflection because both are protected. */
    static void predicates(final SAMFileHeader header, final GenomeLocParser parser, final Path fasta)
            throws Exception {
        final OverhangFixingManager manager = manager(header, parser, fasta, new Collector(), 100, 1, 40, false, false);

        final Method mismatch = OverhangFixingManager.class.getDeclaredMethod(
                "overhangingBasesMismatch", byte[].class, int.class, int.class, byte[].class,
                int.class, int.class);
        mismatch.setAccessible(true);

        final byte[] read = "ACGTACGTAC".getBytes(StandardCharsets.UTF_8);
        final byte[] same = "ACGTACGTAC".getBytes(StandardCharsets.UTF_8);
        final byte[] oneOff = "TCGTACGTAC".getBytes(StandardCharsets.UTF_8);
        final byte[] twoOff = "TTGTACGTAC".getBytes(StandardCharsets.UTF_8);

        // A span of zero and a negative span are refused before anything is compared.
        System.out.printf("mismatch\tspan-0\t%s%n", mismatch.invoke(manager, read, 0, 10, twoOff, 0, 0));
        System.out.printf("mismatch\tspan-negative\t%s%n", mismatch.invoke(manager, read, 0, 10, twoOff, 0, -1));
        // The gate is `spanToTest > readLength / 2`, integer division and a STRICT comparison, so on
        // a ten-base read a span of five is tested and a span of six is refused.
        System.out.printf("mismatch\tspan-6-of-10\t%s%n", mismatch.invoke(manager, read, 0, 10, twoOff, 0, 6));
        System.out.printf("mismatch\tspan-5-of-10\t%s%n", mismatch.invoke(manager, read, 0, 10, twoOff, 0, 5));
        System.out.printf("mismatch\tspan-4-of-10\t%s%n", mismatch.invoke(manager, read, 0, 10, twoOff, 0, 4));
        // An odd read length halves downwards: on nine bases a span of five is already refused.
        System.out.printf("mismatch\tspan-5-of-9\t%s%n", mismatch.invoke(manager, read, 0, 9, twoOff, 0, 5));
        // No mismatch at all, and one mismatch inside the tolerance of one.
        System.out.printf("mismatch\tidentical\t%s%n", mismatch.invoke(manager, read, 0, 10, same, 0, 4));
        System.out.printf("mismatch\tone-of-four\t%s%n", mismatch.invoke(manager, read, 0, 10, oneOff, 0, 4));
        // One mismatch out of two is inside the tolerance and still true, by the half rule.
        System.out.printf("mismatch\tone-of-two\t%s%n", mismatch.invoke(manager, read, 0, 10, oneOff, 0, 2));
        // Two mismatches out of four is above the tolerance of one.
        System.out.printf("mismatch\ttwo-of-four\t%s%n", mismatch.invoke(manager, read, 0, 10, twoOff, 0, 4));
        // A span longer than the tolerated bases, with the same content as a span that passes.
        final OverhangFixingManager narrow =
                manager(header, parser, fasta, new Collector(), 100, 1, 3, false, false);
        System.out.printf("mismatch\tspan-above-max-bases\t%s%n",
                mismatch.invoke(narrow, read, 0, 10, twoOff, 0, 4));
        // And a tolerance of zero, where one mismatch is already too many.
        final OverhangFixingManager strict =
                manager(header, parser, fasta, new Collector(), 100, 0, 40, false, false);
        System.out.printf("mismatch\tzero-tolerance\t%s%n",
                mismatch.invoke(strict, read, 0, 10, oneOff, 0, 4));

        final Method left = OverhangFixingManager.class.getDeclaredMethod(
                "isLeftOverhang", GenomeLoc.class, GenomeLoc.class);
        final Method right = OverhangFixingManager.class.getDeclaredMethod(
                "isRightOverhang", GenomeLoc.class, GenomeLoc.class);
        left.setAccessible(true);
        right.setAccessible(true);

        final GenomeLoc splice = parser.createGenomeLoc("chr1", 50, 60);
        final int[][] reads = {
                // start, stop, and what each one is
                {51, 70},   // starts inside the splice and ends past it: a left overhang
                {50, 70},   // starts exactly at the splice start: not an overhang either way
                {40, 55},   // ends inside the splice and starts before it: a right overhang
                {40, 60},   // ends exactly at the splice stop: not a right overhang
                {40, 70},   // spans the whole splice
                {61, 70},   // entirely after the splice
                {30, 40},   // entirely before it
                {55, 58},   // entirely inside it
        };
        for (final int[] pair : reads) {
            final GenomeLoc loc = parser.createGenomeLoc("chr1", pair[0], pair[1]);
            final boolean isLeft = (Boolean) left.invoke(null, loc, splice);
            final boolean isRight = (Boolean) right.invoke(null, loc, splice);
            System.out.printf("overhang\t%d-%d\t%s%n", pair[0], pair[1],
                    isLeft ? "left" : (isRight ? "right" : "neither"));
        }
    }

    /** makeKey, which is package-private and whose `@` is the point. */
    static void keys() throws Exception {
        final Method makeKey = OverhangFixingManager.class.getDeclaredMethod(
                "makeKey", String.class, boolean.class, int.class);
        makeKey.setAccessible(true);
        System.out.printf("key\tfirst\t%s%n", makeKey.invoke(null, "read1", true, 100));
        System.out.printf("key\tsecond\t%s%n", makeKey.invoke(null, "read1", false, 100));
        System.out.printf("key\tzero-start\t%s%n", makeKey.invoke(null, "read1", true, 0));
    }

    /** Splices: what is new, what is a duplicate, and what a new contig clears. */
    static void lifecycle(final SAMFileHeader header, final GenomeLocParser parser, final Path fasta)
            throws Exception {
        final Collector collector = new Collector();
        final OverhangFixingManager manager =
                manager(header, parser, fasta, collector, 100, 1, 40, false, false);

        splice(manager, "first", "chr1", 50, 60);
        splice(manager, "again", "chr1", 50, 60);
        splice(manager, "second", "chr1", 20, 30);
        // A splice on another contig clears everything seen so far.
        splice(manager, "new-contig", "chr2", 10, 20);
        splice(manager, "back-to-chr1", "chr1", 50, 60);

        // With overhang fixing off, every splice is null and nothing is remembered.
        final OverhangFixingManager off =
                manager(header, parser, fasta, new Collector(), 100, 1, 40, true, false);
        splice(off, "fixing-off", "chr1", 50, 60);

        // Writing can only be activated once.
        manager.activateWriting();
        try {
            manager.activateWriting();
        } catch (final Exception e) {
            System.out.printf("error\tactivate-twice\t%s:%s%n", e.getClass().getName(), e.getMessage());
        }
    }

    /** Reads through the manager: what gets clipped, what gets written, and in what order. */
    static void clipping(final SAMFileHeader header, final GenomeLocParser parser, final Path fasta)
            throws Exception {
        // A left overhang: the read starts inside the splice and runs past its end. Matching the
        // reference across those six bases leaves the read alone; mismatching clips it.
        run(header, parser, fasta, "left-matching", 100, false,
                new String[][] {{"left", "55", "20M", CHR1.substring(54, 74)}},
                new int[][] {{50, 60}});
        run(header, parser, fasta, "left-mismatching", 100, false,
                new String[][] {{"left", "55", "20M", mutate(CHR1.substring(54, 74), 0, 1, 2)}},
                new int[][] {{50, 60}});
        // A right overhang: the read ends inside the splice and starts before it.
        run(header, parser, fasta, "right-mismatching", 100, false,
                new String[][] {{"right", "40", "15M", mutate(CHR1.substring(39, 54), 10, 11, 12)}},
                new int[][] {{50, 60}});
        // The splice arrives after the read, which is the other order the manager has to handle.
        run(header, parser, fasta, "splice-after-read", 100, false,
                new String[][] {{"left", "55", "20M", mutate(CHR1.substring(54, 74), 0, 1, 2)}},
                new int[][] {{50, 60}});
        // Two reads of one family, which are repaired together on the way out.
        run(header, parser, fasta, "family", 100, false,
                new String[][] {
                        {"piece", "45", "10M10S", CHR1.substring(44, 64)},
                        {"piece", "70", "10S10M", CHR1.substring(59, 79)}},
                new int[][] {{50, 60}});
        // A queue smaller than the reads pushed through it, which flushes halfway.
        run(header, parser, fasta, "pressure", 2, false,
                new String[][] {
                        {"a", "10", "10M", CHR1.substring(9, 19)},
                        {"b", "20", "10M", CHR1.substring(19, 29)},
                        {"c", "30", "10M", CHR1.substring(29, 39)},
                        {"d", "40", "10M", CHR1.substring(39, 49)}},
                new int[][] {});
        // A secondary alignment, which is written but not clipped unless asked for.
        run(header, parser, fasta, "secondary", 100, false,
                new String[][] {{"secondary", "55", "20M", mutate(CHR1.substring(54, 74), 0, 1, 2), "256"}},
                new int[][] {{50, 60}});
        run(header, parser, fasta, "secondary-processed", 100, true,
                new String[][] {{"secondary", "55", "20M", mutate(CHR1.substring(54, 74), 0, 1, 2), "256"}},
                new int[][] {{50, 60}});
    }

    /** The mate repair pass: what the first traversal records and the second one applies. */
    static void mates(final SAMFileHeader header, final GenomeLocParser parser, final Path fasta)
            throws Exception {
        final Collector collector = new Collector();
        final OverhangFixingManager manager =
                manager(header, parser, fasta, collector, 100, 1, 40, false, false);

        // First pass: a read whose overhang is clipped, so its mate's information becomes wrong.
        final GATKRead clipped = read(header, "pair", 55, "20M",
                mutate(CHR1.substring(54, 74), 0, 1, 2), 0x1 | 0x40);
        manager.addSplicePosition("chr1", 50, 60);
        manager.addReadGroup(Collections.singletonList(clipped));
        manager.flush();

        // Before writing is active, setPredictedMateInformation does nothing at all.
        System.out.printf("mate\tbefore-activation\t%s\t%s%n",
                manager.setPredictedMateInformation(mate(header, 0x1 | 0x80, 55, true)) ? "edited" : "untouched",
                ReferenceQueryDump.escape(sam(header, mate(header, 0x1 | 0x80, 55, true))));

        manager.activateWriting();

        // The mate of the clipped read, looked up by the key the clipped read stored: its start
        // moves to where the clipped read now begins, and its MC tag becomes the new cigar.
        final GATKRead after = mate(header, 0x1 | 0x80, 55, true);
        System.out.printf("mate\tafter-activation\t%s\t%s%n",
                manager.setPredictedMateInformation(after) ? "edited" : "untouched",
                ReferenceQueryDump.escape(sam(header, after)));

        // A mate with no MC tag: the position is repaired and nothing else is.
        final GATKRead noMc = mate(header, 0x1 | 0x80, 55, false);
        System.out.printf("mate\tno-mc-tag\t%s\t%s%n",
                manager.setPredictedMateInformation(noMc) ? "edited" : "untouched",
                ReferenceQueryDump.escape(sam(header, noMc)));

        // The same mate start with the paired flag clear, which never reaches the lookup.
        final GATKRead unpaired = mate(header, 0, 55, true);
        System.out.printf("mate\tunpaired\t%s\t%s%n",
                manager.setPredictedMateInformation(unpaired) ? "edited" : "untouched",
                ReferenceQueryDump.escape(sam(header, unpaired)));

        // The key carries the mate's own start, so a read pointing at the wrong place misses.
        final GATKRead wrongStart = mate(header, 0x1 | 0x80, 56, true);
        System.out.printf("mate\twrong-mate-start\t%s\t%s%n",
                manager.setPredictedMateInformation(wrongStart) ? "edited" : "untouched",
                ReferenceQueryDump.escape(sam(header, wrongStart)));

        // And the digit in the key is the OTHER end of the pair: a read that claims to be first of
        // pair searches for a key the first-of-pair read stored, which is not this one.
        final GATKRead firstOfPair = mate(header, 0x1 | 0x40, 55, true);
        System.out.printf("mate\tfirst-of-pair\t%s\t%s%n",
                manager.setPredictedMateInformation(firstOfPair) ? "edited" : "untouched",
                ReferenceQueryDump.escape(sam(header, firstOfPair)));
    }

    /** A mate of the clipped read, with its mate fields set directly so the flags stay as given. */
    static GATKRead mate(final SAMFileHeader header, final int flags, final int mateStart,
                         final boolean withMc) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName("pair");
        record.setFlags(flags);
        record.setReferenceName("chr1");
        record.setAlignmentStart(100);
        record.setCigarString("20M");
        record.setReadBases(CHR1.substring(99, 119).getBytes(StandardCharsets.UTF_8));
        final byte[] quals = new byte[20];
        Arrays.fill(quals, (byte) 35);
        record.setBaseQualities(quals);
        record.setMappingQuality(60);
        record.setMateReferenceName("chr1");
        record.setMateAlignmentStart(mateStart);
        record.setAttribute("RG", "rg1");
        if (withMc) {
            record.setAttribute("MC", "20M");
        }
        return new SAMRecordToGATKReadAdapter(record);
    }

    static String sam(final SAMFileHeader header, final GATKRead read) {
        return read.convertToSAMRecord(header).getSAMString().trim();
    }

    /** One manager driven from empty: the splices, then the reads, then the flush. */
    static void run(final SAMFileHeader header, final GenomeLocParser parser, final Path fasta,
                    final String label, final int maxRecords, final boolean processSecondary,
                    final String[][] reads, final int[][] splices) throws Exception {
        final Collector collector = new Collector();
        final OverhangFixingManager manager =
                manager(header, parser, fasta, collector, maxRecords, 1, 40, false, processSecondary);
        manager.activateWriting();

        // "splice-after-read" is the same input with the two calls the other way round.
        final boolean spliceFirst = !label.equals("splice-after-read");
        if (spliceFirst) {
            for (final int[] splice : splices) {
                manager.addSplicePosition("chr1", splice[0], splice[1]);
            }
        }

        final List<GATKRead> family = new ArrayList<>();
        for (final String[] spec : reads) {
            final int flags = spec.length > 4 ? Integer.parseInt(spec[4]) : 0;
            family.add(read(header, spec[0], Integer.parseInt(spec[1]), spec[2], spec[3], flags));
        }
        if (reads.length > 1 && reads[0][0].equals(reads[1][0])) {
            // One family of supplementary pieces, which is what SplitNCigarReads produces.
            manager.addReadGroup(family);
        } else {
            for (final GATKRead read : family) {
                manager.addReadGroup(Collections.singletonList(read));
            }
        }

        if (!spliceFirst) {
            for (final int[] splice : splices) {
                manager.addSplicePosition("chr1", splice[0], splice[1]);
            }
        }

        System.out.printf("queue\t%s\t%d%n", label, collector.written.size());
        manager.flush();
        for (final GATKRead written : collector.written) {
            System.out.printf("written\t%s\t%s%n", label,
                    ReferenceQueryDump.escape(written.convertToSAMRecord(header).getSAMString().trim()));
        }
    }

    static void splice(final OverhangFixingManager manager, final String label, final String contig,
                       final int start, final int stop) {
        final Object added = manager.addSplicePosition(contig, start, stop);
        final List<String> held = new ArrayList<>();
        for (final Object splice : manager.getSplicesForTesting()) {
            held.add(spliceText(splice));
        }
        // `null` covers two different reasons: a splice already seen, and a manager with overhang
        // fixing turned off. The splices held afterwards are what tells them apart.
        System.out.printf("splice\t%s\t%s\t%s%n", label, added == null ? "null" : "new",
                String.join(",", held));
    }

    /** A Splice's location, read off its public `loc` field. */
    static String spliceText(final Object splice) {
        try {
            final java.lang.reflect.Field field = splice.getClass().getDeclaredField("loc");
            field.setAccessible(true);
            final GenomeLoc loc = (GenomeLoc) field.get(splice);
            return loc.getContig() + ":" + loc.getStart() + "-" + loc.getStop();
        } catch (final Exception e) {
            throw new IllegalStateException(e);
        }
    }

    static OverhangFixingManager manager(final SAMFileHeader header, final GenomeLocParser parser,
                                         final Path fasta, final GATKReadWriter writer,
                                         final int maxRecords, final int maxMismatches,
                                         final int maxBases, final boolean doNotFix,
                                         final boolean processSecondary) throws Exception {
        final ReferenceSequenceFile reference = new CachingIndexedFastaSequenceFile(fasta);
        return new OverhangFixingManager(header, writer, parser, reference, maxRecords, maxMismatches,
                maxBases, doNotFix, processSecondary);
    }

    /** The reference bases with the given zero-based offsets changed to a base that is not there. */
    static String mutate(final String bases, final int... offsets) {
        final char[] chars = bases.toCharArray();
        for (final int offset : offsets) {
            chars[offset] = chars[offset] == 'A' ? 'C' : 'A';
        }
        return new String(chars);
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", CHR1.length()),
                new SAMSequenceRecord("chr2", CHR2.length()))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("s1");
        header.addReadGroup(group);
        return header;
    }

    static GATKRead read(final SAMFileHeader header, final String name, final int start,
                         final String cigar, final String bases, final int flags) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setFlags(flags);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setReadBases(bases.getBytes(StandardCharsets.UTF_8));
        final byte[] quals = new byte[bases.length()];
        Arrays.fill(quals, (byte) 35);
        record.setBaseQualities(quals);
        record.setMappingQuality(60);
        record.setAttribute("RG", "rg1");
        // The three tags the family repair clears, so their disappearance is visible.
        record.setAttribute("NM", 1);
        record.setAttribute("MD", "20");
        record.setAttribute("NH", 2);
        return new SAMRecordToGATKReadAdapter(record);
    }

    static Path writeReference(final Path dir) throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        Files.writeString(fasta, ">chr1\n" + CHR1 + "\n>chr2\n" + CHR2 + "\n", StandardCharsets.UTF_8);
        FastaSequenceIndexCreator.create(fasta, true);
        final Path dict = dir.resolve("reference.dict");
        Files.writeString(dict, "@HD\tVN:1.6\tSO:unsorted\n"
                + "@SQ\tSN:chr1\tLN:" + CHR1.length() + "\n"
                + "@SQ\tSN:chr2\tLN:" + CHR2.length() + "\n", StandardCharsets.UTF_8);
        return fasta;
    }

    static void emptyDirectory(final Path dir) throws Exception {
        if (!Files.isDirectory(dir)) {
            return;
        }
        try (final var entries = Files.list(dir)) {
            for (final Path entry : entries.toList()) {
                Files.deleteIfExists(entry);
            }
        }
    }

    /** A writer that keeps what it is given, in the order it is given it. */
    static final class Collector implements GATKReadWriter {
        final List<GATKRead> written = new ArrayList<>();

        @Override
        public void addRead(final GATKRead read) {
            written.add(read);
        }

        @Override
        public void close() {}
    }
}

/*
 * ReadsDataSource's answers, taken from the reference.
 *
 * Every walker reads its reads through ReadsPathDataSource: `-L` becomes an interval query against
 * the .bai, and no `-L` becomes a full traversal. Which records come back is decided in three
 * places, and none of them is "does the read overlap the interval":
 *
 *   - IntervalUtils.convertSimpleIntervalToQueryInterval resolves the contig against the *reads*
 *     dictionary and throws for a contig it does not hold;
 *   - QueryInterval.optimizeIntervals sorts and merges overlapping *and abutting* intervals, so
 *     two adjacent -L arguments return a read spanning the join once rather than twice;
 *   - BAMQueryMultipleIntervalsIteratorFilter is a stateful single-pass filter: its interval index
 *     only advances, unmapped reads carrying their mate's position are special-cased to end=start,
 *     and once every interval is behind the current record the traversal stops.
 *
 * The fixture travels in the golden as base64, index included, so the port queries exactly the
 * bytes the reference queried rather than a reconstruction of them. It is built to make the
 * decisions above visible: reads that abut, reads spanning a linear-index boundary, a read whose
 * deletion makes it span further than its bases, an unmapped read at its mate's position, a mapped
 * read with an empty cigar, a contig with no reads at all, and unplaced unmapped reads at the end.
 *
 * Output:
 *
 *     bam\t<base64 of the BAM>
 *     bai\t<base64 of the .bai>
 *     query\t<label>\t<record>\n<record>...   (\n literal, one segment per returned read)
 *     count\t<label>\t<n>                     (or E if the reference threw)
 *
 * Usage: ReadsDataSourceDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.Iterator;
import java.util.List;

public class ReadsDataSourceDump {

    static final int CONTIG_LENGTH = 100_000;

    public static void main(final String[] args) throws Exception {
        final Path dir = Files.createTempDirectory("readsquery");
        final Path bam = dir.resolve("reads.bam");
        buildFixture(bam.toFile());
        final Path bai = dir.resolve("reads.bai");

        System.out.println("# ReadsDataSourceDump: which records an interval query returns");
        System.out.printf("bam\t%s%n", Base64.getEncoder().encodeToString(Files.readAllBytes(bam)));
        System.out.printf("bai\t%s%n", Base64.getEncoder().encodeToString(Files.readAllBytes(bai)));

        for (final String[] query : QUERIES) {
            final String label = query[0];
            final List<String> intervals = Arrays.asList(query).subList(1, query.length);
            emit(bam, label, intervals, false);
        }

        // queryUnmapped: the unplaced reads at the tail, reached through the last linear bin.
        emit(bam, "unmapped", null, true);
        // A full traversal bounded by intervals *and* unmapped, which is what -L plus
        // --interval-set-rule leaves the engine doing: intervals first, then the unplaced reads.
        emit(bam, "traverse:chr1:1-100000+unmapped", Arrays.asList("chr1:1-100000"), true);
        // A full traversal with no bounds at all.
        emit(bam, "traverse:all", null, false);
    }

    /**
     * The queries. Each row is a label followed by its intervals; several intervals in one row is
     * one query with several -L arguments, which is where optimizeIntervals becomes visible.
     */
    static final String[][] QUERIES = {
            {"chr1:100-100", "chr1:100-100"},
            {"chr1:105-105", "chr1:105-105"},
            {"chr1:99-99", "chr1:99-99"},
            {"chr1:110-110", "chr1:110-110"},
            {"chr1:150-160", "chr1:150-160"},
            {"chr1:200-200", "chr1:200-200"},
            {"chr1:1-100000", "chr1:1-100000"},
            // Abutting: 100-200 and 201-300 merge into one interval, so the read crossing the
            // join is returned once.
            {"abut", "chr1:100-200", "chr1:201-300"},
            // The same interval twice.
            {"duplicate", "chr1:100-200", "chr1:100-200"},
            // Overlapping, merged.
            {"overlap", "chr1:100-160", "chr1:150-250"},
            // Out of order and across contigs: optimizeIntervals sorts by reference index.
            {"unsorted", "chr2:1-1000", "chr1:100-200"},
            // An unmapped read placed at its mate's coordinate.
            {"chr1:300-300", "chr1:300-300"},
            {"chr1:305-305", "chr1:305-305"},
            {"chr1:315-315", "chr1:315-315"},
            // A mapped read with an empty cigar: its alignment end is its start minus one.
            {"chr1:499-499", "chr1:499-499"},
            {"chr1:500-500", "chr1:500-500"},
            {"chr1:495-505", "chr1:495-505"},
            // The 16 kb linear-index boundary.
            {"chr1:16384-16384", "chr1:16384-16384"},
            {"chr1:16380-16395", "chr1:16380-16395"},
            // The end of the contig, and past it.
            {"chr1:99995-100000", "chr1:99995-100000"},
            {"chr1:100001-100010", "chr1:100001-100010"},
            {"chr2:1-1000", "chr2:1-1000"},
            {"chr2:5000-5010", "chr2:5000-5010"},
            // A contig the BAM declares but no read touches.
            {"chr3:1-1000", "chr3:1-1000"},
            // A contig the BAM does not declare at all.
            {"chrX:1-100", "chrX:1-100"},
            // Four intervals, unsorted, overlapping and abutting at once.
            {"messy", "chr1:250-260", "chr2:1-100000", "chr1:100-200", "chr1:100-200"},
    };

    static void emit(final Path bam, final String label, final List<String> intervals,
                     final boolean unmapped) {
        final List<String> rows = new ArrayList<>();
        String count;
        try (final ReadsDataSource source = new ReadsPathDataSource(bam)) {
            final Iterator<GATKRead> iterator;
            if (label.startsWith("traverse:")) {
                source.setTraversalBounds(toIntervals(intervals), unmapped);
                iterator = source.iterator();
            } else if (intervals == null) {
                iterator = source.queryUnmapped();
            } else {
                // A query takes one interval at a time; several -L arguments reach the same
                // filter through setTraversalBounds, so the multi-interval rows go that way.
                if (intervals.size() == 1) {
                    iterator = source.query(new SimpleInterval(intervals.get(0)));
                } else {
                    source.setTraversalBounds(toIntervals(intervals), false);
                    iterator = source.iterator();
                }
            }
            while (iterator.hasNext()) {
                rows.add(describe(iterator.next()));
            }
            count = String.valueOf(rows.size());
        } catch (final Exception | AssertionError e) {
            rows.clear();
            rows.add("E");
            count = "E";
        }
        System.out.printf("query\t%s\t%s%n", label, String.join("\\n", rows));
        System.out.printf("count\t%s\t%s%n", label, count);
    }

    static List<SimpleInterval> toIntervals(final List<String> raw) {
        if (raw == null) {
            // setTraversalBounds(null, false) is an unbounded traversal, not an empty one.
            return null;
        }
        final List<SimpleInterval> out = new ArrayList<>();
        for (final String text : raw) {
            out.add(new SimpleInterval(text));
        }
        return out;
    }

    /** Enough of the read to see which record it is and why it matched. */
    static String describe(final GATKRead read) {
        return String.format("%s|%d|%s|%d|%s|%d",
                read.getName(),
                read.isUnmapped() ? -1 : read.getStart(),
                read.getCigar().toString(),
                read.getFlags(),
                read.getAssignedContig() == null ? "*" : read.getAssignedContig(),
                read.getAssignedStart());
    }

    static void buildFixture(final File bam) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        for (final String contig : new String[] {"chr1", "chr2", "chr3"}) {
            dictionary.addSequence(new SAMSequenceRecord(contig, CONTIG_LENGTH));
        }
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        final List<SAMRecord> records = new ArrayList<>();
        records.add(mapped(header, "r001", "chr1", 100, "10M"));
        records.add(mapped(header, "r002", "chr1", 150, "10M"));
        // Spans the 200/201 join that the abutting query merges over.
        records.add(mapped(header, "r003", "chr1", 195, "10M"));
        records.add(mapped(header, "r004", "chr1", 200, "10M"));
        records.add(mapped(header, "r005", "chr1", 250, "10M"));
        // A deletion makes the reference span longer than the bases: 300-319 from 10 bases.
        records.add(mapped(header, "r006", "chr1", 300, "5M10D5M"));
        // Unmapped, but carrying its mate's coordinate. Its alignment end is 0, so the filter
        // has to special-case it or it is invisible to every query.
        records.add(unmappedAtMate(header, "u001", "chr1", 300));
        // Mapped with no cigar at all: alignment end is start - 1.
        records.add(mapped(header, "m001", "chr1", 500, "*"));
        // Across the first 16 kb linear-index boundary.
        records.add(mapped(header, "r007", "chr1", 16_380, "10M"));
        records.add(mapped(header, "r008", "chr1", 16_390, "10M"));
        // Ends exactly at the contig's last base.
        records.add(mapped(header, "r009", "chr1", 99_995, "6M"));
        records.add(mapped(header, "r101", "chr2", 100, "10M"));
        records.add(mapped(header, "r102", "chr2", 5_000, "10M"));
        // chr3 deliberately has no reads: its linear index is empty, which is what
        // getStartOfLastLinearBin walks past.
        records.add(unplaced(header, "x001"));
        records.add(unplaced(header, "x002"));
        records.add(unplaced(header, "x003"));

        final SAMFileWriterFactory factory = new SAMFileWriterFactory().setCreateIndex(true);
        try (final SAMFileWriter writer = factory.makeBAMWriter(header, true, bam)) {
            for (final SAMRecord record : records) {
                writer.addAlignment(record);
            }
        }
    }

    static SAMRecord mapped(final SAMFileHeader header, final String name, final String contig,
                            final int start, final String cigar) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName(contig);
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setMappingQuality(60);
        record.setReadBases("ACGTACGTAC".getBytes());
        record.setBaseQualities(new byte[] {30, 30, 30, 30, 30, 30, 30, 30, 30, 30});
        return record;
    }

    static SAMRecord unmappedAtMate(final SAMFileHeader header, final String name,
                                    final String contig, final int start) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReadPairedFlag(true);
        record.setFirstOfPairFlag(true);
        record.setReadUnmappedFlag(true);
        record.setReferenceName(contig);
        record.setAlignmentStart(start);
        record.setMateReferenceName(contig);
        record.setMateAlignmentStart(start);
        record.setMappingQuality(0);
        record.setReadBases("ACGTACGTAC".getBytes());
        record.setBaseQualities(new byte[] {30, 30, 30, 30, 30, 30, 30, 30, 30, 30});
        return record;
    }

    static SAMRecord unplaced(final SAMFileHeader header, final String name) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReadUnmappedFlag(true);
        record.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
        record.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
        record.setMappingQuality(0);
        record.setReadBases("ACGTACGTAC".getBytes());
        record.setBaseQualities(new byte[] {30, 30, 30, 30, 30, 30, 30, 30, 30, 30});
        return record;
    }
}

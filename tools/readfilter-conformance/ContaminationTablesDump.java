/*
 * The two contamination output tables, taken from the reference.
 *
 * `CalculateContamination` writes both: a contamination table, which is the answer, and a segments
 * table of minor allele fractions. They are built from the same tsv layer as PileupSummary and yet
 * differ in what they carry, which is what this measures.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE CONTAMINATION TABLE CARRIES NO SAMPLE METADATA AT ALL, so its first line is the header,
 *     while the segments table writes `#<METADATA>SAMPLE=` before its own. The sample of the
 *     contamination table is a COLUMN instead, which is why gathering one is not the other;
 *   - BOTH DOUBLE COLUMNS ARE Double.toString's SPELLING, so an estimate of exactly zero is `0.0`
 *     and a tiny error comes out in exponent form;
 *   - THE SEGMENTS TABLE'S START AND END GO THROUGH SimpleInterval, which validates: a row whose
 *     end is before its start is refused when the record is BUILT, not when the file is parsed, so
 *     the refusal is an IllegalArgumentException from the interval and not the table's BadInput;
 *   - AND A START OF ZERO IS REFUSED THE SAME WAY, because the interval is one-based;
 *   - A ROW READ BACK KEEPS THE INTERVAL'S OWN toString, `contig:start-end`, which is not what the
 *     columns hold;
 *   - THE CONTAMINATION READER TAKES A FILE WITH SAMPLE METADATA WITHOUT COMPLAINT and ignores it,
 *     since the metadata map is filled whether or not anything asks for it;
 *   - AND A MISSING COLUMN IS RAISED WHEN IT IS ASKED FOR, on the first record, so a file with the
 *     right columns in the wrong ORDER is read without complaint while one missing a column is an
 *     IllegalArgumentException naming it.
 *
 * Output:
 *
 *     written\t<label>\t<the whole file, escaped>
 *     read\t<label>\t<sample>\t<row>\t<the fields, comma separated>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: ContaminationTablesDump
 */

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.tools.walkers.contamination.ContaminationRecord;
import org.broadinstitute.hellbender.tools.walkers.contamination.MinorAlleleFractionRecord;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class ContaminationTablesDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("contaminationtables-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# ContaminationTablesDump: the two contamination output tables, from the reference");

        // The contamination table, whose sample is a column and not metadata.
        final List<ContaminationRecord> contamination = new ArrayList<>();
        contamination.add(new ContaminationRecord("sample-a", 0.0, 0.0));
        contamination.add(new ContaminationRecord("sample-a", 0.05, 0.001));
        contamination.add(new ContaminationRecord("sample-a", 1.0 / 3, 1e-7));
        contamination.add(new ContaminationRecord("has\ta tab", 0.5, Double.NaN));
        writeContamination(dir, "contamination", contamination);
        readContamination(dir, "contamination", dir.resolve("contamination.table"));

        // And one with no records at all, which still has its header.
        writeContamination(dir, "contamination-empty", List.of());
        readContamination(dir, "contamination-empty", dir.resolve("contamination-empty.table"));

        // The segments table, whose sample IS metadata.
        final List<MinorAlleleFractionRecord> segments = new ArrayList<>();
        segments.add(new MinorAlleleFractionRecord(new SimpleInterval("chr1", 1, 1000), 0.5));
        segments.add(new MinorAlleleFractionRecord(new SimpleInterval("chr1", 1001, 2000), 0.0));
        segments.add(new MinorAlleleFractionRecord(new SimpleInterval("chr2", 5, 5), 1.0 / 3));
        writeSegments(dir, "segments", "sample-a", segments);
        readSegments(dir, "segments", dir.resolve("segments.table"));

        writeSegments(dir, "segments-empty", "sample-a", List.of());
        readSegments(dir, "segments-empty", dir.resolve("segments-empty.table"));

        // A contamination file carrying metadata nobody asked for.
        handWritten(dir, "contamination-with-metadata",
                "#<METADATA>SAMPLE=sample-a\nsample\tcontamination\terror\nsample-a\t0.1\t0.01\n");
        readContamination(dir, "contamination-with-metadata",
                dir.resolve("contamination-with-metadata.table"));

        // The same columns in another order, which the reader takes by name.
        handWritten(dir, "contamination-reordered",
                "error\tsample\tcontamination\n0.01\tsample-a\t0.1\n");
        readContamination(dir, "contamination-reordered",
                dir.resolve("contamination-reordered.table"));

        // One missing a column the reader asks for.
        handWritten(dir, "contamination-short",
                "sample\tcontamination\nsample-a\t0.1\n");
        readContamination(dir, "contamination-short", dir.resolve("contamination-short.table"));

        // A segment whose end is before its start, and one whose start is zero: the interval
        // validates, not the table.
        handWritten(dir, "segments-backwards",
                "#<METADATA>SAMPLE=sample-a\ncontig\tstart\tend\tminor_allele_fraction\nchr1\t100\t50\t0.5\n");
        readSegments(dir, "segments-backwards", dir.resolve("segments-backwards.table"));
        handWritten(dir, "segments-zero-start",
                "#<METADATA>SAMPLE=sample-a\ncontig\tstart\tend\tminor_allele_fraction\nchr1\t0\t50\t0.5\n");
        readSegments(dir, "segments-zero-start", dir.resolve("segments-zero-start.table"));

        // And a segments file with no sample metadata, which the reader does not mind.
        handWritten(dir, "segments-nameless",
                "contig\tstart\tend\tminor_allele_fraction\nchr1\t1\t10\t0.25\n");
        readSegments(dir, "segments-nameless", dir.resolve("segments-nameless.table"));
    }

    static void writeContamination(final Path dir, final String label,
                                   final List<ContaminationRecord> records) throws Exception {
        final Path file = dir.resolve(label + ".table");
        ContaminationRecord.writeToFile(records, file.toFile());
        System.out.printf("written\t%s\t%s%n", label,
                ReferenceQueryDump.escape(Files.readString(file)));
    }

    static void writeSegments(final Path dir, final String label, final String sample,
                              final List<MinorAlleleFractionRecord> records) throws Exception {
        final Path file = dir.resolve(label + ".table");
        MinorAlleleFractionRecord.writeToFile(sample, records, file.toFile());
        System.out.printf("written\t%s\t%s%n", label,
                ReferenceQueryDump.escape(Files.readString(file)));
    }

    static void handWritten(final Path dir, final String label, final String text) throws Exception {
        final Path file = dir.resolve(label + ".table");
        Files.writeString(file, text, StandardCharsets.UTF_8);
        System.out.printf("written\t%s\t%s%n", label, ReferenceQueryDump.escape(text));
    }

    /** The contamination table has no sample of its own, so the sample column is printed instead. */
    static void readContamination(final Path dir, final String label, final Path file) {
        final List<ContaminationRecord> records;
        try {
            records = ContaminationRecord.readFromFile(file.toFile());
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        if (records.isEmpty()) {
            System.out.printf("read\t%s\t\t\t%n", label);
            return;
        }
        for (int i = 0; i < records.size(); i++) {
            final ContaminationRecord record = records.get(i);
            System.out.printf("read\t%s\t\t%d\t%s,%s,%s%n", label, i,
                    ReferenceQueryDump.escape(record.getSample()),
                    String.valueOf(record.getContamination()),
                    String.valueOf(record.getError()));
        }
    }

    /** The segments table carries a sample, and each record keeps the interval's own toString. */
    static void readSegments(final Path dir, final String label, final Path file) {
        final ImmutablePair<String, List<MinorAlleleFractionRecord>> pair;
        try {
            pair = MinorAlleleFractionRecord.readFromFile(file.toFile());
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        final List<MinorAlleleFractionRecord> records = pair.getRight();
        if (records.isEmpty()) {
            System.out.printf("read\t%s\t%s\t\t%n", label, String.valueOf(pair.getLeft()));
            return;
        }
        for (int i = 0; i < records.size(); i++) {
            final MinorAlleleFractionRecord record = records.get(i);
            System.out.printf("read\t%s\t%s\t%d\t%s,%s,%d,%d,%s%n", label,
                    String.valueOf(pair.getLeft()), i,
                    record.getSegment().toString(),
                    record.getContig(), record.getStart(), record.getEnd(),
                    String.valueOf(record.getMinorAlleleFraction()));
        }
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
}

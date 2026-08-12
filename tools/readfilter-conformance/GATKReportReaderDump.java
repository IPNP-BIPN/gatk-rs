/*
 * Reading a GATKReport, taken from the reference.
 *
 * The writing side is already measured (the gatk-report suite). This is the other direction, and it
 * is the one ApplyBQSR uses: a recalibration report is parsed back into typed values, and which Java
 * type each cell comes back as decides which branch RecalibrationReport's asLong, asDouble and
 * decodeByte take.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE COLUMNS ARE CUT AT THE POSITIONS OF THE HEADER LINE'S WORDS. getWordStarts walks the
 *     column-name line and records every index whose previous character is whitespace and whose own
 *     is not; splitFixedWidth then cuts EVERY data line at exactly those indexes. The widths in the
 *     file are not declared anywhere: they are inferred from where the header's words begin;
 *   - SO A VALUE CONTAINING A SPACE IS SPLIT BY ITS OWN CONTENT if it pushes past the next column's
 *     start, and a value that is merely wider than its header is silently merged with its
 *     neighbour. The writer never produces either, because it pads every column to its widest value,
 *     but a reader given a hand-edited file will;
 *   - EVERY FIELD IS TRIMMED after the cut, so the left-aligned trailing padding disappears, and a
 *     data line SHORTER than the header's last column start is a StringIndexOutOfBoundsException
 *     rather than a short row, because substring is called with the header's indexes and no clamp;
 *   - A `%d` COLUMN PARSES TO Long AND NOT Integer, which is why RecalibrationReport's asLong has an
 *     Integer branch it never takes from a parsed table and decodeByte has a Long branch it always
 *     takes;
 *   - A `%s` COLUMN PARSES TO String, SO `null` COMES BACK AS THE FOUR CHARACTERS and not as a null,
 *     and so does `NaN` in a `%s` column while `NaN` in a `%.2f` column comes back as a Double;
 *   - THE ROW IDS OF A PARSED TABLE ARE ITS INDEXES, and its sorting is DO_NOT_SORT whatever the
 *     table was written with, which is what makes a parse and a second writing reproduce the first;
 *   - AND getReadGroups READS THE RecalTable0 TABLE and returns a SORTED SET, so the read group order
 *     of a parsed report is alphabetical rather than the file's.
 *
 * Output:
 *
 *     report\t<label>\t<version>\t<numTables>
 *     table\t<label>\t<index>\t<name>\t<description>\t<numRows>\t<numColumns>
 *     column\t<label>\t<table>\t<index>\t<name>\t<format>\t<dataType>
 *     cell\t<label>\t<table>\t<row>\t<column>\t<javaClass>\t<value>
 *     starts\t<label>\t<comma separated column start indexes>
 *     split\t<label>\t<field>|<field>|...
 *     readgroups\t<label>\t<comma separated, in the order returned>
 *     roundtrip\t<label>\t<true|false>
 *     error\t<what>\t<exception>\t<message>
 *
 * Usage: GATKReportReaderDump
 */

import org.broadinstitute.hellbender.utils.report.GATKReport;
import org.broadinstitute.hellbender.utils.report.GATKReportTable;
import org.broadinstitute.hellbender.utils.text.TextFormattingUtils;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class GATKReportReaderDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# GATKReportReaderDump: reading a GATKReport");

        splitting();

        for (final Map.Entry<String, String> entry : reports().entrySet()) {
            read(entry.getKey(), entry.getValue());
        }

        readGroups();
        errors();
    }

    /**
     * The fixed-width split itself, which is where a parsed table's columns come from.
     */
    static void splitting() {
        final String[] lines = {
                "Alpha  Beta  Gamma",
                "  Alpha  Beta",
                "Alpha",
                "",
                "   ",
                "A B C",
                "A  B  C  ",
        };
        for (final String line : lines) {
            final List<Integer> starts = TextFormattingUtils.getWordStarts(line);
            System.out.printf("starts\t%s\t%s%n", line.replace(' ', '_'), join(starts));
            System.out.printf("split\t%s\t%s%n", line.replace(' ', '_'), split(line, starts));
        }

        // The header's word starts applied to a data line whose values are wider than the header's,
        // which is what a hand-edited file looks like.
        final String header = "Alpha  Beta  Gamma";
        final List<Integer> starts = TextFormattingUtils.getWordStarts(header);
        for (final String data : new String[] {
                "1      2     3",
                "1234567  2  3",
                "a b    c d   e f",
                "1",
        }) {
            System.out.printf("split\theader:%s\t%s%n", data.replace(' ', '_'),
                    split(data, starts));
        }
    }

    /** The reports to parse, each shaped for a different type or layout. */
    static Map<String, String> reports() {
        final Map<String, String> out = new LinkedHashMap<>();

        // One column per data type, written by the reference so the layout is the reference's.
        final GATKReport types = new GATKReport();
        types.addTable("Types", "one column per data type", 6, GATKReportTable.Sorting.DO_NOT_SORT);
        final GATKReportTable table = types.getTable("Types");
        table.addColumn("Name", "%s");
        table.addColumn("Count", "%d");
        table.addColumn("Rate", "%.4f");
        table.addColumn("Flag", "%b");
        table.addColumn("Letter", "%c");
        table.addColumn("Untyped", "");
        table.set(0, "Name", "short");
        table.set(0, "Count", 1);
        table.set(0, "Rate", 0.5);
        table.set(0, "Flag", true);
        table.set(0, "Letter", 'A');
        table.set(0, "Untyped", 0.5);
        table.set(1, "Name", "null");
        table.set(1, "Count", 1234567);
        table.set(1, "Rate", Double.NaN);
        table.set(1, "Flag", false);
        table.set(1, "Letter", 'z');
        table.set(1, "Untyped", "text");
        out.put("types", render(types));

        // Two tables, which is what a recalibration report is more than one of.
        final GATKReport two = new GATKReport();
        two.addTable("First", "the first table", 1, GATKReportTable.Sorting.DO_NOT_SORT);
        two.getTable("First").addColumn("Argument", "%s");
        two.getTable("First").set(0, "Argument", "value");
        two.addTable("Second", "the second table", 2, GATKReportTable.Sorting.DO_NOT_SORT);
        two.getTable("Second").addColumn("Key", "%s");
        two.getTable("Second").addColumn("Value", "%d");
        two.getTable("Second").set(0, "Key", "k");
        two.getTable("Second").set(0, "Value", 7);
        out.put("two", render(two));

        // A table written SORT_BY_COLUMN, to show a parsed one comes back DO_NOT_SORT.
        final GATKReport sorted = new GATKReport();
        sorted.addTable("Sorted", "rows out of order", 2, GATKReportTable.Sorting.SORT_BY_COLUMN);
        sorted.getTable("Sorted").addColumn("RowKey", "%s");
        sorted.getTable("Sorted").addColumn("Value", "%d");
        for (final String key : new String[] {"bbb", "aaa", "ccc"}) {
            sorted.getTable("Sorted").set(key, "RowKey", key);
            sorted.getTable("Sorted").set(key, "Value", key.length());
        }
        out.put("sorted", render(sorted));

        // A table with no description, which the header may simply not carry.
        out.put("no-description",
                "#:GATKReport.v1.1:1\n"
                        + "#:GATKTable:2:1:%s:%d:;\n"
                        + "#:GATKTable:Bare\n"
                        + "Key  Value\n"
                        + "k        7\n"
                        + "\n");

        // A table with no rows at all.
        out.put("no-rows",
                "#:GATKReport.v1.1:1\n"
                        + "#:GATKTable:2:0:%s:%d:;\n"
                        + "#:GATKTable:Empty:nothing in it\n"
                        + "Key  Value\n"
                        + "\n");

        // A hand-written file whose values are wider than its header's columns, which the writer
        // never produces and the reader will happily cut in the wrong places.
        out.put("ragged",
                "#:GATKReport.v1.1:1\n"
                        + "#:GATKTable:2:2:%s:%d:;\n"
                        + "#:GATKTable:Ragged:values wider than the header\n"
                        + "Key  Value\n"
                        + "kkkkkk  7\n"
                        + "k  8\n"
                        + "\n");

        return out;
    }

    /** One report: its tables, its columns' recovered types, and every cell's Java class. */
    static void read(final String label, final String text) {
        final GATKReport report;
        try {
            report = new GATKReport(new ByteArrayInputStream(text.getBytes(StandardCharsets.UTF_8)));
        } catch (final Exception e) {
            System.out.printf("error\tread@%s\t%s\t%s%n", label, e.getClass().getSimpleName(),
                    e.getMessage());
            return;
        }
        System.out.printf("report\t%s\t%s\t%d%n", label, report.getVersion(),
                report.getTables().size());

        int tableIndex = 0;
        for (final GATKReportTable table : report.getTables()) {
            System.out.printf("table\t%s\t%d\t%s\t%s\t%d\t%d%n", label, tableIndex,
                    table.getTableName(), table.getTableDescription(), table.getNumRows(),
                    table.getNumColumns());
            int columnIndex = 0;
            for (final var column : table.getColumnInfo()) {
                System.out.printf("column\t%s\t%s\t%d\t%s\t%s\t%s%n", label, table.getTableName(),
                        columnIndex, column.getColumnName(), column.getFormat(),
                        column.getDataType());
                columnIndex++;
            }
            for (int row = 0; row < table.getNumRows(); row++) {
                for (final var column : table.getColumnInfo()) {
                    final Object value = table.get(row, column.getColumnName());
                    System.out.printf("cell\t%s\t%s\t%d\t%s\t%s\t%s%n", label, table.getTableName(),
                            row, column.getColumnName(),
                            value == null ? "null" : value.getClass().getSimpleName(),
                            String.valueOf(value));
                }
            }
            tableIndex++;
        }

        // A parse and a second writing, which is what makes a gathered report reproducible.
        System.out.printf("roundtrip\t%s\t%b%n", label, render(report).equals(text));
    }

    /** getReadGroups, which reads one named table and sorts what it finds. */
    static void readGroups() {
        // The read group table of a recalibration report, deliberately out of alphabetical order
        // and with one identifier repeated.
        final String text = "#:GATKReport.v1.1:1\n"
                + "#:GATKTable:3:4:%s:%s:%d:;\n"
                + "#:GATKTable:RecalTable0:\n"
                + "ReadGroup  EventType  EmpiricalQuality\n"
                + "zebra      M                        30\n"
                + "alpha      M                        29\n"
                + "alpha      I                        45\n"
                + "middle     D                        45\n"
                + "\n";
        final GATKReport report =
                new GATKReport(new ByteArrayInputStream(text.getBytes(StandardCharsets.UTF_8)));
        System.out.printf("readgroups\trecal\t%s%n", String.join(",", report.getReadGroups()));
        System.out.printf("table\trecal\t0\t%s\t%s\t%d\t%d%n",
                report.getTable("RecalTable0").getTableName(),
                report.getTable("RecalTable0").getTableDescription(),
                report.getTable("RecalTable0").getNumRows(),
                report.getTable("RecalTable0").getNumColumns());
    }

    /** Every input the reader refuses. */
    static void errors() {
        attempt("empty-stream", "");
        attempt("legacy-version", "#:GATKReport.v0.1:1\n");
        attempt("no-such-version", "#:GATKReport.v9.9:1\n");
        attempt("not-a-report", "hello\n");
        // A header promising more tables than the file holds.
        attempt("too-few-tables",
                "#:GATKReport.v1.1:2\n"
                        + "#:GATKTable:1:1:%s:;\n"
                        + "#:GATKTable:One:\n"
                        + "Key\n"
                        + "k\n"
                        + "\n");
        // A header promising more rows than the file holds.
        attempt("too-few-rows",
                "#:GATKReport.v1.1:1\n"
                        + "#:GATKTable:1:3:%s:;\n"
                        + "#:GATKTable:One:\n"
                        + "Key\n"
                        + "k\n"
                        + "\n");
        // A `%d` column holding something that is not a number.
        attempt("unparseable-integer",
                "#:GATKReport.v1.1:1\n"
                        + "#:GATKTable:1:1:%d:;\n"
                        + "#:GATKTable:One:\n"
                        + "Value\n"
                        + "abc\n"
                        + "\n");
        // And a table asked for by a name the report does not carry.
        try {
            final String text = "#:GATKReport.v1.1:1\n"
                    + "#:GATKTable:1:1:%s:;\n"
                    + "#:GATKTable:One:\n"
                    + "Key\n"
                    + "k\n"
                    + "\n";
            final GATKReport report =
                    new GATKReport(new ByteArrayInputStream(text.getBytes(StandardCharsets.UTF_8)));
            System.out.printf("error\tunknown-table\tnone\t%s%n",
                    String.valueOf(report.getTable("Nonesuch")));
        } catch (final Exception e) {
            System.out.printf("error\tunknown-table\t%s\t%s%n", e.getClass().getSimpleName(),
                    e.getMessage());
        }
    }

    static void attempt(final String what, final String text) {
        try {
            final GATKReport report =
                    new GATKReport(new ByteArrayInputStream(text.getBytes(StandardCharsets.UTF_8)));
            System.out.printf("error\t%s\tnone\t%d tables%n", what, report.getTables().size());
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s\t%s%n", what, e.getClass().getSimpleName(),
                    e.getMessage());
        }
    }

    /**
     * One split, or the exception it is. A data line SHORTER than the header's last column start is
     * a StringIndexOutOfBoundsException rather than a short row: `substring` is called with the
     * header's indexes and no clamping.
     */
    static String split(final String line, final List<Integer> starts) {
        try {
            return String.join("|", TextFormattingUtils.splitFixedWidth(line, starts));
        } catch (final Exception e) {
            return "E:" + e.getClass().getSimpleName() + ":" + e.getMessage();
        }
    }

    static String join(final List<Integer> values) {
        final StringBuilder out = new StringBuilder();
        for (final int value : values) {
            if (out.length() != 0) {
                out.append(',');
            }
            out.append(value);
        }
        return out.toString();
    }

    static String render(final GATKReport report) {
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        try (final PrintStream out = new PrintStream(bytes, true, StandardCharsets.UTF_8)) {
            report.print(out);
        }
        return bytes.toString(StandardCharsets.UTF_8);
    }
}

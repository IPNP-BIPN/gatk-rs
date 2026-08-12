/*
 * GATKReport, taken from the reference.
 *
 * The file format BQSR is written in, measured before either tool that uses it. ApplyBQSR reads a
 * recalibration table and BaseRecalibrator writes one, and both are this format, so the bytes are
 * settled here or they are settled twice.
 *
 * Seven behaviours this is built to catch, and every one of them is about a byte in a text file.
 *
 *   - THE COLUMN WIDTH IS THE WIDEST THING IN THE COLUMN, name included, and it is COMPUTED FROM
 *     THE FORMATTED VALUE rather than from the value. A double written `%.4f` is as wide as its
 *     rendering, not as wide as its `toString`, and the width is frozen the first time
 *     `getColumnFormat()` is called;
 *   - COLUMNS ARE SEPARATED BY EXACTLY TWO SPACES, and the padding is inside the column's own
 *     format rather than between them, so the last column of every row carries its trailing
 *     padding when it is left-aligned and carries none when it is right-aligned;
 *   - ALIGNMENT IS RIGHT UNTIL A VALUE ASKS FOR LEFT. The default is RIGHT, and a value that is
 *     not numeric and is not one of `null`, `NA`, `Infinity`, `-Infinity`, `NaN` turns the whole
 *     column left-aligned;
 *   - A COLUMN WITH AN EMPTY FORMAT IS Unknown, AND A DOUBLE IN IT GETS `%.8f` rather than the
 *     column's `%s`. That is a special case in `writeRow` and nothing in the table declaration
 *     says so;
 *   - A NON-FINITE DOUBLE ESCAPES ITS OWN FORMAT. `Double.isFinite` decides, and `NaN`,
 *     `Infinity` and `-Infinity` are written through `toString()` instead of through the column's
 *     format, so a `%.4f` column can hold a value with no decimal point at all;
 *   - A NULL IS THE FOUR CHARACTERS `null`, whatever the column's type;
 *   - AND THE TABLE HEADER CARRIES THE COLUMN FORMATS BUT NOT THEIR WIDTHS, so a reader recomputes
 *     every width from the values it parses. The dump round-trips a report through the parser to
 *     show whether the second writing equals the first.
 *
 * Output:
 *
 *     version\t<the #:GATKReport line>
 *     report\t<label>\t<escaped whole report text>
 *     roundtrip\t<label>\t<true|false>
 *     width\t<label>\t<column name>\t<width>\t<alignment>
 *     line\t<label>\t<n>\t<the line, with spaces shown as underscores>
 *
 * Usage: GATKReportDump
 */

import org.broadinstitute.hellbender.utils.report.GATKReport;
import org.broadinstitute.hellbender.utils.report.GATKReportTable;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;

public class GATKReportDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# GATKReportDump: GATKReport");

        // Every data type the format knows, plus the two special cases in writeRow.
        final GATKReport types = new GATKReport();
        types.addTable("Types", "one column per data type", 6, GATKReportTable.Sorting.DO_NOT_SORT);
        final GATKReportTable typesTable = types.getTable("Types");
        typesTable.addColumn("Name", "%s");
        typesTable.addColumn("Count", "%d");
        typesTable.addColumn("Rate", "%.4f");
        typesTable.addColumn("Flag", "%b");
        typesTable.addColumn("Letter", "%c");
        // An empty format is the Unknown type, and a double in it is written %.8f.
        typesTable.addColumn("Untyped", "");

        typesTable.set(0, "Name", "short");
        typesTable.set(0, "Count", 1);
        typesTable.set(0, "Rate", 0.5);
        typesTable.set(0, "Flag", true);
        typesTable.set(0, "Letter", 'A');
        typesTable.set(0, "Untyped", 0.5);

        typesTable.set(1, "Name", "a considerably longer value");
        typesTable.set(1, "Count", 1234567);
        typesTable.set(1, "Rate", 0.123456789);
        typesTable.set(1, "Flag", false);
        typesTable.set(1, "Letter", 'z');
        typesTable.set(1, "Untyped", 1.0 / 3.0);

        // The escapes: a null, and the three non-finite doubles that leave their own format.
        typesTable.set(2, "Name", null);
        typesTable.set(2, "Count", 0);
        typesTable.set(2, "Rate", Double.NaN);
        typesTable.set(2, "Flag", true);
        typesTable.set(2, "Letter", 'x');
        typesTable.set(2, "Untyped", Double.POSITIVE_INFINITY);

        typesTable.set(3, "Name", "neg");
        typesTable.set(3, "Count", -42);
        typesTable.set(3, "Rate", Double.NEGATIVE_INFINITY);
        typesTable.set(3, "Flag", false);
        typesTable.set(3, "Letter", 'y');
        typesTable.set(3, "Untyped", Double.NaN);

        emit(types, "types");

        // A table whose every value is numeric, so nothing ever asks for left alignment.
        final GATKReport numeric = new GATKReport();
        numeric.addTable("Numeric", "right aligned throughout", 2,
                GATKReportTable.Sorting.DO_NOT_SORT);
        final GATKReportTable numericTable = numeric.getTable("Numeric");
        numericTable.addColumn("Quality", "%d");
        numericTable.addColumn("EmpiricalQuality", "%.4f");
        for (int i = 0; i < 3; i++) {
            numericTable.set(i, "Quality", 10 + i * 10);
            numericTable.set(i, "EmpiricalQuality", 10.0 + i * 10 + 0.5);
        }
        emit(numeric, "numeric");

        // Two tables in one report, which is what a recalibration report is.
        final GATKReport two = new GATKReport();
        two.addTable("First", "the first table", 1, GATKReportTable.Sorting.DO_NOT_SORT);
        two.getTable("First").addColumn("Argument", "%s");
        two.getTable("First").set(0, "Argument", "value");
        two.addTable("Second", "the second table", 2, GATKReportTable.Sorting.DO_NOT_SORT);
        two.getTable("Second").addColumn("Key", "%s");
        two.getTable("Second").addColumn("Value", "%d");
        two.getTable("Second").set(0, "Key", "k");
        two.getTable("Second").set(0, "Value", 7);
        emit(two, "two");

        // The three sortings, over rows deliberately out of order.
        for (final GATKReportTable.Sorting sorting : GATKReportTable.Sorting.values()) {
            final GATKReport sorted = new GATKReport();
            sorted.addTable("Sorted", "rows added out of order", 2, sorting);
            final GATKReportTable table = sorted.getTable("Sorted");
            table.addColumn("RowKey", "%s");
            table.addColumn("Value", "%d");
            table.set("bbb", "RowKey", "bbb");
            table.set("bbb", "Value", 2);
            table.set("aaa", "RowKey", "aaa");
            table.set("aaa", "Value", 1);
            table.set("ccc", "RowKey", "ccc");
            table.set("ccc", "Value", 3);
            emit(sorted, "sort-" + sorting.name());
        }
    }

    /**
     * One report: its text, its per-column widths, its lines with the spaces made visible, and
     * whether parsing it and writing it again produces the same bytes.
     */
    static void emit(final GATKReport report, final String label) throws Exception {
        final String text = render(report);
        System.out.printf("report\t%s\t%s%n", label, ReferenceQueryDump.escape(text));

        // The version line is the same for every report and is printed once, from the first.
        if (label.equals("types")) {
            System.out.printf("version\t%s%n", text.split("\n")[0]);
        }

        // Spaces decide this format and are invisible in a diff, so they travel as underscores.
        final String[] lines = text.split("\n", -1);
        for (int i = 0; i < lines.length; i++) {
            System.out.printf("line\t%s\t%d\t%s%n", label, i, lines[i].replace(' ', '_'));
        }

        for (final GATKReportTable table : report.getTables()) {
            for (final var column : table.getColumnInfo()) {
                System.out.printf("width\t%s\t%s\t%d\t%s%n", label, column.getColumnName(),
                        column.getColumnFormat().getWidth(),
                        column.getColumnFormat().getAlignment().name());
            }
        }

        // Parse what was written and write it again. The table header carries the formats but not
        // the widths, so the second writing recomputes every width from the values it parsed.
        final GATKReport reparsed =
                new GATKReport(new ByteArrayInputStream(text.getBytes(StandardCharsets.UTF_8)));
        System.out.printf("roundtrip\t%s\t%b%n", label, render(reparsed).equals(text));
    }

    static String render(final GATKReport report) {
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        try (final PrintStream out = new PrintStream(bytes, true, StandardCharsets.UTF_8)) {
            report.print(out);
        }
        return bytes.toString(StandardCharsets.UTF_8);
    }
}

//! Ported from `org.broadinstitute.hellbender.utils.report` (GATK 4.6.2.0): `GATKReport`,
//! `GATKReportTable`, `GATKReportColumn`, `GATKReportColumnFormat` and `GATKReportDataType`.
//!
//! The file format BQSR is written in. `ApplyBQSR` reads a recalibration table and
//! `BaseRecalibrator` writes one, so the bytes are settled here or they are settled twice.
//!
//! # The width is the widest formatted value, and it freezes
//!
//! `GATKReportColumn` starts at the column name's length and grows as values are added, measuring
//! the **formatted** value rather than the value. `getColumnFormat()` caches the result, so the
//! width a column reports after the first call is the width it will always report.
//!
//! # Two spaces between columns, and the padding is inside them
//!
//! `writeRow` prints `"  "` before every column but the first, and each value goes through
//! `%-<w>s` or `%<w>s`. A left-aligned last column therefore carries trailing spaces and a
//! right-aligned one carries none, which is invisible in a diff and decides the bytes.
//!
//! # Alignment is right until a value asks for left
//!
//! The default is [`Alignment::Right`]. A value that is not numeric and is not one of `null`,
//! `NA`, `Infinity`, `-Infinity`, `NaN` turns the whole column left-aligned. Measured: a column of
//! `true`/`false` comes out left-aligned, and a `%.4f` column holding `NaN` stays right-aligned.
//!
//! # Two escapes in `writeRow` that the table declaration does not announce
//!
//! A column whose format is the empty string is [`DataType::Unknown`], and a double in it is
//! written `%.8f` rather than through the column's `%s`. And a **non-finite** double escapes its
//! own format entirely: `Double.isFinite` decides, and `toString()` is used instead, so a `%.4f`
//! column can hold `-Infinity` with no decimal point.
//!
//! # The header carries the formats, not the widths
//!
//! `#:GATKTable:<columns>:<rows>:<format>:...:;` then `#:GATKTable:<name>:<description>`. A reader
//! recomputes every width from the values it parses, which is why a parse and a second writing
//! reproduce the first byte for byte.

use std::fmt::Write as _;

use crate::java_format::format_decimals;

/// `GATKReport.GATKREPORT_HEADER_PREFIX` and the version this port writes.
pub const REPORT_VERSION: &str = "v1.1";

/// `GATKReportColumnFormat.Alignment`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Alignment {
    Left,
    Right,
}

/// `GATKReportDataType`, as far as the writer needs it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DataType {
    /// An empty format string. A double in such a column is written `%.8f`.
    Unknown,
    Boolean,
    Character,
    Decimal,
    Integer,
    String,
}

impl DataType {
    /// `GATKReportDataType.fromFormatString`, matched on the conversion character.
    pub fn from_format(format: &str) -> DataType {
        if format.is_empty() {
            return DataType::Unknown;
        }
        match format.chars().last() {
            Some('b') | Some('B') => DataType::Boolean,
            Some('c') | Some('C') => DataType::Character,
            Some('e') | Some('E') | Some('f') | Some('F') => DataType::Decimal,
            Some('d') | Some('D') => DataType::Integer,
            Some('s') | Some('S') => DataType::String,
            _ => DataType::Unknown,
        }
    }
}

/// One cell. The variants are the Java types `GATKReportTable.set` accepts.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    Null,
    Int(i64),
    Double(f64),
    Bool(bool),
    Char(char),
    Str(String),
}

/// `GATKReportColumn.RIGHT_ALIGN_STRINGS`: the renderings that do not force a left alignment even
/// though they are not numeric.
const RIGHT_ALIGN_STRINGS: [&str; 5] = ["null", "NA", "Infinity", "-Infinity", "NaN"];

/// `GATKReportColumn.isRightAlign`: numeric, or one of the five strings above.
///
/// The value is **not** trimmed first, which the reference says is deliberate: spaces are taken to
/// mean the value is already padded.
pub fn is_right_align(value: &str) -> bool {
    if RIGHT_ALIGN_STRINGS.contains(&value) {
        return true;
    }
    value.parse::<f64>().is_ok()
}

/// `GATKReportColumn`.
#[derive(Debug, Clone)]
pub struct Column {
    pub name: String,
    /// The format as declared. An empty string becomes `%s` with [`DataType::Unknown`].
    pub format: String,
    pub data_type: DataType,
    max_width: usize,
    alignment: Alignment,
}

impl Column {
    pub fn new(name: &str, format: &str) -> Column {
        let (format, data_type) = if format.is_empty() {
            ("%s".to_string(), DataType::Unknown)
        } else {
            (format.to_string(), DataType::from_format(format))
        };
        Column {
            name: name.to_string(),
            format,
            data_type,
            // "this.maxWidth = columnName.length()" before any value is seen.
            max_width: name.chars().count(),
            alignment: Alignment::Right,
        }
    }

    /// `updateFormatting`: the width grows to the rendering, and one non-numeric value is enough to
    /// turn the whole column left-aligned.
    fn observe(&mut self, rendered: &str) {
        self.max_width = self.max_width.max(rendered.chars().count());
        if !is_right_align(rendered) {
            self.alignment = Alignment::Left;
        }
    }

    pub fn width(&self) -> usize {
        self.max_width
    }

    pub fn alignment(&self) -> Alignment {
        self.alignment
    }
}

/// `String.format(format, value)` for the conversions this format uses.
///
/// `%.Nf` goes through [`format_decimals`], which rounds HALF_UP as Java does rather than
/// half-to-even as Rust does.
pub fn format_value(format: &str, data_type: DataType, value: &Value) -> String {
    match value {
        // "if ( obj == null ) value = "null";", before the type is consulted at all.
        Value::Null => "null".to_string(),
        Value::Double(number) => {
            if data_type == DataType::Unknown {
                // The first escape: an untyped column renders a double `%.8f`.
                return format_decimals(*number, 8);
            }
            if !number.is_finite() {
                // The second: a non-finite double leaves its own format for `Double.toString`.
                return format_decimals(*number, 0)
                    .trim_end_matches('.')
                    .to_string();
            }
            apply_format(format, *number)
        }
        Value::Int(number) => number.to_string(),
        Value::Bool(flag) => flag.to_string(),
        Value::Char(character) => character.to_string(),
        Value::Str(text) => text.clone(),
    }
}

/// `GATKReportTable.ROW_COMPARATOR`, which is what `SORT_BY_COLUMN` sorts by.
///
/// **It compares the values, not their renderings.** The reference dispatches on the boxed type and
/// uses `compareTo` for every numeric one, falling back to `toString` only for the rest and for two
/// values of different classes. So a column of integers sorts 2 before 10, where comparing the
/// formatted text would put 10 first. `QuantizationInfo`'s quantization map is 94 rows keyed by
/// quality score and is the table that shows it.
///
/// **`toString` is the value's, not the column's format.** A double written `%.4f` still compares as
/// `Double.toString` when it reaches the fallback, which it does only against a value of a different
/// class.
///
/// Two places where this cannot be exactly the reference and the reason is in the type:
///
///  * Java boxes an `int` to `Integer` and a `long` to `Long`, and a column mixing them takes the
///    different-class branch and compares as text. [`Value::Int`] is one variant, so a mixed column
///    compares numerically here. Nothing in BQSR mixes them in one column;
///  * `String.compareTo` orders by UTF-16 code unit and Rust's `str` ordering is by byte, which
///    disagree only above the basic multilingual plane.
///
/// A null cell is a `NullPointerException` in the reference, because the comparator asks for its
/// class before anything else. Here it sorts as the string `null`, which is what the row would have
/// been written as.
fn compare_values(left: &Value, right: &Value) -> std::cmp::Ordering {
    match (left, right) {
        // `((Integer) a).compareTo((Integer) b)`, numeric.
        (Value::Int(left), Value::Int(right)) => left.cmp(right),
        // `((Double) a).compareTo((Double) b)`, which orders NaN last and -0.0 below 0.0.
        // `total_cmp` is that ordering exactly.
        (Value::Double(left), Value::Double(right)) => left.total_cmp(right),
        // Everything else, and every pair of different classes, is `toString().compareTo(...)`.
        _ => java_to_string(left).cmp(&java_to_string(right)),
    }
}

/// `Object.toString()` on a cell, which is not the column's rendering of it.
fn java_to_string(value: &Value) -> String {
    match value {
        Value::Null => "null".to_string(),
        Value::Int(number) => number.to_string(),
        // `Double.toString`, which the writer also reaches for a non-finite value.
        Value::Double(number) => format_decimals(*number, 0)
            .trim_end_matches('.')
            .to_string(),
        Value::Bool(flag) => flag.to_string(),
        Value::Char(character) => character.to_string(),
        Value::Str(text) => text.clone(),
    }
}

/// `%.Nf` and `%s` over a double, which is every numeric conversion this format uses.
fn apply_format(format: &str, number: f64) -> String {
    if let Some(rest) = format.strip_prefix("%.") {
        if let Some(digits) = rest.strip_suffix('f') {
            if let Ok(places) = digits.parse::<usize>() {
                return format_decimals(number, places);
            }
        }
    }
    // `%s` on a Double is `Double.toString`, which this port only needs for whole values.
    format_decimals(number, 0).trim_end_matches('.').to_string()
}

/// `GATKReportTable.Sorting`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sorting {
    SortByColumn,
    SortByRow,
    DoNotSort,
}

/// `GATKReportTable`.
#[derive(Debug, Clone)]
pub struct Table {
    pub name: String,
    pub description: String,
    pub columns: Vec<Column>,
    /// Parallel to `rows`: the row key each row was set under, for `SORT_BY_ROW`.
    pub row_keys: Vec<String>,
    pub rows: Vec<Vec<Value>>,
    pub sorting: Sorting,
}

impl Table {
    pub fn new(name: &str, description: &str, sorting: Sorting) -> Table {
        Table {
            name: name.to_string(),
            description: description.to_string(),
            columns: Vec::new(),
            row_keys: Vec::new(),
            rows: Vec::new(),
            sorting,
        }
    }

    pub fn add_column(&mut self, name: &str, format: &str) {
        self.columns.push(Column::new(name, format));
        for row in &mut self.rows {
            row.push(Value::Null);
        }
    }

    /// `set(rowKey, columnName, value)`, which creates the row on first use.
    pub fn set(&mut self, row_key: &str, column: &str, value: Value) {
        let index = match self.row_keys.iter().position(|key| key == row_key) {
            Some(index) => index,
            None => {
                self.row_keys.push(row_key.to_string());
                self.rows.push(vec![Value::Null; self.columns.len()]);
                self.rows.len() - 1
            }
        };
        let column_index = self
            .columns
            .iter()
            .position(|c| c.name == column)
            .expect("a declared column");
        let rendered = format_value(
            &self.columns[column_index].format,
            self.columns[column_index].data_type,
            &value,
        );
        self.columns[column_index].observe(&rendered);
        self.rows[index][column_index] = value;
    }

    /// The rows in the order this table's sorting writes them.
    fn ordered(&self) -> Vec<usize> {
        let mut indices: Vec<usize> = (0..self.rows.len()).collect();
        match self.sorting {
            Sorting::DoNotSort => {}
            // `new TreeMap<>(rowIdToIndex)`: the row keys, in their natural order.
            Sorting::SortByRow => {
                indices.sort_by(|a, b| self.row_keys[*a].cmp(&self.row_keys[*b]));
            }
            // `Collections.sort(underlyingData, ROW_COMPARATOR)`: column by column, left to right.
            Sorting::SortByColumn => {
                indices.sort_by(|a, b| {
                    for column in 0..self.columns.len() {
                        match compare_values(&self.rows[*a][column], &self.rows[*b][column]) {
                            std::cmp::Ordering::Equal => continue,
                            other => return other,
                        }
                    }
                    std::cmp::Ordering::Equal
                });
            }
        }
        indices
    }

    /// `GATKReportTable.write`.
    pub fn write(&self, out: &mut String) {
        let _ = write!(
            out,
            "#:GATKTable:{}:{}",
            self.columns.len(),
            self.rows.len()
        );
        for column in &self.columns {
            let _ = write!(out, ":{}", column.format);
        }
        // `ENDLINE` is ":;", printed by `println(ENDLINE)` after the last `:<format>`, so the line
        // ends with a separator the formats do not account for. Measured, not read off the source.
        out.push_str(":;\n");
        let _ = writeln!(out, "#:GATKTable:{}:{}", self.name, self.description);

        let mut first = true;
        for column in &self.columns {
            if !first {
                out.push_str("  ");
            }
            first = false;
            // The name always uses `%-<w>s`, whatever the column's alignment.
            let _ = write!(out, "{:<width$}", column.name, width = column.width());
        }
        out.push('\n');

        for index in self.ordered() {
            let mut first = true;
            for (position, column) in self.columns.iter().enumerate() {
                if !first {
                    out.push_str("  ");
                }
                first = false;
                let value = format_value(
                    &column.format,
                    column.data_type,
                    &self.rows[index][position],
                );
                match column.alignment() {
                    Alignment::Left => {
                        let _ = write!(out, "{:<width$}", value, width = column.width());
                    }
                    Alignment::Right => {
                        let _ = write!(out, "{:>width$}", value, width = column.width());
                    }
                }
            }
            out.push('\n');
        }
        // `out.println()` after the body: every table is followed by a blank line.
        out.push('\n');
    }
}

/// `GATKReport`.
#[derive(Debug, Clone, Default)]
pub struct Report {
    pub tables: Vec<Table>,
}

impl Report {
    pub fn new() -> Report {
        Report::default()
    }

    pub fn add_table(&mut self, table: Table) {
        self.tables.push(table);
    }

    pub fn table(&mut self, name: &str) -> &mut Table {
        self.tables
            .iter_mut()
            .find(|table| table.name == name)
            .expect("a declared table")
    }

    /// `GATKReport.print`: the version line, then every table.
    pub fn write(&self) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "#:GATKReport.{}:{}", REPORT_VERSION, self.tables.len());
        for table in &self.tables {
            table.write(&mut out);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn types_table() -> Report {
        let mut table = Table::new("Types", "one column per data type", Sorting::DoNotSort);
        table.add_column("Name", "%s");
        table.add_column("Count", "%d");
        table.add_column("Rate", "%.4f");
        table.add_column("Flag", "%b");
        table.add_column("Letter", "%c");
        table.add_column("Untyped", "");

        let rows: [(&str, Value, Value, Value, Value, Value, Value); 4] = [
            (
                "0",
                Value::Str("short".into()),
                Value::Int(1),
                Value::Double(0.5),
                Value::Bool(true),
                Value::Char('A'),
                Value::Double(0.5),
            ),
            (
                "1",
                Value::Str("a considerably longer value".into()),
                Value::Int(1234567),
                Value::Double(0.123456789),
                Value::Bool(false),
                Value::Char('z'),
                Value::Double(1.0 / 3.0),
            ),
            (
                "2",
                Value::Null,
                Value::Int(0),
                Value::Double(f64::NAN),
                Value::Bool(true),
                Value::Char('x'),
                Value::Double(f64::INFINITY),
            ),
            (
                "3",
                Value::Str("neg".into()),
                Value::Int(-42),
                Value::Double(f64::NEG_INFINITY),
                Value::Bool(false),
                Value::Char('y'),
                Value::Double(f64::NAN),
            ),
        ];
        for (key, name, count, rate, flag, letter, untyped) in rows {
            table.set(key, "Name", name);
            table.set(key, "Count", count);
            table.set(key, "Rate", rate);
            table.set(key, "Flag", flag);
            table.set(key, "Letter", letter);
            table.set(key, "Untyped", untyped);
        }
        let mut report = Report::new();
        report.add_table(table);
        report
    }

    #[test]
    fn the_widths_are_the_widest_formatted_value() {
        let report = types_table();
        let widths: Vec<(usize, Alignment)> = report.tables[0]
            .columns
            .iter()
            .map(|c| (c.width(), c.alignment()))
            .collect();
        assert_eq!(
            widths,
            vec![
                (27, Alignment::Left),  // "a considerably longer value"
                (7, Alignment::Right),  // "1234567"
                (9, Alignment::Right),  // "-Infinity", which does not force left
                (5, Alignment::Left),   // "false"
                (6, Alignment::Left),   // the name, wider than one character
                (10, Alignment::Right), // "0.50000000"
            ]
        );
    }

    #[test]
    fn the_two_escapes_in_write_row() {
        // An untyped column renders a double %.8f rather than through its own %s.
        assert_eq!(
            format_value("%s", DataType::Unknown, &Value::Double(1.0 / 3.0)),
            "0.33333333"
        );
        // A non-finite double leaves its own format.
        assert_eq!(
            format_value("%.4f", DataType::Decimal, &Value::Double(f64::NEG_INFINITY)),
            "-Infinity"
        );
        assert_eq!(
            format_value("%.4f", DataType::Decimal, &Value::Double(f64::NAN)),
            "NaN"
        );
        // And a finite one does not.
        assert_eq!(
            format_value("%.4f", DataType::Decimal, &Value::Double(0.123456789)),
            "0.1235"
        );
    }

    #[test]
    fn a_null_is_four_characters_whatever_the_column_is() {
        for (format, data_type) in [("%d", DataType::Integer), ("%.4f", DataType::Decimal)] {
            assert_eq!(format_value(format, data_type, &Value::Null), "null");
        }
    }

    #[test]
    fn the_five_strings_that_do_not_force_a_left_alignment() {
        for value in ["null", "NA", "Infinity", "-Infinity", "NaN", "42", "-0.5"] {
            assert!(is_right_align(value), "{value}");
        }
        for value in ["true", "false", "A", "short"] {
            assert!(!is_right_align(value), "{value}");
        }
    }

    #[test]
    fn the_table_header_carries_the_formats_and_not_the_widths() {
        let text = types_table().write();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "#:GATKReport.v1.1:1");
        assert_eq!(lines[1], "#:GATKTable:6:4:%s:%d:%.4f:%b:%c:%s:;");
        assert_eq!(lines[2], "#:GATKTable:Types:one column per data type");
    }

    #[test]
    fn a_row_is_two_spaces_between_columns_and_padding_inside_them() {
        let text = types_table().write();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(
            lines[3].replace(' ', "_"),
            "Name_________________________Count____Rate_______Flag___Letter__Untyped___"
        );
        assert_eq!(
            lines[6].replace(' ', "_"),
            "null_______________________________0________NaN__true___x_________Infinity"
        );
    }

    #[test]
    fn every_table_is_followed_by_a_blank_line() {
        let text = types_table().write();
        assert!(text.ends_with("\n\n"), "{text:?}");
    }
}

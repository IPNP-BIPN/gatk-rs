//! `TableWriter` and `TableReader`, ported from `org.broadinstitute.hellbender.utils.tsv`
//! (GATK 4.6.2.0).
//!
//! GATK's own table format, under every `.table` file the contamination tools read and write. It is
//! not a plain TSV: it is opencsv configured with a tab separator, a quote character and a backslash
//! escape, plus a comment convention that carries metadata.
//!
//! # Metadata is a tagged comment, and the tag is easy to miss
//!
//! ```java
//! if (commentText.startsWith(TableWriter.METADATA_TAG)) {
//!     final String[] keyAndValue = commentText.substring(...).split("=");
//!     metadata.put(keyAndValue[0], keyAndValue[1]);
//! }
//! ```
//!
//! The line is `#<METADATA>key=value`. A hand-written `#sample=s1` looks like metadata and is not.
//! And the map is filled by `processCommentLine` **itself**, so a subclass overriding that method
//! without calling super loses every pair: the same file read twice, once with an override and once
//! without, gives an empty map and `sample=s1`. Nothing in the signature says so, which is why the
//! reader here separates the two: [`Table::parse`] always collects the metadata, and the comments
//! are returned beside it rather than through a hook a caller can swallow.
//!
//! # A value is quoted only when it has to be
//!
//! Measured: a space and a comma pass through bare, a tab and a quote force quotes, a backslash is
//! escaped **and** forces quotes, and a value beginning with the comment prefix is written unquoted.

use std::collections::HashMap;
use std::fmt::Write as _;

/// `TableUtils.COLUMN_SEPARATOR`.
pub const COLUMN_SEPARATOR: char = '\t';
/// `TableUtils.COMMENT_PREFIX`.
pub const COMMENT_PREFIX: &str = "#";
/// `TableWriter.METADATA_TAG`, which follows the comment prefix.
pub const METADATA_TAG: &str = "<METADATA>";

/// What the reader refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TableError {
    /// A row whose column count is not the header's.
    ColumnCount {
        source: String,
        line: u64,
        values: usize,
        columns: usize,
    },
    /// A column the caller asked for and the header does not have.
    NoSuchColumn(String),
    /// A field asked for as a number that is not one. The message names the COLUMN and prints the
    /// offending value unquoted, so an empty field ends the sentence with a space.
    NotAnInteger {
        source: String,
        line: u64,
        column: String,
        value: String,
    },
}

impl TableError {
    /// The message the reference carries.
    pub fn message(&self) -> String {
        match self {
            TableError::ColumnCount {
                source,
                line,
                values,
                columns,
            } => format!(
                "format error in '{source}' at line {line}: mismatch between number of values in line ({values}) and number of columns ({columns})"
            ),
            TableError::NoSuchColumn(name) => format!("there is no such column: {name}"),
            TableError::NotAnInteger {
                source,
                line,
                column,
                value,
            } => format!(
                "format error in '{source}' at line {line}: expected int value for column {column} but found {value}"
            ),
        }
    }

    /// The Java class, which is not the same for all three.
    pub fn java_class(&self) -> &'static str {
        match self {
            TableError::NoSuchColumn(_) => "java.lang.IllegalArgumentException",
            _ => "org.broadinstitute.hellbender.exceptions.UserException$BadInput",
        }
    }
}

/// One table: its columns, its rows, its metadata and the comments that were not metadata.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Table {
    pub columns: Vec<String>,
    pub rows: Vec<Vec<String>>,
    /// Every `#<METADATA>key=value`, in the order they were seen: a repeated key keeps the last.
    pub metadata: HashMap<String, String>,
    /// Every comment line, tagged or not, with the prefix stripped.
    pub comments: Vec<String>,
    /// The line number each row was read from, which the numeric refusals name.
    pub row_lines: Vec<u64>,
}

impl Table {
    /// `TableReader`: the header, the rows and the comments of one file.
    ///
    /// `source` is the name the error messages use, which is the path the reference was given.
    pub fn parse(text: &str, source: &str) -> Result<Table, TableError> {
        let mut table = Table::default();
        let mut header_seen = false;

        for (index, line) in text.lines().enumerate() {
            let number = index as u64 + 1;
            let values = split_line(line);
            // A comment is recognised **after** the line is parsed, on the first field, and its
            // text is rebuilt by joining the fields with tabs. So a metadata value holding a tab,
            // which the writer quotes, is still read as one comment and not as a row.
            if values[0].starts_with(COMMENT_PREFIX) {
                let mut comment = values[0][COMMENT_PREFIX.len()..].to_string();
                for value in &values[1..] {
                    comment.push(COLUMN_SEPARATOR);
                    comment.push_str(value);
                }
                if let Some(pair) = comment.strip_prefix(METADATA_TAG) {
                    // `split("=")` with no limit, so a value containing `=` keeps only its first
                    // piece, and the reference would throw on a line with no `=` at all.
                    let mut pieces = pair.split('=');
                    if let (Some(key), Some(value)) = (pieces.next(), pieces.next()) {
                        table.metadata.insert(key.to_string(), value.to_string());
                    }
                }
                table.comments.push(comment);
                continue;
            }
            if !header_seen {
                table.columns = values;
                header_seen = true;
                continue;
            }
            if values.len() != table.columns.len() {
                return Err(TableError::ColumnCount {
                    source: source.to_string(),
                    line: number,
                    values: values.len(),
                    columns: table.columns.len(),
                });
            }
            table.rows.push(values);
            table.row_lines.push(number);
        }
        Ok(table)
    }

    /// `DataLine.get(column)`.
    pub fn get<'a>(&'a self, row: &'a [String], column: &str) -> Result<&'a str, TableError> {
        let index = self
            .columns
            .iter()
            .position(|name| name == column)
            .ok_or_else(|| TableError::NoSuchColumn(column.to_string()))?;
        Ok(&row[index])
    }

    /// `DataLine.getInt(column)`, whose refusal names the value and not the column.
    pub fn get_int(
        &self,
        row: &[String],
        column: &str,
        source: &str,
        line: u64,
    ) -> Result<i32, TableError> {
        let value = self.get(row, column)?;
        value.parse::<i32>().map_err(|_| TableError::NotAnInteger {
            source: source.to_string(),
            line,
            column: column.to_string(),
            value: value.to_string(),
        })
    }
}

/// `TableWriter`: the metadata, the header and the rows.
///
/// The header is written whether or not any record follows, which is what an empty table with a
/// header shows.
pub fn write_table(columns: &[&str], rows: &[Vec<String>], metadata: &[(&str, &str)]) -> String {
    let mut out = String::new();
    for (key, value) in metadata {
        // `writeComment` hands the whole line to the same csv writer a row goes through, so a
        // value holding a tab gets the **comment line** quoted rather than escaped in place.
        let _ = writeln!(
            out,
            "{}",
            quote_if_needed(&format!("{COMMENT_PREFIX}{METADATA_TAG}{key}={value}"))
        );
    }
    let _ = writeln!(
        out,
        "{}",
        columns
            .iter()
            .map(|column| quote_if_needed(column))
            .collect::<Vec<String>>()
            .join(&COLUMN_SEPARATOR.to_string())
    );
    for row in rows {
        let _ = writeln!(
            out,
            "{}",
            row.iter()
                .map(|value| quote_if_needed(value))
                .collect::<Vec<String>>()
                .join(&COLUMN_SEPARATOR.to_string())
        );
    }
    out
}

/// `DataLine.set(index, double)`, whose short form never survives.
///
/// ```java
/// final long rounded = Math.round(value);
/// if (rounded == value) {
///     set(index, Long.toString(rounded));
/// } else {
///     set(index, Double.toString(value));
/// }
/// return set(index, Double.toString(value));
/// ```
///
/// The `return` writes over whatever the branch above chose, so the rounding is dead: `1.0` is
/// written `1.0` and never `1`. A port that read the intent rather than the code would differ on
/// every integral frequency, which is why this function is only `Double.toString`.
///
/// Java prints the shortest decimal that round-trips, with `E` notation outside `[1e-3, 1e7)` and
/// always a digit after the point. Rust's `{}` and `{:e}` have the same shortest-round-trip
/// property, so the two agree on everything the contamination tables carry; the pre-JDK19
/// `FloatingDecimal` can emit one digit more than needed for a few values needing 16 or 17
/// significant digits, and the golden is what pins the ones this port writes.
pub fn java_double_to_string(value: f64) -> String {
    if value.is_nan() {
        return "NaN".to_string();
    }
    if value.is_infinite() {
        return if value > 0.0 { "Infinity" } else { "-Infinity" }.to_string();
    }
    let magnitude = value.abs();
    if magnitude != 0.0 && !(1e-3..1e7).contains(&magnitude) {
        let rendered = format!("{value:e}");
        let (mantissa, exponent) = rendered
            .split_once('e')
            .expect("an exponent form has its `e`");
        let mantissa = if mantissa.contains('.') {
            mantissa.to_string()
        } else {
            format!("{mantissa}.0")
        };
        return format!("{mantissa}E{exponent}");
    }
    let rendered = format!("{value}");
    if rendered.contains('.') {
        rendered
    } else {
        format!("{rendered}.0")
    }
}

/// opencsv's `processLine` with `applyQuotesToAll = false`.
///
/// A value is quoted when it holds the separator, a quote or a newline; a quote and a backslash are
/// escaped with a backslash. A comma and a space force nothing, which a port that reused a CSV
/// writer would get wrong in the other direction.
pub fn quote_if_needed(value: &str) -> String {
    let needs_quotes = value.contains(COLUMN_SEPARATOR)
        || value.contains('"')
        || value.contains('\n')
        || value.contains('\\');
    if !needs_quotes {
        return value.to_string();
    }
    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for character in value.chars() {
        if character == '"' || character == '\\' {
            out.push('\\');
        }
        out.push(character);
    }
    out.push('"');
    out
}

/// opencsv's parser: tab-separated, with quotes and backslash escapes.
fn split_line(line: &str) -> Vec<String> {
    let mut values = Vec::new();
    let mut current = String::new();
    let mut in_quotes = false;
    let mut characters = line.chars().peekable();

    while let Some(character) = characters.next() {
        match character {
            '\\' if in_quotes => {
                if let Some(next) = characters.next() {
                    current.push(next);
                }
            }
            '"' => in_quotes = !in_quotes,
            c if c == COLUMN_SEPARATOR && !in_quotes => {
                values.push(std::mem::take(&mut current));
            }
            c => current.push(c),
        }
    }
    values.push(current);
    values
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_tab_and_a_quote_force_quotes_and_a_comma_does_not() {
        assert_eq!(quote_if_needed("has a space"), "has a space");
        assert_eq!(quote_if_needed("has a , comma"), "has a , comma");
        assert_eq!(quote_if_needed("has\ta tab"), "\"has\ta tab\"");
        assert_eq!(
            quote_if_needed("has a \"quote\""),
            "\"has a \\\"quote\\\"\""
        );
        assert_eq!(quote_if_needed("#comment"), "#comment");
        assert_eq!(quote_if_needed(""), "");
    }

    #[test]
    fn a_hand_written_comment_is_not_metadata() {
        let tagged = Table::parse("#<METADATA>sample=s1\na\tb\n1\t2\n", "x").expect("parses");
        assert_eq!(
            tagged.metadata.get("sample").map(String::as_str),
            Some("s1")
        );

        let untagged = Table::parse("#sample=s1\na\tb\n1\t2\n", "x").expect("parses");
        assert!(untagged.metadata.is_empty());
        // It is still kept as a comment.
        assert_eq!(untagged.comments, vec!["sample=s1".to_string()]);
    }

    #[test]
    fn a_row_of_the_wrong_width_names_the_line_and_the_counts() {
        let error = Table::parse("a\tb\tc\n1\t2\n", "file.table").unwrap_err();
        assert_eq!(
            error.message(),
            "format error in 'file.table' at line 2: mismatch between number of values in line (2) and number of columns (3)"
        );
    }

    #[test]
    fn a_missing_column_is_a_different_class_from_a_bad_row() {
        let table = Table::parse("a\tb\n1\t2\n", "x").expect("parses");
        let error = table.get(&table.rows[0], "c").unwrap_err();
        assert_eq!(error.java_class(), "java.lang.IllegalArgumentException");
        assert_eq!(error.message(), "there is no such column: c");
    }

    #[test]
    fn a_double_keeps_its_point_and_the_short_form_is_dead_code() {
        assert_eq!(java_double_to_string(1.0), "1.0");
        assert_eq!(java_double_to_string(0.0), "0.0");
        assert_eq!(java_double_to_string(-0.0), "-0.0");
        assert_eq!(java_double_to_string(1.0 / 3.0), "0.3333333333333333");
        assert_eq!(java_double_to_string(1e-7), "1.0E-7");
        assert_eq!(java_double_to_string(1e21), "1.0E21");
        assert_eq!(java_double_to_string(f64::NAN), "NaN");
        assert_eq!(java_double_to_string(f64::INFINITY), "Infinity");
        assert_eq!(java_double_to_string(f64::NEG_INFINITY), "-Infinity");
    }

    #[test]
    fn a_metadata_value_with_a_tab_quotes_the_whole_comment_line() {
        let text = write_table(&["a"], &[], &[("SAMPLE", "has\ta tab")]);
        assert!(
            text.starts_with("\"#<METADATA>SAMPLE=has\ta tab\"\n"),
            "{text}"
        );
        // And it is read back as one comment, not as a row.
        let table = Table::parse(&text, "x").expect("parses");
        assert_eq!(
            table.metadata.get("SAMPLE").map(String::as_str),
            Some("has\ta tab")
        );
        assert!(table.rows.is_empty());
    }

    #[test]
    fn the_header_is_written_even_with_no_rows() {
        let text = write_table(&["a", "b"], &[], &[("sample", "s1")]);
        assert_eq!(text, "#<METADATA>sample=s1\na\tb\n");
    }
}

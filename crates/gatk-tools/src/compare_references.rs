//! `CompareReferences`, ported from
//! `org.broadinstitute.hellbender.tools.reference.CompareReferences`,
//! `ReferenceSequenceTable` and `ReferencePair` (GATK 4.6.2.0).
//!
//! Several references compared by the MD5 of each sequence.
//!
//! # The table is keyed by MD5
//!
//! Two references whose sequences agree base for base share one row however differently they name
//! them, and one name over two different sequences is two rows. That is the whole design: the
//! names are the CELLS, not the keys.
//!
//! # The analysis removes a flag and then adds others
//!
//! ```java
//! analysis = EnumSet.of(Status.EXACT_MATCH);
//! ```
//!
//! Every pair starts as an exact match and loses that on the first disagreement, so a pair can end
//! up carrying two flags at once. And the printing is `EnumSet` order, which is the enum's
//! declaration order rather than the order they were added, so `DIFFER_IN_SEQUENCE_NAMES` always
//! precedes `DIFFER_IN_SEQUENCES_PRESENT`.
//!
//! # Superset and subset are a replacement, not a conclusion
//!
//! `DIFFER_IN_SEQUENCES_PRESENT` is REMOVED and one of the two added, but only when every missing
//! entry points the same way AND no naming discrepancy was found. A pair that both renames and
//! omits keeps both flags, which the golden's three-reference run shows.
//!
//! # The MD5 mode decides what is read
//!
//! `USE_DICT` refuses a dictionary with no `M5` and TRUSTS one that lies, which produces a table
//! with a row per lie. `ALWAYS_RECALCULATE` ignores what the dictionary says.

use std::collections::BTreeSet;

/// `ReferenceSequenceTable.MISSING_ENTRY_DISPLAY_STRING`.
pub const MISSING_ENTRY: &str = "---";

/// One reference, as the table reads it: the file's name and its dictionary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Reference {
    /// `getReferenceColumnName`, which is the file NAME and not the path.
    pub column: String,
    pub sequences: Vec<Sequence>,
}

/// One `@SQ` line, as far as the comparison reads one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sequence {
    pub name: String,
    pub length: i64,
    /// `M5`, absent when the dictionary has none.
    pub md5: Option<String>,
    /// What recalculating would give, which the fixtures carry so the modes can be compared
    /// without a fasta reader here.
    pub calculated_md5: String,
}

/// `CompareReferences.MD5CalculationMode`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Md5Mode {
    UseDict,
    RecalculateIfMissing,
    AlwaysRecalculate,
}

/// What the table refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TableError {
    /// `USE_DICT` over a sequence with no `M5`.
    MissingMd5 { sequence: String },
    /// A name found in more than two rows for one pair.
    DuplicateSequence {
        sequence: String,
        first: String,
        second: String,
    },
}

impl TableError {
    pub fn java_class(&self) -> &str {
        "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
    }

    pub fn message(&self) -> String {
        match self {
            TableError::MissingMd5 { sequence } => format!(
                "Bad input: Running in USE_DICT mode, but MD5 missing for sequence {sequence}. Run \
                 --md5-calculation-mode with a different mode to recalculate MD5."
            ),
            TableError::DuplicateSequence {
                sequence,
                first,
                second,
            } => format!(
                "Bad input: Duplicate of sequence '{sequence}' found in {first} or {second}."
            ),
        }
    }
}

/// One row of the table: an MD5, a length, and one cell per reference.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Row {
    pub md5: String,
    pub length: i64,
    /// The sequence's name in each reference, in the references' own order, or `None` for a
    /// reference that does not have it.
    pub cells: Vec<Option<String>>,
}

impl Row {
    /// The cell as the table prints it.
    pub fn cell(&self, index: usize) -> &str {
        self.cells[index].as_deref().unwrap_or(MISSING_ENTRY)
    }
}

/// The built table, in MD5 insertion order, which is the order the rows are written.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Table {
    pub columns: Vec<String>,
    pub rows: Vec<Row>,
}

/// `getMD5ForRecord`, whose answer depends on the mode alone.
fn md5_of(sequence: &Sequence, mode: Md5Mode) -> Result<String, TableError> {
    match mode {
        Md5Mode::UseDict => match &sequence.md5 {
            Some(md5) if !md5.is_empty() => Ok(md5.clone()),
            _ => Err(TableError::MissingMd5 {
                sequence: sequence.name.clone(),
            }),
        },
        Md5Mode::RecalculateIfMissing => match &sequence.md5 {
            Some(md5) if !md5.is_empty() => Ok(md5.clone()),
            _ => Ok(sequence.calculated_md5.clone()),
        },
        Md5Mode::AlwaysRecalculate => Ok(sequence.calculated_md5.clone()),
    }
}

/// `ReferenceSequenceTable.build`.
pub fn build(references: &[Reference], mode: Md5Mode) -> Result<Table, TableError> {
    let mut columns = vec!["MD5".to_string(), "Length".to_string()];
    for reference in references {
        columns.push(reference.column.clone());
    }
    let mut rows: Vec<Row> = Vec::new();
    for (index, reference) in references.iter().enumerate() {
        for sequence in &reference.sequences {
            let md5 = md5_of(sequence, mode)?;
            let position = match rows.iter().position(|row| row.md5 == md5) {
                Some(position) => position,
                None => {
                    rows.push(Row {
                        md5: md5.clone(),
                        length: sequence.length,
                        cells: vec![None; references.len()],
                    });
                    rows.len() - 1
                }
            };
            rows[position].cells[index] = Some(sequence.name.clone());
        }
    }
    Ok(Table { columns, rows })
}

/// `ReferencePair.Status`, in the enum's own declaration order, which is also `EnumSet`'s.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Status {
    ExactMatch,
    DifferInSequenceNames,
    DifferInSequence,
    DifferInSequencesPresent,
    Superset,
    Subset,
}

impl Status {
    pub fn name(&self) -> &'static str {
        match self {
            Status::ExactMatch => "EXACT_MATCH",
            Status::DifferInSequenceNames => "DIFFER_IN_SEQUENCE_NAMES",
            Status::DifferInSequence => "DIFFER_IN_SEQUENCE",
            Status::DifferInSequencesPresent => "DIFFER_IN_SEQUENCES_PRESENT",
            Status::Superset => "SUPERSET",
            Status::Subset => "SUBSET",
        }
    }
}

/// One compared pair and what the analysis concluded.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Pair {
    pub first: String,
    pub second: String,
    /// A `BTreeSet` so the order is the enum's, which is what `EnumSet` gives the printing.
    pub analysis: BTreeSet<Status>,
}

impl Pair {
    /// `ReferencePair.toString`.
    pub fn rendered(&self) -> String {
        let mut out = format!("REFERENCE PAIR: {}, {}\nStatus:\n", self.first, self.second);
        for status in &self.analysis {
            out.push_str(&format!("\t{}\n", status.name()));
        }
        out
    }
}

/// `compareAllReferences`: every pair in index order, analysed.
pub fn compare_all(table: &Table, references: &[Reference]) -> Result<Vec<Pair>, TableError> {
    let mut pairs = Vec::new();
    for first in 0..references.len() {
        for second in (first + 1)..references.len() {
            pairs.push(analyse(table, references, first, second)?);
        }
    }
    Ok(pairs)
}

/// `analyzeTable` for one pair.
fn analyse(
    table: &Table,
    references: &[Reference],
    first: usize,
    second: usize,
) -> Result<Pair, TableError> {
    let mut analysis: BTreeSet<Status> = BTreeSet::new();
    analysis.insert(Status::ExactMatch);

    for row in &table.rows {
        let left = row.cell(first);
        let right = row.cell(second);
        if left != right {
            analysis.remove(&Status::ExactMatch);
        }
        // Both present and different is the same sequence under two names.
        if left != right && row.cells[first].is_some() && row.cells[second].is_some() {
            analysis.insert(Status::DifferInSequenceNames);
        }
    }

    let mut superset = false;
    let mut subset = false;
    // `tableBySequenceName.keySet()`, which is every name either reference used, in the order the
    // build met them.
    let mut names: Vec<&str> = Vec::new();
    for reference in references {
        for sequence in &reference.sequences {
            if !names.contains(&sequence.name.as_str()) {
                names.push(&sequence.name);
            }
        }
    }
    for name in names {
        let mut found_in_one = 0;
        for row in &table.rows {
            let left = row.cell(first);
            let right = row.cell(second);
            if !(left == name || right == name) {
                continue;
            }
            let left_empty = row.cells[first].is_none();
            let right_empty = row.cells[second].is_none();
            if ((left == name) ^ (right == name)) && (left_empty ^ right_empty) {
                found_in_one += 1;
            }
            if left_empty ^ right_empty {
                if left_empty {
                    subset = true;
                } else {
                    superset = true;
                }
            }
        }
        match found_in_one {
            2 => {
                analysis.insert(Status::DifferInSequence);
            }
            1 => {
                analysis.insert(Status::DifferInSequencesPresent);
            }
            count if count > 2 => {
                return Err(TableError::DuplicateSequence {
                    sequence: name.to_string(),
                    first: references[first].column.clone(),
                    second: references[second].column.clone(),
                })
            }
            _ => {}
        }
    }

    if (superset ^ subset) && !analysis.contains(&Status::DifferInSequenceNames) {
        analysis.remove(&Status::DifferInSequencesPresent);
        if superset {
            analysis.insert(Status::Superset);
        } else {
            analysis.insert(Status::Subset);
        }
    }

    Ok(Pair {
        first: references[first].column.clone(),
        second: references[second].column.clone(),
        analysis,
    })
}

/// `writeTable`: the tab-delimited file, header included.
pub fn write_table(table: &Table) -> String {
    let mut out = table.columns.join("\t");
    out.push('\n');
    for row in &table.rows {
        out.push_str(&row.md5);
        out.push('\t');
        out.push_str(&row.length.to_string());
        for index in 0..row.cells.len() {
            out.push('\t');
            out.push_str(row.cell(index));
        }
        out.push('\n');
    }
    out
}

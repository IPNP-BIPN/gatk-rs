//! The tranches file, ported from
//! `org.broadinstitute.hellbender.tools.walkers.vqsr.Tranche` and `TruthSensitivityTranche`
//! (GATK 4.6.2.0).
//!
//! The file `VariantRecalibrator` writes and `ApplyVQSR` reads: eleven comma-separated columns per
//! tranche, named by a header line rather than by position.
//!
//! # The two optional columns are not optional
//!
//! ```java
//! getOptionalInteger(bindings, "numKnown", -1),
//! ...
//! if ( numKnown < 0 || numNovel < 0) {
//!     throw new GATKException("Invalid tranche " + name + " - no. variants is < 0 : known " + numKnown + " novel " + numNovel);
//! }
//! ```
//!
//! `numKnown` is read with a default of `-1` and the constructor then refuses anything negative, so
//! a header that does not name the column is fatal rather than defaulted. `accessibleTruthSites` and
//! `callsAtTruthSites` default to `-1` in the same way and are never checked, which is what makes
//! `getTruthSensitivity` answer `0.0` for them.
//!
//! # The header names the columns, and there must be eleven of them
//!
//! The bindings are built header-to-value, so the eleven may be in any order, while a header of
//! another length is refused before any row is read and a row of another length against the
//! header's. The two length refusals name the file; the missing-key and bad-value refusals build
//! their `MalformedFile` with no file at all and come out as **`Unknown file is malformed`**.
//!
//! # Two of the columns are not read defensively at all
//!
//! `model` and `filterName` are taken with a bare `bindings.get`, so a header that does not name the
//! model column reaches `Mode.valueOf(null)` and the refusal is Java's own
//! `NullPointerException: Name is null`. And `numNovel` is a `long` field parsed with
//! `Integer.valueOf`, so a count `VariantRecalibrator` could have written cannot be read back past
//! 2^31: it comes out `Invalid value for key numNovel`, exactly as a value that is not a number.

/// `EXPECTED_COLUMN_COUNT`.
pub const EXPECTED_COLUMN_COUNT: usize = 11;
/// `COMMENT_STRING`.
pub const COMMENT_STRING: &str = "#";
/// `VALUE_SEPARATOR`.
pub const VALUE_SEPARATOR: char = ',';

/// `VariantRecalibratorArgumentCollection.Mode`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    Snp,
    Indel,
    Both,
}

impl Mode {
    /// `Mode.valueOf(name)`, whose two failures are two different exceptions.
    pub fn value_of(name: Option<&str>) -> Result<Mode, TrancheError> {
        match name {
            None => Err(TrancheError::NullModel),
            Some("SNP") => Ok(Mode::Snp),
            Some("INDEL") => Ok(Mode::Indel),
            Some("BOTH") => Ok(Mode::Both),
            Some(other) => Err(TrancheError::UnknownModel(other.to_string())),
        }
    }

    /// `toString()`, which is the constant's own name.
    pub fn name(&self) -> &'static str {
        match self {
            Mode::Snp => "SNP",
            Mode::Indel => "INDEL",
            Mode::Both => "BOTH",
        }
    }
}

/// What reading a tranches file refuses with.
#[derive(Debug, Clone, PartialEq)]
pub enum TrancheError {
    /// A header whose length is not eleven. Names the file.
    HeaderLength { file: String, line: String },
    /// A row whose length is not the header's. Names the file.
    RowLength {
        file: String,
        header: usize,
        values: usize,
        line: String,
    },
    /// A key one of the `getRequired*` helpers did not find. Names no file.
    MissingRequiredKey(String),
    /// A value one of the helpers could not parse, which includes a `numNovel` past 2^31.
    InvalidValue(String),
    /// `Mode.valueOf(null)`, from a header that does not name the model column.
    NullModel,
    /// A VQSLOD tranche file whose version line is not six. The message names no file at all: the
    /// reference concatenates an empty string where the path should be, leaving two spaces.
    VqslodVersion { found: String, expected: i32 },
    /// `Mode.valueOf` on something that is not a constant.
    UnknownModel(String),
    /// The `TruthSensitivityTranche` constructor's own range check.
    UnreasonableTargetFdr(f64),
    /// The `Tranche` constructor's check, which the `-1` default of `numKnown` walks straight into.
    NegativeCounts {
        name: String,
        known: i64,
        novel: i64,
    },
}

impl TrancheError {
    /// The exception class the reference throws.
    pub fn class(&self) -> &'static str {
        match self {
            TrancheError::HeaderLength { .. }
            | TrancheError::RowLength { .. }
            | TrancheError::MissingRequiredKey(_)
            | TrancheError::InvalidValue(_) => {
                "org.broadinstitute.hellbender.exceptions.UserException$MalformedFile"
            }
            TrancheError::NullModel => "java.lang.NullPointerException",
            TrancheError::VqslodVersion { .. } => {
                "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
            }
            TrancheError::UnknownModel(_) => "java.lang.IllegalArgumentException",
            TrancheError::UnreasonableTargetFdr(_) | TrancheError::NegativeCounts { .. } => {
                "org.broadinstitute.hellbender.exceptions.GATKException"
            }
        }
    }

    /// The message, whose spacing is the reference's.
    pub fn message(&self) -> String {
        match self {
            TrancheError::HeaderLength { file, line } => format!(
                "File {file} is malformed: Expected 11 elements in header line {line}"
            ),
            TrancheError::RowLength {
                file,
                header,
                values,
                line,
            } => format!(
                "File {file} is malformed: Line had too few/many fields.  Header = {header} vals {values}. The line was: {line}"
            ),
            // Two spaces after the full stop, and no file: the exception was built without one.
            TrancheError::MissingRequiredKey(key) => format!(
                "Unknown file is malformed: Malformed tranches file.  Missing required key {key}"
            ),
            // One space, and the same no-file wording.
            TrancheError::InvalidValue(key) => format!(
                "Unknown file is malformed: Malformed tranches file. Invalid value for key {key}"
            ),
            TrancheError::NullModel => "Name is null".to_string(),
            TrancheError::VqslodVersion { found, expected } => format!(
                "Bad input: The file  contains version {found} tranches, but VQSLOD tranche parsing requires version {expected}"
            ),
            TrancheError::UnknownModel(value) => format!(
                "No enum constant org.broadinstitute.hellbender.tools.walkers.vqsr.VariantRecalibratorArgumentCollection.Mode.{value}"
            ),
            TrancheError::UnreasonableTargetFdr(value) => format!(
                "Target FDR is unreasonable {}",
                crate::tsv_table::java_double_to_string(*value)
            ),
            TrancheError::NegativeCounts {
                name,
                known,
                novel,
            } => format!("Invalid tranche {name} - no. variants is < 0 : known {known} novel {novel}"),
        }
    }
}

/// One `TruthSensitivityTranche`.
#[derive(Debug, Clone, PartialEq)]
pub struct TruthSensitivityTranche {
    pub target_truth_sensitivity: f64,
    pub min_vqslod: f64,
    pub num_known: i64,
    pub known_titv: f64,
    pub num_novel: i64,
    pub novel_titv: f64,
    pub accessible_truth_sites: i32,
    pub calls_at_truth_sites: i32,
    pub model: Mode,
    /// The `filterName` column, which is what an `ApplyVQSR` FILTER line is named after.
    pub name: String,
}

impl TruthSensitivityTranche {
    /// `getTruthSensitivity()`, which is `0.0` rather than a division when nothing is accessible.
    pub fn truth_sensitivity(&self) -> f64 {
        if self.accessible_truth_sites > 0 {
            self.calls_at_truth_sites as f64 / (1.0 * self.accessible_truth_sites as f64)
        } else {
            0.0
        }
    }
}

/// `Double.valueOf`, as the reader's helpers call it.
fn required_double(bindings: &[(String, String)], key: &str) -> Result<f64, TrancheError> {
    match lookup(bindings, key) {
        Some(value) => value
            .parse::<f64>()
            .map_err(|_| TrancheError::InvalidValue(key.to_string())),
        None => Err(TrancheError::MissingRequiredKey(key.to_string())),
    }
}

fn optional_double(
    bindings: &[(String, String)],
    key: &str,
    default: f64,
) -> Result<f64, TrancheError> {
    match lookup(bindings, key) {
        Some(value) => value
            .parse::<f64>()
            .map_err(|_| TrancheError::InvalidValue(key.to_string())),
        None => Ok(default),
    }
}

/// `Integer.valueOf`, which is where a `numNovel` past 2^31 stops.
fn required_integer(bindings: &[(String, String)], key: &str) -> Result<i32, TrancheError> {
    match lookup(bindings, key) {
        Some(value) => value
            .parse::<i32>()
            .map_err(|_| TrancheError::InvalidValue(key.to_string())),
        None => Err(TrancheError::MissingRequiredKey(key.to_string())),
    }
}

fn optional_integer(
    bindings: &[(String, String)],
    key: &str,
    default: i32,
) -> Result<i32, TrancheError> {
    match lookup(bindings, key) {
        Some(value) => value
            .parse::<i32>()
            .map_err(|_| TrancheError::InvalidValue(key.to_string())),
        None => Ok(default),
    }
}

fn lookup<'a>(bindings: &'a [(String, String)], key: &str) -> Option<&'a str> {
    bindings
        .iter()
        .find(|(name, _)| name == key)
        .map(|(_, value)| value.as_str())
}

/// `readTranches(f)`: every non-comment line, the first of them the header.
///
/// The result is sorted by `targetTruthSensitivity`, which is not the order the writer sorts by.
pub fn read_tranches(file: &str, text: &str) -> Result<Vec<TruthSensitivityTranche>, TrancheError> {
    let mut header: Option<Vec<String>> = None;
    let mut tranches: Vec<TruthSensitivityTranche> = Vec::new();

    for line in text.lines() {
        if line.starts_with(COMMENT_STRING) {
            continue;
        }
        let values: Vec<String> = line.split(VALUE_SEPARATOR).map(|v| v.to_string()).collect();
        let Some(columns) = header.as_ref() else {
            if values.len() != EXPECTED_COLUMN_COUNT {
                return Err(TrancheError::HeaderLength {
                    file: file.to_string(),
                    line: line.to_string(),
                });
            }
            header = Some(values);
            continue;
        };
        if columns.len() != values.len() {
            return Err(TrancheError::RowLength {
                file: file.to_string(),
                header: columns.len(),
                values: values.len(),
                line: line.to_string(),
            });
        }
        let bindings: Vec<(String, String)> = columns.iter().cloned().zip(values).collect();
        tranches.push(tranche(&bindings)?);
    }

    // `tranches.sort(TRUTH_SENSITIVITY_ORDER)`, a stable sort on one column.
    tranches.sort_by(|left, right| {
        left.target_truth_sensitivity
            .total_cmp(&right.target_truth_sensitivity)
    });
    Ok(tranches)
}

/// One row, read in the order the reference's argument list evaluates.
fn tranche(bindings: &[(String, String)]) -> Result<TruthSensitivityTranche, TrancheError> {
    let target_truth_sensitivity = required_double(bindings, "targetTruthSensitivity")?;
    let min_vqslod = required_double(bindings, "minVQSLod")?;
    // The "optional" one, whose -1 the constructor then refuses.
    let num_known = optional_integer(bindings, "numKnown", -1)? as i64;
    let known_titv = optional_double(bindings, "knownTiTv", -1.0)?;
    // A long field, parsed as an int.
    let num_novel = required_integer(bindings, "numNovel")? as i64;
    let novel_titv = required_double(bindings, "novelTiTv")?;
    let accessible_truth_sites = optional_integer(bindings, "accessibleTruthSites", -1)?;
    let calls_at_truth_sites = optional_integer(bindings, "callsAtTruthSites", -1)?;
    let model = Mode::value_of(lookup(bindings, "model"))?;
    let name = lookup(bindings, "filterName")
        .unwrap_or_default()
        .to_string();

    // `super(...)` runs before the subclass's own check, so a negative count is reported first.
    if num_known < 0 || num_novel < 0 {
        return Err(TrancheError::NegativeCounts {
            name,
            known: num_known,
            novel: num_novel,
        });
    }
    if !(0.0..=100.0).contains(&target_truth_sensitivity) {
        return Err(TrancheError::UnreasonableTargetFdr(
            target_truth_sensitivity,
        ));
    }

    Ok(TruthSensitivityTranche {
        target_truth_sensitivity,
        min_vqslod,
        num_known,
        known_titv,
        num_novel,
        novel_titv,
        accessible_truth_sites,
        calls_at_truth_sites,
        model,
        name,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const COLUMNS: &str = "targetTruthSensitivity,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,model,accessibleTruthSites,callsAtTruthSites,truthSensitivity";
    const ROW_99: &str =
        "99.00,20,9,2.0000,1.8000,1.5000,VQSRTrancheSNP90.00to99.00,SNP,100,99,0.9900";
    const ROW_90: &str =
        "90.00,10,5,2.1000,1.9000,3.5000,VQSRTrancheSNP0.00to90.00,SNP,100,90,0.9000";

    fn file(body: &str) -> String {
        format!("# Variant quality score tranches file\n# Version number 5\n{body}")
    }

    #[test]
    fn the_rows_come_back_sorted_by_target_truth_sensitivity() {
        let text = file(&format!("{COLUMNS}\n{ROW_99}\n{ROW_90}\n"));
        let tranches = read_tranches("tranches", &text).expect("a good file");
        assert_eq!(
            tranches
                .iter()
                .map(|tranche| tranche.target_truth_sensitivity)
                .collect::<Vec<_>>(),
            vec![90.0, 99.0]
        );
        assert_eq!(tranches[0].name, "VQSRTrancheSNP0.00to90.00");
        assert_eq!(tranches[1].model, Mode::Snp);
    }

    #[test]
    fn the_columns_may_be_in_any_order_and_must_be_eleven() {
        // The header names them, so a permutation reads the same tranche.
        let permuted_columns = "model,targetTruthSensitivity,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,accessibleTruthSites,callsAtTruthSites,truthSensitivity";
        let permuted_row =
            "SNP,99.00,20,9,2.0000,1.8000,1.5000,VQSRTrancheSNP90.00to99.00,100,99,0.9900";
        let permuted = read_tranches("f", &file(&format!("{permuted_columns}\n{permuted_row}\n")))
            .expect("the same eleven columns");
        let plain =
            read_tranches("f", &file(&format!("{COLUMNS}\n{ROW_99}\n"))).expect("a good file");
        assert_eq!(permuted, plain);

        let short = read_tranches("f.tranches", &file("a,b,c\n")).expect_err("three columns");
        assert_eq!(
            short.message(),
            "File f.tranches is malformed: Expected 11 elements in header line a,b,c"
        );
    }

    #[test]
    fn the_optional_num_known_is_not_optional() {
        let text = file(&format!(
            "{}\n{ROW_99}\n",
            COLUMNS.replace("numKnown", "numknown")
        ));
        let error = read_tranches("f", &text).expect_err("the default is -1 and -1 is refused");
        assert_eq!(
            error.class(),
            "org.broadinstitute.hellbender.exceptions.GATKException"
        );
        assert_eq!(
            error.message(),
            "Invalid tranche VQSRTrancheSNP90.00to99.00 - no. variants is < 0 : known -1 novel 9"
        );
        // The other two defaults of -1 are never checked, and make the sensitivity 0.0.
        let unchecked = COLUMNS.replace("accessibleTruthSites", "unread");
        let tranches = read_tranches("f", &file(&format!("{unchecked}\n{ROW_99}\n")))
            .expect("nothing checks this one");
        assert_eq!(tranches[0].accessible_truth_sites, -1);
        assert_eq!(tranches[0].truth_sensitivity(), 0.0);
    }

    #[test]
    fn the_missing_key_refusal_names_no_file() {
        let text = file(&format!(
            "{}\n{ROW_99}\n",
            COLUMNS.replace("minVQSLod", "minVQSLOD")
        ));
        let error = read_tranches("f.tranches", &text).expect_err("a required key");
        assert_eq!(
            error.message(),
            "Unknown file is malformed: Malformed tranches file.  Missing required key minVQSLod"
        );
        // And the length refusal, from the same reader, does name it.
        let short_row = file(&format!("{COLUMNS}\n99.00,20,9\n"));
        let error = read_tranches("f.tranches", &short_row).expect_err("three of eleven");
        assert!(error.message().starts_with(
            "File f.tranches is malformed: Line had too few/many fields.  Header = 11 vals 3."
        ));
    }

    #[test]
    fn a_num_novel_past_an_int_cannot_be_read_back() {
        let text = file(&format!(
            "{COLUMNS}\n{}\n",
            ROW_99.replace(",20,9,", ",20,3000000000,")
        ));
        let error = read_tranches("f", &text).expect_err("a long parsed as an int");
        assert_eq!(
            error.message(),
            "Unknown file is malformed: Malformed tranches file. Invalid value for key numNovel"
        );
        // A value that is not a number at all is the same refusal.
        let text = file(&format!(
            "{COLUMNS}\n{}\n",
            ROW_99.replace(",20,9,", ",20,many,")
        ));
        assert_eq!(
            read_tranches("f", &text)
                .expect_err("not a number")
                .message(),
            "Unknown file is malformed: Malformed tranches file. Invalid value for key numNovel"
        );
    }

    #[test]
    fn the_model_is_read_without_a_guard() {
        let text = file(&format!(
            "{}\n{ROW_99}\n",
            COLUMNS.replace(",model,", ",mdl,")
        ));
        let error = read_tranches("f", &text).expect_err("valueOf(null)");
        assert_eq!(error.class(), "java.lang.NullPointerException");
        assert_eq!(error.message(), "Name is null");

        let text = file(&format!(
            "{COLUMNS}\n{}\n",
            ROW_99.replace(",SNP,", ",GERMLINE,")
        ));
        let error = read_tranches("f", &text).expect_err("no such constant");
        assert_eq!(error.class(), "java.lang.IllegalArgumentException");
        assert!(error.message().ends_with("Mode.GERMLINE"));
    }

    #[test]
    fn a_sensitivity_outside_the_range_is_refused() {
        let text = file(&format!(
            "{COLUMNS}\n{}\n",
            ROW_99.replace("99.00,20", "150.00,20")
        ));
        let error = read_tranches("f", &text).expect_err("over a hundred");
        assert_eq!(error.message(), "Target FDR is unreasonable 150.0");
    }
}

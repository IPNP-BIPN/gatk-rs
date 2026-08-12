//! `CompareBaseQualities` and its `CompareMatrix`, ported from
//! `org.broadinstitute.hellbender.tools.validation` (GATK 4.6.2.0).
//!
//! The first tool here that reads **two** BAMs at once, and the first that is not a walker: it is a
//! `PicardCommandLineProgram` with two positional arguments. It walks both files in lockstep,
//! counts every pair of base qualities into a 94x94 matrix, and prints that matrix twice: once as it
//! stands and once through the static quantization mapping.
//!
//! # No filters, and a strict reader
//!
//! There are no read filters at all, so a duplicate and a vendor failure are counted like anything
//! else. What the tool does skip is a **secondary or supplementary** read, through htsjdk's
//! `SecondaryOrSupplementarySkippingIterator`, and it skips them in each file independently: the two
//! line up on their primary reads and not on their record counts.
//!
//! Its reader is `VALIDATION_STRINGENCY` STRICT rather than the engine's SILENT, so a record that
//! contradicts itself stops the run where a walker would have carried on. That belongs to the
//! reader rather than to this port, which takes the records it is given.
//!
//! # The summary is the matrix collapsed onto its diagonals
//!
//! ```java
//! deltaColumns[(dimension-1) + i-j] += matrix[i][j];
//! ```
//!
//! When every count lands on the main diagonal the summary is one sentence; otherwise it is a table
//! of `diff`, `count` and a percentage formatted `%.4f`. The full matrix that follows prints only
//! non-zero entries, in row-major order, with `diff` = QRead1 - QRead2, so swapping the two inputs
//! flips every sign.
//!
//! The binned half is the same counts through `--static-quantized-quals`, and the two halves can
//! disagree with each other: quantization can map every difference onto one bin, so the same run
//! reports four differences above and "all 8 quality scores are the same" below.

use htsjdk_bam::record::BamRecord;

use gatk_engine::read;

/// `CompareMatrix.dimension`, which is `QualityUtils.MAX_SAM_QUAL_SCORE + 1`.
pub const DIMENSION: usize = 94;

/// What this tool refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CompareError {
    /// Two reads at the same position with different names.
    OutOfOrder(String, String),
    /// One file ran out before the other.
    DifferentCounts,
    /// Two reads whose quality arrays are different lengths.
    RaggedQualities(usize, usize),
    /// `--round-down-quantized` without `--static-quantized-quals`.
    RoundDownAlone,
    /// `--throw-on-diff` and the matrix has an off-diagonal entry.
    QualitiesDiffer,
}

impl CompareError {
    /// The message the reference carries, without any exception prefix.
    pub fn message(&self) -> String {
        match self {
            CompareError::OutOfOrder(first, second) => {
                format!("files do not have the same exact order of reads:{first} vs {second}")
            }
            CompareError::DifferentCounts => {
                "files do not have the same exact number of reads".to_string()
            }
            CompareError::RaggedQualities(first, second) => format!(
                "The length of the quality scores are not the same for read {first},{second}"
            ),
            CompareError::RoundDownAlone => {
                "Argument round-down-quantized has a bad value: true. This option can only be used \
                 if static-quantized-quals is also used."
                    .to_string()
            }
            CompareError::QualitiesDiffer => {
                "Quality scores from the two BAMs do not match".to_string()
            }
        }
    }
}

/// `CompareMatrix`: the raw counts and the binned ones, filled together.
pub struct CompareMatrix {
    bin: Vec<u8>,
    matrix: Vec<Vec<i64>>,
    binned: Vec<Vec<i64>>,
}

impl CompareMatrix {
    /// `new CompareMatrix(binning)`.
    ///
    /// The mapping is cloned, as the reference clones it: a caller that goes on to change its own
    /// array does not change what this matrix bins with.
    pub fn new(binning: &[u8]) -> CompareMatrix {
        CompareMatrix {
            bin: binning.to_vec(),
            matrix: vec![vec![0; DIMENSION]; DIMENSION],
            binned: vec![vec![0; DIMENSION]; DIMENSION],
        }
    }

    /// `add(first, second)`: one pair of reads, quality by quality.
    pub fn add(&mut self, first: &[u8], second: &[u8]) -> Result<(), CompareError> {
        if first.len() != second.len() {
            return Err(CompareError::RaggedQualities(first.len(), second.len()));
        }
        for (&a, &b) in first.iter().zip(second.iter()) {
            self.matrix[a as usize][b as usize] += 1;
            self.binned[self.bin[a as usize] as usize][self.bin[b as usize] as usize] += 1;
        }
        Ok(())
    }

    /// `hasNonDiagonalElements`, which is what decides the tool's exit code.
    pub fn has_non_diagonal_elements(&self) -> bool {
        (0..DIMENSION).any(|i| (0..DIMENSION).any(|j| i != j && self.matrix[i][j] != 0))
    }

    /// `printOutput`: the summary and the matrix, raw then binned.
    pub fn to_text(&self) -> String {
        let mut out = String::new();
        out.push_str(&summary(&self.matrix, "CompareMatrix"));
        out.push('\n');
        out.push_str(&full_matrix(&self.matrix, "CompareMatrix"));
        out.push_str(&summary(&self.binned, "CompareMatrix-binned"));
        out.push('\n');
        out.push_str(&full_matrix(&self.binned, "CompareMatrix-binned"));
        out
    }
}

/// `printSummary`.
fn summary(matrix: &[Vec<i64>], name: &str) -> String {
    let total_size = 2 * (DIMENSION - 1) + 1;
    let mut deltas = vec![0i64; total_size];
    for (i, row) in matrix.iter().enumerate() {
        for (j, count) in row.iter().enumerate() {
            deltas[(DIMENSION - 1) + i - j] += count;
        }
    }

    let mut out = format!("-----------{name} summary------------\n");
    let sum: i64 = deltas.iter().sum();
    if sum == deltas[DIMENSION - 1] {
        out.push_str(&format!("all {sum} quality scores are the same\n"));
        return out;
    }
    out.push_str("diff\tcount\t%total\n");
    for (k, &count) in deltas.iter().enumerate() {
        if count != 0 {
            // `String.format("%.4f")`, which rounds HALF_UP, not `DecimalFormat`'s HALF_EVEN.
            out.push_str(&format!(
                "{}\t{}\t{}\n",
                k as i64 - (DIMENSION as i64 - 1),
                count,
                gatk_engine::java_format::format_decimals(count as f64 * 100.0 / sum as f64, 4),
            ));
        }
    }
    out
}

/// `print`: the non-zero entries, in row-major order.
fn full_matrix(matrix: &[Vec<i64>], name: &str) -> String {
    let mut out = format!("---------{name} full matrix (non-zero entries) ----------\n");
    out.push_str("QRead1\tQRead2\tdiff\tcount\n");
    for (i, row) in matrix.iter().enumerate() {
        for (j, &count) in row.iter().enumerate() {
            if count != 0 {
                out.push_str(&format!("{i}\t{j}\t{}\t{count}\n", i as i64 - j as i64));
            }
        }
    }
    out
}

/// The tool's own arguments.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct CompareArguments {
    /// `--static-quantized-quals`. The mapping is built from it, and the list is sorted in place.
    pub static_quantization_quals: Vec<i32>,
    /// `--round-down-quantized`, which only means anything beside the list above.
    pub round_down: bool,
    /// `--throw-on-diff`.
    pub throw_on_diff: bool,
}

/// What one run produced: the report, and the value the tool returns as its exit code.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CompareResult {
    pub report: String,
    /// `finalMatrix.hasNonDiagonalElements() ? 1 : 0`.
    pub exit_code: i32,
}

/// `CompareBaseQualities.doWork`.
///
/// The two record lists are what the readers handed over; the secondary and supplementary skipping
/// happens here, in each list independently, which is what `SecondaryOrSupplementarySkippingIterator`
/// does to each file.
pub fn compare_base_qualities(
    first: &[BamRecord],
    second: &[BamRecord],
    arguments: &CompareArguments,
) -> Result<CompareResult, CompareError> {
    if arguments.round_down && arguments.static_quantization_quals.is_empty() {
        return Err(CompareError::RoundDownAlone);
    }
    let mut quals = arguments.static_quantization_quals.clone();
    let mapping = gatk_engine::bqsr_transformer::construct_static_quantized_mapping(
        &mut quals,
        arguments.round_down,
    );

    let mut matrix = CompareMatrix::new(&mapping);
    let mut left = first
        .iter()
        .filter(|record| !is_secondary_or_supplementary(record));
    let mut right = second
        .iter()
        .filter(|record| !is_secondary_or_supplementary(record));

    let mut a = left.next();
    let mut b = right.next();
    while let (Some(one), Some(other)) = (a, b) {
        if one.read_name != other.read_name {
            return Err(CompareError::OutOfOrder(
                one.read_name.clone(),
                other.read_name.clone(),
            ));
        }
        matrix.add(&one.base_qualities, &other.base_qualities)?;
        a = left.next();
        b = right.next();
    }
    if a.is_some() || b.is_some() {
        return Err(CompareError::DifferentCounts);
    }

    let differs = matrix.has_non_diagonal_elements();
    if arguments.throw_on_diff && differs {
        return Err(CompareError::QualitiesDiffer);
    }
    Ok(CompareResult {
        report: matrix.to_text(),
        exit_code: if differs { 1 } else { 0 },
    })
}

/// `SAMRecord.isSecondaryOrSupplementary`.
fn is_secondary_or_supplementary(record: &BamRecord) -> bool {
    record.flags & (read::flags::NOT_PRIMARY_ALIGNMENT | read::flags::SUPPLEMENTARY_ALIGNMENT) != 0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(name: &str, quals: &[u8], flags: u16) -> BamRecord {
        BamRecord {
            read_name: name.to_string(),
            flags,
            base_qualities: quals.to_vec(),
            ..Default::default()
        }
    }

    #[test]
    fn identical_files_say_so_in_one_sentence() {
        let a = vec![read("r1", &[30, 30], 0)];
        let b = vec![read("r1", &[30, 30], 0)];
        let result = compare_base_qualities(&a, &b, &CompareArguments::default()).expect("it runs");
        assert!(result.report.contains("all 2 quality scores are the same"));
        assert_eq!(result.exit_code, 0);
    }

    #[test]
    fn the_diff_flips_sign_when_the_inputs_are_swapped() {
        let a = vec![read("r1", &[31], 0)];
        let b = vec![read("r1", &[30], 0)];
        let forward = compare_base_qualities(&a, &b, &CompareArguments::default()).expect("runs");
        let reverse = compare_base_qualities(&b, &a, &CompareArguments::default()).expect("runs");
        assert!(forward.report.contains("31\t30\t1\t1"));
        assert!(reverse.report.contains("30\t31\t-1\t1"));
        assert_eq!(forward.exit_code, 1);
    }

    #[test]
    fn a_secondary_read_is_skipped_in_the_file_that_has_it() {
        let a = vec![read("r1", &[30], 0), read("r2", &[20], 0)];
        let b = vec![
            read("r1", &[30], 0),
            read("sec", &[40], 0x100),
            read("r2", &[20], 0),
        ];
        let result = compare_base_qualities(&a, &b, &CompareArguments::default()).expect("it runs");
        assert_eq!(result.exit_code, 0);
    }

    #[test]
    fn round_down_alone_is_refused_before_anything_is_read() {
        let arguments = CompareArguments {
            round_down: true,
            ..CompareArguments::default()
        };
        assert_eq!(
            compare_base_qualities(&[], &[], &arguments).unwrap_err(),
            CompareError::RoundDownAlone
        );
    }

    #[test]
    fn ragged_qualities_name_the_two_lengths_and_not_the_read() {
        let a = vec![read("r1", &[30, 30, 30, 30], 0)];
        let b = vec![read("r1", &[30, 30, 30], 0)];
        let error = compare_base_qualities(&a, &b, &CompareArguments::default()).unwrap_err();
        assert_eq!(
            error.message(),
            "The length of the quality scores are not the same for read 4,3"
        );
    }
}

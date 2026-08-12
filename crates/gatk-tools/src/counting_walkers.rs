//! `CountBases`, `CountReads` and `FlagStat`, ported from
//! `org.broadinstitute.hellbender.tools` (GATK 4.6.2.0).
//!
//! The first three tools here whose output is a **number** rather than a BAM. All three are
//! `ReadWalker`s whose `apply` is one line, so what they measure is the traversal and the
//! formatting rather than any transform. They take the engine's default filter,
//! `WellformedReadFilter`, which is what decides whether a malformed read is counted at all.
//!
//! # The percentages are computed in `float` and formatted `#0.00`
//!
//! ```java
//! NumberFormat percentFormatter = new DecimalFormat("#0.00");
//! ...
//! percentFormatter.format(((float) mapped / (float) readCount) * 100.0)
//! ```
//!
//! Two things there are not what a reader would assume. The division is in **`float`**, so
//! `11f/12f` is `0.9166667` and not the double `0.9166666666666666`; and `DecimalFormat` rounds
//! **HALF_EVEN**, where `String.format("%.2f")` rounds HALF_UP. The two agree on this corpus and
//! would not on a value landing exactly on a half.
//!
//! With no reads at all the ratio is `0f/0f`, so the percentage is `NaN` and the line reads
//! `0 mapped (NaN%)`. That is a real output, not a failure.
//!
//! # And `read2` is tested before `read1`
//!
//! ```java
//! if ( read.isSecondOfPair() ) { this.read2++; }
//! else if ( read.isFirstOfPair() ) { this.read1++; }
//! ```
//!
//! A read carrying both 0x40 and 0x80 counts as `read2` only. Both tests also go through
//! `isPaired()` first, because `isFirstOfPair` is `isPaired() && flag` and not the bare flag.

use gatk_engine::read;
use htsjdk_bam::record::BamRecord;

/// `CountBases`: the number of **bases**, which is `read.getLength()` and not the span.
///
/// A read carrying a deletion covers more reference than it has bases and one carrying an insertion
/// covers less; neither changes what this counts.
pub fn count_bases(records: &[BamRecord]) -> i64 {
    records
        .iter()
        .map(|record| record.read_bases.len() as i64)
        .sum()
}

/// `CountReads`: how many reads the traversal kept.
pub fn count_reads(records: &[BamRecord]) -> i64 {
    records.len() as i64
}

/// `FlagStat.FlagStatus`: the twelve counters, in the order the output prints them.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct FlagStatus {
    pub read_count: i64,
    pub qc_failure: i64,
    pub duplicates: i64,
    pub mapped: i64,
    pub paired_in_sequencing: i64,
    pub read1: i64,
    pub read2: i64,
    pub properly_paired: i64,
    pub with_itself_and_mate_mapped: i64,
    pub singletons: i64,
    pub with_mate_mapped_to_a_different_chr: i64,
    pub with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5: i64,
}

impl FlagStatus {
    /// `add(read)`.
    ///
    /// Every paired counter is inside the `isPaired()` branch, so an unpaired read contributes to
    /// `readCount`, `QC_failure`, `duplicates` and `mapped` and to nothing else.
    pub fn add(&mut self, record: &BamRecord, contig: Option<&str>, mate_contig: Option<&str>) {
        self.read_count += 1;

        if read::fails_vendor_quality_check(record) {
            self.qc_failure += 1;
        }
        if read::is_duplicate(record) {
            self.duplicates += 1;
        }
        // `isUnmapped` is the three-part test, not the 0x4 flag: see `gatk_engine::read`.
        if !read::is_unmapped(record) {
            self.mapped += 1;
        }
        if read::is_paired(record) {
            self.paired_in_sequencing += 1;

            // Second is tested first, and the two are exclusive, so a read with both flags is read2.
            if read::is_second_of_pair(record) {
                self.read2 += 1;
            } else if read::is_first_of_pair(record) {
                self.read1 += 1;
            }

            if read::is_proper_pair(record) {
                self.properly_paired += 1;
            }

            if !read::is_unmapped(record) && !read::mate_is_unmapped(record) {
                self.with_itself_and_mate_mapped += 1;

                // `read.getContig().equals(read.getMateContig())`, both resolved through the header.
                if contig != mate_contig {
                    self.with_mate_mapped_to_a_different_chr += 1;
                    if record.mapping_quality as i32 >= 5 {
                        self.with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5 += 1;
                    }
                }
            }

            if !read::is_unmapped(record) && read::mate_is_unmapped(record) {
                self.singletons += 1;
            }
        }
    }

    /// `merge(that)`: the counters added pairwise, which is what the Spark version gathers with.
    pub fn merge(&mut self, other: &FlagStatus) {
        self.read_count += other.read_count;
        self.qc_failure += other.qc_failure;
        self.duplicates += other.duplicates;
        self.mapped += other.mapped;
        self.paired_in_sequencing += other.paired_in_sequencing;
        self.read1 += other.read1;
        self.read2 += other.read2;
        self.properly_paired += other.properly_paired;
        self.with_itself_and_mate_mapped += other.with_itself_and_mate_mapped;
        self.singletons += other.singletons;
        self.with_mate_mapped_to_a_different_chr += other.with_mate_mapped_to_a_different_chr;
        self.with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5 +=
            other.with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5;
    }

    /// `toString()`: thirteen lines, the last without a trailing newline.
    pub fn to_text(&self) -> String {
        let percent = |numerator: i64| {
            format_percent((numerator as f32 / self.read_count as f32) as f64 * 100.0)
        };
        format!(
            "{} in total\n\
             {} QC failure\n\
             {} duplicates\n\
             {} mapped ({}%)\n\
             {} paired in sequencing\n\
             {} read1\n\
             {} read2\n\
             {} properly paired ({}%)\n\
             {} with itself and mate mapped\n\
             {} singletons ({}%)\n\
             {} with mate mapped to a different chr\n\
             {} with mate mapped to a different chr (mapQ>=5)",
            self.read_count,
            self.qc_failure,
            self.duplicates,
            self.mapped,
            percent(self.mapped),
            self.paired_in_sequencing,
            self.read1,
            self.read2,
            self.properly_paired,
            percent(self.properly_paired),
            self.with_itself_and_mate_mapped,
            self.singletons,
            percent(self.singletons),
            self.with_mate_mapped_to_a_different_chr,
            self.with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5,
        )
    }
}

/// `new DecimalFormat("#0.00").format(value)`.
///
/// **HALF_EVEN**, which is `DecimalFormat`'s default and not `String.format`'s HALF_UP, and `NaN`
/// for a value that is not a number, which is what an empty file produces.
pub fn format_percent(value: f64) -> String {
    if value.is_nan() {
        // `DecimalFormatSymbols.getNaN()`, which in the root locale is the three letters.
        return "NaN".to_string();
    }
    if value.is_infinite() {
        return if value > 0.0 {
            "∞".to_string()
        } else {
            "-∞".to_string()
        };
    }
    // Two decimal places, ties to even. `f64::round_ties_even` on the scaled value is that rule.
    let scaled = value * 100.0;
    let rounded = round_half_even(scaled) / 100.0;
    // `#0.00` always shows two decimals and at least one integer digit.
    format!("{rounded:.2}")
}

/// Ties to even, which is what `RoundingMode.HALF_EVEN` means.
fn round_half_even(value: f64) -> f64 {
    let floor = value.floor();
    let difference = value - floor;
    match difference.partial_cmp(&0.5) {
        Some(std::cmp::Ordering::Greater) => floor + 1.0,
        Some(std::cmp::Ordering::Less) => floor,
        // An exact tie is the only place this differs from HALF_UP: it goes to the even side.
        _ if (floor as i64) % 2 == 0 => floor,
        _ => floor + 1.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_reads_is_a_nan_percentage_and_not_a_failure() {
        let status = FlagStatus::default();
        let text = status.to_text();
        assert!(text.starts_with("0 in total\n"), "{text}");
        assert!(text.contains("0 mapped (NaN%)"), "{text}");
        // Thirteen lines, and the last carries no newline.
        assert_eq!(text.lines().count(), 12);
        assert!(!text.ends_with('\n'));
    }

    #[test]
    fn the_percentage_is_computed_in_float() {
        let status = FlagStatus {
            read_count: 12,
            mapped: 11,
            ..Default::default()
        };
        // `11f/12f` is 0.9166667, which times a hundred and rounded to two places is 91.67.
        assert!(status.to_text().contains("11 mapped (91.67%)"));
    }

    #[test]
    fn the_rounding_is_half_even_and_not_half_up() {
        // 0.125 to two places: HALF_UP gives 0.13, HALF_EVEN gives 0.12.
        assert_eq!(format_percent(0.125), "0.12");
        assert_eq!(format_percent(0.135), "0.14");
        assert_eq!(format_percent(100.0), "100.00");
        assert_eq!(format_percent(0.0), "0.00");
    }
}

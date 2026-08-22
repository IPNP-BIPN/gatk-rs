//! `MTLowHeteroplasmyFilterTool`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.mutect.filtering.MTLowHeteroplasmyFilterTool`
//! (GATK 4.6.2.0).
//!
//! Low-heteroplasmy mitochondrial calls filtered, but only once there are too many of them.
//!
//! # `--low-het-threshold` does nothing
//!
//! ```java
//! private final double lowHetThreshold = 0.1;
//! ```
//!
//! `private final` with a constant initialiser is a compile-time constant, so javac folds every
//! read of it into the literal and whatever Barclay writes into the field is never looked at
//! again. A run at `--low-het-threshold 0.6` over fractions of 0.05 and 0.5 filters only the 0.05
//! ones, which is what the default does. The control is `maxAllowedLowHets`, declared without
//! `final` in the same class, which works.
//!
//! So this port takes the allowance as an argument and the threshold as a constant, because that
//! is what the reference does whatever its command line says.
//!
//! # The filter is all or nothing across the whole file
//!
//! The first pass counts the UNFILTERED low sites; the second filters every low allele it can
//! find, but only if that count exceeded the allowance. A file with three such sites comes out
//! untouched and a file with four comes out with all four filtered, so one record's fate is
//! decided by the others.
//!
//! A site that is already filtered does not count. `PASS` does: htsjdk's `isFiltered` asks whether
//! the filter SET is non-empty, and `PASS` leaves it empty.
//!
//! # `AF=.` is an absent attribute, not a missing value
//!
//! It takes the `() -> null` default and the first pass throws a `NullPointerException`, exactly as
//! it does for a genotype with no `AF` field at all, because the first pass reads EVERY genotype
//! without a precondition while the second reads only those that have one. The `Double.MAX_VALUE`
//! substitution is for a missing entry WITHIN a multi-valued array: `AF=0.05,.` filters its first
//! alternate and never its second.

use crate::numt_filter::{merged_as_filter_string, NuMTError, Record as FilterRecord};

/// `GATKVCFConstants.LOW_HET_FILTER_NAME`.
pub const FILTER_NAME: &str = "mt_many_low_hets";

/// The threshold, which is a compile-time constant in the reference and therefore one here.
pub const LOW_HET_THRESHOLD: f64 = 0.1;

/// The one argument that works.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Arguments {
    pub max_allowed_low_hets: i32,
}

impl Default for Arguments {
    fn default() -> Self {
        Arguments {
            max_allowed_low_hets: 3,
        }
    }
}

/// One record, reduced to what this filter reads and writes.
#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    /// The alternate alleles, reference excluded: `getAltDataByAllele` keys on those alone.
    pub alternates: Vec<String>,
    /// Per genotype, the `AF` list. `None` is both an absent field and a bare `.`, which the
    /// reference cannot tell apart and throws on either way.
    pub allele_fractions: Vec<Option<Vec<Option<f64>>>>,
    pub filters: Vec<String>,
    pub as_filter_status: Option<String>,
}

/// What the run refuses.
#[derive(Debug, Clone, PartialEq)]
pub enum LowHetError {
    /// The first pass asked a genotype for an allele fraction it does not have.
    NoAlleleFraction,
    /// The merge into `AS_FilterStatus` found nothing to merge into.
    Filter(NuMTError),
}

impl LowHetError {
    pub fn java_class(&self) -> &str {
        match self {
            LowHetError::NoAlleleFraction => "java.lang.NullPointerException",
            LowHetError::Filter(error) => error.java_class(),
        }
    }

    pub fn message(&self) -> String {
        match self {
            LowHetError::NoAlleleFraction => {
                "Cannot read the array length because \"array\" is null".to_string()
            }
            LowHetError::Filter(error) => error.message(),
        }
    }
}

/// `VariantContextGetters.getAttributeAsDoubleArray(g, AF, () -> null, Double.MAX_VALUE)`: a
/// missing entry inside the array becomes `Double.MAX_VALUE`, and an absent attribute is null.
fn fractions(genotype: &Option<Vec<Option<f64>>>) -> Result<Vec<f64>, LowHetError> {
    match genotype {
        None => Err(LowHetError::NoAlleleFraction),
        Some(values) => Ok(values
            .iter()
            .map(|value| value.unwrap_or(f64::MAX))
            .collect()),
    }
}

/// `isSiteLowHeteroplasmy`, which reads EVERY genotype and therefore throws on one without `AF`.
pub fn is_site_low_heteroplasmy(record: &Record) -> Result<bool, LowHetError> {
    for genotype in &record.allele_fractions {
        if fractions(genotype)?
            .iter()
            .any(|fraction| *fraction < LOW_HET_THRESHOLD)
        {
            return Ok(true);
        }
    }
    Ok(false)
}

/// `variant.isNotFiltered()`: the filter SET being empty, which `PASS` also is.
pub fn is_not_filtered(record: &Record) -> bool {
    record.filters.is_empty() || record.filters == ["PASS"]
}

/// `areAllelesArtifacts`: per alternate allele, the maximum across the genotypes that HAVE an
/// allele fraction, compared strictly against the threshold.
pub fn alleles_are_artifacts(record: &Record) -> Vec<bool> {
    let mut by_alternate: Vec<Vec<f64>> = vec![Vec::new(); record.alternates.len()];
    for genotype in record.allele_fractions.iter().flatten() {
        let values: Vec<f64> = genotype
            .iter()
            .map(|value| value.unwrap_or(f64::MAX))
            .collect();
        for (slot, fraction) in by_alternate.iter_mut().zip(values) {
            slot.push(fraction);
        }
    }
    by_alternate
        .iter()
        .map(|fractions| {
            fractions
                .iter()
                .copied()
                .fold(None::<f64>, |best, value| {
                    Some(match best {
                        None => value,
                        Some(best) => best.max(value),
                    })
                })
                .unwrap_or(0.0)
                < LOW_HET_THRESHOLD
        })
        .collect()
}

/// Both passes over a whole file: count, decide, then filter.
pub fn run(records: &[Record], arguments: &Arguments) -> Result<Vec<Record>, LowHetError> {
    let mut unfiltered_low_hets = 0;
    for record in records {
        if is_not_filtered(record) && is_site_low_heteroplasmy(record)? {
            unfiltered_low_hets += 1;
        }
    }
    let failed_low_het = unfiltered_low_hets > arguments.max_allowed_low_hets;

    let mut written = Vec::with_capacity(records.len());
    for record in records {
        let mut out = record.clone();
        if failed_low_het {
            let applied = alleles_are_artifacts(record);
            if !applied.contains(&false) {
                // A record whose filter column was PASS loses it, since the filter set gains a
                // member and PASS is only ever the empty set.
                out.filters.retain(|filter| filter != "PASS");
                out.filters.push(FILTER_NAME.to_string());
            }
            if applied.contains(&true) {
                let filter_record = FilterRecord {
                    alleles: std::iter::once(String::new())
                        .chain(record.alternates.iter().cloned())
                        .collect(),
                    allele_depths: Vec::new(),
                    filters: record.filters.clone(),
                    as_filter_status: record.as_filter_status.clone(),
                };
                out.as_filter_status = Some(
                    merged_as_filter_string(&filter_record, &applied, FILTER_NAME)
                        .map_err(LowHetError::Filter)?,
                );
            }
        }
        written.push(out);
    }
    Ok(written)
}

//! `FilterVcf`, ported from `picard.vcf.filter.FilterVcf` and the filters under it
//! (Picard 3.4.0).
//!
//! Three site filters and two genotype filters applied to every record, written back out with the
//! FILTER column and the genotypes' `FT` set.
//!
//! # A passing genotype is filtered too, and it does not show
//!
//! ```java
//! if (filters == null || filters.isEmpty()) { gtBuilder.filter(PASS_FILTER); }
//! ```
//!
//! Every genotype comes out of the iterator carrying an `FT`, and the writer drops a `PASS` one. So
//! the FORMAT column gains `FT` only on the records where something was filtered, and a port that
//! set `FT` only on filtered genotypes would write the same bytes for the wrong reason.
//!
//! # `AllGtsFiltered` is a site filter the genotypes set
//!
//! A record every one of whose genotypes was filtered is itself filtered, and that replaces the
//! `PASS` the three site filters would otherwise have left.
//!
//! # Two thresholds are tested against a default that means "absent"
//!
//! `QD` defaults to -1 and the test is `qd >= 0 && qd < minimum`, so a record with no QD passes
//! however high the threshold. `FS` defaults to 0 and the test is `fs > max`, so a record with no
//! FS is filtered only when the threshold is negative. `getGQ()` and `getDP()` are -1 when absent,
//! which is BELOW any non-negative threshold, so a genotype carrying neither is filtered by both.
//!
//! # The allele balance filter groups by the genotype's alleles
//!
//! Two samples with the same het call share one tally, and the filter fires on the tally rather
//! than on either sample. A het with no `AD` is skipped, and a record with no het genotype answers
//! nothing at all.

use htsjdk_vcf::encoder::EncodeError;
use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use htsjdk_vcf::reader::read_vcf;
use htsjdk_vcf::variant::{Genotype, Value, VariantContext};
use htsjdk_vcf::vcf_file::write_vcf;

/// The tool's defaults, which filter nothing: every threshold is at the edge of its range, and
/// `MAX_FS` is `Double.MAX_VALUE`. A run with no arguments rewrites the file with `FT` set on the
/// genotypes and changes no verdict.
pub const DEFAULT_MIN_AB: f64 = 0.0;
pub const DEFAULT_MAX_FS: f64 = f64::MAX;
pub const DEFAULT_MIN_QD: f64 = 0.0;
pub const DEFAULT_MIN_GQ: i32 = 0;
pub const DEFAULT_MIN_DP: i32 = 0;

/// The five filter names.
pub const ALL_GTS_FILTERED: &str = "AllGtsFiltered";
pub const ALLELE_BALANCE: &str = "AlleleBalance";
pub const STRAND_BIAS: &str = "StrandBias";
pub const LOW_QD: &str = "LowQD";
pub const LOW_GQ: &str = "LowGQ";
pub const LOW_DP: &str = "LowDP";

/// The thresholds, with the tool's defaults.
#[derive(Debug, Clone, PartialEq)]
pub struct Thresholds {
    pub min_ab: f64,
    pub max_fs: f64,
    pub min_qd: f64,
    pub min_gq: i32,
    pub min_dp: i32,
}

impl Default for Thresholds {
    fn default() -> Self {
        Thresholds {
            min_ab: DEFAULT_MIN_AB,
            max_fs: DEFAULT_MAX_FS,
            min_qd: DEFAULT_MIN_QD,
            min_gq: DEFAULT_MIN_GQ,
            min_dp: DEFAULT_MIN_DP,
        }
    }
}

/// What the tool refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FilterError {
    /// A `.vcf` output whose input header declares no contigs.
    NoSequenceDictionary,
    /// The reader or the writer refused.
    Vcf(String, String),
}

impl FilterError {
    pub fn java_class(&self) -> &str {
        match self {
            FilterError::NoSequenceDictionary => "picard.PicardException",
            FilterError::Vcf(class, _) => class,
        }
    }

    pub fn message(&self) -> String {
        match self {
            FilterError::NoSequenceDictionary => "The input vcf must have a sequence dictionary \
                 in order to create indexed vcf or bcfs."
                .to_string(),
            FilterError::Vcf(_, message) => message.clone(),
        }
    }
}

/// The four FILTER lines and the FT line the tool adds to the header it was given.
pub fn added_header_lines() -> Vec<HeaderLine> {
    let filter = |id: &str, description: &str| HeaderLine::Filter {
        id: id.to_string(),
        description: description.to_string(),
    };
    vec![
        filter(
            ALL_GTS_FILTERED,
            "Site filtered out because all genotypes are filtered out.",
        ),
        HeaderLine::Compound {
            key: "FORMAT".to_string(),
            id: "FT".to_string(),
            number: Cardinality::Unbounded,
            line_type: LineType::String,
            description: "Genotype filters.".to_string(),
            extra: Vec::new(),
        },
        filter(
            ALLELE_BALANCE,
            "Heterozygote allele balance below required threshold.",
        ),
        filter(
            STRAND_BIAS,
            "Site exhibits excessive allele/strand correlation.",
        ),
        filter(LOW_QD, "Site exhibits QD value below a hard limit."),
    ]
}

/// `getAttributeAsDouble(key, default)`, which parses a string and falls back on the default.
fn attribute_as_double(record: &VariantContext, key: &str, default: f64) -> f64 {
    record
        .attributes
        .iter()
        .find(|(name, _)| name == key)
        .and_then(|(_, value)| match value {
            Value::Double(number) => Some(*number),
            Value::Int(number) => Some(*number as f64),
            Value::Str(text) => text.parse().ok(),
            _ => None,
        })
        .unwrap_or(default)
}

/// `AlleleBalanceFilter.filter`.
fn allele_balance(record: &VariantContext, min_ab: f64) -> Option<String> {
    let hets: Vec<&Genotype> = record
        .genotypes
        .iter()
        .filter(|genotype| is_het(genotype))
        .collect();
    if hets.is_empty() {
        return None;
    }
    // `getGenotypesOrderedByName()` and a map keyed by the allele list.
    let mut ordered: Vec<&Genotype> = record.genotypes.iter().collect();
    ordered.sort_by(|left, right| left.sample_name.cmp(&right.sample_name));
    let mut tallies: Vec<(Vec<String>, i64, i64)> = Vec::new();
    for genotype in ordered {
        if !is_het(genotype) || genotype.ad.is_none() {
            continue;
        }
        let depths = genotype.ad.as_ref().expect("a genotype with AD");
        let alleles: Vec<String> = genotype
            .alleles
            .iter()
            .map(|allele| allele.base_string())
            .collect();
        let first = allele_index(record, &alleles[0]);
        let second = allele_index(record, &alleles[1]);
        let (Some(first), Some(second)) = (first, second) else {
            continue;
        };
        let entry = match tallies.iter_mut().find(|(key, _, _)| *key == alleles) {
            Some(entry) => entry,
            None => {
                tallies.push((alleles.clone(), 0, 0));
                tallies.last_mut().expect("just pushed")
            }
        };
        entry.1 += i64::from(*depths.get(first).unwrap_or(&0));
        entry.2 += i64::from(*depths.get(second).unwrap_or(&0));
    }
    for (_, first, second) in &tallies {
        let total = first + second;
        if total > 0 && (*first.min(second) as f64) / (total as f64) < min_ab {
            return Some(ALLELE_BALANCE.to_string());
        }
    }
    None
}

/// `Genotype.isHet()`: called, and two different alleles.
fn is_het(genotype: &Genotype) -> bool {
    genotype.alleles.len() == 2
        && !genotype.alleles.iter().any(|allele| allele.is_no_call())
        && genotype.alleles[0].base_string() != genotype.alleles[1].base_string()
}

fn allele_index(record: &VariantContext, bases: &str) -> Option<usize> {
    record
        .alleles
        .iter()
        .position(|allele| allele.base_string() == bases)
}

/// `FilterApplyingVariantIterator.next`, for one record.
pub fn filter_record(record: &VariantContext, thresholds: &Thresholds) -> VariantContext {
    let mut filters: Vec<String> = Vec::new();
    if let Some(name) = allele_balance(record, thresholds.min_ab) {
        filters.push(name);
    }
    // `getAttributeAsDouble("FS", 0)`, so an absent FS is zero.
    if attribute_as_double(record, "FS", 0.0) > thresholds.max_fs {
        filters.push(STRAND_BIAS.to_string());
    }
    // `getAttributeAsDouble("QD", -1)`, and a negative value means absent.
    let qd = attribute_as_double(record, "QD", -1.0);
    if qd >= 0.0 && qd < thresholds.min_qd {
        filters.push(LOW_QD.to_string());
    }

    let mut genotypes = Vec::new();
    let mut all_filtered = !record.genotypes.is_empty();
    for genotype in &record.genotypes {
        let mut own: Vec<String> = Vec::new();
        // `getGQ()` and `getDP()` are -1 when absent.
        if genotype.gq.unwrap_or(-1) < thresholds.min_gq {
            own.push(LOW_GQ.to_string());
        }
        if genotype.dp.unwrap_or(-1) < thresholds.min_dp {
            own.push(LOW_DP.to_string());
        }
        let mut copy = genotype.clone();
        if own.is_empty() {
            all_filtered = false;
            // The iterator sets `FT=PASS` and htsjdk's writer drops a passing FT, so the field
            // never reaches the file. Carrying `None` here is the same bytes and says why.
            copy.filters = None;
        } else {
            own.sort();
            copy.filters = Some(own.join(";"));
        }
        genotypes.push(copy);
    }
    if all_filtered {
        filters.push(ALL_GTS_FILTERED.to_string());
    }

    let mut out = record.clone();
    out.genotypes = genotypes;
    filters.sort();
    // `passFilters()` for an empty set, which is `Some(empty)` here and prints PASS.
    out.filters = Some(filters);
    out
}

/// `doWork()`: the whole run.
pub fn filter(input: &str, thresholds: &Thresholds) -> Result<String, FilterError> {
    let file = read_vcf(input).map_err(|failure| {
        FilterError::Vcf(failure.error.class().to_string(), failure.error.message())
    })?;
    let has_dictionary = file
        .header
        .lines
        .iter()
        .any(|line| matches!(line, HeaderLine::Contig { .. }));
    if !has_dictionary {
        return Err(FilterError::NoSequenceDictionary);
    }

    // `header.addMetaDataLine(...)` on the reader's own header, which is then written.
    let mut header = VcfHeader {
        lines: file.header.lines.clone(),
        samples: file.header.samples.clone(),
    };
    for line in added_header_lines() {
        if !header.lines.contains(&line) {
            header.lines.push(line);
        }
    }

    let records: Vec<VariantContext> = file
        .records
        .iter()
        .map(|record| filter_record(record, thresholds))
        .collect();

    write_vcf(&header, &records).map_err(|error| match error {
        EncodeError::MissingFromHeader {
            key,
            field,
            contig,
            start,
        } => FilterError::Vcf(
            "java.lang.IllegalStateException".to_string(),
            format!(
                "Key {key} found in VariantContext field {field} at {contig}:{start} but this key \
                 isn't defined in the VCFHeader.  We require all VCFs to have complete VCF headers \
                 by default."
            ),
        ),
        other => FilterError::Vcf(
            "java.lang.IllegalStateException".to_string(),
            format!("{other:?}"),
        ),
    })
}

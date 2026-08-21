//! `CalculateAverageCombinedAnnotations`, ported from
//! `org.broadinstitute.hellbender.tools.CalculateAverageCombinedAnnotations` (GATK 4.6.2.0).
//!
//! Annotations that GenomicsDB summed across samples, divided by the number of samples that were
//! called het or hom-var.
//!
//! # The divisor is two of three counts, read as doubles
//!
//! ```java
//! List<String> genotypeCounts = variant.getAttributeAsStringList(RAW_GENOTYPE_COUNT_KEY, "");
//! double counter = Double.parseDouble(genotypeCounts.get(1)) + Double.parseDouble(genotypeCounts.get(2));
//! ```
//!
//! `RAW_GT_COUNT` is hom-ref, het, hom-var, and the first is ignored. The two that count are
//! parsed as doubles, so a non-integral count divides as written.
//!
//! # A divisor of zero writes the record through untouched
//!
//! Not a zero average, not a missing value: the record is added exactly as it arrived, so a run
//! can produce a file where some records carry `AVERAGE_` fields and others carry none.
//!
//! # The header gains a line per REQUESTED annotation
//!
//! Whether any record carries it or not, and whether the input declared it or not. A run asking
//! for an annotation the file never mentions still declares `AVERAGE_` for it.

use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use htsjdk_vcf::variant::{Value, VariantContext};

/// `GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY`.
pub const RAW_GENOTYPE_COUNT_KEY: &str = "RAW_GT_COUNT";

/// What the run refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AverageError {
    /// A record with no `RAW_GT_COUNT`, which names the site.
    MissingCounts { contig: String, start: i64 },
    /// No annotation given at all. The parser refuses first, so this is unreachable from a
    /// command line and kept because it is the tool's own.
    NoAnnotations,
}

impl AverageError {
    pub fn java_class(&self) -> &str {
        "org.broadinstitute.hellbender.exceptions.UserException"
    }

    pub fn message(&self) -> String {
        match self {
            AverageError::MissingCounts { contig, start } => {
                format!("Need annotation {RAW_GENOTYPE_COUNT_KEY} at site {contig}:{start}")
            }
            AverageError::NoAnnotations => {
                "--summed-annotation-to-divide must be provided.".to_string()
            }
        }
    }
}

/// The INFO line one requested annotation adds, which quotes the source twice.
pub fn average_header_line(annotation: &str) -> HeaderLine {
    HeaderLine::Compound {
        key: "INFO".to_string(),
        id: format!("AVERAGE_{annotation}"),
        number: Cardinality::Fixed(1),
        line_type: LineType::Float,
        description: format!(
            "Average of {annotation} annotation across samples. See {annotation} header line for \
             more information."
        ),
        extra: Vec::new(),
    }
}

/// `onTraversalStart`: the header the writer gets.
pub fn header_with_averages(
    header: &VcfHeader,
    annotations: &[String],
) -> Result<VcfHeader, AverageError> {
    if annotations.is_empty() {
        return Err(AverageError::NoAnnotations);
    }
    let mut out = header.clone();
    for annotation in annotations {
        out.lines.push(average_header_line(annotation));
    }
    Ok(out)
}

/// `getAttributeAsStringList(key, "")`, for the one key this tool reads.
fn counts(record: &VariantContext) -> Option<Vec<String>> {
    record
        .attributes
        .iter()
        .find(|(key, _)| key == RAW_GENOTYPE_COUNT_KEY)
        .map(|(_, value)| match value {
            Value::List(values) => values.iter().map(render).collect(),
            other => vec![render(other)],
        })
}

fn render(value: &Value) -> String {
    match value {
        Value::Int(number) => number.to_string(),
        Value::Double(number) => number.to_string(),
        Value::Str(text) => text.clone(),
        Value::Bool(flag) => flag.to_string(),
        Value::Missing => ".".to_string(),
        Value::List(values) => values.iter().map(render).collect::<Vec<String>>().join(","),
    }
}

/// `getAttributeAsDouble(annot, 0)`.
fn attribute_as_double(record: &VariantContext, key: &str) -> Option<f64> {
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
}

/// `apply` for one record: the record the writer is handed.
pub fn apply(
    record: &VariantContext,
    annotations: &[String],
) -> Result<VariantContext, AverageError> {
    let counts = counts(record).ok_or_else(|| AverageError::MissingCounts {
        contig: record.contig.clone(),
        start: record.start,
    })?;
    let divisor: f64 = counts
        .get(1)
        .and_then(|value| value.parse::<f64>().ok())
        .unwrap_or(0.0)
        + counts
            .get(2)
            .and_then(|value| value.parse::<f64>().ok())
            .unwrap_or(0.0);
    if divisor <= 0.0 {
        // `if (counter > 0)`, and the else branch adds the variant as it stands.
        return Ok(record.clone());
    }
    let mut out = record.clone();
    for annotation in annotations {
        if let Some(value) = attribute_as_double(record, annotation) {
            out.attributes.push((
                format!("AVERAGE_{annotation}"),
                Value::Double(value / divisor),
            ));
        }
    }
    Ok(out)
}

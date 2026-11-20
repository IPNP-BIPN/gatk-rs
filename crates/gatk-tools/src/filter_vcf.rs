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
//! A record whose every called, non-hom-ref genotype was filtered is itself filtered, and that
//! replaces the `PASS` the three site filters would otherwise have left. Hom-ref and no-call
//! genotypes do not count either way, and a record with no variant genotype is never filtered.
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

use std::sync::Arc;

use htsjdk_vcf::encoder::{EncodeError, VcfEncoder};
use htsjdk_vcf::genotype_parse::{parse_genotypes, GenotypeContext};
use htsjdk_vcf::genotype_type::compare_sample_names;
use htsjdk_vcf::genotypes_context::GenotypesContext;
use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use htsjdk_vcf::reader::{read_vcf, ReadError, VcfFile};
use htsjdk_vcf::record_parse::decode_line;
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
    if !record.genotypes.iter().any(is_het) {
        return None;
    }
    // Upstream walks `getGenotypesOrderedByName()` into a map keyed by the allele list. The tallies
    // are integer sums and the verdict is whether any is below the threshold, so neither order can
    // reach the answer, and the genotypes are taken as they come.
    let mut tallies: Vec<(Vec<String>, i64, i64)> = Vec::new();
    for genotype in record.genotypes.iter() {
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

/// `Genotype.isHet()`: called, and two different alleles, compared by `base_string`.
fn is_het(genotype: &Genotype) -> bool {
    let [first, second] = genotype.alleles.as_slice() else {
        return false;
    };
    if first.is_no_call() || second.is_no_call() {
        return false;
    }
    // `base_string` without building it: a symbolic allele's is empty, any other's is its bases.
    let same = match (first.is_symbolic(), second.is_symbolic()) {
        (true, true) => true,
        (true, false) => second.is_empty(),
        (false, true) => first.is_empty(),
        (false, false) => first.equals_ignore_ref(second),
    };
    !same
}

fn allele_index(record: &VariantContext, bases: &str) -> Option<usize> {
    record
        .alleles
        .iter()
        .position(|allele| allele.base_string() == bases)
}

/// `FilterApplyingVariantIterator.next`, for one record.
pub fn filter_record(record: &VariantContext, thresholds: &Thresholds) -> VariantContext {
    filter_owned(record.clone(), thresholds)
}

/// [`filter_record`] on a record the caller hands over, which is changed rather than copied.
pub fn filter_owned(mut record: VariantContext, thresholds: &Thresholds) -> VariantContext {
    let mut filters: Vec<String> = Vec::new();
    if let Some(name) = allele_balance(&record, thresholds.min_ab) {
        filters.push(name);
    }
    // `getAttributeAsDouble("FS", 0)`, so an absent FS is zero.
    if attribute_as_double(&record, "FS", 0.0) > thresholds.max_fs {
        filters.push(STRAND_BIAS.to_string());
    }
    // `getAttributeAsDouble("QD", -1)`, and a negative value means absent.
    let qd = attribute_as_double(&record, "QD", -1.0);
    if qd >= 0.0 && qd < thresholds.min_qd {
        filters.push(LOW_QD.to_string());
    }

    // Taken out of the context, so the record is written from these rather than from its text.
    let mut genotypes = std::mem::take(&mut record.genotypes).into_vec();
    // `variantSamples`: only a called, non-hom-ref genotype counts, and a record with none is
    // never `AllGtsFiltered` however many hom-ref genotypes were filtered.
    let mut any_variant = false;
    let mut all_variants_filtered = true;
    for genotype in &mut genotypes {
        let variant = genotype.is_called() && !genotype.is_hom_ref();
        any_variant |= variant;
        // `getGQ()` and `getDP()` are -1 when absent. `LowDP` sorts before `LowGQ`.
        let low_gq = genotype.gq.unwrap_or(-1) < thresholds.min_gq;
        let low_dp = genotype.dp.unwrap_or(-1) < thresholds.min_dp;
        genotype.filters = match (low_dp, low_gq) {
            (false, false) => {
                all_variants_filtered &= !variant;
                // The iterator sets `FT=PASS` and htsjdk's writer drops a passing FT, so the field
                // never reaches the file. Carrying `None` here is the same bytes and says why.
                None
            }
            (true, true) => Some(format!("{LOW_DP};{LOW_GQ}")),
            (true, false) => Some(LOW_DP.to_string()),
            (false, true) => Some(LOW_GQ.to_string()),
        };
    }
    if any_variant && all_variants_filtered {
        filters.push(ALL_GTS_FILTERED.to_string());
    }

    record.genotypes = genotypes.into();
    filters.sort();
    // `passFilters()` for an empty set, which is `Some(empty)` here and prints PASS.
    record.filters = Some(filters);
    record
}

/// `doWork()`: the whole run, to a `.vcf` or `.bcf` output.
pub fn filter(input: &str, thresholds: &Thresholds) -> Result<String, FilterError> {
    filter_to(input, thresholds, true)
}

/// `doWork()`: the whole run. Only an output named `.vcf` or `.bcf` (`isVcfOrBcf`) needs the
/// input's sequence dictionary; a `.vcf.gz` one does not.
pub fn filter_to(
    input: &str,
    thresholds: &Thresholds,
    vcf_or_bcf_output: bool,
) -> Result<String, FilterError> {
    let mut out = String::new();
    filter_stream(
        input.lines().map(str::to_string),
        thresholds,
        vcf_or_bcf_output,
        |text| {
            out.push_str(text);
            Ok(())
        },
    )?;
    Ok(out)
}

/// [`filter_to`] one record at a time: `lines` is the input without its line ends, and `emit` gets
/// the output in order. What was emitted is the output only if this returns `Ok`.
///
/// The refusals come in the order the whole-file read gives them: a record the reader refuses,
/// anywhere in the file, before a missing dictionary or a record the encoder refuses. So after
/// either of those the rest is still decoded, and only a read refusal can replace it.
pub fn filter_stream<E>(
    lines: impl Iterator<Item = String>,
    thresholds: &Thresholds,
    vcf_or_bcf_output: bool,
    mut emit: E,
) -> Result<(), FilterError>
where
    E: FnMut(&str) -> Result<(), FilterError>,
{
    let mut lines = lines;
    // The `##` lines and the column line, which is all `read_vcf` reads before the first record.
    let mut header_text = String::new();
    let mut complete = false;
    for line in lines.by_ref() {
        header_text.push_str(&line);
        header_text.push('\n');
        if !line.starts_with("##") {
            complete = line.starts_with('#');
            break;
        }
    }
    if !complete {
        // No column line: the header is refused, or the text is all header. Either way the
        // whole-file path answers it.
        for line in lines {
            header_text.push_str(&line);
            header_text.push('\n');
        }
        return emit(&filter_whole(&header_text, thresholds, vcf_or_bcf_output)?);
    }
    let file = read_vcf(&header_text).map_err(|failure| read_error(&failure.error))?;
    let VcfFile {
        header: read_header,
        codec_version,
        ..
    } = file;
    let mut refusal = (vcf_or_bcf_output && !has_dictionary(&read_header))
        .then_some(FilterError::NoSequenceDictionary);
    let header = output_header(&read_header);
    let encoder = VcfEncoder::new(&header);
    if refusal.is_none() {
        emit(&header.write())?;
    }

    // `read_vcf_with`'s eager loop, one line at a time.
    let samples: Arc<[String]> = Arc::from(read_header.samples.clone());
    let samples_were_already_sorted = read_header
        .samples
        .windows(2)
        .all(|pair| compare_sample_names(&pair[0], &pair[1]) != std::cmp::Ordering::Greater);
    let mut line_number = header_text.lines().count();
    let mut out = String::new();
    for line in lines {
        let decoded = match decode_line(&line, &read_header, line_number, codec_version) {
            Ok(Some(decoded)) => decoded,
            Ok(None) => continue,
            Err(error) => return Err(read_error(&ReadError::Record(error))),
        };
        line_number += 1;
        let mut variant = decoded.variant;
        if let Some(block) = decoded.genotype_block {
            let context = GenotypeContext {
                line: &line,
                header: &read_header,
                version: codec_version,
                contig: &variant.contig,
                pos: variant.start,
                line_number,
            };
            let genotypes = parse_genotypes(&block, &variant.alleles, &context)
                .map_err(|error| read_error(&ReadError::Record(error)))?;
            variant.genotypes = GenotypesContext::lazy(genotypes, block, samples.clone());
            if !samples_were_already_sorted {
                variant.genotypes.decode();
            }
        }
        if refusal.is_some() {
            continue;
        }
        out.clear();
        match encoder.encode_into(&filter_owned(variant, thresholds), &mut out) {
            Ok(()) => {
                out.push('\n');
                emit(&out)?;
            }
            Err(error) => refusal = Some(encode_error(error)),
        }
    }
    refusal.map_or(Ok(()), Err)
}

fn read_error(error: &ReadError) -> FilterError {
    FilterError::Vcf(error.class().to_string(), error.message())
}

fn has_dictionary(header: &VcfHeader) -> bool {
    header
        .lines
        .iter()
        .any(|line| matches!(line, HeaderLine::Contig { .. }))
}

/// `header.addMetaDataLine(...)` on the reader's own header, which is then written.
fn output_header(read: &VcfHeader) -> VcfHeader {
    let mut header = VcfHeader {
        lines: read.lines.clone(),
        samples: read.samples.clone(),
    };
    // `VCFHeader.addMetaDataLine` keeps the line already there under the same key and ID, so an
    // input's own `FT` or `LowQD` line, whatever its description, is the one written.
    let id_of = |line: &HeaderLine| match line {
        HeaderLine::Filter { id, .. } => Some(("FILTER".to_string(), id.clone())),
        HeaderLine::Compound { key, id, .. } => Some((key.clone(), id.clone())),
        _ => None,
    };
    for line in added_header_lines() {
        if !header
            .lines
            .iter()
            .any(|existing| id_of(existing) == id_of(&line))
        {
            header.lines.push(line);
        }
    }
    header
}

fn encode_error(error: EncodeError) -> FilterError {
    match error {
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
    }
}

/// The whole-file run, for a text with no column line.
fn filter_whole(
    input: &str,
    thresholds: &Thresholds,
    vcf_or_bcf_output: bool,
) -> Result<String, FilterError> {
    let file = read_vcf(input).map_err(|failure| read_error(&failure.error))?;
    if vcf_or_bcf_output && !has_dictionary(&file.header) {
        return Err(FilterError::NoSequenceDictionary);
    }
    let header = output_header(&file.header);
    let records: Vec<VariantContext> = file
        .records
        .into_iter()
        .map(|record| filter_owned(record, thresholds))
        .collect();
    write_vcf(&header, &records).map_err(encode_error)
}

//! `SVCallRecordUtils.create` and the `SVCallRecord` constructor it calls, ported from GATK 4.6.2.0.
//!
//! Every SV tool of the `sv` package reads its VCF through this one conversion, so it is the brick
//! the runners of `SVStratify`, `SVCluster`, `GroupedSVCluster`, `SVConcordance` and `SVAnnotate`
//! share: the type inferred from `SVTYPE` or from the alternate allele, the second breakpoint from
//! `CHR2`/`END2` or from `END`, the strands and the length resolved per type, and the attributes
//! stripped of the ten keys the record owns.
//!
//! # The checks run in the order the fields are read
//!
//! ```java
//! type = inferStructuralVariantType(variant);
//! cpxSubtype = getComplexSubtype(variant);
//! cpxIntervals = parseComplexIntervals(variant, dictionary);
//! algorithms = getAlgorithms(variant);
//! evidence = getEvidence(variant);
//! ```
//!
//! So a record with no `SVTYPE` and a non-symbolic allele is refused for its allele before its
//! missing `ALGORITHMS` is noticed, and a BND without `CHR2` is refused only after its strands and
//! its evidence were accepted. [`create`] keeps that order, because a record with two problems is
//! reported by the first.
//!
//! # An insertion's second breakpoint is its first
//!
//! ```java
//! if (type == INS) { positionB = positionA; } else { positionB = variant.getEnd(); }
//! ```
//!
//! `END` is ignored for an insertion however far it reaches, and for the other intrachromosomal
//! types it IS the second breakpoint, so a DEL's length is `END - POS + 1` and `SVLEN` is never read
//! for it: the length is passed as `null` and recomputed by the constructor.
//!
//! # Only an insertion and a complex event read `SVLEN`
//!
//! And for those two it is required: a missing field is an `IllegalArgumentException`, `-1` is the
//! undefined length, and any other negative value an `IllegalStateException`, which `Utils.validate`
//! throws where `Utils.validateArg` throws the other.
//!
//! # The strands are the type's own for a DEL and a DUP
//!
//! `STRANDS` is not even read for DEL, DUP and CNV; the constructor then sets `+-` for a DEL, `-+`
//! for a DUP and nothing for a CNV. Every other type keeps what the field says, or nothing.
//!
//! # One refusal cannot be reproduced
//!
//! An unknown `CPX_TYPE` lists the valid subtypes in the iteration order of a `HashBiMap` built from
//! `Map.ofEntries`, whose order is salted per JVM: two runs of the reference print two different
//! messages. [`SvRecordError::InvalidComplexSubtype`] lists them in declaration order, and no
//! fixture should reach it.

use htsjdk_vcf::variant::{Value, VariantContext};

use crate::sv_stratify::SvType;

/// The `INFO` keys the record takes over, which `sanitizeAttributes` removes from the rest.
pub const INVALID_ATTRIBUTES: &[&str] = &[
    "END",
    ALGORITHMS_ATTRIBUTE,
    SVLEN,
    EVIDENCE,
    CONTIG2_ATTRIBUTE,
    END2_ATTRIBUTE,
    STRANDS_ATTRIBUTE,
    SVTYPE,
    CPX_TYPE,
    CPX_INTERVALS,
];

pub const SVTYPE: &str = "SVTYPE";
pub const SVLEN: &str = "SVLEN";
pub const EVIDENCE: &str = "EVIDENCE";
pub const ALGORITHMS_ATTRIBUTE: &str = "ALGORITHMS";
pub const STRANDS_ATTRIBUTE: &str = "STRANDS";
pub const CONTIG2_ATTRIBUTE: &str = "CHR2";
pub const END2_ATTRIBUTE: &str = "END2";
pub const CPX_TYPE: &str = "CPX_TYPE";
pub const CPX_INTERVALS: &str = "CPX_INTERVALS";

/// `SVCallRecord.UNDEFINED_LENGTH`.
pub const UNDEFINED_LENGTH: i32 = -1;

/// The keys of `COMPLEX_VARIANT_SUBTYPE_MAP`, in the enum's declaration order.
pub const COMPLEX_SUBTYPES: &[&str] = &[
    "delINV",
    "INVdel",
    "dupINV",
    "INVdup",
    "delINVdel",
    "dupINVdup",
    "delINVdup",
    "dupINVdel",
    "piDUP_FR",
    "piDUP_RF",
    "dDUP",
    "dDUP_iDEL",
    "INS_iDEL",
    "CTX_PP/QQ",
    "CTX_PQ/QP",
    "CTX_INV",
];

/// `GATKSVVCFConstants.EvidenceTypes`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Evidence {
    Baf,
    Pe,
    Rd,
    Sr,
}

impl Evidence {
    fn value_of(text: &str) -> Option<Evidence> {
        Some(match text {
            "BAF" => Evidence::Baf,
            "PE" => Evidence::Pe,
            "RD" => Evidence::Rd,
            "SR" => Evidence::Sr,
            _ => return None,
        })
    }
}

/// `SVCallRecord.ComplexEventInterval`: a type and a closed interval, `SVTYPE_chr:pos-end`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ComplexEventInterval {
    pub sv_type: SvType,
    pub contig: String,
    pub start: i32,
    pub end: i32,
}

/// What [`create`] refuses, each with the exception class the reference throws.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SvRecordError {
    /// `Utils.validate(!alleles.isEmpty(), ...)`.
    MissingAltAllele {
        id: String,
    },
    /// More than one alternate that is not the `<DEL>`/`<DUP>` pair of a CNV.
    MultiallelicNotSupported {
        id: String,
    },
    /// A single alternate that is not symbolic, which names no record.
    ExpectedSymbolicAllele,
    /// A symbolic alternate whose name is no type.
    NoValidType {
        id: String,
    },
    /// A `CPX_TYPE` outside the map. See the module documentation for why this message cannot match.
    InvalidComplexSubtype {
        subtype: String,
    },
    /// A `CPX_INTERVALS` entry with no `_`.
    ComplexIntervalFormat {
        entry: String,
    },
    /// A `CPX_INTERVALS` entry whose type is no constant.
    ComplexIntervalType {
        name: String,
    },
    /// A `CPX_INTERVALS` entry off the dictionary, printed as `SimpleInterval.toString()`.
    ComplexIntervalOffDictionary {
        interval: String,
    },
    /// A `CPX_INTERVALS` entry `new SimpleInterval(String)` does not parse.
    ComplexIntervalUnparseable {
        message: String,
    },
    MissingAlgorithms {
        id: String,
    },
    UnknownEvidence {
        name: String,
    },
    StrandsLength {
        id: String,
    },
    StartStrand {
        id: String,
    },
    EndStrand {
        id: String,
    },
    MissingLength {
        id: String,
    },
    NegativeLength {
        id: String,
    },
    /// `getAttributeAsInt` on a value `Integer.valueOf` refuses.
    NotAnInteger {
        text: String,
    },
    MissingSecondBreakpoint {
        id: String,
    },
    /// `validatePosition`: a breakpoint on a contig the dictionary does not name.
    ContigNotInDictionary {
        contig: String,
    },
    /// `validatePosition`: a breakpoint outside `[1, length]`.
    InvalidPosition {
        contig: String,
        position: i32,
    },
    /// `validateCoordinates`: the second breakpoint sorts before the first.
    EndPrecedesStart {
        id: String,
    },
}

impl SvRecordError {
    pub fn class(&self) -> &'static str {
        match self {
            SvRecordError::MissingAltAllele { .. }
            | SvRecordError::MultiallelicNotSupported { .. }
            | SvRecordError::ExpectedSymbolicAllele
            | SvRecordError::NegativeLength { .. } => "java.lang.IllegalStateException",
            SvRecordError::NotAnInteger { .. } => "java.lang.NumberFormatException",
            SvRecordError::MissingSecondBreakpoint { .. } => {
                "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
            }
            _ => "java.lang.IllegalArgumentException",
        }
    }

    pub fn is_user(&self) -> bool {
        matches!(self, SvRecordError::MissingSecondBreakpoint { .. })
    }

    pub fn message(&self) -> String {
        match self {
            SvRecordError::MissingAltAllele { id } => format!("Missing alt allele for variant {id}"),
            SvRecordError::MultiallelicNotSupported { id } => {
                format!("Non-CNV multiallelic variants not supported (variant {id})")
            }
            SvRecordError::ExpectedSymbolicAllele => "Expected symbolic alt allele".to_string(),
            SvRecordError::NoValidType { id } => {
                format!("Could not find a valid SV type for variant {id}")
            }
            SvRecordError::InvalidComplexSubtype { subtype } => format!(
                "Invalid CPX subtype: {subtype}, valid values are: {}",
                COMPLEX_SUBTYPES.join(", ")
            ),
            SvRecordError::ComplexIntervalFormat { entry } => format!(
                "Expected complex interval with format \"SVTYPE_chr:pos-end\" but found \"{entry}\""
            ),
            SvRecordError::ComplexIntervalType { name } => format!(
                "No enum constant org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.StructuralVariantAnnotationType.{name}"
            ),
            SvRecordError::ComplexIntervalOffDictionary { interval } => {
                format!("Invalid CPX interval: {interval}")
            }
            SvRecordError::ComplexIntervalUnparseable { message } => message.clone(),
            SvRecordError::MissingAlgorithms { id } => {
                format!("Expected {ALGORITHMS_ATTRIBUTE} field for variant {id}")
            }
            SvRecordError::UnknownEvidence { name } => format!(
                "No enum constant org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.EvidenceTypes.{name}"
            ),
            SvRecordError::StrandsLength { id } => {
                format!("Strands field is not 2 characters long for variant {id}")
            }
            SvRecordError::StartStrand { id } => {
                format!("Valid start strand not found for variant {id}")
            }
            SvRecordError::EndStrand { id } => format!("Valid end strand not found for variant {id}"),
            SvRecordError::MissingLength { id } => {
                format!("Expected {SVLEN} field for variant {id}")
            }
            SvRecordError::NegativeLength { id } => format!(
                "Length must be non-negative or {UNDEFINED_LENGTH} for variant {id}"
            ),
            SvRecordError::NotAnInteger { text } => format!("For input string: \"{text}\""),
            SvRecordError::MissingSecondBreakpoint { id } => format!(
                "Bad input: Attributes {END2_ATTRIBUTE} and {CONTIG2_ATTRIBUTE} are required for BND and CTX records (variant {id})."
            ),
            SvRecordError::ContigNotInDictionary { contig } => {
                format!("Contig {contig} not found in dictionary")
            }
            SvRecordError::InvalidPosition { contig, position } => {
                format!("Invalid position {contig}:{position}")
            }
            SvRecordError::EndPrecedesStart { id } => {
                format!("End precedes start in variant {id}")
            }
        }
    }
}

/// `SVCallRecord`, less the alleles and genotypes, which a caller reads off the variant itself.
#[derive(Debug, Clone, PartialEq)]
pub struct SvCallRecord {
    pub id: String,
    pub contig_a: String,
    pub position_a: i32,
    pub strand_a: Option<bool>,
    pub contig_b: String,
    pub position_b: i32,
    pub strand_b: Option<bool>,
    pub sv_type: SvType,
    pub cpx_subtype: Option<String>,
    /// Sorted by their encoding, as `canonicalizeComplexEventList` sorts them.
    pub cpx_intervals: Vec<ComplexEventInterval>,
    pub length: Option<i32>,
    pub evidence: Vec<Evidence>,
    pub algorithms: Vec<String>,
    /// The variant's attributes without [`INVALID_ATTRIBUTES`].
    pub attributes: Vec<(String, Value)>,
    pub log10_p_error: Option<f64>,
}

impl SvCallRecord {
    /// What [`crate::sv_stratify`] reads off a record.
    pub fn stratify_record(&self) -> crate::sv_stratify::CallRecord {
        crate::sv_stratify::CallRecord {
            id: self.id.clone(),
            sv_type: self.sv_type,
            contig_a: self.contig_a.clone(),
            position_a: self.position_a,
            contig_b: self.contig_b.clone(),
            position_b: self.position_b,
            length: self.length,
        }
    }
}

/// `validateCoordinates`, which the constructor taking a dictionary runs and `create` does not.
///
/// Both breakpoints must lie on a contig of the dictionary and inside it, the second may not sort
/// before the first unless the record is complex, and every complex interval is checked the same
/// way at both of its ends.
pub fn validate_coordinates(
    record: &SvCallRecord,
    dictionary: &[(String, i32)],
) -> Result<(), SvRecordError> {
    let position = |contig: &str, position: i32| -> Result<usize, SvRecordError> {
        let Some(index) = dictionary.iter().position(|(name, _)| name == contig) else {
            return Err(SvRecordError::ContigNotInDictionary {
                contig: contig.to_string(),
            });
        };
        if position <= 0 || position > dictionary[index].1 {
            return Err(SvRecordError::InvalidPosition {
                contig: contig.to_string(),
                position,
            });
        }
        Ok(index)
    };
    let a = position(&record.contig_a, record.position_a)?;
    let b = position(&record.contig_b, record.position_b)?;
    if record.sv_type != SvType::Cpx && (a, record.position_a) > (b, record.position_b) {
        return Err(SvRecordError::EndPrecedesStart {
            id: record.id.clone(),
        });
    }
    for interval in &record.cpx_intervals {
        position(&interval.contig, interval.start)?;
        position(&interval.contig, interval.end)?;
    }
    Ok(())
}

/// `SVCallRecordUtils.compareCalls`: both breakpoints in dictionary order, then the type's
/// ordinal, then the strands and the length, an absent value sorting first.
pub fn compare_calls(
    first: &SvCallRecord,
    second: &SvCallRecord,
    dictionary: &[(String, i32)],
) -> std::cmp::Ordering {
    let index = |contig: &str| {
        dictionary
            .iter()
            .position(|(name, _)| name == contig)
            .unwrap_or(usize::MAX)
    };
    let ordinal = |sv_type: SvType| match sv_type {
        SvType::Del => 0,
        SvType::Dup => 1,
        SvType::Ins => 2,
        SvType::Inv => 3,
        SvType::Cpx => 4,
        SvType::Bnd => 5,
        SvType::Ctx => 6,
        SvType::Cnv => 7,
    };
    (index(&first.contig_a), first.position_a)
        .cmp(&(index(&second.contig_a), second.position_a))
        .then(
            (index(&first.contig_b), first.position_b)
                .cmp(&(index(&second.contig_b), second.position_b)),
        )
        .then(ordinal(first.sv_type).cmp(&ordinal(second.sv_type)))
        .then(first.strand_a.cmp(&second.strand_a))
        .then(first.strand_b.cmp(&second.strand_b))
        .then(first.length.cmp(&second.length))
}

fn evidence_name(evidence: Evidence) -> &'static str {
    match evidence {
        Evidence::Baf => "BAF",
        Evidence::Pe => "PE",
        Evidence::Rd => "RD",
        Evidence::Sr => "SR",
    }
}

/// `SVCallRecordUtils.getVariantBuilder`: the record back into a variant.
///
/// The fields the record owns are written again from the record rather than copied from the input,
/// so what comes out is not what went in: a deletion's `SVLEN` is the length `END - POS + 1`,
/// positive, whatever sign the input used; `END` of a breakend is its own position; `STRANDS` is
/// written only for the types that keep them. The filters are set only when there are some, and
/// `getFilters()` answers the empty set for `PASS` as for `.`, so a passing input comes out
/// UNFILTERED. `alleles` is the variant's own list, reference first; `filters` the input's.
pub fn to_variant(
    record: &SvCallRecord,
    alleles: Vec<htsjdk_vcf::allele::Allele>,
    genotypes: Vec<htsjdk_vcf::variant::Genotype>,
    filters: &[String],
) -> VariantContext {
    let (end, second) = match record.sv_type {
        SvType::Bnd | SvType::Ctx => (
            record.position_a,
            Some((record.position_b, record.contig_b.clone())),
        ),
        _ => (record.position_b, None),
    };
    let mut variant = VariantContext::new(&record.contig_a, record.position_a as i64, alleles);
    variant.id = record.id.clone();
    variant.stop = end as i64;
    let mut attributes = record.attributes.clone();
    let mut put =
        |key: &str, value: Value| match attributes.iter_mut().find(|(name, _)| name == key) {
            Some(slot) => slot.1 = value,
            None => attributes.push((key.to_string(), value)),
        };
    put("END", Value::Int(end as i64));
    put(SVTYPE, Value::Str(type_name(record.sv_type).to_string()));
    put(
        ALGORITHMS_ATTRIBUTE,
        Value::List(
            record
                .algorithms
                .iter()
                .map(|a| Value::Str(a.clone()))
                .collect(),
        ),
    );
    if let Some((position, contig)) = second {
        put(END2_ATTRIBUTE, Value::Int(position as i64));
        put(CONTIG2_ATTRIBUTE, Value::Str(contig));
    }
    if let Some(subtype) = &record.cpx_subtype {
        put(CPX_TYPE, Value::Str(subtype.clone()));
    }
    if !record.cpx_intervals.is_empty() {
        put(
            CPX_INTERVALS,
            Value::List(
                record
                    .cpx_intervals
                    .iter()
                    .map(|interval| {
                        Value::Str(format!(
                            "{}_{}:{}-{}",
                            type_name(interval.sv_type),
                            interval.contig,
                            interval.start,
                            interval.end
                        ))
                    })
                    .collect(),
            ),
        );
    }
    if let Some(length) = record.length {
        put(SVLEN, Value::Int(length as i64));
    }
    if matches!(
        record.sv_type,
        SvType::Bnd | SvType::Inv | SvType::Ins | SvType::Cpx | SvType::Ctx
    ) {
        if let (Some(a), Some(b)) = (record.strand_a, record.strand_b) {
            let strand = |forward: bool| if forward { "+" } else { "-" };
            put(
                STRANDS_ATTRIBUTE,
                Value::Str(format!("{}{}", strand(a), strand(b))),
            );
        }
    }
    if !record.evidence.is_empty() {
        put(
            EVIDENCE,
            Value::List(
                record
                    .evidence
                    .iter()
                    .map(|evidence| Value::Str(evidence_name(*evidence).to_string()))
                    .collect(),
            ),
        );
    }
    variant.attributes = attributes;
    if !filters.is_empty() {
        variant.filters = Some(filters.to_vec());
    }
    if let Some(error) = record.log10_p_error {
        variant.log10_p_error = error;
    }
    // "htsjdk vcf encoder does not allow genotypes to have empty alleles".
    variant.genotypes = genotypes
        .into_iter()
        .map(|mut genotype| {
            if genotype.alleles.is_empty() {
                genotype.alleles = vec![htsjdk_vcf::allele::Allele::no_call()];
            }
            genotype
        })
        .collect::<Vec<_>>()
        .into();
    variant
}

fn attribute<'a>(variant: &'a VariantContext, key: &str) -> Option<&'a Value> {
    variant
        .attributes
        .iter()
        .find(|(name, _)| name == key)
        .map(|(_, value)| value)
}

/// `getAttributeAsString(key, default)`: a list joined at `,`.
fn as_string(value: &Value) -> String {
    match value {
        Value::Str(text) => text.clone(),
        Value::List(items) => items.iter().map(as_string).collect::<Vec<_>>().join(","),
        Value::Bool(flag) => flag.to_string(),
        other => other.format().unwrap_or_default(),
    }
}

/// `getAttributeAsStringList(key, null)`: nothing when absent, one entry per list element.
fn as_string_list(variant: &VariantContext, key: &str) -> Vec<String> {
    match attribute(variant, key) {
        None => Vec::new(),
        Some(Value::List(items)) => items.iter().map(as_string).collect(),
        Some(other) => vec![as_string(other)],
    }
}

/// `getAttributeAsInt(key, default)`.
fn as_int(variant: &VariantContext, key: &str, default: i32) -> Result<i32, SvRecordError> {
    match attribute(variant, key) {
        None | Some(Value::Missing) => Ok(default),
        Some(Value::Int(value)) => Ok(*value as i32),
        Some(Value::Double(value)) => Ok(*value as i32),
        Some(other) => {
            let text = as_string(other);
            text.parse::<i32>()
                .map_err(|_| SvRecordError::NotAnInteger { text })
        }
    }
}

fn has(variant: &VariantContext, key: &str) -> bool {
    attribute(variant, key).is_some()
}

/// `inferStructuralVariantType`: `SVTYPE` when it names a type, else the alternate allele.
pub fn infer_type(variant: &VariantContext) -> Result<SvType, SvRecordError> {
    let named = attribute(variant, SVTYPE)
        .map(as_string)
        .unwrap_or_default();
    if let Some(sv_type) = SvType::parse(&named) {
        return Ok(sv_type);
    }
    let id = variant.id.clone();
    let alternates: Vec<String> = variant
        .alternate_alleles()
        .iter()
        .map(|allele| allele.display_string())
        .collect();
    if alternates.is_empty() {
        return Err(SvRecordError::MissingAltAllele { id });
    }
    if alternates.len() == 2
        && alternates.iter().any(|allele| allele == "<DEL>")
        && alternates.iter().any(|allele| allele == "<DUP>")
    {
        return Ok(SvType::Cnv);
    }
    if alternates.len() != 1 {
        return Err(SvRecordError::MultiallelicNotSupported { id });
    }
    let allele = &variant.alternate_alleles()[0];
    if !allele.is_symbolic() {
        return Err(SvRecordError::ExpectedSymbolicAllele);
    }
    let name = allele.display_string().replace(['<', '>'], "");
    SvType::parse(&name).ok_or(SvRecordError::NoValidType { id })
}

fn type_name(sv_type: SvType) -> &'static str {
    match sv_type {
        SvType::Del => "DEL",
        SvType::Dup => "DUP",
        SvType::Ins => "INS",
        SvType::Inv => "INV",
        SvType::Cnv => "CNV",
        SvType::Cpx => "CPX",
        SvType::Bnd => "BND",
        SvType::Ctx => "CTX",
    }
}

/// `new SimpleInterval(String)`, for the one shape a complex interval takes.
fn simple_interval(text: &str) -> Result<(String, i32, i32), SvRecordError> {
    let unparseable = |message: String| SvRecordError::ComplexIntervalUnparseable { message };
    if text.is_empty() {
        return Err(unparseable("str should not be empty".to_string()));
    }
    let Some(colon) = text.rfind(':') else {
        return Ok((text.to_string(), 1, i32::MAX));
    };
    let contig = text[..colon].to_string();
    let rest = &text[colon + 1..];
    let position = |value: &str| -> Result<i32, SvRecordError> {
        value.replace(',', "").parse::<i32>().map_err(|_| {
            unparseable(format!(
                "Problem parsing start/end value in interval string. Value was: {value}"
            ))
        })
    };
    let (start, end) = match rest.find('-') {
        None if rest.ends_with('+') => (position(&rest[..rest.len() - 1])?, i32::MAX),
        None => {
            let start = position(rest)?;
            (start, start)
        }
        Some(dash) => (position(&rest[..dash])?, position(&rest[dash + 1..])?),
    };
    Ok((contig, start, end))
}

/// `ComplexEventInterval.decode`.
fn decode_complex_interval(
    entry: &str,
    dictionary: &[(String, i32)],
) -> Result<ComplexEventInterval, SvRecordError> {
    let Some((name, rest)) = entry.split_once('_') else {
        return Err(SvRecordError::ComplexIntervalFormat {
            entry: entry.to_string(),
        });
    };
    let (contig, start, end) = simple_interval(rest)?;
    let on_contig = dictionary
        .iter()
        .find(|(sequence, _)| *sequence == contig)
        .is_some_and(|(_, length)| start >= 1 && end <= *length && start <= end);
    if !on_contig {
        return Err(SvRecordError::ComplexIntervalOffDictionary {
            interval: format!("{contig}:{start}-{end}"),
        });
    }
    let sv_type = SvType::parse(name).ok_or_else(|| SvRecordError::ComplexIntervalType {
        name: name.to_string(),
    })?;
    Ok(ComplexEventInterval {
        sv_type,
        contig,
        start,
        end,
    })
}

/// `getStrands`, which only the types without strands of their own reach.
fn strands(variant: &VariantContext) -> Result<Option<(bool, bool)>, SvRecordError> {
    let Some(value) = attribute(variant, STRANDS_ATTRIBUTE) else {
        return Ok(None);
    };
    let text = as_string(value);
    let id = variant.id.clone();
    if text.chars().count() != 2 {
        return Err(SvRecordError::StrandsLength { id });
    }
    let mut chars = text.chars();
    let first = chars.next().unwrap_or(' ');
    let second = chars.next().unwrap_or(' ');
    if first != '+' && first != '-' {
        return Err(SvRecordError::StartStrand { id });
    }
    if second != '+' && second != '-' {
        return Err(SvRecordError::EndStrand { id });
    }
    Ok(Some((first == '+', second == '+')))
}

/// `SVCallRecordUtils.create(variant, true, dictionary)`, the constructor included.
///
/// `dictionary` is the sequence dictionary as name and length, which only `CPX_INTERVALS` reads:
/// the constructor `create` calls is the one that does NOT validate the breakpoints against it.
pub fn create(
    variant: &VariantContext,
    dictionary: &[(String, i32)],
) -> Result<SvCallRecord, SvRecordError> {
    let id = variant.id.clone();
    let sv_type = infer_type(variant)?;

    let cpx_subtype = match attribute(variant, CPX_TYPE) {
        None => None,
        Some(value) => {
            let subtype = as_string(value);
            if !COMPLEX_SUBTYPES.contains(&subtype.as_str()) {
                return Err(SvRecordError::InvalidComplexSubtype { subtype });
            }
            Some(subtype)
        }
    };
    let mut cpx_intervals = as_string_list(variant, CPX_INTERVALS)
        .iter()
        .map(|entry| decode_complex_interval(entry, dictionary))
        .collect::<Result<Vec<_>, _>>()?;
    // `canonicalizeComplexEventList`: sorted by `encode()`, the type's name then the interval.
    cpx_intervals.sort_by_key(|interval| {
        format!(
            "{}_{}:{}-{}",
            type_name(interval.sv_type),
            interval.contig,
            interval.start,
            interval.end
        )
    });

    if !has(variant, ALGORITHMS_ATTRIBUTE) {
        return Err(SvRecordError::MissingAlgorithms { id });
    }
    let algorithms = as_string_list(variant, ALGORITHMS_ATTRIBUTE);

    let evidence = as_string_list(variant, EVIDENCE)
        .into_iter()
        .filter(|value| value != ".")
        .map(|value| {
            Evidence::value_of(&value).ok_or(SvRecordError::UnknownEvidence { name: value })
        })
        .collect::<Result<Vec<_>, _>>()?;

    let input_strands = match sv_type {
        SvType::Del | SvType::Cnv | SvType::Dup => None,
        _ => strands(variant)?,
    };

    let input_length = match sv_type {
        SvType::Bnd | SvType::Del | SvType::Dup | SvType::Cnv | SvType::Inv | SvType::Ctx => None,
        _ => {
            if !has(variant, SVLEN) {
                return Err(SvRecordError::MissingLength { id });
            }
            let length = as_int(variant, SVLEN, UNDEFINED_LENGTH)?;
            if length == UNDEFINED_LENGTH {
                None
            } else if length < 0 {
                return Err(SvRecordError::NegativeLength { id });
            } else {
                Some(length)
            }
        }
    };

    let contig_a = variant.contig.clone();
    let position_a = variant.start as i32;
    let (contig_b, position_b) = match sv_type {
        SvType::Bnd | SvType::Ctx => {
            if !(has(variant, CONTIG2_ATTRIBUTE) && has(variant, END2_ATTRIBUTE)) {
                return Err(SvRecordError::MissingSecondBreakpoint { id });
            }
            let contig = attribute(variant, CONTIG2_ATTRIBUTE)
                .map(as_string)
                .unwrap_or_default();
            (contig, as_int(variant, END2_ATTRIBUTE, 0)?)
        }
        SvType::Ins => (contig_a.clone(), position_a),
        _ => (contig_a.clone(), variant.stop as i32),
    };
    let log10_p_error = variant.has_log10_p_error().then_some(variant.log10_p_error);
    let attributes = variant
        .attributes
        .iter()
        .filter(|(key, _)| !INVALID_ATTRIBUTES.contains(&key.as_str()))
        .cloned()
        .collect();

    // The constructor: `inferLength`, then `inferStrands`. `create` passes a null length for every
    // type the constructor recomputes, so its "does not match" check is never reached from here.
    let length = match sv_type {
        SvType::Cnv | SvType::Del | SvType::Dup | SvType::Inv => Some(position_b - position_a + 1),
        _ => input_length,
    };
    let (strand_a, strand_b) = match sv_type {
        SvType::Cnv => (None, None),
        SvType::Del => (Some(true), Some(false)),
        SvType::Dup => (Some(false), Some(true)),
        _ => match input_strands {
            Some((a, b)) => (Some(a), Some(b)),
            None => (None, None),
        },
    };

    Ok(SvCallRecord {
        id,
        contig_a,
        position_a,
        strand_a,
        contig_b,
        position_b,
        strand_b,
        sv_type,
        cpx_subtype,
        cpx_intervals,
        length,
        evidence,
        algorithms,
        attributes,
        log10_p_error,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_vcf::allele::Allele;

    fn variant(alternate: &str, attributes: &[(&str, &str)]) -> VariantContext {
        let mut record = VariantContext::new(
            "chr1",
            1000,
            vec![
                Allele::create(b"N", true).unwrap(),
                Allele::create(alternate.as_bytes(), false).unwrap(),
            ],
        );
        record.id = "sv1".to_string();
        for (key, value) in attributes {
            let value = if value.contains(',') {
                Value::List(
                    value
                        .split(',')
                        .map(|v| Value::Str(v.to_string()))
                        .collect(),
                )
            } else {
                Value::Str(value.to_string())
            };
            record.attributes.push((key.to_string(), value));
        }
        if let Some((_, Value::Str(end))) = record.attributes.iter().find(|(k, _)| k == "END") {
            record.stop = end.parse().unwrap();
        }
        record
    }

    fn dictionary() -> Vec<(String, i32)> {
        vec![("chr1".to_string(), 100000), ("chr2".to_string(), 50000)]
    }

    #[test]
    fn a_deletion_takes_its_length_from_end_and_its_strands_from_its_type() {
        let record = create(
            &variant(
                "<DEL>",
                &[("END", "1999"), ("ALGORITHMS", "manta"), ("SVLEN", "5")],
            ),
            &dictionary(),
        )
        .unwrap();
        assert_eq!(record.sv_type, SvType::Del);
        assert_eq!(record.position_b, 1999);
        assert_eq!(record.length, Some(1000), "SVLEN is never read for a DEL");
        assert_eq!(
            (record.strand_a, record.strand_b),
            (Some(true), Some(false))
        );
        assert!(
            record.attributes.is_empty(),
            "every key read was a reserved one"
        );
    }

    #[test]
    fn an_insertion_ends_where_it_starts_and_needs_its_length() {
        let record = create(
            &variant(
                "<INS>",
                &[
                    ("END", "1500"),
                    ("ALGORITHMS", "manta"),
                    ("SVLEN", "300"),
                    ("STRANDS", "+-"),
                ],
            ),
            &dictionary(),
        )
        .unwrap();
        assert_eq!(record.position_b, 1000);
        assert_eq!(record.length, Some(300));
        assert_eq!(
            (record.strand_a, record.strand_b),
            (Some(true), Some(false))
        );

        let error =
            create(&variant("<INS>", &[("ALGORITHMS", "manta")]), &dictionary()).unwrap_err();
        assert_eq!(error.message(), "Expected SVLEN field for variant sv1");
        let undefined = create(
            &variant("<INS>", &[("ALGORITHMS", "manta"), ("SVLEN", "-1")]),
            &dictionary(),
        )
        .unwrap();
        assert_eq!(undefined.length, None);
        let negative = create(
            &variant("<INS>", &[("ALGORITHMS", "manta"), ("SVLEN", "-2")]),
            &dictionary(),
        )
        .unwrap_err();
        assert_eq!(negative.class(), "java.lang.IllegalStateException");
    }

    #[test]
    fn a_breakend_needs_both_chr2_and_end2() {
        let error = create(
            &variant("<BND>", &[("ALGORITHMS", "manta"), ("CHR2", "chr2")]),
            &dictionary(),
        )
        .unwrap_err();
        assert!(error.is_user());
        assert_eq!(
            error.message(),
            "Bad input: Attributes END2 and CHR2 are required for BND and CTX records (variant sv1)."
        );
        let record = create(
            &variant(
                "<BND>",
                &[
                    ("ALGORITHMS", "manta"),
                    ("CHR2", "chr2"),
                    ("END2", "700"),
                    ("STRANDS", "-+"),
                ],
            ),
            &dictionary(),
        )
        .unwrap();
        assert_eq!((record.contig_b.as_str(), record.position_b), ("chr2", 700));
        assert_eq!(record.length, None);
    }

    #[test]
    fn svtype_outranks_the_allele_and_the_allele_needs_brackets() {
        let named = create(
            &variant(
                "<DEL>",
                &[("SVTYPE", "INV"), ("END", "1100"), ("ALGORITHMS", "depth")],
            ),
            &dictionary(),
        )
        .unwrap();
        assert_eq!(named.sv_type, SvType::Inv);
        let plain = create(&variant("A", &[("ALGORITHMS", "depth")]), &dictionary()).unwrap_err();
        assert_eq!(plain, SvRecordError::ExpectedSymbolicAllele);
        let unknown = create(&variant("<FOO>", &[]), &dictionary()).unwrap_err();
        assert_eq!(
            unknown.message(),
            "Could not find a valid SV type for variant sv1"
        );
    }

    #[test]
    fn the_type_is_refused_before_the_missing_algorithms() {
        let error = create(&variant("<FOO>", &[]), &dictionary()).unwrap_err();
        assert!(matches!(error, SvRecordError::NoValidType { .. }));
        let error = create(&variant("<DEL>", &[("END", "1100")]), &dictionary()).unwrap_err();
        assert_eq!(error.message(), "Expected ALGORITHMS field for variant sv1");
    }

    #[test]
    fn evidence_skips_the_missing_value_and_refuses_an_unknown_name() {
        let record = create(
            &variant(
                "<DEL>",
                &[
                    ("END", "1100"),
                    ("ALGORITHMS", "depth"),
                    ("EVIDENCE", "RD,.,PE"),
                ],
            ),
            &dictionary(),
        )
        .unwrap();
        assert_eq!(record.evidence, vec![Evidence::Rd, Evidence::Pe]);
        let error = create(
            &variant(
                "<DEL>",
                &[("END", "1100"), ("ALGORITHMS", "depth"), ("EVIDENCE", "XX")],
            ),
            &dictionary(),
        )
        .unwrap_err();
        assert!(error.message().ends_with("EvidenceTypes.XX"));
    }

    #[test]
    fn complex_intervals_are_sorted_by_their_encoding_and_checked_against_the_dictionary() {
        let record = create(
            &variant(
                "<CPX>",
                &[
                    ("ALGORITHMS", "manta"),
                    ("SVLEN", "500"),
                    ("CPX_TYPE", "dupINV"),
                    ("CPX_INTERVALS", "INV_chr1:1000-1500,DUP_chr1:900-1000"),
                ],
            ),
            &dictionary(),
        )
        .unwrap();
        assert_eq!(record.cpx_intervals[0].sv_type, SvType::Dup);
        let error = create(
            &variant(
                "<CPX>",
                &[
                    ("ALGORITHMS", "manta"),
                    ("SVLEN", "500"),
                    ("CPX_INTERVALS", "DUP_chr3:1-10"),
                ],
            ),
            &dictionary(),
        )
        .unwrap_err();
        assert_eq!(error.message(), "Invalid CPX interval: chr3:1-10");
    }

    #[test]
    fn a_strands_field_is_checked_character_by_character() {
        let long = create(
            &variant(
                "<INV>",
                &[("END", "1100"), ("ALGORITHMS", "manta"), ("STRANDS", "+++")],
            ),
            &dictionary(),
        )
        .unwrap_err();
        assert_eq!(
            long.message(),
            "Strands field is not 2 characters long for variant sv1"
        );
        let second = create(
            &variant(
                "<INV>",
                &[("END", "1100"), ("ALGORITHMS", "manta"), ("STRANDS", "+x")],
            ),
            &dictionary(),
        )
        .unwrap_err();
        assert_eq!(
            second.message(),
            "Valid end strand not found for variant sv1"
        );
    }
}

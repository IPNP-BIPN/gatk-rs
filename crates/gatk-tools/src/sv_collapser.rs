//! `CanonicalSVCollapser`, ported from GATK 4.6.2.0: a cluster of SV records flattened into the one
//! record that represents it.
//!
//! [`crate::sv_cluster`] decides which records belong together; this decides what the group is
//! written as. Every field is summarised by a rule of its own, and the rules do not read the same
//! members:
//!
//! * **the breakpoints come from the most precise members only**: when any member has evidence
//!   other than read depth, the depth-only ones are dropped before the interval is summarised;
//! * **the type, the algorithms, the alternate alleles, the genotypes and the filters come from
//!   every member**;
//! * **the ID, the evidence, the complex subtype and the attributes come from ONE representative**,
//!   the precise member closest to the summarised interval, ties to the smallest ID.
//!
//! # The interval is summarised before the representative is chosen
//!
//! So under `MEDIAN_START_MEDIAN_END` the representative is whichever member lies nearest the
//! medians, and its own coordinates are not the output's: the record keeps its ID and loses its
//! positions. Only `REPRESENTATIVE` picks a member first, by a chain of its own (quality, then
//! split-read before paired-end evidence, then carrier count, then total distance to every member).
//!
//! # The reference allele is read again
//!
//! `collapseRefAlleles` asks the reference for the base at the new start, whatever the members'
//! REF said, and every genotype's reference allele is replaced by it.
//!
//! # A genotype is the best of the sample's, not a merge
//!
//! Per sample, one member's genotype wins by a comparator chain (non-reference at all, number
//! called, GQ, number non-reference, CNQ, distance of CN from ECN, deletion over duplication), and
//! it is copied whole. `Stream.max` keeps the FIRST of equal genotypes, which is member order.

use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Genotype, Value};

use crate::sv_call_record::{Evidence, SvCallRecord};
use crate::sv_stratify::SvType;

/// `BreakpointSummaryStrategy`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BreakpointSummary {
    MedianStartMedianEnd,
    MinStartMaxEnd,
    MaxStartMinEnd,
    MeanStartMeanEnd,
    Representative,
}

impl BreakpointSummary {
    pub fn value_of(name: &str) -> Option<BreakpointSummary> {
        Some(match name {
            "MEDIAN_START_MEDIAN_END" => BreakpointSummary::MedianStartMedianEnd,
            "MIN_START_MAX_END" => BreakpointSummary::MinStartMaxEnd,
            "MAX_START_MIN_END" => BreakpointSummary::MaxStartMinEnd,
            "MEAN_START_MEAN_END" => BreakpointSummary::MeanStartMeanEnd,
            "REPRESENTATIVE" => BreakpointSummary::Representative,
            _ => return None,
        })
    }
}

/// `AltAlleleSummaryStrategy`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AltAlleleSummary {
    MostSpecificSubtype,
    CommonSubtype,
}

/// `FlagFieldLogic`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FlagFieldLogic {
    And,
    Or,
    AlwaysFalse,
}

/// The two flag fields collapsed by [`FlagFieldLogic`] rather than copied from the representative.
pub const FLAG_TYPE_INFO_FIELDS: &[&str] = &["BOTHSIDES_SUPPORT", "HIGH_SR_BACKGROUND"];

/// `CLUSTER_MEMBER_IDS_KEY`.
pub const CLUSTER_MEMBER_IDS_KEY: &str = "MEMBERS";

/// What collapsing refuses, with the reference's class.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CollapseError {
    MultipleContigsA,
    MultipleContigsB,
    IncompatibleTypes { types: Vec<String> },
    NonSymbolicAllele { allele: String },
    InvalidReferenceLocus { contig: String, position: i32 },
    NullBreakendStrand { first: bool },
    SubtypedMultiallelic,
    MissingExpectedCopyNumber,
}

impl CollapseError {
    pub fn class(&self) -> &'static str {
        match self {
            CollapseError::MultipleContigsA
            | CollapseError::MultipleContigsB
            | CollapseError::InvalidReferenceLocus { .. }
            | CollapseError::SubtypedMultiallelic => "java.lang.IllegalStateException",
            _ => "java.lang.IllegalArgumentException",
        }
    }

    pub fn message(&self) -> String {
        match self {
            CollapseError::MultipleContigsA => {
                "Cannot collapse intervals with multiple position A contigs".to_string()
            }
            CollapseError::MultipleContigsB => {
                "Cannot collapse intervals with multiple position B contigs".to_string()
            }
            CollapseError::IncompatibleTypes { types } => format!(
                "Incompatible SV types found in cluster: {}",
                types.join(", ")
            ),
            CollapseError::NonSymbolicAllele { allele } => {
                format!("Cannot collapse non-symbolic allele: {allele}")
            }
            CollapseError::InvalidReferenceLocus { contig, position } => {
                format!("Invalid reference locus {contig}:{position}")
            }
            CollapseError::NullBreakendStrand { first } => {
                if *first {
                    "First breakend strand cannot be null".to_string()
                } else {
                    "Second breakend strand cannot be null".to_string()
                }
            }
            CollapseError::SubtypedMultiallelic => {
                "Multi-allelic variants with subtyped alleles are not supported.".to_string()
            }
            CollapseError::MissingExpectedCopyNumber => {
                "Genotype missing required field ECN".to_string()
            }
        }
    }
}

/// One member of a cluster: the converted call and what the collapser reads off the variant.
#[derive(Debug, Clone)]
pub struct Member {
    pub call: SvCallRecord,
    /// The variant's alleles, reference first.
    pub alleles: Vec<Allele>,
    pub genotypes: Vec<Genotype>,
    /// `getFilters()`, empty for `PASS` as for `.`.
    pub filters: Vec<String>,
}

impl Member {
    fn alt_alleles(&self) -> Vec<&Allele> {
        self.alleles
            .iter()
            .filter(|allele| !allele.is_no_call() && !allele.is_reference())
            .collect()
    }

    fn is_depth_only(&self) -> bool {
        self.call.algorithms.len() == 1 && self.call.algorithms[0] == "depth"
    }
}

/// The collapsed record, before the walker fills in missing samples and renames it.
#[derive(Debug, Clone)]
pub struct Collapsed {
    pub call: SvCallRecord,
    pub alleles: Vec<Allele>,
    pub genotypes: Vec<Genotype>,
    pub filters: Vec<String>,
}

/// An integer FORMAT field, `getAttributeAsInt(genotype, key, default)`.
pub fn genotype_int(genotype: &Genotype, key: &str, default: i32) -> i32 {
    if key == "GQ" {
        return genotype.gq.unwrap_or(default);
    }
    genotype
        .extended
        .iter()
        .find(|(name, _)| name == key)
        .and_then(|(_, value)| match value {
            Value::Int(value) => Some(*value as i32),
            Value::Double(value) => Some(*value as i32),
            Value::Missing => None,
            other => other.format().and_then(|text| text.parse().ok()),
        })
        .unwrap_or(default)
}

fn has_extended(genotype: &Genotype, key: &str) -> bool {
    genotype.extended.iter().any(|(name, _)| name == key)
}

/// `SVCallRecord.isCarrier`: whether this sample carries the record's alternate allele.
///
/// A record with no alternate carries nothing; otherwise the genotype must declare `ECN`, a ploidy
/// of zero carries nothing, a called genotype carries when it holds an alternate, and an uncalled
/// one is judged by its copy number for the three copy-number types.
pub fn is_carrier(
    sv_type: SvType,
    has_alt: bool,
    genotype: &Genotype,
) -> Result<bool, CollapseError> {
    if !has_alt {
        return Ok(false);
    }
    if !has_extended(genotype, "ECN") {
        return Err(CollapseError::MissingExpectedCopyNumber);
    }
    let expected = genotype_int(genotype, "ECN", 0);
    if expected == 0 {
        return Ok(false);
    }
    let called = !genotype.alleles.is_empty() && genotype.alleles.iter().all(|a| !a.is_no_call());
    if called {
        return Ok(genotype
            .alleles
            .iter()
            .any(|allele| !allele.is_no_call() && !allele.is_reference()));
    }
    let copy_number = genotype_int(genotype, "CN", genotype_int(genotype, "RD_CN", expected));
    Ok(match sv_type {
        SvType::Del => copy_number < expected,
        SvType::Dup => copy_number > expected,
        SvType::Cnv => copy_number != expected,
        _ => false,
    })
}

fn carrier_count(member: &Member) -> Result<usize, CollapseError> {
    let has_alt = !member.alt_alleles().is_empty();
    let mut count = 0;
    for genotype in &member.genotypes {
        if is_carrier(member.call.sv_type, has_alt, genotype)? {
            count += 1;
        }
    }
    Ok(count)
}

/// `MathUtils.median(sorted, R_1)`: the ceil(n/2)-th smallest.
fn median_r1(sorted: &[i32]) -> i32 {
    let n = sorted.len();
    sorted[n.div_ceil(2) - 1]
}

fn java_round(value: f64) -> i32 {
    (value + 0.5).floor() as i32
}

fn evidence_rank(evidence: &[Evidence]) -> (bool, bool) {
    // Split-read first, then paired-end: `false` sorts first, so the presence is negated.
    (
        !evidence.contains(&Evidence::Sr),
        !evidence.contains(&Evidence::Pe),
    )
}

fn distance_to_all(member: &Member, starts: &[i32], ends: &[i32]) -> i64 {
    let a: i64 = starts
        .iter()
        .map(|start| (i64::from(*start) - i64::from(member.call.position_a)).abs())
        .sum();
    let b: i64 = ends
        .iter()
        .map(|end| (i64::from(*end) - i64::from(member.call.position_b)).abs())
        .sum();
    a + b
}

/// `getRepresentativeIntervalItem`'s chain as one key: quality, breakpoint evidence, carriers,
/// distance, ID, each smaller-is-better.
type RepresentativeKey = (f64, (bool, bool), i64, i64, String);

/// `collapseInterval`.
fn collapse_interval(
    precise: &[&Member],
    strategy: BreakpointSummary,
) -> Result<(i32, i32), CollapseError> {
    let example = &precise[0].call;
    if precise.len() > 1 {
        if precise.iter().any(|m| m.call.contig_a != example.contig_a) {
            return Err(CollapseError::MultipleContigsA);
        }
        if precise.iter().any(|m| m.call.contig_b != example.contig_b) {
            return Err(CollapseError::MultipleContigsB);
        }
    }
    let mut starts: Vec<i32> = precise.iter().map(|m| m.call.position_a).collect();
    let mut ends: Vec<i32> = precise.iter().map(|m| m.call.position_b).collect();
    starts.sort_unstable();
    ends.sort_unstable();
    let (start, end) = match strategy {
        BreakpointSummary::MedianStartMedianEnd => (median_r1(&starts), median_r1(&ends)),
        BreakpointSummary::MinStartMaxEnd => (starts[0], ends[ends.len() - 1]),
        BreakpointSummary::MaxStartMinEnd => (starts[starts.len() - 1], ends[0]),
        BreakpointSummary::MeanStartMeanEnd => {
            let n = starts.len() as f64;
            let sum = |values: &[i32]| values.iter().map(|v| f64::from(*v)).sum::<f64>();
            (java_round(sum(&starts) / n), java_round(sum(&ends) / n))
        }
        BreakpointSummary::Representative => {
            let chosen = if precise.len() == 1 {
                precise[0]
            } else {
                let mut best: Option<(&Member, RepresentativeKey)> = None;
                for member in precise {
                    let key = (
                        member.call.log10_p_error.unwrap_or(0.0),
                        evidence_rank(&member.call.evidence),
                        -(carrier_count(member)? as i64),
                        distance_to_all(member, &starts, &ends),
                        member.call.id.clone(),
                    );
                    let better = match &best {
                        None => true,
                        Some((_, current)) => {
                            key.0
                                .partial_cmp(&current.0)
                                .unwrap_or(std::cmp::Ordering::Equal)
                                .then(key.1.cmp(&current.1))
                                .then(key.2.cmp(&current.2))
                                .then(key.3.cmp(&current.3))
                                .then(key.4.cmp(&current.4))
                                == std::cmp::Ordering::Less
                        }
                    };
                    if better {
                        best = Some((member, key));
                    }
                }
                best.expect("a precise member").0
            };
            (chosen.call.position_a, chosen.call.position_b)
        }
    };
    Ok(if example.sv_type == SvType::Ins {
        (start, start)
    } else if example.contig_a == example.contig_b {
        (start, start.max(end))
    } else {
        (start, end)
    })
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

/// `collapseTypes`. The refusal lists the types in a `HashSet` of enum constants, whose order is
/// their identity hash codes: no fixture should reach it.
fn collapse_types(members: &[Member]) -> Result<SvType, CollapseError> {
    let mut types: Vec<SvType> = Vec::new();
    for member in members {
        if !types.contains(&member.call.sv_type) {
            types.push(member.call.sv_type);
        }
    }
    if types.len() == 1 {
        return Ok(types[0]);
    }
    if types
        .iter()
        .all(|t| matches!(t, SvType::Del | SvType::Dup | SvType::Cnv))
    {
        return Ok(SvType::Cnv);
    }
    Err(CollapseError::IncompatibleTypes {
        types: types.iter().map(|t| type_name(*t).to_string()).collect(),
    })
}

fn symbolic(name: &str) -> Allele {
    Allele::from_str(name, false).expect("a symbolic allele")
}

fn symbols(allele: &Allele) -> Vec<String> {
    allele
        .display_string()
        .replace(['<', '>'], "")
        .split(':')
        .map(str::to_string)
        .collect()
}

fn is_breakpoint_or_single_breakend(allele: &Allele) -> bool {
    let text = allele.display_string();
    let bytes = text.as_bytes();
    bytes.len() > 1
        && (text.contains('[')
            || text.contains(']')
            || bytes[0] == b'.'
            || bytes[bytes.len() - 1] == b'.')
}

/// `collapseAltAlleles`, over every member.
fn collapse_alt_alleles(
    members: &[Member],
    strategy: AltAlleleSummary,
) -> Result<Vec<Allele>, CollapseError> {
    let mut alleles: Vec<Allele> = Vec::new();
    for member in members {
        for allele in member.alt_alleles() {
            if !alleles.contains(allele) {
                alleles.push(allele.clone());
            }
        }
    }
    alleles.sort_by_key(|allele| allele.display_string());
    if alleles.len() <= 1 {
        return Ok(alleles);
    }
    for allele in &alleles {
        if allele.is_reference()
            || !allele.is_symbolic()
            || is_breakpoint_or_single_breakend(allele)
        {
            return Err(CollapseError::NonSymbolicAllele {
                allele: allele.display_string(),
            });
        }
    }
    let (mut cnv, mut multiallelic) = (0, 0);
    for allele in &alleles {
        match allele.display_string().as_str() {
            "<CNV>" => {
                cnv += 1;
                multiallelic += 1;
            }
            "<DUP>" | "<DEL>" => cnv += 1,
            _ => {}
        }
    }
    if alleles.len() == cnv {
        return Ok(if multiallelic > 0 {
            vec![symbolic("<CNV>")]
        } else {
            vec![symbolic("<DEL>"), symbolic("<DUP>")]
        });
    }
    let tokens: Vec<Vec<String>> = alleles.iter().map(symbols).collect();
    let kept: Vec<String> = match strategy {
        AltAlleleSummary::CommonSubtype => {
            let first = &tokens[0];
            let mut size = 0;
            'outer: for (i, subtype) in first.iter().enumerate() {
                for other in &tokens[1..] {
                    if i < other.len() && *subtype == other[i] {
                        size = i + 1;
                    } else {
                        break 'outer;
                    }
                }
            }
            first[..size].to_vec()
        }
        AltAlleleSummary::MostSpecificSubtype => {
            let max = tokens.iter().map(Vec::len).max().unwrap_or(0);
            let (mut index, mut size) = (0, 0);
            'outer: for i in 0..max {
                let mut subtype: Option<&String> = None;
                for (j, other) in tokens.iter().enumerate() {
                    if i < other.len() {
                        match subtype {
                            None => {
                                subtype = Some(&other[i]);
                                index = j;
                                size = i + 1;
                            }
                            Some(seen) if *seen != other[i] => {
                                index = j;
                                size = i;
                                break 'outer;
                            }
                            _ => {}
                        }
                    }
                }
            }
            tokens[index][..size].to_vec()
        }
    };
    Ok(vec![symbolic(&format!("<{}>", kept.join(":")))])
}

/// `getRepresentativeGenotype`'s chain, as one sort key: larger is better.
fn genotype_key(genotype: &Genotype) -> (i64, i64, i32, i64, i32, i64, i32) {
    let non_ref = genotype
        .alleles
        .iter()
        .filter(|a| !a.is_reference() && !a.is_no_call())
        .count() as i64;
    let called = genotype.alleles.iter().filter(|a| !a.is_no_call()).count() as i64;
    let expected = genotype_int(genotype, "ECN", 0);
    let copy_number = genotype_int(genotype, "CN", 0);
    let deletion_expected = genotype_int(genotype, "ECN", 0);
    let is_deletion = genotype_int(genotype, "CN", deletion_expected) < deletion_expected;
    (
        non_ref.min(1),
        called,
        genotype_int(genotype, "GQ", 0),
        // Fewer non-reference alleles is better here: the comparator is reversed.
        -non_ref,
        genotype_int(genotype, "CNQ", 0),
        // A copy number closer to the expected one is better.
        -i64::from((expected - copy_number).abs()),
        // `genotypeDelOverDupComparator` answers 1 for a deletion against a non-deletion, so the
        // deletion is the greater and wins.
        if is_deletion { 1 } else { 0 },
    )
}

/// `collapseAllGenotypes` then `harmonizeAltAlleles`.
fn collapse_genotypes(
    members: &[Member],
    reference: &Allele,
    alt_alleles: &[Allele],
) -> Result<Vec<Genotype>, CollapseError> {
    let mut samples: Vec<String> = Vec::new();
    for member in members {
        for genotype in &member.genotypes {
            if !samples.contains(&genotype.sample_name) {
                samples.push(genotype.sample_name.clone());
            }
        }
    }
    let mut collapsed: Vec<Genotype> = Vec::new();
    for sample in &samples {
        let mut best: Option<&Genotype> = None;
        for member in members {
            for genotype in member.genotypes.iter().filter(|g| g.sample_name == *sample) {
                // `maxBy` keeps the first of equals, so a later one must be strictly greater.
                if best.is_none_or(|current| genotype_key(genotype) > genotype_key(current)) {
                    best = Some(genotype);
                }
            }
        }
        let mut genotype = best.expect("a genotype for the sample").clone();
        genotype.alleles = genotype
            .alleles
            .iter()
            .map(|allele| {
                if allele.is_reference() {
                    reference.clone()
                } else {
                    allele.clone()
                }
            })
            .collect();
        collapsed.push(genotype);
    }
    let genotype_alts: Vec<&Allele> = collapsed
        .iter()
        .flat_map(|g| g.alleles.iter())
        .filter(|a| !a.is_no_call() && !a.is_reference())
        .collect();
    if genotype_alts.iter().all(|a| alt_alleles.contains(a)) {
        return Ok(collapsed);
    }
    if alt_alleles.len() != 1 {
        return Err(CollapseError::SubtypedMultiallelic);
    }
    let new_alt = &alt_alleles[0];
    for genotype in &mut collapsed {
        genotype.alleles = genotype
            .alleles
            .iter()
            .map(|a| {
                if !a.is_no_call() && !a.is_reference() {
                    new_alt.clone()
                } else {
                    a.clone()
                }
            })
            .collect();
    }
    Ok(collapsed)
}

/// `constructBndAllele`.
fn breakend_allele(
    strand_a: Option<bool>,
    strand_b: Option<bool>,
    contig_b: &str,
    position_b: i32,
    reference: &Allele,
) -> Result<Allele, CollapseError> {
    let strand_a = strand_a.ok_or(CollapseError::NullBreakendStrand { first: true })?;
    let strand_b = strand_b.ok_or(CollapseError::NullBreakendStrand { first: false })?;
    let bracket = if strand_b { "]" } else { "[" };
    let base = reference.base_string();
    let text = if strand_a {
        format!("{base}{bracket}{contig_b}:{position_b}{bracket}")
    } else {
        format!("{bracket}{contig_b}:{position_b}{bracket}{base}")
    };
    Ok(Allele::from_str(&text, false).expect("a breakend allele"))
}

/// `CanonicalSVCollapser.collapse`.
///
/// `reference_base` answers the base at a 1-based position, or `None` off the reference, which is
/// `collapseRefAlleles`' refusal.
pub fn collapse(
    members: &[Member],
    breakpoints: BreakpointSummary,
    alternates: AltAlleleSummary,
    flags: FlagFieldLogic,
    reference_base: &mut dyn FnMut(&str, i32) -> Option<u8>,
) -> Result<Collapsed, CollapseError> {
    let all_depth = members.iter().all(Member::is_depth_only);
    let precise: Vec<&Member> = members
        .iter()
        .filter(|m| all_depth || !m.is_depth_only())
        .collect();
    let (start, end) = collapse_interval(&precise, breakpoints)?;
    // `getRepresentativeRecord`: sorted by ID, then the first of the closest.
    let mut by_id = precise.clone();
    by_id.sort_by(|a, b| a.call.id.cmp(&b.call.id));
    let distance = |m: &Member| {
        (i64::from(m.call.position_a) - i64::from(start)).abs()
            + (i64::from(m.call.position_b) - i64::from(end)).abs()
    };
    let representative = by_id
        .iter()
        .min_by_key(|m| distance(m))
        .expect("a precise member");
    let sv_type = collapse_types(members)?;
    let length = match sv_type {
        SvType::Ins | SvType::Cpx => representative.call.length,
        SvType::Bnd | SvType::Ctx => None,
        _ => Some(end - start + 1),
    };
    let mut algorithms: Vec<String> = Vec::new();
    for member in members {
        for algorithm in &member.call.algorithms {
            if !algorithms.contains(algorithm) {
                algorithms.push(algorithm.clone());
            }
        }
    }
    algorithms.sort();

    let mut attributes: Vec<(String, Value)> = representative
        .call
        .attributes
        .iter()
        .filter(|(key, _)| !FLAG_TYPE_INFO_FIELDS.contains(&key.as_str()))
        .cloned()
        .collect();
    let mut ids: Vec<String> = members.iter().map(|m| m.call.id.clone()).collect();
    ids.sort();
    attributes.retain(|(key, _)| key != CLUSTER_MEMBER_IDS_KEY);
    attributes.push((
        CLUSTER_MEMBER_IDS_KEY.to_string(),
        Value::List(ids.into_iter().map(Value::Str).collect()),
    ));
    for key in FLAG_TYPE_INFO_FIELDS {
        let set = |m: &Member| {
            m.call
                .attributes
                .iter()
                .any(|(name, value)| name == key && *value == Value::Bool(true))
        };
        let on = match flags {
            FlagFieldLogic::And => members.iter().all(set),
            FlagFieldLogic::Or => members.iter().any(set),
            FlagFieldLogic::AlwaysFalse => false,
        };
        if on {
            attributes.push((key.to_string(), Value::Bool(true)));
        }
    }

    let (strand_a, strand_b) = if sv_type == SvType::Cnv {
        (None, None)
    } else {
        (representative.call.strand_a, representative.call.strand_b)
    };
    let contig = representative.call.contig_a.clone();
    let base = reference_base(&contig, start).ok_or(CollapseError::InvalidReferenceLocus {
        contig: contig.clone(),
        position: start,
    })?;
    let reference = Allele::create(&[base], true).expect("a reference base");
    let alt_alleles = if sv_type == SvType::Bnd {
        vec![breakend_allele(
            strand_a,
            strand_b,
            &representative.call.contig_b,
            end,
            &reference,
        )?]
    } else {
        collapse_alt_alleles(members, alternates)?
    };
    let genotypes = collapse_genotypes(members, &reference, &alt_alleles)?;
    let mut alleles = vec![reference];
    alleles.extend(alt_alleles);
    let mut filters: Vec<String> = Vec::new();
    for member in members {
        for filter in &member.filters {
            if !filters.contains(filter) {
                filters.push(filter.clone());
            }
        }
    }
    filters.sort();
    let log10_p_error = if members.len() == 1 {
        members[0].call.log10_p_error
    } else {
        None
    };

    let call = SvCallRecord {
        id: representative.call.id.clone(),
        contig_a: contig,
        position_a: start,
        strand_a,
        contig_b: representative.call.contig_b.clone(),
        position_b: end,
        strand_b,
        sv_type,
        cpx_subtype: representative.call.cpx_subtype.clone(),
        cpx_intervals: representative.call.cpx_intervals.clone(),
        length,
        evidence: representative.call.evidence.clone(),
        algorithms,
        attributes,
        log10_p_error,
    };
    Ok(Collapsed {
        call,
        alleles,
        genotypes,
        filters,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn call(id: &str, start: i32, end: i32, algorithms: &[&str]) -> SvCallRecord {
        SvCallRecord {
            id: id.to_string(),
            contig_a: "chr1".to_string(),
            position_a: start,
            strand_a: Some(true),
            contig_b: "chr1".to_string(),
            position_b: end,
            strand_b: Some(false),
            sv_type: SvType::Del,
            cpx_subtype: None,
            cpx_intervals: Vec::new(),
            length: Some(end - start + 1),
            evidence: Vec::new(),
            algorithms: algorithms.iter().map(|a| a.to_string()).collect(),
            attributes: Vec::new(),
            log10_p_error: None,
        }
    }

    fn member(call: SvCallRecord, gt: &[&str]) -> Member {
        let reference = Allele::create(b"N", true).unwrap();
        let del = symbolic("<DEL>");
        let genotype = Genotype::new(
            "s1",
            gt.iter()
                .map(|a| match *a {
                    "0" => reference.clone(),
                    "1" => del.clone(),
                    _ => Allele::no_call(),
                })
                .collect(),
        );
        Member {
            call,
            alleles: vec![reference, del],
            genotypes: vec![genotype],
            filters: Vec::new(),
        }
    }

    #[test]
    fn depth_only_members_do_not_move_the_breakpoints() {
        let members = vec![
            member(call("b", 1000, 2000, &["manta"]), &["0", "1"]),
            member(call("a", 900, 2500, &["depth"]), &["0", "0"]),
        ];
        let collapsed = collapse(
            &members,
            BreakpointSummary::MedianStartMedianEnd,
            AltAlleleSummary::CommonSubtype,
            FlagFieldLogic::Or,
            &mut |_, _| Some(b'a'),
        )
        .unwrap();
        assert_eq!(
            (collapsed.call.position_a, collapsed.call.position_b),
            (1000, 2000)
        );
        assert_eq!(collapsed.call.id, "b");
        assert_eq!(collapsed.call.algorithms, vec!["depth", "manta"]);
        assert_eq!(
            collapsed.alleles[0].display_string(),
            "A",
            "htsjdk upper-cases a base allele"
        );
        // The carrier genotype wins, and its reference allele is the new one.
        assert_eq!(collapsed.genotypes[0].alleles[1].display_string(), "<DEL>");
    }

    #[test]
    fn the_median_is_the_lower_middle_value() {
        assert_eq!(median_r1(&[1, 2, 3, 4]), 2);
        assert_eq!(median_r1(&[1, 2, 3]), 2);
        assert_eq!(median_r1(&[5]), 5);
    }

    #[test]
    fn subtyped_alleles_collapse_to_their_common_prefix() {
        let mut first = member(call("a", 100, 200, &["manta"]), &["0", "1"]);
        first.alleles[1] = symbolic("<INS:ME:ALU>");
        let mut second = member(call("b", 100, 200, &["manta"]), &["0", "1"]);
        second.alleles[1] = symbolic("<INS:ME:LINE1>");
        let common = collapse_alt_alleles(
            &[first.clone(), second.clone()],
            AltAlleleSummary::CommonSubtype,
        )
        .unwrap();
        assert_eq!(common[0].display_string(), "<INS:ME>");
        let specific =
            collapse_alt_alleles(&[first, second], AltAlleleSummary::MostSpecificSubtype).unwrap();
        assert_eq!(specific[0].display_string(), "<INS:ME>");
    }
}

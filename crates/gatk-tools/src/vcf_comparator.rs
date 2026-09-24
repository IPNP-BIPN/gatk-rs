//! `VCFComparator`: what counts as two VCFs disagreeing, and where the first disagreement is found.
//!
//! Ported from `org.broadinstitute.hellbender.tools.walkers.variantutils.VCFComparator` and the
//! grouping of `org.broadinstitute.hellbender.engine.MultiVariantWalkerGroupedByOverlap` in GATK
//! 4.6.2.0.
//!
//! The tool's output is its exceptions. It walks two files merged into one stream, groups the
//! records that overlap, and throws on the first difference inside a group that no argument
//! tolerates. Four things decide what a run answers, and none of them is the comparison alone.
//!
//! # A complaint names the record AFTER the group it is about
//!
//! A group is compared when the next record does not overlap it, and that happens inside the
//! walker's `apply` for the NEXT record. `MultiVariantWalker.traverse` wraps whatever `apply`
//! throws in a `GATKException` naming that record, so a difference at 100 reaches the user as
//! `Exception thrown at chr1:101 [VC expected @ chr1:101-199 ...]`, and the comparison's own
//! message is only the cause. The LAST group is compared after the traversal, outside the wrap,
//! and its complaint is the plain user error. [`Stopped::at`] carries which of the two it was.
//!
//! # The merge breaks a tie by who spoke last
//!
//! The two files are merged by a priority queue over their next records, ordered by contig and
//! start and nothing else. On a tie the input that did NOT just supply a record goes first, and at
//! the very start the input named first on the command line does. Which record of a tied pair
//! comes first decides nothing in a group, but it does decide which record a wrap names.
//!
//! # Two of the tolerances are one-sided or unreachable
//!
//! The allele check runs only when EXPECTED carries an allele actual lacks
//! (`actualHasNewAlleles(expected, actual)`), so an allele added to actual is never checked. And
//! `alleleNumberIsDifferent` is assigned from a comparison that throws rather than returning
//! false, so it is never true and `--mute-acceptable-diffs` mutes only an inbreeding difference or
//! a low GQ.
//!
//! # Trimming a record rewrites its annotations
//!
//! Before comparing, each record is cut down to the alleles its genotypes call, plus `<NON_REF>`
//! unless `--ignore-non-ref-data` says otherwise. When that changes the allele count of a record
//! that is not hom-ref, its INFO map is rebuilt the way `ReblockGVCF` rebuilds one: only the keys
//! an annotation claims survive, and a genotype count is added. The allele-specific annotations
//! are subset allele by allele there, which is not ported, and a record carrying one is
//! [`Failure::Limitation`].

use gatk_engine::interval::SimpleInterval;
use gatk_engine::java_format::format_decimals;
use gatk_engine::java_hash::JavaHashMap;
use gatk_engine::subset_alleles::{subset_alleles, AssignmentMethod, Genotype as EngineGenotype};
use gatk_engine::tsv_table::java_double_to_string;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::genotypes_context::GenotypesContext;
use htsjdk_vcf::variant::{Genotype, Value, VariantContext};

/// What the tool refuses about its inputs, before a record is read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InputError {
    WrongNumberOfInputs,
    NoExpectedInput,
}

impl InputError {
    pub fn message(&self) -> String {
        match self {
            InputError::WrongNumberOfInputs => {
                "Bad input: VCFComparator expects exactly two inputs -- one actual and one \
                 expected."
                    .to_string()
            }
            InputError::NoExpectedInput => {
                "Bad input: Tool requires exactly one expected input file".to_string()
            }
        }
    }
}

/// `onTraversalStart`'s two checks, in its own order: the COUNT first, then the tag.
pub fn check_inputs(names: &[String]) -> Result<(), InputError> {
    if names.len() != 2 {
        return Err(InputError::WrongNumberOfInputs);
    }
    if names.iter().filter(|name| *name == "expected").count() != 1 {
        return Err(InputError::NoExpectedInput);
    }
    Ok(())
}

/// Every argument the comparison reads, with the tool's defaults.
#[derive(Debug, Clone, PartialEq)]
pub struct Options {
    pub warn_on_errors: bool,
    pub finish_before_failing: bool,
    pub default_ploidy: usize,
    pub ignore_quals: bool,
    pub qual_change_allowed: f64,
    pub inbreeding_coeff_change_allowed: f64,
    pub good_qual_threshold: f64,
    pub dp_change_allowed: i32,
    pub ranksum_change_allowed: f64,
    pub likelihood_change_allowed: i32,
    pub ignore_non_ref_data: bool,
    pub ignore_annotations: bool,
    pub ignore_genotype_annotations: bool,
    pub ignore_genotype_phasing: bool,
    pub ignore_filters: bool,
    pub ignore_attributes: Vec<String>,
    pub positions_only: bool,
    pub allow_new_stars: bool,
    pub allow_extra_alleles: bool,
    pub allow_missing_stars: bool,
    pub ignore_star_attributes: bool,
    pub allow_nan_mismatch: bool,
    pub mute_acceptable_diffs: bool,
    pub ignore_hom_ref_attributes: bool,
    pub ignore_dbsnp_ids: bool,
    pub ignore_gq0: bool,
    pub ignore_some_multi_allelics: bool,
    /// `--annotations-to-keep`, and the two arguments that decide which annotations the engine
    /// holds, all three read only when a trim subsets a record's annotations.
    pub annotations_to_keep: Vec<String>,
    pub enable_all_annotations: bool,
    pub disable_tool_default_annotations: bool,
    /// `MultiVariantWalkerGroupedByOverlap`'s two record filters.
    pub ignore_variants_starting_outside_interval: bool,
    pub ignore_reference_blocks: bool,
    /// `--ref-padding`, the window every group's reference context is widened by.
    pub reference_padding: i32,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            warn_on_errors: false,
            finish_before_failing: false,
            default_ploidy: 2,
            ignore_quals: false,
            qual_change_allowed: 0.001,
            inbreeding_coeff_change_allowed: 0.001,
            good_qual_threshold: 100.0,
            dp_change_allowed: 0,
            ranksum_change_allowed: 0.0,
            likelihood_change_allowed: 0,
            ignore_non_ref_data: false,
            ignore_annotations: false,
            ignore_genotype_annotations: false,
            ignore_genotype_phasing: false,
            ignore_filters: false,
            ignore_attributes: Vec::new(),
            positions_only: false,
            allow_new_stars: false,
            allow_extra_alleles: false,
            allow_missing_stars: false,
            ignore_star_attributes: false,
            allow_nan_mismatch: false,
            mute_acceptable_diffs: false,
            ignore_hom_ref_attributes: false,
            ignore_dbsnp_ids: false,
            ignore_gq0: false,
            ignore_some_multi_allelics: false,
            annotations_to_keep: Vec::new(),
            enable_all_annotations: false,
            disable_tool_default_annotations: false,
            ignore_variants_starting_outside_interval: false,
            ignore_reference_blocks: false,
            reference_padding: 1,
        }
    }
}

/// Why a run stopped.
#[derive(Debug, Clone, PartialEq)]
pub enum Failure {
    /// A `UserException`, by its message.
    User(String),
    /// Anything else the reference throws, by class and message: a `NumberFormatException` from an
    /// annotation that does not parse, an index past the end of a short list.
    Runtime {
        class: &'static str,
        message: String,
    },
    /// A path this port does not carry, by its own message.
    Limitation(String),
}

/// A failure, and the record the traversal was on when it happened.
#[derive(Debug, Clone, PartialEq)]
pub struct Stopped {
    pub failure: Failure,
    /// `(input, record)`: the record whose `apply` threw, which the walker names when it wraps the
    /// failure. `None` for a failure after the traversal, which reaches the user unwrapped.
    pub at: Option<(usize, usize)>,
}

/// What a run that did not stop leaves behind.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Finished {
    /// The warnings, in order, as the log prints them.
    pub warnings: Vec<String>,
}

/// One input: the name its tag gave it, which every record carries as its source, the samples of
/// its header and its records.
#[derive(Debug, Clone)]
pub struct Input {
    pub name: String,
    pub samples: Vec<String>,
    pub records: Vec<VariantContext>,
}

/// `Allele.NON_REF_ALLELE`.
const NON_REF: &str = "<NON_REF>";
/// `Allele.SPAN_DEL`.
const SPAN_DEL: &str = "*";

const INBREEDING_COEFFICIENT_KEY: &str = "InbreedingCoeff";
const AS_INBREEDING_COEFFICIENT_KEY: &str = "AS_InbreedingCoeff";
const ALLELE_NUMBER_KEY: &str = "AN";
const ALLELE_COUNT_KEY: &str = "AC";
const MLE_ALLELE_COUNT_KEY: &str = "MLEAC";
const AS_SB_TABLE_KEY: &str = "AS_SB_TABLE";
const RAW_GENOTYPE_COUNT_KEY: &str = "RAW_GT_COUNT";
const RAW_MAPPING_QUALITY_WITH_DEPTH_KEY: &str = "RAW_MQandDP";
const AS_QUAL_BY_DEPTH_KEY: &str = "AS_QD";
const QUAL_BY_DEPTH_KEY: &str = "QD";
const DEPTH_KEY: &str = "DP";
const HAPLOTYPE_CALLER_PHASING_ID_KEY: &str = "PID";
const AS_VARIANT_DEPTH_KEY: &str = "AS_VarDP";
const RAW_QUAL_APPROX_KEY: &str = "QUALapprox";
const RAW_RMS_MAPPING_QUALITY_DEPRECATED: &str = "RAW_MQ";

/// The INFO keys the standard annotations claim, annotation by annotation in class-name order:
/// `BaseQualityRankSumTest`, `ChromosomeCounts`, `Coverage`, `ExcessHet`, `FisherStrand`,
/// `InbreedingCoeff`, `MappingQualityRankSumTest`, `QualByDepth`, `RMSMappingQuality` (with its raw
/// key), `ReadPosRankSumTest` and `StrandOddsRatio`.
const STANDARD_ANNOTATION_KEYS: [&str; 14] = [
    "BaseQRankSum",
    "AN",
    "AC",
    "AF",
    "DP",
    "ExcessHet",
    "FS",
    "InbreedingCoeff",
    "MQRankSum",
    "QD",
    "MQ",
    "RAW_MQandDP",
    "ReadPosRankSum",
    "SOR",
];

/// `ReblockGVCF.infoFieldAnnotationKeyNamesToRemove`.
const REMOVED_ANNOTATION_KEYS: [&str; 8] = [
    "GVCFBlock",
    "HaplotypeScore",
    "InbreedingCoeff",
    "MLEAC",
    "MLEAF",
    "ExcessHet",
    "AS_InbreedingCoeff",
    "DS",
];

/// `annotationsThatVaryWithNoCalls`, which `--mute-acceptable-diffs` would forgive when the allele
/// number differs, and it never does.
const VARY_WITH_NO_CALLS: [&str; 15] = [
    "AN",
    "InbreedingCoeff",
    "AS_InbreedingCoeff",
    "ExcessHet",
    "MLEAC",
    "MLEAF",
    "AF",
    "DP",
    "QD",
    "AS_QD",
    "VQSLOD",
    "AS_VQSLOD",
    "AS_FilterStatus",
    "culprit",
    "AS_culprit",
];

fn is_non_ref(allele: &Allele) -> bool {
    allele.is_symbolic() && !allele.is_reference() && allele.display_string() == NON_REF
}

fn is_span_del(allele: &Allele) -> bool {
    !allele.is_reference() && !allele.is_no_call() && allele.display_string() == SPAN_DEL
}

/// `GATKVariantContextUtils.isConcreteAlt`.
fn is_concrete_alt(allele: &Allele) -> bool {
    !allele.is_reference() && !allele.is_symbolic() && !is_span_del(allele) && !allele.is_no_call()
}

/// `Allele.toString()`: a no-call is `.`, and a reference allele carries a `*`.
fn allele_string(allele: &Allele) -> String {
    let text = if allele.is_no_call() {
        ".".to_string()
    } else {
        allele.display_string()
    };
    if allele.is_reference() {
        format!("{text}*")
    } else {
        text
    }
}

/// `List<Allele>.toString()`.
fn alleles_string(alleles: &[Allele]) -> String {
    format!(
        "[{}]",
        alleles
            .iter()
            .map(allele_string)
            .collect::<Vec<_>>()
            .join(", ")
    )
}

/// `Allele.compareTo`: the reference first, then the bases as written, shorter first on a tie.
fn allele_order(a: &Allele, b: &Allele) -> std::cmp::Ordering {
    use std::cmp::Ordering;
    match (a.is_reference(), b.is_reference()) {
        (true, false) => Ordering::Less,
        (false, true) => Ordering::Greater,
        _ => {
            let bases = |allele: &Allele| {
                if allele.is_no_call() {
                    Vec::new()
                } else {
                    allele.display_string().into_bytes()
                }
            };
            bases(a).cmp(&bases(b))
        }
    }
}

/// `Genotype.getGenotypeString(false)`: sorted unless phased, each allele by its `toString`.
fn genotype_string(genotype: &Genotype) -> String {
    if genotype.alleles.is_empty() {
        return "NA".to_string();
    }
    let separator = if genotype.phased { "|" } else { "/" };
    let mut alleles = genotype.alleles.clone();
    if !genotype.phased {
        alleles.sort_by(allele_order);
    }
    alleles
        .iter()
        .map(allele_string)
        .collect::<Vec<_>>()
        .join(separator)
}

/// `Object.toString()` of a decoded attribute: a String as itself, a list in brackets.
fn value_string(value: &Value) -> String {
    match value {
        Value::Missing => ".".to_string(),
        Value::Bool(flag) => flag.to_string(),
        Value::Int(number) => number.to_string(),
        Value::Double(number) => java_double_to_string(*number),
        Value::Str(text) => text.clone(),
        Value::List(items) => format!(
            "[{}]",
            items
                .iter()
                .map(value_string)
                .collect::<Vec<_>>()
                .join(", ")
        ),
    }
}

/// `int[]` printed through `Arrays.toString`.
fn ints_string(values: &[i32]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(i32::to_string)
            .collect::<Vec<_>>()
            .join(", ")
    )
}

/// `Genotype.toString()`.
fn genotype_to_string(genotype: &Genotype) -> String {
    let join = |values: &[i32]| {
        values
            .iter()
            .map(i32::to_string)
            .collect::<Vec<_>>()
            .join(",")
    };
    let mut text = format!("[{} {}", genotype.sample_name, genotype_string(genotype));
    if let Some(gq) = genotype.gq {
        text.push_str(&format!(" GQ {gq}"));
    }
    if let Some(dp) = genotype.dp {
        text.push_str(&format!(" DP {dp}"));
    }
    if let Some(ad) = &genotype.ad {
        text.push_str(&format!(" AD {}", join(ad)));
    }
    if let Some(pl) = &genotype.pl {
        text.push_str(&format!(" PL {}", join(pl)));
    }
    if let Some(filters) = &genotype.filters {
        text.push_str(&format!(" FT {filters}"));
    }
    let mut extended: Vec<(String, String)> = genotype
        .extended
        .iter()
        .map(|(key, value)| (key.clone(), value_string(value)))
        .collect();
    extended.sort();
    if !extended.is_empty() {
        text.push_str(&format!(
            " {{{}}}",
            extended
                .iter()
                .map(|(key, value)| format!("{key}={value}"))
                .collect::<Vec<_>>()
                .join(", ")
        ));
    }
    text.push(']');
    text
}

/// `getGQ()`, which is -1 when there is none.
fn gq(genotype: &Genotype) -> i32 {
    genotype.gq.unwrap_or(-1)
}

fn number_format(text: &str) -> Failure {
    Failure::Runtime {
        class: "java.lang.NumberFormatException",
        message: if text.trim().is_empty() {
            "empty String".to_string()
        } else {
            format!("For input string: \"{text}\"")
        },
    }
}

/// `Double.parseDouble`: surrounding whitespace ignored, `NaN` and `Infinity` spelled as Java
/// spells them, an optional `d` or `f` suffix, and nothing Rust accepts that Java would not.
fn parse_double(text: &str) -> Result<f64, Failure> {
    let trimmed = text.trim_matches(|c: char| c <= ' ');
    let (sign, body) = match trimmed.strip_prefix('-') {
        Some(rest) => (-1.0, rest),
        None => (1.0, trimmed.strip_prefix('+').unwrap_or(trimmed)),
    };
    match body {
        "NaN" => return Ok(f64::NAN),
        "Infinity" => return Ok(sign * f64::INFINITY),
        _ => {}
    }
    let body = body.strip_suffix(['d', 'D', 'f', 'F']).unwrap_or(body);
    let plausible = !body.is_empty()
        && body.chars().any(|c| c.is_ascii_digit())
        && body
            .chars()
            .all(|c| c.is_ascii_digit() || matches!(c, '.' | 'e' | 'E' | '+' | '-'));
    if !plausible {
        return Err(number_format(text));
    }
    body.parse::<f64>()
        .map(|value| sign * value)
        .map_err(|_| number_format(text))
}

/// `Integer.parseInt`: an optional sign and decimal digits, nothing else.
fn parse_int(text: &str) -> Result<i32, Failure> {
    let digits = text.strip_prefix(['-', '+']).unwrap_or(text);
    if digits.is_empty() || !digits.chars().all(|c| c.is_ascii_digit()) {
        return Err(number_format(text));
    }
    text.parse::<i32>().map_err(|_| number_format(text))
}

fn index_out_of_bounds(index: usize, length: usize) -> Failure {
    Failure::Runtime {
        class: "java.lang.IndexOutOfBoundsException",
        message: format!("Index {index} out of bounds for length {length}"),
    }
}

fn null_pointer() -> Failure {
    Failure::Runtime {
        class: "java.lang.NullPointerException",
        message: String::new(),
    }
}

/// A record as the comparison reads it: the file's record, the input it came from, and its INFO
/// attributes in the order a `HashMap` of them iterates.
#[derive(Debug, Clone)]
struct Record {
    source: String,
    vc: VariantContext,
    attributes: Vec<(String, Value)>,
}

/// The iteration order of a `HashMap` built by putting `entries` in order, at `capacity` buckets.
fn hash_order(
    entries: &[(String, Value)],
    capacity: Option<usize>,
) -> Result<Vec<(String, Value)>, Failure> {
    let mut map: JavaHashMap<String, Value> = match capacity {
        Some(capacity) => JavaHashMap::with_capacity(capacity),
        None => JavaHashMap::new(),
    };
    for (key, value) in entries {
        map.insert(key.clone(), value.clone());
    }
    map.check().map_err(|_| {
        Failure::Limitation(
            "A record with more INFO attributes in one hash bucket than this port has measured \
             the order of. This message is the port's own and not GATK's."
                .to_string(),
        )
    })?;
    Ok(map
        .iter()
        .map(|(key, value)| (key.clone(), value.clone()))
        .collect())
}

/// The bucket count a map copied from `size` entries into an empty map starts at: the entry count
/// over the load factor, plus one, rounded up to a power of two, as measured.
fn copied_capacity(size: usize) -> usize {
    ((size as f32 / 0.75f32) + 1.0f32) as usize
}

impl Record {
    /// A record as the codec decoded it. The INFO field is parsed into one map and handed to the
    /// builder, which copies it into another; every later builder copies it again at the same
    /// size, which keeps the order, so this is the order every comparison iterates in.
    fn decoded(source: &str, vc: &VariantContext) -> Result<Record, Failure> {
        let parsed = hash_order(&vc.attributes, None)?;
        let attributes = hash_order(&parsed, Some(copied_capacity(parsed.len())))?;
        Ok(Record {
            source: source.to_string(),
            vc: vc.clone(),
            attributes,
        })
    }

    fn attribute(&self, key: &str) -> Option<&Value> {
        self.attributes
            .iter()
            .find(|(name, _)| name == key)
            .map(|(_, value)| value)
    }

    fn contig(&self) -> &str {
        &self.vc.contig
    }

    fn start(&self) -> i64 {
        self.vc.start
    }

    fn end(&self) -> i64 {
        self.vc.stop
    }

    fn alleles(&self) -> &[Allele] {
        &self.vc.alleles
    }

    fn reference(&self) -> &Allele {
        &self.vc.alleles[0]
    }

    fn alternates(&self) -> &[Allele] {
        &self.vc.alleles[1..]
    }

    fn genotypes(&self) -> Vec<Genotype> {
        self.vc.genotypes.iter().cloned().collect()
    }

    fn genotype(&self, index: usize) -> Result<Genotype, Failure> {
        let genotypes = self.genotypes();
        genotypes
            .get(index)
            .cloned()
            .ok_or_else(|| index_out_of_bounds(index, genotypes.len()))
    }

    fn qual(&self) -> f64 {
        self.vc.phred_scaled_qual()
    }

    /// `Locatable.overlaps`.
    fn overlaps(&self, other: &Record) -> bool {
        self.contig() == other.contig()
            && self.start() <= other.end()
            && other.start() <= self.end()
    }

    /// `getAttributeAsInt(key, 0)`: the default for an absent key or the missing-value constant.
    fn attribute_as_int(&self, key: &str) -> Result<i32, Failure> {
        match self.attribute(key) {
            None | Some(Value::Missing) => Ok(0),
            Some(Value::List(_)) => Err(Failure::Runtime {
                class: "java.lang.ClassCastException",
                message: String::new(),
            }),
            Some(value) => parse_int(&value_string(value)),
        }
    }
}

/// `GATKVariantContextUtils.isAlleleInList`.
fn is_allele_in_list(
    reference1: &Allele,
    alternate1: &Allele,
    reference2: &Allele,
    alternates: &[Allele],
) -> Result<bool, Failure> {
    if reference1 == reference2 {
        return Ok(alternates.contains(alternate1));
    }
    // `determineReferenceAllele`: the longer one, and two of one length that differ are refused.
    let length = |allele: &Allele| {
        if allele.is_symbolic() {
            0
        } else {
            allele.len()
        }
    };
    if length(reference1) == length(reference2) {
        return Err(Failure::Runtime {
            class: "java.lang.IllegalStateException",
            message: format!(
                "The provided reference alleles do not appear to represent the same position, {} \
                 vs. {}",
                allele_string(reference1),
                allele_string(reference2)
            ),
        });
    }
    let extend = |allele: &Allele, extra: &[u8]| -> Allele {
        let mut bases = allele.display_string().into_bytes();
        bases.extend_from_slice(extra);
        Allele::create(&bases, allele.is_reference()).unwrap_or_else(|_| allele.clone())
    };
    // `createAlleleMapping`: an alternate that can be extended is, the spanning deletion maps to
    // itself, and a symbolic allele is not in the map at all.
    let mapped = |common: &Allele, input_reference: &Allele, input: &[Allele]| -> Vec<Allele> {
        let extra = &common.display_string().into_bytes()[length(input_reference)..];
        input
            .iter()
            .filter_map(|allele| {
                if is_span_del(allele) {
                    Some(allele.clone())
                } else if allele.is_reference() || allele.is_symbolic() {
                    None
                } else {
                    Some(extend(allele, extra))
                }
            })
            .collect()
    };
    if length(reference1) > length(reference2) {
        Ok(mapped(reference1, reference2, alternates).contains(alternate1))
    } else {
        let mapping = mapped(reference2, reference1, std::slice::from_ref(alternate1));
        Ok(mapping
            .first()
            .is_some_and(|allele| alternates.contains(allele)))
    }
}

/// `actualHasNewAlleles(a, b)`: does `a` carry an alternate `b` does not?
fn has_new_alleles(a: &Record, b: &Record) -> Result<bool, Failure> {
    for allele in a.alternates() {
        if !is_allele_in_list(a.reference(), allele, b.reference(), b.alternates())? {
            return Ok(true);
        }
    }
    Ok(false)
}

fn has_span_del(alleles: &[Allele]) -> bool {
    alleles.iter().any(is_span_del)
}

/// `QualByDepth.getDepth(genotypes, null)`.
fn depth_for_qual(genotypes: &[Genotype]) -> i32 {
    let mut depth: i32 = 0;
    let mut restricted: i32 = 0;
    for genotype in genotypes {
        if !genotype.is_het() && !genotype.is_hom_var() {
            continue;
        }
        if let Some(ad) = &genotype.ad {
            let total: i64 = ad.iter().map(|value| i64::from(*value)).sum();
            let total = total as i32;
            if total != 0 {
                if total - ad.first().copied().unwrap_or(0) > 1 {
                    restricted += total;
                }
                depth += total;
                continue;
            }
        }
        if let Some(dp) = genotype.dp {
            depth += dp;
        }
    }
    if restricted > 0 {
        restricted
    } else {
        depth
    }
}

/// `AS_QualByDepth.getAlleleDepths(genotypes)` summed, and 0 where no genotype has an AD.
fn allele_depth_sum(genotypes: &[Genotype]) -> Result<i32, Failure> {
    let Some(alleles) = genotypes
        .iter()
        .find_map(|genotype| genotype.ad.as_ref().map(Vec::len))
    else {
        return Ok(0);
    };
    let mut depths = vec![0i32; alleles];
    for genotype in genotypes {
        if !genotype.is_het() && !genotype.is_hom_var() {
            continue;
        }
        if let Some(ad) = &genotype.ad {
            let total: i64 = ad.iter().map(|value| i64::from(*value)).sum();
            if total as i32 - ad.first().copied().unwrap_or(0) > 1 {
                for (index, value) in ad.iter().enumerate() {
                    let slot = depths
                        .get_mut(index)
                        .ok_or_else(|| index_out_of_bounds(index, alleles))?;
                    *slot += value;
                }
            }
        }
    }
    Ok(depths.iter().sum())
}

/// `qualByDepthWillHaveJitter`. `fixTooHighQD` replaces anything from 35 up with a random draw, so
/// the comparison with the estimate it was given is true there, and true for a NaN.
fn qual_by_depth_will_have_jitter(expected_qual: f64, expected_depth: i32, as_ad: i32) -> bool {
    let estimate = expected_qual / f64::from(expected_depth);
    estimate.partial_cmp(&35.0) != Some(std::cmp::Ordering::Less) || estimate > 34.9 || as_ad == 0
}

fn qual_by_depth_difference_is_acceptable(
    actual: f64,
    expected: f64,
    expected_qual: f64,
    expected_depth: i32,
    as_ad: i32,
) -> bool {
    let difference = (expected - actual).abs();
    let relative = difference / expected;
    expected > 25.0
        || relative < 0.01
        || difference < 0.5
        || qual_by_depth_will_have_jitter(expected_qual, expected_depth, as_ad)
}

/// `AnnotationUtils.decodeAnyASListWithRawDelim`: brackets removed, split at every `|`, and
/// nothing at all for an empty string.
fn decode_raw_list(text: &str) -> Vec<String> {
    let stripped: String = text.chars().filter(|c| *c != '[' && *c != ']').collect();
    if stripped.is_empty() {
        return Vec::new();
    }
    stripped.split('|').map(str::to_string).collect()
}

fn variant_difference(thing: &str, actual: &str, expected: &str) -> Failure {
    Failure::User(format!(
        "Variant contexts have different {thing}: actual has {actual} expected has {expected}"
    ))
}

fn genotype_difference(thing: &str, actual: &str, expected: &str) -> Failure {
    Failure::User(format!(
        "Genotypes have different {thing}: actual has {actual} expected has {expected}"
    ))
}

/// `wrapWithPosition`, which puts the position in FRONT of the message.
fn wrapped(contig: &str, start: i64, failure: Failure) -> Failure {
    match failure {
        Failure::User(message) => Failure::User(format!("At position {contig}:{start} {message}")),
        other => other,
    }
}

/// The comparison, with the three fields the reference keeps between groups.
struct Comparator<'a> {
    options: &'a Options,
    single_sample: bool,
    allele_number_is_different: bool,
    inbreeding_coeff_is_different: bool,
    fail_on_completion: bool,
    warnings: Vec<String>,
}

impl Comparator<'_> {
    fn warn(&mut self, message: String) {
        self.warnings.push(message);
    }

    /// `throwOrWarn`.
    fn throw_or_warn(&mut self, failure: Failure) -> Result<(), Failure> {
        let Failure::User(message) = failure else {
            return Err(failure);
        };
        if self.options.mute_acceptable_diffs
            && (self.allele_number_is_different || self.inbreeding_coeff_is_different)
        {
            return Ok(());
        }
        if self.options.warn_on_errors || self.options.finish_before_failing {
            self.warn(format!("***** {message} *****"));
            if self.options.finish_before_failing {
                self.fail_on_completion = true;
            }
            Ok(())
        } else {
            Err(Failure::User(message))
        }
    }

    /// `hasGoodEvidence`.
    fn has_good_evidence(&self, vc: &Record) -> Result<bool, Failure> {
        let above =
            |value: f64, bound: f64| value.partial_cmp(&bound) == Some(std::cmp::Ordering::Greater);
        if !above(vc.qual(), self.options.good_qual_threshold) {
            return Ok(false);
        }
        let depth = f64::from(vc.attribute_as_int(DEPTH_KEY)?);
        let number = f64::from(vc.attribute_as_int(ALLELE_NUMBER_KEY)?);
        if !above(depth / number, 5.0) {
            return Ok(false);
        }
        let sum: i32 = vc
            .genotypes()
            .iter()
            .filter(|genotype| !genotype.is_hom_ref())
            .map(gq)
            .sum();
        Ok(f64::from(sum) > self.options.good_qual_threshold)
    }

    /// `isHighQuality`.
    fn is_high_quality(&self, vc: &Record) -> Result<bool, Failure> {
        let genotypes = vc.genotypes();
        if genotypes.len() == 1 {
            let genotype = &genotypes[0];
            let Some(pl) = &genotype.pl else {
                return Ok(false);
            };
            if genotype.is_hom_ref()
                || pl.first() == Some(&0)
                || genotype.alleles.iter().any(is_non_ref)
            {
                return Ok(false);
            }
            let max = pl.iter().copied().max().unwrap_or(0);
            return Ok(vc.qual() > 0.01 && max > 0);
        }
        for genotype in &genotypes {
            if genotype.is_hom_ref() {
                continue;
            }
            if passes_gnomad_adj(genotype)? {
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// `apply(List<VariantContext>, ...)`.
    fn apply(&mut self, group: &[Record]) -> Result<(), Failure> {
        if self.single_sample && !group[0].alleles().iter().any(is_non_ref) {
            return Err(Failure::User(
                "Bad input: Single-sample mode expects two GVCFs with <NON_REF> data for \
                 comparison"
                    .to_string(),
            ));
        }

        let gq_zero = |vc: &Record| vc.genotypes().iter().any(|genotype| genotype.gq == Some(0));
        if group.len() == 1 {
            let vc = &group[0];
            let unmatched = || {
                Failure::User(format!(
                    "Unmatched variant in {} at position {}:{}",
                    vc.source,
                    vc.contig(),
                    vc.start()
                ))
            };
            // With the diffs muted, only a site with good evidence is looked at, and one without
            // falls through to the comparison below.
            let looked_at = !self.options.mute_acceptable_diffs || self.has_good_evidence(vc)?;
            if looked_at {
                if !self.options.ignore_gq0 && gq_zero(vc) {
                    self.throw_or_warn(unmatched())?;
                } else {
                    return Ok(());
                }
            }
        }

        for vc in group {
            // Only an expected record is compared, against whatever else starts where it does.
            if vc.source == "actual" {
                continue;
            }
            let matches: Vec<&Record> = group
                .iter()
                .filter(|other| other.start() == vc.start())
                .collect();
            if matches.len() == 1 {
                let spanning = vc
                    .genotypes()
                    .iter()
                    .any(|genotype| has_span_del(&genotype.alleles));
                if self.is_high_quality(vc)? && self.has_good_evidence(vc)? && !spanning {
                    self.throw_or_warn(Failure::User(format!(
                        "Apparent unmatched high quality variant in {} at {}:{}",
                        vc.source,
                        vc.contig(),
                        vc.start()
                    )))?;
                } else {
                    return Ok(());
                }
            } else {
                let skippable = !self.is_high_quality(vc)?
                    || (self.single_sample && !vc.genotype(0)?.alleles.iter().any(is_concrete_alt));
                if skippable && self.options.ignore_hom_ref_attributes {
                    return Ok(());
                }
                let matched = if matches[0].source != vc.source {
                    matches[0]
                } else {
                    matches[1]
                };
                let deletions: Vec<&Record> = group
                    .iter()
                    .filter(|other| other.start() < vc.start() && other.overlaps(vc))
                    .collect();
                let expected_trimmed = self.trim_alleles(vc, &deletions)?;
                let actual_trimmed = self.trim_alleles(matched, &deletions)?;
                let no_calls_nearby = group.iter().any(|other| {
                    other
                        .genotypes()
                        .iter()
                        .any(|genotype| genotype.is_no_call() || genotype.gq == Some(0))
                });
                if self.single_sample && !(self.options.ignore_gq0 && no_calls_nearby) {
                    match self.validate_single_sample_deletions(
                        vc,
                        matched,
                        &expected_trimmed,
                        &actual_trimmed,
                        &deletions,
                        no_calls_nearby,
                    ) {
                        Err(failure @ Failure::User(_)) => self.throw_or_warn(failure)?,
                        Err(other) => return Err(other),
                        Ok(()) => {}
                    }
                }

                if self.options.positions_only {
                    return Ok(());
                }

                match self.check_matching(
                    &actual_trimmed,
                    &expected_trimmed,
                    &deletions,
                    no_calls_nearby,
                ) {
                    Err(failure @ Failure::User(_)) => {
                        let low_quality = expected_trimmed
                            .genotypes()
                            .iter()
                            .any(|genotype| gq(genotype) < 20);
                        if !self.options.mute_acceptable_diffs
                            || !(self.allele_number_is_different
                                || self.inbreeding_coeff_is_different
                                || low_quality)
                        {
                            self.throw_or_warn(failure)?;
                        }
                    }
                    Err(other) => return Err(other),
                    Ok(()) => {}
                }

                if self.allele_number_is_different && !self.options.mute_acceptable_diffs {
                    self.warn(format!(
                        "Observed allele number differed at position {}:{}",
                        vc.contig(),
                        vc.start()
                    ));
                }
                let low_quality = expected_trimmed
                    .genotypes()
                    .iter()
                    .any(|genotype| gq(genotype) < 20);
                if self.inbreeding_coeff_is_different
                    && low_quality
                    && !self.options.mute_acceptable_diffs
                {
                    self.warn(format!(
                        "Low quality genotype may have caused inbreeding coeff differences at \
                         position {}:{}",
                        vc.contig(),
                        vc.start()
                    ));
                }
            }
        }
        Ok(())
    }

    /// `validateSingleSampleDeletions`.
    fn validate_single_sample_deletions(
        &mut self,
        vc: &Record,
        matched: &Record,
        expected_trimmed: &Record,
        actual_trimmed: &Record,
        deletions: &[&Record],
        nearby_gq0: bool,
    ) -> Result<(), Failure> {
        let expected_genotype = expected_trimmed.genotype(0)?;
        let expected_alleles = &expected_genotype.alleles;
        let actual_genotype = actual_trimmed.genotype(0)?;
        let actual_alleles = &actual_genotype.alleles;
        let at = |failure: Failure| wrapped(vc.contig(), vc.start(), failure);

        if self.options.ignore_genotype_phasing
            && !has_new_alleles(expected_trimmed, actual_trimmed)?
        {
            return Ok(());
        }

        let get = |alleles: &Vec<Allele>, index: usize| -> Result<Allele, Failure> {
            alleles
                .get(index)
                .cloned()
                .ok_or_else(|| index_out_of_bounds(index, alleles.len()))
        };
        let first_differs = get(expected_alleles, 0)? != get(actual_alleles, 0)?;
        let second_differs = expected_alleles.len() > 1
            && actual_alleles.len() > 1
            && expected_alleles[1] != actual_alleles[1];
        if !(first_differs || second_differs) {
            return Ok(());
        }
        if expected_genotype.phased
            && get(expected_alleles, 1)? == get(actual_alleles, 0)?
            && get(expected_alleles, 0)? == get(actual_alleles, 1)?
        {
            let phase_set = |genotype: &Genotype| {
                genotype
                    .get(HAPLOTYPE_CALLER_PHASING_ID_KEY)
                    .map(value_string)
                    .unwrap_or_default()
            };
            return Err(at(Failure::User(format!(
                "phasing is swapped. Actual in phaseset {} has {} expected in phaseset {} has {}",
                phase_set(&actual_genotype),
                genotype_string(&actual_genotype),
                phase_set(&expected_genotype),
                genotype_string(&expected_genotype)
            ))));
        }
        if deletions.is_empty() {
            return Err(at(genotype_difference(
                "called genotype alleles",
                &genotype_string(&actual_genotype),
                &genotype_string(&expected_genotype),
            )));
        }

        // Fine if a star was corrected by dropping a hom-ref deletion upstream.
        let mut hom_ref_deletion = false;
        for deletion in deletions {
            if deletion.genotype(0)?.is_hom_ref() {
                hom_ref_deletion = true;
                break;
            }
        }
        if hom_ref_deletion {
            // Note the arguments: the EXPECTED record as actual, and every allele as alternates.
            let dp = depth_for_qual(&vc.genotypes());
            let as_ad = allele_depth_sum(&vc.genotypes())?;
            match self.check_attributes(
                &vc.attributes,
                &matched.attributes,
                vc.alleles(),
                matched.alleles(),
                vc.qual(),
                dp,
                as_ad,
            ) {
                Err(Failure::User(_)) => {
                    return Err(at(Failure::User(format!(
                        "INFO attributes do not match at {}:{}",
                        vc.contig(),
                        vc.start()
                    ))));
                }
                Err(other) => return Err(other),
                Ok(()) => {}
            }
            return self
                .check_genotypes(&expected_genotype, &actual_genotype, nearby_gq0)
                .map_err(at);
        }
        if deletions.iter().any(|deletion| deletion.source == "actual") {
            if !self.options.allow_missing_stars {
                return Err(at(Failure::User(format!(
                    "genotype alleles do not match at spanning deletion site. {} has {} and {} \
                     has {}",
                    vc.source,
                    genotype_to_string(&expected_genotype),
                    matched.source,
                    genotype_to_string(&actual_genotype)
                ))));
            }
            return Ok(());
        }
        // The overlapper's own trim may stop it overlapping, which is fine; if it does not, the
        // reference lets the difference through as well.
        let overlapper = deletions[0];
        let _ = self.trim_alleles(overlapper, deletions)?;
        Ok(())
    }

    /// `checkVariantContextsAreMatching`.
    fn check_matching(
        &mut self,
        actual: &Record,
        expected: &Record,
        deletions: &[&Record],
        nearby_gq0: bool,
    ) -> Result<(), Failure> {
        let at = |failure: Failure| wrapped(expected.contig(), expected.start(), failure);
        if actual.contig() != expected.contig() {
            return Err(at(Failure::User("contigs differ for VCs".to_string())));
        }
        if actual.start() != expected.start() {
            return Err(at(Failure::User(
                "start positions differ for VCs".to_string(),
            )));
        }

        let skip = self.options.ignore_gq0 && nearby_gq0;
        if !skip && has_new_alleles(expected, actual)? {
            self.check_alleles(actual, expected, deletions)
                .map_err(at)?;
        }

        if !self.options.ignore_dbsnp_ids
            && actual.vc.id != expected.vc.id
            && actual.alleles().len() == expected.alleles().len()
        {
            return Err(at(Failure::User("dbsnp IDs differ for VCs".to_string())));
        }

        if self.options.ignore_annotations {
            return Ok(());
        }

        if !self.options.ignore_quals && !skip {
            let difference = (actual.qual() - expected.qual()).abs();
            if difference > self.options.qual_change_allowed {
                return Err(at(Failure::User(format!(
                    "qual scores differ by {}, which is more than {}",
                    java_double_to_string(difference),
                    java_double_to_string(self.options.qual_change_allowed)
                ))));
            }
        }

        if !skip {
            let dp = depth_for_qual(&expected.genotypes());
            let as_ad = allele_depth_sum(&expected.genotypes())?;
            self.check_attributes(
                &actual.attributes,
                &expected.attributes,
                actual.alternates(),
                expected.alternates(),
                expected.qual(),
                dp,
                as_ad,
            )
            .map_err(at)?;
        }

        if !self.allele_number_is_different {
            if !self.options.ignore_filters {
                if actual.vc.filters_were_applied() != expected.vc.filters_were_applied() {
                    return Err(at(Failure::User(
                        " filters were not applied to both variants".to_string(),
                    )));
                }
                let set = |record: &Record| {
                    let mut filters = record.vc.filters.clone().unwrap_or_default();
                    filters.dedup();
                    filters
                };
                let (actual_filters, expected_filters) = (set(actual), set(expected));
                let same = actual_filters.len() == expected_filters.len()
                    && actual_filters
                        .iter()
                        .all(|filter| expected_filters.contains(filter));
                if !same {
                    let expected_text = if expected_filters.is_empty() {
                        "PASS".to_string()
                    } else {
                        format!("[{}]", expected_filters.join(", "))
                    };
                    return Err(at(Failure::User(format!(
                        "variants have different filters: expected has {expected_text} and actual \
                         has [{}]",
                        actual_filters.join(", ")
                    ))));
                }
            }
            // The loop runs once per genotype and compares the FIRST each time.
            for _ in 0..actual.genotypes().len() {
                self.check_genotypes(&actual.genotype(0)?, &expected.genotype(0)?, nearby_gq0)
                    .map_err(at)?;
            }
        }
        Ok(())
    }

    /// `checkAlleles`.
    fn check_alleles(
        &mut self,
        actual: &Record,
        expected: &Record,
        deletions: &[&Record],
    ) -> Result<(), Failure> {
        if self.options.ignore_genotype_phasing && !has_new_alleles(actual, expected)? {
            return Ok(());
        }
        let mismatch = |prefix: &str| {
            Failure::User(format!(
                "{prefix}Alleles are mismatched at {}:{}: actual has {} and expected has {}",
                actual.contig(),
                actual.start(),
                alleles_string(actual.alternates()),
                alleles_string(expected.alternates())
            ))
        };
        let new_star = has_span_del(actual.alleles()) && !has_span_del(expected.alleles());
        let missing_star = has_span_del(expected.alleles()) && !has_span_del(actual.alleles());
        if !self.options.allow_extra_alleles && has_new_alleles(actual, expected)? {
            Err(mismatch(""))
        } else if !self.options.allow_new_stars && new_star {
            if deletions.is_empty() || !deletions.iter().any(|deletion| deletion.source == "actual")
            {
                Err(mismatch("Actual has new unmatched * allele. "))
            } else {
                Ok(())
            }
        } else if !self.options.allow_missing_stars && missing_star {
            let mut remainder: Vec<&Allele> = Vec::new();
            for allele in expected.alleles() {
                if !actual.alleles().contains(allele) && !remainder.contains(&allele) {
                    remainder.push(allele);
                }
            }
            if remainder.len() > 1 || !remainder.iter().any(|allele| is_span_del(allele)) {
                Err(mismatch("Actual missing * allele. "))
            } else {
                Ok(())
            }
        } else {
            Err(mismatch(""))
        }
    }

    /// `isAttributeValueEqual`, which throws rather than answering false.
    fn attribute_value_equal(
        &self,
        key: &str,
        actual: &str,
        expected: &str,
    ) -> Result<bool, Failure> {
        if self.options.allow_nan_mismatch && actual == "." && expected == "NaN" {
            return Ok(true);
        }
        if actual != expected {
            return Err(Failure::User(format!(
                "Variant contexts have different attribute values for {key}: actual has {actual} \
                 and expected has {expected}"
            )));
        }
        Ok(true)
    }

    /// `isAttributeEqualDoubleSmart`.
    fn double_equal(
        key: &str,
        actual: f64,
        expected: f64,
        tolerance: f64,
    ) -> Result<bool, Failure> {
        let difference = (actual - expected).abs();
        if difference > tolerance {
            return Err(Failure::User(format!(
                "Attribute {key} has difference {}, which is larger difference than allowed \
                 delta {}",
                format_decimals(difference, 3),
                java_double_to_string(tolerance)
            )));
        }
        Ok(true)
    }

    /// `isACEqualEnough`.
    fn ac_equal_enough(
        &mut self,
        actual_list: &[String],
        expected_list: &[String],
        actual_alts: &[Allele],
        expected_alts: &[Allele],
    ) -> Result<bool, Failure> {
        if actual_list.len() != expected_list.len() {
            return Ok(false);
        }
        for (index, allele) in actual_alts.iter().enumerate() {
            if is_span_del(allele) {
                continue;
            }
            let expected = expected_alts
                .get(index)
                .ok_or_else(|| index_out_of_bounds(index, expected_alts.len()))?;
            if allele == expected {
                let item = |list: &[String]| {
                    list.get(index)
                        .cloned()
                        .ok_or_else(|| index_out_of_bounds(index, list.len()))
                };
                parse_int(&item(expected_list)?)?;
                parse_int(&item(actual_list)?)?;
            } else {
                self.throw_or_warn(Failure::User(
                    "Alleles are not ordered the same".to_string(),
                ))?;
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// `checkAttributes`.
    #[allow(clippy::too_many_arguments)]
    fn check_attributes(
        &mut self,
        actual: &[(String, Value)],
        expected: &[(String, Value)],
        actual_alts: &[Allele],
        expected_alts: &[Allele],
        expected_qual: f64,
        expected_dp: i32,
        expected_as_ad: i32,
    ) -> Result<(), Failure> {
        let get = |map: &[(String, Value)], key: &str| -> Option<Value> {
            map.iter()
                .find(|(name, _)| name == key)
                .map(|(_, value)| value.clone())
        };
        let text = |value: Option<Value>| -> Result<String, Failure> {
            value.as_ref().map(value_string).ok_or_else(null_pointer)
        };

        // The allele number first, because the rest cannot be expected to match without it. The
        // comparison throws on a difference, so the assignment never sees one.
        if let Some(actual_number) = get(actual, ALLELE_NUMBER_KEY) {
            let expected_number = text(get(expected, ALLELE_NUMBER_KEY))?;
            let equal = self.attribute_value_equal(
                ALLELE_NUMBER_KEY,
                &value_string(&actual_number),
                &expected_number,
            )?;
            self.allele_number_is_different = !equal;
        }

        if self.options.ignore_some_multi_allelics {
            if let Some(actual_ac) = get(actual, ALLELE_COUNT_KEY) {
                let size = |value: Option<&Value>| match value {
                    Some(Value::List(items)) => items.len(),
                    _ => 1,
                };
                let expected_ac = get(expected, ALLELE_COUNT_KEY);
                let actual_size = size(Some(&actual_ac));
                if actual_size != size(expected_ac.as_ref()) && actual_size == actual_alts.len() {
                    return Ok(());
                }
            }
        }

        if let Some(actual_coefficient) = get(actual, INBREEDING_COEFFICIENT_KEY) {
            let result = match get(expected, INBREEDING_COEFFICIENT_KEY) {
                Some(expected_coefficient) => {
                    let a = parse_double(&value_string(&actual_coefficient))?;
                    let e = parse_double(&value_string(&expected_coefficient))?;
                    Self::double_equal(INBREEDING_COEFFICIENT_KEY, a, e, 0.001)
                }
                None => Err(variant_difference(
                    INBREEDING_COEFFICIENT_KEY,
                    &value_string(&actual_coefficient),
                    "missing",
                )),
            };
            match result {
                Ok(equal) => self.inbreeding_coeff_is_different = !equal,
                Err(failure @ Failure::User(_)) => {
                    self.inbreeding_coeff_is_different = true;
                    return Err(failure);
                }
                Err(other) => return Err(other),
            }
        }

        if let Some(value) = get(actual, AS_INBREEDING_COEFFICIENT_KEY) {
            if !matches!(value, Value::List(_)) {
                let a = parse_double(&value_string(&value))?;
                let e = parse_double(&text(get(expected, AS_INBREEDING_COEFFICIENT_KEY))?)?;
                match Self::double_equal(AS_INBREEDING_COEFFICIENT_KEY, a, e, 0.001) {
                    Ok(equal) => self.inbreeding_coeff_is_different = !equal,
                    Err(failure @ Failure::User(_)) => {
                        self.inbreeding_coeff_is_different = true;
                        return Err(failure);
                    }
                    Err(other) => return Err(other),
                }
            }
        }

        let keys = |map: &[(String, Value)]| {
            format!(
                "[{}]",
                map.iter()
                    .map(|(key, _)| key.as_str())
                    .collect::<Vec<_>>()
                    .join(", ")
            )
        };
        for (key, expected_value) in expected {
            let Some(actual_value) = get(actual, key) else {
                continue;
            };
            if key == INBREEDING_COEFFICIENT_KEY
                || key == ALLELE_NUMBER_KEY
                || (key == AS_INBREEDING_COEFFICIENT_KEY
                    && !matches!(expected_value, Value::List(_)))
                || self.options.ignore_attributes.contains(key)
            {
                continue;
            }
            let both_lists =
                matches!(expected_value, Value::List(_)) && matches!(actual_value, Value::List(_));
            if both_lists || key == AS_SB_TABLE_KEY {
                let (expected_list, actual_list): (Vec<String>, Vec<String>) =
                    if (key.contains("AS_") && key.contains("RAW")) || key == AS_SB_TABLE_KEY {
                        (
                            decode_raw_list(&value_string(expected_value)),
                            decode_raw_list(&value_string(&actual_value)),
                        )
                    } else {
                        let items = |value: &Value| match value {
                            Value::List(items) => items.iter().map(value_string).collect(),
                            other => vec![value_string(other)],
                        };
                        (items(expected_value), items(&actual_value))
                    };
                if actual_list.len() != expected_list.len() && !key.contains("AS_") {
                    return Err(variant_difference(
                        "attributes",
                        &keys(actual),
                        &keys(expected),
                    ));
                }
                if (key == ALLELE_COUNT_KEY || key == MLE_ALLELE_COUNT_KEY)
                    && !self.ac_equal_enough(
                        &actual_list,
                        &expected_list,
                        actual_alts,
                        expected_alts,
                    )?
                {
                    return Err(variant_difference(
                        key,
                        &format!("[{}]", actual_list.join(", ")),
                        &format!("[{}]", expected_list.join(", ")),
                    ));
                }
                let end = if key == RAW_GENOTYPE_COUNT_KEY {
                    3
                } else if key == RAW_MAPPING_QUALITY_WITH_DEPTH_KEY {
                    2
                } else if self.options.ignore_non_ref_data && actual_alts.iter().any(is_non_ref) {
                    expected_alts.len() as i64 - 1
                } else {
                    expected_alts.len() as i64
                };
                let mut index: i64 = 0;
                while index < end {
                    let i = index as usize;
                    index += 1;
                    if i >= actual_list.len() || i >= expected_list.len() {
                        return Err(variant_difference(
                            key,
                            &value_string(&actual_value),
                            &value_string(expected_value),
                        ));
                    }
                    let (a, e) = (&actual_list[i], &expected_list[i]);
                    if key == AS_INBREEDING_COEFFICIENT_KEY {
                        let (a, e) = (parse_double(a)?, parse_double(e)?);
                        Self::double_equal(
                            key,
                            a,
                            e,
                            self.options.inbreeding_coeff_change_allowed,
                        )?;
                    }
                    if key.starts_with("AS_") && key.contains("RankSum") && !key.contains("RAW") {
                        if !a.is_empty() && a == "NaN" {
                            let (a, e) = (parse_double(a)?, parse_double(e)?);
                            Self::double_equal(key, a, e, self.options.ranksum_change_allowed)?;
                            continue;
                        } else if !self.options.mute_acceptable_diffs
                            && !self.options.allow_nan_mismatch
                        {
                            self.warn(
                                "GATK version-specific NaN versus empty AS_RAW annotation \
                                 discrepancy"
                                    .to_string(),
                            );
                        }
                    } else if key.starts_with("AS_")
                        && key.contains("RAW")
                        && key.contains("RankSum")
                    {
                        if a.is_empty() && e.is_empty() {
                            continue;
                        }
                        if ((a.is_empty() && e == "NaN") || (a == "NaN" && e.is_empty()))
                            && !self.options.mute_acceptable_diffs
                            && !self.options.allow_nan_mismatch
                        {
                            self.warn(
                                "GATK version-specific NaN versus empty AS_RAW annotation \
                                 discrepancy"
                                    .to_string(),
                            );
                        }
                        if a != e {
                            return Err(variant_difference(key, a, e));
                        }
                    }
                    let checked = if key == AS_QUAL_BY_DEPTH_KEY {
                        let (a, e) = (parse_double(a)?, parse_double(e)?);
                        let depth = parse_int(&text(get(expected, DEPTH_KEY))?)?;
                        Ok(qual_by_depth_difference_is_acceptable(
                            a,
                            e,
                            expected_qual,
                            depth,
                            expected_as_ad,
                        ))
                    } else {
                        self.attribute_value_equal(key, a, e)
                    };
                    match checked {
                        Ok(_) => {}
                        Err(Failure::User(_)) => {
                            if !self.options.ignore_star_attributes {
                                let alternate = actual_alts
                                    .get(i)
                                    .ok_or_else(|| index_out_of_bounds(i, actual_alts.len()))?;
                                if !is_span_del(alternate) {
                                    return Err(variant_difference(key, a, e));
                                }
                            }
                        }
                        Err(other) => return Err(other),
                    }
                }
            } else {
                let compared = if key.contains("RankSum") && !key.contains("AS_") {
                    let a = parse_double(&value_string(&actual_value))?;
                    let e = parse_double(&value_string(expected_value))?;
                    Self::double_equal(key, a, e, self.options.ranksum_change_allowed)
                } else {
                    self.attribute_value_equal(
                        key,
                        &value_string(&actual_value),
                        &value_string(expected_value),
                    )
                };
                let failure = match compared {
                    Ok(_) => continue,
                    Err(failure @ Failure::User(_)) => failure,
                    Err(other) => return Err(other),
                };
                if VARY_WITH_NO_CALLS.contains(&key.as_str()) && self.allele_number_is_different {
                    continue;
                }
                if key == QUAL_BY_DEPTH_KEY || key == AS_QUAL_BY_DEPTH_KEY {
                    let a = parse_double(&value_string(&actual_value))?;
                    let e = parse_double(&value_string(expected_value))?;
                    if qual_by_depth_difference_is_acceptable(
                        a,
                        e,
                        expected_qual,
                        expected_dp,
                        expected_as_ad,
                    ) {
                        if !self.options.mute_acceptable_diffs {
                            self.warn(format!("{key} difference is within expected tolerances"));
                        }
                        continue;
                    }
                    let actual_depth = text(get(actual, DEPTH_KEY))?;
                    let expected_depth = text(get(expected, DEPTH_KEY))?;
                    match self.attribute_value_equal(DEPTH_KEY, &actual_depth, &expected_depth) {
                        Ok(true) => return Err(failure),
                        Ok(false) => {}
                        Err(Failure::User(_)) => {
                            let qd = parse_double(&text(get(actual, QUAL_BY_DEPTH_KEY))?)?;
                            let as_qd = parse_double(&text(get(actual, AS_QUAL_BY_DEPTH_KEY))?)?;
                            if !qual_by_depth_difference_is_acceptable(
                                qd,
                                as_qd,
                                expected_qual,
                                expected_dp,
                                expected_as_ad,
                            ) {
                                self.warn(format!(
                                    "{key} difference (actual = {} versus expected:{} is larger \
                                     than expected, but so is DP (actual={actual_depth} versus \
                                     expected:{expected_depth})",
                                    value_string(&actual_value),
                                    value_string(expected_value)
                                ));
                            }
                            continue;
                        }
                        Err(other) => return Err(other),
                    }
                }
                return Err(failure);
            }
        }
        Ok(())
    }

    /// `checkGenotypes`.
    fn check_genotypes(
        &self,
        actual: &Genotype,
        expected: &Genotype,
        nearby_gq0: bool,
    ) -> Result<(), Failure> {
        let options = self.options;
        if actual.sample_name != expected.sample_name {
            return Err(Failure::User("Sample names do not match".to_string()));
        }
        if options.ignore_gq0 && nearby_gq0 {
            return Ok(());
        }
        let mut sorted_actual = actual.alleles.clone();
        let mut sorted_expected = expected.alleles.clone();
        sorted_actual.sort_by(allele_order);
        sorted_expected.sort_by(allele_order);
        if sorted_actual != sorted_expected {
            return Err(genotype_difference(
                "alleles",
                &alleles_string(&actual.alleles),
                &alleles_string(&expected.alleles),
            ));
        }
        if !options.ignore_genotype_phasing {
            let (a, e) = (genotype_string(actual), genotype_string(expected));
            if a != e {
                return Err(genotype_difference("genotype string", &a, &e));
            }
            if actual.phased != expected.phased {
                return Err(genotype_difference(
                    "phasing status",
                    &actual.phased.to_string(),
                    &expected.phased.to_string(),
                ));
            }
        }
        if options.ignore_hom_ref_attributes && actual.is_hom_ref() {
            return Ok(());
        }
        let presence = |thing: &str, a: bool, e: bool| {
            if a != e {
                Err(genotype_difference(thing, &a.to_string(), &e.to_string()))
            } else {
                Ok(())
            }
        };
        presence("DP presence", actual.dp.is_some(), expected.dp.is_some())?;
        presence("AD presence", actual.ad.is_some(), expected.ad.is_some())?;
        presence("GQ presence", actual.gq.is_some(), expected.gq.is_some())?;
        if !actual.is_hom_ref()
            && actual.gq.is_some()
            && (gq(actual) - gq(expected)).abs() > options.likelihood_change_allowed
            && (gq(actual) < 99 || gq(expected) < 99)
        {
            return Err(genotype_difference(
                "GQ value",
                &gq(actual).to_string(),
                &gq(expected).to_string(),
            ));
        }
        if options.ignore_genotype_annotations {
            return Ok(());
        }
        if let (Some(actual_dp), Some(expected_dp)) = (actual.dp, expected.dp) {
            if actual_dp != expected_dp {
                if options.dp_change_allowed == 0 && expected_dp - actual_dp != 0 {
                    return Err(genotype_difference(
                        "DP value",
                        &actual_dp.to_string(),
                        &expected_dp.to_string(),
                    ));
                }
                if actual_dp - expected_dp > options.dp_change_allowed {
                    return Err(Failure::User(format!(
                        "DP difference exceeds allowable tolerance of {}, actual has {actual_dp} \
                         expected has {expected_dp}",
                        options.dp_change_allowed
                    )));
                }
                if let Some(actual_ad) = &actual.ad {
                    let expected_ad = expected.ad.clone().unwrap_or_default();
                    if *actual_ad != expected_ad {
                        let last = |ad: &[i32]| ad.last().copied();
                        if options.dp_change_allowed == 0 && last(actual_ad) == last(&expected_ad) {
                            return Err(genotype_difference(
                                "AD values",
                                &ints_string(actual_ad),
                                &ints_string(&expected_ad),
                            ));
                        }
                        let compared = actual_ad.len().saturating_sub(1);
                        for (index, actual_value) in actual_ad.iter().take(compared).enumerate() {
                            let expected_value = expected_ad
                                .get(index)
                                .copied()
                                .ok_or_else(|| index_out_of_bounds(index, expected_ad.len()))?;
                            if actual_value - expected_value > options.dp_change_allowed {
                                // `int[].toString()` is an identity hash, which no port can
                                // reproduce; the message is the reference's up to that point.
                                return Err(Failure::User(format!(
                                    "AD difference exceeds allowable tolerance of {}, actual has \
                                     [I@ expected has [I@",
                                    options.dp_change_allowed
                                )));
                            }
                        }
                    }
                }
            }
            presence("PL presence", actual.pl.is_some(), expected.pl.is_some())?;
            if let (Some(actual_pl), Some(expected_pl)) = (&actual.pl, &expected.pl) {
                if actual_pl.len() != expected_pl.len() {
                    return Err(Failure::User(format!(
                        "PL lengths for genotype vary: actual is size {} and expected is {}",
                        actual_pl.len(),
                        expected_pl.len()
                    )));
                }
                if actual_pl
                    .iter()
                    .zip(expected_pl)
                    .any(|(a, e)| (a - e).abs() > options.likelihood_change_allowed)
                {
                    return Err(genotype_difference(
                        "PL value",
                        &ints_string(actual_pl),
                        &ints_string(expected_pl),
                    ));
                }
            }
        }
        Ok(())
    }

    /// `trimAlleles`: the record cut down to the alleles its genotypes call, the genotypes
    /// subset to match, and a record that is not hom-ref trimmed from the right.
    fn trim_alleles(&self, variant: &Record, deletions: &[&Record]) -> Result<Record, Failure> {
        let genotypes = variant.genotypes();
        let reference = variant.reference().clone();
        let mut relevant: Vec<Allele> = Vec::new();
        let add = |allele: &Allele, relevant: &mut Vec<Allele>| {
            if !relevant.contains(allele) {
                relevant.push(allele.clone());
            }
        };
        for genotype in &genotypes {
            add(&reference, &mut relevant);
            if !genotype.alleles.iter().any(Allele::is_no_call) {
                for allele in &genotype.alleles {
                    if !deletions.is_empty() || is_concrete_alt(allele) {
                        add(allele, &mut relevant);
                    }
                }
            }
        }
        if genotypes.len() == 1 && !self.options.ignore_non_ref_data {
            let non_ref = Allele::from_str(NON_REF, false).map_err(|_| null_pointer())?;
            // A list, not a set: a `<NON_REF>` a genotype already called is added twice, and the
            // record built over it refuses the duplicate.
            if relevant.contains(&non_ref) {
                return Err(Failure::Runtime {
                    class: "java.lang.IllegalArgumentException",
                    message: format!("Duplicate allele added to VariantContext: {NON_REF}"),
                });
            }
            relevant.push(non_ref);
        }

        let limitation = || {
            Failure::Limitation(format!(
                "Comparing {}:{} asks for an allele the record does not carry, which this port \
                 does not subset. This message is the port's own and not GATK's.",
                variant.contig(),
                variant.start()
            ))
        };
        let kept: Vec<usize> = relevant
            .iter()
            .map(|allele| variant.alleles().iter().position(|a| a == allele))
            .collect::<Option<Vec<usize>>>()
            .ok_or_else(limitation)?;

        let engine: Vec<EngineGenotype> = genotypes
            .iter()
            .map(|genotype| engine_genotype(genotype, variant.alleles()))
            .collect();
        let subset = subset_alleles(
            &engine,
            self.options.default_ploidy,
            variant.alleles().len(),
            &kept,
            AssignmentMethod::BestMatchToOriginal,
        )
        .map_err(|error| Failure::Runtime {
            class: "java.lang.IllegalArgumentException",
            message: error.message(),
        })?;
        let rebuilt: Vec<Genotype> = genotypes
            .iter()
            .zip(&subset)
            .map(|(original, engine)| vcf_genotype(original, engine, &relevant))
            .collect();

        let mut vc = variant.vc.clone();
        vc.alleles = relevant.clone();
        vc.genotypes = GenotypesContext::new(rebuilt);
        // A hom-ref record is built before the annotations are set, so whatever the subsetting
        // computed for it is thrown away.
        let hom_ref = genotypes.first().is_some_and(Genotype::is_hom_ref);
        if hom_ref {
            return Ok(Record {
                source: variant.source.clone(),
                vc,
                attributes: variant.attributes.clone(),
            });
        }
        let attributes = if relevant.len() != variant.alleles().len() {
            let first = vc.genotypes.iter().next().cloned();
            self.subset_annotations(variant, first.as_ref())?
        } else {
            variant.attributes.clone()
        };
        Ok(Record {
            source: variant.source.clone(),
            vc: reverse_trim(&vc)?,
            attributes,
        })
    }

    /// `ReblockGVCF.subsetAnnotationsIfNecessary` for a record whose allele count changed:
    /// `composeUpdatedAnnotations` with no QUAL approximation.
    ///
    /// The new map holds the raw MQ when the record had none, every key an engine annotation claims
    /// except the eight the reblocking removes, the variant depth and QUAL approximation when
    /// present, a genotype count built from the first genotype, and `--annotations-to-keep`. A key
    /// no annotation claims is DROPPED. Allele-specific annotations are subset per allele, which is
    /// not ported, and neither is the set `--enable-all-annotations` claims beyond the standard
    /// one: a record that reaches either is the port's limitation.
    fn subset_annotations(
        &self,
        variant: &Record,
        genotype: Option<&Genotype>,
    ) -> Result<Vec<(String, Value)>, Failure> {
        let limitation = |what: String| {
            Failure::Limitation(format!(
                "Comparing {}:{} subsets {what} through GATK's annotation engine, which this port \
                 does not carry yet. This message is the port's own and not GATK's.",
                variant.contig(),
                variant.start()
            ))
        };
        let original = &variant.attributes;
        let has = |key: &str| original.iter().any(|(name, _)| name == key);
        for (key, _) in original {
            if key.starts_with("AS_") && key != AS_VARIANT_DEPTH_KEY {
                return Err(limitation(format!("the allele-specific {key}")));
            }
            if key == RAW_RMS_MAPPING_QUALITY_DEPRECATED {
                return Err(limitation(format!("the deprecated {key}")));
            }
            let known = STANDARD_ANNOTATION_KEYS.contains(&key.as_str())
                || REMOVED_ANNOTATION_KEYS.contains(&key.as_str())
                || key == AS_VARIANT_DEPTH_KEY
                || key == RAW_QUAL_APPROX_KEY
                || self.options.annotations_to_keep.contains(key);
            if self.options.enable_all_annotations && !known {
                return Err(limitation(format!("{key} under --enable-all-annotations")));
            }
        }
        let claimed =
            self.options.enable_all_annotations || !self.options.disable_tool_default_annotations;

        let mut destination: Vec<(String, Value)> = Vec::new();
        let put =
            |destination: &mut Vec<(String, Value)>, key: &str, value: Value| match destination
                .iter_mut()
                .find(|(name, _)| name == key)
            {
                Some(slot) => slot.1 = value,
                None => destination.push((key.to_string(), value)),
            };
        // `updateMQAnnotations`: a record with no raw MQ gets one, from MQ (60 when absent) and DP.
        if !has(RAW_MAPPING_QUALITY_WITH_DEPTH_KEY) {
            let attribute_error =
                |error: htsjdk_vcf::attributes::AttributeError| Failure::Runtime {
                    class: error.class(),
                    message: error.message(),
                };
            let find = |key: &str| {
                original
                    .iter()
                    .find(|(name, _)| name == key)
                    .map(|(_, value)| value)
            };
            let mq =
                htsjdk_vcf::attributes::as_double(find("MQ"), 60.0).map_err(attribute_error)?;
            let dp = htsjdk_vcf::attributes::as_int(find(DEPTH_KEY), 0).map_err(attribute_error)?;
            let raw = (mq * mq * f64::from(dp) + 0.5).floor() as i64 as i32;
            put(
                &mut destination,
                RAW_MAPPING_QUALITY_WITH_DEPTH_KEY,
                Value::Str(format!("{raw},{dp}")),
            );
        }
        if claimed {
            for key in STANDARD_ANNOTATION_KEYS {
                if REMOVED_ANNOTATION_KEYS.contains(&key) {
                    continue;
                }
                if let Some((_, value)) = original.iter().find(|(name, _)| name == key) {
                    put(&mut destination, key, value.clone());
                }
            }
        }
        for key in [AS_VARIANT_DEPTH_KEY, RAW_QUAL_APPROX_KEY] {
            if let Some((_, value)) = original.iter().find(|(name, _)| name == key) {
                put(&mut destination, key, value.clone());
            }
        }
        let called_reference = genotype
            .ok_or_else(|| index_out_of_bounds(0, 0))?
            .alleles
            .iter()
            .any(Allele::is_reference);
        let counts: [i64; 3] = if called_reference {
            [0, 1, 0]
        } else {
            [0, 0, 1]
        };
        put(
            &mut destination,
            RAW_GENOTYPE_COUNT_KEY,
            Value::List(counts.iter().map(|count| Value::Int(*count)).collect()),
        );
        for key in &self.options.annotations_to_keep {
            if let Some((_, value)) = original.iter().find(|(name, _)| name == key) {
                put(&mut destination, key, value.clone());
            }
        }
        // The builder copies the new map into one of its own, sized to it.
        hash_order(&destination, Some(copied_capacity(destination.len())))
    }
}

/// `passesGnomadAdjCriteria`.
fn passes_gnomad_adj(genotype: &Genotype) -> Result<bool, Failure> {
    if genotype.gq.is_none_or(|gq| gq < 20) || genotype.dp.is_none_or(|dp| dp < 10) {
        return Ok(false);
    }
    if genotype.is_het() {
        let Some(ad) = &genotype.ad else {
            return Ok(false);
        };
        let at = |index: usize| {
            ad.get(index)
                .map(|value| f64::from(*value))
                .ok_or_else(|| index_out_of_bounds(index, ad.len()))
        };
        let balance = if genotype.is_het_non_ref() {
            at(2)? / (at(1)? + at(2)?)
        } else {
            at(1)? / (at(0)? + at(1)?)
        };
        if (0.2..=0.8).contains(&balance) {
            return Ok(true);
        }
    }
    Ok(false)
}

/// A genotype in the index form the subsetting works in.
fn engine_genotype(genotype: &Genotype, alleles: &[Allele]) -> EngineGenotype {
    EngineGenotype {
        alleles: genotype
            .alleles
            .iter()
            .map(|allele| {
                if allele.is_no_call() {
                    None
                } else {
                    alleles.iter().position(|candidate| candidate == allele)
                }
            })
            .collect(),
        pl: genotype.pl.clone(),
        gq: genotype.gq,
        ad: genotype.ad.clone(),
        dp: genotype.dp,
        attributes: genotype
            .extended
            .iter()
            .map(|(key, value)| (key.clone(), value_string(value)))
            .collect(),
    }
}

/// The subset genotype back in the file's form, over the new allele list. An attribute the
/// subsetting left as it was keeps its decoded value.
fn vcf_genotype(original: &Genotype, subset: &EngineGenotype, alleles: &[Allele]) -> Genotype {
    Genotype {
        sample_name: original.sample_name.clone(),
        alleles: subset
            .alleles
            .iter()
            .map(|index| {
                index
                    .and_then(|index| alleles.get(index).cloned())
                    .unwrap_or_else(Allele::no_call)
            })
            .collect(),
        phased: original.phased,
        gq: subset.gq,
        dp: subset.dp,
        ad: subset.ad.clone(),
        pl: subset.pl.clone(),
        filters: original.filters.clone(),
        extended: subset
            .attributes
            .iter()
            .map(|(key, text)| {
                original
                    .extended
                    .iter()
                    .find(|(name, value)| name == key && value_string(value) == *text)
                    .map(|(name, value)| (name.clone(), value.clone()))
                    .unwrap_or_else(|| (key.clone(), Value::Str(text.clone())))
            })
            .collect(),
    }
}

/// `GATKVariantContextUtils.reverseTrimAlleles`, shared with GenotypeGVCFs.
fn reverse_trim(vc: &VariantContext) -> Result<VariantContext, Failure> {
    crate::variant_trim::reverse_trim_alleles(vc).map_err(|message| Failure::Runtime {
        class: "java.lang.IllegalStateException",
        message,
    })
}

/// The merge of the two inputs: `(input, record)` in the order the priority queue hands them out.
///
/// Ordered by contig, in the merged dictionary's order, then by start. On a tie the input already
/// waiting goes first, and at the start the input named first does.
pub fn merged_order(inputs: &[Input], contigs: &[String]) -> Result<Vec<(usize, usize)>, Failure> {
    let key = |input: usize, record: usize| -> Result<(usize, i64), Failure> {
        let vc = &inputs[input].records[record];
        let contig = contigs
            .iter()
            .position(|name| *name == vc.contig)
            .ok_or_else(null_pointer)?;
        Ok((contig, vc.start))
    };
    // The queue: `(input, next record)`, root first.
    let mut queue: Vec<(usize, usize)> = Vec::new();
    let offer = |queue: &mut Vec<(usize, usize)>, entry: (usize, usize)| -> Result<(), Failure> {
        // A binary heap's sift-up: the new entry moves above a parent only when strictly smaller.
        queue.push(entry);
        let mut child = queue.len() - 1;
        while child > 0 {
            let parent = (child - 1) / 2;
            if key(queue[child].0, queue[child].1)? < key(queue[parent].0, queue[parent].1)? {
                queue.swap(child, parent);
                child = parent;
            } else {
                break;
            }
        }
        Ok(())
    };
    for (index, input) in inputs.iter().enumerate() {
        if !input.records.is_empty() {
            offer(&mut queue, (index, 0))?;
        }
    }
    let mut order = Vec::new();
    let mut last: Option<(usize, i64)> = None;
    while !queue.is_empty() {
        // Poll: the root leaves, the last entry takes its place and sifts down.
        let (input, record) = queue.swap_remove(0);
        let mut parent = 0;
        loop {
            let (left, right) = (2 * parent + 1, 2 * parent + 2);
            if left >= queue.len() {
                break;
            }
            let mut smallest = left;
            if right < queue.len()
                && key(queue[right].0, queue[right].1)? < key(queue[left].0, queue[left].1)?
            {
                smallest = right;
            }
            if key(queue[smallest].0, queue[smallest].1)? < key(queue[parent].0, queue[parent].1)? {
                queue.swap(parent, smallest);
                parent = smallest;
            } else {
                break;
            }
        }
        let this = key(input, record)?;
        if last.is_some_and(|previous| previous > this) {
            return Err(Failure::Runtime {
                class: "java.lang.IllegalStateException",
                message: "The elements of the input Iterators are not sorted according to the \
                          comparator htsjdk.variant.variantcontext.VariantContextComparator"
                    .to_string(),
            });
        }
        last = Some(this);
        order.push((input, record));
        if record + 1 < inputs[input].records.len() {
            offer(&mut queue, (input, record + 1))?;
        }
    }
    Ok(order)
}

/// `makeSpanningReferenceContext`, which runs before every group is compared: the window is
/// widened by the padding, and a padded window's end is trimmed to the reference's contig, which
/// the reference must therefore have. Nothing is read.
fn spanning_reference(
    group: &[Record],
    reference: Option<&[String]>,
    padding: i32,
) -> Result<(), Failure> {
    if padding < 0 {
        return Err(Failure::Runtime {
            class: "org.broadinstitute.hellbender.exceptions.GATKException",
            message: "Reference window starts after the current interval".to_string(),
        });
    }
    if padding == 0 {
        return Ok(());
    }
    let contig = group[0].contig();
    match reference {
        Some(contigs) if !contigs.iter().any(|name| name == contig) => Err(Failure::User(format!(
            "Given reference file does not have data at the requested contig({contig})!"
        ))),
        _ => Ok(()),
    }
}

/// The traversal: the merged records grouped by overlap, each group compared when the next record
/// does not overlap it, the last after the traversal, and the deferred failure at the end.
///
/// `start_intervals` is `-L` for `--ignore-variants-starting-outside-interval`; the records are
/// already the ones the traversal reaches. `reference` is the reference's contigs.
pub fn compare(
    inputs: &[Input],
    contigs: &[String],
    start_intervals: Option<&[SimpleInterval]>,
    reference: Option<&[String]>,
    options: &Options,
) -> Result<Finished, Stopped> {
    let mut samples: Vec<&String> = inputs.iter().flat_map(|input| &input.samples).collect();
    samples.sort();
    samples.dedup();
    let mut comparator = Comparator {
        options,
        single_sample: samples.len() == 1,
        allele_number_is_different: false,
        inbreeding_coeff_is_different: false,
        fail_on_completion: false,
        warnings: Vec::new(),
    };

    let order = merged_order(inputs, contigs).map_err(|failure| Stopped { failure, at: None })?;
    let mut current: Vec<Record> = Vec::new();
    let mut last_contig = String::new();
    let mut last_end: i64 = 0;
    for (input, index) in order {
        let vc = &inputs[input].records[index];
        let stop = |failure: Failure| Stopped {
            failure,
            at: Some((input, index)),
        };
        if options.ignore_variants_starting_outside_interval {
            if let Some(intervals) = start_intervals {
                let start = vc.start as i32;
                if !intervals
                    .iter()
                    .any(|interval| interval.overlaps(&vc.contig, start, start))
                {
                    continue;
                }
            }
        }
        if options.ignore_reference_blocks && vc.alleles.len() == 2 && is_non_ref(&vc.alleles[1]) {
            continue;
        }
        let record = Record::decoded(&inputs[input].name, vc).map_err(stop)?;
        if current.is_empty() {
            last_contig = vc.contig.clone();
        } else if current[0].contig() != vc.contig || vc.start > last_end {
            let group = std::mem::take(&mut current);
            spanning_reference(&group, reference, options.reference_padding).map_err(stop)?;
            comparator.apply(&group).map_err(stop)?;
        }
        current.push(record);
        if vc.stop > last_end || last_contig != vc.contig {
            last_end = vc.stop;
            last_contig = vc.contig.clone();
        }
    }
    if current.is_empty() {
        comparator.warn(
            "Error: The requested interval contained no data in source VCF files".to_string(),
        );
    } else {
        spanning_reference(&current, reference, options.reference_padding)
            .and_then(|()| comparator.apply(&current))
            .map_err(|failure| Stopped { failure, at: None })?;
    }
    if comparator.fail_on_completion {
        return Err(Stopped {
            failure: Failure::User(
                "Some comparisons failed.  See stderr log for details.".to_string(),
            ),
            at: None,
        });
    }
    Ok(Finished {
        warnings: comparator.warnings,
    })
}

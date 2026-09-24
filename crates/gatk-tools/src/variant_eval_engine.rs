//! `VariantEvalEngine`: every stratification crossed with every other, and every evaluation module
//! counting the records that fall in each combination.
//!
//! Ported from `org.broadinstitute.hellbender.tools.walkers.varianteval` (GATK 4.6.2.0). The walker
//! groups the records of every driving input (the evals, `--dbsnp` and the comps) by start; each
//! group is bound per eval track and per sample, subset to the samples under evaluation, and fed
//! to the evaluation context of every combination of states its stratifiers answer.
//!
//! # The report
//!
//! One GATKReport table per evaluation module, sorted by name as `TreeSet<VariantEvaluator>` sorts
//! them, and in each one row per combination of states, sorted by the row key the combination
//! builds (`CompFeatureInput:dbsnpEvalFeatureInput:eval...`). A molten module writes one row per
//! map entry, keyed by the combination and a five-digit counter.
//!
//! # Order does not reach the output
//!
//! The reference keeps its strat tree in `HashMap`s and its evaluators in a `HashSet` of classes,
//! so the order it visits contexts in and the order a context updates its modules in are the
//! JVM's. None of it is visible: every module only counts, the table is sorted by row key, and the
//! modules never read each other except through `MetricsCollection`, which is filled after all of
//! them are finalized.

use std::collections::BTreeMap;

use gatk_engine::gatk_report::{Report, Sorting, Table, Value as ReportValue};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::genotype_type::{
    is_monomorphic_in_samples, is_polymorphic_in_samples, GenotypeType,
};
use htsjdk_vcf::genotypes_context::GenotypesContext;
use htsjdk_vcf::variant::{Value, VariantContext};

/// `VariantEvalArgumentCollection.ALL_SAMPLE_NAME` and `ALL_FAMILY_NAME`.
pub const ALL: &str = "all";
/// `VariantEvalEngine.IS_SINGLETON_KEY`.
pub const IS_SINGLETON_KEY: &str = "ISSINGLETON";

/// What the engine refuses, in the reference's classes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EvalError {
    pub class: String,
    pub message: String,
    /// A `UserException` (printed as a user error) rather than a crash.
    pub user: bool,
}

impl EvalError {
    fn user(message: String) -> EvalError {
        EvalError {
            class: "org.broadinstitute.hellbender.exceptions.UserException".to_string(),
            message,
            user: true,
        }
    }

    fn command_line(message: String) -> EvalError {
        EvalError {
            class: "org.broadinstitute.barclay.argparser.CommandLineException".to_string(),
            message,
            user: true,
        }
    }

    fn bad_argument(argument: &str, message: &str) -> EvalError {
        EvalError::command_line(format!("Argument {argument} has a bad value: {message}"))
    }

    fn gatk(message: String) -> EvalError {
        EvalError {
            class: "org.broadinstitute.hellbender.exceptions.GATKException".to_string(),
            message,
            user: false,
        }
    }
}

// ---------------------------------------------------------------------------------------------
// VariantContext helpers
// ---------------------------------------------------------------------------------------------

/// `VariantContext.Type`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum VcType {
    NoVariation,
    Snp,
    Mnp,
    Indel,
    Symbolic,
    Mixed,
}

impl VcType {
    pub const ALL: [VcType; 6] = [
        VcType::NoVariation,
        VcType::Snp,
        VcType::Mnp,
        VcType::Indel,
        VcType::Symbolic,
        VcType::Mixed,
    ];

    pub fn name(self) -> &'static str {
        match self {
            VcType::NoVariation => "NO_VARIATION",
            VcType::Snp => "SNP",
            VcType::Mnp => "MNP",
            VcType::Indel => "INDEL",
            VcType::Symbolic => "SYMBOLIC",
            VcType::Mixed => "MIXED",
        }
    }
}

/// `VariantContext.getType()`.
pub fn vc_type(vc: &VariantContext) -> VcType {
    if vc.alleles.len() <= 1 {
        return VcType::NoVariation;
    }
    let reference = &vc.alleles[0];
    let mut kind: Option<VcType> = None;
    for allele in &vc.alleles[1..] {
        let pair = if allele.is_symbolic() {
            VcType::Symbolic
        } else if reference.len() == allele.len() {
            if allele.len() == 1 {
                VcType::Snp
            } else {
                VcType::Mnp
            }
        } else {
            VcType::Indel
        };
        match kind {
            None => kind = Some(pair),
            Some(existing) if existing != pair => return VcType::Mixed,
            Some(_) => {}
        }
    }
    kind.unwrap_or(VcType::NoVariation)
}

fn is_biallelic(vc: &VariantContext) -> bool {
    vc.alleles.len() == 2
}

fn is_snp(vc: &VariantContext) -> bool {
    vc_type(vc) == VcType::Snp
}

fn is_indel(vc: &VariantContext) -> bool {
    vc_type(vc) == VcType::Indel
}

fn bases(allele: &Allele) -> Vec<u8> {
    allele.display_string().into_bytes()
}

/// `isSimpleIndel()`.
fn is_simple_indel(vc: &VariantContext) -> bool {
    if vc_type(vc) != VcType::Indel || !is_biallelic(vc) {
        return false;
    }
    let reference = bases(&vc.alleles[0]);
    let alternate = bases(&vc.alleles[1]);
    !reference.is_empty()
        && !alternate.is_empty()
        && reference[0] == alternate[0]
        && (reference.len() == 1 || alternate.len() == 1)
}

fn is_simple_insertion(vc: &VariantContext) -> bool {
    is_simple_indel(vc) && vc.alleles[0].len() == 1
}

fn is_simple_deletion(vc: &VariantContext) -> bool {
    is_simple_indel(vc) && vc.alleles[1].len() == 1
}

fn is_complex_indel(vc: &VariantContext) -> bool {
    is_indel(vc) && !is_simple_deletion(vc) && !is_simple_insertion(vc)
}

/// `getIndelLengths()`, `None` for a record that is neither an indel nor mixed.
fn indel_lengths(vc: &VariantContext) -> Option<Vec<i64>> {
    let kind = vc_type(vc);
    if kind != VcType::Indel && kind != VcType::Mixed {
        return None;
    }
    let reference = vc.alleles[0].len() as i64;
    Some(
        vc.alleles[1..]
            .iter()
            .map(|allele| allele.len() as i64 - reference)
            .collect(),
    )
}

fn is_filtered(vc: &VariantContext) -> bool {
    vc.is_filtered()
}

fn attribute<'a>(vc: &'a VariantContext, key: &str) -> Option<&'a Value> {
    vc.attributes
        .iter()
        .find(|(name, _)| name == key)
        .map(|(_, value)| value)
}

fn value_text(value: &Value) -> String {
    crate::reference_confidence_merger::value_to_string(value)
}

/// `getAttributeAsString(key, default)`.
fn attribute_string(vc: &VariantContext, key: &str) -> Option<String> {
    attribute(vc, key).map(|value| match value {
        Value::List(values) => format!(
            "[{}]",
            values.iter().map(value_text).collect::<Vec<_>>().join(", ")
        ),
        other => value_text(other),
    })
}

/// `getAttributeAsInt(key, 0)`, which throws on text that is not a number.
fn attribute_int(vc: &VariantContext, key: &str) -> Result<i64, String> {
    match attribute(vc, key) {
        None | Some(Value::Missing) => Ok(0),
        Some(Value::Int(value)) => Ok(*value),
        Some(Value::List(_)) => Err("java.lang.ClassCastException".to_string()),
        Some(other) => {
            let text = value_text(other);
            if text == "." {
                return Ok(0);
            }
            text.trim()
                .parse()
                .map_err(|_| format!("For input string: \"{text}\""))
        }
    }
}

/// `getAttributeAsDoubleList(key, 0.0)`.
fn attribute_double_list(vc: &VariantContext, key: &str) -> Result<Vec<f64>, String> {
    let values: Vec<Value> = match attribute(vc, key) {
        None => return Ok(Vec::new()),
        Some(Value::List(values)) => values.clone(),
        Some(other) => vec![other.clone()],
    };
    values
        .iter()
        .map(|value| match value {
            Value::Double(number) => Ok(*number),
            Value::Int(number) => Ok(*number as f64),
            Value::Missing => Ok(0.0),
            other => {
                let text = value_text(other);
                if text == "." {
                    Ok(0.0)
                } else {
                    text.trim()
                        .parse()
                        .map_err(|_| format!("For input string: \"{text}\""))
                }
            }
        })
        .collect()
}

/// `getAttributeAsDouble(key, 0.0)`.
fn attribute_double(vc: &VariantContext, key: &str) -> f64 {
    attribute_double_list(vc, key)
        .ok()
        .and_then(|values| values.first().copied())
        .unwrap_or(0.0)
}

/// `BaseUtils.simpleBaseToBaseIndex`.
fn base_index(base: u8) -> i32 {
    match base {
        b'A' | b'a' | b'*' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => -1,
    }
}

/// `BaseUtils.SNPSubstitutionType(base1, base2) == TRANSITION`.
fn is_transition_bases(first: u8, second: u8) -> bool {
    let (a, b) = (base_index(first), base_index(second));
    (a == 0 && b == 2) || (a == 2 && b == 0) || (a == 1 && b == 3) || (a == 3 && b == 1)
}

/// `GATKVariantContextUtils.isTransition(vc)`.
fn is_transition(vc: &VariantContext) -> bool {
    let reference = bases(&vc.alleles[0]);
    let alternate = bases(&vc.alleles[1]);
    is_transition_bases(reference[0], alternate[0])
}

/// `GenotypesContext.subsetToSamples` then, when `rederive`, the alleles the genotypes carry in
/// their original order, with the reference added back when no genotype holds it.
pub fn sub_context_from_samples(
    vc: &VariantContext,
    samples: &[String],
    rederive: bool,
) -> VariantContext {
    let all_present = vc
        .genotypes
        .iter()
        .all(|genotype| samples.contains(&genotype.sample_name));
    if all_present && !rederive {
        return vc.clone();
    }
    let mut subset = Vec::new();
    for sample in samples {
        if let Some(genotype) = vc.genotype(sample) {
            subset.push(genotype.clone());
        }
    }
    let mut out = vc.clone();
    if rederive {
        let mut from_genotypes: Vec<Allele> = Vec::new();
        let mut added_reference = false;
        for genotype in &subset {
            for allele in &genotype.alleles {
                added_reference = added_reference || allele.is_reference();
                if !allele.is_no_call() && !from_genotypes.contains(allele) {
                    from_genotypes.push(allele.clone());
                }
            }
        }
        if !added_reference && !from_genotypes.contains(&vc.alleles[0]) {
            from_genotypes.push(vc.alleles[0].clone());
        }
        out.alleles = vc
            .alleles
            .iter()
            .filter(|allele| from_genotypes.contains(allele))
            .cloned()
            .collect();
    }
    out.genotypes = GenotypesContext::new(subset);
    out
}

/// `ensureAnnotations(vc, vcsub)`: `ISSINGLETON`, and AC, AF and AN when the subset lacks any.
pub fn ensure_annotations(vc: &VariantContext, sub: &VariantContext) -> VariantContext {
    let allele_count = |record: &VariantContext| -> i64 {
        record
            .genotypes
            .iter()
            .map(|genotype| match genotype.genotype_type() {
                GenotypeType::Het => 1,
                GenotypeType::HomVar => 2,
                _ => 0,
            })
            .sum()
    };
    let original = allele_count(vc);
    let new = allele_count(sub);
    let singleton = original == new && new == 1;
    let has_counts = ["AC", "AF", "AN"]
        .iter()
        .all(|key| attribute(sub, key).is_some());
    if !singleton && has_counts {
        return sub.clone();
    }
    let mut out = sub.clone();
    if singleton {
        put(&mut out.attributes, IS_SINGLETON_KEY, Value::Bool(true));
    }
    if !has_counts {
        let counts = htsjdk_vcf::chromosome_counts::calculate_chromosome_counts(&out, true, &[]);
        let computed = counts.attributes();
        if counts.remove_stale {
            out.attributes
                .retain(|(key, _)| !matches!(key.as_str(), "AC" | "AF" | "AN"));
        } else {
            if computed.iter().all(|(key, _)| key != "AC") {
                out.attributes.retain(|(key, _)| key != "AC" && key != "AF");
            }
            for (key, value) in computed {
                put(&mut out.attributes, &key, value);
            }
        }
    }
    out
}

fn put(attributes: &mut Vec<(String, Value)>, key: &str, value: Value) {
    match attributes.iter_mut().find(|(name, _)| name == key) {
        Some((_, slot)) => *slot = value,
        None => attributes.push((key.to_string(), value)),
    }
}

// ---------------------------------------------------------------------------------------------
// Stratifications
// ---------------------------------------------------------------------------------------------

/// One state of a stratifier: a string, or an integer for the `%d` columns.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum State {
    Str(String),
    Int(i64),
}

impl State {
    fn text(&self) -> String {
        match self {
            State::Str(text) => text.clone(),
            State::Int(number) => number.to_string(),
        }
    }

    fn report(&self) -> ReportValue {
        match self {
            State::Str(text) => ReportValue::Str(text.clone()),
            State::Int(number) => ReportValue::Int(*number),
        }
    }
}

fn s(text: &str) -> State {
    State::Str(text.to_string())
}

/// `AlleleFrequency`'s two scales.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AfScale {
    Linear,
    Logarithmic,
}

/// Which stratifier, with what it needs to answer.
#[derive(Debug, Clone)]
pub enum StratKind {
    AlleleCount { nchrom: i64 },
    AlleleFrequency { scale: AfScale, use_comp: bool },
    CompFeatureInput,
    Contig,
    CpG,
    Degeneracy,
    EvalFeatureInput,
    Family,
    Filter,
    FilterType,
    FunctionalClass,
    IndelSize,
    IntervalStratification,
    JexlExpression,
    Novelty,
    OneBPIndel,
    Sample,
    SnpEffPositionModifier,
    TandemRepeat,
    VariantType,
}

/// `VariantStratifier`: its name, its states in declaration order and its column format.
#[derive(Debug, Clone)]
pub struct Stratifier {
    pub name: &'static str,
    pub states: Vec<State>,
    pub kind: StratKind,
    pub format: &'static str,
}

/// Every stratifier the reference finds by reflection, with whether it is standard or required.
pub const STRATIFIER_NAMES: &[(&str, bool, bool)] = &[
    ("AlleleCount", false, false),
    ("AlleleFrequency", false, false),
    ("CompFeatureInput", false, true),
    ("Contig", false, false),
    ("CpG", false, false),
    ("Degeneracy", false, false),
    ("EvalFeatureInput", false, true),
    ("Family", false, false),
    ("Filter", false, false),
    ("FilterType", false, false),
    ("FunctionalClass", false, false),
    ("IndelSize", false, false),
    ("IntervalStratification", false, false),
    ("JexlExpression", true, false),
    ("Novelty", true, false),
    ("OneBPIndel", false, false),
    ("Sample", false, false),
    ("SnpEffPositionModifier", false, false),
    ("TandemRepeat", false, false),
    ("VariantType", false, false),
];

/// Every evaluator the reference finds, with whether it is standard.
pub const EVALUATOR_NAMES: &[(&str, bool)] = &[
    ("CompOverlap", true),
    ("CountVariants", true),
    ("GenotypeFilterSummary", false),
    ("IndelLengthHistogram", true),
    ("IndelSummary", true),
    ("MendelianViolationEvaluator", false),
    ("MetricsCollection", false),
    ("MultiallelicSummary", true),
    ("PrintMissingComp", false),
    ("ThetaVariantEvaluator", false),
    ("TiTvVariantEvaluator", true),
    ("ValidationReport", true),
    ("VariantAFEvaluator", false),
    ("VariantSummary", true),
];

/// One named JEXL select expression.
pub struct SelectExpression {
    pub name: String,
    pub expression: gatk_engine::jexl::Expression,
}

/// A feature the interval and CNV queries read: its contig and 1-based closed span.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Feature {
    pub contig: String,
    pub start: i64,
    pub end: i64,
}

fn overlapping<'a>(
    features: &'a [Feature],
    contig: &str,
    start: i64,
    end: i64,
) -> Vec<&'a Feature> {
    features
        .iter()
        .filter(|feature| feature.contig == contig && feature.start <= end && start <= feature.end)
        .collect()
}

/// The arguments the engine reads.
pub struct Arguments {
    pub eval_names: Vec<String>,
    /// `comps`: the provided comps' names, then `dbsnp` when given.
    pub comp_names: Vec<String>,
    /// Indices into `comp_names` of the comps that are known.
    pub known: Vec<usize>,
    pub strats_to_use: Vec<String>,
    pub no_standard_strats: bool,
    pub modules_to_use: Vec<String>,
    pub no_standard_modules: bool,
    pub ploidy: i64,
    pub require_strict_allele_match: bool,
    pub keep_ac0: bool,
    pub merge_evals: bool,
    pub num_samples_from_argument: i64,
    pub af_scale: AfScale,
    pub use_comp_af: bool,
    pub samples_for_evaluation: Vec<String>,
    pub all_eval_samples: Vec<String>,
    pub selects: Vec<SelectExpression>,
    pub traversal: Vec<(String, i64, i64)>,
    pub eval_filter_names: Vec<String>,
    pub strat_intervals: Option<Vec<Feature>>,
    /// The `--strat-intervals` path when it has no index: `FeatureDataSource` asks for one at the
    /// first query, which is the first eval `IntervalStratification` sees.
    pub strat_intervals_unindexed: Option<String>,
    pub known_cnvs: Option<Vec<Feature>>,
    /// The `--known-cnvs` path when it has no index, refused at the first query as
    /// [`Arguments::strat_intervals_unindexed`] is.
    pub known_cnvs_unindexed: Option<String>,
    pub gold_standard: Option<Vec<VariantContext>>,
}

impl Arguments {
    fn ignore_ac0(&self) -> bool {
        !self.keep_ac0
    }
}

/// The site's reference, as the stratifiers read it.
pub trait ReferenceBases {
    fn bases(&mut self, contig: &str, start: i64, end: i64) -> Vec<u8>;
}

/// What a stratifier or an evaluator reads beside the records: the site's grouped records and
/// the reference.
struct SiteContext<'a> {
    /// Per comp track, the records starting here.
    comps_here: Vec<Vec<&'a VariantContext>>,
    reference_start_bases: Vec<u8>,
    tandem_context: Option<Vec<u8>>,
}

impl Stratifier {
    #[allow(clippy::too_many_arguments)]
    fn relevant_states(
        &self,
        arguments: &Arguments,
        site: &SiteContext<'_>,
        comp: Option<&VariantContext>,
        comp_name: Option<&str>,
        eval: Option<&VariantContext>,
        eval_name: &str,
        sample: Option<&str>,
        family: Option<&str>,
    ) -> Result<Vec<State>, EvalError> {
        Ok(match &self.kind {
            StratKind::AlleleCount { nchrom } => {
                let Some(eval) = eval else {
                    return Ok(Vec::new());
                };
                let mut ac: i64 = 0;
                if is_biallelic(eval) {
                    if attribute(eval, "MLEAC").is_some() {
                        if let Ok(value) = attribute_int(eval, "MLEAC") {
                            ac = value.min(*nchrom);
                        }
                    } else if attribute(eval, "AC").is_some() {
                        if let Ok(value) = attribute_int(eval, "AC") {
                            ac = value;
                        }
                    }
                }
                if ac == 0 && eval.is_variant() {
                    for allele in &eval.alleles[1..] {
                        let called =
                            htsjdk_vcf::chromosome_counts::called_chr_count_for(eval, allele, &[])
                                as i64;
                        ac = ac.max(called);
                    }
                }
                if ac > *nchrom {
                    return Err(EvalError::user(format!(
                        "The AC value ({ac}) at position {}:{} is larger than the number of \
                         chromosomes over all samples ({nchrom})",
                        eval.contig, eval.start
                    )));
                }
                vec![State::Int(ac)]
            }
            StratKind::AlleleFrequency { scale, use_comp } => {
                let Some(eval) = eval else {
                    return Ok(Vec::new());
                };
                let frequency = (|| -> Result<f64, String> {
                    let mut frequency = max_of(&attribute_double_list(eval, "AF")?)?;
                    if *use_comp {
                        frequency = match comp {
                            Some(comp) => max_of(&attribute_double_list(comp, "AF")?)?,
                            None => 0.0,
                        };
                    }
                    Ok(frequency)
                })();
                let Ok(frequency) = frequency else {
                    return Ok(Vec::new());
                };
                match scale {
                    AfScale::Linear => vec![State::Str(gatk_engine::java_format::format_decimals(
                        5.0 * round_to_n_decimal_places(frequency / 5.0, 3),
                        3,
                    ))],
                    AfScale::Logarithmic => {
                        vec![State::Str(
                            logit_bucket(frequency + 10f64.powi(-6)).to_string(),
                        )]
                    }
                }
            }
            StratKind::CompFeatureInput => vec![s(comp_name.unwrap_or("none"))],
            StratKind::Contig => match eval {
                Some(eval) => vec![s(ALL), State::Str(eval.contig.clone())],
                None => Vec::new(),
            },
            StratKind::CpG => {
                let is_cpg = site.reference_start_bases.starts_with(b"CG");
                vec![s(ALL), s(if is_cpg { "CpG" } else { "non_CpG" })]
            }
            StratKind::Degeneracy => degeneracy_states(eval),
            StratKind::EvalFeatureInput => vec![s(eval_name)],
            StratKind::Family => vec![State::Str(family.unwrap_or("null").to_string())],
            StratKind::Filter => {
                let mut states = vec![s("raw")];
                if let Some(eval) = eval {
                    states.push(s(if is_filtered(eval) {
                        "filtered"
                    } else {
                        "called"
                    }));
                }
                states
            }
            StratKind::FilterType => {
                let Some(eval) = eval else {
                    return Ok(Vec::new());
                };
                if is_filtered(eval) {
                    eval.filters
                        .clone()
                        .unwrap_or_default()
                        .into_iter()
                        .map(State::Str)
                        .collect()
                } else {
                    vec![s("PASS")]
                }
            }
            StratKind::FunctionalClass => functional_class_states(eval),
            StratKind::IndelSize => {
                let Some(eval) = eval else {
                    return Ok(Vec::new());
                };
                if !is_indel(eval) || !is_biallelic(eval) {
                    return Ok(Vec::new());
                }
                let mut length = 0i64;
                if is_simple_insertion(eval) {
                    length = eval.alleles[1].len() as i64;
                } else if is_simple_deletion(eval) {
                    length = -(eval.alleles[0].len() as i64);
                }
                vec![State::Int(length.clamp(-100, 100))]
            }
            StratKind::IntervalStratification => {
                let Some(eval) = eval else {
                    return Ok(Vec::new());
                };
                if let Some(path) = &arguments.strat_intervals_unindexed {
                    return Err(EvalError::user(format!(
                        "Input {path} must support random access to enable queries by interval. \
                         If it's a file, please index it using the bundled tool IndexFeatureFile"
                    )));
                }
                let features = arguments.strat_intervals.as_deref().unwrap_or(&[]);
                if overlapping(features, &eval.contig, eval.start, eval.stop).is_empty() {
                    vec![s(ALL), s("outside.intervals")]
                } else {
                    vec![s(ALL), s("overlaps.intervals")]
                }
            }
            StratKind::JexlExpression => {
                let mut states = vec![s("none")];
                for select in &arguments.selects {
                    if let Some(eval) = eval {
                        if jexl_matches(eval, select)? {
                            states.push(State::Str(select.name.clone()));
                        }
                    }
                }
                states
            }
            StratKind::Novelty => {
                if let Some(eval) = eval {
                    let eval_type = vc_type(eval);
                    for index in &arguments.known {
                        for known in &site.comps_here[*index] {
                            if eval_type == vc_type(known) || eval_type == VcType::NoVariation {
                                return Ok(vec![s(ALL), s("known")]);
                            }
                        }
                    }
                }
                vec![s(ALL), s("novel")]
            }
            StratKind::OneBPIndel => match eval {
                Some(eval) if is_indel(eval) => {
                    let two_plus = indel_lengths(eval)
                        .unwrap_or_default()
                        .iter()
                        .any(|length| length.abs() > 1);
                    if two_plus {
                        vec![s(ALL), s("two.plus.bp")]
                    } else {
                        vec![s(ALL), s("one.bp")]
                    }
                }
                _ => vec![s(ALL), s("one.bp"), s("two.plus.bp")],
            },
            StratKind::Sample => vec![State::Str(sample.unwrap_or("null").to_string())],
            StratKind::SnpEffPositionModifier => snpeff_position_states(eval)?,
            StratKind::TandemRepeat => match eval {
                Some(eval) if is_indel(eval) => {
                    let context = site.tandem_context.clone().unwrap_or_default();
                    if is_tandem_repeat(eval, &context) {
                        vec![s(ALL), s("is.repeat")]
                    } else {
                        vec![s(ALL), s("not.repeat")]
                    }
                }
                _ => vec![s(ALL), s("is.repeat"), s("not.repeat")],
            },
            StratKind::VariantType => match eval {
                Some(eval) => vec![s(vc_type(eval).name())],
                None => Vec::new(),
            },
        })
    }
}

/// `Collections.max`, which throws on an empty list.
fn max_of(values: &[f64]) -> Result<f64, String> {
    let mut best: Option<f64> = None;
    for value in values {
        best = Some(match best {
            None => *value,
            Some(best) if value.total_cmp(&best) == std::cmp::Ordering::Greater => *value,
            Some(best) => best,
        });
    }
    best.ok_or_else(|| "NoSuchElementException".to_string())
}

/// `MathUtils.roundToNDecimalPlaces`: `Math.round((in + Math.ulp(in)) * 10^n) / 10^n`.
fn round_to_n_decimal_places(value: f64, places: i32) -> f64 {
    let multiplier = 10f64.powi(places);
    let ulp = if value == 0.0 {
        f64::from_bits(1)
    } else {
        let bits = value.abs().to_bits();
        f64::from_bits(bits + 1) - value.abs()
    };
    (crate::reference_confidence_merger::java_round((value + ulp) * multiplier)) as f64 / multiplier
}

/// `AlleleFrequency.getLogitBucket`: a float score, rounded as `Math.round(float)` does.
fn logit_bucket(frequency: f64) -> i32 {
    let score = (-10.0 * jmath::math::log10(frequency / (1.0 - frequency))) as f32;
    let rounded = if score.is_nan() {
        0
    } else {
        (score + 0.5).floor() as i32
    };
    rounded.clamp(-30, 30)
}

fn degeneracy_table(amino_acid: &str, frame: i64) -> Option<&'static str> {
    let codons: &[&str] = match amino_acid {
        "Ile" => &["ATT", "ATC", "ATA"],
        "Leu" => &["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"],
        "Val" => &["GTT", "GTC", "GTA", "GTG"],
        "Phe" => &["TTT", "TTC"],
        "Met" => &["ATG"],
        "Cys" => &["TGT", "TGC"],
        "Ala" => &["GCT", "GCC", "GCA", "GCG"],
        "Gly" => &["GGT", "GGC", "GGA", "GGG"],
        "Pro" => &["CCT", "CCC", "CCA", "CCG"],
        "Thr" => &["ACT", "ACC", "ACA", "ACG"],
        "Ser" => &["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
        "Tyr" => &["TAT", "TAC"],
        "Trp" => &["TGG"],
        "Glu" => &["CAA", "CAG"],
        "Asn" => &["AAT", "AAC"],
        "His" => &["CAT", "CAC"],
        "Gln" => &["GAA", "GAG"],
        "Asp" => &["GAT", "GAC"],
        "Lys" => &["AAA", "AAG"],
        "Arg" => &["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "Stop" => &["TAA", "TAG", "TGA"],
        _ => return None,
    };
    if !(0..3).contains(&frame) {
        return None;
    }
    let mut seen: Vec<u8> = Vec::new();
    for codon in codons {
        let base = codon.as_bytes()[frame as usize];
        if !seen.contains(&base) {
            seen.push(base);
        }
    }
    Some(match seen.len() {
        2 => "2-fold",
        3 => "3-fold",
        4 => "4-fold",
        6 => "6-fold",
        _ => "1-fold",
    })
}

fn degeneracy_states(eval: Option<&VariantContext>) -> Vec<State> {
    let mut states = vec![s(ALL)];
    let Some(eval) = eval else { return states };
    if !eval.is_variant() {
        return states;
    }
    let mut kind: Option<String> = None;
    let mut amino_acid: Option<String> = None;
    let mut frame: Option<i64> = None;
    if attribute(eval, "refseq.functionalClass").is_some() {
        amino_acid = attribute_string(eval, "refseq.variantAA");
        frame = Some(attribute_int(eval, "refseq.frame").unwrap_or(0));
    } else if attribute(eval, "refseq.functionalClass_1").is_some() {
        let mut id = 1;
        loop {
            let key = format!("refseq.functionalClass_{id}");
            if let Some(new_kind) = attribute_string(eval, &key) {
                let better = match kind.as_deref() {
                    None => true,
                    Some("silent") => new_kind != "silent",
                    Some("missense") => new_kind == "nonsense",
                    Some(_) => false,
                };
                if better {
                    kind = Some(new_kind);
                    amino_acid = attribute_string(eval, &format!("refseq.variantAA_{id}"));
                    if amino_acid.is_some() {
                        let frame_key = format!("refseq.frame_{id}");
                        if attribute(eval, &frame_key).is_some() {
                            frame = Some(attribute_int(eval, &frame_key).unwrap_or(0));
                        }
                    }
                }
            }
            id += 1;
            if attribute(eval, &key).is_none() {
                break;
            }
        }
    }
    if let (Some(amino_acid), Some(frame)) = (amino_acid, frame) {
        if let Some(fold) = degeneracy_table(&amino_acid, frame) {
            states.push(s(fold));
        } else if [
            "Ile", "Leu", "Val", "Phe", "Met", "Cys", "Ala", "Gly", "Pro", "Thr", "Ser", "Tyr",
            "Trp", "Glu", "Asn", "His", "Gln", "Asp", "Lys", "Arg", "Stop",
        ]
        .contains(&amino_acid.as_str())
        {
            // `degeneracies.get(aa).get(frame)` is null for a frame outside 0..2: a null state.
            states.push(s("null"));
        }
    }
    states
}

fn functional_class_states(eval: Option<&VariantContext>) -> Vec<State> {
    let mut states = vec![s(ALL)];
    let Some(eval) = eval else { return states };
    if !eval.is_variant() {
        return states;
    }
    let rank = |name: &str| match name {
        "silent" => Some(0),
        "missense" => Some(1),
        "nonsense" => Some(2),
        _ => None,
    };
    let mut kind: Option<&'static str> = None;
    let names = ["silent", "missense", "nonsense"];
    if attribute(eval, "refseq.functionalClass").is_some() {
        if let Some(text) = attribute_string(eval, "refseq.functionalClass") {
            kind = rank(&text).map(|index| names[index]);
        }
    } else if attribute(eval, "refseq.functionalClass_1").is_some() {
        let mut id = 1;
        loop {
            let key = format!("refseq.functionalClass_{id}");
            if let Some(text) = attribute_string(eval, &key) {
                if !text.eq_ignore_ascii_case("null") {
                    if let Some(new_index) = rank(&text) {
                        let better = match kind {
                            None => true,
                            Some("silent") => new_index != 0,
                            Some("missense") => new_index == 2,
                            Some(_) => false,
                        };
                        if better {
                            kind = Some(names[new_index]);
                        }
                    }
                }
            }
            id += 1;
            if attribute(eval, &key).is_none() {
                break;
            }
        }
    } else if let Some(text) = attribute_string(eval, "SNPEFF_FUNCTIONAL_CLASS") {
        kind = match text.as_str() {
            "NONSENSE" => Some("nonsense"),
            "MISSENSE" => Some("missense"),
            "SILENT" => Some("silent"),
            _ => None,
        };
    }
    if let Some(kind) = kind {
        states.push(s(kind));
    }
    states
}

/// `SnpEffUtil`'s effect graph: each child's parent.
fn snpeff_parent(effect: &str) -> Option<&'static str> {
    Some(match effect {
        "UPSTREAM" | "DOWNSTREAM" | "INTERGENIC_CONSERVED" => "INTERGENIC",
        "INTRON_CONSERVED" | "SPLICE_SITE_ACCEPTOR" | "SPLICE_SITE_DONOR" => "INTRON",
        "EXON_DELETED" | "SYNONYMOUS_CODING" | "NON_SYNONYMOUS_CODING" => "CDS",
        "SYNONYMOUS_STOP" | "SYNONYMOUS_START" => "SYNONYMOUS_CODING",
        "START_LOST"
        | "STOP_GAINED"
        | "STOP_LOST"
        | "CODON_CHANGE"
        | "CODON_INSERTION"
        | "CODON_DELETION"
        | "CODON_CHANGE_PLUS_CODON_DELETION"
        | "CODON_CHANGE_PLUS_CODON_INSERTION"
        | "FRAME_SHIFT" => "NON_SYNONYMOUS_CODING",
        "UTR_5_DELETED" | "START_GAINED" => "UTR_5_PRIME",
        "UTR_3_DELETED" => "UTR_3_PRIME",
        "UTR_5_PRIME" | "UTR_3_PRIME" | "CDS" => "EXON",
        "INTRON" | "EXON" => "TRANSCRIPT",
        "TRANSCRIPT" | "REGULATION" => "GENE",
        "GENE" | "INTERGENIC" => "CHROMOSOME",
        _ => return None,
    })
}

const SNPEFF_EFFECTS: &[&str] = &[
    "SPLICE_SITE_ACCEPTOR",
    "SPLICE_SITE_DONOR",
    "START_LOST",
    "EXON_DELETED",
    "FRAME_SHIFT",
    "STOP_GAINED",
    "STOP_LOST",
    "NON_SYNONYMOUS_CODING",
    "CODON_CHANGE",
    "CODON_INSERTION",
    "CODON_CHANGE_PLUS_CODON_INSERTION",
    "CODON_DELETION",
    "CODON_CHANGE_PLUS_CODON_DELETION",
    "UTR_5_DELETED",
    "UTR_3_DELETED",
    "SYNONYMOUS_START",
    "NON_SYNONYMOUS_START",
    "START_GAINED",
    "SYNONYMOUS_CODING",
    "SYNONYMOUS_STOP",
    "NON_SYNONYMOUS_STOP",
    "NONE",
    "CHROMOSOME",
    "CUSTOM",
    "CDS",
    "GENE",
    "TRANSCRIPT",
    "EXON",
    "INTRON_CONSERVED",
    "UTR_5_PRIME",
    "UTR_3_PRIME",
    "DOWNSTREAM",
    "INTRAGENIC",
    "INTERGENIC",
    "INTERGENIC_CONSERVED",
    "UPSTREAM",
    "REGULATION",
    "INTRON",
];

fn snpeff_is_subtype(child: &str, parent: &str) -> bool {
    let mut current = Some(child);
    while let Some(effect) = current {
        if effect == parent {
            return true;
        }
        current = snpeff_parent(effect);
    }
    false
}

fn snpeff_position_states(eval: Option<&VariantContext>) -> Result<Vec<State>, EvalError> {
    let mut states = Vec::new();
    let Some(eval) = eval else { return Ok(states) };
    if !eval.is_variant() {
        return Ok(states);
    }
    let Some(effect) = attribute_string(eval, "SNPEFF_EFFECT") else {
        return Ok(states);
    };
    if !SNPEFF_EFFECTS.contains(&effect.as_str()) {
        return Err(EvalError {
            class: "java.lang.IllegalArgumentException".to_string(),
            message: format!(
                "No enum constant org.broadinstitute.hellbender.tools.walkers.varianteval.util.\
                 SnpEffUtil.EffectType.{effect}"
            ),
            user: false,
        });
    }
    if snpeff_is_subtype(&effect, "EXON") {
        states.push(s("GENE"));
    }
    if snpeff_is_subtype(&effect, "CDS") {
        states.push(s("CODING_REGION"));
    }
    if snpeff_is_subtype(&effect, "STOP_GAINED") {
        states.push(s("STOP_GAINED"));
    }
    if snpeff_is_subtype(&effect, "STOP_LOST") {
        states.push(s("STOP_LOST"));
    }
    if snpeff_is_subtype(&effect, "SPLICE_SITE_ACCEPTOR")
        || snpeff_is_subtype(&effect, "SPLICE_SITE_DONOR")
    {
        states.push(s("SPLICE_SITE"));
    }
    Ok(states)
}

/// `GATKVariantContextUtils.isTandemRepeat(vc, refBasesStartingAtVCWithPad)`.
fn is_tandem_repeat(vc: &VariantContext, with_pad: &[u8]) -> bool {
    if !is_indel(vc) {
        return false;
    }
    let without_pad: Vec<u8> = with_pad.iter().skip(1).copied().collect();
    let reference = bases(&vc.alleles[0]);
    for allele in &vc.alleles[1..] {
        let alternate = bases(allele);
        let (short, long) = if reference.len() <= alternate.len() {
            (&reference, &alternate)
        } else {
            (&alternate, &reference)
        };
        if !long.starts_with(short) {
            return false;
        }
        let repeated = if reference.len() > alternate.len() {
            bases_are_repeated(&reference, &alternate, &without_pad, 2)
        } else {
            bases_are_repeated(&alternate, &reference, &without_pad, 1)
        };
        if !repeated {
            return false;
        }
    }
    true
}

fn bases_are_repeated(long: &[u8], short: &[u8], reference: &[u8], matches: usize) -> bool {
    let potential = &long[short.len()..];
    for i in 0..matches {
        let start = i * potential.len();
        let end = (i + 1) * potential.len();
        if reference.len() < end {
            return false;
        }
        if &reference[start..end] != potential {
            return false;
        }
    }
    true
}

/// `VariantContextUtils.match(eval, exp)` with the default treatment of a missing value: false.
fn jexl_matches(vc: &VariantContext, select: &SelectExpression) -> Result<bool, EvalError> {
    use gatk_engine::jexl::{JexlError, Value as JexlValue};
    let mut context = gatk_engine::jexl::Context::new();
    context.insert("CHROM".to_string(), JexlValue::Str(vc.contig.clone()));
    context.insert("POS".to_string(), JexlValue::Int(vc.start as i32));
    context.insert(
        "TYPE".to_string(),
        JexlValue::Str(vc_type(vc).name().to_string()),
    );
    context.insert(
        "QUAL".to_string(),
        JexlValue::Double(-10.0 * vc.log10_p_error),
    );
    context.insert(
        "N_ALLELES".to_string(),
        JexlValue::Int(vc.alleles.len() as i32),
    );
    context.insert(
        "FILTER".to_string(),
        JexlValue::Str(if is_filtered(vc) { "1" } else { "0" }.to_string()),
    );
    for (key, value) in &vc.attributes {
        let converted = match value {
            Value::Int(number) => JexlValue::Int(*number as i32),
            Value::Double(number) => JexlValue::Double(*number),
            Value::Bool(flag) => JexlValue::Bool(*flag),
            other => match other.format() {
                Some(text) => JexlValue::Str(text),
                None => continue,
            },
        };
        context.insert(key.clone(), converted);
    }
    for filter in vc.filters.iter().flatten() {
        context
            .entry(filter.clone())
            .or_insert_with(|| JexlValue::Str("1".to_string()));
    }
    // `JEXLMap.evaluateExpression`: `(Boolean) exp.evaluate(context)`, a null or an unknown name
    // is `TREAT_AS_MISMATCH`, and any other engine complaint becomes an IllegalArgumentException.
    match select.expression.evaluate(&context) {
        Ok(JexlValue::Bool(value)) => Ok(value),
        Ok(JexlValue::Null) | Err(JexlError::UndefinedVariable(_)) => Ok(false),
        Ok(other) => {
            let class = match other {
                JexlValue::Int(_) => "java.lang.Integer",
                JexlValue::Long(_) => "java.lang.Long",
                JexlValue::Float(_) => "java.lang.Float",
                JexlValue::Double(_) => "java.lang.Double",
                _ => "java.lang.String",
            };
            Err(EvalError {
                class: "java.lang.ClassCastException".to_string(),
                message: format!(
                    "class {class} cannot be cast to class java.lang.Boolean ({class} and \
                     java.lang.Boolean are in module java.base of loader 'bootstrap')"
                ),
                user: false,
            })
        }
        Err(JexlError::Unsupported(what)) => Err(EvalError {
            class: "gatk_rs::PortLimitation".to_string(),
            message: format!("JEXL construct not ported: {what}"),
            user: false,
        }),
        Err(_) => Err(EvalError {
            class: "java.lang.IllegalArgumentException".to_string(),
            message: format!("Invalid JEXL expression detected for {}", select.name),
            user: false,
        }),
    }
}

// ---------------------------------------------------------------------------------------------
// Evaluators
// ---------------------------------------------------------------------------------------------

/// `Utils.formattedPercent`.
fn formatted_percent(x: i64, total: i64) -> String {
    if total == 0 {
        "NA".to_string()
    } else {
        gatk_engine::java_format::format_decimals((100.0 * x as f64) / total as f64, 2)
    }
}

/// `Utils.formattedRatio`.
fn formatted_ratio(num: i64, denom: i64) -> String {
    if denom == 0 {
        "NA".to_string()
    } else {
        gatk_engine::java_format::format_decimals(num as f64 / denom as f64, 2)
    }
}

fn rate(n: i64, d: i64) -> f64 {
    n as f64 / (1.0 * d.max(1) as f64)
}

fn inverse_rate(n: i64, d: i64) -> i64 {
    if n == 0 {
        0
    } else {
        d / n.max(1)
    }
}

fn ratio(num: i64, denom: i64) -> f64 {
    num as f64 / denom.max(1) as f64
}

fn was_singleton(vc: &VariantContext) -> bool {
    matches!(attribute(vc, IS_SINGLETON_KEY), Some(Value::Bool(true)))
        || attribute_string(vc, IS_SINGLETON_KEY)
            .is_some_and(|text| text.eq_ignore_ascii_case("true"))
}

/// One cell of a module's row: its column format and the value.
type Cell = (&'static str, &'static str, ReportValue);

/// One evaluation module's state.
#[derive(Debug, Clone)]
pub enum Evaluator {
    CompOverlap {
        n_eval: i64,
        n_at_comp: i64,
        n_concordant: i64,
        comp_rate: f64,
        concordant_rate: f64,
        novel_sites: i64,
    },
    CountVariants(Box<CountVariants>),
    GenotypeFilterSummary {
        called: i64,
        no_call_or_filtered: i64,
    },
    IndelLengthHistogram {
        counts: BTreeMap<i64, i64>,
        n_indels: i64,
        results: Vec<(i64, f64)>,
    },
    IndelSummary(Box<IndelSummary>),
    MetricsCollection(Box<MetricsCollection>),
    MultiallelicSummary(Box<MultiallelicSummary>),
    PrintMissingComp {
        missing: i64,
    },
    ThetaVariantEvaluator(Box<Theta>),
    TiTvVariantEvaluator(Box<TiTv>),
    ValidationReport(Box<ValidationReport>),
    VariantAFEvaluator(Box<VariantAf>),
    VariantSummary(Box<VariantSummary>),
}

#[derive(Debug, Clone, Default)]
pub struct CountVariants {
    processed: i64,
    called: i64,
    reference: i64,
    variant: i64,
    variant_rate: f64,
    variant_rate_per_bp: f64,
    snps: i64,
    mnps: i64,
    insertions: i64,
    deletions: i64,
    complex: i64,
    symbolic: i64,
    mixed: i64,
    no_calls: i64,
    hets: i64,
    hom_ref: i64,
    hom_var: i64,
    singletons: i64,
    hom_derived: i64,
    heterozygosity: f64,
    heterozygosity_per_bp: f64,
    het_hom_ratio: f64,
    indel_rate: f64,
    indel_rate_per_bp: f64,
    insertion_deletion_ratio: f64,
}

#[derive(Debug, Clone, Default)]
pub struct IndelSummary {
    snps: i64,
    singleton_snps: i64,
    indels: i64,
    singleton_indels: i64,
    matching_gold: i64,
    gold_rate: String,
    indel_sites: i64,
    multiallelic_indel_sites: i64,
    percent_multiallelic: String,
    snp_to_indel: String,
    snp_to_indel_singletons: String,
    novel_indels: i64,
    novelty_rate: String,
    insertions: i64,
    deletions: i64,
    insertion_to_deletion: String,
    large_deletions: i64,
    large_insertions: i64,
    insertion_to_deletion_large: String,
    frameshifting: i64,
    in_frame: i64,
    frameshift_rate: String,
    snp_het_to_hom: String,
    indel_het_to_hom: String,
    snp_hets: i64,
    snp_homs: i64,
    indel_hets: i64,
    indel_homs: i64,
    insertion_by_length: [i64; 4],
    deletion_by_length: [i64; 4],
    ratio_insertions: String,
    ratio_deletions: String,
}

#[derive(Debug, Clone, Default)]
pub struct MetricsCollection {
    concordant_rate: f64,
    snps: i64,
    snp_loci: i64,
    indels: i64,
    indel_loci: i64,
    indel_ratio: Option<String>,
    indel_ratio_loci: f64,
    titv: f64,
}

#[derive(Debug, Clone, Default)]
pub struct MultiallelicSummary {
    processed: i64,
    snps: i64,
    multi_snps: i64,
    processed_multi_snp_ratio: f64,
    variant_multi_snp_ratio: f64,
    indels: i64,
    multi_indels: i64,
    processed_multi_indel_ratio: f64,
    variant_multi_indel_ratio: f64,
    ti: i64,
    tv: i64,
    titv: f64,
    known_partial: i64,
    known_complete: i64,
    snp_novelty: String,
}

#[derive(Debug, Clone, Default)]
pub struct Theta {
    avg_het: f64,
    avg_avg_diffs: f64,
    total_het: f64,
    total_avg_diffs: f64,
    theta_region_num_sites: f64,
    num_sites: f64,
}

#[derive(Debug, Clone, Default)]
pub struct TiTv {
    ti: i64,
    tv: i64,
    ratio: f64,
    ti_comp: i64,
    tv_comp: i64,
    ratio_standard: f64,
    ti_derived: i64,
    tv_derived: i64,
    ratio_derived: f64,
}

#[derive(Debug, Clone, Default)]
pub struct ValidationReport {
    counts: [[i64; 4]; 4],
    n_comp: i64,
    tp: i64,
    fp: i64,
    fn_: i64,
    tn: i64,
    sensitivity: f64,
    specificity: f64,
    ppv: f64,
    fdr: f64,
    comp_filtered: i64,
    different_alleles: i64,
}

#[derive(Debug, Clone, Default)]
pub struct VariantAf {
    avg: f64,
    called: i64,
    het: i64,
    hom_var: i64,
    hom_ref: i64,
    sum: f64,
}

/// `VariantSummary`: its per-type, per-sample maps keep the sample order the reference's
/// `HashMap` iterates in, because the ratio averages sum doubles in that order.
#[derive(Debug, Clone)]
pub struct VariantSummary {
    samples_in_hash_order: Vec<String>,
    n_samples: i64,
    processed: i64,
    snps: i64,
    titv: f64,
    snp_novelty: String,
    snps_per_sample: i64,
    titv_per_sample: f64,
    snp_dp_per_sample: f64,
    indels: i64,
    indel_novelty: String,
    indels_per_sample: i64,
    indel_dp_per_sample: f64,
    svs: i64,
    sv_novelty: String,
    svs_per_sample: i64,
    /// Per type (SNP, INDEL, CNV): per key (a sample, or `ALL`), a count. `ALL` is kept apart
    /// because the reference compares it by identity.
    counts_per_sample: [TypeCounts; 3],
    transitions: [TypeCounts; 3],
    transversions: [TypeCounts; 3],
    all_variant_counts: [TypeCounts; 3],
    known_variant_counts: [TypeCounts; 3],
    depth: [TypeCounts; 3],
}

#[derive(Debug, Clone, Default)]
struct TypeCounts {
    all: i64,
    by_sample: Vec<i64>,
}

impl Evaluator {
    fn new(name: &str, samples: &[String]) -> Evaluator {
        match name {
            "CompOverlap" => Evaluator::CompOverlap {
                n_eval: 0,
                n_at_comp: 0,
                n_concordant: 0,
                comp_rate: 0.0,
                concordant_rate: 0.0,
                novel_sites: 0,
            },
            "CountVariants" => Evaluator::CountVariants(Box::default()),
            "GenotypeFilterSummary" => Evaluator::GenotypeFilterSummary {
                called: 0,
                no_call_or_filtered: 0,
            },
            "IndelLengthHistogram" => Evaluator::IndelLengthHistogram {
                counts: (-10..=10).filter(|i| *i != 0).map(|i| (i, 0)).collect(),
                n_indels: 0,
                results: Vec::new(),
            },
            "IndelSummary" => Evaluator::IndelSummary(Box::default()),
            "MetricsCollection" => Evaluator::MetricsCollection(Box::default()),
            "MultiallelicSummary" => {
                Evaluator::MultiallelicSummary(Box::new(MultiallelicSummary {
                    snp_novelty: "NA".to_string(),
                    ..Default::default()
                }))
            }
            "PrintMissingComp" => Evaluator::PrintMissingComp { missing: 0 },
            "ThetaVariantEvaluator" => Evaluator::ThetaVariantEvaluator(Box::default()),
            "TiTvVariantEvaluator" => Evaluator::TiTvVariantEvaluator(Box::default()),
            "ValidationReport" => Evaluator::ValidationReport(Box::default()),
            "VariantAFEvaluator" => Evaluator::VariantAFEvaluator(Box::default()),
            "VariantSummary" => {
                let order = hash_map_sample_order(samples);
                let empty = || TypeCounts {
                    all: 0,
                    by_sample: vec![0; order.len()],
                };
                Evaluator::VariantSummary(Box::new(VariantSummary {
                    n_samples: samples.len() as i64,
                    samples_in_hash_order: order.clone(),
                    processed: 0,
                    snps: 0,
                    titv: 0.0,
                    snp_novelty: "NA".to_string(),
                    snps_per_sample: 0,
                    titv_per_sample: 0.0,
                    snp_dp_per_sample: 0.0,
                    indels: 0,
                    indel_novelty: "NA".to_string(),
                    indels_per_sample: 0,
                    indel_dp_per_sample: 0.0,
                    svs: 0,
                    sv_novelty: "NA".to_string(),
                    svs_per_sample: 0,
                    counts_per_sample: [empty(), empty(), empty()],
                    transitions: [empty(), empty(), empty()],
                    transversions: [empty(), empty(), empty()],
                    all_variant_counts: [empty(), empty(), empty()],
                    known_variant_counts: [empty(), empty(), empty()],
                    depth: [empty(), empty(), empty()],
                }))
            }
            other => unreachable!("no evaluator {other}"),
        }
    }

    fn name(&self) -> &'static str {
        match self {
            Evaluator::CompOverlap { .. } => "CompOverlap",
            Evaluator::CountVariants(_) => "CountVariants",
            Evaluator::GenotypeFilterSummary { .. } => "GenotypeFilterSummary",
            Evaluator::IndelLengthHistogram { .. } => "IndelLengthHistogram",
            Evaluator::IndelSummary(_) => "IndelSummary",
            Evaluator::MetricsCollection(_) => "MetricsCollection",
            Evaluator::MultiallelicSummary(_) => "MultiallelicSummary",
            Evaluator::PrintMissingComp { .. } => "PrintMissingComp",
            Evaluator::ThetaVariantEvaluator(_) => "ThetaVariantEvaluator",
            Evaluator::TiTvVariantEvaluator(_) => "TiTvVariantEvaluator",
            Evaluator::ValidationReport(_) => "ValidationReport",
            Evaluator::VariantAFEvaluator(_) => "VariantAFEvaluator",
            Evaluator::VariantSummary(_) => "VariantSummary",
        }
    }

    fn description(name: &str) -> &'static str {
        match name {
            "CompOverlap" => "The overlap between eval and comp sites",
            "CountVariants" => "Counts different classes of variants in the sample",
            "GenotypeFilterSummary" => "Counts called and filtered genotypes across samples",
            "IndelLengthHistogram" => "Indel length histogram",
            "IndelSummary" => "Evaluation summary for indels",
            "MendelianViolationEvaluator" => "Mendelian Violation Evaluator",
            "MetricsCollection" => "Metrics Collection",
            "MultiallelicSummary" => "Evaluation summary for multi-allelic variants",
            "PrintMissingComp" => "count the number of comp SNP sites that are not in eval",
            "ThetaVariantEvaluator" => {
                "Computes different estimates of theta based on variant sites and genotypes"
            }
            "TiTvVariantEvaluator" => "Ti/Tv Variant Evaluator",
            "ValidationReport" => {
                "Assess site accuracy and sensitivity of callset against follow-up validation assay"
            }
            "VariantAFEvaluator" => {
                "Computes different estimates of theta based on variant sites and genotypes"
            }
            "VariantSummary" => "1000 Genomes Phase I summary of variants table",
            _ => "",
        }
    }

    fn comparison_order(&self) -> u8 {
        match self {
            Evaluator::CountVariants(_)
            | Evaluator::GenotypeFilterSummary { .. }
            | Evaluator::IndelLengthHistogram { .. }
            | Evaluator::ThetaVariantEvaluator(_)
            | Evaluator::VariantAFEvaluator(_) => 1,
            _ => 2,
        }
    }

    fn requires_territory(name: &str) -> bool {
        matches!(
            name,
            "CountVariants" | "MultiallelicSummary" | "VariantSummary"
        )
    }

    fn update1(&mut self, vc: &VariantContext, arguments: &Arguments) -> Result<(), EvalError> {
        match self {
            Evaluator::CountVariants(counts) => {
                counts.called += 1;
                if arguments.ignore_ac0() && is_monomorphic_in_samples(vc) {
                    counts.reference += 1;
                } else {
                    match vc_type(vc) {
                        VcType::NoVariation => {}
                        VcType::Snp => {
                            counts.variant += 1;
                            counts.snps += 1;
                            if was_singleton(vc) {
                                counts.singletons += 1;
                            }
                        }
                        VcType::Mnp => {
                            counts.variant += 1;
                            counts.mnps += 1;
                            if was_singleton(vc) {
                                counts.singletons += 1;
                            }
                        }
                        VcType::Indel => {
                            counts.variant += 1;
                            if is_simple_insertion(vc) {
                                counts.insertions += 1;
                            } else if is_simple_deletion(vc) {
                                counts.deletions += 1;
                            } else {
                                counts.complex += 1;
                            }
                        }
                        VcType::Mixed => {
                            counts.variant += 1;
                            counts.mixed += 1;
                        }
                        VcType::Symbolic => counts.symbolic += 1,
                    }
                }
                let ancestral =
                    attribute_string(vc, "ANCESTRALALLELE").map(|text| text.to_uppercase());
                let reference = ancestral
                    .as_ref()
                    .map(|_| vc.alleles[0].display_string().to_uppercase());
                for genotype in vc.genotypes.iter() {
                    let alternate = vc.alleles.get(1).map(|a| a.display_string().to_uppercase());
                    match genotype.genotype_type() {
                        GenotypeType::NoCall => counts.no_calls += 1,
                        GenotypeType::HomRef => {
                            counts.hom_ref += 1;
                            if let (Some(aa), Some(_), Some(reference)) =
                                (&ancestral, &alternate, &reference)
                            {
                                if !reference.eq_ignore_ascii_case(aa) {
                                    counts.hom_derived += 1;
                                }
                            }
                        }
                        GenotypeType::Het => counts.hets += 1,
                        GenotypeType::HomVar => {
                            counts.hom_var += 1;
                            if let (Some(aa), Some(alternate)) = (&ancestral, &alternate) {
                                if !alternate.eq_ignore_ascii_case(aa) {
                                    counts.hom_derived += 1;
                                }
                            }
                        }
                        GenotypeType::Mixed | GenotypeType::Unavailable => {}
                    }
                }
            }
            Evaluator::GenotypeFilterSummary {
                called,
                no_call_or_filtered,
            } => {
                for genotype in vc.genotypes.iter() {
                    let filtered = genotype.is_filtered();
                    if genotype.is_called() && !filtered {
                        *called += 1;
                    } else if genotype.is_no_call() || filtered {
                        *no_call_or_filtered += 1;
                    }
                }
            }
            Evaluator::IndelLengthHistogram {
                counts, n_indels, ..
            } => {
                if is_indel(vc)
                    && !is_complex_indel(vc)
                    && !(arguments.ignore_ac0() && is_monomorphic_in_samples(vc))
                {
                    for allele in &vc.alleles[1..] {
                        let length = allele.len() as i64 - vc.alleles[0].len() as i64;
                        if length == 0 {
                            return Err(EvalError::gatk(format!(
                                "Allele size not expected to be zero for indel: alt = {} ref = {}*",
                                allele.display_string(),
                                vc.alleles[0].display_string()
                            )));
                        }
                        if length.abs() > 10 {
                            continue;
                        }
                        *n_indels += 1;
                        *counts.get_mut(&length).expect("a histogram bin") += 1;
                    }
                }
            }
            Evaluator::ThetaVariantEvaluator(theta) => {
                if !is_snp(vc) || (arguments.ignore_ac0() && is_monomorphic_in_samples(vc)) {
                    return Ok(());
                }
                let mut allele_counts: Vec<(String, i64)> = Vec::new();
                let (mut hets, mut genotyped, mut individuals) = (0i64, 0i64, 0i64);
                for genotype in vc.genotypes.iter() {
                    individuals += 1;
                    if !genotype.is_no_call() {
                        if genotype.is_het() {
                            hets += 1;
                        }
                        genotyped += 1;
                        for allele in &genotype.alleles {
                            if !allele.is_no_call() {
                                let key = format!(
                                    "{}{}",
                                    allele.display_string(),
                                    if allele.is_reference() { "*" } else { "" }
                                );
                                match allele_counts.iter_mut().find(|(k, _)| *k == key) {
                                    Some((_, count)) => *count += 1,
                                    None => allele_counts.push((key, 1)),
                                }
                            }
                        }
                    }
                }
                if genotyped > 0 {
                    theta.num_sites += 1.0;
                    theta.total_het += hets as f64 / genotyped as f64;
                    let mut harmonic: f32 = 0.0;
                    for i in 1..=individuals {
                        harmonic += (1.0 / i as f64) as f32;
                    }
                    theta.theta_region_num_sites += 1.0 / harmonic as f64;
                    let mut pairwise: f32 = 0.0;
                    let mut diffs: i64 = 0;
                    for (first, first_count) in &allele_counts {
                        for (second, second_count) in &allele_counts {
                            match first.cmp(second) {
                                std::cmp::Ordering::Less => continue,
                                std::cmp::Ordering::Equal => {
                                    pairwise +=
                                        (*first_count as f64 * (*first_count - 1) as f64 * 0.5)
                                            as f32;
                                }
                                std::cmp::Ordering::Greater => {
                                    pairwise += (first_count * second_count) as f32;
                                    diffs += first_count * second_count;
                                }
                            }
                        }
                    }
                    if pairwise > 0.0 {
                        theta.total_avg_diffs += (diffs as f32 / pairwise) as f64;
                    }
                }
            }
            Evaluator::VariantAFEvaluator(af) => {
                if !is_snp(vc) || (arguments.ignore_ac0() && is_monomorphic_in_samples(vc)) {
                    return Ok(());
                }
                for genotype in vc.genotypes.iter() {
                    if !genotype.is_no_call() {
                        if genotype.ploidy() != 2 {
                            return Err(EvalError::user(
                                "Bad input: This tool only works with ploidy 2".to_string(),
                            ));
                        }
                        af.called += 1;
                        let references = genotype
                            .alleles
                            .iter()
                            .filter(|allele| **allele == vc.alleles[0])
                            .count() as i64;
                        af.sum += (2.0 - references as f64) / 2.0;
                        if references == 1 {
                            af.het += 1;
                        }
                        if references == 0 {
                            af.hom_var += 1;
                        }
                        if references == 2 {
                            af.hom_ref += 1;
                        }
                    }
                }
                if vc.genotypes.is_empty() {
                    af.called += 1;
                    af.sum += attribute_double(vc, "AF");
                }
            }
            _ => {}
        }
        Ok(())
    }

    fn update2(
        &mut self,
        eval: Option<&VariantContext>,
        comp: Option<&VariantContext>,
        arguments: &Arguments,
        samples_for_evaluation: &[String],
    ) -> Result<(), EvalError> {
        match self {
            Evaluator::CompOverlap {
                n_eval,
                n_at_comp,
                n_concordant,
                ..
            } => {
                let eval_good = eval.is_some_and(is_polymorphic_in_samples);
                let comp_good = comp.is_some_and(|comp| !is_filtered(comp));
                if eval_good {
                    *n_eval += 1;
                }
                if comp_good && eval_good {
                    *n_at_comp += 1;
                    let (eval, comp) = (eval.expect("eval"), comp.expect("comp"));
                    let discordant = eval.alleles.iter().any(|allele| {
                        !comp.alleles.iter().any(|candidate| {
                            candidate.display_string() == allele.display_string()
                                && candidate.is_reference() == allele.is_reference()
                        })
                    });
                    if !discordant {
                        *n_concordant += 1;
                    }
                }
            }
            Evaluator::IndelSummary(summary) => {
                let Some(eval) = eval else { return Ok(()) };
                if arguments.ignore_ac0() && is_monomorphic_in_samples(eval) {
                    return Ok(());
                }
                match vc_type(eval) {
                    VcType::Snp => {
                        summary.snps += eval.alleles.len() as i64 - 1;
                        if was_singleton(eval) {
                            summary.singleton_snps += 1;
                        }
                        for genotype in eval.genotypes.iter() {
                            if genotype.is_het() {
                                summary.snp_hets += 1;
                            }
                            if genotype.is_hom_var() {
                                summary.snp_homs += 1;
                            }
                        }
                    }
                    VcType::Indel => {
                        let gold = arguments.gold_standard.as_ref().is_some_and(|gold| {
                            gold.iter().any(|record| {
                                record.contig == eval.contig
                                    && record.start <= eval.start
                                    && eval.start <= record.stop
                            })
                        });
                        summary.indel_sites += 1;
                        if !is_biallelic(eval) {
                            summary.multiallelic_indel_sites += 1;
                        }
                        for genotype in eval.genotypes.iter() {
                            if genotype.is_het() {
                                summary.indel_hets += 1;
                            }
                            if genotype.is_hom_var() {
                                summary.indel_homs += 1;
                            }
                        }
                        for allele in &eval.alleles[1..] {
                            summary.indels += 1;
                            if was_singleton(eval) {
                                summary.singleton_indels += 1;
                            }
                            if comp.is_none() {
                                summary.novel_indels += 1;
                            }
                            if gold {
                                summary.matching_gold += 1;
                            }
                            let size = allele.len() as i64 - eval.alleles[0].len() as i64;
                            if size == 0 {
                                return Err(EvalError::gatk(format!(
                                    "Allele size not expected to be zero for indel: alt = {} ref = {}*",
                                    allele.display_string(),
                                    eval.alleles[0].display_string()
                                )));
                            }
                            if size > 0 {
                                summary.insertions += 1;
                            }
                            if size < 0 {
                                summary.deletions += 1;
                            }
                            let biotype = attribute_string(eval, "SNPEFF_GENE_BIOTYPE")
                                .unwrap_or_else(|| "missing".to_string());
                            if biotype == "protein_coding" {
                                let effect = attribute_string(eval, "SNPEFF_EFFECT")
                                    .unwrap_or_else(|| "missing".to_string());
                                if effect == "missing" {
                                    return Err(EvalError::gatk(format!(
                                        "Saw SNPEFF_GENE_BIOTYPE but unexpected no SNPEFF_EFFECT at \
                                         {}:{}",
                                        eval.contig, eval.start
                                    )));
                                }
                                if effect == "FRAME_SHIFT" {
                                    summary.frameshifting += 1;
                                } else if effect.starts_with("CODON") {
                                    summary.in_frame += 1;
                                }
                            }
                            if size > 10 {
                                summary.large_insertions += 1;
                            } else if size < -10 {
                                summary.large_deletions += 1;
                            }
                            let absolute = size.unsigned_abs() as usize;
                            let table = if size < 0 {
                                &mut summary.deletion_by_length
                            } else {
                                &mut summary.insertion_by_length
                            };
                            if absolute < table.len() {
                                table[absolute] += 1;
                            }
                        }
                    }
                    _ => {}
                }
            }
            Evaluator::MultiallelicSummary(summary) => {
                let Some(eval) = eval else { return Ok(()) };
                if arguments.ignore_ac0() && is_monomorphic_in_samples(eval) {
                    return Ok(());
                }
                match vc_type(eval) {
                    VcType::Snp => {
                        summary.snps += 1;
                        if !is_biallelic(eval) {
                            summary.multi_snps += 1;
                            let reference = bases(&eval.alleles[0]);
                            for allele in &eval.alleles[1..] {
                                if is_transition_bases(reference[0], bases(allele)[0]) {
                                    summary.ti += 1;
                                } else {
                                    summary.tv += 1;
                                }
                            }
                            if let Some(comp) = comp {
                                let known = eval.alleles[1..]
                                    .iter()
                                    .filter(|allele| comp.alleles[1..].contains(allele))
                                    .count();
                                if known == eval.alleles.len() - 1 {
                                    summary.known_complete += 1;
                                } else if known > 0 {
                                    summary.known_partial += 1;
                                }
                            }
                        }
                    }
                    VcType::Indel => {
                        summary.indels += 1;
                        if !is_biallelic(eval) {
                            summary.multi_indels += 1;
                        }
                    }
                    _ => {}
                }
            }
            Evaluator::PrintMissingComp { missing } => {
                let comp_good = comp.is_some_and(|comp| !is_filtered(comp) && is_snp(comp));
                let eval_good = eval.is_some_and(is_snp);
                if comp_good && !eval_good {
                    *missing += 1;
                }
            }
            Evaluator::TiTvVariantEvaluator(titv) => {
                if let Some(eval) = eval {
                    titv.update(eval, false);
                }
                if let Some(comp) = comp {
                    titv.update(comp, true);
                }
            }
            Evaluator::ValidationReport(report) => {
                let Some(comp) = comp else { return Ok(()) };
                let eval_status = site_status(eval);
                let do_subset = !comp.genotypes.is_empty()
                    && !samples_for_evaluation.is_empty()
                    && samples_for_evaluation
                        .iter()
                        .all(|sample| comp.genotype(sample).is_some());
                let comp_status = if do_subset {
                    site_status(Some(&sub_context_from_samples(
                        comp,
                        samples_for_evaluation,
                        false,
                    )))
                } else {
                    site_status(Some(comp))
                };
                report.counts[comp_status][eval_status] += 1;
            }
            Evaluator::VariantSummary(summary) => {
                let Some(eval) = eval else { return Ok(()) };
                if arguments.ignore_ac0() && is_monomorphic_in_samples(eval) {
                    return Ok(());
                }
                let kind = match vc_type(eval) {
                    VcType::Snp => 0usize,
                    VcType::Indel => {
                        if indel_lengths(eval)
                            .unwrap_or_default()
                            .iter()
                            .any(|length| length.abs() > 50)
                        {
                            2
                        } else {
                            1
                        }
                    }
                    VcType::Symbolic => 2,
                    _ => return Ok(()),
                };
                if attribute(eval, "DP").is_some() {
                    summary.depth[kind].all += 1;
                }
                summary.all_variant_counts[kind].all += 1;
                let mut titv_table: Option<bool> = None;
                if kind == 0 && is_biallelic(eval) {
                    let transition = is_transition(eval);
                    titv_table = Some(transition);
                    if transition {
                        summary.transitions[kind].all += 1;
                    } else {
                        summary.transversions[kind].all += 1;
                    }
                }
                // `comp != null || (type == CNV && overlapsKnownCNV(eval, context))`: the known
                // CNVs are queried only for a CNV with no comp, and an unindexed file is refused
                // at that first query.
                let known = comp.is_some()
                    || (kind == 2
                        && match &arguments.known_cnvs {
                            None => false,
                            Some(known) => {
                                if let Some(path) = &arguments.known_cnvs_unindexed {
                                    return Err(EvalError::user(format!(
                                        "Input {path} must support random access to enable \
                                         queries by interval. If it's a file, please index it \
                                         using the bundled tool IndexFeatureFile"
                                    )));
                                }
                                overlapping(known, &eval.contig, eval.start, eval.stop)
                                    .iter()
                                    .any(|feature| reciprocal_overlap(eval, feature) > 0.5)
                            }
                        });
                if known {
                    summary.known_variant_counts[kind].all += 1;
                }
                for genotype in eval.genotypes.iter() {
                    if !genotype.is_no_call() && !genotype.is_hom_ref() {
                        let Some(index) = summary
                            .samples_in_hash_order
                            .iter()
                            .position(|sample| *sample == genotype.sample_name)
                        else {
                            return Err(EvalError {
                                class: "java.lang.NullPointerException".to_string(),
                                message: "null".to_string(),
                                user: false,
                            });
                        };
                        summary.counts_per_sample[kind].by_sample[index] += 1;
                        match titv_table {
                            Some(true) => summary.transitions[kind].by_sample[index] += 1,
                            Some(false) => summary.transversions[kind].by_sample[index] += 1,
                            None => {}
                        }
                        if genotype.dp.is_some() {
                            summary.depth[kind].by_sample[index] += 1;
                        }
                    }
                }
            }
            _ => {}
        }
        Ok(())
    }

    fn finalize(&mut self, processed_loci: Option<i64>) {
        match self {
            Evaluator::CompOverlap {
                n_eval,
                n_at_comp,
                n_concordant,
                comp_rate,
                concordant_rate,
                novel_sites,
            } => {
                *comp_rate = 100.0 * rate(*n_at_comp, *n_eval);
                *concordant_rate = 100.0 * rate(*n_concordant, *n_at_comp);
                *novel_sites = *n_eval - *n_at_comp;
            }
            Evaluator::CountVariants(c) => {
                c.processed = processed_loci.unwrap_or(0);
                c.variant_rate = rate(c.variant, c.processed);
                c.variant_rate_per_bp = inverse_rate(c.variant, c.processed) as f64;
                c.heterozygosity = rate(c.hets, c.processed);
                c.heterozygosity_per_bp = inverse_rate(c.hets, c.processed) as f64;
                c.het_hom_ratio = ratio(c.hets, c.hom_var);
                c.indel_rate = rate(c.deletions + c.insertions + c.complex, c.processed);
                c.indel_rate_per_bp =
                    inverse_rate(c.deletions + c.insertions + c.complex, c.processed) as f64;
                c.insertion_deletion_ratio = ratio(c.insertions, c.deletions);
            }
            Evaluator::IndelLengthHistogram {
                counts,
                n_indels,
                results,
            } => {
                *results = counts
                    .iter()
                    .map(|(length, count)| {
                        let value = if *n_indels == 0 {
                            0.0
                        } else {
                            *count as f64 / (1.0 * *n_indels as f64)
                        };
                        (*length, value)
                    })
                    .collect();
            }
            Evaluator::IndelSummary(sm) => {
                sm.percent_multiallelic =
                    formatted_percent(sm.multiallelic_indel_sites, sm.indel_sites);
                sm.snp_to_indel = formatted_ratio(sm.snps, sm.indels);
                sm.snp_to_indel_singletons =
                    formatted_ratio(sm.singleton_snps, sm.singleton_indels);
                sm.gold_rate = formatted_percent(sm.matching_gold, sm.indels);
                sm.novelty_rate = formatted_percent(sm.novel_indels, sm.indels);
                sm.frameshift_rate =
                    formatted_percent(sm.frameshifting, sm.in_frame + sm.frameshifting);
                sm.ratio_deletions = formatted_ratio(
                    sm.deletion_by_length[1] + sm.deletion_by_length[2],
                    sm.deletion_by_length[3],
                );
                sm.ratio_insertions = formatted_ratio(
                    sm.insertion_by_length[1] + sm.insertion_by_length[2],
                    sm.insertion_by_length[3],
                );
                sm.snp_het_to_hom = formatted_ratio(sm.snp_hets, sm.snp_homs);
                sm.indel_het_to_hom = formatted_ratio(sm.indel_hets, sm.indel_homs);
                sm.insertion_to_deletion = formatted_ratio(sm.insertions, sm.deletions);
                sm.insertion_to_deletion_large =
                    formatted_ratio(sm.large_insertions, sm.large_deletions);
            }
            Evaluator::MultiallelicSummary(sm) => {
                sm.processed = processed_loci.unwrap_or(0);
                sm.processed_multi_snp_ratio = sm.multi_snps as f64 / sm.processed as f64;
                sm.variant_multi_snp_ratio = sm.multi_snps as f64 / sm.snps as f64;
                sm.processed_multi_indel_ratio = sm.multi_indels as f64 / sm.processed as f64;
                sm.variant_multi_indel_ratio = sm.multi_indels as f64 / sm.indels as f64;
                sm.titv = sm.ti as f64 / sm.tv as f64;
                let known = sm.known_partial + sm.known_complete;
                sm.snp_novelty = formatted_percent(sm.multi_snps - known, sm.multi_snps);
            }
            Evaluator::ThetaVariantEvaluator(theta) => {
                if theta.num_sites > 0.0 {
                    theta.avg_het = theta.total_het / theta.num_sites;
                    theta.avg_avg_diffs = theta.total_avg_diffs / theta.num_sites;
                }
            }
            Evaluator::TiTvVariantEvaluator(titv) => {
                titv.ratio = rate(titv.ti, titv.tv);
                titv.ratio_derived = rate(titv.ti_derived, titv.tv_derived);
                titv.ratio_standard = rate(titv.ti_comp, titv.tv_comp);
            }
            Evaluator::ValidationReport(r) => {
                // SiteStatus: NO_CALL, FILTERED, MONO, POLY.
                for x in 0..4 {
                    r.comp_filtered += r.counts[1][x];
                }
                let get = |comp: usize, eval: usize| r.counts[comp][eval];
                let (mono_nc, mono_f, mono_m, mono_p) =
                    (get(2, 0), get(2, 1), get(2, 2), get(2, 3));
                let (poly_nc, poly_f, poly_m, poly_p) =
                    (get(3, 0), get(3, 1), get(3, 2), get(3, 3));
                r.tp = poly_p;
                r.fn_ = poly_nc + poly_f + poly_m;
                r.fp = mono_p;
                r.tn = mono_nc + mono_f + mono_m;
                r.n_comp = r.counts.iter().flatten().sum();
                r.sensitivity = (100.0 * r.tp as f64) / (r.tp + r.fn_) as f64;
                r.specificity = if r.tn + r.fp > 0 {
                    (100.0 * r.tn as f64) / (r.tn + r.fp) as f64
                } else {
                    100.0
                };
                r.ppv = (100.0 * r.tp as f64) / (r.tp + r.fp) as f64;
                r.fdr = (100.0 * r.fp as f64) / (r.fp + r.tp) as f64;
            }
            Evaluator::VariantAFEvaluator(af) => {
                af.avg = if af.called == 0 {
                    0.0
                } else {
                    af.sum / af.called as f64
                };
            }
            Evaluator::VariantSummary(sm) => {
                sm.processed = processed_loci.unwrap_or(0);
                sm.snps = sm.all_variant_counts[0].all;
                sm.indels = sm.all_variant_counts[1].all;
                sm.svs = sm.all_variant_counts[2].all;
                sm.titv = ratio(sm.transitions[0].all, sm.transversions[0].all);
                sm.titv_per_sample = {
                    let mut sum = 0.0;
                    let n = sm.samples_in_hash_order.len();
                    for index in 0..n {
                        sum += ratio(
                            sm.transitions[0].by_sample[index],
                            sm.transversions[0].by_sample[index],
                        );
                    }
                    if n > 0 {
                        sum / n as f64
                    } else {
                        0.0
                    }
                };
                let mean = |counts: &TypeCounts| -> i64 {
                    let sum: i64 = counts.by_sample.iter().sum();
                    crate::reference_confidence_merger::java_round(
                        sum as f64 / (1.0 * counts.by_sample.len() as f64),
                    )
                };
                sm.snps_per_sample = mean(&sm.counts_per_sample[0]);
                sm.indels_per_sample = mean(&sm.counts_per_sample[1]);
                sm.svs_per_sample = mean(&sm.counts_per_sample[2]);
                let novelty = |all: i64, known: i64| formatted_percent(all - known, all);
                sm.snp_novelty =
                    novelty(sm.all_variant_counts[0].all, sm.known_variant_counts[0].all);
                sm.indel_novelty =
                    novelty(sm.all_variant_counts[1].all, sm.known_variant_counts[1].all);
                sm.sv_novelty =
                    novelty(sm.all_variant_counts[2].all, sm.known_variant_counts[2].all);
                sm.snp_dp_per_sample = mean(&sm.depth[0]) as f64;
                sm.indel_dp_per_sample = mean(&sm.depth[1]) as f64;
            }
            _ => {}
        }
    }

    /// The `@DataPoint` fields in declaration order, or the molten rows.
    fn cells(&self) -> Result<Vec<Cell>, Vec<(i64, f64)>> {
        use ReportValue::{Double as D, Int as I, Str as S};
        let text = |value: &String| S(value.clone());
        Ok(match self {
            Evaluator::CompOverlap {
                n_eval,
                n_at_comp,
                n_concordant,
                comp_rate,
                concordant_rate,
                novel_sites,
            } => vec![
                ("nEvalVariants", "%d", I(*n_eval)),
                ("novelSites", "%d", I(*novel_sites)),
                ("nVariantsAtComp", "%d", I(*n_at_comp)),
                ("compRate", "%.2f", D(*comp_rate)),
                ("nConcordant", "%d", I(*n_concordant)),
                ("concordantRate", "%.2f", D(*concordant_rate)),
            ],
            Evaluator::CountVariants(c) => vec![
                ("nProcessedLoci", "%d", I(c.processed)),
                ("nCalledLoci", "%d", I(c.called)),
                ("nRefLoci", "%d", I(c.reference)),
                ("nVariantLoci", "%d", I(c.variant)),
                ("variantRate", "%.8f", D(c.variant_rate)),
                ("variantRatePerBp", "%.8f", D(c.variant_rate_per_bp)),
                ("nSNPs", "%d", I(c.snps)),
                ("nMNPs", "%d", I(c.mnps)),
                ("nInsertions", "%d", I(c.insertions)),
                ("nDeletions", "%d", I(c.deletions)),
                ("nComplex", "%d", I(c.complex)),
                ("nSymbolic", "%d", I(c.symbolic)),
                ("nMixed", "%d", I(c.mixed)),
                ("nNoCalls", "%d", I(c.no_calls)),
                ("nHets", "%d", I(c.hets)),
                ("nHomRef", "%d", I(c.hom_ref)),
                ("nHomVar", "%d", I(c.hom_var)),
                ("nSingletons", "%d", I(c.singletons)),
                ("nHomDerived", "%d", I(c.hom_derived)),
                ("heterozygosity", "%.2e", D(c.heterozygosity)),
                ("heterozygosityPerBp", "%.2f", D(c.heterozygosity_per_bp)),
                ("hetHomRatio", "%.2f", D(c.het_hom_ratio)),
                ("indelRate", "%.2e", D(c.indel_rate)),
                ("indelRatePerBp", "%.2f", D(c.indel_rate_per_bp)),
                (
                    "insertionDeletionRatio",
                    "%.2f",
                    D(c.insertion_deletion_ratio),
                ),
            ],
            Evaluator::GenotypeFilterSummary {
                called,
                no_call_or_filtered,
            } => vec![
                ("nCalledNotFiltered", "%d", I(*called)),
                ("nNoCallOrFiltered", "%d", I(*no_call_or_filtered)),
            ],
            Evaluator::IndelLengthHistogram { results, .. } => return Err(results.clone()),
            Evaluator::IndelSummary(sm) => vec![
                ("n_SNPs", "%d", I(sm.snps)),
                ("n_singleton_SNPs", "%d", I(sm.singleton_snps)),
                ("n_indels", "%d", I(sm.indels)),
                ("n_singleton_indels", "%d", I(sm.singleton_indels)),
                ("n_indels_matching_gold_standard", "%d", I(sm.matching_gold)),
                ("gold_standard_matching_rate", "", text(&sm.gold_rate)),
                (
                    "n_multiallelic_indel_sites",
                    "",
                    I(sm.multiallelic_indel_sites),
                ),
                (
                    "percent_of_sites_with_more_than_2_alleles",
                    "",
                    text(&sm.percent_multiallelic),
                ),
                ("SNP_to_indel_ratio", "", text(&sm.snp_to_indel)),
                (
                    "SNP_to_indel_ratio_for_singletons",
                    "",
                    text(&sm.snp_to_indel_singletons),
                ),
                ("n_novel_indels", "%d", I(sm.novel_indels)),
                ("indel_novelty_rate", "", text(&sm.novelty_rate)),
                ("n_insertions", "", I(sm.insertions)),
                ("n_deletions", "", I(sm.deletions)),
                (
                    "insertion_to_deletion_ratio",
                    "",
                    text(&sm.insertion_to_deletion),
                ),
                ("n_large_deletions", "", I(sm.large_deletions)),
                ("n_large_insertions", "", I(sm.large_insertions)),
                (
                    "insertion_to_deletion_ratio_for_large_indels",
                    "",
                    text(&sm.insertion_to_deletion_large),
                ),
                ("n_coding_indels_frameshifting", "", I(sm.frameshifting)),
                ("n_coding_indels_in_frame", "", I(sm.in_frame)),
                (
                    "frameshift_rate_for_coding_indels",
                    "",
                    text(&sm.frameshift_rate),
                ),
                ("SNP_het_to_hom_ratio", "", text(&sm.snp_het_to_hom)),
                ("indel_het_to_hom_ratio", "", text(&sm.indel_het_to_hom)),
                (
                    "ratio_of_1_and_2_to_3_bp_insertions",
                    "",
                    text(&sm.ratio_insertions),
                ),
                (
                    "ratio_of_1_and_2_to_3_bp_deletions",
                    "",
                    text(&sm.ratio_deletions),
                ),
            ],
            Evaluator::MetricsCollection(m) => vec![
                ("concordantRate", "%.2f", D(m.concordant_rate)),
                ("nSNPs", "%d", I(m.snps)),
                ("nSNPloci", "%d", I(m.snp_loci)),
                ("nIndels", "%d", I(m.indels)),
                ("nIndelLoci", "%d", I(m.indel_loci)),
                (
                    "indelRatio",
                    "",
                    match &m.indel_ratio {
                        Some(value) => S(value.clone()),
                        None => ReportValue::Null,
                    },
                ),
                ("indelRatioLociBased", "%.2f", D(m.indel_ratio_loci)),
                ("tiTvRatio", "%.2f", D(m.titv)),
            ],
            Evaluator::MultiallelicSummary(sm) => vec![
                ("nProcessedLoci", "%d", I(sm.processed)),
                ("nSNPs", "%d", I(sm.snps)),
                ("nMultiSNPs", "%d", I(sm.multi_snps)),
                (
                    "processedMultiSnpRatio",
                    "%.5f",
                    D(sm.processed_multi_snp_ratio),
                ),
                (
                    "variantMultiSnpRatio",
                    "%.3f",
                    D(sm.variant_multi_snp_ratio),
                ),
                ("nIndels", "%d", I(sm.indels)),
                ("nMultiIndels", "%d", I(sm.multi_indels)),
                (
                    "processedMultiIndelRatio",
                    "%.5f",
                    D(sm.processed_multi_indel_ratio),
                ),
                (
                    "variantMultiIndelRatio",
                    "%.3f",
                    D(sm.variant_multi_indel_ratio),
                ),
                ("nTi", "%d", I(sm.ti)),
                ("nTv", "%d", I(sm.tv)),
                ("TiTvRatio", "%.2f", D(sm.titv)),
                ("knownSNPsPartial", "%d", I(sm.known_partial)),
                ("knownSNPsComplete", "%d", I(sm.known_complete)),
                ("SNPNoveltyRate", "", text(&sm.snp_novelty)),
            ],
            Evaluator::PrintMissingComp { missing } => vec![("nMissing", "%d", I(*missing))],
            Evaluator::ThetaVariantEvaluator(t) => vec![
                ("avgHet", "%.8f", D(t.avg_het)),
                ("avgAvgDiffs", "%.8f", D(t.avg_avg_diffs)),
                ("totalHet", "%.8f", D(t.total_het)),
                ("totalAvgDiffs", "%.8f", D(t.total_avg_diffs)),
                ("thetaRegionNumSites", "%.8f", D(t.theta_region_num_sites)),
            ],
            Evaluator::TiTvVariantEvaluator(t) => vec![
                ("nTi", "%d", I(t.ti)),
                ("nTv", "%d", I(t.tv)),
                ("tiTvRatio", "%.2f", D(t.ratio)),
                ("nTiInComp", "%d", I(t.ti_comp)),
                ("nTvInComp", "%d", I(t.tv_comp)),
                ("TiTvRatioStandard", "%.2f", D(t.ratio_standard)),
                ("nTiDerived", "%d", I(t.ti_derived)),
                ("nTvDerived", "%d", I(t.tv_derived)),
                ("tiTvDerivedRatio", "%.2f", D(t.ratio_derived)),
            ],
            Evaluator::ValidationReport(r) => {
                let get = |comp: usize, eval: usize| I(r.counts[comp][eval]);
                vec![
                    ("nComp", "%d", I(r.n_comp)),
                    ("TP", "%d", I(r.tp)),
                    ("FP", "%d", I(r.fp)),
                    ("FN", "%d", I(r.fn_)),
                    ("TN", "%d", I(r.tn)),
                    ("sensitivity", "%.2f", D(r.sensitivity)),
                    ("specificity", "%.2f", D(r.specificity)),
                    ("PPV", "%.2f", D(r.ppv)),
                    ("FDR", "%.2f", D(r.fdr)),
                    ("CompMonoEvalNoCall", "%d", get(2, 0)),
                    ("CompMonoEvalFiltered", "%d", get(2, 1)),
                    ("CompMonoEvalMono", "%d", get(2, 2)),
                    ("CompMonoEvalPoly", "%d", get(2, 3)),
                    ("CompPolyEvalNoCall", "%d", get(3, 0)),
                    ("CompPolyEvalFiltered", "%d", get(3, 1)),
                    ("CompPolyEvalMono", "%d", get(3, 2)),
                    ("CompPolyEvalPoly", "%d", get(3, 3)),
                    ("CompFiltered", "%d", I(r.comp_filtered)),
                    ("nDifferentAlleleSites", "%d", I(r.different_alleles)),
                ]
            }
            Evaluator::VariantAFEvaluator(af) => vec![
                ("avgVarAF", "%.8f", D(af.avg)),
                ("totalCalledSites", "%d", I(af.called)),
                ("totalHetSites", "%d", I(af.het)),
                ("totalHomVarSites", "%d", I(af.hom_var)),
                ("totalHomRefSites", "%d", I(af.hom_ref)),
            ],
            Evaluator::VariantSummary(sm) => vec![
                ("nSamples", "%d", I(sm.n_samples)),
                ("nProcessedLoci", "%d", I(sm.processed)),
                ("nSNPs", "%d", I(sm.snps)),
                ("TiTvRatio", "%.2f", D(sm.titv)),
                ("SNPNoveltyRate", "%s", text(&sm.snp_novelty)),
                ("nSNPsPerSample", "%d", I(sm.snps_per_sample)),
                ("TiTvRatioPerSample", "%.2f", D(sm.titv_per_sample)),
                ("SNPDPPerSample", "%.1f", D(sm.snp_dp_per_sample)),
                ("nIndels", "%d", I(sm.indels)),
                ("IndelNoveltyRate", "%s", text(&sm.indel_novelty)),
                ("nIndelsPerSample", "%d", I(sm.indels_per_sample)),
                ("IndelDPPerSample", "%.1f", D(sm.indel_dp_per_sample)),
                ("nSVs", "%d", I(sm.svs)),
                ("SVNoveltyRate", "%s", text(&sm.sv_novelty)),
                ("nSVsPerSample", "%d", I(sm.svs_per_sample)),
            ],
        })
    }
}

impl TiTv {
    fn update(&mut self, vc: &VariantContext, standard: bool) {
        if !(is_snp(vc) && is_biallelic(vc) && is_polymorphic_in_samples(vc)) {
            return;
        }
        let transition = is_transition(vc);
        match (transition, standard) {
            (true, true) => self.ti_comp += 1,
            (true, false) => self.ti += 1,
            (false, true) => self.tv_comp += 1,
            (false, false) => self.tv += 1,
        }
        if let Some(aa) = attribute_string(vc, "ANCESTRALALLELE") {
            let aa = aa.to_uppercase();
            if aa != "." {
                let alternate = bases(&vc.alleles[1]);
                let first = aa.as_bytes().first().copied().unwrap_or(0);
                let (a, b) = (base_index(first), base_index(alternate[0]));
                if a >= 0 && b >= 0 {
                    if is_transition_bases(first, alternate[0]) {
                        self.ti_derived += 1;
                    } else {
                        self.tv_derived += 1;
                    }
                }
            }
        }
    }
}

/// `ValidationReport.calcSiteStatus`: NO_CALL, FILTERED, MONO or POLY, as ordinals.
fn site_status(vc: Option<&VariantContext>) -> usize {
    let Some(vc) = vc else { return 0 };
    if is_filtered(vc) {
        return 1;
    }
    if is_monomorphic_in_samples(vc) {
        return 2;
    }
    if !vc.genotypes.is_empty() {
        return 3;
    }
    if attribute(vc, "AC").is_some() {
        if vc.alleles.len() > 2 {
            return 3;
        }
        return if attribute_int(vc, "AC").unwrap_or(0) > 0 {
            3
        } else {
            2
        };
    }
    3
}

/// `GenomeLoc.reciprocialOverlapFraction`.
fn reciprocal_overlap(vc: &VariantContext, feature: &Feature) -> f64 {
    let start = vc.start.max(feature.start);
    let stop = vc.stop.min(feature.end);
    if start > stop {
        return 0.0;
    }
    let overlap = (stop - start + 1) as f64;
    let first = overlap / (vc.stop - vc.start + 1) as f64;
    let second = overlap / (feature.end - feature.start + 1) as f64;
    first.min(second)
}

/// The order a `HashMap<String, Integer>` holding these samples and `ALL` iterates the samples in.
fn hash_map_sample_order(samples: &[String]) -> Vec<String> {
    let mut keys: Vec<String> = samples.to_vec();
    keys.push("ALL".to_string());
    let ordered = gatk_engine::java_hash::hash_set_order(&keys).unwrap_or(keys);
    ordered.into_iter().filter(|key| key != "ALL").collect()
}

// ---------------------------------------------------------------------------------------------
// The engine
// ---------------------------------------------------------------------------------------------

/// One driving record, with the input it came from: an index into `evals`, then into `comps`.
#[derive(Debug, Clone)]
pub struct Driving {
    pub source: Source,
    pub record: VariantContext,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Source {
    Eval(usize),
    Comp(usize),
}

pub struct Engine {
    pub arguments: Arguments,
    pub strats: Vec<Stratifier>,
    pub evaluator_names: Vec<String>,
    /// One context per combination of states, in mixed-radix order over `strats`.
    contexts: Vec<Vec<Evaluator>>,
    by_filter: bool,
    per_sample: bool,
    sample_names_for_stratification: Vec<String>,
}

impl Engine {
    /// `validateAndInitialize`.
    pub fn new(
        arguments: Arguments,
        contigs_of_dictionary: &[String],
    ) -> Result<Engine, EvalError> {
        // Stratifiers: the required ones, the standard ones unless refused, and the named ones.
        let mut names: Vec<String> = STRATIFIER_NAMES
            .iter()
            .filter(|(_, standard, required)| {
                *required || (*standard && !arguments.no_standard_strats)
            })
            .map(|(name, _, _)| name.to_string())
            .collect();
        for module in &arguments.strats_to_use {
            if !names.contains(module) {
                names.push(module.clone());
            }
        }
        // `stratsToUse` is a HashSet, so an unknown name is found in whatever order it iterates;
        // with one unknown the message is the same.
        for module in &names {
            if !STRATIFIER_NAMES.iter().any(|(name, _, _)| name == module) {
                return Err(EvalError::command_line(format!(
                    "Module {module} could not be found; please check that you have specified the \
                     class name correctly"
                )));
            }
        }
        names.sort();

        let mut evaluator_names: Vec<String> = arguments.modules_to_use.clone();
        if !arguments.no_standard_modules {
            for (name, standard) in EVALUATOR_NAMES {
                if *standard {
                    evaluator_names.push(name.to_string());
                }
            }
        }
        evaluator_names.sort();
        evaluator_names.dedup();
        for module in &evaluator_names {
            if !EVALUATOR_NAMES.iter().any(|(name, _)| name == module) {
                return Err(EvalError::command_line(format!(
                    "Module {module} could not be found; please check that you have specified the \
                     class name correctly"
                )));
            }
        }
        if [
            "CompOverlap",
            "IndelSummary",
            "TiTvVariantEvaluator",
            "CountVariants",
            "MultiallelicSummary",
        ]
        .iter()
        .all(|needed| evaluator_names.iter().any(|name| name == needed))
            && !evaluator_names
                .iter()
                .any(|name| name == "MetricsCollection")
        {
            evaluator_names.push("MetricsCollection".to_string());
            evaluator_names.sort();
        }
        if evaluator_names
            .iter()
            .any(|name| name == "MendelianViolationEvaluator")
        {
            return Err(EvalError::gatk(
                "MendelianViolationEvaluator needs a pedigree, which this port does not carry for \
                 VariantEval yet."
                    .to_string(),
            ));
        }

        let mut sample_names_for_stratification = Vec::new();
        if arguments.strats_to_use.iter().any(|name| name == "Sample") {
            sample_names_for_stratification
                .extend(arguments.samples_for_evaluation.iter().cloned());
        }
        sample_names_for_stratification.push(ALL.to_string());

        let mut strats = Vec::new();
        for name in &names {
            strats.push(make_stratifier(
                name,
                &arguments,
                &sample_names_for_stratification,
                contigs_of_dictionary,
            )?);
        }
        // `checkForIncompatibleEvaluatorsAndStratifiers`.
        for strat in &strats {
            if matches!(strat.name, "AlleleCount" | "Family" | "Sample")
                && evaluator_names.iter().any(|name| name == "VariantSummary")
            {
                return Err(EvalError::bad_argument(
                    "ST and ET",
                    &format!(
                        "The selected stratification {} and evaluator VariantSummary are \
                         incompatible due to combinatorial memory requirements. Please disable one",
                        strat.name
                    ),
                ));
            }
        }
        let by_filter = strats.iter().any(|strat| strat.name == "Filter");
        let per_sample = strats.iter().any(|strat| strat.name == "Sample");
        let per_family = strats.iter().any(|strat| strat.name == "Family");
        if per_sample && per_family {
            return Err(EvalError::bad_argument(
                "ST",
                "Variants cannot be stratified by sample and family at the same time",
            ));
        }
        if per_family {
            return Err(EvalError::bad_argument(
                "ST",
                "Cannot stratify by family without *.ped file",
            ));
        }
        if arguments.strat_intervals.is_some()
            && !strats
                .iter()
                .any(|strat| strat.name == "IntervalStratification")
        {
            return Err(EvalError::bad_argument(
                "ST",
                "stratIntervals argument provided but -ST IntervalStratification not provided",
            ));
        }
        let size: usize = strats
            .iter()
            .map(|strat| distinct_states(&strat.states).len())
            .product();
        let samples = arguments.samples_for_evaluation.clone();
        let contexts = (0..size)
            .map(|_| {
                evaluator_names
                    .iter()
                    .map(|name| Evaluator::new(name, &samples))
                    .collect()
            })
            .collect();
        Ok(Engine {
            arguments,
            strats,
            evaluator_names,
            contexts,
            by_filter,
            per_sample,
            sample_names_for_stratification,
        })
    }

    fn key_of(&self, indices: &[usize]) -> usize {
        let mut key = 0;
        for (strat, index) in self.strats.iter().zip(indices) {
            key = key * distinct_states(&strat.states).len() + index;
        }
        key
    }

    fn states_of_key(&self, mut key: usize) -> Vec<State> {
        let mut out = vec![State::Int(0); self.strats.len()];
        for (position, strat) in self.strats.iter().enumerate().rev() {
            let states = distinct_states(&strat.states);
            out[position] = states[key % states.len()].clone();
            key /= states.len();
        }
        out
    }

    /// `bindVariantContexts` for the evals: per track, per sample, the subset records.
    fn bind_evals<'a>(
        &self,
        by_input: &'a [Vec<&'a VariantContext>],
    ) -> Vec<Vec<(String, Vec<VariantContext>)>> {
        let mut bindings: Vec<Vec<(String, Vec<VariantContext>)>> = Vec::new();
        for (track, records) in by_input.iter().enumerate() {
            let mut mapping: Vec<(String, Vec<VariantContext>)> = Vec::new();
            let mut add = |sample: &str, vc: VariantContext| match mapping
                .iter_mut()
                .find(|(name, _)| name == sample)
            {
                Some((_, list)) => list.push(vc),
                None => mapping.push((sample.to_string(), vec![vc])),
            };
            for vc in records {
                let mut sub = (*vc).clone();
                if !vc.genotypes.is_empty() {
                    sub = self.subset(vc, &self.arguments.samples_for_evaluation);
                }
                if self.by_filter || !is_filtered(&sub) {
                    add(ALL, sub);
                }
                if !vc.genotypes.is_empty() && self.per_sample {
                    for sample in &self.arguments.samples_for_evaluation {
                        let per = self.subset(vc, std::slice::from_ref(sample));
                        if self.by_filter || !is_filtered(&per) {
                            add(sample, per);
                        }
                    }
                }
            }
            if self.arguments.merge_evals && track > 0 && !bindings.is_empty() {
                for (sample, list) in mapping {
                    match bindings[0].iter_mut().find(|(name, _)| *name == sample) {
                        Some((_, existing)) => existing.extend(list),
                        None => bindings[0].push((sample, list)),
                    }
                }
                bindings.push(Vec::new());
            } else {
                bindings.push(mapping);
            }
        }
        bindings
    }

    /// `getSubsetOfVariantContext`.
    fn subset(&self, vc: &VariantContext, samples: &[String]) -> VariantContext {
        let sub = sub_context_from_samples(vc, samples, self.arguments.ignore_ac0());
        ensure_annotations(vc, &sub)
    }

    /// `apply(variantContexts, referenceContext)` for one group of records starting together.
    pub fn apply(
        &mut self,
        group: &[Driving],
        reference: &mut dyn ReferenceBases,
    ) -> Result<(), EvalError> {
        let eval_count = self.arguments.eval_names.len();
        let comp_count = self.arguments.comp_names.len();
        let mut eval_inputs: Vec<Vec<&VariantContext>> = vec![Vec::new(); eval_count];
        let mut comp_inputs: Vec<Vec<&VariantContext>> = vec![Vec::new(); comp_count];
        for driving in group {
            match driving.source {
                Source::Eval(index) => eval_inputs[index].push(&driving.record),
                Source::Comp(index) => comp_inputs[index].push(&driving.record),
            }
        }
        let any_eval = eval_inputs.iter().any(|list| !list.is_empty());
        let first = &group[0].record;
        let start = first.start;
        let contig = first.contig.clone();
        let stop = group
            .iter()
            .map(|driving| driving.record.stop)
            .max()
            .unwrap_or(start);
        let reference_start_bases = reference.bases(&contig, start, start + 1);
        let site = SiteContext {
            comps_here: comp_inputs.clone(),
            reference_start_bases,
            tandem_context: Some(reference.bases(&contig, start, start + 50)),
        };
        let _ = stop;
        let eval_bindings = if any_eval {
            self.bind_evals(&eval_inputs)
        } else {
            vec![Vec::new(); eval_count]
        };
        // Comps are bound with no subsetting and no per-sample split: only `all`.
        let comp_bindings: Vec<Vec<VariantContext>> = comp_inputs
            .iter()
            .map(|records| {
                records
                    .iter()
                    .filter(|vc| self.by_filter || !is_filtered(vc))
                    .map(|vc| (*vc).clone())
                    .collect()
            })
            .collect();

        for (track, eval_name) in self.arguments.eval_names.clone().iter().enumerate() {
            let eval_set = &eval_bindings[track];
            let levels = self.sample_names_for_stratification.clone();
            for level in &levels {
                let by_sample: Vec<Option<VariantContext>> =
                    match eval_set.iter().find(|(name, _)| name == level) {
                        Some((_, list)) => list.iter().cloned().map(Some).collect(),
                        None => vec![None],
                    };
                let present: Vec<&VariantContext> = by_sample.iter().flatten().collect();
                for eval in &by_sample {
                    if comp_count == 0 {
                        self.process_comp(
                            &site,
                            eval.as_ref(),
                            eval_name,
                            None,
                            level,
                            &[],
                            &present,
                        )?;
                    } else {
                        for (comp_index, binding) in
                            comp_bindings.iter().enumerate().take(comp_count)
                        {
                            self.process_comp(
                                &site,
                                eval.as_ref(),
                                eval_name,
                                Some(comp_index),
                                level,
                                binding,
                                &present,
                            )?;
                        }
                    }
                }
            }
            if self.arguments.merge_evals {
                break;
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn process_comp(
        &mut self,
        site: &SiteContext<'_>,
        eval: Option<&VariantContext>,
        eval_name: &str,
        comp_index: Option<usize>,
        level: &str,
        comp_set: &[VariantContext],
        evals_here: &[&VariantContext],
    ) -> Result<(), EvalError> {
        let comp_name = comp_index.map(|index| self.arguments.comp_names[index].clone());
        let comp = self.find_matching_comp(eval, comp_set);
        let family = if level == ALL { Some(ALL) } else { None };
        let keys = self.evaluation_keys(
            site,
            eval,
            eval_name,
            comp,
            comp_name.as_deref(),
            Some(level),
            family,
        )?;
        let samples = self.arguments.samples_for_evaluation.clone();
        for key in keys {
            for evaluator in self.contexts[key].iter_mut() {
                apply_one(evaluator, eval, comp, &self.arguments, &samples)?;
            }
            for other in comp_set {
                let is_matched = comp.is_some_and(|comp| std::ptr::eq(comp, other));
                if !is_matched && !self.comp_has_matching_eval(other, evals_here) {
                    for evaluator in self.contexts[key].iter_mut() {
                        apply_one(evaluator, None, Some(other), &self.arguments, &samples)?;
                    }
                }
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluation_keys(
        &self,
        site: &SiteContext<'_>,
        eval: Option<&VariantContext>,
        eval_name: &str,
        comp: Option<&VariantContext>,
        comp_name: Option<&str>,
        sample: Option<&str>,
        family: Option<&str>,
    ) -> Result<Vec<usize>, EvalError> {
        let mut per_strat: Vec<Vec<usize>> = Vec::new();
        for strat in &self.strats {
            let states = strat.relevant_states(
                &self.arguments,
                site,
                comp,
                comp_name,
                eval,
                eval_name,
                sample,
                family,
            )?;
            let distinct = distinct_states(&strat.states);
            let mut indices = Vec::new();
            for state in states {
                let Some(index) = distinct.iter().position(|candidate| *candidate == state) else {
                    return Err(EvalError::gatk(format!(
                        "Couldn't find state for {} at node StratNode",
                        state.text()
                    )));
                };
                if !indices.contains(&index) {
                    indices.push(index);
                }
            }
            per_strat.push(indices);
        }
        let mut keys = Vec::new();
        let mut current = vec![0usize; per_strat.len()];
        fn walk(
            engine: &Engine,
            per_strat: &[Vec<usize>],
            depth: usize,
            current: &mut Vec<usize>,
            keys: &mut Vec<usize>,
        ) {
            if depth == per_strat.len() {
                let key = engine.key_of(current);
                if !keys.contains(&key) {
                    keys.push(key);
                }
                return;
            }
            for index in &per_strat[depth] {
                current[depth] = *index;
                walk(engine, per_strat, depth + 1, current, keys);
            }
        }
        walk(self, &per_strat, 0, &mut current, &mut keys);
        Ok(keys)
    }

    fn match_type(&self, eval: &VariantContext, comp: &VariantContext) -> u8 {
        // 0 NO_MATCH, 1 STRICT, 2 LENIENT.
        let (eval_type, comp_type) = (vc_type(eval), vc_type(comp));
        if comp_type == VcType::NoVariation || eval_type == VcType::NoVariation {
            return 2;
        }
        if comp_type != eval_type {
            return 0;
        }
        let eval_alt = eval.alleles.get(1);
        let comp_alt = comp.alleles.get(1);
        let strict = match (eval_alt, comp_alt) {
            (None, None) => true,
            (Some(a), Some(b)) => a == b && eval.alleles[0] == comp.alleles[0],
            _ => false,
        };
        if strict {
            1
        } else if self.arguments.require_strict_allele_match {
            0
        } else {
            2
        }
    }

    fn find_matching_comp<'c>(
        &self,
        eval: Option<&VariantContext>,
        comps: &'c [VariantContext],
    ) -> Option<&'c VariantContext> {
        if comps.is_empty() {
            return None;
        }
        let Some(eval) = eval else {
            return comps.first();
        };
        let mut lenient: Option<&VariantContext> = None;
        for comp in comps {
            // `doEvalAndCompMatch(comp, eval, ...)`: the reference passes them the other way round.
            match self.match_type(comp, eval) {
                1 => return Some(comp),
                2 if lenient.is_none() => lenient = Some(comp),
                _ => {}
            }
        }
        lenient
    }

    fn comp_has_matching_eval(&self, comp: &VariantContext, evals: &[&VariantContext]) -> bool {
        evals.iter().any(|eval| self.match_type(comp, eval) != 0)
    }

    /// `finalizeReport`: every module finalized, `MetricsCollection` filled, the report written.
    pub fn finalize_report(&mut self) -> String {
        let processed: Option<i64> = Some(
            self.arguments
                .traversal
                .iter()
                .map(|(_, start, end)| end - start + 1)
                .sum(),
        );
        for context in &mut self.contexts {
            for evaluator in context.iter_mut() {
                evaluator.finalize(processed);
            }
            let mut data = MetricsCollection::default();
            let mut present = [false; 5];
            for evaluator in context.iter() {
                match evaluator {
                    Evaluator::CompOverlap {
                        concordant_rate, ..
                    } => {
                        data.concordant_rate = *concordant_rate;
                        present[0] = true;
                    }
                    Evaluator::IndelSummary(sm) => {
                        data.snps = sm.snps;
                        data.indels = sm.indels;
                        data.indel_ratio = Some(sm.insertion_to_deletion.clone());
                        present[1] = true;
                    }
                    Evaluator::CountVariants(c) => {
                        data.snp_loci = c.snps;
                        data.indel_ratio_loci = c.insertion_deletion_ratio;
                        present[2] = true;
                    }
                    Evaluator::MultiallelicSummary(sm) => {
                        data.indel_loci = sm.indels;
                        present[3] = true;
                    }
                    Evaluator::TiTvVariantEvaluator(t) => {
                        data.titv = t.ratio;
                        present[4] = true;
                    }
                    _ => {}
                }
            }
            if present.iter().all(|seen| *seen) {
                for evaluator in context.iter_mut() {
                    if let Evaluator::MetricsCollection(m) = evaluator {
                        **m = data.clone();
                    }
                }
            }
        }

        let mut report = Report::new();
        for name in &self.evaluator_names {
            let prototype = Evaluator::new(name, &self.arguments.samples_for_evaluation);
            let mut table = Table::new(name, Evaluator::description(name), Sorting::SortByRow);
            table.add_column(name, name);
            for strat in &self.strats {
                table.add_column(strat.name, strat.format);
            }
            match prototype.cells() {
                Ok(cells) => {
                    for (column, format, _) in cells {
                        table.add_column(column, format);
                    }
                }
                Err(_) => {
                    table.add_column("Length", "%d");
                    table.add_column("Freq", "%.2f");
                }
            }
            report.add_table(table);
        }
        for key in 0..self.contexts.len() {
            let states = self.states_of_key(key);
            let key_string: String = self
                .strats
                .iter()
                .zip(&states)
                .map(|(strat, state)| format!("{}:{}", strat.name, state.text()))
                .collect();
            for evaluator in &self.contexts[key] {
                let table = report.table(evaluator.name());
                let set_strats = |table: &mut Table, row: &str| {
                    table.set(
                        row,
                        &table.name.clone(),
                        ReportValue::Str(table.name.clone()),
                    );
                    for (strat, state) in self.strats.iter().zip(&states) {
                        table.set(row, strat.name, state.report());
                    }
                };
                match evaluator.cells() {
                    Ok(cells) => {
                        set_strats(table, &key_string);
                        for (column, _, value) in cells {
                            table.set(&key_string, column, value);
                        }
                    }
                    Err(results) => {
                        for (counter, (length, frequency)) in results.iter().enumerate() {
                            let row = format!("{key_string}{counter:05}");
                            set_strats(table, &row);
                            table.set(&row, "Length", ReportValue::Int(*length));
                            table.set(&row, "Freq", ReportValue::Double(*frequency));
                        }
                    }
                }
            }
        }
        report.write()
    }

    /// The evaluators that need a territory, for `assertThatTerritoryIsSpecifiedIfNecessary`.
    pub fn territory_evaluators(&self) -> Vec<String> {
        self.evaluator_names
            .iter()
            .filter(|name| Evaluator::requires_territory(name))
            .cloned()
            .collect()
    }
}

fn apply_one(
    evaluator: &mut Evaluator,
    eval: Option<&VariantContext>,
    comp: Option<&VariantContext>,
    arguments: &Arguments,
    samples: &[String],
) -> Result<(), EvalError> {
    match evaluator.comparison_order() {
        1 => {
            if let Some(eval) = eval {
                evaluator.update1(eval, arguments)?;
            }
        }
        _ => evaluator.update2(eval, comp, arguments, samples)?,
    }
    Ok(())
}

/// The states a strat node keys on: a `LinkedHashMap` put per state, so a repeated state is one.
fn distinct_states(states: &[State]) -> Vec<State> {
    let mut out: Vec<State> = Vec::new();
    for state in states {
        if !out.contains(state) {
            out.push(state.clone());
        }
    }
    out
}

fn make_stratifier(
    name: &str,
    arguments: &Arguments,
    sample_names_for_stratification: &[String],
    contigs: &[String],
) -> Result<Stratifier, EvalError> {
    let strings = |values: &[&str]| values.iter().map(|value| s(value)).collect::<Vec<_>>();
    let (states, kind, format): (Vec<State>, StratKind, &'static str) = match name {
        "AlleleCount" => {
            if arguments.eval_names.len() != 1 && !arguments.merge_evals {
                return Err(EvalError::bad_argument(
                    "AlleleCount",
                    "AlleleCount stratification only works with a single eval vcf",
                ));
            }
            let samples = if arguments.samples_for_evaluation.is_empty() {
                arguments.num_samples_from_argument
            } else {
                arguments.samples_for_evaluation.len() as i64
            };
            let nchrom = samples * arguments.ploidy;
            if nchrom < 2 {
                return Err(EvalError::bad_argument(
                    "AlleleCount",
                    "AlleleCount stratification requires an eval vcf with at least one sample",
                ));
            }
            (
                (0..=nchrom).map(State::Int).collect(),
                StratKind::AlleleCount { nchrom },
                "%d",
            )
        }
        "AlleleFrequency" => {
            let states = match arguments.af_scale {
                AfScale::Linear => {
                    let mut states = Vec::new();
                    let mut a = 0.000f64;
                    while a <= 1.005 {
                        states.push(State::Str(gatk_engine::java_format::format_decimals(a, 3)));
                        a += 0.005;
                    }
                    states
                }
                AfScale::Logarithmic => {
                    (-30..=30).map(|a: i64| State::Str(a.to_string())).collect()
                }
            };
            (
                states,
                StratKind::AlleleFrequency {
                    scale: arguments.af_scale,
                    use_comp: arguments.use_comp_af,
                },
                "%s",
            )
        }
        "CompFeatureInput" => {
            let mut states: Vec<State> = arguments
                .comp_names
                .iter()
                .map(|name| State::Str(name.clone()))
                .collect();
            if states.is_empty() {
                states.push(s("none"));
            }
            (states, StratKind::CompFeatureInput, "%s")
        }
        "Contig" => {
            let mut sorted: Vec<String> = contigs.to_vec();
            sorted.sort();
            sorted.dedup();
            let mut states: Vec<State> = sorted.into_iter().map(State::Str).collect();
            states.push(s(ALL));
            (states, StratKind::Contig, "%s")
        }
        "CpG" => (strings(&["all", "CpG", "non_CpG"]), StratKind::CpG, "%s"),
        "Degeneracy" => (
            strings(&["1-fold", "2-fold", "3-fold", "4-fold", "6-fold", "all"]),
            StratKind::Degeneracy,
            "%s",
        ),
        "EvalFeatureInput" => {
            let mut states = Vec::new();
            for eval in &arguments.eval_names {
                states.push(State::Str(eval.clone()));
                if arguments.merge_evals {
                    break;
                }
            }
            (states, StratKind::EvalFeatureInput, "%s")
        }
        "Family" => (vec![s(ALL)], StratKind::Family, "%s"),
        "Filter" => (
            strings(&["called", "filtered", "raw"]),
            StratKind::Filter,
            "%s",
        ),
        "FilterType" => {
            let order = gatk_engine::java_hash::hash_set_order(&arguments.eval_filter_names)
                .unwrap_or_else(|_| arguments.eval_filter_names.clone());
            let mut states: Vec<State> = order.into_iter().map(State::Str).collect();
            states.push(s("PASS"));
            (states, StratKind::FilterType, "%s")
        }
        "FunctionalClass" => (
            strings(&["all", "silent", "missense", "nonsense"]),
            StratKind::FunctionalClass,
            "%s",
        ),
        "IndelSize" => (
            (-100..=100).map(State::Int).collect(),
            StratKind::IndelSize,
            "%d",
        ),
        "IntervalStratification" => {
            if arguments.strat_intervals.is_none() {
                return Err(EvalError::command_line(
                    "Argument stratIntervals was missing: Must be provided when \
                     IntervalStratification is enabled"
                        .to_string(),
                ));
            }
            (
                strings(&["all", "overlaps.intervals", "outside.intervals"]),
                StratKind::IntervalStratification,
                "%s",
            )
        }
        "JexlExpression" => {
            let mut states = vec![s("none")];
            let mut names: Vec<String> = arguments
                .selects
                .iter()
                .map(|select| select.name.clone())
                .collect();
            names.sort();
            names.dedup();
            for select in names {
                states.push(State::Str(select));
            }
            (states, StratKind::JexlExpression, "%s")
        }
        "Novelty" => (
            strings(&["all", "known", "novel"]),
            StratKind::Novelty,
            "%s",
        ),
        "OneBPIndel" => (
            strings(&["all", "one.bp", "two.plus.bp"]),
            StratKind::OneBPIndel,
            "%s",
        ),
        "Sample" => (
            sample_names_for_stratification
                .iter()
                .map(|name| State::Str(name.clone()))
                .collect(),
            StratKind::Sample,
            "%s",
        ),
        "SnpEffPositionModifier" => (
            strings(&[
                "GENE",
                "CODING_REGION",
                "SPLICE_SITE",
                "STOP_GAINED",
                "STOP_LOST",
            ]),
            StratKind::SnpEffPositionModifier,
            "%s",
        ),
        "TandemRepeat" => (
            strings(&["all", "is.repeat", "not.repeat"]),
            StratKind::TandemRepeat,
            "%s",
        ),
        "VariantType" => (
            VcType::ALL.iter().map(|kind| s(kind.name())).collect(),
            StratKind::VariantType,
            "%s",
        ),
        other => unreachable!("no stratifier {other}"),
    };
    let name: &'static str = STRATIFIER_NAMES
        .iter()
        .find(|(candidate, _, _)| *candidate == name)
        .map(|(candidate, _, _)| *candidate)
        .expect("a known stratifier");
    Ok(Stratifier {
        name,
        states,
        kind,
        format,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_linear_allele_frequency_bins_step_by_five_thousandths() {
        let arguments_states: Vec<String> = {
            let mut states = Vec::new();
            let mut a = 0.000f64;
            while a <= 1.005 {
                states.push(gatk_engine::java_format::format_decimals(a, 3));
                a += 0.005;
            }
            states
        };
        // `for (double a = 0.000; a <= 1.005; a += 0.005)`: the sum drifts above 1.005 after
        // 1.000, so there are 201 states, which is what the reference's report carries.
        assert_eq!(arguments_states.first().map(String::as_str), Some("0.000"));
        assert_eq!(arguments_states.last().map(String::as_str), Some("1.000"));
        assert_eq!(arguments_states.len(), 201);
    }

    #[test]
    fn a_ratio_or_a_percent_of_nothing_is_na() {
        assert_eq!(formatted_percent(1, 0), "NA");
        assert_eq!(formatted_ratio(3, 2), "1.50");
    }

    #[test]
    fn snpeff_effects_inherit_their_parents() {
        assert!(snpeff_is_subtype("STOP_GAINED", "CDS"));
        assert!(snpeff_is_subtype("CDS", "EXON"));
        assert!(!snpeff_is_subtype("INTRON", "EXON"));
    }
}

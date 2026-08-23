//! `FilterFuncotations`, ported from the tool and its `filtrationRules` package (GATK 4.6.2.0).
//!
//! A Funcotated VCF read back and marked with the clinical-significance rules that match it. Five
//! filters run over the funcotations of each transcript; what comes out is one INFO field naming
//! the filters that matched and one FILTER value when none did.
//!
//! # A filter matches only when every one of its rules does
//!
//! ```java
//! return getRules().stream()
//!         .map(rule -> rule.checkRule(prunedTranscriptFuncotations, variant))
//!         .reduce(Boolean::logicalAnd)
//!         .orElse(false);
//! ```
//!
//! so a filter with no rules is false rather than true, and the funcotations a rule sees have
//! already had their empty values pruned, which is why a missing key and an empty one are the same
//! thing here.
//!
//! # Three kinds of missing ExAC data all read as a frequency of zero
//!
//! An allele number of zero is answered with zero rather than divided by, an allele count that does
//! not parse is caught and answered with zero, and an allele count with no allele number beside it
//! takes the `orElse(0)` and lands on the first case. All three make the variant PASS the frequency
//! rule, which is the permissive direction.
//!
//! # The gnomAD path is not the same code
//!
//! It has a rule ExAC has not: a dataset counts as present only when it carries no `_FILTER`
//! funcotation or when that funcotation says `PASS`, and when neither dataset is present the rule
//! fails whatever the frequencies say. And it catches nothing, so a frequency that does not parse
//! is a `NumberFormatException` out of the tool rather than a zero. [`gnomad_max_maf`] returns a
//! `Result` for that reason and [`exac_max_maf`] does not.
//!
//! # The two autosomal recessive filters share one name
//!
//! Both are registered as `AR`, and the matching names are collected into a `HashSet` before being
//! joined, so a variant matching both contributes one entry. That set is also what decides the
//! order of the joined value, which is not the order the filters were registered in.
//!
//! # And the het rule needs a second pass
//!
//! A gene needs more than one het variant before any of them is compound het, so a lone het variant
//! in an autosomal recessive gene matches nothing. The reference adds a variant to its gene's list
//! ONCE PER TRANSCRIPT that names the gene, so a variant whose funcotations name the same gene on
//! two transcripts reaches a list of two on its own; [`compound_het_variants`] keeps that, though
//! no case in the golden has two transcripts.

use gatk_engine::java_hash::hash_set_order;

/// `FilterFuncotationsConstants`.
pub const CLINSIG_INFO_KEY: &str = "CLINSIG";
pub const CLINSIG_INFO_NOT_SIGNIFICANT: &str = "NONE";
pub const NOT_CLINSIG_FILTER: &str = "NOT_CLINSIG";
pub const FILTER_DELIMITER: &str = ",";

/// `AutosomalRecessiveConstants`.
pub const AR_INFO_VALUE: &str = "AR";
pub const AUTOSOMAL_RECESSIVE_GENES: [&str; 2] = ["ATP7B", "MUTYH"];

pub const CLINVAR_INFO_VALUE: &str = "CLINVAR";
pub const LOF_INFO_VALUE: &str = "LOF";
pub const LMM_INFO_VALUE: &str = "LMM";

const CLINVAR_MAX_MAF: f64 = 0.05;
const LOF_MAX_MAF: f64 = 0.01;

const ACMG_DISEASE_FUNCOTATION: &str = "ACMG_recommendation_Disease_Name";
const CLINVAR_SIGNIFICANCE_FUNCOTATION: &str = "ClinVar_VCF_CLNSIG";
const CLINVAR_SIGNIFICANCE_MATCHING_VALUES: [&str; 3] = [
    "Pathogenic",
    "Likely_pathogenic",
    "Pathogenic/Likely_pathogenic",
];
const LOF_GENE_FUNCOTATION: &str = "ACMGLMMLof_LOF_Mechanism";
const LMM_FLAGGED: &str = "LMMKnown_LMM_FLAGGED";

/// The five `GencodeFuncotation.VariantClassification`s the LOF filter matches.
const LOF_CLASSIFICATIONS: [&str; 5] = [
    "FRAME_SHIFT_DEL",
    "FRAME_SHIFT_INS",
    "NONSENSE",
    "START_CODON_DEL",
    "SPLICE_SITE",
];

const EXAC_SUB_POPULATIONS: [&str; 7] = ["AFR", "AMR", "EAS", "FIN", "NFE", "OTH", "SAS"];

const GNOMAD_DATASETS: [&str; 2] = ["gnomAD_genome", "gnomAD_exome"];

const GNOMAD_SUB_POPULATIONS: [&str; 37] = [
    "afr",
    "afr_female",
    "afr_male",
    "amr",
    "amr_female",
    "amr_male",
    "asj",
    "asj_female",
    "asj_male",
    "eas",
    "eas_female",
    "eas_jpn",
    "eas_kor",
    "eas_male",
    "eas_oea",
    "female",
    "fin",
    "fin_female",
    "fin_male",
    "male",
    "nfe",
    "nfe_bgr",
    "nfe_est",
    "nfe_female",
    "nfe_male",
    "nfe_nwe",
    "nfe_onf",
    "nfe_seu",
    "nfe_swe",
    "oth",
    "oth_female",
    "oth_male",
    "popmax",
    "raw",
    "sas",
    "sas_female",
    "sas_male",
];

/// `--ref-version`, whose gencode number is part of every key the tool looks for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Reference {
    B37,
    Hg19,
    Hg38,
}

impl Reference {
    pub fn gencode_version(&self) -> i32 {
        match self {
            Reference::B37 | Reference::Hg19 => 19,
            Reference::Hg38 => 27,
        }
    }

    fn key(&self, suffix: &str) -> String {
        format!("Gencode_{}_{suffix}", self.gencode_version())
    }
}

/// `--allele-frequency-data-source`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlleleFrequencySource {
    Exac,
    Gnomad,
}

/// One transcript's funcotations, already pruned of empty values as
/// `FilterFuncotationsUtils.getTranscriptFuncotations` prunes them.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Funcotations {
    entries: Vec<(String, String)>,
}

impl Funcotations {
    /// Builds the set, dropping every empty value the way the reference does.
    pub fn new(entries: &[(&str, &str)]) -> Self {
        Funcotations {
            entries: entries
                .iter()
                .filter(|(_, value)| !value.is_empty())
                .map(|(key, value)| (key.to_string(), value.to_string()))
                .collect(),
        }
    }

    pub fn contains_key(&self, key: &str) -> bool {
        self.entries.iter().any(|(name, _)| name == key)
    }

    pub fn get(&self, key: &str) -> Vec<&str> {
        self.entries
            .iter()
            .filter(|(name, _)| name == key)
            .map(|(_, value)| value.as_str())
            .collect()
    }

    /// `matchOnKeyOrDefault`, which hands back the default when nothing matched.
    fn matches_or_default<'a>(&'a self, key: &str, default: &'a str) -> Vec<&'a str> {
        let matched = self.get(key);
        if matched.is_empty() {
            vec![default]
        } else {
            matched
        }
    }
}

/// A variant, reduced to what the rules ask it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub reference_allele: String,
    pub alternate_alleles: Vec<String>,
    pub het_count: i32,
    pub hom_var_count: i32,
}

/// What `variantContextsMatch` compares, which is everything but the attributes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VariantKey {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub reference_allele: String,
    pub alternate_alleles: Vec<String>,
}

impl Variant {
    pub fn key(&self) -> VariantKey {
        VariantKey {
            contig: self.contig.clone(),
            start: self.start,
            end: self.end,
            reference_allele: self.reference_allele.clone(),
            alternate_alleles: self.alternate_alleles.clone(),
        }
    }
}

impl VariantKey {
    /// `variantContextsMatch`, whose alternate-allele test is a size check and a containsAll.
    pub fn matches(&self, other: &VariantKey) -> bool {
        self.contig == other.contig
            && self.start == other.start
            && self.end == other.end
            && self.reference_allele == other.reference_allele
            && self.alternate_alleles.len() == other.alternate_alleles.len()
            && other
                .alternate_alleles
                .iter()
                .all(|allele| self.alternate_alleles.contains(allele))
    }
}

/// A gnomAD frequency that does not parse, which the reference does not catch.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NumberFormatError {
    pub value: String,
}

impl NumberFormatError {
    pub fn java_class(&self) -> &str {
        "java.lang.NumberFormatException"
    }

    pub fn message(&self) -> String {
        format!("For input string: \"{}\"", self.value)
    }
}

/// `AlleleFrequencyExacUtils.getMaxMinorAlleleFreq`, which never fails.
pub fn exac_max_maf(funcotations: &Funcotations) -> f64 {
    let mut max: Option<f64> = None;
    for population in EXAC_SUB_POPULATIONS {
        let count_key = format!("ExAC_AC_{population}");
        if !funcotations.contains_key(&count_key) {
            continue;
        }
        let number_key = format!("ExAC_AN_{population}");
        let frequency = match (
            funcotations
                .get(&count_key)
                .first()
                .map(|v| v.parse::<f64>()),
            funcotations
                .get(&number_key)
                .first()
                .map(|v| v.parse::<i32>()),
        ) {
            // Either value failing to parse is caught and logged, and the whole sub-population
            // answers zero.
            (Some(Err(_)), _) | (_, Some(Err(_))) => 0.0,
            (count, number) => {
                let count = count.map(|value| value.expect("checked")).unwrap_or(0.0);
                let number = number.map(|value| value.expect("checked")).unwrap_or(0);
                if number == 0 {
                    // A variant never seen in ExAC is reported as 0% MAF.
                    0.0
                } else {
                    count / f64::from(number)
                }
            }
        };
        max = Some(match max {
            Some(previous) if previous >= frequency => previous,
            _ => frequency,
        });
    }
    max.unwrap_or(0.0)
}

/// `AlleleFrequencyGnomadUtils.datasetsPresent`.
fn gnomad_datasets_present(funcotations: &Funcotations) -> Vec<&'static str> {
    GNOMAD_DATASETS
        .into_iter()
        .filter(|dataset| {
            let filter_key = format!("{dataset}_FILTER");
            !funcotations.contains_key(&filter_key)
                || funcotations.get(&filter_key).contains(&"PASS")
        })
        .collect()
}

/// `AlleleFrequencyGnomadUtils.allFrequenciesFiltered`.
pub fn gnomad_all_frequencies_filtered(funcotations: &Funcotations) -> bool {
    gnomad_datasets_present(funcotations).is_empty()
}

/// `AlleleFrequencyGnomadUtils.getMaxMinorAlleleFreq`, which parses without catching.
pub fn gnomad_max_maf(funcotations: &Funcotations) -> Result<f64, NumberFormatError> {
    let mut max: Option<f64> = None;
    for dataset in gnomad_datasets_present(funcotations) {
        for population in GNOMAD_SUB_POPULATIONS {
            let key = format!("{dataset}_AF_{population}");
            for value in funcotations.get(&key) {
                let frequency = value.parse::<f64>().map_err(|_| NumberFormatError {
                    value: value.to_string(),
                })?;
                max = Some(match max {
                    Some(previous) if previous >= frequency => previous,
                    _ => frequency,
                });
            }
        }
    }
    Ok(max.unwrap_or(0.0))
}

/// `AlleleFrequencyUtils.buildMaxMafRule`.
pub fn max_maf_rule(
    funcotations: &Funcotations,
    max_maf: f64,
    source: AlleleFrequencySource,
) -> Result<bool, NumberFormatError> {
    match source {
        AlleleFrequencySource::Exac => Ok(exac_max_maf(funcotations) <= max_maf),
        AlleleFrequencySource::Gnomad => Ok(!gnomad_all_frequencies_filtered(funcotations)
            && gnomad_max_maf(funcotations)? <= max_maf),
    }
}

/// `ClinVarFilter`.
pub fn clinvar_filter(
    funcotations: &Funcotations,
    source: AlleleFrequencySource,
) -> Result<bool, NumberFormatError> {
    let on_acmg_list = funcotations.contains_key(ACMG_DISEASE_FUNCOTATION);
    let significance = funcotations.matches_or_default(CLINVAR_SIGNIFICANCE_FUNCOTATION, "");
    // An exact match against three values, so a value carrying anything else does not match.
    let significant = CLINVAR_SIGNIFICANCE_MATCHING_VALUES
        .iter()
        .any(|wanted| significance.iter().any(|value| value == wanted));
    Ok(on_acmg_list && significant && max_maf_rule(funcotations, CLINVAR_MAX_MAF, source)?)
}

/// `LofFilter`.
pub fn lof_filter(
    funcotations: &Funcotations,
    reference: Reference,
    source: AlleleFrequencySource,
) -> Result<bool, NumberFormatError> {
    let classification = reference.key("variantClassification");
    let classified = funcotations
        .matches_or_default(&classification, "")
        .iter()
        .any(|value| LOF_CLASSIFICATIONS.contains(value));
    let mechanism = funcotations
        .matches_or_default(LOF_GENE_FUNCOTATION, "")
        .contains(&"YES");
    Ok(classified && mechanism && max_maf_rule(funcotations, LOF_MAX_MAF, source)?)
}

/// `LmmFilter`, whose flag is read with `Boolean.valueOf`.
pub fn lmm_filter(funcotations: &Funcotations) -> bool {
    funcotations
        .matches_or_default(LMM_FLAGGED, "false")
        .iter()
        .any(|value| value.eq_ignore_ascii_case("true"))
}

/// `ArHomvarFilter`, which asks the variant rather than the funcotations for its hom-var count.
pub fn ar_homvar_filter(
    funcotations: &Funcotations,
    reference: Reference,
    variant: &Variant,
) -> bool {
    let gene_key = reference.key("hugoSymbol");
    let interesting = funcotations
        .get(&gene_key)
        .iter()
        .any(|gene| AUTOSOMAL_RECESSIVE_GENES.contains(gene));
    interesting && variant.hom_var_count > 0
}

/// The first pass of `ArHetvarFilter`: the variants of every autosomal recessive gene that carries
/// more than one het variant.
///
/// The reference appends a variant once per transcript whose funcotations name the gene, so a
/// variant naming the same gene on two transcripts fills a list of two by itself.
pub fn compound_het_variants(
    variants: &[(Variant, Vec<Funcotations>)],
    reference: Reference,
) -> Vec<VariantKey> {
    let gene_key = reference.key("hugoSymbol");
    let mut by_gene: Vec<(String, Vec<VariantKey>)> = Vec::new();
    for (variant, transcripts) in variants {
        for funcotations in transcripts {
            // findFirst, so only the first gene funcotation of a transcript is looked at.
            let Some(gene) = funcotations
                .get(&gene_key)
                .first()
                .map(|gene| gene.to_string())
            else {
                continue;
            };
            if !AUTOSOMAL_RECESSIVE_GENES.contains(&gene.as_str()) || variant.het_count <= 0 {
                continue;
            }
            match by_gene.iter_mut().find(|(name, _)| *name == gene) {
                Some((_, keys)) => keys.push(variant.key()),
                None => by_gene.push((gene, vec![variant.key()])),
            }
        }
    }
    by_gene
        .into_iter()
        .filter(|(_, keys)| keys.len() > 1)
        .flat_map(|(_, keys)| keys)
        .collect()
}

/// `ArHetvarFilter`'s rule.
pub fn ar_hetvar_filter(compound: &[VariantKey], variant: &Variant) -> bool {
    let key = variant.key();
    compound.iter().any(|candidate| candidate.matches(&key))
}

/// The names of every filter matching any transcript of a variant, in the order the `HashSet`
/// they were collected into hands them over.
pub fn matching_filters(
    transcripts: &[Funcotations],
    variant: &Variant,
    reference: Reference,
    source: AlleleFrequencySource,
    compound: &[VariantKey],
) -> Result<Vec<String>, NumberFormatError> {
    let mut matched: Vec<String> = Vec::new();
    for funcotations in transcripts {
        let mut names: Vec<&str> = Vec::new();
        if clinvar_filter(funcotations, source)? {
            names.push(CLINVAR_INFO_VALUE);
        }
        if lof_filter(funcotations, reference, source)? {
            names.push(LOF_INFO_VALUE);
        }
        if lmm_filter(funcotations) {
            names.push(LMM_INFO_VALUE);
        }
        if ar_homvar_filter(funcotations, reference, variant) || ar_hetvar_filter(compound, variant)
        {
            names.push(AR_INFO_VALUE);
        }
        for name in names {
            if !matched.iter().any(|existing| existing == name) {
                matched.push(name.to_string());
            }
        }
    }
    Ok(hash_set_order(&matched).expect("four filter names do not treeify a bucket"))
}

/// `applyFilters`: the CLINSIG value, and whether the variant passes.
pub fn clinsig(matching: &[String]) -> (String, bool) {
    if matching.is_empty() {
        (CLINSIG_INFO_NOT_SIGNIFICANT.to_string(), false)
    } else {
        (matching.join(FILTER_DELIMITER), true)
    }
}

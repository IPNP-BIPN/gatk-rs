//! `VariantAnnotator`: the resolved annotations run over a VCF, with likelihoods built from reads.
//!
//! Ported from `org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotator` and the
//! parts of `VariantAnnotatorEngine` it reaches (GATK 4.6.2.0): `annotateContext` with read
//! likelihoods and no fragment or haplotype likelihoods, the expressions over `--resource`
//! inputs, and the overlap annotator over `--dbsnp` and `--comp`.
//!
//! # The likelihoods are a pileup, not a model
//!
//! `makeLikelihoods` gives every read of a sample a row, then walks the pileup at the record's
//! start and sets, for the element at position `i` of that sample's pileup, the row `i`: minus
//! infinity for every allele and zero for the one `chooseAlleleForRead` picks. The row it sets is
//! the element's rank in the PILEUP, not the read's rank in the sample's read list, so a read that
//! overlaps the record without reaching its first base shifts every later assignment onto the
//! wrong read. That is reproduced, and the rows nothing reaches keep the zeros the matrix was
//! created with.
//!
//! # Which annotations run
//!
//! The info and genotype annotations this port carries with reads. The jumbo annotations need
//! fragment and haplotype likelihoods, which this walker never has, so they contribute their
//! header lines and nothing else, exactly as upstream. An annotation the port does not carry is
//! refused at the first record rather than skipped.
//!
//! # An annotated record is always a decoded one
//!
//! `annotateContext` rebuilds the record with `builder.genotypes(...)`, and `make()` validates
//! genotypes it was handed by iterating them. So every record that reaches the engine is written
//! with its FORMAT keys recomputed and sorted, even when nothing was asked for; only a record the
//! walker skips, over an ambiguous reference base, keeps the file's own genotype text.

use gatk_annotation::catalogue::{Entry, Kind};
use gatk_annotation::info_annotation::{AnnotationValue, InfoFieldAnnotation};
use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::allele_list::{AlleleList, SampleList};
use gatk_engine::java_random::JavaRandom;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::genotypes_context::GenotypesContext;
use htsjdk_vcf::header::Cardinality;
use htsjdk_vcf::variant::{Value, VariantContext};

use crate::genotype_gvcfs::{attribute, attribute_int, put, value_of};
use crate::genotyping_engine::EngineError;

/// `VCFConstants.DBSNP_KEY`.
pub const DBSNP_KEY: &str = "DB";

/// `VariantAnnotatorEngine.VAExpression`: `name.field`, bound to one `--resource`.
#[derive(Debug, Clone)]
pub struct Expression {
    /// The whole expression, which is the INFO key it writes.
    pub full_name: String,
    /// What follows the last dot.
    pub field_name: String,
    /// The index of the `--resource` it names.
    pub binding: usize,
    /// The count of the header line the expression writes, which `sethInfo` records once the
    /// line exists. `None` for `ID`, which never reads it.
    pub count: Option<Cardinality>,
}

/// `VariantAnnotatorEngine.VAExpression(fullExpression, dataSourceList)`.
pub fn expression(text: &str, resource_names: &[String]) -> Result<Expression, EngineError> {
    let Some(dot) = text.rfind('.') else {
        return Err(bad_input(format!(
            "The requested expression '{text}' is invalid, it should be in VCFFile.value format"
        )));
    };
    let binding_name = &text[..dot];
    let Some(binding) = resource_names.iter().position(|name| name == binding_name) else {
        return Err(bad_input(format!(
            "The requested expression '{text}' is invalid, could not find vcf input file"
        )));
    };
    Ok(Expression {
        full_name: text.to_string(),
        field_name: text[dot + 1..].to_string(),
        binding,
        count: None,
    })
}

fn bad_input(message: String) -> EngineError {
    EngineError::Runtime {
        class: "org.broadinstitute.hellbender.exceptions.UserException$BadInput".to_string(),
        message: format!("Bad input: {message}"),
    }
}

fn limitation(what: String) -> EngineError {
    EngineError::Limitation(what)
}

/// What one record is annotated against, beside the record itself.
pub struct Site<'a> {
    /// The read likelihoods `makeLikelihoods` built.
    pub likelihoods: &'a AlleleLikelihoods<BamRecord>,
    /// The reference context's window start and bases: the expanded window when there is a
    /// reference, and the record's own start with no bases when there is none, which is the
    /// empty context the walker still hands over.
    pub window: (i64, &'a [u8]),
    /// Per overlap name, in the engine's order, the records `getValues(input, start)` returns.
    pub overlaps: Vec<Vec<&'a VariantContext>>,
    /// The `--dbsnp` records starting here, when `--dbsnp` was given.
    pub dbsnp: Option<Vec<&'a VariantContext>>,
    /// Per `--resource`, the records starting here.
    pub resources: Vec<Vec<&'a VariantContext>>,
}

/// `VariantAnnotatorEngine` as `VariantAnnotator` builds it.
pub struct Engine {
    /// The resolved annotations, sorted by simple name as the plugin descriptor hands them out.
    pub resolved: Vec<&'static Entry>,
    /// The overlap names: every `--comp`'s, then `DB` when `--dbsnp` was given.
    pub overlap_names: Vec<String>,
    pub expressions: Vec<Expression>,
    /// `--resource-allele-concordance`.
    pub allele_concordance: bool,
}

impl Engine {
    fn info(&self) -> impl Iterator<Item = &'static Entry> + '_ {
        self.resolved
            .iter()
            .copied()
            .filter(|entry| entry.kind == Kind::Info)
    }

    fn genotype(&self) -> impl Iterator<Item = &'static Entry> + '_ {
        self.resolved
            .iter()
            .copied()
            .filter(|entry| entry.kind == Kind::Genotype)
    }

    /// `annotateContext(vc, features, ref, likelihoods, a -> true)`.
    pub fn annotate_context(
        &self,
        vc: &VariantContext,
        site: &Site<'_>,
        random: &mut JavaRandom,
    ) -> Result<VariantContext, EngineError> {
        // `annotateGenotypes`: untouched when no genotype annotation, jumbo or not, was asked
        // for, rebuilt otherwise, each annotation reading the ORIGINAL genotype and writing the
        // builder. A jumbo one never runs here and still costs the file's own genotype text.
        let mut genotype_annotated = vc.clone();
        let genotype_entries: Vec<&Entry> = self.genotype().collect();
        let jumbo_genotype = self
            .resolved
            .iter()
            .any(|entry| entry.kind == Kind::JumboGenotype);
        if !genotype_entries.is_empty() || jumbo_genotype {
            let mut rebuilt = Vec::with_capacity(vc.genotypes.len());
            for genotype in vc.genotypes.iter() {
                let mut built = genotype.clone();
                let called = genotype.alleles.iter().any(|allele| !allele.is_no_call());
                for entry in &genotype_entries {
                    annotate_genotype(entry, vc, genotype, &mut built, called, site.likelihoods)?;
                }
                rebuilt.push(built);
            }
            genotype_annotated.genotypes = GenotypesContext::new(rebuilt);
        }

        // `addInfoAnnotations`: the record's attributes, the expressions, then each annotation.
        let mut attributes = genotype_annotated.attributes.clone();
        self.annotate_expressions(vc, site, &mut attributes)?;
        for entry in self.info() {
            for (key, value) in annotate_info(entry, &genotype_annotated, site, random)? {
                put(&mut attributes, &key, value);
            }
        }
        let mut annotated = genotype_annotated;
        annotated.attributes = attributes;
        // `builder.genotypes(...)` asks `make()` to validate the genotypes, and validating them
        // iterates them: the record is written from its decoded genotypes, never from the file's
        // own text, whatever was annotated.
        let _ = annotated.genotypes.first();

        // `annotateOverlaps(features, annotateRsID(features, annotated))`.
        if let Some(sources) = &site.dbsnp {
            let id = crate::variant_overlap::rs_id(sources, &annotated).map_err(runtime_iae)?;
            if let Some(id) = id {
                if annotated.id == "." {
                    annotated.id = id;
                } else if !annotated.id.contains(&id) {
                    annotated.id = format!("{};{id}", annotated.id);
                }
            }
        }
        for (name, sources) in self.overlap_names.iter().zip(&site.overlaps) {
            let overlaps = crate::variant_overlap::rs_id(sources, &annotated)
                .map_err(runtime_iae)?
                .is_some();
            if overlaps {
                put(&mut annotated.attributes, name, Value::Bool(true));
            }
        }
        Ok(annotated)
    }

    /// `annotateExpressions`.
    fn annotate_expressions(
        &self,
        vc: &VariantContext,
        site: &Site<'_>,
        attributes: &mut Vec<(String, Value)>,
    ) -> Result<(), EngineError> {
        for expression in &self.expressions {
            let Some(source) = site.resources[expression.binding].first() else {
                continue;
            };
            match expression.field_name.as_str() {
                "ID" => {
                    if source.id != "." {
                        put(
                            attributes,
                            &expression.full_name,
                            Value::Str(source.id.clone()),
                        );
                    }
                }
                "ALT" => {
                    // `getAlternateAllele(0)`, which is `alleles.get(1)`.
                    let alternate = source.alleles.get(1).ok_or_else(|| EngineError::Runtime {
                        class: "java.lang.IndexOutOfBoundsException".to_string(),
                        message: format!(
                            "Index 1 out of bounds for length {}",
                            source.alleles.len()
                        ),
                    })?;
                    put(
                        attributes,
                        &expression.full_name,
                        Value::Str(alternate.display_string()),
                    );
                }
                "FILTER" => {
                    // `getFilters()` is the `HashSet` the codec decoded them into, joined in its
                    // own order rather than the file's.
                    let filters = source.filters.clone().unwrap_or_default();
                    let text = if filters.is_empty() {
                        "PASS".to_string()
                    } else {
                        gatk_engine::java_hash::hash_set_order(&filters)
                            .map_err(|error| limitation(format!("{error:?}")))?
                            .join(",")
                    };
                    put(attributes, &expression.full_name, Value::Str(text));
                }
                field => {
                    let Some(value) = attribute(source, field) else {
                        continue;
                    };
                    let count = expression.count;
                    let use_ref_and_alt = count == Some(Cardinality::R);
                    let use_alt = count == Some(Cardinality::A);
                    if use_alt || use_ref_and_alt || self.allele_concordance {
                        let cleaned: String =
                            crate::reference_confidence_merger::value_to_string(value)
                                .chars()
                                .filter(|c| !matches!(c, '[' | ']') && !c.is_whitespace())
                                .collect();
                        let values: Vec<&str> = cleaned.split(',').collect();
                        let mine = min_representation_biallelics(vc)?;
                        let theirs = min_representation_biallelics(source)?;
                        let mut can_annotate = false;
                        let mut annotation_values: Vec<Value> = Vec::new();
                        for biallelic in &mine {
                            let mut concordant = false;
                            let mut i = 0usize;
                            for other in &theirs {
                                if other == biallelic {
                                    if i == 0 && use_ref_and_alt {
                                        annotation_values.push(Value::Str(value_at(&values, i)?));
                                        i += 1;
                                    }
                                    annotation_values.push(Value::Str(value_at(&values, i)?));
                                    concordant = true;
                                    can_annotate = true;
                                    break;
                                }
                                i += 1;
                            }
                            if !concordant {
                                annotation_values.push(Value::Str("0".to_string()));
                            }
                        }
                        if can_annotate {
                            // `VCFEncoder` writes a key declared `Number=0` bare whatever its
                            // value, so a flag carried through this path is still a flag.
                            let value = if count == Some(Cardinality::Fixed(0)) {
                                Value::Bool(true)
                            } else {
                                Value::List(annotation_values)
                            };
                            put(attributes, &expression.full_name, value);
                        }
                    } else {
                        put(attributes, &expression.full_name, value.clone());
                    }
                }
            }
        }
        Ok(())
    }
}

fn value_at(values: &[&str], index: usize) -> Result<String, EngineError> {
    values
        .get(index)
        .map(|value| value.to_string())
        .ok_or_else(|| EngineError::Runtime {
            class: "java.lang.IndexOutOfBoundsException".to_string(),
            message: format!("Index {index} out of bounds for length {}", values.len()),
        })
}

fn runtime_iae(message: String) -> EngineError {
    EngineError::Runtime {
        class: "java.lang.IllegalArgumentException".to_string(),
        message,
    }
}

/// `getMinRepresentationBiallelics`, as the allele lists the concordance test compares: a record
/// with at most two alleles as it stands, and each reference and alternate pair of a
/// multi-allelic one, trimmed at both ends unless both are one base long.
fn min_representation_biallelics(vc: &VariantContext) -> Result<Vec<Vec<Allele>>, EngineError> {
    if vc.alleles.len() <= 2 {
        return Ok(vec![vc.alleles.clone()]);
    }
    let reference = &vc.alleles[0];
    let mut pairs = Vec::new();
    for alternate in &vc.alleles[1..] {
        if java_length(reference) == 1 && java_length(alternate) == 1 {
            pairs.push(vec![reference.clone(), alternate.clone()]);
        } else {
            let (_, trimmed_reference, trimmed_alternate) =
                crate::variant_overlap::trim_pair(vc, alternate).map_err(runtime_iae)?;
            pairs.push(vec![trimmed_reference, trimmed_alternate]);
        }
    }
    Ok(pairs)
}

/// `Allele.length()`, which is zero for a symbolic allele.
fn java_length(allele: &Allele) -> usize {
    if allele.is_symbolic() {
        0
    } else {
        allele.len()
    }
}

fn to_values(pairs: Vec<(String, AnnotationValue)>) -> Vec<(String, Value)> {
    pairs
        .into_iter()
        .map(|(key, value)| (key, value_of(value)))
        .collect()
}

/// One info annotation's `annotate(ref, vc, likelihoods)`.
fn annotate_info(
    entry: &Entry,
    vc: &VariantContext,
    site: &Site<'_>,
    random: &mut JavaRandom,
) -> Result<Vec<(String, Value)>, EngineError> {
    use gatk_annotation as a;
    let likelihoods = Some(site.likelihoods);
    let from_trait = |annotation: &dyn InfoFieldAnnotation| {
        to_values(annotation.annotate(None, vc, likelihoods))
    };
    let string = |key: &str, value: Option<String>| -> Vec<(String, Value)> {
        value
            .map(|value| vec![(key.to_string(), Value::Str(value))])
            .unwrap_or_default()
    };
    Ok(match entry.name {
        "BaseQuality" => per_allele(&a::per_allele::BaseQuality, vc, site)?,
        "BaseQualityHistogram" => from_trait(&a::read_grouping::BaseQualityHistogram),
        "BaseQualityRankSumTest" => from_trait(&a::rank_sum::BaseQualityRankSumTest),
        "ChromosomeCounts" => from_trait(&a::chromosome_counts::ChromosomeCounts),
        "ClippingRankSumTest" => from_trait(&a::rank_sum::ClippingRankSumTest),
        "CountNs" => from_trait(&a::coverage::CountNs),
        "Coverage" => from_trait(&a::coverage::Coverage),
        "ExcessHet" => {
            if !vc.genotypes.is_empty() && vc.is_variant() {
                no_likelihoods(vc)?;
            }
            string("ExcessHet", a::heterozygosity::excess_het(vc))
        }
        "FisherStrand" => from_trait(&a::strand_bias::FisherStrand),
        "FragmentLength" => per_allele(&a::per_allele::FragmentLength, vc, site)?,
        "GenotypeSummaries" => from_trait(&a::site_statistics::GenotypeSummaries),
        "InbreedingCoeff" => {
            if vc.genotypes.len() >= a::heterozygosity::INBREEDING_MIN_SAMPLES && vc.is_variant() {
                no_likelihoods(vc)?;
            }
            string("InbreedingCoeff", a::heterozygosity::inbreeding_coeff(vc))
        }
        "LikelihoodRankSumTest" => to_values(a::rank_sum::annotate(
            &a::site_statistics::LikelihoodRankSumTest,
            None,
            vc,
            likelihoods,
        )),
        "MappingQuality" => per_allele(&a::per_allele::MappingQuality, vc, site)?,
        "MappingQualityRankSumTest" => from_trait(&a::rank_sum::MappingQualityRankSumTest),
        "MappingQualityZero" => from_trait(&a::coverage::MappingQualityZero),
        "OriginalAlignment" => from_trait(&a::original_alignment::OriginalAlignment),
        "QualByDepth" => {
            let raw =
                attribute(vc, "QUALapprox").map(|_| attribute_int(vc, "QUALapprox", 0) as i32);
            string(
                "QD",
                a::site_statistics::qual_by_depth(vc, likelihoods, raw, random),
            )
        }
        "RMSMappingQuality" => from_trait(&a::mapping_quality::RmsMappingQuality),
        "RawGtCount" => from_trait(&a::raw_gt_count::RawGtCount),
        "ReadPosRankSumTest" => from_trait(&a::rank_sum::ReadPosRankSumTest),
        "ReadPosition" => per_allele(&a::per_allele::ReadPosition, vc, site)?,
        "ReferenceBases" => vec![(
            "REF_BASES".to_string(),
            Value::Str(a::read_grouping::ReferenceBases::local_bases(
                site.window.0,
                site.window.1,
                vc,
            )),
        )],
        "SampleList" => from_trait(&a::sample_list::SampleList),
        "StrandOddsRatio" => from_trait(&a::strand_bias::StrandOddsRatio),
        "TandemRepeat" => {
            let (start, bases) = site.window;
            // `Arrays.copyOfRange(refBases, startIndex, refBases.length)`, which refuses a start
            // past the end: an indel annotated with no reference behind it.
            let start_index = vc.start + 1 - start;
            if a::tandem_repeat::is_indel(vc) && start_index > bases.len() as i64 {
                return Err(runtime_iae(format!("{start_index} > {}", bases.len())));
            }
            to_values(a::tandem_repeat::TandemRepeat::local_annotate(
                start, bases, vc,
            ))
        }
        "UniqueAltReadCount" => from_trait(&a::read_grouping::UniqueAltReadCount),
        "AS_FisherStrand" => to_values(a::allele_specific_strand_bias::annotate_direct(
            a::allele_specific_strand_bias::AsStrandBias::Fisher,
            vc,
            likelihoods,
        )),
        "AS_StrandOddsRatio" => to_values(a::allele_specific_strand_bias::annotate_direct(
            a::allele_specific_strand_bias::AsStrandBias::OddsRatio,
            vc,
            likelihoods,
        )),
        "AS_QualByDepth" => crate::genotype_gvcfs::finalize_as_qual_by_depth(vc, vc, random)?
            .into_iter()
            .filter(|(key, _)| key == "AS_QD")
            .collect(),
        "AS_RMSMappingQuality" => vec![(
            a::allele_specific_site_statistics::AS_RMS_MAPPING_QUALITY_KEY.to_string(),
            Value::Str(as_rms_mapping_quality(
                vc,
                &a::allele_specific_site_statistics::as_rms_data(site.likelihoods),
            )?),
        )],
        "AS_InbreedingCoeff" => string(
            "AS_InbreedingCoeff",
            a::allele_specific_site_statistics::as_inbreeding_coefficient(vc),
        ),
        name => {
            return Err(limitation(format!(
                "the {name} annotation is not ported for VariantAnnotator yet."
            )))
        }
    })
}

/// `GenotypeUtils.computeDiploidGenotypeCounts`' refusal of a called genotype with a GQ and no
/// likelihoods, which the ported counts skip.
fn no_likelihoods(vc: &VariantContext) -> Result<(), EngineError> {
    let genotypes: Vec<&htsjdk_vcf::variant::Genotype> = vc.genotypes.iter().collect();
    match gatk_annotation::heterozygosity::genotype_without_likelihoods(&genotypes) {
        Some(genotype) => Err(EngineError::Runtime {
            class: "java.lang.IllegalStateException".to_string(),
            message: format!(
                "Genotype has no likelihoods: {}",
                java_genotype_string(genotype)
            ),
        }),
        None => Ok(()),
    }
}

/// `Genotype.toString()`: the sample, the alleles (the reference starred, sorted unless phased),
/// the four inline fields that are present and the extended ones sorted by key.
pub fn java_genotype_string(genotype: &htsjdk_vcf::variant::Genotype) -> String {
    fn java(value: &Value) -> String {
        match value {
            Value::Missing => "null".to_string(),
            Value::Bool(flag) => flag.to_string(),
            Value::Str(text) => text.clone(),
            Value::List(items) => format!(
                "[{}]",
                items.iter().map(java).collect::<Vec<String>>().join(", ")
            ),
            other => other.format().unwrap_or_default(),
        }
    }
    let allele = |a: &Allele| {
        if a.is_no_call() {
            ".".to_string()
        } else if a.is_reference() {
            format!("{}*", a.display_string())
        } else {
            a.display_string()
        }
    };
    let calls = if genotype.alleles.is_empty() {
        "NA".to_string()
    } else {
        let mut alleles: Vec<&Allele> = genotype.alleles.iter().collect();
        if !genotype.phased {
            alleles.sort_by(|a, b| {
                b.is_reference()
                    .cmp(&a.is_reference())
                    .then(a.display_string().cmp(&b.display_string()))
            });
        }
        alleles
            .into_iter()
            .map(allele)
            .collect::<Vec<String>>()
            .join(if genotype.phased { "|" } else { "/" })
    };
    let int =
        |name: &str, value: Option<i32>| value.map(|v| format!(" {name} {v}")).unwrap_or_default();
    let ints = |name: &str, values: &Option<Vec<i32>>| {
        values
            .as_ref()
            .map(|v| {
                format!(
                    " {name} {}",
                    v.iter()
                        .map(i32::to_string)
                        .collect::<Vec<String>>()
                        .join(",")
                )
            })
            .unwrap_or_default()
    };
    let mut extended: Vec<(String, String)> = genotype
        .extended
        .iter()
        .map(|(key, value)| (key.clone(), java(value)))
        .collect();
    extended.sort();
    let extended = if extended.is_empty() {
        String::new()
    } else {
        format!(
            " {{{}}}",
            extended
                .iter()
                .map(|(key, value)| format!("{key}={value}"))
                .collect::<Vec<String>>()
                .join(", ")
        )
    };
    format!(
        "[{} {}{}{}{}{}{}{}]",
        genotype.sample_name,
        calls,
        int("GQ", genotype.gq),
        int("DP", genotype.dp),
        ints("AD", &genotype.ad),
        ints("PL", &genotype.pl),
        genotype
            .filters
            .as_ref()
            .map(|f| format!(" FT {f}"))
            .unwrap_or_default(),
        extended
    )
}

/// `AS_RMSMappingQuality.makeFinalizedAnnotationString`, as `annotate` reaches it with the
/// per-allele sums of squares the likelihoods give.
///
/// The separator is written BEFORE the allele is looked up, so an alternate no informative read
/// was assigned to leaves a trailing comma rather than nothing, and the depth divided by is the
/// genotypes' own `AD`, whatever the reads said.
fn as_rms_mapping_quality(
    vc: &VariantContext,
    per_allele: &[(Allele, Option<f64>)],
) -> Result<String, EngineError> {
    let counts = gatk_annotation::allele_specific_site_statistics::ad_counts(vc);
    let mut out = String::new();
    for alternate in &vc.alleles[1..] {
        if !out.is_empty() {
            out.push(',');
        }
        match per_allele.iter().find(|(allele, _)| allele == alternate) {
            None => {}
            Some((_, None)) => out.push('.'),
            Some((_, Some(square_sum))) => {
                let counts = counts.as_ref().ok_or_else(|| EngineError::Runtime {
                    class: "java.lang.NullPointerException".to_string(),
                    message:
                        "Cannot invoke \"java.util.Map.get(Object)\" because \"variantADs\" is null"
                            .to_string(),
                })?;
                let depth = counts
                    .iter()
                    .find(|(allele, _)| allele == alternate)
                    .map_or(0, |(_, count)| *count);
                out.push_str(&gatk_engine::java_format::format_decimals(
                    (square_sum / f64::from(depth)).sqrt(),
                    2,
                ));
            }
        }
    }
    Ok(out)
}

fn per_allele<A: gatk_annotation::per_allele::PerAlleleAnnotation>(
    annotation: &A,
    vc: &VariantContext,
    site: &Site<'_>,
) -> Result<Vec<(String, Value)>, EngineError> {
    gatk_annotation::per_allele::annotate(annotation, None, vc, Some(site.likelihoods))
        .map(to_values)
        .map_err(|error| EngineError::Runtime {
            class: error.class().to_string(),
            message: error.message().to_string(),
        })
}

/// One genotype annotation's `annotate(ref, vc, genotype, builder, likelihoods)`.
fn annotate_genotype(
    entry: &Entry,
    vc: &VariantContext,
    genotype: &htsjdk_vcf::variant::Genotype,
    built: &mut htsjdk_vcf::variant::Genotype,
    called: bool,
    likelihoods: &AlleleLikelihoods<BamRecord>,
) -> Result<(), EngineError> {
    use gatk_annotation as a;
    let sample = genotype.sample_name.as_str();
    match entry.name {
        // An `AF` the genotype already carries is kept: the annotation never overwrites one.
        "AlleleFraction" if genotype.extended.iter().any(|(key, _)| key == "AF") => {}
        "AlleleFraction" => {
            if let Some(fractions) = a::depth_per_allele::allele_fractions(
                vc,
                genotype.ad.as_deref(),
                Some(likelihoods),
                sample,
                called,
            ) {
                set_extended(
                    built,
                    "AF",
                    Value::List(fractions.into_iter().map(Value::Double).collect()),
                );
            }
        }
        "DepthPerAlleleBySample" => {
            if let Some(depths) =
                a::depth_per_allele::allele_depths(vc, Some(likelihoods), sample, called)
            {
                built.ad = Some(depths);
            }
        }
        "DepthPerSampleHC" => {
            if let Some(depth) =
                a::depth_per_allele::informative_depth(Some(likelihoods), sample, called)
            {
                built.dp = Some(depth);
            }
        }
        // A genotype that is not called is left alone, likelihoods or not.
        "StrandBiasBySample" if !called => {}
        "StrandBiasBySample" => {
            if let Some(counts) =
                a::strand_bias::StrandBiasBySample.counts(vc, sample, Some(likelihoods))
            {
                set_extended(
                    built,
                    "SB",
                    Value::List(
                        counts
                            .into_iter()
                            .map(|count| Value::Int(count as i64))
                            .collect(),
                    ),
                );
            }
        }
        name => {
            return Err(limitation(format!(
                "the {name} annotation is not ported for VariantAnnotator yet."
            )))
        }
    }
    Ok(())
}

fn set_extended(genotype: &mut htsjdk_vcf::variant::Genotype, key: &str, value: Value) {
    match genotype.extended.iter_mut().find(|(name, _)| name == key) {
        Some(slot) => slot.1 = value,
        None => genotype.extended.push((key.to_string(), value)),
    }
}

/// Why `makeLikelihoods` could not build the matrix.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LikelihoodsFailure {
    /// `splitReadsBySample` found a read whose sample is not one of the VCF's, which with the
    /// startup check in place means a read with no read group: `returnMap.get(null).add(read)`.
    NoSampleSlot,
    /// A negative `--min-base-quality-score`.
    NegativeMinimumBaseQuality,
}

impl LikelihoodsFailure {
    pub fn class(&self) -> &'static str {
        match self {
            LikelihoodsFailure::NoSampleSlot => "java.lang.NullPointerException",
            LikelihoodsFailure::NegativeMinimumBaseQuality => "java.lang.IllegalArgumentException",
        }
    }

    pub fn message(&self) -> &'static str {
        match self {
            LikelihoodsFailure::NoSampleSlot => {
                "Cannot invoke \"java.util.List.add(Object)\" because the return value of \"java.util.Map.get(Object)\" is null"
            }
            LikelihoodsFailure::NegativeMinimumBaseQuality => {
                "minBaseQualityCutoff must be greater than or equal to 0"
            }
        }
    }
}

/// `VariantAnnotator.makeLikelihoods(vc, readsContext)`.
///
/// `reads` are the reads the walker's `ReadsContext` holds for the record, after the read filters,
/// in file order.
pub fn make_likelihoods(
    vc: &VariantContext,
    samples: &[String],
    reads: &[BamRecord],
    header: Option<&SamHeader>,
    minimum_base_quality: i32,
) -> Result<AlleleLikelihoods<BamRecord>, LikelihoodsFailure> {
    let sample_of = |read: &BamRecord| -> Option<String> {
        header.and_then(|header| gatk_engine::read_pileup::sample_name(read, header))
    };
    // `splitReadsBySample(variantSamples, header, reads)`.
    let mut evidence: Vec<Vec<BamRecord>> = vec![Vec::new(); samples.len()];
    for read in reads {
        let index = sample_of(read)
            .and_then(|sample| samples.iter().position(|name| *name == sample))
            .ok_or(LikelihoodsFailure::NoSampleSlot)?;
        evidence[index].push(read.clone());
    }
    let allele_count = vc.alleles.len();
    let mut values: Vec<Vec<Vec<f64>>> = evidence
        .iter()
        .map(|reads| vec![vec![0.0; reads.len()]; allele_count])
        .collect();

    // `new ReadPileup(vc, readsContext)`, then `splitBySample(header, "__UNKNOWN__")`.
    let pileup = gatk_engine::read_pileup::pileup_from_reads(
        &vc.contig,
        vc.start as i32,
        reads,
        |read| read.flags & 0x200 == 0,
        |read| read.flags & 0x400 == 0,
    );
    let engine_allele = |allele: &Allele| {
        gatk_engine::variant_context_utils::Allele::new(
            allele.display_string().as_bytes(),
            allele.is_reference(),
        )
    };
    let reference = engine_allele(&vc.alleles[0]);
    let alternates: Vec<_> = vc.alleles[1..].iter().map(engine_allele).collect();
    let mut next_row: Vec<usize> = vec![0; samples.len()];
    for element in &pileup.elements {
        let Some(sample) = samples
            .iter()
            .position(|name| Some(name.clone()) == sample_of(element.read))
        else {
            // Unreachable: `splitReadsBySample` already refused such a read.
            return Err(LikelihoodsFailure::NoSampleSlot);
        };
        let row = next_row[sample];
        next_row[sample] += 1;
        for allele in values[sample].iter_mut() {
            if let Some(slot) = allele.get_mut(row) {
                *slot = f64::NEG_INFINITY;
            }
        }
        let chosen = gatk_engine::variant_context_utils::choose_allele_for_read(
            element,
            &reference,
            &alternates,
            minimum_base_quality,
        )
        .map_err(|_| LikelihoodsFailure::NegativeMinimumBaseQuality)?;
        if let Some(chosen) = chosen {
            let index = if chosen == reference {
                0
            } else {
                1 + alternates
                    .iter()
                    .position(|alternate| *alternate == chosen)
                    .unwrap_or(0)
            };
            if let Some(slot) = values[sample][index].get_mut(row) {
                *slot = 0.0;
            }
        }
    }
    Ok(AlleleLikelihoods::new(
        SampleList::new(samples),
        AlleleList::new(&vc.alleles),
        evidence,
        values,
    )
    .expect("every row is sized from its sample's evidence"))
}

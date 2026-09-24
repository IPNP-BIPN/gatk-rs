//! `GenotypeGVCFs`: how a GVCF becomes a genotyped VCF.
//!
//! Ported from
//! `org.broadinstitute.hellbender.tools.walkers.GenotypeGVCFs`,
//! `org.broadinstitute.hellbender.tools.walkers.GenotypeGVCFsEngine` and
//! `org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine`
//! (GATK 4.6.2.0), for plain VCF or GVCF input. GenomicsDB and the somatic path are not ported.
//!
//! Each locus is merged by [`crate::reference_confidence_merger`] with `<NON_REF>` removed, the
//! merged record is re-genotyped by [`crate::genotyping_engine`], and the annotations are
//! finalised and then computed again from the called genotypes.
//!
//! # A reference block is written only when output is forced
//!
//! By default the walker visits one record at a time and a record whose only alternate was
//! `<NON_REF>` merges to a record with no alternate at all: it is not a variant, is not
//! re-genotyped, and is dropped. `--include-non-variant-sites` walks locus by locus instead, over
//! every base a record covers, and writes each such base as `0/0` with the GQ moved to `RGQ`.
//!
//! # The annotations are computed twice
//!
//! `finalizeAnnotations` turns the raw data the GVCF carried into final values (`RAW_MQandDP` into
//! `MQ`) and removes the raw keys unless `--keep-combined-raw-annotations` asked to keep them.
//! Then `annotateContext` runs every resolved info annotation over the re-genotyped record, with no
//! reads: `AC`, `AF` and `AN` describe the output, `ExcessHet` and `InbreedingCoeff` the called
//! genotypes, `FS` and `SOR` the samples' `SB`, and `QD` the new QUAL over the variant samples'
//! depth. Annotations that need reads write nothing.
//!
//! # QD above 35 consumes the run's one generator
//!
//! `QD` past 35 is replaced by a Gaussian draw from `Utils.getRandomGenerator()`, one stream per
//! run. [`Engine::call_region`] therefore takes the generator from its caller, so that the value a
//! site gets depends on how many draws the sites before it took, as it does in the reference.

use gatk_annotation::catalogue::{Entry, Kind};
use gatk_annotation::info_annotation::{AnnotationValue, InfoFieldAnnotation};
use gatk_engine::java_random::JavaRandom;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::genotypes_context::GenotypesContext;
use htsjdk_vcf::header::{Cardinality, HeaderLine, VcfHeader};
use htsjdk_vcf::variant::{Genotype, Value, VariantContext};

use crate::genotyping_engine::{
    EngineError, GenotypingEngine, AS_QUAL_KEY, MLE_ALLELE_COUNT_KEY, MLE_ALLELE_FREQUENCY_KEY,
    NUMBER_OF_DISCOVERED_ALLELES_KEY,
};
use crate::reference_confidence_merger::{MergeNotes, Merger};

/// `GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY`.
pub const REFERENCE_GENOTYPE_QUALITY: &str = "RGQ";
/// `GATKVCFConstants.MIN_DP_FORMAT_KEY`.
const MIN_DP_FORMAT_KEY: &str = "MIN_DP";
/// `GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY`.
const STRAND_BIAS_BY_SAMPLE_KEY: &str = "SB";
/// `GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY`.
const HAPLOTYPE_CALLER_PHASING_GT_KEY: &str = "PGT";
/// `GenotypeGVCFs.PHASED_HOM_VAR_STRING`.
const PHASED_HOM_VAR_STRING: &str = "1|1";
/// The raw keys `finalizeAnnotations` removes by hand: `AS_QUAL`, `QUALapprox`, `VarDP` and
/// `RAW_GT_COUNT`.
const MANUAL_RAW_KEYS: [&str; 4] = [AS_QUAL_KEY, "QUALapprox", "VarDP", "RAW_GT_COUNT"];

/// `VariantAnnotatorEngine`, over the resolved annotations, with no reads and no pedigree.
pub struct AnnotationEngine {
    /// The resolved annotations, sorted by simple name as the plugin descriptor hands them out.
    pub resolved: Vec<&'static Entry>,
    /// `--keep-combined-raw-annotations`.
    pub keep_combined: bool,
    /// `rawAnnotationsToKeep`: the raw keys of the annotations `--keep-specific-combined-raw-
    /// annotation` named.
    pub raw_keys_to_keep: Vec<&'static str>,
}

fn value_of(annotation: AnnotationValue) -> Value {
    match annotation {
        AnnotationValue::Int(value) => Value::Int(value as i64),
        AnnotationValue::Long(value) => Value::Int(value),
        AnnotationValue::Double(value) => Value::Double(value),
        AnnotationValue::Str(value) => Value::Str(value),
        AnnotationValue::Flag(value) => Value::Bool(value),
        AnnotationValue::List(values) => Value::List(values.into_iter().map(value_of).collect()),
    }
}

/// `LinkedHashMap.put`: replace in place, or append.
fn put(attributes: &mut Vec<(String, Value)>, key: &str, value: Value) {
    match attributes.iter_mut().find(|(name, _)| name == key) {
        Some(slot) => slot.1 = value,
        None => attributes.push((key.to_string(), value)),
    }
}

fn attribute<'a>(vc: &'a VariantContext, key: &str) -> Option<&'a Value> {
    vc.attributes
        .iter()
        .find(|(name, _)| name == key)
        .map(|(_, value)| value)
}

/// `getAttributeAsString`: a list is the Java list's string, a scalar its own.
fn attribute_string(vc: &VariantContext, key: &str) -> Option<String> {
    attribute(vc, key).map(crate::reference_confidence_merger::value_to_string)
}

/// `getAttributeAsInt(key, default)`.
fn attribute_int(vc: &VariantContext, key: &str, default: i64) -> i64 {
    match attribute(vc, key) {
        Some(Value::Int(value)) => *value,
        Some(Value::Double(value)) => *value as i64,
        Some(Value::Str(text)) => text.trim().parse().unwrap_or(default),
        Some(Value::List(values)) if values.len() == 1 => match &values[0] {
            Value::Int(value) => *value,
            Value::Str(text) => text.trim().parse().unwrap_or(default),
            _ => default,
        },
        _ => default,
    }
}

impl AnnotationEngine {
    /// The resolved info annotations, in order.
    fn info_annotations(&self) -> impl Iterator<Item = &&'static Entry> {
        self.resolved
            .iter()
            .filter(|entry| entry.kind == Kind::Info)
    }

    /// Whether any resolved info annotation is allele-specific, which is what switches on
    /// `AS_QUAL` in the genotyping engine.
    pub fn any_allele_specific(&self) -> bool {
        self.info_annotations()
            .any(|entry| entry.name.starts_with("AS_"))
    }

    /// `finalizeAnnotations(vc, originalVC)`.
    pub fn finalize_annotations(&self, vc: &VariantContext) -> Result<VariantContext, EngineError> {
        let mut attributes = vc.attributes.clone();
        let saved: Vec<(String, Value)> = self
            .raw_keys_to_keep
            .iter()
            .filter_map(|key| attribute(vc, key).map(|value| ((*key).to_string(), value.clone())))
            .collect();
        for entry in self.info_annotations().filter(|entry| entry.is_reducible()) {
            for (key, value) in finalize_raw_data(entry, vc)? {
                put(&mut attributes, &key, value);
            }
            if !self.keep_combined {
                let raw_keys = entry.raw_keys.unwrap_or(&[]);
                attributes.retain(|(key, _)| !raw_keys.contains(&key.as_str()));
            }
        }
        if !self.keep_combined {
            attributes.retain(|(key, _)| !MANUAL_RAW_KEYS.contains(&key.as_str()));
        }
        for (key, value) in saved {
            put(&mut attributes, &key, value);
        }
        let mut out = vc.clone();
        out.attributes = attributes;
        Ok(out)
    }

    /// `annotateContext(vc, features, ref, null, addAnnot)`, with the genotype annotations a
    /// no-op because every one the catalogue carries needs reads.
    ///
    /// `hom_ref_site` is the monomorphic branch, which passes
    /// `GenotypeGVCFsEngine::annotationShouldBeSkippedForHomRefSites` as the predicate. The
    /// predicate says which annotations to ADD, whatever its name says, so at a hom-ref site only
    /// the rank sums and the two MQs run, and with no reads they write nothing: that is why a
    /// forced reference site carries no `AN`.
    pub fn annotate_context(
        &self,
        vc: &VariantContext,
        hom_ref_site: bool,
        random: &mut JavaRandom,
    ) -> Result<VariantContext, EngineError> {
        let mut attributes = vc.attributes.clone();
        for entry in self.info_annotations() {
            if hom_ref_site && !annotation_should_be_skipped_for_hom_ref_sites(entry) {
                continue;
            }
            for (key, value) in annotate(entry, vc, random)? {
                put(&mut attributes, &key, value);
            }
        }
        let mut out = vc.clone();
        out.attributes = attributes;
        Ok(out)
    }
}

/// `annotationShouldBeSkippedForHomRefSites`: the rank sums and the two MQs.
fn annotation_should_be_skipped_for_hom_ref_sites(entry: &Entry) -> bool {
    matches!(
        entry.name,
        "BaseQualityRankSumTest"
            | "MappingQualityRankSumTest"
            | "ReadPosRankSumTest"
            | "ClippingRankSumTest"
            | "LikelihoodRankSumTest"
            | "RMSMappingQuality"
            | "AS_RMSMappingQuality"
    )
}

fn limitation(what: String) -> EngineError {
    EngineError::Limitation(what)
}

/// One reducible annotation's `finalizeRawData(vc, originalVC)`.
fn finalize_raw_data(
    entry: &Entry,
    vc: &VariantContext,
) -> Result<Vec<(String, Value)>, EngineError> {
    match entry.name {
        "RMSMappingQuality" => {
            if let Some(raw) = attribute_string(vc, "RAW_MQandDP") {
                let finalized =
                    gatk_annotation::mapping_quality::RmsMappingQuality::finalize_raw_data(&raw)
                        .map_err(|error| EngineError::Runtime {
                            class:
                                "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
                                    .to_string(),
                            message: format!("{error:?}"),
                        })?;
                Ok(vec![("MQ".to_string(), Value::Str(finalized))])
            } else if attribute(vc, "RAW_MQ").is_some() {
                Err(limitation(
                    "the deprecated RAW_MQ key is refused or reformatted by GATK, which this port \
                     does not carry yet."
                        .to_string(),
                ))
            } else {
                Ok(Vec::new())
            }
        }
        name => {
            // Every other reducible annotation writes nothing unless its raw key is there.
            let raw_keys = entry.raw_keys.unwrap_or(&[]);
            if raw_keys.iter().any(|key| attribute(vc, key).is_some()) {
                Err(limitation(format!(
                    "finalizing {name}'s raw data is not ported yet."
                )))
            } else {
                Ok(Vec::new())
            }
        }
    }
}

/// One info annotation's `annotate(ref, vc, null)`.
fn annotate(
    entry: &Entry,
    vc: &VariantContext,
    random: &mut JavaRandom,
) -> Result<Vec<(String, Value)>, EngineError> {
    use gatk_annotation as a;
    let from_trait = |annotation: &dyn InfoFieldAnnotation| -> Vec<(String, Value)> {
        annotation
            .annotate(None, vc, None)
            .into_iter()
            .map(|(key, value)| (key, value_of(value)))
            .collect()
    };
    Ok(match entry.name {
        "ChromosomeCounts" => from_trait(&a::chromosome_counts::ChromosomeCounts),
        "Coverage" => from_trait(&a::coverage::Coverage),
        "FisherStrand" => from_trait(&a::strand_bias::FisherStrand),
        "StrandOddsRatio" => from_trait(&a::strand_bias::StrandOddsRatio),
        "RMSMappingQuality" => from_trait(&a::mapping_quality::RmsMappingQuality),
        "BaseQualityRankSumTest" => from_trait(&a::rank_sum::BaseQualityRankSumTest),
        "MappingQualityRankSumTest" => from_trait(&a::rank_sum::MappingQualityRankSumTest),
        "ReadPosRankSumTest" => from_trait(&a::rank_sum::ReadPosRankSumTest),
        "ExcessHet" => a::heterozygosity::excess_het(vc)
            .map(|value| vec![("ExcessHet".to_string(), Value::Str(value))])
            .unwrap_or_default(),
        "InbreedingCoeff" => a::heterozygosity::inbreeding_coeff(vc)
            .map(|value| vec![("InbreedingCoeff".to_string(), Value::Str(value))])
            .unwrap_or_default(),
        "QualByDepth" => {
            let raw =
                attribute(vc, "QUALapprox").map(|_| attribute_int(vc, "QUALapprox", 0) as i32);
            a::site_statistics::qual_by_depth(vc, None, raw, random)
                .map(|value| vec![("QD".to_string(), Value::Str(value))])
                .unwrap_or_default()
        }
        name => {
            return Err(limitation(format!(
                "the {name} annotation is not ported for GenotypeGVCFs yet."
            )))
        }
    })
}

/// The walker's arguments that reach the engine.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Arguments {
    /// `--include-non-variant-sites`.
    pub include_non_variants: bool,
    /// Whether `StrandBiasBySample` was resolved, which keeps `SB` in the output.
    pub keep_sb: bool,
}

/// `GenotypeGVCFsEngine`.
pub struct Engine<'a> {
    pub merger: Merger<'a>,
    pub genotyping: GenotypingEngine,
    pub forced: GenotypingEngine,
    pub annotations: AnnotationEngine,
    pub arguments: Arguments,
    /// The output header's FORMAT lines, which type the allele-specific genotype fields.
    pub output_header: VcfHeader,
    pub notes: MergeNotes,
}

/// `getVariantSubsetToProcess`: with non-variants included, a record starting at the locus wins.
fn variant_subset_to_process(
    start: i64,
    variants: &[VariantContext],
    include_non_variants: bool,
) -> Result<Vec<VariantContext>, EngineError> {
    if !include_non_variants {
        return Ok(variants.to_vec());
    }
    let matching: Vec<&VariantContext> = variants.iter().filter(|vc| vc.start == start).collect();
    match matching.len() {
        0 => Ok(variants.to_vec()),
        1 => Ok(vec![matching[0].clone()]),
        _ => Err(EngineError::Runtime {
            class: "java.lang.IllegalStateException".to_string(),
            message: format!(
                "Variant input contains more than one variant starting at location: {}:{}-{}",
                matching[0].contig, matching[0].start, matching[0].stop
            ),
        }),
    }
}

/// `VariantContext.isPolymorphicInSamples`.
fn is_polymorphic_in_samples(vc: &VariantContext) -> bool {
    htsjdk_vcf::genotype_type::is_polymorphic_in_samples(vc)
}

/// `GATKVariantContextUtils.isProperlyPolymorphic`.
fn is_properly_polymorphic(vc: &VariantContext) -> bool {
    let alternates = &vc.alleles[1..];
    if alternates.is_empty() {
        return false;
    }
    let span_del = |allele: &Allele| allele.display_string() == htsjdk_vcf::allele::SPAN_DEL_STRING;
    if alternates.len() == 1 {
        return !(span_del(&alternates[0]) || alternates[0].is_symbolic());
    }
    !(span_del(&alternates[0]) && alternates[1] == crate::reference_confidence_merger::non_ref())
}

/// `excludeFromAnnotations`: a hom-ref or no-call with no depth and a GQ of 0.
fn exclude_from_annotations(genotype: &Genotype) -> bool {
    (crate::reference_confidence_merger::is_hom_ref(genotype)
        || crate::reference_confidence_merger::is_no_call(genotype))
        && genotype.dp.is_none_or(|dp| dp == 0)
        && genotype.gq == Some(0)
}

/// `assignNoCallsAnnotationExcludedGenotypes`.
fn assign_no_calls_annotation_excluded(genotypes: &[Genotype]) -> Vec<Genotype> {
    genotypes
        .iter()
        .map(|genotype| {
            let mut out = genotype.clone();
            if exclude_from_annotations(genotype) {
                out.alleles = vec![Allele::no_call(); genotype.ploidy()];
            }
            out
        })
        .collect()
}

/// `getAnyAttribute(MIN_DP)` read as an int, which is what `parseInt` accepts.
fn int_of(value: &Value) -> Result<i32, EngineError> {
    match value {
        Value::Int(value) => Ok(*value as i32),
        Value::Str(text) => text.trim().parse().map_err(|_| EngineError::Runtime {
            class: "java.lang.NumberFormatException".to_string(),
            message: format!("For input string: \"{text}\""),
        }),
        _ => Err(EngineError::Runtime {
            class: "java.lang.IllegalArgumentException".to_string(),
            message: "Expected a Number or a String but found something else.".to_string(),
        }),
    }
}

/// `cleanupGenotypeAnnotations(vc, createRefGTs, keepSB)`.
pub fn cleanup_genotype_annotations(
    vc: &VariantContext,
    create_ref_gts: bool,
    keep_sb: bool,
) -> Result<Vec<Genotype>, EngineError> {
    let mut recovered = Vec::with_capacity(vc.genotypes.len());
    for old in vc.genotypes.iter() {
        let mut attrs = old.extended.clone();
        let mut built = old.clone();
        let mut depth = old.dp.unwrap_or(0);
        if let Some(position) = attrs.iter().position(|(key, _)| key == MIN_DP_FORMAT_KEY) {
            depth = int_of(&attrs[position].1)?;
            built.dp = Some(depth);
            attrs.remove(position);
        }
        if !keep_sb {
            attrs.retain(|(key, _)| key != STRAND_BIAS_BY_SAMPLE_KEY);
        }
        if old.is_hom_var() {
            if let Some(slot) = attrs
                .iter_mut()
                .find(|(key, _)| key == HAPLOTYPE_CALLER_PHASING_GT_KEY)
            {
                slot.1 = Value::Str(PHASED_HOM_VAR_STRING.to_string());
            }
        }
        if old.ad.is_none() && vc.is_variant() && depth > 0 {
            let mut ad = vec![0; vc.alleles.len()];
            ad[0] = depth;
            built.ad = Some(ad);
        }
        if create_ref_gts {
            if depth > 0 && old.gq.is_some() {
                let gq = old.gq.unwrap_or(0);
                if gq > 0 {
                    built.alleles = vec![vc.reference().clone(); old.ploidy()];
                } else {
                    built.alleles = vec![Allele::no_call(); old.ploidy()];
                }
                built.gq = None;
                attrs.push((
                    REFERENCE_GENOTYPE_QUALITY.to_string(),
                    Value::Int(gq as i64),
                ));
            } else {
                built.alleles = vec![Allele::no_call(); old.ploidy()];
                built.gq = None;
                built.dp = None;
            }
            built.pl = None;
        }
        built.extended = attrs;
        recovered.push(built);
    }
    Ok(recovered)
}

/// `ReferenceConfidenceVariantContextMerger.generateAnnotationValueVector`.
fn annotation_value_vector(number: &Cardinality, values: &[Value], indices: &[usize]) -> Value {
    match number {
        Cardinality::A => Value::List(
            indices
                .iter()
                .skip(1)
                .filter_map(|index| index.checked_sub(1).and_then(|i| values.get(i).cloned()))
                .collect(),
        ),
        Cardinality::R => Value::List(
            indices
                .iter()
                .filter_map(|index| values.get(*index).cloned())
                .collect(),
        ),
        _ => Value::List(values.to_vec()),
    }
}

impl Engine<'_> {
    /// `callRegion(loc, variants, ref, ...)`: merge, then re-genotype.
    pub fn call_region(
        &mut self,
        contig: &str,
        start: i64,
        variants: &[VariantContext],
        reference_base: u8,
        force_output: bool,
        random: &mut JavaRandom,
    ) -> Result<Option<VariantContext>, EngineError> {
        let to_process =
            variant_subset_to_process(start, variants, self.arguments.include_non_variants)?;
        let merged = self
            .merger
            .merge(
                &to_process,
                contig,
                start,
                Some(reference_base),
                true,
                &mut self.notes,
            )
            .map_err(|error| EngineError::Runtime {
                class: error.java_class().to_string(),
                message: error.message(),
            })?
            .expect("a reference base was given");
        self.regenotype(&merged, force_output, random)
    }

    /// `subsetAlleleSpecificFormatFields`.
    fn subset_allele_specific_format_fields(
        &self,
        genotypes: &[Genotype],
        relevant: &[usize],
    ) -> Vec<Genotype> {
        genotypes
            .iter()
            .map(|genotype| {
                let mut out = genotype.clone();
                for (key, value) in out.extended.iter_mut() {
                    let number = self.output_header.lines.iter().find_map(|line| match line {
                        HeaderLine::Compound {
                            key: kind,
                            id,
                            number,
                            ..
                        } if kind == "FORMAT" && id == key => Some(*number),
                        _ => None,
                    });
                    match number {
                        Some(Cardinality::Fixed(1)) | None => {}
                        Some(number) => {
                            let values = match &*value {
                                Value::List(values) => values.clone(),
                                other => vec![other.clone()],
                            };
                            *value = annotation_value_vector(&number, &values, relevant);
                        }
                    }
                }
                out
            })
            .collect()
    }

    /// `regenotypeVC(originalVC, ref, features, includeNonVariants)`.
    fn regenotype(
        &mut self,
        original: &VariantContext,
        include_non_variants: bool,
        random: &mut JavaRandom,
    ) -> Result<Option<VariantContext>, EngineError> {
        let result = if original.is_variant() && attribute_int(original, "DP", 0) > 0 {
            let engine = if include_non_variants {
                &mut self.forced
            } else {
                &mut self.genotyping
            };
            let Some(regenotyped) = engine.calculate_genotypes(original)? else {
                return Ok(None);
            };
            if is_properly_polymorphic(&regenotyped) || include_non_variants {
                let with_genotyping =
                    add_genotyping_annotations(&original.attributes, &regenotyped);
                let with_annotations = self.annotations.finalize_annotations(&with_genotyping)?;
                let relevant: Vec<usize> = regenotyped
                    .alleles
                    .iter()
                    .map(|allele| {
                        original
                            .alleles
                            .iter()
                            .position(|a| a == allele)
                            .unwrap_or(usize::MAX)
                    })
                    .collect();
                let trimmed = crate::variant_trim::reverse_trim_alleles(&with_annotations)
                    .map_err(|message| EngineError::Runtime {
                        class: "java.lang.IllegalStateException".to_string(),
                        message,
                    })?;
                let genotypes =
                    self.subset_allele_specific_format_fields(&trimmed.genotypes, &relevant);
                let mut out = trimmed;
                out.genotypes = GenotypesContext::new(genotypes);
                out
            } else {
                return Ok(None);
            }
        } else {
            original.clone()
        };

        if is_polymorphic_in_samples(&result) && attribute_int(&result, "DP", 0) > 0 {
            let mut prepared = result.clone();
            prepared.genotypes =
                GenotypesContext::new(assign_no_calls_annotation_excluded(&result.genotypes));
            let annotated = self
                .annotations
                .annotate_context(&prepared, false, random)?;
            let cleaned = cleanup_genotype_annotations(&annotated, false, self.arguments.keep_sb)?;
            let mut out = annotated;
            out.genotypes = GenotypesContext::new(cleaned);
            Ok(Some(out))
        } else if include_non_variants {
            let mut preannotated = result.clone();
            preannotated.genotypes =
                GenotypesContext::new(cleanup_genotype_annotations(&result, true, false)?);
            Ok(Some(self.annotations.annotate_context(
                &preannotated,
                true,
                random,
            )?))
        } else {
            Ok(None)
        }
    }
}

/// `addGenotypingAnnotations`: the original record's attributes, and MLEAC, MLEAF, NDA and AS_QUAL
/// from the re-genotyped one. MLEAC and MLEAF are put even when the genotyper wrote none, which is
/// how a forced monomorphic site comes out with `MLEAC=.;MLEAF=.`.
fn add_genotyping_annotations(
    original: &[(String, Value)],
    regenotyped: &VariantContext,
) -> VariantContext {
    let mut attributes = original.to_vec();
    for key in [MLE_ALLELE_COUNT_KEY, MLE_ALLELE_FREQUENCY_KEY] {
        let value = attribute(regenotyped, key)
            .cloned()
            .unwrap_or(Value::Missing);
        put(&mut attributes, key, value);
    }
    for key in [NUMBER_OF_DISCOVERED_ALLELES_KEY, AS_QUAL_KEY] {
        if let Some(value) = attribute(regenotyped, key) {
            put(&mut attributes, key, value.clone());
        }
    }
    let mut out = regenotyped.clone();
    out.attributes = attributes;
    out
}

/// `GATKVariantContextUtils.isSpanningDeletionOnly`: one alternate, and it is `*`.
pub fn is_spanning_deletion_only(vc: &VariantContext) -> bool {
    vc.alleles.len() == 2 && vc.alleles[1].display_string() == htsjdk_vcf::allele::SPAN_DEL_STRING
}

/// A failure inside the traversal, with the record `VariantLocusWalker` names when it wraps an
/// `IllegalStateException`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WalkFailure {
    pub error: EngineError,
    /// The record the wrapper names: the visited one by variant, the first overlapping one by
    /// locus. `None` when the error is not one the wrapper catches.
    pub record: Option<(String, i64, bool)>,
}

/// `VariantLocusWalker.traverse` and `GenotypeGVCFs.apply` over records already restricted to the
/// traversal intervals, in file order.
///
/// `by_locus` is the traversal `--include-non-variant-sites` and `--force-output-intervals`
/// choose: every base some record covers and `traversed` accepts, with all the records covering
/// it. `forced` answers whether a locus lies in `--force-output-intervals`, and `reference_base`
/// the base the merger needs.
pub fn walk(
    engine: &mut Engine<'_>,
    records: &[VariantContext],
    by_locus: bool,
    traversed: &dyn Fn(&str, i64) -> bool,
    forced: &dyn Fn(&str, i64) -> bool,
    reference_base: &mut dyn FnMut(&str, i64) -> u8,
    random: &mut JavaRandom,
) -> Result<Vec<VariantContext>, WalkFailure> {
    let mut out = Vec::new();
    let mut emit = |engine: &mut Engine<'_>,
                    contig: &str,
                    start: i64,
                    overlapping: &[VariantContext],
                    random: &mut JavaRandom,
                    out: &mut Vec<VariantContext>|
     -> Result<(), WalkFailure> {
        let force_output = engine.arguments.include_non_variants || forced(contig, start);
        let base = reference_base(contig, start);
        let called = engine
            .call_region(contig, start, overlapping, base, force_output, random)
            .map_err(|error| {
                let wrapped = matches!(
                    &error,
                    EngineError::Runtime { class, .. } if class == "java.lang.IllegalStateException"
                );
                WalkFailure {
                    record: wrapped.then(|| {
                        (
                            overlapping[0].contig.clone(),
                            overlapping[0].start,
                            by_locus,
                        )
                    }),
                    error,
                }
            })?;
        if let Some(vc) = called {
            if force_output || !is_spanning_deletion_only(&vc) {
                out.push(vc);
            }
        }
        Ok(())
    };
    if !by_locus {
        for record in records {
            emit(
                engine,
                &record.contig,
                record.start,
                std::slice::from_ref(record),
                random,
                &mut out,
            )?;
        }
        return Ok(out);
    }
    // Every base a record covers, contig by contig in the order the records come, with the
    // records covering it in file order.
    let mut contigs: Vec<&str> = Vec::new();
    for record in records {
        if !contigs.contains(&record.contig.as_str()) {
            contigs.push(&record.contig);
        }
    }
    for contig in contigs {
        let on_contig: Vec<&VariantContext> =
            records.iter().filter(|r| r.contig == contig).collect();
        let mut loci: Vec<i64> = on_contig
            .iter()
            .flat_map(|record| record.start..=record.stop.max(record.start))
            .collect();
        loci.sort_unstable();
        loci.dedup();
        for locus in loci.into_iter().filter(|locus| traversed(contig, *locus)) {
            let overlapping: Vec<VariantContext> = on_contig
                .iter()
                .filter(|record| record.start <= locus && locus <= record.stop.max(record.start))
                .map(|record| (*record).clone())
                .collect();
            emit(engine, contig, locus, &overlapping, random, &mut out)?;
        }
    }
    Ok(out)
}

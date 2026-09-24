//! `ReblockGVCF`: a single-sample GVCF's reference blocks re-banded, its weak variants demoted to
//! blocks, and its strong ones trimmed to the alleles they call.
//!
//! Ported from `org.broadinstitute.hellbender.tools.walkers.variantutils.ReblockGVCF` (GATK
//! 4.6.2.0). The writing half is [`crate::gvcf_blocks::ReblockingWriter`]; this is `apply` and
//! everything it calls.
//!
//! # What reaches the writer
//!
//! * A reference block goes to the writer as it is, to be merged into its band there.
//! * A variant whose best genotype is the reference, or calls `<NON_REF>`, or whose `PL[0]` is
//!   below `--rgq-threshold-to-no-call`, is DEMOTED: a GQ0 hom-ref over `<NON_REF>` with
//!   `PL=0,0,0` (or its PLs against the best alternate, for a hom-ref call), and `END` its end.
//! * Anything else keeps the alleles its genotype calls plus `<NON_REF>`, its genotype
//!   re-subset to them, and a fresh set of INFO keys: the engine's annotations copied from the
//!   input, `RAW_MQandDP` rebuilt, `RAW_GT_COUNT`, and under `--do-qual-score-approximation`
//!   `QUALapprox`, `VarDP` and their allele-specific forms.
//!
//! # Two places the reference mutates what it read
//!
//! `removeNonRefADs` zeroes the `<NON_REF>` depth in the genotype's OWN array, so the depths the
//! qual annotations read afterwards from the pre-fix record already have it zeroed. And
//! `addRefBlockIfNecessary`'s "renormalisation" of the new block's PLs zeroes the first entry and
//! then subtracts that zero from the rest, which leaves them as they were.

use gatk_annotation::catalogue::{Entry, Kind};
use gatk_engine::genotype_index::{genotypes_in_canonical_order, subsetted_pl_indices};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::genotypes_context::GenotypesContext;
use htsjdk_vcf::variant::{Genotype, Value, VariantContext, NO_LOG10_PERROR};

use crate::genotyping_engine::{
    make_genotype_call, most_likely_alleles_ensuring_alt, subset_alleles,
    subset_alleles_with_length_annotations, EngineError, GenotypingEngine, SubsetMethod,
};
use crate::gvcf_blocks::{move_start, ReblockingWriter, ReferenceBase};
use crate::reference_confidence_merger::{java_round, value_to_string};

pub const RAW_QUAL_APPROX_KEY: &str = "QUALapprox";
pub const AS_RAW_QUAL_APPROX_KEY: &str = "AS_QUALapprox";
pub const VARIANT_DEPTH_KEY: &str = "VarDP";
pub const AS_VARIANT_DEPTH_KEY: &str = "AS_VarDP";
pub const RAW_GENOTYPE_COUNT_KEY: &str = "RAW_GT_COUNT";
pub const RAW_MAPPING_QUALITY_WITH_DEPTH_KEY: &str = "RAW_MQandDP";
pub const MAPPING_QUALITY_DEPTH_DEPRECATED: &str = "MQ_DP";
pub const RAW_RMS_MAPPING_QUALITY_DEPRECATED: &str = "RAW_MQ";
pub const TREE_SCORE: &str = "TREE_SCORE";

/// `infoFieldAnnotationKeyNamesToRemove`.
pub const INFO_KEYS_TO_REMOVE: &[&str] = &[
    "GVCFBlock",
    "HaplotypeScore",
    "InbreedingCoeff",
    "MLEAC",
    "MLEAF",
    "ExcessHet",
    "AS_InbreedingCoeff",
    "DS",
];

/// The tool's own arguments.
#[derive(Debug, Clone, PartialEq)]
pub struct Arguments {
    pub drop_low_quals: bool,
    pub rgq_threshold: f64,
    pub tree_score_threshold: f64,
    pub annotations_to_keep: Vec<String>,
    pub format_annotations_to_remove: Vec<String>,
    pub do_qual_approx: bool,
    pub allow_missing_hom_ref_data: bool,
    pub keep_all_alts: bool,
    pub keep_filters: bool,
    pub add_filters_to_format_field: bool,
}

fn user_exception(message: String) -> EngineError {
    EngineError::Runtime {
        class: "org.broadinstitute.hellbender.exceptions.UserException".to_string(),
        message,
    }
}

fn non_ref() -> Allele {
    Allele::create(b"<NON_REF>", false).expect("the symbolic <NON_REF>")
}

fn is_span_del(allele: &Allele) -> bool {
    !allele.is_reference() && allele.display_string() == htsjdk_vcf::allele::SPAN_DEL_STRING
}

/// `GATKVariantContextUtils.isConcreteAlt`.
fn is_concrete_alt(allele: &Allele) -> bool {
    !allele.is_reference() && !allele.is_symbolic() && !is_span_del(allele)
}

/// `ReblockGVCF.isHomRefBlock`: one alternate, and it is `<NON_REF>`.
pub fn is_hom_ref_block(vc: &VariantContext) -> bool {
    vc.alleles.len() == 2 && vc.alleles[1] == non_ref()
}

fn attribute<'a>(vc: &'a VariantContext, key: &str) -> Option<&'a Value> {
    vc.attributes
        .iter()
        .find(|(name, _)| name == key)
        .map(|(_, value)| value)
}

/// `getAttributeAsInt(key, default)`.
fn attribute_int(vc: &VariantContext, key: &str, default: i64) -> Result<i64, EngineError> {
    match attribute(vc, key) {
        None | Some(Value::Missing) => Ok(default),
        Some(Value::Int(value)) => Ok(*value),
        Some(other) => {
            let text = value_to_string(other);
            text.trim().parse().map_err(|_| EngineError::Runtime {
                class: "java.lang.NumberFormatException".to_string(),
                message: format!("For input string: \"{text}\""),
            })
        }
    }
}

/// `getAttributeAsDouble(key, default)`.
fn attribute_double(vc: &VariantContext, key: &str, default: f64) -> Result<f64, EngineError> {
    match attribute(vc, key) {
        None | Some(Value::Missing) => Ok(default),
        Some(Value::Int(value)) => Ok(*value as f64),
        Some(Value::Double(value)) => Ok(*value),
        Some(other) => {
            let text = value_to_string(other);
            text.trim().parse().map_err(|_| EngineError::Runtime {
                class: "java.lang.NumberFormatException".to_string(),
                message: format!("For input string: \"{text}\""),
            })
        }
    }
}

/// `MathUtils.minElementIndex`: the first minimum.
fn min_element_index(values: &[i32]) -> usize {
    let mut best = 0;
    for (index, value) in values.iter().enumerate() {
        if *value < values[best] {
            best = index;
        }
    }
    best
}

/// `GenotypesCache.get(ploidy, index).asAlleleList(alleles)`.
fn alleles_of_genotype_index(
    ploidy: usize,
    index: usize,
    alleles: &[Allele],
) -> Result<Vec<Allele>, EngineError> {
    let genotype = genotypes_in_canonical_order(ploidy, alleles.len().max(1))
        .into_iter()
        .nth(index)
        .ok_or_else(|| EngineError::Runtime {
            class: "java.lang.IndexOutOfBoundsException".to_string(),
            message: format!("genotype index {index} out of range for ploidy {ploidy}"),
        })?;
    genotype
        .into_iter()
        .map(|allele| {
            alleles
                .get(allele)
                .cloned()
                .ok_or_else(|| EngineError::Runtime {
                    class: "java.lang.IndexOutOfBoundsException".to_string(),
                    message: format!("Index {allele} out of bounds for length {}", alleles.len()),
                })
        })
        .collect()
}

/// `MathUtils.secondSmallestMinusSmallest(values, default)`.
fn second_smallest_minus_smallest(values: &[i32], default: i32) -> i32 {
    if values.len() <= 1 {
        return default;
    }
    let mut smallest = values[0];
    let mut second = i32::MAX;
    for &value in &values[1..] {
        if value < smallest {
            second = smallest;
            smallest = value;
        } else if value < second {
            second = value;
        }
    }
    second.wrapping_sub(smallest)
}

/// Whether an annotation is an `AlleleSpecificAnnotation`, and the raw value it writes for an
/// allele with none.
fn allele_specific_empty_value(entry: &Entry) -> Option<&'static str> {
    match entry.name {
        "AS_BaseQualityRankSumTest"
        | "AS_MappingQualityRankSumTest"
        | "AS_ReadPosRankSumTest"
        | "AS_ClippingRankSumTest" => Some(""),
        "AS_FisherStrand" | "AS_StrandOddsRatio" => Some("0,0"),
        "AS_RMSMappingQuality" => Some("0.00"),
        "AssemblyComplexity" | "UniqueAltReadCount" => Some("0"),
        name if name.starts_with("AS_") => Some("0"),
        _ => None,
    }
}

/// `remapList(originalList, indexesOfRelevantAlleles, offset, filler)`.
fn remap_list(original: &[String], relevant: &[usize], offset: usize, filler: &str) -> Vec<String> {
    relevant[offset..]
        .iter()
        .map(|old| {
            if *old >= original.len() + offset {
                filler.to_string()
            } else {
                original[*old - offset].clone()
            }
        })
        .collect()
}

/// `getAttributeAsString` then the bracket and whitespace stripping of
/// `getAlleleLengthListOfString[FromRawData]`.
fn allele_list(value: &Value, delimiter: char) -> Vec<String> {
    let text = match value {
        Value::List(values) => format!(
            "[{}]",
            values
                .iter()
                .map(value_to_string)
                .collect::<Vec<_>>()
                .join(", ")
        ),
        other => value_to_string(other),
    };
    let text = if text.starts_with('[') {
        text[1..text.len() - 1]
            .chars()
            .filter(|c| !c.is_whitespace())
            .collect()
    } else {
        text
    };
    text.split(delimiter).map(str::to_string).collect()
}

/// The walker: its arguments, the genotyper `--drop-low-quals` regenotypes with, the engine's
/// annotations, and the writer every record goes through.
pub struct Reblocker {
    pub arguments: Arguments,
    pub genotyping: GenotypingEngine,
    /// `annotationEngine.getInfoAnnotations()` and `getJumboInfoAnnotations()`.
    pub annotations: Vec<&'static Entry>,
    /// The FORMAT keys whose header count is A, R or G.
    pub allele_based_length_annotations: Vec<String>,
    pub writer: ReblockingWriter,
    /// `logger.warn` lines, in order.
    pub warnings: Vec<String>,
}

impl Reblocker {
    pub fn new(
        arguments: Arguments,
        genotyping: GenotypingEngine,
        annotations: &[&'static Entry],
        allele_based_length_annotations: Vec<String>,
        writer: ReblockingWriter,
    ) -> Reblocker {
        Reblocker {
            arguments,
            genotyping,
            annotations: annotations
                .iter()
                .copied()
                .filter(|entry| matches!(entry.kind, Kind::Info | Kind::JumboInfo))
                .collect(),
            allele_based_length_annotations,
            writer,
            warnings: Vec::new(),
        }
    }

    fn vcf_output_end(&self) -> Option<(String, i64)> {
        self.writer.vcf_output_end()
    }

    fn add(
        &mut self,
        vc: VariantContext,
        reference: &mut dyn ReferenceBase,
    ) -> Result<(), EngineError> {
        self.writer.add(vc, reference)
    }

    /// `apply(variant, ...)`.
    pub fn apply(
        &mut self,
        variant: VariantContext,
        reference: &mut dyn ReferenceBase,
    ) -> Result<(), EngineError> {
        if !variant.alleles.contains(&non_ref()) {
            return Err(user_exception(format!(
                "Variant Context at {}:{} does not contain a <NON-REF> allele. This tool is only \
                 intended for use with GVCFs.",
                variant.contig, variant.start
            )));
        }
        let vc = if self.arguments.format_annotations_to_remove.is_empty() {
            variant
        } else {
            self.remove_format_annotations(variant)
        };
        self.regenotype(vc, reference)
    }

    /// `removeVCFFormatAnnotations`.
    fn remove_format_annotations(&self, vc: VariantContext) -> VariantContext {
        let genotype = &vc.genotypes[0];
        if genotype.extended.is_empty() {
            return vc;
        }
        let mut rebuilt = genotype.clone();
        rebuilt
            .extended
            .retain(|(key, _)| !self.arguments.format_annotations_to_remove.contains(key));
        let mut out = vc.clone();
        out.genotypes = GenotypesContext::new(vec![rebuilt]);
        out
    }

    /// `regenotypeVC(originalVC)`.
    fn regenotype(
        &mut self,
        original: VariantContext,
        reference: &mut dyn ReferenceBase,
    ) -> Result<(), EngineError> {
        if is_hom_ref_block(&original) {
            // "if this hom ref block is entirely overlapped by previous VCF output, then drop it"
            if let Some((contig, end)) = self.vcf_output_end() {
                if contig == original.contig && original.stop <= end {
                    return Ok(());
                }
            }
            let genotype = original.genotypes[0].clone();
            if self.arguments.drop_low_quals
                && genotype
                    .gq
                    .is_none_or(|gq| f64::from(gq) < self.arguments.rgq_threshold || gq == 0)
            {
                return Ok(());
            }
            if genotype.pl.is_none() {
                if genotype.gq.is_some() {
                    self.warnings.push(format!(
                        "PL is missing for hom ref genotype at at least one position for sample \
                         {}: {}:{}.  Using GQ to determine quality.",
                        genotype.sample_name, original.contig, original.start
                    ));
                    self.add(original.clone(), reference)?;
                } else {
                    let message = format!(
                        "Homozygous reference genotypes must contain GQ or PL. Both are missing \
                         for hom ref genotype at {}:{}",
                        original.contig, original.start
                    );
                    if !self.arguments.allow_missing_hom_ref_data {
                        return Err(EngineError::Runtime {
                            class:
                                "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
                                    .to_string(),
                            message,
                        });
                    }
                    self.warnings.push(message);
                    let mut filled = genotype.clone();
                    filled.gq = Some(0);
                    filled.pl = Some(vec![0, 0, 0]);
                    let mut built = original.clone();
                    built.genotypes = GenotypesContext::new(vec![filled]);
                    self.add(built, reference)?;
                }
            }
            return self.add(original, reference);
        }

        let mut result = original.clone();
        if self.arguments.drop_low_quals
            && attribute_int(&original, "DP", 0)? > 0
            && !self.is_monomorphic_call_with_alts(&original)
        {
            let Some(regenotyped) = self.genotyping.calculate_genotypes(&original)? else {
                return Ok(());
            };
            let attributes = self.subset_annotations_if_necessary(&original, &regenotyped)?;
            result = regenotyped;
            result.attributes = attributes;
        }

        if self.should_be_reblocked(&result)? {
            if let Some((contig, end)) = self.vcf_output_end() {
                if result.contig == contig && result.stop <= end {
                    return Ok(());
                }
            }
            if let Some(block) = self.low_qual_variant_to_gq0_hom_ref(&result, reference)? {
                self.add(block, reference)?;
            }
        } else if let Some(trimmed) = self.clean_up_high_quality_variant(&result, reference)? {
            self.add(trimmed, reference)?;
        }
        Ok(())
    }

    /// `subsetAnnotationsIfNecessary`.
    fn subset_annotations_if_necessary(
        &self,
        original: &VariantContext,
        regenotyped: &VariantContext,
    ) -> Result<Vec<(String, Value)>, EngineError> {
        if regenotyped.alleles.len() == original.alleles.len() {
            return Ok(original.attributes.clone());
        }
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
        let mut attributes = Vec::new();
        self.compose_updated_annotations(&mut attributes, original, &relevant, regenotyped)?;
        Ok(attributes)
    }

    /// `isMonomorphicCallWithAlts`.
    fn is_monomorphic_call_with_alts(&self, vc: &VariantContext) -> bool {
        let genotype = &vc.genotypes[0];
        let Some(pl) = &genotype.pl else {
            return false;
        };
        let no_called_hom_ref = genotype.is_no_call() && pl.first() == Some(&0);
        (genotype.is_hom_ref() || no_called_hom_ref || min_element_index(pl) == 0)
            && vc.alleles[1..].iter().any(is_concrete_alt)
    }

    /// `shouldBeReblocked`.
    fn should_be_reblocked(&self, vc: &VariantContext) -> Result<bool, EngineError> {
        if vc.genotypes.is_empty() {
            return Err(EngineError::Runtime {
                class: "java.lang.IllegalStateException".to_string(),
                message: "Variant contexts must contain genotypes to be reblocked.".to_string(),
            });
        }
        let genotype = &vc.genotypes[0];
        let Some(pls) = &genotype.pl else {
            return Ok(true);
        };
        let best = min_element_index(pls);
        let called = alleles_of_genotype_index(genotype.ploidy(), best, &vc.alleles)?;
        let tree_score = attribute_double(vc, TREE_SCORE, 0.0)?;
        Ok(f64::from(pls[0]) < self.arguments.rgq_threshold
            || !called.iter().any(is_concrete_alt)
            || called.contains(&non_ref())
            || (genotype.pl.is_none() && genotype.gq.is_none())
            || tree_score < self.arguments.tree_score_threshold)
    }

    /// `lowQualVariantToGQ0HomRef`.
    fn low_qual_variant_to_gq0_hom_ref(
        &mut self,
        low_quality: &VariantContext,
        reference: &mut dyn ReferenceBase,
    ) -> Result<Option<VariantContext>, EngineError> {
        if self.arguments.drop_low_quals
            && (!self.is_monomorphic_call_with_alts(low_quality)
                || !low_quality.genotypes[0].is_called())
        {
            return Ok(None);
        }
        let mut attributes = Vec::new();
        let genotype = self.change_call_to_hom_ref_versus_non_ref(low_quality, &mut attributes)?;
        let mut builder = low_quality.clone();
        builder.alleles = vec![genotype.alleles[0].clone(), non_ref()];
        builder.genotypes = GenotypesContext::new(vec![genotype]);
        if let Some((contig, end)) = self.vcf_output_end() {
            if low_quality.contig == contig && low_quality.start <= end {
                let new_start = end + 1;
                if new_start > low_quality.stop {
                    return Ok(None);
                }
                move_start(&mut builder, new_start, reference);
            }
        }
        builder.filters = None;
        builder.log10_p_error = NO_LOG10_PERROR;
        builder.attributes = attributes;
        Ok(Some(builder))
    }

    /// `changeCallToHomRefVersusNonRef`.
    fn change_call_to_hom_ref_versus_non_ref(
        &self,
        low_quality: &VariantContext,
        attributes: &mut Vec<(String, Value)>,
    ) -> Result<Genotype, EngineError> {
        let genotype = &low_quality.genotypes[0];
        let ploidy = genotype.ploidy();
        let input_reference = low_quality.alleles[0].clone();
        let mut built;
        if genotype.pl.as_ref().is_none_or(|pl| pl.first() != Some(&0)) {
            built = genotype.clone();
            built.pl = Some(vec![0; ploidy + 1]);
            built.gq = Some(0);
            built.ad = None;
            built.alleles = vec![input_reference.clone(); ploidy];
            built.extended.clear();
        } else {
            // "find best ALT so we can use its likelihood for NON_REF"
            let best = most_likely_alleles_ensuring_alt(low_quality, ploidy, 1, true)?;
            let best_alt = best
                .iter()
                .find(|allele| !allele.is_reference())
                .cloned()
                .unwrap_or_else(non_ref);
            let subset = subset_alleles(
                &low_quality.genotypes,
                ploidy,
                &low_quality.alleles,
                &[input_reference.clone(), best_alt],
                SubsetMethod::BestMatchToOriginal,
            )?;
            let subset_genotype = subset[0].clone();
            built = subset_genotype.clone();
            built.extended.clear();
            if subset_genotype.gq.is_none() {
                built.gq = Some(0);
            }
            if subset_genotype.pl.is_none() {
                built.pl = Some(vec![0; ploidy + 1]);
            }
        }
        if attribute(low_quality, "DP").is_some() {
            let depth = attribute_int(low_quality, "DP", 0)? as i32;
            built.dp = Some(depth);
            put_extended(
                &mut built,
                crate::gvcf_blocks::MIN_DP_FORMAT_KEY,
                Value::Int(i64::from(depth)),
            );
        } else if let Some(ad) = &genotype.ad {
            let depth: i32 = ad.iter().sum();
            built.dp = Some(depth);
            put_extended(
                &mut built,
                crate::gvcf_blocks::MIN_DP_FORMAT_KEY,
                Value::Int(i64::from(depth)),
            );
        }
        // "If we're dropping a deletion allele, then we need to trim the reference"
        let output_reference = if input_reference.len() > 1
            || genotype.alleles.iter().any(is_span_del)
            || genotype.alleles.iter().any(Allele::is_no_call)
        {
            let base = input_reference.display_string().as_bytes()[0];
            Allele::create(&[base], true).expect("a reference base")
        } else {
            input_reference
        };
        attributes.push(("END".to_string(), Value::Int(low_quality.stop)));
        built.alleles = vec![output_reference; ploidy];
        Ok(built)
    }

    /// `getCalledGenotype`.
    fn called_genotype(&self, variant: &VariantContext) -> Result<Genotype, EngineError> {
        let original = &variant.genotypes[0];
        let Some(pls) = &original.pl else {
            return Err(EngineError::Runtime {
                class: "java.lang.IllegalStateException".to_string(),
                message: format!(
                    "Cannot verify called genotype without likelihoods or posteriors.  Error at \
                     {}:{}",
                    variant.contig, variant.start
                ),
            });
        };
        let best = min_element_index(pls);
        let called = alleles_of_genotype_index(original.ploidy(), best, &variant.alleles)?;
        let mismatch = !original
            .alleles
            .iter()
            .all(|allele| called.contains(allele));
        if original.is_no_call() || mismatch {
            let mut built = original.clone();
            let likelihoods: Vec<f64> = pls.iter().map(|pl| f64::from(*pl) / -10.0).collect();
            make_genotype_call(
                original.ploidy(),
                &mut built,
                SubsetMethod::UsePlsToAssign,
                Some(&likelihoods),
                &variant.alleles,
                original,
            )?;
            Ok(built)
        } else {
            Ok(original.clone())
        }
    }

    /// `getAllelesToDrop`.
    fn alleles_to_drop(&self, variant: &VariantContext, called: &Genotype) -> Vec<Allele> {
        let mut drop: Vec<Allele> = variant.alleles[1..]
            .iter()
            .filter(|allele| !allele.is_symbolic() && !called.alleles.contains(allele))
            .cloned()
            .collect();
        if called.alleles.iter().any(is_span_del) {
            let before_output = match self.vcf_output_end() {
                None => true,
                Some((_, end)) => end < variant.start,
            };
            if self.writer.site_overlaps_buffer(variant) || before_output {
                drop.push(
                    Allele::create(htsjdk_vcf::allele::SPAN_DEL_STRING.as_bytes(), false)
                        .expect("the spanning deletion"),
                );
            }
        }
        drop
    }

    /// `cleanUpHighQualityVariant`.
    fn clean_up_high_quality_variant(
        &mut self,
        variant: &VariantContext,
        reference: &mut dyn ReferenceBase,
    ) -> Result<Option<VariantContext>, EngineError> {
        let genotype = self.called_genotype(variant)?;
        let mut builder = variant.clone();
        builder.attributes = Vec::new();
        builder.genotypes = GenotypesContext::new(vec![genotype.clone()]);

        let to_drop = self.alleles_to_drop(variant, &genotype);
        // `new int[n]`: zeroes, read only when the alleles were subset and it was overwritten.
        let mut relevant: Vec<usize> = vec![0; variant.alleles.len()];
        let mut untrimmed = variant.alleles.clone();
        if !to_drop.is_empty() && !self.arguments.keep_all_alts {
            untrimmed.retain(|allele| !to_drop.contains(allele));
            let subset = subset_alleles_with_length_annotations(
                &variant.genotypes,
                genotype.ploidy(),
                &variant.alleles,
                &untrimmed,
                SubsetMethod::UsePlsToAssign,
                &self.allele_based_length_annotations,
            )?;
            let first = &subset[0];
            if first.is_hom_ref()
                || first.gq.is_none()
                || first.alleles.iter().any(Allele::is_no_call)
            {
                if self.arguments.drop_low_quals {
                    return Ok(None);
                }
                if let Some(block) = self.low_qual_variant_to_gq0_hom_ref(variant, reference)? {
                    self.add(block, reference)?;
                }
                return Ok(None);
            }
            builder.genotypes = GenotypesContext::new(subset);
            builder.alleles = untrimmed.clone();
            builder.stop = builder.start + builder.alleles[0].len() as i64 - 1;
            let trimmed =
                crate::variant_trim::reverse_trim_alleles(&builder).map_err(|message| {
                    EngineError::Runtime {
                        class: "java.lang.IllegalStateException".to_string(),
                        message,
                    }
                })?;
            builder = trimmed;
            relevant = untrimmed
                .iter()
                .map(|allele| {
                    variant
                        .alleles
                        .iter()
                        .position(|a| a == allele)
                        .unwrap_or(usize::MAX)
                })
                .collect();
            let depth = if attribute(variant, "DP").is_some() {
                attribute_int(variant, "DP", 0)? as i32
            } else {
                genotype.dp.unwrap_or(0)
            };
            self.add_ref_block_if_necessary(variant, &to_drop, &builder, depth, reference)?;
        }

        let mut updated = builder.clone();
        let non_ref_index = updated.alleles.iter().position(|a| *a == non_ref());
        let fixed = remove_non_ref_ads(&mut updated, non_ref_index)?;
        let final_genotype = if self.arguments.add_filters_to_format_field {
            let mut with_filters = fixed;
            let filters: Vec<String> = updated.filters.clone().unwrap_or_default();
            with_filters.filters = match filters.len() {
                0 => None,
                1 => Some(filters[0].clone()),
                _ => {
                    let mut sorted = filters.clone();
                    sorted.sort();
                    Some(sorted.join(";"))
                }
            };
            with_filters
        } else {
            fixed
        };
        builder.genotypes = GenotypesContext::new(vec![final_genotype]);

        let mut attributes = Vec::new();
        self.compose_updated_annotations(&mut attributes, variant, &relevant, &updated)?;
        builder.attributes = attributes;
        if !self.arguments.keep_filters {
            builder.filters = None;
        }
        Ok(Some(builder))
    }

    /// `addRefBlockIfNecessary`.
    fn add_ref_block_if_necessary(
        &mut self,
        original: &VariantContext,
        to_drop: &[Allele],
        trimmed: &VariantContext,
        depth: i32,
        reference: &mut dyn ReferenceBase,
    ) -> Result<(), EngineError> {
        let old_length = original.alleles[0].len() as i64;
        let new_length = trimmed.alleles[0].len() as i64;
        let genotype = &original.genotypes[0];
        let vcf_output_end = self.vcf_output_end().map_or(-1, |(_, end)| end);
        if new_length >= old_length {
            return Ok(());
        }
        let Some(likelihoods) = &genotype.pl else {
            return Ok(());
        };
        // `Stream.min`, which keeps the first of equals.
        let mut shortest: Option<&Allele> = None;
        for allele in to_drop.iter().filter(|allele| !is_span_del(allele)) {
            if shortest
                .is_none_or(|best| allele.display_string().len() < best.display_string().len())
            {
                shortest = Some(allele);
            }
        }
        let Some(shortest) = shortest else {
            return Err(EngineError::Runtime {
                class: "org.broadinstitute.hellbender.exceptions.GATKException".to_string(),
                message: format!(
                    "No shortest ALT at {} across alleles: [{}]",
                    original.start,
                    to_drop
                        .iter()
                        .map(|allele| format!(
                            "{}{}",
                            allele.display_string(),
                            if allele.is_reference() { "*" } else { "" }
                        ))
                        .collect::<Vec<_>>()
                        .join(", ")
                ),
            });
        };
        let shortest_index = original
            .alleles
            .iter()
            .position(|a| a == shortest)
            .expect("a dropped allele of the record");
        let indices =
            subsetted_pl_indices(genotype.ploidy(), &[0, shortest_index]).map_err(|error| {
                EngineError::Runtime {
                    class: error.java_class().to_string(),
                    message: error.message(),
                }
            })?;
        let picked: Vec<i32> = indices.iter().map(|index| likelihoods[*index]).collect();
        let smallest = picked.iter().copied().min().unwrap_or(0);
        let mut normalized: Vec<i32> = picked.iter().map(|pl| pl - smallest).collect();
        if normalized[0] != 0 {
            // The first iteration zeroes the entry every later one subtracts.
            for index in 0..normalized.len() {
                normalized[index] = (normalized[index] - normalized[0]).max(0);
            }
        }
        let ref_start = (original.stop - (old_length - new_length)).max(vcf_output_end) + 1;
        let base = reference.base(&original.contig, ref_start);
        let new_reference = Allele::create(&[base], true).expect("a reference base");
        let mut block_genotype =
            Genotype::new("", vec![new_reference.clone(), new_reference.clone()]);
        block_genotype.gq = Some(second_smallest_minus_smallest(&normalized, 0));
        block_genotype.pl = Some(normalized);
        block_genotype.dp = Some(depth);
        if ref_start > vcf_output_end && original.stop > vcf_output_end {
            let start = ref_start.max(vcf_output_end + 1);
            let mut block =
                VariantContext::new(&original.contig, start, vec![new_reference, non_ref()]);
            block.stop = original.stop;
            block.attributes = vec![("END".to_string(), Value::Int(original.stop))];
            block.genotypes = GenotypesContext::new(vec![block_genotype]);
            self.add(block, reference)?;
        }
        Ok(())
    }

    /// `composeUpdatedAnnotations`.
    fn compose_updated_annotations(
        &self,
        destination: &mut Vec<(String, Value)>,
        variant: &VariantContext,
        relevant: &[usize],
        updated: &VariantContext,
    ) -> Result<(), EngineError> {
        update_mq_annotations(destination, variant)?;
        let needs_subsetting = relevant.len() < variant.alleles.len();
        self.copy_info_annotations(destination, variant, needs_subsetting, relevant);
        let updated_genotype = &updated.genotypes[0];
        if self.arguments.do_qual_approx {
            if updated_genotype.pl.is_some() {
                self.add_qual_annotations(destination, updated)?;
            }
        } else {
            for key in [AS_VARIANT_DEPTH_KEY, RAW_QUAL_APPROX_KEY] {
                if let Some(value) = attribute(variant, key) {
                    put(destination, key, value.clone());
                }
            }
        }
        let counts = if updated_genotype.alleles.iter().any(Allele::is_reference) {
            [0, 1, 0]
        } else {
            [0, 0, 1]
        };
        put(
            destination,
            RAW_GENOTYPE_COUNT_KEY,
            Value::List(counts.iter().map(|count| Value::Int(*count)).collect()),
        );
        for key in &self.arguments.annotations_to_keep {
            if let Some(value) = attribute(variant, key) {
                put(destination, key, value.clone());
            }
        }
        Ok(())
    }

    /// `copyInfoAnnotations`.
    fn copy_info_annotations(
        &self,
        destination: &mut Vec<(String, Value)>,
        source: &VariantContext,
        needs_subsetting: bool,
        relevant: &[usize],
    ) {
        for entry in &self.annotations {
            for key in entry.key_names {
                if INFO_KEYS_TO_REMOVE.contains(key) {
                    continue;
                }
                if let Some(value) = attribute(source, key) {
                    put(destination, key, value.clone());
                }
            }
            let Some(empty) = allele_specific_empty_value(entry) else {
                continue;
            };
            let reducible = entry.is_reducible();
            let keys: &[&str] = if reducible {
                entry.raw_keys.unwrap_or(&[])
            } else {
                entry.key_names
            };
            for raw_key in keys {
                if INFO_KEYS_TO_REMOVE.contains(raw_key) {
                    continue;
                }
                let Some(value) = attribute(source, raw_key) else {
                    continue;
                };
                if needs_subsetting && raw_key.starts_with("AS_") {
                    let values = allele_list(value, if reducible { '|' } else { ',' });
                    let mut subset = if reducible {
                        remap_list(&values, relevant, 0, "")
                    } else {
                        remap_list(&values, relevant, 1, "")
                    };
                    let last = relevant[relevant.len() - 1];
                    if source.alleles.get(last) == Some(&non_ref()) {
                        // "zero out non-ref value, just in case"
                        let end = subset.len() - 1;
                        subset[end] = empty.to_string();
                    }
                    let encoded = if reducible {
                        subset.join("|").replace(['[', ']'], "")
                    } else {
                        subset.join(",")
                    };
                    put(destination, raw_key, Value::Str(encoded));
                } else {
                    put(destination, raw_key, value.clone());
                }
            }
        }
    }

    /// `addQualAnnotations`.
    fn add_qual_annotations(
        &self,
        destination: &mut Vec<(String, Value)>,
        updated: &VariantContext,
    ) -> Result<(), EngineError> {
        let genotype = &updated.genotypes[0];
        let pl0 = genotype
            .pl
            .as_ref()
            .and_then(|pl| pl.first().copied())
            .unwrap_or(0);
        put(destination, RAW_QUAL_APPROX_KEY, Value::Int(i64::from(pl0)));
        let mut depth = gatk_annotation::site_statistics::qual_by_depth_depth(updated, None);
        if depth == 0 {
            // "prevent QD=Infinity": the DP the copy put, or 1.
            depth = match destination.iter().find(|(key, _)| key == "DP") {
                Some((_, value)) => {
                    let text = value_to_string(value);
                    text.trim().parse().map_err(|_| EngineError::Runtime {
                        class: "java.lang.NumberFormatException".to_string(),
                        message: format!("For input string: \"{text}\""),
                    })?
                }
                None => 1,
            };
        }
        put(destination, VARIANT_DEPTH_KEY, Value::Int(i64::from(depth)));
        if self
            .annotations
            .iter()
            .any(|entry| entry.name == "AS_QualByDepth")
        {
            let mut quals = Vec::new();
            for alternate in &updated.alleles[1..] {
                if *alternate == non_ref() || is_span_del(alternate) {
                    quals.push("0".to_string());
                    continue;
                }
                let subset = subset_alleles(
                    &updated.genotypes,
                    genotype.ploidy(),
                    &updated.alleles,
                    &[updated.alleles[0].clone(), alternate.clone()],
                    SubsetMethod::BestMatchToOriginal,
                )?;
                match subset[0].pl.as_ref().and_then(|pl| pl.first()) {
                    Some(pl0) => quals.push(pl0.to_string()),
                    None => quals.push("0".to_string()),
                }
            }
            put(
                destination,
                AS_RAW_QUAL_APPROX_KEY,
                Value::Str(format!("|{}", quals.join("|"))),
            );
            if let Some(depths) =
                gatk_annotation::allele_specific_site_statistics::allele_depths(updated)
            {
                put(
                    destination,
                    AS_VARIANT_DEPTH_KEY,
                    Value::Str(
                        depths
                            .iter()
                            .map(i32::to_string)
                            .collect::<Vec<_>>()
                            .join("|"),
                    ),
                );
            }
        }
        Ok(())
    }
}

fn put(destination: &mut Vec<(String, Value)>, key: &str, value: Value) {
    match destination.iter_mut().find(|(name, _)| name == key) {
        Some((_, slot)) => *slot = value,
        None => destination.push((key.to_string(), value)),
    }
}

fn put_extended(genotype: &mut Genotype, key: &str, value: Value) {
    match genotype.extended.iter_mut().find(|(name, _)| name == key) {
        Some((_, slot)) => *slot = value,
        None => genotype.extended.push((key.to_string(), value)),
    }
}

/// `removeNonRefADs`, which zeroes the `<NON_REF>` depth in the record's own genotype as well as
/// in the one it returns.
fn remove_non_ref_ads(
    updated: &mut VariantContext,
    non_ref_index: Option<usize>,
) -> Result<Genotype, EngineError> {
    let genotype = updated.genotypes[0].clone();
    let (Some(ad), Some(index)) = (&genotype.ad, non_ref_index) else {
        return Ok(genotype);
    };
    if ad.len() < index {
        return Ok(genotype);
    }
    let Some(&count) = ad.get(index) else {
        return Err(EngineError::Runtime {
            class: "java.lang.ArrayIndexOutOfBoundsException".to_string(),
            message: format!("Index {index} out of bounds for length {}", ad.len()),
        });
    };
    if count <= 0 {
        return Ok(genotype);
    }
    let mut new_ad = ad.clone();
    new_ad[index] = 0;
    let mut shared = genotype.clone();
    shared.ad = Some(new_ad.clone());
    updated.genotypes = GenotypesContext::new(vec![shared]);
    let mut fixed = genotype.clone();
    fixed.ad = Some(new_ad.clone());
    fixed.dp = Some(match genotype.dp {
        Some(dp) => dp - count,
        None => new_ad.iter().sum(),
    });
    Ok(fixed)
}

/// `updateMQAnnotations`.
fn update_mq_annotations(
    destination: &mut Vec<(String, Value)>,
    source: &VariantContext,
) -> Result<(), EngineError> {
    if attribute(source, RAW_MAPPING_QUALITY_WITH_DEPTH_KEY).is_some() {
        return Ok(());
    }
    let depth = attribute_int(source, "DP", 0)?;
    let raw = if attribute(source, RAW_RMS_MAPPING_QUALITY_DEPRECATED).is_some() {
        java_round(attribute_double(
            source,
            RAW_RMS_MAPPING_QUALITY_DEPRECATED,
            0.0,
        )?)
    } else {
        let mq = attribute_double(source, "MQ", 60.0)?;
        java_round(mq * mq * depth as f64)
    } as i32;
    put(
        destination,
        RAW_MAPPING_QUALITY_WITH_DEPTH_KEY,
        Value::Str(format!("{raw},{depth}")),
    );
    if attribute(source, RAW_RMS_MAPPING_QUALITY_DEPRECATED).is_some() {
        put(
            destination,
            RAW_RMS_MAPPING_QUALITY_DEPRECATED,
            Value::Double(attribute_double(
                source,
                RAW_RMS_MAPPING_QUALITY_DEPRECATED,
                0.0,
            )?),
        );
        put(
            destination,
            MAPPING_QUALITY_DEPTH_DEPRECATED,
            Value::Int(depth),
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_reducible_list_is_remapped_with_the_reference_and_a_reduced_one_without() {
        let values: Vec<String> = ["", "a", "b", "c"].iter().map(|s| s.to_string()).collect();
        assert_eq!(remap_list(&values, &[0, 2, 3], 0, ""), vec!["", "b", "c"]);
        let alts: Vec<String> = ["a", "b", "c"].iter().map(|s| s.to_string()).collect();
        assert_eq!(remap_list(&alts, &[0, 2, 3], 1, ""), vec!["b", "c"]);
    }

    #[test]
    fn a_bracketed_list_loses_its_brackets_and_spaces() {
        let value = Value::List(vec![Value::Str("1".into()), Value::Str("2".into())]);
        assert_eq!(allele_list(&value, ','), vec!["1", "2"]);
        assert_eq!(
            allele_list(&Value::Str("|1|2".into()), '|'),
            vec!["", "1", "2"]
        );
    }

    #[test]
    fn second_smallest_minus_smallest_is_the_gap_or_the_default() {
        assert_eq!(second_smallest_minus_smallest(&[0, 25, 250], 0), 25);
        assert_eq!(second_smallest_minus_smallest(&[7], 3), 3);
    }
}

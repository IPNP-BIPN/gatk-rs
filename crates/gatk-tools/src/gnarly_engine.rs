//! `GnarlyGenotyperEngine.finalizeGenotype` and the walker's `apply` around it.
//!
//! Ported from `org.broadinstitute.hellbender.tools.walkers.gnarlyGenotyper` (GATK 4.6.2.0). The
//! site's QUAL is not computed but read: `QUALapprox`, the sum of PL[0] over the samples that
//! `ReblockGVCF` wrote, or the sum of `AS_QUALapprox`. It is held against a floor, the confidence
//! threshold less ten times the log of the site prior, 60 for a SNP and about 69.03 for an indel.
//! A site over it gets its MQ, QD and allele-specific annotations finalized, its genotypes called
//! from their own PLs with `<NON_REF>` cut from every array, and AC, AF, AN, ExcessHet, FS and
//! SOR tallied from the calls.
//!
//! # The attributes a genotype keeps
//!
//! `iterateOnGenotypes` removes `MIN_DP` from a COPY of the genotype's attributes and adds the copy
//! back onto a builder that still holds the original, so `MIN_DP` survives. Only a genotype
//! `makeGenotypeCall` stripped (a GQ-0, DP-0 reference call) loses it.

use gatk_annotation::catalogue;
use gatk_engine::genotype_index::genotype_count;
use gatk_engine::java_random::JavaRandom;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::genotypes_context::GenotypesContext;
use htsjdk_vcf::variant::{Genotype, Value, VariantContext};

use crate::genotyping_engine::{make_genotype_call, EngineError, SubsetMethod};
use crate::reference_confidence_merger::value_to_string;

/// `HomoSapiensConstants.SNP_HETEROZYGOSITY` and `INDEL_HETEROZYGOSITY`.
const SNP_HETEROZYGOSITY: f64 = 1e-3;
const INDEL_HETEROZYGOSITY: f64 = 1.25e-4;
/// `GenotypeCalculationArgumentCollection.DEFAULT_STANDARD_CONFIDENCE_FOR_CALLING`.
const DEFAULT_STANDARD_CONFIDENCE_FOR_CALLING: f64 = 30.0;

/// The allele-specific reducible annotations the engine finds by reflection, with the one key
/// each finalizes and trims.
const AS_REDUCIBLE: &[(&str, &str)] = &[
    ("AS_QualByDepth", "AS_QD"),
    ("AS_RMSMappingQuality", "AS_MQ"),
    ("AS_FisherStrand", "AS_FS"),
    ("AS_StrandOddsRatio", "AS_SOR"),
    ("AS_BaseQualityRankSumTest", "AS_BaseQRankSum"),
    ("AS_MappingQualityRankSumTest", "AS_MQRankSum"),
    ("AS_ReadPosRankSumTest", "AS_ReadPosRankSum"),
];

fn non_ref() -> Allele {
    Allele::create(b"<NON_REF>", false).expect("the symbolic <NON_REF>")
}

fn is_span_del(allele: &Allele) -> bool {
    !allele.is_reference() && allele.display_string() == htsjdk_vcf::allele::SPAN_DEL_STRING
}

fn attribute<'a>(vc: &'a VariantContext, key: &str) -> Option<&'a Value> {
    vc.attributes
        .iter()
        .find(|(name, _)| name == key)
        .map(|(_, value)| value)
}

fn put(vc: &mut VariantContext, key: &str, value: Value) {
    match vc.attributes.iter_mut().find(|(name, _)| name == key) {
        Some((_, slot)) => *slot = value,
        None => vc.attributes.push((key.to_string(), value)),
    }
}

fn remove(vc: &mut VariantContext, key: &str) {
    vc.attributes.retain(|(name, _)| name != key);
}

fn number_format(text: &str) -> EngineError {
    EngineError::Runtime {
        class: "java.lang.NumberFormatException".to_string(),
        message: format!("For input string: \"{text}\""),
    }
}

/// `getAttributeAsInt(key, default)`.
fn attribute_int(vc: &VariantContext, key: &str, default: i64) -> Result<i64, EngineError> {
    match attribute(vc, key) {
        None | Some(Value::Missing) => Ok(default),
        Some(Value::Int(value)) => Ok(*value),
        Some(other) => {
            let text = value_to_string(other);
            text.trim().parse().map_err(|_| number_format(&text))
        }
    }
}

/// `getAttributeAsIntList(key, 0)`.
fn attribute_int_list(vc: &VariantContext, key: &str) -> Result<Vec<i64>, EngineError> {
    let values: Vec<String> = match attribute(vc, key) {
        None => Vec::new(),
        Some(Value::List(values)) => values.iter().map(value_to_string).collect(),
        Some(other) => value_to_string(other)
            .split(',')
            .map(str::to_string)
            .collect(),
    };
    values
        .iter()
        .map(|text| {
            if text.trim() == "." {
                Ok(0)
            } else {
                text.trim().parse().map_err(|_| number_format(text))
            }
        })
        .collect()
}

/// `RMSMappingQuality.finalizeRawMQ(vc)`.
pub fn finalize_raw_mq(vc: &VariantContext) -> Result<VariantContext, EngineError> {
    let Some(raw) = attribute(vc, "RAW_MQandDP").map(value_to_string) else {
        if attribute(vc, "MQ_DP").is_some() {
            return Err(EngineError::Limitation(
                "the deprecated MQ_DP and RAW_MQ keys are not ported for GnarlyGenotyper."
                    .to_string(),
            ));
        }
        return Ok(vc.clone());
    };
    let (square_sum, depth) = gatk_annotation::mapping_quality::parse_raw_data_string(&raw)
        .map_err(|error| EngineError::Runtime {
            class: "org.broadinstitute.hellbender.exceptions.UserException$BadInput".to_string(),
            message: format!("{error:?}"),
        })?;
    let mut out = vc.clone();
    remove(&mut out, "RAW_MQ");
    remove(&mut out, "RAW_MQandDP");
    if depth > 0 {
        put(
            &mut out,
            "MQ",
            Value::Str(gatk_annotation::mapping_quality::rms_from_tuple(
                square_sum, depth,
            )),
        );
    }
    Ok(out)
}

/// `AS_QualByDepth.parseQualList(vc)` summed, for a site with `AS_QUALapprox` alone.
fn allele_specific_qual_sum(vc: &VariantContext) -> Result<f64, EngineError> {
    let text = attribute(vc, "AS_QUALapprox")
        .map(value_to_string)
        .unwrap_or_default();
    let values: Vec<&str> = text.split('|').collect();
    if values.len() != vc.alleles.len() {
        return Err(EngineError::Runtime {
            class: "java.lang.IllegalStateException".to_string(),
            message: "Number of AS_QUALapprox values doesn't match the number of alleles in the \
                      variant context."
                .to_string(),
        });
    }
    let mut sum = 0i64;
    for value in &values[1..] {
        if !value.is_empty() {
            sum += value.parse::<i64>().map_err(|_| number_format(value))?;
        }
    }
    Ok(sum as f64)
}

/// `trimASAnnotation`: the entries of `key` at the target alternates' positions among `vc`'s.
fn trim_as_annotation(vc: &VariantContext, targets: &[Allele], key: &str) -> Value {
    let Some(value) = attribute(vc, key) else {
        return Value::Missing;
    };
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
    let entries: Vec<String> = text
        .replace(['[', ']'], "")
        .split(',')
        .map(str::to_string)
        .collect();
    let alternates = &vc.alleles[1..];
    let trimmed: Vec<String> = targets
        .iter()
        .filter(|allele| !allele.is_reference())
        .map(|allele| {
            let index = alternates.iter().position(|a| a == allele);
            match index {
                Some(index) if index < entries.len() => entries[index].clone(),
                _ => ".".to_string(),
            }
        })
        .collect();
    Value::Str(trimmed.join(","))
}

/// `getSBFieldAsIntArray`.
fn sb_counts(genotype: &Genotype) -> Result<Vec<i32>, EngineError> {
    let Some((_, value)) = genotype.extended.iter().find(|(key, _)| key == "SB") else {
        return Ok(Vec::new());
    };
    let text = match value {
        Value::List(values) => values
            .iter()
            .map(value_to_string)
            .collect::<Vec<_>>()
            .join(","),
        other => value_to_string(other),
    };
    text.split(',')
        .map(|count| count.trim().parse::<i32>())
        .collect::<Result<_, _>>()
        .map_err(|_| EngineError::Runtime {
            class: "java.lang.IllegalStateException".to_string(),
            message: "The GnarlyGenotyper tool assumes that input variants have SB FORMAT  \
                      fields as a list of integers separated by commas."
                .to_string(),
        })
}

/// `GnarlyGenotyperEngine`.
pub struct GnarlyEngine {
    pub keep_all_sites: bool,
    pub max_alt_alleles_to_output: usize,
    pub strip_as_annotations: bool,
}

impl GnarlyEngine {
    /// `finalizeGenotype(variant, annotationDBBuilder)`: the site, or `None` under the floor, and
    /// what the annotation database gets when one is written.
    pub fn finalize_genotype(
        &self,
        variant: &VariantContext,
        with_database: bool,
        random: &mut JavaRandom,
    ) -> Result<(Option<VariantContext>, Option<VariantContext>), EngineError> {
        let mut database = with_database.then(|| variant.clone());
        let qual_approx = if attribute(variant, "QUALapprox").is_some() {
            attribute_int(variant, "QUALapprox", 0)? as f64
        } else if attribute(variant, "AS_QUALapprox").is_some() {
            allele_specific_qual_sum(variant)?
        } else {
            0.0
        };
        // "don't count a '*' as a SNP"
        let reference_length = variant.alleles[0].len();
        let has_snp = variant.alleles[1..]
            .iter()
            .any(|allele| !is_span_del(allele) && allele.len() == reference_length);
        let is_indel = !has_snp;
        let prior = if is_indel {
            INDEL_HETEROZYGOSITY
        } else {
            SNP_HETEROZYGOSITY
        };
        let floor = DEFAULT_STANDARD_CONFIDENCE_FOR_CALLING - 10.0 * jmath::math::log10(prior);
        if qual_approx < floor {
            if self.keep_all_sites {
                let mut kept = finalize_raw_mq(variant)?;
                kept.filters = Some(vec!["LowQual".to_string()]);
                put(&mut kept, "AC_adj", Value::Int(0));
                return Ok((Some(kept), database));
            }
            return Ok((None, database));
        }

        let with_mq = finalize_raw_mq(variant)?;
        let mut builder = with_mq.clone();
        let mut modified = with_mq.attributes.clone();
        for (name, _) in AS_REDUCIBLE {
            let entry = catalogue::entry(name).expect("a catalogued annotation");
            let primary = entry.raw_keys.and_then(|keys| keys.first()).copied();
            let Some(primary) = primary else { continue };
            if attribute(variant, primary).is_none() || self.strip_as_annotations {
                continue;
            }
            let finalized =
                crate::genotype_gvcfs::finalize_raw_data(entry, &with_mq, variant, random)?;
            for (key, value) in finalized {
                match modified.iter_mut().find(|(name, _)| *name == key) {
                    Some((_, slot)) => *slot = value,
                    None => modified.push((key, value)),
                }
            }
            if let Some(database) = database.as_mut() {
                let raw = attribute(variant, primary).cloned().expect("the raw key");
                put(database, primary, raw);
            }
        }
        builder.attributes = modified;

        if attribute(variant, "VarDP").is_some() {
            let depth = attribute_int(variant, "VarDP", 0)?;
            put(
                &mut builder,
                "QD",
                Value::Double(qual_approx / depth as f64),
            );
            builder.log10_p_error = qual_approx / -10.0 - jmath::math::log10(prior);
        }
        if !self.keep_all_sites {
            remove(&mut builder, "QUALapprox");
        }

        let mut sb_sum = [0i32; 4];
        let remove_non_ref = variant.alleles.contains(&non_ref());
        let targets: Vec<Allele> = if remove_non_ref {
            variant.alleles[..variant.alleles.len() - 1].to_vec()
        } else {
            variant.alleles.clone()
        };
        let mut counts: Vec<(Allele, i64)> =
            targets.iter().map(|allele| (allele.clone(), 0)).collect();
        let tally_genotypes = attribute(variant, "RAW_GT_COUNT").is_none();
        let mut raw_genotype_counts = [0i64; 3];
        let called = self.iterate_on_genotypes(
            variant,
            &targets,
            &mut counts,
            &mut sb_sum,
            remove_non_ref,
            tally_genotypes.then_some(&mut raw_genotype_counts),
        )?;
        let mut called_alleles = 0i64;
        if !variant.genotypes.is_empty() {
            called_alleles = counts.iter().map(|(_, count)| count).sum();
            let alternates: Vec<i64> = counts
                .iter()
                .filter(|(allele, _)| !allele.is_reference())
                .map(|(_, count)| *count)
                .collect();
            let ac = if alternates.len() == 1 {
                Value::Int(alternates[0])
            } else {
                Value::List(alternates.iter().map(|count| Value::Int(*count)).collect())
            };
            let frequencies: Vec<f64> = alternates
                .iter()
                .map(|count| *count as f64 / called_alleles as f64)
                .collect();
            let af = if frequencies.len() == 1 {
                Value::Double(frequencies[0])
            } else {
                Value::List(frequencies.iter().map(|f| Value::Double(*f)).collect())
            };
            put(&mut builder, "AC", ac.clone());
            put(&mut builder, "AF", af.clone());
            put(&mut builder, "AN", Value::Int(called_alleles));
            if let Some(database) = database.as_mut() {
                put(database, "AC", ac);
                put(database, "AF", af);
                put(database, "AN", Value::Int(called_alleles));
            }
        } else {
            if attribute(variant, "SB_TABLE").is_some() {
                let table = attribute_int_list(variant, "SB_TABLE")?;
                for (slot, value) in sb_sum.iter_mut().zip(table) {
                    *slot = value as i32;
                }
            }
            if let Some(database) = database.as_mut() {
                for key in ["AC", "AF", "AN"] {
                    match attribute(variant, key).cloned() {
                        Some(value) => put(database, key, value),
                        None => put(database, key, Value::Missing),
                    }
                }
            }
        }
        if attribute(variant, "RAW_GT_COUNT").is_some() || !variant.genotypes.is_empty() {
            let mut genotype_counts: Vec<i64> = if attribute(variant, "RAW_GT_COUNT").is_some() {
                attribute_int_list(variant, "RAW_GT_COUNT")?
            } else {
                raw_genotype_counts.to_vec()
            };
            // `gtCounts.get(1)`, then `get(2)`: the first index the list does not reach.
            if genotype_counts.len() < 3 {
                return Err(EngineError::Runtime {
                    class: "java.lang.IndexOutOfBoundsException".to_string(),
                    message: format!(
                        "Index {} out of bounds for length {}",
                        genotype_counts.len().max(1),
                        genotype_counts.len()
                    ),
                });
            }
            let reference_count =
                (called_alleles / 2 - genotype_counts[1] - genotype_counts[2]).max(0);
            genotype_counts[0] = reference_count;
            let excess = gatk_annotation::heterozygosity::excess_het_value(
                gatk_annotation::heterozygosity::GenotypeCounts {
                    refs: genotype_counts[0] as f64,
                    hets: genotype_counts[1] as f64,
                    homs: genotype_counts[2] as f64,
                },
                (called_alleles / 2) as usize,
            )
            .ok_or_else(|| EngineError::Runtime {
                class: "java.lang.IllegalArgumentException".to_string(),
                message: "genotype counts must be non-negative".to_string(),
            })?;
            put(&mut builder, "ExcessHet", Value::Str(excess));
            remove(&mut builder, "RAW_GT_COUNT");
            if let Some(database) = database.as_mut() {
                put(
                    database,
                    "RAW_GT_COUNT",
                    Value::Str(
                        genotype_counts
                            .iter()
                            .map(i64::to_string)
                            .collect::<Vec<_>>()
                            .join(","),
                    ),
                );
            }
        }

        let table = gatk_annotation::strand_bias::decode_sbbs(&sb_sum);
        put(
            &mut builder,
            "FS",
            Value::Str(gatk_annotation::strand_bias::fisher_strand_value(table)),
        );
        put(
            &mut builder,
            "SOR",
            Value::Str(gatk_annotation::strand_bias::strand_odds_ratio_value(table)),
        );
        builder.genotypes = GenotypesContext::new(called);
        if let Some(database) = database.as_mut() {
            put(
                database,
                "SB_TABLE",
                Value::List(sb_sum.iter().map(|v| Value::Int(i64::from(*v))).collect()),
            );
            database.genotypes = GenotypesContext::new(Vec::new());
        }

        for (name, key) in AS_REDUCIBLE {
            let entry = catalogue::entry(name).expect("a catalogued annotation");
            let raw_keys = entry.raw_keys.unwrap_or(&[]);
            let Some(primary) = raw_keys.first() else {
                continue;
            };
            if attribute(variant, primary).is_some() {
                let trimmed = trim_as_annotation(&builder, &targets, key);
                put(&mut builder, key, trimmed);
                if !self.keep_all_sites {
                    remove(&mut builder, primary);
                }
            }
        }
        if !self.keep_all_sites {
            for (name, _) in AS_REDUCIBLE {
                let entry = catalogue::entry(name).expect("a catalogued annotation");
                for raw_key in entry.raw_keys.unwrap_or(&[]) {
                    if attribute(variant, raw_key).is_some() {
                        remove(&mut builder, raw_key);
                    }
                }
            }
        }
        if attribute(variant, "AS_VarDP").is_some() {
            let raw = attribute(variant, "AS_VarDP")
                .map(value_to_string)
                .unwrap_or_default();
            let values: Vec<&str> = raw.split('|').collect();
            let finalized = if values.len() != targets.len() + 1 {
                Value::Missing
            } else {
                Value::Str(values[1..values.len() - 1].join(","))
            };
            put(&mut builder, "AS_AltDP", finalized);
            remove(&mut builder, "AS_VarDP");
        }
        if let Some(database) = database.as_mut() {
            database.alleles = targets.clone();
        }
        builder.alleles = targets;
        Ok((Some(builder), database))
    }

    /// `iterateOnGenotypes`.
    fn iterate_on_genotypes(
        &self,
        vc: &VariantContext,
        targets: &[Allele],
        counts: &mut [(Allele, i64)],
        sb_sum: &mut [i32; 4],
        non_ref_returned: bool,
        mut raw_genotype_counts: Option<&mut [i64; 3]>,
    ) -> Result<Vec<Genotype>, EngineError> {
        let max_alleles_to_output = self.max_alt_alleles_to_output + 1;
        if non_ref_returned && vc.alleles.last() != Some(&non_ref()) {
            return Err(EngineError::Runtime {
                class: "java.lang.IllegalStateException".to_string(),
                message: format!(
                    "This tool assumes that the NON_REF allele is listed last, as in \
                     HaplotypeCaller GVCF output, but that was not the case at position {}:{}.",
                    vc.contig, vc.start
                ),
            });
        }
        let maximum_allele_count = vc.alleles.len();
        let concrete_alternates = maximum_allele_count.saturating_sub(2);
        let mut out = Vec::with_capacity(vc.genotypes.len());
        for g in vc.genotypes.iter() {
            let mut built = g.clone();
            if g.alleles.contains(&non_ref()) {
                built.alleles = vec![Allele::no_call(); g.ploidy()];
                built.gq = None;
            } else if g.pl.is_none() && g.ad.is_none() && g.is_no_call() {
                built.alleles = vec![Allele::no_call(); 2];
                built.gq = None;
            }
            if non_ref_returned {
                if let Some(ad) = &g.ad {
                    built.ad = Some(ad.iter().take(targets.len()).copied().collect());
                }
            }
            if let Some(pl) = &g.pl {
                let size = if maximum_allele_count <= max_alleles_to_output && g.ploidy() == 2 {
                    genotype_count(2, concrete_alternates + 1)
                } else {
                    genotype_count(g.ploidy(), concrete_alternates + 1)
                }
                .map_err(|error| EngineError::Runtime {
                    class: error.java_class().to_string(),
                    message: error.message(),
                })?;
                if pl.len() < size {
                    return Err(EngineError::Runtime {
                        class: "java.lang.ArrayIndexOutOfBoundsException".to_string(),
                        message: format!(
                            "arraycopy: last source index {size} out of bounds for int[{}]",
                            pl.len()
                        ),
                    });
                }
                let pls: Vec<i32> = pl[..size].to_vec();
                built.gq = Some(second_smallest_minus_smallest(&pls, 0));
                built.pl = Some(pls.clone());
                let likelihoods: Vec<f64> = pls.iter().map(|p| f64::from(*p) / -10.0).collect();
                make_genotype_call(
                    g.ploidy(),
                    &mut built,
                    SubsetMethod::UsePlsToAssign,
                    Some(&likelihoods),
                    targets,
                    g,
                )?;
            } else if g.gq == Some(0) {
                make_genotype_call(
                    g.ploidy(),
                    &mut built,
                    SubsetMethod::UsePlsToAssign,
                    None,
                    targets,
                    g,
                )?;
            }
            // `attributes(attrs)` adds the copy without MIN_DP onto what the builder kept.
            for (key, value) in &g.extended {
                if key == crate::gvcf_blocks::MIN_DP_FORMAT_KEY {
                    continue;
                }
                match built.extended.iter_mut().find(|(name, _)| name == key) {
                    Some((_, slot)) => *slot = value.clone(),
                    None => built.extended.push((key.clone(), value.clone())),
                }
            }
            if g.extended.iter().any(|(key, _)| key == "SB") {
                let sb = sb_counts(g)?;
                for (slot, value) in sb_sum.iter_mut().zip(sb) {
                    *slot += value;
                }
            }
            for allele in &built.alleles {
                if allele.is_no_call() {
                    continue;
                }
                if let Some((_, count)) = counts.iter_mut().find(|(a, _)| a == allele) {
                    *count += 1;
                }
            }
            if let Some(raw) = raw_genotype_counts.as_deref_mut() {
                let alternates = g.alleles.iter().filter(|a| !a.is_reference()).count();
                if alternates >= raw.len() {
                    return Err(EngineError::Runtime {
                        class: "java.lang.ArrayIndexOutOfBoundsException".to_string(),
                        message: format!("Index {alternates} out of bounds for length 3"),
                    });
                }
                raw[alternates] += 1;
            }
            out.push(built);
        }
        Ok(out)
    }
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

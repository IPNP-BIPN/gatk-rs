//! `GnarlyGenotyperEngine.finalizeGenotype` and the walker's `apply` around it.
//!
//! Ported from `org.broadinstitute.hellbender.tools.walkers.gnarlyGenotyper` (GATK 4.7.0.0;
//! `GnarlyGenotyperEngine.java` and `GnarlyGenotyper.java` are unchanged since 4.6.2.0). The
//! site's QUAL is not computed but read: `QUALapprox`, the sum of `PL[0]` over the samples that
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
//!
//! # Attributes are a list, not a map
//!
//! The Java's `VariantContext` attributes are a `Map`; here they are an ordered list of
//! `(key, value)` pairs, so `attribute`, `put` and `remove` below are `getAttribute`,
//! `attribute(key, value)` (replace or append) and `rmAttribute`.
//!
//! # Reading this file without reading Rust
//!
//! - `Result<T, E>` is "a value or an error"; [`EngineError`] stands for the exception the Java
//!   throws, and `?` after a call returns its error to the caller at once;
//! - `Option<T>` is "a value or nothing", Java's nullable reference checked by the compiler;
//!   `let Some(x) = y else { .. }` takes the value or the early exit;
//! - `&[T]` is a borrowed, read-only view of a list and `&mut T` a writable borrow; the caller
//!   keeps ownership. `.clone()` makes the copy a Java reference assignment would share;
//! - `.iter().map(..).collect()` and friends are Java's stream operations, `|x| ...` a lambda;
//! - `x as f64` converts a count to a `double`; `cond.then(|| v)` and `cond.then_some(v)` are
//!   `cond ? Optional.of(v) : Optional.empty()`.

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
/// See above.
const INDEL_HETEROZYGOSITY: f64 = 1.25e-4;
/// `GenotypeCalculationArgumentCollection.DEFAULT_STANDARD_CONFIDENCE_FOR_CALLING`.
const DEFAULT_STANDARD_CONFIDENCE_FOR_CALLING: f64 = 30.0;

/// The allele-specific reducible annotations the engine finds by reflection, with the one key
/// each finalizes and trims.
///
/// The Java collects every `InfoFieldAnnotation`, `AS_StrandBiasTest` and `AS_RankSumTest`
/// subtype under `annotator.allelespecific` and keeps those that are both `AS_StandardAnnotation`
/// and `ReducibleAnnotation`; this is that list written out, the key being each one's
/// `getKeyNames().get(0)`.
const AS_REDUCIBLE: &[(&str, &str)] = &[
    ("AS_QualByDepth", "AS_QD"),
    ("AS_RMSMappingQuality", "AS_MQ"),
    ("AS_FisherStrand", "AS_FS"),
    ("AS_StrandOddsRatio", "AS_SOR"),
    ("AS_BaseQualityRankSumTest", "AS_BaseQRankSum"),
    ("AS_MappingQualityRankSumTest", "AS_MQRankSum"),
    ("AS_ReadPosRankSumTest", "AS_ReadPosRankSum"),
];

/// `Allele.NON_REF_ALLELE`.
fn non_ref() -> Allele {
    Allele::create(b"<NON_REF>", false).expect("the symbolic <NON_REF>")
}

/// The `*` allele, `Allele.SPAN_DEL`.
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

/// The `NumberFormatException` `Integer.parseInt` throws on `text`.
fn number_format(text: &str) -> EngineError {
    EngineError::Runtime {
        class: "java.lang.NumberFormatException".to_string(),
        message: format!("For input string: \"{text}\""),
    }
}

/// `getAttributeAsInt(key, default)`.
///
/// An absent or missing value is the default, an integer is itself, anything else is parsed
/// from its text.
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
///
/// A list is read element by element, a single value split at commas, and `.` reads as the
/// default 0. The last `.collect()` gathers a list of `Result`s into one `Result` of a list,
/// failing on the first element that failed.
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
///
/// `RAW_MQandDP` is the sum of squared mapping qualities and the depth; MQ is written only when
/// the depth is positive, and both raw keys are removed. A site carrying only the deprecated
/// `MQ_DP` is refused as a port limitation rather than computed.
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
///
/// The Java splits with `asQuals.split(AnnotationUtils.ALLELE_SPECIFIC_SPLIT_REGEX, -1)`, "allow
/// for empty tokens at the end": the `-1` keeps trailing empty entries, as Rust's `split` does,
/// so the count check sees every `|`. The reference's entry is skipped and empty entries count
/// as nothing.
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
///
/// The Java reads `variant.getAttribute(key).toString()`, which for a list is `[a, b]`, then
/// `decodeAnyASList` strips the brackets and splits at commas; a list is therefore rendered the
/// same way here first. An index past the entries becomes `.`, and the entries are joined with
/// commas (`encodeStringList`).
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

/// `getSBFieldAsIntArray`: the genotype's `SB` as four integers, or the `IllegalStateException`
/// the Java throws when one does not parse.
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

/// `GnarlyGenotyperEngine`: `keepAllSites`, `maxAltAllelesToOutput` and `stripASAnnotations`.
/// `emitPls` is always true here, the three-argument constructor's value.
pub struct GnarlyEngine {
    pub keep_all_sites: bool,
    pub max_alt_alleles_to_output: usize,
    pub strip_as_annotations: bool,
}

impl GnarlyEngine {
    /// `finalizeGenotype(variant, annotationDBBuilder)`: the site, or `None` under the floor, and
    /// what the annotation database gets when one is written.
    ///
    /// The Java mutates the `annotationDBBuilder` it is passed; here the database copy is built
    /// alongside (`with_database`) and returned as the second half of the pair. `random` is the
    /// generator the rank-sum finalizers draw from.
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
        // `allele != Allele.SPAN_DEL && allele.length() == variant.getReference().length()`:
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
        // `INDEL_QUAL_THRESHOLD` or `SNP_QUAL_THRESHOLD`, computed the same way: the default
        // calling confidence less ten times `Math.log10` of the heterozygosity.
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
        // `annotationsToBeModified`: each annotation's
        // `finalizeRawData(vcfBuilder.make(), variant)` reads the MQ-finalized site, and its
        // results are collected separately and set at the end, as the Java's map is.
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

        // "tolerate lack of VarDP annotation": QD and the QUAL (as `log10PError`) only with it.
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

        // `int[] SBsum = {0,0,0,0}`, then the target alleles: `<NON_REF>`, when present, is
        // assumed last and dropped (`subList(0, size - 1)`).
        let mut sb_sum = [0i32; 4];
        let remove_non_ref = variant.alleles.contains(&non_ref());
        let targets: Vec<Allele> = if remove_non_ref {
            variant.alleles[..variant.alleles.len() - 1].to_vec()
        } else {
            variant.alleles.clone()
        };
        let mut counts: Vec<(Allele, i64)> =
            targets.iter().map(|allele| (allele.clone(), 0)).collect();
        // `variant.hasAttribute(RAW_GENOTYPE_COUNT_KEY) ? null : rawGenotypeCounts`: the counts are
        // tallied only when the site does not carry them.
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
            // `targetAlleleCounts.size() == 1 ? targetAlleleCounts.get(0) : targetAlleleCounts`:
            // a single alternate is written as a scalar, several as a list; AF likewise.
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
            // No genotypes: the site's own `SB_TABLE`, and the database copies AC, AF and AN as the
            // site had them (null where absent).
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
            // `Math.max(numCalledAlleles / 2 - gtCounts.get(1) - gtCounts.get(2), 0)`: integer
            // division, and "homRefs don't get counted properly because ref blocks aren't
            // annotated".
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

        // FS and SOR from the summed table, `StrandBiasTest.decodeSBBS(SBsum)`.
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

        // Trim `<NON_REF>` out of each AS value; then, only after every finalized value is in,
        // remove the raw keys: "since AS_FS and AS_SOR share the same raw key, we have to wait".
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
        // `AS_QualByDepth.finalizeRawGVCFVarDPValues(.., targetAlleles.size())`: the reference's
        // and `<NON_REF>`'s entries dropped, or missing when the count is not one more than the
        // targets.
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
    ///
    /// Each genotype: a call on `<NON_REF>`, or a GenomicsDB no-call with neither PL nor AD
    /// (`isGDBnoCall`), becomes a no-call without GQ; AD is trimmed to the targets; PLs are
    /// trimmed to the target genotype count (`trimPLs`, whose `System.arraycopy` throws on a
    /// short array), GQ set from them and the call made from them; a GQ-0 genotype without PLs
    /// is called with null likelihoods. Then SB is summed, the called alleles counted, and the
    /// raw genotype counts tallied by the number of non-reference alleles in the INPUT call.
    ///
    /// `mut raw_genotype_counts: Option<&mut [i64; 3]>` is the Java's nullable `int[]`;
    /// `as_deref_mut` borrows the array out of the `Option` for one genotype at a time.
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
        // `numConcreteAlts = maximumAlleleCount - 2`, "-1 for NON_REF and -1 for ref".
        // `saturating_sub` stops at zero where the Java's `int` would go negative.
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
                // `trimADs`: `take` truncates, where the Java's `System.arraycopy` would throw on
                // an AD shorter than the targets.
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
///
/// One pass keeping the two smallest, the Java's loop line for line; `wrapping_sub` is Java's
/// `int` subtraction.
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

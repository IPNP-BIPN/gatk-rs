//! Between htsjdk's `VariantContext` and the reduced record the ported statics read.
//!
//! # Why there are two models at all
//!
//! `htsjdk_vcf::variant::VariantContext` is the file: an id, a QUAL, applied-or-not filters, typed
//! attributes and genotypes whose alleles are alleles. `gatk_engine::variant_context_utils::Variant`
//! is what the ported GATK statics were written against, and it is deliberately smaller -- "as much
//! of a `VariantContext` as these functions read and write", which is the contig, the span, the
//! alleles, the genotypes and the attributes as strings. Neither is wrong. A tool that only decides
//! whether to keep a record needs the small one; a tool that WRITES the record back needs the file.
//!
//! `SelectVariants` is the first ported tool that needs both, so this is where they meet.
//!
//! # The rule that makes the trip lossless
//!
//! Going down is a rendering: every attribute becomes the string the encoder would have written.
//! Going back up is NOT the inverse rendering, because it cannot be: a decoded INFO flag is
//! `Value::Bool(true)` and prints as `KEY`, while the string it renders to is empty and would print
//! as `KEY=`. So the trip back is a DIFF against the record it came from. A key whose rendered
//! string is unchanged keeps the original `Value` -- flag, list and all -- and only a key the
//! pipeline actually rewrote becomes a `Value::Str`. The three the pipeline rewrites are AC, AN and
//! AF, which `calculateChromosomeCounts` sets as text in the reference as well.
//!
//! Genotype alleles make the same trip through indices, because that is what
//! `gatk_engine::subset_alleles` works in: a call is a position in the record's allele list and a
//! no-call is `None`. The list they index is the record's list AFTER the subset, which is why the
//! write-back reads the alleles it is given rather than the ones the original carried.

use std::collections::HashMap;

use gatk_engine::subset_alleles::Genotype as EngineGenotype;
use gatk_engine::variant_context_utils::{Allele as EngineAllele, Variant};
use gatk_tools::select_variants::{FilterRecord, Record};
use htsjdk_vcf::allele::Allele as VcfAllele;
use htsjdk_vcf::variant::{Genotype as VcfGenotype, Value, VariantContext};

/// One record in both models, with the file's own copy kept for the write-back's diff.
pub struct Bridged {
    /// The record as it was decoded. Nothing writes to it; it is the baseline the diff is against.
    pub original: VariantContext,
    pub record: Record,
    pub filter_record: FilterRecord,
}

/// What `Value.format()` would write, which is what the string model holds.
///
/// `None` is the one thing a string model cannot represent, and only a `false` flag produces it;
/// such a field is dropped by the encoder, so the empty string here is a value no key reaches.
fn rendered(value: &Value) -> String {
    value.format().unwrap_or_default()
}

/// The bases the reduced model measures: the DISPLAY string, so that `*` stays one base and a
/// symbolic allele keeps its angle brackets rather than becoming empty.
fn engine_allele(allele: &VcfAllele) -> EngineAllele {
    EngineAllele {
        bases: allele.display_string().into_bytes(),
        is_reference: allele.is_reference(),
    }
}

/// The index of a genotype's allele in the record's list, which is how the reference stores a call.
///
/// A no-call is `None`. An allele the record does not declare is `None` as well rather than an
/// error: htsjdk's own `GenotypesContext` allows it, and the reference reads such a call as no
/// call when it counts chromosomes.
fn allele_index(alleles: &[VcfAllele], allele: &VcfAllele) -> Option<usize> {
    if allele.is_no_call() {
        return None;
    }
    alleles.iter().position(|candidate| candidate == allele)
}

fn engine_genotype(genotype: &VcfGenotype, alleles: &[VcfAllele]) -> EngineGenotype {
    let mut attributes: Vec<(String, String)> = Vec::new();
    // `FT` is a field of its own in htsjdk and an attribute in the reduced model, and it is the
    // attribute `setFilteredGenotypeToNocall` reads, so it goes first.
    if let Some(filters) = &genotype.filters {
        attributes.push(("FT".to_string(), filters.clone()));
    }
    for (key, value) in &genotype.extended {
        attributes.push((key.clone(), rendered(value)));
    }
    EngineGenotype {
        alleles: genotype
            .alleles
            .iter()
            .map(|allele| allele_index(alleles, allele))
            .collect(),
        pl: genotype.pl.clone(),
        gq: genotype.gq,
        ad: genotype.ad.clone(),
        dp: genotype.dp,
        attributes,
    }
}

/// The record as the ported statics read it, with the file's copy kept beside it.
pub fn to_engine(vc: &VariantContext) -> Bridged {
    let genotypes: Vec<EngineGenotype> = vc
        .genotypes
        .iter()
        .map(|genotype| engine_genotype(genotype, &vc.alleles))
        .collect();
    let variant = Variant {
        contig: vc.contig.clone(),
        start: vc.start as i32,
        stop: vc.stop as i32,
        alleles: vc.alleles.iter().map(engine_allele).collect(),
        genotypes,
        attributes: vc
            .attributes
            .iter()
            .map(|(key, value)| (key.clone(), rendered(value)))
            .collect(),
    };
    // `VariantContextUtils.match`'s view of the record: what a JEXL expression can name, which is
    // the id, the filters, the INFO map and one map per genotype.
    let filter_record = FilterRecord {
        id: vc.id.clone(),
        filters: vc.filters.clone().unwrap_or_default(),
        info: vc
            .attributes
            .iter()
            .map(|(key, value)| (key.clone(), rendered(value)))
            .collect::<HashMap<String, String>>(),
        genotype_fields: vc
            .genotypes
            .iter()
            .map(|genotype| {
                let mut fields: HashMap<String, String> = genotype
                    .extended
                    .iter()
                    .map(|(key, value)| (key.clone(), rendered(value)))
                    .collect();
                if let Some(filters) = &genotype.filters {
                    fields.insert("FT".to_string(), filters.clone());
                }
                if let Some(gq) = genotype.gq {
                    fields.insert("GQ".to_string(), gq.to_string());
                }
                if let Some(dp) = genotype.dp {
                    fields.insert("DP".to_string(), dp.to_string());
                }
                fields
            })
            .collect(),
    };
    Bridged {
        original: vc.clone(),
        record: Record {
            variant,
            samples: vc
                .genotypes
                .iter()
                .map(|genotype| genotype.sample_name.clone())
                .collect(),
        },
        filter_record,
    }
}

/// An engine allele back as the file's, by identity where the record already had it.
///
/// Identity rather than reconstruction, because `Allele::create` is not a total inverse of
/// `display_string`: a symbolic allele's brackets and the star allele go through it, and htsjdk's
/// own equality counts the reference flag. An allele the subset produced that the original did not
/// carry cannot be looked up, and that one is built.
fn vcf_allele(original: &[VcfAllele], allele: &EngineAllele) -> VcfAllele {
    let text = String::from_utf8_lossy(&allele.bases).to_string();
    original
        .iter()
        .find(|candidate| {
            candidate.display_string() == text && candidate.is_reference() == allele.is_reference
        })
        .cloned()
        .unwrap_or_else(|| {
            VcfAllele::from_str(&text, allele.is_reference).unwrap_or_else(|_| VcfAllele::no_call())
        })
}

/// The attributes of the written record: the original's `Value` wherever the string is unchanged,
/// and a `Value::Str` wherever the pipeline rewrote it.
///
/// Order is the engine's, which is the original's order with anything added at the end, because
/// that is where `set_attribute` puts a key the record did not have. The encoder sorts, so the
/// order does not reach the file; it is kept so that a divergence can be read.
fn written_attributes(
    original: &[(String, Value)],
    engine: &[(String, String)],
) -> Vec<(String, Value)> {
    engine
        .iter()
        .map(|(key, text)| {
            let held = original
                .iter()
                .find(|(name, _)| name == key)
                .filter(|(_, value)| rendered(value) == *text);
            match held {
                Some((_, value)) => (key.clone(), value.clone()),
                None => (key.clone(), Value::Str(text.clone())),
            }
        })
        .collect()
}

fn written_genotype(
    original: Option<&VcfGenotype>,
    sample_name: &str,
    genotype: &EngineGenotype,
    alleles: &[VcfAllele],
) -> VcfGenotype {
    let filters = genotype
        .attributes
        .iter()
        .find(|(key, _)| key == "FT")
        .map(|(_, value)| value.clone());
    let extended: Vec<(String, String)> = genotype
        .attributes
        .iter()
        .filter(|(key, _)| key != "FT")
        .cloned()
        .collect();
    VcfGenotype {
        sample_name: sample_name.to_string(),
        alleles: genotype
            .alleles
            .iter()
            .map(|call| match call {
                Some(index) => alleles
                    .get(*index)
                    .cloned()
                    .unwrap_or_else(VcfAllele::no_call),
                None => VcfAllele::no_call(),
            })
            .collect(),
        phased: original.is_some_and(|genotype| genotype.phased),
        gq: genotype.gq,
        dp: genotype.dp,
        ad: genotype.ad.clone(),
        pl: genotype.pl.clone(),
        filters,
        extended: written_attributes(
            original
                .map(|genotype| genotype.extended.as_slice())
                .unwrap_or(&[]),
            &extended,
        ),
    }
}

/// The record to write: the one that was decoded, with what the pipeline changed applied to it.
///
/// The id, the QUAL and the applied-or-not distinction of the filters are carried across
/// untouched, because nothing in the pipeline can reach them: `SelectVariants` selects records and
/// rewrites alleles, and a record it keeps is the record the file had.
pub fn from_engine(original: &VariantContext, record: &Record) -> VariantContext {
    let alleles: Vec<VcfAllele> = record
        .variant
        .alleles
        .iter()
        .map(|allele| vcf_allele(&original.alleles, allele))
        .collect();
    let genotypes = record
        .variant
        .genotypes
        .iter()
        .zip(record.samples.iter())
        .map(|(genotype, sample_name)| {
            let held = original
                .genotypes
                .iter()
                .find(|candidate| candidate.sample_name == *sample_name);
            written_genotype(held, sample_name, genotype, &alleles)
        })
        .collect();
    VariantContext {
        contig: record.variant.contig.clone(),
        start: record.variant.start as i64,
        stop: record.variant.stop as i64,
        id: original.id.clone(),
        alleles,
        log10_p_error: original.log10_p_error,
        filters: original.filters.clone(),
        attributes: written_attributes(&original.attributes, &record.variant.attributes),
        genotypes,
    }
}

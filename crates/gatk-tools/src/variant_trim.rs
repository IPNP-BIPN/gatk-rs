//! `GATKVariantContextUtils.reverseTrimAlleles` over a whole record.
//!
//! Ported from `org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils` (GATK
//! 4.6.2.0). The trim itself is [`gatk_engine::variant_context_utils::trim_alleles`]; this carries
//! it over to an htsjdk record, whose genotypes hold their alleles by value and so have each one
//! replaced by the allele at its index. A symbolic allele or a spanning deletion is never cut.

use gatk_engine::variant_context_utils::{
    trim_alleles, Allele as EngineAllele, Variant as EngineVariant,
};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::genotypes_context::GenotypesContext;
use htsjdk_vcf::variant::{Genotype, VariantContext};

fn is_span_del(allele: &Allele) -> bool {
    !allele.is_reference() && allele.display_string() == htsjdk_vcf::allele::SPAN_DEL_STRING
}

/// `GATKVariantContextUtils.reverseTrimAlleles`: the bases every allele shares at its end, cut.
/// Genotypes hold alleles by value, so each is replaced by the allele at its index.
pub fn reverse_trim_alleles(vc: &VariantContext) -> Result<VariantContext, String> {
    let engine = EngineVariant {
        contig: vc.contig.clone(),
        start: vc.start as i32,
        stop: vc.stop as i32,
        alleles: vc
            .alleles
            .iter()
            .map(|allele| {
                EngineAllele::new(allele.display_string().as_bytes(), allele.is_reference())
            })
            .collect(),
        genotypes: Vec::new(),
        attributes: Vec::new(),
    };
    let trimmed = trim_alleles(&engine, false, true).map_err(|error| format!("{error:?}"))?;
    if trimmed.alleles == engine.alleles && trimmed.stop == engine.stop {
        return Ok(vc.clone());
    }
    let alleles: Vec<Allele> = vc
        .alleles
        .iter()
        .zip(&trimmed.alleles)
        .map(|(original, cut)| {
            if original.is_symbolic() || is_span_del(original) {
                original.clone()
            } else {
                Allele::create(&cut.bases, original.is_reference())
                    .unwrap_or_else(|_| original.clone())
            }
        })
        .collect();
    let genotypes: Vec<Genotype> = vc
        .genotypes
        .iter()
        .map(|genotype| {
            let mut genotype = genotype.clone();
            genotype.alleles = genotype
                .alleles
                .iter()
                .map(|allele| match vc.alleles.iter().position(|a| a == allele) {
                    Some(index) => alleles[index].clone(),
                    None => allele.clone(),
                })
                .collect();
            genotype
        })
        .collect();
    let mut out = vc.clone();
    out.stop = i64::from(trimmed.stop);
    out.alleles = alleles;
    out.genotypes = GenotypesContext::new(genotypes);
    Ok(out)
}

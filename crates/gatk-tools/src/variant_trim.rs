//! `GATKVariantContextUtils.reverseTrimAlleles` over a whole record.
//!
//! Ported from `org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils` (GATK
//! 4.7.0.0; `GATKVariantContextUtils.java` is unchanged since 4.6.2.0). The trim itself is
//! [`gatk_engine::variant_context_utils::trim_alleles`]; this carries it over to an htsjdk
//! record, whose genotypes hold their alleles by value and so have each one replaced by the
//! allele at its index. A symbolic allele or a spanning deletion is never cut.
//!
//! # Reading this file without reading Rust
//!
//! - `Result<T, String>` is "a value or an error message"; `?` after a call returns the error
//!   to the caller at once, as a Java method would let the exception propagate.
//! - `.iter().map(|x| ...).collect()` is Java's `stream().map(...).collect(...)`; `|x| ...` is
//!   a lambda. `.zip(other)` walks two lists in step, as pairs.
//! - `as i32` converts a 64-bit position to Java's `int` width; it truncates rather than
//!   throwing, which is harmless for a VCF position.

use gatk_engine::variant_context_utils::{
    trim_alleles, Allele as EngineAllele, Variant as EngineVariant,
};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::genotypes_context::GenotypesContext;
use htsjdk_vcf::variant::{Genotype, VariantContext};

/// `a.equals(Allele.SPAN_DEL)`: the `*` alternate, which the Java trim leaves untouched.
fn is_span_del(allele: &Allele) -> bool {
    !allele.is_reference() && allele.display_string() == htsjdk_vcf::allele::SPAN_DEL_STRING
}

/// `GATKVariantContextUtils.reverseTrimAlleles`: the bases every allele shares at its end, cut.
/// Genotypes hold alleles by value, so each is replaced by the allele at its index.
///
/// The Java is one line, `return trimAlleles(inputVC, false, true);`, and the work happens in
/// the three-argument `trimAlleles`, which builds each new allele and then remaps the genotypes:
///
/// ```java
/// if (a.isSymbolic() || a.equals(Allele.SPAN_DEL)) {
///     alleles.add(a);
///     originalToTrimmedAlleleMap.put(a, a);
/// } else {
///     // get bases for current allele and create a new one with trimmed bases
///     final byte[] newBases = Arrays.copyOfRange(a.getBases(), fwdTrimEnd+1, a.length()-revTrim);
///     final Allele trimmedAllele = Allele.create(newBases, a.isReference());
/// ```
///
/// The engine's `trim_alleles` computes the cut on a copy of the record without genotypes; this
/// function then applies it to the htsjdk record. When the engine returns the alleles and stop
/// unchanged, the input is returned as it was, as Java returns `inputVC`.
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
    // Forward trim off, reverse trim on: the `(inputVC, false, true)` of `reverseTrimAlleles`.
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
    // `AlleleMapper.remap`: each genotype allele is found by equality among the record's
    // original alleles and replaced by the trimmed allele at the same index.
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

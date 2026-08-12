//! `RemoveNearbyIndels`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.validation.RemoveNearbyIndels` (GATK 4.6.2.0).
//!
//! The first member of the variant-transform archetype: a `VariantWalker` that buffers one indel
//! at a time and drops any **pair** of indels closer than a spacing, keeping whatever non-indels
//! sat between them.
//!
//! # The buffer remembers an indel it threw away
//!
//! ```java
//! } else if (nearby(lastIndel, vc)) {
//!     if (vc.isIndel()) {
//!         emitAllNonIndels(); // throw out {@code lastIndel} and {@code vc}
//! ...
//! lastIndel = vc.isIndel() ? vc : lastIndel;
//! ```
//!
//! The last line runs after the branch that discarded `vc`, so the next indel is measured against
//! one that never reached the output. Three indels in a row lose all three, and with a wide enough
//! spacing a lone indel hundreds of bases past a discarded one is discarded with it.
//!
//! # The last indel survives on an identity, not an equality
//!
//! ```java
//! buffer.stream().filter(vc -> !(vc.isIndel() && nearby(lastIndel, vc)) || vc == lastIndel)
//! ```
//!
//! At the end of the traversal the buffered indel **is** `lastIndel`, and `nearby(lastIndel,
//! lastIndel)` is true for any positive spacing, since `start - end` is at most zero. Only the
//! `vc == lastIndel` reference comparison keeps it. This port buffers indices into the input, so
//! that identity is index equality and not a comparison of contents: two records with the same
//! contents at the same place would not be the same variant here either.
//!
//! # `isIndel` is a type, and a mixed site is not that type
//!
//! `getType()` compares every alternate against the reference and returns `MIXED` as soon as two
//! of them disagree. A record carrying a snp and an insertion side by side is therefore **not** an
//! indel, is never buffered, and cannot pair with anything, while a record whose alternates are
//! all a different length from the reference is an `INDEL` and buffers like one.

use htsjdk_vcf::variant::VariantContext;

/// `VariantContext.Type`, as far as this tool needs it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantType {
    NoVariation,
    Snp,
    Mnp,
    Indel,
    Symbolic,
    Mixed,
}

/// `typeOfBiallelicVariant(ref, allele)`: length decides, and only length.
///
/// The comment in htsjdk is explicit that a prefix test used to be here and was wrong: `REF=CTTA`
/// with `ALT=C,CT,CA` is an indel, and checking prefixes made it mixed.
fn type_of_biallelic(
    reference: &htsjdk_vcf::allele::Allele,
    allele: &htsjdk_vcf::allele::Allele,
) -> VariantType {
    if allele.is_symbolic() {
        return VariantType::Symbolic;
    }
    if reference.len() == allele.len() {
        if allele.len() == 1 {
            VariantType::Snp
        } else {
            VariantType::Mnp
        }
    } else {
        VariantType::Indel
    }
}

/// `determineType()`: one allele is no variation, and disagreeing alternates are mixed.
pub fn variant_type(variant: &VariantContext) -> VariantType {
    if variant.alleles.len() < 2 {
        return VariantType::NoVariation;
    }
    let reference = variant.reference();
    let mut kind: Option<VariantType> = None;
    for allele in variant.alternate_alleles() {
        let this = type_of_biallelic(reference, allele);
        match kind {
            None => kind = Some(this),
            Some(seen) if seen != this => return VariantType::Mixed,
            Some(_) => {}
        }
    }
    kind.unwrap_or(VariantType::NoVariation)
}

/// `isIndel()`, which is `getType() == INDEL` and nothing looser.
pub fn is_indel(variant: &VariantContext) -> bool {
    variant_type(variant) == VariantType::Indel
}

/// The tool: the indices of the variants that reach the output, in the order they are written.
///
/// Indices rather than records, because the reference keeps `lastIndel` by identity and compares
/// it with `==`.
pub fn remove_nearby_indels(variants: &[VariantContext], min_indel_spacing: i32) -> Vec<usize> {
    let mut buffer = VariantBuffer {
        variants,
        min_indel_spacing,
        last_indel: None,
        buffer: Vec::new(),
        emitted: Vec::new(),
    };
    for index in 0..variants.len() {
        buffer.add(index);
    }
    buffer.emit_remaining();
    buffer.emitted
}

struct VariantBuffer<'a> {
    variants: &'a [VariantContext],
    min_indel_spacing: i32,
    last_indel: Option<usize>,
    /// "INVARIANT: this buffer will contain at most one indel", says the reference.
    buffer: Vec<usize>,
    emitted: Vec<usize>,
}

impl VariantBuffer<'_> {
    fn add(&mut self, index: usize) {
        let indel = is_indel(&self.variants[index]);
        if self.last_indel.is_none() && indel {
            // Only ever reached once per traversal: `lastIndel` never returns to null.
            self.buffer.push(index);
            self.last_indel = Some(index);
        } else if self.nearby(self.last_indel, index) {
            if indel {
                // Throws out `lastIndel`, which is in the buffer, and `vc`, which is never added.
                self.emit_all_non_indels();
            } else {
                self.buffer.push(index);
            }
        } else {
            self.emit_all_variants();
            self.buffer.push(index);
        }

        // After the branch above, so an indel that was just discarded becomes the one the next
        // indel is measured against.
        if indel {
            self.last_indel = Some(index);
        }
    }

    /// `nearby(left, right)`: same contig, and END to START, strictly below the spacing.
    fn nearby(&self, left: Option<usize>, right: usize) -> bool {
        match left {
            None => false,
            Some(left) => {
                let (left, right) = (&self.variants[left], &self.variants[right]);
                left.contig == right.contig
                    && (right.start - left.stop) < i64::from(self.min_indel_spacing)
            }
        }
    }

    fn emit_all_non_indels(&mut self) {
        for index in std::mem::take(&mut self.buffer) {
            if !is_indel(&self.variants[index]) {
                self.emitted.push(index);
            }
        }
    }

    fn emit_all_variants(&mut self) {
        self.emitted.append(&mut self.buffer);
    }

    /// The end of the traversal, where the identity comparison is what keeps a trailing indel.
    fn emit_remaining(&mut self) {
        for index in std::mem::take(&mut self.buffer) {
            let keep = !(is_indel(&self.variants[index]) && self.nearby(self.last_indel, index))
                || self.last_indel == Some(index);
            if keep {
                self.emitted.push(index);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_vcf::allele::Allele;

    fn variant(contig: &str, start: i64, reference: &str, alternates: &[&str]) -> VariantContext {
        let mut alleles = vec![Allele::create(reference.as_bytes(), true).expect("a reference")];
        for alternate in alternates {
            alleles.push(Allele::create(alternate.as_bytes(), false).expect("an alternate"));
        }
        VariantContext::new(contig, start, alleles)
    }

    #[test]
    fn a_mixed_site_is_not_an_indel() {
        assert!(is_indel(&variant("chr1", 100, "A", &["AGG"])));
        assert!(is_indel(&variant("chr1", 100, "ACG", &["A", "T"])));
        // A snp and an insertion disagree, so the type is MIXED.
        assert!(!is_indel(&variant("chr1", 100, "A", &["C", "AGG"])));
        assert!(!is_indel(&variant("chr1", 100, "A", &["C", "G"])));
    }

    #[test]
    fn the_distance_is_measured_from_the_end_and_the_test_is_strict() {
        // A two-base reference ends at 101, so an insertion at 106 is five away.
        let variants = vec![
            variant("chr1", 100, "AC", &["A"]),
            variant("chr1", 106, "A", &["AGG"]),
        ];
        assert_eq!(remove_nearby_indels(&variants, 5), vec![0, 1]);
        assert!(remove_nearby_indels(&variants, 6).is_empty());
    }

    #[test]
    fn a_discarded_indel_is_still_what_the_next_one_is_measured_against() {
        let variants = vec![
            variant("chr1", 100, "A", &["AGG"]),
            variant("chr1", 105, "A", &["AGG"]),
            variant("chr1", 110, "A", &["AGG"]),
            variant("chr1", 500, "A", &["C"]),
        ];
        // All three indels go, though the third is only near the second, which was thrown away.
        assert_eq!(remove_nearby_indels(&variants, 20), vec![3]);
    }

    #[test]
    fn the_non_indels_between_a_discarded_pair_are_kept() {
        let variants = vec![
            variant("chr1", 200, "ACG", &["A"]),
            variant("chr1", 205, "A", &["C"]),
            variant("chr1", 210, "A", &["AGG"]),
        ];
        assert_eq!(remove_nearby_indels(&variants, 20), vec![1]);
    }

    #[test]
    fn a_trailing_indel_survives_on_the_identity() {
        let variants = vec![
            variant("chr1", 100, "A", &["C"]),
            variant("chr1", 500, "A", &["AGG"]),
        ];
        assert_eq!(remove_nearby_indels(&variants, 20), vec![0, 1]);
    }

    #[test]
    fn a_contig_boundary_is_never_nearby() {
        let variants = vec![
            variant("chr1", 1000, "A", &["AGG"]),
            variant("chr2", 1, "A", &["AGG"]),
        ];
        assert_eq!(remove_nearby_indels(&variants, 20), vec![0, 1]);
    }
}

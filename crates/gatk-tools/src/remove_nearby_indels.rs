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

/// What the buffer needs of a record: where it is and whether it is an indel. The records
/// themselves stay with the caller, which [`IndelSpacing`] names by their index in the traversal.
#[derive(Debug, Clone)]
pub struct Site {
    pub index: usize,
    pub indel: bool,
    pub contig: String,
    pub start: i64,
    pub stop: i64,
}

impl Site {
    pub fn of(index: usize, variant: &VariantContext) -> Site {
        Site {
            index,
            indel: is_indel(variant),
            contig: variant.contig.clone(),
            start: variant.start,
            stop: variant.stop,
        }
    }
}

/// `RemoveNearbyIndels`' `VariantBuffer`, fed one record at a time. The indices it emits are
/// increasing, so a caller can write each record as it is emitted and drop any record
/// [`buffered`](Self::buffered) no longer names.
pub struct IndelSpacing {
    min_indel_spacing: i32,
    last_indel: Option<Site>,
    /// "INVARIANT: this buffer will contain at most one indel", says the reference.
    buffer: Vec<Site>,
}

impl IndelSpacing {
    pub fn new(min_indel_spacing: i32) -> IndelSpacing {
        IndelSpacing {
            min_indel_spacing,
            last_indel: None,
            buffer: Vec::new(),
        }
    }

    /// The indices still held, which a later call may yet emit.
    pub fn buffered(&self) -> impl Iterator<Item = usize> + '_ {
        self.buffer.iter().map(|site| site.index)
    }

    pub fn add(&mut self, site: Site, emitted: &mut Vec<usize>) {
        let indel = site.indel;
        if self.last_indel.is_none() && indel {
            // Only ever reached once per traversal: `lastIndel` never returns to null.
            self.buffer.push(site.clone());
            self.last_indel = Some(site);
            return;
        }
        if self.nearby(self.last_indel.as_ref(), &site) {
            if indel {
                // Throws out `lastIndel`, which is in the buffer, and `vc`, which is never added.
                for held in std::mem::take(&mut self.buffer) {
                    if !held.indel {
                        emitted.push(held.index);
                    }
                }
            } else {
                self.buffer.push(site.clone());
            }
        } else {
            emitted.extend(self.buffer.drain(..).map(|held| held.index));
            self.buffer.push(site.clone());
        }

        // After the branch above, so an indel that was just discarded becomes the one the next
        // indel is measured against.
        if indel {
            self.last_indel = Some(site);
        }
    }

    /// `nearby(left, right)`: same contig, and END to START, strictly below the spacing.
    fn nearby(&self, left: Option<&Site>, right: &Site) -> bool {
        match left {
            None => false,
            Some(left) => {
                left.contig == right.contig
                    && (right.start - left.stop) < i64::from(self.min_indel_spacing)
            }
        }
    }

    /// The end of the traversal, where the identity comparison is what keeps a trailing indel.
    pub fn finish(mut self, emitted: &mut Vec<usize>) {
        let last = self.last_indel.take();
        for held in std::mem::take(&mut self.buffer) {
            let keep = !(held.indel && self.nearby(last.as_ref(), &held))
                || last.as_ref().map(|site| site.index) == Some(held.index);
            if keep {
                emitted.push(held.index);
            }
        }
    }
}

/// The tool: the indices of the variants that reach the output, in the order they are written.
///
/// Indices rather than records, because the reference keeps `lastIndel` by identity and compares
/// it with `==`.
pub fn remove_nearby_indels(variants: &[VariantContext], min_indel_spacing: i32) -> Vec<usize> {
    let mut spacing = IndelSpacing::new(min_indel_spacing);
    let mut emitted = Vec::new();
    for (index, variant) in variants.iter().enumerate() {
        spacing.add(Site::of(index, variant), &mut emitted);
    }
    spacing.finish(&mut emitted);
    emitted
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

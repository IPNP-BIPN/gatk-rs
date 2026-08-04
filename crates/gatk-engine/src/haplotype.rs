//! Ported from `org.broadinstitute.hellbender.utils.haplotype.Haplotype` (GATK 4.6.2.0), as far as
//! the allele axis of a likelihood matrix needs it.
//!
//! ```java
//! public class Haplotype extends SimpleAllele implements Locatable {
//!     public Haplotype(final byte[] bases, final boolean isRef) {
//!         super(Arrays.copyOf(bases, bases.length), isRef);
//!     }
//! }
//! ```
//!
//! A haplotype **is** an allele in the reference, and that is why `AlleleLikelihoods` can be
//! instantiated over it: the haplotype-typed matrix `HaplotypeFilteringAnnotation` reads is the same
//! class as the allele-typed one every other annotation reads. So the port keeps the `SimpleAllele`
//! it is built on as a field rather than re-deriving base handling: the uppercasing, the
//! symbolic-allele test and the refusals are that constructor's, and a second implementation of them
//! would be a second chance to disagree.
//!
//! # What this slice is, and what it is not
//!
//! Only the parts that are observable through the allele axis are here: the bases, the reference
//! flag, and the uniqueness value that equality reads. The assembly-side state (`cigar`,
//! `alignmentStartHapwrtRef`, `score`, `isCollapsed`, `kmerSize`, the genome location, and above all
//! `EventMap`) is deliberately absent. `getEventMap()` is the assembly event model, which is what
//! blocks `AssemblyComplexity` and puts it in Milestone G3; inventing an empty version of it here
//! would let a G3 annotation compile against something that has never been compared to anything.
//!
//! # `equals` and `hashCode` disagree about what a haplotype is
//!
//! ```java
//! public boolean equals(final Object h) {
//!     return h instanceof Haplotype
//!             && getUniquenessValue() == ((Haplotype) h).getUniquenessValue()
//!             && isReference() == ((Haplotype) h).isReference()
//!             && Arrays.equals(getBases(), ((Haplotype) h).getBases());
//! }
//! public int hashCode() { return Arrays.hashCode(getBases()); }
//! ```
//!
//! Equality is three fields; the hash is one of them. That is legal (equal objects still hash
//! equal) and it is load-bearing twice over.
//!
//! First, it widens what a haplotype list can hold. `IndexedSet` drops a duplicate, and duplicate
//! means `equals`: two haplotypes with the same bases and different uniqueness values are **two**
//! entries, where two alleles with the same bases and the same reference flag are one. So
//! `alleles().size()` on a haplotype list is not "the number of distinct base strings", and
//! `HaplotypeFilteringAnnotation` reports that size as `ASSEMBLED_HAPS`.
//!
//! Second, it collides them in a `HashMap`. Every haplotype with the same bases lands in one
//! bucket whatever its uniqueness value, so the iteration order of a set of near-identical
//! haplotypes is decided by insertion within the bucket rather than by the hash. Nothing in this
//! slice depends on that order, and this note is here so that the first caller that does knows it
//! is not free.

use crate::allele_list::AlleleType;
use htsjdk_vcf::allele::{Allele, AlleleError};

/// `Haplotype`, over the state its allele identity is made of.
///
/// The derived `PartialEq` is **not** used: `Haplotype.equals` is not the field-by-field one, since
/// the `SimpleAllele` it inherits from compares bases and the reference flag while this class also
/// compares the uniqueness value. The impl below is written out for that reason.
#[derive(Debug, Clone)]
pub struct Haplotype {
    /// The `SimpleAllele` super-constructor's result: the bases as it stored them, and the flag.
    allele: Allele,
    /// `uniquenessValue`, "uniquely differentiates the haplotype from others with same ref/bases".
    /// Zero unless the assembler sets it, and part of equality whether or not it was set.
    uniqueness_value: i32,
}

impl Haplotype {
    /// `new Haplotype(bases, isRef)`.
    ///
    /// `Arrays.copyOf` then `super(...)`, which is `SimpleAllele`'s constructor and therefore the
    /// same validation, uppercasing and no-call handling as `Allele.create(bases, isRef)`. The
    /// refusals are that constructor's too, so a null allele or a symbolic reference is an error
    /// here exactly where it is one there.
    pub fn new(bases: &[u8], is_ref: bool) -> Result<Self, AlleleError> {
        Ok(Self {
            allele: Allele::create(bases, is_ref)?,
            uniqueness_value: 0,
        })
    }

    /// `new Haplotype(bases)`, which is the two-argument one with `isRef` false.
    pub fn non_reference(bases: &[u8]) -> Result<Self, AlleleError> {
        Self::new(bases, false)
    }

    /// `getUniquenessValue()`.
    pub fn uniqueness_value(&self) -> i32 {
        self.uniqueness_value
    }

    /// `setUniquenessValue(int)`.
    pub fn set_uniqueness_value(&mut self, value: i32) {
        self.uniqueness_value = value;
    }

    /// `getBases()`, which is `SimpleAllele`'s.
    ///
    /// The allele keeps its bases private and renders them through `getDisplayString`, which is the
    /// same bytes: the substitution `getBaseString` makes for a no-call is that method's, not this
    /// one's.
    pub fn bases(&self) -> Vec<u8> {
        self.allele.display_string().into_bytes()
    }

    /// `isReference()`.
    pub fn is_reference(&self) -> bool {
        self.allele.is_reference()
    }

    /// The haplotype as the allele it is, for a caller that wants the super-class view.
    pub fn as_allele(&self) -> &Allele {
        &self.allele
    }

    /// `hashCode()`: `Arrays.hashCode(getBases())`, and **not** the fields `equals` compares.
    ///
    /// `Arrays.hashCode` is `31 * result + element` from a seed of 1, over the signed byte values,
    /// wrapping like Java's `int`.
    pub fn java_hash_code(&self) -> i32 {
        let mut result: i32 = 1;
        for base in self.bases() {
            result = result.wrapping_mul(31).wrapping_add(i32::from(base as i8));
        }
        result
    }
}

/// `Haplotype.equals`: the uniqueness value, the reference flag, and the bases. Not the allele's
/// own equality, which stops at the last two.
///
/// The last two are asked of the allele in one comparison: `Allele`'s own equality is the bases,
/// the reference flag, and the two flags derived from the bases, so it answers the same question
/// the two clauses do.
impl PartialEq for Haplotype {
    fn eq(&self, other: &Self) -> bool {
        self.uniqueness_value == other.uniqueness_value && self.allele == other.allele
    }
}

impl Eq for Haplotype {}

impl AlleleType for Haplotype {
    fn is_reference(&self) -> bool {
        Haplotype::is_reference(self)
    }
}

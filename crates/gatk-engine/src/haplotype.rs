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
//! `alignmentStartHapwrtRef`, `score`, `kmerSize`, the genome location and the `EventMap`) is in
//! [`crate::haplotype_alignment`], set by the assembler and measured with it; equality never reads
//! it. `isCollapsed`, the flow-based haplotype collapsing flag, is not ported, so no event is ever
//! tagged as collapsed.
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
use crate::haplotype_alignment::HaplotypeAlignment;
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
    /// What the assembler gives a haplotype once it has one: where it is, its CIGAR against the
    /// reference, its path score and the k-mer size of its graph. None of it takes part in
    /// equality; the accessors are in [`crate::haplotype_alignment`].
    pub(crate) alignment: HaplotypeAlignment,
    /// `PartiallyDeterminedHaplotype`'s own state, for a haplotype that is one.
    pd: Option<Box<PartiallyDetermined>>,
}

/// What a `PartiallyDeterminedHaplotype` adds to its haplotype.
///
/// Its bases are the reference with every event of its branch applied except where the event is
/// undetermined: there the bases are the reference's (or the insertion's) and the alternate byte
/// beside each says what else the base may be: `SNP` with the bases it may also be, or the start
/// and end of a deletion that may or may not be there.
#[derive(Debug, Clone, PartialEq)]
pub struct PartiallyDetermined {
    /// `alternateBases`, one flag byte per base.
    pub alternate_bases: Vec<u8>,
    /// `constituentBuiltEvents`, in the order the haplotype was built from them.
    pub constituent_events: Vec<crate::event_map::Event>,
    /// `determinedEvents`: none when the reference is the determined allele, else one.
    pub determined_events: Vec<crate::event_map::Event>,
    /// `determinedPosition`.
    pub determined_position: i32,
    /// `getMaximumExtentOfSiteDeterminedAlleles()`: from the first start to the last end of every
    /// event at the determined locus.
    pub determined_extent: crate::interval::SimpleInterval,
}

/// `PartiallyDeterminedHaplotype.SNP`, `DEL_START`, `DEL_END` and the four bases.
pub mod pd_bytes {
    pub const SNP: u8 = 1;
    pub const DEL_START: u8 = 2;
    pub const DEL_END: u8 = 4;
    pub const A: u8 = 8;
    pub const C: u8 = 16;
    pub const G: u8 = 32;
    pub const T: u8 = 64;
    pub const N: u8 = 128;
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
            // The field initialisers: a null location and CIGAR, an alignment start of 0, a score
            // of `Double.NaN` (which is what a reference haplotype keeps, since the sequence-graph
            // finder never marks a path as the reference) and a k-mer size of 0.
            alignment: HaplotypeAlignment::default(),
            pd: None,
        })
    }

    /// `isPartiallyDetermined()`.
    pub fn is_partially_determined(&self) -> bool {
        self.pd.is_some()
    }

    /// The `PartiallyDeterminedHaplotype` state, for one.
    pub fn partially_determined(&self) -> Option<&PartiallyDetermined> {
        self.pd.as_deref()
    }

    /// Makes this haplotype a `PartiallyDeterminedHaplotype`.
    pub fn set_partially_determined(&mut self, pd: PartiallyDetermined) {
        self.pd = Some(Box::new(pd));
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
        let hash = |bytes: &[u8]| {
            let mut result: i32 = 1;
            for byte in bytes {
                result = result.wrapping_mul(31).wrapping_add(i32::from(*byte as i8));
            }
            result
        };
        let bases = hash(&self.bases());
        match &self.pd {
            None => bases,
            // `Objects.hash(Arrays.hashCode(getBases()), Arrays.hashCode(alternateBases))`.
            Some(pd) => 31i32
                .wrapping_add(bases)
                .wrapping_mul(31)
                .wrapping_add(hash(&pd.alternate_bases)),
        }
    }
}

/// `Haplotype.equals`: the uniqueness value, the reference flag, and the bases, and for a partially
/// determined one its own fields as well. Not the allele's
/// own equality, which stops at the last two.
///
/// The last two are asked of the allele in one comparison: `Allele`'s own equality is the bases,
/// the reference flag, and the two flags derived from the bases, so it answers the same question
/// the two clauses do.
impl PartialEq for Haplotype {
    fn eq(&self, other: &Self) -> bool {
        let same = self.uniqueness_value == other.uniqueness_value && self.allele == other.allele;
        match (&self.pd, &other.pd) {
            // `PartiallyDeterminedHaplotype.equals`: only another one, at the same determined
            // position and with the same alternate bases.
            (Some(this), Some(that)) => {
                same && this.determined_position == that.determined_position
                    && this.alternate_bases == that.alternate_bases
            }
            (Some(_), None) => false,
            // `Haplotype.equals`, which does not ask whether the other is partially determined.
            (None, _) => same,
        }
    }
}

impl Eq for Haplotype {}

impl AlleleType for Haplotype {
    fn is_reference(&self) -> bool {
        Haplotype::is_reference(self)
    }
}

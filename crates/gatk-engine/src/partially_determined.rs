//! `PartiallyDeterminedHaplotype`'s own state, ported from
//! `org.broadinstitute.hellbender.utils.haplotype.PartiallyDeterminedHaplotype` (GATK 4.7.0.0;
//! `PartiallyDeterminedHaplotype.java` is unchanged since 4.6.2.0).
//!
//! The reference makes it a subclass of `Haplotype`; here a [`crate::haplotype::Haplotype`] carries
//! it as an optional field, so that every list of haplotypes stays one type and only the code that
//! cares asks whether a haplotype is partially determined. What the subclass overrides is answered
//! where it matters: `equals` and `hashCode` in [`crate::haplotype`], the event map in
//! [`crate::event_map::build_event_maps_for_haplotypes`], and the allele mapper's skip in
//! [`crate::allele_mapping`].
//!
//! # Reading this file without reading Rust
//!
//! - `Option<&PartiallyDetermined>` is "a borrowed `PartiallyDetermined`, or nothing": the
//!   compiler-checked form of a reference that may be `null`, and here the answer to Java's
//!   `instanceof PartiallyDeterminedHaplotype`;
//! - `let Some(pd) = pd else { return ...; };` unwraps the value or takes the early exit;
//! - `match (this, that) { ... }` tests both values at once, and the compiler refuses it unless
//!   every combination is covered, so the three arms below are the whole truth table;
//! - `u8` is an unsigned byte (0 to 255) where Java's `byte` is signed (-128 to 127), and
//!   `wrapping_mul` / `wrapping_add` are `int` arithmetic that wraps on overflow as Java's does,
//!   where plain `*` and `+` would stop the program in a debug build.

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
    /// `determinedPosition`, the reference-contig position this haplotype is determined at.
    /// An `int` field in the reference; its getter widens it to `long`.
    pub determined_position: i32,
    /// `getMaximumExtentOfSiteDeterminedAlleles()`: from the first start to the last end of every
    /// event at the determined locus.
    pub determined_extent: crate::interval::SimpleInterval,
}

/// `PartiallyDeterminedHaplotype.SNP`, `DEL_START`, `DEL_END` and the four bases (and `N`): the
/// bits of one `alternateBases` byte.
///
/// ```java
/// final public static byte N = (byte) 128;
/// ```
///
/// In Java `N` is therefore the signed byte -128. Here every flag is an unsigned `u8`, so `N` is
/// 128 with the same bit pattern; the one place the sign is observable, the hash, converts back
/// (see `hash_code`).
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

impl crate::haplotype::Haplotype {
    /// `isPartiallyDetermined()`.
    pub fn is_partially_determined(&self) -> bool {
        self.pd.is_some()
    }

    /// The `PartiallyDeterminedHaplotype` state, or `None` for a plain haplotype. Where Java
    /// would test `instanceof` and cast, a caller here asks for this and matches on the answer.
    pub fn partially_determined(&self) -> Option<&PartiallyDetermined> {
        self.pd.as_deref()
    }

    /// Makes this haplotype a `PartiallyDeterminedHaplotype`: the part of the reference
    /// constructor that sets the subclass's fields. The state is boxed (held behind a pointer, as
    /// every Java object is) so that a plain haplotype pays one pointer for it, not its full size.
    pub fn set_partially_determined(&mut self, pd: PartiallyDetermined) {
        self.pd = Some(Box::new(pd));
    }
}

/// `PartiallyDeterminedHaplotype.hashCode()`: `Objects.hash(Arrays.hashCode(getBases()),
/// Arrays.hashCode(alternateBases))`, which is `31 * (31 + bases) + alternate`; any other haplotype
/// keeps the hash of its bases.
///
/// ```java
/// return Objects.hash(Arrays.hashCode(getBases()),Arrays.hashCode(alternateBases));
/// ```
///
/// `bases` is `Arrays.hashCode(getBases())`, already computed by the caller.
pub(crate) fn hash_code(bases: i32, pd: Option<&PartiallyDetermined>) -> i32 {
    let Some(pd) = pd else {
        return bases;
    };
    // `Arrays.hashCode(byte[])`: start at 1, then `31 * h + b` for each byte. Java's `b` is the
    // SIGNED byte, so `as i8` reinterprets the stored `u8` (128 becomes -128, the reference's
    // `N`) before `i32::from` widens it with its sign.
    let mut alternate: i32 = 1;
    for byte in &pd.alternate_bases {
        alternate = alternate
            .wrapping_mul(31)
            .wrapping_add(i32::from(*byte as i8));
    }
    31i32
        .wrapping_add(bases)
        .wrapping_mul(31)
        .wrapping_add(alternate)
}

/// `PartiallyDeterminedHaplotype.equals`: only another partially determined haplotype, at the same
/// determined position and with the same alternate bases. The plain `Haplotype.equals` does not
/// ask whether the other is one, so the relation is not symmetric, as the reference's is not.
///
/// ```java
/// return h instanceof PartiallyDeterminedHaplotype
///         && getUniquenessValue() == ((Haplotype) h).getUniquenessValue()
///         && isReference() == ((Haplotype) h).isReference()
///         && determinedPosition == ((PartiallyDeterminedHaplotype) h).determinedPosition
///         && Arrays.equals(getBases(), ((Haplotype) h).getBases())
///         && Arrays.equals(alternateBases, ((PartiallyDeterminedHaplotype) h).alternateBases);
/// ```
///
/// (The reference's trailing comment on the `determinedPosition` line is left out of the quote.)
///
/// `same` is the three tests the two classes share (uniqueness value, reference flag, bases),
/// which `Haplotype::eq` computes; `this` and `that` are the two sides' subclass state, `None`
/// for a plain haplotype. The arms read: both partially determined, the reference's test above;
/// `this` partially determined and `that` not, the `instanceof` fails; `this` plain, the plain
/// `Haplotype.equals`, which never looks at `that`'s subclass.
pub(crate) fn equals(
    same: bool,
    this: Option<&PartiallyDetermined>,
    that: Option<&PartiallyDetermined>,
) -> bool {
    match (this, that) {
        (Some(this), Some(that)) => {
            same && this.determined_position == that.determined_position
                && this.alternate_bases == that.alternate_bases
        }
        (Some(_), None) => false,
        (None, _) => same,
    }
}

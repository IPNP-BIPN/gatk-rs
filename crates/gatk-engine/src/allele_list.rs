//! Ported from `org.broadinstitute.hellbender.utils.genotyper.AlleleList`, `IndexedAlleleList`,
//! `IndexedSampleList` and `org.broadinstitute.hellbender.utils.collections.IndexedSet`
//! (GATK 4.7.0.0; all four Java files are unchanged since 4.6.2.0).
//!
//! These are the two axes of every likelihood matrix in GATK: a list of samples and a list of
//! alleles, each with a stable index. `AlleleLikelihoods` is defined over them, so their
//! surprises are inherited by everything that reads a likelihood.
//!
//! # An indexed list is a **set**, and it swallows duplicates
//!
//! ```java
//! for (final E value : values) {
//!     if (indexByElement.containsKey(value)) {
//!         continue;
//!     }
//!     indexByElement.put(value, nextIndex++);
//!     elements.add(value);
//! }
//! ```
//!
//! A list built from `[A, C, A]` has two entries, not three, and the surviving `A` keeps the index
//! of its **first** occurrence. Nothing reports the drop. A caller that built its allele list from
//! a variant context with a repeated allele therefore gets a matrix with fewer rows than it asked
//! for, and the indices it was about to use are off by one from that point on.
//!
//! # Membership is `Allele.equals`, which includes the reference flag
//!
//! `indexOfAllele` is a hash lookup on the allele, and two alleles with the same bases and
//! different reference flags are different keys. A list can hold both, and `indexOfReference`
//! returns the first entry whose flag is set, scanning linearly.
//!
//! # A permutation is a **subset** map, and it is directional
//!
//! `permutation(target)` refuses a target longer than the original, and refuses a target holding an
//! allele the original does not: it can drop alleles and reorder them, never invent them. It
//! reports three things that a caller uses to skip work: `is_partial` (the target is shorter),
//! `is_non_permuted` (every kept allele is at the same index it started at) and `is_kept` (per
//! original index).
//!
//! `is_non_permuted` is computed only from the indices, so a target that is the same list in the
//! same order is non-permuted even when it was constructed separately; and `permutation` short
//! circuits on equality before building anything, which is the only path that produces a
//! permutation over a target longer than zero without scanning.
//!
//! # The allele axis is a type parameter, because in the reference it is one
//!
//! ```java
//! public interface AlleleList<A extends Allele> {
//! ```
//!
//! Most of GATK instantiates it at `Allele`, which is why this list was written that way first.
//! `HaplotypeCaller` instantiates it at [`crate::haplotype::Haplotype`], and
//! `HaplotypeFilteringAnnotation` reads that instantiation. The parameter defaults to `Allele`
//! here so that every existing caller keeps its spelling.
//!
//! The Java bound is inheritance: `A` **is** an `Allele`, so every `Allele` member is in scope. The
//! only one this file uses is `isReference()`, so the bound here is [`AlleleType`], a trait with
//! that one method. Anything wider would be claiming a relationship the code does not exercise.
//!
//! # Reading this file without reading Rust
//!
//! - `AlleleList<A: AlleleType = Allele>` is a generic class with a **default** type argument:
//!   plain `AlleleList` means `AlleleList<Allele>`. `A: AlleleType` is the bound, Java's
//!   `A extends ...`, except that a type meets it by implementing the trait (an interface), not
//!   by inheriting from a class.
//! - `#[derive(PartialEq)]` writes `equals` field by field, so two lists are `==` when their
//!   element vectors are: same length, and each pair equal under `A`'s own `==`.
//! - `Option<usize>` is "an index or nothing"; `None` is the reference's `-1`, and `usize` is an
//!   unsigned index type, so a negative index cannot be passed in at all.
//! - `|existing| existing == value` is a closure, Java's lambda. `.iter().position(..)` is the
//!   index of the first element it accepts, `.any(..)` whether any element does.
//! - `x.ok_or(e)?` turns "nothing" into the error `e` and returns it to the caller at once; that
//!   is how this file throws. `.expect(..)` unwraps a value that cannot be missing and aborts the
//!   program if it is.
//! - `.get(i)` on a vector is a bounds-checked read that answers `None` past the end instead of
//!   throwing; `.copied()` turns the borrowed element into a plain value.

use htsjdk_vcf::allele::Allele;

/// The `A extends Allele` bound, narrowed to what an allele list actually asks of its elements.
///
/// `Haplotype` satisfies it by being a `SimpleAllele` in the reference; here it satisfies it by
/// implementing this trait.
pub trait AlleleType: Clone + PartialEq {
    /// `Allele.isReference()`.
    fn is_reference(&self) -> bool;
}

impl AlleleType for Allele {
    fn is_reference(&self) -> bool {
        Allele::is_reference(self)
    }
}

/// `IndexedSet<E>`, which this file uses over alleles and over sample names.
///
/// Insertion-ordered, first occurrence wins, later duplicates dropped without a word.
///
/// **Why a vector scanned with `==` and not a hash map**: the reference keeps a
/// `LinkedHashMap<E, Integer>` beside the list, so membership is `hashCode` and then `equals`.
/// A scan with `==` finds the same entry whenever equal elements hash equally, which holds for
/// the types used here: `SimpleAllele.hashCode` is `Arrays.hashCode(bases) * 31 +
/// Boolean.hashCode(isRef)` over the same fields its `equals` compares, and `Haplotype.hashCode`
/// is `Arrays.hashCode(getBases())`, bases being part of its `equals`. The bound needs only
/// `PartialEq`, which is what `Haplotype` implements.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IndexedSet<T: Clone + PartialEq> {
    elements: Vec<T>,
}

impl<T: Clone + PartialEq> Default for IndexedSet<T> {
    fn default() -> Self {
        Self {
            elements: Vec::new(),
        }
    }
}

impl<T: Clone + PartialEq> IndexedSet<T> {
    /// `IndexedSet(Collection<E>)` and `IndexedSet(E...)`. The reference also refuses a `null`
    /// element, `Utils.nonNull(value, "null element not allowed: index == " + nextIndex);`; a
    /// Rust slice cannot hold one, so there is nothing to check.
    pub fn new(values: &[T]) -> Self {
        // `indexByElement.containsKey(value)`: skip a value already kept, so the survivor keeps
        // the index of its first occurrence.
        let mut elements: Vec<T> = Vec::with_capacity(values.len());
        for value in values {
            if elements.iter().any(|existing| existing == value) {
                continue;
            }
            elements.push(value.clone());
        }
        Self { elements }
    }

    pub fn size(&self) -> usize {
        self.elements.len()
    }

    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }

    /// `indexOf`: the position, or `-1`. The `-1` is kept as `None` here and turned back into a
    /// `-1` only where the reference's own callers compare against it.
    pub fn index_of(&self, value: &T) -> Option<usize> {
        self.elements.iter().position(|existing| existing == value)
    }

    /// `get(index)`, which checks `Utils.validIndex(index, size());` and so throws
    /// `IllegalArgumentException` past the end; here that is `None`.
    pub fn get(&self, index: usize) -> Option<&T> {
        self.elements.get(index)
    }

    pub fn as_slice(&self) -> &[T] {
        &self.elements
    }
}

/// `IndexedAlleleList<A>`: an `IndexedSet<A>` and nothing else. Its constructors document the
/// duplicate rule the set implements: "Repeats in the input collection will be ignored (keeping
/// the first one)."
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlleleList<A: AlleleType = Allele> {
    alleles: IndexedSet<A>,
}

impl<A: AlleleType> Default for AlleleList<A> {
    fn default() -> Self {
        Self::empty()
    }
}

impl<A: AlleleType> AlleleList<A> {
    /// `new IndexedAlleleList<>(alleles)`, which is also what `AlleleList.newList` returns.
    pub fn new(alleles: &[A]) -> Self {
        Self {
            alleles: IndexedSet::new(alleles),
        }
    }

    /// `AlleleList.emptyAlleleList()`.
    pub fn empty() -> Self {
        Self::new(&[])
    }

    pub fn number_of_alleles(&self) -> usize {
        self.alleles.size()
    }

    pub fn is_empty(&self) -> bool {
        self.alleles.is_empty()
    }

    /// `indexOfAllele`, which is a lookup by `Allele.equals`: bases **and** the reference flag.
    ///
    /// At `A = Haplotype` the equality is that class's own, which adds the uniqueness value, so a
    /// haplotype list can hold two entries with identical bases and reference flag where an allele
    /// list cannot.
    ///
    /// The reference takes any allele, `int indexOfAllele(final Allele allele);`; here the
    /// argument has the list's own element type, so a haplotype list is queried with haplotypes.
    pub fn index_of_allele(&self, allele: &A) -> Option<usize> {
        self.alleles.index_of(allele)
    }

    pub fn contains_allele(&self, allele: &A) -> bool {
        self.index_of_allele(allele).is_some()
    }

    /// `getAllele(index)`. The reference throws `IllegalArgumentException` past the end in both
    /// cases, with different messages: the empty list's `"allele index is out of range"`, and
    /// otherwise `Utils.validIndex(index, alleles.size());`, which throws the same type. Both are
    /// `None` here.
    pub fn get_allele(&self, index: usize) -> Option<&A> {
        self.alleles.get(index)
    }

    pub fn as_slice(&self) -> &[A] {
        self.alleles.as_slice()
    }

    /// `indexOfReference()`: the **first** allele whose reference flag is set, found by scanning.
    /// A list can hold more than one and this reports only the first.
    pub fn index_of_reference(&self) -> Option<usize> {
        self.alleles
            .as_slice()
            .iter()
            .position(|allele| allele.is_reference())
    }

    /// `AlleleList.equals(first, second)`: same length, same alleles, same order. Slice `==`
    /// compares the lengths first and then each pair with `A`'s own `==`, which is the reference's
    /// loop over `if (!firstSample.equals(secondSample)) {`.
    pub fn same_alleles(&self, other: &AlleleList<A>) -> bool {
        self.as_slice() == other.as_slice()
    }

    /// `permutation(target)`:
    ///
    /// ```java
    /// if (equals(this, target)) {
    ///     return new NonPermutation<>(this);
    /// } else {
    ///     return new ActualPermutation<>(this, target);
    /// }
    /// ```
    pub fn permutation(
        &self,
        target: &AlleleList<A>,
    ) -> Result<AllelePermutation<A>, PermutationError> {
        AllelePermutation::new(self, target)
    }
}

/// `IndexedSampleList`. The same set semantics, over sample names.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SampleList {
    samples: IndexedSet<String>,
}

impl SampleList {
    pub fn new(samples: &[String]) -> Self {
        Self {
            samples: IndexedSet::new(samples),
        }
    }

    pub fn number_of_samples(&self) -> usize {
        self.samples.size()
    }

    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }

    /// `indexOfSample`, the set's `samples.indexOf(sample)`. The scan is written here rather than
    /// through [`IndexedSet::index_of`] because the argument is a borrowed `&str` and the set
    /// holds `String`s, which compare with `==` directly.
    pub fn index_of_sample(&self, sample: &str) -> Option<usize> {
        self.samples
            .as_slice()
            .iter()
            .position(|existing| existing == sample)
    }

    /// `getSample`, through `IndexedSet.get`: `None` where the reference throws
    /// `IllegalArgumentException`.
    pub fn get_sample(&self, index: usize) -> Option<&String> {
        self.samples.get(index)
    }

    pub fn as_slice(&self) -> &[String] {
        self.samples.as_slice()
    }
}

/// What `ActualPermutation`'s constructor refuses. Both are the same
/// `IllegalArgumentException` with the same message, reached two different ways (`fromSize <
/// toSize`, then `originalIndex < 0`), both thrown with
/// `"target allele list is not a permutation of the original allele list"`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PermutationError {
    /// The target is longer than the original.
    TargetLonger,
    /// The target holds an allele the original does not.
    AlleleNotInOriginal,
}

impl PermutationError {
    pub fn class(&self) -> &'static str {
        "java.lang.IllegalArgumentException"
    }

    pub fn message(&self) -> &'static str {
        "target allele list is not a permutation of the original allele list"
    }
}

/// `AlleleListPermutation`, covering both `NonPermutation` and `ActualPermutation`.
///
/// The reference has two classes: `NonPermutation` holds only the list and answers every index
/// with itself, `ActualPermutation` holds the two lists and the index tables. Here one struct
/// holds the tables in both cases, filled with the identity for the first.
///
/// Two things differ from the reference outside the in-range answers. The reference's
/// permutation is itself an `AlleleList` (its `getAllele` reads the target); this one is not.
/// And an out-of-range index is `None` or `false` here, where `ActualPermutation` indexes its
/// arrays and `NonPermutation` returns the argument unchecked (`isKept` is `return true;`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AllelePermutation<A: AlleleType = Allele> {
    from: AlleleList<A>,
    to: AlleleList<A>,
    /// `fromIndex[toIndex]`: where each target allele came from.
    from_index: Vec<usize>,
    /// `keptFromIndices[fromIndex]`.
    kept: Vec<bool>,
    is_partial: bool,
    is_non_permuted: bool,
}

impl<A: AlleleType> AllelePermutation<A> {
    fn new(original: &AlleleList<A>, target: &AlleleList<A>) -> Result<Self, PermutationError> {
        // `permutation()` short circuits on equality to `NonPermutation`. For equal lists the
        // loop below would build the same tables (identity, all kept, not partial, non-permuted),
        // so the short circuit changes no in-range answer; it only skips the lookups.
        if original.same_alleles(target) {
            return Ok(Self {
                from: original.clone(),
                to: target.clone(),
                from_index: (0..original.number_of_alleles()).collect(),
                kept: vec![true; original.number_of_alleles()],
                is_partial: false,
                is_non_permuted: true,
            });
        }

        let from_size = original.number_of_alleles();
        let to_size = target.number_of_alleles();
        // `if (fromSize < toSize) {`: the first refusal, before any lookup.
        if from_size < to_size {
            return Err(PermutationError::TargetLonger);
        }

        // `vec![false; n]` is `new boolean[n]`.
        let mut kept = vec![false; from_size];
        let mut from_index = Vec::with_capacity(to_size);
        // `nonPermuted` starts as "the sizes agree" and is then ANDed with "every target allele is
        // where it was", so a partial permutation is never non-permuted.
        let mut non_permuted = from_size == to_size;
        for i in 0..to_size {
            let allele = target.get_allele(i).expect("index below the target's size");
            // `if (originalIndex < 0) {`, the second refusal: a target allele the original lacks.
            let original_index = original
                .index_of_allele(allele)
                .ok_or(PermutationError::AlleleNotInOriginal)?;
            kept[original_index] = true;
            from_index.push(original_index);
            non_permuted &= original_index == i;
        }

        Ok(Self {
            from: original.clone(),
            to: target.clone(),
            from_index,
            kept,
            is_partial: from_size != to_size,
            is_non_permuted: non_permuted,
        })
    }

    pub fn is_partial(&self) -> bool {
        self.is_partial
    }

    pub fn is_non_permuted(&self) -> bool {
        self.is_non_permuted
    }

    /// `toIndex(fromIndex)`: where an original allele ended up, or `None` if it was dropped. The
    /// reference computes it by looking the allele up in the target rather than by inverting the
    /// table, so a dropped allele answers `-1` rather than being an error.
    pub fn to_index(&self, from_index: usize) -> Option<usize> {
        let allele = self.from.get_allele(from_index)?;
        self.to.index_of_allele(allele)
    }

    /// `fromIndex(toIndex)`: where a target allele came from, `return fromIndex[toIndex];`.
    pub fn from_index(&self, to_index: usize) -> Option<usize> {
        self.from_index.get(to_index).copied()
    }

    /// `isKept(fromIndex)`, `return keptFromIndices[fromIndex];`: whether the original allele at
    /// that index appears anywhere in the target.
    pub fn is_kept(&self, from_index: usize) -> bool {
        self.kept.get(from_index).copied().unwrap_or(false)
    }

    pub fn from_size(&self) -> usize {
        self.from.number_of_alleles()
    }

    pub fn to_size(&self) -> usize {
        self.to.number_of_alleles()
    }

    pub fn from_list(&self) -> &[A] {
        self.from.as_slice()
    }

    /// `toList()`. On the equal-lists path the reference's `NonPermutation` returns the original's
    /// elements; this returns the target's. They are equal under `A`'s `==`, but a `Haplotype`
    /// carries state outside its equality (its alignment), and that state is the target's here.
    pub fn to_list(&self) -> &[A] {
        self.to.as_slice()
    }
}

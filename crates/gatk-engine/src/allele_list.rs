//! Ported from `org.broadinstitute.hellbender.utils.genotyper.AlleleList`, `IndexedAlleleList`,
//! `IndexedSampleList` and `org.broadinstitute.hellbender.utils.collections.IndexedSet`
//! (GATK 4.6.2.0).
//!
//! These are the two axes of every likelihood matrix in GATK: a list of samples and a list of
//! alleles, each with a stable index. `AlleleLikelihoods` is defined over them, so their
//! surprises are inherited by everything that reads a likelihood.
//!
//! # An indexed list is a **set**, and it swallows duplicates
//!
//! ```java
//! for (final E value : values) {
//!     if (indexByElement.containsKey(value)) { continue; }
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
//! public interface AlleleList<A extends Allele> { ... }
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

/// `IndexedSet`, specialised to the two element types GATK uses it with here.
///
/// Insertion-ordered, first occurrence wins, later duplicates dropped without a word.
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
    pub fn new(values: &[T]) -> Self {
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

    pub fn get(&self, index: usize) -> Option<&T> {
        self.elements.get(index)
    }

    pub fn as_slice(&self) -> &[T] {
        &self.elements
    }
}

/// `IndexedAlleleList<A>`.
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
    /// haplotype list can hold two entries with identical bases where an allele list cannot.
    pub fn index_of_allele(&self, allele: &A) -> Option<usize> {
        self.alleles.index_of(allele)
    }

    pub fn contains_allele(&self, allele: &A) -> bool {
        self.index_of_allele(allele).is_some()
    }

    /// `getAllele(index)`. The reference throws `IllegalArgumentException` past the end for the
    /// empty list and `IndexOutOfBoundsException` for a non-empty one, which is a distinction no
    /// caller depends on; both are `None` here.
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

    /// `AlleleList.equals(first, second)`: same length, same alleles, same order.
    pub fn same_alleles(&self, other: &AlleleList<A>) -> bool {
        self.as_slice() == other.as_slice()
    }

    /// `permutation(target)`.
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

    pub fn index_of_sample(&self, sample: &str) -> Option<usize> {
        self.samples
            .as_slice()
            .iter()
            .position(|existing| existing == sample)
    }

    pub fn get_sample(&self, index: usize) -> Option<&String> {
        self.samples.get(index)
    }

    pub fn as_slice(&self) -> &[String] {
        self.samples.as_slice()
    }
}

/// What `ActualPermutation`'s constructor refuses. Both are the same
/// `IllegalArgumentException` with the same message, reached two different ways.
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
/// The reference has two classes because one of them can answer without any state; here the same
/// state answers both, and the flags say which case was taken.
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
        // `permutation()` short circuits on equality, which is the only way `isPartial` comes back
        // false for a target that was built separately.
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
        if from_size < to_size {
            return Err(PermutationError::TargetLonger);
        }

        let mut kept = vec![false; from_size];
        let mut from_index = Vec::with_capacity(to_size);
        // `nonPermuted` starts as "the sizes agree" and is then ANDed with "every target allele is
        // where it was", so a partial permutation is never non-permuted.
        let mut non_permuted = from_size == to_size;
        for i in 0..to_size {
            let allele = target.get_allele(i).expect("index below the target's size");
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

    /// `fromIndex(toIndex)`: where a target allele came from.
    pub fn from_index(&self, to_index: usize) -> Option<usize> {
        self.from_index.get(to_index).copied()
    }

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

    pub fn to_list(&self) -> &[A] {
        self.to.as_slice()
    }
}

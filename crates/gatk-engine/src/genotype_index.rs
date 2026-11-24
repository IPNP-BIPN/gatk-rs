//! `GenotypeIndexCalculator` and `GenotypeAlleleCounts`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.genotyper` (GATK 4.7.0.0; both `.java` files are
//! unchanged since 4.6.2.0).
//!
//! The combinatorics under every PL array: which genotype each index of a likelihood vector means,
//! and where a genotype lands when the allele list changes.
//!
//! # The order is not the obvious one
//!
//! A diploid site with three alleles has six genotypes, and the order a PL array is indexed by is
//!
//! ```text
//! 0/0  0/1  1/1  0/2  1/2  2/2
//! ```
//!
//! not `0/0 0/1 0/2 1/1 1/2 2/2`. The highest allele of the genotype decides first, which is why
//! adding an allele appends to a PL array rather than interleaving into it, and why a port that
//! iterates the obvious way is wrong from the third entry (index 2) on. The class javadoc states
//! the rule for two genotypes `g` and `h` of N alleles each, sorted ascending:
//!
//! ```java
//!  *    (i)   the order of haploid genotypes is simply the allele ordering
//!  *    (ii)  if g_N < h_N then g < h
//!  *    (iii) if g_N = h_N then the order is that of the first N-1 alleles
//! ```
//!
//! # The index is a sum of binomials
//!
//! `indexOfFirstGenotypeWithAllele` answers `0` for allele 0 (`allele == 0 ? 0 : ...`) and
//! otherwise
//!
//! ```java
//! CombinatoricsUtils.binomialCoefficient(ploidy + allele - 1, allele - 1);
//! ```
//!
//! That is the index of the first genotype containing a given allele at a given ploidy, and the
//! index of a genotype is that function summed over its alleles, sorted, from the highest down.
//! The same function at the allele count is the number of genotypes.
//!
//! Worked example, writing `f(p, a)` for that function, `C(p + a - 1, a - 1)`. `calculateIndex`
//! sorts the genotype and scores its k-th allele from the top (k from 0) at ploidy `p - k`:
//!
//! ```text
//! diploid 1/2:     f(2, 2) + f(1, 1)           = C(3, 1) + C(1, 0)       = 3 + 1     = 4
//! diploid 2/2:     f(2, 2) + f(1, 2)           = C(3, 1) + C(2, 1)       = 3 + 2     = 5
//! triploid 0/1/2:  f(3, 2) + f(2, 1) + f(1, 0) = C(4, 1) + C(2, 0) + 0   = 4 + 1 + 0 = 5
//! genotype count, diploid, three alleles:   f(2, 3) = C(4, 2) = 6
//! ```
//!
//! The diploid rows match the order above; the triploid row is the Java's own "genotype is ABC"
//! example and matches index 5, `0/1/2`, in the `GenotypeAlleleCounts` javadoc table.
//!
//! # Overflow, as the Java has it
//!
//! `indexOfFirstGenotypeWithAllele` returns a `long`, from commons-math3's
//! `CombinatoricsUtils.binomialCoefficient` (no local source, so what that does past a `long` is
//! not restated here). Its two callers then narrow it to `int` in two different ways.
//!
//! `genotypeCount` refuses, and its bound is strict (`MathUtils.LONG_OVERFLOW` is `-1`):
//!
//! ```java
//! Utils.validateArg(result != MathUtils.LONG_OVERFLOW && result < Integer.MAX_VALUE, () ->
//! ```
//!
//! `calculateIndex` never refuses: each term is narrowed by `(int)` and `IndexRange.sumInt` adds
//! them in an `int` (`result += lambda.applyAsInt(i);`), which wraps silently.
//!
//! The port has one rule for both: [`index_of_first_genotype_with_allele`] answers `None` for any
//! binomial above `i32::MAX`. It agrees with the Java wherever every term fits in an `int`, which
//! holds whenever the genotype's alleles are all below an allele count that `genotype_count`
//! accepts (each term is then below that count, and so is their sum, the index). Past that the two
//! differ: `genotype_count` accepts a count of exactly `i32::MAX` (ploidy 1 with 2147483647 alleles
//! is one such input), which the strict `<` refuses, and an index term past `i32::MAX` counts as 0
//! here where the Java narrows it with `(int)`.
//!
//! # Subsetting is a permutation
//!
//! `subsetted_pl_indices` says where each new index takes its likelihood from. Keeping three
//! alleles in the order 0, 2, 1 gives `0,3,5,1,4,2`: nothing is dropped and everything moves.
//!
//! # Reading this file without reading Rust
//!
//! - `Option<usize>` is "an index or nothing": `Some(n)` or `None`. A `?` after an expression that
//!   gives an `Option` returns `None` from the enclosing function at once when it is `None`;
//! - `Result<T, GenotypeIndexError>` is "a `T` or one of the refusals below"; it is how this
//!   codebase models a method that throws `IllegalArgumentException`. `.ok_or(e)` turns a `None`
//!   into the refusal `e`;
//! - `usize` is an unsigned index type Java lacks, so a negative count cannot even be passed, and
//!   `u128` is an unsigned 128-bit integer; `checked_mul` multiplies and gives `None` rather than
//!   wrapping;
//! - `&[usize]` is a borrowed, read-only view of an `int[]`; `&mut [usize]` a writable one, through
//!   which the function changes the caller's array in place;
//! - `|step| expr` is a closure, Java's lambda `step -> expr`, and `.map(..).sum()` over a range is
//!   Java's `IntStream.range(..).map(..).sum()`;
//! - a function's last expression is its return value; `return` is only written for early exits.

/// What the combinatorics refuses, as the reference's `IllegalArgumentException`s.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GenotypeIndexError {
    /// `Utils.validateArg((alleleCountArray.length & 1) == 0, ...)`. Of the reference's checks in
    /// `alleleCountsToIndex` it is the only one this port can fail: `Utils.nonNull` has nothing to
    /// catch (a Rust slice is never null) and `count >= 0` neither (a `usize` is never negative).
    OddLengthCounts,
    /// A binomial too large for the reference's `int`, which it refuses rather than wraps.
    TooManyGenotypes { ploidy: usize, alleles: usize },
}

impl GenotypeIndexError {
    /// The reference's message text. For `TooManyGenotypes` it stops where the Java's format
    /// string goes on to `: approx. %.0f` (a `CombinatoricsUtils.binomialCoefficientDouble`
    /// estimate, `GenotypeIndexCalculator.genotypeCount`), so that tail is not reproduced.
    pub fn message(&self) -> String {
        // `match` is Java's `switch` over the variants; the compiler rejects it if one is missing.
        match self {
            GenotypeIndexError::OddLengthCounts => {
                "the allele counts array cannot have odd length".to_string()
            }
            GenotypeIndexError::TooManyGenotypes { ploidy, alleles } => format!(
                "the number of genotypes is too large for ploidy {ploidy} and {alleles} alleles"
            ),
        }
    }

    /// Both refusals are `Utils.validateArg`, which throws `IllegalArgumentException`.
    pub fn java_class(&self) -> &'static str {
        "java.lang.IllegalArgumentException"
    }
}

/// `C(n, k)` exactly, or `None` once it passes `i32::MAX`: the port's stand-in for commons-math3's
/// `CombinatoricsUtils.binomialCoefficient` together with the `int` bound its callers apply (see
/// the module note on overflow).
///
/// **How**: the multiplicative formula over `min(k, n - k)` steps. After step `s` the accumulator
/// is `C(n, s + 1)`, so each division is exact. `C(n, j)` grows with `j` up to `n / 2`, so the
/// first intermediate past `i32::MAX` proves the answer is past it too, and the early exit loses
/// nothing.
///
/// **Why `u128`**: before its division a step's product is at most `i32::MAX * n`, which fits in
/// 128 bits for any `usize` `n`, so the `checked_mul` never actually fails.
///
/// `k > n` answers 0. [`index_of_first_genotype_with_allele`] never asks it: there
/// `k = allele - 1`, which is at most `n = ploidy + allele - 1`.
fn binomial(n: usize, k: usize) -> Option<usize> {
    if k > n {
        return Some(0);
    }
    // `C(n, k) = C(n, n - k)`. Taking the smaller keeps every step at or below the middle of the
    // row, where `C(n, j)` still grows, which is what makes the early exit below exact.
    let k = k.min(n - k);
    let mut result: u128 = 1;
    for step in 0..k {
        // `(n - step) as u128` widens without loss. The `?` would return `None` on a `u128`
        // overflow, which the bound below rules out.
        result = result.checked_mul((n - step) as u128)? / (step as u128 + 1);
        // The port's single `int` bound: strictly above, so `i32::MAX` itself passes, where
        // `genotypeCount`'s `result < Integer.MAX_VALUE` would refuse it.
        if result > i32::MAX as u128 {
            return None;
        }
    }
    // `try_from` checks the narrowing rather than truncating, and `.ok()` turns its failure into
    // `None`. Under the bound above it always fits.
    usize::try_from(result).ok()
}

/// `GenotypeIndexCalculator.indexOfFirstGenotypeWithAllele(ploidy, allele)`: how many genotypes of
/// this ploidy come before the first one containing `allele`, `C(ploidy + allele - 1, allele - 1)`.
/// The Java's javadoc derives it by counting strings of `ploidy` `x`s and `allele - 1` slashes.
///
/// `None` where the Java would return a `long` above `i32::MAX` (module note on overflow).
pub fn index_of_first_genotype_with_allele(ploidy: usize, allele: usize) -> Option<usize> {
    if allele == 0 {
        return Some(0);
    }
    binomial(ploidy + allele - 1, allele - 1)
}

/// `genotypeCount(ploidy, alleleCount)`, which is the same function at the allele count.
///
/// The first genotype that would hold allele `alleleCount` is the one just past the last genotype
/// over alleles `0..alleleCount`, so its index is their number. Refuses a count above `i32::MAX`;
/// the Java's strict `result < Integer.MAX_VALUE` refuses `i32::MAX` itself as well.
pub fn genotype_count(ploidy: usize, allele_count: usize) -> Result<usize, GenotypeIndexError> {
    index_of_first_genotype_with_allele(ploidy, allele_count).ok_or(
        GenotypeIndexError::TooManyGenotypes {
            ploidy,
            alleles: allele_count,
        },
    )
}

/// Every genotype of one shape, in the canonical order, each as its sorted allele indices.
///
/// The reference walks these with a mutable `GenotypeAlleleCounts` that its iterator hands out
/// **the same object of** each time. Here each genotype is its own value, which is the one liberty
/// taken: nothing observable depends on the sharing, and a caller that kept the reference's object
/// would find it had changed underneath.
///
/// **How**: the reference steps from one genotype to the next with
/// `GenotypeAlleleCounts.increase()`, which looks only at the two lowest alleles. This builds the
/// list from the ordering rule in the module note instead, which produces the same sequence (the
/// first unit test below pins it for a diploid site with three alleles).
///
/// One edge differs: ploidy 0 with no alleles gives one empty genotype here, where the Java's
/// iterator stops at `genotypeCount(0, 0)`, which is 0 (`allele == 0 ? 0 : ...`).
pub fn genotypes_in_canonical_order(ploidy: usize, allele_count: usize) -> Vec<Vec<usize>> {
    // The highest allele decides first, so build by appending it to every genotype of one less
    // ploidy over the alleles up to and including it.
    // (A `fn` declared inside another is visible only there and, unlike a closure, captures
    // nothing; `&mut Vec<..>` lets it append to the caller's list.)
    fn build(ploidy: usize, allele_count: usize, out: &mut Vec<Vec<usize>>) {
        if ploidy == 0 {
            out.push(Vec::new());
            return;
        }
        for highest in 0..allele_count {
            let mut smaller = Vec::new();
            build(ploidy - 1, highest + 1, &mut smaller);
            // `for mut genotype in smaller` takes each list out of `smaller` (which is used up)
            // and lets the body extend it, so no copy is made.
            for mut genotype in smaller {
                genotype.push(highest);
                out.push(genotype);
            }
        }
    }
    let mut out = Vec::new();
    build(ploidy, allele_count, &mut out);
    out
}

/// The allele/count pairs of a genotype given as sorted allele indices, which is what
/// `GenotypeAlleleCounts` stores.
///
/// The pairs come out as the Java's `sortedAlleleCounts` holds them: ascending allele, no zero
/// counts. Pair `i` is what `alleleIndexAt(i)` and `alleleCountAt(i)` read, and the number of
/// pairs is `distinctAlleleCount()`.
pub fn allele_counts_of(genotype: &[usize]) -> Vec<(usize, usize)> {
    let mut out: Vec<(usize, usize)> = Vec::new();
    for allele in genotype {
        // `last_mut()` is a writable reference to the last pair, or `None` for an empty list. The
        // `if` after the pattern is a guard: the arm is taken only when it holds, and `*count += 1`
        // writes through the reference. Any other case, empty list included, starts a new pair.
        match out.last_mut() {
            Some((seen, count)) if seen == allele => *count += 1,
            _ => out.push((*allele, 1)),
        }
    }
    out
}

/// `alleleCountsToIndex(alleleCountArray)`: pairs of allele and count, in any order.
///
/// A genotype is a multiset, so the pairs may come in any order and a count of zero contributes
/// nothing. Only an odd-length array is refused.
pub fn allele_counts_to_index(pairs: &[usize]) -> Result<usize, GenotypeIndexError> {
    // The Java's `(alleleCountArray.length & 1) == 0`.
    if !pairs.len().is_multiple_of(2) {
        return Err(GenotypeIndexError::OddLengthCounts);
    }
    // Expand each pair into `count` copies of its allele, the Java's `alleleContainer`.
    // `chunks(2)` yields consecutive two-element views, the Java's `i += 2` loop.
    let mut alleles = Vec::new();
    for pair in pairs.chunks(2) {
        for _ in 0..pair[1] {
            alleles.push(pair[0]);
        }
    }
    Ok(index_of_sorted_alleles(&mut alleles))
}

/// `calculateIndex`: sort, then sum the first-genotype index of each allele from the highest down.
///
/// ```java
/// Arrays.sort(alleles);
/// return new IndexRange(0, ploidy).sumInt(n ->  {
///     final int allele = alleles[ploidy - n - 1];
///     return (int) indexOfFirstGenotypeWithAllele(ploidy - n, allele);
/// });
/// ```
///
/// The k-th allele from the top is scored at ploidy `ploidy - k`: rule (iii) of the ordering,
/// applied once per allele, each step setting the highest allele aside with one copy of ploidy.
/// See the module note for a worked example.
///
/// `unwrap_or(0)` stands where the Java narrows with `(int)`, and the `usize` sum where `sumInt`
/// adds in an `int`; they agree while every term fits in an `int` (module note on overflow).
/// `sort_unstable` sorts in place and may reorder equal elements, which cannot matter for integers.
fn index_of_sorted_alleles(alleles: &mut [usize]) -> usize {
    alleles.sort_unstable();
    let ploidy = alleles.len();
    (0..ploidy)
        .map(|step| {
            let allele = alleles[ploidy - step - 1];
            index_of_first_genotype_with_allele(ploidy - step, allele).unwrap_or(0)
        })
        .sum()
}

/// `subsettedPLIndices(ploidy, originalAlleles, newAlleles)`.
///
/// `kept` is the index in the original allele list of each allele the new list keeps, in the new
/// list's order. The answer says, for each index of the NEW likelihood vector, which index of the
/// old one it takes its value from, so it is a permutation when nothing is dropped.
///
/// **How**: as `GenotypeIndexCalculator.newToOldGenotypeMap(ploidy, newToOldAlleleMap)` does it,
/// with `kept` as the `newToOldAlleleMap`: walk the NEW genotypes and map each allele back.
/// `AlleleSubsettingUtils.subsettedPLIndices` instead walks every OLD genotype, skips those holding
/// a dropped allele, and records `result[newPLIndex] = oldAlleleCounts.index();`. When `kept` names
/// distinct alleles each new genotype has exactly one old one, so the two fill every slot alike.
/// One difference in what is refused: the Java's walk counts the genotypes over all the original
/// alleles (`GenotypeAlleleCounts.iterable(ploidy, originalAlleles.size())`, which calls
/// `genotypeCount`) and so can refuse a long original list; this counts only the kept ones.
pub fn subsetted_pl_indices(
    ploidy: usize,
    kept: &[usize],
) -> Result<Vec<usize>, GenotypeIndexError> {
    // `?` passes the refusal straight to the caller, Java's `throws` at this call.
    let count = genotype_count(ploidy, kept.len())?;
    let mut out = vec![0; count];
    // `.into_iter().enumerate()` yields (position, genotype) pairs, taking each genotype out of
    // the list, the Java's `newGAC.index()` and `newGAC`.
    for (new_index, genotype) in genotypes_in_canonical_order(ploidy, kept.len())
        .into_iter()
        .enumerate()
    {
        // `kept[*allele]` is `newToOldAlleleMap[newAllele]`; `*` reads the value behind the
        // reference the iterator hands out.
        let mut old: Vec<usize> = genotype.iter().map(|allele| kept[*allele]).collect();
        out[new_index] = index_of_sorted_alleles(&mut old);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_order_is_the_highest_allele_first() {
        let genotypes = genotypes_in_canonical_order(2, 3);
        assert_eq!(
            genotypes,
            vec![
                vec![0, 0],
                vec![0, 1],
                vec![1, 1],
                vec![0, 2],
                vec![1, 2],
                vec![2, 2]
            ]
        );
    }

    #[test]
    fn a_triploid_site_with_four_alleles_has_twenty_genotypes() {
        assert_eq!(genotype_count(3, 4).expect("a count"), 20);
        assert_eq!(genotypes_in_canonical_order(3, 4).len(), 20);
        // And ploidy zero has exactly one, the empty genotype.
        assert_eq!(genotype_count(0, 3).expect("a count"), 1);
        assert_eq!(
            genotypes_in_canonical_order(0, 3),
            vec![Vec::<usize>::new()]
        );
    }

    #[test]
    fn the_pairs_may_come_in_any_order() {
        assert_eq!(allele_counts_to_index(&[1, 1, 2, 1]).expect("an index"), 4);
        assert_eq!(allele_counts_to_index(&[2, 1, 1, 1]).expect("an index"), 4);
        // A count of zero contributes nothing.
        assert_eq!(allele_counts_to_index(&[0, 2, 1, 0]).expect("an index"), 0);
        // And the empty genotype is index 0.
        assert_eq!(allele_counts_to_index(&[]).expect("an index"), 0);
        assert_eq!(
            allele_counts_to_index(&[0, 2, 1]).unwrap_err().message(),
            "the allele counts array cannot have odd length"
        );
    }

    #[test]
    fn subsetting_is_a_permutation_when_nothing_is_dropped() {
        assert_eq!(
            subsetted_pl_indices(2, &[0, 2, 1]).expect("indices"),
            vec![0, 3, 5, 1, 4, 2]
        );
        assert_eq!(
            subsetted_pl_indices(2, &[0, 1, 2]).expect("indices"),
            vec![0, 1, 2, 3, 4, 5]
        );
        // Keeping one alternate of three takes the ref, the het and the hom of THAT allele.
        assert_eq!(
            subsetted_pl_indices(2, &[0, 2]).expect("indices"),
            vec![0, 3, 5]
        );
        assert_eq!(subsetted_pl_indices(2, &[0]).expect("indices"), vec![0]);
    }
}

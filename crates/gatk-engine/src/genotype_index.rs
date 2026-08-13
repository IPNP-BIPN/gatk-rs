//! `GenotypeIndexCalculator` and `GenotypeAlleleCounts`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.genotyper` (GATK 4.6.2.0).
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
//! iterates the obvious way is wrong from the fourth entry on.
//!
//! # The index is a sum of binomials
//!
//! ```java
//! return allele == 0 ? 0 : CombinatoricsUtils.binomialCoefficient(ploidy + allele - 1, allele - 1);
//! ```
//!
//! is the index of the first genotype containing a given allele at a given ploidy, and the index of
//! a genotype is that function summed over its alleles, sorted, from the highest down. The same
//! function at the allele count is the number of genotypes.
//!
//! # Subsetting is a permutation
//!
//! `subsetted_pl_indices` says where each new index takes its likelihood from. Keeping three
//! alleles in the order 0, 2, 1 gives `0,3,5,1,4,2`: nothing is dropped and everything moves.

/// What the combinatorics refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GenotypeIndexError {
    /// `Utils.validateArg((alleleCountArray.length & 1) == 0, ...)`, the only argument check.
    OddLengthCounts,
    /// A binomial too large for the reference's `int`, which it refuses rather than wraps.
    TooManyGenotypes { ploidy: usize, alleles: usize },
}

impl GenotypeIndexError {
    pub fn message(&self) -> String {
        match self {
            GenotypeIndexError::OddLengthCounts => {
                "the allele counts array cannot have odd length".to_string()
            }
            GenotypeIndexError::TooManyGenotypes { ploidy, alleles } => format!(
                "the number of genotypes is too large for ploidy {ploidy} and {alleles} alleles"
            ),
        }
    }

    pub fn java_class(&self) -> &'static str {
        "java.lang.IllegalArgumentException"
    }
}

/// `C(n, k)`, exactly, refusing rather than wrapping.
fn binomial(n: usize, k: usize) -> Option<usize> {
    if k > n {
        return Some(0);
    }
    let k = k.min(n - k);
    let mut result: u128 = 1;
    for step in 0..k {
        result = result.checked_mul((n - step) as u128)? / (step as u128 + 1);
        if result > i32::MAX as u128 {
            return None;
        }
    }
    usize::try_from(result).ok()
}

/// `indexOfFirstGenotypeWithAllele(ploidy, allele)`.
pub fn index_of_first_genotype_with_allele(ploidy: usize, allele: usize) -> Option<usize> {
    if allele == 0 {
        return Some(0);
    }
    binomial(ploidy + allele - 1, allele - 1)
}

/// `genotypeCount(ploidy, alleleCount)`, which is the same function at the allele count.
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
pub fn genotypes_in_canonical_order(ploidy: usize, allele_count: usize) -> Vec<Vec<usize>> {
    // The highest allele decides first, so build by appending it to every genotype of one less
    // ploidy over the alleles up to and including it.
    fn build(ploidy: usize, allele_count: usize, out: &mut Vec<Vec<usize>>) {
        if ploidy == 0 {
            out.push(Vec::new());
            return;
        }
        for highest in 0..allele_count {
            let mut smaller = Vec::new();
            build(ploidy - 1, highest + 1, &mut smaller);
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
pub fn allele_counts_of(genotype: &[usize]) -> Vec<(usize, usize)> {
    let mut out: Vec<(usize, usize)> = Vec::new();
    for allele in genotype {
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
    if !pairs.len().is_multiple_of(2) {
        return Err(GenotypeIndexError::OddLengthCounts);
    }
    let mut alleles = Vec::new();
    for pair in pairs.chunks(2) {
        for _ in 0..pair[1] {
            alleles.push(pair[0]);
        }
    }
    Ok(index_of_sorted_alleles(&mut alleles))
}

/// `calculateIndex`: sort, then sum the first-genotype index of each allele from the highest down.
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
pub fn subsetted_pl_indices(
    ploidy: usize,
    kept: &[usize],
) -> Result<Vec<usize>, GenotypeIndexError> {
    let count = genotype_count(ploidy, kept.len())?;
    let mut out = vec![0; count];
    for (new_index, genotype) in genotypes_in_canonical_order(ploidy, kept.len())
        .into_iter()
        .enumerate()
    {
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

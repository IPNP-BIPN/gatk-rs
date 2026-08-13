//! Conformance for the genotype index combinatorics against GATK 4.6.2.0, compared as every count,
//! every canonical order, every index and every subsetting the golden holds.
//!
//! Golden from `tools/readfilter-conformance/GenotypeIndexDump.java`.
//!
//! # What this suite is for
//!
//!  * **the canonical order**, which is the highest allele first and not the obvious one;
//!  * **the genotype count**, which is a binomial and not a power;
//!  * **a genotype is a multiset**, so the pairs may come in any order;
//!  * **and subsetting is a permutation**, not a filter.

use gatk_corpus as corpus;
use gatk_engine::genotype_index::{
    allele_counts_of, allele_counts_to_index, genotype_count, genotypes_in_canonical_order,
    subsetted_pl_indices,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/genotype_index.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

#[test]
fn every_genotype_count_is_the_reference_s() {
    let text = golden();
    let counts = rows(&text, "count");
    assert!(counts.len() >= 25, "every count is in the golden");
    for row in counts {
        let ploidy: usize = row[0].parse().expect("a ploidy");
        let alleles: usize = row[1].parse().expect("an allele count");
        let expected: usize = row[2].parse().expect("a count");
        assert_eq!(
            genotype_count(ploidy, alleles).expect("a count"),
            expected,
            "count/{ploidy}/{alleles}"
        );
    }
}

/// The order a PL array is indexed by, which is what a rewrite gets wrong from the fourth entry.
#[test]
fn the_canonical_order_is_the_reference_s() {
    let text = golden();
    let mut shapes: Vec<(usize, usize)> = Vec::new();
    for row in rows(&text, "order") {
        let shape = (
            row[0].parse::<usize>().expect("a ploidy"),
            row[1].parse::<usize>().expect("an allele count"),
        );
        if !shapes.contains(&shape) {
            shapes.push(shape);
        }
    }
    assert!(shapes.len() >= 6, "six shapes are in the golden");

    for (ploidy, alleles) in shapes {
        let expected: Vec<String> = rows(&text, "order")
            .into_iter()
            .filter(|row| {
                row[0].parse::<usize>() == Ok(ploidy) && row[1].parse::<usize>() == Ok(alleles)
            })
            .map(|row| row.get(3).copied().unwrap_or("").to_string())
            .collect();

        let ours: Vec<String> = genotypes_in_canonical_order(ploidy, alleles)
            .iter()
            .map(|genotype| {
                allele_counts_of(genotype)
                    .into_iter()
                    .map(|(allele, count)| format!("{allele}:{count}"))
                    .collect::<Vec<_>>()
                    .join(",")
            })
            .collect();
        assert_eq!(ours, expected, "order/{ploidy}/{alleles}");
    }
}

/// Every index the golden holds, including the one whose pairs are reversed.
#[test]
fn every_index_is_the_reference_s() {
    let text = golden();
    for row in rows(&text, "index") {
        let label = row[0];
        let pairs: Vec<usize> = if row[1].is_empty() {
            Vec::new()
        } else {
            row[1]
                .split(',')
                .map(|value| value.parse().expect("a number"))
                .collect()
        };
        let expected: usize = row[2].parse().expect("an index");
        assert_eq!(
            allele_counts_to_index(&pairs).expect("an index"),
            expected,
            "index/{label}"
        );
    }

    // The same genotype with its pairs the other way round is the same index.
    let straight = rows(&text, "index")
        .into_iter()
        .find(|row| row[0] == "diploid-het-non-ref")
        .expect("the row")[2]
        .to_string();
    let reversed = rows(&text, "index")
        .into_iter()
        .find(|row| row[0] == "diploid-het-non-ref-reversed")
        .expect("the row")[2]
        .to_string();
    assert_eq!(straight, reversed);
}

#[test]
fn every_subsetting_is_the_reference_s() {
    let text = golden();
    let subsets = rows(&text, "subset");
    assert!(subsets.len() >= 7, "every subsetting is in the golden");
    for row in subsets {
        let label = row[0];
        let ploidy: usize = row[1].parse().expect("a ploidy");
        let kept: Vec<usize> = row[3]
            .split(',')
            .map(|value| value.parse().expect("an allele index"))
            .collect();
        let expected: Vec<usize> = row[4]
            .split(',')
            .map(|value| value.parse().expect("an index"))
            .collect();
        assert_eq!(
            subsetted_pl_indices(ploidy, &kept).expect("indices"),
            expected,
            "subset/{label}"
        );
    }
}

/// The refusal, which is the only argument check there is.
#[test]
fn an_odd_length_count_array_is_the_reference_s_refusal() {
    let text = golden();
    let expected = rows(&text, "error")
        .into_iter()
        .find(|row| row[0] == "odd-length")
        .expect("the row")[1]
        .to_string();
    let error = allele_counts_to_index(&[0, 2, 1]).unwrap_err();
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        expected
    );
}

//! Conformance for the allele and sample axes, against the oracle.
//!
//! Golden from `tools/genotyper-conformance/AlleleListDump.java`.
//!
//! The rows that decide the port:
//!
//! ```text
//! list  duplicate-separated  2  A*,C          (the third allele repeated the first and vanished)
//! list  ref-flag-pair        2  A*,A          (same bases, different flag, two entries)
//! perm  prefix               true  false      (a subset that kept its order is still permuted)
//! ```

use std::io::Read;

use gatk_engine::allele_list::{AlleleList, PermutationError, SampleList};
use htsjdk_vcf::allele::Allele;

fn golden() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/allele_list.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

fn allele(bases: &str, is_ref: bool) -> Allele {
    Allele::from_str(bases, is_ref).expect("an allele")
}

fn reference() -> Allele {
    allele("A", true)
}

fn reference_bases_as_alt() -> Allele {
    allele("A", false)
}

fn alt1() -> Allele {
    allele("C", false)
}

fn alt2() -> Allele {
    allele("G", false)
}

fn alt3() -> Allele {
    allele("T", false)
}

fn second_reference() -> Allele {
    allele("C", true)
}

/// htsjdk's own rendering: the display string, with a `*` appended for a reference allele.
fn show(allele: &Allele) -> String {
    format!(
        "{}{}",
        allele.display_string(),
        if allele.is_reference() { "*" } else { "" }
    )
}

/// The golden's row for one kind and label, as the remaining tab-separated fields.
fn row(text: &str, kind: &str, label: &str) -> Vec<String> {
    let needle = format!("{kind}\t{label}\t");
    let rest = text
        .lines()
        .find(|line| line.starts_with(&needle))
        .unwrap_or_else(|| panic!("no {kind} row for {label:?}"));
    rest[needle.len()..].split('\t').map(String::from).collect()
}

/// `indexOfAllele` answers `-1` where this port answers `None`.
fn index_or_minus_one(index: Option<usize>) -> i64 {
    index.map_or(-1, |i| i as i64)
}

fn list_cases() -> Vec<(&'static str, Vec<Allele>)> {
    let (r, ru, a1, a2) = (reference(), reference_bases_as_alt(), alt1(), alt2());
    let n = Allele::no_call();
    let span_del = Allele::from_str("*", false).expect("the span-del allele");
    vec![
        ("empty", vec![]),
        ("one-ref", vec![r.clone()]),
        ("ref-and-alt", vec![r.clone(), a1.clone()]),
        ("alt-first", vec![a1.clone(), r.clone()]),
        ("no-reference", vec![a1.clone(), a2]),
        ("duplicate-adjacent", vec![r.clone(), r.clone(), a1.clone()]),
        (
            "duplicate-separated",
            vec![r.clone(), a1.clone(), r.clone()],
        ),
        ("all-duplicates", vec![a1.clone(), a1.clone(), a1.clone()]),
        ("ref-flag-pair", vec![r.clone(), ru]),
        ("two-references", vec![r.clone(), second_reference()]),
        ("with-no-call", vec![r.clone(), n, a1.clone()]),
        ("with-span-del", vec![r, span_del, a1]),
    ]
}

fn sample_cases() -> Vec<(&'static str, Vec<&'static str>)> {
    vec![
        ("empty", vec![]),
        ("one", vec!["s1"]),
        ("two", vec!["s1", "s2"]),
        ("duplicate", vec!["s1", "s2", "s1"]),
        ("out-of-order", vec!["b", "A", "a"]),
    ]
}

fn permutation_cases() -> Vec<(&'static str, Vec<Allele>, Vec<Allele>)> {
    let (r, ru, a1, a2, a3) = (
        reference(),
        reference_bases_as_alt(),
        alt1(),
        alt2(),
        alt3(),
    );
    let long = allele("AA", false);
    let original = vec![r.clone(), a1.clone(), a2.clone(), a3.clone()];
    vec![
        ("identity", original.clone(), original.clone()),
        (
            "reordered",
            original.clone(),
            vec![a3.clone(), a2.clone(), a1.clone(), r.clone()],
        ),
        (
            "swap-two",
            original.clone(),
            vec![r.clone(), a2.clone(), a1.clone(), a3.clone()],
        ),
        (
            "drop-last",
            original.clone(),
            vec![r.clone(), a1.clone(), a2.clone()],
        ),
        (
            "drop-first",
            original.clone(),
            vec![a1.clone(), a2.clone(), a3],
        ),
        ("keep-one", original.clone(), vec![a2.clone()]),
        ("keep-none", original.clone(), vec![]),
        ("prefix", original.clone(), vec![r.clone(), a1.clone()]),
        (
            "longer-target",
            original.clone(),
            vec![r.clone(), a1.clone(), a2.clone(), alt3(), long.clone()],
        ),
        ("unknown-allele", original.clone(), vec![r, long]),
        ("wrong-ref-flag", original.clone(), vec![ru, a1.clone()]),
        ("duplicate-in-target", original, vec![a1.clone(), a1, a2]),
        ("empty-to-empty", vec![], vec![]),
        ("empty-to-one", vec![], vec![reference()]),
    ]
}

#[test]
fn every_allele_list_matches_the_reference() {
    let text = golden();
    for (label, alleles) in list_cases() {
        let list = AlleleList::new(&alleles);
        let expected = row(&text, "list", label);
        let rendered = list
            .as_slice()
            .iter()
            .map(show)
            .collect::<Vec<_>>()
            .join(",");
        assert_eq!(
            list.number_of_alleles().to_string(),
            expected[0],
            "size for {label}"
        );
        assert_eq!(rendered, expected[1], "contents for {label}");
        assert_eq!(
            index_or_minus_one(list.index_of_reference()).to_string(),
            expected[2],
            "indexOfReference for {label}"
        );
    }
}

#[test]
fn every_index_lookup_matches_the_reference() {
    let text = golden();
    let list = AlleleList::new(&[reference(), reference_bases_as_alt(), alt1()]);
    let queries: Vec<(&str, Allele)> = vec![
        ("A*", reference()),
        ("A", reference_bases_as_alt()),
        ("C", alt1()),
        ("G", alt2()),
        ("no-call", Allele::no_call()),
    ];
    for (query, allele) in queries {
        let needle = format!("index\tref-flag-pair\t{query}\t");
        let line = text
            .lines()
            .find(|line| line.starts_with(&needle))
            .unwrap_or_else(|| panic!("no index row for {query}"));
        let fields: Vec<&str> = line[needle.len()..].split('\t').collect();
        assert_eq!(
            index_or_minus_one(list.index_of_allele(&allele)).to_string(),
            fields[0],
            "indexOfAllele({query})"
        );
        assert_eq!(
            list.contains_allele(&allele).to_string(),
            fields[1],
            "containsAllele({query})"
        );
    }
}

#[test]
fn every_sample_list_matches_the_reference() {
    let text = golden();
    for (label, names) in sample_cases() {
        let owned: Vec<String> = names.iter().map(|n| (*n).to_string()).collect();
        let list = SampleList::new(&owned);
        let expected = row(&text, "samples", label);
        assert_eq!(
            list.number_of_samples().to_string(),
            expected[0],
            "size for {label}"
        );
        assert_eq!(
            list.as_slice().join(","),
            expected[1],
            "contents for {label}"
        );
    }
}

#[test]
fn every_permutation_matches_the_reference() {
    let text = golden();
    for (label, from, to) in permutation_cases() {
        let original = AlleleList::new(&from);
        let target = AlleleList::new(&to);
        let expected = row(&text, "perm", label);

        match original.permutation(&target) {
            Err(error) => {
                assert_eq!(
                    format!("E:{}:{}", error.class(), error.message()),
                    expected[0],
                    "the refusal for {label}"
                );
            }
            Ok(permutation) => {
                let from_indices = (0..permutation.to_size())
                    .map(|i| index_or_minus_one(permutation.from_index(i)).to_string())
                    .collect::<Vec<_>>()
                    .join(",");
                let to_indices = (0..permutation.from_size())
                    .map(|i| index_or_minus_one(permutation.to_index(i)).to_string())
                    .collect::<Vec<_>>()
                    .join(",");
                let kept = (0..permutation.from_size())
                    .map(|i| permutation.is_kept(i).to_string())
                    .collect::<Vec<_>>()
                    .join(",");
                let ours = vec![
                    permutation.is_partial().to_string(),
                    permutation.is_non_permuted().to_string(),
                    permutation.from_size().to_string(),
                    permutation.to_size().to_string(),
                    from_indices,
                    to_indices,
                    kept,
                ];
                assert_eq!(ours, expected, "the permutation for {label}");
            }
        }
    }
}

/// The rows a port gets wrong by treating an indexed list as a list, or a subset as an identity.
#[test]
fn the_rows_that_a_list_gets_wrong() {
    let text = golden();

    // A duplicate is dropped and nothing says so, wherever it sits.
    assert_eq!(row(&text, "list", "duplicate-separated")[0], "2");
    assert_eq!(row(&text, "list", "all-duplicates")[0], "1");

    // Same bases, different reference flag: two entries, and only the first is the reference.
    assert_eq!(row(&text, "list", "ref-flag-pair")[0], "2");
    assert_eq!(row(&text, "list", "two-references")[2], "0");

    // A subset that kept its order is partial and NOT non-permuted.
    let prefix = row(&text, "perm", "prefix");
    assert_eq!(prefix[0], "true", "isPartial");
    assert_eq!(prefix[1], "false", "isNonPermuted");
    // While the identity is neither partial nor permuted.
    let identity = row(&text, "perm", "identity");
    assert_eq!(identity[0], "false");
    assert_eq!(identity[1], "true");

    // A target whose duplicate its own constructor collapsed is a permutation, not a refusal.
    assert!(!row(&text, "perm", "duplicate-in-target")[0].starts_with("E:"));

    // And both refusals carry the same message by two different routes.
    assert_eq!(
        row(&text, "perm", "longer-target")[0],
        row(&text, "perm", "unknown-allele")[0]
    );
    assert_eq!(
        PermutationError::TargetLonger.message(),
        PermutationError::AlleleNotInOriginal.message()
    );
}

//! Conformance for `PathSeqBuildReferenceTaxonomy` against GATK 4.6.2.0, compared as the tree the
//! database holds and the map beside it.
//!
//! Golden from `tools/readfilter-conformance/PathSeqBuildReferenceTaxonomyDump.java`.
//!
//! # What this suite is for
//!
//!  * **the map being keyed by the contig name rather than the accession**;
//!  * **the two forms a reference name can place a contig by**, and the fallback when it has
//!    neither;
//!  * **the two catalog formats**, and the blank line that ends one silently;
//!  * **the trim to the taxa the reference holds**, and the totals the nodes keep;
//!  * **and the length filter, which is on the map and not on the tree**.
//!
//! The contig names are the ones htsjdk put in the dictionary, which are the fasta names cut at the
//! first space: that is what the tool saw.

use gatk_corpus as corpus;
use gatk_tools::pathseq_taxonomy::{build, PsTree, TaxonomyError};
use std::collections::BTreeMap;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/pathseq_taxonomy.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn fixture(text: &str, name: &str) -> String {
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(&format!("fixture\t{name}=")))
            .unwrap_or_else(|| panic!("the golden carries the {name} fixture")),
    )
}

/// The tree of one run, as the dump printed it: id, name, parent, rank, length.
fn tree_rows(text: &str, label: &str) -> Vec<(i32, String, i32, String, i64)> {
    let prefix = format!("tree\t{label}\t");
    text.lines()
        .filter_map(|line| line.strip_prefix(prefix.as_str()))
        .map(|rest| {
            let fields: Vec<&str> = rest.split('\t').collect();
            (
                fields[0].parse().expect("a tax id"),
                fields[1].to_string(),
                fields[2].parse().expect("a parent"),
                fields[3].to_string(),
                fields[4].parse().expect("a length"),
            )
        })
        .collect()
}

fn accession_rows(text: &str, label: &str) -> BTreeMap<String, i32> {
    let prefix = format!("accession\t{label}\t");
    text.lines()
        .filter_map(|line| line.strip_prefix(prefix.as_str()))
        .map(|rest| {
            let (name, id) = rest.rsplit_once('=').expect("an accession row");
            (name.to_string(), id.parse().expect("a tax id"))
        })
        .collect()
}

fn refusal(text: &str, label: &str) -> String {
    text.lines()
        .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
        .map(unescape)
        .expect("the golden carries the refusal")
}

fn ours(tree: &PsTree) -> Vec<(i32, String, i32, String, i64)> {
    tree.node_ids()
        .into_iter()
        .map(|id| {
            (
                id,
                tree.name_of(id).unwrap_or("null").to_string(),
                tree.parent_of(id),
                tree.rank_of(id).unwrap_or("null").to_string(),
                tree.length_of(id),
            )
        })
        .collect()
}

/// The reference the dump built, as htsjdk's dictionary holds it: every name cut at the space.
fn contigs() -> Vec<(String, i64)> {
    vec![
        ("ref|NC_VIRUS.1|".to_string(), 300),
        ("ref|NC_BACT.1|".to_string(), 1000),
        ("ref|NC_SHORT.1|".to_string(), 100),
        ("taxid|562|".to_string(), 800),
        ("ACC_PLAIN.1".to_string(), 900),
        ("gi|9|ref|NC_BOTH.1|taxid|11234|".to_string(), 700),
    ]
}

fn check(
    text: &str,
    label: &str,
    refseq: Option<&str>,
    genbank: Option<&str>,
    min_length: i64,
    contigs: &[(String, i64)],
) {
    let (tree, map, _) = build(
        contigs,
        refseq,
        genbank,
        &fixture(text, "names.dmp"),
        &fixture(text, "nodes.dmp"),
        min_length,
    )
    .expect("a database");
    assert_eq!(ours(&tree), tree_rows(text, label), "{label}: the tree");
    assert_eq!(map, accession_rows(text, label), "{label}: the map");
}

#[test]
fn both_catalogs() {
    let text = golden();
    check(
        &text,
        "both-catalogs",
        Some(&fixture(&text, "refseq.catalog")),
        Some(&fixture(&text, "genbank.catalog")),
        0,
        &contigs(),
    );
    // The map is keyed by the whole contig name, not by the accession the catalog was searched
    // with, so no entry reads NC_BACT.1 on its own.
    let keys = accession_rows(&text, "both-catalogs");
    assert!(keys.contains_key("ref|NC_BACT.1|"));
    assert!(!keys.contains_key("NC_BACT.1"));
}

#[test]
fn the_length_filter_is_on_the_map_and_not_the_tree() {
    let text = golden();
    check(
        &text,
        "min-length-500",
        Some(&fixture(&text, "refseq.catalog")),
        Some(&fixture(&text, "genbank.catalog")),
        500,
        &contigs(),
    );
    let unfiltered = accession_rows(&text, "both-catalogs");
    let filtered = accession_rows(&text, "min-length-500");
    // The hundred-base bacterium goes, the three-hundred-base virus stays.
    assert!(unfiltered.contains_key("ref|NC_SHORT.1|"));
    assert!(!filtered.contains_key("ref|NC_SHORT.1|"));
    assert!(filtered.contains_key("ref|NC_VIRUS.1|"));
    // And the tree is untouched: E. coli still counts the contig the map cannot reach.
    let node = tree_rows(&text, "min-length-500")
        .into_iter()
        .find(|row| row.0 == 562)
        .expect("the E. coli node");
    assert_eq!(node.4, 1900);
}

#[test]
fn each_catalog_alone() {
    let text = golden();
    check(
        &text,
        "refseq-only",
        Some(&fixture(&text, "refseq.catalog")),
        None,
        0,
        &contigs(),
    );
    check(
        &text,
        "genbank-only",
        None,
        Some(&fixture(&text, "genbank.catalog")),
        0,
        &contigs(),
    );
}

#[test]
fn a_blank_line_ends_a_catalog() {
    let text = golden();
    check(
        &text,
        "blank-line-truncates",
        Some(&fixture(&text, "truncated.catalog")),
        None,
        0,
        &contigs(),
    );
    // Everything after the blank line is invisible, so the bacterium never lands anywhere.
    let map = accession_rows(&text, "blank-line-truncates");
    assert!(!map.contains_key("ref|NC_BACT.1|"));
    assert!(map.contains_key("ref|NC_VIRUS.1|"));
}

#[test]
fn neither_catalog_is_refused_first() {
    let text = golden();
    let error = build(
        &contigs(),
        None,
        None,
        &fixture(&text, "names.dmp"),
        &fixture(&text, "nodes.dmp"),
        0,
    )
    .expect_err("a refusal");
    assert_eq!(error, TaxonomyError::NoCatalog);
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "no-catalog")
    );
}

#[test]
fn a_narrow_catalog_line_names_genbank_whatever_the_format() {
    let text = golden();
    let error = build(
        &contigs(),
        Some(&fixture(&text, "narrow.catalog")),
        None,
        &fixture(&text, "names.dmp"),
        &fixture(&text, "nodes.dmp"),
        0,
    )
    .expect_err("a refusal");
    assert_eq!(
        error,
        TaxonomyError::TooFewColumns {
            expected: 3,
            found: 2,
            line: 1
        }
    );
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "narrow-catalog")
    );
}

#[test]
fn a_taxon_id_that_is_not_a_number() {
    let text = golden();
    let error = build(
        &[("taxid|abc|".to_string(), 100)],
        Some(&fixture(&text, "refseq.catalog")),
        None,
        &fixture(&text, "names.dmp"),
        &fixture(&text, "nodes.dmp"),
        0,
    )
    .expect_err("a refusal");
    assert_eq!(
        error,
        TaxonomyError::NotAnInteger {
            value: "abc".to_string()
        }
    );
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "bad-taxon-id")
    );
}

#[test]
fn a_reference_no_catalog_knows() {
    let text = golden();
    let error = build(
        &[("ref|NC_NOWHERE.1|".to_string(), 100)],
        Some(&fixture(&text, "refseq.catalog")),
        None,
        &fixture(&text, "names.dmp"),
        &fixture(&text, "nodes.dmp"),
        0,
    )
    .expect_err("a refusal");
    assert_eq!(error, TaxonomyError::NoRelevantTaxa);
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "no-relevant-taxa")
    );
}

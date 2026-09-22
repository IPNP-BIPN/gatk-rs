//! Conformance for the file `PathSeqBuildReferenceTaxonomy` writes, against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/PathSeqTaxonomyKryoDump.java`, which ran the tool in
//! the pinned container and printed the database file as hex, and built one `HashMap` per probe
//! with the constructor the tool uses and printed the order it iterated in.
//!
//! # What this suite is for
//!
//!  * **the probes**: every map shape the file's order passes through -- the default table across
//!    both of its growth points, a table sized at construction down to one bucket,
//!    `new HashSet<>(collection)` and the removals after it, a table grown and then trimmed -- over
//!    integer ids and over contig names, replayed through [`JavaHashMap`];
//!  * **and the runs**: the whole database, byte for byte, for a small taxonomy and for a wide one
//!    built to push the node map, a genus's children, a taxon's contigs and the output map past a
//!    growth point, with a branch trimmed and an unreachable subtree removed.
//!
//! The contig names are the dictionary's, as the dump printed them, which is what the tool read.

use gatk_corpus as corpus;
use gatk_engine::java_hash::{JavaHashCode, JavaHashMap};
use gatk_tools::{pathseq_kryo, pathseq_taxonomy};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/pathseq_taxonomy_kryo.txt.gz"),
    )
}

/// `ReferenceQueryDump.escape` undone: a backslash, a tab or a newline.
fn unescape(text: &str) -> String {
    let mut out = String::new();
    let mut chars = text.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
            Some(other) => out.push(other),
            None => out.push('\\'),
        }
    }
    out
}

fn fixture(text: &str, name: &str) -> String {
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(&format!("fixture\t{name}=")))
            .unwrap_or_else(|| panic!("the golden carries the {name} fixture")),
    )
}

fn value(text: &str, kind: &str, label: &str) -> Option<String> {
    let prefix = format!("{kind}\t{label}=");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .map(str::to_string)
}

/// The reference as the dictionary held it.
fn contigs(text: &str, reference: &str) -> Vec<(String, i64)> {
    let prefix = format!("contig\t{reference}\t");
    let rows: Vec<(String, i64)> = text
        .lines()
        .filter_map(|line| line.strip_prefix(prefix.as_str()))
        .map(|rest| {
            let (_, row) = rest.split_once('=').expect("a contig row");
            let (name, length) = row.rsplit_once('\t').expect("a name and a length");
            (name.to_string(), length.parse().expect("a length"))
        })
        .collect();
    assert!(
        !rows.is_empty(),
        "the golden carries the {reference} reference"
    );
    rows
}

fn split(list: &str) -> Vec<&str> {
    if list.is_empty() {
        Vec::new()
    } else {
        list.split(',').collect()
    }
}

/// One probe, replayed.
fn replay<K: JavaHashCode + PartialEq + Clone + ToString>(
    constructor: &str,
    inserted: &[K],
    removed: &[K],
) -> Vec<String> {
    let mut map: JavaHashMap<K, ()> = match constructor {
        "default" => JavaHashMap::new(),
        "copy" => JavaHashMap::copy_of(inserted.iter()),
        sized => JavaHashMap::with_capacity(
            sized
                .strip_prefix("sized:")
                .and_then(|n| n.parse().ok())
                .unwrap_or_else(|| panic!("a constructor: {sized}")),
        ),
    };
    if constructor != "copy" {
        for key in inserted {
            map.insert(key.clone(), ());
        }
    }
    for key in removed {
        map.remove(key);
    }
    map.check().expect("a measured shape");
    map.keys().map(ToString::to_string).collect()
}

#[test]
fn every_probe_iterates_in_the_reference_s_order() {
    let text = golden();
    let mut seen = 0;
    for line in text.lines().filter(|line| line.starts_with("probe\t")) {
        let (row, order) = line.rsplit_once('=').expect("a probe row");
        let fields: Vec<&str> = row.split('\t').collect();
        let [_, label, kind, constructor, inserted, removed] = fields[..] else {
            panic!("a probe row: {line}");
        };
        let ours = match kind {
            "int" => {
                let parse = |list: &str| -> Vec<i32> {
                    split(list)
                        .iter()
                        .map(|key| key.parse().expect("an int"))
                        .collect()
                };
                replay(constructor, &parse(inserted), &parse(removed))
            }
            "string" => {
                let parse = |list: &str| -> Vec<String> {
                    split(list).iter().map(|key| key.to_string()).collect()
                };
                replay(constructor, &parse(inserted), &parse(removed))
            }
            other => panic!("a probe kind: {other}"),
        };
        assert_eq!(ours, split(order), "{label}");
        seen += 1;
    }
    assert_eq!(seen, 18, "every probe the dump writes");
}

fn check_run(
    text: &str,
    label: &str,
    fixtures: &str,
    refseq: bool,
    genbank: bool,
    min_length: i64,
) {
    let refseq = refseq.then(|| fixture(text, &format!("{fixtures}-refseq.catalog")));
    let genbank = genbank.then(|| fixture(text, &format!("{fixtures}-genbank.catalog")));
    let (tree, map, _) = pathseq_taxonomy::build(
        &contigs(text, fixtures),
        refseq.as_deref(),
        genbank.as_deref(),
        &fixture(text, &format!("{fixtures}-names.dmp")),
        &fixture(text, &format!("{fixtures}-nodes.dmp")),
        min_length,
    )
    .expect("a database");
    let file = pathseq_kryo::taxonomy_database_file(&tree, &map).expect("a measured order");
    let hex: String = file.iter().map(|byte| format!("{byte:02x}")).collect();
    let expected = value(text, "stream", label).unwrap_or_else(|| panic!("{label}: a stream"));
    assert_eq!(
        value(text, "len", label).as_deref(),
        Some(file.len().to_string().as_str()),
        "{label}: the length"
    );
    assert_eq!(hex, expected, "{label}: the bytes");
}

#[test]
fn the_small_taxonomy_is_the_reference_s_file() {
    let text = golden();
    check_run(&text, "small-both-catalogs", "small", true, true, 0);
    check_run(&text, "small-min-length-500", "small", true, true, 500);
    check_run(&text, "small-refseq-only", "small", true, false, 0);
}

#[test]
fn the_wide_taxonomy_is_the_reference_s_file() {
    let text = golden();
    check_run(&text, "wide-both-catalogs", "wide", true, true, 0);
    check_run(&text, "wide-min-length-500", "wide", true, true, 500);
    check_run(&text, "wide-refseq-only", "wide", true, false, 0);
}

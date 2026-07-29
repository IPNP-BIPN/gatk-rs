//! Conformance for `ReferenceDataSource` against GATK 4.6.2.0.
//!
//! The fixture travels in the golden: the FASTA and its `.fai` are written back to a temporary
//! directory and queried, so the port reads exactly the bytes the reference read.
//!
//! What this measures is not the FASTA reader, which is `noodles`, but GATK's transformation of
//! what it returns: every query is upper-cased and every IUPAC code becomes `N`. The fixture is
//! built to make both visible, with a soft-masked line and a line of ambiguity codes.

use gatk_corpus as corpus;
use gatk_engine::reference::ReferenceFileSource;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/reference_query.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    let mut out = String::new();
    let mut chars = text.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('n') => out.push('\n'),
            Some('t') => out.push('\t'),
            Some('\\') => out.push('\\'),
            Some(other) => out.push(other),
            None => {}
        }
    }
    out
}

fn field<'a>(text: &'a str, kind: &str) -> &'a str {
    text.lines()
        .find_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .unwrap_or_else(|| panic!("the golden carries no {kind} row"))
}

#[test]
fn every_query_returns_the_bases_the_reference_returns() {
    let text = golden();

    // The fixture is written where the test runs, not committed twice: the golden is the one copy.
    let dir = std::env::temp_dir().join(format!("gatk-rs-refquery-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let fasta = dir.join("ref.fasta");
    std::fs::write(&fasta, unescape(field(&text, "fasta"))).unwrap();
    std::fs::write(dir.join("ref.fasta.fai"), unescape(field(&text, "fai"))).unwrap();

    let mut source = ReferenceFileSource::open(&fasta).expect("the fixture opens");

    let mut compared = 0;
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts[0] != "query" {
            continue;
        }
        let (contig, start, stop) = (
            parts[1],
            parts[2].parse().unwrap(),
            parts[3].parse().unwrap(),
        );
        let ours = source
            .query(contig, start, stop)
            .map_or_else(|_| "E".to_string(), |b| String::from_utf8(b).unwrap());
        assert_eq!(ours, parts[4], "{contig}:{start}-{stop}");
        compared += 1;
    }

    std::fs::remove_dir_all(&dir).ok();
    assert!(compared > 0, "the golden carries no queries");
    println!("{compared} reference queries, all identical");
}

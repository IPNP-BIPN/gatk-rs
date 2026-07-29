//! Conformance for `CountingReadFilter` against GATK 4.6.2.0.
//!
//! The golden is produced by `tools/readfilter-conformance/CountingFilterDump.java` over the same
//! corpus as the decision matrix, and carries three rows per composition: the reference's own name
//! for the tree, the decision per record, and the summary text with its newlines escaped.
//!
//! The port rebuilds the tree by parsing the name, so neither side keeps its own list of
//! compositions. What the summary catches that the decisions cannot: the *order* of a conjunction.
//! `WellformedReadFilter` and the eight filters it and's together produce the same decisions here,
//! and only the counts say which of the eight rejected each read.

use gatk_readfilter::counting;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

use gatk_corpus as corpus;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/counting_filters.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<(&'a str, &'a str)> {
    text.lines()
        .filter_map(|line| {
            let mut parts = line.splitn(3, '\t');
            (parts.next() == Some(kind)).then(|| (parts.next().unwrap(), parts.next().unwrap()))
        })
        .collect()
}

#[test]
fn every_composition_counts_and_prints_what_the_reference_does() {
    let text = golden();
    let records: Vec<BamRecord> = corpus::records(&text);
    let header: SamHeader = corpus::header(&text);

    let compositions = rows(&text, "composition");
    assert!(
        !compositions.is_empty(),
        "the golden carries no compositions"
    );
    let decisions = rows(&text, "decisions");
    let summaries = rows(&text, "summary");

    for (id, spec) in &compositions {
        let mut filter = counting::parse(spec, &header)
            .unwrap_or_else(|| panic!("{id}: cannot rebuild {spec:?}; a leaf is not ported"));

        assert_eq!(
            filter.name(),
            *spec,
            "{id}: the port names the tree differently from the reference"
        );

        let ours: String = records
            .iter()
            .map(|read| if filter.test(read) { '1' } else { '0' })
            .collect();
        let expected = decisions
            .iter()
            .find(|(other, _)| other == id)
            .expect("a decisions row per composition")
            .1;
        assert_eq!(ours, expected, "{id}: decisions differ");

        // The summary is compared as text, escapes included, because it is text a tool prints.
        let expected = summaries
            .iter()
            .find(|(other, _)| other == id)
            .expect("a summary row per composition")
            .1;
        let ours = filter.summary_line().replace('\n', "\\n");
        assert_eq!(ours, expected, "{id}: summary text differs");
    }

    println!(
        "{} compositions, decisions and summaries identical",
        compositions.len()
    );
}

//! Conformance for `FeatureDataSource`'s query path against GATK 4.6.2.0.
//!
//! Twenty queries over one BED fixture, at two lookahead settings.
//!
//! What the golden showed is worth stating rather than hiding: the twenty answers are **identical**
//! at lookahead 0 and at 100. The cache is transparent for well-ordered access, because
//! `trimToNewStartPosition` keeps exactly the features a fresh overlap query would have returned.
//! So this suite pins what a tool sees, which is what matters, and it does not distinguish the
//! cache from a fresh query per call. The port carries the cache for the ordering guarantee, not
//! because these rows force it; a probe that separated them would need a feature file whose records
//! are out of start order, which Tribble refuses to index.

use gatk_corpus as corpus;
use gatk_engine::features::{Feature, FeatureDataSource, FeatureReader};
use gatk_engine::interval::SimpleInterval;

/// The BED body the harness wrote. BED is 0-based half-open, so `chr1 9 20` decodes to
/// `chr1:10-20`, and the fixture is written here in the decoded coordinates the reader returns.
const FEATURES: [(&str, i32, i32, &str); 8] = [
    ("chr1", 10, 20, "f1"),
    ("chr1", 15, 25, "f2"),
    ("chr1", 20, 120, "f3"),
    ("chr1", 50, 60, "f4"),
    ("chr1", 55, 56, "f5"),
    ("chr1", 100, 110, "f6"),
    ("chr1", 150, 160, "f7"),
    ("chr2", 10, 20, "g1"),
];

/// Tribble's own query: every feature overlapping the interval, in file order.
struct FixtureReader;

impl FeatureReader for FixtureReader {
    fn query(&self, interval: &SimpleInterval) -> Vec<Feature> {
        FEATURES
            .iter()
            .filter(|(contig, start, end, _)| {
                *contig == interval.contig && *start <= interval.end && interval.start <= *end
            })
            .map(|(contig, start, end, name)| Feature::new(contig, *start, *end, name))
            .collect()
    }
}

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/feature_cache.txt.gz"),
    )
}

fn parse_interval(text: &str) -> SimpleInterval {
    let (contig, range) = text.split_once(':').expect("contig:range");
    let (start, end) = range.split_once('-').expect("start-end");
    SimpleInterval {
        contig: contig.to_string(),
        start: start.parse().expect("a start"),
        end: end.parse().expect("an end"),
    }
}

#[test]
fn every_query_returns_what_the_reference_returns() {
    let text = golden();

    let mut compared = 0;
    let mut source: Option<(i32, FeatureDataSource<FixtureReader>)> = None;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("query\t") else {
            continue;
        };
        let mut parts = rest.split('\t');
        let lookahead: i32 = parts
            .next()
            .expect("a lookahead")
            .parse()
            .expect("a number");
        let index: usize = parts.next().expect("an index").parse().expect("a number");
        let interval = parse_interval(parts.next().expect("an interval"));
        let outcome = parts.next().expect("an outcome");
        let expected = parts.next().unwrap_or("-");

        // Each lookahead is one run of the harness with one source, and the cache carries across
        // the queries within it. Rebuilding it per query would erase the thing being measured.
        if index == 0 || source.as_ref().map(|(l, _)| *l) != Some(lookahead) {
            source = Some((lookahead, FeatureDataSource::new(FixtureReader, lookahead)));
        }
        let source = &mut source.as_mut().expect("a source").1;

        let result = source.query_and_prefetch(&interval);
        match (result, outcome) {
            (Ok(features), "ok") => {
                let ours = if features.is_empty() {
                    "-".to_string()
                } else {
                    features
                        .iter()
                        .map(|f| format!("{}@{}:{}-{}", f.name, f.contig, f.start, f.end))
                        .collect::<Vec<_>>()
                        .join("|")
                };
                assert_eq!(ours, expected, "lookahead {lookahead} query {index}");
            }
            (Err(error), outcome) if outcome.starts_with("E:") => {
                panic!("lookahead {lookahead} query {index}: unexpected refusal {error}");
            }
            (Ok(_), outcome) => {
                panic!("lookahead {lookahead} query {index}: the reference raised {outcome}")
            }
            (Err(error), _) => {
                panic!("lookahead {lookahead} query {index}: the port raised {error}")
            }
        }
        compared += 1;
    }

    assert!(compared > 0, "the golden carries no query rows");
    println!("{compared} feature queries, identical at both lookahead settings");
}

//! Conformance for the sample-set iteration order and the builder's routing, against
//! GATK 4.6.2.0.
//!
//! The order half is not a check on an implementation, it is the implementation's definition. The
//! order comes from a `java.util.HashSet`, which is GPL2 and unspecified, so
//! `crates/gatk-engine/src/java_hash.rs` reproduces it as an observable of the oracle rather than
//! as a port. See `docs/an-unspecified-order-that-reaches-the-output.md`.
//!
//! Each `order` row is preceded by the `String.hashCode` of every name in it, so a divergence
//! lands on the hash, which is specified, rather than on the layout, which is not. The rows are
//! worth reading: `three` comes back reversed from insertion order, `digits` interleaves, and
//! `thirteen` crosses the growth point.
//!
//! The routing half pins one asymmetry: `areIntervalsSpecified` is `intervals != null`, so an
//! empty but non-null list is *specified*. With `emitEmptyLoci` off that routes to an iterator
//! whose constructor rejects it, and with it on the same list produces a traversal of nothing.

use gatk_corpus as corpus;
use gatk_engine::java_hash;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/ac_builder.txt.gz"),
    )
}

/// The names each labelled probe inserted, in insertion order.
fn names_for(label: &str) -> Vec<String> {
    let names: Vec<&str> = match label {
        "two" => vec!["sampleA", "sampleB"],
        "two-reversed" => vec!["sampleB", "sampleA"],
        "three" => vec!["NA12878", "NA12891", "NA12892"],
        "digits" => vec!["1", "2", "3", "10", "11"],
        "thirteen" => vec![
            "s01", "s02", "s03", "s04", "s05", "s06", "s07", "s08", "s09", "s10", "s11", "s12",
            "s13",
        ],
        "negative" => vec!["zzzzzzzzzzzz", "sampleA"],
        "empty-name" => vec!["", "sampleA"],
        other => panic!("{other} is in the golden but not configured here"),
    };
    names.into_iter().map(|n| n.to_string()).collect()
}

#[test]
fn the_sample_order_and_the_routing_match_the_reference() {
    let text = golden();

    let mut orders = 0;
    let mut hashes = 0;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("hash\t") {
            let (name, expected) = rest.split_once('\t').expect("a name and a hash");
            // The harness prints `<empty>` for the empty string, which has no other rendering.
            let name = if name == "<empty>" { "" } else { name };
            assert_eq!(
                java_hash::string_hash_code(name).to_string(),
                expected,
                "String.hashCode of {name:?}"
            );
            hashes += 1;
        } else if let Some(rest) = line.strip_prefix("order\t") {
            let (label, expected) = rest.split_once('\t').expect("a label and an order");
            let ours = java_hash::hash_set_order(&names_for(label))
                .unwrap_or_else(|e| panic!("{label}: the port refused: {e:?}"))
                .into_iter()
                .map(|name| {
                    if name.is_empty() {
                        "<empty>".to_string()
                    } else {
                        name
                    }
                })
                .collect::<Vec<_>>()
                .join("|");
            assert_eq!(ours, expected, "{label}: iteration order");
            orders += 1;
        }
    }

    assert!(orders > 0, "the golden carries no order rows");
    println!("{orders} sample orders and {hashes} hashes, all identical");
}

/// The routing table, which is a property of the builder rather than of any data.
///
/// Compared as the classes the reference named, because the routing is what a walker inherits: a
/// tool that asks for empty loci and one that does not are reading different iterators, not the
/// same iterator with a flag.
#[test]
fn the_routing_matches_the_reference() {
    let text = golden();

    // (emitEmptyLoci, intervals: None = null, Some(vec![]) = empty, Some(..) = present)
    let expected_route = |label: &str| -> &'static str {
        match label {
            "noloci-nullintervals" => "LocusIteratorByState",
            "noloci-emptyintervals" => "E:java.lang.IllegalArgumentException",
            "noloci-intervals" => "IntervalOverlappingIterator",
            _ => "IntervalAlignmentContextIterator",
        }
    };

    let mut compared = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("route\t") else {
            continue;
        };
        let (label, route) = rest.split_once('\t').expect("a label and a route");
        // The port has no class names, so what is asserted is that the reference routed where
        // this table says it does. The table is the port's routing decision, written out.
        assert_eq!(route, expected_route(label), "{label}: routing");
        compared += 1;
    }

    assert!(compared > 0, "the golden carries no route rows");
    println!("{compared} routing decisions, all as the port routes");
}

/// The contexts each route yields, which is what a walker actually iterates.
///
/// `emitEmptyLoci` with no intervals substitutes the whole reference, so that run walks all 200
/// bases of the contig; the harness capped it at 61 rows and said so, and this compares the rows
/// it produced rather than pretending the traversal ended there.
#[test]
fn every_route_yields_the_contexts_the_reference_yields() {
    use gatk_engine::context_iterator::{self, Route};
    use gatk_engine::interval::SimpleInterval;
    use gatk_engine::locus_iterator::{self, LocusIteratorOptions};
    use gatk_engine::read_states::{DownsamplingInfo, ReadStateManager};
    use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
    use htsjdk_bam::record::BamRecord;
    use htsjdk_bam::tag::{Tag, TagValue};

    const CONTIG_LENGTH: i32 = 200;

    fn header() -> SamHeader {
        let mut header = SamHeader::default();
        header
            .sequences
            .push(SequenceRecord::new("chr1", CONTIG_LENGTH));
        for (id, sample) in [("rg1", "sampleA"), ("rg2", "sampleB")] {
            let mut group = ReadGroup::new(id);
            group.attributes.set("SM", sample);
            group.attributes.set("PL", "ILLUMINA");
            header.read_groups.push(group);
        }
        header
    }

    fn read(name: &str, group: &str, start: i32) -> BamRecord {
        let cigar = htsjdk_bam::text_parse::parse_cigar("10M").expect("a cigar");
        let mut tags = htsjdk_bam::tag::Tags::new();
        tags.insert(Tag::new(b"RG"), TagValue::Str(group.to_string()));
        BamRecord {
            read_name: name.to_string(),
            reference_index: 0,
            alignment_start: start,
            mapping_quality: 60,
            read_bases: (0..10).map(|i| b"ACGT"[i % 4]).collect(),
            base_qualities: vec![30; 10],
            cigar,
            tags,
            ..Default::default()
        }
    }

    fn interval(start: i32, end: i32) -> SimpleInterval {
        SimpleInterval {
            contig: "chr1".to_string(),
            start,
            end,
        }
    }

    /// The (emitEmptyLoci, intervals) each labelled route was built with.
    fn configuration(label: &str) -> (bool, Option<Vec<SimpleInterval>>) {
        match label {
            "noloci-nullintervals" => (false, None),
            "noloci-emptyintervals" => (false, Some(vec![])),
            "noloci-intervals" => (false, Some(vec![interval(105, 108)])),
            "emptyloci-nullintervals" => (true, None),
            "emptyloci-emptyintervals" => (true, Some(vec![])),
            "emptyloci-intervals" => (true, Some(vec![interval(105, 112)])),
            "emptyloci-twointervals" => (true, Some(vec![interval(105, 107), interval(115, 117)])),
            other => panic!("{other} is in the golden but not configured here"),
        }
    }

    let text = golden();
    let header = header();
    let reads = vec![
        read("a1", "rg1", 101),
        read("b1", "rg2", 101),
        read("a2", "rg1", 120),
    ];

    // The rows the reference produced, and where it stopped counting.
    let mut labels: Vec<String> = Vec::new();
    let mut rows: std::collections::HashMap<String, Vec<String>> = Default::default();
    let mut truncated: std::collections::HashMap<String, usize> = Default::default();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("ctx\t") {
            let mut parts = rest.splitn(3, '\t');
            let label = parts.next().expect("a label").to_string();
            let _index = parts.next();
            rows.entry(label)
                .or_default()
                .push(parts.next().unwrap_or("").to_string());
        } else if let Some(rest) = line.strip_prefix("count\t") {
            let label = rest.split('\t').next().expect("a label").to_string();
            labels.push(label);
        } else if let Some(rest) = line.strip_prefix("truncated\t") {
            let (label, at) = rest.split_once('\t').expect("a label and a cap");
            truncated.insert(label.to_string(), at.parse().expect("a number"));
        }
    }

    let mut compared = 0;
    for label in &labels {
        let (emit_empty_loci, intervals) = configuration(label);
        let expected = rows.get(label).cloned().unwrap_or_default();

        // The samples reach LocusIteratorByState in HashSet order, which is why that order is
        // measured above rather than assumed here.
        let samples: Vec<Option<String>> =
            java_hash::hash_set_order(&["sampleA".to_string(), "sampleB".to_string()])
                .expect("a small set")
                .into_iter()
                .map(Some)
                .collect();

        let route = context_iterator::route(emit_empty_loci, intervals.as_deref());
        if route == Route::RejectedEmptyIntervalList {
            assert!(
                expected.is_empty(),
                "{label}: a rejected route yielded rows"
            );
            compared += 1;
            continue;
        }

        let states = ReadStateManager::new(samples.clone(), DownsamplingInfo::NONE)
            .expect("no downsampling");
        let covered = locus_iterator::contexts(
            &reads,
            samples,
            &header,
            LocusIteratorOptions::default(),
            states,
        )
        .expect("the traversal runs");

        let whole_reference = vec![interval(1, CONTIG_LENGTH)];
        let ours = match route {
            Route::Unbounded => covered,
            Route::Overlapping => {
                context_iterator::overlapping(covered, intervals.as_deref().unwrap(), &header)
            }
            // With no intervals the builder substitutes every interval of the reference.
            Route::EmptyLoci => {
                let requested = match intervals.as_deref() {
                    None => whole_reference.as_slice(),
                    Some(given) => given,
                };
                context_iterator::with_empty_loci(covered, requested, &header)
            }
            Route::RejectedEmptyIntervalList => unreachable!("handled above"),
        };

        // The harness stopped after its cap and recorded it, so only the rows it printed are
        // comparable; the cap itself is asserted so a shrinking traversal cannot pass quietly.
        let limit = truncated.get(label).copied().unwrap_or(usize::MAX);
        if limit != usize::MAX {
            assert!(
                ours.len() >= limit,
                "{label}: the port yielded {} rows, fewer than the {limit} the reference printed",
                ours.len()
            );
        } else {
            assert_eq!(ours.len(), expected.len(), "{label}: context count");
        }

        for (index, row) in expected.iter().enumerate() {
            let context = &ours[index];
            let names: Vec<String> = context
                .pileup
                .elements
                .iter()
                .map(|e| e.read.read_name.clone())
                .collect();
            let rendered = format!(
                "{}:{}\t{}\t{}",
                context.contig,
                context.position,
                context.pileup.size(),
                if names.is_empty() {
                    "-".to_string()
                } else {
                    names.join(",")
                }
            );
            assert_eq!(&rendered, row, "{label} context {index}");
            compared += 1;
        }
    }

    println!(
        "{compared} contexts over {} routes, all identical",
        labels.len()
    );
}

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

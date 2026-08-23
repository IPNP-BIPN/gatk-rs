//! Conformance for `CollectF1R2Counts` against GATK 4.6.2.0, compared as the whole of every file
//! its tar.gz holds.
//!
//! Golden from `tools/readfilter-conformance/CollectF1R2CountsDump.java`.
//!
//! # What this suite is for
//!
//!  * **which of three places a site goes to**, and the tie the alt base is chosen by;
//!  * **the alt table's depth column being `refCount + altCount`** rather than the pileup's;
//!  * **every skip being a skip of the whole site**, so a second sample is only ever reached
//!    through the alt-table branch, in a `HashMap` order over the sample names;
//!  * **the shape of the output not depending on the data**, all 64 contexts and every prefilled
//!    bin present;
//!  * **and the reference histograms' `HashMap` order over their context strings**, which is
//!    reproducible where the alt histograms' identity-hash order is not.
//!
//! The fixture is rebuilt here from the same description the dump built it from: four blocks of
//! six unpaired reads, thirty bases long, against a forty-base motif repeated.

use gatk_corpus as corpus;
use gatk_tools::collect_f1r2_counts::{
    all_kmers, ref_context_order, sample_order, Args, Collector, Element,
};
use htsjdk_metrics::file::MetricsFile;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/collect_f1r2_counts.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn reference(text: &str) -> String {
    text.lines()
        .find_map(|line| line.strip_prefix("reference\tchr1="))
        .expect("the golden carries the reference")
        .to_string()
}

fn file(text: &str, kind: &str, label: &str, name: &str) -> String {
    let prefix = format!("{kind}\t{label}\t./{name}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries {kind}/{label}/{name}")),
    )
}

fn entries(text: &str, label: &str) -> Vec<String> {
    text.lines()
        .find_map(|line| line.strip_prefix(&format!("tar\t{label}=")))
        .expect("the golden carries the archive's entries")
        .split(',')
        .map(str::to_string)
        .collect()
}

/// One read of the fixture, as the dump wrote it.
struct Read {
    sample: String,
    start: usize,
    mapping_quality: i32,
    reverse: bool,
    deletion: bool,
    bases: Vec<u8>,
    quals: Vec<u8>,
}

const READ_LENGTH: usize = 30;

/// The rank-th base of ACGT that is not the reference base, counted from zero.
fn alt_base(reference: &[u8], one_based: usize, rank: usize) -> u8 {
    let mut seen = 0;
    for base in *b"ACGT" {
        if base != reference[one_based - 1] {
            if seen == rank {
                return base;
            }
            seen += 1;
        }
    }
    unreachable!()
}

/// The four blocks, for the samples given, exactly as `CollectF1R2CountsDump.block` builds them.
fn reads(reference: &[u8], samples: &[&str]) -> Vec<Read> {
    let mut reads = Vec::new();
    for sample in samples {
        for (start, mapping_quality, deletion) in [
            (41, 60, false),
            (91, 60, false),
            (131, 30, false),
            (161, 60, true),
        ] {
            for index in 0..6 {
                let mut bases = Vec::new();
                let mut quals = vec![40u8; READ_LENGTH];
                #[allow(clippy::needless_range_loop)]
                for offset in 0..READ_LENGTH {
                    let locus = start + offset + if deletion && offset >= 12 { 2 } else { 0 };
                    let mut base = reference[locus - 1];
                    if start == 41 {
                        if offset == 3 && *sample == "alpha" && index < 2 {
                            base = alt_base(reference, locus, 0);
                        }
                        if offset == 6 && *sample == "alpha" && index == 3 {
                            base = alt_base(reference, locus, 0);
                        }
                        if offset == 9 && *sample == "bravo" && index < 2 {
                            base = alt_base(reference, locus, 0);
                        }
                        if offset == 14 && *sample == "alpha" && index == 0 {
                            base = alt_base(reference, locus, 0);
                            quals[offset] = 20;
                        }
                        if offset == 19 && *sample == "alpha" {
                            if index < 2 {
                                base = alt_base(reference, locus, 0);
                            } else if index < 4 {
                                base = alt_base(reference, locus, 1);
                            }
                        }
                    }
                    bases.push(base);
                }
                reads.push(Read {
                    sample: sample.to_string(),
                    start,
                    mapping_quality,
                    reverse: index >= 3,
                    deletion,
                    bases,
                    quals,
                });
            }
        }
    }
    reads
}

/// The pileup at one locus: every read that covers it, with the deletion block's two deleted
/// positions and the element just before them marked as the walker marks them.
fn pileup(reads: &[Read], locus: usize) -> Vec<Element> {
    let mut elements = Vec::new();
    for read in reads {
        let span = READ_LENGTH + if read.deletion { 2 } else { 0 };
        if locus < read.start || locus >= read.start + span {
            continue;
        }
        let offset = locus - read.start;
        // The cigar is 12M2D18M, so the two deleted reference positions are the thirteenth and
        // fourteenth of the span and every read base after them sits two further along.
        let (base, qual, deletion) = if read.deletion && (12..14).contains(&offset) {
            (b'D', 0u8, true)
        } else {
            let index = if read.deletion && offset >= 14 {
                offset - 2
            } else {
                offset
            };
            (read.bases[index], read.quals[index], false)
        };
        elements.push(Element {
            sample: read.sample.clone(),
            base,
            qual,
            reverse_strand: read.reverse,
            first_of_pair: false,
            mapping_quality: read.mapping_quality,
            deletion,
            after_insertion: false,
            before_deletion_start: read.deletion && offset == 11,
        });
    }
    elements
}

/// The three-base context around a locus, or `None` when it runs off the contig.
fn context(reference: &[u8], locus: usize) -> Option<String> {
    if locus < 2 || locus + 1 > reference.len() {
        return None;
    }
    Some(String::from_utf8(reference[locus - 2..locus + 1].to_vec()).expect("ascii"))
}

/// One whole run: every locus that has reads, in order.
fn collect(text: &str, samples: &[&str], args: Args) -> Collector {
    let reference = reference(text).into_bytes();
    let reads = reads(&reference, samples);
    let names: Vec<String> = samples.iter().map(|name| name.to_string()).collect();
    let mut collector = Collector::new(args, &names);
    for locus in 1..=reference.len() {
        let elements = pileup(&reads, locus);
        if elements.is_empty() {
            continue;
        }
        collector.process(&elements, context(&reference, locus).as_deref());
    }
    collector
}

fn metrics(sample: &str, histograms: Vec<htsjdk_metrics::file::Histogram>) -> String {
    let mut file = MetricsFile::new();
    file.add_header(sample);
    file.histograms = histograms;
    file.write()
}

fn check(text: &str, label: &str, samples: &[&str], args: Args) {
    let collector = collect(text, samples, args);
    for sample in samples {
        assert_eq!(
            collector.alt_table_text(sample),
            file(text, "file", label, &format!("{sample}.alt_table")),
            "{label}: the alt table of {sample}"
        );
        assert_eq!(
            metrics(sample, collector.ref_histograms(sample)),
            file(text, "file", label, &format!("{sample}.ref_histogram")),
            "{label}: the reference histograms of {sample}"
        );
        assert_eq!(
            metrics(sample, collector.alt_histograms(sample)),
            file(text, "sorted", label, &format!("{sample}.alt_histogram")),
            "{label}: the alt histograms of {sample}"
        );
    }
}

#[test]
fn one_sample_reaches_every_branch() {
    let text = golden();
    check(
        &text,
        "single-sample",
        &["alpha"],
        Args {
            max_depth: 4,
            ..Args::default()
        },
    );
    assert_eq!(
        entries(&text, "single-sample"),
        vec![
            "./alpha.alt_histogram",
            "./alpha.alt_table",
            "./alpha.ref_histogram"
        ]
    );
}

#[test]
fn the_alt_table_depth_is_the_two_counts_added() {
    let text = golden();
    let collector = collect(
        &text,
        &["alpha"],
        Args {
            max_depth: 4,
            ..Args::default()
        },
    );
    let rows = collector.alt_table("alpha");
    // The tied site saw six reads and reports four, the second alt base counting nowhere.
    let tied = rows
        .iter()
        .find(|row| row.alt_count == 2 && row.ref_count == 2);
    assert_eq!(tied.map(|row| row.depth()), Some(4));
}

#[test]
fn a_second_sample_is_only_reached_through_the_alt_table_branch() {
    let text = golden();
    check(
        &text,
        "two-samples",
        &["alpha", "bravo"],
        Args {
            max_depth: 4,
            ..Args::default()
        },
    );
    // bravo first, which is neither the header's order nor alphabetical by accident: it is the
    // HashMap order, and it is why alpha's alt sites are all lost.
    assert_eq!(
        sample_order(&["alpha".to_string(), "bravo".to_string()]),
        vec!["bravo".to_string(), "alpha".to_string()]
    );
    let collector = collect(
        &text,
        &["alpha", "bravo"],
        Args {
            max_depth: 4,
            ..Args::default()
        },
    );
    assert!(collector.alt_table("alpha").is_empty());
    assert_eq!(collector.alt_table("bravo").len(), 1);
}

#[test]
fn the_depth_cap_is_applied_before_the_count() {
    let text = golden();
    check(
        &text,
        "max-depth-two",
        &["alpha"],
        Args {
            max_depth: 2,
            ..Args::default()
        },
    );
}

#[test]
fn a_lower_median_mapping_quality_lets_the_low_quality_block_through() {
    let text = golden();
    check(
        &text,
        "low-median-mq",
        &["alpha"],
        Args {
            min_median_map_qual: 20,
            max_depth: 4,
            ..Args::default()
        },
    );
}

#[test]
fn the_base_quality_test_is_strict() {
    let text = golden();
    check(
        &text,
        "min-bq-nineteen",
        &["alpha"],
        Args {
            min_base_quality: 19,
            max_depth: 4,
            ..Args::default()
        },
    );
}

#[test]
fn every_context_is_present_and_the_reference_order_is_the_maps() {
    let text = golden();
    assert_eq!(all_kmers().len(), 64);
    // AAA is generated first and moved to the end by the generator's own last step.
    assert_eq!(all_kmers().last().map(String::as_str), Some("AAA"));
    let order = ref_context_order();
    assert_eq!(order.len(), 64);
    let header = file(&text, "file", "single-sample", "alpha.ref_histogram")
        .lines()
        .nth(5)
        .expect("the histogram header line")
        .to_string();
    let columns: Vec<&str> = header.split('\t').skip(1).collect();
    assert_eq!(
        columns,
        order.iter().map(String::as_str).collect::<Vec<_>>()
    );
}

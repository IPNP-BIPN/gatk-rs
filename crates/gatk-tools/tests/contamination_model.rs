//! Conformance for the contamination model and the tool around it, against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/ContaminationModelDump.java`.
//!
//! # What this suite is for
//!
//!  * **the corpus is asserted before the model is**, so a port that built a different input fails
//!    on the input rather than on the arithmetic;
//!  * **the coverage filter is compared site by site**, because its low threshold is a ratio of the
//!    median and its high threshold a ratio of the mean, and getting one of the two wrong changes
//!    which sites the model ever sees;
//!  * **the segments are compared as intervals and counts**, which is where the decomposition
//!    reaches the answer;
//!  * **the minor allele fractions and the contamination are compared as raw bits**, because they
//!    come out of Brent optimisations whose last bits depend on every sum below them.
//!
//! The short table is the case that costs nothing and says the most: forty sites is below the
//! segmenter's window of fifty, so no changepoint candidate is ever computed and the contig is one
//! segment rather than none.

use gatk_corpus as corpus;
use gatk_engine::contamination_model::ContaminationModel;
use gatk_engine::contamination_segmenter::find_segments;
use gatk_engine::pileup_summary::PileupSummary;
use gatk_tools::calculate_contamination::{
    filter_sites_by_coverage, DEFAULT_HIGH_COVERAGE_RATIO_THRESHOLD,
    DEFAULT_LOW_COVERAGE_RATIO_THRESHOLD,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/contamination_model.txt.gz"),
    )
}

/// The dump's corpus, from the same formula: integer counts throughout, no rounding rule.
///
/// The genotype cycles hom ref, het, het, hom alt, het. With `loh` the second half of the contig
/// has its hets at a quarter rather than a half, which is the loss of heterozygosity the segmenter
/// is meant to find.
fn sites(contig: &str, count: usize, loh: bool) -> Vec<PileupSummary> {
    (0..count)
        .map(|i| {
            let position = 1000 + 100 * i as i32;
            let allele_frequency = 0.05 + (i % 19) as f64 * 0.05;
            let depth = 50 + (i % 7) as i32 * 5;
            let other = if i % 4 == 0 { 1 } else { 0 };
            let alt = match i % 5 {
                0 => 1 + (i % 2) as i32,
                3 => depth - other - (i % 2) as i32,
                _ => {
                    if loh && i >= count / 2 {
                        (depth - other) / 4
                    } else {
                        (depth - other) / 2
                    }
                }
            };
            let reference = depth - alt - other;
            PileupSummary::new(contig, position, reference, alt, other, allele_frequency)
        })
        .collect()
}

fn tumor() -> Vec<PileupSummary> {
    let mut all = sites("chr1", 200, true);
    all.extend(sites("chr2", 150, false));
    all
}

fn normal() -> Vec<PileupSummary> {
    sites("chr1", 200, false)
}

fn short_table() -> Vec<PileupSummary> {
    sites("chr3", 40, false)
}

fn filter(sites: &[PileupSummary]) -> Vec<PileupSummary> {
    filter_sites_by_coverage(
        sites,
        DEFAULT_LOW_COVERAGE_RATIO_THRESHOLD,
        DEFAULT_HIGH_COVERAGE_RATIO_THRESHOLD,
    )
}

/// `%016x` of the raw bits, as the dump prints a double.
fn bits(value: f64) -> String {
    format!("{:016x}", value.to_bits())
}

fn rows<'a>(text: &'a str, kind: &str, label: &str) -> Vec<&'a str> {
    let prefix = format!("{kind}\t{label}\t");
    text.lines()
        .filter_map(|line| line.strip_prefix(prefix.as_str()))
        .collect()
}

#[test]
fn the_corpus_is_the_reference_corpus() {
    let text = golden();
    for (label, built) in [
        ("tumor", tumor()),
        ("normal", normal()),
        ("short", short_table()),
    ] {
        let expected = rows(&text, "sites", label);
        assert_eq!(expected.len(), built.len(), "{label}: site count");
        for (index, site) in built.iter().enumerate() {
            let mine = format!(
                "{index}={},{},{},{},{},{}",
                site.contig,
                site.position,
                site.ref_count,
                site.alt_count,
                site.other_alt_count,
                bits(site.allele_frequency)
            );
            assert_eq!(mine, expected[index], "{label} site {index}");
        }
    }
}

#[test]
fn the_coverage_filter_keeps_the_same_sites() {
    let text = golden();
    for (label, built) in [
        ("tumor", tumor()),
        ("normal", normal()),
        ("short", short_table()),
    ] {
        let expected = rows(&text, "coverage", label);
        assert_eq!(expected.len(), 1, "{label}: one coverage row");
        let kept = filter(&built);
        let mine = if kept.is_empty() {
            "(none)".to_string()
        } else {
            kept.iter()
                .map(|site| format!("{}:{}", site.contig, site.position))
                .collect::<Vec<String>>()
                .join(",")
        };
        assert_eq!(mine, expected[0], "{label}: the sites that survive");
    }
}

#[test]
fn every_segment_matches_the_golden() {
    let text = golden();
    for (label, built) in [
        ("tumor", tumor()),
        ("normal", normal()),
        ("short", short_table()),
    ] {
        let expected = rows(&text, "segments", label);
        let segments = find_segments(&filter(&built));
        assert_eq!(segments.len(), expected.len(), "{label}: segment count");
        for (index, segment) in segments.iter().enumerate() {
            let mine = format!(
                "{index}={},{},{},{}",
                segment[0].contig,
                segment[0].position,
                segment[segment.len() - 1].position,
                segment.len()
            );
            assert_eq!(mine, expected[index], "{label} segment {index}");
        }
    }
}

#[test]
fn every_minor_allele_fraction_and_contamination_matches_the_golden() {
    let text = golden();
    let filtered_tumor = filter(&tumor());
    let filtered_normal = filter(&normal());
    let filtered_short = filter(&short_table());

    let tumor_model = ContaminationModel::new(&filtered_tumor);
    let normal_model = ContaminationModel::new(&filtered_normal);
    let short_model = ContaminationModel::new(&filtered_short);

    for (label, records) in [
        ("tumor-only", tumor_model.segmentation_records()),
        ("matched-normal", normal_model.segmentation_records()),
        ("short", short_model.segmentation_records()),
    ] {
        let expected = rows(&text, "maf", label);
        assert_eq!(records.len(), expected.len(), "{label}: record count");
        for (index, record) in records.iter().enumerate() {
            let mine = format!(
                "{index}={},{},{},{}",
                record.contig,
                record.start,
                record.end,
                bits(record.minor_allele_fraction)
            );
            assert_eq!(mine, expected[index], "{label} record {index}");
        }
    }

    // The tumour genotypes itself, the normal genotypes the tumour, and the short table has only
    // itself. All three answers are compared as raw bits.
    for (label, (estimate, error)) in [
        (
            "tumor-only",
            tumor_model.calculate_contamination_from_homs(&filtered_tumor),
        ),
        (
            "matched-normal",
            normal_model.calculate_contamination_from_homs(&filtered_tumor),
        ),
        (
            "short",
            short_model.calculate_contamination_from_homs(&filtered_short),
        ),
    ] {
        let expected = rows(&text, "contamination", label);
        assert_eq!(expected.len(), 1, "{label}: one contamination row");
        assert_eq!(
            format!("{},{}", bits(estimate), bits(error)),
            expected[0],
            "{label}: the estimate and its standard error"
        );
    }
}

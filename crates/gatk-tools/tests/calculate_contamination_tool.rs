//! Conformance for `CalculateContamination` against GATK 4.6.2.0, compared as the two files it
//! writes from the pileup tables the golden also carries.
//!
//! Golden from `tools/readfilter-conformance/CalculateContaminationDump.java`. The model
//! underneath has a suite of its own (`contamination-model`); this one is the tool around it.
//!
//! # What this suite is for
//!
//!  * **the two ratio arguments being accepted and ignored**;
//!  * **the coverage filter still cutting at three times the mean**;
//!  * **a site at or under `MIN_COVERAGE` moving neither statistic**;
//!  * **a matched normal saying which sites are homozygous while the tumour's counts are read**;
//!  * **a normal with no homozygous-alternate site falling through to the tumour-only answer**;
//!  * **the segmentation being optional and always the tumour's own**;
//!  * **and the sample name coming from the tumour table's metadata line.**

use gatk_corpus as corpus;
use gatk_engine::pileup_summary::PileupSummary;
use gatk_tools::calculate_contamination::{
    filter_sites_by_coverage, run_from_command_line, DEFAULT_HIGH_COVERAGE_RATIO_THRESHOLD,
    DEFAULT_LOW_COVERAGE_RATIO_THRESHOLD, MIN_COVERAGE,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/calculate_contamination_tool.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn field(text: &str, kind: &str, label: &str) -> Option<String> {
    let prefix = format!("{kind}\t{label}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
}

/// One of the pileup tables the dump printed, parsed back.
fn table(text: &str, label: &str) -> Vec<PileupSummary> {
    let prefix = format!("table\t{label}=");
    let body = text
        .lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
        .unwrap_or_else(|| panic!("table/{label}"));
    body.lines()
        .filter(|line| !line.starts_with('#') && !line.starts_with("contig") && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            PileupSummary::new(
                columns[0],
                columns[1].parse().expect("a position"),
                columns[2].parse().expect("a ref count"),
                columns[3].parse().expect("an alt count"),
                columns[4].parse().expect("an other count"),
                columns[5].parse().expect("a frequency"),
            )
        })
        .collect()
}

/// The contamination table a case wrote, as its sample name and two numbers.
fn written(text: &str, label: &str) -> (String, f64, f64) {
    let body = field(text, "contamination", label).unwrap_or_else(|| panic!("{label}"));
    let mut lines = body.lines();
    lines.next().expect("a header");
    let row: Vec<&str> = lines.next().expect("a row").split('\t').collect();
    (
        row[0].to_string(),
        row[1].parse().expect("a contamination"),
        row[2].parse().expect("an error"),
    )
}

/// The port reproduces both numbers, on every case whose inputs the golden carries.
#[test]
fn the_port_recomputes_both_numbers() {
    let text = golden();
    let tumour = table(&text, "tumour");
    let normal = table(&text, "normal");
    let hom_ref = table(&text, "normal-all-hom-ref");
    for (label, sites, matched) in [
        ("plain", &tumour, None),
        (
            "uncovered-sites",
            &table(&text, "with-uncovered-sites"),
            None,
        ),
        ("a-deep-tail", &table(&text, "with-a-deep-tail"), None),
        ("a-shallow-tail", &table(&text, "with-a-shallow-tail"), None),
        ("matched-normal", &tumour, Some(&normal)),
        ("matched-normal-homref", &tumour, Some(&hom_ref)),
    ] {
        let (_, contamination, error) = written(&text, label);
        let ours = run_from_command_line(
            sites,
            matched.map(|m| m.as_slice()),
            false,
            DEFAULT_LOW_COVERAGE_RATIO_THRESHOLD,
            DEFAULT_HIGH_COVERAGE_RATIO_THRESHOLD,
        );
        assert!(
            (ours.contamination - contamination).abs() < 1e-12,
            "{label}: {} against {contamination}",
            ours.contamination
        );
        assert!((ours.error - error).abs() < 1e-12, "{label}");
    }
}

/// The two ratio arguments are accepted and ignored.
#[test]
fn the_ratio_arguments_do_nothing() {
    let text = golden();
    let plain = written(&text, "plain");
    for label in [
        "low-ratio-zero",
        "high-ratio-one",
        "low-ratio-one",
        "low-ratio-ten",
    ] {
        // `low-ratio-zero` and `a-deep-tail` share a table, so compare each against its own base.
        let base = if label == "low-ratio-zero" {
            written(&text, "a-deep-tail")
        } else {
            plain.clone()
        };
        assert_eq!(written(&text, label), base, "{label}");
    }
    // Which is what the port's command-line entry point does with the same numbers: a ratio of ten
    // would drop every site, and does not.
    let tumour = table(&text, "tumour");
    let ignored = run_from_command_line(&tumour, None, false, 10.0, 0.001);
    let (_, contamination, error) = plain;
    assert!((ignored.contamination - contamination).abs() < 1e-12);
    assert!((ignored.error - error).abs() < 1e-12);
    // The filter itself does read thresholds when it is called directly, which is the difference
    // between the two entry points.
    let sites = filter_sites_by_coverage(&tumour, 10.0, 3.0);
    assert!(sites.is_empty());
}

/// The coverage filter still cuts, at three times the mean.
#[test]
fn the_ceiling_still_cuts() {
    let text = golden();
    let deep = written(&text, "a-deep-tail");
    let shallow = written(&text, "a-shallow-tail");
    let plain = written(&text, "plain");
    // The same six homozygous-alternate sites: dropped at depth four hundred, kept at sixty.
    assert_eq!(deep, plain);
    assert_ne!(shallow, plain);
    // And the port drops them for the same reason.
    let with_deep = table(&text, "with-a-deep-tail");
    let kept = filter_sites_by_coverage(
        &with_deep,
        DEFAULT_LOW_COVERAGE_RATIO_THRESHOLD,
        DEFAULT_HIGH_COVERAGE_RATIO_THRESHOLD,
    );
    assert!(kept.iter().all(|site| site.total_count < 400));
    assert_eq!(kept.len(), table(&text, "tumour").len());
}

/// A site at or under `MIN_COVERAGE` moves neither statistic.
#[test]
fn an_uncovered_site_moves_nothing() {
    let text = golden();
    assert_eq!(written(&text, "uncovered-sites"), written(&text, "plain"));
    let with_uncovered = table(&text, "with-uncovered-sites");
    assert!(with_uncovered
        .iter()
        .any(|site| site.total_count <= MIN_COVERAGE));
    let kept = filter_sites_by_coverage(
        &with_uncovered,
        DEFAULT_LOW_COVERAGE_RATIO_THRESHOLD,
        DEFAULT_HIGH_COVERAGE_RATIO_THRESHOLD,
    );
    assert!(kept.iter().all(|site| site.total_count > MIN_COVERAGE));
}

/// A matched normal says which sites are homozygous; the tumour's own counts are read there.
#[test]
fn the_normal_genotypes_and_the_tumour_is_counted() {
    let text = golden();
    let (_, alone, _) = written(&text, "plain");
    let (_, matched, _) = written(&text, "matched-normal");
    assert!((alone - 0.04388459975619667).abs() < 1e-12);
    assert!((matched - 0.9150511103278111).abs() < 1e-12);
    // A normal with no homozygous-alternate site falls through to the strategies that read the
    // tumour and not the model, so the answer returns to the tumour-only one.
    assert_eq!(
        written(&text, "matched-normal-homref"),
        written(&text, "plain")
    );
}

/// The segmentation is optional, and always the tumour's own.
#[test]
fn the_segmentation_is_optional_and_the_tumours() {
    let text = golden();
    assert_eq!(field(&text, "segments", "plain").as_deref(), Some("absent"));
    assert_eq!(
        field(&text, "segments", "matched-normal").as_deref(),
        Some("absent")
    );
    // Asking for it changes neither number in the contamination table.
    assert_eq!(written(&text, "with-segmentation"), written(&text, "plain"));
    assert_eq!(
        written(&text, "matched-normal-segmented"),
        written(&text, "matched-normal")
    );
    // And the table it writes is the same one either way, because it is the tumour's model.
    let alone = field(&text, "segments", "with-segmentation").expect("a segmentation");
    let matched = field(&text, "segments", "matched-normal-segmented").expect("a segmentation");
    assert_eq!(alone, matched);
    // Which is what the port writes from the tumour's own sites.
    let tumour = table(&text, "tumour");
    let ours = run_from_command_line(
        &tumour,
        Some(&table(&text, "normal")),
        true,
        DEFAULT_LOW_COVERAGE_RATIO_THRESHOLD,
        DEFAULT_HIGH_COVERAGE_RATIO_THRESHOLD,
    )
    .segmentation
    .expect("a segmentation");
    let theirs: Vec<&str> = matched.lines().skip(2).collect();
    assert_eq!(ours.len(), theirs.len());
    for (record, row) in ours.iter().zip(theirs) {
        let columns: Vec<&str> = row.split('\t').collect();
        assert_eq!(record.contig, columns[0]);
        assert_eq!(record.start.to_string(), columns[1]);
        assert_eq!(record.end.to_string(), columns[2]);
        let fraction: f64 = columns[3].parse().expect("a fraction");
        assert!((record.minor_allele_fraction - fraction).abs() < 1e-12);
    }
}

/// The sample name is the tumour table's metadata line, and the tool answers `SUCCESS`.
#[test]
fn the_sample_name_is_the_tumours() {
    let text = golden();
    for label in ["plain", "matched-normal"] {
        assert_eq!(written(&text, label).0, "tumour", "{label}");
        assert_eq!(
            field(&text, "returned", label).as_deref(),
            Some("SUCCESS"),
            "{label}"
        );
    }
    // The normal's own metadata line says `normal`, and never reaches either output.
    assert!(text.contains("SAMPLE=normal"));
}

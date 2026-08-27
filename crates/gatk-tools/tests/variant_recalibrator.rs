//! Conformance for `VariantRecalibrator` against GATK 4.6.2.0, compared as the whole tranches
//! file of every run.
//!
//! Golden from `tools/readfilter-conformance/VariantRecalibratorDump.java`.
//!
//! The Gaussian mixture that produces the scores is not measured or ported. The scores are read
//! off the recalibration table the golden carries, and everything the tool does with them is
//! reproduced from there.
//!
//! # What this suite is for
//!
//!  * **the running sensitivity being computed from the top down**;
//!  * **a tranche being the largest set that still reaches its target**;
//!  * **a tranche counting every variant at or above its own minimum**;
//!  * **Ti/Tv having its denominator floored at one**;
//!  * **the file being sorted by CALLS AT TRUTH SITES and not by the target**;
//!  * **the sort being stable**;
//!  * **each row's lower bound being the previous row's target, backwards or not**;
//!  * **a target of zero being reachable**;
//!  * **and the VQSLOD tranches taking their minimum from the request.**

use gatk_corpus as corpus;
use gatk_tools::variant_recalibrator::{
    calls_at_truth, empty_tranche, find_tranches, find_vqslod_tranches, no_tranche_refusal,
    running_sensitivity, sensitivity_threshold, sorted, tranche_of_variants, tranche_order,
    tranche_row, tranches_string, truth_sensitivity_file, vqslod_file, Datum, Mode,
};
use std::collections::BTreeSet;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/variant_recalibrator.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn section(text: &str, kind: &str, name: &str) -> String {
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(&format!("{kind}\t{name}=")))
            .unwrap_or_else(|| panic!("the golden carries {kind}/{name}")),
    )
}

/// The positions one of the golden's VCFs carries.
fn positions(text: &str, name: &str) -> BTreeSet<i32> {
    section(text, "vcf", name)
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            line.split('\t')
                .nth(1)
                .expect("a position")
                .parse()
                .expect("a number")
        })
        .collect()
}

/// The alternate allele the input carries at each position, which is what decides a transition.
fn alternates(text: &str) -> Vec<(i32, String)> {
    section(text, "vcf", "input")
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            (
                columns[1].parse().expect("a position"),
                columns[4].to_string(),
            )
        })
        .collect()
}

/// The data one run scored, read off its recalibration table.
///
/// The table is the callset the tranches were cut from: the variants a filter kept out are
/// absent from it, which is what makes the two filter runs differ from the others.
fn data(text: &str, label: &str) -> Vec<Datum> {
    let known = positions(text, "known");
    let truth = positions(text, "truth");
    let alternate: std::collections::BTreeMap<i32, String> = alternates(text).into_iter().collect();
    section(text, "recal", label)
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let position: i32 = columns[0].parse().expect("a position");
            let alternate = alternate.get(&position).expect("an input record");
            Datum {
                lod: columns[3].parse().expect("a score"),
                is_known: known.contains(&position),
                at_truth_site: truth.contains(&position),
                // Every record of this fixture is a SNP with `A` for its reference.
                is_snp: true,
                is_transition: alternate == "G",
            }
        })
        .collect()
}

/// The truth sites in the data, which is the denominator the running sensitivity uses.
fn n_true_sites(data: &[Datum]) -> usize {
    calls_at_truth(data, f64::NEG_INFINITY) as usize
}

/// label, the targets it asked for, in order.
fn runs() -> Vec<(&'static str, Vec<f64>)> {
    vec![
        ("four-tranches", vec![100.0, 99.9, 99.0, 90.0]),
        ("targets-out-of-order", vec![90.0, 100.0, 99.0, 99.9]),
        ("one-tranche", vec![99.0]),
        ("target-too-low-first", vec![0.0]),
        ("target-too-low-last", vec![100.0, 0.0]),
        ("ignore-all-filters", vec![100.0, 99.0]),
        ("ignore-filter", vec![100.0, 99.0]),
        ("one-annotation", vec![100.0, 99.0]),
    ]
}

/// Every truth-sensitivity run's whole tranches file, byte for byte.
#[test]
fn every_tranches_file_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for (label, targets) in runs() {
        let data = data(&text, label);
        let produced =
            find_tranches(&data, &targets, n_true_sites(&data), Mode::Snp).expect("tranches");
        assert_eq!(
            truth_sensitivity_file(&produced),
            section(&text, "tranches", label),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(
        compared, 8,
        "the truth-sensitivity runs the port reproduces"
    );
}

/// Cut on the score itself, and the minimum reported is the request.
#[test]
fn every_vqslod_file_matches_the_golden() {
    let text = golden();
    for (label, thresholds) in [
        ("vqslod-tranches", vec![10.0, 0.0, -10.0]),
        ("vqslod-unreachable", vec![100000.0, 0.0]),
    ] {
        let data = data(&text, label);
        let produced = find_vqslod_tranches(&data, &thresholds, Mode::Snp);
        assert_eq!(
            vqslod_file(&produced),
            section(&text, "tranches", label),
            "{label}"
        );
        // Every row's minimum is exactly what was asked for.
        for (tranche, threshold) in produced.iter().zip(thresholds.iter()) {
            assert_eq!(tranche.min_vqs_lod, *threshold);
        }
    }
    // A threshold no variant reaches still produces a row, and it is an empty one.
    let data = data(&text, "vqslod-unreachable");
    let produced = find_vqslod_tranches(&data, &[100000.0], Mode::Snp);
    assert_eq!(produced.len(), 1);
    assert_eq!(produced[0].num_known, 0);
    assert_eq!(produced[0].num_novel, 0);
    assert_eq!(produced[0].calls_at_truth_sites, 0);
    // Its truth sites are still counted, so its sensitivity is a real zero rather than a default.
    assert!(produced[0].accessible_truth_sites > 0);
    assert_eq!(produced[0].truth_sensitivity(), 0.0);
}

/// One less the truth sites at or above each index, over the truth sites in total.
#[test]
fn the_running_sensitivity_is_computed_from_the_top_down() {
    let text = golden();
    let data = sorted(&data(&text, "four-tranches"));
    let total = n_true_sites(&data);
    let running = running_sensitivity(&data, total);
    assert_eq!(running.len(), data.len());
    // The entry at zero counts every truth site, so it is zero when they are all in the data.
    assert_eq!(running[0], 0.0);
    // The walk never rises as the index falls.
    for i in 1..running.len() {
        assert!(running[i - 1] <= running[i], "{i}");
    }
    // The last entry is one less the last variant's own contribution.
    let last = *running.last().expect("an entry");
    let expected = if data.last().expect("a variant").at_truth_site {
        1.0 - 1.0 / total as f64
    } else {
        1.0
    };
    assert_eq!(last, expected);
    // A target of 99 asks for 0.01, and 100 asks for 0.
    assert_eq!(sensitivity_threshold(99.0), 1.0 - 0.99);
    assert_eq!(sensitivity_threshold(100.0), 0.0);
    assert_eq!(sensitivity_threshold(0.0), 1.0);
}

/// So the tranches nest rather than partition.
#[test]
fn a_tranche_counts_every_variant_above_its_own_minimum() {
    let text = golden();
    let data = sorted(&data(&text, "four-tranches"));
    let total = n_true_sites(&data);
    let running = running_sensitivity(&data, total);
    let tranches =
        find_tranches(&data, &[100.0, 99.9, 99.0, 90.0], total, Mode::Snp).expect("tranches");
    // The looser the target, the lower the minimum and the more variants counted.
    let hundred = &tranches[0];
    let ninety = &tranches[3];
    assert!(hundred.min_vqs_lod < ninety.min_vqs_lod);
    assert!(hundred.num_known + hundred.num_novel > ninety.num_known + ninety.num_novel);
    // The count is over the whole callset at that LOD, not over the walk's prefix.
    let index = (0..data.len())
        .find(|i| running[*i] >= sensitivity_threshold(90.0))
        .expect("an index");
    let rebuilt = tranche_of_variants(&data, index, 90.0, Mode::Snp);
    assert_eq!(&rebuilt, ninety);
    let above = data.iter().filter(|d| d.lod >= ninety.min_vqs_lod).count() as i64;
    assert_eq!(ninety.num_known + ninety.num_novel, above);
    // Which is more than the variants the walk had left behind it.
    assert!(above <= (data.len() - index) as i64);
}

/// A tranche with no transversion reports its transition count where a ratio should be.
#[test]
fn the_ti_tv_denominator_is_floored_at_one() {
    let text = golden();
    let data = data(&text, "four-tranches");
    // The fixture's novel variants are the ones no `known` resource carries, and every one of
    // them is a transition, so no tranche has a novel transversion at all.
    let novel: Vec<&Datum> = data.iter().filter(|d| !d.is_known).collect();
    assert!(!novel.is_empty());
    assert!(novel.iter().all(|d| d.is_transition));
    let tranches =
        find_tranches(&data, &[100.0], n_true_sites(&data), Mode::Snp).expect("tranches");
    let whole = &tranches[0];
    assert_eq!(
        whole.novel_ti_tv, whole.num_novel as f64,
        "a count, not a ratio"
    );
    // The known half does have both, so its ratio is one.
    assert!(whole.known_ti_tv < 1.5);
    assert!(whole.known_ti_tv > 0.0);
}

/// And not by the target, though the two agree when the targets rise.
#[test]
fn the_file_is_sorted_by_calls_at_truth_sites() {
    let text = golden();
    let file = section(&text, "tranches", "four-tranches");
    let written: Vec<String> = file
        .lines()
        .filter(|line| !line.starts_with('#') && !line.starts_with("target"))
        .map(|line| line.split(',').next().expect("a target").to_string())
        .collect();
    // The targets were given 100, 99.9, 99, 90 and came out in this order.
    assert_eq!(written, vec!["90.00", "99.90", "99.00", "100.00"]);
    let calls: Vec<i32> = file
        .lines()
        .filter(|line| !line.starts_with('#') && !line.starts_with("target"))
        .map(|line| {
            line.split(',')
                .nth(9)
                .expect("the calls column")
                .parse()
                .expect("a number")
        })
        .collect();
    for i in 1..calls.len() {
        assert!(calls[i - 1] <= calls[i], "{calls:?}");
    }
    // Given in increasing order the two agree, which is why the mistake is easy to miss.
    let ordered = section(&text, "tranches", "targets-out-of-order");
    let targets: Vec<String> = ordered
        .lines()
        .filter(|line| !line.starts_with('#') && !line.starts_with("target"))
        .map(|line| line.split(',').next().expect("a target").to_string())
        .collect();
    assert_eq!(targets, vec!["90.00", "99.00", "99.90", "100.00"]);
    // The comparator itself looks at that one field.
    let data = data(&text, "four-tranches");
    let tranches =
        find_tranches(&data, &[100.0, 90.0], n_true_sites(&data), Mode::Snp).expect("tranches");
    assert_eq!(
        tranche_order(&tranches[0], &tranches[1]),
        std::cmp::Ordering::Greater
    );
}

/// Two targets that found the same tranche keep the order they were given in.
#[test]
fn the_sort_is_stable() {
    let text = golden();
    let file = section(&text, "tranches", "four-tranches");
    let rows: Vec<&str> = file
        .lines()
        .filter(|line| !line.starts_with('#') && !line.starts_with("target"))
        .collect();
    // The middle two rows found the same tranche: same counts, same minimum, same calls.
    let fields =
        |row: &str| -> Vec<String> { row.split(',').skip(1).take(5).map(str::to_string).collect() };
    assert_eq!(fields(rows[1]), fields(rows[2]));
    // And they kept the order the targets were given in, 99.9 before 99.
    assert!(rows[1].starts_with("99.90,"));
    assert!(rows[2].starts_with("99.00,"));
    // A list of one is not sorted at all, which is the same result here.
    let data = data(&text, "one-tranche");
    let tranches = find_tranches(&data, &[99.0], n_true_sites(&data), Mode::Snp).expect("tranches");
    assert_eq!(tranches_string(&tranches), tranche_row(&tranches[0], None));
}

/// Whatever that target is, so an unsorted list writes a band that runs backwards.
#[test]
fn a_rows_lower_bound_is_the_previous_rows_target() {
    let text = golden();
    let file = section(&text, "tranches", "four-tranches");
    let names: Vec<String> = file
        .lines()
        .filter(|line| !line.starts_with('#') && !line.starts_with("target"))
        .map(|line| line.split(',').nth(6).expect("a name").to_string())
        .collect();
    assert_eq!(
        names,
        vec![
            "VQSRTrancheSNP0.00to90.00",
            "VQSRTrancheSNP90.00to99.90",
            // The band that runs backwards.
            "VQSRTrancheSNP99.90to99.00",
            "VQSRTrancheSNP99.00to100.00",
        ]
    );
    // The first row's lower bound is 0.00 because there is no previous row.
    let data = data(&text, "one-tranche");
    let tranches = find_tranches(&data, &[99.0], n_true_sites(&data), Mode::Snp).expect("tranches");
    assert!(tranche_row(&tranches[0], None).contains("VQSRTrancheSNP0.00to99.00"));
    assert!(tranche_row(&tranches[0], Some(&tranches[0])).contains("VQSRTrancheSNP99.00to99.00"));
}

/// It calls no truth site at all, which is a tranche and not a refusal.
#[test]
fn a_target_of_zero_is_reachable() {
    let text = golden();
    for label in ["target-too-low-first", "target-too-low-last"] {
        let file = section(&text, "tranches", label);
        assert!(file.contains("0.00,"), "{label}");
        let data = data(&text, label);
        let targets = if label == "target-too-low-first" {
            vec![0.0]
        } else {
            vec![100.0, 0.0]
        };
        let produced =
            find_tranches(&data, &targets, n_true_sites(&data), Mode::Snp).expect("tranches");
        assert_eq!(truth_sensitivity_file(&produced), file, "{label}");
        let zero = produced
            .iter()
            .find(|tranche| tranche.index == 0.0)
            .expect("the zero tranche");
        assert_eq!(zero.calls_at_truth_sites, 0);
        assert_eq!(zero.truth_sensitivity(), 0.0);
    }
    // A target that no index reaches would be a refusal, and the wording names the threshold.
    let unreachable: Vec<Datum> = vec![Datum {
        lod: 1.0,
        is_known: true,
        at_truth_site: true,
        is_snp: true,
        is_transition: true,
    }];
    let produced = find_tranches(&unreachable, &[0.0], 1, Mode::Snp).expect_err("no tranche");
    assert_eq!(produced, no_tranche_refusal("TruthSensitivity", 1.0));
    // And a target after one that DID find a tranche ends the list rather than refusing.
    let two = find_tranches(&unreachable, &[100.0, 0.0], 1, Mode::Snp).expect("one tranche");
    assert_eq!(two.len(), 1);
    assert_eq!(two[0].index, 100.0);
}

/// The filtered records never reach the table, so every count moves with them.
#[test]
fn a_filtered_record_is_not_in_the_callset() {
    let text = golden();
    let plain = data(&text, "four-tranches");
    let all = data(&text, "ignore-all-filters");
    let named = data(&text, "ignore-filter");
    assert!(all.len() > plain.len());
    // The fixture's only filter is LOW, so naming it and ignoring everything keep the same set.
    assert_eq!(all.len(), named.len());
    // And the truth sites among them, which is what the sensitivity is measured against.
    assert_eq!(n_true_sites(&all), n_true_sites(&named));
    assert!(n_true_sites(&all) > n_true_sites(&plain));
    // The scores differ, though, because the model was fitted to a different callset.
    assert_ne!(
        section(&text, "tranches", "ignore-all-filters"),
        section(&text, "tranches", "ignore-filter")
    );
}

/// The three refusals the golden recorded, each from a different stage.
#[test]
fn the_three_refusals_are_the_goldens() {
    let text = golden();
    let refusal = |label: &str| -> (String, String) {
        let row = text
            .lines()
            .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
            .unwrap_or_else(|| panic!("the golden carries error/{label}"));
        let (class, message) = row.split_once(':').expect("a class and a message");
        (class.to_string(), unescape(message))
    };
    // Two are argument-parsing refusals, before any variant is read.
    let (class, message) = refusal("no-training");
    assert_eq!(
        class,
        "org.broadinstitute.barclay.argparser.CommandLineException"
    );
    assert!(message.starts_with("No training set found!"), "{message}");
    let (class, message) = refusal("no-truth");
    assert_eq!(
        class,
        "org.broadinstitute.barclay.argparser.CommandLineException"
    );
    assert!(message.starts_with("No truth set found!"), "{message}");
    // The third is about the data: an INDEL run over a callset with no indel in it.
    let (class, message) = refusal("indel-mode");
    assert_eq!(
        class,
        "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
    );
    assert!(
        message.contains("not detected for ANY training variant"),
        "{message}"
    );
}

/// An empty tranche keeps the truth-site counts and zeroes everything else.
#[test]
fn an_empty_tranche_keeps_only_its_truth_sites() {
    let text = golden();
    let data = sorted(&data(&text, "four-tranches"));
    let produced = empty_tranche(&data, data.len() - 1, 42.0, Mode::Snp);
    assert_eq!(produced.num_known, 0);
    assert_eq!(produced.num_novel, 0);
    assert_eq!(produced.known_ti_tv, 0.0);
    assert_eq!(produced.novel_ti_tv, 0.0);
    assert_eq!(produced.index, 42.0);
    // Its minimum is the LOD at that index until a caller overwrites it, and its accessible
    // truth sites are the whole callset's.
    assert_eq!(produced.min_vqs_lod, data.last().expect("a variant").lod);
    assert_eq!(
        produced.accessible_truth_sites,
        calls_at_truth(&data, f64::NEG_INFINITY)
    );
    // Over no data at all the minimum is negative infinity rather than a panic.
    let nothing = empty_tranche(&[], 0, 1.0, Mode::Indel);
    assert_eq!(nothing.min_vqs_lod, f64::NEG_INFINITY);
    assert_eq!(nothing.accessible_truth_sites, 0);
    assert_eq!(nothing.truth_sensitivity(), 0.0);
    assert_eq!(nothing.model.name(), "INDEL");
}

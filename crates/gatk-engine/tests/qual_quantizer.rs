//! Conformance for `QualQuantizer` and `QuantizationInfo` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/QualQuantizerDump.java`. Six histograms, each shaped
//! for a different merge sequence, at seven level counts and three minimum interesting qualities,
//! with every interval of every final forest and every double as its raw bits.
//!
//! # What this suite is for
//!
//!  * **every leaf carries a fixed quality**, so its error rate is the theoretical one its Phred
//!    score declares and not the one its counts imply, and a quantized quality of zero is reachable;
//!  * **the leaf error count saturates at an `int`** on the way from a `long` count;
//!  * **the merge search keeps the first minimum**, which decides everything on an empty histogram
//!    where every penalty is zero;
//!  * **a leaf at or below the minimum interesting quality is free to merge**;
//!  * **`errorProbToQual` widens before it narrows**, which is why an error rate of zero is 93 and
//!    not 1;
//!  * **the level count counts changes**, not distinct values.

use gatk_corpus as corpus;
use gatk_engine::qual_quantizer::{
    calculate_quantization_levels, error_prob_to_qual, QualQuantizer, QuantizationInfo,
    QuantizerError, MAX_SAM_QUAL_SCORE, MIN_USABLE_Q_SCORE,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/qual_quantizer.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter_map(|line| {
            line.strip_prefix(kind)
                .and_then(|rest| rest.strip_prefix('\t'))
        })
        .map(|rest| rest.split('\t').collect())
        .collect()
}

fn bits(value: f64) -> String {
    format!("{:x}", value.to_bits())
}

/// The same six histograms the dump built, by the same names.
fn histogram(name: &str) -> Vec<i64> {
    match name {
        "flat" => vec![1000; 94],
        "empty" => vec![0; 94],
        "illumina" => {
            let mut counts = vec![0i64; 94];
            counts[2] = 5_000_000;
            for (q, count) in counts.iter_mut().enumerate().take(36).skip(25) {
                *count = 1_000_000 * (11 - (q as i64 - 30).abs());
            }
            counts
        }
        "single" => {
            let mut counts = vec![0i64; 94];
            counts[40] = 1234;
            counts
        }
        "short" => vec![100; 5],
        "overflowing" => {
            let mut counts = vec![0i64; 94];
            counts[0] = 3_000_000_000;
            counts[1] = 3_000_000_000;
            counts
        }
        other => panic!("{other} is in the golden but not configured here"),
    }
}

/// A label like `illumina@16,min40` back into its three arguments.
fn arguments(label: &str) -> (Vec<i64>, i32, i32) {
    let (name, rest) = label.split_once('@').expect("a label is name@levels");
    match rest.split_once(",min") {
        Some((levels, min)) => (
            histogram(name),
            levels.parse().unwrap(),
            min.parse().unwrap(),
        ),
        None => (histogram(name), rest.parse().unwrap(), MIN_USABLE_Q_SCORE),
    }
}

/// `QuantizationInfo` always holds a histogram at the full quality range.
fn padded(name: &str) -> Vec<i64> {
    let mut out = vec![0i64; (MAX_SAM_QUAL_SCORE + 1) as usize];
    for (index, count) in histogram(name).into_iter().enumerate() {
        if index < out.len() {
            out[index] = count;
        }
    }
    out
}

#[test]
fn the_constants_are_the_references() {
    let text = golden();
    let constant = |name: &str| {
        rows(&text, "const")
            .into_iter()
            .find(|row| row[0] == name)
            .unwrap_or_else(|| panic!("no constant {name}"))[1]
            .to_string()
    };
    assert_eq!(
        constant("MIN_USABLE_Q_SCORE"),
        MIN_USABLE_Q_SCORE.to_string()
    );
    assert_eq!(
        constant("MAX_SAM_QUAL_SCORE"),
        MAX_SAM_QUAL_SCORE.to_string()
    );
}

/// The rounding and the two clamps, which decide every merged interval's quality.
#[test]
fn error_prob_to_qual_is_the_reference() {
    let text = golden();
    for row in rows(&text, "errorprob") {
        let rate: f64 = row[0].parse().unwrap();
        assert_eq!(
            error_prob_to_qual(rate).unwrap().to_string(),
            row[1],
            "errorProbToQual({rate})"
        );
    }
    // The three the reference refuses, which are not probabilities.
    for row in rows(&text, "error") {
        if let Some(rate) = row[0].strip_prefix("errorProbToQual@") {
            let rate: f64 = rate.parse().unwrap();
            assert_eq!(error_prob_to_qual(rate), None, "errorProbToQual({rate})");
        }
    }
}

/// Every quantization map, for every histogram at every level count.
#[test]
fn every_quantization_map_is_the_reference() {
    let text = golden();
    let maps = rows(&text, "map");
    assert!(maps.len() >= 50, "six histograms at nine settings each");
    for row in &maps {
        let (histogram, levels, min_interesting) = arguments(row[0]);
        assert_eq!(levels.to_string(), row[1], "{}: levels", row[0]);
        assert_eq!(min_interesting.to_string(), row[2], "{}: min", row[0]);
        let quantizer = QualQuantizer::new(&histogram, levels, min_interesting).unwrap();
        let ours: Vec<String> = quantizer
            .original_to_quantized_map
            .iter()
            .map(|qual| qual.to_string())
            .collect();
        assert_eq!(ours.join(","), row[3], "{}: map", row[0]);
    }
    println!("qual-quantizer: {} maps compared", maps.len());
}

/// Every interval of every final forest, with its counts, its level and its two doubles as bits.
#[test]
fn every_interval_of_every_forest_is_the_reference() {
    let text = golden();
    let intervals = rows(&text, "interval");
    assert!(!intervals.is_empty());

    let mut label = String::new();
    let mut ours: Vec<Vec<String>> = Vec::new();
    let mut compared = 0;
    for row in &intervals {
        if row[0] != label {
            label = row[0].to_string();
            let (histogram, levels, min_interesting) = arguments(&label);
            let quantizer = QualQuantizer::new(&histogram, levels, min_interesting).unwrap();
            ours = quantizer
                .quantized_intervals
                .iter()
                .map(|interval| {
                    vec![
                        interval.name(),
                        interval.n_observations.to_string(),
                        interval.n_errors.to_string(),
                        interval.level.to_string(),
                        interval.fixed_qual.to_string(),
                        interval.qual().to_string(),
                        bits(interval.error_rate()),
                        bits(interval.penalty(min_interesting)),
                    ]
                })
                .collect();
        }
        let theirs: Vec<String> = row[1..].iter().map(|field| field.to_string()).collect();
        let found = ours
            .iter()
            .find(|ours| ours[0] == theirs[0])
            .unwrap_or_else(|| panic!("{label}: no interval {}", theirs[0]));
        assert_eq!(*found, theirs, "{label}: interval {}", theirs[0]);
        compared += 1;
    }
    println!("qual-quantizer: {compared} intervals compared");
}

/// The forests have the same size and the same members, not merely the same intervals somewhere.
#[test]
fn the_forests_are_the_same_size_as_the_references() {
    let text = golden();
    let mut counts: Vec<(String, usize)> = Vec::new();
    for row in rows(&text, "interval") {
        match counts.last_mut() {
            Some((label, count)) if label == row[0] => *count += 1,
            _ => counts.push((row[0].to_string(), 1)),
        }
    }
    for (label, count) in counts {
        let (histogram, levels, min_interesting) = arguments(&label);
        let quantizer = QualQuantizer::new(&histogram, levels, min_interesting).unwrap();
        assert_eq!(quantizer.quantized_intervals.len(), count, "{label}");
    }
}

/// The level count, `noQuantization`, and the level count after it.
#[test]
fn quantization_info_is_the_reference() {
    let text = golden();
    let levels = rows(&text, "levels");
    let value = |label: &str| -> String {
        levels
            .iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no levels row {label}"))[1]
            .to_string()
    };

    for name in [
        "flat",
        "empty",
        "illumina",
        "single",
        "short",
        "overflowing",
    ] {
        let counts = padded(name);
        let quantizer = QualQuantizer::new(&counts, 16, MIN_USABLE_Q_SCORE).unwrap();
        let mut info = QuantizationInfo::new(quantizer.original_to_quantized_map, counts);
        assert_eq!(info.quantization_levels.to_string(), value(name), "{name}");

        info.no_quantization();
        let map: Vec<String> = info
            .quantized_quals
            .iter()
            .map(|qual| qual.to_string())
            .collect();
        let expected = rows(&text, "noquant")
            .into_iter()
            .find(|row| row[0] == name)
            .unwrap_or_else(|| panic!("no noquant row {name}"))[1]
            .to_string();
        assert_eq!(map.join(","), expected, "{name}: after noQuantization");
        assert_eq!(
            info.quantization_levels.to_string(),
            value(&format!("{name}-after-noquant")),
            "{name}: levels after noQuantization"
        );
    }

    // A map that returns to a value it already used: four levels for three distinct values.
    assert_eq!(
        calculate_quantization_levels(&[2, 2, 10, 10, 2, 2, 30]).to_string(),
        value("repeating")
    );
}

/// The report table `QuantizationInfo` writes, which is a recalibration report's third table.
#[test]
fn the_report_table_is_the_reference_space_for_space() {
    use gatk_engine::gatk_report::{Report, Sorting, Table, Value};

    let text = golden();
    for name in [
        "flat",
        "empty",
        "illumina",
        "single",
        "short",
        "overflowing",
    ] {
        let counts = padded(name);
        let quantizer = QualQuantizer::new(&counts, 16, MIN_USABLE_Q_SCORE).unwrap();
        let map = quantizer.original_to_quantized_map;

        let mut table = Table::new(
            "Quantized",
            "Quality quantization map",
            Sorting::SortByColumn,
        );
        table.add_column("QualityScore", "%d");
        table.add_column("Count", "%d");
        table.add_column("QuantizedScore", "%d");
        for (qual, count) in counts.iter().enumerate() {
            let key = qual.to_string();
            table.set(&key, "QualityScore", Value::Int(qual as i64));
            table.set(&key, "Count", Value::Int(*count));
            table.set(&key, "QuantizedScore", Value::Int(map[qual] as i64));
        }
        let mut report = Report::new();
        report.add_table(table);
        let ours: Vec<String> = report
            .write()
            .split('\n')
            .map(|line| line.replace(' ', "_"))
            .collect();

        let mut expected: Vec<(usize, String)> = rows(&text, "table")
            .into_iter()
            .filter(|row| row[0] == name)
            .map(|row| (row[1].parse::<usize>().unwrap(), row[2].to_string()))
            .collect();
        expected.sort_by_key(|(n, _)| *n);
        let expected: Vec<String> = expected.into_iter().map(|(_, line)| line).collect();
        assert_eq!(ours.len(), expected.len(), "{name}: line count");
        for (n, (ours, theirs)) in ours.iter().zip(&expected).enumerate() {
            assert_eq!(ours, theirs, "{name}: line {n}");
        }
    }
}

/// Every argument the quantizer refuses, and the two ends it does not refuse but cannot survive.
#[test]
fn the_refusals_are_worded_like_the_reference() {
    let text = golden();
    let message = |what: &str| -> (String, String) {
        let row = rows(&text, "error")
            .into_iter()
            .find(|row| row[0] == what)
            .unwrap_or_else(|| panic!("no error row {what}"));
        (row[1].to_string(), row[2].to_string())
    };

    let mut negative = vec![1i64; 10];
    negative[3] = -1;
    let (exception, words) = message("negative-counts");
    assert_eq!(exception, "GATKException");
    assert_eq!(
        QualQuantizer::new(&negative, 4, 6).unwrap_err().message(),
        words
    );

    assert_eq!(
        QualQuantizer::new(&[1; 10], -1, 6).unwrap_err().message(),
        message("negative-levels").1
    );
    assert_eq!(
        QualQuantizer::new(&[1; 10], 4, -1).unwrap_err().message(),
        message("negative-min-interesting").1
    );

    // Not refused, and not survivable.
    let (exception, words) = message("zero-levels");
    assert_eq!(exception, "NullPointerException");
    assert_eq!(
        QualQuantizer::new(&[1; 10], 0, 6).unwrap_err(),
        QuantizerError::NoPairToMerge
    );
    assert_eq!(
        QualQuantizer::new(&[1; 10], 0, 6).unwrap_err().message(),
        words
    );

    let (exception, words) = message("empty-histogram");
    assert_eq!(exception, "NoSuchElementException");
    assert_eq!(
        QualQuantizer::new(&[], 4, 6).unwrap_err().message(),
        words,
        "the reference's exception carries no message"
    );

    // One bin and one level needs no merging at all, and the leaf keeps its own quality of zero.
    let (exception, result) = message("one-bin-one-level");
    assert_eq!(exception, "none");
    let quantizer = QualQuantizer::new(&[5], 1, 6).unwrap();
    assert_eq!(
        quantizer
            .original_to_quantized_map
            .iter()
            .map(|qual| qual.to_string())
            .collect::<Vec<_>>()
            .join(","),
        result
    );
}

//! Conformance for `RecalDatum` and `EventType` against GATK 4.6.2.0, compared **bit for bit**.
//!
//! Golden from `tools/readfilter-conformance/RecalDatumDump.java`. Every double in it travels as
//! its raw bit pattern as well as its decimal, because a decimal rendering cannot show the last
//! bits and those are what a recalibration table is made of.
//!
//! # What this suite is for
//!
//! A recalibration table is an array of these, so every number `BaseRecalibrator` writes and every
//! number `ApplyBQSR` reads comes from one:
//!
//!  * **the mismatch count is stored multiplied by 100000**, which is not the number the comment
//!    beside the constant names;
//!  * **the empirical quality is an integer, cached, and cleared by every setter**, so a second
//!    prior gets the first prior's answer back;
//!  * **the prior cache has 41 entries and the search has 61 bins**, and the cast in `getLogPrior`
//!    runs before the absolute value;
//!  * **the smoothing is applied twice with different arithmetic**, once in doubles and once
//!    through a truncating cast after adding a half;
//!  * **`combine` recomputes the reported quality from both sides' expected errors**, which for two
//!    empty datums is NaN, and combining onto that NaN throws.

use gatk_corpus as corpus;
use gatk_engine::recal_datum::{
    bayesian_estimate_of_empirical_quality, get_log_binomial_likelihood, get_log_prior,
    log_prior_cache, EventType, RecalDatum, MAX_GATK_USABLE_Q_SCORE, MAX_REASONABLE_Q_SCORE,
    MAX_RECALIBRATED_Q_SCORE,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/recal_datum.txt.gz"),
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

/// `Long.toHexString(Double.doubleToRawLongBits(v))`, which is how the golden writes a double.
fn bits(value: f64) -> String {
    format!("{:x}", value.to_bits())
}

/// Compare one double against the golden's bit pattern, with **one** exception.
///
/// A NaN is compared as a NaN and not as bits, because the bits of a NaN produced by `0.0 / 0.0`
/// are chosen by the processor and not by the program. x86-64 answers the "floating-point
/// indefinite" pattern, whose **sign bit is set**, and AArch64 answers the default NaN, whose sign
/// bit is clear. The reference and the port make the same NaN the same way and still disagree in
/// that one bit depending on where they run, which is not a difference in the port. See
/// [`the_nan_sign_is_the_processors_and_not_the_programs`], which asserts the golden carries the
/// x86-64 pattern, since that is the architecture goldens are produced on.
fn assert_double(ours: f64, golden: &str, what: &str) {
    let expected = u64::from_str_radix(golden, 16).expect("the golden's bits are hexadecimal");
    if ours.is_nan() {
        assert!(
            f64::from_bits(expected).is_nan(),
            "{what}: the port answered NaN and the reference did not"
        );
    } else {
        assert_eq!(ours.to_bits(), expected, "{what}");
    }
}

fn constant(text: &str, name: &str) -> String {
    rows(text, "const")
        .into_iter()
        .find(|row| row[0] == name)
        .unwrap_or_else(|| panic!("the golden has no constant {name}"))[1]
        .to_string()
}

/// The same datums the dump built, by the same names, in the same order.
fn datum_for(label: &str) -> RecalDatum {
    match label {
        "plain" => RecalDatum::new(1000, 10.0, 30).unwrap(),
        "empty" => RecalDatum::new(0, 0.0, 30).unwrap(),
        "perfect" => RecalDatum::new(1000, 0.0, 30).unwrap(),
        "all-errors" => RecalDatum::new(1000, 1000.0, 30).unwrap(),
        "qual-zero" => RecalDatum::new(1000, 10.0, 0).unwrap(),
        "qual-max" => RecalDatum::new(1000, 10.0, 93).unwrap(),
        "fractional" => RecalDatum::new(1000, 0.1, 30).unwrap(),
        "tiny-fraction" => RecalDatum::new(1000, 1.0e-7, 30).unwrap(),
        "half" => RecalDatum::new(1000, 0.5, 30).unwrap(),
        "just-under-half" => RecalDatum::new(1000, 0.49999, 30).unwrap(),
        "huge" => RecalDatum::new(3_000_000_000, 3000.0, 30).unwrap(),
        "capped" => RecalDatum::new(100_000_000, 0.0, 93).unwrap(),
        // The copy carries the cached quality, because the copy constructor copies every field.
        "copy-of-computed" => {
            let mut source = RecalDatum::new(1000, 10.0, 30).unwrap();
            source.empirical_quality();
            source.clone()
        }
        "after-setters" => {
            let mut datum = RecalDatum::new(1000, 10.0, 30).unwrap();
            datum.set_num_mismatches(0.1).unwrap();
            datum.set_num_observations(7).unwrap();
            datum.set_reported_quality(20.5).unwrap();
            datum
        }
        "forced-empirical" => {
            let mut datum = RecalDatum::new(1000, 10.0, 30).unwrap();
            datum.set_empirical_quality(7).unwrap();
            datum
        }
        "incremented" => {
            let mut datum = RecalDatum::new(0, 0.0, 30).unwrap();
            datum.increment_by_observation(true);
            datum.increment_by_observation(false);
            datum.increment_by_observation(false);
            datum.increment_num_observations(10);
            datum.increment_num_mismatches(0.25);
            datum
        }
        other => panic!("{other} is in the golden but not configured here"),
    }
}

fn combined_for(label: &str) -> RecalDatum {
    let mut datum = match label {
        "different-qualities" | "same-quality" | "one-empty" | "after-computing" => {
            RecalDatum::new(1000, 10.0, 30).unwrap()
        }
        "both-empty" => RecalDatum::new(0, 0.0, 30).unwrap(),
        "with-quality-zero" => RecalDatum::new(1000, 10.0, 40).unwrap(),
        other => panic!("{other} is in the golden but not configured here"),
    };
    // The cache is filled before the combine in exactly one case, to show the combine clears it.
    if label == "after-computing" {
        datum.empirical_quality();
    }
    let other = match label {
        "different-qualities" => RecalDatum::new(1000, 10.0, 20).unwrap(),
        "same-quality" => RecalDatum::new(1000, 10.0, 30).unwrap(),
        "both-empty" => RecalDatum::new(0, 0.0, 30).unwrap(),
        "one-empty" => RecalDatum::new(0, 0.0, 30).unwrap(),
        "with-quality-zero" => RecalDatum::new(1000, 10.0, 0).unwrap(),
        "after-computing" => RecalDatum::new(1000, 500.0, 30).unwrap(),
        other => unreachable!("{other}"),
    };
    datum.combine(&other).unwrap();
    datum
}

#[test]
fn the_constants_are_the_references() {
    let text = golden();
    assert_eq!(
        constant(&text, "MAX_RECALIBRATED_Q_SCORE"),
        MAX_RECALIBRATED_Q_SCORE.to_string()
    );
    assert_eq!(
        constant(&text, "MAX_GATK_USABLE_Q_SCORE"),
        MAX_GATK_USABLE_Q_SCORE.to_string()
    );
    assert_eq!(
        constant(&text, "MAX_REASONABLE_Q_SCORE"),
        MAX_REASONABLE_Q_SCORE.to_string()
    );
    // The one the comment beside it contradicts.
    assert_eq!(constant(&text, "MULTIPLIER"), "100000.0");
    assert_eq!(constant(&text, "SMOOTHING_CONSTANT"), "1");
    assert_eq!(constant(&text, "UNINITIALIZED_EMPIRICAL_QUALITY"), "-1");
    assert_eq!(constant(&text, "MAX_NUMBER_OF_OBSERVATIONS"), "2147483646");
    assert_eq!(
        constant(&text, "logPriorCache.length"),
        log_prior_cache().len().to_string()
    );
}

/// Every entry of the Gaussian prior, bit for bit. This is where `FastMath.log` shows up or does
/// not: the reference's logarithm is table-driven and not correctly rounded, so a port built on the
/// platform's would differ here and nowhere a decimal could show.
#[test]
fn the_prior_cache_is_the_reference_bit_for_bit() {
    let text = golden();
    let entries = rows(&text, "logprior");
    assert_eq!(entries.len(), 41, "the cache is 41 entries");
    for row in entries {
        let index: usize = row[0].parse().unwrap();
        assert_eq!(
            bits(log_prior_cache()[index]),
            row[1],
            "logPriorCache[{index}]"
        );
    }
}

#[test]
fn the_prior_lookup_is_the_reference() {
    let text = golden();
    for row in rows(&text, "getlogprior") {
        let quality: f64 = row[0].parse().unwrap();
        let prior: f64 = row[1].parse().unwrap();
        assert_eq!(
            bits(get_log_prior(quality, prior)),
            row[2],
            "getLogPrior({quality}, {prior})"
        );
    }
}

#[test]
fn the_binomial_likelihood_is_the_reference() {
    let text = golden();
    for row in rows(&text, "loglik") {
        let quality: f64 = row[0].parse().unwrap();
        let observations: i64 = row[1].parse().unwrap();
        let errors: i64 = row[2].parse().unwrap();
        assert_eq!(
            bits(get_log_binomial_likelihood(quality, observations, errors)),
            row[3],
            "getLogBinomialLikelihood({quality}, {observations}, {errors})"
        );
    }
}

#[test]
fn the_posterior_maximum_is_the_reference() {
    let text = golden();
    let cases = rows(&text, "bayes");
    assert!(!cases.is_empty());
    for row in cases {
        let observations: i64 = row[0].parse().unwrap();
        let errors: i64 = row[1].parse().unwrap();
        let prior: f64 = row[2].parse().unwrap();
        assert_eq!(
            bayesian_estimate_of_empirical_quality(observations, errors, prior).to_string(),
            row[3],
            "bayesianEstimateOfEmpiricalQuality({observations}, {errors}, {prior})"
        );
    }
}

/// Every reader of every datum, in the order the dump asked for them, because the order decides
/// what the cache holds.
#[test]
fn every_datum_reads_back_bit_for_bit() {
    let text = golden();
    let texts = rows(&text, "text");
    assert!(!texts.is_empty());

    for row in &texts {
        let label = row[0];
        let mut datum = datum_for(label);
        let expected = |field: &str| -> Vec<String> {
            rows(&text, "datum")
                .into_iter()
                .filter(|r| r[0] == label && r[1] == field)
                .map(|r| r[2].to_string())
                .collect()
        };
        let expected_value = |field: &str| -> Vec<String> {
            rows(&text, "datum")
                .into_iter()
                .filter(|r| r[0] == label && r[1] == field)
                .map(|r| r[3].to_string())
                .collect()
        };

        assert_eq!(
            expected("numObservations"),
            vec![bits(datum.num_observations() as f64)],
            "{label}: numObservations"
        );
        assert_eq!(
            expected("numMismatches"),
            vec![bits(datum.num_mismatches())],
            "{label}: numMismatches"
        );
        assert_eq!(
            expected("reportedQuality"),
            vec![bits(datum.reported_quality())],
            "{label}: reportedQuality"
        );
        assert_eq!(
            expected("empiricalErrorRate"),
            vec![bits(datum.empirical_error_rate())],
            "{label}: empiricalErrorRate"
        );
        assert_eq!(
            expected_value("reportedQualityAsByte"),
            vec![datum.reported_quality_as_byte().to_string()],
            "{label}: reportedQualityAsByte"
        );
        // Before the quality getters, because these two compute and cache it themselves.
        assert_eq!(datum.to_text(), row[1], "{label}: toString");
        assert_eq!(datum.string_for_csv(), row[2], "{label}: stringForCSV");
        assert_eq!(
            expected("empiricalQuality"),
            vec![bits(datum.empirical_quality())],
            "{label}: empiricalQuality"
        );
        assert_eq!(
            expected_value("empiricalQualityAsByte"),
            vec![datum.empirical_quality_as_byte().to_string()],
            "{label}: empiricalQualityAsByte"
        );
    }
    println!("recal-datum: {} datums read back", texts.len());
}

/// The cache and its invalidation, which is the behaviour a reader of the class would get wrong.
#[test]
fn the_empirical_quality_cache_is_the_reference() {
    let text = golden();
    let steps = rows(&text, "cache");
    assert_eq!(steps.len(), 9);
    let value = |label: &str, step: &str| -> String {
        steps
            .iter()
            .find(|row| row[0] == label && row[1] == step)
            .unwrap_or_else(|| panic!("no cache row {label}/{step}"))[2]
            .to_string()
    };
    let decimal = |v: f64| format!("{v:?}");

    let mut datum = RecalDatum::new(1000, 10.0, 30).unwrap();
    assert_eq!(
        decimal(datum.empirical_quality_with_prior(10.0)),
        value("prior-then-prior", "first-with-10")
    );
    assert_eq!(
        decimal(datum.empirical_quality_with_prior(45.0)),
        value("prior-then-prior", "then-with-45")
    );
    datum.set_num_observations(1000).unwrap();
    assert_eq!(
        decimal(datum.empirical_quality_with_prior(45.0)),
        value("prior-then-prior", "after-setter-with-45")
    );

    let mut other = RecalDatum::new(1000, 10.0, 30).unwrap();
    assert_eq!(
        decimal(other.empirical_quality_with_prior(45.0)),
        value("reversed", "first-with-45")
    );
    assert_eq!(
        decimal(other.empirical_quality_with_prior(10.0)),
        value("reversed", "then-with-10")
    );

    let mut implicit = RecalDatum::new(1000, 10.0, 30).unwrap();
    assert_eq!(
        decimal(implicit.empirical_quality()),
        value("implicit", "no-argument")
    );
    assert_eq!(
        decimal(implicit.empirical_quality_with_prior(0.0)),
        value("implicit", "then-with-0")
    );

    let mut forced = RecalDatum::new(1000, 10.0, 30).unwrap();
    forced.set_empirical_quality(3).unwrap();
    assert_eq!(
        decimal(forced.empirical_quality_with_prior(45.0)),
        value("forced", "with-45")
    );
    // Adds nothing, and invalidates anyway.
    forced.increment_num_observations(0);
    assert_eq!(
        decimal(forced.empirical_quality_with_prior(45.0)),
        value("forced", "after-zero-increment")
    );
}

#[test]
fn combining_is_the_reference() {
    let text = golden();
    let all = rows(&text, "combine");
    assert!(!all.is_empty());
    let labels: Vec<String> = {
        let mut seen: Vec<String> = Vec::new();
        for row in &all {
            if !seen.iter().any(|s| s == row[0]) {
                seen.push(row[0].to_string());
            }
        }
        seen
    };
    assert_eq!(labels.len(), 6);

    for label in &labels {
        let mut datum = combined_for(label);
        let expected = |field: &str| -> String {
            all.iter()
                .find(|row| row[0] == label && row[1] == field)
                .unwrap_or_else(|| panic!("no combine row {label}/{field}"))[2]
                .to_string()
        };
        assert_double(
            datum.num_observations() as f64,
            &expected("numObservations"),
            &format!("{label}: numObservations"),
        );
        assert_double(
            datum.num_mismatches(),
            &expected("numMismatches"),
            &format!("{label}: numMismatches"),
        );
        assert_double(
            datum.reported_quality(),
            &expected("reportedQuality"),
            &format!("{label}: reportedQuality"),
        );
        assert_double(
            datum.empirical_quality(),
            &expected("empiricalQuality"),
            &format!("{label}: empiricalQuality"),
        );
    }
}

/// The two narrowing casts, which are the only way a quality score comes back negative.
#[test]
fn the_byte_getters_narrow_like_the_reference() {
    let text = golden();
    let narrow = rows(&text, "narrow");
    let value = |what: &str| -> String {
        narrow
            .iter()
            .find(|row| row[0] == what)
            .unwrap_or_else(|| panic!("no narrow row {what}"))[1]
            .to_string()
    };

    let mut datum = RecalDatum::new(1, 0.0, 30).unwrap();
    for (quality, what) in [
        (200.0, "reported-200"),
        (127.4, "reported-127.4"),
        (127.5, "reported-127.5"),
        (0.5, "reported-0.5"),
        (1.5, "reported-1.5"),
    ] {
        datum.set_reported_quality(quality).unwrap();
        assert_eq!(
            datum.reported_quality_as_byte().to_string(),
            value(what),
            "{what}"
        );
    }

    let mut forced = RecalDatum::new(1, 0.0, 30).unwrap();
    forced.set_empirical_quality(200).unwrap();
    assert_eq!(
        forced.empirical_quality_as_byte().to_string(),
        value("empirical-200")
    );
    forced.set_empirical_quality(93).unwrap();
    assert_eq!(
        forced.empirical_quality_as_byte().to_string(),
        value("empirical-93")
    );
}

/// Every argument the class refuses, with the words it refuses it in.
#[test]
fn the_refusals_are_worded_like_the_reference() {
    let text = golden();
    let errors = rows(&text, "error");
    let message = |what: &str| -> String {
        errors
            .iter()
            .find(|row| row[0] == what)
            .unwrap_or_else(|| panic!("no error row {what}"))[2]
            .to_string()
    };

    assert_eq!(
        RecalDatum::new(-1, 0.0, 30).unwrap_err().message(),
        message("constructor-negative-observations")
    );
    assert_eq!(
        RecalDatum::new(1, -1.0, 30).unwrap_err().message(),
        message("constructor-negative-mismatches")
    );
    assert_eq!(
        RecalDatum::new(1, 0.0, -1).unwrap_err().message(),
        message("constructor-negative-quality")
    );

    let mut datum = RecalDatum::new(1, 0.0, 30).unwrap();
    assert_eq!(
        datum.set_num_observations(-1).unwrap_err().message(),
        message("set-negative-observations")
    );
    assert_eq!(
        datum.set_num_mismatches(-1.0).unwrap_err().message(),
        message("set-negative-mismatches")
    );
    assert_eq!(
        datum.set_reported_quality(-1.0).unwrap_err().message(),
        message("set-negative-reported")
    );
    assert_eq!(
        datum
            .set_reported_quality(f64::INFINITY)
            .unwrap_err()
            .message(),
        message("set-infinite-reported")
    );
    assert_eq!(
        datum.set_reported_quality(f64::NAN).unwrap_err().message(),
        message("set-nan-reported")
    );
    assert_eq!(
        datum.set_empirical_quality(-1).unwrap_err().message(),
        message("set-negative-empirical")
    );

    // Not a refusal: increment checks nothing, so the counts go negative.
    let mut negative = RecalDatum::new(1, 0.0, 30).unwrap();
    negative.increment(-5, -5.0);
    assert_eq!(
        format!(
            "{},{:?}",
            negative.num_observations(),
            negative.num_mismatches()
        ),
        message("increment-negative-is-allowed")
    );
    // Nor is a reported quality above the recalibration ceiling.
    assert_eq!(
        format!(
            "{:?}",
            RecalDatum::new(1, 0.0, 127).unwrap().reported_quality()
        ),
        message("quality-above-max-is-allowed")
    );

    // The reachable end of a NaN reported quality: the guard inside qualToErrorProb.
    let mut nan = RecalDatum::new(0, 0.0, 30).unwrap();
    nan.combine(&RecalDatum::new(0, 0.0, 30).unwrap()).unwrap();
    assert_eq!(
        nan.combine(&RecalDatum::new(1000, 10.0, 30).unwrap())
            .unwrap_err()
            .message(),
        message("combine-onto-nan-reported")
    );

    // And the two ways an event type is refused, which are not the same exception.
    assert_eq!(EventType::from_representation("X"), None);
    assert_eq!(message("event-from-unknown"), "Event X does not exist.");
    assert_eq!(EventType::from_index(3), None);
    assert_eq!(EventType::from_index(-1), None);
}

/// The one double in this suite that is not compared bit for bit, and why.
///
/// `combine` on two empty datums computes `-10 * log10(0.0 / 0.0)`. The value is NaN in every
/// implementation; **which** NaN is the processor's choice. x86-64 produces the floating-point
/// indefinite, `0xFFF8000000000000`, whose sign bit is set, and AArch64 produces the default NaN,
/// `0x7FF8000000000000`, whose sign bit is clear. The reference and the port agree on the
/// architecture they are run on and disagree across architectures, which is why this is asserted
/// about the golden rather than about the port.
///
/// Goldens are produced on real x86-64 (decision 0004), so the pattern below is the one that will
/// always be in the file.
#[test]
fn the_nan_sign_is_the_processors_and_not_the_programs() {
    let text = golden();
    let reported = rows(&text, "combine")
        .into_iter()
        .find(|row| row[0] == "both-empty" && row[1] == "reportedQuality")
        .expect("the golden lost the both-empty combine")[2]
        .to_string();
    assert_eq!(reported, "fff8000000000000", "the x86-64 indefinite");
    assert!(f64::from_bits(u64::from_str_radix(&reported, 16).unwrap()).is_nan());

    let mut empty = RecalDatum::new(0, 0.0, 30).unwrap();
    empty
        .combine(&RecalDatum::new(0, 0.0, 30).unwrap())
        .unwrap();
    assert!(empty.reported_quality().is_nan());
}

#[test]
fn the_event_types_are_the_references() {
    let text = golden();
    let events = rows(&text, "event");
    assert_eq!(events.len(), 3);
    for row in events {
        let ordinal: usize = row[0].parse().unwrap();
        let event = EventType::from_index(ordinal as i32).unwrap();
        assert_eq!(event.name(), row[1]);
        assert_eq!(event.representation(), row[2]);
        assert_eq!(event.pretty_print(), row[3]);
        // The letter round-trips to the same constant, which is how a report's tables are read.
        assert_eq!(
            EventType::from_representation(row[2]).unwrap().name(),
            row[4]
        );
        assert_eq!(event.ordinal(), ordinal);
    }
}

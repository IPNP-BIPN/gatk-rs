//! `VariantRecalibrator`: the tranches a truth sensitivity target cuts out of a scored callset.
//!
//! The Gaussian mixture that produces the scores is not ported. What is ported is everything the
//! tool does with them: the running sensitivity, the search for each target's tranche, the counts
//! a tranche reports, and the two tranche file formats.
//!
//! Ported from `org.broadinstitute.hellbender.tools.walkers.vqsr.TrancheManager`,
//! `org.broadinstitute.hellbender.tools.walkers.vqsr.Tranche`,
//! `org.broadinstitute.hellbender.tools.walkers.vqsr.TruthSensitivityTranche` and
//! `org.broadinstitute.hellbender.tools.walkers.vqsr.VQSLODTranche` in GATK 4.7.0.0; the four
//! `.java` files are unchanged since 4.6.2.0.
//!
//! # Reading this file without reading Rust
//!
//! - `&[Datum]` is a borrowed, read-only view of a list; `.to_vec()` makes an owned copy, which
//!   is how [`sorted`] avoids the Java's in-place `Collections.sort` of the caller's list.
//! - `std::cmp::Ordering` is `Less`, `Equal` or `Greater`, the sign of a Java `compare` result.
//! - `Option<Tranche>` is a tranche or nothing, the Java's `null` return of `findTranche`;
//!   `Result<Vec<Tranche>, String>` is the tranches or the `UserException` text.
//! - `Tranche { min_vqs_lod: threshold, ..other }` builds a copy of `other` with one field
//!   replaced.
//! - `as f64` turns a count into a `double`, as Java's implicit widening in `n / (1.0 * d)`.

/// One scored variant, reduced to what the tranche arithmetic reads off it.
#[derive(Debug, Clone, PartialEq)]
pub struct Datum {
    /// The VQSLOD score, `VariantDatum.lod`.
    pub lod: f64,
    /// Present in a resource tagged `known=true`.
    pub is_known: bool,
    /// Present in a resource tagged `truth=true`.
    pub at_truth_site: bool,
    pub is_snp: bool,
    pub is_transition: bool,
}

/// `SNP`, `INDEL` or `BOTH`, which is only ever written into the filter name and the model column.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    Snp,
    Indel,
    Both,
}

impl Mode {
    /// The enum constant's name, which is what `model.toString()` prints in a tranche row.
    pub fn name(self) -> &'static str {
        match self {
            Mode::Snp => "SNP",
            Mode::Indel => "INDEL",
            Mode::Both => "BOTH",
        }
    }
}

/// `VariantDatumLODComparator`, which is `Double.compare` on the LOD.
///
/// `Double.compare` is a total order: it separates -0.0 from 0.0 and puts NaN last, which a plain
/// comparison would not.
pub fn compare_lod(a: &Datum, b: &Datum) -> std::cmp::Ordering {
    compare_double(a.lod, b.lod)
}

/// `Double.compare`: numeric first, then the signed bit pattern.
///
/// Only equal values and NaN reach the bit comparison. Every NaN maps to the one canonical NaN,
/// as `doubleToLongBits` does, so NaN sorts above `+inf` and all NaNs are equal; `-0.0` has the
/// sign bit set, so as an `i64` it is negative and sorts below `0.0`.
fn compare_double(a: f64, b: f64) -> std::cmp::Ordering {
    if a < b {
        return std::cmp::Ordering::Less;
    }
    if a > b {
        return std::cmp::Ordering::Greater;
    }
    let bits = |value: f64| {
        if value.is_nan() {
            f64::NAN.to_bits() as i64
        } else {
            value.to_bits() as i64
        }
    };
    bits(a).cmp(&bits(b))
}

/// The variants in LOD order, which is what every step below assumes.
///
/// `Collections.sort( data, VariantDatum.VariantDatumLODComparator )`. Both sorts are stable,
/// so equal LODs keep their input order.
pub fn sorted(data: &[Datum]) -> Vec<Datum> {
    let mut sorted = data.to_vec();
    sorted.sort_by(compare_lod);
    sorted
}

/// `TruthSensitivityMetric.calculateRunningMetric`: for each index, one less the truth sites at
/// or above it over the truth sites in total.
///
/// The walk runs from the TOP DOWN, so the entry at zero is one less the whole sensitivity and
/// the last entry is one less the last variant's own contribution. `n_true_sites` is the count
/// the tool was given rather than the count in `data`, which is why it is a parameter.
///
/// ```java
/// for ( int i = data.size() - 1; i >= 0; i-- ) {
///     VariantDatum datum = data.get(i);
///     nCalledAtTruth += datum.atTruthSite ? 1 : 0;
///     runningSensitivity[i] = 1 - nCalledAtTruth / (1.0 * nTrueSites);
/// }
/// ```
pub fn running_sensitivity(data: &[Datum], n_true_sites: usize) -> Vec<f64> {
    let mut running = vec![0.0; data.len()];
    let mut called_at_truth = 0;
    // `(0..n).rev()` counts n-1 down to 0, the Java's `i--` loop.
    for i in (0..data.len()).rev() {
        if data[i].at_truth_site {
            called_at_truth += 1;
        }
        running[i] = 1.0 - called_at_truth as f64 / (1.0 * n_true_sites as f64);
    }
    running
}

/// `TruthSensitivityMetric.getThreshold`: a target of 99 asks for a running value of 0.01.
pub fn sensitivity_threshold(target: f64) -> f64 {
    1.0 - target / 100.0
}

/// `countCallsAtTruth`: the truth sites whose LOD reaches `minimum`.
///
/// `VariantDatum.countCallsAtTruth` counts a `long` and casts it, `//XXX cast to int for
/// compatibility`; `as i32` is that cast.
pub fn calls_at_truth(data: &[Datum], minimum: f64) -> i32 {
    data.iter()
        .filter(|datum| datum.at_truth_site && datum.lod >= minimum)
        .count() as i32
}

/// One row of a tranche file.
#[derive(Debug, Clone, PartialEq)]
pub struct Tranche {
    /// The value the first column carries: the target sensitivity, or the REQUESTED score.
    pub index: f64,
    /// Known variants at or above the minimum LOD; a `long` in the Java.
    pub num_known: i64,
    /// The rest of those variants.
    pub num_novel: i64,
    /// `knownTi / Math.max(1.0 * knownTv, 1.0)`, see [`tranche_of_variants`].
    pub known_ti_tv: f64,
    pub novel_ti_tv: f64,
    /// The LOD that bounds the tranche from below, printed as `minVQSLod`.
    pub min_vqs_lod: f64,
    pub model: Mode,
    /// Truth sites anywhere in the callset (`countCallsAtTruth(data, Double.NEGATIVE_INFINITY)`).
    pub accessible_truth_sites: i32,
    /// Truth sites at or above the minimum LOD.
    pub calls_at_truth_sites: i32,
}

impl Tranche {
    /// `getTruthSensitivity`, which is zero rather than a division when there is no truth site.
    ///
    /// ```java
    /// return accessibleTruthSites > 0 ? callsAtTruthSites / (1.0*accessibleTruthSites) : 0.0;
    /// ```
    pub fn truth_sensitivity(&self) -> f64 {
        if self.accessible_truth_sites > 0 {
            self.calls_at_truth_sites as f64 / (1.0 * self.accessible_truth_sites as f64)
        } else {
            0.0
        }
    }
}

/// `Tranche.trancheOfVariants`: the counts over every variant at or above the LOD at `min_index`.
///
/// The counts are NOT over the variants the walk had reached: they are over the whole callset
/// filtered by that LOD, so a looser target's tranche contains a stricter one's rather than
/// sitting beside it.
///
/// ```java
/// final double knownTiTv = knownTi / Math.max(1.0 * knownTv, 1.0);
/// ```
pub fn tranche_of_variants(data: &[Datum], min_index: usize, index: f64, model: Mode) -> Tranche {
    let min_lod = data[min_index].lod;
    let (mut num_known, mut num_novel) = (0i64, 0i64);
    let (mut known_ti, mut known_tv, mut novel_ti, mut novel_tv) = (0i32, 0i32, 0i32, 0i32);
    for datum in data.iter().filter(|datum| datum.lod >= min_lod) {
        if datum.is_known {
            num_known += 1;
            if datum.is_snp {
                if datum.is_transition {
                    known_ti += 1;
                } else {
                    known_tv += 1;
                }
            }
        } else {
            num_novel += 1;
            if datum.is_snp {
                if datum.is_transition {
                    novel_ti += 1;
                } else {
                    novel_tv += 1;
                }
            }
        }
    }
    Tranche {
        index,
        num_known,
        num_novel,
        // The denominator is floored at one rather than guarded, so a tranche with no
        // transversion reports its transition COUNT where a ratio should be.
        known_ti_tv: known_ti as f64 / (known_tv as f64).max(1.0),
        novel_ti_tv: novel_ti as f64 / (novel_tv as f64).max(1.0),
        min_vqs_lod: min_lod,
        model,
        accessible_truth_sites: calls_at_truth(data, f64::NEG_INFINITY),
        calls_at_truth_sites: calls_at_truth(data, min_lod),
    }
}

/// `Tranche.emptyTranche`: the same truth-site counts, and nothing else.
///
/// The counts are zeroed and both ratios with them, but the LOD at `min_index` is still what the
/// row reports as its minimum until the caller overwrites it.
pub fn empty_tranche(data: &[Datum], min_index: usize, index: f64, model: Mode) -> Tranche {
    let min_lod = if data.is_empty() {
        f64::NEG_INFINITY
    } else {
        data[min_index].lod
    };
    Tranche {
        index,
        num_known: 0,
        num_novel: 0,
        known_ti_tv: 0.0,
        novel_ti_tv: 0.0,
        min_vqs_lod: min_lod,
        model,
        accessible_truth_sites: calls_at_truth(data, f64::NEG_INFINITY),
        calls_at_truth_sites: calls_at_truth(data, min_lod),
    }
}

/// `findTranche`: the FIRST index whose running sensitivity reaches the target's threshold.
///
/// Walking upwards from the worst variant means the tranche found is the LARGEST set that still
/// meets the target. A target no index reaches yields nothing.
pub fn find_tranche(data: &[Datum], running: &[f64], target: f64, model: Mode) -> Option<Tranche> {
    let threshold = sensitivity_threshold(target);
    // `find` yields the first `i` that passes, or `None`; `?` then returns `None` from this
    // function, the Java's `return null` after the loop.
    let index = (0..data.len()).find(|i| running[*i] >= threshold)?;
    Some(tranche_of_variants(data, index, target, model))
}

/// `findVQSLODTranche`: the first index whose LOD reaches the requested score.
///
/// Two things separate it from the sensitivity search. The minimum it reports is the REQUEST
/// rather than the LOD it found, and a request no variant reaches still produces a row: an empty
/// tranche over the last variant, again with the request as its minimum.
///
/// That is `emptyTranche(data, n-1, ...)`. With no data the Java passes -1, which
/// `emptyTranche` never reads because it tests `data.isEmpty()` first; `saturating_sub` stops
/// at 0 instead of wrapping below zero, and is equally unread.
pub fn find_vqslod_tranche(data: &[Datum], threshold: f64, model: Mode) -> Tranche {
    match data.iter().position(|datum| datum.lod >= threshold) {
        Some(index) => Tranche {
            min_vqs_lod: threshold,
            ..tranche_of_variants(data, index, threshold, model)
        },
        None => Tranche {
            min_vqs_lod: threshold,
            ..empty_tranche(data, data.len().saturating_sub(1), threshold, model)
        },
    }
}

/// The refusal a first target that no tranche reaches produces.
///
/// `UserException` text with `%s` the metric's name and `%.2f` its threshold, which
/// `format_decimals` renders as Java's `String.format` does.
pub fn no_tranche_refusal(metric: &str, threshold: f64) -> String {
    format!(
        "Couldn't find any tranche containing variants with a {metric} > {}. Are you \
         sure the truth files contain unfiltered variants which overlap the input data?",
        gatk_engine::java_format::format_decimals(threshold, 2)
    )
}

/// `findTranches`: one tranche per target, in the order the targets were given.
///
/// A target no tranche reaches ENDS the list rather than being skipped, so the result is a prefix
/// of the targets; only a FIRST such target is a refusal.
pub fn find_tranches(
    data: &[Datum],
    targets: &[f64],
    n_true_sites: usize,
    model: Mode,
) -> Result<Vec<Tranche>, String> {
    let data = sorted(data);
    let running = running_sensitivity(&data, n_true_sites);
    let mut tranches = Vec::new();
    for target in targets {
        match find_tranche(&data, &running, *target, model) {
            Some(tranche) => tranches.push(tranche),
            None => {
                if tranches.is_empty() {
                    return Err(no_tranche_refusal(
                        "TruthSensitivity",
                        sensitivity_threshold(*target),
                    ));
                }
                break;
            }
        }
    }
    Ok(tranches)
}

/// `findVQSLODTranches`, which never comes up empty and so never refuses.
pub fn find_vqslod_tranches(data: &[Datum], thresholds: &[f64], model: Mode) -> Vec<Tranche> {
    let data = sorted(data);
    thresholds
        .iter()
        .map(|threshold| find_vqslod_tranche(&data, *threshold, model))
        .collect()
}

// ================================================================================================
// The file.
// ================================================================================================

/// The two headers, which differ in their version number and in their first column.
///
/// `TruthSensitivityTranche.CURRENT_VERSION`.
pub const TRUTH_SENSITIVITY_VERSION: i32 = 5;
/// `VQSLODTranche.CURRENT_VERSION`.
pub const VQSLOD_VERSION: i32 = 6;

/// The three `stream.println` lines `TruthSensitivityTranche.printHeader` writes.
pub fn truth_sensitivity_header() -> String {
    format!(
        "# Variant quality score tranches file\n# Version number {TRUTH_SENSITIVITY_VERSION}\n\
         targetTruthSensitivity,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,\
         model,accessibleTruthSites,callsAtTruthSites,truthSensitivity\n"
    )
}

/// The same three lines from `VQSLODTranche`, whose first column is `requestedVQSLOD`.
pub fn vqslod_header() -> String {
    format!(
        "# Variant quality score tranches file\n# Version number {VQSLOD_VERSION}\n\
         requestedVQSLOD,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,model,\
         accessibleTruthSites,callsAtTruthSites,truthSensitivity\n"
    )
}

/// `TrancheComparator`, which orders on CALLS AT TRUTH SITES and on nothing else.
///
/// It is not the target sensitivity, though the two agree whenever the targets were given in
/// increasing order. Given 100, 99.9, 99 and 90 the file comes out 90, 99.9, 99, 100.
///
/// ```java
/// return Double.compare(tranche1.callsAtTruthSites, tranche2.callsAtTruthSites);
/// ```
///
/// The `int` counts are widened to `double` for that call, which the `as f64` repeats.
pub fn tranche_order(a: &Tranche, b: &Tranche) -> std::cmp::Ordering {
    compare_double(a.calls_at_truth_sites as f64, b.calls_at_truth_sites as f64)
}

/// `getTrancheString`: one row, whose filter name takes its LOWER bound from the previous row.
///
/// The bound is the previous row's own index whatever that index is, so a list that is not in
/// increasing order produces a band that runs backwards. The first row has no previous row and
/// so is bounded at 0.00.
///
/// Every decimal is `String.format`'s, which rounds half-up on `Double.toString`'s digits: a Ti/Tv
/// of 363/160 is `2.26875` there and prints `2.2688`, where Rust's own `{:.4}` prints `2.2687`.
///
/// The Java format string, field for field:
///
/// ```java
/// return String.format("%.2f,%d,%d,%.4f,%.4f,%.4f,VQSRTranche%s%.2fto%.2f,%s,%d,%d,%.4f%n",
/// ```
///
/// `%n` is the platform line separator; this port writes `\n`.
pub fn tranche_row(tranche: &Tranche, previous: Option<&Tranche>) -> String {
    // `use ... as f` imports the function under a short local name for the fields below.
    use gatk_engine::java_format::format_decimals as f;
    format!(
        "{},{},{},{},{},{},VQSRTranche{}{}to{},{},{},{},{}\n",
        f(tranche.index, 2),
        tranche.num_known,
        tranche.num_novel,
        f(tranche.known_ti_tv, 4),
        f(tranche.novel_ti_tv, 4),
        f(tranche.min_vqs_lod, 4),
        tranche.model.name(),
        f(previous.map_or(0.0, |previous| previous.index), 2),
        f(tranche.index, 2),
        tranche.model.name(),
        tranche.accessible_truth_sites,
        tranche.calls_at_truth_sites,
        f(tranche.truth_sensitivity(), 4)
    )
}

/// `tranchesString`: the rows, sorted, each naming the one before it.
///
/// The sort is STABLE, so two targets that found the same tranche keep the order they were given
/// in, and a list of one is not sorted at all.
///
/// The Java sorts the caller's list in place with `TrancheComparator` and prints each row with
/// `stream.printf(t.getTrancheString(prev))`; a row carries no `%`, so using it as a format
/// string changes nothing.
pub fn tranches_string(tranches: &[Tranche]) -> String {
    let mut tranches = tranches.to_vec();
    if tranches.len() > 1 {
        tranches.sort_by(tranche_order);
    }
    let mut text = String::new();
    let mut previous: Option<Tranche> = None;
    for tranche in tranches {
        text.push_str(&tranche_row(&tranche, previous.as_ref()));
        previous = Some(tranche);
    }
    text
}

/// The whole truth-sensitivity file: header, then the sorted rows.
pub fn truth_sensitivity_file(tranches: &[Tranche]) -> String {
    truth_sensitivity_header() + &tranches_string(tranches)
}

/// The whole VQSLOD file.
pub fn vqslod_file(tranches: &[Tranche]) -> String {
    vqslod_header() + &tranches_string(tranches)
}

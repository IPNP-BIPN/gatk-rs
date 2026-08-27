//! `VariantRecalibrator`: the tranches a truth sensitivity target cuts out of a scored callset.
//!
//! The Gaussian mixture that produces the scores is not ported. What is ported is everything the
//! tool does with them: the running sensitivity, the search for each target's tranche, the counts
//! a tranche reports, and the two tranche file formats.
//!
//! Ported from `org.broadinstitute.hellbender.tools.walkers.vqsr.TrancheManager`,
//! `org.broadinstitute.hellbender.tools.walkers.vqsr.Tranche`,
//! `org.broadinstitute.hellbender.tools.walkers.vqsr.TruthSensitivityTranche` and
//! `org.broadinstitute.hellbender.tools.walkers.vqsr.VQSLODTranche` in GATK 4.6.2.0.

/// One scored variant, reduced to what the tranche arithmetic reads off it.
#[derive(Debug, Clone, PartialEq)]
pub struct Datum {
    pub lod: f64,
    /// Present in a resource tagged `known=true`.
    pub is_known: bool,
    /// Present in a resource tagged `truth=true`.
    pub at_truth_site: bool,
    pub is_snp: bool,
    pub is_transition: bool,
}

/// `SNP` or `INDEL`, which is only ever written into the filter name and the model column.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    Snp,
    Indel,
}

impl Mode {
    pub fn name(self) -> &'static str {
        match self {
            Mode::Snp => "SNP",
            Mode::Indel => "INDEL",
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
pub fn running_sensitivity(data: &[Datum], n_true_sites: usize) -> Vec<f64> {
    let mut running = vec![0.0; data.len()];
    let mut called_at_truth = 0;
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
    pub num_known: i64,
    pub num_novel: i64,
    pub known_ti_tv: f64,
    pub novel_ti_tv: f64,
    pub min_vqs_lod: f64,
    pub model: Mode,
    pub accessible_truth_sites: i32,
    pub calls_at_truth_sites: i32,
}

impl Tranche {
    /// `getTruthSensitivity`, which is zero rather than a division when there is no truth site.
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
    let index = (0..data.len()).find(|i| running[*i] >= threshold)?;
    Some(tranche_of_variants(data, index, target, model))
}

/// `findVQSLODTranche`: the first index whose LOD reaches the requested score.
///
/// Two things separate it from the sensitivity search. The minimum it reports is the REQUEST
/// rather than the LOD it found, and a request no variant reaches still produces a row: an empty
/// tranche over the last variant, again with the request as its minimum.
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
pub fn no_tranche_refusal(metric: &str, threshold: f64) -> String {
    format!(
        "Couldn't find any tranche containing variants with a {metric} > {threshold:.2}. Are you \
         sure the truth files contain unfiltered variants which overlap the input data?"
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
pub const TRUTH_SENSITIVITY_VERSION: i32 = 5;
pub const VQSLOD_VERSION: i32 = 6;

pub fn truth_sensitivity_header() -> String {
    format!(
        "# Variant quality score tranches file\n# Version number {TRUTH_SENSITIVITY_VERSION}\n\
         targetTruthSensitivity,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,\
         model,accessibleTruthSites,callsAtTruthSites,truthSensitivity\n"
    )
}

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
pub fn tranche_order(a: &Tranche, b: &Tranche) -> std::cmp::Ordering {
    compare_double(a.calls_at_truth_sites as f64, b.calls_at_truth_sites as f64)
}

/// `getTrancheString`: one row, whose filter name takes its LOWER bound from the previous row.
///
/// The bound is the previous row's own index whatever that index is, so a list that is not in
/// increasing order produces a band that runs backwards. The first row has no previous row and
/// so is bounded at 0.00.
pub fn tranche_row(tranche: &Tranche, previous: Option<&Tranche>) -> String {
    format!(
        "{:.2},{},{},{:.4},{:.4},{:.4},VQSRTranche{}{:.2}to{:.2},{},{},{},{:.4}\n",
        tranche.index,
        tranche.num_known,
        tranche.num_novel,
        tranche.known_ti_tv,
        tranche.novel_ti_tv,
        tranche.min_vqs_lod,
        tranche.model.name(),
        previous.map_or(0.0, |previous| previous.index),
        tranche.index,
        tranche.model.name(),
        tranche.accessible_truth_sites,
        tranche.calls_at_truth_sites,
        tranche.truth_sensitivity()
    )
}

/// `tranchesString`: the rows, sorted, each naming the one before it.
///
/// The sort is STABLE, so two targets that found the same tranche keep the order they were given
/// in, and a list of one is not sorted at all.
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

/// The whole truth-sensitivity file.
pub fn truth_sensitivity_file(tranches: &[Tranche]) -> String {
    truth_sensitivity_header() + &tranches_string(tranches)
}

/// The whole VQSLOD file.
pub fn vqslod_file(tranches: &[Tranche]) -> String {
    vqslod_header() + &tranches_string(tranches)
}

//! `AddFlowSNVQuality`, ported from the tool (GATK 4.6.2.0).
//!
//! The sibling of [`crate::add_flow_base_quality`]: the same enumeration over the ways a flow key
//! could have been misread, but producing a probability for each base that was NOT called. Building
//! the flow matrix is not ported; everything the tool does with it is.
//!
//! # The computed base quality is discarded
//!
//! ```java
//! // at this point, bq becomes trivial (?)
//! baseProbs[ofs] = 1 - snvqProbs[calledIndex][ofs];
//! ```
//!
//! The hmer loop fills `baseProbs` with the same per-hmer work the sibling does, and then the
//! normalisation loop overwrites EVERY entry with one minus the called base's own probability,
//! which is itself one minus the sum of the alternates. The reference's own comment is the
//! question mark.
//!
//! # The phred conversion rounds, where the sibling truncates
//!
//! ```java
//! phred[i] = (byte)Math.round(-10 * Math.log10(p));
//! ```
//!
//! `AddFlowBaseQuality` casts the same expression instead. Two tools, one formula, two answers.
//!
//! # The side walk is bounded by the slice
//!
//! ```java
//! if ( sideFlow < minIndex || sideFlow > maxIndex ) { break; }
//! ```
//!
//! That bound is exactly what the sibling lacks, so a flow order whose cycle is ONE gets past the
//! enumeration here. It dies in the normalisation instead: with a cycle of one only the first base
//! of the order is considered, a read carrying any other leaves `calledIndex` at -1, and
//! `snvqProbs[-1]` is indexed with it. The guard above tests `calledBase < 0`, which an ASCII base
//! never is, so it never fires. [`SnvError`] is that refusal.
//!
//! # The alternates are keyed by base, not by flow
//!
//! `allBaseProb` is a `LinkedHashMap<Byte, Double>`, so two side flows carrying the same base
//! overwrite each other and the LAST one wins.
//!
//! # And --max-phred-score moves the floor as well as the clamp
//!
//! ```java
//! maxQualityScore = (int)maxPhredScore;
//! minLikelihoodProbRate = Math.pow(10, -maxPhredScore / 10.0);
//! ```
//!
//! At 10 the floor is 0.1, so three alternates take 0.3 of the mass and every base comes out at 5:
//! the clamp is never reached.

use crate::add_flow_base_quality::{
    extract_error_prob_bands, RawProbs, ERROR_PROB_BAND_1LESS, ERROR_PROB_BAND_1MORE,
    ERROR_PROB_BAND_KEY, PHRED_ASCII_BASE,
};

/// `AddFlowSNVQuality.minLikelihoodProbRate`.
pub const DEFAULT_MIN_LIKELIHOOD_PROB_RATE: f64 = 1e-6;
/// `AddFlowSNVQuality.maxQualityScore`.
pub const DEFAULT_MAX_QUALITY_SCORE: i32 = 60;

/// `AddFlowSNVQualityArgumentCollection.SnvqModeEnum`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SnvqMode {
    Legacy,
    Optimistic,
    Pessimistic,
    Geometric,
}

/// `getSnvq`: four formulae over the slice's probability and the two flows it touched.
pub fn get_snvq(slice_p: f64, p1: f64, p2: f64, mode: SnvqMode) -> f64 {
    match mode {
        SnvqMode::Legacy => slice_p,
        SnvqMode::Optimistic => p1 * p2,
        SnvqMode::Pessimistic => 1.0 - (1.0 - p1) * (1.0 - p2),
        SnvqMode::Geometric => ((p1 * p2) * (1.0 - (1.0 - p1) * (1.0 - p2))).sqrt(),
    }
}

/// What the tool throws rather than reports.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SnvError {
    /// `snvqProbs[calledIndex]` with `calledIndex` still -1, which a base outside the flow order's
    /// cycle reaches.
    CalledBaseNotInOrder { index: i64, length: usize },
    /// `sliceProbs` given a slice more than one away from the key.
    SliceTooFar { slice: i32, hmer: i32 },
}

impl SnvError {
    pub fn java_class(&self) -> &'static str {
        match self {
            SnvError::CalledBaseNotInOrder { .. } => "java.lang.ArrayIndexOutOfBoundsException",
            SnvError::SliceTooFar { .. } => {
                "org.broadinstitute.hellbender.exceptions.GATKException"
            }
        }
    }

    pub fn message(&self) -> String {
        match self {
            SnvError::CalledBaseNotInOrder { index, length } => {
                format!("Index {index} out of bounds for length {length}")
            }
            SnvError::SliceTooFar { slice, hmer } => {
                format!("slice[i] and hmer are too far apart: {slice} {hmer}")
            }
        }
    }
}

/// `sliceProbs`, which returns the product AND the two probabilities the snvq modes combine.
pub fn slice_probs(
    slice: &[i32],
    min_index: usize,
    key: &[i32],
    bands: &[Vec<f64>],
    flow: usize,
    side_flow: usize,
) -> Result<(f64, f64, f64), SnvError> {
    let mut accumulated = 1.0;
    let mut p1 = 0.0;
    let mut p2 = 0.0;
    for (offset, value) in slice.iter().enumerate() {
        let key_index = min_index + offset;
        let hmer = key[key_index];
        let band = if *value == hmer - 1 {
            ERROR_PROB_BAND_1LESS
        } else if *value == hmer + 1 {
            ERROR_PROB_BAND_1MORE
        } else if *value == hmer {
            ERROR_PROB_BAND_KEY
        } else {
            return Err(SnvError::SliceTooFar {
                slice: *value,
                hmer,
            });
        };
        let probability = bands[band][key_index];
        accumulated *= probability;
        if key_index == flow {
            p1 = probability;
        }
        if key_index == side_flow {
            p2 = probability;
        }
    }
    Ok((accumulated, p1, p2))
}

/// `sliceIsValidForConsideration`, which is the sibling's rule unchanged.
pub fn slice_is_valid(slice: &[i32], flow_order_length: usize) -> bool {
    let mut consecutive_zeros = 0usize;
    for value in slice {
        if *value != 0 {
            consecutive_zeros = 0;
        } else {
            consecutive_zeros += 1;
            if consecutive_zeros + 1 >= flow_order_length {
                return false;
            }
        }
    }
    true
}

/// The alternate probabilities one side walk collected, keyed by BASE.
type AltProbs = Vec<(u8, f64)>;

fn put(alternates: &mut AltProbs, base: u8, probability: f64) {
    match alternates
        .iter_mut()
        .find(|(existing, _)| *existing == base)
    {
        // A LinkedHashMap keeps the first insertion's POSITION and takes the last value.
        Some((_, value)) => *value = probability,
        None => alternates.push((base, probability)),
    }
}

/// `generateSidedHmerBaseErrorProbability`, which also fills the alternate probabilities.
///
/// The eight parameters are the reference's own: splitting them into a struct would hide which of
/// them the walk reads and which it writes.
#[allow(clippy::too_many_arguments)]
fn generate_sided(
    key: &[i32],
    bands: &[Vec<f64>],
    flow: usize,
    side_incr: i32,
    flow_order_length: usize,
    flow_order: &[u8],
    mode: SnvqMode,
    alternates: &mut AltProbs,
) -> Result<f64, SnvError> {
    let min_index = (flow as i64 - (flow_order_length as i64 - 1)).max(0) as usize;
    let max_index = ((flow + flow_order_length - 1) as i64).min(key.len() as i64 - 1) as usize;
    let slice: Vec<i32> = key[min_index..=max_index].to_vec();
    let hmer_length = key[flow];

    let increments: Vec<i32> = if hmer_length != 1 {
        vec![side_incr]
    } else {
        vec![side_incr, -side_incr]
    };
    let mut slices: Vec<(Vec<i32>, u8, usize)> = Vec::new();
    for increment in increments {
        let mut side_flow = flow as i64 + increment as i64;
        while side_flow >= 0 && side_flow < key.len() as i64 {
            // The bound the sibling tool does not have: the walk stops at the slice's edge rather
            // than indexing past it.
            if side_flow < min_index as i64 || side_flow > max_index as i64 {
                break;
            }
            let mut alt = slice.clone();
            alt[(side_flow - min_index as i64) as usize] += 1;
            alt[flow - min_index] -= 1;
            if slice_is_valid(&alt, flow_order_length) {
                slices.push((
                    alt,
                    flow_order[side_flow as usize % flow_order_length],
                    side_flow as usize,
                ));
            }
            if key[side_flow as usize] != 0 {
                break;
            }
            side_flow += increment as i64;
        }
    }

    let (key_probability, _, _) = slice_probs(&slice, min_index, key, bands, flow, flow)?;
    let mut sum = key_probability;
    for (alt, base, side_flow) in &slices {
        let (probability, p1, p2) = slice_probs(alt, min_index, key, bands, flow, *side_flow)?;
        put(alternates, *base, get_snvq(probability, p1, p2, mode));
        sum += probability;
    }
    Ok(1.0 - (key_probability / sum))
}

/// The two arrays one read comes out with.
#[derive(Debug, Clone, PartialEq)]
pub struct ReadProbs {
    pub base_probs: Vec<f64>,
    /// One row per base of the flow order's cycle, each as long as the read.
    pub snvq_probs: Vec<Vec<f64>>,
}

/// `generateFlowReadBaseAndSNVQErrorProbabilities`, including the normalisation that throws the
/// base probabilities away.
pub fn generate_read_probs(
    key: &[i32],
    bands: &[Vec<f64>],
    bases: &[u8],
    flow_order: &[u8],
    flow_order_length: usize,
    min_likelihood_prob_rate: f64,
    mode: SnvqMode,
) -> Result<ReadProbs, SnvError> {
    let mut base_probs = vec![0.0; bases.len()];
    let mut snvq_probs = vec![vec![0.0; bases.len()]; flow_order_length];

    let mut base = 0usize;
    for flow in 0..key.len() {
        if key[flow] == 0 {
            continue;
        }
        let mut first: AltProbs = Vec::new();
        let mut last: AltProbs = Vec::new();
        let flow_i = flow % flow_order_length;
        let hmer_length = key[flow] as usize;

        let left = generate_sided(
            key,
            bands,
            flow,
            -1,
            flow_order_length,
            flow_order,
            mode,
            &mut first,
        )?;
        let right = if key[flow] != 1 {
            generate_sided(
                key,
                bands,
                flow,
                1,
                flow_order_length,
                flow_order,
                mode,
                &mut last,
            )?
        } else {
            0.0
        };

        base_probs[base] = left;
        base += 1;
        for (index, order_base) in flow_order.iter().take(flow_order_length).enumerate() {
            if let Some((_, probability)) = first.iter().find(|(b, _)| b == order_base) {
                snvq_probs[index][base - 1] = *probability;
            } else if index != flow_i {
                snvq_probs[index][base - 1] = min_likelihood_prob_rate;
            }
        }

        if hmer_length > 1 {
            // The middle bases are stepped over, exactly as in the sibling.
            base += hmer_length - 2;
            base_probs[base] = right;
            base += 1;
            for (index, order_base) in flow_order.iter().take(flow_order_length).enumerate() {
                if let Some((_, probability)) = last.iter().find(|(b, _)| b == order_base) {
                    for j in 0..hmer_length - 1 {
                        snvq_probs[index][base - 1 - j] = if j == 0 {
                            *probability
                        } else {
                            min_likelihood_prob_rate
                        };
                    }
                } else if index != flow_i {
                    for j in 0..hmer_length - 1 {
                        snvq_probs[index][base - 1 - j] = min_likelihood_prob_rate;
                    }
                }
            }
        }

        if base == base_probs.len() {
            base_probs[base - 1] = bands[ERROR_PROB_BAND_KEY][flow];
        }
    }

    // And then every base probability computed above is thrown away.
    for (offset, called) in bases.iter().enumerate() {
        let mut alt_p = 0.0;
        let mut called_index: i64 = -1;
        for (index, order_base) in flow_order.iter().take(flow_order_length).enumerate() {
            if called != order_base {
                snvq_probs[index][offset] = snvq_probs[index][offset].max(min_likelihood_prob_rate);
                alt_p += snvq_probs[index][offset];
            } else {
                called_index = index as i64;
            }
        }
        if called_index < 0 {
            // The reference indexes the array with -1 here; the guard above it tests the BASE
            // rather than the index and never fires.
            return Err(SnvError::CalledBaseNotInOrder {
                index: called_index,
                length: flow_order_length,
            });
        }
        snvq_probs[called_index as usize][offset] = (1.0 - alt_p).max(0.0);
        base_probs[offset] = 1.0 - snvq_probs[called_index as usize][offset];
    }

    Ok(ReadProbs {
        base_probs,
        snvq_probs,
    })
}

/// `convertErrorProbToPhred`: a ROUNDING, where the sibling truncates.
pub fn convert_error_prob_to_phred(error_prob: &[f64], max_quality_score: i32) -> Vec<u8> {
    error_prob
        .iter()
        .map(|probability| {
            if *probability == 0.0 {
                max_quality_score as u8
            } else {
                // `Math.round(double)` is floor(x + 0.5), not a half-even rounding.
                ((-10.0 * probability.log10()) + 0.5).floor() as i64 as u8
            }
        })
        .collect()
}

/// `SAMUtils.phredToFastq`.
pub fn phred_to_fastq(phred: &[u8]) -> String {
    phred
        .iter()
        .map(|value| value.wrapping_add(PHRED_ASCII_BASE) as char)
        .collect()
}

/// What one read carries after the run.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadOutput {
    /// The `QUAL` field, unchanged when the base quality went to a tag instead.
    pub qualities: String,
    /// The tag `--output-quality-attribute` named, when it was given.
    pub quality_attribute: Option<String>,
    /// One per base of the flow order's cycle, in that order, named `q` plus the lower-cased base.
    pub snvq_attributes: Vec<(String, String)>,
}

/// `attrNameForNonCalledBase`.
pub fn attr_name_for_non_called_base(base: u8) -> String {
    format!("q{}", (base as char).to_ascii_lowercase())
}

/// `addBaseQuality`, over one read whose key and raw probabilities are already known.
#[allow(clippy::too_many_arguments)]
pub fn add_base_quality(
    key: &[i32],
    probs: &[RawProbs],
    bases: &[u8],
    original_qualities: &str,
    flow_order: &str,
    flow_order_length: usize,
    max_phred_score: Option<f64>,
    mode: SnvqMode,
    output_quality_attribute: Option<&str>,
) -> Result<ReadOutput, SnvError> {
    // The one argument that moves two numbers.
    let (max_quality_score, min_rate) = match max_phred_score {
        Some(score) => (score as i32, 10f64.powf(-score / 10.0)),
        None => (DEFAULT_MAX_QUALITY_SCORE, DEFAULT_MIN_LIKELIHOOD_PROB_RATE),
    };
    let bands = extract_error_prob_bands(probs, min_rate);
    let order = flow_order.as_bytes();
    let read_probs =
        generate_read_probs(key, &bands, bases, order, flow_order_length, min_rate, mode)?;

    let base_quality = phred_to_fastq(&convert_error_prob_to_phred(
        &read_probs.base_probs,
        max_quality_score,
    ));
    let snvq_attributes = read_probs
        .snvq_probs
        .iter()
        .enumerate()
        .map(|(index, row)| {
            (
                attr_name_for_non_called_base(order[index]),
                phred_to_fastq(&convert_error_prob_to_phred(row, max_quality_score)),
            )
        })
        .collect();

    Ok(match output_quality_attribute {
        Some(_) => ReadOutput {
            qualities: original_qualities.to_string(),
            quality_attribute: Some(base_quality),
            snvq_attributes,
        },
        None => ReadOutput {
            qualities: base_quality,
            quality_attribute: None,
            snvq_attributes,
        },
    })
}

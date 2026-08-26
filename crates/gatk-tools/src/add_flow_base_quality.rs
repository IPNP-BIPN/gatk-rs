//! `AddFlowBaseQuality`, ported from the tool (GATK 4.6.2.0).
//!
//! A flow-based read carries the probability that each homopolymer was read one base too short or
//! too long. This turns those into a per-base quality by enumerating the ways the flow key could
//! have been misread. Building the flow matrix is not ported; everything the tool does with it is.
//!
//! # The middle bases of an hmer are never computed
//!
//! ```java
//! result[base++] = hmerBaseErrorProbs[0];
//! if ( hmerLength > 1 ) {
//!     base += (hmerLength - 2);          // skip all but last, leaving zero
//!     result[base++] = hmerBaseErrorProbs[1];
//! }
//! ```
//!
//! The cursor jumps over them and the array keeps the zero it was allocated with. A zero error
//! probability is not a low quality: `convertErrorProbToPhred` special-cases it to the MAXIMAL
//! quality. So `TTTTTGCA` comes out as `!~~~\YY!`, three characters of quality 93 sitting in the
//! middle of the homopolymer.
//!
//! # Both ends of the read are overridden
//!
//! The first base takes the hmer's own key probability rather than the computed one, by a test on
//! the cursor being zero. The last base takes it too, by a second override that fires whenever the
//! cursor has reached the end of the array, which overwrites a value the hmer loop had just
//! written.
//!
//! # An hmer of length one walks both sides
//!
//! ```java
//! final int[] incrs = (hmerLength != 1) ? new int[] { sideIncr } : new int[] { sideIncr, -sideIncr};
//! ```
//!
//! It is also the only hmer whose second probability is never computed: the right-hand call is
//! skipped, so `errorProbs[1]` stays zero and is never read, because an hmer of one has no last
//! base of its own.
//!
//! # The slice window is the flow order's cycle, and a cycle of one throws
//!
//! `calcFlowOrderLength` answers the distance to the SECOND occurrence of the order's first base,
//! so `TGCA` gives four and `TTGCA` gives ONE. With a cycle of one the slice is a single flow, and
//! `altSlice[sideFlow - minIndex] += 1` indexes one past it: measured as
//! `ArrayIndexOutOfBoundsException: Index 1 out of bounds for length 1`, thrown out of
//! `generateSidedHmerBaseErrorProbability`. [`FlowError`] is that refusal rather than a panic.
//!
//! # The floor is applied to the impossible bands too
//!
//! A key of 0 has no shorter neighbour and a key at `maxHmer` has no longer one. Both get
//! `--minimal-error-rate` rather than zero, which is what keeps `sliceProb` from multiplying a
//! whole product away.

/// `ERROR_PROB_BAND_1LESS`.
pub const ERROR_PROB_BAND_1LESS: usize = 0;
/// `ERROR_PROB_BAND_KEY`.
pub const ERROR_PROB_BAND_KEY: usize = 1;
/// `ERROR_PROB_BAND_1MORE`.
pub const ERROR_PROB_BAND_1MORE: usize = 2;
/// `PHRED_ASCII_BASE`.
pub const PHRED_ASCII_BASE: u8 = b'!';
/// The tool's defaults.
pub const DEFAULT_MIN_ERROR_RATE: f64 = 1e-3;
pub const DEFAULT_MAX_QUALITY_SCORE: i32 = 93;

/// `java.lang.ArrayIndexOutOfBoundsException`, which a flow order of cycle one reaches.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FlowError {
    pub index: i64,
    pub length: usize,
}

impl FlowError {
    pub fn java_class(&self) -> &'static str {
        "java.lang.ArrayIndexOutOfBoundsException"
    }

    pub fn message(&self) -> String {
        format!(
            "Index {} out of bounds for length {}",
            self.index, self.length
        )
    }
}

/// `calcFlowOrderLength`: the distance to the second occurrence of the order's first base, or the
/// whole order when there is none.
pub fn calc_flow_order_length(flow_order: &str) -> usize {
    let bytes = flow_order.as_bytes();
    match bytes.iter().skip(1).position(|base| *base == bytes[0]) {
        Some(offset) => offset + 1,
        None => bytes.len(),
    }
}

/// The three probabilities the analyzer reads off the flow matrix for one flow, before the floor.
///
/// `None` is the neighbour that does not exist: a key of 0 has nothing shorter, a key at `maxHmer`
/// nothing longer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RawProbs {
    pub minus: Option<f64>,
    pub key: f64,
    pub plus: Option<f64>,
}

/// `extractErrorProbBands`, which floors every band including the ones that cannot happen.
pub fn extract_error_prob_bands(probs: &[RawProbs], min_value: f64) -> Vec<Vec<f64>> {
    let mut bands = vec![vec![0.0; probs.len()]; 3];
    for (flow, raw) in probs.iter().enumerate() {
        bands[ERROR_PROB_BAND_KEY][flow] = raw.key.max(min_value);
        bands[ERROR_PROB_BAND_1LESS][flow] = match raw.minus {
            Some(value) => value.max(min_value),
            None => min_value,
        };
        bands[ERROR_PROB_BAND_1MORE][flow] = match raw.plus {
            Some(value) => value.max(min_value),
            None => min_value,
        };
    }
    bands
}

/// `sliceIsValid`: `flowOrderLength - 1` consecutive zeros invalidates a slice, so a cycle of one
/// makes every zero fatal and a cycle of four allows at most two.
pub fn slice_is_valid(slice: &[i32], flow_order_length: usize) -> bool {
    let mut consecutive_zeros = 0usize;
    for key in slice {
        if *key != 0 {
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

/// `sliceProb`: the product of one band per flow, chosen by how the slice differs from the key.
pub fn slice_prob(slice: &[i32], min_index: usize, key: &[i32], bands: &[Vec<f64>]) -> f64 {
    let mut product = 1.0;
    for (offset, value) in slice.iter().enumerate() {
        let band = match (*value).cmp(&key[offset + min_index]) {
            std::cmp::Ordering::Less => ERROR_PROB_BAND_1LESS,
            std::cmp::Ordering::Greater => ERROR_PROB_BAND_1MORE,
            std::cmp::Ordering::Equal => ERROR_PROB_BAND_KEY,
        };
        product *= bands[band][offset + min_index];
    }
    product
}

/// `generateSidedHmerBaseErrorProbability`.
fn generate_sided(
    key: &[i32],
    bands: &[Vec<f64>],
    flow: usize,
    side_incr: i32,
    flow_order_length: usize,
) -> Result<f64, FlowError> {
    let min_index = (flow as i64 - flow_order_length as i64 + 1).max(0) as usize;
    let max_index = ((flow + flow_order_length - 1) as i64).min(key.len() as i64 - 1) as usize;
    let slice: Vec<i32> = key[min_index..=max_index].to_vec();
    let hmer_length = key[flow];

    // An hmer of length one is the only one that walks in both directions.
    let increments: Vec<i32> = if hmer_length != 1 {
        vec![side_incr]
    } else {
        vec![side_incr, -side_incr]
    };
    let mut slices: Vec<Vec<i32>> = Vec::new();
    for increment in increments {
        let mut side_flow = flow as i64 + increment as i64;
        while side_flow >= 0 && side_flow < key.len() as i64 {
            let mut alt = slice.clone();
            let index = side_flow - min_index as i64;
            // The reference indexes the slice here with no bound of its own, and a cycle of one
            // makes this the index that leaves the array.
            if index < 0 || index as usize >= alt.len() {
                return Err(FlowError {
                    index,
                    length: alt.len(),
                });
            }
            alt[index as usize] += 1;
            alt[flow - min_index] -= 1;
            if slice_is_valid(&alt, flow_order_length) {
                slices.push(alt);
            }
            if key[side_flow as usize] != 0 {
                break;
            }
            side_flow += increment as i64;
        }
    }

    let key_probability = slice_prob(&slice, min_index, key, bands);
    let mut sum = key_probability;
    for alt in &slices {
        sum += slice_prob(alt, min_index, key, bands);
    }
    Ok(1.0 - (key_probability / sum))
}

/// `generateHmerBaseErrorProbabilities`: one probability for the hmer's first base, one for the
/// rest, and the second left at zero for an hmer of length one.
pub fn generate_hmer_base_error_probabilities(
    key: &[i32],
    bands: &[Vec<f64>],
    flow: usize,
    flow_order_length: usize,
) -> Result<[f64; 2], FlowError> {
    let mut probabilities = [0.0; 2];
    probabilities[0] = generate_sided(key, bands, flow, -1, flow_order_length)?;
    if key[flow] != 1 {
        probabilities[1] = generate_sided(key, bands, flow, 1, flow_order_length)?;
    }
    Ok(probabilities)
}

/// `generateBaseErrorProbability`, over a read of `read_length` bases.
///
/// The entries this never writes stay at zero, which is the whole of the middle-of-the-hmer
/// behaviour.
pub fn generate_base_error_probability(
    key: &[i32],
    bands: &[Vec<f64>],
    read_length: usize,
    flow_order_length: usize,
) -> Result<Vec<f64>, FlowError> {
    let mut result = vec![0.0; read_length];
    let mut base = 0usize;
    for flow in 0..key.len() {
        if key[flow] == 0 {
            continue;
        }
        let hmer_length = key[flow] as usize;
        let hmer_probs =
            generate_hmer_base_error_probabilities(key, bands, flow, flow_order_length)?;

        // The first base of the READ takes the hmer's own key probability, not the computed one.
        result[base] = if base == 0 {
            bands[ERROR_PROB_BAND_KEY][flow]
        } else {
            hmer_probs[0]
        };
        base += 1;

        if hmer_length > 1 {
            base += hmer_length - 2;
            result[base] = hmer_probs[1];
            base += 1;
        }

        // And the last base of the READ takes it too, overwriting whatever was just written.
        if base == result.len() {
            result[base - 1] = bands[ERROR_PROB_BAND_KEY][flow];
        }
    }
    Ok(result)
}

/// `convertErrorProbToPhred`: a truncation, not a rounding, and a probability of exactly zero
/// takes the clamp rather than an infinity.
pub fn convert_error_prob_to_phred(error_prob: &[f64], max_quality_score: i32) -> Vec<u8> {
    error_prob
        .iter()
        .map(|probability| {
            if *probability == 0.0 {
                max_quality_score as u8
            } else {
                std::cmp::min(max_quality_score, (-10.0 * probability.log10()) as i32) as u8
            }
        })
        .collect()
}

/// `convertPhredToString`, which is the SAM quality encoding.
pub fn convert_phred_to_string(phred: &[u8]) -> String {
    phred
        .iter()
        .map(|value| (value.wrapping_add(PHRED_ASCII_BASE)) as char)
        .collect()
}

/// What one read comes out as: the new quality string, and whether it replaces `QUAL` or lands in
/// `XQ`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadOutput {
    /// `XQ`, when the tool is not replacing.
    pub base_quality_attribute: Option<String>,
    /// `OQ`, when it is.
    pub old_quality_attribute: Option<String>,
    /// The quality string the record carries afterwards.
    pub qualities: String,
}

/// `addBaseQuality`, over one read whose key and bands are already known.
#[allow(clippy::too_many_arguments)]
pub fn add_base_quality(
    key: &[i32],
    probs: &[RawProbs],
    read_length: usize,
    original_qualities: &str,
    flow_order: &str,
    min_error_rate: f64,
    max_quality_score: i32,
    replace_quality_mode: bool,
) -> Result<ReadOutput, FlowError> {
    let bands = extract_error_prob_bands(probs, min_error_rate);
    let flow_order_length = calc_flow_order_length(flow_order);
    let error_prob = generate_base_error_probability(key, &bands, read_length, flow_order_length)?;
    let phred = convert_error_prob_to_phred(&error_prob, max_quality_score);
    let written = convert_phred_to_string(&phred);
    Ok(if replace_quality_mode {
        ReadOutput {
            base_quality_attribute: None,
            old_quality_attribute: Some(original_qualities.to_string()),
            qualities: written,
        }
    } else {
        ReadOutput {
            base_quality_attribute: Some(written),
            old_quality_attribute: None,
            qualities: original_qualities.to_string(),
        }
    })
}

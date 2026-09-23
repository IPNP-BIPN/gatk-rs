//! Ported from `org.broadinstitute.hellbender.utils.read.FlowBasedRead` (GATK 4.6.2.0): the flow
//! matrix a flow-based read carries, built from its `tp` and `t0` tags.
//!
//! A flow read's bases are a flow KEY in disguise: the flow order `TGCA` is cycled, and each flow
//! reports how many copies of its base were read, zero included. The matrix has one column per flow
//! and one row per hmer length from zero to `maxHmer`, and a cell is the probability that the flow
//! really held that many bases.
//!
//! # How the `tp` tag fills it
//!
//! ```java
//! final int loc = Math.max(Math.min(flowCall+tp[i], maxHmer),0);
//! if (flowMatrix[loc][flowIdx] == perHmerMinErrorProb) flowMatrix[loc][flowIdx] = probs[i];
//! else                                                 flowMatrix[loc][flowIdx] += probs[i];
//! ```
//!
//! Each base of an hmer carries a `tp` offset and its quality: the quality is the probability of
//! the hmer being `call + tp` long instead. A cell still at the filling value is REPLACED, and one
//! already written is ADDED to, so two bases naming the same length sum. The called cell is then
//! set to `1 - sum`, over the rows below `maxHmer` only, and never lower than 0.1.
//!
//! # The filling value of zero is a different mode
//!
//! `--flow-fill-empty-bins-value 0` does not fill with zero: it estimates a floor from the read's
//! own best quality and divides it by `maxHmer`, which is why the empty-cell probability of the
//! `t0` test and the clipping threshold of `clipProbs` are two different numbers.
//!
//! # Only the production format is ported
//!
//! The vestigial `kr` and `ti` layouts are refused as a port limitation, and so is a read group
//! that is not flow-based, which the reference turns into a `NullPointerException`.

use htsjdk_bam::cigar::Op;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

/// `FlowBasedRead.MAX_CLASS`, the maximal hmer when the read group has no `mc`.
pub const MAX_CLASS: i32 = 12;
/// `MINIMAL_CALL_PROB`, the floor of the called cell.
const MINIMAL_CALL_PROB: f64 = 0.1;
/// `FLOW_MATRIX_TAG_NAME`.
pub const FLOW_MATRIX_TAG: &[u8; 2] = b"tp";
/// `FLOW_MATRIX_T0_TAG_NAME`.
pub const FLOW_MATRIX_T0_TAG: &[u8; 2] = b"t0";
/// The two vestigial layouts, `FLOW_MATRiX_OLD_TAG_KR` and `..._TI`.
const OLD_TAG_KR: &[u8; 2] = b"kr";
const OLD_TAG_TI: &[u8; 2] = b"ti";

/// `FlowBasedArgumentCollection`, the arguments that shape the matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct FlowArguments {
    pub use_t0_tag: bool,
    pub remove_longer_than_one_indels: bool,
    pub remove_one_to_zero_probs: bool,
    pub filling_value: f64,
    pub symmetric_indels: bool,
    pub only_ins_or_del: bool,
    pub disallow_larger_probs: bool,
    pub lump_probs: bool,
    pub retain_max_n_probs: bool,
    pub flow_matrix_mods: Option<String>,
    pub keep_boundary_flows: bool,
}

impl Default for FlowArguments {
    fn default() -> Self {
        FlowArguments {
            use_t0_tag: false,
            remove_longer_than_one_indels: false,
            remove_one_to_zero_probs: false,
            filling_value: 0.001,
            symmetric_indels: false,
            only_ins_or_del: false,
            disallow_larger_probs: false,
            lump_probs: false,
            retain_max_n_probs: false,
            flow_matrix_mods: None,
            keep_boundary_flows: false,
        }
    }
}

/// Why a flow read could not be built: the reference's exception, by class and message.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FlowReadError {
    pub class: &'static str,
    pub message: String,
    /// The reference would have gone on, into a layout the port does not read.
    pub port_limitation: bool,
}

impl FlowReadError {
    fn thrown(class: &'static str, message: String) -> Self {
        FlowReadError {
            class,
            message,
            port_limitation: false,
        }
    }

    fn out_of_bounds(index: i64, length: usize) -> Self {
        Self::thrown(
            "java.lang.ArrayIndexOutOfBoundsException",
            format!("Index {index} out of bounds for length {length}"),
        )
    }

    fn limitation(message: String) -> Self {
        FlowReadError {
            class: "",
            message,
            port_limitation: true,
        }
    }
}

/// `FlowBasedReadUtils.ReadGroupInfo`, for a flow read group: its flow order and maximal class.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadGroupInfo {
    pub flow_order: String,
    pub max_class: i32,
}

/// `NGSPlatform.fromReadGroupPL`, reduced to the two platforms that make a group flow-based.
fn platform(pl: &str) -> Option<&'static str> {
    match pl.to_ascii_uppercase().as_str() {
        "ULTIMA" => Some("ULTIMA"),
        "454" | "LS454" => Some("LS454"),
        _ => None,
    }
}

/// `FlowBasedReadUtils.getReadGroupInfo`, for a read already known to carry flow tags.
///
/// A group that is not flow-based yields a null flow order, which the matrix construction
/// dereferences; the port refuses it instead.
pub fn read_group_info(
    read: &BamRecord,
    header: &SamHeader,
) -> Result<ReadGroupInfo, FlowReadError> {
    let Some(group) = gatk_engine::read_group::resolve(read, header) else {
        return Err(FlowReadError::limitation(
            "a flow read whose read group the header does not declare".to_string(),
        ));
    };
    let flow_order = group.attributes.get("FO");
    let platform = group.attributes.get("PL").and_then(platform);
    match (platform, flow_order) {
        (Some("ULTIMA"), None) => Err(FlowReadError::thrown(
            "java.lang.RuntimeException",
            format!(
                "Malformed Ultima read group identified, aborting: SAMReadGroupRecord{{ID: {}}}",
                group.id
            ),
        )),
        (Some(_), Some(order)) => {
            let max_class = match group.attributes.get("mc") {
                None => MAX_CLASS,
                Some(text) => text.parse::<i32>().map_err(|_| {
                    FlowReadError::thrown(
                        "java.lang.NumberFormatException",
                        format!("For input string: \"{text}\""),
                    )
                })?,
            };
            Ok(ReadGroupInfo {
                flow_order: order.to_string(),
                max_class,
            })
        }
        _ => Err(FlowReadError::limitation(format!(
            "read group {} is not a flow-based read group with a flow order",
            group.id
        ))),
    }
}

/// `FlowBasedReadUtils.hasFlowTags`.
pub fn has_flow_tags(read: &BamRecord) -> bool {
    [FLOW_MATRIX_TAG, OLD_TAG_KR, OLD_TAG_TI]
        .iter()
        .any(|name| read.tags.get(Tag::new(name)).is_some())
}

/// `FlowBasedKeyCodec.baseArrayToKey`, `None` where the period guard trips.
pub fn base_array_to_key(bases: &[u8], flow_order: &str) -> Option<Vec<i32>> {
    let flow = flow_order.as_bytes();
    let period = flow.len();
    let mut result = Vec::new();
    let mut loc = 0usize;
    let mut flow_number = 0usize;
    let mut period_guard = 0usize;
    while loc < bases.len() {
        let flow_base = flow[flow_number % period];
        if bases[loc] != flow_base && bases[loc] != b'N' {
            result.push(0);
            period_guard += 1;
            if period_guard > period {
                return None;
            }
        } else {
            let mut count = 0;
            while loc < bases.len() && (bases[loc] == flow_base || bases[loc] == b'N') {
                loc += 1;
                count += 1;
            }
            result.push(count);
            period_guard = 0;
        }
        flow_number += 1;
    }
    Some(result)
}

/// A flow read's key and matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct FlowRead {
    pub key: Vec<i32>,
    pub max_hmer: i32,
    /// `flowMatrix[hmer][flow]`, `maxHmer + 1` rows of `key.len()` columns.
    pub matrix: Vec<Vec<f64>>,
    per_hmer_min_error_prob: f64,
}

impl FlowRead {
    /// `getProb`: the cell for an hmer of that length, capped at one.
    pub fn prob(&self, flow: usize, hmer: i32) -> f64 {
        let prob = self.matrix[hmer.min(self.max_hmer) as usize][flow];
        if prob <= 1.0 {
            prob
        } else {
            1.0
        }
    }

    /// `new FlowBasedRead(read, flowOrder, maxHmer, fbargs)`.
    pub fn new(
        read: &BamRecord,
        flow_order: &str,
        max_hmer: i32,
        args: &FlowArguments,
    ) -> Result<FlowRead, FlowReadError> {
        if read.tags.get(Tag::new(FLOW_MATRIX_TAG)).is_none() {
            return Err(FlowReadError::limitation(
                "the vestigial kr/ti flow matrix layouts are not ported".to_string(),
            ));
        }
        let mut flow_read = read_flow_matrix(read, flow_order, max_hmer, args)?;
        flow_read.implement_matrix_mods(args.flow_matrix_mods.as_deref())?;
        if !args.keep_boundary_flows {
            let unmapped = read.flags & 0x4 != 0;
            let clipped = |element: Option<&htsjdk_bam::cigar::CigarElement>| {
                element.is_some_and(|element| element.op == Op::H)
            };
            if unmapped || !clipped(read.cigar.elements.first()) {
                flow_read.spread(find_first_non_zero(&flow_read.key));
            }
            if unmapped || !clipped(read.cigar.elements.last()) {
                flow_read.spread(find_last_non_zero(&flow_read.key));
            }
        }
        Ok(flow_read)
    }

    /// `spreadFlowLengthProbsAcrossCountsAtFlow`.
    fn spread(&mut self, flow: Option<usize>) {
        let Some(flow) = flow else {
            return;
        };
        let call = self.key[flow];
        let number_to_fill = self.max_hmer - call + 1;
        let mut total = 0.0;
        for i in call..=self.max_hmer {
            total += self.matrix[i as usize][flow];
        }
        let fill = (total / number_to_fill as f64).max(self.per_hmer_min_error_prob);
        for i in call..=self.max_hmer {
            self.matrix[i as usize][flow] = fill;
        }
    }

    /// `implementMatrixMods` over `getFlowMatrixModsInstructions`.
    fn implement_matrix_mods(&mut self, mods: Option<&str>) -> Result<(), FlowReadError> {
        let Some(mods) = mods else {
            return Ok(());
        };
        let rows = self.max_hmer as usize + 1;
        let mut instructions = vec![0i32; rows];
        let tokens: Vec<&str> = mods.split(',').collect();
        let parse = |text: &str| {
            text.parse::<i32>().map_err(|_| {
                FlowReadError::thrown(
                    "java.lang.NumberFormatException",
                    format!("For input string: \"{text}\""),
                )
            })
        };
        let mut i = 0;
        while i + 1 < tokens.len() {
            let hmer = parse(tokens[i])?;
            if hmer < 0 {
                return Err(FlowReadError::thrown(
                    "java.lang.IllegalArgumentException",
                    format!("the index cannot be negative: {hmer}"),
                ));
            }
            if hmer as usize >= rows {
                return Err(FlowReadError::thrown(
                    "java.lang.IllegalArgumentException",
                    format!(
                        "the index points past the last element of the collection or array: {hmer} > {}",
                        rows - 1
                    ),
                ));
            }
            instructions[hmer as usize] = parse(tokens[i + 1])?;
            i += 2;
        }
        for (hmer, hmer2) in instructions.iter().copied().enumerate() {
            if hmer2 == 0 {
                continue;
            }
            if hmer2 < 0 || hmer2 as usize >= rows {
                return Err(FlowReadError::out_of_bounds(hmer2 as i64, rows));
            }
            let hmer2 = hmer2 as usize;
            for pos in 0..self.key.len() {
                if self.matrix[hmer][pos] > self.matrix[hmer2][pos] {
                    self.matrix[hmer2][pos] = self.matrix[hmer][pos];
                }
                if hmer > hmer2 {
                    self.matrix[hmer][pos] = 0.0;
                }
            }
        }
        Ok(())
    }
}

fn find_first_non_zero(key: &[i32]) -> Option<usize> {
    key.iter().position(|value| *value != 0)
}

fn find_last_non_zero(key: &[i32]) -> Option<usize> {
    key.iter().rposition(|value| *value != 0)
}

/// `Math.pow(10, -q / 10)`, at run time.
fn phred_to_prob(q: f64) -> f64 {
    std::hint::black_box(10.0f64).powf(-q / 10.0)
}

/// `readFlowMatrix`, the production path, then `applyFilteringFlowMatrix`.
fn read_flow_matrix(
    read: &BamRecord,
    flow_order: &str,
    max_hmer: i32,
    args: &FlowArguments,
) -> Result<FlowRead, FlowReadError> {
    let key = base_array_to_key(&read.read_bases, flow_order).ok_or_else(|| {
        FlowReadError::thrown(
            "org.broadinstitute.hellbender.exceptions.GATKException",
            format!(
                "baseArrayToKey periodGuard tripped, on {}, flowOrder: {flow_order} This probably indicates the presence of a base (value) in the sequence that is not included in the provided flow order",
                String::from_utf8_lossy(&read.read_bases)
            ),
        )
    })?;

    let quals: &[u8] = &read.base_qualities;
    let mut per_hmer_min_error_prob = args.filling_value;
    let mut total_min_error_prob = per_hmer_min_error_prob;
    if per_hmer_min_error_prob == 0.0 {
        total_min_error_prob = estimate_min_error_prob(quals);
        per_hmer_min_error_prob = total_min_error_prob / max_hmer as f64;
    }

    let rows = max_hmer as usize + 1;
    let mut matrix = vec![vec![per_hmer_min_error_prob; key.len()]; rows];

    let tp: Vec<i8> = match read.tags.get(Tag::new(FLOW_MATRIX_TAG)) {
        Some(TagValue::ByteArray {
            values,
            unsigned: false,
        }) => values.clone(),
        _ => {
            return Err(FlowReadError::limitation(
                "a tp tag that is not a signed byte array".to_string(),
            ))
        }
    };
    let t0: Option<Vec<i32>> = match read.tags.get(Tag::new(FLOW_MATRIX_T0_TAG)) {
        None => None,
        Some(TagValue::Str(text)) => {
            let mut scores = Vec::with_capacity(text.len());
            for ch in text.chars() {
                if !(33..=126).contains(&(ch as u32)) {
                    return Err(FlowReadError::thrown(
                        "java.lang.IllegalArgumentException",
                        format!("Invalid fastq character: {ch}"),
                    ));
                }
                scores.push(ch as i32 - 33);
            }
            Some(scores)
        }
        Some(_) => {
            return Err(FlowReadError::limitation(
                "a t0 tag that is not a string".to_string(),
            ))
        }
    };

    let mut special_zero = false;
    if let Some(t0) = &t0 {
        if args.use_t0_tag {
            special_zero = true;
            if t0.len() != tp.len() {
                return Err(FlowReadError::thrown(
                    "org.broadinstitute.hellbender.exceptions.GATKException",
                    format!("Illegal read len(t0)!=len(qual): {}", read.read_name),
                ));
            }
        }
    }

    // SAMRecord.getBaseQualities() is signed bytes.
    let mut probs = vec![0.0; quals.len()];
    let mut t0_probs = vec![0.0; quals.len()];
    for i in 0..quals.len() {
        probs[i] = phred_to_prob(f64::from(quals[i] as i8));
        if special_zero {
            let t0 = t0.as_ref().expect("checked above");
            let value = *t0
                .get(i)
                .ok_or_else(|| FlowReadError::out_of_bounds(i as i64, t0.len()))?;
            // SAMUtils.fastqToPhred narrows to a byte.
            t0_probs[i] = phred_to_prob(f64::from(value as i8));
        }
    }

    let mut qual_ofs = 0usize;
    for i in 0..key.len() {
        let run = key[i];
        if run > 0 {
            // parseSingleHmer
            for at in qual_ofs..qual_ofs + run as usize {
                let offset = *tp
                    .get(at)
                    .ok_or_else(|| FlowReadError::out_of_bounds(at as i64, tp.len()))?;
                if offset != 0 {
                    let loc = (run + i32::from(offset)).min(max_hmer).max(0) as usize;
                    let prob = *probs
                        .get(at)
                        .ok_or_else(|| FlowReadError::out_of_bounds(at as i64, probs.len()))?;
                    if matrix[loc][i] == per_hmer_min_error_prob {
                        matrix[loc][i] = prob;
                    } else {
                        matrix[loc][i] += prob;
                    }
                }
            }
        }
        if run == 0 && special_zero {
            // parseZeroQuals
            if qual_ofs != 0 && qual_ofs != t0_probs.len() {
                let mut prob0 = t0_probs[qual_ofs - 1].min(t0_probs[qual_ofs]);
                if prob0 <= total_min_error_prob * 3.0 {
                    prob0 = 0.0;
                }
                matrix[1][i] = matrix[1][i].max(prob0);
            }
        }
        let mut total_error_prob = 0.0;
        for row in matrix.iter().take(max_hmer as usize) {
            total_error_prob += row[i];
        }
        let call_prob = MINIMAL_CALL_PROB.max(1.0 - total_error_prob);
        matrix[run.min(max_hmer) as usize][i] = call_prob;
        qual_ofs += run as usize;
    }

    let mut flow_read = FlowRead {
        key,
        max_hmer,
        matrix,
        per_hmer_min_error_prob,
    };
    flow_read.apply_filtering(args)?;
    Ok(flow_read)
}

/// `estimateMinErrorProb`: the read's best quality as a probability, 40 when it has none.
fn estimate_min_error_prob(quals: &[u8]) -> f64 {
    let mut max_qual = 0.0f64;
    for qual in quals {
        let qual = f64::from(*qual as i8);
        if qual > max_qual {
            max_qual = qual;
        }
    }
    if max_qual == 0.0 {
        max_qual = 40.0;
    }
    phred_to_prob(max_qual)
}

impl FlowRead {
    /// `applyFilteringFlowMatrix`, in the reference's order.
    fn apply_filtering(&mut self, args: &FlowArguments) -> Result<(), FlowReadError> {
        let min = self.per_hmer_min_error_prob;
        let rows = self.max_hmer as usize + 1;
        let flows = self.key.len();
        if args.disallow_larger_probs {
            for row in self.matrix.iter_mut() {
                for cell in row.iter_mut() {
                    if *cell > 1.0 {
                        *cell = 1.0;
                    }
                }
            }
        }
        if args.remove_longer_than_one_indels {
            for i in 0..flows {
                for j in 0..rows {
                    if (j as i32 - self.key[i]).abs() > 1 {
                        self.matrix[j][i] = min;
                    }
                }
            }
        }
        if args.remove_one_to_zero_probs {
            for i in 0..flows {
                if self.key[i] == 0 {
                    for j in 1..rows {
                        self.matrix[j][i] = min;
                    }
                }
            }
        }
        if args.lump_probs {
            for i in 0..self.max_hmer as usize {
                for j in 0..flows {
                    let fkey = self.key[j];
                    if self.matrix[i][j] <= min {
                        continue;
                    }
                    let target = if (i as i32 - fkey) < -1 {
                        fkey - 1
                    } else if (i as i32 - fkey) > 1 {
                        fkey + 1
                    } else {
                        continue;
                    };
                    if target as usize >= rows {
                        return Err(FlowReadError::out_of_bounds(target as i64, rows));
                    }
                    let value = self.matrix[i][j];
                    self.matrix[target as usize][j] += value;
                    self.matrix[i][j] = min;
                }
            }
        }
        // clipProbs
        let threshold = min * 3.0;
        for i in 0..self.max_hmer as usize {
            for j in 0..flows {
                if self.matrix[i][j] <= threshold && self.key[j] != i as i32 {
                    self.matrix[i][j] = min;
                }
            }
        }
        if args.symmetric_indels {
            for i in 0..flows {
                let idx = self.key[i];
                if idx > 1 && idx < self.max_hmer {
                    let idx = idx as usize;
                    let prob = (self.matrix[idx - 1][i] + self.matrix[idx + 1][i]) / 2.0;
                    self.matrix[idx - 1][i] = prob;
                    self.matrix[idx + 1][i] = prob;
                }
            }
        }
        if args.only_ins_or_del {
            for i in 0..flows {
                let idx = self.key[i];
                if idx > 1 && idx < self.max_hmer {
                    let idx = idx as usize;
                    if self.matrix[idx - 1][i] > min && self.matrix[idx + 1][i] > min {
                        let fix = if self.matrix[idx - 1][i] > self.matrix[idx + 1][i] {
                            idx + 1
                        } else {
                            idx - 1
                        };
                        self.matrix[fix][i] = min;
                    }
                }
            }
        }
        if args.retain_max_n_probs {
            let hmers = self.max_hmer as usize;
            for i in 0..flows {
                let column: Vec<f64> = (0..hmers).map(|j| self.matrix[j][i]).collect();
                let k = ((self.key[i] + 1) / 2) as usize;
                let kth = find_kth_largest(&column, k + 1);
                for j in 0..hmers {
                    if self.matrix[j][i] < kth {
                        self.matrix[j][i] = min;
                    }
                }
            }
        }
        Ok(())
    }
}

/// `findKthLargest`: a min-heap of `k` that keeps the largest seen, whose head is then the k-th
/// largest, or the smallest of all when there are fewer than `k`.
fn find_kth_largest(nums: &[f64], k: usize) -> f64 {
    let mut sorted: Vec<f64> = nums.to_vec();
    sorted.sort_by(|a, b| b.total_cmp(a));
    let at = k.min(sorted.len()).saturating_sub(1);
    sorted[at]
}

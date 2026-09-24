//! `FlowPairHMMAlignReadsToHaplotypes`: every read scored against every haplotype, and what the
//! concise format makes of it.
//!
//! The `FlowBased` engine is ported as the tool runs it, with exact matching: the read's flow
//! matrix read off along the haplotype's key, the start stepped by four flows. The `FlowBasedHMM`
//! engine is not. Above all the concise format's reference score is recorded from inside the branch
//! that raises the best score, and so depends on the order of the haplotype FASTA.
//!
//! Ported from
//! `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AlleleLikelihoodWriter`,
//! `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ConciseAlleleLikelihoodWriter`
//! and `org.broadinstitute.hellbender.tools.walkers.featuremapping.FlowPairHMMAlignReadsToHaplotypes`
//! in GATK 4.6.2.0.

/// One haplotype of the FASTA: its name, and whether `--ref-haplotype` named it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Haplotype {
    pub name: String,
    pub is_reference: bool,
}

impl Haplotype {
    /// `new Haplotype(bases, name.equals(refHaplotypeName))`, which is an equality on the NAME.
    ///
    /// A `--ref-haplotype` the FASTA does not name leaves every haplotype non-reference, which is
    /// the same state as not passing the argument at all.
    pub fn new(name: &str, reference_name: Option<&str>) -> Haplotype {
        Haplotype {
            name: name.to_string(),
            is_reference: reference_name == Some(name),
        }
    }
}

/// `%.3f` and `%.03f`, which are the same format written two ways.
///
/// An infinity is printed by name rather than as a number, so the concise format's reference
/// column can read `Infinity` and its score column `-Infinity`.
pub fn three_decimals(value: f64) -> String {
    if value.is_nan() {
        return "NaN".to_string();
    }
    if value.is_infinite() {
        return if value > 0.0 { "Infinity" } else { "-Infinity" }.to_string();
    }
    format!("{value:.3}")
}

// ================================================================================================
// The expanded format.
// ================================================================================================

/// The header of the expanded format: `Read` and then one column per haplotype.
pub fn expanded_header(haplotypes: &[Haplotype]) -> String {
    let mut line = String::from("Read");
    for haplotype in haplotypes {
        line.push('\t');
        line.push_str(&haplotype.name);
    }
    line
}

/// One row of the expanded format: the read's name and its score against each haplotype.
pub fn expanded_row(read: &str, scores: &[f64]) -> String {
    let mut line = String::from(read);
    for score in scores {
        line.push('\t');
        line.push_str(&three_decimals(*score));
    }
    line
}

/// The whole expanded file.
///
/// The header is written once, before the first buffer of reads, rather than once per buffer.
pub fn expanded_file(haplotypes: &[Haplotype], rows: &[(String, Vec<f64>)]) -> String {
    let mut text = expanded_header(haplotypes);
    text.push('\n');
    for (read, scores) in rows {
        text.push_str(&expanded_row(read, scores));
        text.push('\n');
    }
    text
}

// ================================================================================================
// The concise format.
// ================================================================================================

/// The header of the concise format, whose five columns are fixed.
pub const CONCISE_HEADER: &str = "Read\tBest_hap\tBest_score\tDiff_from_second\tDiff_from_ref";

/// What one read's scores reduce to in the concise format.
#[derive(Debug, Clone, PartialEq)]
pub struct Concise {
    /// The best haplotype's name, EMPTY when no score beat negative infinity.
    pub best_haplotype: String,
    pub best_score: f64,
    pub second_best_score: f64,
    /// The reference's score, recorded only while the reference was the best so far.
    pub reference_score: f64,
}

impl Concise {
    pub fn difference_from_second(&self) -> f64 {
        self.best_score - self.second_best_score
    }

    pub fn difference_from_reference(&self) -> f64 {
        self.best_score - self.reference_score
    }
}

/// The concise reduction of one read's row.
///
/// Three things follow from where the assignments sit. The best score is raised only by a STRICTLY
/// greater score, so the first of an exact tie wins. The second-best score is raised either by the
/// score the best one displaces or by a score that beats it without beating the best. And the
/// REFERENCE SCORE IS ASSIGNED INSIDE THE FIRST BRANCH, so it is recorded only while the reference
/// is the best haplotype seen so far: a reference that comes after a better haplotype leaves it at
/// negative infinity, and the difference from it is then an infinity.
pub fn concise(haplotypes: &[Haplotype], scores: &[f64]) -> Concise {
    let mut best_haplotype = String::new();
    let mut best_score = f64::NEG_INFINITY;
    let mut second_best_score = f64::NEG_INFINITY;
    let mut reference_score = f64::NEG_INFINITY;
    for (haplotype, score) in haplotypes.iter().zip(scores.iter()) {
        if *score > best_score {
            second_best_score = best_score;
            best_score = *score;
            best_haplotype = haplotype.name.clone();
            if haplotype.is_reference {
                reference_score = *score;
            }
        } else if *score > second_best_score {
            second_best_score = *score;
        }
    }
    Concise {
        best_haplotype,
        best_score,
        second_best_score,
        reference_score,
    }
}

/// One row of the concise format.
pub fn concise_row(read: &str, concise: &Concise) -> String {
    format!(
        "{read}\t{}\t{}\t{}\t{}",
        concise.best_haplotype,
        three_decimals(concise.best_score),
        three_decimals(concise.difference_from_second()),
        three_decimals(concise.difference_from_reference())
    )
}

/// The whole concise file.
pub fn concise_file(haplotypes: &[Haplotype], rows: &[(String, Vec<f64>)]) -> String {
    let mut text = String::from(CONCISE_HEADER);
    text.push('\n');
    for (read, scores) in rows {
        text.push_str(&concise_row(read, &concise(haplotypes, scores)));
        text.push('\n');
    }
    text
}

// ================================================================================================
// The walk.
// ================================================================================================

/// `BUFFER_SIZE_LIMIT`: how many reads are scored in one go.
///
/// The buffer decides only when the matrix is computed, not what it holds, so a file of fewer
/// reads than this is one buffer and comes out the same as any other.
pub const BUFFER_SIZE_LIMIT: usize = 50;

/// The buffers a run of `count` reads is split into, as their lengths.
pub fn buffers(count: usize) -> Vec<usize> {
    let mut lengths = vec![BUFFER_SIZE_LIMIT; count / BUFFER_SIZE_LIMIT];
    // The last buffer is flushed at the end of the traversal whatever its length, INCLUDING when
    // it is empty: a count that is a multiple of the limit still gets a final empty flush.
    lengths.push(count % BUFFER_SIZE_LIMIT);
    lengths
}

/// The two engines the tool accepts, by the name `-E` takes.
pub const ENGINES: [&str; 2] = ["FlowBased", "FlowBasedHMM"];

/// The refusal any other engine produces.
///
/// It is a bare `RuntimeException` rather than a `UserException`, so it reads as a crash and not
/// as an argument being wrong.
pub const UNKNOWN_ENGINE_MESSAGE: &str = "Accepted engines are FlowBasedHMM or FlowBased";

pub fn is_known_engine(name: &str) -> bool {
    ENGINES.contains(&name)
}

// ================================================================================================
// The FlowBased engine.
// ================================================================================================

/// `FlowBasedAlignmentLikelihoodEngine.ALIGNMENT_UNCERTAINTY`, the step of the start search.
const ALIGNMENT_UNCERTAINTY: usize = 4;
/// `COMMON_PROB_VALUE_1` and `_2`, whose logarithms the optimised computation looks up.
const COMMON_PROB_VALUE_1: f64 = 0.988;
const COMMON_PROB_VALUE_2: f64 = 0.001;

/// `Precision.equals(x, y)`: equal to within one ulp.
fn within_one_ulp(x: f64, y: f64) -> bool {
    if x.is_nan() || y.is_nan() {
        return false;
    }
    let order = |value: f64| {
        let bits = value.to_bits() as i64;
        if bits < 0 {
            i64::MIN - bits
        } else {
            bits
        }
    };
    (order(x) - order(y)).abs() <= 1
}

/// `Math.log10`, at run time.
fn log10(value: f64) -> f64 {
    std::hint::black_box(value).log10()
}

/// `FlowBasedHaplotype`: the haplotype's key and the base each of its flows reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FlowHaplotype {
    pub key: Vec<i32>,
    pub flow_order: Vec<u8>,
}

impl FlowHaplotype {
    /// `None` where `baseArrayToKey`'s period guard trips.
    pub fn new(bases: &[u8], flow_order: &str) -> Option<FlowHaplotype> {
        let key = crate::flow_based_read::base_array_to_key(bases, flow_order)?;
        let order = flow_order.as_bytes();
        let flow_order = (0..key.len()).map(|i| order[i % order.len()]).collect();
        Some(FlowHaplotype { key, flow_order })
    }
}

/// `haplotypeReadMatchingExactLength`: the whole haplotype key, from the first flow reading the
/// read's first base, stepping by four flows, the best sum of the read's log probabilities.
///
/// A haplotype whose key is shorter than the read's never fits, and scores negative infinity.
pub fn exact_length_score(
    haplotype: &FlowHaplotype,
    read: &crate::flow_based_read::FlowRead,
    optimized: bool,
) -> f64 {
    let key_length = haplotype.key.len();
    let read_first = read.flow_order.first().copied().unwrap_or(0);
    let starting_point = haplotype
        .flow_order
        .iter()
        .position(|base| *base == read_first)
        .unwrap_or(0);
    let read_key_length = read.key.len();
    let cap = read.max_hmer + 1;
    let mut best = f64::NEG_INFINITY;
    let mut s = starting_point;
    let mut locations = vec![0i32; read_key_length];
    while s + read_key_length <= key_length {
        for i in s..s + read_key_length {
            locations[i - s] = (haplotype.key[i] & 0xff).min(cap);
        }
        let mut result = 0.0;
        for (i, location) in locations.iter().enumerate() {
            let prob = read.prob(i, *location);
            if optimized {
                result += if within_one_ulp(prob, COMMON_PROB_VALUE_1) {
                    log10(COMMON_PROB_VALUE_1)
                } else if within_one_ulp(prob, COMMON_PROB_VALUE_2) {
                    log10(COMMON_PROB_VALUE_2)
                } else {
                    log10(prob)
                };
                if result < best {
                    break;
                }
            } else {
                result += log10(prob);
            }
        }
        if result > best {
            best = result;
        }
        s += ALIGNMENT_UNCERTAINTY;
    }
    best
}

/// `FastaSequenceFile(path, false)`: every record, its name the whole header line.
pub fn read_haplotype_fasta(text: &str) -> Vec<(String, Vec<u8>)> {
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    for line in text.lines() {
        if let Some(name) = line.strip_prefix('>') {
            records.push((name.trim_end().to_string(), Vec::new()));
        } else if let Some((_, bases)) = records.last_mut() {
            bases.extend(line.trim_end().bytes());
        }
    }
    records
}

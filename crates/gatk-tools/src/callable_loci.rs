//! `CallableLoci`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.coverage.CallableLoci` (GATK 4.6.2.0).
//!
//! Every locus of the requested intervals is given one of six states, the runs of equal state are
//! written as BED, and the six counts are written as a summary.
//!
//! # The run test never compares contigs
//!
//! ```java
//! } else if (callableState.getStart() != currentState.getEnd() + 1 || currentState.getState() != callableState.getState()) {
//! ```
//!
//! A run continues when the next locus starts one past the previous end and carries the same
//! state. Nothing there asks whether the contig changed, and `updateInterval` keeps the FIRST
//! interval's contig, so two stretches on different contigs whose coordinates happen to run on come
//! out as a single BED line under the first contig's name. The `contig-run-on` row of the golden is
//! exactly that: twenty bases of `REF_N` spanning `chr1` and `chr2`, written as `chr1`.
//!
//! This port reproduces it. A tool that repaired the line here would disagree with the reference on
//! any interval list whose contigs run on, which is what `SplitIntervals` hands out.
//!
//! # A deletion is callable whatever its base quality
//!
//! ```java
//! if (e.getMappingQual() >= minMappingQuality && (e.getQual() >= minBaseQuality || e.isDeletion()))
//! ```
//!
//! The deletion is an alternative to the base-quality test, not to the mapping-quality one, so a
//! deleted position can be CALLABLE on reads whose bases would not have counted at all.
//!
//! # The order of the tests is the answer
//!
//! No coverage, then poor mapping quality, then low coverage, then excessive coverage, then
//! callable. A locus whose reads are mostly low-quality mappings is POOR_MAPPING_QUALITY even when
//! its QC depth would have called it LOW_COVERAGE, and EXCESSIVE_COVERAGE is tested on the RAW
//! depth rather than on the QC depth.

/// `CallableLoci.State`, in the order the enum declares them, which is the summary's order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum State {
    RefN,
    Callable,
    NoCoverage,
    LowCoverage,
    ExcessiveCoverage,
    PoorMappingQuality,
}

impl State {
    /// The enum constant's own name, which is what both files carry.
    pub fn name(&self) -> &'static str {
        match self {
            State::RefN => "REF_N",
            State::Callable => "CALLABLE",
            State::NoCoverage => "NO_COVERAGE",
            State::LowCoverage => "LOW_COVERAGE",
            State::ExcessiveCoverage => "EXCESSIVE_COVERAGE",
            State::PoorMappingQuality => "POOR_MAPPING_QUALITY",
        }
    }

    /// `State.values()`.
    pub const ALL: [State; 6] = [
        State::RefN,
        State::Callable,
        State::NoCoverage,
        State::LowCoverage,
        State::ExcessiveCoverage,
        State::PoorMappingQuality,
    ];
}

/// The thresholds, with the reference's defaults.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Arguments {
    pub max_low_mapq: i32,
    pub min_mapping_quality: i32,
    pub min_base_quality: i32,
    pub min_depth: i32,
    pub max_depth: Option<i32>,
    pub min_depth_low_mapq: i32,
    pub max_low_mapq_fraction: f64,
}

impl Default for Arguments {
    fn default() -> Self {
        Arguments {
            max_low_mapq: 1,
            min_mapping_quality: 10,
            min_base_quality: 20,
            min_depth: 4,
            max_depth: None,
            min_depth_low_mapq: 10,
            max_low_mapq_fraction: 0.1,
        }
    }
}

/// As much of a pileup element as the state depends on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Element {
    pub mapping_quality: i32,
    /// `getQual()`, which for a deletion is the constant the pileup carries.
    pub base_quality: i32,
    pub is_deletion: bool,
}

/// `getCurrentState`.
///
/// `reference_base` is the base at the locus; an N there is answered before any depth is counted.
pub fn state_at(reference_base: u8, pileup: &[Element], arguments: &Arguments) -> State {
    if reference_base.eq_ignore_ascii_case(&b'N') {
        return State::RefN;
    }
    let mut raw_depth = 0;
    let mut qc_depth = 0;
    let mut low_mapq_depth = 0;
    for element in pileup {
        raw_depth += 1;
        // `<=` here and `>=` below, so thresholds set equal count a read in both.
        if element.mapping_quality <= arguments.max_low_mapq {
            low_mapq_depth += 1;
        }
        if element.mapping_quality >= arguments.min_mapping_quality
            && (element.base_quality >= arguments.min_base_quality || element.is_deletion)
        {
            qc_depth += 1;
        }
    }
    if raw_depth == 0 {
        return State::NoCoverage;
    }
    if raw_depth >= arguments.min_depth_low_mapq
        && f64::from(low_mapq_depth) / f64::from(raw_depth) >= arguments.max_low_mapq_fraction
    {
        return State::PoorMappingQuality;
    }
    if qc_depth < arguments.min_depth {
        return State::LowCoverage;
    }
    if arguments
        .max_depth
        .is_some_and(|maximum| raw_depth >= maximum)
    {
        return State::ExcessiveCoverage;
    }
    State::Callable
}

/// One BED line's worth of run.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Run {
    pub contig: String,
    /// One-based inclusive, as the tool holds it; the BED line subtracts one from the start.
    pub start: i32,
    pub end: i32,
    pub state: State,
}

impl Run {
    /// `toBedString`: zero-based start, one-based end.
    pub fn bed_line(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}",
            self.contig,
            self.start - 1,
            self.end,
            self.state.name()
        )
    }
}

/// `OutputFormat`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    Bed,
    StatePerBase,
}

/// The whole traversal: the loci in the order they are visited, each with its state.
///
/// Returns the BED file's text and the summary's, which is what the two output streams get.
pub fn write(loci: &[(String, i32, State)], format: OutputFormat) -> (String, String) {
    let mut bed = String::new();
    let mut counts = [0_i64; 6];
    let mut current: Option<Run> = None;
    for (contig, position, state) in loci {
        counts[State::ALL.iter().position(|s| s == state).expect("a state")] += 1;
        let this = Run {
            contig: contig.clone(),
            start: *position,
            end: *position,
            state: *state,
        };
        match format {
            OutputFormat::StatePerBase => {
                bed.push_str(&this.bed_line());
                bed.push('\n');
            }
            OutputFormat::Bed => match &mut current {
                None => current = Some(this),
                // The contig is not part of this test, which is what merges runs across contigs.
                Some(run) if this.start == run.end + 1 && run.state == this.state => {
                    run.end = this.end;
                }
                Some(run) => {
                    bed.push_str(&run.bed_line());
                    bed.push('\n');
                    current = Some(this);
                }
            },
        }
    }
    if format == OutputFormat::Bed {
        if let Some(run) = &current {
            bed.push_str(&run.bed_line());
            bed.push('\n');
        }
    }

    let mut summary = format!("{:>30} {}\n", "state", "nBases");
    for (index, state) in State::ALL.iter().enumerate() {
        summary.push_str(&format!("{:>30} {}\n", state.name(), counts[index]));
    }
    (bed, summary)
}

/// `UserException.BadInput` when the reads carry more than one sample.
pub fn sample_refusal(samples: &[String]) -> String {
    format!(
        "org.broadinstitute.hellbender.exceptions.UserException$BadInput:Bad input: CallableLoci \
         only works for a single sample.  Found {} samples ({}).",
        samples.len(),
        samples.join(", ")
    )
}

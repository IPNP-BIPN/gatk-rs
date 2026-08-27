//! `CombineGVCFs`: how several single-sample GVCFs become one multi-sample one.
//!
//! Every sample's boundaries cut every other sample's blocks, so the output's records are the union
//! of every input's edges rather than any one input's.
//!
//! Reading and writing the GVCFs are not ported, nor are the likelihoods the merge expands. Which
//! records the output has, which samples carry data on each, and what the two band arguments do
//! are.

/// One record of one input GVCF.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub sample: String,
    pub start: i32,
    pub end: i32,
    /// Empty for a reference block, whose only alternate is `<NON_REF>`.
    pub alternates: Vec<String>,
    pub genotype_quality: i32,
}

impl Record {
    pub fn is_reference_block(&self) -> bool {
        self.alternates.is_empty()
    }
}

/// One record of the merged output: a span, and what each sample contributes to it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Merged {
    pub start: i32,
    pub end: i32,
    /// The union of every carrying sample's alternates, in the order the samples are listed.
    pub alternates: Vec<String>,
    /// One entry per sample. `None` where that sample's inputs had run out.
    pub qualities: Vec<Option<i32>>,
}

/// The two band arguments.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct BandArguments {
    /// Every block becomes single bases.
    pub base_pair_resolution: bool,
    /// Blocks are also cut at every multiple of this, when it is above zero.
    pub break_bands_at_multiples_of: i32,
}

impl BandArguments {
    /// The two are NOT mutually exclusive, and base-pair resolution WINS: given together, the grid
    /// is ignored entirely.
    pub fn effective_grid(&self) -> Option<i32> {
        if self.base_pair_resolution {
            Some(1)
        } else if self.break_bands_at_multiples_of > 0 {
            Some(self.break_bands_at_multiples_of)
        } else {
            None
        }
    }
}

/// Every position at which the output must start a new record.
///
/// Each input contributes the start of every record it has, and the base after the end of every
/// record it has. A grid contributes its own multiples. The union of those, clipped to the span
/// the inputs cover, is the output's set of boundaries.
pub fn boundaries(records: &[Record], arguments: &BandArguments) -> Vec<i32> {
    let Some(first) = records.iter().map(|record| record.start).min() else {
        return Vec::new();
    };
    let last = records
        .iter()
        .map(|record| record.end)
        .max()
        .expect("a record");

    let mut out: Vec<i32> = Vec::new();
    for record in records {
        out.push(record.start);
        if record.end < last {
            out.push(record.end + 1);
        }
    }
    if let Some(grid) = arguments.effective_grid() {
        let mut at = (first / grid) * grid;
        if at < first {
            at += grid;
        }
        while at <= last {
            out.push(at);
            at += grid;
        }
    }
    out.retain(|position| *position >= first && *position <= last);
    out.sort();
    out.dedup();
    out
}

/// The merged records.
///
/// A sample whose inputs ran out before a span keeps its column and carries nothing, which the
/// writer renders as `./.` with no fields after it. Nothing is padded with reference.
pub fn combine(records: &[Record], samples: &[String], arguments: &BandArguments) -> Vec<Merged> {
    let edges = boundaries(records, arguments);
    let last = records
        .iter()
        .map(|record| record.end)
        .max()
        .unwrap_or_default();
    let mut out = Vec::new();
    for (index, start) in edges.iter().enumerate() {
        let end = match edges.get(index + 1) {
            Some(next) => next - 1,
            None => last,
        };
        // A variant's record covers its own base alone, whatever the grid says.
        let variant_here = records
            .iter()
            .find(|record| !record.is_reference_block() && record.start == *start);
        let end = match variant_here {
            Some(record) => record.end.min(end),
            None => end,
        };
        let mut alternates: Vec<String> = Vec::new();
        let mut qualities = Vec::new();
        for sample in samples {
            let covering = records.iter().find(|record| {
                record.sample == *sample && record.start <= *start && *start <= record.end
            });
            match covering {
                Some(record) => {
                    for allele in &record.alternates {
                        if !alternates.contains(allele) {
                            alternates.push(allele.clone());
                        }
                    }
                    qualities.push(Some(record.genotype_quality));
                }
                None => qualities.push(None),
            }
        }
        out.push(Merged {
            start: *start,
            end,
            alternates,
            qualities,
        });
    }
    out
}

/// What the tool refuses about its inputs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CombineError {
    /// The same file twice, caught by the feature-input check rather than by anything about
    /// samples.
    DuplicateInput { path: String },
}

impl CombineError {
    pub fn message(&self) -> String {
        match self {
            CombineError::DuplicateInput { path } => format!(
                "Bad input: Feature inputs must be unique, but {path} was specified more than once"
            ),
        }
    }
}

/// The check the walker makes on its driving inputs.
pub fn check_inputs(paths: &[String]) -> Result<(), CombineError> {
    for (index, path) in paths.iter().enumerate() {
        if paths[..index].contains(path) {
            return Err(CombineError::DuplicateInput { path: path.clone() });
        }
    }
    Ok(())
}

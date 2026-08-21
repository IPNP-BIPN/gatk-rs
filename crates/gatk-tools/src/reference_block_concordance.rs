//! `ReferenceBlockConcordance`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.validation.ReferenceBlockConcordance`
//! (GATK 4.6.2.0).
//!
//! Two GVCFs walked side by side by [`gatk_engine::concordance_walker`], and the reference blocks
//! they hold turned into three metrics files.
//!
//! # The histograms are keyed by a Java `Pair`'s `toString`, and sorted as strings
//!
//! The tool's own comment says the key should have been a pair of integers and is a `String`
//! because `MetricsFile` cannot read arbitrary types. `Pair.toString` is `length,GQ` with no space
//! and no brackets, and the metrics file sorts its bins as strings, so the golden's truth
//! histogram runs `1,80`, `100,20`, `50,40`, `50,60`: neither the lengths' order nor the file's.
//!
//! # The concordance histogram counts BASES, not blocks
//!
//! ```java
//! confidenceConcordanceHistogram.increment(pair.toString(), truthInterval.intersect(evalInterval).getLengthOnReference());
//! ```
//!
//! It fires only while both sides currently hold a block, and each side's held block is cleared as
//! soon as a step no longer overlaps it. So the value is the number of bases the two blocks share,
//! and a pair of blocks that merely appear in the same file contributes nothing.
//!
//! # Only the genotype decides what is walked
//!
//! Both filters are `isHomRef` on genotype 0, so a filtered block is walked like any other and a
//! variant site is dropped from both sides. A multi-sample file passes that filter on its first
//! sample and is refused afterwards, by the length extraction.

use htsjdk_metrics::file::{Histogram, MetricsFile};

/// One reference block, as far as this tool reads one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Block {
    pub contig: String,
    pub start: i32,
    /// `getEnd()`, which the `END` attribute sets.
    pub end: i32,
    /// `getGQ()` of genotype 0.
    pub gq: i32,
    /// Whether the record's genotype 0 is hom-ref, which is the only filter the tool has.
    pub is_hom_ref: bool,
    /// How many genotypes the record carries, which the length extraction refuses unless it is 1.
    pub genotypes: usize,
    /// The record as `toStringDecodeGenotypes` renders it, which the refusal quotes.
    pub rendered: String,
}

impl Block {
    /// `getLengthOnReference()`.
    pub fn length(&self) -> i32 {
        self.end - self.start + 1
    }
}

/// The refusal the length extraction raises.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MultiSample {
    pub rendered: String,
}

impl MultiSample {
    pub fn java_class(&self) -> &str {
        "java.lang.IllegalStateException"
    }

    pub fn message(&self) -> String {
        format!(
            "A multisample GVCF file was provided, however, only single sample GVCFs are currently \
             supported. This occurred when reading \"{}\".",
            self.rendered
        )
    }
}

/// A histogram keyed by a string, which is what the tool builds: counts in insertion order, sorted
/// as strings when it is written.
#[derive(Debug, Clone, Default)]
pub struct StringHistogram {
    bins: Vec<(String, f64)>,
}

impl StringHistogram {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn increment(&mut self, key: &str) {
        self.increment_by(key, 1.0);
    }

    pub fn increment_by(&mut self, key: &str, amount: f64) {
        match self.bins.iter_mut().find(|(bin, _)| bin == key) {
            Some((_, value)) => *value += amount,
            None => self.bins.push((key.to_string(), amount)),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.bins.is_empty()
    }

    /// The bins as the metrics file writes them: sorted by the key AS A STRING.
    pub fn sorted(&self) -> Vec<(String, f64)> {
        let mut sorted = self.bins.clone();
        sorted.sort_by(|left, right| left.0.cmp(&right.0));
        sorted
    }
}

/// `new Pair<>(a, b).toString()`, which is what both keys are built from.
pub fn pair_key(first: i32, second: i32) -> String {
    format!("{first},{second}")
}

/// The three histograms one run produces.
#[derive(Debug, Clone, Default)]
pub struct Histograms {
    pub truth_blocks: StringHistogram,
    pub eval_blocks: StringHistogram,
    pub confidence_concordance: StringHistogram,
}

/// `apply` over every step of the walk.
///
/// `steps` are the walker's own, as `(truth index, eval index)`; the two slices hold the blocks
/// that survived the hom-ref filter, in file order.
pub fn accumulate(
    truth: &[Block],
    eval: &[Block],
    steps: &[(Option<usize>, Option<usize>)],
) -> Result<Histograms, MultiSample> {
    let mut histograms = Histograms::default();
    let mut current_truth: Option<&Block> = None;
    let mut current_eval: Option<&Block> = None;

    for (truth_index, eval_index) in steps {
        // The step's own span, which is whichever side it carries: both sides of a step are at the
        // same locus, so either one answers.
        let step = match (truth_index, eval_index) {
            (Some(index), _) => &truth[*index],
            (None, Some(index)) => &eval[*index],
            (None, None) => continue,
        };

        if let Some(index) = truth_index {
            let block = &truth[*index];
            histograms.truth_blocks.increment(&key_of(block)?);
            current_truth = Some(block);
        }
        if let Some(index) = eval_index {
            let block = &eval[*index];
            histograms.eval_blocks.increment(&key_of(block)?);
            current_eval = Some(block);
        }

        // `!currentTruthVariantContext.overlaps(truthVersusEval)`, which drops a held block as soon
        // as the walk has passed it.
        if let Some(block) = current_truth {
            if !overlaps(block, step) {
                current_truth = None;
            }
        }
        if let Some(block) = current_eval {
            if !overlaps(block, step) {
                current_eval = None;
            }
        }

        if let (Some(held_truth), Some(held_eval)) = (current_truth, current_eval) {
            if overlaps(held_truth, held_eval) {
                let bases = intersection_length(held_truth, held_eval);
                histograms
                    .confidence_concordance
                    .increment_by(&pair_key(held_truth.gq, held_eval.gq), f64::from(bases));
            }
        }
    }
    Ok(histograms)
}

/// `extractLengthAndGQ`, which refuses anything but a single sample.
fn key_of(block: &Block) -> Result<String, MultiSample> {
    if block.genotypes != 1 {
        return Err(MultiSample {
            rendered: block.rendered.clone(),
        });
    }
    Ok(pair_key(block.length(), block.gq))
}

fn overlaps(left: &Block, right: &Block) -> bool {
    left.contig == right.contig && left.start <= right.end && right.start <= left.end
}

/// `SimpleInterval.intersect(...).getLengthOnReference()`.
fn intersection_length(left: &Block, right: &Block) -> i32 {
    left.end.min(right.end) - left.start.max(right.start) + 1
}

/// One of the three metrics files, which carry a histogram and no metric rows.
///
/// `headers` are the two the engine adds, which the golden masks.
pub fn write_histogram(histogram: &StringHistogram, headers: &[String]) -> String {
    let mut file = MetricsFile::new();
    for header in headers {
        file.add_header(header);
    }
    file.histograms.push(Histogram {
        bin_label: "BIN".to_string(),
        value_label: "VALUE".to_string(),
        key_class: "java.lang.String".to_string(),
        bins: histogram.sorted(),
    });
    file.write()
}

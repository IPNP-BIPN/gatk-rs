//! `CollectF1R2Counts`, ported from `F1R2CountsCollector`, `DepthOneHistograms`,
//! `F1R2FilterConstants` and `F1R2FilterUtils` (GATK 4.6.2.0).
//!
//! The counts Mutect2's orientation-bias filter is trained on: for every three-base reference
//! context, how deep the reference sites were, and what the alt sites looked like. The traversal is
//! a `LocusWalker`, already ported; what is here is the collector, which decides which pileups
//! count at all and which of three places a site is written to.
//!
//! # A site goes to exactly one of three places
//!
//! A reference site increments a depth histogram for its context. An alt site with exactly one alt
//! read increments a depth-one histogram keyed by context, alt base and orientation. Any other alt
//! site becomes a row of the alt table. The alt base is chosen by
//!
//! ```java
//! final int[] baseCountsCopy = Arrays.copyOf(baseCounts, baseCounts.length);
//! baseCountsCopy[refBase.ordinal()] = -1;
//! final int altBaseIndex = MathUtils.maxElementIndex(baseCountsCopy);
//! final boolean referenceSite = baseCounts[altBaseIndex] == 0;
//! ```
//!
//! so a tie between two alt bases goes to the lower base index, and a pileup with no alt read at
//! all reads as a reference site. The alt table's `depth` column is `refCount + altCount` rather
//! than the pileup's depth, so a site with two different alt bases reports less than it saw.
//!
//! # Every skip is a skip of the whole site
//!
//! ```java
//! for (final Map.Entry<String, ReadPileup> entry : splitPileup.entrySet()) {
//!     ...
//!     if (!isPileupGood(samplePileup)) { return; }
//!     ...
//!     if (referenceSite) { ...; return; }
//!     if (altCount == 1) { ...; return; }
//! ```
//!
//! Every branch but the alt-table one leaves `process` rather than the iteration, so at any site
//! only the samples up to the first one that returns are counted at all. The order they arrive in
//! is a `HashMap` order over the sample names, not the header's: over `alpha` and `bravo` it is
//! `bravo` first, and every one of alpha's alt sites is lost because bravo is a reference site at
//! each of them. [`sample_order`] is that order.
//!
//! # The shape of the output does not depend on the data
//!
//! All 64 contexts are present whether or not they were seen, and every histogram is prefilled with
//! the bins from one to `--f1r2-max-depth`. The reference histograms come out in the `HashMap`
//! order of their 64 context strings, which is reproducible; the alt histograms are keyed by a pair
//! holding an enum, whose `hashCode` is an identity hash, so THEIR order is not reproducible from
//! one JVM to the next and [`alt_histograms`] hands them back sorted by label instead.

use gatk_engine::java_hash::{hash_map_order, string_hash_code};
use htsjdk_metrics::file::Histogram;
use std::collections::HashMap;

/// `F1R2FilterConstants.REF_CONTEXT_PADDING`, which sets the k-mer size to three.
pub const REF_CONTEXT_PADDING: usize = 1;

/// `F1R2FilterConstants.REFERENCE_CONTEXT_SIZE`.
pub const REFERENCE_CONTEXT_SIZE: usize = 2 * REF_CONTEXT_PADDING + 1;

/// `F1R2FilterConstants.DEFAULT_MAX_DEPTH`.
pub const DEFAULT_MAX_DEPTH: i32 = 200;

/// `F1R2FilterConstants.binName`.
pub const BIN_NAME: &str = "depth";

/// The four bases, in the order `BaseUtils` indexes them, which is what settles a tie.
pub const BASES: &[u8; 4] = b"ACGT";

/// `ReadOrientation`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ReadOrientation {
    F1R2,
    F2R1,
}

impl ReadOrientation {
    pub fn label(&self) -> &'static str {
        match self {
            ReadOrientation::F1R2 => "F1R2",
            ReadOrientation::F2R1 => "F2R1",
        }
    }
}

/// `CollectF1R2CountsArgumentCollection`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Args {
    pub min_median_map_qual: i32,
    pub min_base_quality: i32,
    pub max_depth: i32,
}

impl Default for Args {
    fn default() -> Self {
        Args {
            min_median_map_qual: 50,
            min_base_quality: 20,
            max_depth: DEFAULT_MAX_DEPTH,
        }
    }
}

/// `SequenceUtil.generateAllKmers(3)`, mapped to strings.
///
/// The generator is a breadth-first walk over a `LinkedList` that ends by REMOVING the first
/// complete k-mer and adding it back at the end, so the list is lexicographic with `AAA` last
/// rather than first. That is the insertion order into the histogram maps, and therefore part of
/// what decides [`ref_context_order`].
pub fn all_kmers() -> Vec<String> {
    let mut kmers = Vec::with_capacity(64);
    for first in *BASES {
        for second in *BASES {
            for third in *BASES {
                kmers.push(String::from_utf8(vec![first, second, third]).expect("ascii"));
            }
        }
    }
    let first = kmers.remove(0);
    kmers.push(first);
    kmers
}

/// The order `refSiteHistograms.get(sample).values()` walks the 64 contexts in.
pub fn ref_context_order() -> Vec<String> {
    let entries: Vec<(String, i32)> = all_kmers()
        .into_iter()
        .map(|kmer| {
            let hash = string_hash_code(&kmer);
            (kmer, hash)
        })
        .collect();
    hash_map_order(&entries).expect("64 three-base strings do not treeify a bucket")
}

/// The order `splitBySample`'s map hands the samples over in.
pub fn sample_order(samples: &[String]) -> Vec<String> {
    let entries: Vec<(String, i32)> = samples
        .iter()
        .map(|name| (name.clone(), string_hash_code(name)))
        .collect();
    hash_map_order(&entries).expect("a handful of sample names do not treeify a bucket")
}

/// The alt counts of one sample, keyed by context, alt base and orientation.
pub type AltCounts = HashMap<(String, u8, ReadOrientation), HashMap<i32, i64>>;

/// One pileup element, reduced to what the collector reads from it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Element {
    pub sample: String,
    /// `PileupElement.getBase()`, which is `D` for a deletion and counts nowhere.
    pub base: u8,
    pub qual: u8,
    pub reverse_strand: bool,
    pub first_of_pair: bool,
    pub mapping_quality: i32,
    pub deletion: bool,
    pub after_insertion: bool,
    pub before_deletion_start: bool,
}

impl Element {
    /// `ReadUtils.isF1R2`, which for an UNPAIRED read makes a forward read F2R1 and a reverse read
    /// F1R2, `isFirstOfPair` being false when the flag is not set.
    pub fn orientation(&self) -> ReadOrientation {
        if self.reverse_strand != self.first_of_pair {
            ReadOrientation::F1R2
        } else {
            ReadOrientation::F2R1
        }
    }
}

/// One row of the alt table, as `AltSiteRecord` writes it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AltSiteRecord {
    pub context: String,
    pub ref_count: i32,
    pub alt_count: i32,
    pub ref_f1r2: i32,
    pub alt_f1r2: i32,
    pub alt: u8,
}

impl AltSiteRecord {
    /// The `depth` column, which is the two counts added and not the pileup's depth.
    pub fn depth(&self) -> i32 {
        self.ref_count + self.alt_count
    }

    /// The row as the table writer lays it out.
    pub fn line(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.context,
            self.ref_count,
            self.alt_count,
            self.ref_f1r2,
            self.alt_f1r2,
            self.depth(),
            self.alt as char
        )
    }
}

/// The counts of one run.
#[derive(Debug, Clone)]
pub struct Collector {
    args: Args,
    samples: Vec<String>,
    /// Per sample, per context, the counts by depth bin.
    ref_counts: HashMap<String, HashMap<String, HashMap<i32, i64>>>,
    /// Per sample, per (context, alt base, orientation), the counts by depth bin.
    alt_counts: HashMap<String, AltCounts>,
    alt_tables: HashMap<String, Vec<AltSiteRecord>>,
}

impl Collector {
    pub fn new(args: Args, samples: &[String]) -> Self {
        let samples = sample_order(samples);
        let mut collector = Collector {
            args,
            samples: samples.clone(),
            ref_counts: HashMap::new(),
            alt_counts: HashMap::new(),
            alt_tables: HashMap::new(),
        };
        for sample in &samples {
            collector.ref_counts.insert(sample.clone(), HashMap::new());
            collector.alt_counts.insert(sample.clone(), HashMap::new());
            collector.alt_tables.insert(sample.clone(), Vec::new());
        }
        collector
    }

    /// `F1R2CountsCollector.process` over one locus.
    ///
    /// `ref_context` is the three-base k-mer around the locus, or `None` when it ran off the
    /// reference, which is the first thing the method checks.
    pub fn process(&mut self, elements: &[Element], ref_context: Option<&str>) {
        let Some(context) = ref_context else {
            return;
        };
        if context.contains('N') || context.len() != REFERENCE_CONTEXT_SIZE {
            return;
        }
        let ref_base = context.as_bytes()[REF_CONTEXT_PADDING];

        for sample in self.samples.clone() {
            let pileup: Vec<&Element> = elements
                .iter()
                .filter(|element| {
                    element.sample == sample && element.qual as i32 > self.args.min_base_quality
                })
                .collect();
            let base_counts = base_counts(&pileup);
            let depth: i32 = base_counts.iter().sum();

            if !self.is_pileup_good(&pileup, &base_counts, depth) {
                // A return, not a continue: every sample after this one is skipped too.
                return;
            }

            let mut copy = base_counts;
            copy[base_index(ref_base).expect("a reference base")] = -1;
            let alt_index = max_element_index(&copy);
            let reference_site = base_counts[alt_index] == 0;

            if reference_site {
                let capped = depth.min(self.args.max_depth);
                *self
                    .ref_counts
                    .get_mut(&sample)
                    .expect("a known sample")
                    .entry(context.to_string())
                    .or_default()
                    .entry(capped)
                    .or_insert(0) += 1;
                return;
            }

            let alt_base = BASES[alt_index];
            let ref_count = base_counts[base_index(ref_base).expect("a reference base")];
            let alt_count = base_counts[alt_index];
            let ref_f1r2 = count_f1r2(&pileup, ref_base);
            let alt_f1r2 = count_f1r2(&pileup, alt_base);

            if alt_count == 1 {
                let orientation = if alt_f1r2 == 1 {
                    ReadOrientation::F1R2
                } else {
                    ReadOrientation::F2R1
                };
                let capped = depth.min(self.args.max_depth);
                *self
                    .alt_counts
                    .get_mut(&sample)
                    .expect("a known sample")
                    .entry((context.to_string(), alt_base, orientation))
                    .or_default()
                    .entry(capped)
                    .or_insert(0) += 1;
                return;
            }

            self.alt_tables
                .get_mut(&sample)
                .expect("a known sample")
                .push(AltSiteRecord {
                    context: context.to_string(),
                    ref_count,
                    alt_count,
                    ref_f1r2,
                    alt_f1r2,
                    alt: alt_base,
                });
        }
    }

    /// `isPileupGood`, over the already quality-filtered pileup.
    fn is_pileup_good(&self, pileup: &[&Element], base_counts: &[i32; 4], depth: i32) -> bool {
        let mapping_qualities: Vec<f64> = pileup
            .iter()
            .map(|element| element.mapping_quality as f64)
            .collect();
        // A hundredth of the depth, truncated, so at any depth below a hundred a single indel
        // element is enough.
        let indel_threshold = depth / 100;
        let indels = pileup
            .iter()
            .filter(|element| {
                element.deletion || element.after_insertion || element.before_deletion_start
            })
            .count() as i32;
        let mut is_indel = indels > indel_threshold;
        is_indel = is_indel || (depth == 0 && !pileup.is_empty());
        let _ = base_counts;
        depth > 0
            && !is_indel
            && jmath::percentile::median(&mapping_qualities) >= self.args.min_median_map_qual as f64
    }

    pub fn samples(&self) -> &[String] {
        &self.samples
    }

    pub fn alt_table(&self, sample: &str) -> &[AltSiteRecord] {
        self.alt_tables
            .get(sample)
            .map_or(&[], |records| records.as_slice())
    }

    /// The alt table file, header line included.
    pub fn alt_table_text(&self, sample: &str) -> String {
        let mut out = format!("#<METADATA>SAMPLE={sample}\n");
        out.push_str("context\tref_count\talt_count\tref_f1r2\talt_f1r2\tdepth\talt\n");
        for record in self.alt_table(sample) {
            out.push_str(&record.line());
            out.push('\n');
        }
        out
    }

    /// The reference histograms, in the order the reference walks its map.
    pub fn ref_histograms(&self, sample: &str) -> Vec<Histogram> {
        let counts = self.ref_counts.get(sample);
        ref_context_order()
            .into_iter()
            .map(|context| {
                let bins = counts
                    .and_then(|per_context| per_context.get(&context))
                    .cloned()
                    .unwrap_or_default();
                histogram(&context, &bins, self.args.max_depth)
            })
            .collect()
    }

    /// The alt histograms, SORTED BY LABEL rather than in the reference's own order, which is an
    /// identity hash order and is not reproducible.
    pub fn alt_histograms(&self, sample: &str) -> Vec<Histogram> {
        let counts = self.alt_counts.get(sample);
        let mut histograms: Vec<Histogram> = Vec::new();
        for context in all_kmers() {
            let middle = context.as_bytes()[REF_CONTEXT_PADDING];
            for alt in *BASES {
                if alt == middle {
                    continue;
                }
                for orientation in [ReadOrientation::F1R2, ReadOrientation::F2R1] {
                    let key = (context.clone(), alt, orientation);
                    let bins: HashMap<i32, i64> = counts
                        .and_then(|per_key| per_key.get(&key))
                        .cloned()
                        .unwrap_or_default();
                    let label = format!("{context}_{}_{}", alt as char, orientation.label());
                    histograms.push(histogram(&label, &bins, self.args.max_depth));
                }
            }
        }
        histograms.sort_by(|left, right| left.value_label.cmp(&right.value_label));
        histograms
    }
}

/// One prefilled histogram: every bin from one to `max_depth`, and the counts on top.
fn histogram(label: &str, counts: &HashMap<i32, i64>, max_depth: i32) -> Histogram {
    Histogram {
        bin_label: BIN_NAME.to_string(),
        value_label: label.to_string(),
        key_class: "java.lang.Integer".to_string(),
        bins: (1..=max_depth)
            .map(|bin| (bin.to_string(), *counts.get(&bin).unwrap_or(&0) as f64))
            .collect(),
    }
}

/// `ReadPileup.getBaseCounts`, which counts only the four bases.
fn base_counts(pileup: &[&Element]) -> [i32; 4] {
    let mut counts = [0; 4];
    for element in pileup {
        if let Some(index) = base_index(element.base) {
            counts[index] += 1;
        }
    }
    counts
}

/// `BaseUtils.simpleBaseToBaseIndex`, for the upper-case bases this sees.
fn base_index(base: u8) -> Option<usize> {
    BASES.iter().position(|candidate| *candidate == base)
}

/// `MathUtils.maxElementIndex`, which keeps the FIRST maximum and so settles a tie on the lower
/// base index.
fn max_element_index(counts: &[i32; 4]) -> usize {
    let mut best = 0;
    for index in 1..counts.len() {
        if counts[index] > counts[best] {
            best = index;
        }
    }
    best
}

fn count_f1r2(pileup: &[&Element], base: u8) -> i32 {
    pileup
        .iter()
        .filter(|element| element.base == base && element.orientation() == ReadOrientation::F1R2)
        .count() as i32
}

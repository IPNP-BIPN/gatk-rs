//! Mutect2's Permutect datasets, ported from `PermutectDatasetEngine` and `FeaturizedReadSets`
//! (GATK 4.7.0.0; both `.java` files are unchanged since 4.6.2.0).
//!
//! `--permutect-training-dataset` and `--permutect-test-dataset` each get an engine, and every
//! annotated call (before phasing) adds one datum per alternate it labels: the label, the site, the
//! reference bases around it, five variant features, the tumour's featurized reference and
//! alternate reads, the pre-downsampling counts, and the sequencing-error log likelihoods of the
//! tumour and the normal. When the engine closes, `contigs.table` and `read-groups.table` are
//! written beside the dataset.
//!
//! What a reader would not guess:
//!
//!  * the training engine's labels lean on `POPAF`, `TLOD`, the tumour and normal `AD` and the
//!    panel's `PON` flag, and on a queue per variant type that an artifact fills and a variant
//!    drains; a variant's alternate reads are downsampled to `Random(1).nextInt(5, 20)`;
//!  * a read set larger than its downsampling count is cut by `Collections.shuffle` WITHOUT a
//!    generator, so the reference's own output then changes from run to run. The port refuses that
//!    case rather than answer with one of the possible outputs;
//!  * the tumour samples are a `HashSet`, so with two tumours the reads come in that set's order;
//!  * `read-groups.table` lists a `HashMap` of the read group ids, in its iteration order;
//!  * the variant features are written `%.2f` from their `float` values, the log likelihoods
//!    `%.3f` from their doubles.
//!
//! # Reading this file without reading Rust
//!
//! - `Result<T, E>` is "a `T` or the error `E`", how a Java `throw` is modelled: `Ok(..)` returns,
//!   `Err(..)` throws, and a `?` after a call passes an `Err` straight up to the caller.
//!   `map_err` rewrites an error into the Java exception the reference would throw at that call.
//! - `Option<T>` is a nullable value made explicit, `Some(x)` or `None`. A featurized read is a
//!   `Vec<Option<i64>>` because the reference's `List<Integer>` can hold a `null`, and it prints
//!   as "null".
//! - `let Some(x) = y else { .. }` unwraps or runs the block, which must leave (`continue`,
//!   `return`).
//! - `&[T]` is a read-only borrowed view of a list; `&mut T` a writable borrow, so the callee may
//!   change the engine's queues and text while the caller keeps them.
//! - `.iter().map(..).collect()` is Java's stream, map and collect; `|x| ..` is a lambda.
//! - `x as f32`, `as i64` and `as usize` are casts; float to integer truncates toward zero, as
//!   Java's `(int)` does.

use gatk_annotation::assembly_complexity;
use gatk_annotation::per_allele::{BaseQuality, ReadPosition};
use gatk_annotation::read_grouping::ReferenceBases;
use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::fragment::Fragment;
use gatk_engine::haplotype::Haplotype;
use gatk_engine::java_format::format_decimals;
use gatk_engine::java_random::JavaRandom;
use gatk_engine::smith_waterman::{
    SmithWatermanAligner, SmithWatermanJavaAligner, SwOverhangStrategy,
    ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS,
};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Value, VariantContext};
use std::collections::VecDeque;

/// `PermutectDatasetEngine.CAPACITY`, each queue's bound.
const CAPACITY: usize = 100_000;
/// `RARE_POPAF_THRESHOLD`: above this `POPAF` (negative log10 frequency) a low-AF call is an
/// artifact for training.
const RARE_POPAF_THRESHOLD: f64 = 5.9;
/// `COMMON_POPAF_THRESHOLD`: below this `POPAF` a high-AF call is germline for training.
const COMMON_POPAF_THRESHOLD: f64 = 1.0;
/// `TLOD_THRESHOLD`: below this `TLOD` a call is a sequencing error, neither artifact nor germline.
const TLOD_THRESHOLD: f64 = 6.0;
/// `FeaturizedReadSets.DEFAULT_BASE_QUALITY`.
const DEFAULT_BASE_QUALITY: i32 = 25;
/// `FeaturizedReadSets.RANGES`: the distances from the variant the ranged features are counted
/// over, the first range that contains a distance taking it.
const RANGES: [i32; 3] = [5, 10, 20];
/// `FEATURES_PER_RANGE`: inserted bases, deletions, mismatches, very bad and bad qualities.
const FEATURES_PER_RANGE: usize = 5;
/// `VERY_BAD_QUAL_THRESHOLD` and `BAD_QUAL_THRESHOLD`.
const VERY_BAD_QUAL_THRESHOLD: u8 = 10;
const BAD_QUAL_THRESHOLD: u8 = 20;

/// `M2ArgumentCollection.PermutectDatasetMode`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DatasetMode {
    Illumina,
    Ultima,
}

impl DatasetMode {
    /// The enum constant's name, or `None` for anything else.
    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "ILLUMINA" => Some(DatasetMode::Illumina),
            "ULTIMA" => Some(DatasetMode::Ultima),
            _ => None,
        }
    }
}

/// `PermutectDatasetEngine.Label`; `name` is the `toString` the dataset's first line prints.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Label {
    Artifact,
    Variant,
    Unlabeled,
    Ignore,
}

impl Label {
    fn name(self) -> &'static str {
        match self {
            Label::Artifact => "ARTIFACT",
            Label::Variant => "VARIANT",
            Label::Unlabeled => "UNLABELED",
            Label::Ignore => "IGNORE",
        }
    }
}

/// What a datum refuses, in the reference's classes, or the port's own limitation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DatasetError {
    /// An exception the reference throws, by its Java class and message.
    Runtime { class: String, message: String },
    /// An input the port refuses because the reference's answer is not reproducible.
    Limitation(String),
}

/// A `DatasetError::Runtime`, shortened. `impl Into<String>` accepts a `&str` or a `String`.
fn runtime(class: &str, message: impl Into<String>) -> DatasetError {
    DatasetError::Runtime {
        class: class.to_string(),
        message: message.into(),
    }
}

/// The arguments every engine shares.
#[derive(Debug, Clone)]
pub struct DatasetArguments {
    /// `--permutect-ref-downsample` and `--permutect-alt-downsample`.
    pub max_ref_count: usize,
    pub max_alt_count: usize,
    /// `--permutect-non-artifact-ratio`.
    pub non_artifact_ratio: i32,
    /// `--permutect-dataset-mode`.
    pub mode: DatasetMode,
}

/// One `PermutectDatasetEngine`: its dataset's text so far.
///
/// `queues` is `unmatchedArtifactAltCounts`, one `ArrayBlockingQueue` per `VariantType` (SNV,
/// insertion, deletion, the order `kind` indexes them in). `random` is the engine's
/// `new Random(1)`, so each engine draws its own stream.
#[derive(Debug, Clone)]
pub struct DatasetEngine {
    pub training: bool,
    /// `--permutect-training-truth` or `--permutect-test-truth`, read whole.
    pub truth: Option<Vec<VariantContext>>,
    queues: [VecDeque<i32>; 3],
    random: JavaRandom,
    pub text: String,
}

impl DatasetEngine {
    /// The constructor's state: empty queues, `Random(1)`, no text yet.
    pub fn new(training: bool, truth: Option<Vec<VariantContext>>) -> Self {
        DatasetEngine {
            training,
            truth,
            queues: [VecDeque::new(), VecDeque::new(), VecDeque::new()],
            random: JavaRandom::new(1),
            text: String::new(),
        }
    }
}

/// Both engines, and what they share.
#[derive(Debug, Clone)]
pub struct Datasets {
    pub arguments: DatasetArguments,
    /// The training engine, then the test engine, each that was asked for.
    pub engines: Vec<DatasetEngine>,
    pub normal_samples: Vec<String>,
    /// The header's read groups, by index.
    pub read_groups: Vec<String>,
    /// The sequence dictionary's contigs, in order.
    pub contigs: Vec<String>,
}

impl Datasets {
    /// `close()`'s two tables. `contigs.table` is each dictionary contig and its index, one per
    /// line, in dictionary order.
    pub fn contigs_table(&self) -> String {
        self.contigs
            .iter()
            .enumerate()
            .map(|(index, contig)| format!("{contig}\t{index}\n"))
            .collect()
    }

    /// `read-groups.table`: `readGroupIndices`, a `new HashMap<>()` of id to header position,
    /// printed in its iteration order.
    ///
    /// ```java
    /// for (final Map.Entry<String, Integer> entry : readGroupIndices.entrySet()) {
    /// ```
    ///
    /// `hash_set_order` gives that order for the ids; if it cannot (an order it has not
    /// measured), the header's order is used instead.
    pub fn read_groups_table(&self) -> String {
        let order = gatk_engine::java_hash::hash_set_order(&self.read_groups)
            .unwrap_or_else(|_| self.read_groups.clone());
        order
            .iter()
            .map(|id| {
                let index = self.read_groups.iter().position(|g| g == id).unwrap_or(0);
                format!("{id}\t{index}\n")
            })
            .collect()
    }

    /// Every engine's `addData`, the training engine first.
    ///
    /// `vc` is the annotated call, `numbers` the same call before annotation (where `POPAF`,
    /// `TLOD` and `NALOD` still hold the doubles the reference reads), `reads` the read-allele
    /// likelihoods and `fragments` the fragment-haplotype ones.
    #[allow(clippy::too_many_arguments)]
    pub fn add_data(
        &mut self,
        reference_window: (i64, &[u8]),
        vc: &VariantContext,
        numbers: &VariantContext,
        reads: &AlleleLikelihoods<BamRecord, Allele>,
        fragments: &AlleleLikelihoods<Fragment, Haplotype>,
        haplotypes: &[&Haplotype],
        contig_index: usize,
    ) -> Result<(), DatasetError> {
        // Copies, so the loop can borrow each engine writably while these are read.
        let arguments = self.arguments.clone();
        let normal_samples = self.normal_samples.clone();
        let read_groups = self.read_groups.clone();
        for engine in &mut self.engines {
            // `featureContext.getValues(truth, mergedVC.getStart())`.
            let truth: Option<Vec<VariantContext>> = engine.truth.as_ref().map(|records| {
                records
                    .iter()
                    .filter(|t| t.contig == vc.contig && t.start <= vc.start && t.stop >= vc.start)
                    .cloned()
                    .collect()
            });
            let truth_refs: Option<Vec<&VariantContext>> =
                truth.as_ref().map(|records| records.iter().collect());
            add_datum(
                engine,
                &arguments,
                &normal_samples,
                &read_groups,
                reference_window,
                vc,
                numbers,
                truth_refs.as_deref(),
                reads,
                fragments,
                haplotypes,
                contig_index,
            )?;
        }
        Ok(())
    }
}

/// `getAttributeAsDoubleArray` / `getAttributeAsDoubleList(key, 0)`.
///
/// A missing key is an empty list; a string that does not parse is 0.0. The callers then read a
/// missing index as 0.0 where the reference's array or list would throw.
fn doubles(vc: &VariantContext, key: &str) -> Vec<f64> {
    vc.attributes
        .iter()
        .find(|(name, _)| name == key)
        .map(|(_, value)| match value {
            Value::List(values) => values.iter().map(double_of).collect(),
            other => vec![double_of(other)],
        })
        .unwrap_or_default()
}

/// One attribute value as a double, 0.0 for anything that is not a number.
fn double_of(value: &Value) -> f64 {
    match value {
        Value::Double(d) => *d,
        Value::Int(i) => *i as f64,
        Value::Str(s) => s.parse().unwrap_or(0.0),
        _ => 0.0,
    }
}

/// `sumADsOverSamples`: each allele's `AD` summed over the named samples.
///
/// ```java
/// .forEach(ad -> new IndexRange(0, vc.getNAlleles()).forEach(n -> ADs[n] += ad[n]));
/// ```
///
/// A genotype without `AD` is skipped and a short one adds what it has, where the reference would
/// throw on the `null` or the index.
fn sum_ads(vc: &VariantContext, samples: &[String]) -> Vec<i32> {
    let mut ads = vec![0; vc.alleles.len()];
    for genotype in vc.genotypes.iter() {
        if !samples.contains(&genotype.sample_name) {
            continue;
        }
        if let Some(ad) = &genotype.ad {
            for (n, value) in ad.iter().enumerate().take(ads.len()) {
                ads[n] += *value;
            }
        }
    }
    ads
}

/// `Random.nextInt(origin, bound)`: `RandomSupport.boundedNextInt`.
///
/// A power-of-two range masks one draw; any other range takes the draw's top 31 bits modulo the
/// range and redraws while that lands in the final partial block. The overflow test is Java's
/// `int` arithmetic, so it uses `wrapping_add` and `wrapping_sub`, which wrap instead of
/// panicking. `(r as u32) >> 1` is Java's `r >>> 1`.
fn bounded_next_int(random: &mut JavaRandom, origin: i32, bound: i32) -> i32 {
    let mut r = random.next_int();
    let n = bound - origin;
    let m = n - 1;
    if n & m == 0 {
        return (r & m) + origin;
    }
    let mut u = ((r as u32) >> 1) as i32;
    loop {
        r = u % n;
        if u.wrapping_add(m).wrapping_sub(r) >= 0 {
            break;
        }
        u = ((random.next_int() as u32) >> 1) as i32;
    }
    r + origin
}

/// `PermutectDatasetEngine.addData` for one engine: label every alternate, then write one datum
/// for each alternate not labelled `IGNORE`.
///
/// The training labels follow the reference's chain (`PermutectDatasetEngine.java`):
///
/// ```java
/// if  (probableArtifact) {
///     if (unmatchedQueue.size() > 0.9 * CAPACITY) { // this should rarely come up
/// ```
///
/// An artifact pushes `nonArtifactPerArtifact` copies of its alt count onto its type's queue; a
/// true variant is labelled only while that queue is non-empty, pops one, and draws its alt
/// downsampling count from the engine's `Random`.
#[allow(clippy::too_many_arguments)]
fn add_datum(
    engine: &mut DatasetEngine,
    arguments: &DatasetArguments,
    normal_samples: &[String],
    read_groups: &[String],
    reference_window: (i64, &[u8]),
    vc: &VariantContext,
    numbers: &VariantContext,
    truth: Option<&[&VariantContext]>,
    reads: &AlleleLikelihoods<BamRecord, Allele>,
    fragments: &AlleleLikelihoods<Fragment, Haplotype>,
    haplotypes: &[&Haplotype],
    contig_index: usize,
) -> Result<(), DatasetError> {
    let ref_bases = ReferenceBases::local_bases(reference_window.0, reference_window.1, vc);
    let ref_allele = vc.alleles[0].base_string();
    let position = vc.start;
    // `Collectors.toSet()` over the likelihoods' samples.
    let likelihood_samples: Vec<String> = (0..reads.number_of_samples())
        .filter_map(|s| reads.get_sample(s).cloned())
        .collect();
    let tumor_in_order: Vec<String> = likelihood_samples
        .iter()
        .filter(|s| !normal_samples.contains(s))
        .cloned()
        .collect();
    // That set's order, which is the order the tumour reads are gathered in; if it cannot be
    // reproduced the likelihoods' order stands.
    let tumor_samples =
        gatk_engine::java_hash::hash_set_order(&tumor_in_order).unwrap_or(tumor_in_order);
    let num_alt = vc.alleles.len() - 1;
    // The doubles the reference's call still holds, which the annotated record here carries as
    // written: read from the call before its annotations.
    let popafs = doubles(numbers, "POPAF");
    let tumor_lods = doubles(numbers, "TLOD");
    let blacklisted = vc.attributes.iter().any(|(name, _)| name == "PON");
    let tumor_ads = sum_ads(vc, &tumor_samples);
    let normal_ads = sum_ads(vc, normal_samples);
    let tumor_depth: i32 = tumor_ads.iter().sum();
    let normal_depth: i32 = normal_ads.iter().sum();
    let has_normal = normal_depth > 0;

    // The longest reference among the call's and the truth records', the first of a length:
    // `sorted(Comparator.comparingInt(Allele::length).reversed()).findFirst()` is a stable sort,
    // so among equal lengths the earliest wins, which strict `>` reproduces.
    let mut longest = &vc.alleles[0];
    for record in truth.unwrap_or(&[]) {
        if record.alleles[0].len() > longest.len() {
            longest = &record.alleles[0];
        }
    }
    // `ReferenceConfidenceVariantContextMerger.remapAlleles(record, longestRef)` as a closure,
    // its refusal turned into the exception it would be.
    let remap = |record: &VariantContext| {
        crate::reference_confidence_merger::remap_alleles(record, longest)
            .map_err(|e| runtime("java.lang.IllegalStateException", format!("{e:?}")))
    };
    let remapped_alts: Vec<Allele> = remap(vc)?.into_iter().skip(1).collect();
    // `truthAlleles`, a `Set` upstream, only asked `contains`, so a list does the same.
    let mut truth_alleles: Vec<Allele> = Vec::new();
    for record in truth.unwrap_or(&[]) {
        if record.is_filtered() {
            continue;
        }
        truth_alleles.extend(remap(record)?);
    }

    let mut labels = Vec::with_capacity(num_alt);
    // `altDownsampleMap`, a `HashMap` only looked up by allele, so a list of pairs serves.
    let mut alt_downsample: Vec<(Allele, usize)> = Vec::new();
    for n in 0..num_alt {
        let tumor_af = f64::from(tumor_ads[n + 1]) / f64::from(tumor_depth);
        let normal_af = if has_normal {
            f64::from(normal_ads[n + 1]) / f64::from(normal_depth)
        } else {
            0.0
        };
        let alt = &vc.alleles[n + 1];
        let remapped = &remapped_alts[n];
        // `VariantType`: SNV for equal lengths (any same-length change), insertion when the
        // alternate is longer, deletion when shorter.
        let diff = alt.base_string().len() as i64 - ref_allele.len() as i64;
        let kind = match diff.signum() {
            0 => 0,
            1 => 1,
            _ => 2,
        };
        let tlod = tumor_lods.get(n).copied().unwrap_or(0.0);
        let popaf = popafs.get(n).copied().unwrap_or(0.0);
        if engine.training {
            let likely_seq_error = tlod < TLOD_THRESHOLD;
            let definite_germline = !likely_seq_error
                && popaf < COMMON_POPAF_THRESHOLD
                && tumor_af > 0.35
                && (!has_normal || normal_af > 0.35);
            let true_variant = if truth.is_some() {
                truth_alleles.contains(remapped)
            } else {
                definite_germline
            };
            let probable_artifact = !likely_seq_error
                && if truth.is_some() {
                    !truth_alleles.contains(remapped)
                } else {
                    tumor_af < 0.2 && popaf > RARE_POPAF_THRESHOLD
                };
            let queue = &mut engine.queues[kind];
            if probable_artifact {
                if queue.len() as f64 > 0.9 * CAPACITY as f64 {
                    labels.push(Label::Ignore);
                } else if tumor_af < 0.3 {
                    labels.push(Label::Artifact);
                    // `addAll` into a bounded queue throws once it is full.
                    let copies = arguments.non_artifact_ratio.max(0) as usize;
                    if queue.len() + copies > CAPACITY {
                        return Err(runtime("java.lang.IllegalStateException", "Queue full"));
                    }
                    for _ in 0..copies {
                        queue.push_back(tumor_ads[n + 1]);
                    }
                } else {
                    labels.push(Label::Unlabeled);
                }
            } else if true_variant && !queue.is_empty() {
                labels.push(Label::Variant);
                // `unmatchedQueue.poll(); // just to pop off a value`, then the draw.
                queue.pop_front();
                let count = bounded_next_int(&mut engine.random, 5, 20);
                alt_downsample.push((alt.clone(), count as usize));
            } else if tlod > 5.0 && blacklisted {
                labels.push(Label::Artifact);
            } else if tlod > 4.0 && tumor_af < 0.3 {
                labels.push(Label::Unlabeled);
            } else {
                labels.push(Label::Ignore);
            }
        } else if truth.is_some() {
            labels.push(if truth_alleles.contains(remapped) {
                Label::Variant
            } else {
                Label::Artifact
            });
        } else {
            labels.push(Label::Unlabeled);
        }
    }
    // Nothing to write: return before the assembly complexity and the read vectors are built.
    if engine.training && labels.iter().all(|l| *l == Label::Ignore) {
        return Ok(());
    }

    let complexity = assembly_complexity::annotate_with_event_maps(
        vc, fragments, haplotypes, false,
    )
    .map_err(|e| match e {
        assembly_complexity::ComplexityError::NoCarrier => {
            runtime("java.util.NoSuchElementException", "No value present")
        }
        assembly_complexity::ComplexityError::NullPointer => {
            runtime("java.lang.NullPointerException", "null")
        }
        assembly_complexity::ComplexityError::NoReference => runtime(
            "java.lang.IndexOutOfBoundsException",
            "no reference haplotype",
        ),
    })?;
    // Only the tumour's reads: the reference also builds the normal's, then writes `0 0` in
    // their place and never prints them. So a refusal only a normal read would raise (an
    // unpaired normal read in Illumina mode, at `getMateStart`) is not raised here.
    let tumor_reads = read_vectors(
        vc,
        &tumor_samples,
        reads,
        fragments,
        arguments,
        &alt_downsample,
        read_groups,
    )?;

    for n in 0..num_alt {
        if labels[n] == Label::Ignore {
            continue;
        }
        let tlod = tumor_lods.get(n).copied().unwrap_or(0.0);
        let alt = vc.alleles[n + 1].base_string();
        let total: i32 = complexity.equivalence_counts.iter().sum();
        let counts = &complexity.equivalence_counts;
        let features = [
            if counts.len() < 2 {
                0.0
            } else {
                f64::from(counts[1]) / f64::from(total)
            },
            if counts.len() < 3 {
                0.0
            } else {
                f64::from(counts[2]) / f64::from(total)
            },
            f64::from(complexity.edit_distances[n]),
            complexity.dominance[n],
            if tumor_ads[n + 1] == 0 {
                0.0
            } else {
                tlod / f64::from(tumor_ads[n + 1])
            },
        ];
        let text = &mut engine.text;
        text.push_str(labels[n].name());
        text.push('\n');
        text.push_str(&format!("{contig_index}:{position},{ref_allele}->{alt}\n"));
        text.push_str(&ref_bases);
        text.push('\n');
        // `String.format("%.2f", x.floatValue())`: each feature is rounded to `float` first.
        let featured: Vec<String> = features
            .iter()
            .map(|x| format_decimals(f64::from(*x as f32), 2))
            .collect();
        text.push_str(&featured.join(" "));
        text.push('\n');
        let refs = &tumor_reads[0];
        let alts = &tumor_reads[n + 1];
        text.push_str(&format!("{} {} 0 0\n", refs.len(), alts.len()));
        for read in refs.iter().chain(alts.iter()) {
            let fields: Vec<String> = read
                .iter()
                .map(|v| v.map_or("null".to_string(), |v| v.to_string()))
                .collect();
            text.push_str(&fields.join(" "));
            text.push('\n');
        }
        text.push_str(&format!(
            "{} {} {} {}\n",
            tumor_depth,
            tumor_ads[n + 1],
            normal_depth,
            normal_ads[n + 1]
        ));
        // `-MathUtils.log10ToLog(tlod)` minus the natural log of `tumorDepth + 1`, with
        // `jmath::math::log` standing for Java's `Math.log`.
        let log10_to_log = gatk_engine::allele_likelihoods::log10_to_log;
        let seq_error = -log10_to_log(tlod) - jmath::math::log(f64::from(tumor_depth) + 1.0);
        text.push_str(&format_decimals(seq_error, 3));
        text.push('\n');
        let nalod = if normal_samples.is_empty() {
            0.0
        } else {
            doubles(numbers, "NALOD").get(n).copied().unwrap_or(0.0)
        };
        let normal_seq_error =
            -log10_to_log(nalod) - jmath::math::log(f64::from(normal_depth) + 1.0);
        text.push_str(&format_decimals(normal_seq_error, 3));
        text.push('\n');
    }
    Ok(())
}

/// One featurized read: `None` where the reference adds a null.
///
/// The read group index (`readGroupIndices.get` of an unknown group) and Ultima's `si` are the
/// two entries that can be `null`.
type ReadVector = Vec<Option<i64>>;

/// `FeaturizedReadSets.getReadVectors` for the tumour: in the call's allele order, each allele's
/// informative best reads over the samples, downsampled, then featurized.
///
/// The downsampling is where the port stops:
///
/// ```java
/// if (readsByAllele.get(allele).size() > downsample) {
///     Collections.shuffle(readsByAllele.get(allele));
/// ```
///
/// `Collections.shuffle` with no `Random` draws from a generator nobody seeds, so any set over its
/// count refuses as a `Limitation` instead of answering one of its outputs.
fn read_vectors(
    vc: &VariantContext,
    samples: &[String],
    reads: &AlleleLikelihoods<BamRecord, Allele>,
    fragments: &AlleleLikelihoods<Fragment, Haplotype>,
    arguments: &DatasetArguments,
    alt_downsample: &[(Allele, usize)],
    read_groups: &[String],
) -> Result<Vec<Vec<ReadVector>>, DatasetError> {
    let allele_count = reads.number_of_alleles();
    let mut by_allele: Vec<Vec<BamRecord>> = vec![Vec::new(); allele_count];
    for sample in samples {
        let Some(s) = reads.index_of_sample(sample) else {
            continue;
        };
        for best in reads.best_alleles_breaking_ties_for_sample(s, None) {
            if !best.is_informative() {
                continue;
            }
            let Some(index) = best.allele.as_ref().and_then(|a| reads.index_of_allele(a)) else {
                continue;
            };
            let evidence = reads.sample_evidence(s).unwrap_or(&[]);
            by_allele[index].push(evidence[best.evidence_index].clone());
        }
    }
    for (index, reads_of_allele) in by_allele.iter().enumerate() {
        let allele = reads.get_allele(index).expect("an allele");
        let downsample = if allele.is_reference() {
            arguments.max_ref_count
        } else {
            alt_downsample
                .iter()
                .find(|(a, _)| a == allele)
                .map_or(arguments.max_alt_count, |(_, n)| *n)
        };
        if reads_of_allele.len() > downsample {
            return Err(DatasetError::Limitation(format!(
                "a Permutect read set of {} reads above its downsampling count of {downsample} is \
                 cut by Collections.shuffle without a generator, which the reference does not \
                 reproduce from run to run.",
                reads_of_allele.len()
            )));
        }
    }
    // `bestHaplotypes`: each fragment's reads, the later put winning.
    //
    // ```java
    // .forEach(ba -> ba.evidence.getReads().forEach(read -> bestHaplotypes.put(read, ba.allele)));
    // ```
    //
    // A `HashMap<GATKRead, Haplotype>` only looked up, so a list searched by equality serves.
    let mut best_haplotypes: Vec<(BamRecord, Haplotype)> = Vec::new();
    for sample in samples {
        let Some(s) = fragments.index_of_sample(sample) else {
            continue;
        };
        let evidence = fragments.sample_evidence(s).unwrap_or(&[]);
        for best in fragments.best_alleles_breaking_ties_for_sample(s, None) {
            let Some(haplotype) = best.allele.clone() else {
                continue;
            };
            for read in &evidence[best.evidence_index].reads {
                match best_haplotypes.iter_mut().find(|(r, _)| r == read) {
                    Some(entry) => entry.1 = haplotype.clone(),
                    None => best_haplotypes.push((read.clone(), haplotype.clone())),
                }
            }
        }
    }
    // In the call's allele order, which need not be the likelihoods'.
    let mut out = Vec::new();
    for allele in &vc.alleles {
        let Some(index) = reads.index_of_allele(allele) else {
            out.push(Vec::new());
            continue;
        };
        let mut vectors = Vec::new();
        for read in &by_allele[index] {
            vectors.push(featurize(
                read,
                vc,
                &best_haplotypes,
                arguments.mode,
                read_groups,
            )?);
        }
        out.push(vectors);
    }
    Ok(out)
}

/// A read's tag by its two-letter name. `'r` says the answer borrows from the read.
fn tag<'r>(read: &'r BamRecord, name: &[u8; 2]) -> Option<&'r htsjdk_bam::tag::TagValue> {
    read.tags.get(htsjdk_bam::tag::Tag::new(name))
}

/// `FeaturizedReadSets.featurize`.
///
/// In order: read group index, mapping quality, base quality at the variant (25 when absent),
/// first of pair, reverse strand, the variant's position in the read and the distance to the
/// read's end (hard clips counted), the absolute fragment length; for Illumina the distances from
/// the fragment's ends, for Ultima `si` and `1000 * rq`; then fifteen ranged features from the read
/// re-aligned to its best haplotype, or fifteen zeros when the read has none.
fn featurize(
    read: &BamRecord,
    vc: &VariantContext,
    best_haplotypes: &[(BamRecord, Haplotype)],
    mode: DatasetMode,
    read_groups: &[String],
) -> Result<ReadVector, DatasetError> {
    use gatk_engine::read as gread;
    let mut result: ReadVector = Vec::new();
    let group = match tag(read, b"RG") {
        Some(htsjdk_bam::tag::TagValue::Str(id)) => read_groups.iter().position(|g| g == id),
        _ => None,
    };
    result.push(group.map(|g| g as i64));
    result.push(Some(i64::from(read.mapping_quality)));
    result.push(Some(i64::from(
        BaseQuality::base_quality(read, vc).unwrap_or(DEFAULT_BASE_QUALITY),
    )));
    result.push(Some(i64::from(gread::is_first_of_pair(read))));
    result.push(Some(i64::from(read.flags & 0x10 != 0)));
    let position = ReadPosition::position(read, vc).unwrap_or(0);
    result.push(Some(i64::from(position)));
    let hard_clips: u32 = read
        .cigar
        .elements
        .iter()
        .filter(|e| e.op == htsjdk_bam::cigar::Op::H)
        .map(|e| e.length)
        .sum();
    let actual_length = read.read_bases.len() as i64 + i64::from(hard_clips);
    result.push(Some(actual_length - i64::from(position)));
    let fragment_length = i64::from(gread::fragment_length(read)).abs();
    result.push(Some(fragment_length));
    if mode == DatasetMode::Illumina {
        // `getMateStart()`, which an unpaired read refuses.
        if !gread::is_paired(read) {
            return Err(runtime(
                "java.lang.IllegalStateException",
                "Cannot get mate information for an unpaired read",
            ));
        }
        let mate_start = i64::from(gatk_engine::read_utils::mate_start(read));
        let fragment_start =
            mate_start.min(i64::from(gatk_engine::read_utils::unclipped_start(read)));
        let fragment_end = fragment_start + fragment_length;
        result.push(Some(vc.start - fragment_start));
        result.push(Some(fragment_end - vc.stop));
    }
    if mode == DatasetMode::Ultima {
        result.push(match tag(read, b"si") {
            Some(htsjdk_bam::tag::TagValue::Int(v)) => Some(*v),
            Some(htsjdk_bam::tag::TagValue::Str(s)) => s.trim().parse().ok(),
            _ => None,
        });
        let rq = match tag(read, b"rq") {
            Some(htsjdk_bam::tag::TagValue::Float(v)) => *v,
            Some(htsjdk_bam::tag::TagValue::Str(s)) => s.trim().parse().unwrap_or(0.0),
            // Java 17's helpful message for the unboxing of a missing tag.
            _ => {
                return Err(runtime(
                    "java.lang.NullPointerException",
                    "Cannot invoke \"java.lang.Float.floatValue()\" because the return value of \
                     \"org.broadinstitute.hellbender.utils.read.GATKRead.getAttributeAsFloat(String)\" \
                     is null",
                ))
            }
        };
        result.push(Some((1000.0 * rq) as i64));
    }
    let Some((_, haplotype)) = best_haplotypes.iter().find(|(r, _)| r == read) else {
        result.extend(std::iter::repeat_n(
            Some(0),
            FEATURES_PER_RANGE * RANGES.len(),
        ));
        return Ok(result);
    };
    let haplotype_bases = haplotype.bases();
    // `aligner.align(haplotypeBases, read.getBases(), ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS,
    // SWOverhangStrategy.SOFTCLIP)`, on the Java aligner the reference names.
    let alignment = SmithWatermanJavaAligner
        .align(
            &haplotype_bases,
            &read.read_bases,
            &ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS,
            SwOverhangStrategy::SoftClip,
        )
        .map_err(|e| runtime("java.lang.IllegalArgumentException", format!("{e:?}")))?;
    let mut copy = read.clone();
    copy.cigar = alignment.cigar.clone();
    let start_in_haplotype = alignment.alignment_offset;
    let mut machine = gatk_engine::alignment_state::AlignmentStateMachine::new(&copy);
    let malformed = |e: gatk_engine::alignment_state::MalformedRead| {
        runtime("htsjdk.samtools.SAMFormatException", format!("{e:?}"))
    };
    // Walk the copy, re-cigared against the haplotype, one genome position at a time.
    machine.step_forward_on_genome().map_err(malformed)?;
    let mut ranged = [[0i64; FEATURES_PER_RANGE]; 3];
    while !machine.is_right_edge() {
        let Some(pe) = gatk_engine::pileup::PileupElement::from_state(&copy, &machine) else {
            break;
        };
        // `IntStream.range(0, RANGES.size()).filter(n -> distanceFromVariant <= RANGES.get(n))
        // .findFirst()`: `position` is the first range that holds the distance, `None` past 20.
        let distance = (machine.read_offset() - position).abs();
        if let Some(range) = RANGES.iter().position(|&r| distance <= r) {
            let features = &mut ranged[range];
            if pe.is_before_insertion() {
                features[0] += i64::from(pe.length_of_immediately_following_indel());
            }
            if pe.is_deletion() {
                features[1] += 1;
            } else {
                let base = pe.base();
                let qual = pe.qual();
                // `haplotypeBases[asm.getGenomeOffset() + readStartInHaplotype]`, whose
                // out-of-range index is the reference's `ArrayIndexOutOfBoundsException`.
                let index = (machine.genome_offset() + start_in_haplotype) as usize;
                let haplotype_base = *haplotype_bases.get(index).ok_or_else(|| {
                    runtime(
                        "java.lang.ArrayIndexOutOfBoundsException",
                        format!(
                            "Index {index} out of bounds for length {}",
                            haplotype_bases.len()
                        ),
                    )
                })?;
                if base != haplotype_base {
                    features[2] += 1;
                }
                if qual < VERY_BAD_QUAL_THRESHOLD {
                    features[3] += 1;
                } else if qual < BAD_QUAL_THRESHOLD {
                    features[4] += 1;
                }
            }
        }
        machine.step_forward_on_genome().map_err(malformed)?;
    }
    for features in ranged {
        result.extend(features.iter().map(|v| Some(*v)));
    }
    Ok(result)
}

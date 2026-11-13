//! `CanonicalSVLinkage` and `SVClusterEngine`, ported from GATK 4.7.0.0.
//!
//! Which structural variants are judged the same event. Reading the VCF and collapsing a cluster
//! into a representative record are not ported; deciding which records belong together is.
//!
//! # The parameter set is chosen by the pair, not by the user
//!
//! ```java
//! if (evidenceParams.isValidPair(a, b))      { ... evidenceParams ... }
//! else if (mixedParams.isValidPair(a, b))    { ... mixedParams ... }
//! else if (depthOnlyParams.isValidPair(a, b)){ ... depthOnlyParams ... }
//! ```
//!
//! The predicates are `!a.depthOnly && !b.depthOnly`, `a.depthOnly != b.depthOnly` and
//! `a.depthOnly && b.depthOnly`, and a record is depth-only when its `ALGORITHMS` is exactly
//! `["depth"]`. So one run applies three different thresholds to different pairs of the same file,
//! and no argument changes which.
//!
//! # The class documentation and the factory disagree
//!
//! The javadoc says a depth-only pair needs "only one of interval overlap or break-end proximity".
//! Every factory passes `true`:
//!
//! ```java
//! public static ClusteringParameters createDepthParameters(...) {
//!     return new ClusteringParameters(..., true, (a,b) -> a.isDepthOnly() && b.isDepthOnly());
//! }
//! ```
//!
//! `requiresOverlapAndProximity` is therefore true for all three sets, and the `||` branch of
//! `clusterTogetherWithParams` is unreachable through this linkage. Measured: widening the PESR
//! breakend window twentyfold changes nothing, because proximity alone never suffices.
//!
//! # Single linkage and max clique disagree on a chain
//!
//! And the chain need not be geometric. A deletion and a duplication at one locus do not cluster
//! with each other; both cluster with a CNV between them. Single linkage returns one group of
//! three; max clique returns two overlapping pairs that each contain the CNV. `--enable-cnv` is
//! only observable under max clique, because under single linkage the CNV already chains them.
//!
//! # An insertion is given two different assumed lengths
//!
//! ```java
//! public static final int INSERTION_ASSUMED_LENGTH_FOR_OVERLAP = 50;
//! public static final int INSERTION_ASSUMED_LENGTH_FOR_SIZE_SIMILARITY = 1;
//! ```
//!
//! Fifty for the reciprocal overlap and one for the size similarity, so the same record is a
//! different size depending on which test is asking.

use crate::sv_stratify::SvType;

/// `CanonicalSVLinkage.INSERTION_ASSUMED_LENGTH_FOR_OVERLAP`.
pub const INSERTION_ASSUMED_LENGTH_FOR_OVERLAP: i32 = 50;
/// `CanonicalSVLinkage.INSERTION_ASSUMED_LENGTH_FOR_SIZE_SIMILARITY`.
pub const INSERTION_ASSUMED_LENGTH_FOR_SIZE_SIMILARITY: i32 = 1;
/// `GATKSVVCFConstants.DEPTH_ALGORITHM`.
pub const DEPTH_ALGORITHM: &str = "depth";

/// Which pairs a parameter set applies to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairKind {
    /// Both records depth-only.
    Depth,
    /// Exactly one of them depth-only.
    Mixed,
    /// Neither depth-only.
    Pesr,
}

/// `ClusteringParameters`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClusteringParameters {
    pub reciprocal_overlap: f64,
    pub size_similarity: f64,
    pub window: i32,
    pub sample_overlap: f64,
    /// True in every factory, whatever the class documentation says.
    pub requires_overlap_and_proximity: bool,
    pub kind: PairKind,
}

impl ClusteringParameters {
    pub fn depth(
        reciprocal_overlap: f64,
        size_similarity: f64,
        window: i32,
        sample_overlap: f64,
    ) -> Self {
        ClusteringParameters {
            reciprocal_overlap,
            size_similarity,
            window,
            sample_overlap,
            requires_overlap_and_proximity: true,
            kind: PairKind::Depth,
        }
    }

    pub fn mixed(
        reciprocal_overlap: f64,
        size_similarity: f64,
        window: i32,
        sample_overlap: f64,
    ) -> Self {
        ClusteringParameters {
            reciprocal_overlap,
            size_similarity,
            window,
            sample_overlap,
            requires_overlap_and_proximity: true,
            kind: PairKind::Mixed,
        }
    }

    pub fn pesr(
        reciprocal_overlap: f64,
        size_similarity: f64,
        window: i32,
        sample_overlap: f64,
    ) -> Self {
        ClusteringParameters {
            reciprocal_overlap,
            size_similarity,
            window,
            sample_overlap,
            requires_overlap_and_proximity: true,
            kind: PairKind::Pesr,
        }
    }

    /// `isValidPair`.
    pub fn is_valid_pair(&self, a: &CallRecord, b: &CallRecord) -> bool {
        match self.kind {
            PairKind::Depth => a.is_depth_only() && b.is_depth_only(),
            PairKind::Mixed => a.is_depth_only() != b.is_depth_only(),
            PairKind::Pesr => !a.is_depth_only() && !b.is_depth_only(),
        }
    }
}

/// The three defaults, as `CanonicalSVLinkage` builds them.
pub fn default_depth_parameters() -> ClusteringParameters {
    ClusteringParameters::depth(0.8, 0.0, 10_000_000, 0.0)
}

pub fn default_mixed_parameters() -> ClusteringParameters {
    ClusteringParameters::mixed(0.8, 0.0, 1000, 0.0)
}

pub fn default_pesr_parameters() -> ClusteringParameters {
    ClusteringParameters::pesr(0.5, 0.0, 500, 0.0)
}

/// One record, reduced to what the linkage reads off it.
#[derive(Debug, Clone, PartialEq)]
pub struct CallRecord {
    pub id: String,
    pub sv_type: SvType,
    pub contig_a: String,
    pub position_a: i32,
    pub contig_b: String,
    pub position_b: i32,
    pub strand_a: Option<bool>,
    pub strand_b: Option<bool>,
    /// `getLength()`, absent for the types that have none.
    pub length: Option<i32>,
    pub algorithms: Vec<String>,
    /// The samples carrying a non-reference allele, which sample overlap is a fraction of.
    pub carriers: Vec<String>,
    /// `getComplexSubtype()`, which two records of one type must also share.
    pub cpx_subtype: Option<String>,
    /// `getComplexEventIntervals()`, which a pair of complex events is judged on instead.
    pub cpx_intervals: Vec<crate::sv_call_record::ComplexEventInterval>,
}

impl CallRecord {
    /// `isDepthOnly`: exactly one algorithm, and it is `depth`.
    pub fn is_depth_only(&self) -> bool {
        self.algorithms.len() == 1 && self.algorithms[0] == DEPTH_ALGORITHM
    }

    /// `isSimpleCNV`: a DEL, a DUP or a CNV.
    pub fn is_simple_cnv(&self) -> bool {
        matches!(self.sv_type, SvType::Del | SvType::Dup | SvType::Cnv)
    }

    pub fn null_strands(&self) -> bool {
        self.strand_a.is_none() && self.strand_b.is_none()
    }

    pub fn is_intrachromosomal(&self) -> bool {
        self.contig_a == self.contig_b
    }

    /// `getLength(record, lengthIfMissing)`, per type.
    ///
    /// A breakend spans its two positions, at least 1 since 4.7.0.0 allows the empty one, and an
    /// interchromosomal one is 0. Every other type takes its length, at least 1, and the assumed
    /// one only when it has none; the two tests assume different ones. Measured on
    /// `SVConcordance`'s array, where two insertions 50 bp apart with SVLEN 400 and 380 paired in
    /// the reference and not here, since an assumed 50 bp made them disjoint.
    pub fn length_for(&self, assumed_if_missing: i32) -> i32 {
        match self.sv_type {
            SvType::Bnd if self.is_intrachromosomal() => (self.position_b - self.position_a).max(1),
            SvType::Bnd => 0,
            _ => self.length.unwrap_or(assumed_if_missing).max(1),
        }
    }
}

/// `typesMatch`, which lets CNVs stand in for both deletions and duplications. Two records of one
/// type must also share their complex subtype, which is null for every type but CPX.
pub fn types_match(a: &CallRecord, b: &CallRecord, cluster_del_with_dup: bool) -> bool {
    if a.sv_type == b.sv_type && a.cpx_subtype == b.cpx_subtype {
        return true;
    }
    if a.is_simple_cnv() && b.is_simple_cnv() {
        // A DEL and a DUP only meet through a CNV, unless the flag opens the door.
        if cluster_del_with_dup || a.sv_type == SvType::Cnv || b.sv_type == SvType::Cnv {
            return true;
        }
    }
    false
}

/// `strandsMatch`, where a record with no strands at all matches anything.
pub fn strands_match(a: &CallRecord, b: &CallRecord) -> bool {
    if a.null_strands() || b.null_strands() {
        return true;
    }
    a.strand_a == b.strand_a && a.strand_b == b.strand_b
}

/// `computeReciprocalOverlap(SimpleInterval, SimpleInterval)`: the shared bases over the LONGER
/// interval, and 0 for two that do not overlap.
pub fn reciprocal_overlap(start_a: i32, end_a: i32, start_b: i32, end_b: i32) -> f64 {
    let overlap = end_a.min(end_b) - start_a.max(start_b) + 1;
    if overlap <= 0 {
        return 0.0;
    }
    f64::from(overlap) / f64::from((end_a - start_a + 1).max(end_b - start_b + 1))
}

/// `testReciprocalOverlap` over [`reciprocal_overlap`]. Since 4.7.0.0 a threshold of 0 passes two
/// disjoint intervals, which `IntervalUtils.isReciprocalOverlap` did not.
pub fn is_reciprocal_overlap(
    start_a: i32,
    end_a: i32,
    start_b: i32,
    end_b: i32,
    threshold: f64,
) -> bool {
    reciprocal_overlap(start_a, end_a, start_b, end_b) >= threshold
}

/// `computeSizeSimilarity`: the smaller over the larger.
pub fn size_similarity(length_a: i32, length_b: i32) -> f64 {
    f64::from(length_a.min(length_b)) / f64::from(length_a.max(length_b))
}

/// `testSizeSimilarity` over [`size_similarity`].
pub fn test_size_similarity(length_a: i32, length_b: i32, threshold: f64) -> bool {
    size_similarity(length_a, length_b) >= threshold
}

/// `testBreakendProximity`: both distances measured, which needs both contigs shared, and both
/// within the window.
fn test_breakend_proximity(distance_1: Option<i32>, distance_2: Option<i32>, window: i32) -> bool {
    matches!((distance_1, distance_2), (Some(d1), Some(d2)) if d1 <= window && d2 <= window)
}

/// `getFirstBreakpointProximity` and `getSecondBreakpointProximity`: null across contigs.
fn distance(contig_a: &str, position_a: i32, contig_b: &str, position_b: i32) -> Option<i32> {
    (contig_a == contig_b).then(|| (position_a - position_b).abs())
}

/// `getSampleSetOverlap` and `testSampleOverlap`: the shared carriers over the LARGER set, two
/// empty sets overlapping perfectly. A threshold of 0 is never computed at all.
pub fn has_sample_overlap(a: &CallRecord, b: &CallRecord, threshold: f64) -> bool {
    if threshold <= 0.0 {
        return true;
    }
    let shared = a
        .carriers
        .iter()
        .filter(|sample| b.carriers.contains(sample))
        .count();
    let larger = a.carriers.len().max(b.carriers.len());
    if larger == 0 {
        return true;
    }
    shared as f64 / larger as f64 >= threshold
}

/// `CanonicalSVLinkage.CanonicalLinkageResult`: the decision and what was measured on the way to
/// it, which `SVConcordance` writes as its `TRUTH_*` fields. A pair refused before anything was
/// measured, and every complex pair, carries no measurement at all.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct LinkageResult {
    pub result: bool,
    pub reciprocal_overlap: Option<f64>,
    pub size_similarity: Option<f64>,
    pub breakpoint_distance_1: Option<i32>,
    pub breakpoint_distance_2: Option<i32>,
}

impl LinkageResult {
    fn refused() -> LinkageResult {
        LinkageResult::default()
    }

    fn bare(result: bool) -> LinkageResult {
        LinkageResult {
            result,
            ..LinkageResult::default()
        }
    }
}

/// The three parameter sets one run was given.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Linkage {
    pub depth: ClusteringParameters,
    pub mixed: ClusteringParameters,
    pub pesr: ClusteringParameters,
    pub cluster_del_with_dup: bool,
}

impl Default for Linkage {
    fn default() -> Self {
        Linkage {
            depth: default_depth_parameters(),
            mixed: default_mixed_parameters(),
            pesr: default_pesr_parameters(),
            cluster_del_with_dup: false,
        }
    }
}

impl Linkage {
    /// `areClusterable`.
    pub fn are_clusterable(&self, a: &CallRecord, b: &CallRecord) -> bool {
        self.linkage(a, b).result
    }

    /// `areClusterable` with what it measured.
    pub fn linkage(&self, a: &CallRecord, b: &CallRecord) -> LinkageResult {
        if !types_match(a, b, self.cluster_del_with_dup) {
            return LinkageResult::refused();
        }
        // Only these two types care about strands.
        if matches!(a.sv_type, SvType::Bnd | SvType::Inv) && !strands_match(a, b) {
            return LinkageResult::refused();
        }
        // The order is PESR, then mixed, then depth: the first set the pair is valid for wins.
        let parameters = if self.pesr.is_valid_pair(a, b) {
            self.pesr
        } else if self.mixed.is_valid_pair(a, b) {
            self.mixed
        } else if self.depth.is_valid_pair(a, b) {
            self.depth
        } else {
            return LinkageResult::refused();
        };
        cluster_together_with(a, b, &parameters)
    }

    /// `getMaxClusterableStartingPosition`: the furthest start any record linked to this one can
    /// have, under its own parameter set or the mixed one, whichever reaches further.
    pub fn max_clusterable_start(&self, record: &CallRecord, contig_length: i32) -> i32 {
        let own = if record.is_depth_only() {
            self.depth
        } else {
            self.pesr
        };
        max_clusterable_start_with(record, &own, contig_length).max(max_clusterable_start_with(
            record,
            &self.mixed,
            contig_length,
        ))
    }
}

/// `getMaxClusterableStartingPositionWithParams`.
fn max_clusterable_start_with(
    record: &CallRecord,
    parameters: &ClusteringParameters,
    contig_length: i32,
) -> i32 {
    let by_window = record
        .position_a
        .wrapping_add(parameters.window)
        .min(contig_length);
    if !record.is_intrachromosomal() {
        return by_window;
    }
    let assumed = record.length_for(INSERTION_ASSUMED_LENGTH_FOR_OVERLAP);
    let by_overlap = ((f64::from(record.position_a)
        + (1.0 - parameters.reciprocal_overlap) * f64::from(assumed)) as i32)
        .min(contig_length);
    if parameters.requires_overlap_and_proximity {
        by_overlap.min(by_window)
    } else {
        by_overlap.max(by_window)
    }
}

/// `clusterTogetherWithParams`, as 4.7.0.0 has it: every measurement is taken first, then judged.
///
/// Size similarity is now required in both branches, and an interchromosomal pair measures neither
/// overlap nor size, so each passes and the pair is judged on proximity alone.
pub fn cluster_together_with(
    a: &CallRecord,
    b: &CallRecord,
    parameters: &ClusteringParameters,
) -> LinkageResult {
    if a.contig_a != b.contig_a || a.contig_b != b.contig_b {
        return LinkageResult::refused();
    }
    if a.sv_type == SvType::Cpx && b.sv_type == SvType::Cpx {
        return test_complex_intervals(a, b, parameters);
    }

    let distance_1 = distance(&a.contig_a, a.position_a, &b.contig_a, b.position_a);
    let distance_2 = distance(&a.contig_b, a.position_b, &b.contig_b, b.position_b);
    let overlap = (a.is_intrachromosomal() && b.is_intrachromosomal()).then(|| {
        let length_a = a.length_for(INSERTION_ASSUMED_LENGTH_FOR_OVERLAP);
        let length_b = b.length_for(INSERTION_ASSUMED_LENGTH_FOR_OVERLAP);
        reciprocal_overlap(
            a.position_a,
            a.position_a + length_a - 1,
            b.position_a,
            b.position_a + length_b - 1,
        )
    });
    let size = (a.is_intrachromosomal() && b.is_intrachromosomal()).then(|| {
        size_similarity(
            a.length_for(INSERTION_ASSUMED_LENGTH_FOR_SIZE_SIMILARITY),
            b.length_for(INSERTION_ASSUMED_LENGTH_FOR_SIZE_SIMILARITY),
        )
    });

    let proximity = test_breakend_proximity(distance_1, distance_2, parameters.window);
    let has_overlap = overlap.is_none_or(|value| value >= parameters.reciprocal_overlap);
    let has_size = size.is_none_or(|value| value >= parameters.size_similarity);
    let passes = if parameters.requires_overlap_and_proximity {
        proximity && has_overlap
    } else {
        // Unreachable through CanonicalSVLinkage: every factory sets the flag.
        proximity || has_overlap
    };
    LinkageResult {
        result: passes && has_size && has_sample_overlap(a, b, parameters.sample_overlap),
        reciprocal_overlap: overlap,
        size_similarity: size,
        breakpoint_distance_1: distance_1,
        breakpoint_distance_2: distance_2,
    }
}

/// `testComplexIntervals`: the whole events' proximity, then each pair of intervals in order by
/// type, overlap, size and proximity, then the samples. It returns the bare decision, so a matched
/// complex event writes none of the `TRUTH_*` measurements.
fn test_complex_intervals(
    a: &CallRecord,
    b: &CallRecord,
    parameters: &ClusteringParameters,
) -> LinkageResult {
    let overall_1 = distance(&a.contig_a, a.position_a, &b.contig_a, b.position_a);
    let overall_2 = distance(&a.contig_b, a.position_b, &b.contig_b, b.position_b);
    if !test_breakend_proximity(overall_1, overall_2, parameters.window)
        || a.cpx_intervals.len() != b.cpx_intervals.len()
    {
        return LinkageResult::refused();
    }
    for (interval_a, interval_b) in a.cpx_intervals.iter().zip(&b.cpx_intervals) {
        if interval_a.sv_type != interval_b.sv_type {
            return LinkageResult::refused();
        }
        let distance_1 = distance(
            &interval_a.contig,
            interval_a.start,
            &interval_b.contig,
            interval_b.start,
        );
        let distance_2 = distance(
            &interval_a.contig,
            interval_a.end,
            &interval_b.contig,
            interval_b.end,
        );
        // `SimpleInterval.overlaps` asks for the same contig first.
        let overlap = if interval_a.contig == interval_b.contig {
            reciprocal_overlap(
                interval_a.start,
                interval_a.end,
                interval_b.start,
                interval_b.end,
            )
        } else {
            0.0
        };
        let size = size_similarity(
            interval_a.end - interval_a.start + 1,
            interval_b.end - interval_b.start + 1,
        );
        if !(overlap >= parameters.reciprocal_overlap
            && size >= parameters.size_similarity
            && test_breakend_proximity(distance_1, distance_2, parameters.window))
        {
            return LinkageResult::refused();
        }
    }
    LinkageResult::bare(has_sample_overlap(a, b, parameters.sample_overlap))
}

/// `SVClusterEngine.CLUSTERING_TYPE`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Algorithm {
    SingleLinkage,
    MaxClique,
}

/// The clusters one run produced, each a list of member ids in the order the items were added.
pub fn cluster(
    records: &[CallRecord],
    linkage: &Linkage,
    algorithm: Algorithm,
) -> Vec<Vec<String>> {
    cluster_indices(records, linkage, algorithm)
        .into_iter()
        .map(|members| {
            let mut ids: Vec<String> = members
                .into_iter()
                .map(|index| records[index].id.clone())
                .collect();
            // The output writes them sorted, which is what MEMBERS holds.
            ids.sort();
            ids
        })
        .collect()
}

/// [`cluster`] as positions in `records`, each cluster in its members' order, which is what a
/// collapser reads: the first of two equal genotypes is the one kept.
pub fn cluster_indices(
    records: &[CallRecord],
    linkage: &Linkage,
    algorithm: Algorithm,
) -> Vec<Vec<usize>> {
    let mut clusters: Vec<Vec<usize>> = Vec::new();

    for (index, item) in records.iter().enumerate() {
        let linked: Vec<usize> = (0..index)
            .filter(|other| linkage.are_clusterable(item, &records[*other]))
            .collect();

        let mut to_augment: Vec<usize> = Vec::new();
        let mut to_seed: Vec<Vec<usize>> = Vec::new();
        for (cluster_index, members) in clusters.iter().enumerate() {
            match algorithm {
                Algorithm::MaxClique => {
                    let linked_here: Vec<usize> = members
                        .iter()
                        .filter(|member| linked.contains(member))
                        .cloned()
                        .collect();
                    if linked_here.len() == members.len() {
                        to_augment.push(cluster_index);
                    } else if !linked_here.is_empty() && !to_seed.contains(&linked_here) {
                        to_seed.push(linked_here);
                    }
                }
                Algorithm::SingleLinkage => {
                    if members.iter().any(|member| linked.contains(member)) {
                        to_augment.push(cluster_index);
                    }
                }
            }
        }

        // Max clique seeds new clusters from the subsets that linked, dropping any seed that is
        // contained in another one being created.
        let mut seeded: Vec<Vec<usize>> = Vec::new();
        if !to_seed.is_empty() {
            let mut candidates: Vec<Vec<usize>> = to_seed.clone();
            candidates.extend(to_augment.iter().map(|index| clusters[*index].clone()));
            // A stable sort by size, as `Comparator.comparingInt(Set::size)` is.
            candidates.sort_by_key(Vec::len);
            for (position, seed) in candidates.iter().enumerate() {
                let is_subset = candidates[position + 1..]
                    .iter()
                    .any(|other| seed.iter().all(|member| other.contains(member)));
                if !is_subset {
                    let mut members = seed.clone();
                    members.push(index);
                    seeded.push(members);
                }
            }
        }

        match algorithm {
            Algorithm::SingleLinkage => {
                if !to_augment.is_empty() {
                    // Every matching cluster is merged into one, with the item added.
                    let mut merged: Vec<usize> = Vec::new();
                    for cluster_index in to_augment.iter().rev() {
                        let removed = clusters.remove(*cluster_index);
                        for member in removed {
                            if !merged.contains(&member) {
                                merged.push(member);
                            }
                        }
                    }
                    merged.sort_unstable();
                    merged.push(index);
                    clusters.push(merged);
                }
            }
            Algorithm::MaxClique => {
                for cluster_index in &to_augment {
                    clusters[*cluster_index].push(index);
                }
            }
        }
        clusters.extend(seeded);

        if to_augment.is_empty() && to_seed.is_empty() {
            clusters.push(vec![index]);
        }
    }

    clusters
}

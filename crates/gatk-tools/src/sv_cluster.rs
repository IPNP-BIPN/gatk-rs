//! `CanonicalSVLinkage` and `SVClusterEngine`, ported from GATK 4.6.2.0.
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

    /// `getLength(record, assumed)`: an insertion has no length of its own, so one is assumed, and
    /// the two tests assume different ones.
    pub fn length_for(&self, assumed_for_insertion: i32) -> i32 {
        match self.sv_type {
            SvType::Ins => assumed_for_insertion,
            _ => self.length.unwrap_or(assumed_for_insertion),
        }
    }
}

/// `typesMatch`, which lets CNVs stand in for both deletions and duplications.
pub fn types_match(a: &CallRecord, b: &CallRecord, cluster_del_with_dup: bool) -> bool {
    if a.sv_type == b.sv_type {
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

/// `IntervalUtils.isReciprocalOverlap`: each interval must cover the threshold fraction of the
/// other.
pub fn is_reciprocal_overlap(
    start_a: i32,
    end_a: i32,
    start_b: i32,
    end_b: i32,
    threshold: f64,
) -> bool {
    let overlap = (end_a.min(end_b) - start_a.max(start_b) + 1).max(0);
    if overlap == 0 {
        return false;
    }
    let length_a = end_a - start_a + 1;
    let length_b = end_b - start_b + 1;
    f64::from(overlap) / f64::from(length_a) >= threshold
        && f64::from(overlap) / f64::from(length_b) >= threshold
}

/// `testSizeSimilarity`: the smaller over the larger.
pub fn test_size_similarity(length_a: i32, length_b: i32, threshold: f64) -> bool {
    f64::from(length_a.min(length_b)) / f64::from(length_a.max(length_b)) >= threshold
}

/// `testBreakendProximity`: both ends within the window of each other, on the same contigs.
pub fn test_breakend_proximity(a: &CallRecord, b: &CallRecord, window: i32) -> bool {
    a.contig_a == b.contig_a
        && a.contig_b == b.contig_b
        && (a.position_a - b.position_a).abs() <= window
        && (a.position_b - b.position_b).abs() <= window
}

/// `hasSampleOverlap`: a fraction of the SMALLER carrier set.
pub fn has_sample_overlap(a: &CallRecord, b: &CallRecord, threshold: f64) -> bool {
    if threshold <= 0.0 {
        return true;
    }
    let shared = a
        .carriers
        .iter()
        .filter(|sample| b.carriers.contains(sample))
        .count();
    let smaller = a.carriers.len().min(b.carriers.len());
    if smaller == 0 {
        return false;
    }
    shared as f64 / smaller as f64 >= threshold
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
        if !types_match(a, b, self.cluster_del_with_dup) {
            return false;
        }
        // Only these two types care about strands.
        if matches!(a.sv_type, SvType::Bnd | SvType::Inv) && !strands_match(a, b) {
            return false;
        }
        // The order is PESR, then mixed, then depth: the first set the pair is valid for wins.
        let parameters = if self.pesr.is_valid_pair(a, b) {
            self.pesr
        } else if self.mixed.is_valid_pair(a, b) {
            self.mixed
        } else if self.depth.is_valid_pair(a, b) {
            self.depth
        } else {
            return false;
        };
        cluster_together_with(a, b, &parameters)
    }
}

/// `clusterTogetherWithParams`.
pub fn cluster_together_with(
    a: &CallRecord,
    b: &CallRecord,
    parameters: &ClusteringParameters,
) -> bool {
    if a.contig_a != b.contig_a || a.contig_b != b.contig_b {
        return false;
    }

    // Overlap is skipped entirely for an interchromosomal pair, which is then judged on proximity
    // alone whatever `requires_overlap_and_proximity` says.
    let overlap_and_size = if a.is_intrachromosomal() {
        let length_a = a.length_for(INSERTION_ASSUMED_LENGTH_FOR_OVERLAP);
        let length_b = b.length_for(INSERTION_ASSUMED_LENGTH_FOR_OVERLAP);
        let overlap = is_reciprocal_overlap(
            a.position_a,
            a.position_a + length_a - 1,
            b.position_a,
            b.position_a + length_b - 1,
            parameters.reciprocal_overlap,
        );
        let size = test_size_similarity(
            a.length_for(INSERTION_ASSUMED_LENGTH_FOR_SIZE_SIMILARITY),
            b.length_for(INSERTION_ASSUMED_LENGTH_FOR_SIZE_SIMILARITY),
            parameters.size_similarity,
        );
        let both = overlap && size;
        if parameters.requires_overlap_and_proximity && !both {
            return false;
        }
        Some(both)
    } else {
        None
    };

    let proximity = test_breakend_proximity(a, b, parameters.window);
    match overlap_and_size {
        None => proximity && has_sample_overlap(a, b, parameters.sample_overlap),
        Some(both) => {
            if parameters.requires_overlap_and_proximity {
                both && proximity && has_sample_overlap(a, b, parameters.sample_overlap)
            } else {
                // Unreachable through CanonicalSVLinkage: every factory sets the flag.
                (both || proximity) && has_sample_overlap(a, b, parameters.sample_overlap)
            }
        }
    }
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

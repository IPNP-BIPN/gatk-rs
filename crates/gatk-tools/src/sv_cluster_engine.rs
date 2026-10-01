//! `SVClusterEngine` as a stream, and the two CNV linkages, ported from GATK 4.6.2.0.
//!
//! [`crate::sv_cluster`] answers which records of a whole contig belong together, which is what
//! `SVCluster` needs because it collapses at the end of a contig. `JointGermlineCNVSegmentation`
//! drives two engines one record at a time and writes what each call to `addAndFlush` returns, so
//! the ORDER in which clusters leave the engine is part of the output. That order is the one
//! `idToClusterMap`, a `HashMap<Integer, Cluster>`, iterates in, and the members of a merged or
//! seeded cluster come out of `HashSet<Integer>` copies; both are [`JavaHashMap`] here, whose
//! layout is the measured one.
//!
//! # A cluster is closed by a position, not by a contig
//!
//! Every cluster carries the largest start an item may have and still join it. A new item past
//! that position closes the cluster, which is collapsed there and then, before the item is placed.
//! A new contig closes everything.
//!
//! # Two linkages beside the canonical one
//!
//! `CNVLinkage` pads each record by a fraction of its OWN length and joins depth-only copy-number
//! records whose padded intervals overlap; one carrier on both sides must also have lost or gained
//! the same number of copies. `BinnedCNVLinkage` pads by whole bins of the model's intervals
//! instead, and lets any item on the contig join, since its bound is the contig's length.

use gatk_engine::java_hash::{JavaHashCode, JavaHashMap};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::Genotype;

use crate::sv_cluster::{CallRecord, Linkage};
use crate::sv_collapser::{genotype_int, is_carrier, Member};
use crate::sv_stratify::SvType;

/// A refusal raised inside an engine or a linkage, as the Java class and message.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EngineError {
    pub class: &'static str,
    pub message: String,
}

impl EngineError {
    fn illegal_state(message: String) -> Self {
        EngineError {
            class: "java.lang.IllegalStateException",
            message,
        }
    }

    fn illegal_argument(message: String) -> Self {
        EngineError {
            class: "java.lang.IllegalArgumentException",
            message,
        }
    }
}

/// `SVClusterLinkage`.
pub trait EngineLinkage {
    fn are_clusterable(&self, a: &Member, b: &Member) -> Result<bool, EngineError>;
    fn max_clusterable_start(&self, item: &Member) -> Result<i32, EngineError>;
}

/// `SVClusterEngine.CLUSTERING_TYPE`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClusteringType {
    SingleLinkage,
    MaxClique,
}

/// A list of item ids as a hash key, `List.hashCode` being the one its Javadoc specifies.
#[derive(Debug, Clone, PartialEq)]
struct IdList(Vec<i32>);

impl JavaHashCode for IdList {
    fn java_hash_code(&self) -> i32 {
        let mut hash: i32 = 1;
        for id in &self.0 {
            hash = hash.wrapping_mul(31).wrapping_add(*id);
        }
        hash
    }
}

#[derive(Debug, Clone)]
struct Cluster {
    max_clusterable_start: i32,
    items: Vec<i32>,
}

/// `new HashSet<>(collection)`'s iteration order over item ids.
fn hash_set_order(items: &[i32]) -> Vec<i32> {
    JavaHashMap::<i32, ()>::copy_of(items.iter())
        .keys()
        .copied()
        .collect()
}

/// `SVClusterEngine`, one item at a time. What leaves it is each closed cluster's members, in the
/// order the engine collapses them; the caller collapses.
pub struct Engine<L: EngineLinkage> {
    clustering: ClusteringType,
    pub linkage: L,
    clusters: JavaHashMap<i32, Cluster>,
    items: Vec<Member>,
    current_contig: Option<String>,
    last_start: i32,
    next_cluster: i32,
}

impl<L: EngineLinkage> Engine<L> {
    pub fn new(clustering: ClusteringType, linkage: L) -> Self {
        Engine {
            clustering,
            linkage,
            clusters: JavaHashMap::new(),
            items: Vec::new(),
            current_contig: None,
            last_start: 0,
            next_cluster: 0,
        }
    }

    /// `addAndFlush`.
    pub fn add_and_flush(&mut self, item: Member) -> Result<Vec<Vec<Member>>, EngineError> {
        if self.current_contig.as_deref() != Some(item.call.contig_a.as_str()) {
            let result = self.flush();
            self.current_contig = Some(item.call.contig_a.clone());
            self.last_start = 0;
            let id = self.register(item)?;
            self.seed(vec![id])?;
            return Ok(result);
        }
        let id = self.register(item)?;
        let to_process = self.cluster(id)?;
        Ok(to_process
            .into_iter()
            .map(|cluster| self.process(cluster))
            .collect())
    }

    /// `flush`: every open cluster, in the map's order.
    pub fn flush(&mut self) -> Vec<Vec<Member>> {
        let ids: Vec<i32> = self.clusters.keys().copied().collect();
        let result = ids.into_iter().map(|id| self.process(id)).collect();
        self.items.clear();
        self.next_cluster = 0;
        result
    }

    fn register(&mut self, item: Member) -> Result<i32, EngineError> {
        if item.call.position_a < self.last_start {
            return Err(EngineError::illegal_state(
                "Items must be added in order of increasing start coordinate".to_string(),
            ));
        }
        self.last_start = item.call.position_a;
        self.items.push(item);
        Ok(self.items.len() as i32 - 1)
    }

    fn item(&self, id: i32) -> &Member {
        &self.items[id as usize]
    }

    fn max_start_of(&self, ids: &[i32]) -> Result<i32, EngineError> {
        let mut contigs: Vec<&str> = Vec::new();
        for id in ids {
            let contig = self.item(*id).call.contig_a.as_str();
            if !contigs.contains(&contig) {
                contigs.push(contig);
            }
        }
        if contigs.len() > 1 {
            return Err(EngineError::illegal_argument(
                "Items start on multiple contigs".to_string(),
            ));
        }
        let mut max = i32::MIN;
        for id in ids {
            max = max.max(self.linkage.max_clusterable_start(self.item(*id))?);
        }
        Ok(max)
    }

    fn seed(&mut self, items: Vec<i32>) -> Result<(), EngineError> {
        let max_clusterable_start = self.max_start_of(&items)?;
        let id = self.next_cluster;
        self.next_cluster += 1;
        self.clusters.insert(
            id,
            Cluster {
                max_clusterable_start,
                items,
            },
        );
        Ok(())
    }

    fn cluster(&mut self, id: i32) -> Result<Vec<i32>, EngineError> {
        let mut linked: Vec<i32> = Vec::new();
        for (_, cluster) in self.clusters.iter() {
            for other in &cluster.items {
                if *other != id
                    && !linked.contains(other)
                    && self
                        .linkage
                        .are_clusterable(self.item(id), self.item(*other))?
                {
                    linked.push(*other);
                }
            }
        }
        let position = self.item(id).call.position_a;
        let mut to_process = Vec::new();
        let mut to_augment = Vec::new();
        let mut to_seed: JavaHashMap<IdList, ()> = JavaHashMap::new();
        for (cluster_id, cluster) in self.clusters.iter() {
            if position > cluster.max_clusterable_start {
                to_process.push(*cluster_id);
                continue;
            }
            match self.clustering {
                ClusteringType::MaxClique => {
                    let linked_here: Vec<i32> = cluster
                        .items
                        .iter()
                        .copied()
                        .filter(|item| linked.contains(item))
                        .collect();
                    if linked_here.len() == cluster.items.len() {
                        to_augment.push(*cluster_id);
                    } else if !linked_here.is_empty() {
                        to_seed.insert(IdList(linked_here), ());
                    }
                }
                ClusteringType::SingleLinkage => {
                    if cluster.items.iter().any(|item| linked.contains(item)) {
                        to_augment.push(*cluster_id);
                    }
                }
            }
        }

        if !to_seed.is_empty() {
            let mut triggered: Vec<Vec<i32>> =
                to_seed.keys().map(|list| hash_set_order(&list.0)).collect();
            for cluster_id in &to_augment {
                let items = &self.clusters.get(cluster_id).expect("a cluster").items;
                triggered.push(hash_set_order(items));
            }
            triggered.sort_by_key(Vec::len);
            for index in 0..triggered.len() {
                let seed = &triggered[index];
                let is_subset = triggered[index + 1..]
                    .iter()
                    .any(|other| seed.iter().all(|item| other.contains(item)));
                if !is_subset {
                    let mut items = seed.clone();
                    items.push(id);
                    self.seed(items)?;
                }
            }
        }

        match self.clustering {
            ClusteringType::SingleLinkage => {
                if !to_augment.is_empty() {
                    let mut items: Vec<i32> = Vec::new();
                    for cluster_id in &to_augment {
                        let removed = self.clusters.remove(cluster_id).expect("a cluster");
                        for item in removed.items {
                            if !items.contains(&item) {
                                items.push(item);
                            }
                        }
                    }
                    items.push(id);
                    self.seed(items)?;
                }
            }
            ClusteringType::MaxClique => {
                let start = self.linkage.max_clusterable_start(self.item(id))?;
                for cluster_id in &to_augment {
                    let cluster = self.clusters.get_mut(cluster_id).expect("a cluster");
                    cluster.items.push(id);
                    cluster.max_clusterable_start = cluster.max_clusterable_start.max(start);
                }
            }
        }
        if to_augment.is_empty() && to_seed.is_empty() {
            self.seed(vec![id])?;
        }
        Ok(to_process)
    }

    fn process(&mut self, cluster_id: i32) -> Vec<Member> {
        let cluster = self.clusters.remove(&cluster_id).expect("a cluster");
        cluster
            .items
            .iter()
            .map(|id| self.item(*id).clone())
            .collect()
    }
}

/// `SVCallRecord.isDepthOnly`.
pub fn is_depth_only(member: &Member) -> bool {
    member.call.algorithms.len() == 1
        && member.call.algorithms[0] == crate::sv_cluster::DEPTH_ALGORITHM
}

fn is_simple_cnv(member: &Member) -> bool {
    matches!(member.call.sv_type, SvType::Del | SvType::Dup | SvType::Cnv)
}

fn has_alt(member: &Member) -> bool {
    member
        .alleles
        .iter()
        .any(|allele| !allele.is_no_call() && !allele.is_reference())
}

/// `getCarrierSampleSet`, in genotype order.
fn carriers(member: &Member) -> Result<Vec<String>, EngineError> {
    let alt = has_alt(member);
    let mut out = Vec::new();
    for genotype in &member.genotypes {
        let carrier = is_carrier(member.call.sv_type, alt, genotype)
            .map_err(|error| EngineError::illegal_argument(error.message()))?;
        if carrier && !out.contains(&genotype.sample_name) {
            out.push(genotype.sample_name.clone());
        }
    }
    Ok(out)
}

fn has_extended(genotype: &Genotype, key: &str) -> bool {
    genotype.extended.iter().any(|(name, _)| name == key)
}

/// `SVClusterLinkage.getCopyState`.
fn copy_state(genotype: Option<&Genotype>, matched: Option<&Genotype>) -> Result<i32, EngineError> {
    match genotype {
        None => match matched {
            Some(matched) => Ok(genotype_int(matched, "ECN", -1)),
            None => Err(EngineError::illegal_argument(
                "Both genotypes are null".to_string(),
            )),
        },
        Some(genotype) => Ok(genotype_int(
            genotype,
            "CN",
            genotype_int(genotype, "RD_CN", -1),
        )),
    }
}

/// `SVClusterLinkage.hasSampleOverlap`.
pub fn has_sample_overlap(a: &Member, b: &Member, minimum: f64) -> Result<bool, EngineError> {
    if minimum <= 0.0 {
        return Ok(true);
    }
    if a.call.sv_type == SvType::Cnv || b.call.sv_type == SvType::Cnv {
        let mut samples: Vec<&str> = Vec::new();
        for genotype in a.genotypes.iter().chain(b.genotypes.iter()) {
            if !samples.contains(&genotype.sample_name.as_str()) {
                samples.push(&genotype.sample_name);
            }
        }
        if samples.is_empty() {
            return Ok(true);
        }
        let mut matches = 0;
        for sample in &samples {
            let genotype_a = a.genotypes.iter().find(|g| g.sample_name == *sample);
            let genotype_b = b.genotypes.iter().find(|g| g.sample_name == *sample);
            if copy_state(genotype_a, genotype_b)? == copy_state(genotype_b, genotype_a)? {
                matches += 1;
            }
        }
        return Ok(f64::from(matches) / samples.len() as f64 >= minimum);
    }
    let samples_a = carriers(a)?;
    let samples_b = carriers(b)?;
    let denominator = samples_a.len().max(samples_b.len());
    if denominator == 0 {
        return Ok(true);
    }
    let shared = samples_a.iter().filter(|s| samples_b.contains(s)).count();
    Ok(shared as f64 / denominator as f64 >= minimum)
}

/// `SVCallRecordUtils.sortAlleles`: by display string.
fn sorted_alleles(alleles: &[Allele]) -> Vec<String> {
    let mut out: Vec<String> = alleles.iter().map(Allele::display_string).collect();
    out.sort();
    out
}

/// A contig's length in the dictionary.
fn contig_length(dictionary: &[(String, i32)], contig: &str) -> i32 {
    dictionary
        .iter()
        .find(|(name, _)| name == contig)
        .map(|(_, length)| *length)
        .unwrap_or(0)
}

/// How a CNV linkage pads a record: by a fraction of its length, or by whole model bins.
#[derive(Debug, Clone)]
pub enum Padding {
    Fraction,
    /// The model's call intervals as (contig index, start, end), in their list order.
    Bins(Vec<(usize, i32, i32)>),
}

/// `CNVLinkage` and `BinnedCNVLinkage`.
#[derive(Debug, Clone)]
pub struct CnvLinkage {
    pub dictionary: Vec<(String, i32)>,
    pub padding_fraction: f64,
    pub min_sample_overlap: f64,
    pub padding: Padding,
}

impl CnvLinkage {
    fn contig_index(&self, contig: &str) -> usize {
        self.dictionary
            .iter()
            .position(|(name, _)| name == contig)
            .unwrap_or(usize::MAX)
    }

    /// `getPaddedRecordInterval`, as (start, end), or `None` for an interval that trims away.
    fn padded(
        &self,
        contig: &str,
        start: i32,
        end: i32,
    ) -> Result<Option<(i32, i32)>, EngineError> {
        let length = contig_length(&self.dictionary, contig);
        match &self.padding {
            Padding::Fraction => {
                let padding = (self.padding_fraction * f64::from(end - start + 1)) as i32;
                Ok(Some((
                    (start - padding).max(1),
                    (end + padding).min(length),
                )))
            }
            Padding::Bins(bins) => {
                let index = self.contig_index(contig);
                let key = |bin: &(usize, i32, i32)| (bin.0, bin.1, bin.2);
                // `TreeMap<GenomeLoc, Integer>`: the map is sorted by locus, and each entry's value
                // is its position in the list.
                let mut sorted: Vec<(usize, &(usize, i32, i32))> =
                    bins.iter().enumerate().collect();
                sorted.sort_by_key(|(_, bin)| key(bin));
                let start_bin = sorted
                    .iter()
                    .find(|(_, bin)| key(bin) >= (index, start, start))
                    .map(|(at, _)| *at);
                let Some(start_index) = start_bin else {
                    return Err(EngineError::illegal_state(format!(
                        "Call start {contig}:{start} for  call at {contig}:{start}-{end} not found in model call intervals."
                    )));
                };
                let end_bin = sorted
                    .iter()
                    .rev()
                    .find(|(_, bin)| key(bin) <= (index, end, end))
                    .map(|(at, _)| *at);
                let Some(end_index) = end_bin else {
                    return Err(EngineError::illegal_state(format!(
                        "Call end {contig}:{end} for call at {contig}:{start}-{end} not found in model call intervals."
                    )));
                };
                let bins_in_call = end_index as i64 - start_index as i64 + 1;
                if bins_in_call <= 0 {
                    return Err(EngineError::illegal_state(format!(
                        "Copy number call at {contig}:{start}-{end} does not align with supplied model calling intervals. Use the filtered intervals input from GermlineCNVCaller for this cohort/model."
                    )));
                }
                let pad = java_round(bins_in_call as f64 * self.padding_fraction);
                let padded_start_index = (start_index as i64 - pad).max(0) as usize;
                let padded_start = if bins[padded_start_index].0 == index {
                    bins[padded_start_index].1
                } else {
                    start
                };
                let padded_end_index =
                    ((end_index as i64 + pad).min(bins.len() as i64 - 1)) as usize;
                let padded_end = if bins[padded_end_index].0 == index {
                    bins[padded_end_index].2
                } else {
                    end
                };
                // `IntervalUtils.trimIntervalToContig`.
                let trimmed_start = padded_start.max(1);
                let trimmed_end = padded_end.min(length);
                if trimmed_start > trimmed_end {
                    return Ok(None);
                }
                Ok(Some((trimmed_start, trimmed_end)))
            }
        }
    }
}

/// `Math.round(double)` as its Javadoc specifies it: the closest long, ties toward positive infinity.
fn java_round(value: f64) -> i64 {
    (value + 0.5).floor() as i64
}

impl EngineLinkage for CnvLinkage {
    fn are_clusterable(&self, a: &Member, b: &Member) -> Result<bool, EngineError> {
        if !is_depth_only(a) || !is_depth_only(b) {
            return Ok(false);
        }
        if !is_simple_cnv(a) || !is_simple_cnv(b) {
            return Ok(false);
        }
        if a.call.contig_a != a.call.contig_b {
            return Err(EngineError::illegal_state(
                "Variant A is a CNV but interchromosomal".to_string(),
            ));
        }
        if b.call.contig_a != b.call.contig_b {
            return Err(EngineError::illegal_state(
                "Variant B is a CNV but interchromosomal".to_string(),
            ));
        }
        if a.call.sv_type != b.call.sv_type {
            return Ok(false);
        }
        let interval_a = self.padded(&a.call.contig_a, a.call.position_a, a.call.position_b)?;
        let interval_b = self.padded(&b.call.contig_a, b.call.position_a, b.call.position_b)?;
        let (Some(interval_a), Some(interval_b)) = (interval_a, interval_b) else {
            return Err(EngineError {
                class: "java.lang.IllegalArgumentException",
                message: "Invalid interval".to_string(),
            });
        };
        let overlaps = a.call.contig_a == b.call.contig_a
            && interval_a.0 <= interval_b.1
            && interval_b.0 <= interval_a.1;
        if !overlaps {
            return Ok(false);
        }
        if !has_sample_overlap(a, b, self.min_sample_overlap)? {
            return Ok(false);
        }
        let carriers_a = carriers(a)?;
        let carriers_b = carriers(b)?;
        let mut sorted_a = carriers_a.clone();
        let mut sorted_b = carriers_b.clone();
        sorted_a.sort();
        sorted_b.sort();
        if carriers_a.len() == 1 && sorted_a == sorted_b {
            let sample = &carriers_a[0];
            let genotype_a = a.genotypes.iter().find(|g| g.sample_name == *sample);
            let genotype_b = b.genotypes.iter().find(|g| g.sample_name == *sample);
            if let (Some(genotype_a), Some(genotype_b)) = (genotype_a, genotype_b) {
                if has_extended(genotype_a, "CN") && has_extended(genotype_b, "CN") {
                    let delta_a = genotype_a.ploidy() as i32 - genotype_int(genotype_a, "CN", 0);
                    let delta_b = genotype_b.ploidy() as i32 - genotype_int(genotype_b, "CN", 0);
                    if delta_a != delta_b {
                        return Ok(false);
                    }
                } else if sorted_alleles(&genotype_a.alleles) != sorted_alleles(&genotype_b.alleles)
                {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }

    fn max_clusterable_start(&self, item: &Member) -> Result<i32, EngineError> {
        let length = contig_length(&self.dictionary, &item.call.contig_a);
        if let Padding::Bins(_) = self.padding {
            return Ok(length);
        }
        if !is_simple_cnv(item) {
            return Ok(0);
        }
        let call_length = item.call.length.unwrap_or(0);
        let theoretical = ((f64::from(item.call.position_b)
            + self.padding_fraction * f64::from(call_length + length))
            / (1.0 + self.padding_fraction))
            .floor() as i32;
        Ok(theoretical.min(length))
    }
}

/// `CanonicalSVLinkage` over collapser members, with every sample overlap at zero so that no
/// carrier is ever asked for.
#[derive(Debug, Clone)]
pub struct CanonicalLinkage {
    pub dictionary: Vec<(String, i32)>,
    pub linkage: Linkage,
}

impl CanonicalLinkage {
    fn record(member: &Member) -> CallRecord {
        CallRecord {
            id: member.call.id.clone(),
            sv_type: member.call.sv_type,
            contig_a: member.call.contig_a.clone(),
            position_a: member.call.position_a,
            contig_b: member.call.contig_b.clone(),
            position_b: member.call.position_b,
            strand_a: member.call.strand_a,
            strand_b: member.call.strand_b,
            length: member.call.length,
            algorithms: member.call.algorithms.clone(),
            carriers: Vec::new(),
        }
    }

    fn max_with(
        &self,
        record: &CallRecord,
        parameters: &crate::sv_cluster::ClusteringParameters,
    ) -> i32 {
        let length = contig_length(&self.dictionary, &record.contig_a);
        let by_window = record
            .position_a
            .wrapping_add(parameters.window)
            .min(length);
        if !record.is_intrachromosomal() {
            return by_window;
        }
        let assumed = record.length_for(crate::sv_cluster::INSERTION_ASSUMED_LENGTH_FOR_OVERLAP);
        let by_overlap = ((f64::from(record.position_a)
            + (1.0 - parameters.reciprocal_overlap) * f64::from(assumed))
            as i32)
            .min(length);
        if parameters.requires_overlap_and_proximity {
            by_overlap.min(by_window)
        } else {
            by_overlap.max(by_window)
        }
    }
}

impl EngineLinkage for CanonicalLinkage {
    fn are_clusterable(&self, a: &Member, b: &Member) -> Result<bool, EngineError> {
        Ok(self
            .linkage
            .are_clusterable(&Self::record(a), &Self::record(b)))
    }

    fn max_clusterable_start(&self, item: &Member) -> Result<i32, EngineError> {
        let record = Self::record(item);
        let own = if record.is_depth_only() {
            self.linkage.depth
        } else {
            self.linkage.pesr
        };
        Ok(self
            .max_with(&record, &own)
            .max(self.max_with(&record, &self.linkage.mixed)))
    }
}

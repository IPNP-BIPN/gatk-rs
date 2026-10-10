//! Ported from `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.
//! ReadThreadingAssembler` (GATK 4.6.2.0), in its default sequence-graph mode: the local assembly
//! of one region into candidate haplotypes.
//!
//! Every graph stage is its own module (`read_threading_graph`, its dangling-end recovery,
//! `seq_graph`, `kbest_haplotype`); this is the order they run in and what is kept of them:
//!
//!  * the k-mer sizes are tried in increasing order, and only if none of them gives a graph are up
//!    to six more tried, each ten larger than the last, the sixth allowed a low-complexity graph and
//!    non-unique reference k-mers;
//!  * `createGraph` refuses a k-mer size whose reference has a repeated k-mer, whose graph has a
//!    cycle once pruned, or whose graph is of low complexity, and answers `FAILED` for a reference
//!    shorter than k (which, being a result, stops the retries);
//!  * the k best paths of every graph become haplotypes, aligned to the reference haplotype by
//!    `CigarUtils.calculateCigar`; one whose CIGAR has an N or spans fewer than 30 reference bases
//!    is dropped, and one whose alignment does not span the reference is dropped if the INDEL
//!    strategy would span it, else refused;
//!  * the haplotypes are a `LinkedHashSet`: after each graph the reference bases found as an
//!    alternative are removed and the reference haplotype added at the end, where it stays for
//!    the graphs that follow, and the result set is replaced by that set.
//!
//! The linked de Bruijn graph mode, the read error corrector, the flow-based haplotype collapsing
//! and the debugging graph output are not ported.

use htsjdk_bam::cigar::{Cigar, Op};
use htsjdk_bam::header::SamHeader;
use htsjdk_vcf::allele::AlleleError;

use crate::alignment_utils::AlignmentError;
use crate::assembly_region::AssemblyRegion;
use crate::assembly_result_set::{
    AssemblyResult, AssemblyResultSet, AssemblyStatus, ResultSetError,
};
use crate::cigar_utils::calculate_cigar;
use crate::clipping::{hard_clip_soft_clipped_bases, ClipError};
use crate::haplotype::Haplotype;
use crate::haplotype_alignment::HaplotypeCigarError;
use crate::interval::SimpleInterval;
use crate::kbest_haplotype::find_best_haplotypes;
use crate::read_pileup::sample_name;
use crate::read_threading_graph::{
    prune_adaptive, prune_low_weight_chains, AdaptivePruning, GraphError, ReadThreadingGraph,
};
use crate::seq_graph::{SeqGraph, SeqGraphError};
use crate::smith_waterman::{SmithWatermanAligner, SwOverhangStrategy, SwParameters};

/// `DEFAULT_NUM_PATHS_PER_GRAPH`.
pub const DEFAULT_NUM_PATHS_PER_GRAPH: usize = 128;
/// `KMER_SIZE_ITERATION_INCREASE`.
const KMER_SIZE_ITERATION_INCREASE: usize = 10;
/// `MAX_KMER_ITERATIONS_TO_ATTEMPT`.
const MAX_KMER_ITERATIONS_TO_ATTEMPT: usize = 6;
/// `DEFAULT_MIN_BASE_QUALITY_TO_USE`.
const DEFAULT_MIN_BASE_QUALITY_TO_USE: u8 = 10;
/// `MIN_HAPLOTYPE_REFERENCE_LENGTH`.
const MIN_HAPLOTYPE_REFERENCE_LENGTH: u32 = 30;
/// The sample a read without one is filed under: the reference files it under `null`, a key no
/// read group can produce.
const NULL_SAMPLE: &str = "\u{0}null";

/// `ChainPruner`: `LowWeightChainPruner(pruneFactor)` or `AdaptiveChainPruner`.
#[derive(Debug, Clone, Copy)]
pub enum ChainPruner {
    LowWeight(i32),
    Adaptive(AdaptivePruning),
}

/// What the assembler refuses.
#[derive(Debug, Clone, PartialEq)]
pub enum AssemblerError {
    /// `Utils.validateArg` and `ParamUtils` on the constructor's and the entry point's arguments.
    IllegalArgument(String),
    /// `IllegalStateException`: a graph without its reference ends, a reference path that is not
    /// the reference haplotype, or a Smith-Waterman alignment that cannot be explained.
    IllegalState(String),
    /// `NullPointerException`: the INDEL realignment the reference dereferences unchecked.
    NullPointer(&'static str),
    Graph(GraphError),
    SeqGraph(SeqGraphError),
    Alignment(AlignmentError),
    Clip(ClipError),
    Cigar(HaplotypeCigarError),
    Allele(AlleleError),
    ResultSet(ResultSetError),
}

impl From<GraphError> for AssemblerError {
    fn from(e: GraphError) -> Self {
        AssemblerError::Graph(e)
    }
}

impl From<SeqGraphError> for AssemblerError {
    fn from(e: SeqGraphError) -> Self {
        AssemblerError::SeqGraph(e)
    }
}

impl From<AlignmentError> for AssemblerError {
    fn from(e: AlignmentError) -> Self {
        AssemblerError::Alignment(e)
    }
}

impl From<ResultSetError> for AssemblerError {
    fn from(e: ResultSetError) -> Self {
        AssemblerError::ResultSet(e)
    }
}

/// `ReadThreadingAssembler`.
#[derive(Debug, Clone)]
pub struct ReadThreadingAssembler {
    kmer_sizes: Vec<usize>,
    dont_increase_kmer_sizes_for_cycles: bool,
    allow_non_unique_kmers_in_ref: bool,
    num_pruning_samples: usize,
    num_best_haplotypes_per_graph: usize,
    prune_before_cycle_counting: bool,
    remove_paths_not_connected_to_ref: bool,
    recover_dangling_branches: bool,
    recover_all_dangling_branches: bool,
    min_dangling_branch_length: i32,
    min_base_quality_to_use_in_assembly: u8,
    prune_factor: i32,
    chain_pruner: ChainPruner,
    min_matching_bases_to_dangling_end_recovery: i32,
}

/// The constructor's arguments, in its order, without `useLinkedDebruijnGraphs`.
#[derive(Debug, Clone)]
pub struct AssemblerSettings {
    pub max_allowed_paths: usize,
    pub kmer_sizes: Vec<usize>,
    pub dont_increase_kmer_sizes_for_cycles: bool,
    pub allow_non_unique_kmers_in_ref: bool,
    pub num_pruning_samples: usize,
    pub prune_factor: i32,
    pub use_adaptive_pruning: bool,
    pub initial_error_rate_for_pruning: f64,
    pub pruning_log_odds_threshold: f64,
    pub pruning_seeding_log_odds_threshold: f64,
    pub max_unpruned_variants: i32,
    pub enable_legacy_graph_cycle_detection: bool,
    pub min_matching_bases_to_dangling_end_recovery: i32,
}

impl ReadThreadingAssembler {
    /// The fourteen-argument constructor, in the sequence-graph mode.
    pub fn new(settings: &AssemblerSettings) -> Result<ReadThreadingAssembler, AssemblerError> {
        if settings.max_allowed_paths < 1 {
            return Err(AssemblerError::IllegalArgument(format!(
                "numBestHaplotypesPerGraph should be >= 1 but got {}",
                settings.max_allowed_paths
            )));
        }
        let mut kmer_sizes = settings.kmer_sizes.clone();
        kmer_sizes.sort_unstable();
        let chain_pruner = if settings.use_adaptive_pruning {
            ChainPruner::Adaptive(AdaptivePruning {
                initial_error_probability: settings.initial_error_rate_for_pruning,
                log_odds_threshold: settings.pruning_log_odds_threshold,
                seeding_log_odds_threshold: settings.pruning_seeding_log_odds_threshold,
                max_unpruned_variants: settings.max_unpruned_variants,
            })
        } else {
            ChainPruner::LowWeight(settings.prune_factor)
        };
        Ok(ReadThreadingAssembler {
            kmer_sizes,
            dont_increase_kmer_sizes_for_cycles: settings.dont_increase_kmer_sizes_for_cycles,
            allow_non_unique_kmers_in_ref: settings.allow_non_unique_kmers_in_ref,
            num_pruning_samples: settings.num_pruning_samples,
            num_best_haplotypes_per_graph: settings.max_allowed_paths,
            prune_before_cycle_counting: !settings.enable_legacy_graph_cycle_detection,
            remove_paths_not_connected_to_ref: true,
            recover_dangling_branches: true,
            recover_all_dangling_branches: false,
            min_dangling_branch_length: 0,
            min_base_quality_to_use_in_assembly: DEFAULT_MIN_BASE_QUALITY_TO_USE,
            prune_factor: settings.prune_factor,
            chain_pruner,
            min_matching_bases_to_dangling_end_recovery: settings
                .min_matching_bases_to_dangling_end_recovery,
        })
    }

    /// `setRecoverDanglingBranches`.
    pub fn set_recover_dangling_branches(&mut self, value: bool) {
        self.recover_dangling_branches = value;
    }

    /// `setRecoverAllDanglingBranches`, which also turns recovery on.
    pub fn set_recover_all_dangling_branches(&mut self, value: bool) {
        self.recover_all_dangling_branches = value;
        self.recover_dangling_branches = true;
    }

    /// `setMinDanglingBranchLength`.
    pub fn set_min_dangling_branch_length(&mut self, value: i32) {
        self.min_dangling_branch_length = value;
    }

    /// `setMinBaseQualityToUseInAssembly`.
    pub fn set_min_base_quality_to_use_in_assembly(&mut self, value: u8) {
        self.min_base_quality_to_use_in_assembly = value;
    }

    /// `setRemovePathsNotConnectedToRef`.
    pub fn set_remove_paths_not_connected_to_ref(&mut self, value: bool) {
        self.remove_paths_not_connected_to_ref = value;
    }

    /// `runLocalAssembly`, without a read error corrector or a haplotype collapsing engine.
    #[allow(clippy::too_many_arguments)]
    pub fn run_local_assembly(
        &self,
        region: &AssemblyRegion,
        mut ref_haplotype: Haplotype,
        full_reference_with_padding: &[u8],
        ref_loc: &SimpleInterval,
        header: &SamHeader,
        aligner: &dyn SmithWatermanAligner,
        dangling_end_sw_parameters: &SwParameters,
        haplotype_to_reference_sw_parameters: &SwParameters,
    ) -> Result<AssemblyResultSet, AssemblerError> {
        let ref_loc_size = (ref_loc.end - ref_loc.start + 1) as usize;
        if full_reference_with_padding.len() != ref_loc_size {
            return Err(AssemblerError::IllegalArgument(
                "Reference bases and reference loc must be the same size.".to_string(),
            ));
        }
        if self.prune_factor < 0 {
            return Err(AssemblerError::IllegalArgument(format!(
                "Pruning factor cannot be negative ({})",
                self.prune_factor
            )));
        }
        let reads = self.hard_clipped_reads(region, header)?;

        let mut result_set = AssemblyResultSet::new();
        result_set.set_region_for_genotyping(region.clone());
        result_set.set_full_reference_with_padding(full_reference_with_padding.to_vec());
        result_set.set_padded_reference_loc(ref_loc.clone());
        let active_region_extended_location = region.padded_span().clone();
        ref_haplotype.set_genome_location(active_region_extended_location.clone());
        result_set.add(ref_haplotype.clone())?;

        // `assembleKmerGraphsAndHaplotypeCall`.
        let mut graphs = Vec::new();
        for result in self.assemble(&reads, &ref_haplotype, aligner, dangling_end_sw_parameters)? {
            if result.status == AssemblyStatus::AssembledSomeVariation {
                let graph = result.seq_graph.expect("a graph that assembled variation");
                sanity_check_reference_graph(&graph, &ref_haplotype)?;
                graphs.push((graph, result.kmer_size.expect("a graph has a k-mer size")));
            }
        }
        self.find_best_paths(
            &graphs,
            &ref_haplotype,
            &active_region_extended_location,
            &mut result_set,
            aligner,
            haplotype_to_reference_sw_parameters,
        )?;
        Ok(result_set)
    }

    /// `ReadClipper.hardClipSoftClippedBases` over the region's reads, each with its sample: what
    /// `runLocalAssembly` hands to `assemble`.
    pub fn hard_clipped_reads(
        &self,
        region: &AssemblyRegion,
        header: &SamHeader,
    ) -> Result<Vec<SampleRead>, AssemblerError> {
        region
            .reads()
            .iter()
            .map(|read| {
                let clipped = hard_clip_soft_clipped_bases(read, Some(header), 0)
                    .map_err(AssemblerError::Clip)?;
                Ok(SampleRead {
                    sample: sample_name(read, header).unwrap_or_else(|| NULL_SAMPLE.to_string()),
                    bases: clipped.read_bases,
                    qualities: clipped.base_qualities,
                })
            })
            .collect()
    }

    /// `findBestPaths`, writing into the result set.
    fn find_best_paths(
        &self,
        graphs: &[(SeqGraph, usize)],
        ref_haplotype: &Haplotype,
        active_region_window: &SimpleInterval,
        result_set: &mut AssemblyResultSet,
        aligner: &dyn SmithWatermanAligner,
        haplotype_to_reference_sw_parameters: &SwParameters,
    ) -> Result<(), AssemblerError> {
        let mut return_haplotypes: Vec<Haplotype> = Vec::new();
        let active_region_start = ref_haplotype.alignment_start_hap_wrt_ref();
        let ref_bases = ref_haplotype.bases();
        let ref_reference_length = ref_haplotype
            .cigar()
            .ok_or(AssemblerError::NullPointer("refHaplotype.getCigar()"))?
            .reference_length();
        for (graph, kmer_size) in graphs {
            let (Some(source), Some(sink)) = (
                graph.reference_source_vertex(),
                graph.reference_sink_vertex(),
            ) else {
                return Err(AssemblerError::IllegalArgument(
                    "Both source and sink cannot be null".to_string(),
                ));
            };
            let best = find_best_haplotypes(
                graph,
                &[source],
                &[sink],
                self.num_best_haplotypes_per_graph,
            )?;
            for path in best {
                // `KBestHaplotype.haplotype()`: the flag the sequence-graph finder never sets.
                let mut h = Haplotype::new(path.bases(), path.is_reference)
                    .map_err(AssemblerError::Allele)?;
                h.set_score(path.score);
                h.set_kmer_size(*kmer_size as i32);
                if return_haplotypes.contains(&h) {
                    continue;
                }
                let h_bases = h.bases();
                let Some(cigar) = calculate_cigar(
                    &ref_bases,
                    &h_bases,
                    aligner,
                    haplotype_to_reference_sw_parameters,
                    SwOverhangStrategy::SoftClip,
                )?
                else {
                    // A failed cigar: the haplotype is ignored.
                    continue;
                };
                if cigar.is_empty() {
                    return Err(AssemblerError::IllegalState(format!(
                        "Smith-Waterman alignment failure. Cigar = {} with reference length {} but \
                         expecting reference length of {}",
                        cigar.to_text(),
                        cigar.reference_length(),
                        ref_reference_length
                    )));
                } else if path_is_too_divergent_from_reference(&cigar)
                    || cigar.reference_length() < MIN_HAPLOTYPE_REFERENCE_LENGTH
                {
                    continue;
                } else if cigar.reference_length() != ref_reference_length {
                    let with_indel_strategy = calculate_cigar(
                        &ref_bases,
                        &h_bases,
                        aligner,
                        haplotype_to_reference_sw_parameters,
                        SwOverhangStrategy::Indel,
                    )?
                    .ok_or(AssemblerError::NullPointer("cigarWithIndelStrategy"))?;
                    if with_indel_strategy.reference_length() == ref_reference_length {
                        continue;
                    }
                    return Err(AssemblerError::IllegalState(format!(
                        "Smith-Waterman alignment failure. Cigar = {} with reference length {} but \
                         expecting reference length of {}",
                        cigar.to_text(),
                        cigar.reference_length(),
                        ref_reference_length
                    )));
                }
                h.set_cigar(&cigar).map_err(AssemblerError::Cigar)?;
                h.set_alignment_start_hap_wrt_ref(active_region_start);
                h.set_genome_location(active_region_window.clone());
                return_haplotypes.push(h.clone());
                result_set.add(h)?;
            }

            // The reference bases found as an alternative give way to the reference haplotype,
            // which is added at the end unless an earlier graph already put it there.
            if !return_haplotypes.is_empty() {
                let tmp_ref = Haplotype::new(&ref_bases, false).map_err(AssemblerError::Allele)?;
                return_haplotypes.retain(|h| *h != tmp_ref);
                if !return_haplotypes.contains(ref_haplotype) {
                    return_haplotypes.push(ref_haplotype.clone());
                }
                result_set.replace_all_haplotypes(&return_haplotypes)?;
            }
        }
        Ok(())
    }

    /// `assemble`: the configured k-mer sizes, then the larger ones if none gave a result.
    pub fn assemble(
        &self,
        reads: &[SampleRead],
        ref_haplotype: &Haplotype,
        aligner: &dyn SmithWatermanAligner,
        dangling_end_sw_parameters: &SwParameters,
    ) -> Result<Vec<AssemblyResult>, AssemblerError> {
        let mut results = Vec::new();
        for &kmer_size in &self.kmer_sizes {
            results.extend(self.create_graph(
                reads,
                ref_haplotype,
                kmer_size,
                self.dont_increase_kmer_sizes_for_cycles,
                self.allow_non_unique_kmers_in_ref,
                aligner,
                dangling_end_sw_parameters,
            )?);
        }
        if results.is_empty() && !self.dont_increase_kmer_sizes_for_cycles {
            let mut kmer_size = self.max_kmer_size()? + KMER_SIZE_ITERATION_INCREASE;
            let mut iterations = 1;
            while results.is_empty() && iterations <= MAX_KMER_ITERATIONS_TO_ATTEMPT {
                let last_attempt = iterations == MAX_KMER_ITERATIONS_TO_ATTEMPT;
                results.extend(self.create_graph(
                    reads,
                    ref_haplotype,
                    kmer_size,
                    last_attempt,
                    last_attempt,
                    aligner,
                    dangling_end_sw_parameters,
                )?);
                kmer_size += KMER_SIZE_ITERATION_INCREASE;
                iterations += 1;
            }
        }
        Ok(results)
    }

    /// `arrayMaxInt(kmerSizes)`.
    fn max_kmer_size(&self) -> Result<usize, AssemblerError> {
        self.kmer_sizes
            .iter()
            .copied()
            .max()
            .ok_or_else(|| AssemblerError::IllegalArgument("Array size cannot be 0!".to_string()))
    }

    /// `createGraph`: `None` where the reference returns null.
    #[allow(clippy::too_many_arguments)]
    fn create_graph(
        &self,
        reads: &[SampleRead],
        ref_haplotype: &Haplotype,
        kmer_size: usize,
        allow_low_complexity_graphs: bool,
        allow_non_unique_kmers_in_ref: bool,
        aligner: &dyn SmithWatermanAligner,
        dangling_end_sw_parameters: &SwParameters,
    ) -> Result<Option<AssemblyResult>, AssemblerError> {
        let ref_bases = ref_haplotype.bases();
        if ref_bases.len() < kmer_size {
            return Ok(Some(AssemblyResult {
                status: AssemblyStatus::Failed,
                seq_graph: None,
                kmer_size: None,
                discovered_haplotypes: Vec::new(),
            }));
        }
        if !allow_non_unique_kmers_in_ref && has_non_unique_kmers(&ref_bases, kmer_size) {
            return Ok(None);
        }

        let mut graph = ReadThreadingGraph::new(
            kmer_size,
            self.min_base_quality_to_use_in_assembly,
            self.num_pruning_samples,
        );
        graph.set_min_matching_bases_to_dangling_end_recovery(
            self.min_matching_bases_to_dangling_end_recovery,
        );
        graph.set_threading_start_only_at_existing_vertex(!self.recover_dangling_branches);
        graph.add_sequence(&ref_bases, true)?;
        for read in reads {
            graph.add_read(&read.sample, &read.bases, &read.qualities)?;
        }
        graph.build_graph_if_necessary()?;

        if self.prune_before_cycle_counting {
            self.prune(&mut graph)?;
        }
        if graph.has_cycles() {
            return Ok(None);
        }
        if !allow_low_complexity_graphs && graph.is_low_quality_graph() {
            return Ok(None);
        }
        let result = self.get_assembly_result(&mut graph, aligner, dangling_end_sw_parameters)?;
        if self.recover_all_dangling_branches && graph.has_cycles() {
            return Ok(None);
        }
        Ok(Some(result))
    }

    /// `chainPruner.pruneLowWeightChains(graph)`.
    fn prune(&self, graph: &mut ReadThreadingGraph) -> Result<(), AssemblerError> {
        match &self.chain_pruner {
            ChainPruner::LowWeight(factor) => prune_low_weight_chains(graph, *factor),
            ChainPruner::Adaptive(params) => prune_adaptive(graph, params)?,
        }
        Ok(())
    }

    /// `getAssemblyResult`.
    fn get_assembly_result(
        &self,
        graph: &mut ReadThreadingGraph,
        aligner: &dyn SmithWatermanAligner,
        dangling_end_sw_parameters: &SwParameters,
    ) -> Result<AssemblyResult, AssemblerError> {
        if !self.prune_before_cycle_counting {
            self.prune(graph)?;
        }
        if self.recover_dangling_branches {
            graph.recover_dangling_tails(
                self.prune_factor,
                self.min_dangling_branch_length,
                self.recover_all_dangling_branches,
                aligner,
                dangling_end_sw_parameters,
            )?;
            graph.recover_dangling_heads(
                self.prune_factor,
                self.min_dangling_branch_length,
                self.recover_all_dangling_branches,
                aligner,
                dangling_end_sw_parameters,
            )?;
        }
        if self.remove_paths_not_connected_to_ref {
            graph.remove_paths_not_connected_to_ref()?;
        }
        let mut seq_graph = SeqGraph::from_read_threading_graph(graph);
        seq_graph.clean_non_ref_paths();
        let status = cleanup_seq_graph(&mut seq_graph)?;
        Ok(AssemblyResult {
            status,
            seq_graph: Some(seq_graph),
            kmer_size: Some(graph.kmer_size()),
            discovered_haplotypes: Vec::new(),
        })
    }
}

/// A read as the graph takes it: its sample, and its bases and qualities once soft clips are hard
/// clipped.
#[derive(Debug, Clone)]
pub struct SampleRead {
    pub sample: String,
    pub bases: Vec<u8>,
    pub qualities: Vec<u8>,
}

/// `cleanupSeqGraph`.
fn cleanup_seq_graph(graph: &mut SeqGraph) -> Result<AssemblyStatus, AssemblerError> {
    graph.zip_linear_chains();
    graph.remove_singleton_orphan_vertices();
    graph.remove_vertices_not_connected_to_ref_regardless_of_edge_direction();
    graph.simplify_graph()?;
    if graph.reference_source_vertex().is_none() || graph.reference_sink_vertex().is_none() {
        return Ok(AssemblyStatus::JustAssembledReference);
    }
    graph.remove_paths_not_connected_to_ref()?;
    graph.simplify_graph()?;
    if graph.vertex_count() == 1 {
        // A graph of one vertex is given an empty one after it, so that it has a path.
        let complete = graph.vertex_ids().next().expect("one vertex");
        let dummy = graph.add_vertex(b"");
        graph.add_edge(complete, dummy, true, 0);
    }
    Ok(AssemblyStatus::AssembledSomeVariation)
}

/// `ReadThreadingGraph.determineNonUniqueKmers` over the reference is not empty.
fn has_non_unique_kmers(sequence: &[u8], kmer_size: usize) -> bool {
    let mut seen = std::collections::HashSet::new();
    sequence.windows(kmer_size).any(|kmer| !seen.insert(kmer))
}

/// `pathIsTooDivergentFromReference`: the CIGAR has an N.
fn path_is_too_divergent_from_reference(cigar: &Cigar) -> bool {
    cigar.elements.iter().any(|e| e.op == Op::N)
}

/// `sanityCheckReferenceGraph`: the graph has both reference ends, and the reference path between
/// them spells the reference haplotype.
fn sanity_check_reference_graph(
    graph: &SeqGraph,
    ref_haplotype: &Haplotype,
) -> Result<(), AssemblerError> {
    let Some(source) = graph.reference_source_vertex() else {
        return Err(AssemblerError::IllegalState(
            "All reference graphs must have a reference source vertex.".to_string(),
        ));
    };
    let Some(sink) = graph.reference_sink_vertex() else {
        return Err(AssemblerError::IllegalState(
            "All reference graphs must have a reference sink vertex.".to_string(),
        ));
    };
    if reference_bytes(graph, source, sink) != ref_haplotype.bases() {
        return Err(AssemblerError::IllegalState(
            "Mismatch between the reference haplotype and the reference assembly graph path."
                .to_string(),
        ));
    }
    Ok(())
}

/// `getReferenceBytes(source, sink, true, true)`: each vertex's sequence along the first reference
/// edge out of the last.
fn reference_bytes(graph: &SeqGraph, source: usize, sink: usize) -> Vec<u8> {
    let next = |v: usize| {
        graph
            .outgoing_edges(v)
            .iter()
            .map(|&e| graph.edge(e))
            .find(|e| e.is_ref)
            .map(|e| e.target)
    };
    let mut bytes = graph.sequence(source).to_vec();
    let mut v = next(source);
    while let Some(vertex) = v {
        if vertex == sink {
            break;
        }
        bytes.extend_from_slice(graph.sequence(vertex));
        v = next(vertex);
    }
    if v == Some(sink) {
        bytes.extend_from_slice(graph.sequence(sink));
    }
    bytes
}

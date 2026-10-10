//! Ported from `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading`
//! (`AbstractReadThreadingGraph`: `recoverDanglingTails`, `recoverDanglingHeads` and their
//! helpers), `haplotypecaller.graphs.BaseGraph` (`isReferenceNode`, `getNextReferenceVertex`,
//! `getPrevReferenceVertex`) and `utils.read.AlignmentUtils.removeTrailingDeletions`, GATK 4.6.2.0.
//!
//! A dangling tail is a vertex with no outgoing edge that is not the reference sink; a dangling
//! head, one with no incoming edge that is not the reference source. Each is walked back to where
//! it leaves the reference, the branch is aligned against the reference from that vertex with
//! Smith-Waterman (`LEADING_INDEL`, trailing deletions dropped), and a CIGAR of at most three
//! elements that starts (heads) or ends (tails) with an M lets one edge, multiplicity one, join the
//! branch back to the reference.
//!
//! What a reader would not guess:
//!
//! * tails are recovered while the vertex set is iterated, heads from a list taken first, and the
//!   head merge can add vertices (`extendDanglingPathAgainstReference`) that never enter the k-mer
//!   map, while the head it replaces stays in the graph with no edge at all;
//! * a recovered edge goes through JGraphT's `addEdge`, which refuses a second edge between the
//!   same two vertices, and the merge still counts as recovered;
//! * a head's bases are read backwards from the merge point: every vertex adds its last base,
//!   except a source, which adds its whole k-mer reversed;
//! * the legacy head merge (`numDanglingMatchingPrefixBases` of -1, the default) merges at the
//!   LAST mismatch it tolerates, not at the end of the match;
//! * the extension copies a new vertex's bases with `Arrays.copyOfRange`, so the walk never
//!   reads past the extended sequence, and its edges take the source edge's multiplicity but keep
//!   a pruning multiplicity of one.

use htsjdk_bam::cigar::{Cigar, CigarElement, Op};

use super::{Edge, GraphError, ReadThreadingGraph, Vertex};
use crate::smith_waterman::{SmithWatermanAligner, SwOverhangStrategy, SwParameters};

/// `MAX_CIGAR_COMPLEXITY`.
const MAX_CIGAR_COMPLEXITY: usize = 3;

/// `DanglingChainMergeHelper`.
struct DanglingChainMergeHelper {
    dangling_path: Vec<usize>,
    reference_path: Vec<usize>,
    dangling_path_string: Vec<u8>,
    reference_path_string: Vec<u8>,
    cigar: Cigar,
}

/// `List.get` and an array read: Java's bounds check, as an error.
fn at<T: Copy>(items: &[T], index: i32) -> Result<T, GraphError> {
    usize::try_from(index)
        .ok()
        .and_then(|i| items.get(i).copied())
        .ok_or(GraphError::IndexOutOfBounds)
}

/// `cigarIsOkayToMerge`.
fn cigar_is_okay_to_merge(cigar: &Cigar, require_first_m: bool, require_last_m: bool) -> bool {
    let elements = &cigar.elements;
    if elements.is_empty() || elements.len() > MAX_CIGAR_COMPLEXITY {
        return false;
    }
    if require_first_m && elements[0].op != Op::M {
        return false;
    }
    !(require_last_m && elements[elements.len() - 1].op != Op::M)
}

/// `longestSuffixMatch`: how many bases of `kmer`'s end match `seq` ending at `seq_start`.
fn longest_suffix_match(seq: &[u8], kmer: &[u8], seq_start: i32) -> Result<i32, GraphError> {
    for len in 1..=kmer.len() as i32 {
        let seq_i = seq_start - len + 1;
        let kmer_i = kmer.len() as i32 - len;
        if seq_i < 0 || at(seq, seq_i)? != at(kmer, kmer_i)? {
            return Ok(len - 1);
        }
    }
    Ok(kmer.len() as i32)
}

/// `AlignmentUtils.removeTrailingDeletions`.
fn remove_trailing_deletions(cigar: Cigar) -> Result<Cigar, GraphError> {
    let last = cigar.elements.last().ok_or(GraphError::IndexOutOfBounds)?;
    if last.op != Op::D {
        return Ok(cigar);
    }
    let mut elements = cigar.elements;
    elements.pop();
    Ok(Cigar::new(elements))
}

fn reference_bases(elements: &[CigarElement]) -> i32 {
    elements
        .iter()
        .filter(|e| e.op.consumes_reference_bases())
        .map(|e| e.length as i32)
        .sum()
}

impl ReadThreadingGraph {
    /// `recoverDanglingTails(pruneFactor, minDanglingBranchLength, recoverAll, aligner,
    /// danglingTailSWParameters)`.
    pub fn recover_dangling_tails(
        &mut self,
        prune_factor: i32,
        min_dangling_branch_length: i32,
        recover_all: bool,
        aligner: &dyn SmithWatermanAligner,
        parameters: &SwParameters,
    ) -> Result<(), GraphError> {
        self.validate_recovery(prune_factor, min_dangling_branch_length)?;
        // Only edges are added, so the vertex set iterated is the one the loop started with.
        for vertex in 0..self.vertices.len() {
            if self.vertex_alive[vertex]
                && self.out_degree(vertex) == 0
                && !self.is_ref_sink(vertex)
            {
                self.recover_dangling_tail(
                    vertex,
                    prune_factor,
                    min_dangling_branch_length,
                    recover_all,
                    aligner,
                    parameters,
                )?;
            }
        }
        Ok(())
    }

    /// `recoverDanglingHeads`, over the heads listed before any is merged.
    pub fn recover_dangling_heads(
        &mut self,
        prune_factor: i32,
        min_dangling_branch_length: i32,
        recover_all: bool,
        aligner: &dyn SmithWatermanAligner,
        parameters: &SwParameters,
    ) -> Result<(), GraphError> {
        self.validate_recovery(prune_factor, min_dangling_branch_length)?;
        let heads: Vec<usize> = self
            .vertex_ids()
            .filter(|&v| self.in_degree(v) == 0 && !self.is_ref_source(v))
            .collect();
        for vertex in heads {
            self.recover_dangling_head(
                vertex,
                prune_factor,
                min_dangling_branch_length,
                recover_all,
                aligner,
                parameters,
            )?;
        }
        Ok(())
    }

    fn validate_recovery(&self, prune_factor: i32, min_branch: i32) -> Result<(), GraphError> {
        if prune_factor < 0 {
            return Err(GraphError::IllegalArgument(
                "pruneFactor must be non-negative",
            ));
        }
        if min_branch < 0 {
            return Err(GraphError::IllegalArgument(
                "minDanglingBranchLength must be non-negative",
            ));
        }
        if !self.already_built {
            return Err(GraphError::IllegalState(
                "recovering dangling ends requires the graph be already built",
            ));
        }
        Ok(())
    }

    /// `recoverDanglingTail`: 1 when an edge was added (or refused as a duplicate), 0 otherwise.
    fn recover_dangling_tail(
        &mut self,
        vertex: usize,
        prune_factor: i32,
        min_branch: i32,
        recover_all: bool,
        aligner: &dyn SmithWatermanAligner,
        parameters: &SwParameters,
    ) -> Result<i32, GraphError> {
        if self.out_degree(vertex) != 0 {
            return Err(GraphError::IllegalState(
                "dangling tail with out-degree > 0",
            ));
        }
        let merge = self.generate_cigar_against_downwards_reference_path(
            vertex,
            prune_factor,
            min_branch,
            recover_all,
            aligner,
            parameters,
        )?;
        match merge {
            Some(merge) if cigar_is_okay_to_merge(&merge.cigar, false, true) => {
                self.merge_dangling_tail(&merge)
            }
            _ => Ok(0),
        }
    }

    /// `recoverDanglingHead`.
    fn recover_dangling_head(
        &mut self,
        vertex: usize,
        prune_factor: i32,
        min_branch: i32,
        recover_all: bool,
        aligner: &dyn SmithWatermanAligner,
        parameters: &SwParameters,
    ) -> Result<i32, GraphError> {
        if self.in_degree(vertex) != 0 {
            return Err(GraphError::IllegalState("dangling head with in-degree > 0"));
        }
        let merge = self.generate_cigar_against_upwards_reference_path(
            vertex,
            prune_factor,
            min_branch,
            recover_all,
            aligner,
            parameters,
        )?;
        match merge {
            Some(mut merge) if cigar_is_okay_to_merge(&merge.cigar, true, false) => {
                if self.min_matching_bases >= 0 {
                    self.merge_dangling_head(&mut merge)
                } else {
                    self.merge_dangling_head_legacy(&mut merge)
                }
            }
            _ => Ok(0),
        }
    }

    /// `mergeDanglingTail`.
    fn merge_dangling_tail(&mut self, merge: &DanglingChainMergeHelper) -> Result<i32, GraphError> {
        let elements = &merge.cigar.elements;
        let last = elements[elements.len() - 1];
        if last.op != Op::M {
            return Err(GraphError::IllegalArgument(
                "The last Cigar element must be an M",
            ));
        }
        let last_ref_index = merge.cigar.reference_length() as i32 - 1;
        let matching_suffix = longest_suffix_match(
            &merge.reference_path_string,
            &merge.dangling_path_string,
            last_ref_index,
        )?
        .min(last.length as i32);
        let too_short = if self.min_matching_bases >= 0 {
            matching_suffix < self.min_matching_bases
        } else {
            matching_suffix == 0
        };
        if too_short {
            return Ok(0);
        }

        let alt_index = (merge.cigar.read_length() as i32 - matching_suffix - 1).max(0);
        // A deletion left-aligned onto the common ancestor, with the rest a perfect suffix match:
        // merge one further down so the deletion keeps its base.
        let leading_deletion = elements[0].op == Op::D
            && elements[0].length as i32 + matching_suffix == last_ref_index + 1;
        let ref_index = last_ref_index - matching_suffix + 1 + i32::from(leading_deletion);
        // The whole tail in an insertion would merge back onto the ancestor: a cycle.
        if ref_index == 0 {
            return Ok(0);
        }
        let source = at(&merge.dangling_path, alt_index)?;
        let target = at(&merge.reference_path, ref_index)?;
        self.add_edge_unless_present(source, target);
        Ok(1)
    }

    /// `mergeDanglingHeadLegacy`: the old merge, which does not handle indels.
    fn merge_dangling_head_legacy(
        &mut self,
        merge: &mut DanglingChainMergeHelper,
    ) -> Result<i32, GraphError> {
        let first = merge.cigar.elements[0];
        if first.op != Op::M {
            return Err(GraphError::IllegalArgument(
                "The first Cigar element must be an M",
            ));
        }
        let index = self.best_prefix_match_legacy(
            &merge.reference_path_string,
            &merge.dangling_path_string,
            first.length as i32,
        )?;
        if index <= 0 || index >= merge.reference_path.len() as i32 - 1 {
            return Ok(0);
        }
        if index >= merge.dangling_path.len() as i32 {
            let extend = index - merge.dangling_path.len() as i32 + 2;
            let elements = merge.cigar.elements.clone();
            if !self.extend_dangling_path_against_reference(merge, extend, &elements)? {
                return Ok(0);
            }
        }
        let source = at(&merge.reference_path, index + 1)?;
        let target = at(&merge.dangling_path, index)?;
        self.add_edge_unless_present(source, target);
        Ok(1)
    }

    /// `mergeDanglingHead`.
    fn merge_dangling_head(
        &mut self,
        merge: &mut DanglingChainMergeHelper,
    ) -> Result<i32, GraphError> {
        let elements = merge.cigar.elements.clone();
        if elements[0].op != Op::M {
            return Err(GraphError::IllegalArgument(
                "The first Cigar element must be an M",
            ));
        }
        let (ref_index, read_index) = self.best_prefix_match(
            &elements,
            &merge.reference_path_string,
            &merge.dangling_path_string,
        )?;
        if ref_index <= 0 || read_index <= 0 {
            return Ok(0);
        }
        if ref_index >= merge.reference_path.len() as i32 - 1 {
            return Ok(0);
        }
        if read_index >= merge.dangling_path.len() as i32 {
            let extend = read_index - merge.dangling_path.len() as i32 + 2;
            if !self.extend_dangling_path_against_reference(merge, extend, &elements)? {
                return Ok(0);
            }
        }
        let source = at(&merge.reference_path, ref_index + 1)?;
        let target = at(&merge.dangling_path, read_index)?;
        self.add_edge_unless_present(source, target);
        Ok(1)
    }

    /// `bestPrefixMatch`: walking back from the alignment's end through its trailing M elements,
    /// the offsets of the first mismatch, or `(-1, -1)` with fewer matches than the minimum.
    fn best_prefix_match(
        &self,
        elements: &[CigarElement],
        path1: &[u8],
        path2: &[u8],
    ) -> Result<(i32, i32), GraphError> {
        let mut ref_index = reference_bases(elements) - 1;
        let mut read_index = path2.len() as i32 - 1;
        'cigar: for element in elements.iter().rev() {
            if !(element.op.consumes_read_bases() && element.op.consumes_reference_bases()) {
                break;
            }
            for _ in 0..element.length {
                if at(path1, ref_index)? != at(path2, read_index)? {
                    break 'cigar;
                }
                ref_index -= 1;
                read_index -= 1;
            }
        }
        let matches = path2.len() as i32 - 1 - read_index;
        Ok(if matches < self.min_matching_bases {
            (-1, -1)
        } else {
            (ref_index, read_index)
        })
    }

    /// `bestPrefixMatchLegacy`: the index of the last mismatch before `max_index`, or -1 past
    /// `getMaxMismatchesLegacy` of them (the branch length over the k-mer size, at least one; the
    /// test-only override is never set).
    fn best_prefix_match_legacy(
        &self,
        path1: &[u8],
        path2: &[u8],
        max_index: i32,
    ) -> Result<i32, GraphError> {
        let max_mismatches = (max_index / self.kmer_size as i32).max(1);
        let mut mismatches = 0;
        let mut last_good_index = -1;
        for index in 0..max_index {
            if at(path1, index)? != at(path2, index)? {
                mismatches += 1;
                if mismatches > max_mismatches {
                    return Ok(-1);
                }
                last_good_index = index;
            }
        }
        Ok(last_good_index)
    }

    /// `extendDanglingPathAgainstReference`: replace the head's source by `num_nodes` new
    /// vertices that prepend the reference's bases to it.
    fn extend_dangling_path_against_reference(
        &mut self,
        merge: &mut DanglingChainMergeHelper,
        num_nodes: i32,
        elements: &[CigarElement],
    ) -> Result<bool, GraphError> {
        let last = merge.dangling_path.len() as i32 - 1;
        let offset: i32 = elements
            .iter()
            .map(|e| {
                let reference = if e.op.consumes_reference_bases() {
                    e.length as i32
                } else {
                    0
                };
                let read = if e.op.consumes_read_bases() {
                    e.length as i32
                } else {
                    0
                };
                reference - read
            })
            .sum();
        let ref_node = last + offset + num_nodes;
        if ref_node >= merge.reference_path.len() as i32 {
            return Ok(false);
        }
        let ref_vertex = at(&merge.reference_path, ref_node)?;
        let source = at(&merge.dangling_path, last)?;
        merge.dangling_path.remove(last as usize);
        let ref_sequence = &self.vertices[ref_vertex].sequence;
        let mut sequence = Vec::with_capacity(num_nodes as usize + self.kmer_size);
        for i in 0..num_nodes {
            sequence.push(at(ref_sequence, i)?);
        }
        sequence.extend_from_slice(&self.vertices[source].sequence);

        let source_edge = self.heaviest_outgoing_edge(source)?;
        let multiplicity = self.edges[source_edge].multiplicity;
        let mut previous = self.edges[source_edge].target;
        self.remove_edge(source_edge);
        for i in (1..=num_nodes as usize).rev() {
            // `Arrays.copyOfRange`, which pads with zeros past the end.
            let mut kmer = vec![0u8; self.kmer_size];
            let available = sequence.len().saturating_sub(i).min(self.kmer_size);
            kmer[..available].copy_from_slice(&sequence[i..i + available]);
            let vertex = self.add_vertex_outside_kmer_map(kmer);
            let edge = self.edges.len();
            self.add_edge(vertex, previous, false, 1);
            self.edges[edge].multiplicity = multiplicity;
            merge.dangling_path.push(vertex);
            previous = vertex;
        }
        Ok(true)
    }

    /// `generateCigarAgainstDownwardsReferencePath`.
    fn generate_cigar_against_downwards_reference_path(
        &self,
        vertex: usize,
        prune_factor: i32,
        min_branch: i32,
        recover_all: bool,
        aligner: &dyn SmithWatermanAligner,
        parameters: &SwParameters,
    ) -> Result<Option<DanglingChainMergeHelper>, GraphError> {
        // Heads can be zero long, tails cannot.
        let min_tail = min_branch.max(1);
        let Some(alt_path) =
            self.find_path_upwards_to_lowest_common_ancestor(vertex, prune_factor, !recover_all)?
        else {
            return Ok(None);
        };
        if self.is_ref_source(alt_path[0]) || (alt_path.len() as i32) < min_tail + 1 {
            return Ok(None);
        }
        let blacklisted = self.heaviest_incoming_edge(alt_path[1])?;
        let ref_path = self.reference_path_downwards(alt_path[0], Some(blacklisted));
        let ref_bases = self.bases_for_path(&ref_path, false);
        let alt_bases = self.bases_for_path(&alt_path, false);
        let alignment = aligner
            .align(
                &ref_bases,
                &alt_bases,
                parameters,
                SwOverhangStrategy::LeadingIndel,
            )
            .map_err(GraphError::Alignment)?;
        Ok(Some(DanglingChainMergeHelper {
            dangling_path: alt_path,
            reference_path: ref_path,
            dangling_path_string: alt_bases,
            reference_path_string: ref_bases,
            cigar: remove_trailing_deletions(alignment.cigar)?,
        }))
    }

    /// `generateCigarAgainstUpwardsReferencePath`.
    fn generate_cigar_against_upwards_reference_path(
        &self,
        vertex: usize,
        prune_factor: i32,
        min_branch: i32,
        recover_all: bool,
        aligner: &dyn SmithWatermanAligner,
        parameters: &SwParameters,
    ) -> Result<Option<DanglingChainMergeHelper>, GraphError> {
        let Some(alt_path) = self.find_path_downwards_to_highest_common_descendant_of_reference(
            vertex,
            prune_factor,
            !recover_all,
        )?
        else {
            return Ok(None);
        };
        if self.is_ref_sink(alt_path[0]) || (alt_path.len() as i32) < min_branch + 1 {
            return Ok(None);
        }
        let ref_path = self.reference_path_upwards(alt_path[0]);
        let ref_bases = self.bases_for_path(&ref_path, true);
        let alt_bases = self.bases_for_path(&alt_path, true);
        let alignment = aligner
            .align(
                &ref_bases,
                &alt_bases,
                parameters,
                SwOverhangStrategy::LeadingIndel,
            )
            .map_err(GraphError::Alignment)?;
        Ok(Some(DanglingChainMergeHelper {
            dangling_path: alt_path,
            reference_path: ref_path,
            dangling_path_string: alt_bases,
            reference_path_string: ref_bases,
            cigar: remove_trailing_deletions(alignment.cigar)?,
        }))
    }

    /// `findPathUpwardsToLowestCommonAncestor`.
    fn find_path_upwards_to_lowest_common_ancestor(
        &self,
        vertex: usize,
        prune_factor: i32,
        give_up_at_branch: bool,
    ) -> Result<Option<Vec<usize>>, GraphError> {
        if give_up_at_branch {
            self.find_path(
                vertex,
                prune_factor,
                |g, v| g.in_degree(v) != 1 || g.out_degree(v) >= 2,
                |g, v| g.out_degree(v) > 1,
                |g, v| g.singleton_edge(g.incoming_edges(v)),
                true,
            )
        } else {
            self.find_path(
                vertex,
                prune_factor,
                |g, v| g.has_incident_ref_edge(v) || g.in_degree(v) == 0,
                |g, v| g.out_degree(v) > 1 && g.has_incident_ref_edge(v),
                |g, v| g.heaviest_incoming_edge(v),
                true,
            )
        }
    }

    /// `findPathDownwardsToHighestCommonDescendantOfReference`.
    fn find_path_downwards_to_highest_common_descendant_of_reference(
        &self,
        vertex: usize,
        prune_factor: i32,
        give_up_at_branch: bool,
    ) -> Result<Option<Vec<usize>>, GraphError> {
        if give_up_at_branch {
            self.find_path(
                vertex,
                prune_factor,
                |g, v| g.is_reference_node(v) || g.out_degree(v) != 1,
                |g, v| g.is_reference_node(v),
                |g, v| g.singleton_edge(g.outgoing_edges(v)),
                false,
            )
        } else {
            self.find_path(
                vertex,
                prune_factor,
                |g, v| g.is_reference_node(v) || g.out_degree(v) == 0,
                |g, v| g.is_reference_node(v),
                |g, v| g.heaviest_outgoing_edge(v),
                false,
            )
        }
    }

    /// `findPath`: step from `vertex` until `done`; an edge under the prune factor drops the path
    /// so far (its vertices still count as visited), and a vertex seen twice gives up. The path
    /// is built front first, so it starts at the vertex where the walk stopped.
    fn find_path(
        &self,
        vertex: usize,
        prune_factor: i32,
        done: impl Fn(&Self, usize) -> bool,
        return_path: impl Fn(&Self, usize) -> bool,
        next_edge: impl Fn(&Self, usize) -> Result<usize, GraphError>,
        upwards: bool,
    ) -> Result<Option<Vec<usize>>, GraphError> {
        let mut path: Vec<usize> = Vec::new();
        let mut visited = vec![false; self.vertices.len()];
        let mut v = vertex;
        while !done(self, v) {
            let edge = next_edge(self, v)?;
            if self.edges[edge].pruning_multiplicity() < prune_factor {
                for &seen in &path {
                    visited[seen] = true;
                }
                path.clear();
            } else {
                path.push(v);
            }
            let Edge { source, target, .. } = self.edges[edge];
            v = if upwards { source } else { target };
            // The reference prints "Dangling End recovery killed because of a loop" to stderr.
            if path.contains(&v) || visited[v] {
                return Ok(None);
            }
        }
        path.push(v);
        path.reverse();
        Ok(if return_path(self, v) {
            Some(path)
        } else {
            None
        })
    }

    /// `getSingletonEdge`: the one edge, `validateArg` past one. None is a `null` the caller then
    /// dereferences, which `findPath`'s stopping rules never reach.
    fn singleton_edge(&self, edges: &[usize]) -> Result<usize, GraphError> {
        match edges {
            [edge] => Ok(*edge),
            [] => Err(GraphError::NoSuchElement),
            _ => Err(GraphError::IllegalArgument(
                "Cannot get a single incoming edge for a vertex with multiple incoming edges",
            )),
        }
    }

    /// `hasIncidentRefEdge`: an INCOMING reference edge.
    fn has_incident_ref_edge(&self, vertex: usize) -> bool {
        self.incoming[vertex].iter().any(|&e| self.edges[e].is_ref)
    }

    /// `getHeaviestIncomingEdge`: `Stream.max`, which keeps the first of equal multiplicities.
    fn heaviest_incoming_edge(&self, vertex: usize) -> Result<usize, GraphError> {
        self.heaviest(&self.incoming[vertex])
    }

    /// `getHeaviestOutgoingEdge`.
    fn heaviest_outgoing_edge(&self, vertex: usize) -> Result<usize, GraphError> {
        self.heaviest(&self.outgoing[vertex])
    }

    fn heaviest(&self, edges: &[usize]) -> Result<usize, GraphError> {
        let mut best: Option<usize> = None;
        for &edge in edges {
            if best.is_none_or(|b| self.edges[edge].multiplicity > self.edges[b].multiplicity) {
                best = Some(edge);
            }
        }
        best.ok_or(GraphError::NoSuchElement)
    }

    /// `BaseGraph.isReferenceNode`: a reference edge in or out, or the only vertex.
    fn is_reference_node(&self, vertex: usize) -> bool {
        self.incoming[vertex]
            .iter()
            .chain(self.outgoing[vertex].iter())
            .any(|&e| self.edges[e].is_ref)
            || self.vertex_count() == 1
    }

    /// `getReferencePath(start, downwards, blacklistedEdge)`: `getNextReferenceVertex(v, true,
    /// blacklistedEdge)` until there is none.
    fn reference_path_downwards(&self, start: usize, blacklisted: Option<usize>) -> Vec<usize> {
        let mut path = Vec::new();
        let mut v = Some(start);
        while let Some(vertex) = v {
            path.push(vertex);
            v = self.next_reference_vertex(vertex, blacklisted);
        }
        path
    }

    /// `getReferencePath(start, upwards, empty)`: `getPrevReferenceVertex` until there is none.
    fn reference_path_upwards(&self, start: usize) -> Vec<usize> {
        let mut path = Vec::new();
        let mut v = Some(start);
        while let Some(vertex) = v {
            path.push(vertex);
            v = self.incoming[vertex]
                .iter()
                .map(|&e| self.edges[e].source)
                .find(|&source| self.is_reference_node(source));
        }
        path
    }

    /// `getNextReferenceVertex(v, true, blacklistedEdge)`: the target of the first reference edge
    /// out, else of the only edge out that is not blacklisted.
    fn next_reference_vertex(&self, vertex: usize, blacklisted: Option<usize>) -> Option<usize> {
        let outgoing = &self.outgoing[vertex];
        if let Some(&edge) = outgoing.iter().find(|&&e| self.edges[e].is_ref) {
            return Some(self.edges[edge].target);
        }
        let mut allowed = outgoing.iter().filter(|&&e| Some(e) != blacklisted);
        match (allowed.next(), allowed.next()) {
            (Some(&edge), None) => Some(self.edges[edge].target),
            _ => None,
        }
    }

    /// `getBasesForPath`: each vertex's last base, or, with `expand_source`, a source's whole
    /// k-mer reversed.
    fn bases_for_path(&self, path: &[usize], expand_source: bool) -> Vec<u8> {
        let mut bases = Vec::new();
        for &v in path {
            if expand_source && self.in_degree(v) == 0 {
                bases.extend(self.vertices[v].sequence.iter().rev());
            } else {
                bases.push(self.vertices[v].suffix());
            }
        }
        bases
    }

    /// `addVertex(new MultiDeBruijnVertex(bases))`: a vertex the k-mer map never hears of.
    fn add_vertex_outside_kmer_map(&mut self, sequence: Vec<u8>) -> usize {
        let vertex = self.vertices.len();
        self.vertices.push(Vertex {
            sequence,
            additional_info: String::new(),
        });
        self.vertex_alive.push(true);
        self.outgoing.push(Vec::new());
        self.incoming.push(Vec::new());
        vertex
    }

    /// `addEdge(source, target, createEdge(false, 1))`: JGraphT refuses a second edge between the
    /// same two vertices.
    fn add_edge_unless_present(&mut self, source: usize, target: usize) {
        if self.outgoing[source]
            .iter()
            .any(|&e| self.edges[e].target == target)
        {
            return;
        }
        self.add_edge(source, target, false, 1);
    }
}

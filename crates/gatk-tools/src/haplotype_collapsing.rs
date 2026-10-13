//! Ported from `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.
//! LongHomopolymerHaplotypeCollapsingEngine` (GATK 4.6.2.0): `--flow-assembly-collapse-hmer-size`,
//! the haplotypes put back to the reference's hmer lengths once the reads (flow reads, whose hmers
//! are capped) have been scored against them.
//!
//! An engine exists for a region only when its reference haplotype holds a run longer than the
//! threshold. Each non-reference haplotype is then aligned to the reference under its span, both
//! ways round (reverse-complemented for the second), and wherever the alignment deletes reference
//! bases next to an hmer of the threshold's length in the haplotype, on an hmer as long in the
//! reference, the deleted bases are put back (in partial mode only while they are the hmer's
//! base); the longer of the two results is kept and marked collapsed when anything was put back.
//! Every haplotype is realigned to the reference haplotype, and haplotypes that became identical
//! are merged, their likelihoods the best of the group's, the reference's group under it.

use gatk_engine::haplotype::Haplotype;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::smith_waterman::{SmithWatermanAligner, SwOverhangStrategy, SwParameters};
use htsjdk_bam::cigar::Op;

/// `DETERMINE_COLLAPSE_THRESHOLD`: the threshold is the largest maximal class of the read groups.
pub const DETERMINE_COLLAPSE_THRESHOLD: i32 = -1;

/// `LongHomopolymerHaplotypeCollapsingEngine`.
#[derive(Debug, Clone)]
pub struct HaplotypeCollapsing {
    pub hmer_size_threshold: i32,
    pub partial_mode: bool,
    pub full_ref: Vec<u8>,
    pub ref_loc: SimpleInterval,
    pub parameters: SwParameters,
}

struct Uncollapsed {
    bases: Vec<u8>,
    offset: i32,
    collapsed: bool,
}

/// `needsCollapsing(bases, threshold)`: a run of more than `threshold` of one base.
pub fn needs_collapsing(bases: &[u8], threshold: i32) -> bool {
    let mut last = 0u8;
    let mut same = 0i32;
    for &base in bases {
        if base == last {
            same += 1;
            if same >= threshold {
                return true;
            }
        } else {
            last = base;
            same = 0;
        }
    }
    false
}

/// `SequenceUtil.reverseComplement`, in place.
fn reverse_complement(bases: &mut [u8]) {
    bases.reverse();
    for base in bases.iter_mut() {
        *base = match *base {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            other => other,
        };
    }
}

const OUT_OF_BOUNDS: &str = "java.lang.ArrayIndexOutOfBoundsException";

impl HaplotypeCollapsing {
    /// `uncollapseHmersInHaplotypes(haplotypes, false, null)`.
    pub fn uncollapse_hmers_in_haplotypes(
        &self,
        haplotypes: &[Haplotype],
        aligner: &dyn SmithWatermanAligner,
    ) -> Result<Vec<Haplotype>, String> {
        let mut reference: Option<(Vec<u8>, i32)> = None;
        for h in haplotypes {
            if h.is_reference() && reference.is_none() {
                reference = Some((h.bases(), h.alignment_start_hap_wrt_ref()));
            }
        }
        let mut result: Vec<Haplotype> = Vec::with_capacity(haplotypes.len());
        for h in haplotypes {
            let mut aligned = self.uncollapse_single(h, aligner)?;
            aligned.set_uniqueness_value(result.len() as i32);
            result.push(aligned);
        }
        if let Some((ref_bases, start)) = reference {
            for h in result.iter_mut().filter(|h| !h.is_reference()) {
                let alignment = aligner
                    .align(
                        &ref_bases,
                        &h.bases(),
                        &self.parameters,
                        SwOverhangStrategy::Indel,
                    )
                    .map_err(|e| format!("{e:?}"))?;
                h.set_cigar(&alignment.cigar)
                    .map_err(|e| format!("{e:?}"))?;
                h.set_alignment_start_hap_wrt_ref(alignment.alignment_offset + start);
            }
        }
        Ok(result)
    }

    /// `uncollapseSingleHaplotype`.
    fn uncollapse_single(
        &self,
        h: &Haplotype,
        aligner: &dyn SmithWatermanAligner,
    ) -> Result<Haplotype, String> {
        if h.is_reference() {
            return Ok(h.clone());
        }
        let location = h
            .genome_location()
            .ok_or("java.lang.NullPointerException")?
            .clone();
        let reference = self.uncollapsed_partial_ref(&location)?;
        let bases = h.bases();
        let forward = self.uncollapse_by_ref(&bases, &reference, false, aligner)?;
        let reverse = self.uncollapse_by_ref(&bases, &reference, true, aligner)?;
        let chosen = if reverse.bases.len() > forward.bases.len() {
            reverse
        } else {
            forward
        };
        let mut aligned =
            Haplotype::new(&chosen.bases, h.is_reference()).map_err(|e| format!("{e:?}"))?;
        aligned.set_score(h.score());
        aligned.set_genome_location(location);
        if let Some(map) = h.event_map() {
            aligned.set_event_map(map.clone());
        }
        aligned.set_alignment_start_hap_wrt_ref(chosen.offset);
        aligned.set_collapsed(chosen.collapsed);
        Ok(aligned)
    }

    /// `uncollapsedPartialRef`: the full reference under the haplotype's span.
    fn uncollapsed_partial_ref(&self, location: &SimpleInterval) -> Result<Vec<u8>, String> {
        let offset = location.start - self.ref_loc.start;
        let size = location.end - location.start + 1;
        if offset < 0 || offset as usize > self.full_ref.len() {
            return Err(OUT_OF_BOUNDS.to_string());
        }
        let from = offset as usize;
        let to = from + size.max(0) as usize;
        // `Arrays.copyOfRange` pads past the end with zeros.
        let mut out = self.full_ref[from..to.min(self.full_ref.len())].to_vec();
        out.resize(to - from, 0);
        Ok(out)
    }

    /// `uncollapseByRef(bases, ref, rev)`.
    fn uncollapse_by_ref(
        &self,
        bases_arg: &[u8],
        ref_arg: &[u8],
        rev: bool,
        aligner: &dyn SmithWatermanAligner,
    ) -> Result<Uncollapsed, String> {
        let mut bases = bases_arg.to_vec();
        let mut reference = ref_arg.to_vec();
        if rev {
            reverse_complement(&mut bases);
            reverse_complement(&mut reference);
        }
        let alignment = aligner
            .align(
                &reference,
                &bases,
                &self.parameters,
                SwOverhangStrategy::Indel,
            )
            .map_err(|e| format!("{e:?}"))?;
        let threshold = self.hmer_size_threshold;
        let deleted: u32 = alignment
            .cigar
            .elements
            .iter()
            .filter(|e| e.op == Op::D)
            .map(|e| e.length)
            .sum();
        let capacity = bases.len() + deleted as usize;
        let mut result: Vec<u8> = Vec::with_capacity(capacity);
        let mut bases_offset = alignment.alignment_offset;
        let mut ref_offset: i32 = 0;
        let mut collapsed = false;
        let at = |array: &[u8], i: i32| -> Result<u8, String> {
            if i < 0 {
                return Err(OUT_OF_BOUNDS.to_string());
            }
            array
                .get(i as usize)
                .copied()
                .ok_or_else(|| OUT_OF_BOUNDS.to_string())
        };
        for element in &alignment.cigar.elements {
            let length = element.length as i32;
            if element.op != Op::D {
                if element.op.consumes_read_bases() {
                    let from = bases_offset.max(0) as usize;
                    let to = from + length as usize;
                    if bases_offset < 0 || to > bases.len() || result.len() + (to - from) > capacity
                    {
                        return Err(OUT_OF_BOUNDS.to_string());
                    }
                    result.extend_from_slice(&bases[from..to]);
                    bases_offset += length;
                }
            } else {
                let fwd_end = (bases_offset + threshold).min(bases.len() as i32);
                let forward: &[u8] =
                    if bases_offset as usize <= bases.len() && bases_offset <= fwd_end {
                        &bases[bases_offset.max(0) as usize..fwd_end.max(0) as usize]
                    } else {
                        &[]
                    };
                let bck_start = (bases_offset - threshold).max(0);
                let backward: &[u8] =
                    &bases[bck_start as usize..(bases_offset.max(0) as usize).min(bases.len())];
                if needs_collapsing(forward, threshold - 1)
                    || needs_collapsing(backward, threshold - 1)
                {
                    let first = at(&reference, ref_offset)?;
                    if self.on_homopolymer(&reference, ref_offset - threshold, first, threshold) {
                        for size in 0..length {
                            if self.partial_mode && at(&reference, ref_offset + size)? != first {
                                break;
                            }
                            result.push(first);
                        }
                        collapsed = true;
                    } else {
                        let last = at(&reference, ref_offset + length - 1)?;
                        if self.on_homopolymer(&reference, ref_offset + length, last, threshold) {
                            for size in 0..length {
                                if self.partial_mode
                                    && at(&reference, ref_offset + length - 1 - size)? != last
                                {
                                    break;
                                }
                                result.push(last);
                            }
                            collapsed = true;
                        }
                    }
                }
            }
            if element.op.consumes_reference_bases() {
                ref_offset += length;
            }
        }
        if rev {
            reverse_complement(&mut result);
        }
        Ok(Uncollapsed {
            bases: result,
            offset: alignment.alignment_offset,
            collapsed,
        })
    }

    /// `onHomoPolymer`: a run of `length` of the base starting at any of the `threshold` offsets
    /// from `offset`.
    fn on_homopolymer(&self, bases: &[u8], offset: i32, base: u8, length: i32) -> bool {
        (0..self.hmer_size_threshold).any(|tick| same_base(bases, offset + tick, base, length))
    }
}

/// `sameBase`: `length` of the base from `offset`, an index outside the array answering false.
fn same_base(bases: &[u8], offset: i32, base: u8, length: i32) -> bool {
    if offset + length > bases.len() as i32 {
        return false;
    }
    if offset < 0 {
        return false;
    }
    (0..length).all(|i| bases[(offset + i) as usize] == base)
}

/// `identicalBySequence`: the haplotypes grouped by their bases, each group under its first, and
/// the reference's group moved under the reference (to the end of the map) when it is not its
/// first.
pub fn identical_by_sequence(haplotypes: &[Haplotype]) -> Result<Vec<(usize, Vec<usize>)>, String> {
    let mut groups: Vec<(Vec<u8>, Vec<usize>)> = Vec::new();
    for (index, h) in haplotypes.iter().enumerate() {
        let bases = h.bases();
        match groups.iter_mut().find(|(b, _)| *b == bases) {
            Some((_, members)) => members.push(index),
            None => groups.push((bases, vec![index])),
        }
    }
    let mut result: Vec<(usize, Vec<usize>)> = groups
        .into_iter()
        .map(|(_, members)| (members[0], members))
        .collect();
    let reference = haplotypes.iter().position(Haplotype::is_reference).ok_or(
        "java.lang.IllegalArgumentException: Reference haplotype missing from the list of alleles",
    )?;
    if !result.iter().any(|(key, _)| *key == reference) {
        if let Some(at) = result
            .iter()
            .position(|(_, members)| members.contains(&reference))
        {
            let (_, members) = result.remove(at);
            result.push((reference, members));
        }
    }
    Ok(result)
}

/// `determineFlowAssemblyColapseHmer`: the largest maximal class over the header's read groups, a
/// flow group's `mc` or twelve, any other group's zero.
pub fn determine_collapse_threshold(
    header: &htsjdk_bam::header::SamHeader,
) -> Result<i32, (&'static str, String)> {
    let mut result = 0;
    for group in &header.read_groups {
        let platform = group
            .attributes
            .get("PL")
            .and_then(|pl| crate::flow_based_read::platform(pl));
        let max_class = match platform {
            Some("ULTIMA") if group.attributes.get("FO").is_none() => {
                return Err((
                    "java.lang.RuntimeException",
                    format!(
                        "Malformed Ultima read group identified, aborting: SAMReadGroupRecord{{ID: {}}}",
                        group.id
                    ),
                ))
            }
            Some(_) => match group.attributes.get("mc") {
                None => crate::flow_based_read::MAX_CLASS,
                Some(text) => text.parse::<i32>().map_err(|_| {
                    (
                        "java.lang.NumberFormatException",
                        format!("For input string: \"{text}\""),
                    )
                })?,
            },
            None => 0,
        };
        if max_class >= result {
            result = max_class;
        }
    }
    Ok(result)
}

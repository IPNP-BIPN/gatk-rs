//! `FilterAlignmentArtifacts`: whether a variant is an artefact of alignment.
//!
//! The tool has three steps: reassemble the reads supporting a variant into unitigs, realign the
//! unitigs, and filter the variant if they map somewhere else just as well. The middle step is
//! BWA's and the first is the Mutect2 assembler's; neither is ported. The two rules the tool owns
//! are, along with the filter decision they feed.

/// One read, reduced to what `supportsVariant` reads off it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Read {
    pub name: String,
    pub start: i32,
    /// Cigar as (operator, length).
    pub cigar: Vec<(char, i32)>,
    pub bases: Vec<u8>,
}

/// One variant, reduced to the same.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant {
    pub start: i32,
    pub reference: Vec<u8>,
    pub alternates: Vec<Vec<u8>>,
}

impl Variant {
    /// `isSNP`: every alternate the same single base length as the reference.
    pub fn is_snp(&self) -> bool {
        self.reference.len() == 1
            && !self.alternates.is_empty()
            && self.alternates.iter().all(|allele| allele.len() == 1)
    }
}

/// What `getReadIndexForReferenceCoordinate` answers: an index into the read's bases, and the
/// operator the coordinate landed in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReadIndex {
    pub offset: usize,
    pub operator: char,
}

/// `READ_INDEX_NOT_FOUND` is returned as `None` here.
///
/// Three things about this walk are not what they look like. It starts from the SOFT start, which
/// is the alignment start less any leading soft clip. A SOFT CLIP CONSUMES REFERENCE here, though
/// it consumes none anywhere else, so a coordinate past the last aligned base still lands inside a
/// trailing clip. And an INSERTION consumes none, so it can never bracket a coordinate: the
/// coordinate after an insertion falls in the following match, whose first read position is
/// already past the inserted bases. That last one is why an insertion cannot be supported by an
/// `I` at a tolerance of zero.
pub fn read_index_for_reference_coordinate(read: &Read, coordinate: i32) -> Option<ReadIndex> {
    let leading_clip = match read.cigar.first() {
        Some(('S', length)) => *length,
        _ => 0,
    };
    let soft_start = read.start - leading_clip;
    if coordinate < soft_start {
        return None;
    }
    let mut last_read_position = 0usize;
    let mut last_reference_position = soft_start;
    for (operator, length) in &read.cigar {
        let first_read_position = last_read_position;
        let first_reference_position = last_reference_position;
        let consumes_read = matches!(operator, 'M' | 'I' | 'S' | '=' | 'X');
        let consumes_reference = matches!(operator, 'M' | 'D' | 'N' | '=' | 'X' | 'S');
        if consumes_read {
            last_read_position += *length as usize;
        }
        if consumes_reference {
            last_reference_position += length;
        }
        if first_reference_position <= coordinate && coordinate < last_reference_position {
            let offset = first_read_position
                + if consumes_read {
                    (coordinate - first_reference_position) as usize
                } else {
                    0
                };
            return Some(ReadIndex {
                offset,
                operator: *operator,
            });
        }
    }
    None
}

fn might_support_deletion(operator: char) -> bool {
    operator == 'D' || operator == 'S'
}

fn might_support_insertion(operator: char) -> bool {
    operator == 'I' || operator == 'S'
}

/// `RealignmentEngine.supportsVariant`.
///
/// A SNP is matched by BASES from the offset the coordinate maps to, and a coordinate that fell in
/// a deletion never supports one. An indel is matched by CIGAR OPERATOR instead, because indel
/// representation is not unique.
///
/// The indel walk keeps a running sum of element lengths and compares it against the variant's
/// offset. The sum is advanced ONLY when the element is outside the tolerance, so as soon as one
/// element is within it the sum freezes and every element after it is treated as being at the
/// variant, however far away it really is.
pub fn supports_variant(read: &Read, variant: &Variant, indel_start_tolerance: i32) -> bool {
    let Some(index) = read_index_for_reference_coordinate(read, variant.start) else {
        return false;
    };
    if variant.is_snp() && index.operator == 'D' {
        return false;
    }
    let variant_position = index.offset;

    for allele in &variant.alternates {
        let reference_length = variant.reference.len();
        if allele.len() == reference_length {
            // A SNP or an MNP: the read's bases from the offset, truncated at its end.
            let end = (variant_position + allele.len()).min(read.bases.len());
            if variant_position <= end && read.bases[variant_position..end] == allele[..] {
                return true;
            }
        } else {
            let is_deletion = allele.len() < reference_length;
            let mut read_position: i32 = 0;
            for (operator, length) in &read.cigar {
                if (read_position - variant_position as i32).abs() <= indel_start_tolerance {
                    if (is_deletion && might_support_deletion(*operator))
                        || (!is_deletion && might_support_insertion(*operator))
                    {
                        return true;
                    }
                    // The sum is NOT advanced here, which is the whole of the quirk.
                } else {
                    read_position += length;
                }
            }
        }
    }
    false
}

/// One realignment of one unitig, with only the fields the joint rule reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    pub reference_id: i32,
    pub reference_start: i32,
    pub reference_end: i32,
    pub score: i32,
    pub mismatches: i32,
    pub reverse_strand: bool,
}

impl Alignment {
    /// `convertToInterval`, as the pair the overlap detector compares.
    fn interval(&self) -> (i32, i32, i32) {
        (self.reference_id, self.reference_start, self.reference_end)
    }
}

/// Two intervals overlap once each has been padded by half the maximum fragment length on each
/// side. The reference builds the detector with a NEGATIVE padding of that size, which widens the
/// query rather than narrowing it.
fn overlaps(a: &Alignment, b: &Alignment, max_reasonable_fragment_length: i32) -> bool {
    let padding = max_reasonable_fragment_length / 2;
    let (contig_a, start_a, end_a) = a.interval();
    let (contig_b, start_b, end_b) = b.interval();
    contig_a == contig_b && start_a - padding <= end_b && start_b - padding <= end_a
}

/// `findJointAlignments`.
///
/// With no unitigs there is nothing; with one, every alignment is its own joint alignment. With
/// more, a joint alignment is one in which EVERY unitig has a same-strand alignment overlapping
/// the first unitig's, and the one kept per unitig is the best-scoring of those.
///
/// The reference collects the result in a `HashSet<List<BwaMemAlignment>>` whose elements define no
/// `equals`, so its order is identity-hash order and is not reproducible. The order here is the
/// first unitig's, and any comparison against the reference has to sort.
pub fn find_joint_alignments(
    unitig_alignments: &[Vec<Alignment>],
    max_reasonable_fragment_length: i32,
) -> Vec<Vec<Alignment>> {
    if unitig_alignments.is_empty() {
        return Vec::new();
    }
    if unitig_alignments.len() == 1 {
        return unitig_alignments[0]
            .iter()
            .map(|alignment| vec![alignment.clone()])
            .collect();
    }

    // It is cheaper to start from the unitig with the fewest alignments, and the reference sorts
    // its argument in place to do so. The sort is stable, so equal sizes keep their order.
    let mut ordered: Vec<Vec<Alignment>> = unitig_alignments.to_vec();
    ordered.sort_by_key(|alignments| alignments.len());

    let mut out: Vec<Vec<Alignment>> = Vec::new();
    for alignment in &ordered[0] {
        let reaches_all = ordered.iter().all(|unitig| {
            unitig.iter().any(|other| {
                other.reverse_strand == alignment.reverse_strand
                    && overlaps(alignment, other, max_reasonable_fragment_length)
            })
        });
        if !reaches_all {
            continue;
        }
        let group: Vec<Alignment> = ordered
            .iter()
            .map(|unitig| {
                unitig
                    .iter()
                    .filter(|other| {
                        other.reverse_strand == alignment.reverse_strand
                            && overlaps(alignment, other, max_reasonable_fragment_length)
                    })
                    .max_by_key(|other| other.score)
                    .expect("a reachable alignment")
                    .clone()
            })
            .collect();
        if !out.contains(&group) {
            out.push(group);
        }
    }
    out
}

/// `jointAlignmentScore`: the sum over the group.
pub fn joint_alignment_score(group: &[Alignment]) -> i32 {
    group.iter().map(|alignment| alignment.score).sum()
}

/// `totalMismatches`: the sum over the group.
pub fn total_mismatches(group: &[Alignment]) -> i32 {
    group.iter().map(|alignment| alignment.mismatches).sum()
}

/// What one variant's decision produced.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Decision {
    pub filtered: bool,
    /// Written only when there was a second joint alignment to compare against.
    pub score_difference: Option<i32>,
    pub joint_alignment_count: usize,
}

/// The filter itself.
///
/// The best joint alignment landing on ANOTHER contig filters outright. Otherwise a second joint
/// alignment filters only when BOTH the score difference and the mismatch difference, each per
/// base of unitig, fall below their thresholds: either one on its own leaves the variant alone.
pub fn decide(
    joint_alignments: &[Vec<Alignment>],
    variant_contig_id: i32,
    total_unitig_bases: i32,
    min_aligner_score_difference_per_base: f64,
    min_mismatch_difference_per_base: f64,
) -> Decision {
    let mut sorted: Vec<&Vec<Alignment>> = joint_alignments.iter().collect();
    sorted.sort_by_key(|group| -joint_alignment_score(group));

    let count = joint_alignments.len();
    if let Some(best) = sorted.first() {
        if best[0].reference_id != variant_contig_id {
            return Decision {
                filtered: true,
                score_difference: None,
                joint_alignment_count: count,
            };
        }
    }
    if sorted.len() > 1 {
        let score_difference = joint_alignment_score(sorted[0]) - joint_alignment_score(sorted[1]);
        let mismatch_difference = total_mismatches(sorted[1]) - total_mismatches(sorted[0]);
        let multimapping = (score_difference as f64) / (total_unitig_bases as f64)
            < min_aligner_score_difference_per_base
            && (mismatch_difference as f64) / (total_unitig_bases as f64)
                < min_mismatch_difference_per_base;
        return Decision {
            filtered: multimapping,
            score_difference: Some(score_difference),
            joint_alignment_count: count,
        };
    }
    Decision {
        filtered: false,
        score_difference: None,
        joint_alignment_count: count,
    }
}

/// The filter's own name in the output VCF.
pub const ALIGNMENT_ARTIFACT_FILTER_NAME: &str = "alignment_artifact";

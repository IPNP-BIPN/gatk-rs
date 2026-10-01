//! What a contig that aligns in exactly two pieces says about the reference: where the two
//! breakpoints are, what is inserted or repeated between them, and which structural variant
//! records that makes.
//!
//! Ported from `org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SimpleChimera`,
//! `org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.TypeInferredFromSimpleChimera`,
//! `org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.BreakpointComplications`,
//! `org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.BreakpointsInference`,
//! `org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.NovelAdjacencyAndAltHaplotype`,
//! `org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType`,
//! `org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType` and
//! `org.broadinstitute.hellbender.tools.spark.sv.discovery.BreakEndVariantType`.

use crate::sv_contig_alignments::{
    has_incomplete_picture_from_two_alignments, overlap_on_ref_span, AlignmentInterval, Cigar,
    CigarOp, Dictionary, Interval, SvError, GATK_EXCEPTION, ILLEGAL_ARGUMENT,
    SHOULD_NEVER_REACH_HERE, STRUCTURAL_VARIANT_SIZE_LOWER_BOUND, UNSUPPORTED_OPERATION,
};
use std::cmp::Ordering;

/// `StrandSwitch`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrandSwitch {
    NoSwitch,
    ForwardToReverse,
    ReverseToForward,
}

impl StrandSwitch {
    fn name(self) -> &'static str {
        match self {
            StrandSwitch::NoSwitch => "NO_SWITCH",
            StrandSwitch::ForwardToReverse => "FORWARD_TO_REVERSE",
            StrandSwitch::ReverseToForward => "REVERSE_TO_FORWARD",
        }
    }
}

/// `TypeInferredFromSimpleChimera`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InferredType {
    SimpleDel,
    DelDupContraction,
    SimpleIns,
    Rpl,
    SmallDupExpansion,
    SmallDupCpx,
    IntraChrStrandSwitch55,
    IntraChrStrandSwitch33,
    IntraChrRefOrderSwap,
    InterChrStrandSwitch55,
    InterChrStrandSwitch33,
    InterChrNoSsWithLeftMateFirstInPartner,
    InterChrNoSsWithLeftMateSecondInPartner,
}

/// `SequenceUtil.reverseComplement`, which complements the four bases in either case and leaves
/// anything else as it is.
pub fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            other => *other,
        })
        .collect()
}

fn text(bases: &[u8]) -> String {
    String::from_utf8_lossy(bases).into_owned()
}

/// `SimpleChimera.SPLIT_PAIR_MIN_ALIGNMENT_MQ` and `SPLIT_PAIR_MIN_ALIGNMENT_LENGTH`.
const SPLIT_PAIR_MIN_ALIGNMENT_MQ: i32 = 20;
const SPLIT_PAIR_MIN_ALIGNMENT_LENGTH: i32 = 30;

/// `SimpleChimera.splitPairStrongEnoughEvidenceForCA`.
pub fn split_pair_strong_enough(one: &AlignmentInterval, two: &AlignmentInterval) -> bool {
    if one.map_qual < SPLIT_PAIR_MIN_ALIGNMENT_MQ || two.map_qual < SPLIT_PAIR_MIN_ALIGNMENT_MQ {
        return false;
    }
    one.size_on_read().min(two.size_on_read()) >= SPLIT_PAIR_MIN_ALIGNMENT_LENGTH
}

/// `SimpleChimera.DistancesBetweenAlignmentsOnRefAndOnRead`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Distances {
    gap_on_ref: i32,
    gap_on_contig: i32,
    left_ref_end: i32,
    right_ref_start: i32,
    first_contig_end: i32,
    second_contig_start: i32,
}

/// `SimpleChimera`: the two alignments of a contig, in contig order, with the strand switch
/// between them and whether the contig reads along the reference's forward strand.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SimpleChimera {
    pub source_contig_name: String,
    pub lower: AlignmentInterval,
    pub higher: AlignmentInterval,
    pub strand_switch: StrandSwitch,
    pub forward_representation: bool,
    pub insertion_mappings: Vec<String>,
    pub non_canonical_sa_tag: String,
}

impl SimpleChimera {
    pub fn new(
        lower: AlignmentInterval,
        higher: AlignmentInterval,
        insertion_mappings: Vec<String>,
        source_contig_name: &str,
        non_canonical_sa_tag: &str,
        dictionary: &Dictionary,
    ) -> Result<SimpleChimera, SvError> {
        let strand_switch = if lower.forward_strand == higher.forward_strand {
            StrandSwitch::NoSwitch
        } else if lower.forward_strand {
            StrandSwitch::ForwardToReverse
        } else {
            StrandSwitch::ReverseToForward
        };
        // `isForwardStrandRepresentation`.
        let forward_representation = if lower.reference_span.contig == higher.reference_span.contig
        {
            match strand_switch {
                StrandSwitch::NoSwitch => lower.forward_strand,
                StrandSwitch::ForwardToReverse => lower.reference_span.end < higher.reference_span.end,
                StrandSwitch::ReverseToForward => {
                    lower.reference_span.start < higher.reference_span.start
                }
            }
        } else if strand_switch == StrandSwitch::NoSwitch {
            lower.forward_strand
        } else {
            dictionary.compare_contigs(&lower.reference_span.contig, &higher.reference_span.contig)?
                == Ordering::Less
        };
        Ok(SimpleChimera {
            source_contig_name: source_contig_name.to_string(),
            lower,
            higher,
            strand_switch,
            forward_representation,
            insertion_mappings,
            non_canonical_sa_tag: non_canonical_sa_tag.to_string(),
        })
    }

    /// `toString`, which is what an exception about a chimera prints.
    pub fn to_java_string(&self) -> String {
        format!(
            "SimpleChimera{{sourceContigName='{}', regionWithLowerCoordOnContig={}, regionWithHigherCoordOnContig={}, strandSwitch={}, isForwardStrandRepresentation={}, insertionMappings=[{}], goodNonCanonicalMappingSATag='{}'}}",
            self.source_contig_name,
            self.lower.to_sa_tag_string(),
            self.higher.to_sa_tag_string(),
            self.strand_switch.name(),
            self.forward_representation,
            self.insertion_mappings.join(", "),
            self.non_canonical_sa_tag
        )
    }

    /// `firstContigRegionRefSpanAfterSecond`.
    fn first_after_second(&self, dictionary: &Dictionary) -> Result<bool, SvError> {
        Ok(dictionary
            .compare_locatables(&self.lower.reference_span, &self.higher.reference_span)?
            == Ordering::Greater)
    }

    /// `getCoordinateSortedRefSpans`.
    fn coordinate_sorted_spans(
        &self,
        dictionary: &Dictionary,
    ) -> Result<(Interval, Interval), SvError> {
        Ok(if self.first_after_second(dictionary)? {
            (
                self.higher.reference_span.clone(),
                self.lower.reference_span.clone(),
            )
        } else {
            (
                self.lower.reference_span.clone(),
                self.higher.reference_span.clone(),
            )
        })
    }

    fn has_incomplete_picture(&self) -> bool {
        has_incomplete_picture_from_two_alignments(&self.lower, &self.higher)
    }

    /// `isCandidateSimpleTranslocation`.
    fn is_candidate_simple_translocation(&self) -> bool {
        if self.has_incomplete_picture() {
            return false;
        }
        if self.lower.reference_span.contig != self.higher.reference_span.contig {
            return true;
        }
        if self.strand_switch != StrandSwitch::NoSwitch {
            return false;
        }
        let (one, two) = (&self.lower.reference_span, &self.higher.reference_span);
        if self.lower.forward_strand {
            one.start > two.end
        } else {
            two.start > one.end
        }
    }

    /// `isCandidateInvertedDuplication`.
    fn is_candidate_inverted_duplication(&self) -> bool {
        if self.lower.forward_strand == self.higher.forward_strand {
            return false;
        }
        2 * overlap_on_ref_span(&self.lower, &self.higher)
            > self.lower.size_on_read().min(self.higher.size_on_read())
    }

    /// `inferType`.
    pub fn infer_type(&self, dictionary: &Dictionary) -> Result<InferredType, SvError> {
        if self.is_candidate_simple_translocation() {
            if self.lower.reference_span.contig == self.higher.reference_span.contig {
                return Ok(InferredType::IntraChrRefOrderSwap);
            }
            return Ok(match self.strand_switch {
                StrandSwitch::ForwardToReverse => InferredType::InterChrStrandSwitch55,
                StrandSwitch::ReverseToForward => InferredType::InterChrStrandSwitch33,
                StrandSwitch::NoSwitch => {
                    if self.forward_representation != self.first_after_second(dictionary)? {
                        InferredType::InterChrNoSsWithLeftMateFirstInPartner
                    } else {
                        InferredType::InterChrNoSsWithLeftMateSecondInPartner
                    }
                }
            });
        }
        match self.strand_switch {
            StrandSwitch::ForwardToReverse => return Ok(InferredType::IntraChrStrandSwitch55),
            StrandSwitch::ReverseToForward => return Ok(InferredType::IntraChrStrandSwitch33),
            StrandSwitch::NoSwitch => {}
        }
        let distances = self.distances()?;
        let (on_ref, on_contig) = (distances.gap_on_ref, distances.gap_on_contig);
        Ok(match on_ref.cmp(&0) {
            Ordering::Greater => {
                if on_contig <= 0 {
                    InferredType::SimpleDel
                } else {
                    InferredType::Rpl
                }
            }
            Ordering::Less => {
                if on_contig >= 0 {
                    InferredType::SmallDupExpansion
                } else {
                    InferredType::SmallDupCpx
                }
            }
            Ordering::Equal => match on_contig.cmp(&0) {
                Ordering::Greater => InferredType::SimpleIns,
                Ordering::Less => InferredType::DelDupContraction,
                Ordering::Equal => {
                    return Err(SvError::new(
                        SHOULD_NEVER_REACH_HERE,
                        format!(
                            "Detected badly parsed chimeric alignment for identifying SV breakpoints; no rearrangement found: {}",
                            self.to_java_string()
                        ),
                    ))
                }
            },
        })
    }

    /// `getDistancesBetweenAlignmentsOnRefAndOnRead`.
    fn distances(&self) -> Result<Distances, SvError> {
        let neither = !self.has_incomplete_picture() && !self.is_candidate_simple_translocation();
        if !(neither && self.strand_switch == StrandSwitch::NoSwitch) {
            return Err(SvError::new(
                UNSUPPORTED_OPERATION,
                format!(
                    "Assumption that the simple chimera is neither incomplete picture nor simple translocation is violated.\n{}",
                    self.to_java_string()
                ),
            ));
        }
        let (left, right) = if self.forward_representation {
            (&self.lower.reference_span, &self.higher.reference_span)
        } else {
            (&self.higher.reference_span, &self.lower.reference_span)
        };
        let (r1e, r2b) = (left.end, right.start);
        let (c1e, c2b) = (self.lower.end_in_contig, self.higher.start_in_contig);
        Ok(Distances {
            gap_on_ref: r2b - r1e - 1,
            gap_on_contig: c2b - c1e - 1,
            left_ref_end: r1e,
            right_ref_start: r2b,
            first_contig_end: c1e,
            second_contig_start: c2b,
        })
    }
}

/// The complications at a pair of breakpoints: `BreakpointComplications` and its subclasses, one
/// variant each. Equality is the subclass and every field, as the Java `equals` compares them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Complication {
    SimpleInsDel {
        homology: String,
        inserted: String,
    },
    SmallDupPrecise {
        homology: String,
        inserted: String,
        repeat_unit: Interval,
        repeats_on_ref: i32,
        repeats_on_contig: i32,
        strands_on_ref: Vec<bool>,
        strands_on_contig: Vec<bool>,
        cigars: Vec<String>,
    },
    SmallDupImprecise {
        homology: String,
        inserted: String,
        repeat_unit: Interval,
        repeats_on_ref: i32,
        repeats_on_contig: i32,
        strands_on_ref: Vec<bool>,
        strands_on_contig: Vec<bool>,
        affected: Interval,
    },
    IntraChrStrandSwitch {
        homology: String,
        inserted: String,
    },
    InvertedDuplication {
        homology: String,
        inserted: String,
        repeat_unit: Interval,
        strands_on_contig: Vec<bool>,
        inverted_trans_insertion: Option<Interval>,
    },
    IntraChrRefOrderSwap {
        homology: String,
        inserted: String,
    },
    InterChromosome {
        homology: String,
        inserted: String,
    },
}

/// `Strand.toString` over a list, joined as `DUP_ORIENTATIONS` writes it.
fn strands(list: &[bool]) -> String {
    list.iter().map(|&plus| if plus { '+' } else { '-' }).collect()
}

impl Complication {
    pub fn homology(&self) -> &str {
        match self {
            Complication::SimpleInsDel { homology, .. }
            | Complication::SmallDupPrecise { homology, .. }
            | Complication::SmallDupImprecise { homology, .. }
            | Complication::IntraChrStrandSwitch { homology, .. }
            | Complication::InvertedDuplication { homology, .. }
            | Complication::IntraChrRefOrderSwap { homology, .. }
            | Complication::InterChromosome { homology, .. } => homology,
        }
    }

    pub fn inserted(&self) -> &str {
        match self {
            Complication::SimpleInsDel { inserted, .. }
            | Complication::SmallDupPrecise { inserted, .. }
            | Complication::SmallDupImprecise { inserted, .. }
            | Complication::IntraChrStrandSwitch { inserted, .. }
            | Complication::InvertedDuplication { inserted, .. }
            | Complication::IntraChrRefOrderSwap { inserted, .. }
            | Complication::InterChromosome { inserted, .. } => inserted,
        }
    }

    /// `hasDuplicationAnnotation`, which the inverted duplication does NOT override: only the two
    /// small-duplication classes answer true.
    pub fn has_duplication_annotation(&self) -> bool {
        matches!(
            self,
            Complication::SmallDupPrecise { .. } | Complication::SmallDupImprecise { .. }
        )
    }

    /// The small-duplication fields: repeat unit, copies on the reference and on the contig.
    fn duplication(&self) -> Option<(&Interval, i32, i32)> {
        match self {
            Complication::SmallDupPrecise {
                repeat_unit,
                repeats_on_ref,
                repeats_on_contig,
                ..
            }
            | Complication::SmallDupImprecise {
                repeat_unit,
                repeats_on_ref,
                repeats_on_contig,
                ..
            } => Some((repeat_unit, *repeats_on_ref, *repeats_on_contig)),
            _ => None,
        }
    }

    fn is_dup_contraction(&self) -> bool {
        self.duplication()
            .is_some_and(|(_, on_ref, on_contig)| on_ref > on_contig)
    }

    /// `toVariantAttributes`, each value as the VCF writes it and a flag as the empty string.
    pub fn to_variant_attributes(&self) -> Vec<(&'static str, String)> {
        let mut attributes = Vec::new();
        if !self.inserted().is_empty() {
            attributes.push(("INSSEQ", self.inserted().to_string()));
            attributes.push(("INSLEN", self.inserted().len().to_string()));
        }
        if !self.homology().is_empty() {
            attributes.push(("HOMSEQ", self.homology().to_string()));
            attributes.push(("HOMLEN", self.homology().len().to_string()));
        }
        let small = match self {
            Complication::SmallDupPrecise {
                repeat_unit,
                repeats_on_ref,
                repeats_on_contig,
                strands_on_contig,
                ..
            }
            | Complication::SmallDupImprecise {
                repeat_unit,
                repeats_on_ref,
                repeats_on_contig,
                strands_on_contig,
                ..
            } => Some((repeat_unit, repeats_on_ref, repeats_on_contig, strands_on_contig)),
            _ => None,
        };
        if let Some((unit, on_ref, on_contig, orientations)) = small {
            attributes.push(("DUP_REPEAT_UNIT_REF_SPAN", unit.to_string()));
            attributes.push(("DUP_NUM", format!("{on_ref},{on_contig}")));
            attributes.push(("DUP_ORIENTATIONS", strands(orientations)));
            attributes.push((
                if on_ref < on_contig {
                    "EXPANSION"
                } else {
                    "CONTRACTION"
                },
                String::new(),
            ));
        }
        match self {
            Complication::SmallDupPrecise { cigars, .. } if !cigars.is_empty() => {
                attributes.push(("DUP_SEQ_CIGARS", cigars.join(",")));
            }
            Complication::SmallDupImprecise { affected, .. } => {
                attributes.push(("DUP_ANNOTATIONS_IMPRECISE", String::new()));
                attributes.push(("DUP_IMPRECISE_AFFECTED_RANGE", affected.to_string()));
            }
            _ => {}
        }
        attributes
    }
}

/// `reverseComplementIfNecessary`: the sequence as the forward strand of the reference reads it.
fn forward_strand_representation(
    bases: &[u8],
    first: &AlignmentInterval,
    second: &AlignmentInterval,
    first_after_second: bool,
) -> Vec<u8> {
    let flip = if first.forward_strand == second.forward_strand {
        !first.forward_strand
    } else {
        first_after_second == first.forward_strand
    };
    if flip {
        reverse_complement(bases)
    } else {
        bases.to_vec()
    }
}

fn slice(sequence: &[u8], from: i32, to: i32) -> Result<&[u8], SvError> {
    if from < 0 || to < from || to as usize > sequence.len() {
        return Err(SvError::new(
            "java.lang.ArrayIndexOutOfBoundsException",
            format!(
                "Range [{from}, {to}) out of bounds for length {}",
                sequence.len()
            ),
        ));
    }
    Ok(&sequence[from as usize..to as usize])
}

/// `BreakpointComplications.inferHomology`: the contig bases both alignments claim.
fn infer_homology(
    first: &AlignmentInterval,
    second: &AlignmentInterval,
    contig: &[u8],
    first_after_second: bool,
) -> Result<String, SvError> {
    if first.end_in_contig >= second.start_in_contig {
        let bytes = slice(contig, second.start_in_contig - 1, first.end_in_contig)?;
        Ok(text(&forward_strand_representation(
            bytes,
            first,
            second,
            first_after_second,
        )))
    } else {
        Ok(String::new())
    }
}

/// `BreakpointComplications.inferInsertedSequence`: the contig bases neither alignment claims.
fn infer_inserted_sequence(
    first: &AlignmentInterval,
    second: &AlignmentInterval,
    contig: &[u8],
    first_after_second: bool,
) -> Result<String, SvError> {
    if first.end_in_contig < second.start_in_contig - 1 {
        let bytes = slice(contig, first.end_in_contig, second.start_in_contig - 1)?;
        Ok(text(&forward_strand_representation(
            bytes,
            first,
            second,
            first_after_second,
        )))
    } else {
        Ok(String::new())
    }
}

/// `SmallDuplicationWithPreciseDupRangeBreakpointComplications.extractCigarForTandupExpansion`:
/// the part of an alignment's cigar that lies over the duplicated reference span.
fn extract_cigar_for_tandup_expansion(
    region: &AlignmentInterval,
    one_end: i32,
    two_begin: i32,
) -> Cigar {
    let mut result = Vec::new();
    let forward = region.forward_strand;
    let mut initiated = false;
    let mut ref_pos = if forward {
        region.reference_span.start
    } else {
        region.reference_span.end
    };
    for &(length, op) in &region.cigar.0 {
        if op.is_clipping() {
            continue;
        }
        if op.consumes_reference_bases() {
            ref_pos += if forward { length } else { -length };
        }
        let offset_into_repeat = if forward {
            ref_pos - two_begin
        } else {
            one_end - ref_pos
        };
        let overshoot = if forward {
            ref_pos - one_end - 1
        } else {
            two_begin - ref_pos - 1
        };
        if offset_into_repeat > 0 {
            if overshoot <= 0 {
                result.push(if initiated {
                    (length, op)
                } else {
                    (offset_into_repeat, op)
                });
                initiated = true;
            } else {
                result.push((length - overshoot, op));
                break;
            }
        }
    }
    Cigar(result)
}

/// `SmallDuplicationWithImpreciseDupRangeBreakpointComplications.TandemRepeatStructure`: the
/// repeat-unit length and copy numbers that best explain two negative gaps.
fn tandem_repeat_structure(on_ref: i32, on_contig: i32) -> (i32, i32, i32, i32) {
    const MAX_LOWER_CN: i32 = 10;
    let expansion = on_ref < on_contig;
    let (lower_overlap, higher_overlap) = if expansion {
        (on_ref.abs(), on_contig.abs())
    } else {
        (on_contig.abs(), on_ref.abs())
    };
    let (mut higher, mut lower, mut unit, mut pseudo) = (0, 0, 0, 0);
    let mut error = f64::MAX;
    for cn2 in 1..MAX_LOWER_CN {
        for cn1 in cn2 + 1..=2 * cn2 {
            let bound = if cn1 == 2 * cn2 {
                lower_overlap
            } else {
                higher_overlap
            };
            for l in 2..=bound {
                for lambda in 0..l {
                    let d1 = (2 * cn2 - cn1) * l + lambda;
                    let d2 = cn2 * l + lambda;
                    let new_error =
                        ((higher_overlap - d1).abs() + (lower_overlap - d2).abs()) as f64;
                    if new_error < error {
                        error = new_error;
                        higher = cn1;
                        lower = cn2;
                        unit = l;
                        pseudo = lambda;
                    }
                    if error < 1.0 {
                        return (lower, higher, unit, pseudo);
                    }
                }
            }
        }
    }
    (lower, higher, unit, pseudo)
}

/// What `BreakpointsInference` concludes: the two left-justified breakpoints, the complications
/// and the alternate haplotype.
struct Inference {
    upstream: Interval,
    downstream: Interval,
    complication: Complication,
    alt_haplotype: Vec<u8>,
}

fn point(contig: &str, position: i32) -> Result<Interval, SvError> {
    Interval::new(contig, position, position)
}

/// `BreakpointsInference.validateInferredLocations`.
fn validate_locations(
    left: &Interval,
    right: &Interval,
    dictionary: &Dictionary,
) -> Result<(), SvError> {
    let order = dictionary.compare_contigs(&right.contig, &left.contig)?;
    if order == Ordering::Less || (order == Ordering::Equal && right.end < left.start) {
        return Err(SvError::new(
            SHOULD_NEVER_REACH_HERE,
            "Inferred novel adjacency reference locations have left location after right location.",
        ));
    }
    for (one, space) in [(left, ""), (right, " ")] {
        if dictionary.length(&one.contig).is_some_and(|length| one.end > length) {
            return Err(SvError::new(
                SHOULD_NEVER_REACH_HERE,
                format!("Inferred breakpoint beyond reference sequence length.{space}"),
            ));
        }
    }
    Ok(())
}

/// `BreakpointsInference.getInferenceClass` and the constructor of each class it picks.
fn infer_breakpoints(
    chimera: &SimpleChimera,
    contig: &[u8],
    dictionary: &Dictionary,
) -> Result<Inference, SvError> {
    let first_after_second = chimera.first_after_second(dictionary)?;
    let (lower, higher) = (&chimera.lower, &chimera.higher);
    let lower_contig = lower.reference_span.contig.clone();
    match chimera.infer_type(dictionary)? {
        InferredType::SimpleDel | InferredType::Rpl | InferredType::SimpleIns => {
            // `SimpleInsDelOrReplacementBreakpointComplications`.
            let distances = chimera.distances()?;
            let (mut homology, mut inserted) = (String::new(), String::new());
            if distances.gap_on_ref > 0 {
                if distances.gap_on_contig >= 0 {
                    inserted = infer_inserted_sequence(lower, higher, contig, first_after_second)?;
                } else {
                    homology = infer_homology(lower, higher, contig, first_after_second)?;
                }
            } else if distances.gap_on_ref == 0 && distances.gap_on_contig > 0 {
                inserted = infer_inserted_sequence(lower, higher, contig, first_after_second)?;
            } else {
                return Err(SvError::new(
                    SHOULD_NEVER_REACH_HERE,
                    format!(
                        "Inferring breakpoint complications with the wrong unit: using simple ins-del unit for simple chimera:\n{}",
                        chimera.to_java_string()
                    ),
                ));
            }
            if distances.gap_on_contig > 0 && inserted.is_empty() {
                return Err(SvError::new(
                    SHOULD_NEVER_REACH_HERE,
                    format!(
                        "An identified breakpoint pair seem to suggest insertion but the inserted sequence is empty: {}",
                        chimera.to_java_string()
                    ),
                ));
            }
            let (left, right) = chimera.coordinate_sorted_spans(dictionary)?;
            let upstream = point(&lower_contig, left.end - homology.len() as i32)?;
            let downstream = point(&lower_contig, right.start - 1)?;
            let alt_haplotype = inserted.as_bytes().to_vec();
            Ok(Inference {
                upstream,
                downstream,
                complication: Complication::SimpleInsDel { homology, inserted },
                alt_haplotype,
            })
        }
        InferredType::SmallDupExpansion | InferredType::DelDupContraction => {
            let distances = chimera.distances()?;
            let left_span = if chimera.forward_representation {
                &lower.reference_span
            } else {
                &higher.reference_span
            };
            let complication;
            let (mut homology, mut inserted) = (String::new(), String::new());
            if distances.gap_on_ref > 0 {
                return Err(SvError::new(
                    SHOULD_NEVER_REACH_HERE,
                    format!(
                        "Simple chimera being sent down the wrong path, where the signature indicates a simple deletion but complication being resolve for small duplication. \n{}",
                        chimera.to_java_string()
                    ),
                ));
            } else if distances.gap_on_ref < 0 {
                // `resolveComplicationForSimpleTandupExpansion`.
                if distances.gap_on_contig < 0 {
                    return Err(SvError::new(
                        SHOULD_NEVER_REACH_HERE,
                        format!(
                            "Simple chimera being sent down the wrong path, where the signature indicates complex duplication but complication being resolved for simple small duplication. \n{}",
                            chimera.to_java_string()
                        ),
                    ));
                }
                if distances.gap_on_contig != 0 {
                    inserted = infer_inserted_sequence(lower, higher, contig, first_after_second)?;
                }
                let mut cigars = if lower.forward_strand {
                    vec![
                        extract_cigar_for_tandup_expansion(
                            lower,
                            distances.left_ref_end,
                            distances.right_ref_start,
                        )
                        .to_string(),
                        extract_cigar_for_tandup_expansion(
                            higher,
                            distances.left_ref_end,
                            distances.right_ref_start,
                        )
                        .to_string(),
                    ]
                } else {
                    vec![
                        extract_cigar_for_tandup_expansion(
                            lower,
                            distances.left_ref_end,
                            distances.right_ref_start,
                        )
                        .inverted()
                        .to_string(),
                        extract_cigar_for_tandup_expansion(
                            higher,
                            distances.left_ref_end,
                            distances.right_ref_start,
                        )
                        .inverted()
                        .to_string(),
                    ]
                };
                if !lower.forward_strand {
                    cigars.reverse();
                }
                complication = Complication::SmallDupPrecise {
                    homology: String::new(),
                    inserted: inserted.clone(),
                    repeat_unit: Interval::new(
                        &left_span.contig,
                        distances.right_ref_start,
                        distances.left_ref_end,
                    )?,
                    repeats_on_ref: 1,
                    repeats_on_contig: 2,
                    strands_on_ref: vec![true],
                    strands_on_contig: vec![true, true],
                    cigars,
                };
            } else if distances.gap_on_contig < 0 {
                // `resolveComplicationForSimpleTandupContraction`.
                homology = infer_homology(lower, higher, contig, first_after_second)?;
                complication = Complication::SmallDupPrecise {
                    homology: homology.clone(),
                    inserted: String::new(),
                    repeat_unit: Interval::new(
                        &left_span.contig,
                        distances.left_ref_end
                            - (distances.first_contig_end - distances.second_contig_start),
                        distances.left_ref_end,
                    )?,
                    repeats_on_ref: 2,
                    repeats_on_contig: 1,
                    strands_on_ref: vec![true, true],
                    strands_on_contig: vec![true],
                    cigars: Vec::new(),
                };
            } else if distances.gap_on_contig > 0 {
                return Err(SvError::new(
                    SHOULD_NEVER_REACH_HERE,
                    format!(
                        "Simple chimera being sent down the wrong path, where the signature indicates an insertion but complication being resolve for small duplication. \n{}",
                        chimera.to_java_string()
                    ),
                ));
            } else {
                return Err(SvError::new(
                    SHOULD_NEVER_REACH_HERE,
                    format!(
                        "Detected badly parsed chimeric alignment for identifying SV breakpoints; no rearrangement found: {}",
                        chimera.to_java_string()
                    ),
                ));
            }
            let _ = inserted;
            let (left, right) = chimera.coordinate_sorted_spans(dictionary)?;
            let homology_length = homology.len() as i32;
            let (upstream_position, alt_haplotype) = if complication.is_dup_contraction() {
                (left.end - homology_length, Vec::new())
            } else {
                let (unit, on_ref, on_contig) = complication
                    .duplication()
                    .expect("a small duplication carries its repeat unit");
                let position = left.end - homology_length - (on_contig - on_ref) * unit.size();
                let cigars = match &complication {
                    Complication::SmallDupPrecise { cigars, .. } => cigars.clone(),
                    _ => Vec::new(),
                };
                let parse = |text: &str| Cigar::parse(text).unwrap_or_default();
                let (first_copy, second_copy) = if chimera.forward_representation {
                    (parse(&cigars[0]), parse(&cigars[1]))
                } else {
                    (parse(&cigars[1]), parse(&cigars[0]))
                };
                let start = distances.first_contig_end - first_copy.read_length();
                let end = distances.second_contig_start + second_copy.read_length() - 1;
                let mut bases = slice(contig, start, end)?.to_vec();
                if !chimera.forward_representation {
                    bases = reverse_complement(&bases);
                }
                (position, bases)
            };
            let upstream = point(&lower_contig, upstream_position)?;
            let downstream = point(&lower_contig, right.start - 1)?;
            validate_locations(&upstream, &downstream, dictionary)?;
            Ok(Inference {
                upstream,
                downstream,
                complication,
                alt_haplotype,
            })
        }
        InferredType::SmallDupCpx => {
            let distances = chimera.distances()?;
            if distances.gap_on_ref > 0 {
                return Err(SvError::new(
                    SHOULD_NEVER_REACH_HERE,
                    format!(
                        "Simple chimera being sent down the wrong path, where the signature indicates a simple deletion but complication being resolve for small duplication. \n{}",
                        chimera.to_java_string()
                    ),
                ));
            }
            if !(distances.gap_on_ref < 0 && distances.gap_on_contig < 0) {
                return Err(SvError::new(
                    SHOULD_NEVER_REACH_HERE,
                    format!(
                        "Simple chimera being sent down the wrong path, where the signature indicates simple duplication but complication being resolved for complex small duplication. \n{}",
                        chimera.to_java_string()
                    ),
                ));
            }
            let (lower_cn, higher_cn, unit_length, pseudo_homology) =
                tandem_repeat_structure(distances.gap_on_ref, distances.gap_on_contig);
            let expansion = distances.gap_on_ref < distances.gap_on_contig;
            let unit_start =
                distances.left_ref_end - pseudo_homology - unit_length * lower_cn + 1;
            let homology = infer_homology(lower, higher, contig, first_after_second)?;
            let (on_ref, on_contig) = if expansion {
                (lower_cn, higher_cn)
            } else {
                (higher_cn, lower_cn)
            };
            let (left_ref, right_ref) = if chimera.forward_representation {
                (&lower.reference_span, &higher.reference_span)
            } else {
                (&higher.reference_span, &lower.reference_span)
            };
            let affected = if expansion {
                Interval::new(
                    &left_ref.contig,
                    distances.right_ref_start,
                    distances.left_ref_end,
                )?
            } else {
                Interval::new(
                    &left_ref.contig,
                    right_ref.start - unit_length,
                    left_ref.end + unit_length,
                )?
            };
            let repeat_unit = Interval::new(
                &lower_contig,
                unit_start,
                unit_start + unit_length - 1,
            )?;
            let complication = Complication::SmallDupImprecise {
                homology: homology.clone(),
                inserted: String::new(),
                repeat_unit: repeat_unit.clone(),
                repeats_on_ref: on_ref,
                repeats_on_contig: on_contig,
                strands_on_ref: vec![true; on_ref.max(0) as usize],
                strands_on_contig: vec![true; on_contig.max(0) as usize],
                affected,
            };
            let (left, right) = chimera.coordinate_sorted_spans(dictionary)?;
            let homology_length = homology.len() as i32;
            let (upstream_position, alt_haplotype) = if on_ref > on_contig {
                let mut bases = slice(
                    contig,
                    distances.second_contig_start - 1,
                    distances.first_contig_end,
                )?
                .to_vec();
                if !chimera.forward_representation {
                    bases = reverse_complement(&bases);
                }
                (left.end - homology_length, bases)
            } else {
                let position =
                    left.end - homology_length - (on_contig - on_ref) * repeat_unit.size();
                let hard_clip_offset = |cigar: &Cigar| match cigar.0.first() {
                    Some((length, CigarOp::H)) => *length,
                    _ => 0,
                };
                let distance = lower.cigar.associated_distance_on_read(
                    distances.first_contig_end - hard_clip_offset(&lower.cigar),
                    -distances.gap_on_ref,
                    true,
                )?;
                let start = distances.first_contig_end - distance;
                let distance = higher.cigar.associated_distance_on_read(
                    distances.second_contig_start - hard_clip_offset(&higher.cigar),
                    -distances.gap_on_ref,
                    false,
                )?;
                let end = distances.second_contig_start + distance - 1;
                let mut bases = slice(contig, start, end)?.to_vec();
                if !chimera.forward_representation {
                    bases = reverse_complement(&bases);
                }
                (position, bases)
            };
            let upstream = point(&lower_contig, upstream_position)?;
            let downstream = point(&lower_contig, right.start - 1)?;
            validate_locations(&upstream, &downstream, dictionary)?;
            Ok(Inference {
                upstream,
                downstream,
                complication,
                alt_haplotype,
            })
        }
        InferredType::IntraChrStrandSwitch55 | InferredType::IntraChrStrandSwitch33 => {
            if chimera.is_candidate_inverted_duplication() {
                inverted_duplication(chimera, contig, dictionary, first_after_second)
            } else {
                let homology = infer_homology(lower, higher, contig, first_after_second)?;
                let inserted = infer_inserted_sequence(lower, higher, contig, first_after_second)?;
                let (left, right) = chimera.coordinate_sorted_spans(dictionary)?;
                let homology_length = homology.len() as i32;
                let (up, down) = if chimera.strand_switch == StrandSwitch::ForwardToReverse {
                    (left.end - homology_length, right.end)
                } else {
                    (left.start, right.start + homology_length)
                };
                let upstream = point(&lower_contig, up)?;
                let downstream = point(&lower_contig, down)?;
                validate_locations(&upstream, &downstream, dictionary)?;
                Ok(Inference {
                    upstream,
                    downstream,
                    complication: Complication::IntraChrStrandSwitch { homology, inserted },
                    alt_haplotype: Vec::new(),
                })
            }
        }
        InferredType::IntraChrRefOrderSwap => {
            let homology = infer_homology(lower, higher, contig, first_after_second)?;
            let inserted = infer_inserted_sequence(lower, higher, contig, first_after_second)?;
            let (left, right) = chimera.coordinate_sorted_spans(dictionary)?;
            let upstream = point(&lower_contig, left.start)?;
            let downstream = point(&lower_contig, right.end - homology.len() as i32)?;
            validate_locations(&upstream, &downstream, dictionary)?;
            Ok(Inference {
                upstream,
                downstream,
                complication: Complication::IntraChrRefOrderSwap { homology, inserted },
                alt_haplotype: Vec::new(),
            })
        }
        InferredType::InterChrStrandSwitch55
        | InferredType::InterChrStrandSwitch33
        | InferredType::InterChrNoSsWithLeftMateFirstInPartner
        | InferredType::InterChrNoSsWithLeftMateSecondInPartner => {
            let homology = infer_homology(lower, higher, contig, first_after_second)?;
            let inserted = infer_inserted_sequence(lower, higher, contig, first_after_second)?;
            let homology_length = homology.len() as i32;
            // `isFirstInPartner`.
            let first_in_partner = match chimera.strand_switch {
                StrandSwitch::NoSwitch => {
                    dictionary
                        .compare_contigs(&lower.reference_span.contig, &higher.reference_span.contig)?
                        == Ordering::Less
                }
                _ => chimera.forward_representation,
            };
            let (lo, hi) = (&lower.reference_span, &higher.reference_span);
            let (up_contig, down_contig) = if first_in_partner {
                (&lo.contig, &hi.contig)
            } else {
                (&hi.contig, &lo.contig)
            };
            let (up, down) = match (first_in_partner, chimera.strand_switch) {
                (true, StrandSwitch::NoSwitch) => {
                    if chimera.forward_representation {
                        (lo.end - homology_length, hi.start)
                    } else {
                        (lo.start, hi.end - homology_length)
                    }
                }
                (true, StrandSwitch::ForwardToReverse) => (lo.end - homology_length, hi.end),
                (true, StrandSwitch::ReverseToForward) => (lo.start, hi.start + homology_length),
                (false, StrandSwitch::NoSwitch) => {
                    if chimera.forward_representation {
                        (hi.start, lo.end - homology_length)
                    } else {
                        (hi.end - homology_length, lo.start)
                    }
                }
                (false, StrandSwitch::ForwardToReverse) => (hi.end - homology_length, lo.end),
                (false, StrandSwitch::ReverseToForward) => (hi.start, lo.start + homology_length),
            };
            let upstream = point(up_contig, up)?;
            let downstream = point(down_contig, down)?;
            validate_locations(&upstream, &downstream, dictionary)?;
            Ok(Inference {
                upstream,
                downstream,
                complication: Complication::InterChromosome { homology, inserted },
                alt_haplotype: Vec::new(),
            })
        }
    }
}

/// `InvertedDuplicationBreakpointsInference` with its complications.
fn inverted_duplication(
    chimera: &SimpleChimera,
    contig: &[u8],
    dictionary: &Dictionary,
    first_after_second: bool,
) -> Result<Inference, SvError> {
    let (first, second) = (&chimera.lower, &chimera.higher);
    let inserted = infer_inserted_sequence(first, second, contig, first_after_second)?;
    let jump_start = if first.forward_strand {
        first.reference_span.end
    } else {
        first.reference_span.start
    };
    let jump_landing = if second.forward_strand {
        second.reference_span.start
    } else {
        second.reference_span.end
    };
    let chromosome = &first.reference_span.contig;
    let (repeat_unit, inverted_trans_insertion, strands_on_contig) = if first.forward_strand {
        let (alpha, omega) = (first.reference_span.start, second.reference_span.start);
        let unit = Interval::new(chromosome, alpha.max(omega), jump_start.min(jump_landing))?;
        let trans = ((alpha <= omega && jump_start < jump_landing)
            || (alpha > omega && jump_landing < jump_start))
            .then(|| {
                Interval::new(
                    chromosome,
                    jump_start.min(jump_landing) + 1,
                    jump_start.max(jump_landing),
                )
            })
            .transpose()?;
        (unit, trans, vec![true, false])
    } else {
        let (alpha, omega) = (first.reference_span.end, second.reference_span.end);
        let unit = Interval::new(chromosome, jump_start.max(jump_landing), alpha.min(omega))?;
        let trans = ((alpha >= omega && jump_landing < jump_start)
            || (alpha < omega && jump_start < jump_landing))
            .then(|| {
                Interval::new(
                    chromosome,
                    jump_start.min(jump_landing) + 1,
                    jump_start.max(jump_landing),
                )
            })
            .transpose()?;
        (unit, trans, vec![false, true])
    };
    // `extractAltHaplotypeForInvDup`.
    let (start, end, reverse) = if first.forward_strand {
        let (alpha, omega) = (first.reference_span.start, second.reference_span.start);
        if alpha <= omega {
            let walk = if alpha == omega {
                0
            } else {
                first
                    .cigar
                    .associated_distance_on_read(first.start_in_contig, omega - alpha, false)?
            };
            (
                first.start_in_contig + walk - 1,
                second.end_in_contig,
                false,
            )
        } else {
            let walk =
                second
                    .cigar
                    .associated_distance_on_read(second.end_in_contig, alpha - omega, true)?;
            (first.start_in_contig - 1, second.end_in_contig - walk, true)
        }
    } else {
        let (alpha, omega) = (first.reference_span.end, second.reference_span.end);
        if alpha >= omega {
            let walk = if alpha == omega {
                0
            } else {
                first
                    .cigar
                    .associated_distance_on_read(first.start_in_contig, alpha - omega, false)?
            };
            (first.start_in_contig + walk - 1, second.end_in_contig, true)
        } else {
            let walk =
                second
                    .cigar
                    .associated_distance_on_read(second.end_in_contig, omega - alpha, true)?;
            (first.start_in_contig - 1, second.end_in_contig - walk, false)
        }
    };
    let mut alt_haplotype = slice(contig, start, end)?.to_vec();
    if reverse {
        alt_haplotype = reverse_complement(&alt_haplotype);
    }
    let lower_contig = &first.reference_span.contig;
    let upstream = point(lower_contig, repeat_unit.start - 1)?;
    let downstream = point(lower_contig, repeat_unit.end)?;
    validate_locations(&upstream, &downstream, dictionary)?;
    Ok(Inference {
        upstream,
        downstream,
        complication: Complication::InvertedDuplication {
            homology: String::new(),
            inserted,
            repeat_unit,
            strands_on_contig,
            inverted_trans_insertion,
        },
        alt_haplotype,
    })
}

/// `NovelAdjacencyAndAltHaplotype`: the two left-justified breakpoints, the strand switch, the
/// complications and the alternate haplotype. Equality ignores the inferred type, as the Java
/// `equals` does, so two contigs that infer it differently still merge.
#[derive(Debug, Clone)]
pub struct NovelAdjacency {
    pub left: Interval,
    pub right: Interval,
    pub strand_switch: StrandSwitch,
    pub complication: Complication,
    pub inferred_type: InferredType,
    pub alt_haplotype: Vec<u8>,
}

impl PartialEq for NovelAdjacency {
    fn eq(&self, other: &NovelAdjacency) -> bool {
        self.left == other.left
            && self.right == other.right
            && self.strand_switch == other.strand_switch
            && self.complication == other.complication
            && self.alt_haplotype == other.alt_haplotype
    }
}

impl NovelAdjacency {
    /// The constructor from a chimera, which wraps an `IllegalArgumentException` in a
    /// `GATKException` naming the chimera.
    pub fn new(
        chimera: &SimpleChimera,
        contig: &[u8],
        dictionary: &Dictionary,
    ) -> Result<NovelAdjacency, SvError> {
        let wrap = |error: SvError| {
            if error.class == ILLEGAL_ARGUMENT {
                SvError::new(
                    GATK_EXCEPTION,
                    format!(
                        "Erred when inferring breakpoint location and event type from chimeric alignment:\n{}",
                        chimera.to_java_string()
                    ),
                )
            } else {
                error
            }
        };
        let inference = infer_breakpoints(chimera, contig, dictionary).map_err(wrap)?;
        let inferred_type = chimera.infer_type(dictionary).map_err(wrap)?;
        Ok(NovelAdjacency {
            left: inference.upstream,
            right: inference.downstream,
            strand_switch: chimera.strand_switch,
            complication: inference.complication,
            inferred_type,
            alt_haplotype: inference.alt_haplotype,
        })
    }

    /// `getDistanceBetweenNovelAdjacencies`.
    fn distance(&self) -> i32 {
        if self.left.contig == self.right.contig {
            self.right.end - self.left.start
        } else {
            -1
        }
    }

    /// `isCandidateForFatInsertion`.
    fn is_fat_insertion(&self) -> bool {
        self.inferred_type == InferredType::Rpl
            && self.right.end - self.left.start < STRUCTURAL_VARIANT_SIZE_LOWER_BOUND
    }

    /// `getLengthForDupTandemExpansion`.
    fn length_for_dup_tandem_expansion(&self) -> Result<i32, SvError> {
        let Some((unit, on_ref, on_contig)) = self.complication.duplication() else {
            return Err(SvError::new(
                "java.lang.ClassCastException",
                "complication is not a small duplication",
            ));
        };
        if on_ref > on_contig {
            return Err(SvError::new(
                UNSUPPORTED_OPERATION,
                "Trying to extract length from a duplication contraction",
            ));
        }
        Ok(self.complication.inserted().len() as i32 + (on_contig - on_ref) * unit.size())
    }

    /// `makeLocationString(narl)`.
    fn location(&self) -> String {
        location_string(
            &self.left.contig,
            self.left.start,
            &self.right.contig,
            self.right.end,
        )
    }

    /// `toSimpleOrBNDTypes`: the one record, or the linked pair, this adjacency is written as.
    pub fn to_simple_or_bnd_types(
        &self,
        reference: &mut dyn FnMut(&Interval) -> Result<Vec<u8>, SvError>,
    ) -> Result<Vec<SvRecordType>, SvError> {
        use InferredType::*;
        Ok(match self.inferred_type {
            InterChrStrandSwitch55
            | InterChrStrandSwitch33
            | InterChrNoSsWithLeftMateFirstInPartner
            | InterChrNoSsWithLeftMateSecondInPartner => vec![
                self.inter_chromosome_breakend(true, reference)?,
                self.inter_chromosome_breakend(false, reference)?,
            ],
            IntraChrRefOrderSwap => vec![
                self.ref_order_swap_breakend(true, reference)?,
                self.ref_order_swap_breakend(false, reference)?,
            ],
            IntraChrStrandSwitch55 | IntraChrStrandSwitch33 => {
                if self.complication.has_duplication_annotation() {
                    return Err(SvError::new(
                        "java.lang.ClassCastException",
                        "an inverted duplication is never annotated as a small duplication",
                    ));
                }
                let five_prime = self.strand_switch == StrandSwitch::ForwardToReverse;
                vec![
                    self.strand_switch_breakend(true, five_prime, reference)?,
                    self.strand_switch_breakend(false, five_prime, reference)?,
                ]
            }
            SimpleDel | DelDupContraction => vec![self.deletion(reference)?],
            Rpl => {
                if self.is_fat_insertion() {
                    vec![self.insertion(reference)?]
                } else {
                    let deletion = self.deletion(reference)?;
                    if (self.complication.inserted().len() as i32)
                        < STRUCTURAL_VARIANT_SIZE_LOWER_BOUND
                    {
                        vec![deletion]
                    } else {
                        vec![deletion, self.insertion(reference)?]
                    }
                }
            }
            SimpleIns => vec![self.insertion(reference)?],
            SmallDupExpansion => {
                let (unit, ..) = self.complication.duplication().ok_or_else(|| {
                    SvError::new("java.lang.ClassCastException", "not a small duplication")
                })?;
                if unit.size() < STRUCTURAL_VARIANT_SIZE_LOWER_BOUND {
                    vec![self.insertion(reference)?]
                } else {
                    vec![self.duplication_tandem(reference)?]
                }
            }
            SmallDupCpx => {
                let (unit, ..) = self.complication.duplication().ok_or_else(|| {
                    SvError::new("java.lang.ClassCastException", "not a small duplication")
                })?;
                if self.complication.is_dup_contraction() {
                    vec![self.deletion(reference)?]
                } else if unit.size() < STRUCTURAL_VARIANT_SIZE_LOWER_BOUND {
                    vec![self.insertion(reference)?]
                } else {
                    vec![self.duplication_tandem(reference)?]
                }
            }
        })
    }

    /// `SimpleSVType.Deletion`.
    fn deletion(
        &self,
        reference: &mut dyn FnMut(&Interval) -> Result<Vec<u8>, SvError>,
    ) -> Result<SvRecordType, SvError> {
        let duplication = self.complication.has_duplication_annotation();
        Ok(SvRecordType {
            kind: SvKind::Deletion,
            chromosome: self.left.contig.clone(),
            start: self.left.start,
            stop: self.right.end,
            id: format!(
                "{}_{}",
                if duplication {
                    "DEL-DUPLICATION-TANDEM-CONTRACTION"
                } else {
                    "DEL"
                },
                self.location()
            ),
            reference_allele: text(&reference(&self.left)?),
            alternate_allele: "<DEL>".to_string(),
            sv_len: -self.distance(),
            extra: if duplication {
                vec![("CONTRACTION", String::new())]
            } else {
                Vec::new()
            },
        })
    }

    /// `SimpleSVType.Insertion`, a point at the left breakpoint unless it is a fat insertion.
    fn insertion(
        &self,
        reference: &mut dyn FnMut(&Interval) -> Result<Vec<u8>, SvError>,
    ) -> Result<SvRecordType, SvError> {
        let fat = self.is_fat_insertion();
        let stop = if fat { self.right.end } else { self.left.start };
        let window = if fat {
            Interval::new(&self.left.contig, self.left.start, self.right.end)?
        } else {
            self.left.clone()
        };
        let reference_allele = text(&reference(&window)?);
        let sv_len = if self.complication.has_duplication_annotation() {
            self.length_for_dup_tandem_expansion()?
        } else {
            self.complication.inserted().len() as i32
        };
        let id = if fat {
            format!("INS_{}", self.location())
        } else {
            format!(
                "INS_{}",
                location_string(
                    &self.left.contig,
                    self.left.start,
                    &self.left.contig,
                    self.left.start
                )
            )
        };
        Ok(SvRecordType {
            kind: SvKind::Insertion,
            chromosome: self.left.contig.clone(),
            start: self.left.start,
            stop,
            id,
            reference_allele,
            alternate_allele: "<INS>".to_string(),
            sv_len,
            extra: Vec::new(),
        })
    }

    /// `SimpleSVType.DuplicationTandem`, named by the repeat unit.
    fn duplication_tandem(
        &self,
        reference: &mut dyn FnMut(&Interval) -> Result<Vec<u8>, SvError>,
    ) -> Result<SvRecordType, SvError> {
        let (unit, ..) = self
            .complication
            .duplication()
            .expect("a tandem duplication carries its repeat unit");
        let id = format!(
            "INS-DUPLICATION-TANDEM-EXPANSION_{}",
            location_string(&unit.contig, unit.start, &unit.contig, unit.end)
        );
        Ok(SvRecordType {
            kind: SvKind::DuplicationTandem,
            chromosome: self.left.contig.clone(),
            start: self.left.start,
            stop: self.left.start,
            id,
            reference_allele: text(&reference(&self.left)?),
            alternate_allele: "<DUP>".to_string(),
            sv_len: self.length_for_dup_tandem_expansion()?,
            extra: vec![("EXPANSION", String::new())],
        })
    }

    /// `BreakEndVariantType.getIDString`.
    fn breakend_id(&self, upstream: bool) -> String {
        let kind = if self.strand_switch == StrandSwitch::NoSwitch
            || self.left.contig != self.right.contig
        {
            ""
        } else if self.strand_switch == StrandSwitch::ForwardToReverse {
            "INV55_"
        } else {
            "INV33_"
        };
        format!(
            "BND_{kind}{}_{}",
            self.location(),
            if upstream { "1" } else { "2" }
        )
    }

    fn breakend_base(
        &self,
        upstream: bool,
        reference: &mut dyn FnMut(&Interval) -> Result<Vec<u8>, SvError>,
    ) -> Result<String, SvError> {
        Ok(text(&reference(if upstream {
            &self.left
        } else {
            &self.right
        })?))
    }

    fn breakend(&self, upstream: bool, base: String, alternate: String) -> SvRecordType {
        let (chromosome, start) = if upstream {
            (self.left.contig.clone(), self.left.start)
        } else {
            (self.right.contig.clone(), self.right.end)
        };
        SvRecordType {
            kind: SvKind::BreakEnd,
            chromosome,
            start,
            stop: -1,
            id: self.breakend_id(upstream),
            reference_allele: base,
            alternate_allele: alternate,
            sv_len: -1,
            extra: Vec::new(),
        }
    }

    /// `IntraChromosomalStrandSwitch55BreakEnd` and `...33BreakEnd`.
    fn strand_switch_breakend(
        &self,
        upstream: bool,
        five_prime: bool,
        reference: &mut dyn FnMut(&Interval) -> Result<Vec<u8>, SvError>,
    ) -> Result<SvRecordType, SvError> {
        let base = self.breakend_base(upstream, reference)?;
        let inserted = self.complication.inserted();
        let inserted = if upstream {
            inserted.to_string()
        } else {
            text(&reverse_complement(inserted.as_bytes()))
        };
        let mate = if upstream { &self.right } else { &self.left };
        let alternate = if five_prime {
            format!("{base}{inserted}]{}:{}]", mate.contig, mate.end)
        } else {
            format!("[{}:{}[{inserted}{base}", mate.contig, mate.end)
        };
        let mut record = self.breakend(upstream, base, alternate);
        record.extra = vec![(if five_prime { "INV55" } else { "INV33" }, String::new())];
        Ok(record)
    }

    /// `IntraChromosomeRefOrderSwap`.
    fn ref_order_swap_breakend(
        &self,
        upstream: bool,
        reference: &mut dyn FnMut(&Interval) -> Result<Vec<u8>, SvError>,
    ) -> Result<SvRecordType, SvError> {
        let base = self.breakend_base(upstream, reference)?;
        let inserted = self.complication.inserted();
        let alternate = if upstream {
            format!("]{}:{}]{inserted}{base}", self.right.contig, self.right.end)
        } else {
            format!("{base}{inserted}[{}:{}[", self.left.contig, self.left.end)
        };
        Ok(self.breakend(upstream, base, alternate))
    }

    /// `InterChromosomeBreakend`.
    fn inter_chromosome_breakend(
        &self,
        upstream: bool,
        reference: &mut dyn FnMut(&Interval) -> Result<Vec<u8>, SvError>,
    ) -> Result<SvRecordType, SvError> {
        let base = self.breakend_base(upstream, reference)?;
        let inserted = self.complication.inserted();
        let inserted = if inserted.is_empty() || self.strand_switch == StrandSwitch::NoSwitch {
            inserted.to_string()
        } else if upstream == (self.strand_switch == StrandSwitch::ForwardToReverse) {
            inserted.to_string()
        } else {
            text(&reverse_complement(inserted.as_bytes()))
        };
        let mate = if upstream { &self.right } else { &self.left };
        let first_in_partner =
            self.inferred_type == InferredType::InterChrNoSsWithLeftMateFirstInPartner;
        let alternate = match self.strand_switch {
            StrandSwitch::NoSwitch => {
                if upstream == first_in_partner {
                    format!("{base}{inserted}[{}:{}[", mate.contig, mate.end)
                } else {
                    format!("]{}:{}]{inserted}{base}", mate.contig, mate.start)
                }
            }
            StrandSwitch::ForwardToReverse => {
                format!("{base}{inserted}]{}:{}]", mate.contig, mate.end)
            }
            StrandSwitch::ReverseToForward => {
                format!("[{}:{}[{inserted}{base}", mate.contig, mate.end)
            }
        };
        Ok(self.breakend(upstream, base, alternate))
    }
}

/// `SvType.makeLocationString`: the second contig is written only when it differs.
fn location_string(first: &str, one: i32, second: &str, two: i32) -> String {
    if first == second {
        format!("{first}_{one}_{two}")
    } else {
        format!("{first}_{one}_{second}_{two}")
    }
}

/// The subclass of `SvType` a record was made by.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SvKind {
    Deletion,
    Insertion,
    DuplicationTandem,
    BreakEnd,
}

/// `SvType`: what `getBasicInformation` builds a record from.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SvRecordType {
    pub kind: SvKind,
    pub chromosome: String,
    pub start: i32,
    /// `NO_APPLICABLE_END` (-1) for a breakend.
    pub stop: i32,
    pub id: String,
    pub reference_allele: String,
    pub alternate_allele: String,
    pub sv_len: i32,
    pub extra: Vec<(&'static str, String)>,
}

impl SvRecordType {
    /// `toString`, which is the `SVTYPE` value.
    pub fn sv_type(&self) -> &'static str {
        match self.kind {
            SvKind::Deletion => "DEL",
            SvKind::Insertion => "INS",
            SvKind::DuplicationTandem => "DUP",
            SvKind::BreakEnd => "BND",
        }
    }
}

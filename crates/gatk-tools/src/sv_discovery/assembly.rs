//! `SvDiscoverFromLocalAssemblyContigAlignmentsSpark`: each contig's best alignment configuration,
//! its classification, the simple and breakend calls, the complex calls, and the simple calls
//! derived from those.
//!
//! The classes this is ported from are `AssemblyContigAlignmentsRDDProcessor`, `AlignedContig`'s
//! configuration picking, `AssemblyContigWithFineTunedAlignments`,
//! `SimpleNovelAdjacencyInterpreter`, `NovelAdjacencyAndAltHaplotype.toSimpleOrBNDTypes`,
//! `SimpleSVType`, `BreakEndVariantType`, `CpxVariantInterpreter`,
//! `CpxVariantInducingAssemblyContig`, `CpxVariantCanonicalRepresentation`, `CpxVariantType`,
//! `ContigAlignmentsModifier.clipAlignmentInterval`, `SegmentedCpxVariantSimpleVariantExtractor` and
//! `AnnotatedVariantProducer.filterMergedVariantList`.
//!
//! Three things in them are not what they look like:
//!
//! * **the configurations come out of a `HashSet` of alignments**, whose order rests on a cigar
//!   operator's identity hash and differs between runs. It cannot be visible: configurations that
//!   tie make the contig ambiguous, and an ambiguous contig reaches nothing but its reason;
//! * **an inverted duplication cannot be typed from here**: two alignments overlapping across a
//!   strand switch are an incomplete picture and never reach `toSimpleOrBNDTypes`;
//! * **`getInsFromOneEnd`'s rear insertion of a segment that opens the arrangement carries the
//!   front anchor**, which is the base before the segment, so its REF is not the base at its POS.

use super::*;
use gatk_engine::cigar_builder::CigarBuilder;
use gatk_engine::java_hash::hash_map_order;

/// `AlignedContig.ALIGNMENT_MQ_THRESHOLD`.
const ALIGNMENT_MQ_THRESHOLD: i32 = 20;
/// `ALIGNMENT_MQ_THRESHOLD_FOR_SPEED_BOOST`.
const SPEED_BOOST_MQ_THRESHOLD: i32 = 10;
/// `ALIGNMENT_LOW_READ_UNIQUENESS_THRESHOLD`.
const LOW_READ_UNIQUENESS: i32 = 10;
/// `SECONDARY_CONFIGURATION_MQ_FILTER_THRESHOLD`.
const SECONDARY_CONFIGURATION_MQ: i32 = 0;
/// `COVERAGE_MQ_NORMALIZATION_CONST`.
const MQ_NORMALIZATION: f64 = 60.0;
/// `SimpleChimera.SPLIT_PAIR_MIN_ALIGNMENT_LENGTH`.
const SPLIT_PAIR_MIN_LENGTH: i32 = 30;
/// `SimpleChimera.SPLIT_PAIR_MIN_ALIGNMENT_MQ`.
const SPLIT_PAIR_MIN_MQ: i32 = 20;
/// `CpxVariantInterpreter.MIN_READ_SPAN_AFTER_DEOVERLAP`.
const MIN_READ_SPAN_AFTER_DEOVERLAP: i32 = 2;
/// `STRUCTURAL_VARIANT_SIZE_LOWER_BOUND`.
const SIZE_LOWER_BOUND: i32 = 50;
/// `SegmentedCpxVariantSimpleVariantExtractor.EVENT_SIZE_THRESHOLD`.
const EVENT_SIZE_THRESHOLD: i32 = SIZE_LOWER_BOUND - 1;
/// `AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME`.
const NO_NON_CANONICAL: &str = "NONE";

fn size_on_read(alignment: &Alignment) -> i32 {
    alignment.end_on_contig - alignment.start_on_contig + 1
}

fn contains_on_read(one: &Alignment, two: &Alignment) -> bool {
    one.start_on_contig <= two.start_on_contig && one.end_on_contig >= two.end_on_contig
}

/// `containsGapOfEqualOrLargerSize(GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY)`.
fn has_large_gap(alignment: &Alignment) -> bool {
    alignment
        .cigar
        .iter()
        .any(|element| is_indel(element.op) && element.length as i32 >= GAP_SENSITIVITY)
}

/// `CigarUtils.countUnclippedReadBases`.
fn unclipped_read_bases(cigar: &[CigarElement]) -> i32 {
    cigar
        .iter()
        .filter(|element| is_clipping(element.op) || consumes_read(element.op))
        .map(|element| element.length as i32)
        .sum()
}

/// `getAlignmentIntervalComparator`, in a stable sort.
fn sort_alignments(alignments: &mut [Alignment]) {
    alignments.sort_by(|one, two| {
        one.start_on_contig
            .cmp(&two.start_on_contig)
            .then_with(|| one.span.contig.cmp(&two.span.contig))
            .then(one.span.start.cmp(&two.span.start))
    });
}

/// `GoodAndBadMappings`.
#[derive(Clone)]
struct Mappings {
    good: Vec<Alignment>,
    bad: Vec<Alignment>,
    non_canonical: Option<Alignment>,
}

/// `AssemblyContigWithFineTunedAlignments`.
#[derive(Clone)]
pub struct FineTuned {
    contig: Contig,
    insertion_mappings: Vec<String>,
    ambiguous: bool,
    non_canonical: String,
}

impl FineTuned {
    fn head(&self) -> &Alignment {
        &self.contig.alignments[0]
    }

    fn tail(&self) -> &Alignment {
        &self.contig.alignments[self.contig.alignments.len() - 1]
    }

    /// `hasIncompletePictureFromMultipleAlignments` and its two-alignment form.
    fn incomplete(&self) -> bool {
        let alignments = &self.contig.alignments;
        if alignments.len() < 2 {
            return false;
        }
        let (head, tail) = (self.head(), self.tail());
        if alignments.len() == 2 {
            return incomplete_from_two(head, tail);
        }
        if head.span.contig != tail.span.contig || head.forward != tail.forward {
            return true;
        }
        if head.span.contains(&tail.span) || tail.span.contains(&head.span) {
            return true;
        }
        let valid = Interval {
            contig: head.span.contig.clone(),
            start: head.span.start.min(tail.span.start),
            end: head.span.end.max(tail.span.end),
        };
        if alignments[1..alignments.len() - 1]
            .iter()
            .any(|middle| middle.span.overlaps(&valid) && !valid.contains(&middle.span))
        {
            return true;
        }
        if head.forward {
            head.span.start >= tail.span.start
        } else {
            head.span.end <= tail.span.end
        }
    }

    /// `getAlignmentSignatureBasicType`, as `AssemblyContigsClassifiedByAlignmentSignatures` sorts it.
    fn class(&self) -> Signature {
        if self.ambiguous || self.incomplete() {
            Signature::Unknown
        } else if self.contig.alignments.len() < 2 {
            if self.contig.alignments.len() == 1 && self.insertion_mappings.is_empty() {
                Signature::Normal
            } else {
                Signature::Unknown
            }
        } else if self.contig.alignments.len() == 2 {
            Signature::Simple
        } else {
            Signature::Complex
        }
    }

    /// `getReasonForAlignmentClassificationFailure`.
    fn reason(&self) -> Result<usize, Failed> {
        if self.ambiguous {
            Ok(0)
        } else if self.incomplete() {
            Ok(1)
        } else if !(self.insertion_mappings.is_empty() || self.contig.alignments.len() > 1) {
            Ok(2)
        } else {
            Err(format!(
                "operating on contig without a suspicious alignment signature, contig: {}",
                self.text()
            ))
        }
    }

    /// `toString`.
    fn text(&self) -> String {
        format!(
            "AssemblyContigWithFineTunedAlignments{{sourceTig=({}, [{}]), insertionMappings=[{}], \
             hasEquallyGoodAlnConfigurations={}, saTAGForGoodMappingToNonCanonicalChromosome='{}'}}",
            self.contig.name,
            self.contig
                .alignments
                .iter()
                .map(Alignment::packed)
                .collect::<Vec<_>>()
                .join(", "),
            self.insertion_mappings.join(", "),
            self.ambiguous,
            self.non_canonical
        )
    }
}

/// `AlignmentSignatureBasicType`.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Signature {
    Normal,
    Simple,
    Complex,
    Unknown,
}

/// `AlignedContig.hasGoodMQ`.
fn has_good_mq(contig: &Contig) -> bool {
    let alignments = &contig.alignments;
    if alignments.len() < 2 {
        return alignments
            .first()
            .is_some_and(|only| only.map_quality > ALIGNMENT_MQ_THRESHOLD);
    }
    let mut not_bad = 0;
    for alignment in alignments {
        if alignment.map_quality > ALIGNMENT_MQ_THRESHOLD {
            if has_large_gap(alignment) {
                return true;
            }
            not_bad += 1;
        }
    }
    not_bad > 1
}

/// `computeScoreOfConfiguration`.
fn configuration_score(
    configuration: &[Alignment],
    canonical: &[String],
    maximum_canonical_score: i32,
) -> f64 {
    let mut explained = 0.0;
    for alignment in configuration {
        let quality = f64::from(alignment.map_quality) / MQ_NORMALIZATION;
        let weight = if canonical.contains(&alignment.span.contig) {
            quality
        } else {
            quality.max(if alignment.score > maximum_canonical_score {
                1.0
            } else {
                0.0
            })
        };
        explained += weight * f64::from(size_on_read(alignment));
    }
    let mut redundancy = 0;
    for (index, one) in configuration.iter().enumerate() {
        for two in &configuration[index + 1..] {
            redundancy += overlap_on_contig(one, two);
        }
    }
    explained - f64::from(redundancy)
}

/// `Math.ulp`, of a value that is never infinite here.
fn ulp(value: f64) -> f64 {
    let magnitude = value.abs();
    if magnitude == 0.0 {
        return f64::from_bits(1);
    }
    f64::from_bits(magnitude.to_bits() + 1) - magnitude
}

fn maximum_canonical_score(alignments: &[Alignment], canonical: &[String]) -> i32 {
    alignments
        .iter()
        .filter(|alignment| canonical.contains(&alignment.span.contig))
        .map(|alignment| alignment.score)
        .max()
        .unwrap_or(0)
}

/// `pickBestConfigurations`.
fn pick_best_configurations(
    contig: &Contig,
    canonical: &[String],
) -> Result<Vec<Mappings>, Failed> {
    let alignments = &contig.alignments;
    if alignments.len() == 1 {
        return Ok(vec![Mappings {
            good: alignments.clone(),
            bad: Vec::new(),
            non_canonical: None,
        }]);
    }
    let maximum = maximum_canonical_score(alignments, canonical);
    // `heuristicSpeedUpWhenFacingManyMappings`.
    let (mut good, bad): (Vec<Alignment>, Vec<Alignment>) = if alignments.len() > 10 {
        alignments.iter().cloned().partition(|alignment| {
            (!canonical.contains(&alignment.span.contig) && alignment.score > maximum)
                || alignment.map_quality > SPEED_BOOST_MQ_THRESHOLD
        })
    } else {
        (alignments.clone(), Vec::new())
    };
    let maximum = maximum_canonical_score(&good, canonical);

    // `getBetterNonCanonicalMapping`.
    let (canonical_mappings, non_canonical_mappings): (Vec<Alignment>, Vec<Alignment>) = good
        .iter()
        .cloned()
        .partition(|alignment| canonical.contains(&alignment.span.contig));
    let mut non_canonical = None;
    if !canonical_mappings.is_empty()
        && non_canonical_mappings.len() == 1
        && (canonical_mappings.len() > 1 || has_large_gap(&canonical_mappings[0]))
        && configuration_score(&canonical_mappings, canonical, maximum)
            <= configuration_score(&non_canonical_mappings, canonical, maximum)
    {
        non_canonical = Some(non_canonical_mappings[0].clone());
    }
    if let Some(taken) = &non_canonical {
        if let Some(at) = good.iter().position(|alignment| alignment == taken) {
            good.remove(at);
        }
    }

    // `generateScoreAndPickConfigurations`.
    let mut distinct: Vec<Alignment> = Vec::with_capacity(good.len());
    for alignment in &good {
        if !distinct.contains(alignment) {
            distinct.push(alignment.clone());
        }
    }
    if distinct.len() > 30 {
        return Err(format!(
            "Too many elements to create power set: {} > 30",
            distinct.len()
        ));
    }
    let subsets = 1u64 << distinct.len();
    let configuration_of = |mask: u64| {
        let mut configuration: Vec<Alignment> = distinct
            .iter()
            .enumerate()
            .filter(|(bit, _)| mask >> bit & 1 == 1)
            .map(|(_, alignment)| alignment.clone())
            .collect();
        sort_alignments(&mut configuration);
        configuration
    };
    let mut scores = Vec::with_capacity(subsets as usize);
    for mask in 0..subsets {
        scores.push(configuration_score(
            &configuration_of(mask),
            canonical,
            maximum,
        ));
    }
    let best = scores.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut picked = Vec::new();
    for (mask, score) in scores.iter().enumerate() {
        if *score >= best || best - score <= ulp(*score).max(0.0) {
            let configuration = configuration_of(mask as u64);
            let mut left: Vec<Alignment> = good
                .iter()
                .filter(|alignment| !configuration.contains(alignment))
                .cloned()
                .collect();
            left.extend(bad.iter().cloned());
            picked.push(Mappings {
                good: configuration,
                bad: left,
                non_canonical: non_canonical.clone(),
            });
        }
    }
    Ok(picked)
}

/// `filterSecondaryConfigurationsByMappingQualityThreshold`.
fn filter_secondary_configurations(configurations: Vec<Mappings>) -> Vec<Mappings> {
    if configurations.len() == 1 {
        return configurations;
    }
    let above: Vec<&Mappings> = configurations
        .iter()
        .filter(|mappings| {
            mappings
                .good
                .iter()
                .map(|alignment| alignment.map_quality)
                .min()
                .unwrap_or(SECONDARY_CONFIGURATION_MQ)
                > SECONDARY_CONFIGURATION_MQ
        })
        .collect();
    if above.len() == 1 {
        let only = above[0].clone();
        vec![only]
    } else {
        configurations
    }
}

/// `splitGapsAndDropAlignmentContainedByOtherOnRead`.
fn split_gaps_and_drop_contained(configuration: Mappings) -> Result<Mappings, Failed> {
    let mut split = Vec::with_capacity(configuration.good.len());
    for alignment in configuration.good {
        if has_large_gap(&alignment) {
            let unclipped = unclipped_read_bases(&alignment.cigar);
            split.extend(split_gaps(alignment, unclipped)?);
        } else {
            split.push(alignment);
        }
    }
    let mut bad = configuration.bad;
    for (index, one) in split.iter().enumerate() {
        for two in &split[index + 1..] {
            if contains_on_read(one, two) {
                bad.push(two.clone());
            } else if contains_on_read(two, one) {
                bad.push(one.clone());
            }
        }
    }
    split.retain(|alignment| !bad.contains(alignment));
    sort_alignments(&mut split);
    Ok(Mappings {
        good: split,
        bad,
        non_canonical: configuration.non_canonical,
    })
}

/// `removeNonUniqueMappings`, the map `getMaxOverlapPairs` returns being read back by position:
/// two equal alignments cannot both be here, having been a set.
fn remove_non_unique_mappings(mappings: Mappings) -> Mappings {
    if mappings.good.len() <= 2 {
        return mappings;
    }
    let mut bad = mappings.bad;
    let mut selected = Vec::with_capacity(mappings.good.len());
    for alignment in mappings.good {
        if alignment.map_quality >= ALIGNMENT_MQ_THRESHOLD {
            selected.push(alignment);
        } else {
            bad.push(alignment);
        }
    }
    let count = selected.len();
    let mut front = vec![-1; count];
    let mut rear = vec![-1; count];
    for index in 0..count.saturating_sub(1) {
        let mut maximum = -1;
        let mut target = usize::MAX;
        for other in index + 1..count {
            let overlap = overlap_on_contig(&selected[index], &selected[other]);
            if overlap > maximum {
                maximum = overlap;
                target = other;
            } else {
                break;
            }
        }
        if maximum > 0 {
            rear[index] = maximum;
            if front[target] < maximum {
                front[target] = maximum;
            }
        }
    }
    let mut kept = Vec::with_capacity(count);
    for (index, alignment) in selected.into_iter().enumerate() {
        let unique = size_on_read(&alignment) - front[index].max(0) - rear[index].max(0);
        if unique < LOW_READ_UNIQUENESS {
            bad.push(alignment);
        } else {
            kept.push(alignment);
        }
    }
    Mappings {
        good: kept,
        bad,
        non_canonical: mappings.non_canonical,
    }
}

/// `alignmentShouldNotBeStitchedTogether`.
fn should_not_be_stitched(alignments: &[Alignment]) -> Result<bool, Failed> {
    Ok(alignments.len() != 2 || !stitchable(&alignments[0], &alignments[1])?)
}

/// `createContigGivenClassifiedAlignments`.
fn fine_tuned(contig: &Contig, mappings: Mappings, ambiguous: bool) -> FineTuned {
    let mut alignments = mappings.good;
    sort_alignments(&mut alignments);
    FineTuned {
        contig: Contig {
            name: contig.name.clone(),
            sequence: contig.sequence.clone(),
            alignments,
        },
        insertion_mappings: mappings.bad.iter().map(Alignment::packed).collect(),
        ambiguous,
        non_canonical: mappings
            .non_canonical
            .as_ref()
            .map_or_else(|| NO_NON_CANONICAL.to_string(), Alignment::sa_text),
    }
}

/// `reconstructContigFromBestConfiguration`.
fn reconstruct(contig: &Contig, canonical: &[String]) -> Result<Vec<FineTuned>, Failed> {
    let mut configurations =
        filter_secondary_configurations(pick_best_configurations(contig, canonical)?);
    if configurations.len() > 1 {
        let mut contigs = Vec::with_capacity(configurations.len());
        for configuration in configurations {
            let mappings =
                remove_non_unique_mappings(split_gaps_and_drop_contained(configuration)?);
            if should_not_be_stitched(&mappings.good)? {
                contigs.push(fine_tuned(contig, mappings, true));
            }
        }
        // `getConfigurationComparator`.
        contigs.sort_by_key(|tig| {
            (
                tig.contig.alignments.len(),
                tig.contig
                    .alignments
                    .iter()
                    .map(|alignment| alignment.mismatches)
                    .sum::<i32>(),
            )
        });
        return Ok(contigs);
    }
    let only = configurations.remove(0);
    let mappings = remove_non_unique_mappings(split_gaps_and_drop_contained(only)?);
    Ok(if should_not_be_stitched(&mappings.good)? {
        vec![fine_tuned(contig, mappings, false)]
    } else {
        Vec::new()
    })
}

/// `AssemblyContigsClassifiedByAlignmentSignatures`.
pub struct Classified {
    unknown: Vec<FineTuned>,
    simple: Vec<FineTuned>,
    complex: Vec<FineTuned>,
}

/// `SvDiscoverFromLocalAssemblyContigAlignmentsSpark.preprocess`.
pub fn preprocess(
    records: &[BamRecord],
    header: &SamHeader,
    canonical: &[String],
) -> Result<Classified, Failed> {
    let by_name = append_only_group(
        records
            .iter()
            .filter(|record| record.flags & 0x100 == 0)
            .map(|record| (record.read_name.clone(), record))
            .collect(),
        |name: &String| string_hash_code(name),
    );
    let mut parsed = Vec::with_capacity(by_name.len());
    for (_, reads) in &by_name {
        parsed.push(parse_contig(reads, header, false)?);
    }
    let mut classified = Classified {
        unknown: Vec::new(),
        simple: Vec::new(),
        complex: Vec::new(),
    };
    for contig in parsed.iter().filter(|contig| has_good_mq(contig)) {
        for tuned in reconstruct(contig, canonical)? {
            match tuned.class() {
                Signature::Unknown => classified.unknown.push(tuned),
                Signature::Simple => classified.simple.push(tuned),
                Signature::Complex => classified.complex.push(tuned),
                Signature::Normal => {}
            }
        }
    }
    Ok(classified)
}

impl Classified {
    /// `writeSAMfilesForUnknown`: the raw records of each reason in `AMBIGUOUS`, `INCOMPLETE` and
    /// `UNINFORMATIVE` order, sorted by name as a stable sort leaves them.
    pub fn unknown_reads<'a>(
        &self,
        records: &'a [BamRecord],
    ) -> Result<[Vec<&'a BamRecord>; 3], Failed> {
        let mut reasons: Vec<(&str, usize)> = Vec::with_capacity(self.unknown.len());
        for tig in &self.unknown {
            reasons.push((tig.contig.name.as_str(), tig.reason()?));
        }
        let mut sorted: Vec<&BamRecord> = records
            .iter()
            .filter(|record| reasons.iter().any(|(name, _)| *name == record.read_name))
            .collect();
        sorted.sort_by(|one, two| one.read_name.cmp(&two.read_name));
        let mut sets: [Vec<&BamRecord>; 3] = [Vec::new(), Vec::new(), Vec::new()];
        for record in sorted {
            // `collectAsMap`, where the last of a name wins.
            let reason = reasons
                .iter()
                .rev()
                .find(|(name, _)| *name == record.read_name)
                .map_or(0, |(_, reason)| *reason);
            sets[reason].push(record);
        }
        Ok(sets)
    }
}

/// `SvType`, as far as `toVariantContexts` tells them apart.
struct TypedCall {
    variant: VariantContext,
    breakend: bool,
    deletion: bool,
}

fn allele(bases: &[u8], is_reference: bool) -> Result<Allele, Failed> {
    Allele::create(bases, is_reference).map_err(|error| format!("{error:?}"))
}

fn symbolic(symbol: &str) -> Result<Allele, Failed> {
    Allele::from_str(symbol, false).map_err(|error| format!("{error:?}"))
}

/// `SvType.getBasicInformation`, of a type with an end and a length.
fn simple_call(
    contig: &str,
    position: i32,
    stop: i32,
    id: String,
    alleles: Vec<Allele>,
    svtype: &str,
    length: i32,
) -> VariantContext {
    let mut variant = VariantContext::new(contig, i64::from(position), alleles);
    variant.stop = i64::from(stop);
    variant.id = id;
    set(&mut variant, "END", Value::Int(i64::from(stop)));
    set(&mut variant, "SVTYPE", Value::Str(svtype.to_string()));
    set(&mut variant, "SVLEN", Value::Int(i64::from(length)));
    variant
}

/// `SimpleSVType.Deletion`.
fn deletion_call(
    adjacency: &NovelAdjacency,
    reference: &mut dyn Reference,
) -> Result<TypedCall, Failed> {
    let contig = &adjacency.left.contig;
    let position = adjacency.left.start;
    let duplication = adjacency.complication.has_duplication();
    let id = format!(
        "{}_{}",
        if duplication {
            "DEL-DUPLICATION-TANDEM-CONTRACTION"
        } else {
            "DEL"
        },
        adjacency.location()
    );
    let bases = reference.bases(contig, position, position)?;
    let mut variant = simple_call(
        contig,
        position,
        adjacency.right.end,
        id,
        vec![allele(&bases, true)?, symbolic("<DEL>")?],
        "DEL",
        -(adjacency.right.end - adjacency.left.start),
    );
    if duplication {
        set(&mut variant, "CONTRACTION", Value::Bool(true));
    }
    Ok(TypedCall {
        variant,
        breakend: false,
        deletion: true,
    })
}

/// `NovelAdjacencyAndAltHaplotype.isCandidateForFatInsertion`.
fn fat_insertion(adjacency: &NovelAdjacency) -> bool {
    adjacency.kind == ChimeraType::Replacement
        && adjacency.right.end - adjacency.left.start < SIZE_LOWER_BOUND
}

/// `SimpleSVType.Insertion`.
fn insertion_call(
    adjacency: &NovelAdjacency,
    reference: &mut dyn Reference,
) -> Result<TypedCall, Failed> {
    let contig = &adjacency.left.contig;
    let position = adjacency.left.start;
    let length = if adjacency.complication.has_duplication() {
        adjacency.expansion_length()?
    } else {
        adjacency.complication.inserted.len() as i32
    };
    let (stop, id, bases);
    if fat_insertion(adjacency) {
        stop = adjacency.right.end;
        id = format!("INS_{}", adjacency.location());
        bases = reference.bases(contig, position, stop)?;
    } else {
        stop = position;
        id = format!(
            "INS_{}",
            location_string(contig, position, contig, position)
        );
        bases = reference.bases(contig, position, position)?;
    }
    let variant = simple_call(
        contig,
        position,
        stop,
        id,
        vec![allele(&bases, true)?, symbolic("<INS>")?],
        "INS",
        length,
    );
    Ok(TypedCall {
        variant,
        breakend: false,
        deletion: false,
    })
}

/// `SimpleSVType.DuplicationTandem`.
fn tandem_call(
    adjacency: &NovelAdjacency,
    reference: &mut dyn Reference,
) -> Result<TypedCall, Failed> {
    let contig = &adjacency.left.contig;
    let position = adjacency.left.start;
    let unit = adjacency
        .complication
        .unit
        .as_ref()
        .expect("a duplication has its unit");
    let id = format!(
        "INS-DUPLICATION-TANDEM-EXPANSION_{}",
        location_string(&unit.contig, unit.start, &unit.contig, unit.end)
    );
    let bases = reference.bases(contig, position, position)?;
    let mut variant = simple_call(
        contig,
        position,
        position,
        id,
        vec![allele(&bases, true)?, symbolic("<DUP>")?],
        "DUP",
        adjacency.expansion_length()?,
    );
    set(&mut variant, "EXPANSION", Value::Bool(true));
    Ok(TypedCall {
        variant,
        breakend: false,
        deletion: false,
    })
}

/// The mates of `BreakEndVariantType`'s subclasses.
fn breakend_calls(
    adjacency: &NovelAdjacency,
    reference: &mut dyn Reference,
) -> Result<Vec<TypedCall>, Failed> {
    let (left, right) = (&adjacency.left, &adjacency.right);
    let same_contig = left.contig == right.contig;
    let tag = if adjacency.switch == StrandSwitch::None || !same_contig {
        ""
    } else if adjacency.switch == StrandSwitch::ForwardToReverse {
        "INV55_"
    } else {
        "INV33_"
    };
    let inserted = &adjacency.complication.inserted;
    let reverse_inserted = || text_of(reverse_complement(inserted.as_bytes()));
    let mut calls = Vec::with_capacity(2);
    for upstream in [true, false] {
        let (own, mate) = if upstream {
            (left, right)
        } else {
            (right, left)
        };
        let bases = reference.bases(&own.contig, own.start, own.end)?;
        let base = String::from_utf8_lossy(&bases).into_owned();
        let at_end = format!("{}:{}", mate.contig, mate.end);
        let mut flag = None;
        let alternate = match adjacency.kind {
            ChimeraType::InterChromosomeStrandSwitch
            | ChimeraType::InterChromosomeLeftMateFirst
            | ChimeraType::InterChromosomeLeftMateSecond => {
                // `InterChromosomeBreakend.extractInsertedSequence`.
                let sequence = if inserted.is_empty()
                    || adjacency.switch == StrandSwitch::None
                    || upstream == (adjacency.switch == StrandSwitch::ForwardToReverse)
                {
                    inserted.clone()
                } else {
                    reverse_inserted()
                };
                let left_first = adjacency.kind == ChimeraType::InterChromosomeLeftMateFirst;
                match adjacency.switch {
                    StrandSwitch::None if upstream == left_first => {
                        format!("{base}{sequence}[{at_end}[")
                    }
                    StrandSwitch::None => {
                        format!("]{}:{}]{sequence}{base}", mate.contig, mate.start)
                    }
                    StrandSwitch::ForwardToReverse => format!("{base}{sequence}]{at_end}]"),
                    StrandSwitch::ReverseToForward => format!("[{at_end}[{sequence}{base}"),
                }
            }
            ChimeraType::ReferenceOrderSwap => {
                if upstream {
                    format!("]{at_end}]{inserted}{base}")
                } else {
                    format!("{base}{inserted}[{at_end}[")
                }
            }
            _ => {
                let sequence = if upstream {
                    inserted.clone()
                } else {
                    reverse_inserted()
                };
                if adjacency.switch == StrandSwitch::ForwardToReverse {
                    flag = Some("INV55");
                    format!("{base}{sequence}]{at_end}]")
                } else {
                    flag = Some("INV33");
                    format!("[{at_end}[{sequence}{base}")
                }
            }
        };
        let position = if upstream { own.start } else { own.end };
        let mut variant = VariantContext::new(
            &own.contig,
            i64::from(position),
            vec![allele(&bases, true)?, symbolic(&alternate)?],
        );
        variant.stop = i64::from(position);
        variant.id = format!(
            "BND_{tag}{}_{}",
            adjacency.location(),
            if upstream { "1" } else { "2" }
        );
        set(&mut variant, "SVTYPE", Value::Str("BND".to_string()));
        if let Some(flag) = flag {
            set(&mut variant, flag, Value::Bool(true));
        }
        calls.push(TypedCall {
            variant,
            breakend: true,
            deletion: false,
        });
    }
    Ok(calls)
}

/// `NovelAdjacencyAndAltHaplotype.toSimpleOrBNDTypes`.
fn typed_calls(
    adjacency: &NovelAdjacency,
    reference: &mut dyn Reference,
) -> Result<Vec<TypedCall>, Failed> {
    let unit_is_short = || {
        adjacency
            .complication
            .unit
            .as_ref()
            .map_or(0, Interval::size)
            < SIZE_LOWER_BOUND
    };
    Ok(match adjacency.kind {
        ChimeraType::InterChromosomeStrandSwitch
        | ChimeraType::InterChromosomeLeftMateFirst
        | ChimeraType::InterChromosomeLeftMateSecond
        | ChimeraType::ReferenceOrderSwap
        | ChimeraType::StrandSwitch => breakend_calls(adjacency, reference)?,
        ChimeraType::SimpleDeletion | ChimeraType::DuplicationContraction => {
            vec![deletion_call(adjacency, reference)?]
        }
        ChimeraType::Replacement => {
            if fat_insertion(adjacency) {
                vec![insertion_call(adjacency, reference)?]
            } else if (adjacency.complication.inserted.len() as i32) < SIZE_LOWER_BOUND {
                vec![deletion_call(adjacency, reference)?]
            } else {
                vec![
                    deletion_call(adjacency, reference)?,
                    insertion_call(adjacency, reference)?,
                ]
            }
        }
        ChimeraType::SimpleInsertion => vec![insertion_call(adjacency, reference)?],
        ChimeraType::DuplicationExpansion => {
            if unit_is_short() {
                vec![insertion_call(adjacency, reference)?]
            } else {
                vec![tandem_call(adjacency, reference)?]
            }
        }
        ChimeraType::DuplicationComplex => {
            if adjacency.complication.contraction() {
                vec![deletion_call(adjacency, reference)?]
            } else if unit_is_short() {
                vec![insertion_call(adjacency, reference)?]
            } else {
                vec![tandem_call(adjacency, reference)?]
            }
        }
    })
}

/// What `produceAnnotatedVcFromAssemblyEvidence` makes of one type.
fn annotated_call(
    call: &TypedCall,
    adjacency: &NovelAdjacency,
    evidence: &[Chimera],
    header: &SamHeader,
    cnv_calls: Option<&[CnvCall]>,
) -> Result<VariantContext, Failed> {
    let mut variant = call.variant.clone();
    let stop = if call.breakend {
        NO_VALUE
    } else {
        variant.stop as i32
    };
    annotate_from_evidence(&mut variant, adjacency, evidence, header, cnv_calls, stop)?;
    if call.breakend {
        remove(&mut variant, "SEQ_ALT_HAPLOTYPE");
    }
    Ok(variant)
}

fn remove(variant: &mut VariantContext, key: &str) {
    variant.attributes.retain(|(name, _)| name != key);
}

/// The novel adjacencies of the simple chimeras, each with its evidence and its types.
pub struct SimpleCalls(Vec<(NovelAdjacency, Vec<Chimera>, Vec<TypedCall>)>);

/// `SimpleNovelAdjacencyInterpreter.inferTypeFromSingleContigSimpleChimera`.
pub fn simple_calls(
    classified: &Classified,
    header: &SamHeader,
    reference: &mut dyn Reference,
) -> Result<SimpleCalls, Failed> {
    let mut found = Vec::new();
    for tig in &classified.simple {
        let (head, tail) = (tig.head(), tig.tail());
        // `splitPairStrongEnoughEvidenceForCA`.
        if head.map_quality < SPLIT_PAIR_MIN_MQ
            || tail.map_quality < SPLIT_PAIR_MIN_MQ
            || size_on_read(head).min(size_on_read(tail)) < SPLIT_PAIR_MIN_LENGTH
        {
            continue;
        }
        let chimera = Chimera::new(
            head,
            tail,
            &tig.insertion_mappings,
            &tig.contig.name,
            &tig.non_canonical,
            header,
        )?;
        let adjacency = novel_adjacency(&chimera, &tig.contig.sequence, header)?;
        found.push((adjacency, chimera));
    }
    let mut calls = Vec::new();
    for (adjacency, evidence) in append_only_group(found, NovelAdjacency::hash_code) {
        let typed = typed_calls(&adjacency, reference)?;
        calls.push((adjacency, evidence, typed));
    }
    Ok(SimpleCalls(calls))
}

impl SimpleCalls {
    /// `readBreakpointsFromNarls`, as in [`Adjacencies::check_padded_breakpoints`].
    pub fn check_padded_breakpoints(&self, padding: i32) -> Result<(), Failed> {
        for (adjacency, _, _) in &self.0 {
            let padded = padding + adjacency.complication.homology.len() as i32;
            for location in [&adjacency.left, &adjacency.right] {
                check_padded_start(location.start - padded)?;
            }
        }
        Ok(())
    }

    /// `SimpleNovelAdjacencyAndChimericAlignmentEvidence.toVariantContexts`, unsorted.
    pub fn variants(
        &self,
        header: &SamHeader,
        cnv_calls: Option<&[CnvCall]>,
    ) -> Result<Vec<VariantContext>, Failed> {
        let mut variants = Vec::with_capacity(self.0.len());
        for (adjacency, evidence, typed) in &self.0 {
            if let [only] = typed.as_slice() {
                variants.push(annotated_call(
                    only, adjacency, evidence, header, cnv_calls,
                )?);
                continue;
            }
            // `produceLinkedAssemblyBasedVariants`.
            let (one, two) = (&typed[0], &typed[1]);
            let key = if one.breakend { "MATEID" } else { "LINK" };
            let mut first = annotated_call(one, adjacency, evidence, header, cnv_calls)?;
            let mut second = annotated_call(two, adjacency, evidence, header, cnv_calls)?;
            let (first_id, second_id) = (first.id.clone(), second.id.clone());
            set(&mut first, key, Value::Str(second_id));
            set(&mut second, key, Value::Str(first_id));
            for (call, variant) in [(one, &mut first), (two, &mut second)] {
                if call.deletion {
                    remove(variant, "INSSEQ");
                    remove(variant, "INSLEN");
                    remove(variant, "SEQ_ALT_HAPLOTYPE");
                }
            }
            variants.push(first);
            variants.push(second);
        }
        Ok(variants)
    }
}

/// [`sort_variants`], for the runner's writes.
pub fn sorted(variants: &[VariantContext], header: &SamHeader) -> Vec<VariantContext> {
    let mut copy = variants.to_vec();
    sort_variants(&mut copy, header);
    copy
}

/// `CigarOperator.isPadding` or `N`, which `computeNewRefSpanAndCigar` refuses.
fn clip_alignment(input: &Alignment, clip: i32, from_3_prime: bool) -> Result<Alignment, Failed> {
    if clip < 0 {
        return Err(format!(
            "requesting negative clip length: {clip} on {}",
            input.packed()
        ));
    }
    if clip >= size_on_read(input) {
        return Err(format!(
            "input alignment to be clipped away: {}\twith clip length: {clip}",
            input.packed()
        ));
    }
    if input
        .cigar
        .iter()
        .any(|element| matches!(element.op, Op::N | Op::P))
    {
        return Err(format!(
            "Input alignment contains padding or skip operations, which is currently unsupported: \
             {}",
            input.packed()
        ));
    }
    // `splitCigarByLeftAndRightClipping`.
    let (mut left, mut middle, mut right) = (Vec::new(), Vec::new(), Vec::new());
    let mut consumed = 0;
    for element in &input.cigar {
        if consumed < input.start_on_contig - 1 {
            left.push(*element);
        } else if consumed < input.end_on_contig {
            middle.push(*element);
        } else {
            right.push(*element);
        }
        if consumes_read(element.op) || element.op == Op::H {
            consumed += element.length as i32;
        }
    }
    if middle.is_empty() {
        return Err(format!(
            "Logic error: cigar elements corresponding to alignment block is empty. {}",
            input.packed()
        ));
    }
    if from_3_prime {
        middle.reverse();
    }
    let element = |length: i32, op: Op| CigarElement {
        length: length as u32,
        op,
    };
    let (mut read_consumed, mut reference_consumed) = (0, 0);
    let mut section: Vec<CigarElement> = Vec::with_capacity(middle.len());
    for (index, current) in middle.iter().enumerate() {
        let length = current.length as i32;
        if consumes_read(current.op) {
            if read_consumed + length < clip {
                read_consumed += length;
            } else {
                if !is_alignment(current.op) && current.op != Op::I {
                    return Err(format!(
                        "Logic error, should not reach here: operation consumes read bases but is \
                         neither alignment nor insertion.\n Original cigar({})\t\
                         clipLengthOnRead({clip})\t{}",
                        cigar_text(&input.cigar),
                        if from_3_prime {
                            "clipFrom3PrimeEnd"
                        } else {
                            "clipFrom5PrimeEnd"
                        }
                    ));
                }
                section.push(element(clip, Op::S));
                let left_over = read_consumed + length - clip;
                if left_over != 0 {
                    section.push(element(
                        left_over,
                        if is_alignment(current.op) {
                            Op::M
                        } else {
                            Op::S
                        },
                    ));
                }
                section.extend_from_slice(&middle[index + 1..]);
                if is_alignment(current.op) {
                    reference_consumed += clip - read_consumed;
                }
                break;
            }
        }
        if consumes_reference(current.op) {
            reference_consumed += length;
        }
    }
    let when_erred = format!(
        "{}\tclip length: {clip}\tclip from end: {} of read",
        input.packed(),
        if from_3_prime { "3" } else { "5" }
    );
    if section.len() < 2 {
        return Err(format!(
            "The input alignment after clipping contains no or only one operation. This indicates \
             a logic error or too large a clipping length (whole alignment being clipped away) or \
             invalid input cigar.\n{when_erred}"
        ));
    }
    let (first, second) = (section[0], section[1]);
    if first.op == Op::S && is_indel(second.op) {
        if section.len() < 3 {
            return Err(format!(
                "After clipping, the new cigar would contain NO aligned bases (i.e. either clipped \
                 or gap or mix), This indicates a logic error or too large a clipping length \
                 (whole alignment being clipped away) or invalid input cigar.\n{when_erred}"
            ));
        }
        if consumes_reference(second.op) {
            reference_consumed += second.length as i32;
        }
        if second.op == Op::D {
            section.remove(1);
        } else {
            section.remove(0);
            section[0] = element((first.length + second.length) as i32, Op::S);
        }
    }
    if from_3_prime {
        section.reverse();
    }
    let mut builder = CigarBuilder::new(true);
    for current in left.iter().chain(&section).chain(&right) {
        builder
            .add(*current)
            .map_err(|error| format!("{error:?}"))?;
    }
    let cigar = builder
        .make(false)
        .map_err(|error| format!("{error:?}"))?
        .elements;
    let span = if from_3_prime == input.forward {
        Interval::new(
            &input.span.contig,
            input.span.start,
            input.span.end - reference_consumed,
        )?
    } else {
        Interval::new(
            &input.span.contig,
            input.span.start + reference_consumed,
            input.span.end,
        )?
    };
    // `computeNewReadSpan`.
    let (start_on_contig, end_on_contig) = if from_3_prime {
        (
            input.start_on_contig,
            (input.end_on_contig - clip)
                .min(unclipped_read_bases(&cigar) - clipped_bases(&cigar, false)?),
        )
    } else {
        (
            (input.start_on_contig + clip).max(clipped_bases(&cigar, true)? + 1),
            input.end_on_contig,
        )
    };
    Ok(Alignment {
        span,
        start_on_contig,
        end_on_contig,
        cigar,
        forward: input.forward,
        map_quality: input.map_quality,
        mismatches: NO_VALUE,
        score: NO_VALUE,
        modification: 'H',
    })
}

/// `CpxVariantInterpreter.yieldOverlapToAlignmentTwo`.
fn yield_overlap_to_two(one: &Alignment, two: &Alignment, header: &SamHeader) -> bool {
    if one.span.contig == two.span.contig {
        if one.forward != two.forward {
            (one.span.end <= two.span.start) == one.forward
        } else {
            one.forward
        }
    } else {
        contig_index(header, &one.span.contig) > contig_index(header, &two.span.contig)
    }
}

/// `CpxVariantInterpreter.deOverlapAlignments`.
fn de_overlap_alignments(
    original: &[Alignment],
    header: &SamHeader,
) -> Result<Vec<Alignment>, Failed> {
    let mut result = Vec::with_capacity(original.len());
    let mut one = original[0].clone();
    for two in &original[1..] {
        if two.modification == 'S' && contains_on_read(&one, two) {
            continue;
        }
        let overlap = overlap_on_contig(&one, two);
        if overlap == 0 {
            if size_on_read(&one) >= MIN_READ_SPAN_AFTER_DEOVERLAP {
                result.push(one);
            }
            one = two.clone();
            continue;
        }
        if contains_on_read(&one, two) || contains_on_read(two, &one) {
            return Err(format!(
                "Two input alignments' overlap on read consumes completely one of them.\t{}\t{}",
                one.packed(),
                two.packed()
            ));
        }
        // `removeOverlap`.
        let yields = if one == original[0] {
            false
        } else if *two == original[original.len() - 1] {
            true
        } else {
            yield_overlap_to_two(&one, two, header)
        };
        let (first, second) = if yields {
            (clip_alignment(&one, overlap, true)?, two.clone())
        } else {
            (one, clip_alignment(two, overlap, false)?)
        };
        if size_on_read(&first) >= MIN_READ_SPAN_AFTER_DEOVERLAP {
            result.push(first);
        }
        one = second;
    }
    result.push(one);
    if overlap_on_contig(&result[0], &result[result.len() - 1]) > 0 {
        return Err(format!(
            "After alignment preprocessing, head and tail alignments should not overlap on read; \
             otherwise it would mean some middle alignments are contained in head/tail which \
             should not happen.[{}]",
            original
                .iter()
                .map(Alignment::packed)
                .collect::<Vec<_>>()
                .join(", ")
        ));
    }
    Ok(result)
}

/// `CpxVariantInducingAssemblyContig.Jump`.
struct Jump {
    start: Interval,
    landing: Interval,
    gap: i32,
}

fn jump(one: &Alignment, two: &Alignment) -> Result<Jump, Failed> {
    if overlap_on_contig(one, two) > 0 {
        return Err("assumption that input alignments DO NOT overlap is violated.".to_string());
    }
    let start_at = if one.forward {
        one.span.end
    } else {
        one.span.start
    };
    let landing_at = if two.forward {
        two.span.start
    } else {
        two.span.end
    };
    Ok(Jump {
        start: Interval::new(&one.span.contig, start_at, start_at)?,
        landing: Interval::new(&two.span.contig, landing_at, landing_at)?,
        gap: (two.start_on_contig - one.end_on_contig - 1).max(0),
    })
}

fn intersect(one: &Interval, two: &Interval) -> Interval {
    Interval {
        contig: one.contig.clone(),
        start: one.start.max(two.start),
        end: one.end.min(two.end),
    }
}

/// `AlignmentInterval.readIntervalAlignedToRefSpan`.
fn read_interval_aligned_to(alignment: &Alignment, other: &Interval) -> Result<(i32, i32), Failed> {
    if !other.overlaps(&alignment.span) {
        return Ok((-1, -1));
    }
    if other.contains(&alignment.span) {
        return Ok((alignment.start_on_contig, alignment.end_on_contig));
    }
    let hard_offset = match alignment.cigar.first() {
        Some(first) if first.op == Op::H => first.length as i32,
        _ => 0,
    };
    let lower = (other.start - alignment.span.start).max(0);
    let upper = (alignment.span.end - other.end).max(0);
    let (from_start, from_end) = if alignment.forward {
        (lower, upper)
    } else {
        (upper, lower)
    };
    let walk_start = if from_start != 0 {
        associated_distance_on_read(
            &alignment.cigar,
            alignment.start_on_contig - hard_offset,
            from_start,
            false,
        )?
    } else {
        0
    };
    let walk_end = if from_end != 0 {
        associated_distance_on_read(
            &alignment.cigar,
            alignment.end_on_contig - hard_offset,
            from_end,
            true,
        )?
    } else {
        0
    };
    Ok((
        alignment.start_on_contig + walk_start,
        alignment.end_on_contig - walk_end,
    ))
}

/// `CpxVariantCanonicalRepresentation`.
struct Canonical {
    affected: Interval,
    segments: Vec<Interval>,
    descriptions: Vec<String>,
    alternate: Vec<u8>,
}

impl PartialEq for Canonical {
    fn eq(&self, other: &Self) -> bool {
        self.segments == other.segments
            && self.affected == other.affected
            && self.descriptions == other.descriptions
            && self.alternate == other.alternate
    }
}

/// `List.hashCode`.
fn list_hash(values: impl Iterator<Item = i32>) -> i32 {
    values.fold(1, fold)
}

impl Canonical {
    fn hash_code(&self) -> i32 {
        let mut result = list_hash(self.segments.iter().map(Interval::hash_code));
        result = fold(result, self.affected.hash_code());
        result = fold(
            result,
            list_hash(self.descriptions.iter().map(|text| string_hash_code(text))),
        );
        fold(result, byte_array_hash_code(&self.alternate))
    }
}

/// The insertion arrangement of `descriptionsForSingleSharedBaseByHeadTail`.
fn shared_base_arrangement(alignments: &[Alignment], jumps: &[Jump], forward: bool) -> Vec<String> {
    let mut result = vec!["1".to_string()];
    if jumps[0].gap > 0 {
        result.push(format!("UINS-{}", jumps[0].gap));
    }
    for index in 1..alignments.len() - 1 {
        let alignment = &alignments[index];
        result.push(if forward == alignment.forward {
            alignment.span.text()
        } else {
            format!("-{}", alignment.span.text())
        });
        if jumps[index].gap > 0 {
            result.push(format!("UINS-{}", jumps[index].gap));
        }
    }
    result.push("1".to_string());
    result
}

/// `extractAltArrangements`.
fn alt_arrangements(
    region: &Interval,
    forward: bool,
    alignments: &[Alignment],
    jumps: &[Jump],
    segments: &[Interval],
) -> Vec<String> {
    if segments.len() == 1 && segments[0].size() == 1 {
        return shared_base_arrangement(alignments, jumps, forward);
    }
    let walk: Vec<&Alignment> = if forward {
        alignments.iter().collect()
    } else {
        alignments.iter().rev().collect()
    };
    let gaps: Vec<i32> = if forward {
        jumps.iter().map(|jump| jump.gap).collect()
    } else {
        jumps.iter().rev().map(|jump| jump.gap).collect()
    };
    let mut gap_at = 0;
    let mut current_gap = gaps[gap_at];
    let mut jump_is_last = false;
    let mut descriptions: Vec<String> = Vec::with_capacity(2 * segments.len());
    let mut disjoint: Vec<(String, i32)> = Vec::new();
    for alignment in walk {
        if !alignment.span.overlaps(region) {
            let at = descriptions.len() as i32;
            disjoint.push((
                alignment.span.text(),
                if forward == alignment.forward {
                    at
                } else {
                    -at
                },
            ));
            if current_gap > 0 {
                descriptions.push(format!("UINS-{current_gap}"));
            }
        } else {
            let order: Vec<usize> = if forward == alignment.forward {
                (0..segments.len()).collect()
            } else {
                (0..segments.len()).rev().collect()
            };
            for index in order {
                if alignment.span.contains(&segments[index]) {
                    let number = index as i32 + 1;
                    let positive = forward == alignment.forward;
                    descriptions.push(if positive {
                        number.to_string()
                    } else {
                        (-number).to_string()
                    });
                }
            }
            if current_gap > 0 && !jump_is_last {
                descriptions.push(format!("UINS-{current_gap}"));
            }
        }
        gap_at += 1;
        if gap_at < gaps.len() {
            current_gap = gaps[gap_at];
        } else {
            jump_is_last = true;
        }
    }
    for (span, index) in disjoint.into_iter().rev() {
        if index > 0 {
            descriptions.insert(index as usize, span);
        } else {
            descriptions.insert((-index) as usize, format!("-{span}"));
        }
    }
    descriptions
}

/// `new CpxVariantInducingAssemblyContig` and `new CpxVariantCanonicalRepresentation`.
fn canonical_representation(tig: &FineTuned, header: &SamHeader) -> Result<Canonical, Failed> {
    let alignments = &tig.contig.alignments;
    let (head, tail) = (tig.head(), tig.tail());
    let chromosome = head.span.contig.as_str();
    let forward = head.forward;
    let (alpha, omega) = if forward {
        (head.span.start, tail.span.end)
    } else {
        (tail.span.start, head.span.end)
    };
    let mut jumps = Vec::with_capacity(alignments.len() - 1);
    for pair in alignments.windows(2) {
        jumps.push(jump(&pair[0], &pair[1])?);
    }
    let region = Interval::new(chromosome, alpha, omega)?;
    let mut locations: Vec<Interval> = jumps
        .iter()
        .flat_map(|jump| [jump.start.clone(), jump.landing.clone()])
        .filter(|location| location.overlaps(&region))
        .collect();
    locations.sort_by(|one, two| {
        compare_locatables(
            header,
            (&one.contig, one.start, one.end),
            (&two.contig, two.start, two.end),
        )
    });
    locations.dedup();
    let invalid_location =
        || "List index out of bounds: no segmenting location was found.".to_string();

    // `checkBoundedBySharedSingleBase`.
    if locations.len() == 1 {
        let unseen = |what: &str| {
            format!(
                "run into unseen case where only one reference segmenting location is found{what}"
            )
        };
        if locations[0].size() != 1 {
            return Err(unseen(" but its size is not 1:\n"));
        }
        if alignments[1..alignments.len() - 1]
            .iter()
            .any(|middle| middle.span.overlaps(&region))
        {
            return Err(unseen(&format!(
                " but some middle alignments are overlapping alpha-omega region:\t{}\n",
                region.text()
            )));
        }
    }
    let first_location = locations.first().ok_or_else(invalid_location)?;
    let last_location = &locations[locations.len() - 1];
    let affected = Interval::new(
        &first_location.contig,
        first_location.start,
        last_location.end,
    )?;

    // `getTwoBaseBoundaries` and `extractRefSegments`.
    let mut two_base = Vec::with_capacity(2 * alignments.len());
    for alignment in alignments {
        if region.contains(&alignment.span) {
            two_base.push(Interval::new(
                chromosome,
                alignment.span.start,
                alignment.span.start + 1,
            )?);
            two_base.push(Interval::new(
                chromosome,
                alignment.span.end - 1,
                alignment.span.end,
            )?);
        }
    }
    let segments: Vec<Interval> = if locations.len() == 1 {
        locations.clone()
    } else {
        let mut segments = Vec::with_capacity(locations.len() - 1);
        for pair in locations.windows(2) {
            let distance = pair[1].start - pair[0].end;
            if distance >= 1 {
                let segment = Interval::new(chromosome, pair[0].end, pair[1].start)?;
                if distance > 1 || two_base.contains(&segment) {
                    segments.push(segment);
                }
            }
        }
        segments
    };
    let descriptions = alt_arrangements(&region, forward, alignments, &jumps, &segments);

    // `extractAltHaplotypeSeq`.
    let sequence = &tig.contig.sequence;
    let (start, end);
    if segments.is_empty() {
        start = head.end_on_contig;
        end = tail.start_on_contig;
    } else {
        let (first, last) = if forward {
            (&segments[0], &segments[segments.len() - 1])
        } else {
            (&segments[segments.len() - 1], &segments[0])
        };
        start = if !first.overlaps(&head.span) {
            let neighbors = if forward {
                first.start - head.span.end == 1
            } else {
                head.span.start - first.end == 1
            };
            if !neighbors {
                return Err(format!(
                    "1st segment is not overlapping with head alignment but it is not immediately \
                     before/after the head alignment either\n{}\nSegments:\t[{}]",
                    tig.text(),
                    segments
                        .iter()
                        .map(Interval::text)
                        .collect::<Vec<_>>()
                        .join(", ")
                ));
            }
            head.end_on_contig
        } else {
            let common = intersect(first, &head.span);
            if common.size() == 1 {
                head.end_on_contig
            } else {
                read_interval_aligned_to(head, &common)?.0
            }
        };
        end = if !last.overlaps(&tail.span) {
            let expected = if forward {
                tail.span.start - last.end == 1
            } else {
                last.start - tail.span.end == 1
            };
            if !expected {
                return Err(tig.text());
            }
            tail.start_on_contig
        } else {
            let common = intersect(last, &tail.span);
            if common.size() == 1 {
                tail.start_on_contig
            } else {
                read_interval_aligned_to(tail, &common)?.1
            }
        };
    }
    let mut alternate = copy_of_range(sequence, start - 1, end)?;
    if !forward {
        alternate = reverse_complement(&alternate);
    }
    Ok(Canonical {
        affected,
        segments,
        descriptions,
        alternate,
    })
}

/// A complex call, as much as the extraction reads back from its `VariantContext`.
pub struct ComplexCall {
    variant: VariantContext,
    segments: Vec<Interval>,
    arrangement: Vec<String>,
    alternate_length: i32,
    contigs: Vec<String>,
    qualities: Vec<String>,
    maximum_length: i32,
}

/// `CpxVariantInterpreter.makeInterpretation`.
pub fn complex_calls(
    classified: &Classified,
    header: &SamHeader,
    reference: &mut dyn Reference,
) -> Result<Vec<ComplexCall>, Failed> {
    let mut found = Vec::with_capacity(classified.complex.len());
    for tig in &classified.complex {
        // `furtherPreprocess`.
        let mut alignments = de_overlap_alignments(&tig.contig.alignments, header)?;
        sort_alignments(&mut alignments);
        let processed = FineTuned {
            contig: Contig {
                name: tig.contig.name.clone(),
                sequence: tig.contig.sequence.clone(),
                alignments,
            },
            insertion_mappings: tig.insertion_mappings.clone(),
            ambiguous: tig.ambiguous,
            non_canonical: tig.non_canonical.clone(),
        };
        found.push((canonical_representation(&processed, header)?, processed));
    }
    let mut calls = Vec::new();
    for (canonical, mut evidence) in append_only_group(found, Canonical::hash_code) {
        let affected = &canonical.affected;
        let bases = reference.bases(&affected.contig, affected.start, affected.start)?;
        let mut variant = simple_call(
            &affected.contig,
            affected.start,
            affected.end,
            format!("CPX_{}", affected.text()),
            vec![allele(&bases, true)?, symbolic("<CPX>")?],
            "CPX",
            canonical.alternate.len() as i32 - affected.size(),
        );
        set(
            &mut variant,
            "ALT_ARRANGEMENT",
            Value::Str(canonical.descriptions.join(",")),
        );
        if !canonical.segments.is_empty() {
            set(
                &mut variant,
                "SEGMENTS",
                Value::Str(
                    canonical
                        .segments
                        .iter()
                        .map(Interval::text)
                        .collect::<Vec<_>>()
                        .join(","),
                ),
            );
        }
        set(
            &mut variant,
            "SEQ_ALT_HAPLOTYPE",
            Value::Str(text_of(canonical.alternate.clone())),
        );

        evidence.sort_by(|one, two| one.contig.name.cmp(&two.contig.name));
        let contigs: Vec<String> = evidence.iter().map(|tig| tig.contig.name.clone()).collect();
        let tags: Vec<String> = evidence
            .iter()
            .map(|tig| {
                if tig.non_canonical == NO_NON_CANONICAL {
                    ".".to_string()
                } else {
                    tig.non_canonical.clone()
                }
            })
            .collect();
        let minimum_qualities: Vec<i32> = evidence
            .iter()
            .map(|tig| tig.head().map_quality.min(tig.tail().map_quality))
            .collect();
        let qualities: Vec<String> = minimum_qualities.iter().map(i32::to_string).collect();
        let lengths: Vec<i32> = evidence
            .iter()
            .map(|tig| size_on_read(tig.head()).min(size_on_read(tig.tail())))
            .collect();
        let maximum_length = lengths.iter().copied().max().unwrap_or(0);
        set(
            &mut variant,
            "TOTAL_MAPPINGS",
            Value::Str(contigs.len().to_string()),
        );
        set(&mut variant, "CTG_NAMES", Value::Str(contigs.join(",")));
        if tags.iter().any(|tag| tag != ".") {
            set(
                &mut variant,
                "CTG_GOOD_NONCANONICAL_MAPPING",
                Value::Str(tags.join(",")),
            );
        }
        set(
            &mut variant,
            "HQ_MAPPINGS",
            Value::Str(
                minimum_qualities
                    .iter()
                    .filter(|quality| **quality >= HIGH_MQ)
                    .count()
                    .to_string(),
            ),
        );
        set(
            &mut variant,
            "MAPPING_QUALITIES",
            Value::Str(qualities.join(",")),
        );
        set(
            &mut variant,
            "ALIGN_LENGTHS",
            Value::Str(
                lengths
                    .iter()
                    .map(i32::to_string)
                    .collect::<Vec<_>>()
                    .join(","),
            ),
        );
        set(
            &mut variant,
            "MAX_ALIGN_LENGTH",
            Value::Str(maximum_length.to_string()),
        );
        calls.push(ComplexCall {
            variant,
            alternate_length: canonical.alternate.len() as i32,
            segments: canonical.segments,
            arrangement: canonical.descriptions,
            contigs,
            qualities,
            maximum_length,
        });
    }
    Ok(calls)
}

/// `SVUtils.getAttributeAsStringList`: nothing when the key is absent, and every comma a separator.
fn attribute_strings(variant: &VariantContext, key: &str) -> Vec<String> {
    let Some((_, value)) = variant.attributes.iter().find(|(name, _)| name == key) else {
        return Vec::new();
    };
    let items = match value {
        Value::List(items) => items.iter().filter_map(Value::format).collect::<Vec<_>>(),
        other => other.format().into_iter().collect(),
    };
    items
        .iter()
        .flat_map(|item| item.split(','))
        .map(str::to_string)
        .collect()
}

/// The records of a complex VCF as `SegmentedCpxVariantSimpleVariantExtractor` reads them back: the
/// zero- and one-segment calls' pass comes first, and whichever record it cannot read first fails.
pub fn complex_calls_from_variants(
    variants: Vec<VariantContext>,
) -> Result<Vec<ComplexCall>, Failed> {
    let mut calls: Vec<Option<ComplexCall>> = variants.iter().map(|_| None).collect();
    for few in [true, false] {
        for (at, variant) in variants.iter().enumerate() {
            let texts = attribute_strings(variant, "SEGMENTS");
            if (texts.len() < 2) != few {
                continue;
            }
            let segments = texts
                .iter()
                .map(|text| parse_interval(text))
                .collect::<Result<Vec<_>, _>>()?;
            let maximum_length = match attribute_strings(variant, "MAX_ALIGN_LENGTH").first() {
                Some(text) => text
                    .parse()
                    .map_err(|_| format!("For input string: \"{text}\""))?,
                None => 0,
            };
            calls[at] = Some(ComplexCall {
                variant: variant.clone(),
                segments,
                arrangement: attribute_strings(variant, "ALT_ARRANGEMENT"),
                alternate_length: attribute_strings(variant, "SEQ_ALT_HAPLOTYPE")
                    .first()
                    .map_or(0, |text| text.len() as i32),
                contigs: attribute_strings(variant, "CTG_NAMES"),
                qualities: attribute_strings(variant, "MAPPING_QUALITIES"),
                maximum_length,
            });
        }
    }
    Ok(calls.into_iter().flatten().collect())
}

/// The complex calls' `VariantContext`s, for `Complex.vcf`.
pub fn complex_variants(calls: &[ComplexCall]) -> Vec<VariantContext> {
    calls.iter().map(|call| call.variant.clone()).collect()
}

/// `makeID`.
fn derived_id(kind: &str, contig: &str, start: i32, stop: i32) -> String {
    format!("{kind}_{contig}_{start}_{stop}_CPX_DERIVED")
}

/// `getAnchorBaseRefAllele`.
fn anchor(reference: &mut dyn Reference, contig: &str, position: i32) -> Result<Allele, Failed> {
    Interval::new(contig, position, position)?;
    allele(&reference.bases(contig, position, position)?, true)
}

/// `makeDeletion`.
fn make_deletion(range: &Interval, reference_allele: Allele) -> Result<VariantContext, Failed> {
    Ok(simple_call(
        &range.contig,
        range.start,
        range.end,
        derived_id("DEL", &range.contig, range.start, range.end),
        vec![reference_allele, symbolic("<DEL>")?],
        "DEL",
        -range.size() + 1,
    ))
}

/// `makeInsertion`.
fn make_insertion(
    contig: &str,
    position: i32,
    end: i32,
    length: i32,
    reference_allele: Allele,
) -> Result<VariantContext, Failed> {
    Ok(simple_call(
        contig,
        position,
        end,
        derived_id("INS", contig, position, end),
        vec![reference_allele, symbolic("<INS>")?],
        "INS",
        length,
    ))
}

/// `makeInversion`.
fn make_inversion(region: &Interval, reference_allele: Allele) -> Result<VariantContext, Failed> {
    Ok(simple_call(
        &region.contig,
        region.start - 1,
        region.end,
        derived_id("INV", &region.contig, region.start - 1, region.end),
        vec![reference_allele, symbolic("<INV>")?],
        "INV",
        0,
    ))
}

/// `new SimpleInterval(String)`, of what this tool printed.
fn parse_interval(text: &str) -> Result<Interval, Failed> {
    // `parsePosition`, whose `UserException` a Spark failure reports by its cause.
    let position = |digits: &str| -> Result<i32, Failed> {
        digits
            .replace(',', "")
            .parse()
            .map_err(|_| format!("For input string: \"{}\"", digits.replace(',', "")))
    };
    let Some((contig, range)) = text.rsplit_once(':') else {
        return Interval::new(text, 1, i32::MAX);
    };
    match range.split_once('-') {
        Some((start, end)) => Interval::new(contig, position(start)?, position(end)?),
        None => match range.strip_suffix('+') {
            Some(start) => Interval::new(contig, position(start)?, i32::MAX),
            None => {
                let start = position(range)?;
                Interval::new(contig, start, start)
            }
        },
    }
}

/// `getInsLen`.
fn insertion_length(description: &str, segment_lengths: &[i32]) -> Result<i32, Failed> {
    if let Some(length) = description.strip_prefix("UINS-") {
        return length.parse().map_err(|_| description.to_string());
    }
    let bare = description.strip_prefix('-').unwrap_or(description);
    if let Ok(number) = bare.parse::<usize>() {
        return Ok(segment_lengths[number - 1]);
    }
    Ok(parse_interval(bare)?.size())
}

/// `getInsFromOneEnd`, always told to increase the length by one.
fn insertion_from_one_end(
    from_front: bool,
    index: usize,
    contig: &str,
    position: i32,
    reference_allele: Allele,
    segment_lengths: &[i32],
    arrangement: &[String],
) -> Result<Option<VariantContext>, Failed> {
    // A rear walk starts after the last segment, which is before the first entry if there is none.
    let walk = if from_front {
        &arrangement[..index]
    } else {
        &arrangement[index..]
    };
    let mut length = 1;
    for description in walk {
        length += insertion_length(description, segment_lengths)?;
    }
    if length > EVENT_SIZE_THRESHOLD {
        Ok(Some(make_insertion(
            contig,
            position,
            position,
            length,
            reference_allele,
        )?))
    } else {
        Ok(None)
    }
}

/// What both extractors add to each of their records.
fn with_source(mut variant: VariantContext, call: &ComplexCall) -> VariantContext {
    // A file's record may name none, and an empty list is written as a missing value.
    let list = |items: &[String]| {
        if items.is_empty() {
            Value::Missing
        } else {
            Value::Str(items.join(","))
        }
    };
    let qualities = list(&call.qualities);
    set(
        &mut variant,
        "CPX_EVENT",
        Value::Str(call.variant.id.clone()),
    );
    set(&mut variant, "CTG_NAMES", list(&call.contigs));
    set(&mut variant, "MAPPING_QUALITIES", qualities);
    set(
        &mut variant,
        "MAX_ALIGN_LENGTH",
        Value::Int(i64::from(call.maximum_length)),
    );
    variant
}

/// `ZeroAndOneSegmentCpxVariantExtractor.extract`.
fn extract_zero_or_one_segment(
    call: &ComplexCall,
    reference: &mut dyn Reference,
) -> Result<Vec<VariantContext>, Failed> {
    let contig = call.variant.contig.as_str();
    let start = call.variant.start as i32;
    let Some(segment) = call.segments.first() else {
        let insertion = make_insertion(
            contig,
            start,
            start,
            call.alternate_length - 2,
            anchor(reference, contig, start)?,
        )?;
        return Ok(vec![with_source(insertion, call)]);
    };
    let arrangement = &call.arrangement;
    let lengths = [segment.size()];
    let mut result = Vec::new();
    let as_is = arrangement.iter().position(|text| text == "1");
    let inverted = arrangement.iter().position(|text| text == "-1");
    let front_and_rear = |index: usize,
                          front: Allele,
                          rear: Allele,
                          result: &mut Vec<VariantContext>|
     -> Result<(), Failed> {
        result.extend(insertion_from_one_end(
            true,
            index,
            contig,
            segment.start - 1,
            front,
            &lengths,
            arrangement,
        )?);
        result.extend(insertion_from_one_end(
            false,
            index + 1,
            contig,
            segment.end,
            rear,
            &lengths,
            arrangement,
        )?);
        Ok(())
    };
    match (inverted, as_is) {
        (Some(index), _) if segment.size() > EVENT_SIZE_THRESHOLD => {
            let inversion_anchor = anchor(reference, contig, segment.start)?;
            result.push(make_inversion(segment, inversion_anchor)?);
            let front = anchor(reference, contig, segment.start - 1)?;
            let rear = anchor(reference, contig, segment.end)?;
            front_and_rear(index, front, rear, &mut result)?;
        }
        (_, Some(index)) => {
            if call.alternate_length - segment.size() > EVENT_SIZE_THRESHOLD {
                let front = anchor(reference, contig, segment.start - 1)?;
                let rear = anchor(reference, contig, segment.end)?;
                let length = call.alternate_length - segment.size();
                if arrangement[arrangement.len() - 1] == "1" {
                    result.push(make_insertion(
                        contig,
                        segment.start - 1,
                        segment.start - 1,
                        length,
                        front,
                    )?);
                } else if arrangement[0] == "1" {
                    result.push(make_insertion(
                        contig,
                        segment.end,
                        segment.end,
                        length,
                        front,
                    )?);
                } else {
                    front_and_rear(index, front, rear, &mut result)?;
                }
            }
        }
        _ => {
            if segment.size() > EVENT_SIZE_THRESHOLD {
                let front = anchor(reference, contig, segment.start)?;
                result.push(make_deletion(
                    &Interval::new(contig, segment.start, segment.end - 1)?,
                    front.clone(),
                )?);
                if call.alternate_length - 2 > EVENT_SIZE_THRESHOLD {
                    result.push(make_insertion(
                        contig,
                        segment.start,
                        segment.start,
                        call.alternate_length,
                        front,
                    )?);
                }
            } else if call.alternate_length - 2 > EVENT_SIZE_THRESHOLD {
                Interval::new(contig, segment.start, segment.end - 1)?;
                let bases = reference.bases(contig, segment.start, segment.end - 1)?;
                result.push(make_insertion(
                    contig,
                    segment.start,
                    segment.end - 1,
                    call.alternate_length - segment.size(),
                    allele(&bases, true)?,
                )?);
            }
        }
    }
    Ok(result
        .into_iter()
        .map(|variant| with_source(variant, call))
        .collect())
}

/// The zero- and one-segment complex calls' derived simple variants.
pub fn zero_or_one_segment_variants(
    calls: &[ComplexCall],
    reference: &mut dyn Reference,
) -> Result<Vec<VariantContext>, Failed> {
    let mut variants = Vec::new();
    for call in calls.iter().filter(|call| call.segments.len() < 2) {
        variants.extend(extract_zero_or_one_segment(call, reference)?);
    }
    Ok(variants)
}

/// The missing segments, then the indices of the present ones and of the inverted ones.
type SegmentAppearances = (Vec<Interval>, Vec<i32>, Vec<i32>);

/// `getMissingAndPresentAndInvertedSegments`, the missing ones in `HashSet` order.
fn segment_appearances(call: &ComplexCall) -> Result<SegmentAppearances, Failed> {
    let mut inverted = Vec::new();
    let mut present = Vec::new();
    for description in &call.arrangement {
        if description.contains(':') || description.starts_with("UINS") {
            continue;
        }
        if let Some(number) = description.strip_prefix('-') {
            inverted.push(number.parse::<i32>().map_err(|_| description.clone())?);
        } else {
            let number = description
                .parse::<i32>()
                .map_err(|_| description.clone())?;
            if !present.contains(&number) {
                present.push(number);
            }
        }
    }
    let mut missing = Vec::new();
    for (index, segment) in call.segments.iter().enumerate() {
        let number = index as i32 + 1;
        if !present.contains(&number) && !inverted.contains(&number) {
            missing.push((segment.clone(), segment.hash_code()));
        }
    }
    let missing = hash_map_order(&missing).map_err(|error| format!("{error:?}"))?;
    Ok((missing, present, inverted))
}

/// `findAllSegments`.
fn find_all_segments(arrangement: &[String], count: usize) -> Option<usize> {
    let mut looking_for = count;
    let last = count.to_string();
    for (index, description) in arrangement.iter().enumerate().rev() {
        if *description == looking_for.to_string() {
            if looking_for == 1 {
                return Some(index);
            }
            looking_for -= 1;
        } else {
            looking_for = if *description == last {
                count - 1
            } else {
                count
            };
        }
    }
    None
}

/// `compactifyMissingSegments`.
fn compactify_missing(missing: &[Interval]) -> Vec<Interval> {
    if missing.len() == 1 {
        return missing.to_vec();
    }
    let mut sorted = missing.to_vec();
    sorted.sort_by_key(|segment| segment.start);
    let mut result = Vec::with_capacity(sorted.len());
    let mut walk = sorted.into_iter();
    let mut current = walk.next().expect("missing segments were checked");
    for next in walk {
        if current.contig == next.contig
            && current.start <= next.end + 1
            && next.start <= current.end + 1
        {
            current.end = next.end;
        } else {
            result.push(current);
            current = next;
        }
    }
    result.push(current);
    result
}

/// `descriptionIndicatesInsertion`.
fn indicates_insertion(description: &str) -> bool {
    description.starts_with("UINS") || description.parse::<i32>().is_err()
}

/// `MultiSegmentsCpxVariantExtractor.extract`.
fn extract_multi_segment(
    call: &ComplexCall,
    reference: &mut dyn Reference,
) -> Result<Vec<VariantContext>, Failed> {
    let contig = call.variant.contig.as_str();
    let (start, end) = (call.variant.start as i32, call.variant.stop as i32);
    let arrangement = &call.arrangement;
    let lengths: Vec<i32> = call.segments.iter().map(Interval::size).collect();
    let (missing, present, inverted) = segment_appearances(call)?;
    let mut result = Vec::new();
    if let Some(index) = find_all_segments(arrangement, call.segments.len()) {
        if index != 0 {
            let front = anchor(reference, contig, start - 1)?;
            result.extend(insertion_from_one_end(
                true,
                index,
                contig,
                start - 1,
                front,
                &lengths,
                arrangement,
            )?);
        }
        let last = index + call.segments.len() - 1;
        if last < arrangement.len() - 1 {
            let rear = anchor(reference, contig, end)?;
            result.extend(insertion_from_one_end(
                false,
                last + 1,
                contig,
                end,
                rear,
                &lengths,
                arrangement,
            )?);
        }
    } else {
        for number in &inverted {
            let segment = &call.segments[*number as usize - 1];
            if segment.size() > EVENT_SIZE_THRESHOLD && !present.contains(number) {
                let bases = anchor(reference, contig, segment.start)?;
                result.push(make_inversion(segment, bases)?);
            }
        }
        if !missing.is_empty() {
            for gone in compactify_missing(&missing) {
                if gone.size() > EVENT_SIZE_THRESHOLD {
                    let bases = anchor(reference, &gone.contig, gone.start)?;
                    result.push(make_deletion(
                        &Interval::new(&gone.contig, gone.start, gone.end - 1)?,
                        bases,
                    )?);
                }
            }
        }
        let leading = arrangement
            .iter()
            .take_while(|text| indicates_insertion(text))
            .count();
        if leading > 0 {
            let front = anchor(reference, contig, start)?;
            result.extend(insertion_from_one_end(
                true,
                leading,
                contig,
                start,
                front,
                &lengths,
                arrangement,
            )?);
        }
        let trailing = arrangement
            .iter()
            .rev()
            .take_while(|text| indicates_insertion(text))
            .count();
        if trailing > 0 {
            let index = arrangement.len() - trailing;
            let rear = anchor(reference, contig, end)?;
            result.extend(insertion_from_one_end(
                false,
                index,
                contig,
                end,
                rear,
                &lengths,
                arrangement,
            )?);
        }
    }
    Ok(result
        .into_iter()
        .map(|variant| with_source(variant, call))
        .collect())
}

/// The multi-segment complex calls' records from the segments alone.
pub fn multi_segment_variants(
    calls: &[ComplexCall],
    reference: &mut dyn Reference,
) -> Result<Vec<VariantContext>, Failed> {
    let mut variants = Vec::new();
    for call in calls.iter().filter(|call| call.segments.len() > 1) {
        variants.extend(extract_multi_segment(call, reference)?);
    }
    Ok(variants)
}

/// `reInterpretMultiSegmentComplexVarThroughAlignmentPairIteration`'s first half: the contigs of the
/// multi-segment calls, preprocessed anew, to the novel adjacencies of `parseOneContig`.
pub fn pair_iteration_adjacencies(
    calls: &[ComplexCall],
    records: &[BamRecord],
    header: &SamHeader,
    canonical: &[String],
) -> Result<Adjacencies, Failed> {
    let relevant: Vec<&String> = calls
        .iter()
        .filter(|call| call.segments.len() > 1)
        .flat_map(|call| call.contigs.iter())
        .collect();
    let chosen: Vec<BamRecord> = records
        .iter()
        .filter(|record| relevant.contains(&&record.read_name))
        .cloned()
        .collect();
    let classified = preprocess(&chosen, header, canonical)?;
    let contigs: Vec<&Contig> = classified
        .complex
        .iter()
        .map(|tig| &tig.contig)
        .filter(|contig| contig.alignments.len() > 1)
        .collect();
    adjacencies_of(&contigs, header)
}

fn attribute_text(variant: &VariantContext, key: &str) -> Option<String> {
    variant
        .attributes
        .iter()
        .find(|(name, _)| name == key)
        .and_then(|(_, value)| match value {
            Value::Str(text) => Some(text.clone()),
            Value::Int(number) => Some(number.to_string()),
            _ => None,
        })
}

fn attribute_int(variant: &VariantContext, key: &str) -> i32 {
    variant
        .attributes
        .iter()
        .find(|(name, _)| name == key)
        .and_then(|(_, value)| match value {
            Value::Int(number) => Some(*number as i32),
            Value::Str(text) => text.parse().ok(),
            _ => None,
        })
        .unwrap_or(0)
}

/// `deletionConsistencyCheck`.
fn deletion_is_consistent(variant: &VariantContext, missing: &[Interval]) -> Result<bool, Failed> {
    if missing.is_empty() {
        return Ok(false);
    }
    let deleted = Interval::new(
        &variant.contig,
        variant.start as i32 + 1,
        variant.stop as i32,
    )?;
    for segment in missing {
        if !segment.overlaps(&deleted) {
            return Ok(false);
        }
        if (segment.size() - deleted.size()).abs() > 2 {
            return Ok(false);
        }
        // The half-open `SVInterval`s, on one dummy contig.
        let overlap =
            (segment.end.min(deleted.end) - (segment.start - 1).max(deleted.start - 1)).max(0);
        if (segment.size().min(deleted.size()) - overlap).abs() <= 2 {
            return Ok(true);
        }
    }
    Ok(false)
}

/// `isConsistentWithCPX`.
fn consistent_with_complex(variant: &VariantContext, call: &ComplexCall) -> Result<bool, Failed> {
    match attribute_text(variant, "SVTYPE").as_deref() {
        Some("DEL") => deletion_is_consistent(variant, &segment_appearances(call)?.0),
        Some("INV") => Ok(false),
        _ => Ok(true),
    }
}

/// `postProcessConvertReplacementToFatInsOrInsAndDel`.
fn convert_replacement(
    simple: VariantContext,
    reference: &mut dyn Reference,
) -> Result<Vec<VariantContext>, Failed> {
    if attribute_text(&simple, "SVTYPE").as_deref() != Some("DEL") {
        return Ok(vec![simple]);
    }
    let deletion_length = -attribute_int(&simple, "SVLEN");
    let insertion_length = attribute_int(&simple, "INSLEN");
    let (start, stop) = (simple.start as i32, simple.stop as i32);
    let contig = simple.contig.clone();
    let adopt = |target: &mut VariantContext, without: &[&str]| {
        for (key, value) in &simple.attributes {
            if !without.contains(&key.as_str()) {
                set(target, key, value.clone());
            }
        }
    };
    if insertion_length > EVENT_SIZE_THRESHOLD && deletion_length > EVENT_SIZE_THRESHOLD {
        let dropped = ["INSSEQ_MAP", "SVLEN", "SVTYPE"];
        let reference_allele = simple.alleles[0].clone();
        let mut insertion = make_insertion(
            &contig,
            start,
            start,
            insertion_length,
            reference_allele.clone(),
        )?;
        adopt(&mut insertion, &dropped);
        remove(&mut insertion, "HOMSEQ");
        remove(&mut insertion, "HOMLEN");
        set(&mut insertion, "END", Value::Int(i64::from(start)));
        let mut deletion = make_deletion(&Interval::new(&contig, start, stop)?, reference_allele)?;
        adopt(&mut deletion, &dropped);
        remove(&mut deletion, "INSSEQ");
        remove(&mut deletion, "INSLEN");
        remove(&mut deletion, "SEQ_ALT_HAPLOTYPE");
        set(
            &mut insertion,
            "LINK",
            Value::Str(derived_id("DEL", &contig, start, stop)),
        );
        set(
            &mut deletion,
            "LINK",
            Value::Str(derived_id("INS", &contig, start, start)),
        );
        Ok(vec![deletion, insertion])
    } else if insertion_length > EVENT_SIZE_THRESHOLD {
        let bases = reference.bases(&contig, start, stop)?;
        let mut fat = make_insertion(
            &contig,
            start,
            stop,
            insertion_length,
            allele(&bases, true)?,
        )?;
        adopt(
            &mut fat,
            &["INSSEQ_MAP", "HOMLEN", "HOMSEQ", "SVLEN", "SVTYPE"],
        );
        fat.id = simple.id.replace("DEL", "INS");
        Ok(vec![fat])
    } else if deletion_length > EVENT_SIZE_THRESHOLD {
        Ok(vec![simple])
    } else {
        Ok(Vec::new())
    }
}

/// `reInterpretMultiSegmentComplexVarThroughAlignmentPairIteration`'s second half.
pub fn pair_iteration_variants(
    adjacencies: &Adjacencies,
    calls: &[ComplexCall],
    header: &SamHeader,
    reference: &mut dyn Reference,
    cnv_calls: Option<&[CnvCall]>,
) -> Result<Vec<VariantContext>, Failed> {
    let mut interpreted = Vec::with_capacity(adjacencies.0.len());
    for (adjacency, evidence) in &adjacencies.0 {
        interpreted.push(annotated_variant(
            adjacency, evidence, reference, header, cnv_calls,
        )?);
    }
    let multi: Vec<&ComplexCall> = calls
        .iter()
        .filter(|call| call.segments.len() > 1)
        .collect();
    let mut result = Vec::with_capacity(interpreted.len());
    for mut variant in interpreted {
        let names = attribute_text(&variant, "CTG_NAMES").unwrap_or_default();
        let mut consistent = Vec::new();
        for name in names.split(',') {
            // `collectAsMap`: a name two calls share is the later call's.
            if let Some(call) = multi
                .iter()
                .rfind(|call| call.contigs.iter().any(|contig| contig == name))
            {
                if consistent_with_complex(&variant, call)? {
                    consistent.push(call.variant.id.clone());
                }
            }
        }
        if consistent.is_empty() {
            continue;
        }
        variant.id.push_str("_CPX_DERIVED");
        set(&mut variant, "CPX_EVENT", Value::Str(consistent.join(",")));
        // `postProcessConvertShortDupToIns`.
        if attribute_text(&variant, "SVTYPE").as_deref() == Some("DUP") {
            let unit = attribute_text(&variant, "DUP_REPEAT_UNIT_REF_SPAN").unwrap_or_default();
            if parse_interval(&unit)?.size() <= EVENT_SIZE_THRESHOLD {
                variant.alleles[1] = symbolic("<INS>")?;
                set(&mut variant, "SVTYPE", Value::Str("INS".to_string()));
            }
        }
        result.extend(convert_replacement(variant, reference)?);
    }
    Ok(result)
}

/// `AnnotatedInterval`.
#[derive(PartialEq)]
struct AnnotatedKey {
    interval: Interval,
    id: String,
    kind: String,
    length: i32,
    alleles: Vec<Allele>,
}

impl AnnotatedKey {
    fn of(variant: &VariantContext) -> Self {
        Self {
            interval: Interval {
                contig: variant.contig.clone(),
                start: variant.start as i32,
                end: variant.stop as i32,
            },
            id: variant.id.clone(),
            kind: attribute_text(variant, "SVTYPE").unwrap_or_default(),
            length: attribute_int(variant, "SVLEN"),
            alleles: variant.alleles.clone(),
        }
    }

    fn hash_code(&self) -> i32 {
        let mut result = self.interval.hash_code();
        result = fold(result, string_hash_code(&self.id));
        result = fold(result, string_hash_code(&self.kind));
        result = fold(result, self.length);
        // `SimpleAllele.hashCode`.
        fold(
            result,
            list_hash(self.alleles.iter().map(|allele| {
                fold(
                    byte_array_hash_code(allele.base_string().as_bytes()),
                    if allele.is_reference() { 1231 } else { 1237 },
                )
            })),
        )
    }
}

/// `getAnnotatedIntervalToSourceCpxIDsAndContigNames`: the `HashSet` that keeps a key's first
/// record, then the `HashMap` filled from it, in that map's order.
fn annotated_order(
    variants: Vec<VariantContext>,
) -> Result<Vec<(AnnotatedKey, VariantContext)>, Failed> {
    let mut distinct: Vec<(AnnotatedKey, VariantContext)> = Vec::with_capacity(variants.len());
    for variant in variants {
        let key = AnnotatedKey::of(&variant);
        if !distinct.iter().any(|(other, _)| *other == key) {
            distinct.push((key, variant));
        }
    }
    let hashes: Vec<i32> = distinct.iter().map(|(key, _)| key.hash_code()).collect();
    let in_set: Vec<(usize, i32)> = hashes.iter().copied().enumerate().collect();
    let refused = |error: gatk_engine::java_hash::HashOrderError| format!("{error:?}");
    let set_order = hash_map_order(&in_set).map_err(refused)?;
    let in_map: Vec<(usize, i32)> = set_order.iter().map(|at| (*at, hashes[*at])).collect();
    let map_order = hash_map_order(&in_map).map_err(refused)?;
    let mut taken: Vec<Option<(AnnotatedKey, VariantContext)>> =
        distinct.into_iter().map(Some).collect();
    Ok(map_order
        .into_iter()
        .map(|at| taken[at].take().expect("each index once"))
        .collect())
}

/// `removeDuplicates`.
pub fn remove_duplicates(
    with_less: Vec<VariantContext>,
    with_more: Vec<VariantContext>,
) -> Result<Vec<VariantContext>, Failed> {
    let mut less = annotated_order(with_less)?;
    let more = annotated_order(with_more)?;
    let mut result = Vec::with_capacity(less.len() + more.len());
    let union = |one: &VariantContext, two: &VariantContext, key: &str| {
        let mut all: Vec<String> = [one, two]
            .iter()
            .filter_map(|variant| attribute_text(variant, key))
            .flat_map(|text| text.split(',').map(str::to_string).collect::<Vec<String>>())
            .collect();
        all.sort();
        all.dedup();
        all.join(",")
    };
    for (key, mut variant) in more {
        if let Some(at) = less.iter().position(|(other, _)| *other == key) {
            let (_, other) = less.remove(at);
            let events = union(&variant, &other, "CPX_EVENT");
            let names = union(&variant, &other, "CTG_NAMES");
            set(&mut variant, "CPX_EVENT", Value::Str(events));
            set(&mut variant, "CTG_NAMES", Value::Str(names));
        }
        result.push(variant);
    }
    result.extend(less.into_iter().map(|(_, variant)| variant));
    Ok(result)
}

/// `AnnotatedVariantProducer.filterMergedVariantList`.
pub fn filter_merged(
    variants: Vec<VariantContext>,
    minimum_mq: i32,
    minimum_length: i32,
) -> Vec<VariantContext> {
    let mut kept = Vec::with_capacity(variants.len());
    for mut variant in variants {
        let kind = attribute_text(&variant, "SVTYPE").unwrap_or_default();
        if matches!(kind.as_str(), "DEL" | "INS" | "DUP")
            && attribute_int(&variant, "SVLEN").abs() < SIZE_LOWER_BOUND
        {
            continue;
        }
        if attribute_text(&variant, "CTG_NAMES").is_some() {
            let mut failed = Vec::new();
            let qualities = attribute_text(&variant, "MAPPING_QUALITIES").unwrap_or_default();
            let maximum = qualities
                .split(',')
                .filter_map(|text| text.parse::<i32>().ok())
                .fold(0, i32::max);
            if maximum < minimum_mq {
                failed.push("LOW_MQ".to_string());
            }
            if attribute_int(&variant, "MAX_ALIGN_LENGTH") < minimum_length {
                failed.push("SHORT_ALN".to_string());
            }
            if !failed.is_empty() {
                variant.filters = Some(failed);
            }
        }
        kept.push(variant);
    }
    kept
}

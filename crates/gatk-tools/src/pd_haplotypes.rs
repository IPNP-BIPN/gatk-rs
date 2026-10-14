//! The partially determined haplotypes, ported from
//! `tools.walkers.haplotypecaller.PartiallyDeterminedHaplotypeComputationEngine`,
//! `utils.haplotype.PartiallyDeterminedHaplotype` and `utils.SmallBitSet` (GATK 4.6.2.0).
//!
//! DRAGEN's alternative to the assembled haplotypes: every event the assembly found, less the bad
//! pileup events and with the good ones, is grouped with the events it cannot share a haplotype
//! with (overlapping ones, and pairs or trios that Smith-Waterman realigns into another of the
//! events). Then each locus with events is in turn the determined one, and each of its alleles,
//! the reference included, the determined allele: one haplotype per branch of compatible events
//! elsewhere, each with that allele fixed and every other event left undetermined, which the
//! PDHMM scores as either.
//!
//! # Orders
//!
//! The events are sorted by `HAPLOTYPE_SNP_FIRST_COMPARATOR` (start, reference length, alternate
//! length, alternate allele). The groups come out of JGraphT's connectivity inspector in the order
//! of each one's first event. The haplotypes are a `LinkedHashSet`, sorted by their bases at the
//! end; ties keep the order they were made in.
//!
//! # Giving up
//!
//! A group of more than seventeen events, more than 128 branches at a locus, or more than 512
//! haplotypes (or 128 fully determined ones from a branch) abandons the construction and leaves the
//! assembled haplotypes as they were, without marking the set as overwritten.

use gatk_engine::assembly_result_set::AssemblyResultSet;
use gatk_engine::cigar_builder::CigarBuilder;
use gatk_engine::cigar_utils::calculate_cigar;
use gatk_engine::event_map::{build_event_maps_for_haplotypes, Event};
use gatk_engine::haplotype::{pd_bytes, Haplotype, PartiallyDetermined};
use gatk_engine::interval::SimpleInterval;
use gatk_engine::smith_waterman::{SmithWatermanAligner, SwOverhangStrategy, SwParameters};
use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_vcf::allele::Allele;
use std::cmp::Ordering;

/// `MAX_PD_HAPS_TO_GENERATE`.
const MAX_PD_HAPS_TO_GENERATE: usize = 256 * 2;
/// `MAX_BRANCH_PD_HAPS`.
const MAX_BRANCH_PD_HAPS: usize = 128;
/// `MAX_VAR_IN_EVENT_GROUP`.
const MAX_VAR_IN_EVENT_GROUP: usize = 17;

/// What `generatePDHaplotypes` reads of the arguments.
#[derive(Debug, Clone, PartialEq)]
pub struct PdArguments {
    /// `--make-determined-haps-from-pd-code`.
    pub use_determined_haplotypes: bool,
    /// `--pileup-detection-snp-adjacent-to-assembled-indel-range`.
    pub snp_adjacent_to_assembly_indel: i32,
}

type Failure = String;

fn illegal_state(message: String) -> Failure {
    format!("java.lang.IllegalStateException: {message}")
}

/// `Allele.compareTo` between two alternate alleles: their bases as strings.
fn compare_alleles(a: &Allele, b: &Allele) -> Ordering {
    b.is_reference()
        .cmp(&a.is_reference())
        .then_with(|| a.display_string().cmp(&b.display_string()))
}

/// `HAPLOTYPE_SNP_FIRST_COMPARATOR`.
pub fn snp_first(a: &Event, b: &Event) -> Ordering {
    a.start()
        .cmp(&b.start())
        .then(a.ref_allele().len().cmp(&b.ref_allele().len()))
        .then(a.alt_allele().len().cmp(&b.alt_allele().len()))
        .then_with(|| compare_alleles(a.alt_allele(), b.alt_allele()))
}

/// `dragenStart`.
fn dragen_start(event: &Event) -> f64 {
    f64::from(event.start())
        + if event.is_indel() {
            if event.is_simple_deletion() {
                1.0
            } else {
                0.5
            }
        } else {
            0.0
        }
}

/// `eventsOverlapForPDHapsCode`.
fn events_overlap(e1: &Event, e2: &Event) -> bool {
    if e1.contig() != e2.contig() {
        return false;
    }
    let end = |e: &Event| f64::from(e.end()) + if e.is_simple_insertion() { 0.5 } else { 0.0 };
    !(dragen_start(e1) > end(e2) || dragen_start(e2) > end(e1))
}

/// `Locatable.withinDistanceOf`.
fn within_distance_of(event: &Event, other: &Event, distance: i32) -> bool {
    event.contig() == other.contig()
        && event.start() <= other.end() + distance
        && other.start() - distance <= event.end()
}

fn bases(allele: &Allele) -> Vec<u8> {
    allele.display_string().into_bytes()
}

/// `getPDBytesForHaplotypes(refAllele, altAllele)`.
fn pd_bytes_for(reference: &Allele, alt: &Allele) -> Result<Vec<u8>, Failure> {
    let length = if alt.len() == reference.len() {
        reference.len()
    } else {
        reference.len() - 1
    };
    let mut output = vec![0u8; length];
    if alt.len() == reference.len() {
        let flag = match bases(alt)[0] {
            b'A' => pd_bytes::A,
            b'C' => pd_bytes::C,
            b'T' => pd_bytes::T,
            b'G' => pd_bytes::G,
            _ => {
                return Err(
                    "java.lang.RuntimeException: Found unexpected base in alt alleles".to_string(),
                )
            }
        };
        output[0] = output[0].wrapping_add(pd_bytes::SNP).wrapping_add(flag);
    } else {
        output[0] = output[0].wrapping_add(pd_bytes::DEL_START);
        let last = output.len() - 1;
        output[last] = output[last].wrapping_add(pd_bytes::DEL_END);
    }
    Ok(output)
}

/// `actualStartExcludingInitialIndelBase`.
fn actual_start(event: &Event) -> i32 {
    event.start() + i32::from(event.is_indel())
}

fn make_cigar(builder: &mut CigarBuilder) -> Result<Cigar, Failure> {
    builder.make(false).map_err(|e| format!("{e:?}"))
}

fn add(builder: &mut CigarBuilder, length: i32, op: Op) -> Result<(), Failure> {
    builder
        .add(CigarElement {
            length: length as u32,
            op,
        })
        .map_err(|e| format!("{e:?}"))
}

fn reference_start(reference: &Haplotype) -> Result<i32, Failure> {
    reference
        .genome_location()
        .map(|l| l.start)
        .ok_or_else(|| "java.lang.NullPointerException".to_string())
}

fn validate_reference(reference: &Haplotype) -> Result<(), Failure> {
    let single = reference.cigar().is_some_and(|c| c.elements.len() == 1);
    if !(reference.is_reference() && single) {
        return Err(illegal_state(
            "This is not a valid base haplotype for construction".to_string(),
        ));
    }
    Ok(())
}

/// `constructHaplotypeFromEvents(refHap, events, setEventMap)`.
pub fn construct_haplotype_from_events(
    reference: &Haplotype,
    events: &[Event],
    set_event_map: bool,
) -> Result<Haplotype, Failure> {
    validate_reference(reference)?;
    let location = reference
        .genome_location()
        .cloned()
        .ok_or_else(|| "java.lang.NullPointerException".to_string())?;
    for event in events {
        if !(event.contig() == location.contig
            && location.start <= event.start()
            && event.end() <= location.end)
        {
            return Err(illegal_state(format!(
                "Provided Variant Context{event:?}doesn't overlap haplotype"
            )));
        }
    }
    let ref_start = location.start;
    let mut last = ref_start;
    let ref_bases = reference.bases();
    let mut cigar = CigarBuilder::new(true);
    for pair in events.windows(2) {
        if pair[0].end() >= actual_start(&pair[1]) {
            return Err(illegal_state(format!(
                "PD event list: {events:?} is out of order."
            )));
        }
    }
    let mut out: Vec<u8> = Vec::new();
    for event in events {
        let start = actual_start(event);
        add(&mut cigar, start - last, Op::M)?;
        let difference = event.alt_allele().len() as i32 - event.ref_allele().len() as i32;
        if difference == 0 {
            add(&mut cigar, event.ref_allele().len() as i32, Op::X)?;
        } else {
            add(
                &mut cigar,
                difference.abs(),
                if difference < 0 { Op::D } else { Op::I },
            )?;
        }
        out.extend_from_slice(
            &ref_bases[(last - ref_start) as usize..(start - ref_start) as usize],
        );
        let alt = bases(event.alt_allele());
        out.extend_from_slice(if difference == 0 { &alt } else { &alt[1..] });
        last = event.end() + 1;
    }
    let index = (last - ref_start) as usize;
    out.extend_from_slice(&ref_bases[index..]);
    add(&mut cigar, (ref_bases.len() - index) as i32, Op::M)?;
    let mut result = Haplotype::non_reference(&out).map_err(|e| format!("{e:?}"))?;
    result.set_genome_location(location.clone());
    result
        .set_cigar(&make_cigar(&mut cigar)?)
        .map_err(|e| format!("{e:?}"))?;
    if set_event_map {
        let mut one = vec![result];
        build_event_maps_for_haplotypes(&mut one, &ref_bases, &location, 0)
            .map_err(|e| e.message())?;
        result = one.pop().expect("one haplotype");
        result.set_alignment_start_hap_wrt_ref(reference.alignment_start_hap_wrt_ref());
    }
    Ok(result)
}

/// `constructArtificialHaplotypeAndTestEquivalentEvents`: whether the events, applied to the
/// reference and realigned, give the reference back or an event of the list not among them.
fn equivalent_events(
    reference: &Haplotype,
    aligner: &dyn SmithWatermanAligner,
    sw: &SwParameters,
    all: &[Event],
    tested: &[Event],
) -> Result<bool, Failure> {
    let mut realigned = construct_haplotype_from_events(reference, tested, false)?;
    let ref_bases = reference.bases();
    if realigned.bases() == ref_bases {
        return Ok(true);
    }
    let cigar = calculate_cigar(
        &ref_bases,
        &realigned.bases(),
        aligner,
        sw,
        SwOverhangStrategy::Indel,
    )
    .map_err(|e| format!("{e:?}"))?
    .ok_or_else(|| "java.lang.NullPointerException".to_string())?;
    realigned.set_cigar(&cigar).map_err(|e| format!("{e:?}"))?;
    let location = reference
        .genome_location()
        .cloned()
        .ok_or_else(|| "java.lang.NullPointerException".to_string())?;
    let mut one = vec![realigned];
    build_event_maps_for_haplotypes(&mut one, &ref_bases, &location, 0).map_err(|e| e.message())?;
    let map = one[0].event_map().expect("just built");
    let equivalent = map
        .events()
        .filter(|event| !tested.contains(event))
        .any(|event| all.contains(event));
    Ok(equivalent)
}

/// `smithWatermanRealignPairsOfVariantsForEquivalentEvents`: the pairs, then the trios, that
/// realign into another event of the list.
fn disallowed_combinations(
    reference: &Haplotype,
    aligner: &dyn SmithWatermanAligner,
    sw: &SwParameters,
    events: &[Event],
) -> Result<Vec<Vec<Event>>, Failure> {
    let mut disallowed: Vec<Vec<Event>> = Vec::new();
    let sorted = |mut list: Vec<Event>| {
        list.sort_by(snp_first);
        list
    };
    for (i, first) in events.iter().enumerate() {
        if !first.is_indel() {
            continue;
        }
        for (j, second) in events.iter().enumerate() {
            if j != i && !events_overlap(first, second) && (!second.is_indel() || j > i) {
                let pair = sorted(vec![first.clone(), second.clone()]);
                if equivalent_events(reference, aligner, sw, events, &pair)? {
                    disallowed.push(pair);
                }
            }
        }
    }
    let contains_both = |list: &[Vec<Event>], a: &Event, b: &Event| {
        list.iter().any(|p| p.contains(a) && p.contains(b))
    };
    for (i, first) in events.iter().enumerate() {
        if !first.is_indel() {
            continue;
        }
        for (j, second) in events.iter().enumerate() {
            if !(j != i && !events_overlap(first, second) && (!second.is_indel() || j > i)) {
                continue;
            }
            if contains_both(&disallowed, first, second) {
                continue;
            }
            for (k, third) in events.iter().enumerate().skip(j + 1) {
                if k != i && !events_overlap(third, first) && !events_overlap(third, second) {
                    if contains_both(&disallowed, first, third)
                        || contains_both(&disallowed, second, third)
                    {
                        continue;
                    }
                    let trio = sorted(vec![first.clone(), second.clone(), third.clone()]);
                    if equivalent_events(reference, aligner, sw, events, &trio)? {
                        disallowed.push(trio);
                    }
                }
            }
        }
    }
    Ok(disallowed)
}

/// `EventGroup`: a connected set of mutually exclusive events and which subsets of it may share a
/// haplotype.
struct EventGroup {
    events: Vec<Event>,
    /// `allowedSubsets`, indexed by the subset's bits.
    allowed: Vec<bool>,
    cached: Option<Vec<Vec<Event>>>,
}

impl EventGroup {
    fn new(component: Vec<Event>, mutexes: &[Vec<Event>]) -> Result<EventGroup, Failure> {
        let mut events = component;
        events.sort_by(snp_first);
        let index = |e: &Event| events.iter().position(|x| x == e);
        let overlapping: Vec<&Vec<Event>> = mutexes
            .iter()
            .filter(|m| {
                !(m.len() == 2 && m[0].start() == m[1].start() && m[0].is_snp() && m[1].is_snp())
            })
            .filter(|m| m.iter().any(|e| index(e).is_some()))
            .collect();
        for mutex in &overlapping {
            if !mutex.iter().all(|e| index(e).is_some()) {
                return Err(illegal_state(format!(
                    "Mutex group {mutex:?} only partially overlaps event group"
                )));
            }
        }
        let size = events.len();
        let mut allowed = vec![false; 1 << size];
        if size >= 2 {
            for slot in allowed.iter_mut().skip(1) {
                *slot = true;
            }
            let bits: Vec<usize> = overlapping
                .iter()
                .map(|m| {
                    m.iter()
                        .map(|e| 1usize << index(e).expect("checked"))
                        .fold(0, |a, b| a | b)
                })
                .collect();
            if !bits.is_empty() {
                for subset in 1..(1usize << size) {
                    if bits.iter().any(|&m| subset & m == m) {
                        allowed[subset] = false;
                    }
                }
            }
        }
        Ok(EventGroup {
            events,
            allowed,
            cached: None,
        })
    }

    fn overlap_set(&self, events: &[Event]) -> usize {
        events
            .iter()
            .filter_map(|e| self.events.iter().position(|x| x == e))
            .fold(0, |a, i| a | (1 << i))
    }

    /// `eventSetsForPDHaplotypes(determinedEvents, locusEvents)`.
    fn event_sets(&mut self, determined: &[Event], locus: &[Event]) -> Vec<Vec<Event>> {
        let locus_set = self.overlap_set(locus);
        let determined_set = self.overlap_set(determined);
        if self.events.len() == 1 {
            return if locus_set != 0 && determined.is_empty() {
                vec![Vec::new()]
            } else {
                vec![vec![self.events[0].clone()]]
            };
        }
        let cached_case = locus_set == 0 && determined_set == 0;
        if cached_case {
            if let Some(cached) = &self.cached {
                return cached.clone();
            }
        }
        let mut maximal: Vec<usize> = Vec::new();
        let mut subset = (1usize << self.events.len()) - 1;
        while subset != 0 {
            if self.allowed[subset]
                && subset & locus_set == determined_set
                && !maximal.iter().any(|&group| group & subset == subset)
            {
                maximal.push(subset);
            }
            subset -= 1;
        }
        let output: Vec<Vec<Event>> = maximal
            .iter()
            .map(|&bits| {
                (0..self.events.len())
                    .filter(|i| bits & (1 << i) != 0)
                    .map(|i| self.events[i].clone())
                    .collect()
            })
            .collect();
        if cached_case {
            self.cached = Some(output.clone());
        }
        output
    }
}

/// `getEventGroupClusters`: none when a group is too large.
fn event_groups(
    events: &[Event],
    disallowed: &[Vec<Event>],
) -> Result<Option<Vec<EventGroup>>, Failure> {
    let mut mutexes: Vec<Vec<Event>> = disallowed.to_vec();
    for (e1, first) in events.iter().enumerate() {
        for second in events.iter().skip(e1 + 1) {
            if second.start() > first.end() + 1 {
                break;
            }
            if events_overlap(first, second) {
                mutexes.push(vec![first.clone(), second.clone()]);
            }
        }
    }
    // The connected components, in the order of their first event.
    let mut parent: Vec<usize> = (0..events.len()).collect();
    fn find(parent: &mut [usize], i: usize) -> usize {
        let mut root = i;
        while parent[root] != root {
            root = parent[root];
        }
        let mut node = i;
        while parent[node] != root {
            let next = parent[node];
            parent[node] = root;
            node = next;
        }
        root
    }
    let position = |e: &Event| {
        events
            .iter()
            .position(|x| x == e)
            .expect("an event of the list")
    };
    for mutex in &mutexes {
        for pair in mutex.windows(2) {
            let (a, b) = (
                find(&mut parent, position(&pair[0])),
                find(&mut parent, position(&pair[1])),
            );
            if a != b {
                let (low, high) = if a < b { (a, b) } else { (b, a) };
                parent[high] = low;
            }
        }
    }
    let mut components: Vec<(usize, Vec<Event>)> = Vec::new();
    for (i, event) in events.iter().enumerate() {
        let root = find(&mut parent, i);
        match components.iter_mut().find(|(r, _)| *r == root) {
            Some((_, list)) => list.push(event.clone()),
            None => components.push((root, vec![event.clone()])),
        }
    }
    if components
        .iter()
        .any(|(_, c)| c.len() > MAX_VAR_IN_EVENT_GROUP)
    {
        return Ok(None);
    }
    components
        .into_iter()
        .map(|(_, c)| EventGroup::new(c, &mutexes))
        .collect::<Result<Vec<_>, _>>()
        .map(Some)
}

/// `computeBranches`: none past the branch limit.
fn compute_branches(
    groups: &mut [EventGroup],
    determined: &[Event],
    locus: &[Event],
) -> Option<Vec<Vec<Event>>> {
    let union = |a: &[Event], b: &[Event]| {
        let mut out = a.to_vec();
        for e in b {
            if !out.contains(e) {
                out.push(e.clone());
            }
        }
        out
    };
    let mut branches: Vec<Vec<Event>> = vec![Vec::new()];
    for group in groups.iter_mut() {
        let sets = group.event_sets(determined, locus);
        let extra: Vec<Vec<Event>> = if sets.len() < 2 {
            Vec::new()
        } else {
            branches
                .iter()
                .flat_map(|branch| sets.iter().skip(1).map(|set| union(branch, set)))
                .collect()
        };
        if let Some(first) = sets.first() {
            for branch in branches.iter_mut() {
                *branch = union(branch, first);
            }
        }
        branches.extend(extra);
        if branches.len() > MAX_BRANCH_PD_HAPS {
            return None;
        }
    }
    Some(branches)
}

/// `createNewPDHaplotypeFromEvents(refHap, determinedEvents, determinedLocus, constituentEvents,
/// allEventsAtDeterminedLocus)`.
fn new_pd_haplotype(
    reference: &Haplotype,
    determined: &[Event],
    determined_locus: i32,
    constituent: &[Event],
    locus_events: &[Event],
) -> Result<Haplotype, Failure> {
    validate_reference(reference)?;
    let ref_is_determined = determined.is_empty();
    let ref_start = reference_start(reference)?;
    let mut last = ref_start;
    let ref_bases = reference.bases();
    let mut cigar = CigarBuilder::new(false);
    let mut out: Vec<u8> = Vec::new();
    let mut pd: Vec<u8> = Vec::new();
    let mut last_was_snp = false;
    for event in constituent {
        let start = actual_start(event);
        let before = start - last;
        if before == -1 && event.is_snp() && last_was_snp {
            let this = pd_bytes_for(event.ref_allele(), event.alt_allele())?[0];
            let at = pd.len() - 1;
            pd[at] |= this;
            continue;
        } else if event.start() == determined_locus && ref_is_determined {
            continue;
        }
        if before < 0 {
            return Err(illegal_state(format!(
                "Event {event:?} is out of order in PD event list: {constituent:?}."
            )));
        }
        let reference_allele = event.ref_allele();
        let alt_allele = event.alt_allele();
        let difference = alt_allele.len() as i32 - reference_allele.len() as i32;
        let is_insertion = difference > 0;
        let is_event = determined.contains(event);
        add(&mut cigar, before, Op::M)?;
        if event.is_snp() {
            add(
                &mut cigar,
                reference_allele.len() as i32,
                if ref_is_determined || !is_event {
                    Op::M
                } else {
                    Op::X
                },
            )?;
        } else if is_event {
            let length = if is_insertion && ref_is_determined {
                0
            } else {
                reference_allele.len().max(alt_allele.len()) as i32 - 1
            };
            let op = if is_insertion {
                Op::I
            } else if ref_is_determined {
                Op::M
            } else {
                Op::D
            };
            add(&mut cigar, length, op)?;
        } else {
            add(
                &mut cigar,
                difference.abs(),
                if difference > 0 { Op::I } else { Op::M },
            )?;
        }
        out.extend_from_slice(
            &ref_bases[(last - ref_start) as usize..(start - ref_start) as usize],
        );
        pd.extend(std::iter::repeat_n(0u8, (start - last) as usize));
        let use_ref = (is_event && ref_is_determined) || (!is_event && difference <= 0);
        let (to_use, other) = if use_ref {
            (reference_allele, alt_allele)
        } else {
            (alt_allele, reference_allele)
        };
        let all = bases(to_use);
        let added: &[u8] = if event.is_indel() { &all[1..] } else { &all };
        out.extend_from_slice(added);
        if is_event {
            pd.extend(std::iter::repeat_n(0u8, added.len()));
        } else {
            pd.extend(pd_bytes_for(to_use, other)?);
        }
        last = event.end() + 1;
        last_was_snp = event.is_snp();
    }
    let index = (last - ref_start) as usize;
    out.extend_from_slice(&ref_bases[index..]);
    pd.extend(std::iter::repeat_n(0u8, ref_bases.len() - index));
    add(&mut cigar, (ref_bases.len() - index) as i32, Op::M)?;
    let made = make_cigar(&mut cigar)?;

    let mut haplotype = Haplotype::non_reference(&out).map_err(|e| format!("{e:?}"))?;
    let location = reference
        .genome_location()
        .cloned()
        .ok_or_else(|| "java.lang.NullPointerException".to_string())?;
    if out.len() != pd.len() {
        return Err(
            "java.lang.IllegalArgumentException: pdBytes array must have same length as base haplotype."
                .to_string(),
        );
    }
    haplotype.set_cigar(&made).map_err(|e| format!("{e:?}"))?;
    haplotype.set_genome_location(location.clone());
    haplotype.set_alignment_start_hap_wrt_ref(reference.alignment_start_hap_wrt_ref());
    let start = locus_events
        .iter()
        .map(Event::start)
        .min()
        .unwrap_or(determined_locus);
    let end = locus_events
        .iter()
        .map(Event::end)
        .max()
        .unwrap_or(determined_locus);
    haplotype.set_partially_determined(PartiallyDetermined {
        alternate_bases: pd,
        constituent_events: constituent.to_vec(),
        determined_events: determined.to_vec(),
        determined_position: determined_locus,
        determined_extent: SimpleInterval::new(&location.contig, start, end)
            .ok_or_else(|| illegal_state("an empty determined extent".to_string()))?,
    });
    Ok(haplotype)
}

/// `makeFinalListOfEventsInOrder`: the assembled events less the bad pileup events, the good ones
/// added unless a SNP within range of an assembled indel, in `HAPLOTYPE_SNP_FIRST_COMPARATOR` order.
fn final_events(
    assembled: &[Event],
    bad: &[Event],
    good: &[Event],
    arguments: &PdArguments,
) -> Vec<Event> {
    let mut passing: Vec<Event> = assembled
        .iter()
        .filter(|e| !bad.contains(e))
        .cloned()
        .collect();
    let indels: Vec<Event> = passing.iter().filter(|e| e.is_indel()).cloned().collect();
    for event in good {
        if passing.contains(event) {
            continue;
        }
        if event.is_indel()
            || indels.iter().all(|indel| {
                !within_distance_of(event, indel, arguments.snp_adjacent_to_assembly_indel)
            })
        {
            passing.push(event.clone());
        }
    }
    passing.sort_by(snp_first);
    passing.dedup();
    passing
}

/// `generatePDHaplotypes(sourceSet, badPileupEvents, goodPileupEvents, aligner, args)`, in place.
pub fn generate_pd_haplotypes(
    set: &mut AssemblyResultSet,
    bad: &[Event],
    good: &[Event],
    aligner: &dyn SmithWatermanAligner,
    sw: &SwParameters,
    arguments: &PdArguments,
) -> Result<(), Failure> {
    let reference = set
        .reference_haplotype()
        .cloned()
        .ok_or_else(|| "java.lang.NullPointerException".to_string())?;
    let calling_span = set
        .region_for_genotyping()
        .ok_or_else(|| "java.lang.NullPointerException".to_string())?
        .span()
        .clone();
    let assembled = set.variation_events(0).map_err(|e| format!("{e:?}"))?;
    let events = final_events(&assembled, bad, good, arguments);

    let disallowed = disallowed_combinations(&reference, aligner, sw, &events)?;
    let Some(mut groups) = event_groups(&events, &disallowed)? else {
        return Ok(());
    };

    let mut loci: Vec<(i32, Vec<Event>)> = Vec::new();
    for event in &events {
        match loci.iter_mut().find(|(start, _)| *start == event.start()) {
            Some((_, list)) => list.push(event.clone()),
            None => loci.push((event.start(), vec![event.clone()])),
        }
    }
    loci.sort_by_key(|(start, _)| *start);

    let determined_mode = arguments.use_determined_haplotypes;
    let mut output: Vec<Haplotype> = Vec::new();
    let push = |output: &mut Vec<Haplotype>, haplotype: Haplotype| {
        if !output.contains(&haplotype) {
            output.push(haplotype);
        }
    };
    if determined_mode {
        output.push(reference.clone());
    }
    for (locus, here) in &loci {
        if !(calling_span.start <= *locus && *locus <= calling_span.end) {
            continue;
        }
        let first: i64 = if determined_mode { 0 } else { -1 };
        for index in first..here.len() as i64 {
            let determined: Vec<Event> = if index < 0 {
                Vec::new()
            } else {
                vec![here[index as usize].clone()]
            };
            let Some(branches) = compute_branches(&mut groups, &determined, here) else {
                return Ok(());
            };
            for branch in &branches {
                if !determined_mode {
                    let mut in_order = branch.clone();
                    in_order.sort_by(snp_first);
                    let haplotype =
                        new_pd_haplotype(&reference, &determined, *locus, &in_order, here)?;
                    push(&mut output, haplotype);
                } else {
                    let mut full: Vec<Vec<Event>> = vec![determined.clone()];
                    for (other, at_other) in &loci {
                        if locus < other {
                            let mut children: Vec<Vec<Event>> = Vec::new();
                            for event in at_other.iter().filter(|e| branch.contains(e)) {
                                for group in &full {
                                    let mut child = group.clone();
                                    child.push(event.clone());
                                    child.sort_by(snp_first);
                                    children.push(child);
                                }
                            }
                            full.extend(children);
                            if full.len() > MAX_BRANCH_PD_HAPS {
                                return Ok(());
                            }
                        }
                    }
                    for events in &full {
                        let haplotype = construct_haplotype_from_events(&reference, events, true)?;
                        push(&mut output, haplotype);
                    }
                }
                if output.len() > MAX_PD_HAPS_TO_GENERATE {
                    return Ok(());
                }
            }
        }
    }
    set.store_assembly_haplotypes();
    output.sort_by_key(|h| h.bases());
    set.replace_all_haplotypes(&output)
        .map_err(|e| format!("{e:?}"))?;
    set.set_partially_determined_mode(!determined_mode);
    Ok(())
}

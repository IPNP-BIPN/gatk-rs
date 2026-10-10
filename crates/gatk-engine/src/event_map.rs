//! Ported from `org.broadinstitute.hellbender.utils.haplotype.Event` and `EventMap` (GATK 4.6.2.0):
//! the variants a haplotype carries, read off its CIGAR against the reference.
//!
//! What a reader would not guess:
//!
//!  * an insertion is only an event if it is neither the first nor the last CIGAR element and does
//!    not start the reference window, and a deletion only if it does not start the window; both
//!    are padded with the reference base before them;
//!  * mismatches inside one M, = or X element are grouped into an MNP while the gap from the last
//!    one is at most `maxMnpDistance`, so with 1 the substitutions at 10, 11, 12, 14, 15 and 17 are
//!    an MNP at 10-12, one at 14-15 and a SNP at 17;
//!  * "regular" is `BaseUtils.isRegularBase`: A, C, G and T in either case, **and the wildcard
//!    `*`**, which the base index table maps to A;
//!  * events at one start are compounded as they arrive (a SNP with an insertion or a deletion, an
//!    insertion with a deletion), and two SNPs there are refused;
//!  * every `Event` is cut to its minimal representation, dropping the trailing bases its two
//!    alleles share unless either is a single base. Only a compound of an insertion and a deletion
//!    can need it, and no canonical CIGAR makes one: the CIGAR builder puts a deletion before an
//!    adjacent insertion, which moves the insertion's start.

use std::collections::BTreeMap;

use htsjdk_bam::cigar::Op;
use htsjdk_vcf::allele::{Allele, AlleleError};

use crate::base_utils::simple_base_to_base_index;
use crate::haplotype::Haplotype;
use crate::interval::SimpleInterval;

/// What the event map refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EventError {
    /// `IllegalArgumentException` from `Utils.validateArg` or `ParamUtils`.
    IllegalArgument(&'static str),
    /// `GATKException("Unsupported cigar operator created during SW alignment: <op>")`.
    UnsupportedOperator(char),
    /// `NullPointerException`: a haplotype without a CIGAR.
    NoCigar,
    /// `Utils.nonNull(h.getEventMap(), "Haplotype event map has not been set")`.
    NoEventMap,
    /// `ArrayIndexOutOfBoundsException`: a CIGAR that walks past the reference or the haplotype.
    IndexOutOfBounds,
    /// The allele constructor's own refusal.
    Allele(AlleleError),
}

impl EventError {
    /// The simple name of the exception the reference throws.
    pub fn class(&self) -> &'static str {
        match self {
            EventError::IllegalArgument(_) | EventError::Allele(_) => "IllegalArgumentException",
            EventError::UnsupportedOperator(_) => "GATKException",
            EventError::NoCigar | EventError::NoEventMap => "NullPointerException",
            EventError::IndexOutOfBounds => "ArrayIndexOutOfBoundsException",
        }
    }

    /// Its message, where the reference gives one.
    pub fn message(&self) -> String {
        match self {
            EventError::IllegalArgument(message) => message.to_string(),
            EventError::UnsupportedOperator(op) => {
                format!("Unsupported cigar operator created during SW alignment: {op}")
            }
            EventError::NoCigar => "null".to_string(),
            EventError::NoEventMap => "Haplotype event map has not been set".to_string(),
            EventError::IndexOutOfBounds => "index out of bounds".to_string(),
            EventError::Allele(e) => e.to_string(),
        }
    }
}

impl From<AlleleError> for EventError {
    fn from(e: AlleleError) -> Self {
        EventError::Allele(e)
    }
}

/// `Event`: a reference and an alternate allele at a start, in minimal representation.
#[derive(Debug, Clone)]
pub struct Event {
    contig: String,
    start: i32,
    stop: i32,
    ref_allele: Allele,
    alt_allele: Allele,
}

impl Event {
    /// `new Event(contig, start, ref, alt)`.
    pub fn new(
        contig: &str,
        start: i32,
        reference: Allele,
        alt: Allele,
    ) -> Result<Event, EventError> {
        if !reference.is_reference() {
            return Err(EventError::IllegalArgument("ref is not ref"));
        }
        let (ref_allele, alt_allele) = make_minimal_representation(reference, alt)?;
        let stop = start + ref_allele.len() as i32 - 1;
        Ok(Event {
            contig: contig.to_string(),
            start,
            stop,
            ref_allele,
            alt_allele,
        })
    }

    pub fn contig(&self) -> &str {
        &self.contig
    }

    pub fn start(&self) -> i32 {
        self.start
    }

    pub fn end(&self) -> i32 {
        self.stop
    }

    pub fn ref_allele(&self) -> &Allele {
        &self.ref_allele
    }

    pub fn alt_allele(&self) -> &Allele {
        &self.alt_allele
    }

    /// `isSNP()`.
    pub fn is_snp(&self) -> bool {
        self.ref_allele.len() == 1 && self.ref_allele.len() == self.alt_allele.len()
    }

    /// `isIndel()`.
    pub fn is_indel(&self) -> bool {
        self.ref_allele.len() != self.alt_allele.len() && !self.alt_allele.is_symbolic()
    }

    /// `isSimpleInsertion()`.
    pub fn is_simple_insertion(&self) -> bool {
        self.ref_allele.len() == 1 && self.alt_allele.len() > 1
    }

    /// `isSimpleDeletion()`.
    pub fn is_simple_deletion(&self) -> bool {
        self.ref_allele.len() > 1 && self.alt_allele.len() == 1
    }

    /// `isMNP()`.
    pub fn is_mnp(&self) -> bool {
        self.ref_allele.len() > 1 && self.ref_allele.len() == self.alt_allele.len()
    }
}

/// `Event.equals`: the start and the two alleles, not the contig.
impl PartialEq for Event {
    fn eq(&self, other: &Self) -> bool {
        self.start == other.start
            && self.ref_allele == other.ref_allele
            && self.alt_allele == other.alt_allele
    }
}

/// `makeMinimalRepresentation`: the trailing bases the alleles share are dropped, unless either
/// is one base long or their last bases differ.
fn make_minimal_representation(
    reference: Allele,
    alt: Allele,
) -> Result<(Allele, Allele), EventError> {
    let ref_bases = reference.display_string().into_bytes();
    let alt_bases = alt.display_string().into_bytes();
    let different_last_base =
        ref_bases.is_empty() || alt_bases.is_empty() || ref_bases.last() != alt_bases.last();
    if reference.len() == 1 || alt.len() == 1 || different_last_base {
        return Ok((reference, alt));
    }
    if ref_bases == alt_bases {
        return Err(EventError::IllegalArgument(
            "ref and alt alleles are identical",
        ));
    }
    let min_len = ref_bases.len().min(alt_bases.len());
    let mut overlap = 0;
    while overlap < min_len
        && ref_bases[ref_bases.len() - 1 - overlap] == alt_bases[alt_bases.len() - 1 - overlap]
    {
        overlap += 1;
    }
    Ok((
        Allele::create(&ref_bases[..ref_bases.len() - overlap], true)?,
        Allele::create(&alt_bases[..alt_bases.len() - overlap], false)?,
    ))
}

/// `BaseUtils.isRegularBase`.
fn is_regular_base(base: u8) -> bool {
    simple_base_to_base_index(base) != -1
}

/// `EventMap`: the events of a haplotype by start, a `TreeMap<Integer, Event>`.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct EventMap {
    events: BTreeMap<i32, Event>,
}

impl EventMap {
    /// `new EventMap(events)`: each added in turn, compounded with the one already at its start.
    pub fn new(events: Vec<Event>) -> Result<EventMap, EventError> {
        let mut map = EventMap::default();
        for event in events {
            map.add_event(event)?;
        }
        Ok(map)
    }

    /// `fromHaplotype(haplotype, ref, refLoc, maxMnpDistance)`.
    pub fn from_haplotype(
        haplotype: &Haplotype,
        reference: &[u8],
        ref_loc: &SimpleInterval,
        max_mnp_distance: i32,
    ) -> Result<EventMap, EventError> {
        EventMap::new(get_events(haplotype, reference, ref_loc, max_mnp_distance)?)
    }

    /// `addEvent`: `computeIfPresent` with the compound, then `putIfAbsent`.
    fn add_event(&mut self, event: Event) -> Result<(), EventError> {
        let merged = match self.events.get(&event.start) {
            Some(old) => make_compound_events(old, &event)?,
            None => event,
        };
        self.events.insert(merged.start, merged);
        Ok(())
    }

    /// `getStartPositions()`.
    pub fn start_positions(&self) -> impl Iterator<Item = i32> + '_ {
        self.events.keys().copied()
    }

    /// `getEvents()`, in start order.
    pub fn events(&self) -> impl Iterator<Item = &Event> + '_ {
        self.events.values()
    }

    /// `get(start)`, the `TreeMap` lookup.
    pub fn get(&self, start: i32) -> Option<&Event> {
        self.events.get(&start)
    }

    /// `getNumberOfEvents()`.
    pub fn number_of_events(&self) -> usize {
        self.events.len()
    }

    /// `getOverlappingEvents(loc)`: the events starting at or before `loc` and ending at or after
    /// it, without a deletion that ends at `loc` when an insertion is among them.
    pub fn overlapping_events(&self, loc: i32) -> Vec<&Event> {
        let overlapping: Vec<&Event> = self
            .events
            .range(..=loc)
            .map(|(_, e)| e)
            .filter(|e| e.end() >= loc)
            .collect();
        if overlapping.iter().any(|e| e.is_simple_insertion()) {
            overlapping
                .into_iter()
                .filter(|e| !(e.is_simple_deletion() && e.end() == loc))
                .collect()
        } else {
            overlapping
        }
    }
}

/// `makeCompoundEvents`.
pub fn make_compound_events(e1: &Event, e2: &Event) -> Result<Event, EventError> {
    if e1.start != e2.start {
        return Err(EventError::IllegalArgument(
            "e1 and e2 must have the same start",
        ));
    }
    if e1.is_snp() || e2.is_snp() {
        if e1.is_snp() && e2.is_snp() {
            return Err(EventError::IllegalArgument(
                "Trying to put two overlapping SNPs in one EventMap.  This could be a CIGAR bug.",
            ));
        }
        let (snp, indel) = if e1.is_snp() { (e1, e2) } else { (e2, e1) };
        if snp.ref_allele == indel.ref_allele {
            // SNP plus insertion: A to G and A to CT make A to GT.
            let alt = format!(
                "{}{}",
                snp.alt_allele.display_string(),
                &indel.alt_allele.display_string()[1..]
            );
            Event::new(
                &snp.contig,
                snp.start,
                snp.ref_allele.clone(),
                Allele::create(alt.as_bytes(), false)?,
            )
        } else {
            // SNP plus deletion: A to T and AC to A make AC to T.
            Event::new(
                &snp.contig,
                snp.start,
                indel.ref_allele.clone(),
                snp.alt_allele.clone(),
            )
        }
    } else {
        // Insertion plus deletion: AC to A and A to AGT make AC to AGT.
        if !((e1.is_simple_deletion() && e2.is_simple_insertion())
            || (e1.is_simple_insertion() && e2.is_simple_deletion()))
        {
            return Err(EventError::IllegalArgument(
                "Can only merge single insertion with deletion (or vice versa)",
            ));
        }
        let (insertion, deletion) = if e1.is_simple_insertion() {
            (e1, e2)
        } else {
            (e2, e1)
        };
        Event::new(
            &e1.contig,
            e1.start,
            deletion.ref_allele.clone(),
            insertion.alt_allele.clone(),
        )
    }
}

/// `getEvents`: the CIGAR walked from `alignmentStartHapwrtRef`.
fn get_events(
    haplotype: &Haplotype,
    reference: &[u8],
    ref_loc: &SimpleInterval,
    max_mnp_distance: i32,
) -> Result<Vec<Event>, EventError> {
    if max_mnp_distance < 0 {
        return Err(EventError::IllegalArgument(
            "maxMnpDistance may not be negative.",
        ));
    }
    let cigar = haplotype.cigar().ok_or(EventError::NoCigar)?;
    let alignment = haplotype.bases();
    if haplotype.alignment_start_hap_wrt_ref() < 0 {
        // Protection against Smith-Waterman failures.
        return Ok(Vec::new());
    }
    let mut ref_pos = haplotype.alignment_start_hap_wrt_ref() as usize;
    // `Arrays.copyOfRange`: past the end it pads with zeros, which are not regular bases.
    let slice = |bases: &[u8], from: usize, to: usize| -> Result<Vec<u8>, EventError> {
        if from > bases.len() {
            return Err(EventError::IndexOutOfBounds);
        }
        let mut copy = bases[from..to.min(bases.len())].to_vec();
        copy.resize(to - from, 0);
        Ok(copy)
    };
    let at = |bases: &[u8], i: usize| bases.get(i).copied().ok_or(EventError::IndexOutOfBounds);

    let mut proposed = Vec::new();
    let mut alignment_pos = 0usize;
    let count = cigar.elements.len();
    for (index, element) in cigar.elements.iter().enumerate() {
        let length = element.length as usize;
        match element.op {
            Op::I => {
                // No insertion at the start of the window, or not resolved within the haplotype.
                if ref_pos > 0 && index > 0 && index < count - 1 {
                    let insertion_start = ref_loc.start + ref_pos as i32 - 1;
                    let ref_byte = at(reference, ref_pos - 1)?;
                    let mut insertion = vec![ref_byte];
                    insertion.extend(slice(&alignment, alignment_pos, alignment_pos + length)?);
                    if insertion.iter().all(|&b| is_regular_base(b)) {
                        proposed.push(Event::new(
                            &ref_loc.contig,
                            insertion_start,
                            Allele::create(&[ref_byte], true)?,
                            Allele::create(&insertion, false)?,
                        )?);
                    }
                }
                alignment_pos += length;
            }
            Op::S => alignment_pos += length,
            Op::D => {
                // No deletion at the start of the window.
                if ref_pos > 0 {
                    let deletion = slice(reference, ref_pos - 1, ref_pos + length)?;
                    let deletion_start = ref_loc.start + ref_pos as i32 - 1;
                    let ref_byte = at(reference, ref_pos - 1)?;
                    if is_regular_base(ref_byte) && deletion.iter().all(|&b| is_regular_base(b)) {
                        proposed.push(Event::new(
                            &ref_loc.contig,
                            deletion_start,
                            Allele::create(&deletion, true)?,
                            Allele::create(&[ref_byte], false)?,
                        )?);
                    }
                }
                ref_pos += length;
            }
            Op::M | Op::Eq | Op::X => {
                let mut mismatches = std::collections::VecDeque::new();
                for offset in 0..length {
                    let ref_byte = at(reference, ref_pos + offset)?;
                    let alt_byte = at(&alignment, alignment_pos + offset)?;
                    if ref_byte != alt_byte
                        && is_regular_base(ref_byte)
                        && is_regular_base(alt_byte)
                    {
                        mismatches.push_back(offset);
                    }
                }
                while let Some(start) = mismatches.pop_front() {
                    let mut end = start;
                    while let Some(&next) = mismatches.front() {
                        if (next - end) as i32 > max_mnp_distance {
                            break;
                        }
                        end = next;
                        mismatches.pop_front();
                    }
                    let ref_allele = Allele::create(
                        &slice(reference, ref_pos + start, ref_pos + end + 1)?,
                        true,
                    )?;
                    let alt_allele = Allele::create(
                        &slice(&alignment, alignment_pos + start, alignment_pos + end + 1)?,
                        false,
                    )?;
                    proposed.push(Event::new(
                        &ref_loc.contig,
                        ref_loc.start + (ref_pos + start) as i32,
                        ref_allele,
                        alt_allele,
                    )?);
                }
                ref_pos += length;
                alignment_pos += length;
            }
            other => return Err(EventError::UnsupportedOperator(other.to_char() as char)),
        }
    }
    Ok(proposed)
}

/// `buildEventMapsForHaplotypes`: each haplotype's map, set on it.
pub fn build_event_maps_for_haplotypes(
    haplotypes: &mut [Haplotype],
    reference: &[u8],
    ref_loc: &SimpleInterval,
    max_mnp_distance: i32,
) -> Result<(), EventError> {
    if max_mnp_distance < 0 {
        return Err(EventError::IllegalArgument(
            "maxMnpDistance may not be negative.",
        ));
    }
    for haplotype in haplotypes.iter_mut() {
        let map = EventMap::from_haplotype(haplotype, reference, ref_loc, max_mnp_distance)?;
        haplotype.set_event_map(map);
    }
    Ok(())
}

/// `getEventStartPositions`: every start of every haplotype's map, sorted and once each.
pub fn event_start_positions(haplotypes: &[Haplotype]) -> Result<Vec<i32>, EventError> {
    let mut starts = std::collections::BTreeSet::new();
    for haplotype in haplotypes {
        let map = haplotype.event_map().ok_or(EventError::NoEventMap)?;
        starts.extend(map.start_positions());
    }
    Ok(starts.into_iter().collect())
}

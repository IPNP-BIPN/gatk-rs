//! `HaplotypeBasedVariantRecaller`: every allele a haplotype carries, scored against every read
//! that haplotype spans.
//!
//! The PairHMM that produces the likelihoods is not ported. What is ported is everything around
//! it: which haplotype group a variant is scored against, how a matrix line is built and sorted,
//! and the ways a line is dropped or comes out wrong.
//!
//! Ported from
//! `org.broadinstitute.hellbender.tools.walkers.variantrecalling.HaplotypeRegionWalker` and
//! `org.broadinstitute.hellbender.tools.walkers.variantrecalling.VariantRecallerResultWriter`
//! in GATK 4.6.2.0.

/// A half-open-free interval, closed at both ends, which is what every span here is.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Span {
    pub contig: String,
    pub start: i32,
    pub end: i32,
}

impl Span {
    pub fn new(contig: &str, start: i32, end: i32) -> Span {
        Span {
            contig: contig.to_string(),
            start,
            end,
        }
    }

    /// `SimpleInterval.contains`, which needs the whole of the other interval.
    pub fn contains(&self, other: &Span) -> bool {
        self.contig == other.contig && self.start <= other.start && self.end >= other.end
    }

    pub fn overlaps(&self, other: &Span) -> bool {
        self.contig == other.contig && self.start <= other.end && self.end >= other.start
    }
}

impl std::fmt::Display for Span {
    /// `SimpleInterval.toString`, which is what the matrix's header line carries.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}-{}", self.contig, self.start, self.end)
    }
}

// ================================================================================================
// The haplotype groups.
// ================================================================================================

/// The prefix that makes a record in the haplotype BAM a haplotype.
///
/// Any other record is passed over however well it fits, so the file may hold anything else
/// alongside them.
pub const HAPLOTYPE_NAME_PREFIX: &str = "HC_";

pub fn is_haplotype_record(name: &str) -> bool {
    name.starts_with(HAPLOTYPE_NAME_PREFIX)
}

/// One record of the haplotype BAM, reduced to what the walk reads off it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HaplotypeRecord {
    pub name: String,
    pub span: Span,
}

/// The groups a query yields, in the order the reader hands them over.
///
/// A record whose span differs from the group's first CLOSES the group and opens a new one, so
/// two runs of the same span separated by a third are two groups rather than one.
pub fn groups(records: &[HaplotypeRecord]) -> Vec<Vec<HaplotypeRecord>> {
    let mut groups: Vec<Vec<HaplotypeRecord>> = Vec::new();
    let mut current: Vec<HaplotypeRecord> = Vec::new();
    for record in records
        .iter()
        .filter(|record| is_haplotype_record(&record.name))
    {
        if let Some(first) = current.first() {
            if first.span != record.span {
                groups.push(std::mem::take(&mut current));
            }
        }
        current.push(record.clone());
    }
    if !current.is_empty() {
        groups.push(current);
    }
    groups
}

/// `fitnessScore`: one less twice the distance of the variant from the group's halfway point.
///
/// Both gaps are floored at one before the ratio is taken, so a variant flush against either end
/// does not divide by zero and does not score zero either. An empty group scores zero.
pub fn fitness_score(location: &Span, group: &[HaplotypeRecord]) -> f64 {
    let Some(first) = group.first() else {
        return 0.0;
    };
    let before = std::cmp::max(1, location.start - first.span.start) as f64;
    let after = std::cmp::max(1, first.span.end - location.end) as f64;
    1.0 - 2.0 * (0.5 - before / (before + after)).abs()
}

/// `forBest`: the group with the highest fitness, ties going to the FIRST.
///
/// The comparison is strict, so a later group has to beat the one held rather than match it.
pub fn best_group<'a>(
    location: &Span,
    groups: &'a [Vec<HaplotypeRecord>],
) -> Option<&'a Vec<HaplotypeRecord>> {
    let mut best: Option<&Vec<HaplotypeRecord>> = None;
    for group in groups {
        match best {
            None => best = Some(group),
            Some(held) if fitness_score(location, group) > fitness_score(location, held) => {
                best = Some(group)
            }
            Some(_) => {}
        }
    }
    best.filter(|group| !group.is_empty())
}

// ================================================================================================
// The cigar walk.
// ================================================================================================

/// One cigar element: an operator and a length.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarElement {
    pub operator: char,
    pub length: i32,
}

impl CigarElement {
    pub fn consumes_read_bases(self) -> bool {
        matches!(self.operator, 'M' | 'I' | 'S' | '=' | 'X')
    }

    pub fn consumes_reference_bases(self) -> bool {
        matches!(self.operator, 'M' | 'D' | 'N' | '=' | 'X')
    }
}

/// A cigar string, as its elements.
pub fn parse_cigar(text: &str) -> Vec<CigarElement> {
    let mut elements = Vec::new();
    let mut length = 0i32;
    for character in text.chars() {
        if let Some(digit) = character.to_digit(10) {
            length = length * 10 + digit as i32;
        } else {
            elements.push(CigarElement {
                operator: character,
                length,
            });
            length = 0;
        }
    }
    elements
}

/// `getOffsetOnRead`: the offset in the read of the base at `offset` positions into the alignment.
///
/// The walk has a hole in it. A read-consuming element returns as soon as the remaining offset
/// fits inside it, and the reference-consuming subtraction happens AFTERWARDS, so a deletion
/// drives the offset negative and the very next match element then returns
/// `read_offset + offset` with a negative offset: an index that many bases too early rather than
/// the refusal a variant inside a deletion should be. Only an offset that runs off the end of the
/// read returns nothing.
pub fn offset_on_read(cigar: &[CigarElement], initial: i32) -> Option<i32> {
    let mut read_offset = 0i32;
    let mut offset = initial;
    for element in cigar {
        if element.consumes_read_bases() {
            if offset < element.length {
                return Some(read_offset + offset);
            }
            read_offset += element.length;
        }
        if element.consumes_reference_bases() {
            offset -= element.length;
        }
    }
    None
}

// ================================================================================================
// The matrix lines.
// ================================================================================================

/// One read, reduced to what a matrix line carries.
#[derive(Debug, Clone, PartialEq)]
pub struct Read {
    pub name: String,
    pub span: Span,
    pub cigar: Vec<CigarElement>,
    pub bases: Vec<u8>,
    pub is_duplicate: bool,
    pub is_reverse: bool,
    pub mapping_quality: i32,
    /// The flow key length, which is zero for a read that is not flow-based.
    pub key_length: i32,
    pub sample: String,
    pub unclipped_start: i32,
    pub unclipped_end: i32,
}

/// `Double.toString` for the likelihood columns, which is Rust's `{:?}` for these values.
fn java_double(value: f64) -> String {
    gatk_engine::tsv_table::java_double_to_string(value)
}

/// The variant's type, which decides only whether its end is written.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantKind {
    Mixed,
    Other,
}

/// The header line one variant gets: its position, the haplotype span and the alleles.
///
/// The end is omitted for a one-base variant AND for a MIXED one however long it is, which is the
/// one place the variant's type is consulted at all.
pub fn header_line(
    contig: &str,
    start: i32,
    end: i32,
    kind: VariantKind,
    haplotype_span: &Span,
    alleles: &[String],
) -> String {
    let mut line = format!("#{contig}:{start}");
    if kind != VariantKind::Mixed && end != start {
        line.push_str(&format!("-{end}"));
    }
    line.push_str(&format!(" {haplotype_span}"));
    for allele in alleles {
        line.push(' ');
        line.push_str(allele);
    }
    line
}

/// One matrix line, and the key it is sorted by.
#[derive(Debug, Clone, PartialEq)]
pub struct MatrixLine {
    pub sort_key: f64,
    pub text: String,
}

/// The bases a read carries over the variant, and the unclipped offset of the first of them.
///
/// Nothing comes back when the read does not span the whole variant, and nothing comes back when
/// the cigar walk runs off the end of the read: both are lines that are never added.
pub fn variant_bases(read: &Read, variant: &Span) -> Option<(String, i32)> {
    if !read.span.contains(variant) {
        return None;
    }
    let offset = variant.start - read.span.start;
    let length = variant.end - variant.start + 1;
    let mut bases = String::new();
    let mut first_unclipped = 0;
    for i in 0..length {
        let read_offset = offset_on_read(&read.cigar, offset + i)?;
        bases.push(*read.bases.get(read_offset as usize)? as char);
        first_unclipped = if read.is_reverse {
            (read.bases.len() as i32 - read_offset - 1) + (read.unclipped_end - read.span.end)
        } else {
            read_offset + (read.span.start - read.unclipped_start)
        };
    }
    Some((bases, first_unclipped))
}

/// One line of the matrix, or nothing.
///
/// Two things drop a line. Every likelihood being negative infinity is an alignment failure, and
/// the read not yielding bases over the variant leaves the line unfinished. The SORT KEY is the
/// LAST allele's likelihood rather than the best of them: the loop that collects the columns
/// overwrites it each time round.
pub fn matrix_line(read: &Read, variant: &Span, likelihoods: &[f64]) -> Option<MatrixLine> {
    if likelihoods.iter().all(|value| *value == f64::NEG_INFINITY) {
        return None;
    }
    let sort_key = *likelihoods.last().unwrap_or(&f64::NEG_INFINITY);
    let (bases, first_unclipped) = variant_bases(read, variant)?;
    if bases.is_empty() {
        return None;
    }
    let columns: Vec<String> = likelihoods
        .iter()
        .map(|value| java_double(*value))
        .collect();
    Some(MatrixLine {
        sort_key,
        text: format!(
            "{} {} {} {} {} {} {} {} {}",
            read.name,
            read.key_length,
            if read.is_duplicate { 1 } else { 0 },
            if read.is_reverse { 1 } else { 0 },
            read.mapping_quality,
            columns.join(" "),
            bases,
            first_unclipped,
            read.sample
        ),
    })
}

/// The lines of one variant, sorted by their key, descending.
///
/// The sort is `-Double.compare(a, b)`, which is stable, so two reads with the same last column
/// keep the order they were evidenced in.
pub fn sorted_lines(lines: &[MatrixLine]) -> Vec<String> {
    let mut lines = lines.to_vec();
    lines.sort_by(|a, b| compare_double(b.sort_key, a.sort_key));
    lines.into_iter().map(|line| line.text).collect()
}

/// `Double.compare`: numeric first, then the signed bit pattern.
fn compare_double(a: f64, b: f64) -> std::cmp::Ordering {
    if a < b {
        return std::cmp::Ordering::Less;
    }
    if a > b {
        return std::cmp::Ordering::Greater;
    }
    let bits = |value: f64| {
        if value.is_nan() {
            f64::NAN.to_bits() as i64
        } else {
            value.to_bits() as i64
        }
    };
    bits(a).cmp(&bits(b))
}

/// One variant's whole block: its header line, then its sorted matrix lines.
pub fn variant_block(header: &str, lines: &[MatrixLine]) -> String {
    let mut text = String::from(header);
    text.push('\n');
    for line in sorted_lines(lines) {
        text.push_str(&line);
        text.push('\n');
    }
    text
}

// ================================================================================================
// The events a haplotype carries, and the alleles they merge into at one position.
// ================================================================================================

/// An allele as the likelihoods are keyed by it: bases and whether it is the reference. The
/// spanning deletion is the one-base `*`, which is not symbolic.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Allele {
    pub bases: Vec<u8>,
    pub is_reference: bool,
}

impl Allele {
    pub fn new(bases: &[u8], is_reference: bool) -> Allele {
        Allele {
            bases: bases.to_vec(),
            is_reference,
        }
    }

    pub fn span_del() -> Allele {
        Allele::new(b"*", false)
    }

    /// `Allele.toString`: the bases, and a `*` after them for the reference.
    pub fn text(&self) -> String {
        let mut text = String::from_utf8_lossy(&self.bases).into_owned();
        if self.is_reference {
            text.push('*');
        }
        text
    }
}

/// `Event`: one change a haplotype makes, as a reference and an alternate allele at a start.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Event {
    pub start: i32,
    pub reference: Vec<u8>,
    pub alternate: Vec<u8>,
}

impl Event {
    pub fn end(&self) -> i32 {
        self.start + self.reference.len() as i32 - 1
    }

    fn is_snp(&self) -> bool {
        self.reference.len() == 1 && self.alternate.len() == 1
    }

    fn is_simple_insertion(&self) -> bool {
        self.reference.len() == 1 && self.alternate.len() > 1
    }

    fn is_simple_deletion(&self) -> bool {
        self.alternate.len() == 1 && self.reference.len() > 1
    }
}

fn is_regular_base(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
}

/// `EventMap.getEvents` followed by the map's own merging: the events a haplotype's cigar makes
/// against the reference window it starts at, one per start, an indel and a SNP at the same start
/// compounded and an insertion and a deletion there combined.
///
/// `Err` is the `GATKException` for a cigar operator the walk does not take, or the
/// `IllegalArgumentException` two events that cannot share a start raise.
pub fn event_map(
    haplotype: &[u8],
    cigar: &[CigarElement],
    reference: &[u8],
    reference_start: i32,
    max_mnp_distance: i32,
) -> Result<std::collections::BTreeMap<i32, Event>, (String, String)> {
    let mut proposed: Vec<Event> = Vec::new();
    let mut ref_pos: usize = 0;
    let mut alignment_pos: usize = 0;
    let count = cigar.len();
    for (index, element) in cigar.iter().enumerate() {
        let length = element.length as usize;
        match element.operator {
            'I' => {
                if ref_pos > 0 && index > 0 && index < count - 1 {
                    let ref_byte = reference[ref_pos - 1];
                    let mut bases = vec![ref_byte];
                    bases.extend_from_slice(&haplotype[alignment_pos..alignment_pos + length]);
                    if bases.iter().all(|b| is_regular_base(*b)) {
                        proposed.push(Event {
                            start: reference_start + ref_pos as i32 - 1,
                            reference: vec![ref_byte],
                            alternate: bases,
                        });
                    }
                }
                alignment_pos += length;
            }
            'S' => alignment_pos += length,
            'D' => {
                if ref_pos > 0 {
                    let deleted = reference[ref_pos - 1..ref_pos + length].to_vec();
                    let ref_byte = reference[ref_pos - 1];
                    if is_regular_base(ref_byte) && deleted.iter().all(|b| is_regular_base(*b)) {
                        proposed.push(Event {
                            start: reference_start + ref_pos as i32 - 1,
                            reference: deleted,
                            alternate: vec![ref_byte],
                        });
                    }
                }
                ref_pos += length;
            }
            'M' | '=' | 'X' => {
                let mut mismatches: std::collections::VecDeque<usize> =
                    std::collections::VecDeque::new();
                for offset in 0..length {
                    let ref_byte = reference[ref_pos + offset];
                    let alt_byte = haplotype[alignment_pos + offset];
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
                        if (next - end) as i32 <= max_mnp_distance {
                            end = next;
                            mismatches.pop_front();
                        } else {
                            break;
                        }
                    }
                    proposed.push(Event {
                        start: reference_start + (ref_pos + start) as i32,
                        reference: reference[ref_pos + start..ref_pos + end + 1].to_vec(),
                        alternate: haplotype[alignment_pos + start..alignment_pos + end + 1]
                            .to_vec(),
                    });
                }
                ref_pos += length;
                alignment_pos += length;
            }
            other => {
                return Err((
                    "org.broadinstitute.hellbender.exceptions.GATKException".to_string(),
                    format!("Unsupported cigar operator created during SW alignment: {other}"),
                ))
            }
        }
    }
    let mut map: std::collections::BTreeMap<i32, Event> = std::collections::BTreeMap::new();
    for event in proposed {
        let merged = match map.get(&event.start) {
            Some(old) => compound(old, &event)?,
            None => event.clone(),
        };
        map.insert(event.start, merged);
    }
    Ok(map)
}

/// `EventMap.makeCompoundEvents`.
fn compound(e1: &Event, e2: &Event) -> Result<Event, (String, String)> {
    let refusal = |message: &str| {
        (
            "java.lang.IllegalArgumentException".to_string(),
            message.to_string(),
        )
    };
    if e1.is_snp() || e2.is_snp() {
        if e1.is_snp() && e2.is_snp() {
            return Err(refusal(
                "Trying to put two overlapping SNPs in one EventMap.  This could be a CIGAR bug.",
            ));
        }
        let (snp, indel) = if e1.is_snp() { (e1, e2) } else { (e2, e1) };
        if snp.reference == indel.reference {
            let mut alternate = snp.alternate.clone();
            alternate.extend_from_slice(&indel.alternate[1..]);
            Ok(Event {
                start: snp.start,
                reference: snp.reference.clone(),
                alternate,
            })
        } else {
            Ok(Event {
                start: snp.start,
                reference: indel.reference.clone(),
                alternate: snp.alternate.clone(),
            })
        }
    } else {
        if !((e1.is_simple_deletion() && e2.is_simple_insertion())
            || (e1.is_simple_insertion() && e2.is_simple_deletion()))
        {
            return Err(refusal(
                "Can only merge single insertion with deletion (or vice versa)",
            ));
        }
        let (insertion, deletion) = if e1.is_simple_insertion() {
            (e1, e2)
        } else {
            (e2, e1)
        };
        Ok(Event {
            start: e1.start,
            reference: deletion.reference.clone(),
            alternate: insertion.alternate.clone(),
        })
    }
}

/// `EventMap.getOverlappingEvents`: the events starting at or before `loc` that reach it, a
/// deletion ending exactly there dropped when an insertion is among them.
pub fn overlapping_events(map: &std::collections::BTreeMap<i32, Event>, loc: i32) -> Vec<Event> {
    let overlapping: Vec<Event> = map
        .range(..=loc)
        .map(|(_, event)| event.clone())
        .filter(|event| event.end() >= loc)
        .collect();
    if overlapping.iter().any(Event::is_simple_insertion) {
        overlapping
            .into_iter()
            .filter(|event| !(event.is_simple_deletion() && event.end() == loc))
            .collect()
    } else {
        overlapping
    }
}

/// The merged alleles at `loc` over a group's event maps, in the order `simpleMerge` leaves them:
/// the longest reference, then every alternate in the order the haplotypes first carry it, each
/// extended to that reference; an event starting before `loc` counts as the spanning deletion.
/// `None` where no haplotype has an event over `loc`.
pub fn merged_alleles(
    maps: &[std::collections::BTreeMap<i32, Event>],
    loc: i32,
    reference_base: u8,
) -> Result<Option<Vec<Allele>>, (String, String)> {
    // getVariantsFromActiveHaplotypes, each event once, then replaceSpanDels.
    let mut seen: Vec<Event> = Vec::new();
    let mut vcs: Vec<(Allele, Allele)> = Vec::new();
    for map in maps {
        for event in overlapping_events(map, loc) {
            if seen.contains(&event) {
                continue;
            }
            seen.push(event.clone());
            if event.start == loc {
                vcs.push((
                    Allele::new(&event.reference, true),
                    Allele::new(&event.alternate, false),
                ));
            } else {
                vcs.push((Allele::new(&[reference_base], true), Allele::span_del()));
            }
        }
    }
    if vcs.is_empty() {
        return Ok(None);
    }
    // determineReferenceAllele: the longest, two of one length having to agree.
    let mut reference: Option<Allele> = None;
    for (vc_ref, _) in &vcs {
        reference = Some(match reference {
            None => vc_ref.clone(),
            Some(held) => {
                if held.bases.len() < vc_ref.bases.len() {
                    vc_ref.clone()
                } else if vc_ref.bases.len() < held.bases.len() {
                    held
                } else if held != *vc_ref {
                    return Err((
                        "java.lang.IllegalStateException".to_string(),
                        format!(
                            "The provided variant file(s) have inconsistent references for the same position(s) at chr:{loc}, {} vs. {}",
                            held.text(),
                            vc_ref.text()
                        ),
                    ));
                } else {
                    held
                }
            }
        });
    }
    let reference = reference.expect("at least one event");
    let mut alleles: Vec<Allele> = Vec::new();
    let add = |allele: Allele, alleles: &mut Vec<Allele>| {
        if !alleles.contains(&allele) {
            alleles.push(allele);
        }
    };
    for (vc_ref, alt) in &vcs {
        if *vc_ref == reference {
            add(vc_ref.clone(), &mut alleles);
            add(alt.clone(), &mut alleles);
        } else {
            // createAlleleMapping's values, then the reference put LAST.
            let extra = &reference.bases[vc_ref.bases.len()..];
            if *alt == Allele::span_del() {
                add(alt.clone(), &mut alleles);
            } else {
                let mut bases = alt.bases.clone();
                bases.extend_from_slice(extra);
                add(Allele::new(&bases, false), &mut alleles);
            }
            add(reference.clone(), &mut alleles);
        }
    }
    // `VariantContext.makeAlleles` moves the reference to the front.
    let position = alleles
        .iter()
        .position(|allele| allele.is_reference)
        .expect("the reference allele");
    let reference = alleles.remove(position);
    alleles.insert(0, reference);
    Ok(Some(alleles))
}

/// `createAlleleMapper`: for each merged allele, the haplotypes that carry it at `loc`, keyed in
/// the merged order with the spanning deletion added when a haplotype needs it and the merged
/// alleles do not already hold it.
pub fn allele_mapper(
    merged: &[Allele],
    maps: &[std::collections::BTreeMap<i32, Event>],
    loc: i32,
    emit_spanning_dels: bool,
) -> Vec<(Allele, Vec<usize>)> {
    let reference = &merged[0];
    let mut result: Vec<(Allele, Vec<usize>)> = vec![(reference.clone(), Vec::new())];
    for alt in &merged[1..] {
        if !result.iter().any(|(key, _)| key == alt) {
            result.push((alt.clone(), Vec::new()));
        }
    }
    let index_of = |result: &Vec<(Allele, Vec<usize>)>, allele: &Allele| {
        result.iter().position(|(key, _)| key == allele)
    };
    for (haplotype, map) in maps.iter().enumerate() {
        let overlapping = overlapping_events(map, loc);
        if overlapping.is_empty() {
            result[0].1.push(haplotype);
            continue;
        }
        for event in overlapping {
            if event.start == loc {
                let alt = if event.reference.len() == reference.bases.len() {
                    Allele::new(&event.alternate, false)
                } else if event.reference.len() < reference.bases.len() {
                    let mut bases = event.alternate.clone();
                    bases.extend_from_slice(&reference.bases[event.reference.len()..]);
                    Allele::new(&bases, false)
                } else {
                    continue;
                };
                if let Some(at) = index_of(&result, &alt) {
                    result[at].1.push(haplotype);
                }
            } else if emit_spanning_dels {
                let span = Allele::span_del();
                let at = match index_of(&result, &span) {
                    Some(at) => at,
                    None => {
                        result.push((span, Vec::new()));
                        result.len() - 1
                    }
                };
                result[at].1.push(haplotype);
                // The first spanning event ends the haplotype's walk, whatever follows it.
                break;
            } else {
                result[0].1.push(haplotype);
                break;
            }
        }
    }
    result
}

// ================================================================================================
// The PairHMM engine: the reads' qualities prepared, every read against every haplotype, and the
// poorly modelled reads taken away.
// ================================================================================================

/// `PCRErrorModel`'s rate factors.
pub fn pcr_rate_factor(model: &str) -> f64 {
    match model {
        "HOSTILE" => 1.0,
        "AGGRESSIVE" => 2.0,
        "CONSERVATIVE" => 3.0,
        _ => 0.0,
    }
}

/// `ReadLikelihoodCalculationEngine.MAX_REPEAT_LENGTH` and `MAX_STR_UNIT_LENGTH`.
const MAX_REPEAT_LENGTH: i32 = 20;
const MAX_STR_UNIT_LENGTH: usize = 8;

/// `findTandemRepeatUnits(readBases, offset).getRight()`: the repeat length at `offset`, counted
/// backwards and, where the next unit is the same, forwards too, capped at twenty.
pub fn tandem_repeat_length(bases: &[u8], offset: usize) -> i32 {
    use gatk_annotation::tandem_repeat::find_number_of_repetitions as repetitions;
    let mut max_bw = 0;
    let mut best_bw: Vec<u8> = vec![bases[offset]];
    for unit in 1..=MAX_STR_UNIT_LENGTH {
        if offset + 1 < unit {
            break;
        }
        let candidate = &bases[offset + 1 - unit..offset + 1];
        max_bw = repetitions(candidate, &bases[..offset + 1], false).unwrap_or(0);
        if max_bw > 1 {
            best_bw = candidate.to_vec();
            break;
        }
    }
    let mut max_rl = max_bw;
    if offset + 1 < bases.len() {
        let mut best_fw: Vec<u8> = vec![bases[offset + 1]];
        let mut max_fw = 0;
        for unit in 1..=MAX_STR_UNIT_LENGTH {
            if offset + unit + 1 > bases.len() {
                break;
            }
            let candidate = &bases[offset + 1..offset + unit + 1];
            max_fw = repetitions(candidate, &bases[offset + 1..], true).unwrap_or(0);
            if max_fw > 1 {
                best_fw = candidate.to_vec();
                break;
            }
        }
        if best_fw == best_bw {
            max_rl = max_bw + max_fw;
        } else {
            let bw = repetitions(&best_fw, &bases[..offset + 1], false).unwrap_or(0);
            max_rl = max_fw + bw;
        }
    }
    max_rl.min(MAX_REPEAT_LENGTH)
}

/// The engine's settings, as the likelihood arguments set them.
#[derive(Debug, Clone, PartialEq)]
pub struct EngineSettings {
    pub pcr_error_model: String,
    pub base_quality_score_threshold: u8,
    pub gap_continuation_penalty: u8,
    pub disable_cap_read_qualities_to_map_q: bool,
    pub dynamic_disqualification: bool,
    pub read_disqualification_scale: f64,
    pub expected_error_rate_per_base: f64,
}

/// One read as the engine takes it: already trimmed to the haplotypes.
#[derive(Debug, Clone, PartialEq)]
pub struct EngineRead {
    pub bases: Vec<u8>,
    pub quals: Vec<u8>,
    pub insertion_quals: Option<Vec<u8>>,
    pub deletion_quals: Option<Vec<u8>>,
    pub mapping_quality: u8,
}

/// `QualityUtils.MIN_USABLE_Q_SCORE`.
const MIN_USABLE_Q_SCORE: u8 = 6;
/// `ReadUtils.DEFAULT_INSERTION_DELETION_QUAL`.
const DEFAULT_INDEL_QUAL: u8 = 45;

/// `dynamicReadQualThreshLookupTable`, the mean and variance per base quality from 1 to 40.
const DYNAMIC_TABLE: [(f64, f64); 40] = [
    (5.996842844, 0.196616587),
    (5.870018422, 1.388545569),
    (5.401558531, 5.641990128),
    (4.818940919, 10.33176216),
    (4.218758304, 14.25799688),
    (3.646319832, 17.02880749),
    (3.122346753, 18.64537883),
    (2.654731979, 19.27521677),
    (2.244479156, 19.13584613),
    (1.88893867, 18.43922003),
    (1.583645342, 17.36842261),
    (1.3233807, 16.07088712),
    (1.102785365, 14.65952563),
    (0.916703025, 13.21718577),
    (0.760361881, 11.80207947),
    (0.629457387, 10.45304833),
    (0.520175654, 9.194183767),
    (0.42918208, 8.038657241),
    (0.353590663, 6.991779595),
    (0.290923699, 6.053379213),
    (0.23906788, 5.219610436),
    (0.196230431, 4.484302033),
    (0.160897421, 3.839943445),
    (0.131795374, 3.27839108),
    (0.1078567, 2.791361596),
    (0.088189063, 2.370765375),
    (0.072048567, 2.008921719),
    (0.058816518, 1.698687797),
    (0.047979438, 1.433525748),
    (0.039111985, 1.207526336),
    (0.031862437, 1.015402928),
    (0.025940415, 0.852465956),
    (0.021106532, 0.714585285),
    (0.017163711, 0.598145851),
    (0.013949904, 0.500000349),
    (0.011332027, 0.41742159),
    (0.009200898, 0.348056286),
    (0.007467036, 0.289881373),
    (0.006057179, 0.241163527),
    (0.004911394, 0.200422214),
];

/// What the engine answers: the likelihood of every kept read against every haplotype, as
/// `[haplotype][read]`, and which of the input reads were kept.
#[derive(Debug, Clone, PartialEq)]
pub struct EngineResult {
    pub likelihoods: Vec<Vec<f64>>,
    pub kept: Vec<usize>,
}

/// `PairHMMLikelihoodCalculationEngine.computeReadLikelihoods` with the Java `LOGLESS_CACHING`
/// PairHMM, no likelihood cap (the tool's mismapping rate is off), and
/// `filterPoorlyModeledEvidence` after.
pub fn compute_likelihoods(
    haplotypes: &[Vec<u8>],
    reads: &[EngineRead],
    settings: &EngineSettings,
) -> EngineResult {
    let rate = pcr_rate_factor(&settings.pcr_error_model);
    let pcr_cache: Vec<u8> = (0..=MAX_REPEAT_LENGTH)
        .map(|i| {
            let value = 40.0
                - (std::hint::black_box(f64::from(i)) / (rate * std::f64::consts::PI)).exp()
                + 1.0;
            gatk_engine::pair_hmm::fast_round(value).max(10) as u8
        })
        .collect();
    let mut matrix = vec![Vec::with_capacity(reads.len()); haplotypes.len()];
    let mut hmm_quals: Vec<Vec<u8>> = Vec::with_capacity(reads.len());
    for read in reads {
        let length = read.bases.len();
        let mut quals = read.quals.clone();
        let mut ins = read
            .insertion_quals
            .clone()
            .unwrap_or_else(|| vec![DEFAULT_INDEL_QUAL; length]);
        let mut del = read
            .deletion_quals
            .clone()
            .unwrap_or_else(|| vec![DEFAULT_INDEL_QUAL; length]);
        if rate != 0.0 {
            for i in 1..length {
                let repeat = tandem_repeat_length(&read.bases, i - 1) as usize;
                ins[i - 1] = ins[i - 1].min(pcr_cache[repeat]);
                del[i - 1] = del[i - 1].min(pcr_cache[repeat]);
            }
        }
        for i in 0..quals.len() {
            if !settings.disable_cap_read_qualities_to_map_q {
                quals[i] = quals[i].min(read.mapping_quality);
            }
            // `setToFixedValueIfTooLow` compares SIGNED bytes.
            let low = |value: u8, min: u8| (value as i8) < (min as i8);
            if low(quals[i], settings.base_quality_score_threshold) {
                quals[i] = MIN_USABLE_Q_SCORE;
            }
            if low(ins[i], MIN_USABLE_Q_SCORE) {
                ins[i] = MIN_USABLE_Q_SCORE;
            }
            if low(del[i], MIN_USABLE_Q_SCORE) {
                del[i] = MIN_USABLE_Q_SCORE;
            }
        }
        let gcp = vec![settings.gap_continuation_penalty; length];
        for (h, haplotype) in haplotypes.iter().enumerate() {
            let value = gatk_engine::pair_hmm::read_likelihood_given_haplotype_log10(
                haplotype,
                &read.bases,
                &quals,
                &ins,
                &del,
                &gcp,
            );
            matrix[h].push(value);
        }
        hmm_quals.push(quals);
    }
    // filterPoorlyModeledEvidence.
    let mut kept = Vec::new();
    for (r, quals) in hmm_quals.iter().enumerate() {
        let best = matrix
            .iter()
            .map(|row| row[r])
            .fold(f64::NEG_INFINITY, |a, b| if b > a { b } else { a });
        let length = quals.len() as f64;
        let standard = |cap: bool| {
            let errors = (length * settings.expected_error_rate_per_base).ceil();
            let errors = if cap { errors.min(2.0) } else { errors };
            errors * -4.0
        };
        let threshold = if settings.dynamic_disqualification {
            let mut mean = 0.0;
            let mut variance = 0.0;
            for &q in quals {
                let entry = if q <= 1 { 0 } else { (q.min(40) - 1) as usize };
                mean += DYNAMIC_TABLE[entry].0;
                variance += DYNAMIC_TABLE[entry].1;
            }
            let dynamic = (mean + settings.read_disqualification_scale * variance.sqrt()) * -0.1;
            let plain = standard(false);
            if dynamic < plain {
                dynamic
            } else {
                plain
            }
        } else {
            standard(true)
        };
        // Removed only when strictly below: a NaN maximum is kept, as `<` keeps it.
        if best.partial_cmp(&threshold) != Some(std::cmp::Ordering::Less) {
            kept.push(r);
        }
    }
    let likelihoods = matrix
        .into_iter()
        .map(|row| kept.iter().map(|r| row[*r]).collect())
        .collect();
    EngineResult { likelihoods, kept }
}

//! Pileup detection, ported from `utils.pileup.PileupBasedAlleles`, the pileup half of
//! `AssemblyResultSet` (`removeHaplotypesWithBadAlleles`, `injectPileupEvents`,
//! `makeHaplotypeWithInsertedEvent`) and `AssemblyBasedCallerUtils.getKmerReadCounts` and
//! `filterPileupHaplotypes` (GATK 4.6.2.0).
//!
//! # Where the pileups come from
//!
//! The reference keeps every locus the region iterator walked, with the score `isActive` gave it,
//! and hands a region the ones inside its active span. What it reads from them is only counts, so
//! the port keeps a summary of each locus instead of its pileup: the depth, the reads below the
//! base quality floor, and for each element what it could add (a substitution, the inserted bases,
//! the deletion length) and whether its read is "bad". A read's badness is a function of the read
//! alone, the `MZ` attribute the reference caches on it being its mismatch rate.
//!
//! # Orders
//!
//! The good events are a `HashSet<Event>`, sorted by start before injection, so two events at one
//! start keep the set's order, and the set's order is its hash layout. The per-locus maps
//! (`HashMap<Byte, ..>`, `HashMap<String, ..>`, `HashMap<Integer, ..>`, each turned into a
//! `HashMap<Event, ..>`) decide the insertion order into it, which matters for events sharing a
//! bucket. All of them are modelled with [`JavaHashMap`].
//!
//! # A removal that usually removes nothing
//!
//! `removeHaplotypesWithBadAlleles` matches events against each haplotype's event map, and both
//! callers run it on a freshly trimmed set, whose haplotypes have none: `Haplotype.trim` does not
//! carry the map over. The bad events then act only in HaplotypeCaller's trimmer. The exception is
//! the fallback after the PD engine gave up, which built the maps on its way.

use gatk_engine::assembly_region::AssemblyRegion;
use gatk_engine::assembly_result_set::{AssemblyResultSet, ResultSetError};
use gatk_engine::cigar_utils::calculate_cigar;
use gatk_engine::event_map::Event;
use gatk_engine::haplotype::Haplotype;
use gatk_engine::java_hash::{byte_array_hash_code, JavaHashCode, JavaHashMap};
use gatk_engine::locus_iterator::AlignmentContext;
use gatk_engine::read::flags;
use gatk_engine::read_utils::{read_index_for_reference_coordinate, READ_INDEX_NOT_FOUND};
use gatk_engine::smith_waterman::{SmithWatermanAligner, SwOverhangStrategy, SwParameters};
use htsjdk_bam::cigar::Op;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_vcf::allele::Allele;

/// `MISMATCH_BASES_PERCENTAGE_TAG`.
pub const MISMATCH_BASES_PERCENTAGE_TAG: [u8; 2] = *b"MZ";
/// `MISMATCH_BASES_PERCENTAGE_ADJUSMTENT`.
const MISMATCH_BASES_PERCENTAGE_ADJUSTMENT: f64 = 1000.0;
/// `AssemblyBasedCallerUtils.NUM_HAPLOTYPES_TO_INJECT_FORCE_CALLING_ALLELES_INTO`.
const NUM_HAPLOTYPES_TO_INJECT_INTO: usize = 5;

/// `PileupDetectionArgumentCollection`, the detection half.
#[derive(Debug, Clone, PartialEq)]
pub struct PileupDetectionArguments {
    pub use_pileup_detection: bool,
    pub detect_indels: bool,
    pub active_region_phred_threshold: f64,
    pub num_haplotypes_to_iterate: i32,
    pub filtering_kmer_size: i32,
    pub snp_threshold: f64,
    pub indel_threshold: f64,
    pub pileup_absolute_depth: f64,
    pub snp_adjacent_to_assembly_indel: i32,
    pub quality_for_snps_in_pileup_detection: i32,
    pub bad_read_threshold: f64,
    pub bad_read_proper_pair: bool,
    pub bad_read_edit_distance: f64,
    pub bad_read_secondary_or_supplementary: bool,
    pub template_length_mean: f64,
    pub template_length_std: f64,
    pub assembly_bad_read_threshold: f64,
    pub assembly_bad_read_edit_distance: f64,
    /// `--print-pileupcalling-status`: the haplotype building printed as it goes.
    pub debug_stdout: bool,
}

impl Default for PileupDetectionArguments {
    fn default() -> Self {
        PileupDetectionArguments {
            use_pileup_detection: false,
            detect_indels: false,
            active_region_phred_threshold: 0.0,
            num_haplotypes_to_iterate: 5,
            filtering_kmer_size: 10,
            snp_threshold: 0.1,
            indel_threshold: 0.1,
            pileup_absolute_depth: 0.0,
            snp_adjacent_to_assembly_indel: 5,
            quality_for_snps_in_pileup_detection: 12,
            bad_read_threshold: 0.0,
            bad_read_proper_pair: true,
            bad_read_edit_distance: 0.08,
            bad_read_secondary_or_supplementary: true,
            template_length_mean: 0.0,
            template_length_std: 0.0,
            assembly_bad_read_threshold: 0.0,
            assembly_bad_read_edit_distance: 0.12,
            debug_stdout: false,
        }
    }
}

/// What one pileup element can contribute.
#[derive(Debug, Clone, PartialEq)]
struct ElementSummary {
    /// The base, when it is a substitution of high enough quality.
    snp: Option<u8>,
    /// The inserted bases, when the element is before an insertion.
    insertion: Option<Vec<u8>>,
    /// The deletion length, when the element is before a deletion.
    deletion: Option<u32>,
    bad_pileup: bool,
    bad_assembly: bool,
}

/// One locus of `region.getAlignmentData()`: its score and what its pileup counts.
#[derive(Debug, Clone, PartialEq)]
pub struct LocusSummary {
    pub contig: String,
    pub position: i32,
    /// `getActivityScore()`, the `originalActiveProb` of the locus's activity state.
    pub activity_score: f64,
    ref_base: u8,
    /// `alignmentContext.size()` less the elements below the base quality floor.
    depth: i32,
    elements: Vec<ElementSummary>,
}

impl LocusSummary {
    /// A locus no read covers.
    pub fn empty(contig: &str, position: i32, activity_score: f64, ref_base: u8) -> LocusSummary {
        LocusSummary {
            contig: contig.to_string(),
            position,
            activity_score,
            ref_base,
            depth: 0,
            elements: Vec::new(),
        }
    }

    /// The counting half of `goodAndBadPileupEvents` for one locus, over `contig_bases` (the whole
    /// contig, 0-based) for the reads' mismatch rates.
    pub fn new(
        context: &AlignmentContext<'_>,
        activity_score: f64,
        ref_base: u8,
        contig_bases: &[u8],
        min_base_quality_score: i32,
        arguments: &PileupDetectionArguments,
    ) -> LocusSummary {
        let mut depth = context.pileup.elements.len() as i32;
        let mut elements = Vec::new();
        for element in &context.pileup.elements {
            let base = element.base();
            let qual = element.qual() as i8 as i32;
            if qual < min_base_quality_score {
                depth -= 1;
            }
            let snp = (ref_base != base
                && base != b'D'
                && qual > arguments.quality_for_snps_in_pileup_detection)
                .then_some(base);
            let insertion = if arguments.detect_indels && element.is_before_insertion() {
                element.bases_of_immediately_following_insertion()
            } else {
                None
            };
            let deletion = (arguments.detect_indels && element.is_before_deletion_start())
                .then(|| element.length_of_immediately_following_indel());
            if snp.is_none() && insertion.is_none() && deletion.is_none() {
                continue;
            }
            let mismatch = mismatch_percentage(element.read, contig_bases);
            elements.push(ElementSummary {
                snp,
                insertion,
                deletion,
                bad_pileup: is_bad_pileup_read(element.read, mismatch, arguments),
                bad_assembly: is_bad_assembly_read(mismatch, arguments),
            });
        }
        LocusSummary {
            contig: context.contig.clone(),
            position: context.position,
            activity_score,
            ref_base,
            depth,
            elements,
        }
    }
}

/// `addMismatchPercentageToRead`: an `MZ` the read already carries, else the mismatches (the `NM`
/// tag, or `calculateSamNmTag` over the reference) less the indel bases, per thousand aligned
/// bases.
pub fn mismatch_percentage(read: &BamRecord, contig_bases: &[u8]) -> i32 {
    if let Some(TagValue::Int(value)) = read.tags.get(Tag::new(&MISMATCH_BASES_PERCENTAGE_TAG)) {
        return *value as i32;
    }
    let nm = match read.tags.get(Tag::new(b"NM")) {
        Some(TagValue::Int(value)) => *value as i32,
        _ => {
            let blocks =
                htsjdk_bam::alignment_block::alignment_blocks(&read.cigar, read.alignment_start);
            htsjdk_bam::sequence::calculate_sam_nm_tag(
                &read.read_bases,
                &blocks,
                &read.cigar,
                contig_bases,
                0,
                read.flags & flags::READ_REVERSE_STRAND != 0,
                false,
            )
        }
    };
    let sum = |keep: fn(Op) -> bool| -> i32 {
        read.cigar
            .elements
            .iter()
            .filter(|e| keep(e.op))
            .map(|e| e.length as i32)
            .sum()
    };
    let indel = sum(|op| matches!(op, Op::I | Op::D));
    let aligned = sum(|op| matches!(op, Op::M | Op::Eq | Op::X));
    // `(int) (1000.0 * adjusted / aligned)`, a saturating cast as Java's is.
    (MISMATCH_BASES_PERCENTAGE_ADJUSTMENT * f64::from(nm - indel) / f64::from(aligned)) as i32
}

/// A read as `isActive`'s pileups hold it: what identifies it among the walker's reads.
pub type ReadKey = (String, u16, i32);

pub fn read_key(read: &BamRecord) -> ReadKey {
    (read.read_name.clone(), read.flags, read.alignment_start)
}

/// `addMismatchPercentageToRead` as `isActive` calls it on every read of a non-empty pileup with
/// pileup detection: each of the region's reads the walker saw gains `MZ`, computed as it was then.
pub fn tag_mismatch_percentage(
    region: &mut gatk_engine::assembly_region::AssemblyRegion,
    seen: &std::collections::HashSet<ReadKey>,
    contig_bases: &[u8],
) {
    for read in region.reads_mut() {
        if read
            .tags
            .get(Tag::new(&MISMATCH_BASES_PERCENTAGE_TAG))
            .is_some()
            || !seen.contains(&read_key(read))
        {
            continue;
        }
        let value = mismatch_percentage(read, contig_bases);
        read.tags.insert(
            Tag::new(&MISMATCH_BASES_PERCENTAGE_TAG),
            TagValue::Int(i64::from(value)),
        );
    }
}

/// `isBadPileupRead`.
fn is_bad_pileup_read(
    read: &BamRecord,
    mismatch: i32,
    arguments: &PileupDetectionArguments,
) -> bool {
    if arguments.bad_read_threshold <= 0.0 {
        return false;
    }
    let paired = read.flags & flags::READ_PAIRED != 0;
    let proper = paired && read.flags & flags::PROPER_PAIR != 0;
    if arguments.bad_read_proper_pair && !proper {
        return true;
    }
    if arguments.bad_read_secondary_or_supplementary
        && (read.flags & flags::NOT_PRIMARY_ALIGNMENT != 0
            || read.tags.get(Tag::new(b"SA")).is_some())
    {
        return true;
    }
    if f64::from(mismatch) / MISMATCH_BASES_PERCENTAGE_ADJUSTMENT > arguments.bad_read_edit_distance
    {
        return true;
    }
    if arguments.template_length_std > 0.0 && arguments.template_length_mean > 0.0 {
        let length = f64::from(read.inferred_insert_size);
        if length < arguments.template_length_mean - 2.25 * arguments.template_length_std
            || length > arguments.template_length_mean + 2.25 * arguments.template_length_std
        {
            return true;
        }
    }
    false
}

/// `isBadAssemblyRead`.
fn is_bad_assembly_read(mismatch: i32, arguments: &PileupDetectionArguments) -> bool {
    arguments.assembly_bad_read_threshold > 0.0
        && f64::from(mismatch) / MISMATCH_BASES_PERCENTAGE_ADJUSTMENT
            > arguments.assembly_bad_read_edit_distance
}

/// An `Event` as a hash key: `HashCodeBuilder` over the start and the two alleles.
#[derive(Debug, Clone, PartialEq)]
pub struct EventKey(pub Event);

/// `SimpleAllele.hashCode`.
fn allele_hash(allele: &Allele) -> i32 {
    let bases = allele.display_string();
    byte_array_hash_code(bases.as_bytes())
        .wrapping_mul(31)
        .wrapping_add(if allele.is_reference() { 1231 } else { 1237 })
}

impl JavaHashCode for EventKey {
    fn java_hash_code(&self) -> i32 {
        let mut total: i32 = 17;
        total = total.wrapping_mul(37).wrapping_add(self.0.start());
        total = total
            .wrapping_mul(37)
            .wrapping_add(allele_hash(self.0.ref_allele()));
        total
            .wrapping_mul(37)
            .wrapping_add(allele_hash(self.0.alt_allele()))
    }
}

/// `Byte` as a key: its value, signed.
#[derive(Debug, Clone, Copy, PartialEq)]
struct ByteKey(u8);

impl JavaHashCode for ByteKey {
    fn java_hash_code(&self) -> i32 {
        self.0 as i8 as i32
    }
}

/// `int[]{count, badCount, assemblyBadCount}`.
type Counts = [i32; 3];

fn increment<K: JavaHashCode + PartialEq + Clone>(
    map: &mut JavaHashMap<K, Counts>,
    key: K,
    element: &ElementSummary,
) {
    if map.get(&key).is_none() {
        map.insert(key.clone(), [0; 3]);
    }
    let values = map.get_mut(&key).expect("just inserted");
    values[0] += 1;
    values[1] += i32::from(element.bad_pileup);
    values[2] += i32::from(element.bad_assembly);
}

/// The events of one map, through `Collectors.toMap` into a `HashMap<Event, int[]>`.
fn to_event_map<K: JavaHashCode + PartialEq + Clone>(
    map: &JavaHashMap<K, Counts>,
    make: impl Fn(&K) -> Result<Event, String>,
) -> Result<JavaHashMap<EventKey, Counts>, String> {
    let mut out = JavaHashMap::new();
    for (key, counts) in map.iter() {
        out.insert(EventKey(make(key)?), *counts);
    }
    Ok(out)
}

/// `goodAndBadPileupEvents(alignmentData, args, header, minBaseQualityScore)`: the good and the bad
/// events, each in its `HashSet`'s iteration order. `contig_bases` is the whole contig, for the
/// deletions' reference alleles.
pub fn good_and_bad_pileup_events(
    loci: &[LocusSummary],
    arguments: &PileupDetectionArguments,
    contig_bases: &[u8],
) -> Result<(Vec<Event>, Vec<Event>), String> {
    let mut good: JavaHashMap<EventKey, ()> = JavaHashMap::new();
    let mut bad: JavaHashMap<EventKey, ()> = JavaHashMap::new();
    if !arguments.use_pileup_detection {
        return Ok((Vec::new(), Vec::new()));
    }
    let threshold = arguments.active_region_phred_threshold;
    for (i, locus) in loci.iter().enumerate() {
        let mut only_track_deletions = false;
        if threshold > 0.0 && locus.activity_score < threshold {
            if arguments.detect_indels
                && i + 1 < loci.len()
                && loci[i + 1].activity_score > threshold
            {
                only_track_deletions = true;
            } else {
                continue;
            }
        }
        let mut snps: JavaHashMap<ByteKey, Counts> = JavaHashMap::new();
        let mut insertions: JavaHashMap<String, Counts> = JavaHashMap::new();
        let mut deletions: JavaHashMap<i32, Counts> = JavaHashMap::new();
        for element in &locus.elements {
            let snp = element.snp.filter(|_| !only_track_deletions);
            let insertion = element.insertion.as_ref().filter(|_| !only_track_deletions);
            if snp.is_none() && insertion.is_none() && element.deletion.is_none() {
                continue;
            }
            if let Some(base) = snp {
                increment(&mut snps, ByteKey(base), element);
            }
            if let Some(bases) = insertion {
                increment(
                    &mut insertions,
                    String::from_utf8_lossy(bases).into_owned(),
                    element,
                );
            }
            if let Some(length) = element.deletion {
                increment(&mut deletions, length as i32, element);
            }
        }
        let reference = || Allele::create(&[locus.ref_base], true).map_err(|e| format!("{e:?}"));
        let event = |r: Allele, a: Allele| {
            Event::new(&locus.contig, locus.position, r, a).map_err(|e| e.message())
        };
        let snp_events = to_event_map(&snps, |base| {
            event(
                reference()?,
                Allele::create(&[base.0], false).map_err(|e| format!("{e:?}"))?,
            )
        })?;
        let insertion_events = to_event_map(&insertions, |bases| {
            let mut alt = vec![locus.ref_base];
            alt.extend_from_slice(bases.as_bytes());
            event(
                reference()?,
                Allele::create(&alt, false).map_err(|e| format!("{e:?}"))?,
            )
        })?;
        let deletion_events = to_event_map(&deletions, |length| {
            let from = locus.position as usize - 1;
            let to = (from + 1 + *length as usize).min(contig_bases.len());
            event(
                Allele::create(&contig_bases[from..to], true).map_err(|e| format!("{e:?}"))?,
                Allele::create(&[locus.ref_base], false).map_err(|e| format!("{e:?}"))?,
            )
        })?;
        for map in [&snp_events, &insertion_events, &deletion_events] {
            for (key, counts) in map.iter() {
                if passes_pileup_filters(arguments, counts, locus.depth, key.0.is_indel()) {
                    good.insert(key.clone(), ());
                }
                if fails_assembly_filters(arguments, counts) {
                    bad.insert(key.clone(), ());
                }
            }
        }
    }
    let events = |set: &JavaHashMap<EventKey, ()>| set.keys().map(|k| k.0.clone()).collect();
    Ok((events(&good), events(&bad)))
}

/// `passesPileupFilters`.
fn passes_pileup_filters(
    arguments: &PileupDetectionArguments,
    counts: &Counts,
    depth: i32,
    is_indel: bool,
) -> bool {
    let threshold = if is_indel {
        arguments.indel_threshold
    } else {
        arguments.snp_threshold
    };
    f64::from(counts[0]) / f64::from(depth) > threshold
        && f64::from(depth) >= arguments.pileup_absolute_depth
        && (arguments.bad_read_threshold <= 0.0
            || f64::from(counts[1]) / f64::from(counts[0]) <= arguments.bad_read_threshold)
}

/// `failsAssemblyFilters`.
fn fails_assembly_filters(arguments: &PileupDetectionArguments, counts: &Counts) -> bool {
    arguments.assembly_bad_read_threshold > 0.0
        && f64::from(counts[2]) / f64::from(counts[0]) >= arguments.assembly_bad_read_threshold
}

/// The loci of `region.getAlignmentData()`: those of the region's active span.
pub fn region_loci<'a>(loci: &'a [LocusSummary], region: &AssemblyRegion) -> &'a [LocusSummary] {
    let span = region.span();
    let from = loci.partition_point(|l| l.position < span.start);
    let to = loci.partition_point(|l| l.position <= span.end);
    &loci[from..to.max(from)]
}

/// The trimmer's events in `HAPLOTYPE_EVENT_COMPARATOR` order: the assembled ones less the bad
/// pileup events, then the good ones, each added only when the comparator has no equal yet.
pub fn events_for_trimming(assembled: &[Event], good: &[Event], bad: &[Event]) -> Vec<Event> {
    let mut events: Vec<Event> = assembled
        .iter()
        .filter(|event| !bad.contains(event))
        .cloned()
        .collect();
    events.extend(good.iter().cloned());
    let key = |e: &Event| {
        (
            e.start(),
            e.ref_allele().len(),
            e.alt_allele().display_string(),
        )
    };
    events.sort_by_key(key);
    events.dedup_by(|a, b| key(a) == key(b));
    events
}

/// `Locatable.withinDistanceOf`.
fn within_distance_of(event: &Event, other: &Event, distance: i32) -> bool {
    event.contig() == other.contig()
        && event.start() <= other.end() + distance
        && other.start() - distance <= event.end()
}

/// `Allele.compareTo`, which the `TreeSet<Haplotype>` orders by: the reference first, then the
/// bases as a string.
fn compare_haplotypes(a: &Haplotype, b: &Haplotype) -> std::cmp::Ordering {
    b.is_reference()
        .cmp(&a.is_reference())
        .then_with(|| a.bases().cmp(&b.bases()))
}

/// A `TreeSet<Haplotype>`'s `add`: in order, and nothing the comparator calls equal.
fn tree_set_add(set: &mut Vec<Haplotype>, haplotype: Haplotype) {
    match set.binary_search_by(|probe| compare_haplotypes(probe, &haplotype)) {
        Ok(_) => {}
        Err(at) => set.insert(at, haplotype),
    }
}

/// `removeHaplotypesWithBadAlleles`: every haplotype whose event map holds a bad event. A
/// haplotype with no map is kept, which is every haplotype of a freshly trimmed set.
pub fn remove_haplotypes_with_bad_alleles(
    assembly: &mut AssemblyResultSet,
    bad: &[Event],
    debug: bool,
) {
    if bad.is_empty() {
        return;
    }
    // Per bad event, every haplotype carrying it, so one carrying two is listed twice.
    let mut flagged: Vec<Haplotype> = Vec::new();
    for event in bad {
        for haplotype in assembly.haplotype_list() {
            if haplotype
                .event_map()
                .is_some_and(|map| map.events().any(|e| e == event))
            {
                flagged.push(haplotype.clone());
            }
        }
    }
    if !flagged.is_empty() && debug {
        crate::tool_stdout::println(&format!(
            "Found Assembly Haps with filtered Variants:\n{}",
            flagged
                .iter()
                .map(crate::pd_haplotypes::haplotype_string)
                .collect::<Vec<_>>()
                .join("\n")
        ));
    }
    for haplotype in &flagged {
        assembly.remove_haplotype(haplotype);
    }
}

/// `injectPileupEvents(region, args, aligner, goodPileupEvents)` on a trimmed set.
#[allow(clippy::too_many_arguments)]
pub fn inject_pileup_events(
    assembly: &mut AssemblyResultSet,
    region: &AssemblyRegion,
    arguments: &PileupDetectionArguments,
    aligner: &dyn SmithWatermanAligner,
    haplotype_to_reference_sw: &SwParameters,
    good: &[Event],
    max_mnp_distance: i32,
) -> Result<(), String> {
    if good.is_empty() {
        return Ok(());
    }
    let error = |e: ResultSetError| format!("{e:?}");
    let reference = assembly
        .reference_haplotype()
        .cloned()
        .ok_or_else(|| "no reference haplotype".to_string())?;
    let assembled = assembly.variation_events(max_mnp_distance).map_err(error)?;
    let assembled_indels: Vec<&Event> = assembled.iter().filter(|e| e.is_indel()).collect();

    // The reference first, then the best scores, five of them, into a `TreeSet`.
    let mut ranked: Vec<&Haplotype> = assembly.haplotype_list().iter().collect();
    ranked.sort_by(|a, b| {
        b.is_reference()
            .cmp(&a.is_reference())
            .then_with(|| b.score().total_cmp(&a.score()))
    });
    let mut base_haplotypes: Vec<Haplotype> = Vec::new();
    for haplotype in ranked.into_iter().take(NUM_HAPLOTYPES_TO_INJECT_INTO) {
        tree_set_add(&mut base_haplotypes, haplotype.clone());
    }

    let kmer_size = arguments.filtering_kmer_size.max(0) as usize;
    let read_kmers = kmer_read_counts(region.hard_clipped_pileup_reads(), kmer_size);

    let mut ordered: Vec<&Event> = good.iter().collect();
    ordered.sort_by_key(|e| e.start());
    let debug = arguments.debug_stdout;
    for event in ordered {
        if debug {
            crate::tool_stdout::println(&format!(
                "Processing new Haplotypes for Pileup Allele that was not in the assembly: {}",
                crate::pd_haplotypes::event_string(event)
            ));
        }
        if !event.is_indel()
            && assembled_indels.iter().any(|indel| {
                within_distance_of(event, indel, arguments.snp_adjacent_to_assembly_indel)
            })
        {
            continue;
        }
        if assembled
            .iter()
            .any(|e| e.start() == event.start() && e == event)
            || is_symbolic(event.alt_allele())
        {
            continue;
        }
        let mut new_haplotypes: Vec<Haplotype> = Vec::new();
        for base in &base_haplotypes {
            if let Some(inserted) = make_haplotype_with_inserted_event(
                base,
                &reference,
                event,
                aligner,
                haplotype_to_reference_sw,
            )? {
                if !new_haplotypes.contains(&inserted) {
                    new_haplotypes.push(inserted);
                }
            }
        }
        let kept = filter_pileup_haplotypes(
            new_haplotypes,
            &read_kmers,
            arguments.num_haplotypes_to_iterate,
            kmer_size,
        );
        if debug {
            crate::tool_stdout::println(&format!(
                "Constructed the following new Pileup Haplotypes after filtering:\n{}",
                kept.iter()
                    .map(crate::pd_haplotypes::haplotype_string)
                    .collect::<Vec<_>>()
                    .join("\n")
            ));
        }
        for haplotype in kept {
            tree_set_add(&mut base_haplotypes, haplotype);
        }
    }
    for haplotype in base_haplotypes {
        assembly.add(haplotype).map_err(error)?;
    }
    assembly.variation_events(max_mnp_distance).map_err(error)?;
    Ok(())
}

/// `isSymbolic`: a no-call, `-`, `*`, or a symbolic allele.
pub(crate) fn is_symbolic(allele: &Allele) -> bool {
    let text = allele.display_string();
    allele.is_no_call() || text == "-" || text == "*" || allele.is_symbolic()
}

/// `makeHaplotypeWithInsertedEvent`: none when the base haplotype already has an event over this
/// one (a haplotype made here has no event map, so it never does) or `insertAllele` refuses.
pub(crate) fn make_haplotype_with_inserted_event(
    base: &Haplotype,
    reference: &Haplotype,
    event: &Event,
    aligner: &dyn SmithWatermanAligner,
    haplotype_to_reference_sw: &SwParameters,
) -> Result<Option<Haplotype>, String> {
    if let Some(map) = base.event_map() {
        if map.events().any(|e| within_distance_of(event, e, 0)) {
            return Ok(None);
        }
    }
    let Some(bases) = insert_allele(base, event)? else {
        return Ok(None);
    };
    let mut inserted = Haplotype::non_reference(&bases).map_err(|e| format!("{e:?}"))?;
    let cigar = calculate_cigar(
        &reference.bases(),
        &bases,
        aligner,
        haplotype_to_reference_sw,
        SwOverhangStrategy::Indel,
    )
    .map_err(|e| format!("{e:?}"))?
    .ok_or_else(|| "java.lang.NullPointerException".to_string())?;
    inserted.set_cigar(&cigar).map_err(|e| format!("{e:?}"))?;
    if let Some(location) = reference.genome_location() {
        inserted.set_genome_location(location.clone());
    }
    inserted.set_alignment_start_hap_wrt_ref(reference.alignment_start_hap_wrt_ref());
    Ok(Some(inserted))
}

/// `Haplotype.insertAllele(ref, alt, start)`: the new bases, or none.
fn insert_allele(base: &Haplotype, event: &Event) -> Result<Option<Vec<u8>>, String> {
    let bases = base.bases();
    let reference = event.ref_allele();
    let alt = event.alt_allele().display_string().into_bytes();
    let start = base
        .genome_location()
        .ok_or_else(|| "java.lang.NullPointerException".to_string())?
        .start;
    let location = event.start();
    // The deletion of the zeroth base, whose common base is outside the contig.
    if reference.len() > alt.len() && location == start - 1 {
        let size = reference.len() - alt.len();
        if size > bases.len() {
            return Ok(None);
        }
        return Ok(Some(bases[size..].to_vec()));
    }
    let cigar = base
        .cigar()
        .ok_or_else(|| "java.lang.NullPointerException".to_string())?;
    let (index, op) = read_index_for_reference_coordinate(start, cigar, location);
    if index == READ_INDEX_NOT_FOUND || !op.is_some_and(|op| op.consumes_read_bases()) {
        return Ok(None);
    }
    let index = index as usize;
    if index + reference.len() > bases.len() {
        return Ok(None);
    }
    let mut out = bases[..index].to_vec();
    out.extend_from_slice(&alt);
    out.extend_from_slice(&bases[index + reference.len()..]);
    Ok(Some(out))
}

/// `getKmerReadCounts`, as far as anything reads it: which k-mers some read holds.
fn kmer_read_counts(reads: &[BamRecord], kmer_size: usize) -> std::collections::HashSet<Vec<u8>> {
    let mut kmers = std::collections::HashSet::new();
    for read in reads {
        let bases = &read.read_bases;
        if kmer_size == 0 || bases.len() < kmer_size {
            continue;
        }
        for i in 0..=bases.len() - kmer_size {
            kmers.insert(bases[i..i + kmer_size].to_vec());
        }
    }
    kmers
}

/// `filterPileupHaplotypes`: at most `count`, by the number of their distinct k-mers some read
/// holds, ties kept as good as the last place and then ordered by bases.
fn filter_pileup_haplotypes(
    haplotypes: Vec<Haplotype>,
    read_kmers: &std::collections::HashSet<Vec<u8>>,
    count: i32,
    kmer_size: usize,
) -> Vec<Haplotype> {
    if haplotypes.len() as i64 <= i64::from(count) {
        return haplotypes;
    }
    let scores: Vec<i64> = haplotypes
        .iter()
        .map(|h| {
            let bases = h.bases();
            let mut seen = std::collections::HashSet::new();
            if kmer_size <= bases.len() {
                for i in 0..=bases.len() - kmer_size {
                    seen.insert(&bases[i..i + kmer_size]);
                }
            }
            seen.iter().filter(|k| read_kmers.contains(**k)).count() as i64
        })
        .collect();
    let mut sorted = scores.clone();
    sorted.sort_unstable_by(|a, b| b.cmp(a));
    let minimum = sorted[(count - 1).max(0) as usize];
    let mut kept: Vec<(i64, Haplotype)> = haplotypes
        .into_iter()
        .zip(scores)
        .filter(|(_, score)| *score >= minimum)
        .map(|(h, score)| (score, h))
        .collect();
    kept.sort_by(|a, b| b.0.cmp(&a.0).then_with(|| a.1.bases().cmp(&b.1.bases())));
    kept.into_iter()
        .take(count.max(0) as usize)
        .map(|(_, h)| h)
        .collect()
}

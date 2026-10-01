//! `StructuralVariantDiscoverer`: the alignments of locally assembled contigs, read as structural
//! variant breakpoints and written as a VCF.
//!
//! The tool walks a queryname-sorted file, gathers the consecutive lines of each contig, keeps
//! the configuration of alignments that best explains it (see [`crate::sv_contig_alignments`]),
//! and reads a contig left with exactly two alignments as a novel adjacency (see
//! [`crate::sv_novel_adjacency`]). Contigs that imply the same adjacency are merged into one
//! record carrying every contig as evidence. At the end the records are annotated, filtered by
//! mapping quality and alignment length, sorted, and written by `SVVCFWriter`, which honours
//! none of the engine's VCF writer arguments: no index, no MD5, no compression, whatever the
//! output's name.
//!
//! Contigs left with three or more alignments are complex events, which go through
//! `CpxVariantInterpreter`. That path is not carried, and a run that reaches it is refused as the
//! port's own limitation rather than answered without the records it would write.
//!
//! Ported from `org.broadinstitute.hellbender.tools.StructuralVariantDiscoverer`,
//! `org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SimpleNovelAdjacencyAndChimericAlignmentEvidence`,
//! `org.broadinstitute.hellbender.tools.spark.sv.discovery.AnnotatedVariantProducer`,
//! `org.broadinstitute.hellbender.tools.spark.sv.discovery.SVMappingQualityFilter`,
//! `org.broadinstitute.hellbender.tools.spark.sv.discovery.SVAlignmentLengthFilter`,
//! `org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils`,
//! `org.broadinstitute.hellbender.tools.spark.sv.utils.CNVInputReader`,
//! `org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter` and
//! `org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines`.

use crate::sv_contig_alignments::{
    overlap_on_contig, AlignedContig, AlignmentInterval, AlignmentSignature, ContigRead,
    Dictionary, Interval, SvError, ILLEGAL_ARGUMENT, ILLEGAL_STATE, USER_EXCEPTION,
};
use crate::sv_novel_adjacency::{
    reverse_complement, split_pair_strong_enough, NovelAdjacency, SimpleChimera, SvKind,
    SvRecordType,
};
use std::collections::BTreeMap;

/// `gatk_rs::PortLimitation`, the class a refusal of the complex-event path is reported under.
pub const PORT_LIMITATION: &str = crate::main_entry::PORT_LIMITATION;

/// `StructuralVariationDiscoveryArgumentCollection.STRUCTURAL_VARIANT_SIZE_LOWER_BOUND`.
const SIZE_LOWER_BOUND: i32 = 50;
/// `CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD`, which `HQ_MAPPINGS` counts.
const HIGH_MQ: i32 = 60;
/// `StructuralVariantDiscoverer.SCORE_DIFF_TOLERANCE`.
const SCORE_DIFF_TOLERANCE: f64 = 0.0;

/// The message the tool refuses anything but a queryname-sorted file with.
pub const NOT_QUERYNAME_SORTED: &str = "This tool requires a queryname-sorted source of reads.";

/// The sort order a reads header declares.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortOrder {
    Queryname,
    Coordinate,
    Unsorted,
}

impl SortOrder {
    /// `SAMFileHeader.getSortOrder`, from the `@HD` `SO` value: absent is unsorted.
    pub fn from_header_value(value: Option<&str>) -> SortOrder {
        match value {
            Some("queryname") => SortOrder::Queryname,
            Some("coordinate") => SortOrder::Coordinate,
            _ => SortOrder::Unsorted,
        }
    }
}

/// `onTraversalStart`'s first refusal: the tool gathers a contig by walking CONSECUTIVE lines of
/// one name, which only a queryname-sorted file keeps together.
pub fn check_sort_order(order: SortOrder) -> Result<(), SvError> {
    if order == SortOrder::Queryname {
        Ok(())
    } else {
        Err(SvError::new(USER_EXCEPTION, NOT_QUERYNAME_SORTED))
    }
}

/// The tool's own read filters, `MAPPED` and `NOT_SECONDARY_ALIGNMENT`, as a predicate over the
/// two flags they read.
pub fn passes_default_read_filters(unmapped: bool, secondary: bool) -> bool {
    !unmapped && !secondary
}

/// A `HashSet<String>`'s iteration order: by bucket in the final table, then by insertion.
///
/// `SVUtils.getSampleId` prints the set of samples it found when there is not exactly one, and
/// that set is a `HashSet`, so the order it prints is `String.hashCode`'s rather than the header's.
fn java_hash_set_order(values: &[Option<String>]) -> Vec<Option<String>> {
    let mut unique: Vec<Option<String>> = Vec::new();
    for value in values {
        if !unique.contains(value) {
            unique.push(value.clone());
        }
    }
    let mut capacity = 16usize;
    while unique.len() as f64 > capacity as f64 * 0.75 {
        capacity *= 2;
    }
    let hash = |value: &Option<String>| -> u32 {
        let h = match value {
            None => 0i32,
            Some(text) => text
                .encode_utf16()
                .fold(0i32, |h, unit| h.wrapping_mul(31).wrapping_add(unit as i32)),
        } as u32;
        h ^ (h >> 16)
    };
    let mut ordered: Vec<(usize, usize, Option<String>)> = unique
        .into_iter()
        .enumerate()
        .map(|(position, value)| ((hash(&value) as usize) & (capacity - 1), position, value))
        .collect();
    ordered.sort_by_key(|(bucket, position, _)| (*bucket, *position));
    ordered.into_iter().map(|(_, _, value)| value).collect()
}

/// `SVUtils.getSampleId`: the one sample the read groups carry, which `Utils.validate` insists on.
pub fn sample_id(read_group_samples: &[Option<String>]) -> Result<Option<String>, SvError> {
    let set = java_hash_set_order(read_group_samples);
    if set.len() != 1 {
        let listed: Vec<String> = set
            .iter()
            .map(|sample| sample.clone().unwrap_or_else(|| "null".to_string()))
            .collect();
        return Err(SvError::new(
            ILLEGAL_STATE,
            format!(
                "Read groups must contain reads from one and only one sample, but we are finding the following ones in the given header: \t[{}]",
                listed.join(", ")
            ),
        ));
    }
    Ok(set.into_iter().next().flatten())
}

/// `SVUtils.getCanonicalChromosomes`: every contig of the dictionary, less the names a file lists
/// one per line.
pub fn canonical_chromosomes(
    dictionary: &Dictionary,
    non_canonical_file: Option<(&str, Option<&str>)>,
) -> Result<Vec<String>, SvError> {
    let mut contigs = dictionary.names.clone();
    let Some((path, text)) = non_canonical_file else {
        return Ok(contigs);
    };
    let Some(text) = text else {
        return Err(SvError::new(
            USER_EXCEPTION,
            format!("Can't read nonCanonicalContigNamesFile file {path}"),
        ));
    };
    // `Files.lines` splits on any of the three line terminators.
    for line in text.split('\n').flat_map(|line| line.split('\r')) {
        contigs.retain(|name| name != line);
    }
    Ok(contigs)
}

/// One record of the external copy-number calls, reduced to what the annotation reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CnvCall {
    pub contig_index: i32,
    pub start: i32,
    pub end: i32,
    pub id: String,
    pub copy_number: Option<String>,
    pub copy_number_quality: Option<String>,
}

/// `CNVInputReader.loadCNVCalls`: a single-sample VCF over the reads' own dictionary, held in an
/// interval tree keyed on `[start, end)` where a second record with the same key replaces the
/// first.
pub fn load_cnv_calls(
    text: &str,
    sample: Option<&str>,
    dictionary: &Dictionary,
) -> Result<Vec<CnvCall>, SvError> {
    let mut samples: Vec<String> = Vec::new();
    let mut contigs: Vec<(String, Option<i32>)> = Vec::new();
    let mut body = Vec::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("##contig=<") {
            let fields = rest.trim_end_matches('>');
            let mut id = None;
            let mut length = None;
            for field in fields.split(',') {
                if let Some(value) = field.strip_prefix("ID=") {
                    id = Some(value.to_string());
                } else if let Some(value) = field.strip_prefix("length=") {
                    length = value.parse().ok();
                }
            }
            if let Some(id) = id {
                contigs.push((id, length));
            }
        } else if line.starts_with("##") {
            continue;
        } else if let Some(header) = line.strip_prefix('#') {
            samples = header.split('\t').skip(9).map(str::to_string).collect();
        } else if !line.is_empty() {
            body.push(line);
        }
    }
    if samples.len() != 1 {
        return Err(SvError::new(
            ILLEGAL_STATE,
            "CNV call VCF should be single sample",
        ));
    }
    if sample != Some(samples[0].as_str()) {
        return Err(SvError::new(
            ILLEGAL_STATE,
            format!(
                "CNV call VCF does not contain calls for sample {}",
                sample.unwrap_or("null")
            ),
        ));
    }
    if contigs.is_empty() {
        return Err(SvError::new(
            ILLEGAL_STATE,
            "CNV calls file does not have a valid sequence dictionary",
        ));
    }
    let same = contigs.len() == dictionary.names.len()
        && contigs
            .iter()
            .zip(dictionary.names.iter().zip(&dictionary.lengths))
            .all(|((name, length), (other, other_length))| {
                name == other && length.is_none_or(|length| length == *other_length)
            });
    if !same {
        return Err(SvError::new(
            ILLEGAL_STATE,
            "CNV calls file does not have the same sequence dictionary as the read evidence",
        ));
    }
    let mut calls: Vec<CnvCall> = Vec::new();
    for line in body {
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 10 {
            continue;
        }
        let start: i32 = columns[1].parse().unwrap_or(0);
        let end = columns[7]
            .split(';')
            .find_map(|field| field.strip_prefix("END="))
            .and_then(|value| value.parse().ok())
            .unwrap_or(start + columns[3].len() as i32 - 1);
        let keys: Vec<&str> = columns[8].split(':').collect();
        let values: Vec<&str> = columns[9].split(':').collect();
        let field = |key: &str| {
            keys.iter()
                .position(|candidate| *candidate == key)
                .and_then(|index| values.get(index))
                .map(|value| value.to_string())
        };
        let call = CnvCall {
            contig_index: dictionary.index(columns[0]),
            start,
            end,
            id: columns[2].to_string(),
            copy_number: field("CN"),
            copy_number_quality: field("CNQ"),
        };
        if call.contig_index < 0 || call.start < 0 || call.end < 0 {
            return Err(SvError::new(
                ILLEGAL_ARGUMENT,
                if call.contig_index < 0 {
                    format!("provided contig is negative: {}", call.contig_index)
                } else if call.start < 0 {
                    format!("provided start is negative: {}", call.start)
                } else {
                    format!("provided end is negative: {}", call.end)
                },
            ));
        }
        match calls.iter_mut().find(|other| {
            (other.contig_index, other.start, other.end)
                == (call.contig_index, call.start, call.end)
        }) {
            Some(other) => *other = call,
            None => calls.push(call),
        }
    }
    calls.sort_by_key(|call| (call.contig_index, call.start, call.end));
    Ok(calls)
}

/// The arguments the tool reads for itself.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiscoveryArguments {
    /// `--min-mq`.
    pub min_mq: i32,
    /// `--min-align-length`.
    pub min_align_length: i32,
}

impl Default for DiscoveryArguments {
    fn default() -> Self {
        DiscoveryArguments {
            min_mq: 30,
            min_align_length: 50,
        }
    }
}

/// A `VariantContext` as the writer sees it: the eight columns, with the INFO fields keyed for the
/// `TreeMap` the encoder sorts them into and a flag held as the empty string.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VariantRecord {
    pub contig: String,
    pub start: i32,
    pub stop: i32,
    pub id: String,
    pub reference: String,
    pub alternate: String,
    pub filters: Vec<String>,
    pub attributes: BTreeMap<String, String>,
}

impl VariantRecord {
    /// `SvType.getBasicInformation`.
    fn from_type(sv: &SvRecordType) -> VariantRecord {
        let mut attributes = BTreeMap::new();
        let stop = if sv.kind == SvKind::BreakEnd {
            sv.start
        } else {
            attributes.insert("END".to_string(), sv.stop.to_string());
            attributes.insert("SVLEN".to_string(), sv.sv_len.to_string());
            sv.stop
        };
        attributes.insert("SVTYPE".to_string(), sv.sv_type().to_string());
        for (key, value) in &sv.extra {
            attributes.insert(key.to_string(), value.clone());
        }
        VariantRecord {
            contig: sv.chromosome.clone(),
            start: sv.start,
            stop,
            id: sv.id.clone(),
            reference: sv.reference_allele.clone(),
            alternate: sv.alternate_allele.clone(),
            filters: Vec::new(),
            attributes,
        }
    }

    /// `SVUtils.getAttributeAsStringList`: the value split at every comma.
    fn attribute_list(&self, key: &str) -> Vec<String> {
        self.attributes
            .get(key)
            .map(|value| value.split(',').map(str::to_string).collect())
            .unwrap_or_default()
    }

    /// `VariantContext.getAttributeAsInt`, with the default where the key is absent.
    fn attribute_as_int(&self, key: &str, default: i32) -> i32 {
        self.attributes
            .get(key)
            .and_then(|value| value.parse().ok())
            .unwrap_or(default)
    }

    /// One line of the VCF body, as `VCFEncoder` writes a site with no samples.
    fn to_line(&self) -> String {
        let info: Vec<String> = self
            .attributes
            .iter()
            .map(|(key, value)| {
                if value.is_empty() {
                    key.clone()
                } else {
                    format!("{key}={value}")
                }
            })
            .collect();
        format!(
            "{}\t{}\t{}\t{}\t{}\t.\t{}\t{}",
            self.contig,
            self.start,
            self.id,
            self.reference,
            self.alternate,
            if self.filters.is_empty() {
                ".".to_string()
            } else {
                self.filters.join(";")
            },
            if info.is_empty() {
                ".".to_string()
            } else {
                info.join(";")
            }
        )
    }
}

/// `getAssemblyEvidenceRelatedAnnotations`: one entry per contig, in contig-name order.
fn evidence_annotations(evidence: &[SimpleChimera]) -> Vec<(&'static str, String)> {
    let mut sorted: Vec<&SimpleChimera> = evidence.iter().collect();
    sorted.sort_by(|one, two| one.source_contig_name.cmp(&two.source_contig_name));
    // `ChimericContigAlignmentEvidenceAnnotations`.
    let annotations: Vec<(i32, i32, &SimpleChimera)> = sorted
        .iter()
        .map(|chimera| {
            let min_mq = chimera.lower.map_qual.min(chimera.higher.map_qual);
            let min_al = chimera
                .lower
                .reference_span
                .size()
                .min(chimera.higher.reference_span.size())
                - overlap_on_contig(&chimera.lower, &chimera.higher);
            (min_mq, min_al, *chimera)
        })
        .collect();
    let join = |values: Vec<String>| values.join(",");
    let mut attributes = vec![
        ("TOTAL_MAPPINGS", annotations.len().to_string()),
        (
            "HQ_MAPPINGS",
            annotations
                .iter()
                .filter(|(mq, ..)| *mq == HIGH_MQ)
                .count()
                .to_string(),
        ),
        (
            "MAPPING_QUALITIES",
            join(annotations.iter().map(|(mq, ..)| mq.to_string()).collect()),
        ),
        (
            "ALIGN_LENGTHS",
            join(
                annotations
                    .iter()
                    .map(|(_, al, _)| al.to_string())
                    .collect(),
            ),
        ),
        (
            "MAX_ALIGN_LENGTH",
            annotations
                .iter()
                .map(|(_, al, _)| *al)
                .max()
                .unwrap_or(0)
                .to_string(),
        ),
        (
            "CTG_NAMES",
            join(
                annotations
                    .iter()
                    .map(|(.., chimera)| chimera.source_contig_name.clone())
                    .collect(),
            ),
        ),
    ];
    let mut insertion_mappings: Vec<String> = annotations
        .iter()
        .flat_map(|(.., chimera)| chimera.insertion_mappings.iter().cloned())
        .collect();
    insertion_mappings.sort();
    if !insertion_mappings.is_empty() {
        attributes.push(("INSSEQ_MAP", join(insertion_mappings)));
    }
    let non_canonical: Vec<String> = annotations
        .iter()
        .map(|(.., chimera)| chimera.non_canonical_sa_tag.clone())
        .filter(|tag| {
            tag != crate::sv_contig_alignments::NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME
        })
        .collect();
    if !non_canonical.is_empty() {
        attributes.push(("CTG_GOOD_NONCANONICAL_MAPPING", join(non_canonical)));
    }
    attributes
}

/// `annotateWithExternalCNVCalls`: the calls whose half-open span overlaps the record's.
fn annotate_with_cnv_calls(
    record: &mut VariantRecord,
    sv: &SvRecordType,
    dictionary: &Dictionary,
    cnv_calls: Option<&[CnvCall]>,
) -> Result<(), SvError> {
    let Some(calls) = cnv_calls else {
        return Ok(());
    };
    // `new SVInterval(contig, pos, end)` refuses a negative coordinate, and a breakend's stop is
    // `NO_APPLICABLE_END`, so a breakend always refuses here.
    let contig = dictionary.index(&sv.chromosome);
    if contig < 0 {
        return Err(SvError::new(
            ILLEGAL_ARGUMENT,
            format!("provided contig is negative: {contig}"),
        ));
    }
    if sv.start < 0 {
        return Err(SvError::new(
            ILLEGAL_ARGUMENT,
            format!("provided start is negative: {}", sv.start),
        ));
    }
    if sv.stop < 0 {
        return Err(SvError::new(
            ILLEGAL_ARGUMENT,
            format!("provided end is negative: {}", sv.stop),
        ));
    }
    let annotation: Vec<String> = calls
        .iter()
        .filter(|call| call.contig_index == contig && call.start < sv.stop && sv.start < call.end)
        .map(|call| {
            format!(
                "{}:{}:{}",
                call.id,
                call.copy_number.as_deref().unwrap_or("null"),
                call.copy_number_quality.as_deref().unwrap_or("null")
            )
        })
        .collect();
    if !annotation.is_empty() {
        record
            .attributes
            .insert("EXTERNAL_CNV_CALLS".to_string(), annotation.join(","));
    }
    Ok(())
}

/// `produceAnnotatedVcFromAssemblyEvidence`.
fn annotated_record(
    sv: &SvRecordType,
    adjacency: &NovelAdjacency,
    evidence: &[SimpleChimera],
    dictionary: &Dictionary,
    cnv_calls: Option<&[CnvCall]>,
) -> Result<VariantRecord, SvError> {
    let mut record = VariantRecord::from_type(sv);
    for (key, value) in adjacency.complication.to_variant_attributes() {
        record.attributes.insert(key.to_string(), value);
    }
    for (key, value) in evidence_annotations(evidence) {
        record.attributes.insert(key.to_string(), value);
    }
    if sv.kind != SvKind::BreakEnd && !adjacency.alt_haplotype.is_empty() {
        record.attributes.insert(
            "SEQ_ALT_HAPLOTYPE".to_string(),
            String::from_utf8_lossy(&adjacency.alt_haplotype).into_owned(),
        );
    }
    annotate_with_cnv_calls(&mut record, sv, dictionary, cnv_calls)?;
    Ok(record)
}

/// `toVariantContexts`: one record, or a pair linked by `MATEID` (breakends) or `LINK`.
fn to_variant_records(
    types: &[SvRecordType],
    adjacency: &NovelAdjacency,
    evidence: &[SimpleChimera],
    dictionary: &Dictionary,
    cnv_calls: Option<&[CnvCall]>,
) -> Result<Vec<VariantRecord>, SvError> {
    if types.is_empty() || types.len() > 2 {
        return Err(SvError::new(
            crate::sv_contig_alignments::GATK_EXCEPTION,
            "Wrong number of variants sent for analysis",
        ));
    }
    if types.len() == 1 {
        return Ok(vec![annotated_record(
            &types[0], adjacency, evidence, dictionary, cnv_calls,
        )?]);
    }
    let link = if types[0].kind == SvKind::BreakEnd {
        "MATEID"
    } else {
        "LINK"
    };
    let mut first = annotated_record(&types[0], adjacency, evidence, dictionary, cnv_calls)?;
    let mut second = annotated_record(&types[1], adjacency, evidence, dictionary, cnv_calls)?;
    first.attributes.insert(link.to_string(), second.id.clone());
    second.attributes.insert(link.to_string(), first.id.clone());
    let strip = |record: &mut VariantRecord| {
        for key in ["INSSEQ", "INSLEN", "SEQ_ALT_HAPLOTYPE"] {
            record.attributes.remove(key);
        }
    };
    if types[0].kind == SvKind::Deletion {
        strip(&mut first);
    } else if types[1].kind == SvKind::Deletion {
        strip(&mut second);
    }
    Ok(vec![first, second])
}

/// `AnnotatedVariantProducer.filterMergedVariantList`: small simple variants dropped, and the two
/// assembly filters applied to what is left.
fn filter_merged_variant_list(
    variants: Vec<VariantRecord>,
    arguments: &DiscoveryArguments,
) -> Result<Vec<VariantRecord>, SvError> {
    let mut kept = Vec::with_capacity(variants.len());
    for mut variant in variants {
        let sv_type = variant
            .attributes
            .get("SVTYPE")
            .cloned()
            .unwrap_or_default();
        if matches!(sv_type.as_str(), "DEL" | "INS" | "DUP")
            && variant.attribute_as_int("SVLEN", 0).abs() < SIZE_LOWER_BOUND
        {
            continue;
        }
        let mut applied = Vec::new();
        if variant.attributes.contains_key("CTG_NAMES") {
            // `SVMappingQualityFilter`: the best of the evidence's mapping qualities.
            let mut max_mq = 0;
            for value in variant.attribute_list("MAPPING_QUALITIES") {
                let quality: i32 = value.parse().map_err(|_| {
                    SvError::new(
                        "java.lang.NumberFormatException",
                        format!("For input string: \"{value}\""),
                    )
                })?;
                max_mq = max_mq.max(quality);
            }
            if max_mq < arguments.min_mq {
                applied.push("LOW_MQ".to_string());
            }
            // `SVAlignmentLengthFilter`.
            if variant.attribute_as_int("MAX_ALIGN_LENGTH", 0) < arguments.min_align_length {
                applied.push("SHORT_ALN".to_string());
            }
        }
        // The encoder sorts the filters it writes.
        applied.sort();
        variant.filters = applied;
        kept.push(variant);
    }
    Ok(kept)
}

/// `processContigAlignments`: one contig's lines, interpreted, and what they add to the evidence.
fn process_contig(
    reads: &[ContigRead],
    dictionary: &Dictionary,
    canonical: &[String],
    evidence: &mut Vec<(NovelAdjacency, Vec<SimpleChimera>)>,
) -> Result<(), SvError> {
    let mut alignments = Vec::with_capacity(reads.len());
    let mut name = String::new();
    let mut sequence: Option<Vec<u8>> = None;
    for read in reads {
        name = read.name.clone();
        if !read.supplementary {
            sequence = Some(if read.reverse_strand {
                reverse_complement(&read.bases)
            } else {
                read.bases.clone()
            });
        }
        alignments.push(AlignmentInterval::from_read(read)?);
    }
    let Some(sequence) = sequence else {
        return Err(SvError::new(
            USER_EXCEPTION,
            format!("No primary line for {name}"),
        ));
    };
    let contig = AlignedContig::new(&name, sequence, alignments);
    if !contig.has_good_mq() {
        return Ok(());
    }
    for fine_tuned in contig.reconstruct_from_best_configuration(canonical, SCORE_DIFF_TOLERANCE)? {
        match fine_tuned.signature() {
            AlignmentSignature::SimpleChimera => {
                let head = &fine_tuned.contig.alignments[0];
                let tail = &fine_tuned.contig.alignments[fine_tuned.contig.alignments.len() - 1];
                if !split_pair_strong_enough(head, tail) {
                    continue;
                }
                let chimera = SimpleChimera::new(
                    fine_tuned.contig.alignments[0].clone(),
                    fine_tuned.contig.alignments[1].clone(),
                    fine_tuned.insertion_mappings.clone(),
                    &fine_tuned.contig.name,
                    &fine_tuned.non_canonical_sa_tag,
                    dictionary,
                )?;
                let adjacency =
                    NovelAdjacency::new(&chimera, &fine_tuned.contig.sequence, dictionary)?;
                match evidence.iter_mut().find(|(key, _)| *key == adjacency) {
                    Some((_, chimeras)) => chimeras.push(chimera),
                    None => evidence.push((adjacency, vec![chimera])),
                }
            }
            AlignmentSignature::Complex => {
                return Err(SvError::new(
                    PORT_LIMITATION,
                    format!(
                        "Contig {} aligns as a complex event, which StructuralVariantDiscoverer interprets through CpxVariantInterpreter; that path is a GATK feature this port does not carry yet. This message is the port's own and not GATK's.",
                        fine_tuned.contig.name
                    ),
                ));
            }
            AlignmentSignature::Normal | AlignmentSignature::Unknown => {}
        }
    }
    Ok(())
}

/// The whole tool from the reads that reached `apply` to the records `SVVCFWriter` is handed,
/// in the order it writes them.
///
/// `reference` answers `ReferenceContext.getBases(window)`, trimming included. The evidence map
/// is walked in insertion order where the reference walks a `HashMap`; the records are sorted
/// before they are written, so only which of two refusals comes first could tell the two apart.
pub fn discover(
    reads: &[ContigRead],
    dictionary: &Dictionary,
    canonical: &[String],
    cnv_calls: Option<&[CnvCall]>,
    arguments: &DiscoveryArguments,
    reference: &mut dyn FnMut(&Interval) -> Result<Vec<u8>, SvError>,
) -> Result<Vec<VariantRecord>, SvError> {
    let mut evidence: Vec<(NovelAdjacency, Vec<SimpleChimera>)> = Vec::new();
    let mut start = 0;
    while start < reads.len() {
        let mut end = start + 1;
        while end < reads.len() && reads[end].name == reads[start].name {
            end += 1;
        }
        process_contig(&reads[start..end], dictionary, canonical, &mut evidence)?;
        start = end;
    }
    let mut variants = Vec::new();
    for (adjacency, chimeras) in &evidence {
        let types = adjacency.to_simple_or_bnd_types(reference)?;
        variants.extend(to_variant_records(
            &types, adjacency, chimeras, dictionary, cnv_calls,
        )?);
    }
    let mut kept = filter_merged_variant_list(variants, arguments)?;
    sort_variants(&mut kept, dictionary)?;
    Ok(kept)
}

/// `SVVCFWriter.sortVariantsByCoordinate`: position, then inserted sequence, then ID, stably.
fn sort_variants(variants: &mut [VariantRecord], dictionary: &Dictionary) -> Result<(), SvError> {
    for variant in variants.iter() {
        if dictionary.index(&variant.contig) < 0 {
            return Err(SvError::new(
                ILLEGAL_ARGUMENT,
                "Can't do comparison because Locatables' contigs not found in sequence dictionary",
            ));
        }
    }
    variants.sort_by(|one, two| {
        dictionary
            .index(&one.contig)
            .cmp(&dictionary.index(&two.contig))
            .then(one.start.cmp(&two.start))
            .then(one.stop.cmp(&two.stop))
            .then_with(|| {
                let inserted = |record: &VariantRecord| {
                    record.attributes.get("INSSEQ").cloned().unwrap_or_default()
                };
                inserted(one).cmp(&inserted(two))
            })
            .then_with(|| one.id.cmp(&two.id))
    });
    Ok(())
}

/// One structured header line of `GATKSVVCFHeaderLines`: key, ID, number and type where the key
/// has them, and the description.
struct HeaderLine {
    key: &'static str,
    id: &'static str,
    number: Option<&'static str>,
    kind: Option<&'static str>,
    description: &'static str,
}

const fn info(
    id: &'static str,
    number: &'static str,
    kind: &'static str,
    description: &'static str,
) -> HeaderLine {
    HeaderLine {
        key: "INFO",
        id,
        number: Some(number),
        kind: Some(kind),
        description,
    }
}

const fn simple(key: &'static str, id: &'static str, description: &'static str) -> HeaderLine {
    HeaderLine {
        key,
        id,
        number: None,
        kind: None,
        description,
    }
}

/// `GATKSVVCFHeaderLines`: the symbolic alleles, INFO, FORMAT and FILTER lines, and the standard
/// `END` line `SVVCFWriter.getVcfHeader` adds. `LINK` is declared twice there and the second
/// description is the one the map keeps.
const HEADER_LINES: &[HeaderLine] = &[
    simple("ALT", "INV", "Inversion of reference sequence"),
    simple("ALT", "DEL", "Deletion relative to the reference"),
    simple("ALT", "INS", "Insertion of novel sequence relative to the reference"),
    simple("ALT", "DUP", "Region of elevated copy number relative to the reference"),
    simple("ALT", "DUP:INV", "Region of elevated copy number relative to the reference, with some copies inverted"),
    simple("ALT", "CPX", "Complex rearrangement of reference sequence"),
    info("SVTYPE", "1", "String", "Type of structural variant"),
    info("SVLEN", ".", "Integer", "Difference in length between REF and ALT alleles"),
    info("SEQ_ALT_HAPLOTYPE", "A", "Character", "Alt haplotype sequence, one per alt allele"),
    info("INSSEQ", ".", "String", "Inserted sequence at the breakpoint"),
    info("INSLEN", "A", "Integer", "Length of inserted sequence (note for duplication records, this does not count the extra copies of the duplicated sequence)"),
    info("INSSEQ_MAP", ".", "String", "Alignments of inserted sequence"),
    info("HOMSEQ", ".", "String", "Sequence of base pair identical micro-homology at event breakpoints"),
    info("HOMLEN", "1", "Integer", "Length of base pair identical micro-homology at event breakpoints"),
    info("READ_PAIR_SUPPORT", "1", "Integer", "Number of discordant read pairs supporting the variant"),
    info("SPLIT_READ_SUPPORT", "1", "Integer", "Number of split read supplementary mappings supporting the variant"),
    info("EXTERNAL_CNV_CALLS", "1", "String", "Comma-delimited list of external copy number calls that overlap with this variant in format ID:CN:CNQ"),
    info("CTG_NAMES", ".", "String", "Name of local assembly contigs supporting this variant, formatted as \"asm%06d:tig%05d\""),
    info("TOTAL_MAPPINGS", "1", "Integer", "Number of local assembly contigs supporting the variant, i.e. number of entries in CTG_NAMES"),
    info("MAPPING_QUALITIES", ".", "Integer", "Mapping qualities of the contig alignments that support the variant"),
    info("HQ_MAPPINGS", "1", "Integer", "Number of high-quality contig alignments that support the variant"),
    info("ALIGN_LENGTHS", ".", "Integer", "Minimum lengths of the flanking aligned region from each contig alignment"),
    info("MAX_ALIGN_LENGTH", "1", "Integer", "Maximum of the values listed in ALIGN_LENGTHS"),
    info("CTG_GOOD_NONCANONICAL_MAPPING", ".", "String", "Good mapping of evidence contigs as listed in CTG_NAMES to non-canonical reference chromosomes that could potentially offer better explanation of the assembly contig without the SV record. One for each evidence assembly contig, if available; otherwise a \".\". If no evidence contig has such mapping, this annotation is omitted for the record."),
    info("LINK", ".", "String", "ID(s) of other variants that are linked to this record, i.e. they constitute a larger more complex variant"),
    info("MATEID", ".", "String", "ID(s) for mate(s) of a BND record"),
    info("IMPRECISE", "0", "Flag", "Imprecise structural variation"),
    info("CIPOS", "2", "Integer", "Confidence interval around POS for imprecise variants"),
    info("CIEND", "2", "Integer", "Confidence interval around END for imprecise variants"),
    info("INV33", "0", "Flag", "Whether the event represents a 3' to 5' breakpoint"),
    info("INV55", "0", "Flag", "Whether the event represents a 5' to 3' breakpoint"),
    info("DUP_REPEAT_UNIT_REF_SPAN", "1", "String", "Reference span of the suspected repeated unit in a tandem duplication"),
    info("DUP_SEQ_CIGARS", ".", "String", "CIGARs of the repeated sequence on the locally-assembled contigs when aligned to DUP_REPEAT_UNIT_REF_SPAN (currently only available for repeats when DUP_ANNOTATIONS_IMPRECISE is false)"),
    info("DUP_NUM", "R", "Integer", "Number of times the sequence is duplicated on reference and on the alternate alleles"),
    info("DUP_ANNOTATIONS_IMPRECISE", "0", "Flag", "Whether the duplication annotations are from an experimental optimization procedure"),
    info("DUP_IMPRECISE_AFFECTED_RANGE", ".", "String", "Affected reference range for duplications annotated with the flag DUP_ANNOTATIONS_IMPRECISE"),
    info("DUP_ORIENTATIONS", "A", "String", "Orientations of the duplicated sequence on alt allele relative to the copy on ref; one group for each alt allele (currently only available for inverted duplication variants)"),
    info("CONTRACTION", "0", "Flag", "Tandem repeats contraction compared to reference"),
    info("EXPANSION", "0", "Flag", "Tandem repeats expansion compared to reference"),
    info("ALT_ARRANGEMENT", ".", "String", "For CPX variants only; specifies how reference segments given in SEGMENTS are re-arranged"),
    info("SEGMENTS", ".", "String", "For CPX variants only; segments of reference that are rearranged"),
    info("CPX_EVENT", ".", "String", "ID(s) of CPX(s) events from which current simple variant record is extracted"),
    info("END", "1", "Integer", "Stop position of the interval"),
    HeaderLine {
        key: "FORMAT",
        id: "CN",
        number: Some("1"),
        kind: Some("Integer"),
        description: "Copy number genotype for imprecise events",
    },
    HeaderLine {
        key: "FORMAT",
        id: "CNQ",
        number: Some("1"),
        kind: Some("Float"),
        description: "Copy number genotype quality for imprecise events",
    },
    simple("FILTER", "LOW_MQ", "Assembly evidence based record that whose maximum value specified in MAPPING_QUALITIES is lower than user specified threshold"),
    simple("FILTER", "SHORT_ALN", "Assembly evidence based record that whose MAPPING_QUALITIES value is lower than user specified threshold"),
    simple("FILTER", "LOW_QS", "Depth-only copy number record whose QS value is lower than user specified threshold"),
    simple("FILTER", "FREQ", "Depth-only copy number record whose AF value is higher than user specified threshold"),
];

/// `SVVCFWriter.writeVCF`: the header, sorted as `VCFWriter` sorts it, and the records.
///
/// `getMetaDataInSortedOrder` orders the lines by their text, except that contig lines compare
/// among themselves by dictionary index. `default_lines` are the tool's own `##source` and
/// `##GATKCommandLine` lines, already rendered, which join the sort.
pub fn write_vcf(
    variants: &[VariantRecord],
    dictionary: &Dictionary,
    default_lines: &[String],
) -> String {
    // (text without the leading `##`, contig index for a contig line)
    let mut lines: Vec<(String, Option<usize>)> = HEADER_LINES
        .iter()
        .map(|line| {
            let description = line.description.replace('"', "\\\"");
            let text = match (line.number, line.kind) {
                (Some(number), Some(kind)) => format!(
                    "{}=<ID={},Number={},Type={},Description=\"{}\">",
                    line.key, line.id, number, kind, description
                ),
                _ => format!(
                    "{}=<ID={},Description=\"{}\">",
                    line.key, line.id, description
                ),
            };
            (text, None)
        })
        .collect();
    for (index, (name, length)) in dictionary.names.iter().zip(&dictionary.lengths).enumerate() {
        let assembly = dictionary
            .assemblies
            .get(index)
            .cloned()
            .flatten()
            .map(|assembly| format!(",assembly={assembly}"))
            .unwrap_or_default();
        lines.push((
            format!("contig=<ID={name},length={length}{assembly}>"),
            Some(index),
        ));
    }
    lines.extend(
        default_lines
            .iter()
            .map(|line| (line.trim_start_matches("##").to_string(), None)),
    );
    // `VCFHeaderLine.compareTo` is the text, and `VCFContigHeaderLine`'s the index when both
    // sides are contig lines.
    lines.sort_by(|one, two| match (one.1, two.1) {
        (Some(first), Some(second)) => first.cmp(&second),
        _ => one.0.cmp(&two.0),
    });
    let mut out = String::from("##fileformat=VCFv4.2\n");
    for (line, _) in lines {
        out.push_str("##");
        out.push_str(&line);
        out.push('\n');
    }
    out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for variant in variants {
        out.push_str(&variant.to_line());
        out.push('\n');
    }
    out
}

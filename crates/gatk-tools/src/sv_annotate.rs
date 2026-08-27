//! `SVAnnotate`: what a structural variant is predicted to do to a gene.
//!
//! The answer is one word per transcript, and which word it is depends on the SV type, on which
//! features the two breakpoints land in, and above all on whether the variant SPANS a feature or
//! merely overlaps it. The duplication rule is built entirely out of that difference.
//!
//! Reading the GTF, the BED and the VCF is not ported, nor is the interval tree that makes the
//! lookup fast: the transcripts a variant reaches are given to the rules directly. Every rule that
//! turns a variant and a transcript into a word is.

/// The consequence words, as `GATKSVVCFConstants` spells them in the output.
pub const LOF: &str = "PREDICTED_LOF";
pub const INT_EXON_DUP: &str = "PREDICTED_INTRAGENIC_EXON_DUP";
pub const DUP_PARTIAL: &str = "PREDICTED_DUP_PARTIAL";
pub const PARTIAL_EXON_DUP: &str = "PREDICTED_PARTIAL_EXON_DUP";
pub const COPY_GAIN: &str = "PREDICTED_COPY_GAIN";
pub const TSS_DUP: &str = "PREDICTED_TSS_DUP";
pub const INV_SPAN: &str = "PREDICTED_INV_SPAN";
pub const MSV_EXON_OVERLAP: &str = "PREDICTED_MSV_EXON_OVERLAP";
pub const UTR: &str = "PREDICTED_UTR";
pub const INTRONIC: &str = "PREDICTED_INTRONIC";
pub const BREAKEND_EXON: &str = "PREDICTED_BREAKEND_EXONIC";
pub const PARTIAL_DISPERSED_DUP: &str = "PREDICTED_PARTIAL_DISPERSED_DUP";
pub const PROMOTER: &str = "PREDICTED_PROMOTER";
pub const NEAREST_TSS: &str = "PREDICTED_NEAREST_TSS";
pub const NONCODING_SPAN: &str = "PREDICTED_NONCODING_SPAN";
pub const NONCODING_BREAKPOINT: &str = "PREDICTED_NONCODING_BREAKPOINT";
pub const INTERGENIC: &str = "PREDICTED_INTERGENIC";

/// The six duplication answers a multiallelic CNV has reclassified to `MSV_EXON_OVERLAP`.
pub const MSV_EXON_OVERLAP_CLASSIFICATIONS: &[&str] = &[
    LOF,
    INT_EXON_DUP,
    DUP_PARTIAL,
    PARTIAL_EXON_DUP,
    COPY_GAIN,
    TSS_DUP,
];

/// The consequences that make a variant NOT intergenic. `PARTIAL_DISPERSED_DUP` is deliberately
/// absent: a dispersed duplication over a gene leaves the variant intergenic.
pub const PROTEIN_CODING_CONSEQUENCES: &[&str] = &[
    LOF,
    INT_EXON_DUP,
    DUP_PARTIAL,
    PARTIAL_EXON_DUP,
    COPY_GAIN,
    TSS_DUP,
    INV_SPAN,
    MSV_EXON_OVERLAP,
    UTR,
    INTRONIC,
    BREAKEND_EXON,
];

/// `GATKSVVCFConstants.StructuralVariantAnnotationType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SvType {
    Del,
    Dup,
    Cnv,
    Ins,
    Inv,
    Bnd,
    Ctx,
    Cpx,
}

/// The complex subtypes the annotation rules ask about by name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ComplexSubtype {
    DDup,
    DupInv,
    InvDup,
    DupInvDup,
    DupInvDel,
    DelInvDup,
    DDupIDel,
    InsIDel,
    DelInv,
    InvDel,
    DelInvDel,
    CtxPpQq,
    CtxPqQp,
    CtxInv,
}

/// The subtypes every DUP segment of which is dispersed rather than tandem.
pub const COMPLEX_SUBTYPES_WITH_DISPERSED_DUP: &[ComplexSubtype] = &[
    ComplexSubtype::DDup,
    ComplexSubtype::DupInv,
    ComplexSubtype::InvDup,
    ComplexSubtype::DupInvDup,
    ComplexSubtype::DupInvDel,
    ComplexSubtype::DelInvDup,
    ComplexSubtype::DDupIDel,
];

/// `SimpleInterval`, closed on both ends.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interval {
    pub contig: String,
    pub start: i32,
    pub end: i32,
}

impl Interval {
    /// The constructor's own validation, which is where two of this tool's arguments end up.
    pub fn new(contig: &str, start: i32, end: i32) -> Result<Interval, AnnotateError> {
        if start > end || start < 1 {
            return Err(AnnotateError::InvalidInterval {
                contig: contig.to_string(),
                start,
                end,
            });
        }
        Ok(Interval {
            contig: contig.to_string(),
            start,
            end,
        })
    }

    pub fn overlaps(&self, other: &Interval) -> bool {
        self.contig == other.contig && self.start <= other.end && other.start <= self.end
    }

    pub fn contains(&self, other: &Interval) -> bool {
        self.contig == other.contig && self.start <= other.start && other.end <= self.end
    }
}

/// What this tool refuses, and what it accepts and then crashes on.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AnnotateError {
    CpxWithoutIntervals,
    CpxWithoutType,
    CtxWithoutContig2,
    /// `SimpleInterval`'s own message, which is what both the zero promoter window and the
    /// negative breakend length end up producing.
    InvalidInterval {
        contig: String,
        start: i32,
        end: i32,
    },
    /// `Integer.parseInt` on a BED header row, which the argument documentation asks for.
    NumberFormat {
        text: String,
    },
}

impl AnnotateError {
    pub fn message(&self) -> String {
        match self {
            AnnotateError::CpxWithoutIntervals => {
                "Complex (CPX) variant must contain CPX_INTERVALS INFO field".to_string()
            }
            AnnotateError::CpxWithoutType => {
                "Complex (CPX) variant must contain CPX_TYPE INFO field".to_string()
            }
            AnnotateError::CtxWithoutContig2 => {
                "Translocation (CTX) variant represented as a single record must contain CHR2 \
                 INFO field"
                    .to_string()
            }
            AnnotateError::InvalidInterval { contig, start, end } => {
                format!("Invalid interval. Contig:{contig} start:{start} end:{end}")
            }
            AnnotateError::NumberFormat { text } => format!("For input string: \"{text}\""),
        }
    }
}

/// The GTF feature types the rules read. Everything else in the file is walked past.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FeatureKind {
    Transcript,
    Exon,
    Cds,
    StartCodon,
    StopCodon,
    Utr,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Feature {
    pub kind: FeatureKind,
    pub start: i32,
    pub end: i32,
}

/// One transcript, with its features in the order the GTF listed them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Transcript {
    pub gene_name: String,
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub negative_strand: bool,
    /// `getAllFeatures`, which includes the transcript record itself.
    pub features: Vec<Feature>,
}

impl Transcript {
    pub fn interval(&self) -> Interval {
        Interval {
            contig: self.contig.clone(),
            start: self.start,
            end: self.end,
        }
    }

    /// `getTranscriptionStartSite`: the GTF codec always reassigns start to the earlier coordinate,
    /// so on the minus strand the start site is the transcript's END.
    pub fn transcription_start_site(&self) -> i32 {
        if self.negative_strand {
            self.end
        } else {
            self.start
        }
    }

    fn feature_interval(&self, feature: &Feature) -> Interval {
        Interval {
            contig: self.contig.clone(),
            start: feature.start,
            end: feature.end,
        }
    }
}

/// `getPromoterInterval`: the window upstream of the start site, which on the plus strand ends one
/// base BEFORE it. A window of zero therefore asks for an interval that ends before it begins, and
/// `--promoter-window-length` declares `minValue = 0`, so the argument accepts it.
pub fn promoter_interval(
    transcript: &Transcript,
    promoter_window: i32,
) -> Result<Interval, AnnotateError> {
    let tss = transcript.transcription_start_site();
    let (left, right) = if transcript.negative_strand {
        (tss + 1, tss + promoter_window)
    } else {
        ((tss - promoter_window).max(1), tss - 1)
    };
    Interval::new(&transcript.contig, left, right)
}

/// `variantSpansFeature`.
pub fn variant_spans_feature(variant: &Interval, feature: &Interval) -> bool {
    variant.contains(feature)
}

/// `countBreakendsInsideFeature`: 0, 1 or 2, and a variant that CONTAINS the feature has none
/// inside it however much they overlap.
pub fn count_breakends_inside_feature(variant: &Interval, feature: &Interval) -> i32 {
    if !feature.overlaps(variant) || variant.contains(feature) {
        0
    } else if feature.contains(variant) {
        2
    } else {
        1
    }
}

fn variant_overlaps_transcription_start_site(variant: &Interval, transcript: &Transcript) -> bool {
    let tss = transcript.transcription_start_site();
    variant.overlaps(&Interval {
        contig: transcript.contig.clone(),
        start: tss,
        end: tss,
    })
}

/// `getSimpleConsequence`: coding sequence wins outright, a UTR is remembered and can still be
/// overwritten by a later coding feature, and anything else leaves it intronic.
fn simple_consequence(variant: &Interval, transcript: &Transcript) -> &'static str {
    let mut consequence = INTRONIC;
    for feature in &transcript.features {
        if !variant.overlaps(&transcript.feature_interval(feature)) {
            continue;
        }
        match feature.kind {
            FeatureKind::Cds => return LOF,
            FeatureKind::Utr => consequence = UTR,
            _ => {}
        }
    }
    consequence
}

pub fn annotate_insertion(variant: &Interval, transcript: &Transcript) -> &'static str {
    simple_consequence(variant, transcript)
}

/// `annotateDeletion`: the start site is checked FIRST, so a deletion over it is LOF whatever else
/// it lands in.
pub fn annotate_deletion(variant: &Interval, transcript: &Transcript) -> &'static str {
    if variant_overlaps_transcription_start_site(variant, transcript) {
        LOF
    } else {
        simple_consequence(variant, transcript)
    }
}

/// `annotateDuplication`, which is the whole spanning-versus-overlapping question in one function.
pub fn annotate_duplication(
    variant: &Interval,
    transcript: &Transcript,
    is_dispersed_duplication: bool,
) -> &'static str {
    let transcript_interval = transcript.interval();
    if variant_spans_feature(variant, &transcript_interval) {
        // The same answer whether the duplication is tandem or dispersed, so it is returned before
        // the dispersed case is ever asked.
        return COPY_GAIN;
    }
    if is_dispersed_duplication {
        return PARTIAL_DISPERSED_DUP;
    }
    if variant_overlaps_transcription_start_site(variant, transcript) {
        return TSS_DUP;
    }
    if !transcript_interval.contains(variant) {
        // One breakpoint inside the transcript and one past its end.
        return DUP_PARTIAL;
    }

    // Both breakpoints inside the transcript: the answer is a count of where they landed.
    let mut breakpoints_in_cds = 0;
    let mut breakpoints_in_utr = 0;
    let mut cds_spanned = 0;
    let mut utr_spanned = 0;
    for feature in &transcript.features {
        let interval = transcript.feature_interval(feature);
        if !variant.overlaps(&interval) {
            continue;
        }
        match feature.kind {
            FeatureKind::Cds => {
                if variant_spans_feature(variant, &interval) {
                    cds_spanned += 1;
                } else {
                    breakpoints_in_cds += count_breakends_inside_feature(variant, &interval);
                }
            }
            FeatureKind::Utr => {
                if variant_spans_feature(variant, &interval) {
                    utr_spanned += 1;
                } else {
                    breakpoints_in_utr += count_breakends_inside_feature(variant, &interval);
                }
            }
            _ => {}
        }
    }
    if breakpoints_in_cds == 2 || (breakpoints_in_cds == 1 && breakpoints_in_utr == 1) {
        // The only place the UTR count is read for anything but the UTR answer.
        LOF
    } else if breakpoints_in_cds == 1 {
        PARTIAL_EXON_DUP
    } else if cds_spanned > 0 {
        INT_EXON_DUP
    } else if breakpoints_in_utr > 0 || utr_spanned > 0 {
        UTR
    } else {
        INTRONIC
    }
}

/// `annotateCopyNumberVariant`: annotated as a duplication and then reclassified, because the
/// consequence of a multiallelic CNV depends on the individual's copy number and cannot be decided
/// at the site.
pub fn annotate_copy_number_variant(variant: &Interval, transcript: &Transcript) -> &'static str {
    let consequence = annotate_duplication(variant, transcript, false);
    if MSV_EXON_OVERLAP_CLASSIFICATIONS.contains(&consequence) {
        MSV_EXON_OVERLAP
    } else {
        consequence
    }
}

/// `annotateInversion`: the deletion rule plus the spanning case.
pub fn annotate_inversion(variant: &Interval, transcript: &Transcript) -> &'static str {
    if variant_spans_feature(variant, &transcript.interval()) {
        INV_SPAN
    } else {
        annotate_deletion(variant, transcript)
    }
}

/// `annotateTranslocation`: called only with transcripts the variant reaches, and breaking a gene
/// at all is predicted to break it, so nothing else is looked at.
pub fn annotate_translocation(_variant: &Interval, _transcript: &Transcript) -> &'static str {
    LOF
}

/// `annotateBreakend`: the SIMPLE consequence, not the deletion one, and then LOF downgraded
/// because a low-confidence breakend should not be called loss of function.
pub fn annotate_breakend(variant: &Interval, transcript: &Transcript) -> &'static str {
    let consequence = simple_consequence(variant, transcript);
    if consequence == LOF {
        BREAKEND_EXON
    } else {
        consequence
    }
}

/// `annotateTranscript`. An insertion, deletion, duplication, CNV, inversion, translocation or
/// breakend gets a word; a complex variant gets none at this level, its segments having been
/// resolved into the types above already.
pub fn annotate_transcript(
    variant: &Interval,
    sv_type: SvType,
    includes_dispersed_duplication: bool,
    transcript: &Transcript,
) -> Option<&'static str> {
    match sv_type {
        SvType::Del => Some(annotate_deletion(variant, transcript)),
        SvType::Ins => Some(annotate_insertion(variant, transcript)),
        SvType::Dup => Some(annotate_duplication(
            variant,
            transcript,
            includes_dispersed_duplication,
        )),
        SvType::Cnv => Some(annotate_copy_number_variant(variant, transcript)),
        SvType::Inv => Some(annotate_inversion(variant, transcript)),
        SvType::Ctx => Some(annotate_translocation(variant, transcript)),
        SvType::Bnd => Some(annotate_breakend(variant, transcript)),
        SvType::Cpx => None,
    }
}

/// One piece of a variant, with the type it is annotated as rather than the record's own.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SvSegment {
    pub sv_type: SvType,
    pub interval: Interval,
}

/// The fields of one record the segment rules read.
#[derive(Debug, Clone, PartialEq)]
pub struct Variant {
    pub id: String,
    pub contig: String,
    pub position: i32,
    pub end: i32,
    pub sv_type: SvType,
    pub sv_length: i32,
    pub contig2: Option<String>,
    pub end2: Option<i32>,
    pub strands: Option<String>,
    pub complex_type: Option<ComplexSubtype>,
    pub complex_intervals: Vec<String>,
}

/// `BND_DELETION_STRANDS` and `BND_DUPLICATION_STRANDS`.
pub const BND_DELETION_STRANDS: &str = "+-";
pub const BND_DUPLICATION_STRANDS: &str = "-+";

/// `getAnnotationTypeForBreakend`.
///
/// The length test is `svLen <= maxBreakendLen`, so the conventional SVLEN of -1 is judged to be
/// under ANY maximum, and the interval built from it below ends one base before it starts.
pub fn annotation_type_for_breakend(variant: &Variant, max_breakend_len: i32) -> SvType {
    if matches!(
        variant.complex_type,
        Some(ComplexSubtype::CtxPpQq) | Some(ComplexSubtype::CtxPqQp)
    ) {
        return SvType::Ctx;
    }
    let same_contig = variant.contig2.as_deref() == Some(variant.contig.as_str());
    if max_breakend_len > 0 && same_contig && variant.sv_length <= max_breakend_len {
        match variant.strands.as_deref() {
            // Not enough information to say which, so it stays a breakend.
            None => return SvType::Bnd,
            Some(BND_DELETION_STRANDS) => return SvType::Del,
            Some(BND_DUPLICATION_STRANDS) => return SvType::Dup,
            Some(_) => {}
        }
    }
    SvType::Bnd
}

/// `getSVSegments`. A breakend yields TWO segments, which is why one record can carry two
/// different consequences for the same gene.
pub fn sv_segments(
    variant: &Variant,
    max_breakend_len: i32,
) -> Result<Vec<SvSegment>, AnnotateError> {
    let contig = variant.contig.as_str();
    let pos = variant.position;
    let end2 = variant.end2.unwrap_or(pos);
    match variant.sv_type {
        SvType::Cpx => {
            if variant.complex_intervals.is_empty() {
                return Err(AnnotateError::CpxWithoutIntervals);
            }
            if variant.complex_type.is_none() {
                return Err(AnnotateError::CpxWithoutType);
            }
            let mut out = complex_annotation_intervals(
                &parse_complex_intervals(&variant.complex_intervals)?,
                variant.complex_type,
            );
            if variant.complex_type == Some(ComplexSubtype::DDup) {
                // The sink site, which a dDUP carries in CHROM and POS rather than in the
                // intervals field.
                out.push(SvSegment {
                    sv_type: SvType::Ins,
                    interval: Interval::new(contig, pos, pos + 1)?,
                });
            }
            Ok(out)
        }
        SvType::Ctx => {
            let mut out = if variant.complex_type == Some(ComplexSubtype::CtxInv)
                && !variant.complex_intervals.is_empty()
            {
                complex_annotation_intervals(
                    &parse_complex_intervals(&variant.complex_intervals)?,
                    variant.complex_type,
                )
            } else {
                Vec::new()
            };
            // POS and END separately, in case END is not POS + 1.
            out.push(SvSegment {
                sv_type: SvType::Ctx,
                interval: Interval::new(contig, pos, pos)?,
            });
            out.push(SvSegment {
                sv_type: SvType::Ctx,
                interval: Interval::new(contig, variant.end, variant.end)?,
            });
            let Some(contig2) = variant.contig2.as_deref() else {
                return Err(AnnotateError::CtxWithoutContig2);
            };
            out.push(SvSegment {
                sv_type: SvType::Ctx,
                interval: Interval::new(contig2, end2, end2)?,
            });
            out.push(SvSegment {
                sv_type: SvType::Ctx,
                interval: Interval::new(contig2, end2 + 1, end2 + 1)?,
            });
            Ok(out)
        }
        SvType::Bnd => {
            let annotate_as = annotation_type_for_breakend(variant, max_breakend_len);
            if matches!(annotate_as, SvType::Del | SvType::Dup) {
                return Ok(vec![SvSegment {
                    sv_type: annotate_as,
                    // With SVLEN of -1 this is pos..pos-1, which the interval refuses.
                    interval: Interval::new(contig, pos, pos + variant.sv_length)?,
                }]);
            }
            let mut out = vec![SvSegment {
                sv_type: annotate_as,
                interval: Interval::new(contig, pos, pos)?,
            }];
            match variant.contig2.as_deref() {
                Some(contig2) if contig2 == contig => {
                    // Whichever of SVLEN and END2 the record carries.
                    if variant.sv_length > 0 {
                        out.push(SvSegment {
                            sv_type: annotate_as,
                            interval: Interval::new(
                                contig,
                                pos + variant.sv_length,
                                pos + variant.sv_length,
                            )?,
                        });
                    } else if end2 != pos {
                        out.push(SvSegment {
                            sv_type: annotate_as,
                            interval: Interval::new(contig, end2, end2)?,
                        });
                    }
                }
                Some(contig2) => out.push(SvSegment {
                    sv_type: annotate_as,
                    interval: Interval::new(contig2, end2, end2)?,
                }),
                None => {}
            }
            Ok(out)
        }
        SvType::Ins => Ok(vec![SvSegment {
            sv_type: SvType::Ins,
            interval: Interval::new(contig, pos, pos + 1)?,
        }]),
        other => Ok(vec![SvSegment {
            sv_type: other,
            interval: Interval::new(contig, pos, variant.end)?,
        }]),
    }
}

/// `parseComplexIntervals`: each item reads `SVTYPE_CHROM:POS-END`.
pub fn parse_complex_intervals(items: &[String]) -> Result<Vec<SvSegment>, AnnotateError> {
    let mut out = Vec::new();
    for item in items {
        let (type_text, locus) = item.split_once('_').expect("a type and a locus");
        let (contig, range) = locus.split_once(':').expect("a contig and a range");
        let (start, end) = range.split_once('-').expect("a start and an end");
        let parse = |text: &str| {
            text.parse::<i32>()
                .map_err(|_| AnnotateError::NumberFormat {
                    text: text.to_string(),
                })
        };
        out.push(SvSegment {
            sv_type: match type_text {
                "DEL" => SvType::Del,
                "DUP" => SvType::Dup,
                "INV" => SvType::Inv,
                "INS" => SvType::Ins,
                "CNV" => SvType::Cnv,
                "CTX" => SvType::Ctx,
                _ => SvType::Bnd,
            },
            interval: Interval::new(contig, parse(start)?, parse(end)?)?,
        });
    }
    Ok(out)
}

/// `getComplexAnnotationIntervals`, reduced to what the segments carry here: the intervals are
/// taken as they were parsed.
pub fn complex_annotation_intervals(
    segments: &[SvSegment],
    _complex_type: Option<ComplexSubtype>,
) -> Vec<SvSegment> {
    segments.to_vec()
}

/// `getSegmentsForNonCodingAnnotations`: a complex variant's DUP segments are always dispersed, so
/// they are not considered for the promoter or the non-coding elements.
pub fn segments_for_non_coding(segments: &[SvSegment], is_complex: bool) -> Vec<SvSegment> {
    if is_complex {
        segments
            .iter()
            .filter(|segment| segment.sv_type != SvType::Dup)
            .cloned()
            .collect()
    } else {
        segments.to_vec()
    }
}

/// `isIntergenic`: true when nothing protein-coding was found. `PARTIAL_DISPERSED_DUP` is not in
/// that set, so a dispersed duplication over a gene leaves the variant intergenic.
pub fn is_intergenic(consequences: &[(String, Vec<String>)]) -> bool {
    !consequences
        .iter()
        .any(|(consequence, _)| PROTEIN_CODING_CONSEQUENCES.contains(&consequence.as_str()))
}

pub fn includes_dispersed_duplication(complex_type: Option<ComplexSubtype>) -> bool {
    complex_type.is_some_and(|kind| COMPLEX_SUBTYPES_WITH_DISPERSED_DUP.contains(&kind))
}

/// One non-coding element.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NonCodingElement {
    pub name: String,
    pub interval: Interval,
}

/// The running consequence map, kept as a sorted association list so the output order is the
/// output order.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Consequences {
    entries: Vec<(String, Vec<String>)>,
}

impl Consequences {
    /// `updateVariantConsequenceDict`: a set per consequence, so a name is never repeated.
    pub fn add(&mut self, consequence: &str, name: &str) {
        match self.entries.iter_mut().find(|(key, _)| key == consequence) {
            Some((_, names)) => {
                if !names.iter().any(|existing| existing == name) {
                    names.push(name.to_string());
                }
            }
            None => self
                .entries
                .push((consequence.to_string(), vec![name.to_string()])),
        }
    }

    pub fn contains(&self, consequence: &str) -> bool {
        self.entries.iter().any(|(key, _)| key == consequence)
    }

    /// Every gene named by any consequence, which is what the promoter rule checks against.
    fn named_genes(&self) -> Vec<String> {
        self.entries
            .iter()
            .flat_map(|(_, names)| names.iter().cloned())
            .collect()
    }

    /// `sortVariantConsequenceDict` plus the writer's own key order: the names are sorted inside
    /// each consequence, and the consequences come out alphabetically because the VCF writer sorts
    /// the INFO keys.
    pub fn sorted(&self) -> Vec<(String, Vec<String>)> {
        let mut out: Vec<(String, Vec<String>)> = self
            .entries
            .iter()
            .map(|(key, names)| {
                let mut sorted = names.clone();
                sorted.sort();
                (key.clone(), sorted)
            })
            .collect();
        out.sort();
        out
    }
}

/// What one annotated record carries.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Annotation {
    pub id: String,
    pub consequences: Vec<(String, Vec<String>)>,
    /// Absent when no GTF was given, in which case the flag is not written at all.
    pub intergenic: Option<bool>,
}

/// `annotateStructuralVariant`, in the order it works: the gene overlaps first, then whether it is
/// intergenic, then the promoter and the non-coding elements, then the nearest start site.
pub fn annotate_structural_variant(
    variant: &Variant,
    transcripts: &[Transcript],
    non_coding: &[NonCodingElement],
    have_gtf: bool,
    have_bed: bool,
    promoter_window: i32,
    max_breakend_len: i32,
) -> Result<Annotation, AnnotateError> {
    let mut consequences = Consequences::default();
    let complex_type = variant.complex_type;
    let dispersed = includes_dispersed_duplication(complex_type);
    let segments = sv_segments(variant, max_breakend_len)?;

    if have_gtf {
        for segment in &segments {
            for transcript in transcripts {
                if !segment.interval.overlaps(&transcript.interval()) {
                    continue;
                }
                if let Some(consequence) =
                    annotate_transcript(&segment.interval, segment.sv_type, dispersed, transcript)
                {
                    consequences.add(consequence, &transcript.gene_name);
                }
            }
        }
    }

    // Read BEFORE the promoter and the non-coding elements are added, so neither can change it.
    let intergenic = is_intergenic(&consequences.sorted());

    let non_coding_segments = segments_for_non_coding(&segments, dispersed);

    if have_gtf {
        // The promoter is inferred, and only when the gene has no coding consequence already.
        for segment in &non_coding_segments {
            let named = consequences.named_genes();
            for transcript in transcripts {
                let promoter = promoter_interval(transcript, promoter_window)?;
                if !segment.interval.overlaps(&promoter) {
                    continue;
                }
                if !named.contains(&transcript.gene_name) {
                    consequences.add(PROMOTER, &transcript.gene_name);
                }
            }
        }
    }

    if have_bed {
        for segment in &non_coding_segments {
            for element in non_coding {
                if !segment.interval.overlaps(&element.interval) {
                    continue;
                }
                let consequence = if variant_spans_feature(&segment.interval, &element.interval) {
                    NONCODING_SPAN
                } else {
                    NONCODING_BREAKPOINT
                };
                consequences.add(consequence, &element.name);
            }
        }
    }

    // The nearest start site, for an intergenic variant that reached no promoter.
    if have_gtf && !consequences.contains(PROMOTER) && intergenic {
        for segment in segments_for_nearest_tss(&non_coding_segments, complex_type) {
            if let Some(gene) = nearest_transcription_start_site(&segment.interval, transcripts) {
                consequences.add(NEAREST_TSS, &gene);
            }
        }
    }

    Ok(Annotation {
        id: variant.id.clone(),
        consequences: consequences.sorted(),
        intergenic: have_gtf.then_some(intergenic),
    })
}

/// `getSegmentForNearestTSS`: the subtypes whose remaining segments are merged into one, so that a
/// complex event gets a single nearest start site from its outer breakpoints.
pub fn segments_for_nearest_tss(
    segments: &[SvSegment],
    complex_type: Option<ComplexSubtype>,
) -> Vec<SvSegment> {
    let merged = matches!(
        complex_type,
        Some(ComplexSubtype::InsIDel)
            | Some(ComplexSubtype::DDupIDel)
            | Some(ComplexSubtype::DelInv)
            | Some(ComplexSubtype::InvDel)
            | Some(ComplexSubtype::DupInvDel)
            | Some(ComplexSubtype::DelInvDup)
            | Some(ComplexSubtype::DelInvDel)
    );
    if !merged || segments.is_empty() {
        return segments.to_vec();
    }
    let mut span = segments[0].interval.clone();
    for segment in &segments[1..] {
        span.start = span.start.min(segment.interval.start);
        span.end = span.end.max(segment.interval.end);
    }
    vec![SvSegment {
        sv_type: SvType::Del,
        interval: span,
    }]
}

/// `annotateNearestTranscriptionStartSite`: the nearest start site on the SAME contig.
///
/// The reference asks the interval tree for the greatest site at or below the variant and the
/// least site at or above it, then keeps the one with the smaller gap, comparing with a strict
/// `<`. That comparison decides a tie in favour of the site AFTER the variant, which is what the
/// `before` half of the sort key reproduces here.
///
/// The measured fixture separates its two candidates by more than ten thousand bases either way,
/// so what the golden pins is the ORDERING, not the exact gap: an off-by-one in the arithmetic
/// below would not show up in it.
pub fn nearest_transcription_start_site(
    variant: &Interval,
    transcripts: &[Transcript],
) -> Option<String> {
    let mut best: Option<(i32, bool, String)> = None;
    for transcript in transcripts {
        if transcript.contig != variant.contig {
            continue;
        }
        let tss = transcript.transcription_start_site();
        // The site is held as the half-open interval [tss, tss + 1), so the gap to a variant
        // starting after it is measured from tss + 1.
        let (gap, before) = if tss < variant.start {
            (variant.start - (tss + 1), true)
        } else if variant.end < tss {
            (tss - variant.end - 1, false)
        } else {
            (0, false)
        };
        let candidate = (gap, before, transcript.gene_name.clone());
        best = match best {
            Some(current) if (current.0, current.1) <= (candidate.0, candidate.1) => Some(current),
            _ => Some(candidate),
        };
    }
    best.map(|(_, _, gene)| gene)
}

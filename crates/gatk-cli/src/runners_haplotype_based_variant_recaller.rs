//! `HaplotypeBasedVariantRecaller`, in a file of its own.
//!
//! The groups, the events, the merged alleles, the PairHMM engine and the matrix lines are
//! [`gatk_tools::haplotype_based_variant_recaller`]; what is here is the walk over the alleles
//! file and the reads and haplotypes each record reaches.
use super::*;

/// A read's unclipped start and end, soft and hard clips both counted.
fn unclipped_span(read: &BamRecord) -> (i32, i32) {
    use htsjdk_bam::cigar::Op;
    let trailing: i32 = read
        .cigar
        .elements
        .iter()
        .rev()
        .take_while(|e| matches!(e.op, Op::S | Op::H))
        .map(|e| e.length as i32)
        .sum();
    (
        gatk_engine::read_utils::unclipped_start(read),
        gatk_engine::read_utils::end(read) + trailing,
    )
}

/// `HaplotypeBasedVariantRecaller`: every allele of a VCF, as the haplotypes around it carry it,
/// scored against every read that spans it.
///
/// Per VCF record in each interval: the haplotype group that centres it best, the reads that
/// contain it hard clipped to the group's span, the PairHMM likelihood of each read against each
/// haplotype (the Java `LOGLESS_CACHING` kernel; the native ones are not ported), the reads the
/// engine disqualifies taken away, the haplotypes' events merged at each start the span holds,
/// the likelihoods marginalized onto the merged alleles, and one block per start that is the
/// record's own.
pub fn haplotype_based_variant_recaller(parser: &Parser) -> Outcome {
    use gatk_engine::interval::SimpleInterval;
    use gatk_tools::haplotype_based_variant_recaller as hbvr;
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "HaplotypeBasedVariantRecaller")?;

    // traverse(): the engine first.
    let engine =
        scalar(parser, "likelihood-calculation-engine").unwrap_or_else(|| "PairHMM".to_string());
    let threshold = number_or(parser, "base-quality-score-threshold", 18);
    if engine == "PairHMM" {
        match scalar(parser, "pair-hmm-implementation").as_deref() {
            Some("LOGLESS_CACHING") => {}
            _ => {
                return Err(Thrown::non_user(
                    PORT_LIMITATION,
                    "only the LOGLESS_CACHING PairHMM is ported. This message is the port's own and not GATK's."
                        .to_string(),
                ))
            }
        }
        if threshold < 6 {
            return Err(Thrown::non_user(
                "java.lang.IllegalArgumentException",
                "baseQualityScoreThreshold must be greater than or equal to 6 (QualityUtils.MIN_USABLE_Q_SCORE)"
                    .to_string(),
            ));
        }
    }
    let Some(reference_path) = argument(parser, "reference") else {
        return Err(Thrown::non_user(
            "java.lang.IllegalArgumentException",
            "Null object is not allowed here.".to_string(),
        ));
    };
    let mut reference =
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&reference_path))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let output = scalar(parser, "matrix-file-csv").unwrap_or_default();
    write_file(&output, b"")?;
    let alleles_path = scalar(parser, "alleles-file-vcf").unwrap_or_default();
    let (variants, _) = feature_variants(&alleles_path)?;
    let haplotypes_path = scalar(parser, "haplotypes-file-bam").unwrap_or_default();
    let haplotype_source = {
        let path = std::path::Path::new(&haplotypes_path);
        match htsjdk_bam::sam_files::find_index(path) {
            Some(index) => ReadsDataSource::open(path, &index),
            None => ReadsDataSource::open_unindexed(path),
        }
        .map_err(|error| Thrown::user(format!("{error:?}")))?
    };
    let filter = read_filter(parser, &filters, &header)?;
    let intervals: Vec<SimpleInterval> = if intervals.is_empty() {
        header
            .sequences
            .iter()
            .map(|sequence| SimpleInterval {
                contig: sequence.name.clone(),
                start: 1,
                end: sequence.length,
            })
            .collect()
    } else {
        intervals
    };
    let settings = hbvr::EngineSettings {
        pcr_error_model: scalar(parser, "pcr-indel-model")
            .unwrap_or_else(|| "CONSERVATIVE".to_string()),
        base_quality_score_threshold: threshold as u8,
        gap_continuation_penalty: number_or(parser, "pair-hmm-gap-continuation-penalty", 10) as u8,
        disable_cap_read_qualities_to_map_q: flag(
            parser,
            "disable-cap-base-qualities-to-map-quality",
        ),
        dynamic_disqualification: flag(
            parser,
            "enable-dynamic-read-disqualification-for-genotyping",
        ),
        read_disqualification_scale: double_or(
            parser,
            "dynamic-read-disqualification-threshold",
            1.0,
        ),
        expected_error_rate_per_base: double_or(
            parser,
            "expected-mismatch-rate-for-read-disqualification",
            0.02,
        ),
    };
    let max_mnp_distance = number_or(parser, "max-mnp-distance", 0);
    let emit_spanning_dels = !flag(parser, "disable-spanning-event-genotyping");
    let alleles_indexed = has_feature_index(&alleles_path);
    let sample_of = |read: &BamRecord| -> String {
        gatk_engine::read_group::resolve(read, &header)
            .and_then(|group| group.attributes.get("SM").map(str::to_string))
            .unwrap_or_else(|| "null".to_string())
    };
    let cigar_of = |read: &BamRecord| -> Vec<hbvr::CigarElement> {
        read.cigar
            .elements
            .iter()
            .map(|element| hbvr::CigarElement {
                operator: element.op.to_char() as char,
                length: element.length as i32,
            })
            .collect()
    };
    let read_end =
        |read: &BamRecord| read.alignment_start + read.cigar.reference_length() as i32 - 1;

    let mut text = String::new();
    for region in &intervals {
        if !alleles_indexed {
            return Err(Thrown::user(format!(
                "Input {alleles_path} must support random access to enable queries by interval. If \
                 it's a file, please index it using the bundled tool IndexFeatureFile"
            )));
        }
        for vc in variants.iter().filter(|vc| {
            vc.contig == region.contig
                && vc.start as i32 <= region.end
                && vc.stop as i32 >= region.start
        }) {
            let vc_loc = hbvr::Span::new(&vc.contig, vc.start as i32, vc.stop as i32);
            // forBest: the haplotype records over the variant, grouped by span.
            let records = haplotype_source
                .query(&[SimpleInterval {
                    contig: vc_loc.contig.clone(),
                    start: vc_loc.start,
                    end: vc_loc.end,
                }])
                .map_err(reads_traversal_error)?;
            let haplotype_records: Vec<&BamRecord> = records
                .iter()
                .filter(|record| hbvr::is_haplotype_record(&record.read_name))
                .collect();
            let listed: Vec<hbvr::HaplotypeRecord> = haplotype_records
                .iter()
                .map(|record| hbvr::HaplotypeRecord {
                    name: record.read_name.clone(),
                    span: hbvr::Span::new(&vc_loc.contig, record.alignment_start, read_end(record)),
                })
                .collect();
            let groups = hbvr::groups(&listed);
            let Some(best) = hbvr::best_group(&vc_loc, &groups) else {
                continue;
            };
            // The group's records, in the order the query returned them.
            let first_index = {
                let mut at = 0usize;
                for group in &groups {
                    if std::ptr::eq(group, best) {
                        break;
                    }
                    at += group.len();
                }
                at
            };
            let members: Vec<&BamRecord> =
                haplotype_records[first_index..first_index + best.len()].to_vec();
            let span = best[0].span.clone();
            let span_interval = SimpleInterval {
                contig: span.contig.clone(),
                start: span.start,
                end: span.end,
            };
            let ref_bases = reference
                .query(&span.contig, span.start, span.end)
                .map_err(|error| Thrown::user(format!("{error:?}")))?;

            // The reads that contain the variant, filtered, then clipped to the span.
            let candidates = source
                .query(std::slice::from_ref(&span_interval))
                .map_err(reads_traversal_error)?;
            let mut reads: Vec<BamRecord> = Vec::new();
            for record in &candidates {
                if record.flags & 0x4 != 0 {
                    continue;
                }
                if !(record.alignment_start <= vc_loc.start && read_end(record) >= vc_loc.end) {
                    continue;
                }
                if !filter(record) {
                    continue;
                }
                let clipped =
                    gatk_engine::clipping::hard_clip_soft_clipped_bases(record, Some(&header), 0)
                        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{error:?}")))?;
                let clipped = gatk_engine::clipping::hard_clip_to_region(
                    &clipped,
                    Some(&header),
                    span.start,
                    span.end,
                )
                .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{error:?}")))?;
                if clipped.flags & 0x4 != 0 || clipped.cigar.elements.is_empty() {
                    continue;
                }
                reads.push(clipped);
            }
            if engine != "PairHMM" {
                if let Some(read) = reads.first() {
                    return Err(Thrown::non_user(
                        "java.lang.IllegalArgumentException",
                        format!(
                            "read must be flow based: {} {}:{}-{}",
                            read.read_name,
                            span.contig,
                            read.alignment_start,
                            read_end(read)
                        ),
                    ));
                }
            }

            // The haplotypes, deduplicated by bases as the allele list keeps them.
            let mut distinct: Vec<Vec<u8>> = Vec::new();
            let mut index_of_member: Vec<usize> = Vec::new();
            for member in &members {
                let at = match distinct
                    .iter()
                    .position(|bases| *bases == member.read_bases)
                {
                    Some(at) => at,
                    None => {
                        distinct.push(member.read_bases.clone());
                        distinct.len() - 1
                    }
                };
                index_of_member.push(at);
            }
            let engine_reads: Vec<hbvr::EngineRead> = reads
                .iter()
                .map(|read| {
                    let tag = |name: &[u8; 2]| match read.tags.get(htsjdk_bam::tag::Tag::new(name))
                    {
                        Some(htsjdk_bam::tag::TagValue::Str(text)) => {
                            Some(text.bytes().map(|b| b.wrapping_sub(33)).collect())
                        }
                        _ => None,
                    };
                    hbvr::EngineRead {
                        bases: read.read_bases.clone(),
                        quals: read.base_qualities.clone(),
                        insertion_quals: tag(b"BI"),
                        deletion_quals: tag(b"BD"),
                        mapping_quality: read.mapping_quality,
                    }
                })
                .collect();
            let result = hbvr::compute_likelihoods(&distinct, &engine_reads, &settings);

            // simplifiedAssignGenotypeLikelihood.
            let mut maps = Vec::with_capacity(members.len());
            for member in &members {
                maps.push(
                    hbvr::event_map(
                        &member.read_bases,
                        &cigar_of(member),
                        &ref_bases,
                        span.start,
                        max_mnp_distance,
                    )
                    .map_err(|(class, message)| {
                        Thrown::non_user(Box::leak(class.into_boxed_str()), message)
                    })?,
                );
            }
            let starts: std::collections::BTreeSet<i32> =
                maps.iter().flat_map(|map| map.keys().copied()).collect();
            let kind = {
                let reference_length = vc.reference().len();
                let types: std::collections::BTreeSet<u8> = vc
                    .alternate_alleles()
                    .iter()
                    .map(|alt| {
                        if alt.is_symbolic() || vc.reference().is_symbolic() {
                            0
                        } else if alt.len() == reference_length {
                            if reference_length == 1 {
                                1
                            } else {
                                2
                            }
                        } else {
                            3
                        }
                    })
                    .collect();
                if types.len() > 1 {
                    hbvr::VariantKind::Mixed
                } else {
                    hbvr::VariantKind::Other
                }
            };
            for loc in starts {
                if !(span.start <= loc && loc <= span.end) {
                    continue;
                }
                let reference_base = ref_bases[(loc - span.start) as usize];
                let Some(merged) = hbvr::merged_alleles(&maps, loc, reference_base).map_err(
                    |(class, message)| Thrown::non_user(Box::leak(class.into_boxed_str()), message),
                )?
                else {
                    continue;
                };
                if loc != vc.start as i32 {
                    continue;
                }
                let mapper = hbvr::allele_mapper(&merged, &maps, loc, emit_spanning_dels);
                let header_line = hbvr::header_line(
                    &vc.contig,
                    vc.start as i32,
                    vc.stop as i32,
                    kind,
                    &span,
                    &mapper
                        .iter()
                        .map(|(allele, _)| allele.text())
                        .collect::<Vec<_>>(),
                );
                let mut lines = Vec::new();
                for (column, &r) in result.kept.iter().enumerate() {
                    let values: Vec<f64> = mapper
                        .iter()
                        .map(|(_, members_of)| {
                            members_of
                                .iter()
                                .map(|m| result.likelihoods[index_of_member[*m]][column])
                                .fold(f64::NEG_INFINITY, f64::max)
                        })
                        .collect();
                    let read = &reads[r];
                    let (unclipped_start, unclipped_end) = unclipped_span(read);
                    let line_read = hbvr::Read {
                        name: read.read_name.clone(),
                        span: hbvr::Span::new(&vc.contig, read.alignment_start, read_end(read)),
                        cigar: cigar_of(read),
                        bases: read.read_bases.clone(),
                        is_duplicate: read.flags & 0x400 != 0,
                        is_reverse: read.flags & 0x10 != 0,
                        mapping_quality: i32::from(read.mapping_quality),
                        key_length: 0,
                        sample: sample_of(read),
                        unclipped_start,
                        unclipped_end,
                    };
                    if let Some(line) = hbvr::matrix_line(&line_read, &vc_loc, &values) {
                        lines.push(line);
                    }
                }
                text.push_str(&hbvr::variant_block(&header_line, &lines));
            }
        }
    }
    write_file(&output, text.as_bytes())?;
    Ok(None)
}

//! `GroundTruthReadsBuilder`: every read scored against the haplotypes its two ancestral
//! references give it, and a CSV row per read that passes.
//!
//! The translation to the ancestral contigs and the row's constants are
//! [`gatk_tools::ground_truth_reads_builder`]; the score is `FlowFeatureMapper`'s
//! `computeLikelihoodLocal`, as for `GroundTruthScorer`. What is here is the `PartialReadWalker`
//! around them, in the reference's order:
//!
//! * **the traversal stops** at the first read after `--max-output-reads` reads were counted, and
//!   a read is counted before `--gt-no-output` decides whether it is written;
//! * **the cheap filters come first** (mapping quality, supplementary alignments, the soft-clip
//!   rule), then the subsampling draw, then the read quality, which is the sum over the flows of
//!   the probability of the maximal hmer;
//! * **a read whose span collapses on either ancestor is skipped**, counted as a translation error,
//!   and a position before a translation table's first row indexes the table at minus one;
//! * **a reverse read is scored in synthesis direction**: its matrix, key and flow order flipped
//!   and its zero edge flows dropped, against a haplotype that is the reverse complement of the
//!   reference in the read group's own flow order;
//! * **the CSV is written by `closeTool`**, which runs after a failure too.
//!
//! `--subsampling-ratio` draws from an unseeded `java.util.Random`, so a ratio strictly between
//! nought and one is the port's limitation: the reference itself answers differently on each run.
use super::*;
use gatk_tools::flow_based_read::FlowRead;
use gatk_tools::flow_pairhmm_align_reads_to_haplotypes::FlowHaplotype;
use gatk_tools::ground_truth_reads_builder as gtrb;
use htsjdk_bam::cigar::Op;
use htsjdk_bam::record::BamRecord;

/// One ancestor's reference, haplotype and score.
struct Scored {
    contig: String,
    start: i32,
    end: i32,
    haplotype: Vec<u8>,
    haplotype_length: i32,
    clipped: Option<Vec<u8>>,
    unclipped: Vec<u8>,
    softclip_front_fill: i64,
    score: f64,
}

fn thrown(error: gtrb::Thrown) -> Thrown {
    Thrown::non_user(error.class, error.message)
}

fn first_op(read: &BamRecord) -> Option<Op> {
    read.cigar.elements.first().map(|element| element.op)
}

fn last_op(read: &BamRecord) -> Option<Op> {
    read.cigar.elements.last().map(|element| element.op)
}

fn reverse(read: &BamRecord) -> bool {
    read.flags & 0x10 != 0
}

/// `isEndSoftclipped`: the end the read was synthesised towards.
fn end_softclipped(read: &BamRecord) -> bool {
    if reverse(read) {
        first_op(read) == Some(Op::S)
    } else {
        last_op(read) == Some(Op::S)
    }
}

/// `isStartSoftclipped`.
fn start_softclipped(read: &BamRecord) -> bool {
    if reverse(read) {
        last_op(read) == Some(Op::S)
    } else {
        first_op(read) == Some(Op::S)
    }
}

/// `isEndPolyTSoftclipped`, which for a forward read tests the length of the FIRST element over
/// the first bases, whatever the clip at the end is.
fn end_polyt_softclipped(read: &BamRecord) -> bool {
    if !end_softclipped(read) {
        return false;
    }
    let bases = &read.read_bases;
    if !reverse(read) {
        let length = read.cigar.elements.first().map_or(0, |e| e.length as usize);
        (0..length).all(|n| bases.get(n) == Some(&b'T'))
    } else {
        let length = read.cigar.elements.last().map_or(0, |e| e.length as usize);
        (0..length).all(|n| bases.len() > n && bases[bases.len() - n - 1] == b'A')
    }
}

/// `getSoftclippedBases`: the bases between the soft clips.
fn softclipped_bases(read: &BamRecord) -> Vec<u8> {
    let start = match read.cigar.elements.first() {
        Some(e) if e.op == Op::S => e.length as usize,
        _ => 0,
    };
    let end = match read.cigar.elements.last() {
        Some(e) if e.op == Op::S => e.length as usize,
        _ => 0,
    };
    let bases = &read.read_bases;
    bases[start.min(bases.len())..bases.len().saturating_sub(end).max(start.min(bases.len()))]
        .to_vec()
}

/// `getUnclippedEnd`: the end plus the soft and hard clips after it.
fn unclipped_end(read: &BamRecord) -> i32 {
    let trailing: i32 = read
        .cigar
        .elements
        .iter()
        .rev()
        .take_while(|e| matches!(e.op, Op::S | Op::H))
        .map(|e| e.length as i32)
        .sum();
    gatk_engine::read_utils::end(read) + trailing
}

/// `buildReferenceHaplotype`: the reference oriented as the read was synthesised, with the false
/// SNP compensation's skip taken off its front. Returns the bases and the length on reference.
fn reference_haplotype(
    reference: &[u8],
    length: i32,
    read: &BamRecord,
    false_snp_compensation: bool,
) -> Result<(Vec<u8>, i32), Thrown> {
    let mut bases = gtrb::oriented(reference, reverse(read));
    let read_bases = gtrb::oriented(&softclipped_bases(read), reverse(read));
    let mut length = length;
    if false_snp_compensation {
        let first_h = *bases.first().ok_or_else(|| {
            Thrown::non_user(
                "java.lang.ArrayIndexOutOfBoundsException",
                "Index 0 out of bounds for length 0".to_string(),
            )
        })?;
        let first_r = *read_bases.first().ok_or_else(|| {
            Thrown::non_user(
                "java.lang.ArrayIndexOutOfBoundsException",
                "Index 0 out of bounds for length 0".to_string(),
            )
        })?;
        if first_h != first_r {
            let skip = gtrb::detect_false_snp(&bases, &read_bases).map_err(thrown)?;
            if skip != 0 {
                bases = bases[skip..].to_vec();
                length -= skip as i32;
            }
        }
    }
    gtrb::check_haplotype(&bases).map_err(thrown)?;
    Ok((bases, length))
}

/// `scoreReadAgainstHaplotype` and `scoreReadAgainstReference`, which are the same computation.
fn score(
    read: &BamRecord,
    haplotype: &[u8],
    info: &gatk_tools::flow_based_read::ReadGroupInfo,
    flow: &gatk_tools::flow_based_read::FlowArguments,
) -> Result<f64, Thrown> {
    let flow_haplotype = FlowHaplotype::new(haplotype, &info.flow_order).ok_or_else(|| {
        Thrown::non_user(
            "org.broadinstitute.hellbender.exceptions.GATKException",
            format!(
                "baseArrayToKey periodGuard tripped, on {}, flowOrder: {} This probably indicates the presence of a base (value) in the sequence that is not included in the provided flow order",
                String::from_utf8_lossy(haplotype),
                info.flow_order
            ),
        )
    })?;
    let mut flow_read =
        FlowRead::new(read, &info.flow_order, info.max_class, flow).map_err(flow_refusal)?;
    if reverse(read) {
        flow_read.flip_to_synthesis();
        flow_read.apply_alignment().map_err(flow_refusal)?;
    }
    if !flow_read.valid {
        return Ok(-1.0);
    }
    gatk_tools::flow_feature_mapper::compute_likelihood_local(
        &flow_read,
        &flow_haplotype,
        flow_haplotype.key.len(),
    )
    .map_err(flow_refusal)
}

/// A span on an ancestral contig: its name, its start and its end.
type Span = (String, i32, i32);

/// The translation tables, read the first time an ancestor and a contig need one.
struct Translators {
    base: String,
    tables: Vec<(String, gtrb::Translator)>,
}

impl Translators {
    fn position(&mut self, ancestor: &str, contig: &str, from: i32) -> Result<i32, Thrown> {
        let key = format!("{ancestor}.{contig}.csv");
        if !self.tables.iter().any(|(k, _)| *k == key) {
            let path = format!("{}{key}", self.base);
            let text = std::fs::read_to_string(&path).map_err(|_| {
                Thrown::non_user(
                    PORT_LIMITATION,
                    format!("{path}: an unreadable translation table is not ported. This message is the port's own and not GATK's."),
                )
            })?;
            self.tables
                .push((key.clone(), gtrb::Translator::parse(&text)));
        }
        let table = &self
            .tables
            .iter()
            .find(|(k, _)| *k == key)
            .expect("a table")
            .1;
        table.translate(from).ok_or_else(|| {
            Thrown::non_user(
                "java.lang.ArrayIndexOutOfBoundsException",
                format!(
                    "Index -1 out of bounds for length {}",
                    table.positions.len()
                ),
            )
        })
    }

    /// Both ancestors, maternal first; `None` for a span that collapses on either.
    fn translate(
        &mut self,
        contig: &str,
        start: i32,
        end: i32,
    ) -> Result<Option<[Span; 2]>, Thrown> {
        let mut spans = Vec::new();
        for ancestor in [gtrb::MATERNAL, gtrb::PATERNAL] {
            let s = self.position(ancestor, contig, start)?;
            let e = self.position(ancestor, contig, end)?;
            if e <= s {
                return Ok(None);
            }
            spans.push((format!("{contig}_{ancestor}"), s, e));
        }
        let paternal = spans.pop().expect("two spans");
        let maternal = spans.pop().expect("two spans");
        Ok(Some([maternal, paternal]))
    }
}

pub fn ground_truth_reads_builder(parser: &Parser) -> Outcome {
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "GroundTruthReadsBuilder")?;
    let filter = read_filter(parser, &filters, &header)?;
    let required = |name: &str| {
        argument(parser, name).ok_or_else(|| {
            Thrown::command_line(format!(
                "Argument {name} was missing: Argument '{name}' is required"
            ))
        })
    };
    let output = required("output-csv")?;
    let open = |path: &str| {
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(path))
            .map_err(|error| Thrown::user(format!("{error:?}")))
    };
    let mut maternal_reference = open(&required("maternal-ref")?)?;
    let mut paternal_reference = open(&required("paternal-ref")?)?;
    let mut translators = Translators {
        base: required("ancestral-translators-base-path")?,
        tables: Vec::new(),
    };
    let engine =
        scalar(parser, "likelihood-calculation-engine").unwrap_or_else(|| "PairHMM".to_string());
    if engine != "FlowBased" {
        return Err(Thrown::non_user(
            "org.broadinstitute.hellbender.exceptions.GATKException",
            "must use a flow based likelihood calculation engine".to_string(),
        ));
    }
    let mut csv = gtrb::header();
    csv.push('\n');
    let result = ground_truth_reads_builder_traverse(
        parser,
        &source,
        &header,
        &intervals,
        &filter,
        [&mut maternal_reference, &mut paternal_reference],
        &mut translators,
        &mut csv,
    );
    let bytes = if output.ends_with(".gz") {
        java_gzip(csv.as_bytes(), 6)
    } else {
        csv.into_bytes()
    };
    write_file(&output, &bytes)?;
    result?;
    Ok(None)
}

#[allow(clippy::too_many_arguments)]
fn ground_truth_reads_builder_traverse(
    parser: &Parser,
    source: &ReadsDataSource,
    header: &SamHeader,
    intervals: &[gatk_engine::interval::SimpleInterval],
    filter: &Filter,
    ancestors: [&mut gatk_engine::reference::ReferenceFileSource; 2],
    translators: &mut Translators,
    csv: &mut String,
) -> Result<(), Thrown> {
    use gatk_engine::read_utils;
    let [maternal_reference, paternal_reference] = ancestors;
    let flow = flow_arguments(parser);
    let ratio = double_or(parser, "subsampling-ratio", 1.0);
    let max_output_reads = number_or(parser, "max-output-reads", 20_000_000);
    let output_flow_length = number_or(parser, "output-flow-length", 0);
    let prepend = scalar(parser, "prepend-sequence").filter(|v| v != "null");
    let append = scalar(parser, "append-sequence").filter(|v| v != "null");
    let min_mq = double_or(parser, "min-mq", 0.0);
    let max_rq = double_or(parser, "max-rq", 0.0);
    let include_supp = flag(parser, "include-supp-align");
    let min_score = double_or(parser, "min-haplotype-score", 0.0);
    let min_delta = double_or(parser, "min-haplotype-score-delta", 0.0);
    let padding = number_or(parser, "haplotype-output-padding-size", 8);
    let discard_softclipped = flag(parser, "discard-non-polyt-softclipped-reads");
    let fill_q = flag(parser, "fill-trimmed-reads-Q");
    let fill_z = flag(parser, "fill-trimmed-reads-Z");
    let fill_trimmed = flag(parser, "fill-trimmed-reads");
    let fill_softclipped = flag(parser, "fill-softclipped-reads");
    let false_snp = flag(parser, "false-snp-compensation");
    let no_output = flag(parser, "gt-no-output");
    if ratio > 0.0 && ratio < 1.0 {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "--subsampling-ratio between nought and one draws from an unseeded java.util.Random. This message is the port's own and not GATK's."
                .to_string(),
        ));
    }
    let mut reference = match argument(parser, "reference") {
        Some(path) => Some(
            gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&path))
                .map_err(|error| Thrown::user(format!("{error:?}")))?,
        ),
        None => None,
    };
    let reads = gatk_tools::read_walker::traverse(source, intervals, filter)
        .map_err(reads_traversal_error)?;
    let query = |source: &mut gatk_engine::reference::ReferenceFileSource,
                 contig: &str,
                 start: i32,
                 end: i32| {
        source
            .query(contig, start, end)
            .map_err(|error| Thrown::user(format!("{error:?}")))
    };
    let mut count = 0i32;
    for read in &reads {
        if max_output_reads != 0 && count >= max_output_reads {
            break;
        }
        if min_mq != 0.0 && f64::from(read.mapping_quality) < min_mq {
            continue;
        }
        if read.flags & 0x800 != 0 && !include_supp {
            continue;
        }
        if discard_softclipped && end_softclipped(read) && !end_polyt_softclipped(read) {
            continue;
        }
        // `random.nextDouble() > subsamplingRatio`, which a ratio of one never passes and a
        // negative one always does.
        if ratio < 0.0 {
            continue;
        }
        let info =
            gatk_tools::flow_based_read::read_group_info(read, header).map_err(flow_refusal)?;
        let build =
            || FlowRead::new(read, &info.flow_order, info.max_class, &flow).map_err(flow_refusal);
        let mut flow_read: Option<FlowRead> = None;
        if max_rq != 0.0 {
            let built = build()?;
            let quality: f64 = (0..built.key.len())
                .map(|n| built.prob(n, built.max_hmer))
                .sum();
            if quality > max_rq {
                continue;
            }
            flow_read = Some(built);
        }
        let flow_read = match flow_read {
            Some(built) => built,
            None => build()?,
        };

        let contig = usize::try_from(read.reference_index)
            .ok()
            .and_then(|at| header.sequences.get(at))
            .map(|sequence| sequence.name.clone())
            .unwrap_or_default();
        let (start, end) = (read_utils::start(read), read_utils::end(read));
        let Some(spans) = translators.translate(&contig, start, end)? else {
            continue;
        };
        let rev = reverse(read);
        let mut scored: Vec<Scored> = Vec::new();
        for (index, (ancestor_contig, s, e)) in spans.iter().enumerate() {
            let source: &mut gatk_engine::reference::ReferenceFileSource = if index == 0 {
                maternal_reference
            } else {
                paternal_reference
            };
            let bases = query(source, ancestor_contig, *s, *e)?;
            let (haplotype, haplotype_length) =
                reference_haplotype(&bases, e - s + 1, read, false_snp)?;
            scored.push(Scored {
                contig: ancestor_contig.clone(),
                start: *s,
                end: *e,
                haplotype,
                haplotype_length,
                clipped: None,
                unclipped: Vec::new(),
                softclip_front_fill: 0,
                score: 0.0,
            });
        }
        // `buildExtendedRef`, maternal then paternal.
        for (index, sh) in scored.iter_mut().enumerate() {
            let source: &mut gatk_engine::reference::ReferenceFileSource = if index == 0 {
                maternal_reference
            } else {
                paternal_reference
            };
            let (mut extend_start, mut extend_end) = (0i32, 0i32);
            if fill_softclipped {
                let element = if !rev {
                    read.cigar.elements.last()
                } else {
                    read.cigar.elements.first()
                };
                if let Some(element) = element.filter(|e| e.op == Op::S) {
                    if !rev {
                        extend_end += element.length as i32;
                    } else {
                        extend_start += element.length as i32;
                    }
                }
            }
            if !rev {
                extend_end += padding;
            } else {
                extend_start += padding;
            }
            if output_flow_length != 0
                && should_fill_from_haplotype(read, fill_softclipped, fill_trimmed, fill_q, fill_z)
            {
                let length = (sh.end + extend_end) - (sh.start - extend_start);
                let delta = (output_flow_length - length).max(0) + 50;
                if !rev {
                    extend_end += delta;
                } else {
                    extend_start += delta;
                }
            }
            let delta = (sh.end - sh.start + 1) - sh.haplotype_length;
            if delta != 0 {
                if !rev {
                    extend_start -= delta;
                } else {
                    extend_end -= delta;
                }
            }
            let max_end = source.sequence_length(&sh.contig).unwrap_or(0) as i32;
            if start_softclipped(read) {
                let from = 1.max(sh.start - extend_start);
                let to = max_end.min(sh.end + extend_end);
                sh.clipped = Some(query(source, &sh.contig, from, to)?);
            }
            let front = if !rev {
                read.cigar.elements.first()
            } else {
                read.cigar.elements.last()
            };
            if let Some(element) = front.filter(|e| e.op == Op::S) {
                if !rev {
                    extend_start += element.length as i32;
                } else {
                    extend_end += element.length as i32;
                }
            }
            let from = 1.max(sh.start - extend_start);
            let to = max_end.min(sh.end + extend_end);
            sh.unclipped = query(source, &sh.contig, from, to)?;
        }

        // The reference's own score, then each ancestor's unless its haplotype IS the reference.
        let reference_bases = match reference.as_mut() {
            Some(source) => query(source, &contig, start, end)?,
            None => Vec::new(),
        };
        let (reference_haplotype_bases, _) =
            reference_haplotype(&reference_bases, end - start + 1, read, false_snp)?;
        let reference_score = score(read, &reference_haplotype_bases, &info, &flow)?;
        let oriented_reference = gtrb::oriented(&reference_bases, rev);
        for sh in scored.iter_mut() {
            sh.score = if sh.haplotype == oriented_reference {
                reference_score
            } else {
                score(read, &sh.haplotype, &info, &flow)?
            };
        }
        let (maternal_score, paternal_score) = (scored[0].score, scored[1].score);
        if min_score != 0.0 && maternal_score.min(paternal_score) > min_score {
            continue;
        }
        if min_delta != 0.0 && (maternal_score - paternal_score).abs() > min_delta {
            continue;
        }
        count += 1;

        // `emit`.
        let tm =
            match read.tags.get(htsjdk_bam::tag::Tag::new(b"tm")) {
                Some(htsjdk_bam::tag::TagValue::Str(text)) => Some(text.clone()),
                Some(htsjdk_bam::tag::TagValue::Char(c)) => Some((*c as char).to_string()),
                Some(_) => return Err(Thrown::non_user(
                    PORT_LIMITATION,
                    "a tm tag that is not a string. This message is the port's own and not GATK's."
                        .to_string(),
                )),
                None => None,
            };
        let has = |c: char| tm.as_deref().is_some_and(|tm| tm.contains(c));
        let fill_value = if end_softclipped(read) {
            gtrb::SOFTCLIP_FILL_VALUE
        } else if has('Q') || has('Z') {
            if has('A') {
                gtrb::UNKNOWN_FILL_VALUE
            } else {
                gtrb::NONREF_FILL_VALUE
            }
        } else {
            gtrb::DEFAULT_FILL_VALUE
        };
        // Paternal first, as the reference builds them.
        let mut keys: [Vec<i32>; 2] = [Vec::new(), Vec::new()];
        for index in [1usize, 0] {
            let sh = &mut scored[index];
            let key = gtrb::haplotype_key(&sh.unclipped, &info.flow_order, rev).map_err(thrown)?;
            sh.softclip_front_fill = if start_softclipped(read) {
                let clipped = gtrb::haplotype_key(
                    sh.clipped.as_deref().unwrap_or_default(),
                    &info.flow_order,
                    rev,
                )
                .map_err(thrown)?;
                key.len() as i64 - clipped.len() as i64
            } else {
                0
            };
            if output_flow_length < 0 {
                return Err(Thrown::non_user(
                    "java.lang.NegativeArraySizeException",
                    output_flow_length.to_string(),
                ));
            }
            let length = if output_flow_length != 0 {
                output_flow_length as usize
            } else {
                key.len()
            };
            let mut out = vec![fill_value; length];
            let copied = length.min(key.len());
            out[..copied].copy_from_slice(&key[..copied]);
            keys[index] = out;
        }
        let mut sequences: [String; 2] = [String::new(), String::new()];
        for index in [1usize, 0] {
            let sh = &scored[index];
            let oriented = gtrb::oriented(&sh.unclipped, rev);
            let count = gtrb::key_bases(&keys[index]);
            if count > oriented.len() {
                return Err(Thrown::non_user(
                    "java.lang.StringIndexOutOfBoundsException",
                    format!("begin 0, end {count}, length {}", oriented.len()),
                ));
            }
            let mut text = prepend.clone().unwrap_or_default();
            text.push_str(&String::from_utf8_lossy(&oriented[..count]));
            if let Some(append) = &append {
                text.push_str(append);
            }
            sequences[index] = text;
        }
        if !fill_softclipped {
            for index in [1usize, 0] {
                let limit = scored[index]
                    .softclip_front_fill
                    .min(keys[index].len() as i64);
                for value in keys[index].iter_mut().take(limit.max(0) as usize) {
                    *value = gtrb::SOFTCLIP_FILL_VALUE;
                }
            }
        }
        let same = sequences[0] == sequences[1];
        let best = if paternal_score > maternal_score {
            1
        } else {
            0
        };
        let consensus = gtrb::consensus_key(&keys[1], &keys[0]);
        let read_sequence = gtrb::oriented(&read.read_bases, rev);
        let mut read_key = flow_read.key.clone();
        if rev {
            read_key.reverse();
        }
        let read_flow_order = gtrb::oriented(info.flow_order.as_bytes(), rev);
        let interval = |sh: &Scored| format!("{}:{}-{}", sh.contig, sh.start, sh.end);
        let columns = [
            read.read_name.clone(),
            contig.clone(),
            start.to_string(),
            end.to_string(),
            gatk_engine::java_format::format_decimals(paternal_score, 6),
            gatk_engine::java_format::format_decimals(maternal_score, 6),
            gatk_engine::java_format::format_decimals(reference_score, 6),
            gtrb::read_key_csv(&read_key, &read_sequence, &read_flow_order).map_err(thrown)?,
            if same {
                gtrb::key_csv(&consensus)
            } else {
                gtrb::key_csv(&keys[best])
            },
            gtrb::key_csv(&consensus),
            tm.clone().unwrap_or_default(),
            read.mapping_quality.to_string(),
            read.flags.to_string(),
            read.cigar.to_text(),
            String::from_utf8_lossy(&read_sequence).into_owned(),
            sequences[1].clone(),
            sequences[0].clone(),
            sequences[best].clone(),
            read_utils::unclipped_start(read).to_string(),
            unclipped_end(read).to_string(),
            interval(&scored[1]),
            interval(&scored[0]),
        ];
        if !no_output {
            csv.push_str(&columns.join(","));
            csv.push('\n');
        }
    }
    Ok(())
}

/// `shouldFillFromHaplotype`: a soft clip decides first, then the `tm` tag's letters.
fn should_fill_from_haplotype(
    read: &BamRecord,
    fill_softclipped: bool,
    fill_trimmed: bool,
    fill_q: bool,
    fill_z: bool,
) -> bool {
    if end_softclipped(read) {
        return fill_softclipped;
    }
    let tm = match read.tags.get(htsjdk_bam::tag::Tag::new(b"tm")) {
        Some(htsjdk_bam::tag::TagValue::Str(text)) => text.clone(),
        Some(htsjdk_bam::tag::TagValue::Char(c)) => (*c as char).to_string(),
        _ => return true,
    };
    if tm.contains('A') {
        false
    } else if tm.contains('Z') && (fill_trimmed || fill_z) {
        true
    } else {
        tm.contains('Q') && (fill_trimmed || fill_q)
    }
}

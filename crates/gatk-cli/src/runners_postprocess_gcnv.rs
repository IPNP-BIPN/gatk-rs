//! `PostprocessGermlineCNVCalls` past its import check: `onTraversalStart` and the three outputs
//! of `onTraversalSuccess`, in the reference's order.
//!
//! * **the intervals VCF** (`GermlineCNVIntervalVariantComposer`) and **the denoised table**
//!   (`concatenateDenoisedCopyRatioFiles`) are Java over the call shards;
//! * **the segments VCF** is `segment_gcnv_calls.py`, ported in [`gatk_tools::gcnv_segments`] as
//!   the frozen reference computes it, then `GermlineCNVSegmentVariantComposer` over the qualities
//!   as the script prints them, two decimals, and Java reads them back.
//!
//! Not ported, each refused in the port's words: `--clustered-breakpoints` and
//! `--input-intervals-vcf` (the script then reads the joint call set with PyVCF instead of running
//! Viterbi), and a `--reference`, whose only use is the REF base of each segment.

use super::*;
use gatk_tools::gcnv_segments;
use gatk_tools::model_segments_models::Span;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use htsjdk_vcf::variant::{Genotype, Value, VariantContext};

const TOOL: &str = "PostprocessGermlineCNVCalls";

fn illegal_state(message: impl Into<String>) -> Thrown {
    Thrown::non_user("java.lang.IllegalStateException", message.into())
}

/// `CopyNumberArgumentValidationUtils.validateInputs(File...)`: a file must be readable, a
/// directory is taken as it is.
fn validate_input(path: &str) -> Result<(), Thrown> {
    if std::path::Path::new(path).is_file() {
        tables::can_read_file(path)?;
    }
    Ok(())
}

/// `CopyNumberArgumentValidationUtils.getCanonicalPath`, which is what the script is handed.
fn canonical(path: &str) -> String {
    std::fs::canonicalize(path)
        .map(|p| p.to_string_lossy().into_owned())
        .unwrap_or_else(|_| java_absolute_path(path))
}

/// `FastMath.round(double)`: half up, exactly.
fn fast_math_round(value: f64) -> i64 {
    let floor = value.floor();
    (if value - floor >= 0.5 {
        floor + 1.0
    } else {
        floor
    }) as i64
}

fn compound(
    key: &str,
    id: &str,
    number: Cardinality,
    line_type: LineType,
    text: &str,
) -> HeaderLine {
    HeaderLine::Compound {
        key: key.to_string(),
        id: id.to_string(),
        number,
        line_type,
        description: text.to_string(),
        extra: Vec::new(),
    }
}

/// `VCFHeader.setSequenceDictionary`: one contig line per sequence.
fn contig_lines(dictionary: &[SequenceRecord]) -> Vec<HeaderLine> {
    dictionary
        .iter()
        .enumerate()
        .map(|(index, sequence)| {
            let mut fields = vec![
                ("ID".to_string(), sequence.name.clone()),
                ("length".to_string(), sequence.length.to_string()),
            ];
            if let Some(assembly) = sequence.attributes.get("AS") {
                fields.push(("assembly".to_string(), assembly.to_string()));
            }
            HeaderLine::Contig {
                index: index as i32,
                fields,
            }
        })
        .collect()
}

fn ref_allele() -> Allele {
    Allele::from_str("N", true).expect("a base")
}

/// What `onTraversalStart` settles.
struct Start {
    dictionary: Vec<SequenceRecord>,
    calls: Vec<String>,
    models: Vec<String>,
    intervals: Vec<Vec<Span>>,
    sample_index: usize,
    sample: String,
    allosomal: Vec<String>,
    autosomal_ref: i64,
}

/// `IntervalUtils.compareLocatables` on each shard's first interval: the shard order.
fn shard_order(spans: &[Vec<Span>], dictionary: &[SequenceRecord]) -> Vec<usize> {
    let index_of = |contig: &str| dictionary.iter().position(|s| s.name == contig);
    let mut order: Vec<usize> = (0..spans.len()).collect();
    order.sort_by_key(|&i| {
        let first = &spans[i][0];
        (index_of(&first.contig), first.start, first.end)
    });
    order
}

fn on_traversal_start(parser: &Parser, dictionary: &[SequenceRecord]) -> Result<Start, Thrown> {
    let calls = arguments(parser, "calls-shard-path");
    let models = arguments(parser, "model-shard-path");
    if calls.len() != models.len() {
        return Err(illegal(
            "The number of input call shards must match the number of input model shards.",
        ));
    }
    for path in calls.iter().chain(&models) {
        validate_input(path)?;
    }
    validate_input(&require(parser, "contig-ploidy-calls"))?;

    let read = |shards: &[String]| -> Result<Vec<Vec<Span>>, Thrown> {
        shards
            .iter()
            .map(|shard| Ok(tables::read_interval_spans(&format!("{shard}/interval_list.tsv"))?.1))
            .collect()
    };
    let call_spans = read(&calls)?;
    let model_spans = read(&models)?;
    let call_order = shard_order(&call_spans, dictionary);
    let model_order = shard_order(&model_spans, dictionary);
    let sorted_calls: Vec<Vec<Span>> = call_order.iter().map(|&i| call_spans[i].clone()).collect();
    let sorted_models: Vec<Vec<Span>> = model_order
        .iter()
        .map(|&i| model_spans[i].clone())
        .collect();
    if sorted_calls != sorted_models {
        return Err(illegal(
            "The interval lists found in model and call shards do not match. Make sure that the \
             calls and model paths are provided in matching order.",
        ));
    }
    // `checkForSingletonInterval`. shortcut: with several singleton contigs the reference names
    // the first of a `HashMap`, here the first in genomic order; upgrade if a test meets two.
    let all: Vec<&Span> = sorted_calls.iter().flatten().collect();
    for span in &all {
        if all
            .iter()
            .filter(|other| other.contig == span.contig)
            .count()
            == 1
        {
            return Err(illegal(format!(
                "Records contain a singleton interval on contig ({}). Please run FilterIntervals \
                 tool first.",
                span.contig
            )));
        }
    }

    let allosomal = arguments(parser, "allosomal-contig");
    let names: Vec<String> = dictionary.iter().map(|s| s.name.clone()).collect();
    if !allosomal.iter().all(|contig| names.contains(contig)) {
        let order = |list: &[String]| {
            gatk_engine::java_hash::hash_set_order(list)
                .unwrap_or_else(|_| list.to_vec())
                .join(", ")
        };
        return Err(illegal(format!(
            "The specified allosomal contigs must be contained in the SAM sequence dictionary of \
             the call-set (specified allosomal contigs: [{}], all contigs: [{}])",
            order(&allosomal),
            order(&names)
        )));
    }

    let sample_index: usize = scalar(parser, "sample-index")
        .and_then(|v| v.parse().ok())
        .unwrap_or(0);
    let sorted_call_paths: Vec<String> = call_order.iter().map(|&i| calls[i].clone()).collect();
    let sample_name = |shard: &str| -> Result<String, Thrown> {
        let path = format!("{shard}/SAMPLE_{sample_index}/sample_name.txt");
        let text = std::fs::read_to_string(&path).map_err(|_| {
            bad_input(format!(
                "Could not read the sample name text file at {}.",
                java_absolute_path(&path)
            ))
        })?;
        Ok(text
            .split('\n')
            .next()
            .unwrap_or("")
            .trim_end_matches('\r')
            .to_string())
    };
    let sample = sample_name(&sorted_call_paths[0])?;
    for shard in &sorted_call_paths[1..] {
        if sample_name(shard)? != sample {
            return Err(illegal_state(
                "The sample name is not the same for all of the shards.",
            ));
        }
    }
    Ok(Start {
        dictionary: dictionary.to_vec(),
        calls: sorted_call_paths,
        models: model_order.iter().map(|&i| models[i].clone()).collect(),
        intervals: sorted_calls,
        sample_index,
        sample,
        allosomal,
        autosomal_ref: i64::from(number_or(parser, "autosomal-ref-copy-number", 2)),
    })
}

/// A sample table of one call shard, its sample name checked against the run's.
fn shard_rows(
    start: &Start,
    shard: usize,
    file: &str,
    what: &str,
) -> Result<(Vec<String>, Vec<Vec<String>>), Thrown> {
    let path = format!(
        "{}/SAMPLE_{}/{file}",
        start.calls[shard], start.sample_index
    );
    let (sample, columns, rows) = tables::read_sample_rows(&path)?;
    if sample != start.sample {
        return Err(illegal_state(format!(
            "Sample name found in the header of {what} for shard {shard} {}different from the \
             expected sample name (found: {sample}, expected: {}).",
            if what.starts_with("denoised") {
                "is "
            } else {
                ""
            },
            start.sample
        )));
    }
    Ok((columns, rows))
}

fn java_double(text: &str) -> Result<f64, Thrown> {
    gatk_tools::variant_recalibrator_model::java_parse_double(text)
        .ok_or_else(|| bad_input(format!("Expected a double value but found {text}")))
}

/// `generateIntervalsVCFFileFromAllShards`.
fn intervals_vcf(parser: &Parser, start: &Start) -> Result<String, Thrown> {
    let mut header = VcfHeader::new();
    header.samples = vec![start.sample.clone()];
    header.lines.extend(contig_lines(&start.dictionary));
    header
        .lines
        .extend(default_tool_vcf_header_lines(parser, TOOL));
    header.lines.extend([
        compound("FORMAT", "GT", Cardinality::Fixed(1), LineType::String, "Genotype"),
        compound("FORMAT", "CN", Cardinality::Fixed(1), LineType::Integer, "Copy number maximum a posteriori value"),
        compound("FORMAT", "CNLP", Cardinality::Unbounded, LineType::Integer, "Copy number log posterior (in Phred-scale) rounded down"),
        compound("FORMAT", "CNQ", Cardinality::Fixed(1), LineType::Integer, "Genotype call quality as the difference between the best and second best phred-scaled log posterior scores"),
        compound("INFO", "END", Cardinality::Fixed(1), LineType::Integer, "End coordinate of the variant"),
    ]);
    let all_alleles = vec![
        ref_allele(),
        super::joint_germline_cnv::del_allele(),
        super::joint_germline_cnv::dup_allele(),
    ];
    let mut records = Vec::new();
    for (shard, spans) in start.intervals.iter().enumerate() {
        let (columns, posteriors) =
            shard_rows(start, shard, "log_q_c_tc.tsv", "copy-number posterior file")?;
        let (_, baselines) = shard_rows(
            start,
            shard,
            "baseline_copy_number_t.tsv",
            "baseline copy-number file",
        )?;
        if spans.len() != posteriors.len() {
            return Err(illegal_state(format!(
                "The number of entries in the copy-number posterior file for shard {shard} does not \
                 match the number of entries in the shard interval list (posterior list size: {}, \
                 interval list size: {})",
                posteriors.len(),
                spans.len()
            )));
        }
        if spans.len() != baselines.len() {
            return Err(illegal_state(format!(
                "The number of entries in the baseline copy-number file for shard {shard} does not \
                 match the number of entries in the shard interval list (baseline copy-number list \
                 size: {}, interval list size: {})",
                baselines.len(),
                spans.len()
            )));
        }
        let states: Vec<i64> = columns
            .iter()
            .map(|c| {
                c.strip_prefix("COPY_NUMBER_")
                    .and_then(|n| n.parse().ok())
                    .unwrap_or(0)
            })
            .collect();
        for ((span, row), baseline) in spans.iter().zip(&posteriors).zip(&baselines) {
            let values: Vec<f64> = row
                .iter()
                .map(|v| java_double(v))
                .collect::<Result<_, _>>()?;
            // `calculateCopyNumberPLVector`: `DoubleStream.max`, NaN if any is NaN.
            let largest = values.iter().copied().fold(f64::NEG_INFINITY, |a, b| {
                if a.is_nan() || b.is_nan() {
                    f64::NAN
                } else {
                    a.max(b)
                }
            });
            let pl: Vec<i64> = values
                .iter()
                .map(|v| (-10.0 * (v - largest) * std::f64::consts::LOG10_E).floor() as i32 as i64)
                .collect();
            // `Stream.max(comparingDouble)`: the first of equal maxima.
            let mut best = 0;
            for (i, v) in values.iter().enumerate() {
                if v.total_cmp(&values[best]) == std::cmp::Ordering::Greater {
                    best = i;
                }
            }
            let map = states[best];
            let mut sorted = pl.clone();
            sorted.sort_unstable();
            if sorted[0] != 0 {
                return Err(illegal_state(
                    "Something went wrong. Smallest copy number posterior score must be 0",
                ));
            }
            let reference = if start.allosomal.contains(&span.contig) {
                baseline[0].parse::<i64>().map_err(|_| {
                    bad_input(format!("Expected an int value but found {}", baseline[0]))
                })?
            } else {
                start.autosomal_ref
            };
            let allele = if map > reference {
                super::joint_germline_cnv::dup_allele()
            } else if map < reference {
                super::joint_germline_cnv::del_allele()
            } else {
                ref_allele()
            };
            let mut genotype = Genotype::new(&start.sample, vec![allele]);
            genotype.extended = vec![
                ("CN".to_string(), Value::Int(map)),
                (
                    "CNLP".to_string(),
                    Value::List(pl.iter().map(|&p| Value::Int(p)).collect()),
                ),
                ("CNQ".to_string(), Value::Int(sorted[1])),
            ];
            let mut record =
                VariantContext::new(&span.contig, i64::from(span.start), all_alleles.clone());
            record.stop = i64::from(span.end);
            record.id = format!("CNV_{}_{}_{}", span.contig, span.start, span.end);
            record.attributes = vec![("END".to_string(), Value::Int(i64::from(span.end)))];
            record.genotypes = vec![genotype].into();
            records.push(record);
        }
    }
    write_vcf_honouring_lenient(parser, &header, &records)
}

/// `generateSegmentsVCFFileFromAllShards`, segmentation included.
fn segments_vcf(parser: &Parser, start: &Start) -> Result<String, Thrown> {
    let models: Vec<String> = start.models.iter().map(|p| canonical(p)).collect();
    let calls: Vec<String> = start.calls.iter().map(|p| canonical(p)).collect();
    let segments = gcnv_segments::read_sample_inputs(
        &models,
        &calls,
        &canonical(&require(parser, "contig-ploidy-calls")),
        start.sample_index,
    )
    .and_then(|inputs| gcnv_segments::segments(&inputs))
    // Whatever stops the script surfaces in Java as its exit status.
    .map_err(|_| Thrown::user("Python return code was non-zero."))?;

    let mut header = VcfHeader::new();
    header.samples = vec![start.sample.clone()];
    header
        .lines
        .extend(default_tool_vcf_header_lines(parser, TOOL));
    header.lines.push(HeaderLine::Unstructured {
        key: "gcnvVcfSchemaVersion".to_string(),
        value: "2.0".to_string(),
    });
    header.lines.extend(contig_lines(&start.dictionary));
    header.lines.extend([
        compound("INFO", "SVTYPE", Cardinality::Fixed(1), LineType::String, "Type of structural variant"),
        compound("INFO", "SVLEN", Cardinality::Unbounded, LineType::Integer, "Difference in length between REF and ALT alleles"),
        compound("INFO", "AC_Orig", Cardinality::A, LineType::Integer, "Original AC"),
        compound("INFO", "AF_Orig", Cardinality::A, LineType::Float, "Original AF"),
        compound("INFO", "AN_Orig", Cardinality::Fixed(1), LineType::Integer, "Original AN"),
        compound("FORMAT", "GT", Cardinality::Fixed(1), LineType::String, "Segment genotype"),
        compound("FORMAT", "CN", Cardinality::Fixed(1), LineType::Integer, "Segment most-likely copy-number call"),
        compound("FORMAT", "NP", Cardinality::Fixed(1), LineType::Integer, "Number of points (i.e. targets or bins) in the segment"),
        compound("FORMAT", "QS", Cardinality::Fixed(1), LineType::Integer, "Complementary Phred-scaled probability that at least one point (i.e. target or bin) in the segment agrees with the segment copy-number call"),
        compound("FORMAT", "QA", Cardinality::Fixed(1), LineType::Integer, "Complementary Phred-scaled probability that all points (i.e. targets or bins) in the segment agree with the segment copy-number call"),
        compound("FORMAT", "QSS", Cardinality::Fixed(1), LineType::Integer, "Complementary Phred-scaled probability that the segment start position is a genuine copy-number changepoint"),
        compound("FORMAT", "QSE", Cardinality::Fixed(1), LineType::Integer, "Complementary Phred-scaled probability that the segment end position is a genuine copy-number changepoint"),
        compound("INFO", "END", Cardinality::Fixed(1), LineType::Integer, "End coordinate of the variant"),
        HeaderLine::Filter {
            id: "LOW_QS".to_string(),
            description: "Depth-only copy number record whose QS value is lower than user specified threshold".to_string(),
        },
        HeaderLine::Filter {
            id: "FREQ".to_string(),
            description: "Depth-only copy number record whose AF value is higher than user specified threshold".to_string(),
        },
    ]);

    let mut records = Vec::new();
    for segment in segments {
        // The script prints each quality with two decimals and Java parses that text.
        let [qs, qa, qss, qse] = segment.qualities.map(|q| {
            gcnv_segments::python_fixed2(q)
                .parse::<f64>()
                .unwrap_or(f64::NAN)
        });
        let call = segment.call as i64;
        let reference = if start.allosomal.contains(&segment.contig) {
            segment.baseline as i64
        } else {
            start.autosomal_ref
        };
        let refa = ref_allele();
        let mut alleles = vec![refa.clone()];
        if call > reference {
            alleles.push(super::joint_germline_cnv::dup_allele());
        } else if call < reference {
            alleles.push(super::joint_germline_cnv::del_allele());
        }
        let genotype_alleles = super::joint_germline_cnv::alleles_from_copy_number(
            call as i32,
            reference as i32,
            &refa,
        );
        let mut genotype = Genotype::new(&start.sample, genotype_alleles);
        genotype.extended = vec![
            ("CN".to_string(), Value::Int(call)),
            ("NP".to_string(), Value::Int(segment.num_points as i64)),
            ("QS".to_string(), Value::Int(fast_math_round(qs))),
            ("QA".to_string(), Value::Int(fast_math_round(qa))),
            ("QSS".to_string(), Value::Int(fast_math_round(qss))),
            ("QSE".to_string(), Value::Int(fast_math_round(qse))),
        ];
        let mut record = VariantContext::new(&segment.contig, segment.start, alleles);
        record.stop = segment.end;
        record.id = format!("CNV_{}_{}_{}", segment.contig, segment.start, segment.end);
        record.attributes = vec![("END".to_string(), Value::Int(segment.end))];
        record.genotypes = vec![genotype].into();
        record.log10_p_error = qs / -10.0;
        records.push(record);
    }
    write_vcf_honouring_lenient(parser, &header, &records)
}

/// `concatenateDenoisedCopyRatioFiles`.
fn denoised_table(start: &Start) -> Result<String, Thrown> {
    let mut text = tables::header(&start.sample, Some(&start.dictionary));
    text.push_str("CONTIG\tSTART\tEND\tLINEAR_COPY_RATIO\n");
    for (shard, spans) in start.intervals.iter().enumerate() {
        let (_, rows) = shard_rows(
            start,
            shard,
            "mu_denoised_copy_ratio_t.tsv",
            "denoised copy ratio file",
        )?;
        if spans.len() != rows.len() {
            return Err(illegal_state(format!(
                "The number of entries in denoised copy ratio file for shard {shard} does not match \
                 the number of entries in the shard interval list (copy ratio list size: {}, \
                 interval list size: {})",
                rows.len(),
                spans.len()
            )));
        }
        for (span, row) in spans.iter().zip(&rows) {
            let value = java_double(&row[0])?;
            // `DataLine.append(double)`: a whole value as a long.
            let rounded = fast_math_round(value);
            let printed = if rounded as f64 == value {
                rounded.to_string()
            } else {
                jmath::to_string::double_to_string(value)
            };
            text.push_str(&format!(
                "{}\t{}\t{}\t{printed}\n",
                span.contig, span.start, span.end
            ));
        }
    }
    Ok(text)
}

pub(super) fn run(parser: &Parser, dictionary: &[SequenceRecord]) -> Outcome {
    if argument(parser, "clustered-breakpoints").is_some()
        || argument(parser, "input-intervals-vcf").is_some()
    {
        return Err(port_limit(
            "Segmenting along --clustered-breakpoints and --input-intervals-vcf",
        ));
    }
    if argument(parser, "reference").is_some() {
        return Err(port_limit(
            "Reading each segment's REF base from --reference",
        ));
    }
    let start = on_traversal_start(parser, dictionary)?;
    let intervals = intervals_vcf(parser, &start)?;
    write_variant_output(
        parser,
        &require(parser, "output-genotyped-intervals"),
        &intervals,
    )?;
    let segments = segments_vcf(parser, &start)?;
    write_variant_output(
        parser,
        &require(parser, "output-genotyped-segments"),
        &segments,
    )?;
    let denoised = denoised_table(&start)?;
    write_file(
        &require(parser, "output-denoised-copy-ratios"),
        denoised.as_bytes(),
    )?;
    Ok(None)
}

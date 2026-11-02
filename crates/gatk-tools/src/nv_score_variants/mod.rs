//! `NVScoreVariants`: the Python half of the tool, without Python.
//!
//! `nvscorevariants.py` encodes each variant, runs a pickled PyTorch model on batches of them, writes
//! a score for each to a temporary file with `%.3f`, and has `create_output_vcf.py` put that text
//! into a copy of the VCF through pysam. This is that, in pure Rust.
//!
//! What it does and does not reproduce is measured in `docs/nvscorevariants.md`: the printed scores
//! are what the conformance suite compares. The reference's own intermediate bits move with its
//! thread count and batch size, so they are a target only at one thread and one batch size.

pub mod encode;
pub mod kernels;
pub mod model;
pub mod vcf;
pub mod weights;

use self::encode::ReadTensorEncoder;
use self::model::Trace;
use self::vcf::{Header, Record};
use self::weights::Weights;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::reads::ReadsDataSource;
use std::fmt::Write as _;
use std::path::Path;

/// The arguments the Java half hands to the script.
#[derive(Clone, Debug)]
pub struct Arguments {
    pub output: String,
    pub variant: String,
    pub reference: String,
    pub input: Option<String>,
    pub read_tensor: bool,
    pub batch_size: i32,
    pub tmp_file: Option<String>,
    pub accelerator: String,
    pub weights_directory: String,
}

/// How the script ends other than well.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Stop {
    /// An exception or `sys.exit` in the reference: the Java half reports `Tool returned: 1`.
    Python(String),
    /// Something this port does not do, in its own words.
    Limitation(String),
}

impl From<String> for Stop {
    fn from(message: String) -> Stop {
        Stop::Python(message)
    }
}

/// The description `create_output_vcf.py` gives, with its source line's continuation in it.
const SCORE_DESCRIPTION: &str = "Log odds of being a true variant versus                     \
                                 being false under the trained Convolutional Neural Network";

/// `'%.3f' % float`, which is a correctly rounded decimal and spells the non-numbers its own way.
pub fn format_score(score: f32) -> String {
    if score.is_nan() {
        // Python prints no sign on a NaN.
        "nan".to_string()
    } else if score.is_infinite() {
        if score > 0.0 { "inf" } else { "-inf" }.to_string()
    } else {
        format!("{:.3}", f64::from(score))
    }
}

/// `float(text)`, into the float the record holds.
fn parse_score(text: &str) -> f32 {
    match text {
        "nan" => f32::NAN,
        "inf" => f32::INFINITY,
        "-inf" => f32::NEG_INFINITY,
        // What `format_score` wrote always parses.
        _ => text.parse::<f64>().unwrap_or(f64::NAN) as f32,
    }
}

/// What a VCF is once its header has been read.
struct Variants {
    header: Header,
    records: Vec<Record>,
}

/// `VariantFile(path)`, and a pass over all of it.
fn read_variants(path: &str) -> Result<Variants, Stop> {
    let bytes = std::fs::read(path).map_err(|e| Stop::Python(format!("OSError: {path}: {e}")))?;
    if bytes.starts_with(b"BCF") {
        return Err(Stop::Limitation(format!(
            "{path} is BCF, and this port reads VCF"
        )));
    }
    let bytes = if bytes.starts_with(&[0x1f, 0x8b]) {
        htsjdk_bgzf::read::decompress_all(&bytes)
            .map_err(|e| Stop::Python(format!("OSError: {path}: {e:?}")))?
    } else {
        bytes
    };
    let text = String::from_utf8_lossy(&bytes);
    let mut meta = Vec::new();
    let mut columns = None;
    let mut body = Vec::new();
    for line in text.lines() {
        if columns.is_some() {
            if !line.is_empty() {
                body.push(line);
            }
        } else if line.starts_with("#CHROM") {
            columns = Some(line);
        } else if line.starts_with("##") {
            meta.push(line);
        }
    }
    let Some(columns) = columns else {
        return Err(Stop::Python(format!(
            "ValueError: file `{path}` does not have valid header (mode=`r`) - is it VCF/BCF format?"
        )));
    };
    let header = Header::parse(&meta, columns);
    let mut records = Vec::with_capacity(body.len());
    for line in body {
        records.push(Record::parse(&header, line).map_err(Stop::Python)?);
    }
    Ok(Variants { header, records })
}

/// `htslib`'s two names for a BAM's index, in the order it looks.
fn find_index(bam: &str) -> Option<String> {
    let beside = format!("{bam}.bai");
    if Path::new(&beside).exists() {
        return Some(beside);
    }
    let without = bam.strip_suffix(".bam").map(|stem| format!("{stem}.bai"));
    without.filter(|path| Path::new(path).exists())
}

/// `AlignmentFile.fetch(contig, start, stop)`.
fn fetch(
    reads: &ReadsDataSource,
    record: &Record,
) -> Result<Vec<htsjdk_bam::record::BamRecord>, Stop> {
    let (begin, end) = (record.pos - 1, record.stop());
    let contig_exists = reads
        .header()
        .sequences
        .iter()
        .any(|s| s.name == record.contig);
    if !contig_exists {
        return Err(Stop::Python(format!(
            "ValueError: invalid contig `{}`",
            record.contig
        )));
    }
    // One position wider than asked: `bam_endpos` gives a read with nothing on the reference one
    // position and htsjdk gives it none, and `overlaps` is the test that decides.
    let interval = SimpleInterval {
        contig: record.contig.clone(),
        start: begin.max(1) as i32,
        end: end as i32,
    };
    let fetched = reads
        .query(std::slice::from_ref(&interval))
        .map_err(|e| Stop::Python(format!("OSError: {e:?}")))?;
    Ok(fetched
        .into_iter()
        .filter(|read| encode::overlaps(read, begin, end))
        .collect())
}

/// One batch's tensors.
pub struct Batch {
    pub types: Vec<i64>,
    pub best_practices: Vec<f32>,
    pub tensor: Vec<f32>,
}

/// The two ways to encode a variant.
enum Encoder {
    Reference,
    Reads {
        reads: Option<Box<ReadsDataSource>>,
        encoder: Box<ReadTensorEncoder>,
    },
}

/// The logits, then the softmax and the score, for one batch. With `trace.keep`, every layer.
pub fn score_batch(
    weights: &Weights,
    read_tensor: bool,
    batch: &Batch,
    threads: usize,
    trace: &mut Trace,
) -> Vec<f32> {
    let logits = if read_tensor {
        model::forward_2d(
            weights,
            &batch.tensor,
            &batch.best_practices,
            threads,
            trace,
        )
    } else {
        model::forward_1d(
            weights,
            &batch.tensor,
            &batch.best_practices,
            threads,
            trace,
        )
    };
    let probabilities = kernels::softmax(&logits);
    if trace.keep {
        trace
            .layers
            .insert("softmax".to_string(), probabilities.clone());
    }
    let scores = kernels::scores(&probabilities, &batch.types);
    if trace.keep {
        trace.layers.insert("score".to_string(), scores.clone());
    }
    scores
}

/// `variant_representation` for each record of a batch.
fn encode_batch(
    encoder: &mut Encoder,
    contigs: &std::collections::HashMap<String, Vec<u8>>,
    header: &Header,
    records: &[Record],
) -> Result<Batch, Stop> {
    let mut batch = Batch {
        types: Vec::with_capacity(records.len()),
        best_practices: Vec::with_capacity(records.len() * 7),
        tensor: Vec::new(),
    };
    for record in records {
        batch.types.push(encode::variant_type(record));
        batch
            .best_practices
            .extend_from_slice(&encode::annotations(record, &header.info)?);
        match encoder {
            Encoder::Reference => {
                batch
                    .tensor
                    .extend(encode::reference_tensor(contigs, record)?);
            }
            Encoder::Reads { reads, encoder } => {
                let Some(reads) = reads else {
                    return Err(Stop::Python(
                        "ValueError: fetch called on bamfile without index".to_string(),
                    ));
                };
                let fetched = fetch(reads, record)?;
                batch
                    .tensor
                    .extend(encoder.encode(contigs, record, fetched)?);
            }
        }
    }
    Ok(batch)
}

/// `main` of `nvscorevariants.py`.
pub fn run(arguments: &Arguments) -> Result<(), Stop> {
    let output = arguments.output.as_str();
    if output.ends_with(".gz") || output.ends_with(".bcf") {
        return Err(Stop::Limitation(format!(
            "{output}: this port writes the plain VCF the reference writes for other names, and has \
             not measured these"
        )));
    }

    // `TensorReader.from_files`, in its order.
    let variants = read_variants(&arguments.variant)?;
    let contigs = encode::read_fasta(&arguments.reference)?;
    let mut encoder = if arguments.read_tensor {
        let input = arguments.input.clone().unwrap_or_default();
        if !Path::new(&input).exists() {
            return Err(Stop::Python(format!(
                "FileNotFoundError: [Errno 2] could not open alignment file `{input}`"
            )));
        }
        let reads = match find_index(&input) {
            Some(index) => Some(Box::new(
                ReadsDataSource::open(Path::new(&input), Path::new(&index))
                    .map_err(|e| Stop::Python(format!("OSError: {input}: {e:?}")))?,
            )),
            None => None,
        };
        Encoder::Reads {
            reads,
            encoder: Box::new(ReadTensorEncoder::new()),
        }
    } else {
        Encoder::Reference
    };

    // `get_model`, which is where the reference would load its pickle.
    let pinned = if arguments.read_tensor {
        weights::TWO_D
    } else {
        weights::ONE_D
    };
    let weights = Weights::load(&arguments.weights_directory, pinned).map_err(Stop::Limitation)?;

    // `pl.Trainer(accelerator=...)`. The GPU backends are not ported, and on the machines this is
    // measured on they are not there for the reference either.
    if !matches!(arguments.accelerator.as_str(), "auto" | "cpu") {
        return Err(Stop::Python(format!(
            "MisconfigurationException: accelerator {} is not available",
            arguments.accelerator
        )));
    }
    // `DataLoader(batch_size=...)`.
    if arguments.batch_size <= 0 {
        return Err(Stop::Python(format!(
            "ValueError: batch_size should be a positive integer value, but got batch_size={}",
            arguments.batch_size
        )));
    }

    let threads = std::thread::available_parallelism().map_or(1, |n| n.get());
    let mut scored = String::new();
    let mut printed = Vec::with_capacity(variants.records.len());
    let tmp_file = arguments.tmp_file.clone();
    let flush = |text: &str| {
        if let Some(path) = &tmp_file {
            let _ = std::fs::write(path, text);
        }
    };
    let batches = variants.records.chunks(arguments.batch_size as usize);
    for records in batches {
        let batch = match encode_batch(&mut encoder, &contigs, &variants.header, records) {
            Ok(batch) => batch,
            Err(stop) => {
                flush(&scored);
                return Err(stop);
            }
        };
        let mut trace = Trace::default();
        let scores = score_batch(&weights, arguments.read_tensor, &batch, threads, &mut trace);
        for (record, score) in records.iter().zip(&scores) {
            let text = format_score(*score);
            let _ = writeln!(
                scored,
                "{}\t{}\t{}\t[{}]\t{}",
                record.contig,
                record.pos,
                record.reference,
                record.alternates.join(", "),
                text
            );
            printed.push(text);
        }
    }
    flush(&scored);

    create_output_vcf(arguments, variants, &printed)
}

/// `create_output_vcf.py`.
fn create_output_vcf(
    arguments: &Arguments,
    variants: Variants,
    printed: &[String],
) -> Result<(), Stop> {
    let label = if arguments.read_tensor {
        "CNN_2D"
    } else {
        "CNN_1D"
    };
    let Variants {
        mut header,
        records,
    } = variants;
    header.add_info(label, SCORE_DESCRIPTION);
    let mut text = header.format_text();
    let mut stop = None;
    if records.is_empty() {
        // `next(scoredVariants)` on nothing, and the header is all that was written.
        stop = Some("StopIteration".to_string());
    }
    for (mut record, score) in records.into_iter().zip(printed) {
        // `re.sub('[\[\]]', '', ...)` takes a breakend's brackets out of one side of the test.
        if record.alternates.iter().any(|alt| alt.contains(['[', ']'])) {
            stop = Some(format!(
                "Score file out of sync with original VCF at {}:{}",
                record.contig, record.pos
            ));
            break;
        }
        record.set_float(label, parse_score(score));
        match record.format_line(&header) {
            Ok(line) => text.push_str(&line),
            Err(why) => {
                stop = Some(format!("OSError: {why}"));
                break;
            }
        }
    }
    std::fs::write(&arguments.output, text)
        .map_err(|e| Stop::Python(format!("OSError: {}: {e}", arguments.output)))?;
    match stop {
        Some(why) => Err(Stop::Python(why)),
        None => Ok(()),
    }
}

/// The directory the weights are in: the variable, or beside the executable.
pub fn weights_directory() -> Option<String> {
    if let Some(directory) = std::env::var_os(weights::DIRECTORY_VARIABLE) {
        return directory.into_string().ok();
    }
    let executable = std::env::current_exe().ok()?;
    let directory = executable.parent()?;
    directory
        .join(weights::ONE_D.0)
        .exists()
        .then(|| directory.to_string_lossy().into_owned())
}

/// `encode_batch`, for the comparison with the reference's dump.
pub fn encode_for_comparison(arguments: &Arguments) -> Result<Vec<(Vec<Record>, Batch)>, Stop> {
    let variants = read_variants(&arguments.variant)?;
    let contigs = encode::read_fasta(&arguments.reference)?;
    let mut encoder = match (&arguments.input, arguments.read_tensor) {
        (Some(input), true) => {
            let index = find_index(input).ok_or_else(|| Stop::Python("no index".to_string()))?;
            Encoder::Reads {
                reads: Some(Box::new(
                    ReadsDataSource::open(Path::new(input), Path::new(&index))
                        .map_err(|e| Stop::Python(format!("{e:?}")))?,
                )),
                encoder: Box::new(ReadTensorEncoder::new()),
            }
        }
        _ => Encoder::Reference,
    };
    let mut out = Vec::new();
    for records in variants.records.chunks(arguments.batch_size as usize) {
        let batch = encode_batch(&mut encoder, &contigs, &variants.header, records)?;
        out.push((records.to_vec(), batch));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_score_prints_as_percent_three_f_does() {
        assert_eq!(format_score(-14.168865), "-14.169");
        assert_eq!(format_score(0.0625), "0.062");
        assert_eq!(format_score(0.1875), "0.188");
        assert_eq!(format_score(-0.0001), "-0.000");
        assert_eq!(format_score(f32::NAN), "nan");
        assert_eq!(format_score(f32::INFINITY), "inf");
        assert_eq!(format_score(f32::NEG_INFINITY), "-inf");
    }
}

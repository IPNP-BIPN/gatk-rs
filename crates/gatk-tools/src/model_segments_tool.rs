//! `ModelSegments.doWork`, end to end: the inputs read and checked, the hets genotyped, the kernel
//! segmentation or the given segments, the two models fit and smoothed, and every output file.
//!
//! Ported from `org.broadinstitute.hellbender.tools.copynumber.ModelSegments` (GATK 4.6.2.0), with
//! `MultisampleMultidimensionalKernelSegmenter`, `NaiveHeterozygousPileupGenotypingUtils` and the
//! copy-number collections' readers and writers it goes through.
//!
//! # Two run modes, three data modes
//!
//! More than one file in `--denoised-copy-ratios` or `--allelic-counts` is the multi-sample mode,
//! which segments jointly and writes one Picard interval list. One sample is modelled: its hets
//! written, segmented (or `--segments` read), the copy-ratio and allele-fraction chains fit, the
//! segments smoothed, and ten files written. Either input may be absent, and the missing one is
//! imputed as an empty collection carrying the other's metadata, which is why the copy-ratio-only
//! run still writes allele-fraction files and the allele-fraction-only run copy-ratio ones.
//!
//! # The order is the reference's
//!
//! `setModesAndValidateArguments` runs first, output directory included, then every input is read,
//! then the cross-input checks, then the genotyping. A refused row is compared by its first line, so
//! the check that fires first is the answer.
//!
//! # One `HashSet` reaches a byte
//!
//! The `.cr.seg` mean of a segment averages the copy ratios `OverlapDetector.getOverlaps` returns,
//! a `HashSet`, with a compensated stream sum. The order of the terms is that set's iteration
//! order, so it is taken from [`gatk_engine::java_hash::hash_map_order`] over each ratio's
//! `CopyRatio.hashCode`, inserted in the interval tree's order, which is by start.

use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};

use gatk_engine::kernel_segmenter::{find_changepoints_of, ChangepointSortOrder};

use crate::main_entry::{Failure, Thrown};
use crate::model_segments_models::{
    illegal, ChainLengths, Het, ModeledSegment, MultidimensionalModeller, ParameterDeciles, Span,
};

/// Everything `ModelSegments` reads from its command line.
#[derive(Debug, Clone)]
pub struct Options {
    pub denoised_copy_ratios: Vec<String>,
    pub allelic_counts: Vec<String>,
    pub normal_allelic_counts: Option<String>,
    pub segments: Option<String>,
    pub output_prefix: String,
    pub output_dir: String,
    pub minimum_total_allele_count_case: i32,
    pub minimum_total_allele_count_normal: i32,
    pub genotyping_homozygous_log_ratio_threshold: f64,
    pub genotyping_base_error_rate: f64,
    pub maximum_number_of_segments_per_chromosome: i32,
    pub kernel_variance_copy_ratio: f64,
    pub kernel_variance_allele_fraction: f64,
    pub kernel_scaling_allele_fraction: f64,
    pub kernel_approximation_dimension: i32,
    pub window_sizes: Vec<i32>,
    pub number_of_changepoints_penalty_factor: f64,
    pub minor_allele_fraction_prior_alpha: f64,
    pub number_of_samples_copy_ratio: i32,
    pub number_of_burn_in_samples_copy_ratio: i32,
    pub number_of_samples_allele_fraction: i32,
    pub number_of_burn_in_samples_allele_fraction: i32,
    pub smoothing_credible_interval_threshold_copy_ratio: f64,
    pub smoothing_credible_interval_threshold_allele_fraction: f64,
    pub maximum_number_of_smoothing_iterations: i32,
    pub number_of_smoothing_iterations_per_fit: i32,
}

/// `SampleLocatableMetadata`: a sample name and a sequence dictionary.
#[derive(Debug, Clone, PartialEq)]
struct Metadata {
    sample: String,
    dictionary: Vec<SequenceRecord>,
}

/// A copy-ratio file.
#[derive(Debug, Clone)]
struct CopyRatios {
    metadata: Metadata,
    records: Vec<(Span, f64)>,
}

/// An allelic-count file, with the two nucleotides kept as `Nucleotide.name()` prints them.
#[derive(Debug, Clone)]
struct AllelicCounts {
    metadata: Metadata,
    records: Vec<(Het, String, String)>,
}

impl AllelicCounts {
    fn sites(&self) -> Vec<(String, i32)> {
        self.records
            .iter()
            .map(|(het, _, _)| (het.contig.clone(), het.position))
            .collect()
    }

    fn hets(&self) -> Vec<Het> {
        self.records.iter().map(|(het, _, _)| het.clone()).collect()
    }

    fn filtered(&self, keep: impl Fn(&Het) -> bool) -> AllelicCounts {
        AllelicCounts {
            metadata: self.metadata.clone(),
            records: self
                .records
                .iter()
                .filter(|(het, _, _)| keep(het))
                .cloned()
                .collect(),
        }
    }
}

const READ_GROUP_ID: &str = "GATKCopyNumber";

/// `metadata.toHeader().getSAMString()` for a sample and, when there is one, a dictionary.
fn header(sample: &str, dictionary: Option<&[SequenceRecord]>) -> String {
    let mut header = SamHeader::default();
    if let Some(sequences) = dictionary {
        header.sequences = sequences.to_vec();
    }
    let mut read_group = ReadGroup::new(READ_GROUP_ID);
    read_group.attributes.set("SM", sample);
    header.read_groups.push(read_group);
    header.encode()
}

/// `CopyNumberFormatsUtils.formatDouble`, `%.6f`.
fn format_double(value: f64) -> String {
    gatk_engine::java_format::format_decimals(value, 6)
}

/// A table row, tab-separated, with the quoting opencsv applies where a value needs it.
fn row(values: &[String]) -> String {
    let quoted: Vec<String> = values
        .iter()
        .map(|value| gatk_engine::tsv_table::quote_if_needed(value))
        .collect();
    format!("{}\n", quoted.join("\t"))
}

fn write(path: &str, text: &str) -> Result<(), Thrown> {
    std::fs::write(path, text).map_err(|error| Thrown {
        failure: Failure::User,
        exception: crate::main_entry::USER_EXCEPTION,
        message: Some(format!("Couldn't write file {path} because exception {error}")),
    })
}

/// A copy-number table read the way `AbstractRecordCollection` reads one: the SAM header, the
/// column line, and the rows under it.
struct Table {
    header: SamHeader,
    columns: Vec<String>,
    rows: Vec<Vec<String>>,
}

fn read_table(path: &str) -> Result<Table, Thrown> {
    let text = std::fs::read_to_string(path).map_err(|_| {
        Thrown::user(format!(
            "Couldn't read file {path}. Error was: The input file does not exist."
        ))
    })?;
    let header = htsjdk_bam::reader::parse_header_text(
        &text
            .lines()
            .take_while(|line| line.starts_with('@'))
            .map(|line| format!("{line}\n"))
            .collect::<String>(),
    );
    let mut lines = text
        .lines()
        .filter(|line| !line.starts_with('@') && !line.is_empty());
    let columns: Vec<String> = lines
        .next()
        .unwrap_or_default()
        .split('\t')
        .map(str::to_string)
        .collect();
    let rows = lines
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect();
    Ok(Table {
        header,
        columns,
        rows,
    })
}

impl Table {
    fn column(&self, name: &str) -> Result<usize, Thrown> {
        self.columns
            .iter()
            .position(|column| column == name)
            .ok_or_else(|| {
                Thrown::user(format!(
                    "Bad input: Missing mandatory column {name} in the input file."
                ))
            })
    }

    /// `MetadataUtils.fromHeader(header, SAMPLE_LOCATABLE)`.
    fn metadata(&self) -> Result<Metadata, Thrown> {
        if self.header.read_groups.is_empty() {
            return Err(illegal(
                "The input header does not contain any read groups.  Cannot determine a sample name.",
            ));
        }
        let mut samples: Vec<String> = Vec::new();
        for read_group in &self.header.read_groups {
            let sample = read_group.attributes.get("SM").unwrap_or("null").to_string();
            if !samples.contains(&sample) {
                samples.push(sample);
            }
        }
        if samples.len() > 1 {
            return Err(illegal(format!(
                "The input header contains more than one unique sample name: {}",
                samples.join(", ")
            )));
        }
        Ok(Metadata {
            sample: samples.remove(0),
            dictionary: self.header.sequences.clone(),
        })
    }
}

fn parse_int(text: &str) -> Result<i32, Thrown> {
    text.trim()
        .parse()
        .map_err(|_| Thrown::user(format!("Bad input: Expected an integer value but found {text}")))
}

fn parse_double(text: &str) -> Result<f64, Thrown> {
    let trimmed = text.trim();
    match trimmed {
        "NaN" => Ok(f64::NAN),
        "Infinity" | "+Infinity" => Ok(f64::INFINITY),
        "-Infinity" => Ok(f64::NEG_INFINITY),
        _ => trimmed
            .parse()
            .map_err(|_| Thrown::user(format!("Bad input: Expected a double value but found {text}"))),
    }
}

/// `new CopyRatioCollection(file)`.
fn read_copy_ratios(path: &str) -> Result<CopyRatios, Thrown> {
    let table = read_table(path)?;
    let metadata = table.metadata()?;
    let (contig, start, end, value) = (
        table.column("CONTIG")?,
        table.column("START")?,
        table.column("END")?,
        table.column("LOG2_COPY_RATIO")?,
    );
    let mut records = Vec::with_capacity(table.rows.len());
    for fields in &table.rows {
        records.push((
            Span {
                contig: fields[contig].clone(),
                start: parse_int(&fields[start])?,
                end: parse_int(&fields[end])?,
            },
            parse_double(&fields[value])?,
        ));
    }
    Ok(CopyRatios { metadata, records })
}

/// `Nucleotide.decode(char).name()`.
fn nucleotide(text: &str) -> String {
    match text.chars().next().map(|c| c.to_ascii_uppercase()) {
        Some(
            c @ ('A' | 'C' | 'G' | 'T' | 'R' | 'Y' | 'S' | 'W' | 'K' | 'M' | 'B' | 'D' | 'H' | 'V'
            | 'N' | 'X' | 'U'),
        ) => c.to_string(),
        _ => "X".to_string(),
    }
}

/// `new AllelicCountCollection(file)`.
fn read_allelic_counts(path: &str) -> Result<AllelicCounts, Thrown> {
    let table = read_table(path)?;
    let metadata = table.metadata()?;
    let (contig, position, reference, alternate, reference_base, alternate_base) = (
        table.column("CONTIG")?,
        table.column("POSITION")?,
        table.column("REF_COUNT")?,
        table.column("ALT_COUNT")?,
        table.column("REF_NUCLEOTIDE")?,
        table.column("ALT_NUCLEOTIDE")?,
    );
    let mut records = Vec::with_capacity(table.rows.len());
    for fields in &table.rows {
        records.push((
            Het {
                contig: fields[contig].clone(),
                position: parse_int(&fields[position])?,
                ref_count: parse_int(&fields[reference])?,
                alt_count: parse_int(&fields[alternate])?,
            },
            nucleotide(&fields[reference_base]),
            nucleotide(&fields[alternate_base]),
        ));
    }
    Ok(AllelicCounts { metadata, records })
}

/// `IntervalList.fromFile`: the header's dictionary and the intervals under it.
fn read_interval_list(path: &str) -> Result<(Vec<SequenceRecord>, Vec<Span>), Thrown> {
    let text = std::fs::read_to_string(path).map_err(|_| {
        Thrown::user(format!(
            "Couldn't read file {path}. Error was: The input file does not exist."
        ))
    })?;
    let header = htsjdk_bam::reader::parse_header_text(
        &text
            .lines()
            .take_while(|line| line.starts_with('@'))
            .map(|line| format!("{line}\n"))
            .collect::<String>(),
    );
    let mut intervals = Vec::new();
    for line in text.lines().filter(|line| !line.starts_with('@') && !line.is_empty()) {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }
        intervals.push(Span {
            contig: fields[0].to_string(),
            start: parse_int(fields[1])?,
            end: parse_int(fields[2])?,
        });
    }
    Ok((header.sequences, intervals))
}

/// `IntervalList.write`: the dictionary header, then `contig start end + .` per interval.
fn interval_list_text(dictionary: &[SequenceRecord], intervals: &[Span]) -> String {
    let mut header = SamHeader::default();
    header.sequences = dictionary.to_vec();
    let mut text = header.encode_replacing_version();
    for interval in intervals {
        text.push_str(&format!(
            "{}\t{}\t{}\t+\t.\n",
            interval.contig, interval.start, interval.end
        ));
    }
    text
}

/// `AllelicCountCollection.write`.
fn allelic_counts_text(counts: &AllelicCounts) -> String {
    let mut text = header(&counts.metadata.sample, Some(&counts.metadata.dictionary));
    text.push_str("CONTIG\tPOSITION\tREF_COUNT\tALT_COUNT\tREF_NUCLEOTIDE\tALT_NUCLEOTIDE\n");
    for (het, reference, alternate) in &counts.records {
        text.push_str(&row(&[
            het.contig.clone(),
            het.position.to_string(),
            het.ref_count.to_string(),
            het.alt_count.to_string(),
            reference.clone(),
            alternate.clone(),
        ]));
    }
    text
}

/// `ModeledSegmentCollection.write`.
fn modeled_segments_text(metadata: &Metadata, segments: &[ModeledSegment]) -> String {
    let mut text = header(&metadata.sample, Some(&metadata.dictionary));
    text.push_str(
        "CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tNUM_POINTS_ALLELE_FRACTION\t\
         LOG2_COPY_RATIO_POSTERIOR_10\tLOG2_COPY_RATIO_POSTERIOR_50\tLOG2_COPY_RATIO_POSTERIOR_90\t\
         MINOR_ALLELE_FRACTION_POSTERIOR_10\tMINOR_ALLELE_FRACTION_POSTERIOR_50\t\
         MINOR_ALLELE_FRACTION_POSTERIOR_90\n",
    );
    for segment in segments {
        text.push_str(&row(&[
            segment.span.contig.clone(),
            segment.span.start.to_string(),
            segment.span.end.to_string(),
            segment.num_points_copy_ratio.to_string(),
            segment.num_points_allele_fraction.to_string(),
            format_double(segment.log2_copy_ratio.decile10),
            format_double(segment.log2_copy_ratio.decile50),
            format_double(segment.log2_copy_ratio.decile90),
            format_double(segment.minor_allele_fraction.decile10),
            format_double(segment.minor_allele_fraction.decile50),
            format_double(segment.minor_allele_fraction.decile90),
        ]));
    }
    text
}

/// `ParameterDecileCollection.write`, whose metadata is the sample alone.
fn parameters_text(sample: &str, deciles: &ParameterDeciles) -> String {
    let mut text = header(sample, None);
    text.push_str(
        "PARAMETER_NAME\tPOSTERIOR_10\tPOSTERIOR_20\tPOSTERIOR_30\tPOSTERIOR_40\tPOSTERIOR_50\t\
         POSTERIOR_60\tPOSTERIOR_70\tPOSTERIOR_80\tPOSTERIOR_90\n",
    );
    for (name, values) in deciles {
        let mut fields = vec![name.to_string()];
        fields.extend(values.iter().map(|value| format_double(*value)));
        text.push_str(&row(&fields));
    }
    text
}

/// `LegacySegmentCollection.write`, which has no header at all.
fn legacy_text(sample: &str, segments: &[(Span, i32, f64)]) -> String {
    let mut text = String::from("Sample\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n");
    for (span, num_probes, mean) in segments {
        text.push_str(&row(&[
            sample.to_string(),
            span.contig.clone(),
            span.start.to_string(),
            span.end.to_string(),
            num_probes.to_string(),
            format_double(*mean),
        ]));
    }
    text
}

/// `SimpleInterval.hashCode`.
fn interval_hash_code(contig: &str, start: i32, end: i32) -> i32 {
    let mut result = start;
    result = result.wrapping_mul(31).wrapping_add(end);
    result
        .wrapping_mul(31)
        .wrapping_add(gatk_engine::java_hash::string_hash_code(contig))
}

/// `CopyRatio.hashCode` of a midpoint ratio.
fn midpoint_hash_code(span: &Span, value: f64) -> i32 {
    let midpoint = (span.start + span.end) / 2;
    let bits = if value.is_nan() {
        0x7ff8_0000_0000_0000_u64
    } else {
        value.to_bits()
    };
    let folded = (bits ^ (bits >> 32)) as u32 as i32;
    interval_hash_code(&span.contig, midpoint, midpoint)
        .wrapping_mul(31)
        .wrapping_add(folded)
}

/// `CopyRatioSegmentCollection.write`, one segment per modelled segment, whose mean is over the
/// copy ratios with a midpoint inside it in the order their `HashSet` iterates.
fn copy_ratio_segments_text(
    metadata: &Metadata,
    segments: &[ModeledSegment],
    copy_ratios: &[(Span, f64)],
) -> Result<String, Thrown> {
    let mut text = header(&metadata.sample, Some(&metadata.dictionary));
    text.push_str("CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tMEAN_LOG2_COPY_RATIO\n");
    for segment in segments {
        let inside: Vec<(usize, i32)> = copy_ratios
            .iter()
            .enumerate()
            .filter(|(_, (span, _))| {
                let midpoint = (span.start + span.end) / 2;
                span.contig == segment.span.contig
                    && segment.span.start <= midpoint
                    && midpoint <= segment.span.end
            })
            .map(|(index, (span, value))| (index, midpoint_hash_code(span, *value)))
            .collect();
        let order = gatk_engine::java_hash::hash_map_order(&inside).map_err(|_| Thrown {
            failure: Failure::Other,
            exception: crate::main_entry::PORT_LIMITATION,
            message: Some(
                "a copy-ratio segment whose HashSet crowds a bucket past the measured layout"
                    .to_string(),
            ),
        })?;
        let values: Vec<f64> = order.iter().map(|index| copy_ratios[*index].1).collect();
        let mean = gatk_engine::copy_number_mcmc::double_stream_average(&values).unwrap_or(f64::NAN);
        text.push_str(&row(&[
            segment.span.contig.clone(),
            segment.span.start.to_string(),
            segment.span.end.to_string(),
            values.len().to_string(),
            format_double(mean),
        ]));
    }
    Ok(text)
}

/// `NaiveHeterozygousPileupGenotypingUtils.calculateHomozygousLogRatio` against the threshold.
fn is_heterozygous(het: &Het, threshold: f64, base_error_rate: f64) -> bool {
    let count = crate::model_segments::AllelicCount {
        position: het.position,
        reference_count: het.ref_count,
        alternate_count: het.alt_count,
    };
    crate::model_segments::homozygous_log_ratio(count, base_error_rate).unwrap_or(f64::NAN)
        < threshold
}

/// `filterByOverlap` against copy-ratio intervals: a site inside any interval.
fn overlaps_any(het: &Het, intervals: &[Span]) -> bool {
    intervals
        .iter()
        .any(|interval| interval.contains_position(&het.contig, het.position))
}

/// `genotypeHets`: the case samples' hets, and the normal's when there is one.
fn genotype_hets(
    options: &Options,
    cases: &[AllelicCounts],
    normal: Option<&AllelicCounts>,
    copy_ratio_intervals: &[Span],
) -> (Vec<AllelicCounts>, Option<AllelicCounts>) {
    let threshold = options.genotyping_homozygous_log_ratio_threshold;
    let error_rate = options.genotyping_base_error_rate;
    let het_normal = normal.map(|normal| {
        let minimum = options.minimum_total_allele_count_normal;
        let mut filtered = if minimum == 0 {
            normal.clone()
        } else {
            normal.filtered(|het| het.ref_count + het.alt_count >= minimum)
        };
        if !copy_ratio_intervals.is_empty() {
            filtered = filtered.filtered(|het| overlaps_any(het, copy_ratio_intervals));
        }
        filtered.filtered(|het| is_heterozygous(het, threshold, error_rate))
    });
    let mut hets: Vec<AllelicCounts> = cases
        .iter()
        .map(|counts| {
            let minimum = options.minimum_total_allele_count_case;
            let mut filtered = if minimum == 0 {
                counts.clone()
            } else {
                counts.filtered(|het| het.ref_count + het.alt_count >= minimum)
            };
            if !copy_ratio_intervals.is_empty() {
                filtered = filtered.filtered(|het| overlaps_any(het, copy_ratio_intervals));
            }
            match &het_normal {
                Some(normal_hets) => {
                    let sites = normal_hets.sites();
                    filtered.filtered(|het| {
                        !sites.is_empty() && sites.contains(&(het.contig.clone(), het.position))
                    })
                }
                None => filtered.filtered(|het| is_heterozygous(het, threshold, error_rate)),
            }
        })
        .collect();
    if hets.len() > 1 {
        let first = hets[0].sites();
        let common: Vec<(String, i32)> = first
            .into_iter()
            .filter(|site| hets.iter().all(|counts| counts.sites().contains(site)))
            .collect();
        hets = hets
            .iter()
            .map(|counts| counts.filtered(|het| common.contains(&(het.contig.clone(), het.position))))
            .collect();
    }
    (hets, het_normal)
}

/// One point of `MultisampleMultidimensionalKernelSegmenter`: the interval, and a copy ratio and
/// an alternate-allele fraction per sample.
struct Point {
    span: Span,
    copy_ratios: Vec<f64>,
    fractions: Vec<f64>,
}

/// `AllelicCount.getAlternateAlleleFraction`.
fn alternate_fraction(het: &Het) -> f64 {
    let total = het.ref_count + het.alt_count;
    if total == 0 {
        0.0
    } else {
        f64::from(het.alt_count) / f64::from(total)
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum SegmenterMode {
    CopyRatioOnly,
    AlleleFractionOnly,
    Both,
}

/// `KERNEL.apply(sd).apply(x, y)`: linear at zero, else a normal density centred on `x`.
fn kernel_value(standard_deviation: f64, x: f64, y: f64) -> f64 {
    if standard_deviation == 0.0 {
        x * y
    } else {
        gatk_engine::copy_number_mcmc::normal_density(x, standard_deviation, y)
    }
}

/// `MultisampleMultidimensionalKernelSegmenter.findSegmentation`.
fn find_segmentation(
    options: &Options,
    copy_ratios: &[Vec<(Span, f64)>],
    hets: &[Vec<Het>],
) -> Result<Vec<Span>, Thrown> {
    if options.window_sizes.iter().any(|size| *size <= 0) {
        return Err(illegal("Window sizes must all be positive."));
    }
    let num_copy_ratio = copy_ratios[0].len();
    let num_allele_fraction = hets[0].len();
    let (mode, points): (SegmenterMode, Vec<Point>) = if num_allele_fraction == 0 {
        (
            SegmenterMode::CopyRatioOnly,
            (0..num_copy_ratio)
                .map(|i| Point {
                    span: copy_ratios[0][i].0.clone(),
                    copy_ratios: copy_ratios.iter().map(|sample| sample[i].1).collect(),
                    fractions: Vec::new(),
                })
                .collect(),
        )
    } else if num_copy_ratio == 0 {
        (
            SegmenterMode::AlleleFractionOnly,
            (0..num_allele_fraction)
                .map(|i| Point {
                    span: Span {
                        contig: hets[0][i].contig.clone(),
                        start: hets[0][i].position,
                        end: hets[0][i].position,
                    },
                    copy_ratios: Vec::new(),
                    fractions: hets.iter().map(|sample| alternate_fraction(&sample[i])).collect(),
                })
                .collect(),
        )
    } else {
        // The first het inside each copy-ratio interval, or a balanced count where there is none.
        let site_of: Vec<Option<usize>> = copy_ratios[0]
            .iter()
            .map(|(span, _)| {
                hets[0].iter().position(|het| {
                    het.contig == span.contig && span.start <= het.position && het.position <= span.end
                })
            })
            .collect();
        (
            SegmenterMode::Both,
            (0..num_copy_ratio)
                .map(|i| Point {
                    span: copy_ratios[0][i].0.clone(),
                    copy_ratios: copy_ratios.iter().map(|sample| sample[i].1).collect(),
                    fractions: hets
                        .iter()
                        .map(|sample| match site_of[i] {
                            Some(index) => alternate_fraction(&sample[index]),
                            None => 0.5,
                        })
                        .collect(),
                })
                .collect(),
        )
    };

    let num_samples = copy_ratios.len();
    let sd_copy_ratio = options.kernel_variance_copy_ratio.sqrt();
    let sd_allele_fraction = options.kernel_variance_allele_fraction.sqrt();
    let scaling = options.kernel_scaling_allele_fraction;
    let kernel = move |p1: &Point, p2: &Point| -> f64 {
        let mut sum = 0.0;
        for sample in 0..num_samples {
            match mode {
                SegmenterMode::CopyRatioOnly => {
                    sum += kernel_value(sd_copy_ratio, p1.copy_ratios[sample], p2.copy_ratios[sample]);
                }
                SegmenterMode::AlleleFractionOnly => {
                    sum += kernel_value(sd_allele_fraction, p1.fractions[sample], p2.fractions[sample]);
                }
                SegmenterMode::Both => {
                    sum += kernel_value(sd_copy_ratio, p1.copy_ratios[sample], p2.copy_ratios[sample])
                        + scaling
                            * kernel_value(
                                sd_allele_fraction,
                                p1.fractions[sample],
                                p2.fractions[sample],
                            );
                }
            }
        }
        sum
    };

    let max_changepoints = (options.maximum_number_of_segments_per_chromosome - 1).max(0) as usize;
    let window_sizes: Vec<usize> = crate::model_segments::window_sizes(&options.window_sizes)
        .into_iter()
        .map(|size| size as usize)
        .collect();

    // `Collectors.groupingBy(contig, LinkedHashMap::new, toList())`: contigs in first appearance.
    let mut contigs: Vec<String> = Vec::new();
    for point in &points {
        if !contigs.contains(&point.span.contig) {
            contigs.push(point.span.contig.clone());
        }
    }
    let mut segments = Vec::new();
    for contig in contigs {
        let chromosome: Vec<&Point> = points.iter().filter(|p| p.span.contig == contig).collect();
        let n = chromosome.len();
        if n < crate::model_segments::MINIMUM_POINTS_REQUIRED_PER_CHROMOSOME {
            segments.push(Span {
                contig: contig.clone(),
                start: chromosome[0].span.start,
                end: chromosome[n - 1].span.end,
            });
            continue;
        }
        if window_sizes.is_empty() {
            return Err(illegal("At least one window size must be provided."));
        }
        let changepoints = find_changepoints_of(
            &chromosome,
            max_changepoints,
            |a: &&Point, b: &&Point| kernel(a, b),
            options.kernel_approximation_dimension.max(1) as usize,
            &window_sizes,
            options.number_of_changepoints_penalty_factor,
            options.number_of_changepoints_penalty_factor,
            ChangepointSortOrder::Index,
        );
        let bounds: Vec<(i32, i32)> = chromosome
            .iter()
            .map(|point| (point.span.start, point.span.end))
            .collect();
        for (start, end) in crate::model_segments::segments_from_changepoints(&bounds, &changepoints) {
            segments.push(Span {
                contig: contig.clone(),
                start,
                end,
            });
        }
    }
    Ok(segments)
}

/// The dictionary-order check `AbstractLocatableCollection` makes of the segments it is given,
/// after sorting them, which is the one a segments file can fail.
fn sorted_segments(dictionary: &[SequenceRecord], mut segments: Vec<Span>) -> Result<Vec<Span>, Thrown> {
    let index_of = |contig: &str| dictionary.iter().position(|s| s.name == contig);
    if segments.iter().any(|s| {
        index_of(&s.contig).is_none_or(|i| s.start < 1 || s.end > dictionary[i].length || s.start > s.end)
    }) {
        return Err(illegal(
            "Records contained at least one interval that did not validate against the sequence dictionary.",
        ));
    }
    segments.sort_by_key(|s| (index_of(&s.contig), s.start, s.end));
    for pair in segments.windows(2) {
        if pair[0] == pair[1] {
            return Err(illegal("Records were not strictly sorted in dictionary order."));
        }
        if pair[0].contig == pair[1].contig && pair[1].start <= pair[0].end {
            return Err(illegal(format!(
                "Records contain at least two overlapping intervals: {} and {}",
                pair[0].render(),
                pair[1].render()
            )));
        }
    }
    Ok(segments)
}

/// `ModelSegments.doWork`.
pub fn run(options: &Options) -> Result<(), Thrown> {
    // `setModesAndValidateArguments`.
    if options.output_prefix.is_empty() {
        return Err(illegal("The collection is empty: null"));
    }
    let output_dir = options.output_dir.trim_end_matches('/').to_string();
    let output_dir = if output_dir.is_empty() {
        "/".to_string()
    } else {
        output_dir
    };
    if std::fs::metadata(&output_dir).is_err() && std::fs::create_dir_all(&output_dir).is_err() {
        return Err(Thrown::user(format!(
            "Couldn't write file {output_dir} because : The output directory does not exist and could not be created."
        )));
    }
    let has_copy_ratios = !options.denoised_copy_ratios.is_empty();
    let has_allelic_counts = !options.allelic_counts.is_empty();
    if !has_copy_ratios && !has_allelic_counts {
        return Err(illegal(
            "Must provide at least one denoised-copy-ratios file or allelic-counts file.",
        ));
    }
    if !has_allelic_counts && options.normal_allelic_counts.is_some() {
        return Err(illegal(
            "Must provide an allelic-counts file for the case sample to run in matched-normal mode.",
        ));
    }
    if options.normal_allelic_counts.is_some() && options.minimum_total_allele_count_case > 0 {
        return Err(illegal(
            "The minimum total count for filtering allelic counts in case samples must be set to zero in matched-normal mode. \
             If the effect of statistical noise due to low depth in case samples on segmentation is a concern, \
             consider using only denoised copy ratios or externally preprocessing allelic-count files \
             to remove sites that are poorly covered across all samples.",
        ));
    }
    let multiple = options.denoised_copy_ratios.len() > 1 || options.allelic_counts.len() > 1;
    if multiple {
        if has_copy_ratios
            && has_allelic_counts
            && options.denoised_copy_ratios.len() != options.allelic_counts.len()
        {
            return Err(illegal(
                "Number of denoised-copy-ratios and allelic-counts files for the case samples must be equal \
                 if both input types are specified in multisample mode.",
            ));
        }
        if options.segments.is_some() {
            return Err(illegal("Segments file cannot be specified in multisample mode."));
        }
    }
    if options.number_of_samples_copy_ratio <= options.number_of_burn_in_samples_copy_ratio {
        return Err(illegal(
            "Number of copy-ratio samples must be greater than number of copy-ratio burn-in samples.",
        ));
    }
    if options.number_of_samples_allele_fraction <= options.number_of_burn_in_samples_allele_fraction {
        return Err(illegal(
            "Number of allele-fraction samples must be greater than number of allele-fraction burn-in samples.",
        ));
    }

    // `ModelSegmentsData`: read what there is and impute the rest.
    let (copy_ratios, allelic_counts): (Vec<CopyRatios>, Vec<AllelicCounts>) =
        match (has_copy_ratios, has_allelic_counts) {
            (true, false) => {
                let ratios = options
                    .denoised_copy_ratios
                    .iter()
                    .map(|path| read_copy_ratios(path))
                    .collect::<Result<Vec<_>, _>>()?;
                let counts = ratios
                    .iter()
                    .map(|r| AllelicCounts {
                        metadata: r.metadata.clone(),
                        records: Vec::new(),
                    })
                    .collect();
                (ratios, counts)
            }
            (false, true) => {
                let counts = options
                    .allelic_counts
                    .iter()
                    .map(|path| read_allelic_counts(path))
                    .collect::<Result<Vec<_>, _>>()?;
                let ratios = counts
                    .iter()
                    .map(|c| CopyRatios {
                        metadata: c.metadata.clone(),
                        records: Vec::new(),
                    })
                    .collect();
                (ratios, counts)
            }
            _ => {
                let ratios = options
                    .denoised_copy_ratios
                    .iter()
                    .map(|path| read_copy_ratios(path))
                    .collect::<Result<Vec<_>, _>>()?;
                let counts = options
                    .allelic_counts
                    .iter()
                    .map(|path| read_allelic_counts(path))
                    .collect::<Result<Vec<_>, _>>()?;
                for (ratio, count) in ratios.iter().zip(&counts) {
                    if ratio.metadata != count.metadata {
                        return Err(illegal("Metadata do not match."));
                    }
                }
                (ratios, counts)
            }
        };
    let normal = match &options.normal_allelic_counts {
        Some(path) => Some(read_allelic_counts(path)?),
        None => None,
    };
    if let Some(path) = &options.segments {
        let (dictionary, intervals) = read_interval_list(path)?;
        sorted_segments(&dictionary, intervals)?;
    }

    let first_intervals: Vec<Span> = copy_ratios[0].records.iter().map(|r| r.0.clone()).collect();
    if copy_ratios
        .iter()
        .any(|r| r.records.iter().map(|x| &x.0).ne(first_intervals.iter()))
    {
        return Err(illegal("Copy-ratio intervals must be identical across all case samples."));
    }
    let first_sites = allelic_counts[0].sites();
    if allelic_counts
        .iter()
        .chain(normal.iter())
        .any(|counts| counts.sites() != first_sites)
    {
        return Err(illegal("Allelic-count sites must be identical across all samples."));
    }

    let (hets, het_normal) = genotype_hets(options, &allelic_counts, normal.as_ref(), &first_intervals);
    let path = |suffix: &str| format!("{output_dir}/{}{suffix}", options.output_prefix);

    if multiple {
        let segments = find_segmentation(
            options,
            &copy_ratios.iter().map(|r| r.records.clone()).collect::<Vec<_>>(),
            &hets.iter().map(AllelicCounts::hets).collect::<Vec<_>>(),
        )?;
        write(
            &path(".interval_list"),
            &interval_list_text(&copy_ratios[0].metadata.dictionary, &segments),
        )?;
        return Ok(());
    }

    let denoised = &copy_ratios[0];
    let het_counts = &hets[0];
    let metadata = &denoised.metadata;
    if has_allelic_counts {
        if let Some(normal_hets) = &het_normal {
            write(&path(".hets.normal.tsv"), &allelic_counts_text(normal_hets))?;
        }
        write(&path(".hets.tsv"), &allelic_counts_text(het_counts))?;
    }

    let segments = match &options.segments {
        None => find_segmentation(
            options,
            std::slice::from_ref(&denoised.records),
            &[het_counts.hets()],
        )?,
        Some(file) => {
            let (_, intervals) = read_interval_list(file)?;
            if intervals.is_empty() {
                return Err(illegal("Segments file must contain at least one segment."));
            }
            sorted_segments(&metadata.dictionary, intervals)?
        }
    };

    let het_records = het_counts.hets();
    let mut modeller = MultidimensionalModeller::new(
        segments,
        &denoised.records,
        &het_records,
        options.minor_allele_fraction_prior_alpha,
        ChainLengths {
            num_samples_copy_ratio: options.number_of_samples_copy_ratio as usize,
            num_burn_in_copy_ratio: options.number_of_burn_in_samples_copy_ratio as usize,
            num_samples_allele_fraction: options.number_of_samples_allele_fraction as usize,
            num_burn_in_allele_fraction: options.number_of_burn_in_samples_allele_fraction as usize,
        },
    )?;
    write_fit(&path, metadata, &mut modeller, ".modelBegin")?;
    modeller.smooth_segments(
        options.maximum_number_of_smoothing_iterations,
        options.number_of_smoothing_iterations_per_fit,
        options.smoothing_credible_interval_threshold_copy_ratio,
        options.smoothing_credible_interval_threshold_allele_fraction,
    )?;
    write_fit(&path, metadata, &mut modeller, ".modelFinal")?;

    let modeled = modeller.modeled_segments.clone();
    write(
        &path(".cr.seg"),
        &copy_ratio_segments_text(metadata, &modeled, &denoised.records)?,
    )?;
    let copy_ratio_legacy: Vec<(Span, i32, f64)> = modeled
        .iter()
        .map(|s| (s.span.clone(), s.num_points_copy_ratio, s.log2_copy_ratio.decile50))
        .collect();
    let allele_fraction_legacy: Vec<(Span, i32, f64)> = modeled
        .iter()
        .map(|s| {
            (
                s.span.clone(),
                s.num_points_allele_fraction,
                s.minor_allele_fraction.decile50,
            )
        })
        .collect();
    write(&path(".cr.igv.seg"), &legacy_text(&metadata.sample, &copy_ratio_legacy))?;
    write(
        &path(".af.igv.seg"),
        &legacy_text(&metadata.sample, &allele_fraction_legacy),
    )?;
    Ok(())
}

/// `writeModeledSegmentsAndParameterFiles`: the segments, then the two parameter files.
fn write_fit(
    path: &dyn Fn(&str) -> String,
    metadata: &Metadata,
    modeller: &mut MultidimensionalModeller,
    tag: &str,
) -> Result<(), Thrown> {
    write(
        &path(&format!("{tag}.seg")),
        &modeled_segments_text(metadata, &modeller.modeled_segments),
    )?;
    let (copy_ratio, allele_fraction) = modeller.parameter_deciles()?;
    write(
        &path(&format!("{tag}.cr.param")),
        &parameters_text(&metadata.sample, &copy_ratio),
    )?;
    write(
        &path(&format!("{tag}.af.param")),
        &parameters_text(&metadata.sample, &allele_fraction),
    )?;
    Ok(())
}

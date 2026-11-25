//! `FilterVcf.doWork`.

use gatk_barclay::Parser;
use gatk_tools::filter_vcf::{filter_to, Thresholds};
use gatk_tools::main_entry::{Thrown, PORT_LIMITATION};
use htsjdk_vcf::reader::{read_vcf, ReadError, ReadFailure};

use super::vcf_format_converter::output_type;
use super::{assert_readable, assert_writable, thrown, vcf_text, write_vcf};
use crate::runners::{argument, flag, scalar};

pub(super) fn run(parser: &Parser, _stderr: &mut String) -> Result<i32, Thrown> {
    let input = argument(parser, "INPUT").unwrap_or_default();
    let output = argument(parser, "OUTPUT").unwrap_or_default();
    assert_readable(&input)?;
    assert_writable(&output)?;
    let defaults = Thresholds::default();
    let number = |name: &str, default: f64| {
        scalar(parser, name)
            .and_then(|value| value.parse::<f64>().ok())
            .unwrap_or(default)
    };
    let integer = |name: &str, default: i32| {
        scalar(parser, name)
            .and_then(|value| value.parse::<i32>().ok())
            .unwrap_or(default)
    };
    let thresholds = Thresholds {
        min_ab: number("MIN_AB", defaults.min_ab),
        max_fs: number("MAX_FS", defaults.max_fs),
        min_qd: number("MIN_QD", defaults.min_qd),
        min_gq: integer("MIN_GQ", defaults.min_gq),
        min_dp: integer("MIN_DP", defaults.min_dp),
    };
    let file_name = output.rsplit('/').next().unwrap_or_default();
    let vcf_or_bcf = file_name.ends_with(".vcf") || file_name.ends_with(".bcf");
    let text = vcf_text(&input)?;
    // The script filter is built once `new VCFFileReader` has read the header, before the first
    // record is read and before the dictionary check: only a header failure comes first.
    if argument(parser, "JAVASCRIPT_FILE").is_some_and(|file| file != "null")
        && !matches!(
            read_vcf(&text),
            Err(ReadFailure {
                error: ReadError::Header(_) | ReadError::HeaderLine(_),
                ..
            })
        )
    {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "FilterVcf's JAVASCRIPT_FILE runs a Nashorn script, which this port does not. This \
             message is the port's own and not GATK's.",
        ));
    }
    let filtered = filter_to(&text, &thresholds, vcf_or_bcf)
        .map_err(|error| thrown(error.java_class(), error.message()))?;
    output_type(&output)?;
    // Picard 3.5.0 sets or unsets `INDEX_ON_THE_FLY` from `CREATE_INDEX` (b976447cf); 3.4.0 kept
    // the builder's default and indexed whatever it said.
    if crate::bcf::is_bcf(&input) && crate::bcf::is_bcf(&output) {
        write_bcf(&input, &output, &filtered, flag(parser, "CREATE_INDEX"))?;
    } else {
        write_vcf(parser, &output, &filtered, flag(parser, "CREATE_INDEX"))?;
    }
    Ok(0)
}

/// A BCF filtered into a BCF, as `FilterApplyingVariantIterator.next` builds each record: the BCF
/// codec's own, with its INFO order and its 32-bit floats, under the filters the text gives it.
///
/// The text does not say how Java holds the site's filters, and a BCF codec record writes them in
/// the order it holds them: a `HashSet`, or `passFilters()`'s empty set. A genotype's `PASS` is no
/// filter to `BCF2Writer` as it is to the VCF codec, so the text's genotype filters are the
/// record's. A VCF input is the VCF codec's records, which the text path writes as they are.
fn write_bcf(input: &str, output: &str, filtered: &str, index: bool) -> Result<(), Thrown> {
    use htsjdk_vcf::bcf2::java_hash::{order_through, Fill};
    use htsjdk_vcf::bcf2::writer::Origin;

    let file = read_vcf(filtered)
        .map_err(|failure| thrown(failure.error.class(), failure.error.message()))?;
    // FilterVcf keeps every record, in order, so the n-th filtered line is the n-th record read.
    let sources = crate::bcf::read(input)?.records;
    let mut records = Vec::with_capacity(sources.len());
    for (mut record, out) in sources.into_iter().zip(&file.records) {
        record.filters = out.filters.as_ref().map(|filters| {
            let order = order_through(filters, &[Fill::Put]);
            order.into_iter().map(|at| filters[at].clone()).collect()
        });
        let genotypes = record
            .genotypes
            .try_get_mut()
            .map_err(|error| thrown(error.class, error.message))?;
        for genotype in genotypes {
            genotype.filters = out
                .genotype(&genotype.sample_name)
                .and_then(|filtered| filtered.filters.clone());
        }
        records.push((record, Origin::Bcf2Codec));
    }
    let dictionary: Vec<(String, i32)> = file
        .header
        .lines
        .iter()
        .filter_map(|line| match line {
            htsjdk_vcf::header::HeaderLine::Contig { fields, .. } => {
                let field = |key: &str| fields.iter().find(|(k, _)| k == key).map(|(_, v)| v);
                Some((
                    field("ID")?.clone(),
                    field("length").and_then(|v| v.parse().ok()).unwrap_or(0),
                ))
            }
            _ => None,
        })
        .collect();
    if index && dictionary.is_empty() {
        return Err(thrown(
            "java.lang.IllegalArgumentException",
            "A reference dictionary is required for creating Tribble indices on the fly"
                .to_string(),
        ));
    }
    crate::bcf::write(
        output,
        &file.header,
        &records,
        index.then_some(dictionary.as_slice()),
        false,
    )
}

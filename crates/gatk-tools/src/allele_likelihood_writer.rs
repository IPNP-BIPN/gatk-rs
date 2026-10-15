//! `--alm-path` and `--alm-interval`, ported from `tools.walkers.haplotypecaller.AlleleLikelihoodWriter`
//! (GATK 4.6.2.0).
//!
//! Each region's read-by-haplotype matrix, as HaplotypeCaller computes it (the ramped engine
//! before its uncollapsing, HaplotypeCaller after), when the interval contains the first
//! haplotype: the haplotypes, then per sample the reads, the matrix row by row, and every
//! read-haplotype pair again with its likelihood. Values are `Double.toString`. The file is opened
//! when the engine starts and never closed, only flushed after each region.

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::haplotype::Haplotype;
use gatk_engine::tsv_table::java_double_to_string;
use htsjdk_bam::record::BamRecord;

/// `new SimpleInterval(String)`: `contig`, `contig:start`, `contig:start+` or
/// `contig:start-end`, commas allowed in the numbers.
pub fn parse_simple_interval(text: Option<&str>) -> Result<(String, i32, i32), String> {
    let text = text.ok_or_else(|| {
        "java.lang.IllegalArgumentException: Null object is not allowed here.".to_string()
    })?;
    if text.is_empty() {
        return Err("java.lang.IllegalArgumentException: str should not be empty".to_string());
    }
    let position = |p: &str| -> Result<i32, String> {
        p.replace(',', "")
            .parse::<i32>()
            .map_err(|_| format!("java.lang.NumberFormatException: For input string: \"{p}\""))
    };
    let Some(colon) = text.rfind(':') else {
        return Ok((text.to_string(), 1, i32::MAX));
    };
    let contig = text[..colon].to_string();
    let rest = &text[colon + 1..];
    match rest.find('-') {
        None if rest.ends_with('+') => Ok((contig, position(&rest[..rest.len() - 1])?, i32::MAX)),
        None => {
            let start = position(rest)?;
            Ok((contig, start, start))
        }
        Some(dash) => Ok((
            contig,
            position(&rest[..dash])?,
            position(&rest[dash + 1..])?,
        )),
    }
}

/// The writer: its interval and what it has written.
#[derive(Debug, Clone)]
pub struct AlleleLikelihoodWriter {
    interval: (String, i32, i32),
    pub text: String,
}

impl AlleleLikelihoodWriter {
    pub fn new(interval: (String, i32, i32)) -> AlleleLikelihoodWriter {
        AlleleLikelihoodWriter {
            interval,
            text: String::new(),
        }
    }

    /// `writeAlleleLikelihoods(likelihoods)`.
    pub fn write(&mut self, likelihoods: &AlleleLikelihoods<BamRecord, Haplotype>) {
        let haplotypes: Vec<&Haplotype> = (0..likelihoods.number_of_alleles())
            .filter_map(|a| likelihoods.get_allele(a))
            .collect();
        let Some(first) = haplotypes.first().and_then(|h| h.genome_location()) else {
            return;
        };
        let (contig, start, end) = &self.interval;
        if !(first.contig == *contig && *start <= first.start && first.end <= *end) {
            return;
        }
        let out = &mut self.text;
        out.push_str(&format!(
            "> Location {}:{}-{}\n>> Haplotypes\n",
            first.contig, first.start, first.end
        ));
        let names: Vec<String> = haplotypes
            .iter()
            .map(|h| crate::pd_haplotypes::haplotype_string(h))
            .collect();
        for (i, name) in names.iter().enumerate() {
            out.push_str(&format!("{i:04}\t{name}\n"));
        }
        for s in 0..likelihoods.number_of_samples() {
            out.push_str(&format!(
                ">> Sample {}\n>>> Reads\n",
                likelihoods.get_sample(s).map_or("", |n| n.as_str())
            ));
            let reads = likelihoods.sample_evidence(s).unwrap_or(&[]);
            for (i, read) in reads.iter().enumerate() {
                out.push_str(&format!("{i:04}\t{}\n", read.read_name));
            }
            out.push_str(">>> Matrix\n");
            for a in 0..haplotypes.len() {
                for r in 0..reads.len() {
                    out.push_str(&java_double_to_string(likelihoods.value(s, a, r)));
                    out.push(' ');
                }
                out.push('\n');
            }
            out.push_str(">>> Read->Haplotype in Full\n");
            for (a, name) in names.iter().enumerate() {
                for (r, read) in reads.iter().enumerate() {
                    out.push_str(&format!(
                        "{r:04}\t{}\t{}\t{a:04}\t{name}\n",
                        read.read_name,
                        java_double_to_string(likelihoods.value(s, a, r))
                    ));
                }
                out.push('\n');
            }
        }
    }
}

//! `--alm-path` and `--alm-interval`, ported from
//! `tools.walkers.haplotypecaller.AlleleLikelihoodWriter` (GATK 4.7.0.0;
//! `AlleleLikelihoodWriter.java` and `SimpleInterval.java` are unchanged since 4.6.2.0).
//!
//! Each region's read-by-haplotype matrix, as HaplotypeCaller computes it (the ramped engine
//! before its uncollapsing, HaplotypeCaller after), when the interval contains the first
//! haplotype: the haplotypes, then per sample the reads, the matrix row by row, and every
//! read-haplotype pair again with its likelihood. Values are `Double.toString`. The file is opened
//! when the engine starts and never closed, only flushed after each region.
//!
//! The Java tests `outputInterval == null || outputInterval.contains(first_hap)`, but
//! `createAlleleLikelihoodWriter` always builds the interval with `new SimpleInterval(..)`, which
//! throws on a null string. So the null branch cannot be reached from the command line, and
//! `--alm-path` without `--alm-interval` fails at start-up: the golden
//! `crates/gatk-tools/tests/data/hc_debug_outputs.txt.gz` has the row `error hc-alm-no-interval
//! java.lang.IllegalArgumentException:Null object is not allowed here.`. That is why the writer
//! here holds a plain interval and the parser takes an `Option`.
//!
//! # Reading this file without reading Rust
//!
//! - `Result<T, String>` is "a `T` or an error message"; `Err(..)` is the throw, `Ok(..)` the
//!   return, and `?` after a call passes its error straight up to the caller, Java's `throws`.
//! - `let Some(x) = y else { return ..; };` unwraps `y` or takes the early exit.
//! - `format!("{i:04}")` is Java's `String.format("%04d", i)`: zero-padded to four digits.
//! - `&mut self` is a method that may change the writer; `&self` one that only reads it.

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::haplotype::Haplotype;
use gatk_engine::tsv_table::java_double_to_string;
use htsjdk_bam::record::BamRecord;

/// `new SimpleInterval(String)`: `contig`, `contig:start`, `contig:start+` or
/// `contig:start-end`, commas allowed in the numbers.
///
/// The split is on the LAST colon and then the first dash after it, as the Java does:
///
/// ```java
/// final int colonIndex = str.lastIndexOf(CONTIG_SEPARATOR);
/// ...
/// final int dashIndex = str.indexOf(START_END_SEPARATOR, colonIndex);
/// ```
///
/// so a contig name may itself contain colons. A name with no colon at all is the whole contig,
/// `1` to `Integer.MAX_VALUE`.
///
/// Two parts of the constructor are not ported here. `parsePosition` wraps a bad number in a
/// `UserException` ("Problem parsing start/end value in interval string. Value was: ..."), where
/// this returns the bare `NumberFormatException` text; and `validatePositions` (start at least 1,
/// end at least start) is not applied.
pub fn parse_simple_interval(text: Option<&str>) -> Result<(String, i32, i32), String> {
    // `Utils.nonNull(str)` and `Utils.validateArg(!str.isEmpty(), "str should not be empty")`.
    // `ok_or_else` turns `None` into the error, and `?` returns it.
    let text = text.ok_or_else(|| {
        "java.lang.IllegalArgumentException: Null object is not allowed here.".to_string()
    })?;
    if text.is_empty() {
        return Err("java.lang.IllegalArgumentException: str should not be empty".to_string());
    }
    // `parsePositionThrowOnFailure`: `Integer.parseInt(pos.replaceAll(",", ""))`. A closure,
    // so the three call sites below share it.
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
    // `match` must cover every case. The first arm has a guard (`if ..`): no dash and a trailing
    // `+` is `start` to the end of the contig. `&rest[..rest.len() - 1]` drops that `+`.
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

/// `AlleleLikelihoodWriter`: its interval and what it has written.
///
/// The Java holds a `FileWriter`; this holds the text, which the caller writes out. `text` is
/// `pub` so the caller can read it.
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
    ///
    /// The layout, line by line:
    ///
    /// ```text
    /// > Location <contig>:<start>-<end>        (of the first haplotype)
    /// >> Haplotypes
    /// <i, %04d>\t<haplotype>                   (one per haplotype)
    /// >> Sample <name>                         (then, per sample:)
    /// >>> Reads
    /// <i, %04d>\t<read name>                   (one per read)
    /// >>> Matrix
    /// <v> <v> ... <v> \n                       (one row per haplotype, a space after EVERY value)
    /// >>> Read->Haplotype in Full
    /// <r, %04d>\t<read>\t<v>\t<a, %04d>\t<haplotype>   (one per read, then an empty line,
    ///                                                  per haplotype)
    /// ```
    ///
    /// Every value is `Double.toString`, so `-0.5`, `-1.0E-5`, `-Infinity`. The haplotype is its
    /// `toString()`, which for a plain `Haplotype` is `getDisplayString()`; a partially determined
    /// one has its own (see `crate::pd_haplotypes::haplotype_string`).
    pub fn write(&mut self, likelihoods: &AlleleLikelihoods<BamRecord, Haplotype>) {
        // `likelihoods.alleles()`. `filter_map` keeps the `Some` results of `get_allele`, which is
        // every index in range.
        let haplotypes: Vec<&Haplotype> = (0..likelihoods.number_of_alleles())
            .filter_map(|a| likelihoods.get_allele(a))
            .collect();
        // `haplotypes.get(0)` and its location. The Java would throw on an empty list or a
        // haplotype without a location; this writes nothing.
        let Some(first) = haplotypes.first().and_then(|h| h.genome_location()) else {
            return;
        };
        // `outputInterval.contains(first_hap)`: same contig and the haplotype wholly inside,
        // htsjdk's `CoordMath.encloses`. Anything else writes nothing for this region.
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
            // `&[]` is an empty borrowed list, standing in for a sample with no evidence.
            let reads = likelihoods.sample_evidence(s).unwrap_or(&[]);
            for (i, read) in reads.iter().enumerate() {
                out.push_str(&format!("{i:04}\t{}\n", read.read_name));
            }
            // Haplotypes down, reads across: `sampleMatrix(s).get(allele, read)`, each value
            // followed by `output.write(" ")`, so every row ends in a space before its newline.
            out.push_str(">>> Matrix\n");
            for a in 0..haplotypes.len() {
                for r in 0..reads.len() {
                    out.push_str(&java_double_to_string(likelihoods.value(s, a, r)));
                    out.push(' ');
                }
                out.push('\n');
            }
            // The same values again, one line each, haplotype-major, with the Java's
            // `output.write("\n")` after each haplotype's reads leaving an empty line.
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

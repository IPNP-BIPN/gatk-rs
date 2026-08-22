//! `SiteDepthtoBAF`, ported from `org.broadinstitute.hellbender.tools.sv.SiteDepthtoBAF`,
//! `BafEvidence` and `BafEvidenceCodec` (GATK 4.6.2.0).
//!
//! Per-sample allele depths at a set of sites turned into B-allele fractions.
//!
//! # Which of the four counts is "ref" is decided by another file
//!
//! The depth record carries four counts, A C G T. Which one is the reference and which the
//! alternate comes from the sites VCF walked in lockstep beside it, so the same depths produce
//! different fractions under different sites files. The two must agree on the locus, and the VCF
//! must not run out first; either is a `UserException` naming both positions.
//!
//! # The value written is rarely the fraction measured
//!
//! ```java
//! if ( nBafs <= 1 ) {
//!     if ( nBafs == 1 ) { writer.write(new BafEvidence(beList.get(0), .5)); }
//!     return;
//! }
//! ```
//!
//! A locus with exactly one surviving sample is written as `0.5` **whatever it measured**: the
//! value is replaced, not adjusted. A locus with two or more is shifted so its median lands on
//! 0.5, so every value moves. Only the ordering within a locus survives untouched.
//!
//! # The deviation is the sample deviation, and it is Welford's
//!
//! `MathUtils.RunningAverage.stddev()` is `sqrt(s / (n - 1))`, accumulated by Welford's
//! incremental update rather than by a second pass. Both matter: `n - 1` rather than `n`, and the
//! order the values arrive in, which decides the last bits.
//!
//! # The chi-squared test
//!
//! One degree of freedom, against an expectation of half the total, computed from the ref and alt
//! counts alone while the total that feeds the expectation is the sum of **all four**. So a site
//! whose other two counts are large has a high expectation and fails a test its two named alleles
//! would have passed.

use jmath::gamma::regularized_gamma_p;

/// One row of a `.sd.txt` file, one-based as the tool holds it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SiteDepth {
    pub contig: String,
    pub position: i32,
    pub sample: String,
    /// A, C, G, T.
    pub counts: [i32; 4],
}

impl SiteDepth {
    /// `getTotalDepth()`: all four, not the two the alleles name.
    pub fn total_depth(&self) -> i32 {
        self.counts.iter().sum()
    }
}

/// One row of a `.baf.txt` file.
#[derive(Debug, Clone, PartialEq)]
pub struct BafEvidence {
    pub sample: String,
    pub contig: String,
    pub position: i32,
    pub value: f64,
}

/// One line of the sites VCF, reduced to what the tool reads from it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Site {
    pub contig: String,
    pub position: i32,
    pub reference: Vec<u8>,
    pub alternate: Vec<u8>,
}

impl Site {
    /// `BAFSiteIterator.advance`: a SNP, biallelic, and at a locus after the last one taken.
    fn usable(&self, last: Option<&Site>) -> bool {
        let snp = self.reference.len() == 1 && self.alternate.len() == 1;
        match last {
            _ if !snp => false,
            None => true,
            Some(last) => last.contig != self.contig || last.position < self.position,
        }
    }
}

/// The arguments, with the tool's own defaults.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Arguments {
    pub max_std_dev: f64,
    pub min_total_depth: i32,
    pub min_het_probability: f64,
}

impl Default for Arguments {
    fn default() -> Self {
        Arguments {
            max_std_dev: 0.2,
            min_total_depth: 10,
            min_het_probability: 0.5,
        }
    }
}

/// What the run refuses, all of it while walking.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BafError {
    /// The sites VCF ran out while depth records remained.
    SitesExhausted,
    /// The sites VCF and the depths disagree about where they are.
    LocusMismatch {
        expected_contig: String,
        expected_position: i32,
        found_contig: String,
        found_position: i32,
    },
    /// The reference base is not one of A, C, G, T.
    ReferenceNotAcgt { contig: String, position: i32 },
    /// The alternate base is not one of A, C, G, T.
    AlternateNotAcgt { contig: String, position: i32 },
}

impl BafError {
    pub fn java_class(&self) -> &str {
        "org.broadinstitute.hellbender.exceptions.UserException"
    }

    pub fn message(&self) -> String {
        match self {
            BafError::SitesExhausted => {
                "baf sites vcf exhausted before site depth data".to_string()
            }
            BafError::LocusMismatch {
                expected_contig,
                expected_position,
                found_contig,
                found_position,
            } => format!(
                "expecting locus {expected_contig}:{expected_position}, but found locus \
                 {found_contig}:{found_position} in baf sites vcf"
            ),
            BafError::ReferenceNotAcgt { contig, position } => {
                format!("ref call is not [ACGT] in vcf at {contig}:{position}")
            }
            BafError::AlternateNotAcgt { contig, position } => {
                format!("alt call is not [ACGT] in vcf at {contig}:{position}")
            }
        }
    }
}

/// `Nucleotide.decode(base).ordinal()`, for the four the tool accepts. Anything else answers a
/// value above three, which is how both refusals are phrased.
pub fn nucleotide_ordinal(base: u8) -> usize {
    match base.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

/// `new ChiSquaredDistribution(1.).cumulativeProbability(x)`.
///
/// Apache's chi-squared with one degree of freedom is a gamma with shape `1/2` and scale `2`, and
/// a gamma's CDF is `Gamma.regularizedGammaP(shape, x / scale)` at the default epsilon of `1e-14`
/// and an iteration cap of `Integer.MAX_VALUE`.
pub fn chi_squared_cdf_one_dof(x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    regularized_gamma_p(0.5, x / 2.0, 1e-14, i32::MAX).unwrap_or(f64::NAN)
}

/// `calcBAF`: the depth floor, then the chi-squared fit, then the alt fraction.
pub fn calc_baf(
    depth: &SiteDepth,
    reference_index: usize,
    alternate_index: usize,
    arguments: &Arguments,
) -> Option<BafEvidence> {
    let total_depth = depth.total_depth();
    if total_depth < arguments.min_total_depth {
        return None;
    }
    let expect_ref_alt = f64::from(total_depth) / 2.0;
    let alt_depth = f64::from(depth.counts[alternate_index]);
    let ref_diff = f64::from(depth.counts[reference_index]) - expect_ref_alt;
    let alt_diff = alt_depth - expect_ref_alt;
    let chi_sq = (ref_diff * ref_diff + alt_diff * alt_diff) / expect_ref_alt;
    let fit_probability = 1.0 - chi_squared_cdf_one_dof(chi_sq);
    if fit_probability < arguments.min_het_probability {
        return None;
    }
    Some(BafEvidence {
        sample: depth.sample.clone(),
        contig: depth.contig.clone(),
        position: depth.position,
        value: alt_depth / f64::from(total_depth),
    })
}

/// `MathUtils.RunningAverage.stddev()`: Welford's incremental variance, over `n - 1`.
///
/// Written as the accumulation rather than as a two-pass formula because the two do not agree in
/// the last bits, and the order the values arrive in is part of the answer.
pub fn running_stddev(values: &[f64]) -> f64 {
    let mut mean = 0.0;
    let mut s = 0.0;
    let mut count = 0i64;
    for value in values {
        count += 1;
        let old_mean = mean;
        mean += (value - mean) / count as f64;
        s += (value - old_mean) * (value - mean);
    }
    (s / (count - 1) as f64).sqrt()
}

/// The median of the sorted values, the even case being the mean of the two in the middle.
fn median(values: &mut [f64]) -> f64 {
    values.sort_by(|left, right| left.partial_cmp(right).expect("a comparable baf"));
    let middle = values.len() / 2;
    if !values.len().is_multiple_of(2) {
        values[middle]
    } else {
        (values[middle] + values[middle - 1]) / 2.0
    }
}

/// `processBuffer`: one locus of depth records against one site, and what it writes.
pub fn process_locus(
    buffer: &[SiteDepth],
    site: &Site,
    arguments: &Arguments,
) -> Result<Vec<BafEvidence>, BafError> {
    let held = buffer.first().expect("a non-empty locus");
    if held.contig != site.contig || held.position != site.position {
        return Err(BafError::LocusMismatch {
            expected_contig: held.contig.clone(),
            expected_position: held.position,
            found_contig: site.contig.clone(),
            found_position: site.position,
        });
    }
    let reference_index = nucleotide_ordinal(site.reference[0]);
    if reference_index > 3 {
        return Err(BafError::ReferenceNotAcgt {
            contig: site.contig.clone(),
            position: site.position,
        });
    }
    let alternate_index = nucleotide_ordinal(site.alternate[0]);
    if alternate_index > 3 {
        return Err(BafError::AlternateNotAcgt {
            contig: site.contig.clone(),
            position: site.position,
        });
    }
    let kept: Vec<BafEvidence> = buffer
        .iter()
        .filter_map(|depth| calc_baf(depth, reference_index, alternate_index, arguments))
        .collect();
    if kept.len() <= 1 {
        // One survivor is written as a half whatever it measured, and none is written at all.
        return Ok(kept
            .into_iter()
            .map(|evidence| BafEvidence {
                value: 0.5,
                ..evidence
            })
            .collect());
    }
    let values: Vec<f64> = kept.iter().map(|evidence| evidence.value).collect();
    if running_stddev(&values) > arguments.max_std_dev {
        // The whole locus goes, not the outlying sample, and nothing in the file says so.
        return Ok(Vec::new());
    }
    let adjustment = 0.5 - median(&mut values.clone());
    Ok(kept
        .into_iter()
        .map(|evidence| BafEvidence {
            value: evidence.value + adjustment,
            ..evidence
        })
        .collect())
}

/// The whole traversal: depth records grouped by locus, each group against the next usable site.
pub fn run(
    depths: &[SiteDepth],
    sites: &[Site],
    arguments: &Arguments,
) -> Result<Vec<BafEvidence>, BafError> {
    let mut written = Vec::new();
    let mut buffer: Vec<SiteDepth> = Vec::new();
    let mut remaining = sites.iter();
    let mut last: Option<Site> = None;
    let mut next_site = |last: &mut Option<Site>| -> Result<Site, BafError> {
        for site in remaining.by_ref() {
            if site.usable(last.as_ref()) {
                *last = Some(site.clone());
                return Ok(site.clone());
            }
        }
        Err(BafError::SitesExhausted)
    };
    for depth in depths {
        let same = buffer
            .first()
            .is_none_or(|held| held.position == depth.position && held.contig == depth.contig);
        if !same {
            written.extend(process_locus(&buffer, &next_site(&mut last)?, arguments)?);
            buffer.clear();
        }
        buffer.push(depth.clone());
    }
    if !buffer.is_empty() {
        written.extend(process_locus(&buffer, &next_site(&mut last)?, arguments)?);
    }
    Ok(written)
}

/// `new DecimalFormat("#.00").format(value)`.
///
/// Two decimals always, and no integer digit when the integer part is zero, so a half comes out
/// `.50` rather than `0.50`. Ties go to even, which is `DecimalFormat`'s default and not
/// `String.format`'s.
pub fn format_baf(value: f64) -> String {
    let scaled = value * 100.0;
    let floor = scaled.floor();
    let rounded = match (scaled - floor).partial_cmp(&0.5) {
        Some(std::cmp::Ordering::Greater) => floor + 1.0,
        Some(std::cmp::Ordering::Less) => floor,
        _ if (floor as i64) % 2 == 0 => floor,
        _ => floor + 1.0,
    } / 100.0;
    let text = format!("{rounded:.2}");
    // `#` before the point means the zero is not printed at all.
    if let Some(stripped) = text.strip_prefix("0.") {
        return format!(".{stripped}");
    }
    if let Some(stripped) = text.strip_prefix("-0.") {
        return format!("-.{stripped}");
    }
    text
}

/// `BafEvidenceCodec.encode` over a whole file: zero-based on disk, one-based inside.
pub fn write(records: &[BafEvidence]) -> String {
    let mut out = String::new();
    for record in records {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\n",
            record.contig,
            record.position - 1,
            format_baf(record.value),
            record.sample
        ));
    }
    out
}

/// `SiteDepthCodec.decode` over a whole file.
pub fn read_depths(text: &str) -> Vec<SiteDepth> {
    text.lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            SiteDepth {
                contig: columns[0].to_string(),
                // Zero-based on disk, one-based here.
                position: columns[1].parse::<i32>().expect("a start") + 1,
                sample: columns[2].to_string(),
                counts: [
                    columns[3].parse().expect("an A count"),
                    columns[4].parse().expect("a C count"),
                    columns[5].parse().expect("a G count"),
                    columns[6].parse().expect("a T count"),
                ],
            }
        })
        .collect()
}

/// The sites VCF, reduced to the four fields the tool reads.
pub fn read_sites(text: &str) -> Vec<Site> {
    text.lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            Site {
                contig: columns[0].to_string(),
                position: columns[1].parse().expect("a position"),
                reference: columns[3].as_bytes().to_vec(),
                alternate: columns[4].as_bytes().to_vec(),
            }
        })
        .collect()
}

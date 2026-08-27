//! `GnarlyGenotyper`: the engine's own rules, which decide whether a combined site is worth
//! calling and what its annotations become.
//!
//! The genotyping arithmetic that rewrites the PLs is not ported. What is ported is the decision
//! in front of it: the quality floor, which floor applies, where QUALapprox comes from, and the
//! two numbers a called site carries.
//!
//! Ported from
//! `org.broadinstitute.hellbender.tools.walkers.gnarlyGenotyper.GnarlyGenotyperEngine` and
//! `org.broadinstitute.hellbender.utils.variant.HomeSapiensConstants` in GATK 4.6.2.0.

/// `HomoSapiensConstants.SNP_HETEROZYGOSITY` and its indel companion.
pub const SNP_HETEROZYGOSITY: f64 = 1e-3;
pub const INDEL_HETEROZYGOSITY: f64 = 1.25e-4;

/// `DEFAULT_STANDARD_CONFIDENCE_FOR_CALLING`, which is NOT the floor.
pub const DEFAULT_STANDARD_CONFIDENCE_FOR_CALLING: f64 = 30.0;

/// The quality floor a site has to reach.
///
/// It is the confidence argument LESS ten times the logarithm of the site prior, so the SNP floor
/// is exactly 60 and the indel floor about 69.03. A reader who takes the argument for the floor
/// reads every site under 30 as called and every site between 30 and 60 as called too, which is
/// what makes the confusion silent.
pub fn quality_floor(prior: f64) -> f64 {
    DEFAULT_STANDARD_CONFIDENCE_FOR_CALLING - 10.0 * prior.log10()
}

pub fn snp_quality_floor() -> f64 {
    quality_floor(SNP_HETEROZYGOSITY)
}

pub fn indel_quality_floor() -> f64 {
    quality_floor(INDEL_HETEROZYGOSITY)
}

/// The spanning-deletion allele, which does not count as a SNP however it is written.
pub const SPAN_DEL: &str = "*";
/// The symbolic allele a combined GVCF carries, which a called site loses.
pub const NON_REF: &str = "<NON_REF>";

/// Whether any alternate is the reference's own length, which is what decides the floor.
///
/// The spanning deletion is excluded by name before the length is looked at, so a site whose only
/// alternate is `*` is judged as an indel however long the reference is.
pub fn has_snp_allele(reference: &str, alternates: &[String]) -> bool {
    alternates
        .iter()
        .any(|allele| allele != SPAN_DEL && allele.len() == reference.len())
}

/// The prior a site is judged against.
pub fn site_prior(reference: &str, alternates: &[String]) -> f64 {
    if has_snp_allele(reference, alternates) {
        SNP_HETEROZYGOSITY
    } else {
        INDEL_HETEROZYGOSITY
    }
}

/// `RAW_QUAL_APPROX_KEY` and its allele-specific companion.
pub const QUAL_APPROX_KEY: &str = "QUALapprox";
pub const AS_QUAL_APPROX_KEY: &str = "AS_QUALapprox";

/// Where QUALapprox comes from.
///
/// The plain key is preferred, the allele-specific list is SUMMED when there is no plain key, and
/// a site with neither scores zero, which is under every floor.
pub fn qual_approx(plain: Option<i32>, allele_specific: Option<&str>) -> f64 {
    if let Some(plain) = plain {
        return plain as f64;
    }
    if let Some(list) = allele_specific {
        return parse_qual_list(list).into_iter().sum::<i32>() as f64;
    }
    0.0
}

/// `AS_QualByDepth.parseQualList`, which splits on `|` and skips the empty leading field.
pub fn parse_qual_list(list: &str) -> Vec<i32> {
    list.split('|')
        .filter(|part| !part.is_empty())
        .filter_map(|part| part.parse().ok())
        .collect()
}

/// Whether a site clears its floor.
pub fn clears_the_floor(qual_approx: f64, reference: &str, alternates: &[String]) -> bool {
    let is_indel = !has_snp_allele(reference, alternates);
    if is_indel {
        qual_approx >= indel_quality_floor()
    } else {
        qual_approx >= snp_quality_floor()
    }
}

/// The filter a site under its floor is given when `--keep-all-sites` asks for it back.
pub const LOW_QUAL_FILTER_NAME: &str = "LowQual";
/// The key that site carries instead of its allele counts.
pub const AC_ADJUSTED_KEY: &str = "AC_adj";

/// What the engine does with a site.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Outcome {
    /// Under its floor and dropped: the engine returns nothing at all.
    Dropped,
    /// Under its floor and kept, filtered, with an adjusted count of zero.
    LowQual,
    /// Called.
    Called,
}

pub fn outcome(
    qual_approx: f64,
    reference: &str,
    alternates: &[String],
    keep_all_sites: bool,
) -> Outcome {
    if clears_the_floor(qual_approx, reference, alternates) {
        return Outcome::Called;
    }
    if keep_all_sites {
        Outcome::LowQual
    } else {
        Outcome::Dropped
    }
}

/// `QD`: QUALapprox over the variant depth.
pub fn quality_by_depth(qual_approx: f64, variant_depth: i32) -> f64 {
    qual_approx / variant_depth as f64
}

/// The site's own QUAL, in Phred.
///
/// The engine sets `log10PError` to `QUALapprox / -10 - log10(prior)`, so the Phred-scaled quality
/// the writer prints is ten times the negation of that: `QUALapprox + 10 * log10(prior)`. For a
/// QUALapprox of 900 at the SNP prior that is 870, which is neither QUALapprox nor QD.
pub fn phred_quality(qual_approx: f64, prior: f64) -> f64 {
    qual_approx + 10.0 * prior.log10()
}

/// The alternates a called site keeps: everything but `<NON_REF>`.
pub fn called_alternates(alternates: &[String]) -> Vec<String> {
    alternates
        .iter()
        .filter(|allele| *allele != NON_REF)
        .cloned()
        .collect()
}

/// The annotations a `LowQual` site carries, which is not what a called one carries.
///
/// It keeps QUALapprox, VarDP and the finalized MQ, gains `AC_adj=0`, and has none of AC, AF, AN,
/// QD, FS, SOR or ExcessHet. Its alternates keep their `<NON_REF>`.
pub fn low_qual_keys() -> Vec<&'static str> {
    vec![AC_ADJUSTED_KEY, "MQ", QUAL_APPROX_KEY, "VarDP"]
}

/// The keys a called site gains.
pub fn called_keys() -> Vec<&'static str> {
    vec![
        "AC",
        "AF",
        "AN",
        "ExcessHet",
        "FS",
        "MQ",
        "QD",
        "SOR",
        "VarDP",
    ]
}

/// `RAW_MQandDP`, which is a sum of squared mapping qualities and a depth.
pub const RAW_MQ_AND_DP_KEY: &str = "RAW_MQandDP";

/// `MQ`, finalized from the raw sum and its depth.
///
/// A depth of zero yields nothing, and the key is then absent from the output rather than written
/// as zero or as a NaN, which is the same absence a site with no raw key at all shows.
pub fn finalize_mq(raw: &str) -> Option<f64> {
    let (sum, depth) = raw.split_once(',')?;
    let sum: f64 = sum.trim().parse().ok()?;
    let depth: f64 = depth.trim().parse().ok()?;
    if depth == 0.0 {
        return None;
    }
    Some((sum / depth).sqrt())
}

/// The allele-specific keys `--strip-allele-specific-annotations` acts on.
///
/// It does NOT remove them: the finalizer writes each of them as `null`, so the output carries
/// `AS_QD=null` whether the argument was given or not, and the argument is told apart by which
/// of the other AS keys survive.
pub const ALLELE_SPECIFIC_PREFIX: &str = "AS_";

pub fn is_allele_specific(key: &str) -> bool {
    key.starts_with(ALLELE_SPECIFIC_PREFIX)
}

/// The refusal an allele-specific list of the wrong length produces.
pub const ALLELE_COUNT_MISMATCH_MESSAGE: &str =
    "Number of alleles and number of allele-specific entries do not match.";

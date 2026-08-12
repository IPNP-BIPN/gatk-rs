//! `BAQ`, ported from `org.broadinstitute.hellbender.utils.baq` (GATK 4.6.2.0).
//!
//! Base Alignment Quality: a hidden Markov model that caps each base's quality by how confidently
//! the aligner could have placed it. `BaseRecalibrator` runs it over every read before counting, so
//! the recalibration table rests on these bytes.
//!
//! The reference's own banner says the kernel is kept in sync with samtools and that changes belong
//! upstream, so the transcription below is deliberately literal: the same loop bounds, the same
//! rescaling, the same fall-through.
//!
//! # Confidence is about placement, not about matching
//!
//! An exact match against a reference with no internal repeat scores up to 91. The same kind of
//! exact match against `ACACACAC...` scores 4 everywhere, because the model cannot tell the
//! placements apart. That contrast is the whole point of the algorithm, and both cases are in the
//! golden.
//!
//! # And most of the emission table is one
//!
//! ```java
//! for ( int i = 0; i < 256; i++ )
//!     for ( int j = 0; j < 256; j++ )
//!         for ( int q = 0; q <= SAMUtils.MAX_PHRED_SCORE; q++ )
//!             EPSILONS[i][j][q] = 1.0;
//! for ( char b1 : "ACGTacgt".toCharArray() ) { ... }
//! ```
//!
//! Only the sixteen ACGT-by-acgt pairs are ever filled. **An `N` against an `A` therefore has an
//! emission probability of one**, which is not a special case in the algorithm: it is what the
//! uninitialised table holds. The golden shows the consequence: a read with an `N` scores
//! `43,46,78,83,86,77,62,45` where the same read with a real mismatch drops to
//! `43,46,60,49,46,46,46,42`. The `N` costs almost nothing and the mismatch costs a lot.
//!
//! The floor at `minBaseQual` is applied to the **emission's** quality only, leaving the read's own
//! byte alone: `qual2prob[q < minBaseQual ? minBaseQual : q]`.

use std::sync::LazyLock;

use htsjdk_bam::cigar::Op;
use htsjdk_bam::record::BamRecord;

use crate::interval::SimpleInterval;
use crate::math_utils::pow10;
use crate::read_utils;

/// `BAQ.BAQ_TAG`.
pub const BAQ_TAG: &str = "BQ";

/// `BAQ.DEFAULT_GOP`, **Phred scaled** since 2011 and converted in the constructor.
pub const DEFAULT_GOP: f64 = 40.0;

/// `BAQ.DEFAULT_BANDWIDTH`.
pub const DEFAULT_BANDWIDTH: i32 = 7;

/// `SAMUtils.MAX_PHRED_SCORE`, the last quality the emission table is filled for.
const MAX_PHRED_SCORE: usize = 93;

/// The two emission constants, named as the reference names them.
const EM: f64 = 0.33333333333;
const EI: f64 = 0.25;

/// `BAQ.qual2prob`: `10^(-q/10)` for every byte value.
static QUAL_TO_PROB: LazyLock<[f64; 256]> = LazyLock::new(|| {
    let mut table = [0.0; 256];
    for (i, slot) in table.iter_mut().enumerate() {
        *slot = pow10(-(i as f64) / 10.0);
    }
    table
});

/// The quality-to-probability cache, for the suite that compares it entry by entry.
pub fn qual_to_prob_table() -> &'static [f64; 256] {
    &QUAL_TO_PROB
}

/// `convertFromPhredScale`.
fn convert_from_phred_scale(x: f64) -> f64 {
    pow10(-x / 10.0)
}

/// What the model refuses. Every one of these is a `GATKException` in the reference, and every one
/// of them says "BUG" because they describe a caller that is already wrong.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BaqError {
    /// `query.length != _iqual.length`.
    QualityLengthMismatch,
    /// `l_query < 1`.
    QueryTooShort(i32),
    /// `qstart < 0`.
    NegativeQueryStart(i32),
    /// `calcBAQFromTag` with no tag and `useRawQualsIfNoBAQTag` false. An `IllegalStateException`,
    /// not a `GATKException`.
    MissingTag(String),
    /// `UserException.MalformedRead`: the tag's difference takes the quality below zero. Only a
    /// hand-written tag can reach this, because `encode_bq_tag` cannot produce one.
    BaqLargerThanQuality { name: String, locus: String },
}

impl BaqError {
    pub fn message(&self) -> String {
        match self {
            BaqError::QualityLengthMismatch => {
                "BUG: read sequence length != qual length".to_string()
            }
            BaqError::QueryTooShort(length) => {
                format!("BUG: length of query sequence < 0: {length}")
            }
            BaqError::NegativeQueryStart(start) => {
                format!("BUG: query sequence start < 0: {start}")
            }
            BaqError::MissingTag(name) => {
                format!("Required BAQ tag to be present, but none was on read {name}")
            }
            // `UserException.MalformedRead` prefixes the read's name and locus before the reason.
            BaqError::BaqLargerThanQuality { name, locus } => format!(
                "Read {name} {locus} is malformed: BAQ tag error: the BAQ value is larger than the \
                 base quality"
            ),
        }
    }
}

/// `BAQ`: the model with its parameters and its emission table.
#[derive(Debug, Clone)]
pub struct Baq {
    /// `cd`, the gap **open** probability, already converted out of Phred scale.
    gap_open_prob: f64,
    /// `ce`, the gap extension probability.
    gap_extension_prob: f64,
    /// `cb`, the band width.
    band_width: i32,
    /// Every quality below this is raised to it, **for the emission only**.
    min_base_qual: u8,
    /// `EPSILONS`, flattened. See the module note: everything outside the sixteen ACGT pairs is one.
    epsilons: Vec<f64>,
}

impl Default for Baq {
    /// `new BAQ()`, which is `new BAQ(DEFAULT_GOP)`.
    fn default() -> Self {
        Baq::new(DEFAULT_GOP)
    }
}

impl Baq {
    /// `new BAQ(gapOpenPenalty)`, whose argument is Phred scaled.
    pub fn new(gap_open_penalty: f64) -> Baq {
        Baq::with_parameters(
            convert_from_phred_scale(gap_open_penalty),
            0.1,
            DEFAULT_BANDWIDTH,
            4,
        )
    }

    /// `new BAQ(d, e, b, minBaseQual)`, whose gap open probability is **not** Phred scaled.
    pub fn with_parameters(
        gap_open_prob: f64,
        gap_extension_prob: f64,
        band_width: i32,
        min_base_qual: u8,
    ) -> Baq {
        let mut baq = Baq {
            gap_open_prob,
            gap_extension_prob,
            band_width,
            min_base_qual,
            epsilons: Vec::new(),
        };
        baq.initialize_cached_data();
        baq
    }

    pub fn min_base_qual(&self) -> u8 {
        self.min_base_qual
    }

    pub fn gap_open_prob(&self) -> f64 {
        self.gap_open_prob
    }

    pub fn gap_extension_prob(&self) -> f64 {
        self.gap_extension_prob
    }

    pub fn band_width(&self) -> i32 {
        self.band_width
    }

    /// `initializeCachedData`: fill everything with one, then overwrite the sixteen base pairs.
    fn initialize_cached_data(&mut self) {
        let qualities = MAX_PHRED_SCORE + 1;
        self.epsilons = vec![1.0; 256 * 256 * qualities];
        for first in b"ACGTacgt" {
            for second in b"ACGTacgt" {
                for q in 0..qualities {
                    let floored = QUAL_TO_PROB[q.max(self.min_base_qual as usize)];
                    let same = first.eq_ignore_ascii_case(second);
                    let e = if same { 1.0 - floored } else { floored * EM };
                    self.epsilons[(*first as usize) * 256 * qualities
                        + (*second as usize) * qualities
                        + q] = e;
                }
            }
        }
    }

    /// `calcEpsilon(ref, read, qual)`.
    ///
    /// The three indexes are Java `byte`s widened to `int`, so a base above 127 would index
    /// negatively and throw. Nothing produces one: bases come out of a BAM and are ASCII.
    pub fn calc_epsilon(&self, reference: u8, read: u8, qual: u8) -> f64 {
        let qualities = MAX_PHRED_SCORE + 1;
        self.epsilons
            [(reference as usize) * 256 * qualities + (read as usize) * qualities + qual as usize]
    }

    /// `hmm_glocal(ref, query, qstart, l_query, iqual, state, q)`.
    ///
    /// A forward-backward pass over a banded matrix, three states a cell. Written out rather than
    /// factored, because the reference's own banner says this code is kept in sync with samtools
    /// and any rearrangement would make the two harder to compare.
    ///
    /// **The forward pass rescales and the backward pass divides by the same factors**, which is why
    /// `s` is carried between them and why the posterior at the end is a ratio of differently scaled
    /// numbers.
    ///
    /// Fills `state` and `q` in place from `qstart`, and answers 0 as the reference does.
    #[allow(clippy::too_many_arguments)]
    pub fn hmm_glocal(
        &self,
        reference: &[u8],
        query: &[u8],
        qstart: i32,
        l_query: i32,
        iqual: &[u8],
        state: &mut [i32],
        q: &mut [u8],
    ) -> Result<i32, BaqError> {
        if query.len() != iqual.len() {
            return Err(BaqError::QualityLengthMismatch);
        }
        if l_query < 1 {
            return Err(BaqError::QueryTooShort(l_query));
        }
        if qstart < 0 {
            return Err(BaqError::NegativeQueryStart(qstart));
        }

        let l_ref = reference.len() as i32;
        let qstart = qstart as usize;
        let l_query_usize = l_query as usize;

        // The band width, narrowed and then widened again by three separate tests in this order.
        let mut bw = if l_ref > l_query { l_ref } else { l_query };
        if self.band_width < (l_ref - l_query).abs() {
            bw = (l_ref - l_query).abs() + 3;
        }
        if bw > self.band_width {
            bw = self.band_width;
        }
        if bw < (l_ref - l_query).abs() {
            bw = (l_ref - l_query).abs();
        }
        let bw2 = bw * 2 + 1;

        let width = (bw2 * 3 + 6) as usize;
        let mut f = vec![vec![0.0f64; width]; l_query_usize + 1];
        let mut b = vec![vec![0.0f64; width]; l_query_usize + 1];
        let mut s = vec![0.0f64; l_query_usize + 2];

        let s_m = 1.0 / (2.0 * l_query as f64 + 2.0);
        let s_i = s_m;
        let b_m = (1.0 - self.gap_open_prob) / l_ref as f64;
        let b_i = self.gap_open_prob / l_ref as f64;

        let cd = self.gap_open_prob;
        let ce = self.gap_extension_prob;
        let mut m = [0.0f64; 9];
        m[0] = (1.0 - cd - cd) * (1.0 - s_m);
        m[1] = cd * (1.0 - s_m);
        m[2] = cd * (1.0 - s_m);
        m[3] = (1.0 - ce) * (1.0 - s_i);
        m[4] = ce * (1.0 - s_i);
        m[5] = 0.0;
        m[6] = 1.0 - ce;
        m[7] = 0.0;
        m[8] = ce;

        // f[0]
        s[0] = 1.0;
        f[0][set_u(bw, 0, 0) as usize] = 1.0;

        // f[1]
        {
            let beg = 1;
            let end = if l_ref < bw + 1 { l_ref } else { bw + 1 };
            let mut sum = 0.0;
            for k in beg..=end {
                let e =
                    self.calc_epsilon(reference[(k - 1) as usize], query[qstart], iqual[qstart]);
                let u = set_u(bw, 1, k) as usize;
                f[1][u] = e * b_m;
                f[1][u + 1] = EI * b_i;
                sum += f[1][u] + f[1][u + 1];
            }
            s[1] = sum;
            let begin = set_u(bw, 1, beg) as usize;
            let stop = set_u(bw, 1, end) as usize + 2;
            // The rescale is over a RANGE OF CELLS of one row, not over its values: the bounds come
            // from the band and the row is longer than them.
            for slot in f[1][begin..=stop].iter_mut() {
                *slot /= sum;
            }
        }

        // f[2..l_query]
        for i in 2..=l_query {
            let qyi = query[qstart + i as usize - 1];
            let mut beg = 1;
            let mut end = l_ref;
            let x = i - bw;
            beg = if beg > x { beg } else { x };
            let x = i + bw;
            end = if end < x { end } else { x };
            let mut sum = 0.0;
            for k in beg..=end {
                let e = self.calc_epsilon(
                    reference[(k - 1) as usize],
                    qyi,
                    iqual[qstart + i as usize - 1],
                );
                let u = set_u(bw, i, k) as usize;
                let v11 = set_u(bw, i - 1, k - 1) as usize;
                let v10 = set_u(bw, i - 1, k) as usize;
                let v01 = set_u(bw, i, k - 1) as usize;
                let (previous, current) = f.split_at_mut(i as usize);
                let fi1 = &previous[i as usize - 1];
                let fi = &mut current[0];
                fi[u] = e * (m[0] * fi1[v11] + m[3] * fi1[v11 + 1] + m[6] * fi1[v11 + 2]);
                fi[u + 1] = EI * (m[1] * fi1[v10] + m[4] * fi1[v10 + 1]);
                fi[u + 2] = m[2] * fi[v01] + m[8] * fi[v01 + 2];
                sum += fi[u] + fi[u + 1] + fi[u + 2];
            }
            s[i as usize] = sum;
            let begin = set_u(bw, i, beg) as usize;
            let stop = set_u(bw, i, end) as usize + 2;
            // `sum = 1./sum` then a multiply, not a divide: the two differ in the last bits.
            let scale = 1.0 / sum;
            for slot in f[i as usize][begin..=stop].iter_mut() {
                *slot *= scale;
            }
        }

        // f[l_query+1]
        {
            let mut sum = 0.0;
            for k in 1..=l_ref {
                let u = set_u(bw, l_query, k);
                if u < 3 || u >= bw2 * 3 + 3 {
                    continue;
                }
                let u = u as usize;
                sum += f[l_query_usize][u] * s_m + f[l_query_usize][u + 1] * s_i;
            }
            s[l_query_usize + 1] = sum;
        }

        // b[l_query]
        for k in 1..=l_ref {
            let u = set_u(bw, l_query, k);
            if u < 3 || u >= bw2 * 3 + 3 {
                continue;
            }
            let u = u as usize;
            b[l_query_usize][u] = s_m / s[l_query_usize] / s[l_query_usize + 1];
            b[l_query_usize][u + 1] = s_i / s[l_query_usize] / s[l_query_usize + 1];
        }

        // b[l_query-1..1]
        for i in (1..=l_query - 1).rev() {
            let mut beg = 1;
            let mut end = l_ref;
            // `y` is zero on the first row, which is what stops the deletion state propagating
            // past the start.
            let y = if i > 1 { 1.0 } else { 0.0 };
            let qyi1 = query[qstart + i as usize];
            let x = i - bw;
            beg = if beg > x { beg } else { x };
            let x = i + bw;
            end = if end < x { end } else { x };
            for k in (beg..=end).rev() {
                let u = set_u(bw, i, k) as usize;
                let v11 = set_u(bw, i + 1, k + 1) as usize;
                let v10 = set_u(bw, i + 1, k) as usize;
                let v01 = set_u(bw, i, k + 1) as usize;
                let epsilon = if k >= l_ref {
                    0.0
                } else {
                    self.calc_epsilon(reference[k as usize], qyi1, iqual[qstart + i as usize])
                };
                let (current, next) = b.split_at_mut(i as usize + 1);
                let bi = &mut current[i as usize];
                let bi1 = &next[0];
                // `bi1[v11]` is folded into `e`, exactly as the reference's comment says.
                let e = epsilon * bi1[v11];
                bi[u] = e * m[0] + EI * m[1] * bi1[v10 + 1] + m[2] * bi[v01 + 2];
                bi[u + 1] = e * m[3] + EI * m[4] * bi1[v10 + 1];
                bi[u + 2] = (e * m[6] + m[8] * bi[v01 + 2]) * y;
            }
            let begin = set_u(bw, i, beg) as usize;
            let stop = set_u(bw, i, end) as usize + 2;
            let y = 1.0 / s[i as usize];
            for slot in b[i as usize][begin..=stop].iter_mut() {
                *slot *= y;
            }
        }

        // b[0], whose only use is the `pb` the reference prints in its debug line.
        {
            let beg = 1;
            let end = if l_ref < bw + 1 { l_ref } else { bw + 1 };
            let mut sum = 0.0;
            for k in (beg..=end).rev() {
                let u = set_u(bw, 1, k);
                let e =
                    self.calc_epsilon(reference[(k - 1) as usize], query[qstart], iqual[qstart]);
                // The bounds test comes AFTER the epsilon lookup, which is only observable if the
                // lookup could fail. It cannot: k is inside the reference.
                if u < 3 || u >= bw2 * 3 + 3 {
                    continue;
                }
                let u = u as usize;
                sum += e * b[1][u] * b_m + EI * b[1][u + 1] * b_i;
            }
            b[0][set_u(bw, 0, 0) as usize] = sum / s[0];
        }

        // MAP
        for i in 1..=l_query {
            let mut sum = 0.0;
            let mut max = 0.0;
            let mut beg = 1;
            let mut end = l_ref;
            let mut max_k: i32 = -1;
            let x = i - bw;
            beg = if beg > x { beg } else { x };
            let x = i + bw;
            end = if end < x { end } else { x };
            for k in beg..=end {
                let u = set_u(bw, i, k) as usize;
                let z = f[i as usize][u] * b[i as usize][u];
                sum += z;
                if z > max {
                    max = z;
                    // `| 0` in the reference, kept as a bare shift because the match state is
                    // zero and the insertion state below is one.
                    max_k = (k - 1) << 2;
                }
                let z = f[i as usize][u + 1] * b[i as usize][u + 1];
                sum += z;
                if z > max {
                    max = z;
                    max_k = ((k - 1) << 2) | 1;
                }
            }
            max /= sum;
            state[qstart + i as usize - 1] = max_k;
            // `-4.343 * log(1 - max) + .499` is `10*log10(1-max)` with the rounding folded in, and
            // the cast truncates rather than rounding.
            let k = (-4.343 * jmath::math::log(1.0 - max) + 0.499) as i32;
            q[qstart + i as usize - 1] = if k > 100 {
                99
            } else if k < self.min_base_qual as i32 {
                self.min_base_qual
            } else {
                k as u8
            };
        }

        Ok(0)
    }
}

/// `set_u(b, i, k)`: the offset of cell `(i, k)` in the banded matrix, three states a cell.
fn set_u(b: i32, i: i32, k: i32) -> i32 {
    let x = i - b;
    let x = if x > 0 { x } else { 0 };
    (k + 1 - x) * 3
}

/// `stateIsIndel(state)`: the low two bits.
pub fn state_is_indel(state: i32) -> bool {
    (state & 3) != 0
}

/// `stateAlignedPosition(state)`: everything above them.
pub fn state_aligned_position(state: i32) -> i32 {
    state >> 2
}

/// `encodeBQTag(read, baq)`: the **difference** between the quality and the BAQ, not the BAQ.
///
/// `(char)(quality - baq + 64)`, so a read whose BAQ equals its quality carries a tag of all `@`,
/// and a BAQ above the quality encodes below `@` rather than being refused.
pub fn encode_bq_tag(qualities: &[u8], baq: &[u8]) -> String {
    qualities
        .iter()
        .zip(baq)
        .map(|(quality, baq)| {
            // The arithmetic is on Java `int`s, so it can leave the byte range before the cast.
            let bq = 64 + *quality as i32 - *baq as i32;
            char::from_u32(bq as u32).unwrap_or('\u{0}')
        })
        .collect()
}

/// `calcBAQFromTag(read, overwriteOriginalQuals, useRawQualsIfNoBAQTag)`: the tag decoded back into
/// qualities.
///
/// The tag holds a difference, so the quality is `min(rawQual, rawQual - (tag - 64))`. Where there
/// is no tag, `useRawQualsIfNoBAQTag` decides between the raw qualities and a refusal.
pub fn calc_baq_from_tag(
    name: &str,
    locus: &str,
    qualities: &[u8],
    tag: Option<&str>,
    use_raw_quals_if_no_tag: bool,
) -> Result<Vec<u8>, BaqError> {
    match tag {
        Some(tag) => {
            let mut out = Vec::with_capacity(qualities.len());
            for (quality, encoded) in qualities.iter().zip(tag.bytes()) {
                // `rawQuals[i] - baq[i] + 64`, with no clamp at either end: a BAQ above the
                // quality is allowed and comes back above it, and only a negative result is
                // refused.
                let value = *quality as i32 - encoded as i32 + 64;
                if value < 0 {
                    return Err(BaqError::BaqLargerThanQuality {
                        name: name.to_string(),
                        locus: locus.to_string(),
                    });
                }
                out.push(value as u8);
            }
            Ok(out)
        }
        None if use_raw_quals_if_no_tag => Ok(qualities.to_vec()),
        None => Err(BaqError::MissingTag(name.to_string())),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn most_of_the_emission_table_is_one() {
        let baq = Baq::default();
        // An N against an A was never filled in, so it emits with probability one and costs the
        // model nothing.
        assert_eq!(baq.calc_epsilon(b'A', b'N', 20), 1.0);
        assert_eq!(baq.calc_epsilon(b'-', b'A', 20), 1.0);
        // A real match and a real mismatch are not one.
        assert!(baq.calc_epsilon(b'A', b'A', 20) < 1.0);
        assert!(baq.calc_epsilon(b'A', b'C', 20) < 1.0);
        // Case is ignored on both sides.
        assert_eq!(
            baq.calc_epsilon(b'a', b'A', 20),
            baq.calc_epsilon(b'A', b'a', 20)
        );
    }

    #[test]
    fn the_emission_quality_is_floored_and_the_reads_is_not() {
        let baq = Baq::default();
        // Everything below four emits as four.
        for q in 0..4u8 {
            assert_eq!(
                baq.calc_epsilon(b'A', b'A', q),
                baq.calc_epsilon(b'A', b'A', 4)
            );
        }
        assert_ne!(
            baq.calc_epsilon(b'A', b'A', 5),
            baq.calc_epsilon(b'A', b'A', 4)
        );
    }

    #[test]
    fn confidence_is_about_placement_and_not_about_matching() {
        let baq = Baq::default();
        let run = |reference: &[u8], query: &[u8]| -> Vec<u8> {
            let quals = vec![40u8; query.len()];
            let mut state = vec![0i32; query.len()];
            let mut q = vec![0u8; query.len()];
            baq.hmm_glocal(
                reference,
                query,
                0,
                query.len() as i32,
                &quals,
                &mut state,
                &mut q,
            )
            .unwrap();
            q
        };
        // A unique placement is confident.
        let unique = run(b"GATTACAGGCTCTAGCAT", b"TTACAGGC");
        assert!(unique.iter().any(|q| *q > 50), "{unique:?}");
        // An ambiguous one is floored everywhere.
        let repeated = run(b"ACACACACACACACACACAC", b"ACACAC");
        assert!(repeated.iter().all(|q| *q == 4), "{repeated:?}");
    }

    #[test]
    fn the_tag_is_the_difference_and_not_the_value() {
        // Equal quality and BAQ is the zero of the encoding.
        assert_eq!(encode_bq_tag(&[40, 40, 40, 40], &[40, 40, 40, 40]), "@@@@");
        // A BAQ above the quality encodes below `@`.
        assert_eq!(encode_bq_tag(&[40], &[41]), "?");
        assert_eq!(
            calc_baq_from_tag("r", "chr1:1-4", &[40, 40, 40, 40], Some("@@@@"), false).unwrap(),
            vec![40, 40, 40, 40]
        );
        // A BAQ above the quality comes back above it: there is no clamp on that side.
        assert_eq!(
            calc_baq_from_tag("r", "chr1:1-1", &[40], Some("?"), false).unwrap(),
            vec![41]
        );
        assert_eq!(
            calc_baq_from_tag("r", "chr1:1-1", &[40], None, false)
                .unwrap_err()
                .message(),
            "Required BAQ tag to be present, but none was on read r"
        );
        assert_eq!(
            calc_baq_from_tag("r", "chr1:1-1", &[40], None, true).unwrap(),
            vec![40]
        );
    }

    #[test]
    fn the_state_is_bit_packed() {
        assert!(!state_is_indel(0));
        assert!(state_is_indel(1));
        assert_eq!(state_aligned_position(400), 100);
        assert_eq!(state_aligned_position(401), 100);
        assert!(state_is_indel(401));
    }
}

// ---------------------------------------------------------------------------------------------
// The per-read calculation
// ---------------------------------------------------------------------------------------------

/// `ReadUtils.getFirstInsertionOffset`, which indexes the first cigar element without checking
/// there is one.
fn first_insertion_offset(record: &BamRecord) -> Option<i32> {
    let first = record.cigar.elements.first()?;
    Some(if first.op == Op::I {
        first.length as i32
    } else {
        0
    })
}

/// `getReferenceWindowForRead(read, bandWidth)`.
///
/// **Not the read's span.** The start is the alignment start minus half the band width minus the
/// **first** insertion's offset, clamped at one; the stop is the end plus half the band width plus
/// the **last** insertion's offset, unclamped. A ten-base read at 50 needs 47-62; the same read
/// with a leading `3I` needs 44-59.
///
/// The contig is a parameter because a `BamRecord` carries a reference **index** and the reference's
/// `read.getContig()` resolves it through the header.
pub fn reference_window_for_read(
    record: &BamRecord,
    contig: &str,
    band_width: i32,
) -> Option<SimpleInterval> {
    let offset = band_width / 2;
    let start = (read_utils::start(record) - offset - first_insertion_offset(record)?).max(1);
    let stop = read_utils::end(record) + offset + read_utils::last_insertion_offset(record)?;
    Some(SimpleInterval {
        contig: contig.to_string(),
        start,
        end: stop,
    })
}

/// `calculateQueryRange(read)`: the first and last read offsets the model runs over.
///
/// Soft clips move the offset but are **not** included in the range, and a cigar with an `N` or a
/// read clipped entirely away answers nothing at all, which is how those reads get no BAQ.
fn calculate_query_range(record: &BamRecord) -> Option<(i32, i32)> {
    let mut query_start: i32 = -1;
    let mut query_stop: i32 = -1;
    let mut read_i: i32 = 0;
    for element in &record.cigar.elements {
        match element.op {
            Op::N => return None,
            Op::H | Op::P | Op::D => {}
            Op::I | Op::S | Op::M | Op::Eq | Op::X => {
                let previous = read_i;
                read_i += element.length as i32;
                if element.op != Op::S {
                    if query_start == -1 {
                        query_start = previous;
                    }
                    query_stop = read_i;
                }
            }
        }
    }
    if query_stop == query_start {
        return None;
    }
    Some((query_start, query_stop))
}

/// `BAQCalculationResult`: what one read's pass produced.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaqCalculationResult {
    pub raw_quals: Vec<u8>,
    pub bq: Vec<u8>,
    pub state: Vec<i32>,
}

impl Baq {
    /// `calcBAQFromHMM(ref, query, quals, queryStart, queryEnd)`.
    pub fn calc_baq_from_hmm_bases(
        &self,
        reference: &[u8],
        query: &[u8],
        quals: &[u8],
        query_start: i32,
        query_end: i32,
    ) -> Result<BaqCalculationResult, BaqError> {
        let mut result = BaqCalculationResult {
            raw_quals: quals.to_vec(),
            bq: vec![0; quals.len()],
            state: vec![0; quals.len()],
        };
        let query_len = query_end - query_start;
        self.hmm_glocal(
            reference,
            query,
            query_start,
            query_len,
            &result.raw_quals.clone(),
            &mut result.state,
            &mut result.bq,
        )?;
        Ok(result)
    }

    /// `calcBAQFromHMM(read, ref, refOffset)`: the model over one read, then the capping walk.
    ///
    /// Two things in the walk are the reference's rather than the obvious implementation's:
    ///
    ///  * **an insertion or a soft clip keeps its raw quality**, because the loop copies `rawQuals`
    ///    over `bq` for those elements, and `case S:` **falls through** into `case I:` after moving
    ///    the reference. The reference's own comment asks whether that is really intended;
    ///  * **a base that did not align where it was expected is floored**, whatever the model said:
    ///    `capBaseByBAQ` takes `minBaseQual` for an indel state or a position mismatch, and the
    ///    minimum of the raw and the model's quality otherwise.
    ///
    /// A cigar whose lengths do not cover the read replaces the whole answer with the raw
    /// qualities, which is the reference's "odd cigar string" line.
    pub fn calc_baq_from_hmm(
        &self,
        record: &BamRecord,
        reference: &[u8],
        ref_offset: i32,
    ) -> Result<Option<BaqCalculationResult>, BaqError> {
        let Some((query_start, query_end)) = calculate_query_range(record) else {
            return Ok(None);
        };

        let mut result = self.calc_baq_from_hmm_bases(
            reference,
            &record.read_bases,
            &record.base_qualities,
            query_start,
            query_end,
        )?;

        let mut read_i: i32 = 0;
        let mut ref_i: i32 = 0;
        for element in &record.cigar.elements {
            let l = element.length as i32;
            match element.op {
                Op::N => return Ok(None),
                Op::H | Op::P => {}
                // `case S:` moves the reference and then falls through into `case I:`.
                Op::S | Op::I => {
                    if element.op == Op::S {
                        ref_i += l;
                    }
                    for i in read_i..read_i + l {
                        result.bq[i as usize] = result.raw_quals[i as usize];
                    }
                    read_i += l;
                }
                Op::D => ref_i += l,
                Op::M | Op::Eq | Op::X => {
                    for i in read_i..read_i + l {
                        let expected = ref_i - ref_offset + (i - read_i);
                        result.bq[i as usize] = self.cap_base_by_baq(
                            result.raw_quals[i as usize],
                            result.bq[i as usize],
                            result.state[i as usize],
                            expected,
                        );
                    }
                    read_i += l;
                    ref_i += l;
                }
            }
        }
        if read_i != record.read_bases.len() as i32 {
            result.bq.copy_from_slice(&result.raw_quals);
        }
        Ok(Some(result))
    }

    /// `capBaseByBAQ(oq, bq, state, expectedPos)`.
    pub fn cap_base_by_baq(&self, oq: u8, bq: u8, state: i32, expected_pos: i32) -> u8 {
        if state_is_indel(state) || state_aligned_position(state) != expected_pos {
            self.min_base_qual
        } else if bq < oq {
            bq
        } else {
            oq
        }
    }
}

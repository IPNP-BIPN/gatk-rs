//! `ComposeSTRTableFile`, ported from the tool and `STRDecimationTable` (GATK 4.6.2.0).
//!
//! A reference scanned for short tandem repeats, every site reported with the period and the number
//! of repeats that fit it best. The zip the tool writes is not ported; the sites in it are.
//!
//! # The scan does not look at a position twice, and the sites still overlap
//!
//! ```java
//! for (long pos = seqStart; pos <= seqEnd; pos++) {
//!     final BestPeriodRepeat best = findBestPeriodRepeatCombination(...);
//!     if (best != null) { pos = best.end; emitOrDecimateSTR(...); }
//! }
//! ```
//!
//! The loop jumps past the site it found, so no position starts a search twice. But the search
//! itself reaches BACKWARDS from the position that started it, so the next site can begin at a base
//! the previous one ended on: a homopolymer ending at 9 and a dinucleotide repeat beginning at 9
//! are both reported.
//!
//! # The best period is the one with the most repeats
//!
//! ```java
//! if (newRepeats > repeats || (newRepeats == repeats && newPeriod < period)) { ... }
//! ```
//!
//! Ties go to the shorter period, and the repeat count is an integer division of the span by the
//! period, so trailing bases that do not complete a unit sit inside the interval and count for
//! nothing.
//!
//! # The mask starts at the contig's index
//!
//! ```java
//! for (final int[] masks : nextMasks[i]) { Arrays.fill(masks, i); }
//! ```
//!
//! The counter that decides decimation is per contig, per period and per capped repeat, and it
//! starts at the contig's INDEX rather than at zero. So the first site of the second contig carries
//! mask 1, and under the default table that is what removes it while the same repeat on the first
//! contig is kept.
//!
//! # And the cap changes the mask
//!
//! The counter is indexed by `min(maxRepeat, repeats)`, so lowering `--max-repeat` makes distinct
//! sites share a counter. The masks change, and with them which sites decimation removes, though
//! the repeat REPORTED is never capped.

/// `STRDecimationTable.DEFAULT_DECIMATION_MATRIX`.
pub const DEFAULT_DECIMATION_MATRIX: &[&[i32]] = &[
    &[0],
    &[0, 10, 10, 9, 8, 7, 5, 3, 1, 0],
    &[0, 0, 9, 6, 3, 0],
    &[0, 0, 8, 4, 1, 0],
    &[0, 0, 6, 0],
    &[0, 0, 5, 0],
    &[0, 0, 4, 0],
    &[0, 0, 1, 0],
    &[0],
];

/// `STRDecimationTable`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DecimationTable {
    /// `(1 << bits) - 1` per entry, computed as an int shift exactly as the reference does.
    masks: Vec<Vec<i64>>,
}

impl DecimationTable {
    /// `STRDecimationTable.DEFAULT`.
    pub fn default_table() -> Self {
        DecimationTable::from_matrix(DEFAULT_DECIMATION_MATRIX)
    }

    /// `STRDecimationTable.NONE`, which keeps everything.
    pub fn none() -> Self {
        DecimationTable { masks: Vec::new() }
    }

    pub fn from_matrix(matrix: &[&[i32]]) -> Self {
        DecimationTable {
            masks: matrix
                .iter()
                .map(|row| {
                    row.iter()
                        .map(|bits| i64::from((1i32 << bits) - 1))
                        .collect()
                })
                .collect(),
        }
    }

    /// `decimate`: a bit test, so it keeps one site in every 2^n rather than a fraction.
    ///
    /// A period or a repeat past the end of the table is never decimated.
    pub fn decimate(&self, mask: i64, period: usize, repeats: usize) -> bool {
        let Some(row) = self.masks.get(period) else {
            return false;
        };
        let Some(right) = row.get(repeats) else {
            return false;
        };
        ((mask as i32) & (*right as i32)) != 0 || ((mask >> 32) & (right >> 32)) != 0
    }
}

/// `Nucleotide.same`, which compares the decoded values and is false for anything undecodable.
fn same(left: u8, right: u8) -> bool {
    decode(left).is_some() && decode(left) == decode(right)
}

/// The IUPAC letters the reference decodes, upper or lower case.
fn decode(base: u8) -> Option<u8> {
    let upper = base.to_ascii_uppercase();
    if b"ACGTUMRWSYKVHDBN".contains(&upper) {
        // U decodes to the same value as T.
        Some(if upper == b'U' { b'T' } else { upper })
    } else {
        None
    }
}

/// `Nucleotide.isStandard`.
fn is_standard(base: u8) -> bool {
    matches!(
        decode(base),
        Some(b'A') | Some(b'C') | Some(b'G') | Some(b'T')
    )
}

/// `BestPeriodRepeat`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Best {
    pub period: usize,
    pub repeats: usize,
    /// One-based and closed, as the reference keeps them.
    pub start: i64,
    pub end: i64,
}

impl Best {
    fn new(period: usize, start: i64, end: i64) -> Self {
        Best {
            period,
            repeats: ((end - start + 1) / period as i64) as usize,
            start,
            end,
        }
    }

    fn update_if_better(&mut self, period: usize, start: i64, end: i64) {
        let repeats = ((end - start + 1) / period as i64) as usize;
        if repeats > self.repeats || (repeats == self.repeats && period < self.period) {
            *self = Best {
                period,
                repeats,
                start,
                end,
            };
        }
    }
}

/// `findBestPeriodRepeatCombination`, over a one-based position.
pub fn find_best(bases: &[u8], pos: i64, max_period: usize) -> Option<Best> {
    let length = bases.len() as i64;
    let at = |index: i64| bases[(index - 1) as usize];
    // copyBytesAt gives back how many bases it could copy, which is short near the end.
    let max_period_at_pos = std::cmp::min(max_period as i64, length - pos + 1) as usize;
    let first = at(pos);
    if !is_standard(first) {
        return None;
    }
    let mut beg = pos - 1;
    while beg >= 1 && same(at(beg), first) {
        beg -= 1;
    }
    beg += 1;
    let mut end = pos + 1;
    while end <= length && same(at(end), first) {
        end += 1;
    }
    end -= 1;
    let mut best = Best::new(1, beg, end);

    for period in 2..=max_period_at_pos {
        let unit: Vec<u8> = (0..period).map(|offset| at(pos + offset as i64)).collect();
        // The search stops at the first period whose unit reaches a base that is not ACGT.
        if !is_standard(unit[period - 1]) {
            break;
        }
        let mut beg = pos - 1;
        let mut cmp = period - 1;
        while beg >= 1 && same(at(beg), unit[cmp]) {
            beg -= 1;
            if cmp == 0 {
                cmp = period - 1;
            } else {
                cmp -= 1;
            }
        }
        beg += 1;
        let mut cmp = 0usize;
        let mut end = pos + period as i64;
        while end <= length && same(at(end), unit[cmp]) {
            end += 1;
            cmp += 1;
            if cmp == period {
                cmp = 0;
            }
        }
        end -= 1;
        best.update_if_better(period, beg, end);
    }
    Some(best)
}

/// One emitted site.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Locus {
    pub contig_index: usize,
    pub start: i64,
    pub end: i64,
    pub period: usize,
    pub repeats: usize,
    pub mask: i64,
}

/// The counters `initializeMasks` fills, one per contig, period and capped repeat.
#[derive(Debug, Clone)]
pub struct Masks {
    counters: Vec<Vec<Vec<i32>>>,
}

impl Masks {
    pub fn new(contigs: usize, max_period: usize, max_repeat: usize) -> Self {
        Masks {
            counters: (0..contigs)
                .map(|index| vec![vec![index as i32; max_repeat + 1]; max_period + 1])
                .collect(),
        }
    }

    fn next(&mut self, contig: usize, period: usize, repeat: usize) -> i32 {
        let slot = &mut self.counters[contig][period][repeat];
        let value = *slot;
        *slot += 1;
        value
    }
}

/// What one contig's scan produced: the sites kept, and the ones decimation removed.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Scan {
    pub emitted: Vec<Locus>,
    /// The period and repeat of every site decimation removed, in the order they were found.
    pub decimated: Vec<(usize, usize)>,
}

/// The settings one scan runs under.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Settings {
    pub max_period: usize,
    pub max_repeat: usize,
}

/// `traverseInterval` plus `emitOrDecimateSTR`, over one contig.
///
/// `intervals` are one-based and closed. They decide where the scan STARTS, not how far a site may
/// reach, so a site can begin before an interval and end after it.
pub fn scan_contig(
    bases: &[u8],
    contig_index: usize,
    intervals: &[(i64, i64)],
    settings: Settings,
    table: &DecimationTable,
    masks: &mut Masks,
    into: &mut Scan,
) {
    let Settings {
        max_period,
        max_repeat,
    } = settings;
    for (start, stop) in intervals {
        let mut pos = *start;
        while pos <= *stop {
            if let Some(best) = find_best(bases, pos, max_period) {
                pos = best.end;
                let effective = std::cmp::min(max_repeat, best.repeats);
                let mask = i64::from(masks.next(contig_index, best.period, effective));
                if table.decimate(mask, best.period, best.repeats) {
                    into.decimated.push((best.period, best.repeats));
                } else {
                    into.emitted.push(Locus {
                        contig_index,
                        start: best.start,
                        end: best.end,
                        period: best.period,
                        repeats: best.repeats,
                        mask,
                    });
                }
            }
            pos += 1;
        }
    }
}

/// The whole traversal, contig by contig in the dictionary's order.
pub fn scan(
    contigs: &[(String, Vec<u8>)],
    intervals: &[(String, i64, i64)],
    settings: Settings,
    table: &DecimationTable,
) -> Scan {
    let mut masks = Masks::new(contigs.len(), settings.max_period, settings.max_repeat);
    let mut scan = Scan::default();
    for (index, (name, bases)) in contigs.iter().enumerate() {
        let chosen: Vec<(i64, i64)> = if intervals.is_empty() {
            vec![(1, bases.len() as i64)]
        } else {
            intervals
                .iter()
                .filter(|(contig, _, _)| contig == name)
                .map(|(_, start, stop)| (*start, *stop))
                .collect()
        };
        if chosen.is_empty() {
            continue;
        }
        scan_contig(
            bases, index, &chosen, settings, table, &mut masks, &mut scan,
        );
    }
    scan
}

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
    /// The bits per period and repeat, as given; what `print` writes and `decimationBit` reads.
    matrix: Vec<Vec<i32>>,
    /// `(1 << bits) - 1` per entry, computed as an int shift exactly as the reference does.
    masks: Vec<Vec<i64>>,
}

/// What reading a decimation file refuses, with the reference's own wording.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DecimationError {
    /// `UserException.BadInput`.
    BadInput(String),
}

impl DecimationTable {
    /// `STRDecimationTable.DEFAULT`.
    pub fn default_table() -> Self {
        DecimationTable::from_matrix(DEFAULT_DECIMATION_MATRIX)
    }

    /// `STRDecimationTable.NONE`, whose matrix is a single zero, so nothing past period zero is
    /// ever decimated.
    pub fn none() -> Self {
        DecimationTable::from_rows(vec![vec![0]])
    }

    pub fn from_matrix(matrix: &[&[i32]]) -> Self {
        DecimationTable::from_rows(matrix.iter().map(|row| row.to_vec()).collect())
    }

    pub fn from_rows(matrix: Vec<Vec<i32>>) -> Self {
        let masks = matrix
            .iter()
            .map(|row| {
                row.iter()
                    .map(|bits| i64::from(1i32.wrapping_shl(*bits as u32).wrapping_sub(1)))
                    .collect()
            })
            .collect();
        DecimationTable { matrix, masks }
    }

    /// `new STRDecimationTable(spec)` over a file's text: the lines that are neither blank nor a
    /// `#` comment, each split on runs of whitespace.
    ///
    /// `path` is only what the messages quote.
    pub fn parse(text: &str, path: &str) -> Result<Self, DecimationError> {
        let rows: Vec<Vec<&str>> = text
            .lines()
            .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
            .map(split_on_whitespace)
            .collect();
        if rows.is_empty() {
            return Ok(DecimationTable::from_rows(Vec::new()));
        }
        let mut matrix = Vec::new();
        let mut total = 0;
        for (i, row) in rows.iter().enumerate() {
            let mut values = Vec::new();
            for (j, cell) in row.iter().enumerate() {
                let bad = |details: &str| {
                    DecimationError::BadInput(format!(
                        "bad decimation value found in {path} for period and repeats ({i}, {j}) \
                         with string ({cell}): {details}"
                    ))
                };
                let value: i32 = cell
                    .parse()
                    .map_err(|_| bad("not a valid double literal"))?;
                if value < 0 {
                    return Err(bad("negatives are not allowed"));
                }
                values.push(value);
                total += 1;
            }
            matrix.push(values);
        }
        if total == 0 {
            return Err(DecimationError::BadInput(format!(
                "the input decimation matrix does contain any values:{path}"
            )));
        }
        Ok(DecimationTable::from_rows(matrix))
    }

    /// `print`: each row's values joined by tabs, one row per line.
    pub fn print(&self) -> String {
        let mut out = String::new();
        for row in &self.matrix {
            let cells: Vec<String> = row.iter().map(|value| value.to_string()).collect();
            out.push_str(&cells.join("\t"));
            out.push('\n');
        }
        out
    }

    /// `decimationBit`, zero past the end of the table.
    pub fn decimation_bit(&self, period: usize, repeats: usize) -> i32 {
        self.matrix
            .get(period)
            .and_then(|row| row.get(repeats))
            .copied()
            .unwrap_or(0)
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

/// `String.split("\\s+")`: a leading run of whitespace leaves an empty first field, a trailing one
/// leaves nothing.
fn split_on_whitespace(line: &str) -> Vec<&str> {
    let mut fields: Vec<&str> = line.split(char::is_whitespace).collect();
    let leading_empty = fields.first().is_some_and(|field| field.is_empty());
    fields = std::iter::once(if leading_empty { Some("") } else { None })
        .flatten()
        .chain(fields.into_iter().filter(|field| !field.is_empty()))
        .collect();
    fields
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

/// `DragstrLocusUtils.INDEX_BYTE_INTERVAL`: an index entry at least every 64 KB of sites.
pub const INDEX_BYTE_INTERVAL: u64 = 1 << 16;

/// `sites.bin` and `sites.idx`, as `DragstrLocusUtils.binaryWriter` lays them out.
///
/// A site is 23 big-endian bytes: the contig index as an int, the start as a long, the period as a
/// byte, the length in bases as a short and the mask as a long. The index holds an entry, the
/// contig index and the start as ints and the byte offset as a long, where a contig begins and then
/// wherever 64 KB have passed since the last entry.
pub fn sites_binary(loci: &[Locus]) -> (Vec<u8>, Vec<u8>) {
    let mut sites = Vec::with_capacity(loci.len() * 23);
    let mut index = Vec::new();
    let mut last_contig: i64 = -1;
    let mut last_entry_offset: u64 = 0;
    for locus in loci {
        let offset = sites.len() as u64;
        let contig = locus.contig_index as i64;
        if contig != last_contig || offset - last_entry_offset >= INDEX_BYTE_INTERVAL {
            last_contig = contig;
            index.extend_from_slice(&(locus.contig_index as i32).to_be_bytes());
            index.extend_from_slice(&(locus.start as i32).to_be_bytes());
            index.extend_from_slice(&(offset as i64).to_be_bytes());
            last_entry_offset = offset;
        }
        sites.extend_from_slice(&(locus.contig_index as i32).to_be_bytes());
        sites.extend_from_slice(&locus.start.to_be_bytes());
        sites.push(locus.period as u8);
        sites.extend_from_slice(&(locus.length() as i16).to_be_bytes());
        sites.extend_from_slice(&locus.mask.to_be_bytes());
    }
    (sites, index)
}

/// `sites.txt`, `DragstrLocusUtils.textWriter`'s table.
pub fn sites_text(loci: &[Locus], contig_names: &[String]) -> String {
    let mut out =
        String::from("chridx\tchrid\tstart\tend\tperiod\tmask\tmask_bin\tlength_bp\tlength_rp\n");
    for locus in loci {
        let length = locus.length() as i16;
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{:b}\t{}\t{}\n",
            locus.contig_index,
            contig_names
                .get(locus.contig_index)
                .map(String::as_str)
                .unwrap_or(""),
            locus.start,
            locus.start + i64::from(length) - 1,
            locus.period,
            locus.mask,
            locus.mask,
            length,
            i32::from(length) / locus.period as i32,
        ));
    }
    out
}

impl Locus {
    /// `getLength`, the site's span in bases, which the writer narrows to a short.
    pub fn length(&self) -> i64 {
        self.end - self.start + 1
    }
}

/// `STRTableFileBuilder.writeSummary`.
///
/// Every site the scan found is counted in `total`, and the kept ones in `emitted` too, both by the
/// period and the repeat CAPPED at the maxima. `double_to_string` is `Double.toString`, which the
/// actual decimation column is rendered with.
pub fn summary(
    scan: &Scan,
    settings: Settings,
    annotations: &[(String, String)],
    table: &DecimationTable,
    double_to_string: impl Fn(f64) -> String,
) -> String {
    let Settings {
        max_period,
        max_repeat,
    } = settings;
    let mut total = vec![vec![0i64; max_repeat + 1]; max_period + 1];
    let mut emitted = vec![vec![0i64; max_repeat + 1]; max_period + 1];
    for locus in &scan.emitted {
        let period = locus.period.min(max_period);
        let repeats = locus.repeats.min(max_repeat);
        total[period][repeats] += 1;
        emitted[period][repeats] += 1;
    }
    for (period, repeats) in &scan.decimated {
        total[(*period).min(max_period)][(*repeats).min(max_repeat)] += 1;
    }
    let rule = "##########################################################################################\n";
    let mut out = String::new();
    out.push_str(rule);
    out.push_str("# STRTableSummary\n");
    out.push_str("# ---------------------------------------\n");
    out.push_str(&format!("# maxPeriod = {max_period}\n"));
    out.push_str(&format!("# maxRepeatLength = {max_repeat}\n"));
    for (name, value) in annotations {
        out.push_str(&format!("# {name} = {value}\n"));
    }
    out.push_str(rule);
    out.push_str(
        "period\trepeatLength\ttotalCounts\temittedCounts\tintendedDecimation\tactualDecimation\n",
    );
    for period in 1..=max_period {
        let first = if period == 1 { 1 } else { 2 };
        for repeats in first..=max_repeat {
            let all = total[period][repeats];
            let kept = emitted[period][repeats];
            let actual = if all > 0 {
                std::f64::consts::LOG2_E * ((all as f64).ln() - (kept as f64).ln())
            } else {
                0.0
            };
            // `Math.round(x * 100) / 100.0`: a long, saturating, divided back.
            let rounded = java_round(actual * 100.0) as f64 / 100.0;
            out.push_str(&format!(
                "{period}\t{repeats}\t{all}\t{kept}\t{}\t{}\n",
                table.decimation_bit(period, repeats),
                double_to_string(rounded)
            ));
        }
    }
    out
}

/// `Math.round(double)`: `floor(x + 0.5)` saturated to a long, NaN to zero.
fn java_round(value: f64) -> i64 {
    if value.is_nan() {
        return 0;
    }
    let floored = (value + 0.5).floor();
    if floored >= i64::MAX as f64 {
        i64::MAX
    } else if floored <= i64::MIN as f64 {
        i64::MIN
    } else {
        floored as i64
    }
}

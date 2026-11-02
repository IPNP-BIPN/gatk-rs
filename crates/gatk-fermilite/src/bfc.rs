//! `bfc.c`: k-mer counting, one-read error correction, and the unique-k-mer filter.
//!
//! The same function, `fml_correct_core`, serves twice: once to correct with `ec_k`, once with
//! `min_asm_ovlp` (33) to trim every read to its longest streak of counted k-mers, or to empty it.

use crate::htab::{bfc_kmer_append, bfc_kmer_change, bfc_kmer_hash, BfcCh, BfcKmer, BFC_KMER_NULL};
use crate::ksort::{ks_heapdown, ks_heapup};
use crate::{cvttsd2si, Bseq, FmlOpt, SEQ_NT6_TABLE};

const BFC_MAX_PATHS: usize = 4;
const BFC_EC_HIST: usize = 5;
const BFC_EC_HIST_HIGH: usize = 2;
const BFC_EC_MIN_COV_COEF: f64 = 0.1;

/// `bfc_opt_t`, as `bfc_opt_init` leaves it.
struct BfcOpt {
    q: i32,
    k: i32,
    min_cov: i32,
    max_end_ext: i32,
    win_multi_ec: i32,
    min_trim_frac: f32,
    w_ec: i32,
    w_ec_high: i32,
    w_absent: i32,
    w_absent_high: i32,
    max_path_diff: i32,
    max_heap: i32,
}

impl BfcOpt {
    /// `bfc_opt_init`.
    fn init() -> BfcOpt {
        BfcOpt {
            q: 20,
            k: -1,
            min_cov: 4, // in BFC, this defaults to 3 because it has Bloom pre-filter
            max_end_ext: 5,
            win_multi_ec: 10,
            min_trim_frac: 0.8,
            w_ec: 1,
            w_ec_high: 7,
            w_absent: 3,
            w_absent_high: 1,
            max_path_diff: 15,
            max_heap: 100,
        }
    }
}

/// The signed `char` of `s->qual[i] - 33 >= q`.
#[inline]
fn qual_at_least(qual: u8, q: i32) -> bool {
    (qual as i8 as i32) - 33 >= q
}

/// `worker_count`, through `bfc_kmer_insert`.
fn worker_count(ch: &mut BfcCh, s: &Bseq, k: i32, q: i32) {
    let mut x = BFC_KMER_NULL;
    let mut qmer = 0u64;
    let mask = (1u64 << k) - 1;
    let mut l = 0;
    for i in 0..s.l_seq as usize {
        let c = SEQ_NT6_TABLE[s.seq[i] as usize] as i32 - 1;
        if c < 4 {
            bfc_kmer_append(k, &mut x, c);
            qmer = ((qmer << 1) | u64::from(qual_at_least(s.qual[i], q))) & mask;
            l += 1;
            if l >= k {
                ch.insert(&bfc_kmer_hash(k, &x), qmer == mask);
            }
        } else {
            l = 0;
            qmer = 0;
            x = BFC_KMER_NULL;
        }
    }
}

/// `fml_count`.
fn fml_count(seq: &[Bseq], k: i32, q: i32, l_pre: i32) -> BfcCh {
    let mut ch = BfcCh::new(k, l_pre);
    for s in seq {
        worker_count(&mut ch, s, k, q);
    }
    ch
}

/// `ecbase_t`, its `ec`, `absent` and `i` left out: they are written and never read.
#[derive(Clone, Copy, Default)]
struct EcBase {
    b: u8,  // :3
    q: u8,  // :1
    ob: u8, // :3
    oq: u8, // :1
    lcov: u8,
    hcov: u8,
    solid_end: bool,
    high_end: bool,
}

/// `bfc_seq_conv`. `seq_nt6_table[(int)s[i]]` with a signed `char`: see [`crate::seq_nt6_signed`].
fn bfc_seq_conv(s: &[u8], q: &[u8], qthres: i32, seq: &mut Vec<EcBase>) {
    seq.clear();
    for i in 0..s.len() {
        let t = crate::seq_nt6_signed(s[i]);
        let b = t.wrapping_sub(1) & 7;
        let mut qb = u8::from(qual_at_least(q[i], qthres));
        if b > 3 {
            qb = 0;
        }
        seq.push(EcBase {
            b,
            ob: b,
            q: qb,
            oq: qb,
            ..EcBase::default()
        });
    }
}

/// `ecbase_comp`.
#[inline]
fn ecbase_comp(b: &EcBase) -> EcBase {
    let mut r = *b;
    r.b = if b.b < 4 { 3 - b.b } else { 4 };
    r.ob = if b.ob < 4 { 3 - b.ob } else { 4 };
    r
}

/// `bfc_seq_revcomp`.
fn bfc_seq_revcomp(seq: &mut [EcBase]) {
    let n = seq.len();
    for i in 0..n >> 1 {
        let tmp = ecbase_comp(&seq[i]);
        seq[i] = ecbase_comp(&seq[n - 1 - i]);
        seq[n - 1 - i] = tmp;
    }
    if n & 1 != 0 {
        seq[n >> 1] = ecbase_comp(&seq[n >> 1]);
    }
}

/// `bfc_ec_greedy_k`.
fn bfc_ec_greedy_k(k: i32, mode: i32, x: &BfcKmer, ch: &BfcCh) -> i32 {
    let (mut max, mut max_ec, mut max2) = (0i32, -1i32, 0i32);
    for i in 0..k {
        let c = (((x[1] >> i) & 1) << 1 | ((x[0] >> i) & 1)) as i32;
        for j in 0..4 {
            let mut y = *x;
            if j == c {
                continue;
            }
            bfc_kmer_change(k, &mut y, i, j);
            let ret = ch.kmer_occ(&y);
            if ret < 0 {
                continue;
            }
            if (max & 0xff) < (ret & 0xff) {
                max2 = max;
                max = ret;
                max_ec = (i << 2) | j;
            } else if (max2 & 0xff) < (ret & 0xff) {
                max2 = ret;
            }
        }
    }
    if (max & 0xff) * 3 > mode && (max2 & 0xff) < 3 {
        max_ec
    } else {
        -1
    }
}

/// `bfc_ec_first_kmer`.
fn bfc_ec_first_kmer(k: i32, s: &[EcBase], start: i32, x: &mut BfcKmer) -> i32 {
    *x = BFC_KMER_NULL;
    let mut i = start;
    let mut l = 0;
    while (i as usize) < s.len() {
        let c = &s[i as usize];
        if c.b < 4 {
            bfc_kmer_append(k, x, c.b as i32);
            l += 1;
            if l == k {
                break;
            }
        } else {
            l = 0;
            *x = BFC_KMER_NULL;
        }
        i += 1;
    }
    i
}

/// `bfc_ec_kcov`.
fn bfc_ec_kcov(k: i32, min_occ: i32, s: &mut [EcBase], ch: &BfcCh) {
    let mut x = BFC_KMER_NULL;
    let mut l = 0;
    for i in 0..s.len() {
        let c = &mut s[i];
        c.high_end = false;
        c.solid_end = false;
        c.lcov = 0;
        c.hcov = 0;
        if c.b < 4 {
            bfc_kmer_append(k, &mut x, c.b as i32);
            l += 1;
            if l >= k {
                let r = ch.kmer_occ(&x);
                if r >= 0 {
                    if ((r >> 8) & 0x3f) >= min_occ + 1 {
                        c.high_end = true;
                    }
                    if (r & 0xff) >= min_occ {
                        c.solid_end = true;
                        let high = u8::from(c.high_end);
                        for sj in &mut s[i + 1 - k as usize..=i] {
                            sj.lcov = (sj.lcov + 1) & 0x3f;
                            sj.hcov = (sj.hcov + high) & 0x3f;
                        }
                    }
                }
            }
        } else {
            l = 0;
            x = BFC_KMER_NULL;
        }
    }
}

/// `bfc_ec_best_island`.
fn bfc_ec_best_island(k: i32, s: &[EcBase]) -> u64 {
    let (mut l, mut max, mut max_i) = (0i32, 0i32, -1i32);
    let mut i = k - 1;
    while (i as usize) < s.len() {
        if !s[i as usize].solid_end {
            if l > max {
                max = l;
                max_i = i;
            }
            l = 0;
        } else {
            l += 1;
        }
        i += 1;
    }
    if l > max {
        max = l;
        max_i = i;
    }
    if max > 0 {
        (((max_i - max - k + 1) as i64 as u64) << 32) | (max_i as i64 as u64)
    } else {
        0
    }
}

/// `bfc_penalty_t`.
#[derive(Clone, Copy, Default)]
struct Penalty {
    ec: u8,
    ec_high: u8,
    absent: u8,
    absent_high: u8,
    b: u8,
}

/// `echeap1_t`.
#[derive(Clone, Copy)]
struct EcHeap1 {
    tot_pen: i32,
    i: i32,
    k: i32,
    ecpos_high: [i32; BFC_EC_HIST_HIGH],
    ecpos: [i32; BFC_EC_HIST],
    x: BfcKmer,
}

/// `ecstack1_t`, its `cnt` left out: written and never read.
#[derive(Clone, Copy)]
struct EcStack1 {
    parent: i32,
    i: i32,
    tot_pen: i32,
    b: u8,
}

/// `bfc_ec1buf_t`.
struct Ec1Buf<'a> {
    opt: &'a BfcOpt,
    ch: &'a BfcCh,
    heap: Vec<EcHeap1>,
    stack: Vec<EcStack1>,
    seq: Vec<EcBase>,
    ec: [Vec<EcBase>; 2],
    mode: i32,
}

/// `heap_lt`.
#[inline]
fn heap_lt(a: &EcHeap1, b: &EcHeap1) -> bool {
    a.tot_pen > b.tot_pen
}

/// `weighted_penalty`.
#[inline]
fn weighted_penalty(o: &BfcOpt, p: &Penalty) -> i32 {
    o.w_ec * p.ec as i32
        + o.w_ec_high * p.ec_high as i32
        + o.w_absent * p.absent as i32
        + o.w_absent_high * p.absent_high as i32
}

/// `buf_update`.
fn buf_update(
    opt: &BfcOpt,
    heap: &mut Vec<EcHeap1>,
    stack: &mut Vec<EcStack1>,
    prev: &EcHeap1,
    pen: Penalty,
) {
    let tot_pen = prev.tot_pen + weighted_penalty(opt, &pen);
    stack.push(EcStack1 {
        parent: prev.k,
        i: prev.i,
        tot_pen,
        b: pen.b,
    });
    let mut r = EcHeap1 {
        tot_pen,
        i: prev.i + 1,
        k: stack.len() as i32 - 1,
        ecpos_high: prev.ecpos_high,
        ecpos: prev.ecpos,
        x: prev.x,
    };
    if pen.ec_high != 0 {
        r.ecpos_high[1] = prev.ecpos_high[0];
        r.ecpos_high[0] = prev.i;
    }
    if pen.ec != 0 {
        r.ecpos[1..BFC_EC_HIST].copy_from_slice(&prev.ecpos[..BFC_EC_HIST - 1]);
        r.ecpos[0] = prev.i;
    }
    bfc_kmer_append(opt.k, &mut r.x, pen.b as i32);
    heap.push(r);
    ks_heapup(heap, heap_lt);
}

/// `buf_backtrack`, its `ec`, `absent` and return value left out: nothing reads them but the sign,
/// and the sign is never negative.
fn buf_backtrack(s: &[EcStack1], mut end: i32, seq_n: usize, path: &mut [EcBase]) {
    while end >= 0 {
        let e = &s[end as usize];
        if (e.i as usize) < seq_n {
            path[e.i as usize].b = e.b;
        }
        end = e.parent;
    }
}

/// `bfc_ec1dir`, on `e.seq`, into `e.ec[ec_i]`: negative when there is no path.
fn bfc_ec1dir(e: &mut Ec1Buf, ec_i: usize, start: i32, end: i32) -> i32 {
    let o = e.opt;
    let seq = &e.seq;
    let n = seq.len() as i32;
    c_assert!(
        end <= n && end - start >= o.k,
        "bfc_ec1dir: end <= seq->n && end - start >= k"
    );
    e.heap.clear();
    e.stack.clear();
    let mut rv = -1;
    let mut path = [0i32; BFC_MAX_PATHS];
    let mut n_paths = 0usize;
    let mut min_path = -1i32;
    let mut min_path_pen = i32::MAX;
    let mut n_failures = 0;
    let ec = &mut e.ec[ec_i];
    ec.resize(seq.len(), EcBase::default());
    let mut z = EcHeap1 {
        tot_pen: 0,
        i: start,
        k: 0,
        ecpos_high: [-1; BFC_EC_HIST_HIGH],
        ecpos: [-1; BFC_EC_HIST],
        x: BFC_KMER_NULL,
    };
    let mut l = 0;
    while z.i < end {
        let c = seq[z.i as usize].b as i32;
        if c < 4 {
            l += 1;
            if l == o.k {
                break;
            }
            bfc_kmer_append(o.k, &mut z.x, c);
        } else {
            l = 0;
            z.x = BFC_KMER_NULL;
        }
        z.i += 1;
    }
    c_assert!(z.i < end, "bfc_ec1dir: z.i < end");
    z.k = -1;
    e.heap.push(z);
    for i in 0..seq.len() {
        ec[i].b = seq[i].b;
        ec[i].ob = seq[i].ob;
    }
    // exhaustive error correction
    loop {
        let mut stop = false;
        if e.heap.is_empty() {
            // may happen when there is an uncorrectable "N"
            rv = -2;
            break;
        }
        let z = e.heap[0];
        let last = e.heap.pop().unwrap();
        if !e.heap.is_empty() {
            e.heap[0] = last;
        }
        ks_heapdown(0, &mut e.heap, heap_lt);
        if min_path >= 0 && z.tot_pen > min_path_pen + o.max_path_diff {
            break;
        }
        if z.i - end > o.max_end_ext {
            stop = true;
        }
        if !stop {
            let c = if z.i < n {
                Some(&seq[z.i as usize])
            } else {
                None
            };
            let mut os = -1;
            let mut fixed = false;
            let mut other_ext = 0;
            let mut added = [Penalty::default(); 4];
            let mut n_added = 0usize;
            // test if the read extension alone is enough
            if z.i > end {
                fixed = true;
            }
            if let Some(c) = c {
                if c.b < 4 {
                    let mut x = z.x;
                    bfc_kmer_append(o.k, &mut x, c.b as i32);
                    os = e.ch.kmer_occ(&x);
                    if c.q != 0 && (os & 0xff) >= o.min_cov + 1 && c.lcov as i32 >= o.min_cov + 1 {
                        fixed = true;
                    } else if c.hcov as f64 > o.k as f64 * 0.75 {
                        fixed = true;
                    }
                }
            }
            // extension
            for b in 0..4i32 {
                if let (true, Some(c)) = (fixed, c) {
                    if b != c.b as i32 {
                        continue;
                    }
                }
                if c.is_none() || b != c.unwrap().b as i32 {
                    let mut x = z.x;
                    if let Some(c) = c {
                        // not over the end
                        if c.q != 0
                            && z.ecpos_high[BFC_EC_HIST_HIGH - 1] >= 0
                            && z.i - z.ecpos_high[BFC_EC_HIST_HIGH - 1] < o.win_multi_ec
                        {
                            continue; // no close highQ corrections
                        }
                        if z.ecpos[BFC_EC_HIST - 1] >= 0
                            && z.i - z.ecpos[BFC_EC_HIST - 1] < o.win_multi_ec
                        {
                            continue; // no clustered corrections
                        }
                    }
                    bfc_kmer_append(o.k, &mut x, b);
                    let s = e.ch.kmer_occ(&x);
                    if s < 0 || (s & 0xff) < o.min_cov {
                        continue; // not solid
                    }
                    let is_ec = u8::from(matches!(c, Some(c) if c.b < 4));
                    added[n_added] = Penalty {
                        ec: is_ec,
                        ec_high: if is_ec != 0 { c.unwrap().oq } else { 0 },
                        absent: 0,
                        absent_high: u8::from(((s >> 8) & 0xff) < o.min_cov),
                        b: b as u8,
                    };
                    n_added += 1;
                    other_ext += 1;
                } else {
                    added[n_added] = Penalty {
                        ec: 0,
                        ec_high: 0,
                        absent: u8::from(os < 0 || (os & 0xff) < o.min_cov),
                        absent_high: u8::from(os < 0 || ((os >> 8) & 0xff) < o.min_cov),
                        b: b as u8,
                    };
                    n_added += 1;
                }
            }
            if !fixed && other_ext == 0 {
                n_failures += 1;
            }
            if n_failures > n * 2 {
                rv = -3;
                break;
            }
            if c.is_some() || n_added == 1 {
                if n_added > 1 && e.heap.len() as i32 > o.max_heap {
                    // to prevent heap explosion
                    let mut min_b = -1i32;
                    let mut min = i32::MAX;
                    for (b, pen) in added.iter().enumerate().take(n_added) {
                        let t = weighted_penalty(o, pen);
                        if min > t {
                            min = t;
                            min_b = b as i32;
                        }
                    }
                    buf_update(o, &mut e.heap, &mut e.stack, &z, added[min_b as usize]);
                } else {
                    for pen in added.iter().take(n_added) {
                        buf_update(o, &mut e.heap, &mut e.stack, &z, *pen);
                    }
                }
            } else {
                if n_added == 0 {
                    e.stack[z.k as usize].tot_pen += o.w_absent * (o.max_end_ext - (z.i - end));
                }
                stop = true;
            }
        }
        if stop {
            let tp = e.stack[z.k as usize].tot_pen;
            if tp < min_path_pen {
                min_path_pen = tp;
                min_path = n_paths as i32;
            }
            path[n_paths] = z.k;
            n_paths += 1;
            if n_paths == BFC_MAX_PATHS {
                break;
            }
        }
    }
    // backtrack
    if n_paths == 0 {
        return rv;
    }
    c_assert!(
        min_path >= 0
            && (min_path as usize) < n_paths
            && e.stack[path[min_path as usize] as usize].tot_pen == min_path_pen,
        "bfc_ec1dir: min_path"
    );
    buf_backtrack(&e.stack, path[min_path as usize], seq.len(), ec);
    for (i, a) in ec.iter_mut().enumerate() {
        // mask out uncorrected regions
        if (i as i32) < start + o.k || (i as i32) >= end {
            a.b = 4;
        }
    }
    0
}

/// `bfc_ec1`: corrects `seq` and `qual` in place, or leaves both alone.
fn bfc_ec1(e: &mut Ec1Buf, seq: &mut [u8], qual: &mut [u8]) {
    let k = e.opt.k;
    bfc_seq_conv(seq, qual, e.opt.q, &mut e.seq);
    let n = e.seq.len() as i32;
    let mut n_n = 0;
    for c in &e.seq {
        if c.ob > 3 {
            n_n += 1;
        }
    }
    if n_n as f64 > n as f64 * 0.05 {
        return; // ECCODE_MANY_N
    }
    bfc_ec_kcov(k, e.opt.min_cov, &mut e.seq, e.ch);
    let r = bfc_ec_best_island(k, &e.seq);
    let mut start: i32 = 0;
    let mut end: i32;
    if r == 0 {
        // no solid k-mer
        let mut x = BFC_KMER_NULL;
        let mut ec = -1;
        loop {
            end = bfc_ec_first_kmer(k, &e.seq, start, &mut x);
            if end >= n {
                break;
            }
            ec = bfc_ec_greedy_k(k, e.mode, &x, e.ch);
            if ec >= 0 {
                break;
            }
            if end + (k >> 1) >= n {
                break;
            }
            start = end - (k >> 1);
        }
        if ec >= 0 {
            e.seq[(end - (ec >> 2)) as usize].b = (ec & 3) as u8;
            end += 1;
            start = end - k;
        } else {
            return; // ECCODE_NO_SOLID
        }
    } else {
        start = (r >> 32) as i32;
        end = r as u32 as i32;
    }
    if bfc_ec1dir(e, 0, start, n) < 0 {
        return;
    }
    bfc_seq_revcomp(&mut e.seq);
    if bfc_ec1dir(e, 1, n - end, n) < 0 {
        return;
    }
    bfc_seq_revcomp(&mut e.ec[1]);
    bfc_seq_revcomp(&mut e.seq);
    for i in 0..e.seq.len() {
        let (b0, b1) = (e.ec[0][i].b, e.ec[1][i].b);
        let c = &mut e.seq[i];
        if b0 == b1 {
            c.b = if b0 > 3 { c.b } else { b0 };
        } else if b1 > 3 {
            c.b = b0;
        } else if b0 > 3 {
            c.b = b1;
        } else {
            c.b = c.ob;
        }
    }
    for i in 0..e.seq.len() {
        let c = &e.seq[i];
        let is_diff = c.b != c.ob;
        seq[i] = (if is_diff { b"acgtn" } else { b"ACGTN" })[c.b as usize];
        qual[i] = if is_diff {
            34 + c.ob
        } else {
            b"+?"[c.q as usize]
        };
    }
}

/// `max_streak`.
fn max_streak(k: i32, ch: &BfcCh, s: &Bseq) -> u64 {
    let (mut max, mut t) = (0u64, 0u64);
    let mut x = BFC_KMER_NULL;
    let mut l = 0;
    for i in 0..s.l_seq as usize {
        let c = SEQ_NT6_TABLE[s.seq[i] as usize] as i32 - 1;
        if c < 4 {
            bfc_kmer_append(k, &mut x, c);
            l += 1;
            if l >= k {
                if ch.kmer_occ(&x) > 0 {
                    t += 1u64 << 32;
                } else {
                    t = i as u64 + 1;
                }
            } else {
                t = i as u64 + 1;
            }
        } else {
            l = 0;
            x = BFC_KMER_NULL;
            t = i as u64 + 1;
        }
        max = if max > t { max } else { t };
    }
    max
}

/// The filter's branch of `worker_ec`.
fn worker_fltuniq(o: &BfcOpt, ch: &BfcCh, s: &mut Bseq) {
    let max = max_streak(o.k, ch, s);
    if max >> 32 != 0
        && ((max >> 32) + o.k as u64 - 1) as f64 / s.l_seq as f64 > o.min_trim_frac as f64
    {
        let mut start = max as u32 as i32;
        let end = start.wrapping_add((max >> 32) as i32);
        start -= o.k - 1;
        c_assert!(
            start >= 0 && end <= s.l_seq,
            "worker_ec: start >= 0 && end <= s->l_seq"
        );
        s.seq.copy_within(start as usize..end as usize, 0);
        s.l_seq = end - start;
        s.seq.truncate(s.l_seq as usize);
        s.qual.copy_within(start as usize..end as usize, 0);
        s.qual.truncate(s.l_seq as usize);
    } else {
        s.l_seq = 0;
        s.seq = Vec::new();
        s.qual = Vec::new();
    }
}

/// `fml_correct_core`: the k-mer coverage.
pub(crate) fn fml_correct_core(opt: &FmlOpt, flt_uniq: bool, seq: &mut [Bseq]) -> f32 {
    let mut bfc_opt = BfcOpt::init();
    bfc_opt.k = if flt_uniq { opt.min_asm_ovlp } else { opt.ec_k };
    let mut tot_len = 0u64;
    for s in seq.iter() {
        tot_len = tot_len.wrapping_add(s.l_seq as i64 as u64);
    }
    // Unsigned, so a tot_len under 8 wraps to something enormous and picks 20.
    let t = tot_len.wrapping_sub(8);
    let l_pre = if t < 20 { t as i32 } else { 20 };

    let ch = fml_count(seq, bfc_opt.k, bfc_opt.q, l_pre);
    let (hist, mode) = ch.hist();
    let (mut sum_k, mut tot_k) = (0u64, 0u64);
    // `opt->min_cnt` is 4 here and always a valid index.
    for (i, &h) in hist.iter().enumerate().skip(opt.min_cnt as usize) {
        sum_k = sum_k.wrapping_add(h);
        tot_k = tot_k.wrapping_add((i as u64).wrapping_mul(h));
    }
    // NaN when no k-mer reaches min_cnt, which is ordinary at low depth; see cvttsd2si.
    let kcov = tot_k as f32 / sum_k as f32;
    bfc_opt.min_cov = cvttsd2si(BFC_EC_MIN_COV_COEF * kcov as f64 + 0.499);
    bfc_opt.min_cov = if bfc_opt.min_cov < opt.max_cnt {
        bfc_opt.min_cov
    } else {
        opt.max_cnt
    };
    bfc_opt.min_cov = if bfc_opt.min_cov > opt.min_cnt {
        bfc_opt.min_cov
    } else {
        opt.min_cnt
    };

    if flt_uniq {
        for s in seq.iter_mut() {
            worker_fltuniq(&bfc_opt, &ch, s);
        }
    } else {
        let mut e = Ec1Buf {
            opt: &bfc_opt,
            ch: &ch,
            heap: Vec::new(),
            stack: Vec::new(),
            seq: Vec::new(),
            ec: [Vec::new(), Vec::new()],
            mode,
        };
        for s in seq.iter_mut() {
            bfc_ec1(&mut e, &mut s.seq, &mut s.qual);
        }
    }
    kcov
}

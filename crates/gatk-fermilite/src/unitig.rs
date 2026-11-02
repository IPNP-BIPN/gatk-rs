//! `unitig.c`: unitigs read off the FM-index.
//!
//! Seeds are visited in the order `prime * i % n`, so a contig's index, its orientation and which
//! read wins a tie are all decided by that walk. `kstring_t` becomes a `Vec<u8>` whose length is
//! `s->l`: nothing here reads past `s->l`, the bytes a shrink leaves behind are only ever
//! overwritten.

use crate::khash::Hash64;
use crate::ksort::ks_introsort;
use crate::mag::{mag_g_amend, mag_g_build_hash, Ku128, Mag, Magv};
use crate::rld::{Rld, RldIntv};
use crate::{seq_revcomp6, seq_reverse};

const UTG_PRIMES: [u64; 5] = [123457, 234571, 345679, 456791, 567899];

/// `fm6_comp`.
#[inline]
fn fm6_comp(a: u8) -> u8 {
    if (1..=4).contains(&a) {
        5 - a
    } else {
        a
    }
}

/// `fm6_set_intv`.
#[inline]
fn fm6_set_intv(e: &Rld, c: u8, ik: &mut RldIntv) {
    let cu = c as usize;
    ik.x[0] = e.cnt[cu];
    ik.x[2] = e.cnt[cu + 1].wrapping_sub(e.cnt[cu]);
    ik.x[1] = e.cnt[fm6_comp(c) as usize];
    ik.info = 0;
}

/// `rld_extend0`: only the sentinel's interval.
#[inline]
fn rld_extend0(e: &Rld, ik: &RldIntv, ok0: &mut RldIntv, is_back: bool) {
    let (fw, bw) = (usize::from(!is_back), usize::from(is_back));
    let (mut tk, mut tl) = ([0u64; 6], [0u64; 6]);
    e.rank2a(ik.x[fw], ik.x[fw].wrapping_add(ik.x[2]), &mut tk, &mut tl);
    ok0.x[fw] = tk[0];
    ok0.x[bw] = ik.x[bw];
    ok0.x[2] = tl[0].wrapping_sub(tk[0]);
}

/// `fm6_retrieve`: `k`; the bases go to `s`, the interval to `k2`, the containment to `contained`.
fn fm6_retrieve(e: &Rld, x: u64, s: &mut Vec<u8>, k2: &mut RldIntv, contained: &mut i32) -> u64 {
    let mut k = x;
    let mut ok = [0u64; 6];
    let mut ok2 = [RldIntv::default(); 6];
    s.clear();
    *contained = 0;
    loop {
        let c = e.rank1a(k + 1, &mut ok);
        let cu = c as usize;
        k = e.cnt[cu].wrapping_add(ok[cu]).wrapping_sub(1);
        if c == 0 {
            break;
        }
        if !s.is_empty() {
            if k2.x[2] == 1 {
                k2.x[0] = k;
            } else {
                e.extend(k2, &mut ok2, true);
                *k2 = ok2[cu];
            }
        } else {
            fm6_set_intv(e, c as u8, k2);
        }
        s.push(c as u8);
    }
    if k2.x[2] != 1 {
        e.extend(k2, &mut ok2, true);
        if ok2[0].x[2] != k2.x[2] {
            *contained |= 1; // left contained
        }
        *k2 = ok2[0];
    } else {
        k2.x[0] = k;
    }
    e.extend(k2, &mut ok2, false);
    if ok2[0].x[2] != k2.x[2] {
        *contained |= 2; // right contained
    }
    *k2 = ok2[0];
    k
}

/// `set_bit`.
#[inline]
fn set_bit(bits: &mut [u64], x: u64) {
    bits[(x >> 6) as usize] |= 1u64 << (x & 0x3f);
}

/// `set_bits`.
#[inline]
fn set_bits(bits: &mut [u64], p: &RldIntv) {
    for k in 0..p.x[2] {
        set_bit(bits, p.x[0].wrapping_add(k));
        set_bit(bits, p.x[1].wrapping_add(k));
    }
}

/// `overlap_intv`; `seq[j]` matches the end of a read.
fn overlap_intv(
    e: &Rld,
    seq: &[u8],
    min: i32,
    mut j: i32,
    at5: bool,
    p: &mut Vec<RldIntv>,
    inc_sentinel: bool,
) -> RldIntv {
    let len = seq.len() as i32;
    let mut ik = RldIntv::default();
    let mut ok = [RldIntv::default(); 6];
    p.clear();
    let dir = if at5 { 1 } else { -1 }; // at5 is true iff we start from the 5'-end of a read
    let end = if at5 { len } else { -1 };
    fm6_set_intv(e, seq[j as usize], &mut ik);
    let mut depth = 1;
    j += dir;
    while j != end {
        let c = if at5 {
            fm6_comp(seq[j as usize])
        } else {
            seq[j as usize]
        } as usize;
        e.extend(&ik, &mut ok, !at5);
        if ok[c].x[2] == 0 {
            break; // cannot be extended
        }
        if depth >= min && ok[0].x[2] != 0 {
            if inc_sentinel {
                ok[0].info = (j - dir) as i64 as u64;
                p.push(ok[0]);
            } else {
                ik.info = (j - dir) as i64 as u64;
                p.push(ik);
            }
        }
        ik = ok[c];
        j += dir;
        depth += 1;
    }
    p.reverse(); // such that the smallest interval comes first
    ik
}

/// `aux_t`, without what is never read (`sum`, `sum2`, `n`, `unpaired`).
struct Aux<'a> {
    e: &'a Rld,
    min_match: i32,
    min_merge_len: i32,
    a: [Vec<RldIntv>; 2],
    nei: Vec<RldIntv>,
    cat: Vec<i32>,
    used: Vec<u64>,
    bend: Vec<u64>,
    str_: Vec<u8>,
}

/// `fm6_is_contained`, whose return value has no reader.
fn fm6_is_contained(
    e: &Rld,
    min_match: i32,
    s: &[u8],
    intv: &mut RldIntv,
    ovlp: &mut Vec<RldIntv>,
) {
    let mut ok = [RldIntv::default(); 6];
    c_assert!(
        s.len() as i32 > min_match,
        "fm6_is_contained: s->l > min_match"
    );
    let mut ik = overlap_intv(e, s, min_match, s.len() as i32 - 1, false, ovlp, false);
    e.extend(&ik, &mut ok, true);
    c_assert!(ok[0].x[2] != 0, "fm6_is_contained: ok[0].x[2]");
    ik = ok[0];
    e.extend(&ik, &mut ok, false);
    c_assert!(ok[0].x[2] != 0, "fm6_is_contained: ok[0].x[2]");
    *intv = ok[0];
}

/// `fm6_get_nei`, through `try_right`; `a.a[0]` and `a.a[1]` are `prev` and `curr` on entry, and
/// the swaps inside stay inside.
fn try_right(a: &mut Aux, beg: i32, s: &mut Vec<u8>) -> i32 {
    let e = a.e;
    let min_match = a.min_match;
    let ori_l = s.len() as i32;
    let mut is_forked = false;
    let mut ok = [RldIntv::default(); 6];
    let mut ok0 = RldIntv::default();
    let (pa, pb) = a.a.split_at_mut(1);
    let mut prev: &mut Vec<RldIntv> = &mut pa[0];
    let mut curr: &mut Vec<RldIntv> = &mut pb[0];
    let nei = &mut a.nei;
    let cat = &mut a.cat;
    let used = &mut a.used;

    curr.clear();
    nei.clear();
    if prev.is_empty() {
        // when this routine is called for the seed, prev may be filled by fm6_is_contained()
        overlap_intv(
            e,
            &s[beg as usize..],
            min_match,
            s.len() as i32 - beg - 1,
            false,
            prev,
            false,
        );
        if prev.is_empty() {
            return -1; // no overlap
        }
        for p in prev.iter_mut() {
            p.info = p.info.wrapping_add(beg as i64 as u64);
        }
    }
    cat.clear();
    cat.resize(prev.len(), 0); // only one interval; all point to 0
    while !prev.is_empty() {
        curr.clear();
        for j in 0..prev.len() {
            let p = prev[j];
            if cat[j] < 0 {
                continue;
            }
            e.extend(&p, &mut ok, false); // forward extension
            if ok[0].x[2] != 0 && ori_l != s.len() as i32 {
                // some (partial) reads end here
                rld_extend0(e, &ok[0], &mut ok0, true); // backward extension to look for sentinels
                if ok0.x[2] != 0 {
                    // the match is bounded by sentinels - a full-length match
                    if ok[0].x[2] == p.x[2] && p.x[2] == ok0.x[2] {
                        // never consider a read contained in another read
                        let cat0 = cat[j]; // a category approximately corresponds to one neighbor
                        c_assert!(
                            j == 0 || cat[j] > cat[j - 1],
                            "fm6_get_nei: j == 0 || cat->a[j] > cat->a[j-1]"
                        );
                        ok0.info = (ori_l as i64 as u64).wrapping_sub(p.info & 0xffff_ffff);
                        let mut i = j;
                        while i < prev.len() && cat[i] == cat0 {
                            cat[i] = -1; // mask out other intervals of the same cat
                            i += 1;
                        }
                        nei.push(ok0);
                        continue; // do NOT set "used" as this neighbor may be rejected later
                    } else {
                        set_bits(used, &ok0); // the read is contained in another; mark it as used
                    }
                }
            }
            if cat[j] < 0 {
                continue; // no need to proceed if we have finished this path
            }
            for c in 1..5usize {
                // collect extensible intervals
                if ok[c].x[2] != 0 {
                    rld_extend0(e, &ok[c], &mut ok0, true);
                    if ok0.x[2] != 0 {
                        // do not extend intervals whose left end is not bounded by a sentinel
                        ok[c].info = (p.info & 0xffff_fff0_ffff_ffff) | ((c as u64) << 32);
                        curr.push(ok[c]);
                    }
                }
            }
        }
        if !curr.is_empty() {
            // update category
            cat.resize(curr.len(), 0);
            let c = ((curr[0].info >> 32) & 0xf) as u8;
            s.push(fm6_comp(c));
            ks_introsort(curr, |a, b| a.info < b.info);
            let mut last = (curr[0].info >> 32) as u32;
            cat[0] = 0;
            curr[0].info &= 0xffff_ffff;
            let mut cat0 = 0u32;
            for j in 1..curr.len() {
                // this loop recalculates cat
                if curr[j].info >> 32 != last as u64 {
                    last = (curr[j].info >> 32) as u32;
                    cat0 = j as u32;
                }
                cat[j] = cat0 as i32;
                curr[j].info = (curr[j].info & 0xffff_ffff) | ((cat0 as u64) << 36);
            }
            if cat0 != 0 {
                is_forked = true;
            }
        }
        std::mem::swap(&mut curr, &mut prev);
    }
    if nei.is_empty() {
        return -1; // no overlap
    }
    let rbeg = ori_l.wrapping_sub(nei[0].info as u32 as i32);
    if nei.len() == 1 && is_forked {
        // this may happen if there are contained reads; fix this
        fm6_set_intv(e, 0, &mut ok0);
        let mut i = rbeg;
        while i < ori_l {
            e.extend(&ok0, &mut ok, false);
            ok0 = ok[fm6_comp(s[i as usize]) as usize];
            i += 1;
        }
        let mut i = ori_l;
        while i < s.len() as i32 {
            let mut c0 = -1i32;
            e.extend(&ok0, &mut ok, false);
            let mut j = 0;
            for c in 1..5usize {
                if ok[c].x[2] != 0
                    && ok[c].x[0] <= nei[0].x[0]
                    && ok[c].x[0].wrapping_add(ok[c].x[2]) >= nei[0].x[0].wrapping_add(nei[0].x[2])
                {
                    j += 1;
                    c0 = c as i32;
                }
            }
            if j == 0 && ok[0].x[2] != 0 {
                break;
            }
            c_assert!(j == 1, "fm6_get_nei: j == 1");
            s[i as usize] = fm6_comp(c0 as u8);
            ok0 = ok[c0 as usize];
            i += 1;
        }
        s.truncate(i as usize);
    }
    if nei.len() > 1 {
        s.truncate(ori_l as usize);
    }
    rbeg
}

/// `check_left_simple`.
fn check_left_simple(a: &mut Aux, beg: i32, rbeg: i32, s: &[u8]) -> i32 {
    let mut ok = [RldIntv::default(); 6];
    let (pa, pb) = a.a.split_at_mut(1);
    let mut prev: &mut Vec<RldIntv> = &mut pa[0];
    let mut curr: &mut Vec<RldIntv> = &mut pb[0];
    overlap_intv(a.e, s, a.min_match, rbeg, true, prev, true);
    let mut i = rbeg - 1;
    while i >= beg {
        curr.clear();
        let si = s[i as usize] as usize;
        for j in 0..prev.len() {
            let p = prev[j];
            a.e.extend(&p, &mut ok, true);
            if ok[0].x[2] != 0 {
                set_bits(&mut a.used, &ok[0]); // they must be contained in a longer read
            }
            if ok[0].x[2].wrapping_add(ok[si].x[2]) != p.x[2] {
                return -1; // potential backward bifurcation
            }
            curr.push(ok[si]);
        }
        std::mem::swap(&mut curr, &mut prev);
        i -= 1;
    }
    0
}

/// `check_left`.
fn check_left(a: &mut Aux, beg: i32, rbeg: i32, s: &[u8]) -> i32 {
    c_assert!(a.nei.len() == 1, "check_left: a->nei.n == 1");
    let ret = check_left_simple(a, beg, rbeg, s);
    if ret == 0 {
        return 0;
    }
    // when ret<0, the back fork may be caused by a contained read
    let tmp = a.nei[0]; // try_right() will overwrite it
    a.a[0].clear();
    a.a[1].clear();
    a.nei.clear();
    let mut str_ = std::mem::take(&mut a.str_);
    str_.clear();
    let mut i = s.len() as i32 - 1;
    while i >= rbeg {
        str_.push(fm6_comp(s[i as usize]));
        i -= 1;
    }
    try_right(a, 0, &mut str_);
    a.str_ = str_;
    c_assert!(!a.nei.is_empty(), "check_left: a->nei.n >= 1");
    let ret = if a.nei.len() > 1 { -1 } else { 0 };
    a.nei.clear();
    a.nei.push(tmp); // recover the original neighbour
    ret
}

/// `unitig_unidir`: the number of reads taken in.
fn unitig_unidir(
    a: &mut Aux,
    s: &mut Vec<u8>,
    cov: &mut Vec<u8>,
    beg0: i32,
    k0: u64,
    end: &mut u64,
    is_loop: &mut bool,
) -> i32 {
    let mut beg = beg0;
    let mut ori_l = s.len() as i32;
    let mut n_reads = 0;
    *is_loop = false;
    loop {
        let rbeg = try_right(a, beg, s);
        if rbeg < 0 {
            break; // loop if there is at least one overlap
        }
        if a.nei.len() > 1 {
            // forward bifurcation
            set_bit(&mut a.bend, *end);
            break;
        }
        let k = a.nei[0].x[0];
        if k == *end {
            break; // a loop like b>>c>>a><a; keep the link but stop extension
        }
        if (a.bend[(k >> 6) as usize] >> (k & 0x3f)) & 1 != 0 || check_left(a, beg, rbeg, s) < 0 {
            // backward bifurcation
            set_bit(&mut a.bend, k);
            break;
        }
        if k == k0 {
            // a loop like a>>b>>c>>a
            *is_loop = true;
            break;
        }
        if a.nei[0].x[1] == *end {
            // a loop like b>>c>>a>>a; cut the last link
            a.nei.clear();
            break;
        }
        if (a.nei[0].info as i32) < a.min_merge_len {
            break; // the overlap is not long enough
        }
        *end = a.nei[0].x[1];
        let n0 = a.nei[0];
        set_bits(&mut a.used, &n0); // successful extension
        n_reads += 1;
        cov.resize(s.len(), b'"');
        let mut i = rbeg;
        while i < ori_l {
            // update the coverage string
            let c = &mut cov[i as usize];
            if *c != b'~' {
                *c += 1;
            }
            i += 1;
        }
        for c in &mut cov[ori_l as usize..] {
            *c = b'"';
        }
        beg = rbeg;
        ori_l = s.len() as i32;
        a.a[0].clear(); // prepare for the next round of loop
        a.a[1].clear();
    }
    s.truncate(ori_l as usize);
    cov.truncate(ori_l as usize);
    n_reads
}

/// `copy_nei`.
fn copy_nei(dst: &mut Vec<Ku128>, src: &[RldIntv]) {
    for x in src {
        dst.push(Ku128 {
            x: x.x[0],
            y: x.info,
        });
    }
}

/// `unitig1`.
fn unitig1(
    a: &mut Aux,
    seed: i64,
    s: &mut Vec<u8>,
    cov: &mut Vec<u8>,
    end: &mut [u64; 2],
    nei: &mut [Vec<Ku128>; 2],
    n_reads: &mut i32,
) -> i32 {
    let mut intv0 = RldIntv::default();
    let mut contained = 0;
    let mut is_loop = false;
    *n_reads = 0;
    nei[0].clear();
    nei[1].clear();
    if (a.used[(seed >> 6) as usize] >> (seed & 0x3f)) & 1 != 0 {
        return -2; // used
    }
    // retrieve the sequence pointed by seed
    let k = fm6_retrieve(a.e, seed as u64, s, &mut intv0, &mut contained) as i64;
    seq_reverse(s);
    let seed_len = s.len() as i32;
    // check contained status
    if intv0.x[2] > 1 && k as u64 != intv0.x[0] {
        return -3; // duplicated, but not the first
    }
    set_bits(&mut a.used, &intv0);
    if contained != 0 {
        return -3; // contained
    }
    // check length, containment and if used before
    if s.len() as i32 <= a.min_match {
        return -1; // too short
    }
    {
        let (e, mm) = (a.e, a.min_match);
        fm6_is_contained(e, mm, s, &mut intv0, &mut a.a[0]);
    }
    *n_reads = 1;
    // initialize the coverage string
    cov.clear();
    cov.resize(s.len(), b'"');
    // left-wards extension
    end[0] = intv0.x[1];
    end[1] = intv0.x[0];
    if !a.a[0].is_empty() {
        // no need to extend to the right if there is no overlap
        *n_reads += unitig_unidir(a, s, cov, 0, intv0.x[0], &mut end[0], &mut is_loop);
        copy_nei(&mut nei[0], &a.nei);
        if is_loop {
            nei[1].push(Ku128 {
                x: end[0],
                y: a.nei[0].info,
            });
            return 0;
        }
    }
    // right-wards extension
    a.a[0].clear();
    a.a[1].clear();
    a.nei.clear();
    seq_revcomp6(s); // reverse complement for extension in the other direction
    seq_reverse(cov); // reverse the coverage
    let beg = s.len() as i32 - seed_len;
    *n_reads += unitig_unidir(a, s, cov, beg, intv0.x[1], &mut end[1], &mut is_loop);
    copy_nei(&mut nei[1], &a.nei);
    0
}

/// `fml_fmi2mag_core`, with one thread, through `worker`.
pub(crate) fn fml_fmi2mag_core(e: &Rld, min_match: i32, min_merge_len: i32) -> Mag {
    let mcnt1 = e.mcnt[1];
    let words = (mcnt1.wrapping_add(63) / 64) as usize;
    let mut visited = vec![0u64; words];
    c_assert!(mcnt1 >= 2, "fml_fmi2mag_core: e->mcnt[1] >= n_threads * 2");
    let mut prime = 0u64;
    for &p in &UTG_PRIMES {
        if mcnt1 % p != 0 {
            prime = p;
            break;
        }
    }
    c_assert!(prime != 0, "fml_fmi2mag_core: w.prime");
    let mut a = Aux {
        e,
        min_match,
        min_merge_len,
        a: [Vec::new(), Vec::new()],
        nei: Vec::new(),
        cat: Vec::new(),
        used: vec![0u64; words],
        bend: vec![0u64; words],
        str_: Vec::new(),
    };
    let mut v: Vec<Magv> = Vec::new();
    let mut s = Vec::new();
    let mut cov = Vec::new();
    let mut k = [0u64; 2];
    let mut nei = [Vec::new(), Vec::new()];
    let mut nsr = 0;
    for ii in 0..mcnt1 {
        let i = prime.wrapping_mul(ii) % mcnt1;
        if unitig1(
            &mut a, i as i64, &mut s, &mut cov, &mut k, &mut nei, &mut nsr,
        ) >= 0
        {
            // then we keep the unitig, unless an end has been seen; || stops at the first
            let mut seen = false;
            for &kj in &k {
                let (w, bit) = ((kj >> 6) as usize, 1u64 << (kj & 0x3f));
                let was = visited[w] & bit != 0;
                visited[w] |= bit;
                if was {
                    seen = true;
                    break;
                }
            }
            if seen {
                continue;
            }
            v.push(Magv {
                len: s.len() as i32,
                nsr,
                k,
                nei: nei.clone(),
                seq: s.clone(),
                cov: cov.clone(),
                ptr: -1,
            });
        }
    }
    let mut g = Mag {
        v,
        h: Hash64::new(),
        min_ovlp: 0,
    };
    mag_g_build_hash(&mut g);
    mag_g_amend(&mut g);
    // g->rdist = mag_cal_rdist(g): only verbose messages read rdist, and it changes nothing.
    g
}

//! `mag.c`: the overlap graph and its cleaning.
//!
//! A vertex's edges are read through its index on every use, as C re-reads `r->n` through the
//! pointer: an edge list can change under the loop walking it when a vertex neighbours itself.

use crate::bubble::{mag_g_pop_open, mag_g_pop_simple, mag_g_simplify_bubble};
use crate::khash::Hash64;
use crate::ksort::ks_introsort;
use crate::{seq_revcomp6, seq_reverse, MagOpt, MAG_F_AGGRESSIVE, MAG_F_NO_SIMPL, MAG_F_POPOPEN};

/// `ku128_t`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct Ku128 {
    pub(crate) x: u64,
    pub(crate) y: u64,
}

/// `magv_t`, `max_len` left out: `seq` and `cov` have the length `len`. `ptr` is an index into the
/// pool of `mag_vh_simplify_bubble`, -1 for NULL.
#[derive(Clone, Default)]
pub(crate) struct Magv {
    pub(crate) len: i32,
    pub(crate) nsr: i32,
    pub(crate) k: [u64; 2],
    pub(crate) nei: [Vec<Ku128>; 2],
    pub(crate) seq: Vec<u8>,
    pub(crate) cov: Vec<u8>,
    pub(crate) ptr: i32,
}

/// `mag_t`, `rdist` left out: no reader.
pub(crate) struct Mag {
    pub(crate) v: Vec<Magv>,
    pub(crate) h: Hash64,
    pub(crate) min_ovlp: i32,
}

/// `edge_mark_del`.
#[inline]
pub(crate) fn edge_mark_del(e: &mut Ku128) {
    e.x = -2i64 as u64;
    e.y = 0;
}

/// `edge_is_del`.
#[inline]
pub(crate) fn edge_is_del(e: &Ku128) -> bool {
    e.x == -2i64 as u64 || e.y == 0
}

/// `ku128_xlt`.
#[inline]
fn ku128_xlt(a: &Ku128, b: &Ku128) -> bool {
    a.x < b.x || (a.x == b.x && a.y > b.y)
}

/// `v128_clean`.
pub(crate) fn v128_clean(r: &mut Vec<Ku128>) {
    let mut j = 0;
    for i in 0..r.len() {
        if !edge_is_del(&r[i]) {
            if j != i {
                r[j] = r[i];
            }
            j += 1;
        }
    }
    r.truncate(j);
}

/// `v128_rmdup`.
fn v128_rmdup(r: &mut Vec<Ku128>) {
    if r.len() > 1 {
        ks_introsort(r, ku128_xlt);
    }
    let mut l = 0;
    let mut cnt = 0;
    while l < r.len() {
        // jump to the first node to be retained
        if edge_is_del(&r[l]) {
            cnt += 1;
        } else {
            break;
        }
        l += 1;
    }
    if l == r.len() {
        r.clear(); // no good arcs
        return;
    }
    let mut x = r[l].x;
    l += 1;
    while l < r.len() {
        // mark duplicated node
        if edge_is_del(&r[l]) || r[l].x == x {
            edge_mark_del(&mut r[l]);
            cnt += 1;
        } else {
            x = r[l].x;
        }
        l += 1;
    }
    if cnt != 0 {
        v128_clean(r);
    }
}

/// `mag_g_build_hash`. A repeated terminal would map to -1, but `worker` in `unitig.c` keeps one
/// end from being added twice, even for a vertex whose two ends are the same.
pub(crate) fn mag_g_build_hash(g: &mut Mag) {
    for i in 0..g.v.len() {
        for j in 0..2 {
            let (k, ret) = g.h.put(g.v[i].k[j]);
            g.h.vals[k as usize] = if ret == 0 {
                -1i64 as u64
            } else {
                ((i as u64) << 1) | j as u64
            };
        }
    }
}

/// `tid2idd`.
#[inline]
pub(crate) fn tid2idd(h: &Hash64, tid: u64) -> u64 {
    let k = h.get(tid);
    c_assert!(k != h.end(), "tid2idd: k != kh_end(h)");
    h.vals[k as usize]
}

/// `mag_g_amend`.
pub(crate) fn mag_g_amend(g: &mut Mag) {
    for i in 0..g.v.len() {
        for j in 0..2 {
            let mut l = 0;
            while l < g.v[i].nei[j].len() {
                let x = g.v[i].nei[j][l].x;
                let k = g.h.get(x);
                if k == g.h.end() {
                    // neighbor is not in the hash table; likely due to tip removal
                    edge_mark_del(&mut g.v[i].nei[j][l]);
                    l += 1;
                    continue;
                }
                let z = g.h.vals[k as usize];
                let pk = g.v[i].k[j];
                let found = g.v[(z >> 1) as usize].nei[(z & 1) as usize]
                    .iter()
                    .any(|e| e.x == pk);
                if !found {
                    // not in neighbor's neighor
                    edge_mark_del(&mut g.v[i].nei[j][l]);
                }
                l += 1;
            }
            v128_rmdup(&mut g.v[i].nei[j]);
        }
    }
}

/// `mag_v_destroy`.
fn mag_v_destroy(v: &mut Magv) {
    *v = Magv {
        len: -1,
        ptr: -1,
        ..Magv::default()
    };
}

/// `mag_eh_add`: add `v` to `u`.
fn mag_eh_add(g: &mut Mag, u: u64, v: u64, ovlp: i32) {
    if (u as i64) < 0 {
        return;
    }
    let idd = tid2idd(&g.h, u);
    let r = &mut g.v[(idd >> 1) as usize].nei[(idd & 1) as usize];
    if r.iter().any(|e| e.x == v) {
        return; // no multi-edges
    }
    r.push(Ku128 {
        x: v,
        y: ovlp as i64 as u64,
    });
}

/// `mag_eh_markdel`: mark deletion of `v` from `u`.
pub(crate) fn mag_eh_markdel(g: &mut Mag, u: u64, v: u64) {
    if (u as i64) < 0 {
        return;
    }
    let idd = tid2idd(&g.h, u);
    for e in &mut g.v[(idd >> 1) as usize].nei[(idd & 1) as usize] {
        if e.x == v {
            edge_mark_del(e);
        }
    }
}

/// `mag_v_del`.
pub(crate) fn mag_v_del(g: &mut Mag, p: usize) {
    if g.v[p].len < 0 {
        return;
    }
    let (k0, k1) = (g.v[p].k[0], g.v[p].k[1]);
    for i in 0..2 {
        let mut j = 0;
        while j < g.v[p].nei[i].len() {
            let e = g.v[p].nei[i][j];
            if !edge_is_del(&e) && e.x != k0 && e.x != k1 {
                mag_eh_markdel(g, e.x, g.v[p].k[i]);
            }
            j += 1;
        }
    }
    for i in 0..2 {
        let k = g.h.get(g.v[p].k[i]);
        g.h.del(k);
    }
    mag_v_destroy(&mut g.v[p]);
}

/// `mag_v_transdel`.
fn mag_v_transdel(g: &mut Mag, p: usize, min_ovlp: i32) {
    if !g.v[p].nei[0].is_empty() && !g.v[p].nei[1].is_empty() {
        let (k0, k1) = (g.v[p].k[0], g.v[p].k[1]);
        let mut i = 0;
        while i < g.v[p].nei[0].len() {
            let a = g.v[p].nei[0][i];
            if !(edge_is_del(&a) || a.x == k0 || a.x == k1) {
                let mut j = 0;
                while j < g.v[p].nei[1].len() {
                    let b = g.v[p].nei[1][j];
                    if !(edge_is_del(&b) || b.x == k0 || b.x == k1) {
                        let ovlp = (a.y.wrapping_add(b.y) as i32).wrapping_sub(g.v[p].len);
                        if ovlp >= min_ovlp {
                            mag_eh_add(g, a.x, b.x, ovlp);
                            mag_eh_add(g, b.x, a.x, ovlp);
                        }
                    }
                    j += 1;
                }
            }
            i += 1;
        }
    }
    mag_v_del(g, p);
}

/// `mag_v_flip`.
fn mag_v_flip(g: &mut Mag, p: usize) {
    let v = &mut g.v[p];
    seq_revcomp6(&mut v.seq);
    seq_reverse(&mut v.cov);
    v.k.swap(0, 1);
    v.nei.swap(0, 1);
    let (k0, k1) = (v.k[0], v.k[1]);
    let k = g.h.get(k0);
    c_assert!(k != g.h.end(), "mag_v_flip: k != kh_end(h)");
    g.h.vals[k as usize] ^= 1;
    let k = g.h.get(k1);
    c_assert!(k != g.h.end(), "mag_v_flip: k != kh_end(h)");
    g.h.vals[k as usize] ^= 1;
}

/// `mag_vh_merge_try`: merge p's neighbor to the right end of p, 0 if it did.
fn mag_vh_merge_try(g: &mut Mag, p: usize, min_merge_len: i32) -> i32 {
    // check if an unambiguous merge can be performed
    if g.v[p].nei[1].len() != 1 {
        return -1; // multiple or no neighbor; do not merge
    }
    let e = g.v[p].nei[1][0];
    if (e.x as i64) < 0 {
        return -2; // deleted neighbor
    }
    if (e.y as i32) < min_merge_len {
        return -5;
    }
    let kq = g.h.get(e.x);
    c_assert!(kq != g.h.end(), "mag_vh_merge_try: kq != kh_end(h)");
    let qv = g.h.vals[kq as usize];
    let q = (qv >> 1) as usize;
    if p == q {
        return -3; // we have a loop p->p. We cannot merge in this case
    }
    if g.v[q].nei[(qv & 1) as usize].len() != 1 {
        return -4; // the neighbor q has multiple neighbors
    }
    // we can perform a merge; do further consistency check (mostly check bugs)
    if qv & 1 != 0 {
        mag_v_flip(g, q); // a "><" bidirectional arc; flip q
    }
    let kp = g.h.get(g.v[p].k[1]);
    c_assert!(kp != g.h.end(), "mag_vh_merge_try: kp != kh_end(h)");
    g.h.del(kp);
    g.h.del(kq); // remove the two ends of the arc in the hash table
    c_assert!(
        g.v[p].k[1] == g.v[q].nei[0][0].x && g.v[q].k[0] == e.x,
        "mag_vh_merge_try: p->k[1] == q->nei[0].a[0].x && q->k[0] == p->nei[1].a[0].x"
    );
    c_assert!(
        e.y == g.v[q].nei[0][0].y,
        "mag_vh_merge_try: same overlap length"
    );
    c_assert!(
        g.v[p].len as i64 as u64 >= e.y && g.v[q].len as i64 as u64 >= e.y,
        "mag_vh_merge_try: p->len >= y && q->len >= y"
    );
    // update the read count and sequence length; q goes away below
    let mut qt = std::mem::take(&mut g.v[q]);
    let pv = &mut g.v[p];
    pv.nsr = pv.nsr.wrapping_add(qt.nsr);
    let y = e.y as i32;
    let new_l = pv.len + qt.len - y;
    pv.seq.resize(new_l as usize, 0);
    pv.cov.resize(new_l as usize, 0);
    // merge seq and cov
    let mut i = pv.len - y;
    for j in 0..qt.len as usize {
        let iu = i as usize;
        pv.seq[iu] = qt.seq[j];
        if i < pv.len {
            let sum = pv.cov[iu] as i8 as i32 + (qt.cov[j] as i8 as i32 - 33);
            if sum > 126 {
                pv.cov[iu] = 126;
            } else {
                pv.cov[iu] =
                    (pv.cov[iu] as i8).wrapping_add((qt.cov[j] as i8).wrapping_sub(33)) as u8;
            }
        } else {
            pv.cov[iu] = qt.cov[j];
        }
        i += 1;
    }
    pv.len = new_l;
    // merge neighbors
    pv.nei[1] = std::mem::take(&mut qt.nei[1]);
    pv.k[1] = qt.k[1];
    // update the hash table for the right end of p
    let k1 = pv.k[1];
    let kp = g.h.get(k1);
    c_assert!(kp != g.h.end(), "mag_vh_merge_try: kp != kh_end(h)");
    g.h.vals[kp as usize] = ((p as u64) << 1) | 1;
    // clean up q
    mag_v_destroy(&mut g.v[q]);
    0
}

/// `mag_g_merge`. Every live vertex comes out flipped, merged or not.
pub(crate) fn mag_g_merge(g: &mut Mag, rmdup: bool, min_merge_len: i32) {
    for v in &mut g.v {
        // remove multiedges
        if rmdup {
            v128_rmdup(&mut v.nei[0]);
            v128_rmdup(&mut v.nei[1]);
        } else {
            v128_clean(&mut v.nei[0]);
            v128_clean(&mut v.nei[1]);
        }
    }
    for i in 0..g.v.len() {
        if g.v[i].len < 0 {
            continue;
        }
        while mag_vh_merge_try(g, i, min_merge_len) == 0 {}
        mag_v_flip(g, i);
        while mag_vh_merge_try(g, i, min_merge_len) == 0 {}
    }
}

/// The `vlt1` introsort over vertex indices, `mag_vlt1`.
fn sort_vlt1(g: &Mag, a: &mut [u32]) {
    ks_introsort(a, |&x, &y| {
        let (p, q) = (&g.v[x as usize], &g.v[y as usize]);
        p.nsr < q.nsr || (p.nsr == q.nsr && p.len < q.len)
    });
}

/// `mag_g_rm_vext`.
fn mag_g_rm_vext(g: &mut Mag, min_len: i32, min_nsr: i32) {
    let mut a: Vec<u32> = Vec::new();
    for (i, p) in g.v.iter().enumerate() {
        if p.len < 0 || (!p.nei[0].is_empty() && !p.nei[1].is_empty()) {
            continue;
        }
        if p.len >= min_len || p.nsr >= min_nsr {
            continue;
        }
        a.push(i as u32);
    }
    sort_vlt1(g, &mut a);
    for &i in &a {
        mag_v_del(g, i as usize);
    }
}

/// `mag_g_rm_vint`.
fn mag_g_rm_vint(g: &mut Mag, min_len: i32, min_nsr: i32, min_ovlp: i32) {
    let mut a: Vec<u32> = Vec::new();
    for (i, p) in g.v.iter().enumerate() {
        if p.len >= 0 && p.len < min_len && p.nsr < min_nsr {
            a.push(i as u32);
        }
    }
    sort_vlt1(g, &mut a);
    for &i in &a {
        mag_v_transdel(g, i as usize, min_ovlp);
    }
}

/// `mag_g_rm_edge`.
fn mag_g_rm_edge(g: &mut Mag, min_ovlp: i32, min_ratio: f64, min_len: i32, min_nsr: i32) {
    let mut a: Vec<u32> = Vec::new();
    for (i, p) in g.v.iter().enumerate() {
        if p.len < 0 {
            continue;
        }
        if (p.nei[0].is_empty() || p.nei[1].is_empty()) && p.len < min_len && p.nsr < min_nsr {
            continue; // skip tips
        }
        a.push(i as u32);
    }
    sort_vlt1(g, &mut a);
    for &pi in a.iter().rev() {
        let p = pi as usize;
        for j in 0..2 {
            if g.v[p].nei[j].is_empty() {
                continue; // no overlapping reads
            }
            let mut max_ovlp = min_ovlp;
            let mut max_k = -1i32;
            for (k, e) in g.v[p].nei[j].iter().enumerate() {
                // get the max overlap length
                if (max_ovlp as i64 as u64) < e.y {
                    max_ovlp = e.y as i32;
                    max_k = k as i32;
                }
            }
            if max_k >= 0 {
                // test if max_k is a tip
                let x = tid2idd(&g.h, g.v[p].nei[j][max_k as usize].x);
                let q = &g.v[(x >> 1) as usize];
                if q.len >= 0
                    && (q.nei[0].is_empty() || q.nei[1].is_empty())
                    && q.len < min_len
                    && q.nsr < min_nsr
                {
                    max_ovlp = min_ovlp;
                }
            }
            let mut k = 0;
            while k < g.v[p].nei[j].len() {
                let e = g.v[p].nei[j][k];
                if !edge_is_del(&e)
                    && (e.y < min_ovlp as i64 as u64
                        || (e.y as f64) / (max_ovlp as f64) < min_ratio)
                {
                    let pk = g.v[p].k[j];
                    mag_eh_markdel(g, e.x, pk); // FIXME upstream: should we check if r->a[k] is p?
                    edge_mark_del(&mut g.v[p].nei[j][k]);
                }
                k += 1;
            }
        }
    }
}

/// `mag_g_clean`.
pub(crate) fn mag_g_clean(g: &mut Mag, opt: &MagOpt) {
    if g.min_ovlp < opt.min_ovlp {
        g.min_ovlp = opt.min_ovlp;
    }
    for j in 2..=opt.min_ensr {
        mag_g_rm_vext(g, opt.min_elen, j);
    }
    mag_g_merge(g, false, opt.min_merge_len);
    let m = g.min_ovlp;
    mag_g_rm_edge(g, m, opt.min_dratio1 as f64, opt.min_elen, opt.min_ensr);
    mag_g_merge(g, true, opt.min_merge_len);
    for j in 2..=opt.min_ensr {
        mag_g_rm_vext(g, opt.min_elen, j);
    }
    mag_g_merge(g, false, opt.min_merge_len);
    if (opt.flag & MAG_F_AGGRESSIVE) != 0 || (opt.flag & MAG_F_POPOPEN) != 0 {
        mag_g_pop_open(g, opt.min_elen);
    }
    if opt.flag & MAG_F_NO_SIMPL == 0 {
        mag_g_simplify_bubble(g, opt.max_bvtx, opt.max_bdist);
    }
    mag_g_pop_simple(
        g,
        opt.max_bcov,
        opt.max_bfrac,
        opt.min_merge_len,
        opt.max_bdiff,
        opt.flag & MAG_F_AGGRESSIVE != 0,
    );
    let m = g.min_ovlp;
    mag_g_rm_vint(g, opt.min_elen, opt.min_insr, m);
    mag_g_rm_edge(g, m, opt.min_dratio1 as f64, opt.min_elen, opt.min_ensr);
    mag_g_merge(g, true, opt.min_merge_len);
    mag_g_rm_vext(g, opt.min_elen, opt.min_ensr);
    mag_g_merge(g, false, opt.min_merge_len);
    if (opt.flag & MAG_F_AGGRESSIVE) != 0 || (opt.flag & MAG_F_POPOPEN) != 0 {
        mag_g_pop_open(g, opt.min_elen);
    }
    mag_g_rm_vext(g, opt.min_elen, opt.min_ensr);
    mag_g_merge(g, false, opt.min_merge_len);
}

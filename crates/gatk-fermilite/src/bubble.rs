//! `bubble.c`: closed bubbles (`mag_g_simplify_bubble`, `mag_g_pop_simple`) and open ones
//! (`mag_g_pop_open`). Which of these runs is decided by `mag_opt.flag`: the default flag skips the
//! first, and GATK's `setCleaningFlag(0x60)` clears `MAG_F_NO_SIMPL` and so runs it.

use crate::khash::Hash64;
use crate::ksw::{ksw_i16, ksw_qinit, Kswq};
use crate::mag::{edge_is_del, edge_mark_del, mag_g_merge, mag_v_del, tid2idd, v128_clean, Mag};
use crate::seq_revcomp6;

const MAX_N_DIFF: f64 = 2.01; // for evaluating alignment after SW
const MAX_R_DIFF: f64 = 0.1;
const L_DIFF_COEF: f64 = 0.2; // n_diff=|l_0 - l_1|*L_DIFF_COEF

/// `trinfo_t`, starting from `g_trinull`.
#[derive(Clone, Copy)]
struct TrInfo {
    id: u64,
    cnt: [i32; 2],
    n: [[i32; 2]; 2],
    d: [[i32; 2]; 2],
    v: [[u64; 2]; 2],
}

/// `tip_alloc`, on `g_trinull`.
fn tip_alloc(pool: &mut Vec<TrInfo>, id: u32) -> i32 {
    pool.push(TrInfo {
        id: id as u64,
        cnt: [0, 0],
        n: [[i32::MIN; 2]; 2],
        d: [[i32::MIN; 2]; 2],
        v: [[u64::MAX; 2]; 2],
    });
    pool.len() as i32 - 1
}

/// `backtrace`.
fn backtrace(g: &Mag, pool: &[TrInfo], mut end: u64, start: u64, h: &mut Hash64) {
    while end >> 32 != start {
        h.put(end >> 33);
        let p = g.v[(end >> 33) as usize].ptr as usize;
        end = pool[p].v[(((end >> 32) ^ 1) & 1) as usize][(end & 1) as usize];
    }
}

/// `mag_vh_simplify_bubble`.
fn mag_vh_simplify_bubble(
    g: &mut Mag,
    idd: u64,
    max_vtx: i32,
    max_dist: i32,
    pool: &mut Vec<TrInfo>,
    stack: &mut Vec<u64>,
) {
    let p0 = (idd >> 1) as usize;
    if g.v[p0].len < 0 || g.v[p0].nei[(idd & 1) as usize].len() < 2 {
        return; // stop if p is deleted or it has 0 or 1 neighbor
    }
    // reset aux data
    stack.clear();
    pool.clear();
    let mut h = Hash64::new();
    let mut n_pending = 0i32;
    // add the initial vertex
    let t = tip_alloc(pool, (idd >> 1) as u32);
    g.v[p0].ptr = t;
    pool[t as usize].d[((idd & 1) ^ 1) as usize][0] = g.v[p0].len.wrapping_neg();
    pool[t as usize].n[((idd & 1) ^ 1) as usize][0] = g.v[p0].nsr.wrapping_neg();
    stack.push(idd ^ 1);
    // essentially a topological sorting
    while !stack.is_empty() {
        if stack.len() == 1 && stack[0] != (idd ^ 1) && n_pending == 0 {
            break; // found the other end of the bubble
        }
        let x = stack.pop().unwrap();
        let xu = (x & 1) as usize;
        let p = (x >> 1) as usize;
        let rd = xu ^ 1; // we will look at the neighbors from the other end of the unitig
        let tp = g.v[p].ptr as usize;
        if pool.len() as i32 > max_vtx
            || pool[tp].d[xu][0] > max_dist
            || pool[tp].d[xu][1] > max_dist
            || g.v[p].nei[rd].is_empty()
        {
            break; // we failed
        }
        // set the distance to p's neighbors
        let mut i = 0usize;
        while i < g.v[p].nei[rd].len() {
            let e = g.v[p].nei[rd][i];
            if (e.x as i64) < 0 {
                i += 1;
                continue;
            }
            let y = tid2idd(&g.h, e.x);
            if y == (idd ^ 1) {
                // there is a loop involving the initial vertex
                stack.clear();
                break; // not a bubble; this will jump out of the while() loop
            }
            let q = (y >> 1) as usize;
            let yu = (y & 1) as usize;
            if g.v[q].ptr < 0 {
                // has not been attempted
                g.v[q].ptr = tip_alloc(pool, (y >> 1) as u32);
                n_pending += 1;
                v128_clean(&mut g.v[q].nei[yu]); // make sure there are no deleted edges
            }
            let tq = g.v[q].ptr as usize;
            let tp = g.v[p].ptr as usize;
            let (plen, pnsr) = (g.v[p].len, g.v[p].nsr);
            let (mut nsr, mut dist, mut which);
            nsr = pool[tp].n[xu][0].wrapping_add(pnsr);
            which = 0u64;
            dist = pool[tp].d[xu][0]
                .wrapping_add(plen)
                .wrapping_sub(e.y as i32);
            // test and possibly update the tentative distance
            let link = |which: u64| ((x ^ 1) << 32) | ((i as u64) << 1) | which;
            if nsr > pool[tq].n[yu][0] {
                // then move the best to the 2nd best and update the best
                pool[tq].n[yu][1] = pool[tq].n[yu][0];
                pool[tq].n[yu][0] = nsr;
                pool[tq].v[yu][1] = pool[tq].v[yu][0];
                pool[tq].v[yu][0] = link(which);
                pool[tq].d[yu][1] = pool[tq].d[yu][0];
                pool[tq].d[yu][0] = dist;
                nsr = pool[tp].n[xu][1].wrapping_add(pnsr); // now nsr is the 2nd best
                which = 1;
                dist = pool[tp].d[xu][1]
                    .wrapping_add(plen)
                    .wrapping_sub(e.y as i32);
            }
            if nsr > pool[tq].n[yu][1] {
                // update the 2nd best
                pool[tq].n[yu][1] = nsr;
                pool[tq].v[yu][1] = link(which);
                pool[tq].d[yu][1] = dist;
            }
            pool[tq].cnt[yu] += 1;
            if pool[tq].cnt[yu] as usize == g.v[q].nei[yu].len() {
                // all q's predecessors have been processed; then push
                stack.push(y);
                n_pending -= 1;
            }
            i += 1;
        }
    }
    if n_pending == 0 && stack.len() == 1 {
        // found a bubble
        let x = stack[0];
        let tp = g.v[(x >> 1) as usize].ptr as usize;
        let (a, b) = (
            pool[tp].v[(x & 1) as usize][0],
            pool[tp].v[(x & 1) as usize][1],
        );
        backtrace(g, pool, a, idd, &mut h);
        backtrace(g, pool, b, idd, &mut h);
    }
    for t in pool.iter() {
        // reset p->ptr
        g.v[t.id as usize].ptr = -1;
    }
    if h.size() != 0 {
        // bubble detected; then remove verticies not in the top two paths
        for t in pool.iter().skip(1) {
            // the first one is the initial vertex which we want to exclude
            let id = t.id;
            if id != stack[0] >> 1 && h.get(id) == h.end() {
                mag_v_del(g, id as usize);
            }
        }
    }
}

/// `mag_g_simplify_bubble`.
pub(crate) fn mag_g_simplify_bubble(g: &mut Mag, max_vtx: i32, max_dist: i32) {
    let mut pool = Vec::new();
    let mut stack = Vec::new();
    for i in 0..g.v.len() as u64 {
        mag_vh_simplify_bubble(g, i << 1, max_vtx, max_dist, &mut pool, &mut stack);
        mag_vh_simplify_bubble(g, (i << 1) | 1, max_vtx, max_dist, &mut pool, &mut stack);
    }
    mag_g_merge(g, false, 0);
}

/// `abs(l[0] - l[1])`, whose two sides are `int`.
#[inline]
fn c_abs(x: i32) -> i32 {
    x.wrapping_abs()
}

/// One step of `mag_vh_pop_simple`: `q->cov[i]`, where the terminating NUL at `len` can be read.
#[inline]
fn cov_at(cov: &[u8], i: i32) -> f32 {
    let b = if i as usize == cov.len() {
        0
    } else {
        cov[i as usize]
    };
    (b as i8 as i32 - 33) as f32
}

/// `mag_vh_pop_simple`.
fn mag_vh_pop_simple(
    g: &mut Mag,
    idd: u64,
    max_cov: f32,
    max_frac: f32,
    max_bdiff: i32,
    aggressive: bool,
) {
    let p = (idd >> 1) as usize;
    let side = (idd & 1) as usize;
    let max_n_diff: f32 = if aggressive {
        (MAX_N_DIFF * 2.) as f32
    } else {
        MAX_N_DIFF as f32
    };
    if g.v[p].len < 0 || g.v[p].nei[side].len() != 2 {
        return; // deleted or no bubble
    }
    let mut dir = [0usize; 2];
    let mut q = [0usize; 2];
    let mut l = [0i32; 2];
    for j in 0..2 {
        let e = g.v[p].nei[side][j];
        if (e.x as i64) < 0 {
            return;
        }
        let x = tid2idd(&g.h, e.x);
        dir[j] = (x & 1) as usize;
        q[j] = (x >> 1) as usize;
        let qv = &g.v[q[j]];
        if qv.nei[0].len() != 1 || qv.nei[1].len() != 1 {
            return; // no bubble
        }
        l[j] = qv
            .len
            .wrapping_sub(qv.nei[0][0].y.wrapping_add(qv.nei[1][0].y) as i32);
    }
    if g.v[q[0]].nei[dir[0] ^ 1][0].x != g.v[q[1]].nei[dir[1] ^ 1][0].x {
        return; // no bubble
    }
    if l[0] - l[1] > max_bdiff || l[1] - l[0] > max_bdiff {
        return; // huge bubble differences
    }
    let mut seq: [Vec<u8>; 2] = [Vec::new(), Vec::new()];
    let mut avg = [0f32; 2];
    for j in 0..2 {
        let qv = &g.v[q[j]];
        let y0 = qv.nei[0][0].y as i32;
        if l[j] > 0 {
            let lj = l[j] as usize;
            let mut s = qv.seq[y0 as usize..y0 as usize + lj].to_vec();
            let mut c = qv.cov[y0 as usize..y0 as usize + lj].to_vec();
            if dir[j] != 0 {
                seq_revcomp6(&mut s);
                c.reverse();
            }
            let mut a = 0f32;
            for i in 0..lj {
                s[i] = s[i].wrapping_sub(1); // change DNA6 encoding to DNA4 for SW below
                a += (c[i] as i8 as i32 - 33) as f32;
            }
            avg[j] = a / l[j] as f32;
            seq[j] = s;
        } else {
            // l[j] <= 0; this may happen around a tandem repeat
            let mut beg = y0;
            let mut end = qv.len.wrapping_sub(qv.nei[1][0].y as i32);
            if beg > end {
                std::mem::swap(&mut beg, &mut end);
            }
            if beg < end {
                let mut a = 0f32;
                for i in beg..end {
                    a += cov_at(&qv.cov, i);
                }
                avg[j] = a / (end - beg) as f32;
            } else {
                avg[j] = cov_at(&qv.cov, beg); // FIXME upstream: when q[j] is contained...
            }
        }
    }
    let n_diff: f32;
    let r_diff: f32;
    if l[0] > 0 && l[1] > 0 {
        // then do SW to compute n_diff and r_diff
        let mut mat = [0i8; 16];
        let mut k = 0;
        for i in 0..4 {
            for j in 0..4 {
                mat[k] = if i == j { 5 } else { -4 };
                k += 1;
            }
        }
        let qry = ksw_qinit(&seq[0], &mat);
        let score = ksw_i16(&qry, &seq[1], 5, 2);
        // 5: matching score; -4: mismatchig score
        n_diff = ((l[0].min(l[1]) as f64 * 5. - score as f64) / (5. + 4.)) as f32;
        r_diff = (n_diff as f64 / ((l[0] + l[1]) as f64 / 2.)) as f32;
    } else {
        n_diff = (c_abs(l[0] - l[1]) as f64 * L_DIFF_COEF) as f32;
        r_diff = 1.;
    }
    if n_diff < max_n_diff || (r_diff as f64) < MAX_R_DIFF {
        let j = if avg[0] < avg[1] { 0 } else { 1 };
        if aggressive || (avg[j] < max_cov && avg[j] / (avg[j ^ 1] + avg[j]) < max_frac) {
            mag_v_del(g, q[j]);
        }
    }
}

/// `mag_g_pop_simple`.
pub(crate) fn mag_g_pop_simple(
    g: &mut Mag,
    max_cov: f32,
    max_frac: f32,
    min_merge_len: i32,
    max_bdiff: i32,
    aggressive: bool,
) {
    for i in 0..g.v.len() as u64 {
        mag_vh_pop_simple(g, i << 1, max_cov, max_frac, max_bdiff, aggressive);
        mag_vh_pop_simple(g, (i << 1) | 1, max_cov, max_frac, max_bdiff, aggressive);
    }
    mag_g_merge(g, false, min_merge_len);
}

/// `mag_v_pop_open`.
fn mag_v_pop_open(g: &mut Mag, p: usize, min_elen: i32) {
    if g.v[p].len < 0 || g.v[p].len >= min_elen {
        return;
    }
    if g.v[p].nei[0].len() + g.v[p].nei[1].len() != 1 {
        return;
    }
    let dir = if !g.v[p].nei[0].is_empty() { 0 } else { 1 };
    // initialize the scoring system
    let mut mat = [0i8; 16];
    let mut m = 0;
    for i in 0..4 {
        for j in 0..4 {
            mat[m] = if i == j { 5 } else { -4 };
            m += 1;
        }
    }
    let mut l = 0;
    while l < g.v[p].nei[dir].len() {
        let sl = g.v[p].nei[dir][l];
        if (sl.x as i64) < 0 {
            l += 1;
            continue;
        }
        let v = tid2idd(&g.h, sl.x);
        let q = (v >> 1) as usize;
        let vu = (v & 1) as usize;
        if q == p || g.v[q].nei[vu].len() == 1 {
            l += 1;
            continue;
        }
        // get the query ready
        let plen = g.v[p].len;
        let y = sl.y as i32;
        let max_l = (plen - y).wrapping_mul(2);
        let mut query = Vec::with_capacity((plen - y).max(0) as usize);
        if dir == 0 {
            // forward strand
            for j in y..plen {
                query.push(g.v[p].seq[j as usize].wrapping_sub(1));
            }
        } else {
            // reverse
            let mut j = plen - y - 1;
            while j >= 0 {
                query.push(4u8.wrapping_sub(g.v[p].seq[j as usize]));
                j -= 1;
            }
        }
        let l_qry = query.len() as i32;
        let qry: Kswq = ksw_qinit(&query, &mat);
        let pk = g.v[p].k[dir];
        let mut target: Vec<u8> = Vec::new();
        let mut hit = false;
        for i in 0..g.v[q].nei[vu].len() {
            let ri = g.v[q].nei[vu][i];
            if ri.x == pk || (ri.x as i64) < 0 {
                continue;
            }
            let w = tid2idd(&g.h, ri.x);
            // get the target sequence
            let t = &g.v[(w >> 1) as usize];
            target.clear();
            let ry = ri.y as i32;
            if w & 1 != 0 {
                // reverse strand
                let mut j = t.len - ry - 1;
                while j >= 0 && (target.len() as i32) < max_l {
                    target.push(4u8.wrapping_sub(t.seq[j as usize]));
                    j -= 1;
                }
            } else {
                let mut j = ry;
                while j < t.len && (target.len() as i32) < max_l {
                    target.push(t.seq[j as usize].wrapping_sub(1));
                    j += 1;
                }
            }
            let score = ksw_i16(&qry, &target, 5, 2);
            if score >= l_qry * 5 / 2 {
                // 5: matching score; -4: mismatchig score
                let n_diff = (l_qry as f64 * 5. - score as f64) / (5. + 4.);
                let r_diff = n_diff / l_qry as f64;
                if n_diff < MAX_N_DIFF || r_diff < MAX_R_DIFF {
                    hit = true;
                    break;
                }
            }
        }
        if hit {
            // mark delete in p and delete in q
            edge_mark_del(&mut g.v[p].nei[dir][l]);
            for e in &mut g.v[q].nei[vu] {
                if e.x == pk {
                    edge_mark_del(e);
                }
            }
        }
        l += 1;
    }
    if g.v[p].nei[dir].iter().all(edge_is_del) {
        mag_v_del(g, p); // p is not connected to any other vertices
    }
}

/// `mag_g_pop_open`.
pub(crate) fn mag_g_pop_open(g: &mut Mag, min_elen: i32) {
    for i in 0..g.v.len() {
        mag_v_pop_open(g, i, min_elen);
    }
    mag_g_merge(g, false, 0);
}

//! `ksw.c`: the striped Smith-Waterman, as the 16-bit lanes it computes on.
//!
//! fermi-lite calls `ksw_align` with `xtra` 0, so the query is always 16-bit (`ksw_u8` is never
//! reached), no second pass runs, and only `.score` is read. Those lanes saturate and subtract
//! unsigned (`_mm_adds_epi16` over `_mm_subs_epu16`), and the lazy-F loop stops after 16 rounds, so
//! this is not an ordinary Smith-Waterman; each SSE2 intrinsic is spelled out on eight `i16`.

type V8 = [i16; 8];
const ZERO: V8 = [0; 8];

/// `_mm_adds_epi16`.
#[inline]
fn adds(a: V8, b: V8) -> V8 {
    std::array::from_fn(|i| a[i].saturating_add(b[i]))
}

/// `_mm_max_epi16`.
#[inline]
fn max(a: V8, b: V8) -> V8 {
    std::array::from_fn(|i| a[i].max(b[i]))
}

/// `_mm_subs_epu16`, with `b` the same in each lane.
#[inline]
fn subs_u(a: V8, b: u16) -> V8 {
    std::array::from_fn(|i| (a[i] as u16).saturating_sub(b) as i16)
}

/// `_mm_slli_si128(x, 2)`.
#[inline]
fn slli2(a: V8) -> V8 {
    [0, a[0], a[1], a[2], a[3], a[4], a[5], a[6]]
}

/// `kswq_t` for `size` 2 and `m` 4: only the query profile matters, the rest is scratch.
pub(crate) struct Kswq {
    slen: usize,
    qp: Vec<V8>,
}

/// `ksw_qinit(2, qlen, query, 4, mat)`.
pub(crate) fn ksw_qinit(query: &[u8], mat: &[i8; 16]) -> Kswq {
    const M: usize = 4;
    const P: usize = 8; // # values per __m128i
    let qlen = query.len();
    let slen = (qlen + P - 1) / P; // segmented length
    let mut qp = Vec::with_capacity(slen * M);
    for a in 0..M {
        let ma = &mat[a * M..];
        for i in 0..slen {
            let mut v = ZERO;
            for (lane, x) in v.iter_mut().enumerate() {
                let k = i + lane * slen;
                *x = if k >= qlen {
                    0
                } else {
                    ma[query[k] as usize] as i16
                };
            }
            qp.push(v);
        }
    }
    Kswq { slen, qp }
}

/// `ksw_i16`: `.score`. A query of length 0 has `slen` 0, the core loop does nothing and C's load
/// of `H0[-1]`, inside its own struct, goes nowhere: the score is 0.
pub(crate) fn ksw_i16(q: &Kswq, target: &[u8], gapo: i32, gape: i32) -> i32 {
    let slen = q.slen;
    if slen == 0 {
        return 0;
    }
    let gapoe = (gapo + gape) as u16;
    let gape = gape as u16;
    let mut h0 = vec![ZERO; slen];
    let mut h1 = vec![ZERO; slen];
    let mut e_ = vec![ZERO; slen];
    let mut gmax = 0i32;
    for &t in target {
        let s = &q.qp[t as usize * slen..t as usize * slen + slen];
        let mut f = ZERO;
        let mut mx = ZERO;
        let mut h = slli2(h0[slen - 1]);
        for j in 0..slen {
            h = adds(h, s[j]);
            let mut e = e_[j];
            h = max(h, e);
            h = max(h, f);
            mx = max(mx, h);
            h1[j] = h;
            h = subs_u(h, gapoe);
            e = subs_u(e, gape);
            e = max(e, h);
            e_[j] = e;
            f = subs_u(f, gape);
            f = max(f, h);
            h = h0[j];
        }
        'end_loop8: for _ in 0..16 {
            f = slli2(f);
            for h1j in h1.iter_mut() {
                let mut h = *h1j;
                h = max(h, f);
                *h1j = h;
                h = subs_u(h, gapoe);
                f = subs_u(f, gape);
                // !_mm_movemask_epi8(_mm_cmpgt_epi16(f, h))
                if !f.iter().zip(h.iter()).any(|(a, b)| a > b) {
                    break 'end_loop8;
                }
            }
        }
        // __max_8: mx starts at zero, and _mm_extract_epi16 zero-extends.
        let imax = *mx.iter().max().unwrap() as u16 as i32;
        if imax > gmax {
            gmax = imax;
        }
        std::mem::swap(&mut h0, &mut h1);
    }
    gmax
}

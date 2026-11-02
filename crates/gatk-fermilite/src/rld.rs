//! `rld0.c`: the run-length delta-coded FM-index.
//!
//! `rld0.h` defines `_DNA_ONLY`, so the rank functions go through `rld_dec0_fast_dna`, except for
//! the second half of `rld_rank2a`, which uses `rld_dec0`. Both are here because they stop at
//! different places when a block runs out.

const RLD_LBITS: u32 = 23;
const RLD_LSIZE: usize = 1 << RLD_LBITS;
const RLD_LMASK: u64 = RLD_LSIZE as u64 - 1;
const RLD_IBITS_PLUS: i32 = 4;

/// `LogTable256`.
fn log_table_256(v: u32) -> i32 {
    if v == 0 {
        -1
    } else {
        31 - v.leading_zeros() as i32
    }
}

/// `ilog2_32`.
#[inline]
fn ilog2_32(v: u32) -> i32 {
    let tt = v >> 16;
    if tt != 0 {
        let t = tt >> 8;
        return if t != 0 {
            24 + log_table_256(t)
        } else {
            16 + log_table_256(tt)
        };
    }
    let t = v >> 8;
    if t != 0 {
        8 + log_table_256(t)
    } else {
        log_table_256(v)
    }
}

/// `ilog2`.
#[inline]
fn ilog2(v: u64) -> i32 {
    if v >> 32 != 0 {
        32 + ilog2_32((v >> 32) as u32)
    } else {
        ilog2_32(v as u32)
    }
}

/// `rld_delta_enc1`.
#[inline]
fn rld_delta_enc1(x: u64) -> (u64, i32) {
    let y = ilog2(x);
    let z = ilog2_32((y + 1) as u32);
    let width = (z << 1) + 1 + y;
    ((x ^ (1u64 << y)) | (((y + 1) as u64) << y), width)
}

/// `rlditr_t`: `p`, `shead` and `stail` are word offsets into chunk `i`.
#[derive(Clone, Copy)]
pub(crate) struct RldItr {
    r: i32,
    c: i32,
    l: i64,
    p: usize,
    shead: usize,
    stail: usize,
    i: usize,
    pc: usize,
}

/// `rld_t`, for `asize` 6 and `bbits` 3, which is the only `rld_init` in the library.
pub(crate) struct Rld {
    abits: i32,
    ibits: i32,
    offset0: [usize; 3],
    n_bytes: u64,
    z: Vec<Vec<u64>>,
    /// After `enc_finish`, `cnt[c]` is the number of symbols below `c`.
    pub(crate) cnt: [u64; 7],
    /// After `enc_finish`, `mcnt[c + 1]` is the number of `c`, and `mcnt[0]` the total.
    pub(crate) mcnt: [u64; 7],
    n_frames: u64,
    frame: Vec<u64>,
}

const ASIZE: usize = 6;
const ASIZE1: usize = 7;
const SSIZE: usize = 8;

/// `rld_block_type`.
#[inline]
fn block_type(x: u64) -> usize {
    (x >> 62) as usize
}

impl Rld {
    /// `rld_init(6, 3)` and `rld_itr_init(e, itr, 0)`.
    pub(crate) fn new() -> (Rld, RldItr) {
        let asize1 = ASIZE1;
        let e = Rld {
            abits: ilog2(ASIZE as u64) + 1,
            ibits: 0,
            offset0: [(asize1 * 16 + 63) / 64, (asize1 * 32 + 63) / 64, asize1],
            n_bytes: 0,
            z: vec![vec![0u64; RLD_LSIZE]],
            cnt: [0; 7],
            mcnt: [0; 7],
            n_frames: 0,
            frame: Vec::new(),
        };
        let itr = RldItr {
            i: 0,
            shead: 0,
            stail: e.stail(0),
            p: e.offset0[0],
            r: 64,
            c: -1,
            l: 0,
            pc: 0,
        };
        (e, itr)
    }

    /// `rld_get_stail`.
    #[inline]
    fn stail(&self, shead: usize) -> usize {
        shead + SSIZE - if shead + SSIZE == RLD_LSIZE { 2 } else { 1 }
    }

    /// `enc_next_block`: a block's header holds the counts of the block before it.
    fn enc_next_block(&mut self, itr: &mut RldItr) {
        if itr.stail + 2 == RLD_LSIZE {
            self.z.push(vec![0u64; RLD_LSIZE]);
            itr.i = self.z.len() - 1;
            itr.shead = 0;
        } else {
            itr.shead += SSIZE;
        }
        let d: [u64; 7] = std::array::from_fn(|i| self.cnt[i].wrapping_sub(self.mcnt[i]));
        let head = &mut self.z[itr.i][itr.shead..];
        let ty: usize;
        if d[0] < 0x4000 {
            for (i, &x) in d.iter().enumerate() {
                let (w, s) = (i / 4, 16 * (i % 4));
                head[w] = (head[w] & !(0xffffu64 << s)) | ((x & 0xffff) << s);
            }
            ty = 0;
        } else if d[0] < 0x4000_0000 {
            for (i, &x) in d.iter().enumerate() {
                let (w, s) = (i / 2, 32 * (i % 2));
                head[w] = (head[w] & !(0xffff_ffffu64 << s)) | ((x & 0xffff_ffff) << s);
            }
            ty = 1;
        } else {
            head[..7].copy_from_slice(&d);
            ty = 2;
        }
        head[0] |= (ty as u64) << 62;
        itr.p = itr.shead + self.offset0[ty];
        itr.stail = self.stail(itr.shead);
        itr.r = 64;
        self.mcnt = self.cnt;
    }

    /// `rld_enc1`.
    fn enc1(&mut self, itr: &mut RldItr, l: i64, c: u8) {
        let (d, mut w) = rld_delta_enc1(l as u64);
        let x = (d << self.abits) | c as u64;
        w += self.abits;
        if w >= itr.r && itr.p == itr.stail {
            self.enc_next_block(itr);
        }
        let z = &mut self.z[itr.i];
        if w > itr.r {
            w -= itr.r;
            z[itr.p] |= x >> w;
            itr.p += 1;
            itr.r = 64 - w;
            z[itr.p] = x << itr.r;
        } else {
            itr.r -= w;
            z[itr.p] |= x << itr.r;
        }
        self.cnt[0] = self.cnt[0].wrapping_add(l as u64);
        self.cnt[c as usize + 1] = self.cnt[c as usize + 1].wrapping_add(l as u64);
    }

    /// `rld_enc`.
    pub(crate) fn enc(&mut self, itr: &mut RldItr, l: i64, c: u8) {
        if l == 0 {
            return;
        }
        if itr.c != c as i32 {
            if itr.l != 0 {
                let (pl, pc) = (itr.l, itr.c as u8);
                self.enc1(itr, pl, pc);
            }
            itr.l = l;
            itr.c = c as i32;
        } else {
            itr.l += l;
        }
    }

    /// `rld_seek_blk` and the marginal counts of the block there.
    fn block_counts(&self, k: u64) -> [u64; ASIZE] {
        let q = &self.z[(k >> RLD_LBITS) as usize][(k & RLD_LMASK) as usize..];
        let mut out = [0u64; ASIZE];
        match block_type(q[0]) {
            0 => {
                for (j, o) in out.iter_mut().enumerate() {
                    let i = j + 1;
                    *o = (q[i / 4] >> (16 * (i % 4))) & 0xffff;
                }
            }
            1 => {
                for (j, o) in out.iter_mut().enumerate() {
                    let i = j + 1;
                    *o = (q[i / 2] >> (32 * (i % 2))) & 0x3fff_ffff;
                }
            }
            _ => out.copy_from_slice(&q[1..7]),
        }
        out
    }

    /// `rld_rank_index`.
    fn rank_index(&mut self) {
        let n_blks = self.n_bytes * 8 / 64 / SSIZE as u64 + 1;
        let last = (self.n_bytes >> 3 >> 3) << 3;
        self.ibits = ilog2(self.mcnt[0] / n_blks) + RLD_IBITS_PLUS;
        self.n_frames = ((self.mcnt[0] + (1u64 << self.ibits) - 1) >> self.ibits) + 1;
        self.frame = vec![0u64; self.n_frames as usize * ASIZE1];
        let mut cnt = [0u64; ASIZE];
        let mut i = SSIZE as u64;
        let mut k = 1u64;
        while i <= last {
            let bc = self.block_counts(i);
            for j in 0..ASIZE {
                cnt[j] = cnt[j].wrapping_add(bc[j]);
            }
            let sum = cnt.iter().fold(0u64, |a, b| a.wrapping_add(*b));
            while sum >= k << self.ibits {
                k += 1;
            }
            if k < self.n_frames {
                let x = k as usize * ASIZE1;
                self.frame[x] = i;
                self.frame[x + 1..x + 1 + ASIZE].copy_from_slice(&cnt);
            }
            i += SSIZE as u64;
        }
        c_assert!(
            k >= self.n_frames - 1,
            "rld_rank_index: k >= e->n_frames - 1"
        );
        for k in 1..self.n_frames as usize {
            // fill zero cells
            let x = k * ASIZE1;
            if self.frame[x] == 0 {
                self.frame.copy_within(x - ASIZE1..x, x);
            }
        }
    }

    /// `rld_enc_finish`.
    pub(crate) fn enc_finish(&mut self, itr: &mut RldItr) {
        if itr.l != 0 {
            let (pl, pc) = (itr.l, itr.c as u8);
            self.enc1(itr, pl, pc);
        }
        self.enc_next_block(itr);
        self.n_bytes = ((self.z.len() as u64 - 1) * RLD_LSIZE as u64 + itr.p as u64) * 8;
        self.cnt[0] = 0;
        for i in 1..=ASIZE {
            self.cnt[i] = self.cnt[i].wrapping_add(self.cnt[i - 1]);
        }
        self.rank_index();
    }

    /// `rld_locate_blk`: fills `cnt` and `sum`, and returns the end of the block.
    ///
    /// When the block after the located one opens a new chunk, C's `itr->i` has moved on while
    /// `itr->p` stays behind, and `rld_get_stail` subtracts one from the other. `pc` keeps the
    /// chunk `p` is in. It takes a 64 MB index to get there.
    #[inline]
    fn locate_blk(&self, itr: &mut RldItr, k: u64, cnt: &mut [u64; 6], sum: &mut u64) -> u64 {
        let zi = (k >> self.ibits) as usize * ASIZE1;
        let z = &self.frame[zi..zi + ASIZE1];
        itr.i = (z[0] >> RLD_LBITS) as usize;
        itr.pc = itr.i;
        let mut q = (z[0] & RLD_LMASK) as usize;
        itr.p = q;
        *sum = 0;
        for j in 1..ASIZE1 {
            cnt[j - 1] = z[j];
            *sum = sum.wrapping_add(z[j]);
        }
        let mut c;
        loop {
            // seek to the small block
            q += SSIZE;
            if q == RLD_LSIZE {
                itr.i += 1;
                q = 0;
            }
            let chunk = &self.z[itr.i];
            let ty = block_type(chunk[q]);
            c = if ty == 2 {
                chunk[q] & 0x3fff_ffff_ffff_ffff
            } else if ty == 1 {
                chunk[q] & 0xffff_ffff
            } else {
                chunk[q] & 0xffff
            };
            if sum.wrapping_add(c) > k {
                break;
            }
            match ty {
                0 => {
                    for (j, x) in cnt.iter_mut().enumerate() {
                        let i = j + 1;
                        *x = x.wrapping_add((chunk[q + i / 4] >> (16 * (i % 4))) & 0xffff);
                    }
                }
                1 => {
                    for (j, x) in cnt.iter_mut().enumerate() {
                        let i = j + 1;
                        *x = x.wrapping_add((chunk[q + i / 2] >> (32 * (i % 2))) & 0x3fff_ffff);
                    }
                }
                _ => {
                    for (j, x) in cnt.iter_mut().enumerate() {
                        *x = x.wrapping_add(chunk[q + j + 1]);
                    }
                }
            }
            *sum = sum.wrapping_add(c);
            itr.p = q;
            itr.pc = itr.i;
        }
        itr.shead = itr.p;
        itr.stail = if itr.pc == itr.i {
            self.stail(itr.shead)
        } else {
            itr.shead + SSIZE - 1
        };
        let ty = block_type(self.z[itr.pc][itr.shead]);
        itr.p += self.offset0[ty];
        itr.r = 64;
        c.wrapping_add(*sum)
    }

    /// `rld_dec0_fast_dna`. Does not check the block boundary; the callers stop first.
    #[inline]
    fn dec0_fast_dna(&self, itr: &mut RldItr, c: &mut i32) -> i64 {
        let z = &self.z[itr.pc];
        let x = if itr.r == 64 {
            z[itr.p]
        } else {
            (z[itr.p] << (64 - itr.r)) | (z[itr.p + 1] >> itr.r)
        };
        if x >> 63 == 0 {
            let mut w = ((0x3333_3333_5555_779bu64 >> ((x >> 59) << 2)) & 0xf) as i32;
            let l = ((x >> (64 - w)) as i32) - 1;
            let y = ((x << w) >> (64 - l)) | (1u32.wrapping_shl(l as u32) as u64);
            w += l;
            *c = ((x << w) >> 61) as i32;
            w += 3;
            itr.r -= w;
            if itr.r <= 0 {
                itr.p += 1;
                itr.r += 64;
            }
            y as i64
        } else {
            *c = ((x << 1) >> 61) as i32;
            itr.r -= 4;
            if itr.r <= 0 {
                itr.p += 1;
                itr.r += 64;
            }
            1
        }
    }

    /// `rld_dec0`.
    #[inline]
    fn dec0(&self, itr: &mut RldItr, c: &mut i32) -> i64 {
        let z = &self.z[itr.pc];
        let tail = if itr.p != itr.stail && itr.r != 64 {
            z[itr.p + 1] >> itr.r
        } else {
            0
        };
        let x = (z[itr.p] << (64 - itr.r)) | tail;
        let mut w: i32;
        let y: i64;
        if x >> 63 == 0 {
            w = ((0x3333_3333_5555_779bu64 >> ((x >> 59) << 2)) & 0xf) as i32;
            if w == 0xb && x >> 58 == 0 {
                return 0;
            }
            let l = ((x >> (64 - w)) as i32) - 1;
            y = (((x << w) >> (64 - l)) | (1u32.wrapping_shl(l as u32) as u64)) as i64;
            w += l;
        } else {
            w = 1;
            y = 1;
        }
        *c = ((x << w) >> (64 - self.abits)) as i32;
        w += self.abits;
        if itr.r > w {
            itr.r -= w;
        } else {
            itr.p += 1;
            itr.r = 64 + itr.r - w;
        }
        y
    }

    /// `rld_rank1a`: the symbol at `k - 1`, -1 when `k` is 0.
    pub(crate) fn rank1a(&self, k: u64, ok: &mut [u64; 6]) -> i32 {
        if k == 0 {
            *ok = [0; 6];
            return -1;
        }
        let mut itr = self.blank_itr();
        let mut z = 0u64;
        let mut a = -1;
        self.locate_blk(&mut itr, k - 1, ok, &mut z);
        loop {
            let l = self.dec0_fast_dna(&mut itr, &mut a) as u64;
            if z + l >= k {
                break;
            }
            z += l;
            ok[a as usize] += l;
        }
        ok[a as usize] += k - z;
        a
    }

    fn blank_itr(&self) -> RldItr {
        RldItr {
            r: 0,
            c: 0,
            l: 0,
            p: 0,
            shead: 0,
            stail: 0,
            i: 0,
            pc: 0,
        }
    }

    /// `rld_rank2a`.
    pub(crate) fn rank2a(&self, k: u64, l: u64, ok: &mut [u64; 6], ol: &mut [u64; 6]) {
        if k == 0 {
            *ok = [0; 6];
            self.rank1a(l, ol);
            return;
        }
        let mut itr = self.blank_itr();
        let mut z = 0u64;
        let mut a = -1;
        let mut len;
        let y = self.locate_blk(&mut itr, k - 1, ok, &mut z); // the block bracketing k
        loop {
            // compute ok[]
            len = self.dec0_fast_dna(&mut itr, &mut a) as u64;
            if z + len >= k {
                break;
            }
            z += len;
            ok[a as usize] += len;
        }
        if y > l {
            // we do not need to decode other blocks
            *ol = *ok;
            ok[a as usize] += k - z;
            if z + len < l {
                // we need to decode the next run
                z += len;
                ol[a as usize] += len;
                loop {
                    len = self.dec0(&mut itr, &mut a) as u64;
                    if z + len >= l {
                        break;
                    }
                    z += len;
                    ol[a as usize] += len;
                }
            }
            ol[a as usize] = ol[a as usize].wrapping_add(l.wrapping_sub(z));
        } else {
            // we have to decode other blocks
            ok[a as usize] += k - z;
            self.rank1a(l, ol);
        }
    }

    /// `rld_extend`.
    pub(crate) fn extend(&self, ik: &RldIntv, ok: &mut [RldIntv; 6], is_back: bool) {
        let (fw, bw) = (usize::from(!is_back), usize::from(is_back));
        let (mut tk, mut tl) = ([0u64; 6], [0u64; 6]);
        self.rank2a(ik.x[fw], ik.x[fw].wrapping_add(ik.x[2]), &mut tk, &mut tl);
        for i in 0..6 {
            ok[i].x[fw] = self.cnt[i].wrapping_add(tk[i]);
            tl[i] = tl[i].wrapping_sub(tk[i]);
            ok[i].x[2] = tl[i];
        }
        ok[0].x[bw] = ik.x[bw];
        ok[4].x[bw] = ok[0].x[bw].wrapping_add(tl[0]);
        ok[3].x[bw] = ok[4].x[bw].wrapping_add(tl[4]);
        ok[2].x[bw] = ok[3].x[bw].wrapping_add(tl[3]);
        ok[1].x[bw] = ok[2].x[bw].wrapping_add(tl[2]);
        ok[5].x[bw] = ok[1].x[bw].wrapping_add(tl[1]);
    }
}

/// `rldintv_t`.
#[derive(Clone, Copy, Default)]
pub(crate) struct RldIntv {
    pub(crate) x: [u64; 3],
    pub(crate) info: u64,
}

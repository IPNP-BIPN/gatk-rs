//! `rle.c`, `rope.c` and `mrope.c`: the BWT under construction, a B+ tree of run-length blocks per
//! preceding symbol, built by inserting all strings a column at a time (`mr_insert_multi`).
//!
//! The multi-threaded branch of `mr_insert_multi` (`is_thr`) is not here: GATK has `n_threads` at 1.

const ROPE_MAX_DEPTH: usize = 80;
const ROPE_DEF_MAX_NODES: i32 = 64;
const ROPE_DEF_BLOCK_LEN: i32 = 512;
const RLE_MIN_SPACE: i32 = 18;
const NIL: usize = usize::MAX;

/// `rle_auxtab`.
const RLE_AUXTAB: [u8; 8] = [0x01, 0x11, 0x21, 0x31, 0x03, 0x13, 0x07, 0x17];

/// `*rle_nptr(block)`, a little-endian `uint16_t` as on x86-64.
#[inline]
pub(crate) fn rle_nptr(block: &[u8]) -> usize {
    u16::from_le_bytes([block[0], block[1]]) as usize
}

#[inline]
fn set_rle_nptr(block: &mut [u8], n: usize) {
    let b = (n as u16).to_le_bytes();
    block[0] = b[0];
    block[1] = b[1];
}

/// `rle_dec1`: moves `p` over one run and gives `(c, l)`.
#[inline]
pub(crate) fn rle_dec1(b: &[u8], p: &mut usize) -> (usize, i64) {
    let v = b[*p];
    let c = (v & 7) as usize;
    let l;
    if v & 0x80 == 0 {
        l = (v >> 3) as i64;
        *p += 1;
    } else if v >> 5 == 6 {
        l = ((v as i64 & 0x18) << 3) | (b[*p + 1] as i64 & 0x3f);
        *p += 2;
    } else {
        let n = (((v & 0x10) >> 2) + 4) as usize;
        let mut x = ((v >> 3) & 1) as i64;
        *p += 1;
        for _ in 1..n {
            x = (x << 6) | (b[*p] as i64 & 0x3f);
            *p += 1;
        }
        l = x;
    }
    (c, l)
}

/// `rle_enc1`.
fn rle_enc1(p: &mut [u8], c: i32, l: i64) -> usize {
    if l < 1 << 4 {
        p[0] = ((l << 3) | c as i64) as u8;
        1
    } else if l < 1 << 8 {
        p[0] = (0xC0 | ((l >> 6) << 3) | c as i64) as u8;
        p[1] = (0x80 | (l & 0x3f)) as u8;
        2
    } else if l < 1 << 19 {
        p[0] = (0xE0 | ((l >> 18) << 3) | c as i64) as u8;
        p[1] = (0x80 | ((l >> 12) & 0x3f)) as u8;
        p[2] = (0x80 | ((l >> 6) & 0x3f)) as u8;
        p[3] = (0x80 | (l & 0x3f)) as u8;
        4
    } else {
        p[0] = (0xF0 | ((l >> 42) << 3) | c as i64) as u8;
        let mut shift = 36;
        for q in p.iter_mut().take(8).skip(1) {
            *q = (0x80 | ((l >> shift) & 0x3f)) as u8;
            shift -= 6;
        }
        8
    }
}

/// `rle_insert_cached`: insert `a` after `x` symbols, marginal counts into `cnt`, the block's new
/// byte count back. `beg` and `bc` remember where the last insertion happened in the block.
#[allow(clippy::too_many_arguments)]
fn rle_insert_cached(
    b: &mut [u8],
    x: i64,
    a: i32,
    rl: i64,
    cnt: &mut [i64; 6],
    ec: &[i64; 6],
    beg: &mut i32,
    bc: &mut [i64; 6],
) -> i32 {
    let nptr = rle_nptr(b);
    let diff: i32;
    if nptr == 0 {
        *cnt = [0; 6];
        diff = rle_enc1(&mut b[2..], a, rl) as i32;
    } else {
        let end = 2 + nptr;
        let mut l = 0i64;
        let mut c: i32 = -1;
        let mut tmp = [0u8; 24];
        let mut beg_l: i64 = bc.iter().sum();
        let tot: i64 = ec.iter().sum();
        if x < beg_l {
            beg_l = 0;
            *beg = 0;
            *bc = [0; 6];
        }
        let mut p;
        let mut q;
        let mut z;
        if x == beg_l {
            p = 2 + *beg as usize;
            q = p;
            z = beg_l;
            *cnt = *bc;
        } else if x - beg_l <= ((tot - beg_l) >> 1) + ((tot - beg_l) >> 3) {
            // forward
            z = beg_l;
            p = 2 + *beg as usize;
            *cnt = *bc;
            while z < x {
                let (cc, ll) = rle_dec1(b, &mut p);
                c = cc as i32;
                l = ll;
                z += l;
                cnt[c as usize] += l;
            }
            q = p - 1;
            while b[q] >> 6 == 2 {
                q -= 1;
            }
        } else {
            // backward
            *cnt = *ec;
            z = tot;
            p = end;
            let mut t = 0;
            while z >= x {
                p -= 1;
                if b[p] >> 6 != 2 {
                    l |= if b[p] >> 7 != 0 {
                        ((RLE_AUXTAB[((b[p] >> 3) & 7) as usize] >> 4) as i64) << t
                    } else {
                        (b[p] >> 3) as i64
                    };
                    z -= l;
                    cnt[(b[p] & 7) as usize] -= l;
                    l = 0;
                    t = 0;
                } else {
                    l |= (b[p] as i64 & 0x3f) << t;
                    t += 6;
                }
            }
            q = p;
            let (cc, ll) = rle_dec1(b, &mut p);
            c = cc as i32;
            l = ll;
            z += l;
            cnt[c as usize] += l;
        }
        *beg = (q - 2) as i32;
        *bc = *cnt;
        if c >= 0 {
            // `bc[c] -= l` with `c` still -1 and `l` 0 writes nothing that matters: bc[-1] is
            // the stack word before it, and the amount is 0.
            bc[c as usize] -= l;
        }
        let mut n_bytes = p - q;
        if x == z && a != c && p < end {
            // then try the next run
            q = p;
            let (tc, tl) = rle_dec1(b, &mut q);
            if a == tc as i32 {
                c = tc as i32;
                n_bytes = q - p;
                l = tl;
                z += l;
                p = q;
                cnt[tc] += tl;
            }
        }
        if z != x {
            cnt[c as usize] -= z - x;
        }
        let pre = x - (z - l);
        p -= n_bytes;
        let mut n_bytes2;
        if a == c {
            // insert to the same run
            n_bytes2 = rle_enc1(&mut tmp, c, l + rl);
        } else if x == z {
            // at the end; append to the existing run
            p += n_bytes;
            n_bytes = 0;
            n_bytes2 = rle_enc1(&mut tmp, a, rl);
        } else {
            // break the current run
            n_bytes2 = rle_enc1(&mut tmp, c, pre);
            n_bytes2 += rle_enc1(&mut tmp[n_bytes2..], a, rl);
            n_bytes2 += rle_enc1(&mut tmp[n_bytes2..], c, l - pre);
        }
        if n_bytes != n_bytes2 && end != p + n_bytes {
            b.copy_within(p + n_bytes..end, p + n_bytes2);
        }
        b[p..p + n_bytes2].copy_from_slice(&tmp[..n_bytes2]);
        diff = n_bytes2 as i32 - n_bytes as i32;
    }
    let n = (nptr as i32 + diff) as u16;
    set_rle_nptr(b, n as usize);
    n as i32
}

/// `rle_split`.
fn rle_split(block: &[u8], new_block: &mut [u8]) -> usize {
    let n = rle_nptr(block);
    let end = 2 + n;
    let mut q = 2 + (n >> 1);
    while block[q] >> 6 == 2 {
        q -= 1;
    }
    new_block[2..2 + end - q].copy_from_slice(&block[q..end]);
    set_rle_nptr(new_block, end - q);
    q - 2
}

/// `rle_count`.
fn rle_count(block: &[u8], cnt: &mut [i64; 6]) {
    let mut q = 2;
    let end = q + rle_nptr(block);
    while q < end {
        let (c, l) = rle_dec1(block, &mut q);
        cnt[c] += l;
    }
}

/// `move_backward`, inside `rle_rank2a`.
#[inline]
fn move_backward(
    block: &[u8],
    x: i64,
    p: &mut usize,
    z: &mut i64,
    l: &mut i64,
    t: &mut i32,
    cnt: &mut [i64; 6],
) {
    while *z >= x {
        *p -= 1;
        let v = block[*p];
        if v >> 6 != 2 {
            *l |= if v >> 7 != 0 {
                ((RLE_AUXTAB[((v >> 3) & 7) as usize] >> 4) as i64) << *t
            } else {
                (v >> 3) as i64
            };
            *z -= *l;
            cnt[(v & 7) as usize] -= *l;
            *l = 0;
            *t = 0;
        } else {
            *l |= (v as i64 & 0x3f) << *t;
            *t += 6;
        }
    }
}

/// `rle_rank2a`; `rle_rank1a` with `cy` absent.
fn rle_rank2a(
    block: &[u8],
    x: i64,
    y: i64,
    cx: &mut [i64; 6],
    cy: Option<&mut [i64; 6]>,
    ec: &[i64; 6],
) {
    let y = if y >= x { y } else { x };
    let tot: i64 = ec.iter().sum();
    if tot == 0 {
        return;
    }
    if x <= (tot - y) + (tot >> 3) {
        let mut c = 0usize;
        let mut z = 0i64;
        let mut cnt = [0i64; 6];
        let mut p = 2usize;
        while z < x {
            let (cc, l) = rle_dec1(block, &mut p);
            c = cc;
            z += l;
            cnt[c] += l;
        }
        for a in 0..6 {
            cx[a] += cnt[a];
        }
        cx[c] -= z - x;
        if let Some(cy) = cy {
            while z < y {
                let (cc, l) = rle_dec1(block, &mut p);
                c = cc;
                z += l;
                cnt[c] += l;
            }
            for a in 0..6 {
                cy[a] += cnt[a];
            }
            cy[c] -= z - y;
        }
    } else {
        let mut t = 0;
        let mut l = 0i64;
        let mut z = tot;
        let mut cnt = *ec;
        let mut p = 2 + rle_nptr(block);
        if let Some(cy) = cy {
            move_backward(block, y, &mut p, &mut z, &mut l, &mut t, &mut cnt);
            for a in 0..6 {
                cy[a] += cnt[a];
            }
            cy[(block[p] & 7) as usize] += y - z;
        }
        move_backward(block, x, &mut p, &mut z, &mut l, &mut t, &mut cnt);
        for a in 0..6 {
            cx[a] += cnt[a];
        }
        cx[(block[p] & 7) as usize] += x - z;
    }
}

/// `rpnode_t`, in a flat pool; `n` and `is_bottom` count only on the first node of a bucket.
#[derive(Clone, Copy, Default)]
struct RpNode {
    p: usize,
    l: i64,
    n: i32,
    is_bottom: bool,
    c: [i64; 6],
}

/// `rpcache_t`.
struct RpCache {
    beg: i32,
    bc: [i64; 6],
    p: usize,
}

impl RpCache {
    fn new() -> RpCache {
        RpCache {
            beg: 0,
            bc: [0; 6],
            p: NIL,
        }
    }
    fn reset(&mut self) {
        *self = RpCache::new();
    }
}

/// `rope_t`, its two memory pools as flat vectors.
pub(crate) struct Rope {
    max_nodes: i32,
    block_len: i32,
    pub(crate) c: [i64; 6],
    root: usize,
    nodes: Vec<RpNode>,
    leaf: Vec<u8>,
}

impl Rope {
    /// `rope_init`.
    fn new(max_nodes: i32, mut block_len: i32) -> Rope {
        if block_len < 32 {
            block_len = 32;
        }
        let mut r = Rope {
            max_nodes: ((max_nodes + 1) >> 1) << 1,
            block_len: ((block_len + 7) >> 3) << 3,
            c: [0; 6],
            root: 0,
            nodes: Vec::new(),
            leaf: Vec::new(),
        };
        r.root = r.alloc_node();
        r.nodes[r.root].n = 1;
        r.nodes[r.root].is_bottom = true;
        r.nodes[r.root].p = r.alloc_leaf();
        r
    }

    fn alloc_node(&mut self) -> usize {
        let at = self.nodes.len();
        self.nodes
            .resize(at + self.max_nodes as usize, RpNode::default());
        at
    }

    fn alloc_leaf(&mut self) -> usize {
        let at = self.leaf.len();
        self.leaf.resize(at + self.block_len as usize, 0);
        at
    }

    /// `split_node`: `u` is the first node of the bucket of `v`, both `NIL` at the root.
    fn split_node(&mut self, mut u: usize, mut v: usize) -> usize {
        let i = if u == NIL { 0 } else { v - u };
        if u == NIL {
            // only happens at the root; add a new root
            v = self.alloc_node();
            u = v;
            let c = self.c;
            self.nodes[v].n = 1;
            self.nodes[v].p = self.root;
            self.nodes[v].c = c;
            self.nodes[v].l += c.iter().sum::<i64>();
            self.root = v;
        }
        let un = self.nodes[u].n as usize;
        if i != un - 1 {
            // then make room for a new node
            self.nodes.copy_within(v + 1..v + 1 + (un - i - 1), v + 2);
        }
        self.nodes[u].n += 1;
        let w = v + 1;
        self.nodes[w] = RpNode::default();
        let bottom = self.nodes[u].is_bottom;
        self.nodes[w].p = if bottom {
            self.alloc_leaf()
        } else {
            self.alloc_node()
        };
        if bottom {
            let (p, q) = (self.nodes[v].p, self.nodes[w].p);
            let bl = self.block_len as usize;
            // The new leaf is always later in the pool than the one being split.
            let (left, right) = self.leaf.split_at_mut(q);
            let block = &mut left[p..p + bl];
            let new_n = rle_split(block, &mut right[..bl]);
            set_rle_nptr(block, new_n);
            let mut c = [0i64; 6];
            rle_count(&right[..bl], &mut c);
            self.nodes[w].c = c;
        } else {
            let (p, q) = (self.nodes[v].p, self.nodes[w].p);
            let half = (self.max_nodes >> 1) as usize;
            self.nodes[p].n -= half as i32;
            let pn = self.nodes[p].n as usize;
            self.nodes.copy_within(p + pn..p + pn + half, q);
            self.nodes[q].n = half as i32;
            self.nodes[q].is_bottom = self.nodes[p].is_bottom;
            let mut c = [0i64; 6];
            for node in &self.nodes[q..q + half] {
                for (cj, nc) in c.iter_mut().zip(node.c.iter()) {
                    *cj += nc;
                }
            }
            self.nodes[w].c = c;
        }
        for j in 0..6 {
            let wc = self.nodes[w].c[j];
            self.nodes[w].l += wc;
            self.nodes[v].c[j] -= wc;
        }
        let wl = self.nodes[w].l;
        self.nodes[v].l -= wl;
        v
    }

    /// `rope_insert_run`: insert `rl` of `a` after `x` symbols, and return rank(a, x).
    fn insert_run(&mut self, x: i64, a: i32, rl: i64, cache: &mut RpCache) -> i64 {
        let (mut u, mut v, mut p) = (NIL, NIL, self.root);
        let (mut y, mut z) = (0i64, 0i64);
        let au = a as usize;
        loop {
            // top-down update; searching and node splitting are done together in one pass
            if self.nodes[p].n == self.max_nodes {
                v = self.split_node(u, v);
                if y + self.nodes[v].l < x {
                    y += self.nodes[v].l;
                    z += self.nodes[v].c[au];
                    v += 1;
                    p = self.nodes[v].p;
                }
            }
            u = p;
            if v != NIL && x - y > self.nodes[v].l >> 1 {
                // search backwardly for the right node to descend
                p += self.nodes[p].n as usize - 1;
                y += self.nodes[v].l;
                z += self.nodes[v].c[au];
                while y >= x {
                    y -= self.nodes[p].l;
                    z -= self.nodes[p].c[au];
                    p = p.wrapping_sub(1);
                }
                p = p.wrapping_add(1);
            } else {
                // search forwardly
                while y + self.nodes[p].l < x {
                    y += self.nodes[p].l;
                    z += self.nodes[p].c[au];
                    p += 1;
                }
            }
            c_assert!(
                ((p - u) as i32) < self.nodes[u].n,
                "rope_insert_run: p - u < u->n"
            );
            if v != NIL {
                self.nodes[v].c[au] += rl;
                self.nodes[v].l += rl;
            }
            v = p;
            p = self.nodes[p].p;
            if self.nodes[u].is_bottom {
                break;
            }
        }
        self.c[au] += rl;
        if cache.p != p {
            cache.reset();
        }
        let mut cnt = [0i64; 6];
        let bl = self.block_len as usize;
        let ec = self.nodes[v].c;
        let n_runs = rle_insert_cached(
            &mut self.leaf[p..p + bl],
            x - y,
            a,
            rl,
            &mut cnt,
            &ec,
            &mut cache.beg,
            &mut cache.bc,
        );
        cache.p = p;
        z += cnt[au];
        self.nodes[v].c[au] += rl;
        self.nodes[v].l += rl;
        if n_runs + RLE_MIN_SPACE > self.block_len {
            self.split_node(u, v);
            cache.reset();
        }
        z
    }

    /// `rope_count_to_leaf`: the bottom node, and `rest`.
    fn count_to_leaf(&self, x: i64, cx: &mut [i64; 6]) -> (usize, i64) {
        let (mut v, mut p) = (NIL, self.root);
        let mut y = 0i64;
        *cx = [0; 6];
        loop {
            let u = p;
            if v != NIL && x - y > self.nodes[v].l >> 1 {
                p += self.nodes[p].n as usize - 1;
                y += self.nodes[v].l;
                for (a, c) in cx.iter_mut().enumerate() {
                    *c += self.nodes[v].c[a];
                }
                while y >= x {
                    y -= self.nodes[p].l;
                    for (a, c) in cx.iter_mut().enumerate() {
                        *c -= self.nodes[p].c[a];
                    }
                    p = p.wrapping_sub(1);
                }
                p = p.wrapping_add(1);
            } else {
                while y + self.nodes[p].l < x {
                    y += self.nodes[p].l;
                    for (a, c) in cx.iter_mut().enumerate() {
                        *c += self.nodes[p].c[a];
                    }
                    p += 1;
                }
            }
            v = p;
            p = self.nodes[p].p;
            if self.nodes[u].is_bottom {
                break;
            }
        }
        (v, x - y)
    }

    /// `rope_rank2a`, always with `y >= x` and `cy` here.
    fn rank2a(&self, x: i64, y: i64, cx: &mut [i64; 6], cy: &mut [i64; 6]) {
        let bl = self.block_len as usize;
        let (mut v, mut rest) = self.count_to_leaf(x, cx);
        let n = self.nodes[v];
        let block = &self.leaf[n.p..n.p + bl];
        if y < x {
            rle_rank2a(block, rest, -1, cx, None, &n.c);
        } else if rest + (y - x) <= n.l {
            *cy = *cx;
            rle_rank2a(block, rest, rest + (y - x), cx, Some(cy), &n.c);
        } else {
            rle_rank2a(block, rest, -1, cx, None, &n.c);
            (v, rest) = self.count_to_leaf(y, cy);
            let n = self.nodes[v];
            rle_rank2a(&self.leaf[n.p..n.p + bl], rest, -1, cy, None, &n.c);
        }
    }

    /// `rope_itr_first` and `rope_itr_next_block`, a run of callbacks: each leaf block in order.
    fn for_each_block(&self, mut f: impl FnMut(&[u8])) {
        let bl = self.block_len as usize;
        let mut pa = [0usize; ROPE_MAX_DEPTH];
        let mut ia = [0i32; ROPE_MAX_DEPTH];
        let mut d: i32 = 0;
        pa[0] = self.root;
        while !self.nodes[pa[d as usize]].is_bottom {
            d += 1;
            pa[d as usize] = self.nodes[pa[d as usize - 1]].p;
        }
        loop {
            c_assert!(
                (d as usize) < ROPE_MAX_DEPTH || d < 0,
                "rope_itr_next_block: i->d < ROPE_MAX_DEPTH"
            );
            if d < 0 {
                return;
            }
            let du = d as usize;
            let ret = self.nodes[pa[du] + ia[du] as usize].p;
            while d >= 0 {
                let du = d as usize;
                ia[du] += 1;
                if ia[du] != self.nodes[pa[du]].n {
                    break;
                }
                ia[du] = 0;
                d -= 1;
            }
            if d >= 0 {
                while !self.nodes[pa[d as usize]].is_bottom {
                    d += 1;
                    let du = d as usize;
                    pa[du] = self.nodes[pa[du - 1] + ia[du - 1] as usize].p;
                }
            }
            f(&self.leaf[ret..ret + bl]);
        }
    }
}

/// `triple64_t`.
#[derive(Clone, Copy, Default)]
struct Triple {
    l: u64,
    u: u64,
    c: u8,
    p: usize,
}

/// `rope_comp6` is never needed: `is_comp` is written out in `mr_insert_multi_aux`.
fn mr_insert_multi_aux(rope: &mut Rope, a: &mut [Triple], s: &[u8], is_comp: bool) {
    let m = a.len();
    let mut cache = RpCache::new();
    for t in a.iter_mut() {
        t.c = s[t.p];
        t.p += 1;
    }
    let mut beg = 0usize;
    for k in 1..=m {
        if k == m || a[k].u != a[k - 1].u {
            let (l, u) = (a[beg].l as i64, a[beg].u as i64);
            let (mut tl, mut tu) = ([0i64; 6], [0i64; 6]);
            if l == u && k == beg + 1 {
                // special case; still works without the following block
                let r = rope.insert_run(l, a[beg].c as i32, 1, &mut cache) as u64;
                a[beg].l = r;
                a[beg].u = r;
                beg = k;
                continue;
            } else if l != u {
                rope.rank2a(l, u, &mut tl, &mut tu);
            }
            let mut c = [0i64; 6];
            for t in &a[beg..k] {
                c[t.c as usize] += 1;
            }
            // insert sentinel
            if c[0] != 0 {
                rope.insert_run(l, 0, c[0], &mut cache);
            }
            // insert A/C/G/T
            let mut x = l + c[0] + (tu[0] - tl[0]);
            let (start, end, step): (i32, i32, i32) = if is_comp { (4, 0, -1) } else { (1, 5, 1) };
            let mut b = start;
            while b != end {
                let bu = b as usize;
                let size = tu[bu] - tl[bu];
                if c[bu] != 0 {
                    tl[bu] = rope.insert_run(x, b, c[bu], &mut cache);
                    tu[bu] = tl[bu] + size;
                }
                x += c[bu] + size;
                b += step;
            }
            // insert N
            if c[5] != 0 {
                tu[5] -= tl[5];
                tl[5] = rope.insert_run(x, 5, c[5], &mut cache);
                tu[5] += tl[5];
            }
            // update a[]
            for t in &mut a[beg..k] {
                t.l = tl[t.c as usize] as u64;
                t.u = tu[t.c as usize] as u64;
            }
            beg = k;
        }
    }
}

/// `mrope_t`: six ropes, by the symbol the suffixes in them start with.
pub(crate) struct Mrope {
    r: [Rope; 6],
}

impl Mrope {
    /// `mr_init(ROPE_DEF_MAX_NODES, ROPE_DEF_BLOCK_LEN, MR_SO_RCLO)`.
    pub(crate) fn new() -> Mrope {
        Mrope {
            r: std::array::from_fn(|_| Rope::new(ROPE_DEF_MAX_NODES, ROPE_DEF_BLOCK_LEN)),
        }
    }

    /// `mr_insert_multi(mr, len, s, 0)` under `MR_SO_RCLO`: the strings end up in reverse
    /// complement lexicographic order, not in the order they were given.
    pub(crate) fn insert_multi(&mut self, s: &[u8]) {
        let is_comp = true;
        let len = s.len();
        c_assert!(
            len > 0 && s[len - 1] == 0,
            "mr_insert_multi: len > 0 && s[len-1] == 0"
        );
        let mut m = 0usize;
        for &c in s {
            if c == 0 {
                m += 1;
            }
        }
        let mut curr = vec![Triple::default(); m];
        let mut prev = vec![Triple::default(); m];
        let mut k = 0;
        let mut q = 0;
        for (p, &c) in s.iter().enumerate() {
            if c == 0 {
                prev[k].p = q;
                k += 1;
                q = p + 1;
            }
        }
        // the first (actually the last) column; n0, the strings already in, is 0
        mr_insert_multi_aux(&mut self.r[0], &mut prev, s, is_comp);

        let mut n0 = 0usize; // the number of inserted strings
        loop {
            let mut c = [0i64; 6];
            for t in &prev[n0..m] {
                c[t.c as usize] += 1;
            }
            let mut q = [0usize; 6];
            q[0] = n0;
            for b in 1..6 {
                q[b] = q[b - 1] + c[b - 1] as usize;
            }
            if n0 + (c[0] as usize) < m {
                for t in &prev[n0..m] {
                    curr[q[t.c as usize]] = *t;
                    q[t.c as usize] += 1;
                }
                for b in 0..6 {
                    q[b] -= c[b] as usize;
                }
            }
            n0 += c[0] as usize;
            for b in 1..6 {
                if c[b] != 0 {
                    let (qb, cb) = (q[b], c[b] as usize);
                    mr_insert_multi_aux(&mut self.r[b], &mut curr[qb..qb + cb], s, is_comp);
                }
            }
            if n0 == m {
                break;
            }
            let mut ac = [0i64; 6];
            for b in 1..6 {
                // update the intervals to account for buckets ahead
                for (a, acc) in ac.iter_mut().enumerate() {
                    *acc += self.r[b - 1].c[a];
                }
                for t in &mut curr[q[b]..q[b] + c[b] as usize] {
                    t.l = t.l.wrapping_add(ac[t.c as usize] as u64);
                    t.u = t.u.wrapping_add(ac[t.c as usize] as u64);
                }
            }
            std::mem::swap(&mut curr, &mut prev);
        }
    }

    /// `mr_itr_first` and `mr_itr_next_block`.
    pub(crate) fn for_each_block(&self, mut f: impl FnMut(&[u8])) {
        for r in &self.r {
            r.for_each_block(&mut f);
        }
    }
}

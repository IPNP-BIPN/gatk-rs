//! `kmer.h` and `htab.c`: the k-mer arithmetic of BFC and its counting hash.

use crate::khash::{KhOps, Khash};

const BFC_CH_KEYBITS: i32 = 50;
const BFC_CH_MAXPRE: i32 = 20;

/// `bfc_kmer_t`: the k-mer as two bit planes, and its reverse complement as two more.
pub(crate) type BfcKmer = [u64; 4];

/// `bfc_kmer_null`.
pub(crate) const BFC_KMER_NULL: BfcKmer = [0, 0, 0, 0];

/// `bfc_kmer_append`.
#[inline]
pub(crate) fn bfc_kmer_append(k: i32, x: &mut BfcKmer, c: i32) {
    let mask = (1u64 << k) - 1;
    x[0] = ((x[0] << 1) | (c & 1) as u64) & mask;
    x[1] = ((x[1] << 1) | (c >> 1) as u64) & mask;
    x[2] = (x[2] >> 1) | ((1u64 ^ (c & 1) as u64) << (k - 1));
    x[3] = (x[3] >> 1) | ((1u64 ^ (c >> 1) as u64) << (k - 1));
}

/// `bfc_kmer_change`.
#[inline]
pub(crate) fn bfc_kmer_change(k: i32, x: &mut BfcKmer, d: i32, c: i32) {
    let mut t = !(1u64 << d);
    x[0] = (((c & 1) as u64) << d) | (x[0] & t);
    x[1] = (((c >> 1) as u64) << d) | (x[1] & t);
    t = !(1u64 << (k - 1 - d));
    x[2] = ((1u64 ^ (c & 1) as u64) << (k - 1 - d)) | (x[2] & t);
    x[3] = ((1u64 ^ (c >> 1) as u64) << (k - 1 - d)) | (x[3] & t);
}

/// `bfc_hash_64`.
#[inline]
fn bfc_hash_64(mut key: u64, mask: u64) -> u64 {
    key = (!key).wrapping_add(key << 21) & mask;
    key ^= key >> 24;
    key = key.wrapping_add(key << 3).wrapping_add(key << 8) & mask;
    key ^= key >> 14;
    key = key.wrapping_add(key << 2).wrapping_add(key << 4) & mask;
    key ^= key >> 28;
    key = key.wrapping_add(key << 31) & mask;
    key
}

/// `bfc_kmer_hash`, the `h` it fills.
#[inline]
pub(crate) fn bfc_kmer_hash(k: i32, x: &BfcKmer) -> [u64; 2] {
    let t = k >> 1;
    let u = usize::from(((x[1] >> t) & 1) > ((x[3] >> t) & 1)); // the middle base always differs
    let mask = (1u64 << k) - 1;
    let h0 = bfc_hash_64(x[u << 1].wrapping_add(x[(u << 1) | 1]) & mask, mask);
    let h1 = bfc_hash_64(h0 ^ x[(u << 1) | 1], mask);
    [h0.wrapping_add(h1) & mask, h1]
}

/// `_cnt_hash` and `_cnt_eq`: the low 14 bits of a key are its counts.
struct CntOps;

impl KhOps for CntOps {
    fn hash(key: u64) -> u32 {
        (key >> 14) as u32
    }
    fn equal(a: u64, b: u64) -> bool {
        a >> 14 == b >> 14
    }
}

type CntHash = Khash<CntOps, ()>;

/// `bfc_ch_t`. A million empty sub-tables are allocated per pass in C; here they arrive on first
/// insertion, which nothing can see.
pub(crate) struct BfcCh {
    k: i32,
    l_pre: i32,
    h: Vec<Option<Box<CntHash>>>,
}

impl BfcCh {
    /// `bfc_ch_init`.
    pub(crate) fn new(k: i32, mut l_pre: i32) -> BfcCh {
        c_assert!(k <= 63, "bfc_ch_init: k <= 63");
        if k * 2 - l_pre > BFC_CH_KEYBITS {
            l_pre = k * 2 - BFC_CH_KEYBITS;
        }
        if l_pre > BFC_CH_MAXPRE {
            l_pre = BFC_CH_MAXPRE;
        }
        c_assert!(
            k - l_pre < BFC_CH_KEYBITS,
            "bfc_ch_init: k - l_pre < BFC_CH_KEYBITS"
        );
        let mut h = Vec::new();
        h.resize_with(1usize << l_pre, || None);
        BfcCh { k, l_pre, h }
    }

    /// `get_subhash`.
    #[inline]
    fn subhash(&self, x: &[u64; 2]) -> (usize, u64) {
        if self.k <= 32 {
            let t = self.k * 2 - self.l_pre;
            let z = (x[0] << self.k) | x[1];
            let key = ((z & ((1u64 << t) - 1)) << 14) | 1;
            ((z >> t) as usize, key)
        } else {
            let t = self.k - self.l_pre;
            let shift = if t + self.k < BFC_CH_KEYBITS {
                self.k
            } else {
                BFC_CH_KEYBITS - t
            };
            let key = ((((x[0] & ((1u64 << t) - 1)) << shift) ^ x[1]) << 14) | 1;
            ((x[0] >> t) as usize, key)
        }
    }

    /// `bfc_ch_insert`. With one thread the sub-table lock never fails, so the return value and
    /// `forced` have no case left.
    pub(crate) fn insert(&mut self, x: &[u64; 2], is_high: bool) {
        let (i, key) = self.subhash(x);
        let h = self.h[i].get_or_insert_with(|| Box::new(CntHash::new()));
        let (k, ret) = h.put(key);
        let slot = &mut h.keys[k as usize];
        if ret != 0 {
            if is_high {
                *slot |= 1 << 8;
            }
        } else {
            if (*slot & 0xff) != 0xff {
                *slot += 1;
            }
            if is_high && ((*slot >> 8) & 0x3f) != 0x3f {
                *slot += 1 << 8;
            }
        }
    }

    /// `bfc_ch_kmer_occ`, through `bfc_ch_get`.
    #[inline]
    pub(crate) fn kmer_occ(&self, z: &BfcKmer) -> i32 {
        let x = bfc_kmer_hash(self.k, z);
        let (i, key) = self.subhash(&x);
        let Some(h) = &self.h[i] else {
            return -1;
        };
        let itr = h.get(key);
        if itr == h.end() {
            -1
        } else {
            (h.keys[itr as usize] & 0x3fff) as i32
        }
    }

    /// `bfc_ch_hist`: `cnt`, and the mode the function returns.
    pub(crate) fn hist(&self) -> ([u64; 256], i32) {
        let mut cnt = [0u64; 256];
        for h in self.h.iter().flatten() {
            for k in 0..h.end() {
                if h.exist(k) {
                    cnt[(h.keys[k as usize] & 0xff) as usize] += 1;
                }
            }
        }
        let mut max = 0u64;
        let mut max_i = -1;
        for (i, &c) in cnt.iter().enumerate().skip(3) {
            if c > max {
                max = c;
                max_i = i as i32;
            }
        }
        (cnt, max_i)
    }
}

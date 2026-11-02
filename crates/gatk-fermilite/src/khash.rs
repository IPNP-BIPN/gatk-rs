//! `khash.h` 0.2.8: quadratic probing, two flag bits per bucket, kick-out rehash.
//!
//! Nothing in fermi-lite iterates a table where order could show: `bfc_ch_hist` only adds up, and
//! every other table is looked up by key. The hash and probing are ported regardless, so the keys
//! alias exactly when they do in C (`_cnt_eq` compares only the top 50 bits).

/// `__hash_func` and `__hash_equal` of one `KHASH_INIT`.
pub(crate) trait KhOps {
    fn hash(key: u64) -> u32;
    fn equal(a: u64, b: u64) -> bool;
}

/// `kh_int64_hash_func` and `kh_int64_hash_equal`.
pub(crate) struct Int64Ops;

impl KhOps for Int64Ops {
    fn hash(key: u64) -> u32 {
        ((key >> 33) ^ key ^ (key << 11)) as u32
    }
    fn equal(a: u64, b: u64) -> bool {
        a == b
    }
}

fn isempty(f: &[u32], i: u32) -> bool {
    (f[(i >> 4) as usize] >> ((i & 0xf) << 1)) & 2 != 0
}
fn isdel(f: &[u32], i: u32) -> bool {
    (f[(i >> 4) as usize] >> ((i & 0xf) << 1)) & 1 != 0
}
fn iseither(f: &[u32], i: u32) -> bool {
    (f[(i >> 4) as usize] >> ((i & 0xf) << 1)) & 3 != 0
}
fn set_isempty_false(f: &mut [u32], i: u32) {
    f[(i >> 4) as usize] &= !(2u32 << ((i & 0xf) << 1));
}
fn set_isboth_false(f: &mut [u32], i: u32) {
    f[(i >> 4) as usize] &= !(3u32 << ((i & 0xf) << 1));
}
fn set_isdel_true(f: &mut [u32], i: u32) {
    f[(i >> 4) as usize] |= 1u32 << ((i & 0xf) << 1);
}
fn fsize(m: u32) -> usize {
    if m < 16 {
        1
    } else {
        (m >> 4) as usize
    }
}

/// `kh_##name##_t` with 64-bit keys; `V = ()` is a set.
pub(crate) struct Khash<O: KhOps, V: Copy + Default> {
    n_buckets: u32,
    size: u32,
    n_occupied: u32,
    flags: Vec<u32>,
    pub(crate) keys: Vec<u64>,
    pub(crate) vals: Vec<V>,
    _ops: std::marker::PhantomData<O>,
}

impl<O: KhOps, V: Copy + Default> Khash<O, V> {
    /// `kh_init`.
    pub(crate) fn new() -> Self {
        Khash {
            n_buckets: 0,
            size: 0,
            n_occupied: 0,
            flags: Vec::new(),
            keys: Vec::new(),
            vals: Vec::new(),
            _ops: std::marker::PhantomData,
        }
    }

    /// `kh_end`.
    pub(crate) fn end(&self) -> u32 {
        self.n_buckets
    }

    /// `kh_size`.
    pub(crate) fn size(&self) -> u32 {
        self.size
    }

    /// `kh_exist`.
    pub(crate) fn exist(&self, x: u32) -> bool {
        !iseither(&self.flags, x)
    }

    /// `kh_get_##name`.
    pub(crate) fn get(&self, key: u64) -> u32 {
        if self.n_buckets != 0 {
            let mask = self.n_buckets - 1;
            let mut step = 0u32;
            let k = O::hash(key);
            let mut i = k & mask;
            let last = i;
            while !isempty(&self.flags, i)
                && (isdel(&self.flags, i) || !O::equal(self.keys[i as usize], key))
            {
                step = step.wrapping_add(1);
                i = i.wrapping_add(step) & mask;
                if i == last {
                    return self.n_buckets;
                }
            }
            if iseither(&self.flags, i) {
                self.n_buckets
            } else {
                i
            }
        } else {
            0
        }
    }

    /// `kh_resize_##name`.
    fn resize(&mut self, mut new_n_buckets: u32) {
        // kroundup32
        new_n_buckets = new_n_buckets.wrapping_sub(1);
        new_n_buckets |= new_n_buckets >> 1;
        new_n_buckets |= new_n_buckets >> 2;
        new_n_buckets |= new_n_buckets >> 4;
        new_n_buckets |= new_n_buckets >> 8;
        new_n_buckets |= new_n_buckets >> 16;
        new_n_buckets = new_n_buckets.wrapping_add(1);
        if new_n_buckets < 4 {
            new_n_buckets = 4;
        }
        if self.size >= (new_n_buckets >> 1) + (new_n_buckets >> 2) {
            return; // the requested size is too small
        }
        let mut new_flags = vec![0xaaaa_aaaau32; fsize(new_n_buckets)];
        if self.n_buckets < new_n_buckets {
            self.keys.resize(new_n_buckets as usize, 0);
            self.vals.resize(new_n_buckets as usize, V::default());
        }
        let new_mask = new_n_buckets - 1;
        for j in 0..self.n_buckets {
            if !iseither(&self.flags, j) {
                let mut key = self.keys[j as usize];
                let mut val = self.vals[j as usize];
                set_isdel_true(&mut self.flags, j);
                loop {
                    // the kick-out process
                    let mut step = 0u32;
                    let k = O::hash(key);
                    let mut i = k & new_mask;
                    while !isempty(&new_flags, i) {
                        step = step.wrapping_add(1);
                        i = i.wrapping_add(step) & new_mask;
                    }
                    set_isempty_false(&mut new_flags, i);
                    if i < self.n_buckets && !iseither(&self.flags, i) {
                        std::mem::swap(&mut self.keys[i as usize], &mut key);
                        std::mem::swap(&mut self.vals[i as usize], &mut val);
                        set_isdel_true(&mut self.flags, i);
                    } else {
                        self.keys[i as usize] = key;
                        self.vals[i as usize] = val;
                        break;
                    }
                }
            }
        }
        if self.n_buckets > new_n_buckets {
            self.keys.truncate(new_n_buckets as usize);
            self.vals.truncate(new_n_buckets as usize);
        }
        self.flags = new_flags;
        self.n_buckets = new_n_buckets;
        self.n_occupied = self.size;
    }

    /// `kh_put_##name`: the bucket, and `*ret` (0 when the key was already present).
    pub(crate) fn put(&mut self, key: u64) -> (u32, i32) {
        if self.n_occupied >= (self.n_buckets >> 2) + (self.n_buckets >> 1) {
            if self.n_buckets > (self.size << 1) {
                self.resize(self.n_buckets - 1); // clear the "deleted" elements
            } else {
                self.resize(self.n_buckets + 1); // expand
            }
        }
        let mask = self.n_buckets - 1;
        let mut step = 0u32;
        let mut x = self.n_buckets;
        let mut site = self.n_buckets;
        let k = O::hash(key);
        let mut i = k & mask;
        if isempty(&self.flags, i) {
            x = i;
        } else {
            let last = i;
            while !isempty(&self.flags, i)
                && (isdel(&self.flags, i) || !O::equal(self.keys[i as usize], key))
            {
                if isdel(&self.flags, i) {
                    site = i;
                }
                step = step.wrapping_add(1);
                i = i.wrapping_add(step) & mask;
                if i == last {
                    x = site;
                    break;
                }
            }
            if x == self.n_buckets {
                if isempty(&self.flags, i) && site != self.n_buckets {
                    x = site;
                } else {
                    x = i;
                }
            }
        }
        let ret;
        if isempty(&self.flags, x) {
            self.keys[x as usize] = key;
            set_isboth_false(&mut self.flags, x);
            self.size += 1;
            self.n_occupied += 1;
            ret = 1;
        } else if isdel(&self.flags, x) {
            self.keys[x as usize] = key;
            set_isboth_false(&mut self.flags, x);
            self.size += 1;
            ret = 2;
        } else {
            ret = 0;
        }
        (x, ret)
    }

    /// `kh_del_##name`.
    pub(crate) fn del(&mut self, x: u32) {
        if x != self.n_buckets && !iseither(&self.flags, x) {
            set_isdel_true(&mut self.flags, x);
            self.size -= 1;
        }
    }
}

/// `KHASH_INIT2(64,, khint64_t, uint64_t, 1, kh_int64_hash_func, kh_int64_hash_equal)`.
pub(crate) type Hash64 = Khash<Int64Ops, u64>;

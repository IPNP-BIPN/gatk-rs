//! `ksort.h`: the sorts and heap steps fermi-lite instantiates, comparison for comparison.
//!
//! Every comparator here has ties (`vlt1` on nsr and length, `infocmp` on one field, the error
//! correction heap on the total penalty), so the arrangement of equal keys is output. These are the
//! macros, not something merely equivalent; `slice::sort` is not. `ks_shuffle` and `ks_sample` are
//! the only callers of `drand48` in the library and nothing calls them, so they are not here.

/// `ks_heapdown_##name(i, n, l)`.
pub(crate) fn ks_heapdown<T: Copy>(mut i: usize, l: &mut [T], lt: impl Fn(&T, &T) -> bool) {
    let n = l.len();
    if n == 0 {
        // C reads and writes back the slot just past the live heap, which changes nothing.
        return;
    }
    let tmp = l[i];
    let mut k = i;
    loop {
        k = (k << 1) + 1;
        if k >= n {
            break;
        }
        if k != n - 1 && lt(&l[k], &l[k + 1]) {
            k += 1;
        }
        if lt(&l[k], &tmp) {
            break;
        }
        l[i] = l[k];
        i = k;
    }
    l[i] = tmp;
}

/// `ks_heapup_##name(n, l)`.
pub(crate) fn ks_heapup<T: Copy>(l: &mut [T], lt: impl Fn(&T, &T) -> bool) {
    let mut k = l.len() - 1;
    let tmp = l[k];
    while k != 0 {
        let i = (k - 1) >> 1;
        if lt(&tmp, &l[i]) {
            break;
        }
        l[k] = l[i];
        k = i;
    }
    l[k] = tmp;
}

/// `__ks_insertsort_##name(s, t)`.
fn ks_insertsort<T: Copy>(a: &mut [T], lt: &impl Fn(&T, &T) -> bool) {
    for i in 1..a.len() {
        let mut j = i;
        while j > 0 && lt(&a[j], &a[j - 1]) {
            a.swap(j, j - 1);
            j -= 1;
        }
    }
}

/// `ks_combsort_##name(n, a)`.
fn ks_combsort<T: Copy>(a: &mut [T], lt: &impl Fn(&T, &T) -> bool) {
    const SHRINK_FACTOR: f64 = 1.2473309501039786540366528676643;
    let n = a.len();
    let mut gap = n;
    loop {
        if gap > 2 {
            gap = (gap as f64 / SHRINK_FACTOR) as usize;
            if gap == 9 || gap == 10 {
                gap = 11;
            }
        }
        let mut do_swap = false;
        for i in 0..n - gap {
            let j = i + gap;
            if lt(&a[j], &a[i]) {
                a.swap(i, j);
                do_swap = true;
            }
        }
        if !(do_swap || gap > 2) {
            break;
        }
    }
    if gap != 1 {
        ks_insertsort(a, lt);
    }
}

/// `ks_introsort_##name(n, a)`: median of three over positions `s`, `s + (t-s)/2 + 1` and `t`, the
/// larger side pushed, a side of at most 16 left to the final insertion sort, combsort at depth 0.
pub(crate) fn ks_introsort<T: Copy>(a: &mut [T], lt: impl Fn(&T, &T) -> bool) {
    let n = a.len();
    if n < 1 {
        return;
    } else if n == 2 {
        if lt(&a[1], &a[0]) {
            a.swap(0, 1);
        }
        return;
    }
    let mut d: i32 = 2;
    while (1usize << d) < n {
        d += 1;
    }
    let mut stack: Vec<(isize, isize, i32)> = Vec::new();
    let mut s: isize = 0;
    let mut t: isize = n as isize - 1;
    d <<= 1;
    loop {
        if s < t {
            d -= 1;
            if d == 0 {
                ks_combsort(&mut a[s as usize..=t as usize], &lt);
                t = s;
                continue;
            }
            let mut i = s;
            let mut j = t;
            let mut k = i + ((j - i) >> 1) + 1;
            if lt(&a[k as usize], &a[i as usize]) {
                if lt(&a[k as usize], &a[j as usize]) {
                    k = j;
                }
            } else {
                k = if lt(&a[j as usize], &a[i as usize]) {
                    i
                } else {
                    j
                };
            }
            let rp = a[k as usize];
            if k != t {
                a.swap(k as usize, t as usize);
            }
            loop {
                loop {
                    i += 1;
                    if !lt(&a[i as usize], &rp) {
                        break;
                    }
                }
                loop {
                    j -= 1;
                    if !(i <= j && lt(&rp, &a[j as usize])) {
                        break;
                    }
                }
                if j <= i {
                    break;
                }
                a.swap(i as usize, j as usize);
            }
            a.swap(i as usize, t as usize);
            if i - s > t - i {
                if i - s > 16 {
                    stack.push((s, i - 1, d));
                }
                s = if t - i > 16 { i + 1 } else { t };
            } else {
                if t - i > 16 {
                    stack.push((i + 1, t, d));
                }
                t = if i - s > 16 { i - 1 } else { s };
            }
        } else {
            match stack.pop() {
                None => {
                    ks_insertsort(a, &lt);
                    return;
                }
                Some(top) => {
                    s = top.0;
                    t = top.1;
                    d = top.2;
                }
            }
        }
    }
}

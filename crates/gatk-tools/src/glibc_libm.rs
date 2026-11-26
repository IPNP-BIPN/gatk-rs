//! glibc 2.35's `exp`, `log`, `log1p` and `expm1` as an x86-64 without FMA or AVX2 runs them:
//! the SSE2 builds that `GLIBC_TUNABLES=glibc.cpu.hwcaps=-AVX2,-FMA,...` leaves the ifunc to
//! pick, which is what the frozen gCNV reference runs (`docs/gcnv-reproducible-oracle.md`).
//! numpy's float64 `exp`, `log`, `log1p`, `expm1` and pytensor's generated C both call these.
//!
//! `exp` and `log` are ARM's optimized-routines (MIT) as glibc compiles them without
//! `__FP_FAST_FMA`; `log1p` and `expm1` are fdlibm's, which have no ifunc variant. None is
//! correctly rounded, so each is the algorithm, operation for operation, with the tables glibc
//! carries. Checked against numpy in the image on 1.2 million points each
//! (`tools/python-tools/probe_gcnv_libm.py`).
//!
//! # Reading this file without reading Rust
//!
//! - `x.to_bits()` is Java's `Double.doubleToRawLongBits(x)` and `f64::from_bits(b)` is
//!   `Double.longBitsToDouble(b)`; the C reaches the same bits through a union. A constant
//!   written `f64::from_bits(0x...)` is the C's hexadecimal float literal, bit for bit;
//! - `u32`/`u64` are unsigned 32/64-bit integers, `i32`/`i64` signed ones. `a as u32` keeps the
//!   low 32 bits, and `u32 as i32` reinterprets the same bits as signed, as the C's
//!   `int32_t`/`uint32_t` casts do;
//! - `wrapping_sub` / `wrapping_add` are C's unsigned arithmetic modulo 2^32 or 2^64. Plain `-` on
//!   an unsigned Rust integer aborts on underflow in a debug build, so wherever the C relies on
//!   the wrap, the Rust says so;
//! - `f64 as i32` truncates toward zero like the C cast, and saturates where C is undefined; the
//!   one use here (`expm1`'s `k`) is bounded well inside `i32`;
//! - `|i| log_data(2 + i)` is a closure, Java's `i -> logData(2 + i)`;
//! - `0x8000_0000`: the underscore is a digit separator and changes nothing.

// fdlibm's constants as glibc's source spells them.
#![allow(clippy::excessive_precision)]

use crate::libm_pow_two::EXP_TABLE;

/// The top 12 bits of `x`: the sign and the 11-bit exponent, as glibc's `top12` helper.
fn top12(x: f64) -> u32 {
    (x.to_bits() >> 52) as u32
}

/// `__exp`, `e_exp.c` with `TOINT_INTRINSICS` and `EXP_USE_TOINT_NARROW` off.
///
/// **How**: with N = 128, `k = round(x * N / ln 2)` and `r = x - k * ln 2 / N`, so
/// `exp(x) = 2^(k/N) * exp(r)`. `EXP_TABLE` holds, for each `k mod N`, a pair: a small `tail`
/// correction and the bits of `2^(i/N)`; adding `k << 45` to those bits puts `k / N` into the
/// exponent field. `exp(r) - 1` is a degree-5 polynomial, and the result is
/// `scale + scale * tmp`.
///
/// **Why** each operation is kept in its place: the function is not correctly rounded, so only
/// the same operations in the same order give numpy's bits.
pub fn exp(x: f64) -> f64 {
    const INV_LN2_N: f64 = f64::from_bits(0x40671547652b82fe);
    const NEG_LN2_HI_N: f64 = f64::from_bits(0xbf762e42fefa0000);
    const NEG_LN2_LO_N: f64 = f64::from_bits(0xbd0cf79abc9e3b3a);
    const SHIFT: f64 = f64::from_bits(0x4338000000000000);
    const C2: f64 = f64::from_bits(0x3fdffffffffffdbd);
    const C3: f64 = f64::from_bits(0x3fc555555555543c);
    const C4: f64 = f64::from_bits(0x3fa55555cf172b91);
    const C5: f64 = f64::from_bits(0x3f81111167a4d017);

    // `abstop` is the biased exponent of |x|. One unsigned compare, through the wrap of
    // `wrapping_sub`, catches both |x| < 2^-54 (0x3c9 is that exponent) and |x| >= 512.
    let mut abstop = top12(x) & 0x7ff;
    if abstop.wrapping_sub(top12(f64::from_bits(0x3c90000000000000)))
        >= top12(512.0) - top12(f64::from_bits(0x3c90000000000000))
    {
        // The subtraction wrapped: |x| is tiny and the answer is `1 + x`.
        if abstop.wrapping_sub(top12(f64::from_bits(0x3c90000000000000))) >= 0x8000_0000 {
            return 1.0 + x;
        }
        // |x| >= 1024, infinite or NaN: -inf gives 0, inf and NaN give `1 + x`, a finite
        // value overflows to inf or underflows to 0. glibc reaches the last two through its
        // overflow and underflow helpers, which also set `errno`; the value is the same.
        if abstop >= top12(1024.0) {
            if x.to_bits() == f64::NEG_INFINITY.to_bits() {
                return 0.0;
            }
            if abstop >= top12(f64::INFINITY) {
                return 1.0 + x;
            }
            return if x.to_bits() >> 63 != 0 {
                0.0
            } else {
                f64::INFINITY
            };
        }
        // 512 <= |x| < 1024: the result may leave the normal range, so `abstop = 0` routes the
        // end of the computation through `exp_special_case`.
        abstop = 0;
    }
    // Adding `SHIFT` (1.5 * 2^52) rounds `z` to an integer held in the low bits of `kd`, so
    // `ki` is `k` and `kd - SHIFT` is `k` as a double. `ln 2 / N` is split in two (`HI` + `LO`)
    // so that `r` keeps its low bits.
    let z = INV_LN2_N * x;
    let kd = z + SHIFT;
    let ki = kd.to_bits();
    let kd = kd - SHIFT;
    let r = x + kd * NEG_LN2_HI_N + kd * NEG_LN2_LO_N;
    // `ki % 128` is `k mod N` (the table pair), and the high bits of `ki` shifted by 45 land in
    // the exponent field; `wrapping_add` is the C's `uint64_t` addition.
    let idx = 2 * (ki % 128) as usize;
    let top = ki << 45;
    let tail = f64::from_bits(EXP_TABLE[idx]);
    let sbits = EXP_TABLE[idx + 1].wrapping_add(top);
    let r2 = r * r;
    let tmp = tail + r + r2 * (C2 + r * C3) + r2 * r2 * (C4 + r * C5);
    if abstop == 0 {
        return exp_special_case(tmp, sbits, ki);
    }
    let scale = f64::from_bits(sbits);
    scale + scale * tmp
}

/// glibc's `specialcase` for `exp`, when `2^(k/N)` alone would overflow or underflow.
///
/// **How**: a positive `k` (bit 31 of `ki` clear) scales down by 2^1009 and multiplies the result
/// back by 2^1009 (`0x7f00000000000000`). A negative one scales up by 2^1022 and multiplies by
/// 2^-1022 (`0x0010000000000000`) at the end. When the scaled result is below one, the final
/// multiply lands in the subnormal range and would round twice, so the sum is first redone
/// around 1.0 to round once.
fn exp_special_case(tmp: f64, sbits: u64, ki: u64) -> f64 {
    if ki & 0x8000_0000 == 0 {
        let scale = f64::from_bits(sbits.wrapping_sub(1009u64 << 52));
        return f64::from_bits(0x7f00000000000000) * (scale + scale * tmp);
    }
    let scale = f64::from_bits(sbits.wrapping_add(1022u64 << 52));
    let mut y = scale + scale * tmp;
    if y < 1.0 {
        let lo = scale - y + scale * tmp;
        let hi = 1.0 + y;
        let lo = 1.0 - hi + y + lo;
        y = (hi + lo) - 1.0;
        // `-0.0 == 0.0` is true, so this turns a negative zero into a positive one.
        if y == 0.0 {
            y = 0.0;
        }
    }
    f64::from_bits(0x0010000000000000) * y
}

// Pastes `LOG_DATA`, glibc's `__log_data` table, into this file (C's `#include`).
include!("glibc_log_data.rs");

/// Entry `i` of `__log_data` as a double; the index layout is described on `LOG_DATA`.
fn log_data(i: usize) -> f64 {
    f64::from_bits(LOG_DATA[i])
}

/// `__log`, `e_log.c` without `__FP_FAST_FMA`: the reduction goes through `tab2`.
///
/// **How**: close to 1 a degree-11 polynomial (`poly1`, here `b`) in `r = x - 1`. Elsewhere
/// `x = 2^k * z` with `z` in about [0.69, 1.38), the table entry `i` from the top mantissa bits
/// of `z` gives `1/c` and `log(c)`, and `log(x) = k ln 2 + log(c) + log1p(r)` with
/// `r = z/c - 1` and a degree-5 polynomial (`poly`, here `a`).
///
/// **Why `tab2`**: with an FMA, glibc computes `r` as one fused `z * invc - 1`. Without it, as
/// in the reference image, `r = (z - chi - clo) * invc`, with `chi` and `clo` from `tab2`; that
/// is the form here.
pub fn log(x: f64) -> f64 {
    const OFF: u64 = 0x3fe6000000000000;
    let ln2hi = log_data(0);
    let ln2lo = log_data(1);
    let a = |i: usize| log_data(2 + i);
    let b = |i: usize| log_data(7 + i);
    // The near-one window: `x` in [1 - 2^-4, 1 + 0x1.09p-4), tested with one unsigned compare
    // on the bit patterns.
    let mut ix = x.to_bits();
    let top = (ix >> 48) as u32;
    let lo_bound = (1.0 - f64::from_bits(0x3fb0000000000000)).to_bits();
    let hi_bound = (1.0 + f64::from_bits(0x3fb0900000000000)).to_bits();
    if ix.wrapping_sub(lo_bound) < hi_bound - lo_bound {
        if ix == 1.0f64.to_bits() {
            return 0.0;
        }
        let r = x - 1.0;
        let r2 = r * r;
        let r3 = r * r2;
        let mut y = r3
            * (b(1)
                + r * b(2)
                + r2 * b(3)
                + r3 * (b(4)
                    + r * b(5)
                    + r2 * b(6)
                    + r3 * (b(7) + r * b(8) + r2 * b(9) + r3 * b(10))));
        // Splits `r` into `rhi`, its high half (multiply by 2^27, add, subtract), and `rlo`, so
        // `rhi * rhi` is exact without an FMA.
        let w = r * f64::from_bits(0x41a0000000000000);
        let rhi = r + w - w;
        let rlo = r - rhi;
        let w = rhi * rhi * b(0);
        let hi = r + w;
        let mut lo = r - hi + w;
        lo += b(0) * rlo * (rhi + r);
        y += lo;
        y += hi;
        return y;
    }
    // One unsigned compare catches a subnormal, zero, negative, infinite or NaN `x`.
    if top.wrapping_sub(0x0010) >= 0x7ff0 - 0x0010 {
        if ix << 1 == 0 {
            return f64::NEG_INFINITY;
        }
        if ix == f64::INFINITY.to_bits() {
            return x;
        }
        if top & 0x8000 != 0 || top & 0x7ff0 == 0x7ff0 {
            return f64::NAN;
        }
        // A subnormal: multiply by 2^52 to normalise it and take 52 back off the exponent.
        ix = (x * f64::from_bits(0x4330000000000000)).to_bits();
        ix = ix.wrapping_sub(52u64 << 52);
    }
    // `tmp` is `x`'s bits less those of about 0.69: bits 45 to 51 pick the table entry, and the
    // arithmetic shift `(tmp as i64) >> 52` (signed, so it keeps the sign) gives `k`.
    let tmp = ix.wrapping_sub(OFF);
    let i = ((tmp >> 45) % 128) as usize;
    let k = (tmp as i64) >> 52;
    let iz = ix.wrapping_sub(tmp & (0xfffu64 << 52));
    let invc = log_data(18 + 2 * i);
    let logc = log_data(18 + 2 * i + 1);
    let chi = log_data(274 + 2 * i);
    let clo = log_data(274 + 2 * i + 1);
    let z = f64::from_bits(iz);
    let r = (z - chi - clo) * invc;
    let kd = k as f64;
    let w = kd * ln2hi + logc;
    let hi = w + r;
    let lo = w - hi + r + kd * ln2lo;
    let r2 = r * r;
    lo + r2 * a(0) + r * r2 * (a(1) + r * a(2) + r2 * (a(3) + r * a(4))) + hi
}

/// fdlibm's `GET_HIGH_WORD`: the sign, exponent and top 20 mantissa bits of `x`.
fn high_word(x: f64) -> u32 {
    (x.to_bits() >> 32) as u32
}

/// fdlibm's `SET_HIGH_WORD`: `x` with its upper 32 bits replaced by `high`.
fn with_high_word(x: f64, high: u32) -> f64 {
    f64::from_bits((u64::from(high) << 32) | (x.to_bits() & 0xffff_ffff))
}

/// `__log1p`, glibc's `s_log1p.c` (fdlibm with its polynomial split four ways).
///
/// **How**: write `1 + x = 2^k * (1 + f)` with `1 + f` in [sqrt(2)/2, sqrt(2)), keep `c`, the
/// rounding error of forming `1 + x`, and approximate `log(1 + f)` through `s = f / (2 + f)` and
/// an even polynomial `R(z)` in `z = s * s`. For `x` already in about (-0.29, 0.41) there is no
/// reduction (`k = 0`, `f = x`).
///
/// **Why the split**: fdlibm writes `R` as one Horner chain; glibc evaluates it as four partial
/// sums `r1 + z2 * r2 + z4 * r3 + z6 * r4`, which rounds differently, and that is the form numpy
/// gets in the image.
pub fn log1p(x: f64) -> f64 {
    const LN2_HI: f64 = 6.93147180369123816490e-01;
    const LN2_LO: f64 = 1.90821492927058770002e-10;
    const LP: [f64; 8] = [
        0.0,
        6.666666666666735130e-01,
        3.999999999940941908e-01,
        2.857142874366239149e-01,
        2.222219843214978396e-01,
        1.818357216161805012e-01,
        1.531383769920937332e-01,
        1.479819860511658591e-01,
    ];
    let hx = high_word(x) as i32;
    let ax = hx & 0x7fff_ffff;
    let mut k = 1i32;
    let mut f = 0.0;
    let mut c = 0.0;
    let mut hu: i32 = 0;
    // `hx` below the high word of sqrt(2) - 1, which includes every negative `x`.
    if hx < 0x3FDA827A {
        // x <= -1: -1 gives -inf, anything below gives NaN.
        if ax >= 0x3ff00000 {
            return if x == -1.0 {
                f64::NEG_INFINITY
            } else {
                f64::NAN
            };
        }
        // |x| < 2^-29: `x` itself below 2^-54, else `x - x^2/2`.
        if ax < 0x3e200000 {
            if ax < 0x3c900000 {
                return x;
            }
            return x - x * x * 0.5;
        }
        // fdlibm's `(int)0xbfd2bec3`: as signed high words, a negative `x` of smaller magnitude
        // has the smaller word, so this keeps x > 0 and -0.2929 < x < 0 unreduced. `hu = 1` only
        // marks that the `hu == 0` shortcut below does not apply.
        if hx > 0 || hx <= 0xbfd2bec3u32 as i32 {
            k = 0;
            f = x;
            hu = 1;
        }
    } else if hx >= 0x7ff00000 {
        // Infinite or NaN.
        return x + x;
    }
    if k != 0 {
        let mut u;
        // Below 2^53, `1 + x` rounds, and `c / u` is the part it lost; above, `1 + x` is `x`.
        if hx < 0x43400000 {
            u = 1.0 + x;
            hu = high_word(u) as i32;
            k = (hu >> 20) - 1023;
            c = if k > 0 { 1.0 - (u - x) } else { x - (u - 1.0) };
            c /= u;
        } else {
            u = x;
            hu = high_word(u) as i32;
            k = (hu >> 20) - 1023;
            c = 0.0;
        }
        // Keep the mantissa's top 20 bits and move `u` into [sqrt(2)/2, sqrt(2)) by setting its
        // exponent (0x6a09e is sqrt(2)'s top mantissa bits). In the upper case `hu` becomes a
        // value that is zero only when `f` is tiny, which the `hu == 0` test below reads.
        hu &= 0x000fffff;
        if hu < 0x6a09e {
            u = with_high_word(u, (hu | 0x3ff00000) as u32);
        } else {
            k += 1;
            u = with_high_word(u, (hu | 0x3fe00000) as u32);
            hu = (0x00100000 - hu) >> 2;
        }
        f = u - 1.0;
    }
    let hfsq = 0.5 * f * f;
    let kd = f64::from(k);
    // |f| < 2^-20: a short series is enough.
    if hu == 0 {
        if f == 0.0 {
            if k == 0 {
                return 0.0;
            }
            c += kd * LN2_LO;
            return kd * LN2_HI + c;
        }
        let r = hfsq * (1.0 - 0.66666666666666666 * f);
        if k == 0 {
            return f - r;
        }
        return kd * LN2_HI - ((r - (kd * LN2_LO + c)) - f);
    }
    let s = f / (2.0 + f);
    let z = s * s;
    let r1 = z * LP[1];
    let z2 = z * z;
    let r2 = LP[2] + z * LP[3];
    let z4 = z2 * z2;
    let r3 = LP[4] + z * LP[5];
    let z6 = z4 * z2;
    let r4 = LP[6] + z * LP[7];
    let r = r1 + z2 * r2 + z4 * r3 + z6 * r4;
    if k == 0 {
        return f - (hfsq - s * (hfsq + r));
    }
    kd * LN2_HI - ((hfsq - (s * (hfsq + r) + (kd * LN2_LO + c))) - f)
}

/// `__expm1`, glibc's `s_expm1.c` (fdlibm with its polynomial split three ways).
///
/// **How**: `x = k ln 2 + r` with |r| <= ln 2 / 2, `r` carried as `hi - lo` with the rounding
/// correction `c`. `exp(r) - 1` comes from a rational approximation whose polynomial in
/// `r^2 / 2` is evaluated as three partial sums (`r1 + h2 * r2 + h4 * r3`) rather than fdlibm's
/// single Horner chain. The result is then scaled by 2^k by adding `k` to the exponent field,
/// with separate formulas for `k` of 0, -1, 1, large and small, as fdlibm has them.
pub fn expm1(x: f64) -> f64 {
    const O_THRESHOLD: f64 = 7.09782712893383973096e+02;
    const LN2_HI: f64 = 6.93147180369123816490e-01;
    const LN2_LO: f64 = 1.90821492927058770002e-10;
    const INVLN2: f64 = std::f64::consts::LOG2_E;
    const Q: [f64; 6] = [
        1.0,
        -3.33333333333331316428e-02,
        1.58730158725481460165e-03,
        -7.93650757867487942473e-05,
        4.00821782732936239552e-06,
        -2.01099218183624371326e-07,
    ];
    let mut x = x;
    let mut hx = high_word(x);
    let negative = hx & 0x8000_0000 != 0;
    hx &= 0x7fff_ffff;
    // |x| >= 56 ln 2: the overflow and the NaN/infinity cases, and a negative `x` whose answer
    // is -1 to double precision.
    if hx >= 0x4043687A {
        // |x| >= 709.78, about ln(DBL_MAX).
        if hx >= 0x40862E42 {
            if hx >= 0x7ff00000 {
                // A NaN (non-zero mantissa) propagates; -inf gives -1 and +inf itself.
                if (hx & 0xfffff) | (x.to_bits() as u32) != 0 {
                    return x + x;
                }
                return if negative { -1.0 } else { x };
            }
            if x > O_THRESHOLD {
                return f64::INFINITY;
            }
        }
        // fdlibm's `tiny - one`: -1, with the inexact flag raised.
        if negative {
            return 1.0e-300 - 1.0;
        }
    }
    // `let k: i32;` declares `k` without a value; the compiler checks every path below gives it
    // one before it is read. The same holds for `let (hi, lo);`.
    let k: i32;
    let mut c = 0.0;
    // |x| > ln 2 / 2: reduce. Below 1.5 ln 2, `k` is +1 or -1; otherwise it is `x / ln 2`
    // rounded half away from zero.
    if hx > 0x3fd62e42 {
        let (hi, lo);
        if hx < 0x3FF0A2B2 {
            if !negative {
                hi = x - LN2_HI;
                lo = LN2_LO;
                k = 1;
            } else {
                hi = x + LN2_HI;
                lo = -LN2_LO;
                k = -1;
            }
        } else {
            k = (INVLN2 * x + if negative { -0.5 } else { 0.5 }) as i32;
            let t = f64::from(k);
            hi = x - t * LN2_HI;
            lo = t * LN2_LO;
        }
        x = hi - lo;
        c = (hi - x) - lo;
    } else if hx < 0x3c900000 {
        // |x| < 2^-54.
        return x;
    } else {
        k = 0;
    }
    let hfx = 0.5 * x;
    let hxs = x * hfx;
    let r1 = 1.0 + hxs * Q[1];
    let h2 = hxs * hxs;
    let r2 = Q[2] + hxs * Q[3];
    let h4 = h2 * h2;
    let r3 = Q[4] + hxs * Q[5];
    let r1 = r1 + h2 * r2 + h4 * r3;
    let t = 3.0 - r1 * hfx;
    let mut e = hxs * ((r1 - t) / (6.0 - x * t));
    if k == 0 {
        return x - (x * e - hxs);
    }
    e = x * (e - c) - c;
    e -= hxs;
    if k == -1 {
        return 0.5 * (x - e) - 0.5;
    }
    if k == 1 {
        if x < -0.25 {
            return -2.0 * (e - (x + 0.5));
        }
        return 1.0 + 2.0 * (x - e);
    }
    // For these `k`, `1 - (e - x)` loses nothing worth keeping, and the scaling by 2^k goes
    // into the exponent. `k == 1024` would overflow the exponent field, so it multiplies by
    // 2 * 2^1023 instead. `(k << 20) as u32` keeps a negative `k`'s two's complement bits, and
    // `wrapping_add` adds it as the C's unsigned high word does.
    if k <= -2 || k > 56 {
        let mut y = 1.0 - (e - x);
        if k == 1024 {
            y = y * 2.0 * f64::from_bits(0x7fe0000000000000);
        } else {
            y = with_high_word(y, high_word(y).wrapping_add((k << 20) as u32));
        }
        return y - 1.0;
    }
    // 2 <= k <= 56 (k = 1 and k = -1 returned above): `t` is `1 - 2^-k` below 20 and `2^-k`
    // from 20 on, built directly in the exponent field.
    let v = if k < 20 {
        let t = with_high_word(1.0, 0x3ff00000 - (0x200000 >> k));
        t - (e - x)
    } else {
        let t = with_high_word(1.0, ((0x3ff - k) << 20) as u32);
        x - (e + t) + 1.0
    };
    with_high_word(v, high_word(v).wrapping_add((k << 20) as u32))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// numpy's answers in the image, written by `probe_gcnv_libm.py` to `$GCNV_LIBM_VECTORS`.
    ///
    /// Without the variable the test passes having checked nothing (`let Ok(..) else`, the
    /// early return): the vectors are not committed, so this runs only where they were made.
    fn against_numpy(name: &str, f: fn(f64) -> f64) {
        let Ok(dir) = std::env::var("GCNV_LIBM_VECTORS") else {
            return;
        };
        let read = |suffix: &str| -> Vec<f64> {
            std::fs::read(format!("{dir}/{name}.{suffix}"))
                .unwrap()
                .chunks(8)
                .map(|b| f64::from_le_bytes(b.try_into().unwrap()))
                .collect()
        };
        let (inputs, outputs) = (read("in"), read("out"));
        let wrong: Vec<_> = inputs
            .iter()
            .zip(&outputs)
            .filter(|(x, y)| f(**x).to_bits() != y.to_bits())
            .take(5)
            .collect();
        assert!(
            wrong.is_empty(),
            "{name}: {} points, first {wrong:?}",
            inputs.len()
        );
    }

    #[test]
    fn exp_is_glibcs() {
        against_numpy("exp", exp);
    }

    #[test]
    fn log_is_glibcs() {
        against_numpy("log", log);
    }

    #[test]
    fn log1p_is_glibcs() {
        against_numpy("log1p", log1p);
    }

    #[test]
    fn expm1_is_glibcs() {
        against_numpy("expm1", expm1);
    }
}

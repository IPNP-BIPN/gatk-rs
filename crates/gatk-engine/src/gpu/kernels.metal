// The device half of Milestone GPU (docs/gpu/eligibility.md).
//
// Apple GPUs have no binary64 unit, and LoglessPairHMM is double precision throughout. A float
// kernel would therefore be a different function, so every double here is a `ulong` holding its
// IEEE 754 bits and every operation is integer arithmetic: correctly rounded to nearest-even,
// subnormals kept, no fused multiply-add. Integer work has one answer on every device, which is
// what makes this kernel eligible without a per-architecture argument. tools/audit/gpu_kernels.py
// refuses a floating-point type anywhere in this file.

#include <metal_stdlib>
using namespace metal;

constant ulong SIGN = 0x8000000000000000ul;
constant ulong EXP_MASK = 0x7FF0000000000000ul;
constant ulong FRAC_MASK = 0x000FFFFFFFFFFFFFul;
constant ulong QUIET_BIT = 0x0008000000000000ul;
constant ulong HIDDEN_BIT = 0x0010000000000000ul;
// The NaN an invalid operation produces on AArch64 with FPCR.DN clear, the host this kernel is
// compared against. x86-64 produces the negative one instead, so the CPU path is itself
// host-dependent here; PairHMM refuses a NaN likelihood, so neither reaches an output.
constant ulong DEFAULT_NAN = 0x7FF8000000000000ul;

static inline bool is_nan(ulong a) { return (a & ~SIGN) > EXP_MASK; }

// AArch64's propagation order: a signalling operand first, then the first quiet one.
static inline ulong propagate_nan(ulong a, ulong b) {
    bool a_signalling = is_nan(a) && (a & QUIET_BIT) == 0;
    bool b_signalling = is_nan(b) && (b & QUIET_BIT) == 0;
    if (a_signalling) return a | QUIET_BIT;
    if (b_signalling) return b | QUIET_BIT;
    return is_nan(a) ? a : b;
}

static inline uint clz64(ulong x) {
    uint hi = uint(x >> 32);
    return hi != 0 ? clz(hi) : 32 + clz(uint(x));
}

// Shift right, OR-ing every bit shifted out into bit zero, so the rounding still sees them.
static inline ulong shift_right_sticky(ulong x, uint shift) {
    if (shift == 0) return x;
    if (shift >= 64) return x != 0 ? 1ul : 0ul;
    return (x >> shift) | ((x << (64 - shift)) != 0 ? 1ul : 0ul);
}

// The value is sig * 2^(exponent - 1085), sig nonzero and below 2^63. Normalised to bit 62, the
// 53 bits from 62 down to 10 are kept and bits 9 to 0 decide the rounding.
static ulong round_and_pack(ulong sign, int exponent, ulong sig) {
    if (sig == 0) return sign;
    int shift = int(clz64(sig)) - 1;
    sig <<= shift;
    exponent -= shift;
    if (exponent >= 0x7FF) return sign | EXP_MASK;
    if (exponent <= 0) {
        sig = shift_right_sticky(sig, uint(1 - exponent));
        exponent = 0;
    }
    ulong kept = sig >> 10;
    ulong rest = sig & 0x3FFul;
    if (rest > 0x200ul || (rest == 0x200ul && (kept & 1ul) != 0)) kept += 1;
    if (exponent == 0) return sign | kept;
    // `kept` carries the hidden bit, so it lands one exponent up; a carry out of the significand
    // moves the exponent again, and reaching 0x7FF this way is the overflow to infinity.
    return sign | ((ulong(exponent - 1) << 52) + kept);
}

static inline void unpack(ulong a, thread int &exponent, thread ulong &sig) {
    int field = int((a & EXP_MASK) >> 52);
    sig = a & FRAC_MASK;
    if (field == 0) {
        // A subnormal stays unnormalised: the addition aligns it as it is.
        exponent = 1;
    } else {
        exponent = field;
        sig |= HIDDEN_BIT;
    }
}

// a + b, as an IEEE 754 binary64 addition rounded to nearest-even.
static ulong f64_add(ulong a, ulong b) {
    if (is_nan(a) || is_nan(b)) return propagate_nan(a, b);
    ulong a_mag = a & ~SIGN;
    ulong b_mag = b & ~SIGN;
    if (a_mag == EXP_MASK || b_mag == EXP_MASK) {
        if (a_mag == EXP_MASK && b_mag == EXP_MASK && (a & SIGN) != (b & SIGN)) return DEFAULT_NAN;
        return a_mag == EXP_MASK ? a : b;
    }
    if (a_mag == 0 && b_mag == 0) return (a & SIGN) == (b & SIGN) ? a : 0ul;
    if (a_mag == 0) return b;
    if (b_mag == 0) return a;
    // Order by magnitude: the result takes the larger operand's sign.
    if (b_mag > a_mag) {
        ulong t = a;
        a = b;
        b = t;
    }
    int a_exp, b_exp;
    ulong a_sig, b_sig;
    unpack(a, a_exp, a_sig);
    unpack(b, b_exp, b_sig);
    a_sig <<= 9;
    b_sig = shift_right_sticky(b_sig << 9, uint(a_exp - b_exp));
    ulong sum = (a & SIGN) == (b & SIGN) ? a_sig + b_sig : a_sig - b_sig;
    // Exact cancellation is +0 under round-to-nearest.
    if (sum == 0) return 0ul;
    return round_and_pack(a & SIGN, a_exp + 1, sum);
}

// a * b, as an IEEE 754 binary64 multiplication rounded to nearest-even.
static ulong f64_mul(ulong a, ulong b) {
    if (is_nan(a) || is_nan(b)) return propagate_nan(a, b);
    ulong sign = (a ^ b) & SIGN;
    ulong a_mag = a & ~SIGN;
    ulong b_mag = b & ~SIGN;
    if (a_mag == EXP_MASK || b_mag == EXP_MASK) {
        if (a_mag == 0 || b_mag == 0) return DEFAULT_NAN;
        return sign | EXP_MASK;
    }
    if (a_mag == 0 || b_mag == 0) return sign;
    int a_exp, b_exp;
    ulong a_sig, b_sig;
    unpack(a, a_exp, a_sig);
    unpack(b, b_exp, b_sig);
    // A subnormal is normalised first, or the product below would keep too few of its bits.
    uint a_shift = clz64(a_sig) - 11, b_shift = clz64(b_sig) - 11;
    a_sig <<= a_shift;
    a_exp -= int(a_shift);
    b_sig <<= b_shift;
    b_exp -= int(b_shift);
    // The 106-bit product from 32-bit limbs.
    ulong a0 = a_sig & 0xFFFFFFFFul, a1 = a_sig >> 32;
    ulong b0 = b_sig & 0xFFFFFFFFul, b1 = b_sig >> 32;
    ulong p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
    ulong middle = (p00 >> 32) + (p01 & 0xFFFFFFFFul) + (p10 & 0xFFFFFFFFul);
    ulong low = (p00 & 0xFFFFFFFFul) | (middle << 32);
    ulong high = p11 + (p01 >> 32) + (p10 >> 32) + (middle >> 32);
    // Both significands are in [2^52, 2^53), so the product is in [2^104, 2^106). Keep its top
    // bits, with everything below folded into the sticky bit.
    ulong sig = (high << 21) | (low >> 43) | ((low & ((1ul << 43) - 1)) != 0 ? 1ul : 0ul);
    // sig * 2^43 * 2^(a_exp + b_exp - 2150) is sig * 2^(exponent - 1085).
    return round_and_pack(sign, a_exp + b_exp - 1022, sig);
}

// One operation per thread over two operand arrays, for the test that holds the arithmetic above
// to the host's own: 0 is a + b, 1 is a * b.
kernel void f64_ops(device const ulong *a [[buffer(0)]],
                    device const ulong *b [[buffer(1)]],
                    device ulong *out [[buffer(2)]],
                    constant uint &operation [[buffer(3)]],
                    constant uint &count [[buffer(4)]],
                    uint id [[thread_position_in_grid]]) {
    if (id >= count) return;
    out[id] = operation == 0 ? f64_add(a[id], b[id]) : f64_mul(a[id], b[id]);
}

struct PairJob {
    uint read;
    uint haplotype;
};

// `LoglessPairHMM.computeReadLikelihoodGivenHaplotypeLog10` up to the sum of the last row, one
// read/haplotype pair per thread. The parallelism is between pairs, never inside one, so each
// pair's operations run in the CPU path's order (crates/gatk-engine/src/pair_hmm.rs). The host
// supplies, per read base, the six transitions and the two priors, and takes the logarithm.
//
// Three rolling rows replace the three matrices: cell (i, j) reads (i-1, j-1), (i-1, j) and
// (i, j-1), so the diagonal is carried in registers and a row is overwritten left to right.
kernel void pair_hmm(device const uchar *read_bases [[buffer(0)]],
                     device const uint2 *read_spans [[buffer(1)]],
                     device const ulong *read_rows [[buffer(2)]],
                     device const uchar *haplotype_bases [[buffer(3)]],
                     device const uint2 *haplotype_spans [[buffer(4)]],
                     device const ulong *haplotype_initial [[buffer(5)]],
                     device const PairJob *jobs [[buffer(6)]],
                     device ulong *scratch [[buffer(7)]],
                     device ulong *out [[buffer(8)]],
                     constant uint &job_count [[buffer(9)]],
                     constant uint &stride [[buffer(10)]],
                     uint id [[thread_position_in_grid]]) {
    if (id >= job_count) return;
    PairJob job = jobs[id];
    uint2 read = read_spans[job.read];
    uint2 haplotype = haplotype_spans[job.haplotype];
    uint width = haplotype.y + 1;
    device ulong *m = scratch + ulong(id) * ulong(stride);
    device ulong *ins = m + width;
    device ulong *del = ins + width;
    ulong initial = haplotype_initial[job.haplotype];
    for (uint j = 0; j < width; j++) {
        m[j] = 0ul;
        ins[j] = 0ul;
        del[j] = initial;
    }
    for (uint i = 0; i < read.y; i++) {
        device const ulong *t = read_rows + ulong(read.x + i) * 8;
        uchar read_base = read_bases[read.x + i];
        ulong diagonal_m = m[0], diagonal_i = ins[0], diagonal_d = del[0];
        m[0] = 0ul;
        ins[0] = 0ul;
        del[0] = 0ul;
        ulong left_m = 0ul, left_d = 0ul;
        for (uint j = 1; j < width; j++) {
            uchar haplotype_base = haplotype_bases[haplotype.x + j - 1];
            bool same = read_base == haplotype_base || read_base == 'N' || haplotype_base == 'N';
            ulong prior = same ? t[6] : t[7];
            ulong up_m = m[j], up_i = ins[j], up_d = del[j];
            ulong next_m = f64_mul(
                prior,
                f64_add(f64_add(f64_mul(diagonal_m, t[0]), f64_mul(diagonal_i, t[1])),
                        f64_mul(diagonal_d, t[1])));
            ulong next_i = f64_add(f64_mul(up_m, t[2]), f64_mul(up_i, t[3]));
            ulong next_d = f64_add(f64_mul(left_m, t[4]), f64_mul(left_d, t[5]));
            m[j] = next_m;
            ins[j] = next_i;
            del[j] = next_d;
            diagonal_m = up_m;
            diagonal_i = up_i;
            diagonal_d = up_d;
            left_m = next_m;
            left_d = next_d;
        }
    }
    ulong total = 0ul;
    for (uint j = 1; j < width; j++) total = f64_add(total, f64_add(m[j], ins[j]));
    out[id] = total;
}

// ---- binary32, for NVScoreVariants' convolutions (crates/gatk-tools/src/nv_score_variants) ----
//
// A float is carried as the `uint` holding its bits, and every operation goes through the binary64
// arithmetic above, so this half is integer work too.

constant ulong F32_MIN_NORMAL_AS_F64 = 0x3810000000000000ul; // 2^-126

// `f64::from(f32)`: exact. A signalling NaN comes back quiet, as `fcvt` returns it.
static ulong f32_to_f64(uint a) {
    ulong sign = ulong(a & 0x80000000u) << 32;
    uint field = (a >> 23) & 0xFFu;
    ulong frac = ulong(a & 0x7FFFFFu);
    if (field == 0xFFu) return sign | EXP_MASK | (frac != 0 ? QUIET_BIT : 0ul) | (frac << 29);
    if (field == 0u) {
        if (frac == 0) return sign;
        // A subnormal float is a normal double: bring its leading bit up to the hidden bit.
        uint shift = clz64(frac) - 40;
        frac = (frac << shift) & 0x7FFFFFul;
        return sign | (ulong(897 - shift) << 52) | (frac << 29);
    }
    return sign | (ulong(field + 896) << 52) | (frac << 29);
}

// `f64 as f32`: rounded to nearest-even, subnormals kept, overflow to infinity, NaN quietened.
static uint f64_to_f32(ulong a) {
    uint sign = uint(a >> 32) & 0x80000000u;
    int field = int((a & EXP_MASK) >> 52);
    ulong frac = a & FRAC_MASK;
    if (field == 0x7FF) return frac != 0 ? sign | 0x7FC00000u | uint(frac >> 29) : sign | 0x7F800000u;
    // A subnormal double is below half the smallest subnormal float.
    if (field == 0) return sign;
    int exponent = field - 896;
    ulong sig = frac | HIDDEN_BIT;
    uint shift = 29;
    if (exponent <= 0) {
        if (30 - exponent >= 64) return sign;
        shift = uint(30 - exponent);
        exponent = 0;
    }
    ulong kept = sig >> shift;
    ulong rest = sig & ((1ul << shift) - 1ul);
    ulong midpoint = 1ul << (shift - 1);
    if (rest > midpoint || (rest == midpoint && (kept & 1ul) != 0)) kept += 1;
    if (exponent == 0) return sign | uint(kept);
    // As in round_and_pack: the hidden bit carries the exponent up by one, and so does a carry.
    ulong bits = (ulong(exponent - 1) << 23) + kept;
    return bits >= 0x7F800000ul ? sign | 0x7F800000u : sign | uint(bits);
}

// a - b as the CPU's `fsub`: a NaN `b` propagates with its own sign.
static inline ulong f64_sub(ulong a, ulong b) { return f64_add(a, is_nan(b) ? b : b ^ SIGN); }

static inline bool f64_is_zero(ulong a) { return (a & ~SIGN) == 0; }
static inline bool f64_positive(ulong a) { return !is_nan(a) && (a & SIGN) == 0 && !f64_is_zero(a); }

// f32 + f32. The double holds the exact sum to 53 bits, and 53 >= 2 * 24 + 2 makes rounding it
// again to a float the same as rounding the exact sum once.
static uint f32_add(uint a, uint b) {
    // NaNs first, in binary32: widening would quieten a signalling one before it could take
    // precedence, as AArch64's order gives it.
    bool a_nan = (a & 0x7FFFFFFFu) > 0x7F800000u, b_nan = (b & 0x7FFFFFFFu) > 0x7F800000u;
    if (a_nan && (a & 0x00400000u) == 0) return a | 0x00400000u;
    if (b_nan && (b & 0x00400000u) == 0) return b | 0x00400000u;
    if (a_nan) return a;
    if (b_nan) return b;
    return f64_to_f32(f64_add(f32_to_f64(a), f32_to_f64(b)));
}

// `kernels::fused`, operation for operation: `vfmadd231ps` on one lane, a * b + c rounded once.
static uint f32_fused(uint a, uint b, uint c) {
    ulong product = f64_mul(f32_to_f64(a), f32_to_f64(b));
    ulong addend = f32_to_f64(c);
    ulong sum = f64_add(product, addend);
    if ((sum & 0x1FFFFFFFul) != 0x10000000ul && !is_nan(sum) && (sum & ~SIGN) >= F32_MIN_NORMAL_AS_F64)
        return f64_to_f32(sum);
    // `rounded_once`: round to odd in the double.
    ulong from_product = f64_sub(sum, addend);
    ulong error = f64_add(f64_sub(product, from_product), f64_sub(addend, f64_sub(sum, from_product)));
    if (f64_is_zero(error) || (sum & 1ul) != 0) return f64_to_f32(sum);
    bool outward = f64_positive(error) == f64_positive(sum);
    return f64_to_f32(outward ? sum + 1ul : sum - 1ul);
}

// One binary32 operation per thread, for the test that holds the arithmetic above to the host's:
// 0 is `fused(a, b, c)`, 1 is a + b.
kernel void f32_ops(device const uint *a [[buffer(0)]],
                    device const uint *b [[buffer(1)]],
                    device const uint *c [[buffer(2)]],
                    device uint *out [[buffer(3)]],
                    constant uint &operation [[buffer(4)]],
                    constant uint &count [[buffer(5)]],
                    uint id [[thread_position_in_grid]]) {
    if (id >= count) return;
    out[id] = operation == 0 ? f32_fused(a[id], b[id], c[id]) : f32_add(a[id], b[id]);
}

// The words of `ConvShape`, in order.
struct ConvShape {
    uint batch, inputs, height, width;
    uint outputs, kernel_height, kernel_width;
    uint block, bias_last, padded, first;
};

// `kernels::conv2d`, one output element per thread: every accumulator is one thread's own chain of
// `fused`, in the CPU path's order, and no sum is split between threads.
kernel void conv2d_fused(device const uint *src [[buffer(0)]],
                         device const uint *weight [[buffer(1)]],
                         device const uint *bias [[buffer(2)]],
                         device uint *dst [[buffer(3)]],
                         device const ConvShape *shape [[buffer(4)]],
                         constant uint &count [[buffer(5)]],
                         uint thread_id [[thread_position_in_grid]]) {
    if (thread_id >= count) return;
    ConvShape s = *shape;
    uint id = s.first + thread_id;
    uint out_height = s.height + 1 - s.kernel_height;
    uint out_width = s.width + 1 - s.kernel_width;
    uint plane = out_height * out_width;
    uint batch_channel = id / plane;
    uint y = (id % plane) / out_width;
    uint x = (id % plane) % out_width;
    uint channel = batch_channel % s.outputs;
    ulong span = ulong(s.height) * s.width;
    uint span_w = s.kernel_height * s.kernel_width;
    device const uint *image = src + ulong(batch_channel / s.outputs) * s.inputs * span;
    device const uint *w = weight + ulong(channel) * s.inputs * span_w;
    uint sum = s.bias_last != 0 ? 0u : bias[channel];
    for (uint first = 0; first < s.inputs; first += s.block) {
        uint last = min(first + s.block, s.inputs);
        for (uint ky = 0; ky < s.kernel_height; ky++) {
            for (uint kx = 0; kx < s.kernel_width; kx++) {
                uint at_k = ky * s.kernel_width + kx;
                ulong at_i = ulong(y + ky) * s.width + x + kx;
                for (uint c = first; c < last; c++) {
                    sum = f32_fused(image[c * span + at_i], w[c * span_w + at_k], sum);
                }
            }
        }
    }
    if (s.padded != 0) sum = f32_add(sum, 0u);
    if (s.bias_last != 0) sum = f32_add(sum, bias[channel]);
    dst[thread_id] = sum;
}

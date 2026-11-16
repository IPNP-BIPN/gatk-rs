//! Milestone GPU: the device path against the CPU path, bit for bit.
//!
//! Without `GATK_RS_GPU` these pass vacuously, which is the CPU run. With it they refuse to pass
//! unless the device actually ran, so a GPU job cannot go green on a host that fell back.

use gatk_engine::gpu;
use gatk_engine::pair_hmm::read_likelihood_given_haplotype_log10;

/// SplitMix64, so the cases are the same on every run without a dependency.
struct Random(u64);
impl Random {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

fn device_ran<T>(answer: Option<T>) -> Option<T> {
    if gpu::requested() {
        assert!(
            answer.is_some(),
            "GATK_RS_GPU is set but no device path ran"
        );
    }
    answer
}

/// Operands that reach every branch: zeros, subnormals, the normal extremes, infinities, NaNs,
/// values one ulp apart, and pairs whose sum or product rounds at a tie.
fn operands() -> (Vec<u64>, Vec<u64>) {
    let special = [
        0u64,
        1u64 << 63,
        1,
        0x000F_FFFF_FFFF_FFFF,
        0x0010_0000_0000_0000,
        0x7FEF_FFFF_FFFF_FFFF,
        0x7FF0_0000_0000_0000,
        0xFFF0_0000_0000_0000,
        0x7FF8_0000_0000_0000,
        0x7FF0_0000_0000_0001,
        1.0f64.to_bits(),
        1.5f64.to_bits(),
        f64::EPSILON.to_bits(),
        (1.0f64 + f64::EPSILON).to_bits(),
        0x7FB0_0000_0000_0000, // 2^1020, the initial condition
        0x0030_0000_0000_0000, // 2^-1020
        0.1f64.to_bits(),
        (1.0f64 / 3.0).to_bits(),
    ];
    let (mut a, mut b) = (Vec::new(), Vec::new());
    for &x in &special {
        for &y in &special {
            for sign in [0, 1u64 << 63] {
                a.push(x);
                b.push(y ^ sign);
            }
        }
    }
    let mut random = Random(0x6A09_E667_F3BC_C908);
    for _ in 0..400_000 {
        let x = random.next();
        // Half the second operands sit near the first, where cancellation and ties live.
        let y = match random.below(4) {
            0 => random.next(),
            1 => x ^ random.below(1 << 12),
            2 => (x & !(0x7FFu64 << 52)) | ((random.below(0x7FF)) << 52),
            _ => (x ^ (1u64 << 63))
                .wrapping_add(random.below(3))
                .wrapping_sub(1),
        };
        a.push(x);
        b.push(y);
    }
    (a, b)
}

fn same(expected: f64, actual: u64) -> bool {
    // An invalid operation's NaN is host-dependent on the CPU itself (positive on AArch64,
    // negative on x86-64), and PairHMM refuses a NaN, so only "is a NaN" is held there.
    if expected.is_nan() {
        return f64::from_bits(actual).is_nan();
    }
    expected.to_bits() == actual
}

#[test]
fn the_device_adds_and_multiplies_as_the_host_does() {
    let (a, b) = operands();
    for (operation, name) in [(0u32, "+"), (1, "*")] {
        let Some(device) = device_ran(gpu::f64_ops(operation, &a, &b)) else {
            return;
        };
        let mut wrong = 0;
        for i in 0..a.len() {
            let (x, y) = (f64::from_bits(a[i]), f64::from_bits(b[i]));
            let host = if operation == 0 { x + y } else { x * y };
            if !same(host, device[i]) {
                if wrong < 10 {
                    eprintln!(
                        "{:016x} {name} {:016x}: host {:016x}, device {:016x}",
                        a[i],
                        b[i],
                        host.to_bits(),
                        device[i]
                    );
                }
                wrong += 1;
            }
        }
        assert_eq!(wrong, 0, "{wrong} of {} results differ for {name}", a.len());
    }
}

struct Read {
    bases: Vec<u8>,
    qualities: Vec<u8>,
    insertion: Vec<u8>,
    deletion: Vec<u8>,
    gap: Vec<u8>,
}

fn bases(random: &mut Random, length: usize) -> Vec<u8> {
    (0..length)
        .map(|_| {
            let alphabet = if random.below(20) == 0 { 5 } else { 4 };
            b"ACGTN"[random.below(alphabet) as usize]
        })
        .collect()
}

#[test]
fn the_device_pair_hmm_is_the_cpu_pair_hmm() {
    let mut random = Random(0xBB67_AE85_84CA_A73B);
    let reference = bases(&mut random, 600);
    // Haplotypes as an assembler makes them: windows of one reference, some edited.
    let haplotypes: Vec<Vec<u8>> = (0..12)
        .map(|_| {
            let start = random.below(200) as usize;
            let mut h = reference[start..start + 40 + random.below(300) as usize].to_vec();
            for _ in 0..random.below(4) {
                let at = random.below(h.len() as u64) as usize;
                match random.below(3) {
                    0 => h[at] = b"ACGT"[random.below(4) as usize],
                    1 => h.insert(at, b"ACGT"[random.below(4) as usize]),
                    _ => {
                        h.remove(at);
                    }
                }
            }
            h
        })
        .collect();
    let mut reads: Vec<Read> = (0..60)
        .map(|_| {
            let length = 1 + random.below(150) as usize;
            let source = &haplotypes[random.below(haplotypes.len() as u64) as usize];
            let mut b = if length < source.len() {
                let start = random.below((source.len() - length) as u64) as usize;
                source[start..start + length].to_vec()
            } else {
                bases(&mut random, length)
            };
            for _ in 0..random.below(5) {
                let at = random.below(b.len() as u64) as usize;
                b[at] = b"ACGTN"[random.below(5) as usize];
            }
            let quality = |r: &mut Random, low: u64, high: u64| -> Vec<u8> {
                (0..b.len())
                    .map(|_| (low + r.below(high - low)) as u8)
                    .collect()
            };
            Read {
                qualities: quality(&mut random, 2, 60),
                insertion: quality(&mut random, 10, 70),
                deletion: quality(&mut random, 10, 70),
                gap: quality(&mut random, 5, 15),
                bases: b,
            }
        })
        .collect();
    // An empty read, whose likelihood is the log of zero.
    reads.push(Read {
        bases: vec![],
        qualities: vec![],
        insertion: vec![],
        deletion: vec![],
        gap: vec![],
    });
    let inputs: Vec<gpu::ReadInputs> = reads
        .iter()
        .map(|r| gpu::ReadInputs {
            bases: &r.bases,
            qualities: &r.qualities,
            insertion: &r.insertion,
            deletion: &r.deletion,
            gap_continuation: &r.gap,
        })
        .collect();
    let Some(device) = device_ran(gpu::pair_hmm_log10(&haplotypes, &inputs)) else {
        return;
    };
    assert_eq!(device.len(), reads.len() * haplotypes.len());
    for (r, read) in reads.iter().enumerate() {
        for (h, haplotype) in haplotypes.iter().enumerate() {
            let cpu = read_likelihood_given_haplotype_log10(
                haplotype,
                &read.bases,
                &read.qualities,
                &read.insertion,
                &read.deletion,
                &read.gap,
            );
            let gpu = device[r * haplotypes.len() + h];
            assert_eq!(
                cpu.to_bits(),
                gpu.to_bits(),
                "read {r} haplotype {h}: cpu {cpu:e}, gpu {gpu:e}"
            );
        }
    }
}

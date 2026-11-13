//! Milestone GPU: NVScoreVariants' convolutions on the device against `kernels::conv2d_on_host`,
//! bit for bit, and the device's `fused` against the host's.
//!
//! Without `GATK_RS_GPU` these pass vacuously, which is the CPU run. With it they refuse to pass
//! unless the device actually ran. No weights are needed: the host path is already held to the
//! reference by `nv_score_variants_layers`, so random tensors at the models' shapes say enough.

use gatk_engine::gpu;
use gatk_tools::nv_score_variants::kernels::{
    conv2d_on_device, conv2d_on_host, fused, ConvWeights, Shape,
};

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
    /// A float of the size activations and weights have, sometimes zero, either sign.
    fn value(&mut self) -> f32 {
        match self.below(8) {
            0 => 0.0,
            1 => -0.0,
            _ => {
                let magnitude = (self.below(1 << 24) as f32 / (1 << 24) as f32) * 4.0;
                if self.below(2) == 0 {
                    magnitude
                } else {
                    -magnitude
                }
            }
        }
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

/// Triples that reach every branch of `fused`: zeros, subnormals, the extremes, infinities, NaNs,
/// products that cancel the addend, and sums that land on a midpoint between two floats.
fn operands() -> (Vec<u32>, Vec<u32>, Vec<u32>) {
    let special = [
        0u32,
        1 << 31,
        1,
        0x007F_FFFF,
        0x0080_0000,
        0x7F7F_FFFF,
        0x7F80_0000,
        0x7FC0_0000,
        0x7F80_0001,
        1.0f32.to_bits(),
        1.5f32.to_bits(),
        f32::EPSILON.to_bits(),
        (1.0f32 + f32::EPSILON).to_bits(),
        0.1f32.to_bits(),
        (1.0f32 / 3.0).to_bits(),
        0x2000_0000, // 2^-63, whose square is subnormal
    ];
    let (mut a, mut b, mut c) = (Vec::new(), Vec::new(), Vec::new());
    for &x in &special {
        for &y in &special {
            for &z in &special {
                for sign in [0, 1u32 << 31] {
                    a.push(x);
                    b.push(y);
                    c.push(z ^ sign);
                }
            }
        }
    }
    let mut random = Random(0x6A09_E667_F3BC_C908);
    for _ in 0..400_000 {
        let x = random.next() as u32;
        let y = random.next() as u32;
        let product = f32::from_bits(x) * f32::from_bits(y);
        let z = match random.below(4) {
            0 => random.next() as u32,
            // Near minus the product, where the sum cancels and the single rounding shows.
            1 => (-product)
                .to_bits()
                .wrapping_add(random.below(5) as u32)
                .wrapping_sub(2),
            // Near the product's own size, where ties are likeliest.
            2 => product.to_bits() ^ (random.below(1 << 8) as u32),
            _ => random.value().to_bits(),
        };
        // Half the triples are of the size the convolutions see, the rest are any bits at all.
        if random.below(2) == 0 {
            a.push(random.value().to_bits());
            b.push(random.value().to_bits());
            c.push(random.value().to_bits());
        } else {
            a.push(x);
            b.push(y);
            c.push(z);
        }
    }
    (a, b, c)
}

#[test]
fn the_device_fuses_and_adds_as_the_host_does() {
    let (a, b, c) = operands();
    for (operation, name) in [(0, "fused"), (1, "add")] {
        let Some(device) = device_ran(gpu::f32_ops(operation, &a, &b, &c)) else {
            return;
        };
        let mut wrong = 0;
        for i in 0..a.len() {
            let (x, y, z) = (
                f32::from_bits(a[i]),
                f32::from_bits(b[i]),
                f32::from_bits(c[i]),
            );
            let host = if operation == 0 {
                fused(x, y, z)
            } else {
                x + y
            }
            .to_bits();
            if host != device[i] {
                if wrong < 10 {
                    eprintln!(
                        "{name}({:08x}, {:08x}, {:08x}): host {host:08x}, device {:08x}",
                        a[i], b[i], c[i], device[i]
                    );
                }
                wrong += 1;
            }
        }
        assert_eq!(wrong, 0, "{name}: {wrong} of {} differ", a.len());
    }
}

/// One convolution at a shape a model has, with random weights and input.
fn convolve(input: Shape, outputs: usize, kernel: (usize, usize), bias_last: bool, seed: u64) {
    let mut random = Random(seed);
    let weight: Vec<f32> = (0..outputs * input.channels * kernel.0 * kernel.1)
        .map(|_| random.value() / 16.0)
        .collect();
    let bias: Vec<f32> = (0..outputs).map(|_| random.value()).collect();
    let src: Vec<f32> = (0..input.len()).map(|_| random.value()).collect();
    let conv = ConvWeights {
        outputs,
        inputs: input.channels,
        height: kernel.0,
        width: kernel.1,
        weight: &weight,
        bias: &bias,
        bias_last,
    };
    let Some(device) = device_ran(conv2d_on_device(&src, input, &conv)) else {
        return;
    };
    let (host, _) = conv2d_on_host(&src, input, &conv, 4);
    assert_eq!(host.len(), device.len());
    let wrong = (0..host.len())
        .filter(|&i| host[i].to_bits() != device[i].to_bits())
        .count();
    assert_eq!(
        wrong,
        0,
        "{input:?} -> {outputs}, kernel {kernel:?}: {wrong} of {} differ",
        host.len()
    );
}

#[test]
fn the_device_convolves_as_the_host_does() {
    let shape = |batch, channels, height, width| Shape {
        batch,
        channels,
        height,
        width,
    };
    // The 1D model: `jit:avx2`, a flat first layer of 4 channels, then blocks of 8.
    convolve(shape(3, 4, 1, 128), 32, (1, 12), false, 1);
    convolve(shape(2, 256, 1, 117), 16, (1, 12), false, 2);
    // A partial last block, which gets the `+ 0.0` of its zero padding.
    convolve(shape(2, 20, 1, 40), 8, (1, 12), false, 3);
    // The 2D model: `brgconv:avx2`, from zero over all the channels, the bias last.
    convolve(shape(2, 15, 32, 32), 16, (3, 3), true, 4);
    convolve(shape(1, 64, 20, 24), 8, (2, 5), true, 5);
}

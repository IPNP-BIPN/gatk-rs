# The targets, each with its bytes and its speed

Milestone GPU, GPU.2 (#77) and GPU.3 (#78). A kernel is listed with the test that holds it to the
CPU path's bytes and the measurement that says whether it pays. One that does not pay is deleted
rather than kept, and a target that cannot pay is recorded here instead of built.

## 1. `PairHMM`: shipped

**What runs on the device.** `LoglessPairHMM`'s lattice and its last-row sum, one read/haplotype
pair per thread, in `crates/gatk-engine/src/gpu/kernels.metal`. The transitions, the priors and the
initial condition are computed on the host by the CPU path's own functions and sent down as bits;
the sum comes back as bits and the host takes `jmath`'s `log10`. No transcendental runs on the
device.

**Why it cannot move a byte.** Each cell is the CPU expression, operation for operation, in
software binary64 on integers (docs/gpu/eligibility.md); the parallelism is between pairs, never
inside one, so no sum is regrouped. `PairHMMLikelihoodCalculationEngine.computeReadLikelihoods`
takes the device's values when `GATK_RS_GPU` is set and its own otherwise, and raises a DRAGstr
imputation refusal at the same read, after the same results-file lines, either way.

**The tests that prove it.** `crates/gatk-engine/tests/gpu.rs`: 800,648 additions and
multiplications against the host bit for bit, and 732 read/haplotype pairs (reads of 1 to 150
bases with `N`s and an empty read, edited haplotypes of 40 to 340) against
`read_likelihood_given_haplotype_log10`. Flipping the tie-breaking of the rounding fails both,
on 50,357 additions and on the 40th read. Then every suite of the workspace runs again on the
device path (`GATK_RS_GPU=1 cargo test --workspace --release --features gatk-engine/gpu`), the
HaplotypeCaller and Mutect2 end-to-end goldens among them.

**Speed.** `crates/gatk-engine/examples/gpu_speed.rs` on an Apple M4 Max (16 cores, 40 GPU
cores), second of two runs, with the device's compilation charged to a warm-up. The engine runs the
CPU path on one thread, so that is the column the gate reads; all cores are shown for honesty. The
host had other work on it (load average about 24), so the CPU columns are pessimistic and the
ratios are a floor rather than a figure to quote to two digits.

| Workload | pairs | cells | CPU, 1 thread | CPU, 16 threads | device | vs 1 thread | vs 16 threads | same bits |
|---|---:|---:|---:|---:|---:|---:|---:|:---:|
| 4 haplotypes of 300, 100 reads of 150 | 400 | 18 M | 0.270 s | 0.062 s | 0.158 s | 1.7x | 0.39x | yes |
| 16 of 400, 500 reads | 8,000 | 480 M | 6.30 s | 1.07 s | 0.297 s | 21x | 3.6x | yes |
| 32 of 500, 2,000 reads | 64,000 | 4.8 G | 58.9 s | 9.65 s | 2.63 s | 22x | 3.7x | yes |

The kernel stays: it is faster than the path it replaces at every size measured, and from a few
thousand pairs it is faster than the CPU path spread over every core. At the smallest region the
fixed cost of a dispatch is most of its time, which is the case a size threshold would route back
to the CPU if a profile ever shows small regions dominating.

To run it: build with `--features gatk-engine/gpu` on macOS and set `GATK_RS_GPU=1`.

## 1b. NVScoreVariants' convolutions: shipped (GPU.4, #1021)

**Whether a device could pay at all.** #1021 asked first how much of the tool is inference. With
`crates/gatk-tools/examples/nvscore_profile.rs` on the `make_fixture.py` corpus (71 variants), the
encoding, FASTA, VCF and BAM reading included, takes 0.001 s in 1D and 0.03 s in 2D; inference is
86% of the whole tool in both modes, and the rest is the tool's own start-up around it. The bar was
10%.

**What runs on the device.** Every convolution of both models, one output element per thread, in
`conv2d_fused` (`crates/gatk-engine/src/gpu/kernels.metal`). The batch norm, the pooling, the dense
layers and the softmax stay on the host. Each accumulator is `kernels::conv2d_on_host`'s chain of
`fused`, oneDNN's order, in the same blocks; `fused` is transliterated operation for operation onto
the kernel's software binary64, and floats travel as their bits, so the kernel is of class
`integer` like the PairHMM. #1021 expected that "a GPU path cannot be the byte-identical path";
that holds for a GPU's own float units, not for this.

**The tests that prove it.** `crates/gatk-tools/tests/nv_score_variants_gpu.rs`: 408,192 `fused`
and as many binary32 additions against the host bit for bit, specials (zeros, subnormals,
infinities, both kinds of NaN) and products that cancel their addend among them; then five
convolutions at the two models' shapes, flat first layer, blocks of 8, a partial block and
`brgconv`'s bias-last, against `conv2d_on_host`. Rounding ties to odd, dropping the round-to-odd
step and summing in one block instead of 8 each fail it. Then, locally with the reference's dumps
made under Docker as the `nvscorevariants` job makes them, the device path against PyTorch itself:
`nv_score_variants_layers` finds every convolution of both models bit-identical to the reference's
layer dump (0 of 60,489,728 values at the 2D model's `conv.0`, 0 of 2,126,592 at the 1D model's,
and 0 at every other layer), and `nv_score_variants_output` finds 31 cases, 337 score rows, 0
differing parts. Those two need the weights and run in CI only on the CPU path, in the
`nvscorevariants` job; the `gpu` job runs the weight-free test above.

**Speed.** The same Apple M4 Max as above, best of two runs, each mode on the CPU path then the
device path. The host was loaded (load average 18 to 33 on 16 cores), so the CPU columns are
pessimistic and the ratios a floor.

| Mode | variants | tool, CPU, 16 threads | tool, device | inference, CPU, 1 thread | inference, device | vs 16 threads | vs 1 thread | same bytes |
|---|---:|---:|---:|---:|---:|---:|---:|:---:|
| 1D | 71 | 3.07 s | 0.84 s | 27.2 s | 0.97 s | 3.7x | 28x | yes |
| 2D | 71 | 38.8 s | 9.9 s | 385 s | 9.7 s | 3.9x | 40x | yes |

"Same bytes" is the output VCF and the scores file of each mode, compared with `cmp` between the two
paths. The kernel stays: it beats the CPU path spread over every core in both modes.

## 2. The read filters and the pileup floor: not built, because they cannot pay

A filter is a handful of loads and compares on fields of a decoded `BamRecord`. To run it on the
device those same fields must first be gathered on the host into a packed buffer, which reads
exactly what the CPU predicate reads; the device then adds a copy, a dispatch and a read-back. The
host's share alone is at least the CPU path's whole cost, so no kernel can beat it, on any device,
however fast. The same holds for the pileup's base and quality floors, which are one compare per
element of an array the host has just built. GPU.3 deletes an accelerated path that is not faster;
this one is recorded instead of written. It becomes a target only if decoding itself moves to the
device, which is the BGZF entry below.

## 3. BGZF inflate and deflate: open, in htsjdk-rs

`htsjdk-bgzf` and the GKL-exact deflate (`gkl-deflate`) live in IPNP-BIPN/htsjdk-rs, pinned here by
revision. Inflate is decidable (any correct decoder gives the bytes) and block-independent; deflate
must reproduce GKL's choice per level (zlib's for most, ISA-L igzip's at level 2, which is
GATKConfig's), so the device kernel would be a transliteration of each matcher, one block per
thread. Neither is built yet.

## 4. The assembly graph and genotyping: open

Pointer-chasing over a de Bruijn graph and floating point throughout (G3). Every reduction would
need the CPU path's tree; not attempted.

## 5. CRAM codecs: open, in htsjdk-rs

rANS and the other codecs are integral with defined output, and live in `htsjdk-cram`. The only
CRAM a ported tool touches today is `SplitCRAM`, which copies containers without decoding them, so
no tool here would take a device decoder yet.

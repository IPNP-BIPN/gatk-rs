# The targets, each with its bytes and its speed

Milestone GPU, GPU.2 (#77) and GPU.3 (#78). A kernel is listed with the test that holds it to the
CPU path's bytes and the measurement that says whether it pays. One that does not pay is deleted
rather than kept, and a target that cannot pay is recorded here instead of built.

## 1. `PairHMM`: built, measured, deleted because it cannot pay

**What ran on the device.** `LoglessPairHMM`'s lattice and its last-row sum, one read/haplotype pair
per thread, in software binary64 on `ulong` (class `integer`). It was bit-identical to the CPU path
on every test and on every HaplotypeCaller and Mutect2 golden run on the device path. Removed on
2026-11-26. Before that it lived in `kernels.metal` (`pair_hmm`), `gpu::pair_hmm_log10`, the
`PairHMMLikelihoodCalculationEngine` call site, `tests/gpu.rs` and `examples/gpu_speed.rs`; see git
history.

**Why it was deleted: PairHMM is 3.3% of HaplotypeCaller.** HaplotypeCaller on chr20:10–11 Mb of
NA12878 30x (Apple M4 Max, 16 cores, load average 57–90 from other work). The port is built with
rustc 1.98.0, `--features gatk-engine/gpu`. The VCF is the same on both paths (records md5
`f0281f9b…`, 1,629 records), and no device call fell back.

| path | CPU user+sys | peak RSS | wall | PairHMM time over all calls |
|---|---:|---:|---:|---:|
| CPU | 56.7 / 59.6 s | 19.1 / 18.0 GB | 64.9 / 78.0 s | 1.89 s (0.68 ns/cell) |
| device (`GATK_RS_GPU=1`) | 59.0 / 57.7 s | 18.8 / 17.3 GB | 108.7 / 128.9 s | 47.4 s (17.1 ns/cell) |

- **The ceiling.** The tool made 1,237 PairHMM calls totalling 2.77 G cells. The CPU path spent 1.89 s
  on them, 3.3% of the tool's CPU. A free device could take at most that off, which is under the
  10% bar §1b sets.
- **The device loses besides.** A call costs about 38 ms on the device against 1.5 ms on the CPU.
  Even the 839 calls under 1 M cells cost 11.5 ms each, so the dispatch's fixed cost dominates.
  Regions are small: the mean call is 2.2 M cells, none reaches 10^8, and only 66 reach 10^7.
- **A size threshold cannot rescue it.** In `gpu_speed` the device was still 0.19× one CPU thread at
  18 M cells, and 1.06–1.17× sixteen threads at 480 M to 4.8 G. Real regions sit far below where it
  pays, so a threshold would route almost every call back to the CPU.
- **GitHub's paravirtualised macOS GPU** (run 45304648109): 0.12–0.47× the CPU path, and the
  4.8 G-cell workload hung the GPU (`kIOGPUCommandBufferCallbackErrorHang`).
- The CPU path got faster since this kernel was first measured: 16 threads on the 480 M-cell workload
  went from 1.07 s to 0.31 s, while the device stayed at 0.29 s.

**Scope of the verdict.** One 30x germline slice on a loaded host. Mutect2 and deep or high-ploidy
data, where regions are larger, were not measured. A kernel for them would have to clear the same
ceiling first.

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

## 3. BGZF inflate: not built, because it cannot pay here; deflate: open, in htsjdk-rs

`htsjdk-bgzf` and the GKL-exact deflate (`gkl-deflate`) live in IPNP-BIPN/htsjdk-rs, pinned here by
revision.

**Inflate** is decidable (any correct decoder gives the bytes) and block-independent, so a kernel
would be easy to hold to the CPU path. The question is what it could save, and that can be bounded
without writing one: run a read-bound tool on a BAM, then on the same BAM with stored blocks,
whose inflate is a copy. The difference is the most any inflater, on any device, at any speed and
with no transfer cost, could take off the tool. `tools/speed/bgzf_inflate_bound.py` does that.
With CountReads, the most read-bound tool there is, on 531,540 reads (`snvq_metrics_test.bam`'s
records 20 times, 894 MB inflated, 193 MB at level 5), on the M4 Max with other work on it (load
average 17), best of three:

| BAM | 1 thread, wall | 1 thread, CPU | 16 threads, wall | 16 threads, CPU |
|---|---:|---:|---:|---:|
| level 5 | 5.54 s | 4.91 s | 1.77 s | 4.07 s |
| stored | 1.73 s | 1.26 s | 1.69 s | 1.34 s |

On one thread, inflate is 69% of the tool. But the port already inflates on every core
(`decompress_all` on rayon), and at 16 threads it is 0.08 s of 1.77 s, 4.5%: the record decoding
that follows is what the wall clock waits for. A device inflater would save at most that 4.5%, under
the 10% #1021 set as the bar for a GPU path. So it is recorded instead of built, like the read
filters. It becomes worth measuring on a host with few cores and a large GPU, or if record decoding
also moves to the device.

One CPU-side finding on the way: inside the tool, inflate costs 2.6 s of CPU where one pass of
`decompress_all` over the same bytes costs 0.87 s. The inflater is not the cause: a read walker
inflates the whole BAM three times, once to sniff it at startup (`read_walker_startup_with` in
`crates/gatk-cli/src/runners.rs`), once for the header (`ReadsDataSource::open_unindexed` in
`crates/gatk-engine/src/reads.rs`) and once for the records (`iter_all`), and 3 x 0.87 s is the
gap. Inflating once and handing the buffer on belongs to Milestone S; it also shrinks the bound
above, since a single pass is a third of what stored blocks saved here.

**Deflate** must reproduce GKL's choice per level (zlib's for most, ISA-L igzip's at level 2, which
is GATKConfig's), so a kernel would transliterate each matcher, one block per thread. The pure-Rust
igzip port (htsjdk-rs#281, pinned here since 2026-11-14) made level 2 reachable on ARM, so the bound
can now be taken: `tools/speed/bgzf_inflate_bound.py --deflate` has PrintReads read the stored BAM
above and write it at level 2, then at level 0, whose deflate is a copy, and checks both hold the
same records. On the M4 Max, best of three, with the host at a load average of 145 to 165 on 16
cores:

| level | output | 1 thread, wall | 1 thread, CPU | 16 threads, wall | 16 threads, CPU |
|---:|---:|---:|---:|---:|---:|
| 2 | 220 MB | 7.86 s | 4.63 s | 29.5 s | 5.91 s |
| 0 | 894 MB | 3.73 s | 2.38 s | 10.3 s | 2.86 s |

Deflate is about half of PrintReads' CPU (2.25 s of 4.63 s on one thread, 3.05 s of 5.91 s on 16),
far over the 10% bar, so unlike inflate it stays a target. What it is worth in wall time on every
core cannot be read from this run: at that load the 16-thread wall is slower than the 1-thread one.
That figure, on an idle host, decides whether the kernel is built; the kernel itself would be
igzip's level-2 matcher and zlib's per-level matchers transliterated onto integers, which is open.

## 4. The assembly graph and genotyping: open

Pointer-chasing over a de Bruijn graph and floating point throughout (G3). Every reduction would
need the CPU path's tree; not attempted.

## 5. CRAM codecs: open, in htsjdk-rs

rANS and the other codecs are integral with defined output, and live in `htsjdk-cram`. The only
CRAM a ported tool touches today is `SplitCRAM`, which copies containers without decoding them, so
no tool here would take a device decoder yet.

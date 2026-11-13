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
imputation refusal at the same read, after the same results-file lines, either way. A Metal command
that ends in error (in run 43382302494, on GitHub's paravirtualised GPU, one sample came back with no
reads, which is what an output left at zero gives; the host did not look then) makes the call answer `None`, so the CPU path computes it, and says so on stderr;
the `gpu` job counts those calls.

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

One CPU-side finding on the way: inside the tool, inflate costs 3.6 s of CPU (4.91 - 1.26 s), where
`htsjdk-bgzf`'s `inflate_speed` example inflates the same bytes in about 0.95 s. That gap is the
reader's, not the decoder's, and belongs to Milestone S.

**Deflate** must reproduce GKL's choice per level (zlib's for most, ISA-L igzip's at level 2, which
is GATKConfig's), so a kernel would transliterate each matcher, one block per thread. Deflate is
several times inflate's cost, so it is the side that could pay. It cannot be measured on this host
yet: level 2 needs ISA-L igzip, which has no ARM build here (`PrintReads` refuses with "GKL routes
level 2 through ISA-L igzip, which this build cannot reach"). The pure-Rust igzip port in
htsjdk-rs is what would make it measurable.

## 4. The assembly graph and genotyping: open

Pointer-chasing over a de Bruijn graph and floating point throughout (G3). Every reduction would
need the CPU path's tree; not attempted.

## 5. CRAM codecs: open, in htsjdk-rs

rANS and the other codecs are integral with defined output, and live in `htsjdk-cram`. The only
CRAM a ported tool touches today is `SplitCRAM`, which copies containers without decoding them, so
no tool here would take a device decoder yet.

# NVScoreVariants: what the port reproduces, in numbers

`NVScoreVariants` runs in pure Rust in `crates/gatk-tools/src/nv_score_variants/`, with no PyTorch, no
MKL and no host libm. Every number below was measured against `broadinstitute/gatk:4.6.2.0` **under
Rosetta on an arm64 Mac, where the image dispatches AVX2**. None comes from real x86-64 yet.

`CNNScoreVariants` does not exist in 4.6.2.0: `DeprecatedToolsRegistry` retired it as of 4.6.1.0, so
this is the only ML tool. Stack of the image: torch 2.1.0 (GCC 12.3), MKL 2022.2, oneDNN 3.1.1.

The fixture is `tools/nvscorevariants/make_fixture.py`: 71 records, windows clipped at both ends of a
contig, IUPAC reference, a pile deeper than the reservoir, and 15 more VCFs for the corners.

## The reference is not bit-stable with itself

Its VCF and its score file are identical across two runs and at 1, 2, 4 and 16 threads, and its score
file at batch sizes 1, 7 and 64. Its bits are not, on all 71 rows of the 1D model, against one thread
and batch 64:

| run | conv, normalize, dense_1 | dense_2.0 | dense_3 | softmax | score |
|---|---|---|---|---|---|
| 4 threads, batch 64 | identical | 2047/2840, 12288 ulp | 135/284, 64 ulp | 111/284, 38 ulp | 20/71, 5 ulp |
| 16 threads, batch 64 | identical | 2077/2840, 10240 ulp | 148/284, 32 ulp | 114/284, 56 ulp | 19/71, 7 ulp |
| 1 thread, batch 1 | conv.0 1429647/2126592, 2949120 ulp | | | | 43/71, 16 ulp |
| 16 threads, batch 7 | | | | | 28/71, 9 ulp |

So the port's bit target is one thread at batch 64, and the printed `%.3f` is the only thing the tool
itself keeps constant. That is what the suite compares.

## The port against one thread, batch 64

Values differing, of all values, on all 71 rows, with the largest distance:

| layer | 1D | 2D |
|---|---|---|
| encoded input | 0/36352 | 0/17448960 |
| annotations | 0/497 | 0/497 |
| convolutions | 0 of 4916608 (3 layers) | 0 of 119798016 (4 layers) |
| batch norm | 0/497 | 0/497 |
| first dense | 19/2840, 8 ulp | 23/4544, 20 ulp |
| second dense | 29/2840, 256 ulp | 13/1704, 64 ulp |
| logits | 3/284, 3 ulp | 2/284, 2 ulp |
| softmax | 4/284, 13 ulp | 3/284, 4 ulp |
| score | 2/71, 2 ulp, at most 9.5e-7 | 0/71 |
| **score at `%.3f`** | **0/71** | **0/71** |

Split by batch, every difference above but one is in the file's short second batch of 7 rows. **The
first batch of 64 rows is bit-identical through every layer of both models.** The one exception is its
1D score: 1/64, 1 ulp, the closed-source `log` below.

The port's own `gatk NVScoreVariants` writes VCFs with the same MD5 as the reference's measured in
the stability run: 1D `fbc9a12a55a1dca37995555de16dcd24`, 2D `26becf88371be2f2dcd16fb81efe57c1`.

The orders behind the zeros:

* **1D convolutions** are `jit:avx2` and **2D convolutions** are `brgconv:avx2`, per `ONEDNN_VERBOSE`.
  The two sum differently: the bias first over channel blocks of 8, or the bias last over all channels.
* **Dense layers** are MKL `SGEMM`, whose source is closed. `tools/nvscorevariants/probe_gemm.py` reads
  the order from outside: blocks of 256 terms with the last block and a bit split in two, and either
  one accumulator or an even and an odd one. Which one depends on the position, mapped for outputs 4,
  24, 40, 64 and batches 1 to 69 and 127 to 129.
* **Softmax** is SLEEF's `expf8_u10` under torch's own loop. **Batch norm** carries GCC's contraction.

## What is not reached

* **The tail column of a short batch.** Batch of 7, row 6: the probe's tree matches the even and odd
  kernel, and the bits do not on 19 of its 40 outputs. All 19 dense_1.0 differences of the 1D model
  are there, and the 2D model's dense_1.1 has the same shape. That kernel is not decoded.
* **A batch of one** takes torch's native convolution and a different `SGEMM` path. Neither is ported:
  the port uses its batch-64 orders there, and a file with 1 record over a multiple of the batch size
  lands on it. Against the reference at batch 1 the port's 1D score is 33/71 different, 18 ulp, at most
  2.0e-6, and 0/71 at `%.3f`.
* **`torch.log`** is MKL `vsLn`, closed source. The port's is correctly rounded. It differs by 1 ulp on
  1 of 71 rows of the softmax measured from the reference's own bits.
* Thread counts above one are not modelled, since the reference itself moves there.

## Whole-tool comparison

`tools/nvscorevariants/NVScoreVariantsDump.py` against `crates/gatk-tools/tests/nv_score_variants_output.rs`:
16 VCFs by 2 models, **32 runs, 337 score rows, every line of the score files and VCFs equal, and the
same `Tool returned` value**. Those include what `create_output_vcf.py` does on the way out: Tool
returned 1 with a truncated VCF on a symbolic allele without `END`, on an undefined INFO, FILTER, FORMAT
or contig, on a breakend and on an empty file, and six-digit floats through htslib everywhere else.

## What #73 can claim, and what it cannot

It can claim, under Rosetta at AVX2: encoding and all convolutions bit-identical at one thread, the
printed score and the VCF byte-identical on all fixture rows of both models at the default settings.

It cannot claim: bit-identical dense layers in a short batch, a bit-identical score, or anything about
AVX-512. It is a **bio-identical at printed precision** result, labelled as such, on this fixture.

## What needs real x86-64

The `nvscorevariants` job in CI, which pulls a second image next to the pinned oracle by digest:

* a host with AVX-512 gets `jit:avx512` or another `brgconv`, other MKL kernels and torch's AVX-512
  `Vec`. So all of the orders above, the class map, block size and `vsLn` are claims about AVX2 only;
* Rosetta does not implement AVX-512 and reports no `xsave`, so the AVX2 path itself should be
  confirmed on a real AVX2 host, which the job's `model name` line identifies;
* the suite's golden may only come from there. It is `golden-pending`, with `golden: null`.

The job prints the layer table and asserts only the inputs and the printed scores and VCFs.

## On real x86-64, GATK 4.7.0.0

The `nvscorevariants` job of run 44799995145 landed on an Intel Xeon 6975P-C. oneDNN reported AVX-512
with AMX and ran the 1D convolutions as `jit:avx512_core`. Values differing, port against the
reference at one thread and batch 64:

| layer | 1D | 2D |
|---|---|---|
| encoded input, annotations | 0 | 0 |
| convolutions | 1217823/2126592 (conv.0), every layer | **0 of 119798016**, all four |
| batch norm | 0/497 | 0/497 |
| first dense | 293/2840 | first batch 0/4096, short batch 30/448 |
| second dense on, softmax, score | most values, up to 8192 ulp | most values, up to 5120 ulp |
| **score at `%.3f`, every VCF** | **0/71; 32 cases, 337 rows equal** | **0/71** |

The weights the job extracted have the two SHA-256 that `weights.rs` pins; any other file is refused.

So under AVX-512 the reference's own 1D convolutions take another order than under AVX2: its bits
there depend on the CPU it runs on, not only on the image. The port takes the AVX2 order on every
host. The layer test now asserts the convolutions under a dispatch where they were measured exact
(`jit:avx2` for 1D; `brgconv:avx2` and AVX-512 for 2D), and it tries ten candidate orders for each 1D
convolution from the reference's own previous layer. Under AVX2 only `jit:avx2`'s order misses
nothing, and each of the others misses 25% to 95% of the values. The next run on an AVX-512 host
shows which order that dispatch takes, if any of the ten.

## The limit

These bits are not a property of the reference, because the reference does not hold them still:

* **From the first `SGEMM` on** (the dense layers, the logits, the softmax, the score's bits): they
  move with the thread count (2047 of 2840 values of dense_2.0 between one thread and four), with
  the batch size and with a row's place in its batch. MKL's kernel is closed and chosen per shape and
  per ISA. The port takes the order of one thread at batch 64 under AVX2, and it matches the first
  batch there exactly.
* **The convolutions across CPUs**: AVX2 and AVX-512 hosts give different bits from the same image.

What the tool itself keeps constant is the printed `%.3f` score and the VCF. The port matches those on
every run measured: Rosetta at AVX2 and a real AVX-512 host, both models, 32 cases.

## Running it

The weights are 8 MB and not committed. `tools/nvscorevariants/extract_weights.py` takes them from the
jar once, and `weights.rs` pins each by SHA-256 and refuses any other. Point `GATK_RS_NVSCORE_WEIGHTS`
at the directory, or put them beside the executable. Without them the tool stops where
`interpreter-tools` says it does. `.vcf.gz` and `.bcf` output, BCF input and GPU accelerators are
refused or fail in the port's own words.

For #1021 and #250: no `tch` or runtime is needed, so that evaluation can stop at this result. The cost
is speed, because the sums are scalar software FMA: on this machine, under load from other jobs and so
only roughly, the port took 3.5 s for the 71 rows in 1D and 45 s in 2D. No exclusive timing of port and
reference has been taken. That is the next thing #1021 would want.

#1021 has since been measured (`docs/gpu/targets.md`, §1b): encoding is under 0.1% of the tool and
inference 86%, and the convolutions now also run on Metal, bit-identical to the CPU path above and
so to the reference's layers, at about 4 times the CPU path on every core.

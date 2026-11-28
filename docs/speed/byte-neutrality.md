# A change made for speed may not move a byte

Milestone S, S.2 (#109). The rule is one sentence; this page is what makes it a gate rather than
an intention, what that gate was shown to catch and to miss, and the list of paths that are slow
because of the claim rather than in spite of it.

## The gate has two halves

**The suites.** Every pull request re-runs every conformance suite and every golden-backed test:
`ci.yml` triggers on `pull_request` with no path filter, and `CI complete` is the context `main`
requires. A performance change is not exempted for being one. A golden that moves fails the PR.

**The arithmetic guard**, `tools/audit/arithmetic.py`, in the same CI job as the provenance guard
and in `tools/preflight.py`. A suite sees only the orders and the inputs its corpus reaches, so
the guard refuses the constructs that make reordering possible before a suite has to notice:

| Refused | Where | Why |
|---|---|---|
| `mul_add` | Rust source | one rounding where the reference rounds twice |
| `fadd_fast` and family, `algebraic_*`, `core_intrinsics` | Rust source | they license reassociation |
| `std::simd` | Rust source | a vector reduction is a different summation tree |
| `std::arch`, `target_feature`, `is_*_feature_detected` | Rust source, outside the files `tools/audit/simd_kernels.json` declares | explicit SIMD is admitted only lane-wise, under classification (below) |
| fused, horizontal and approximate intrinsics (`fmadd`, `vfma`, `hadd`, `_reduce_`, `vaddv`, `_rcp`, ...) and enabling `fma` | the declared SIMD files | each rounds once, sums across lanes or approximates |
| `#[global_allocator]` | Rust source, outside `gatk-cli/src/main.rs` | admitted only in the shipped binary (below) |
| `rustflags`, `RUSTFLAGS`, `target-cpu`, `llvm-args`, `CARGO_PROFILE_*`, any `[profile]` but the one below | `.cargo/config*`, every `Cargo.toml`, the workflows | LLVM already folded `powf` of a constant differently between a debug and a release test run once, so a profile is admitted only by name and key, and only once every suite runs on it |
| `-ffast-math`, `-Ofast`, `-ffp-contract=fast` | `build.rs` | a C dependency compiled with fast-math reorders too |
| a new host transcendental (`exp`, `ln`, `log10`, `powf`, `powi`, ...) | Rust source, per file | the platform libm is not the JVM's |

The last row is a ratchet, not a ban: `tools/audit/arithmetic.json` records 155 call sites in 61
files today, some of them deliberate (`pow10` measured closer to `Math.pow` with `powf` than with
`strict_math::pow`; see `docs/numeric-functions-a-ported-call-site-reaches.md`). A file may not
gain one unless `--record` rewrites the count in the same change, where a reviewer sees it.

`arithmetic.py --self-test` plants each refused construct, one at a time, in a scratch copy of the
tree and fails unless every one is reported, and fails if the clean tree is not clean. CI runs it
next to the check, so the guard is verified in both directions on every PR, as the provenance
guard was (#82).

## What is admitted under classification (the S.2 amendment, W12)

The rule above refused every profile, every allocator and every intrinsic. The user's call for
extreme performance amended it: three things are admitted, each by name, each with the argument
for why it cannot move a byte and the suite that would catch it if the argument were wrong. The
guard enforces the names; the arguments are here.

**One profile, `[profile.dist]`** in the workspace `Cargo.toml`, exactly `inherits = "release"`,
`lto = "fat"`, `codegen-units = 1`; any other table or key is refused. Neither setting gives LLVM
a licence it does not already have: Rust emits no fast-math flags and never contracts
`a * b + c`, so inlining across crates and allocating registers over one unit cannot reorder,
reassociate or fuse an operation. What LTO does widen is constant folding: a host transcendental
(`powf`, `ln_1p`, ... , the 155 ratcheted call sites) whose argument becomes a constant after
cross-crate inlining may be folded at compile time with the build host's libm instead of called at
run time. That is the failure mode recorded above for debug against release, and the reason the
port already hides its constant bases behind `std::hint::black_box` (`pair_hmm::jacobian_entry`).
The proof is therefore not the argument but the suites run on the `dist` binary: every job that
runs the port against the reference (the covering arrays, the fuzzer replay, thread independence,
`Speed`) builds it with `--profile dist`; `Speed at scale` (#1493) still builds `release` and must
switch, so a folded constant that moves a byte fails the
array that reaches it. `release` stays the default profile, because the 391 test binaries would
each pay a fat link, and the in-process suites keep running on it: both code generations are held
to the same goldens.

**A global allocator, in the shipped binary only** (`crates/gatk-cli/src/main.rs`). An allocator
decides addresses and nothing else. An address reaches an output only where the code hashes,
orders or prints a pointer; address-space layout randomisation already moves every address on
every run, so such an output would already fail the suites run to run, and the one known case
(`BwaMemIndexImageCreator`, [pointers-that-reach-the-output.md](../pointers-that-reach-the-output.md))
is masked by layout. What a new allocator can change that ASLR cannot is the RELATIVE order of two
allocations; that is checked by the `determinism` job and by every end-to-end test, which spawn the
binary and therefore run its allocator.

**`std::arch`, in the files `tools/audit/simd_kernels.json` declares**, the way
`gpu_kernels.json` declares device kernels: each entry says what the code does and names the test
that holds every one of its paths, each forced in turn, to the scalar function. The only admitted
shape is lane-wise: one independent problem per lane, each lane doing the scalar operations in the
scalar order. A lane-wise IEEE 754 multiply or add is the scalar multiply or add, bit for bit, so
such a kernel is identical by construction, whatever width the host dispatches; that is the
eligibility rule of `docs/gpu/eligibility.md` applied to the CPU. Inside a declared file the guard
still refuses fused multiply-adds, horizontal operations (a reduction across lanes sums in its own
tree), approximate reciprocals and enabling the `fma` feature. `target-cpu` stays refused: a
declared kernel dispatches at run time, so one binary runs everywhere and every path it can take
is tested.

## What the suites caught, and what they did not

Three deliberately byte-moving "optimisations", each run against the suite that covers it
(`cargo test --release`, macOS on Apple silicon, 2026-10-09):

| Change | Suite | Result |
|---|---|---|
| `somatic_likelihoods::sum` iterated in reverse | `somatic-likelihoods` | **passed**. Every sum the corpus reaches has two terms, and `a + b` is `b + a` exactly. |
| `effective_counts` accumulating over reads in reverse (`sumArrayFunction`'s order) | `somatic-likelihoods` | **failed**: `ten-reads-split[0]: 2 ulp apart`. |
| `natural_log_utils` calling the host's `f64::exp` instead of `jmath::strict_math::exp` | `natural-log-utils`, `somatic-likelihoods` | **passed** on this host. Apple's libm agrees with FDLIBM on every point the corpus visits; glibc on the runner need not. |

The first and third are why the guard exists. A reordering the corpus cannot observe is still a
reordering, and a libm substitution that happens to pass on one host is the failure decision 0007
of htsjdk-rs describes for `pow` across x86 CPUs. The guard refuses the third outright now (a new
`.exp(` in `natural_log_utils.rs` fails the ratchet); the first has no syntactic signature, which is
what the suites are for, and the honest statement is that a two-term corpus does not pin a sum's
order. Widening a corpus is the fix for that, not a rule.

A benchmark that cannot be made byte-neutral is not merged with a quarantine. It is not merged.

## Slow on purpose

These paths cost what they cost because of the claim. They are not optimisation targets, and
each names what pins it, so nobody spends a week rediscovering that the slow loop is slow on
purpose.

| Path | The cost | Pinned by |
|---|---|---|
| `somatic_likelihoods::sum`, `effective_counts`, and every port of `MathUtils.sum` / `sumArrayFunction` | accumulate in index order, one term at a time: no vectorised reduction, no pairwise tree | suite `somatic-likelihoods` (the 2-ulp failure above) |
| `jmath::strict_math::exp`, `pow`, and the `jmath::math` functions | FDLIBM in software where the JVM uses a hardware-assisted intrinsic | suites `natural-log-utils`, `somatic-likelihoods`, `mutect-engine-arithmetic`; htsjdk-rs decisions 0014, 0025, 0027 |
| `gatk-annotation/src/decimal_format.rs`, `exact_decimal` | the double's exact decimal expansion, up to 1,074 digits, so `DecimalFormat`'s rounding is decided on the value and not on a shortest representation | suite `decimal-format`; htsjdk-rs decision 0017 |
| `gatk-engine/src/java_hash.rs` (`JavaHashMap`) | a hash table kept in the JVM's bucket order so iteration order matches, instead of the fastest map | `docs/an-unspecified-order-that-reaches-the-output.md`; suites `annotation-depth-per-allele`, `annotation-fragment-counts`, `trim-alleles`, `acbuilder` |
| `gatk-engine/src/java_random.rs`, `well19937c.rs` | the reference's generators, sequence for sequence, instead of a faster one | suites `javarandom`, `well19937c`, `reservoir` |

What is **not** on this list belongs to S.4 (#111): costs the port pays that the claim does not
require, such as `AllelePseudoDepth` cloning rows the reference aliases. Those are candidates;
these are not.

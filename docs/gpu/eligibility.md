# What makes a kernel eligible

Milestone GPU, GPU.1 (#76). Enforced per kernel by `tools/audit/gpu_kernels.py`, which reads its
record from `tools/audit/gpu_kernels.json`.

GATK 4.6.2.0 has no GPU path, so a device kernel here has no vendor output to agree with. Its only
reference is the CPU path of this port, which is already byte-identical to GATK. A kernel ships if
it gives that path's bytes, on every device it runs on, and not otherwise. "Close", "concordant"
and "within tolerance" are not grades this milestone hands out.

## The rule

A kernel belongs to exactly one class.

**`integer`: eligible as it stands.** Integer and byte work with a fixed traversal order. Integer
addition and multiplication wrap identically on every device, so the order of independent
operations cannot change a result, and no rounding mode, denormal mode or vendor library is
involved. The guard refuses every floating-point type in such a kernel, including in a helper it
calls.

**`float-pinned`: eligible only when all of these hold**, and the record says how:

- the operation order is the CPU path's, operation for operation;
- no fast math, and no contraction of `a * b + c` into a fused multiply-add unless the CPU path
  contracts the same expression (Rust never contracts on its own; the port never asks it to, see
  `tools/audit/arithmetic.py`);
- no TF32, tensor-core or reduced-precision substitution;
- every reduction has a fixed tree that does not depend on scheduling, atomics order, SIMD width
  or occupancy;
- every transcendental the CPU path takes from `jmath` (`StrictMath`'s and `Math`'s functions) is
  computed in software on the device, because a vendor libm is not the one `jmath` reproduces, or
  stays on the host.

**`quarantined`: never ships.** Anything that cannot meet the above. It is reported as
bio-identical, exactly as a CPU quarantine is, and no tool's byte-identity claim may rest on it.

On every class the guard also refuses `fma`, `mad`, the `fast::` and `precise::` namespaces,
SIMD-group and quad-group operations, atomics and threadgroup memory, unless the record lists the
construct for that source with a reason, and it refuses host code that compiles a kernel source
without turning fast math off.

## Why no kernel here is floating point

Apple GPUs, the first hardware this milestone runs on, have no binary64 unit, and `LoglessPairHMM`
is binary64 throughout. A `float` kernel would be a different function. So `kernels.metal` carries
every double as the `ulong` holding its IEEE 754 bits and implements addition and multiplication
on integers: rounded to nearest-even, subnormals kept, NaNs propagated in AArch64's order. That
turns the one floating-point target of GPU.2 into an `integer` kernel, which is the strongest of
the three classes rather than the weakest: there is no floating-point unit left for a second
device to disagree with.

The arithmetic is held to the host's by `crates/gatk-engine/tests/gpu.rs`: four hundred thousand
additions and as many multiplications over random bit patterns, operands one ulp apart, opposite
signs that cancel, and every pair of zeros, subnormals, extremes, infinities and NaNs. The single
exception is the bit pattern of a NaN produced by an invalid operation, which is not portable on
the CPU either (positive on AArch64, negative on x86-64); PairHMM refuses a NaN likelihood, so that
NaN never reaches an output.

## Determinism across architectures

GPU.1 asks for a determinism gate across two different GPU architectures, because a reduction that
is deterministic at one SIMD width or occupancy can stop being so at another. For the kernels of
class `integer` the argument is structural rather than empirical: the guard refuses every
width-dependent construct (SIMD-group operations, atomics, threadgroup memory), each thread reads
only its own pair, and integer arithmetic has no mode to differ in. A second architecture can still
expose a compiler bug, so the gate is run, not assumed: the CI job runs the same tests on GitHub's
macOS runner, and `docs/gpu/devices.md` records each device the bit-equality was observed on.
Until two architectures appear there, that box stays open.

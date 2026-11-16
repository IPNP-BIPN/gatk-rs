# The devices the bytes were observed on

Milestone GPU, GPU.1's determinism gate (#76). One row per device on which
`crates/gatk-engine/tests/gpu.rs` and the whole workspace were run on the device path, with what
`tools/provenance/record.py` records for it. GPU.1 asks for two different architectures before the
gate is closed; a row is added from a run, never from a datasheet.

| Device | Architecture | Driver (OS build) | Bindings | Kernel tests | Workspace on the device path | Date |
|---|---|---|---|---|---|---|
| Apple M4 Max, 40 cores | Apple GPU, Metal 4 | macOS 27.2 (26B123) | `metal` 0.33.0 | 800,648 operations and 732 pairs bit-identical | 391 test binaries; 2,085 pass and the same 3 host-bound ones fail as on the CPU path | 2026-11-05 |

The second row is meant to come from the `gpu` job of `.github/workflows/ci.yml` on GitHub's
macOS runner. That runner is a virtual machine on Apple silicon whose GPU is paravirtualised, so
even when it runs, it is a second driver stack rather than a second architecture: the gate stays
open until a device of another family (AMD or Intel on an Intel Mac, or a non-Apple GPU behind a
second backend) has a row.

## The likeliest second backend: CUDA from Rust

NVIDIA's CUDA Rust (github.com/NVIDIA/cuda-rust, announced 2026-09-08) compiles Rust kernels to
PTX: `cuda-oxide` per thread, through a custom `rustc` backend, and `cutile-rs` per tile. The
per-thread track fits these kernels, which are one pair or one output element per thread, and
would let the device code be the CPU path's own Rust rather than a Metal transliteration. What it
needs, none of which is here:

- an NVIDIA GPU of compute capability 8.0 or later on Linux. GitHub's hosted GPU runners are T4s
  (7.5), below that line, so the row would come from a self-hosted or rented A10, L4 or A100;
- a second toolchain: `cuda-oxide` pins its own `rustc` (1.99.1 at the time of writing) with
  `rustc-dev`, against this workspace's 1.97.1, so it would build as a separate crate;
- an eligibility check on the emitted PTX, because the project says nothing yet about contraction
  or fast math: refuse `fma.` and any `add`, `sub` or `mul` on `f32` or `f64` without the `.rn`
  modifier, which is what stops `ptxas` from contracting. With that, CUDA's native binary64 would
  replace the software one, and the kernels would become `float-pinned` rather than `integer`.

It describes itself as an early-stage research project, so it is a candidate for the gate's second
row, not a dependency of anything that ships.

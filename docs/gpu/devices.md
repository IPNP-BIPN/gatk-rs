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

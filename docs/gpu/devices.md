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

## The second architecture: an NVIDIA T4, through CUDA

`crates/gatk-engine/src/gpu/cuda.rs` is the Linux backend. Its kernels are not a second
implementation: `tools/gpu/metal_to_cuda.py` derives `kernels.cu` from `kernels.metal` by rewriting
only the address-space qualifiers, the parameter attributes and `clz`, so the arithmetic is the
Metal file's line for line, and the `guard` job refuses a `kernels.cu` that has drifted from it.
NVRTC compiles it at run time with `--fmad=false` and no fast math; the driver and NVRTC are
loaded with `dlopen`, so a host without them runs the CPU path. The kernels are integer only
(the PTX for sm_75 holds no floating-point instruction, only integer `mad.lo.s32`), so they
need nothing newer than Turing: a T4 (sm_75) is enough. sm_80 was cuda-oxide's floor, not theirs.

`.github/workflows/gpu-gate.yml` runs `gpu_speed`, then every suite on the CPU and the CUDA path
on the organisation's T4 runner, as the `gpu` job does on macOS. It runs when a kernel, its host or
the workflow changes, and by hand, because that runner is billed per minute. Setting it up is an
organisation setting with billing, so it is done by an owner:

1. The organisation needs GitHub Team or Enterprise Cloud: larger runners are not on the Free plan.
2. Organisation settings, Actions, Runners, New runner, New GitHub-hosted runner:
   - name `gpu-t4` (the workflow's `runs-on` label is the runner's name);
   - platform Linux x64, image: Partner tab, "NVIDIA GPU-Optimized Image for AI and HPC" (it
     carries the driver, and the GPU sizes appear only once it is chosen);
   - size GPU-powered, 4 cores, 28 GB, one Tesla T4;
   - maximum concurrency 1, which caps the spend at one job at a time;
   - runner group Default, with "Allow public repositories" ticked, since gatk-rs is public.
3. Run the workflow once by hand (Actions, GPU gate, Run workflow) and add its row here.

It costs $0.07 a minute, always billed, included minutes do not apply. The macOS `gpu` job takes
26 minutes on three cores; on four, with the cache warm, a run should take 30 to 45 minutes, so
about $2 to $3, and up to about $4 for the first, cold one.

#!/usr/bin/env python3
"""Derive crates/gatk-engine/src/gpu/kernels.cu from kernels.metal, so the CUDA backend runs the
same integer arithmetic, line for line, as the Metal one (Milestone GPU.1's second architecture).

    metal_to_cuda.py           write kernels.cu
    metal_to_cuda.py --check   exit 1 if kernels.cu is not what kernels.metal gives
"""

import pathlib
import re
import sys

GPU = pathlib.Path(__file__).resolve().parents[2] / "crates/gatk-engine/src/gpu"
PRELUDE = """\
// Generated from kernels.metal by tools/gpu/metal_to_cuda.py; edit that file, not this one.
// Compiled by NVRTC with --fmad=false (crates/gatk-engine/src/gpu/cuda.rs).

typedef unsigned long ulong; // 64 bits: CUDA follows the LP64 host
typedef unsigned int uint;
typedef unsigned char uchar;
__device__ static inline uint min_u32(uint a, uint b) { return a < b ? a : b; }
"""


def convert(metal: str) -> str:
    text = metal.replace("#include <metal_stdlib>\nusing namespace metal;\n", PRELUDE)
    text = re.sub(r"^constant ulong ", "__device__ constexpr ulong ", text, flags=re.M)
    text = re.sub(r"^static ", "static __device__ ", text, flags=re.M)
    text = text.replace("thread int &", "int &").replace("thread ulong &", "ulong &")
    text = re.sub(r"\bclz\(", "__clz((int)", text)
    text = re.sub(r"\bmin\(", "min_u32(", text)
    # Kernel parameters: buffers become pointers, `constant T &` words are passed by value, and
    # the thread index is computed from the grid, each thread still checking its own bound.
    text = re.sub(r"device (const )?(\w+) \*(\w+) \[\[buffer\(\d+\)\]\]", r"\1\2 *\3", text)
    text = re.sub(r"constant uint &(\w+) \[\[buffer\(\d+\)\]\]", r"uint \1", text)
    text = re.sub(
        r",\s*uint (\w+) \[\[thread_position_in_grid\]\]\) \{",
        r") {\n    uint \1 = blockIdx.x * blockDim.x + threadIdx.x;",
        text,
    )
    text = re.sub(r"^kernel void ", 'extern "C" __global__ void ', text, flags=re.M)
    text = re.sub(r"\bdevice (const )?(\w+) \*", r"\1\2 *", text)
    code = re.sub(r"//[^\n]*", "", text)
    leftover = re.findall(r"\[\[|\bthread\b|\bdevice\b|\bconstant\b|\bmetal\b", code)
    assert not leftover, f"Metal left in the CUDA source: {sorted(set(leftover))}"
    return text


def main() -> int:
    source = convert((GPU / "kernels.metal").read_text())
    target = GPU / "kernels.cu"
    if sys.argv[1:] == ["--check"]:
        if not target.exists() or target.read_text() != source:
            print(f"{target} is stale: run tools/gpu/metal_to_cuda.py")
            return 1
        return 0
    target.write_text(source)
    print(f"wrote {target}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

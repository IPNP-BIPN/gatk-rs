#!/usr/bin/env python
"""Prints the CPU as comments, then `file<TAB>path<TAB>size<TAB>sha256` for every file of a run,
then `input<TAB>path<TAB>hex` for every shard file, which crates/gatk-cli/tests/
gcnv_postprocess_end_to_end.rs writes back to replay PostprocessGermlineCNVCalls.

The pytensor cache is left out: it holds the compiled modules, whose keys name the host.
Usage: dump_gcnv.py <directory measure_gcnv.sh wrote>
"""
import hashlib
import os
import sys

import numpy as np

root = sys.argv[1]
model = next((line.split(":", 1)[1].strip() for line in open("/proc/cpuinfo") if line.startswith("model name")), "?")
flags = next((line.split(":", 1)[1].split() for line in open("/proc/cpuinfo") if line.startswith("flags")), [])
print(f"# cpu {model}; " + " ".join(f for f in ("avx", "fma", "avx2", "avx512f", "avx512_vnni", "amx_tile") if f in flags))
print("# numpy " + " ".join(k for k, v in np.core._multiarray_umath.__cpu_features__.items() if v))
for name in ("PYTENSOR_FLAGS", "PYTENSOR_EXTRA", "MKL_CBWR", "MKL_NUM_THREADS", "NPY_DISABLE_CPU_FEATURES", "GLIBC_TUNABLES", "LD_PRELOAD"):
    print(f"# env {name}={os.environ.get(name, '')}")
paths = []
for directory, subdirs, names in os.walk(root):
    subdirs[:] = [d for d in subdirs if d != ".pytensor"]
    paths += [os.path.join(directory, n) for n in names]
for path in sorted(paths):
    with open(path, "rb") as f:
        content = f.read()
    print(f"file\t{os.path.relpath(path, root)}\t{len(content)}\t{hashlib.sha256(content).hexdigest()}")
for path in sorted(p for p in paths if os.path.relpath(p, root).split(os.sep)[0] in ("gcnv", "ploidy")):
    with open(path, "rb") as f:
        print(f"input\t{os.path.relpath(path, root)}\t{f.read().hex()}")

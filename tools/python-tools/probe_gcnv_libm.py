"""Inside broadinstitute/gatk:4.7.0.0 under gcnv_frozen_env.sh: locate glibc's __log_data, test
numpy's sum order, and write libm test vectors (inputs, numpy's outputs) as raw little-endian
float64, which crates/gatk-tools/src/glibc_libm.rs's tests read from $GCNV_LIBM_VECTORS.

Usage: probe_gcnv_libm.py <output directory>
"""
import struct
import sys

import numpy as np

out = sys.argv[1]
lib = open("/lib/x86_64-linux-gnu/libm.so.6", "rb").read()
pat = struct.pack("<Q", 0x3fe62e42fefa3800)
pos = 0
with open(f"{out}/log_data_candidates.txt", "w") as f:
    while True:
        pos = lib.find(pat, pos)
        if pos < 0:
            break
        vals = struct.unpack_from("<20d", lib, pos)
        f.write(f"offset {pos:#x}: " + " ".join(v.hex() for v in vals) + "\n")
        pos += 8
open(f"{out}/libm.so.6", "wb").write(lib)


def pairwise(a):
    n = len(a)
    if n < 8:
        res = 0.0
        for x in a:
            res += x
        return res
    if n <= 128:
        r = list(a[:8])
        i = 8
        while i < n - (n % 8):
            for j in range(8):
                r[j] += a[i + j]
            i += 8
        res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]))
        while i < n:
            res += a[i]
            i += 1
        return res
    n2 = n // 2
    n2 -= n2 % 8
    return pairwise(a[:n2]) + pairwise(a[n2:])


rng = np.random.default_rng(5)
counts = {"all": 0, "first": 0, "neither": 0}
for n in list(range(1, 300)) * 20:
    a = rng.random(n) * 10.0 ** rng.integers(-5, 5, n)
    s = float(np.sum(a))
    full = pairwise([float(x) for x in a])
    first = float(a[0]) + pairwise([float(x) for x in a[1:]])
    if s == full:
        counts["all"] += 1
    if s == first:
        counts["first"] += 1
    if s != full and s != first:
        counts["neither"] += 1
print("sum order", counts)

m = 400000
inputs = {
    "exp": np.concatenate([rng.uniform(-745, 709, m), rng.uniform(-40, 1, m), rng.uniform(-1e-3, 1e-3, m // 4)]),
    "log": np.concatenate([np.exp(rng.uniform(-700, 700, m)), rng.uniform(0.9, 1.1, m), rng.uniform(0, 64, m)]),
    "log1p": np.concatenate([-np.exp(rng.uniform(-745, 0, m)), np.exp(rng.uniform(-745, 0, m)), rng.uniform(-0.5, 2, m)]),
    "expm1": np.concatenate([rng.uniform(-40, 0, m), -np.exp(rng.uniform(-60, 0, m)), rng.uniform(-1, 1, m)]),
}
for name, x in inputs.items():
    y = getattr(np, name)(x)
    x.astype("<f8").tofile(f"{out}/{name}.in")
    y.astype("<f8").tofile(f"{out}/{name}.out")
a = rng.uniform(-30, 0, m)
b = rng.uniform(-30, 0, m)
b[: m // 10] = a[: m // 10]
np.stack([a, b, np.logaddexp(a, b)], axis=1).astype("<f8").tofile(f"{out}/logaddexp.bin")
print("done")

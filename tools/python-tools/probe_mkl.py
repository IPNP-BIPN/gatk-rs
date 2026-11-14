#!/usr/bin/env python
"""Does MKL's answer depend on where numpy put the operands? Distinct hashes per product.

The same operands, and the output, are placed in buffers shifted by 0 to 7 doubles from a 64-byte
boundary, and dgemm, dgemv and ddot run on each placement. One hash per product means alignment
does not matter.
"""
import hashlib

import numpy as np

rng = np.random.RandomState(1984)
a0 = rng.normal(size=(300, 257))
b0 = rng.normal(size=(257, 301))
x0 = rng.normal(size=100003)


def placed(arr, shift):
    raw = np.empty(arr.size + 64, dtype=np.float64)
    start = (-(raw.ctypes.data // 8) % 8) + shift
    out = raw[start:start + arr.size].reshape(arr.shape)
    out[...] = arr
    return out


def h(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()[:12]


seen = {"dgemm_out": set(), "dgemm": set(), "dgemv": set(), "ddot": set(), "dgemm_small": set()}
for shift in range(8):
    for other in range(8):
        a, b = placed(a0, shift), placed(b0, other)
        seen["dgemm"].add(h(a @ b))
        seen["dgemv"].add(h(a @ placed(b0[:, 0].copy(), other)))
        seen["dgemm_small"].add(h(placed(a0[:, :5].copy(), shift) @ placed(b0[:5, :].copy(), other)))
    out = placed(np.zeros((300, 301)), shift)
    np.matmul(a0, b0, out=out)
    seen["dgemm_out"].add(h(out))
    x = placed(x0, shift)
    seen["ddot"].add(h(np.array([x @ x])))
for k, v in seen.items():
    print(k, len(v), " ".join(sorted(v)))

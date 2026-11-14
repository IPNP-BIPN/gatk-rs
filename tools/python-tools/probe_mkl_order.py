#!/usr/bin/env python
"""In which order does MKL add the products of a dgemm, for the shapes gCNV hands it?

For each (m, k, n), A @ B through numpy (MKL), against candidate orders computed exactly with numpy
elementwise operations: one accumulator from k=0 up, from k=K-1 down, p accumulators over k mod p
summed in a few ways, and an FMA chain (exact via fractions). Prints the share of entries each
candidate reproduces. Run under gcnv_frozen_env.sh.
"""
from fractions import Fraction

import numpy as np

rng = np.random.RandomState(20261105)


def chain(a, b, ks):
    s = np.zeros((a.shape[0], b.shape[1]))
    for k in ks:
        s = s + a[:, k:k + 1] * b[k:k + 1, :]
    return s


def split(a, b, p, combine):
    K = a.shape[1]
    parts = [chain(a, b, range(r, K, p)) for r in range(p)]
    return combine(parts)


def seq(parts):
    s = parts[0]
    for q in parts[1:]:
        s = s + q
    return s


def tree(parts):
    while len(parts) > 1:
        parts = [parts[i] + parts[i + 1] if i + 1 < len(parts) else parts[i] for i in range(0, len(parts), 2)]
    return parts[0]


def fma_chain(a, b):
    out = np.zeros((a.shape[0], b.shape[1]))
    for i in range(a.shape[0]):
        for j in range(b.shape[1]):
            s = 0.0
            for k in range(a.shape[1]):
                s = float(Fraction(a[i, k]) * Fraction(b[k, j]) + Fraction(s))
            out[i, j] = s
    return out


for m, k, n in [(6, 3, 3), (180, 5, 6), (6, 5, 180), (5, 6, 180), (180, 6, 5), (5, 180, 6), (6, 180, 5), (1, 180, 1), (180, 1, 6)]:
    a = rng.normal(size=(m, k))
    b = rng.normal(size=(k, n))
    got = a @ b
    cands = {
        "up": chain(a, b, range(k)),
        "down": chain(a, b, range(k - 1, -1, -1)),
        "fma-up": fma_chain(a, b) if m * n * k <= 20000 else None,
    }
    for p in (2, 4, 8):
        if k > p:
            cands[f"split{p}-seq"] = split(a, b, p, seq)
            cands[f"split{p}-tree"] = split(a, b, p, tree)
    line = " ".join(f"{name}={np.mean(c == got):.3f}" for name, c in cands.items() if c is not None)
    print(f"m={m} k={k} n={n}: {line}")

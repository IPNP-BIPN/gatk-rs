#!/usr/bin/env python
"""Hashes of what the knobs of gcnv_frozen_env.sh are meant to pin, to see each knob take effect.

Run twice inside broadinstitute/gatk:4.7.0.0, with and without the knob, and compare the lines.
libm: glibc's exp, log, pow, log1p, expm1 through ctypes; numpy's own ufuncs; MKL through numpy's
matmul (libblas.so.3 is MKL) and through a pytensor function compiled with the session's flags.
"""
import ctypes
import hashlib
import os
import sys

import numpy as np

n = 200000
rng = np.random.RandomState(1984)
x = rng.uniform(-30.0, 30.0, n)
y = rng.uniform(0.01, 5.0, n)
libm = ctypes.CDLL("libm.so.6")


def h(a):
    return hashlib.sha256(np.ascontiguousarray(a, dtype=np.float64).tobytes()).hexdigest()[:16]


def via_libm(name, *args):
    f = getattr(libm, name)
    f.restype = ctypes.c_double
    f.argtypes = [ctypes.c_double] * len(args)
    return np.array([f(*v) for v in zip(*args)])


print("# cpu " + " ".join(k for k, v in np.core._multiarray_umath.__cpu_features__.items() if v))
m = 20000
print("libm.exp", h(via_libm("exp", x[:m])))
print("libm.log", h(via_libm("log", y[:m])))
print("libm.pow", h(via_libm("pow", y[:m], x[:m] / 7)))
print("libm.log1p", h(via_libm("log1p", y[:m])))
print("libm.expm1", h(via_libm("expm1", x[:m] / 10)))
print("numpy.exp", h(np.exp(x)))
print("numpy.log", h(np.log(y)))
print("numpy.sum", h(np.add.reduce(x.reshape(400, 500), axis=1)))
a = rng.normal(size=(300, 257))
b = rng.normal(size=(257, 301))
print("mkl.dgemm", h(a @ b))
print("mkl.dgemv", h(a @ b[:, 0]))
print("mkl.ddot", h(np.array([x @ y])))
if len(sys.argv) > 1:
    import pytensor
    import pytensor.tensor as pt

    u, v = pt.dvectors("u", "v")
    w = pt.dmatrix("w")
    f = pytensor.function([u, v, w], [u * v + pt.exp(u / 10) * v - pt.log1p(v), pt.dot(w, v[:257]),
                                      pt.logsumexp(pt.outer(u[:300], v[:300]), axis=1)])
    for name, out in zip(["pt.elemwise", "pt.gemv", "pt.logsumexp"], f(x, y, a)):
        print(name, h(out))
    from pytensor.link.c.cmodule import GCC_compiler

    print("# compile_args " + " ".join(GCC_compiler.compile_args()), "openmp", pytensor.config.openmp)

"""Runs segment_gcnv_calls.py for one sample with every compiled-function call and every segment quality
traced as raw float64 bits, which crates/gatk-tools/tests/gcnv_segments_trace.rs replays
(docs/gcnv-postprocess.md). Inside broadinstitute/gatk:4.6.2.0, after gcnv_frozen_env.sh.

Usage: pp_trace.py <segment_gcnv_calls.py> <trace out> <graph out> <the script's arguments...>
"""
import os
import sys

script, trace_path, graph_path = sys.argv[1:4]
user = os.environ.get("PYTENSOR_FLAGS")
default = "device=cpu,floatX=float64,optimizer=fast_run,compute_test_value=ignore," + \
          "openmp=true,blas__ldflags=-lmkl_rt,openmp_elemwise_minsize=10,exception_verbosity=high"
os.environ["PYTENSOR_FLAGS"] = default + ("" if user is None else "," + user)

import numpy as np
import pytensor
from pytensor.compile.function import types as ftypes

trace = open(trace_path, "w")
names = {}


def bits(value):
    a = np.asarray(value)
    if a.dtype.kind == "f":
        return a.shape, " ".join("%016x" % v for v in a.astype(np.float64).ravel().view(np.uint64))
    return a.shape, " ".join(str(int(v)) for v in a.ravel())


def emit(key, value):
    shape, text = bits(value)
    trace.write(f"{key}\t{','.join(map(str, shape))}\t{text}\n")


original_call = ftypes.Function.__call__


def call(self, *args, **kwargs):
    if id(self) not in names:
        names[id(self)] = len(names)
        with open(graph_path, "a") as g:
            g.write(f"=== function {names[id(self)]}\n")
            pytensor.printing.debugprint(self, file=g, print_type=True)
    k = names[id(self)]
    out = original_call(self, *args, **kwargs)
    for i, a in enumerate(args):
        emit(f"f{k}.in{i}", a)
    for i, o in enumerate(out if isinstance(out, (list, tuple)) else [out]):
        emit(f"f{k}.out{i}", o)
    return out


ftypes.Function.__call__ = call

from gcnvkernel.postprocess import viterbi_segmentation as vs
from gcnvkernel.postprocess import segment_quality_utils as squ

for method in ["get_segment_quality_some_called", "get_segment_quality_all_called",
               "get_segment_quality_start", "get_segment_quality_end"]:
    def wrap(m=method, f=getattr(squ.HMMSegmentationQualityCalculator, method)):
        def w(self, *a, **kw):
            r = f(self, *a, **kw)
            emit(f"q.{m}.{','.join(map(str, a))}", np.float64(r))
            return r
        return w
    setattr(squ.HMMSegmentationQualityCalculator, method, wrap())

original_init = vs.ViterbiSegmentationEngine.__init__


def init(self, *a, **kw):
    original_init(self, *a, **kw)
    emit("log_q_tau_tk", self.log_q_tau_tk)
    for contig in self.ordered_contig_list:
        emit(f"stay.{contig}", self.cnv_stay_prob_t_j[contig])


vs.ViterbiSegmentationEngine.__init__ = init

if user is None:
    del os.environ["PYTENSOR_FLAGS"]
else:
    os.environ["PYTENSOR_FLAGS"] = user
sys.argv = [script] + sys.argv[4:]
import runpy
runpy.run_path(script, run_name="__main__")
trace.close()

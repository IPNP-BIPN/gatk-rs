"""On PYTHONPATH: logs every pytensor function compiled, node by node, to $GCNV_GRAPH_LOG.

The patch waits for pytensor.compile.function.types to be imported, because GATK's scripts set
PYTENSOR_FLAGS before importing gcnvkernel and an early import would read the config too soon.
One line per node: pid, function number, op, how it runs (c, py, or the thunk's type).
"""
import importlib.abc
import importlib.util
import os
import sys

LOG = os.environ.get("GCNV_GRAPH_LOG")
TARGET = "pytensor.compile.function.types"


def describe(op):
    scalar = getattr(op, "scalar_op", None)
    name = type(op).__name__
    if scalar is not None:
        return f"{name}{{{scalar}}}"[:400]
    return f"{name}{{{op}}}"[:400] if name in ("Scan",) else name


def patch(module):
    original = module.FunctionMaker.create
    count = [0]

    def create(self, *args, **kwargs):
        fn = original(self, *args, **kwargs)
        count[0] += 1
        vm = getattr(fn, "vm", None) or getattr(fn, "fn", None)
        nodes = getattr(vm, "nodes", None) or self.fgraph.toposort()
        thunks = getattr(vm, "thunks", None) or [None] * len(nodes)
        with open(LOG, "a") as out:
            for node, thunk in zip(nodes, thunks):
                if thunk is None:
                    how = "?"
                elif hasattr(thunk, "cthunk"):
                    how = "c"
                elif hasattr(thunk, "perform"):
                    how = "py"
                else:
                    how = type(thunk).__name__
                shapes = ",".join(str(getattr(o.type, "ndim", "?")) for o in node.outputs)
                out.write(f"{os.getpid()}\t{count[0]}\t{describe(node.op)}\t{how}\t{shapes}\n")
        return fn

    module.FunctionMaker.create = create


class Hook(importlib.abc.MetaPathFinder):
    def find_spec(self, name, path, target=None):
        if name != TARGET:
            return None
        sys.meta_path.remove(self)
        spec = importlib.util.find_spec(name)
        run = spec.loader.exec_module

        def exec_module(module):
            run(module)
            patch(module)

        spec.loader.exec_module = exec_module
        return spec


if LOG:
    sys.meta_path.insert(0, Hook())

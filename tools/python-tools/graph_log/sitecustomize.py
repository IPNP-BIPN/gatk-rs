"""Log every pytensor function compiled, node by node, to $GCNV_GRAPH_LOG, when on PYTHONPATH.

The patch waits for pytensor.compile.function.types to be imported, because GATK's scripts set
PYTENSOR_FLAGS before importing gcnvkernel and an early import would read the config too soon.
One line per node: pid, function number, op, how it runs (c, py, or the thunk's type), and the
number of dimensions of each output. Nothing is patched when $GCNV_GRAPH_LOG is unset.
"""
import importlib.abc
import importlib.util
import os
import sys

LOG = os.environ.get("GCNV_GRAPH_LOG")
TARGET = "pytensor.compile.function.types"


def describe(op):
    """Return a short name for `op`: its type, with its scalar op, or its own text for a Scan.

    The text in braces is cut to 400 characters along with the name.
    """
    scalar = getattr(op, "scalar_op", None)
    name = type(op).__name__
    if scalar is not None:
        return f"{name}{{{scalar}}}"[:400]
    return f"{name}{{{op}}}"[:400] if name in ("Scan",) else name


def patch(module):
    """Wrap `module.FunctionMaker.create` so that every function it makes is logged."""
    original = module.FunctionMaker.create
    count = [0]

    def create(self, *args, **kwargs):
        """Make the function as before, then append one line per node of it to the log."""
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
    """An import hook that patches the pytensor module named by `TARGET` once it has loaded."""

    def find_spec(self, name, path, target=None):
        """Return the target's spec with a loader that patches the module after running it.

        Any other name is left to the other finders. The hook removes itself first, so it acts
        once and the real finder can be asked for the spec.
        """
        if name != TARGET:
            return None
        sys.meta_path.remove(self)
        spec = importlib.util.find_spec(name)
        run = spec.loader.exec_module

        def exec_module(module):
            """Run the module as the original loader would, then patch it."""
            run(module)
            patch(module)

        spec.loader.exec_module = exec_module
        return spec


if LOG:
    sys.meta_path.insert(0, Hook())

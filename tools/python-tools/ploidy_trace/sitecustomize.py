"""On PYTHONPATH with PLOIDY_TRACE=<dir>: traces cohort_determine_ploidy_and_depth.py for the port
(crates/gatk-tools/src/gcnv_ploidy.rs). Every compiled pytensor function's final graph goes to
graph.txt the first time it runs; every call's outputs and every shared variable it reads, after
the call, go to trace.tsv as raw float64 bits; reseed_rngs' seeds and the PCG64 states it sets too.

The patches wait for their modules to be imported: GATK's script sets PYTENSOR_FLAGS first.
"""
import importlib.abc
import importlib.util
import os
import shutil
import sys

OUT = os.environ.get("PLOIDY_TRACE")


def bits(value):
    import numpy as np
    a = np.asarray(value)
    if a.dtype.kind == "f":
        return a.shape, " ".join("%016x" % v for v in a.astype(np.float64).ravel().view(np.uint64))
    if a.dtype.kind in "iub":
        return a.shape, " ".join(str(int(v)) for v in a.ravel())
    return None, repr(value)


def emit(key, value):
    shape, text = bits(value)
    shape = "?" if shape is None else ",".join(map(str, shape))
    with open(os.path.join(OUT, "trace.tsv"), "a") as t:
        t.write(f"{key}\t{shape}\t{text}\n")


def state(value):
    import numpy as np
    if isinstance(value, np.random.Generator):
        s = value.bit_generator.state
        return "pcg64 %x %x %d %d" % (s["state"]["state"], s["state"]["inc"], s["has_uint32"], s["uinteger"])
    return None


def describe_nodes(nodes, out, prefix):
    """One line per node: its index, op, inputs (n<node>.<output>, a name, or a constant), and the
    C template of a fused composite; a Scan's inner graph follows under its own prefix."""
    from pytensor.graph.basic import Constant
    ids = {}
    for i, node in enumerate(nodes):
        for j, o in enumerate(node.outputs):
            ids[o] = f"{prefix}n{i}.{j}"
    for i, node in enumerate(nodes):
        args = []
        for v in node.inputs:
            if v in ids:
                args.append(ids[v])
            elif isinstance(v, Constant):
                args.append("const:" + repr(v.data).replace("\n", " "))
            else:
                args.append(f"in:{v.name}:{v.type}")
        scalar = getattr(node.op, "scalar_op", None)
        out.write(f"{prefix}n{i}\t{node.op}\t{' '.join(args)}\t{[str(o.type) for o in node.outputs]}\n")
        template = getattr(scalar, "c_code_template", None)
        if template is not None:
            inner = [str(x) for x in scalar.fgraph.inputs]
            out.write(f"  C inputs {inner}\n" + "".join("  " + l + "\n" for l in template.splitlines()))
        if type(node.op).__name__ == "Scan":
            fg = node.op.fgraph
            out.write(f"  scan inputs {[str(x) for x in fg.inputs]} outputs {[str(x) for x in fg.outputs]}\n")
            describe_nodes(fg.toposort(), out, prefix + f"n{i}/")
            ids_inner = {}
            for k, o in enumerate(fg.outputs):
                out.write(f"  scan out {k} = {o.owner and fg.toposort().index(o.owner)}\n")


def patch_function(module):
    import pytensor
    from pytensor.compile.sharedvalue import SharedVariable
    original = module.Function.__call__
    names = {}
    calls = {}

    def call(self, *args, **kwargs):
        first = id(self) not in names
        if first:
            names[id(self)] = len(names)
            if names[id(self)] == 0:
                with open(os.path.join(OUT, "argv.txt"), "w") as a:
                    a.write("\n".join(sys.argv) + "\n")
                for arg in sys.argv:
                    path = arg.split("=", 1)[-1]
                    if os.path.isfile(path):
                        shutil.copy(path, os.path.join(OUT, "in." + os.path.basename(path)))
            with open(os.path.join(OUT, "graph.txt"), "a") as g:
                g.write(f"=== function {names[id(self)]}\n")
                pytensor.printing.debugprint(self, file=g, print_type=True, print_destroy_map=True,
                                             print_view_map=True)
            with open(os.path.join(OUT, "nodes.txt"), "a") as g:
                g.write(f"=== function {names[id(self)]} inputs {[str(i.variable) for i in self.maker.inputs]}"
                        f" outputs {[str(o) for o in self.maker.fgraph.outputs]}\n")
                describe_nodes(list(self.vm.nodes), g, "")
        k = names[id(self)]
        shared = [i for i in self.maker.inputs if isinstance(i.variable, SharedVariable)]
        for i in shared:
            s = state(i.variable.get_value(borrow=True))
            if s is not None:
                emit(f"f{k}.rng.{i.variable.name or id(i.variable)}", s)
        out = original(self, *args, **kwargs)
        calls[k] = calls.get(k, 0) + 1
        if calls[k] <= 2:
            for i, (node, thunk) in enumerate(zip(self.vm.nodes, self.vm.thunks)):
                for j, cell in enumerate(getattr(thunk, "outputs", [])):
                    try:
                        emit(f"f{k}.c{calls[k]}.n{i}.{j}", cell[0])
                    except Exception:
                        pass
        for i, a in enumerate(args):
            emit(f"f{k}.in{i}", a)
        for i, o in enumerate(out if isinstance(out, (list, tuple)) else [out]):
            emit(f"f{k}.out{i}", o)
        for n, i in enumerate(shared):
            v = i.variable.get_value(borrow=True)
            if state(v) is None and i.update is not None:
                emit(f"f{k}.upd.{n}.{i.variable.name}", v)
        return out

    module.Function.__call__ = call


def patch_pytensorf(module):
    original = module.reseed_rngs

    def reseed(rngs, seed):
        original(rngs, seed)
        emit("reseed", f"{seed} {len(rngs)} " + " ".join(state(r.get_value(borrow=True)) or "rs" for r in rngs))

    module.reseed_rngs = reseed


TARGETS = {"pytensor.compile.function.types": patch_function, "pymc.pytensorf": patch_pytensorf}


class Hook(importlib.abc.MetaPathFinder):
    def find_spec(self, name, path, target=None):
        if name not in TARGETS:
            return None
        sys.meta_path.remove(self)
        try:
            spec = importlib.util.find_spec(name)
        finally:
            sys.meta_path.insert(0, self)
        run = spec.loader.exec_module

        def exec_module(module):
            run(module)
            TARGETS[name](module)

        spec.loader.exec_module = exec_module
        return spec


if OUT:
    sys.meta_path.insert(0, Hook())

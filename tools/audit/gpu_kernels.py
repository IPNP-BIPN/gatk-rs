#!/usr/bin/env python3
"""Hold every device kernel to the eligibility rule of Milestone GPU (docs/gpu/eligibility.md).

A kernel is a second implementation of a path that is already byte-identical, so it may only do
what cannot change a byte across devices. Per kernel, `gpu_kernels.json` records its class and
the test that holds it to the CPU path; this script refuses:

* a `.metal` source, or a `kernel` in one, that the record does not declare, and a declared one
  that no longer exists;
* in a kernel of class `integer`: any floating-point type (`float`, `half`, `double`, `bfloat`
  and their vectors), because a device's float unit is not the host's;
* in any kernel: fused or fast-math operations (`fma`, `mad`, `fast::`, `precise::`), SIMD-group
  and threadgroup reductions and atomics, whose order is the scheduler's, unless the record lists
  the construct for that file with a reason;
* host code that compiles a kernel source without turning fast math off, or turning it on;
* a record naming a CPU-equality test that does not exist.

Class `float-pinned` admits floating point only under the rule's conditions and must say how it
meets them; `quarantined` kernels never ship and are reported as bio-identical. Neither exists yet.

Usage:
    gpu_kernels.py              check the tree (exit 1 on any violation)
    gpu_kernels.py --self-test  plant each violation in a scratch tree and check it is caught
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import shutil
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[2]
RECORD = "tools/audit/gpu_kernels.json"
CLASSES = {"integer", "float-pinned", "quarantined"}

FLOAT_TYPE = re.compile(r"\b(float|half|double|bfloat)\d*(x\d+)?\b")
FORBIDDEN = {
    r"\bfma\s*\(|\bmad\s*\(": "a fused multiply-add rounds once where the CPU path rounds twice",
    r"\b(fast|precise)::": "a math namespace other than the default picks its own precision",
    r"\bsimd_(sum|product|prefix\w*|min|max|and|or|xor|shuffle\w*|broadcast\w*|ballot)\b":
        "a SIMD-group operation depends on the execution width",
    r"\bquad_\w+\b": "a quad-group operation depends on the execution width",
    r"\batomic\w*\b": "an atomic's order is the scheduler's",
    r"\bthreadgroup\b": "threadgroup memory shares state across threads",
}
KERNEL = re.compile(r"\bkernel\s+void\s+(\w+)")
INCLUDE = re.compile(r'include_str!\s*\(\s*"([^"]+\.metal)"\s*\)')


def strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    return re.sub(r"//[^\n]*", "", text)


def check(root: pathlib.Path) -> list[str]:
    errors = []
    record_path = root / RECORD
    record = json.loads(record_path.read_text()) if record_path.exists() else {"sources": {}}
    sources = record["sources"]
    found = {p.relative_to(root).as_posix() for p in (root / "crates").rglob("*.metal")
             if "target" not in p.parts}
    for path in sorted(found - set(sources)):
        errors.append(f"{path}: a kernel source {RECORD} does not declare")
    for path in sorted(set(sources) - found):
        errors.append(f"{path}: declared in {RECORD} and not in the tree")
    for path in sorted(found & set(sources)):
        entry = sources[path]
        code = strip_comments((root / path).read_text())
        declared = entry.get("kernels", {})
        present = set(KERNEL.findall(code))
        for name in sorted(present - set(declared)):
            errors.append(f"{path}: kernel `{name}` is not declared")
        for name in sorted(set(declared) - present):
            errors.append(f"{path}: kernel `{name}` is declared and not defined")
        classes = set()
        for name, kernel in declared.items():
            if kernel.get("class") not in CLASSES:
                errors.append(f"{path}: kernel `{name}` has no class among {sorted(CLASSES)}")
            classes.add(kernel.get("class"))
            test = kernel.get("test", "")
            if not test or not (root / test).exists():
                errors.append(f"{path}: kernel `{name}` names no existing CPU-equality test")
        allowed = entry.get("allowed", {})
        for pattern, reason in FORBIDDEN.items():
            for match in re.finditer(pattern, code):
                if match.group(0).rstrip(" (") not in allowed:
                    line = code[: match.start()].count("\n") + 1
                    errors.append(f"{path}:{line}: `{match.group(0).strip()}`: {reason}")
        if classes <= {"integer"}:
            for match in FLOAT_TYPE.finditer(code):
                line = code[: match.start()].count("\n") + 1
                errors.append(f"{path}:{line}: `{match.group(0)}` in an integer kernel: a "
                              "device's float unit is not the host's")
    # The host side: whatever compiles a kernel source must turn fast math off, and never on.
    for rust in (root / "crates").rglob("*.rs"):
        if "target" in rust.parts:
            continue
        text = rust.read_text()
        if not INCLUDE.search(text):
            continue
        relative = rust.relative_to(root).as_posix()
        if "set_fast_math_enabled(false)" not in text:
            errors.append(f"{relative}: compiles a kernel without `set_fast_math_enabled(false)`")
        if re.search(r"set_fast_math_enabled\s*\(\s*true", text):
            errors.append(f"{relative}: turns fast math on")
    return errors


def self_test() -> int:
    failures = []
    if check(ROOT):
        failures.append("the clean tree is not clean")
    source = next(iter(json.loads((ROOT / RECORD).read_text())["sources"]))
    host = "crates/gatk-engine/src/gpu/metal.rs"
    plants = [
        (source, "kernel void stray(device uint *o) { o[0] = 1; }\n"),
        (source, "static float f(float a) { return a; }\n"),
        (source, "static ulong f(ulong a) { return fma(a, a, a); }\n"),
        (source, "static uint f(uint a) { return fast::max(a, a); }\n"),
        (source, "static uint f(uint a) { return simd_sum(a); }\n"),
        (source, "static void f(device atomic_uint *a) { }\n"),
        (host, "fn f(o: &CompileOptions) { o.set_fast_math_enabled(true); }\n"),
        ("crates/planted/src/k.metal", "kernel void k(device uint *o) { o[0] = 1; }\n"),
    ]
    for relative, text in plants:
        with tempfile.TemporaryDirectory(prefix="gpu-kernels-self-test-") as tmp:
            scratch = pathlib.Path(tmp) / "tree"
            shutil.copytree(ROOT / "crates", scratch / "crates",
                            ignore=shutil.ignore_patterns("target", "data"))
            (scratch / "tools/audit").mkdir(parents=True)
            shutil.copy(ROOT / RECORD, scratch / RECORD)
            target = scratch / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            with target.open("a") as f:
                f.write(text)
            if not check(scratch):
                failures.append(f"not caught: {text.strip()} in {relative}")
    for failure in failures:
        print(failure)
    print("self-test:", "FAILED" if failures else f"all {len(plants)} plants caught, clean tree clean")
    return 1 if failures else 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    if parser.parse_args().self_test:
        return self_test()
    errors = check(ROOT)
    for error in errors:
        print(error)
    record = json.loads((ROOT / RECORD).read_text())
    kernels = sum(len(e["kernels"]) for e in record["sources"].values())
    print(f"{len(errors)} violations in {len(record['sources'])} sources, {kernels} kernels"
          if errors else f"{kernels} kernels in {len(record['sources'])} sources meet the rule")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())

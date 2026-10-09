#!/usr/bin/env python3
"""Keep a change made for speed from moving a byte. Milestone S, S.2 (#109).

The conformance suites are the first half of the gate: every pull request re-runs all of them,
and a golden that moves fails it, whatever the change was for. They are not enough on their own,
because a suite only sees the orders and the inputs its corpus reaches. Reversing a two-element
sum is invisible to a corpus of pairs, and swapping `jmath`'s `exp` for the host's passes on a
host whose libm happens to agree on the points the corpus visits. `docs/speed/byte-neutrality.md`
records both experiments. This script is the second half: it refuses the constructs that make
arithmetic reordering possible before any suite has to notice them.

Three kinds of rule:

* **Forbidden in Rust source.** Fused multiply-add (`mul_add`), the fast-math intrinsics
  (`fadd_fast` and family, the `algebraic_*` methods that license reassociation), the unstable
  intrinsics feature that exposes them, and explicit SIMD (`std::simd`, `std::arch`), whose
  reductions are a different tree. None occurs in the workspace today.
* **Forbidden in build configuration.** Any `rustflags`, `RUSTFLAGS`, `[profile]` section or
  `CARGO_PROFILE_*` override, in `.cargo/config*`, any `Cargo.toml` or a workflow, and fast-math
  flags in a `build.rs`. The goldens were measured with the default release profile; LLVM folds
  `powf` of a constant differently between profiles, which is a byte that moved once already
  between a debug and a release test run.
* **A ratchet on the host's transcendentals.** `exp`, `ln`, `log10`, `powf`, `powi` and the rest
  are the platform libm, not the JVM's. Some call sites are deliberate (`pow10` measured closer
  to `Math.pow` with `powf` than with `strict_math::pow`, see
  `docs/numeric-functions-a-ported-call-site-reaches.md`), so they are counted per file in
  `arithmetic.json` rather than forbidden, and a file may not gain one without the record being
  rewritten in the same change, where a reviewer sees it. `sqrt` is correctly rounded by IEEE 754
  and is not counted.

Usage:
    arithmetic.py              check the tree (exit 1 on any violation)
    arithmetic.py --record     rewrite the ratchet from the current tree (review the diff)
    arithmetic.py --self-test  plant every kind of violation in a scratch tree and check each is
                               caught, and that the clean tree passes: the guard verified in both
                               directions, as the provenance guard was (#82)
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
RATCHET = pathlib.Path(__file__).resolve().parent / "arithmetic.json"

FORBIDDEN_SOURCE = {
    r"\bmul_add\s*\(": "fused multiply-add rounds once where the reference rounds twice",
    r"\b(fadd|fsub|fmul|fdiv|frem)_fast\b": "a fast-math intrinsic licenses reassociation",
    r"\b(fadd|fsub|fmul|fdiv|frem)_algebraic\b|\balgebraic_(add|sub|mul|div|rem)\b":
        "an algebraic float operation licenses reassociation",
    r"feature\s*\(\s*core_intrinsics": "the intrinsics feature exposes the fast-math operations",
    r"\b(std|core)::simd\b|portable_simd": "a SIMD reduction sums in a different tree",
    r"\b(std|core)::arch\b": "explicit SIMD intrinsics reorder what they vectorise",
}

FORBIDDEN_CONFIG = {
    r"rustflags|RUSTFLAGS": "compiler flags are not part of the measured build",
    r"^\s*\[profile": "the goldens were measured with the default release profile",
    r"CARGO_PROFILE_": "the goldens were measured with the default release profile",
    r"target-cpu|target-feature|llvm-args": "a code-generation flag outside the measured build",
}

FORBIDDEN_BUILD_SCRIPT = {
    r"-ffast-math|-Ofast|-ffp-contract=fast|-funsafe-math-optimizations|-fassociative-math":
        "a C fast-math flag reorders arithmetic in what the build script compiles",
}

TRANSCENDENTAL = re.compile(
    r"\.(exp|exp2|exp_m1|ln|ln_1p|log|log2|log10|powf|powi|cbrt|sin|cos|tan|asin|acos|atan|atan2"
    r"|sinh|cosh|tanh|asinh|acosh|atanh|hypot)\("
)


def rust_sources(root: pathlib.Path) -> list[pathlib.Path]:
    return sorted(p for p in (root / "crates").rglob("*.rs") if "target" not in p.parts)


def config_files(root: pathlib.Path) -> list[pathlib.Path]:
    found = list((root / ".cargo").glob("config*")) if (root / ".cargo").is_dir() else []
    found += [p for p in root.rglob(".cargo/config*") if "target" not in p.parts]
    found += [p for p in root.rglob("Cargo.toml") if "target" not in p.parts]
    workflows = root / ".github" / "workflows"
    if workflows.is_dir():
        found += sorted(workflows.glob("*.yml")) + sorted(workflows.glob("*.yaml"))
    template = root / "tools" / "conformance" / "ci.template.yml"
    if template.exists():
        found.append(template)
    return sorted(set(found))


def scan(path: pathlib.Path, rules: dict[str, str], root: pathlib.Path) -> list[str]:
    """Every line of `path` matching a rule. Every line is read: none is skipped as a comment."""
    found = []
    text = path.read_text(encoding="utf-8", errors="replace")
    for number, line in enumerate(text.splitlines(), 1):
        for pattern, reason in rules.items():
            if re.search(pattern, line):
                found.append(f"{path.relative_to(root)}:{number}: {line.strip()}\n    -> {reason}")
    return found


def transcendentals(root: pathlib.Path) -> dict[str, int]:
    counts = {}
    for path in rust_sources(root):
        n = len(TRANSCENDENTAL.findall(path.read_text(encoding="utf-8", errors="replace")))
        if n:
            counts[str(path.relative_to(root))] = n
    return counts


def check(root: pathlib.Path, ratchet: dict[str, int]) -> list[str]:
    violations = []
    for path in rust_sources(root):
        violations += scan(path, FORBIDDEN_SOURCE, root)
        if path.name == "build.rs":
            violations += scan(path, FORBIDDEN_BUILD_SCRIPT, root)
    for path in config_files(root):
        violations += scan(path, FORBIDDEN_CONFIG, root)
    for path, count in sorted(transcendentals(root).items()):
        allowed = ratchet.get(path, 0)
        if count > allowed:
            violations.append(
                f"{path}: {count} host transcendental call(s), {allowed} recorded\n"
                "    -> the platform libm is not the JVM's; use jmath, or record the call site "
                "with --record and say why in the change"
            )
    return violations


def self_test() -> int:
    """Plant one violation of each kind in a copy of the tree, and check each is reported."""
    ratchet = json.loads(RATCHET.read_text())
    failures = []
    if check(ROOT, ratchet):
        failures.append("the clean tree is not clean")
    plants = [
        ("crates/planted/src/lib.rs", "pub fn f(a: f64, b: f64, c: f64) -> f64 { a.mul_add(b, c) }\n"),
        ("crates/planted/src/lib.rs", "pub fn f(a: f64, b: f64) -> f64 { unsafe { fadd_fast(a, b) } }\n"),
        ("crates/planted/src/lib.rs", "pub fn f(a: f64, b: f64) -> f64 { a.algebraic_add(b) }\n"),
        ("crates/planted/src/lib.rs", "#![feature(core_intrinsics)]\n"),
        ("crates/planted/src/lib.rs", "use std::simd::f64x4;\n"),
        ("crates/planted/src/lib.rs", "use std::arch::x86_64::_mm256_add_pd;\n"),
        ("crates/planted/src/lib.rs", "pub fn f(x: f64) -> f64 { x.exp() }\n"),
        ("crates/planted/build.rs", 'fn main() { cc::Build::new().flag("-ffast-math"); }\n'),
        (".cargo/config.toml", '[build]\nrustflags = ["-C", "target-cpu=native"]\n'),
        ("Cargo.toml", "\n[profile.release]\nlto = true\n"),
        (".github/workflows/fast.yml", "env:\n  RUSTFLAGS: -Ctarget-cpu=native\n"),
    ]
    with tempfile.TemporaryDirectory(prefix="arithmetic-self-test-") as tmp:
        for relative, text in plants:
            scratch = pathlib.Path(tmp) / "tree"
            if scratch.exists():
                shutil.rmtree(scratch)
            scratch.mkdir()
            shutil.copytree(ROOT / "crates", scratch / "crates", ignore=shutil.ignore_patterns("target"))
            shutil.copy(ROOT / "Cargo.toml", scratch / "Cargo.toml")
            target = scratch / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            with open(target, "a", encoding="utf-8") as fh:
                fh.write(text)
            if not check(scratch, ratchet):
                failures.append(f"not caught: {relative}: {text.strip()}")
    for failure in failures:
        print(failure)
    print(f"self-test: {len(plants)} planted violations, {len(failures)} problem(s)")
    return 1 if failures else 0


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--record", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args(argv)
    if args.self_test:
        return self_test()
    if args.record:
        counts = transcendentals(ROOT)
        RATCHET.write_text(json.dumps(counts, indent=1, sort_keys=True) + "\n")
        print(f"recorded {sum(counts.values())} host transcendental call(s) in {len(counts)} file(s)")
        return 0
    violations = check(ROOT, json.loads(RATCHET.read_text()))
    for violation in violations:
        print(violation)
    if violations:
        print(f"\n{len(violations)} arithmetic violation(s). See docs/speed/byte-neutrality.md.")
        return 1
    print("no construct that reorders arithmetic, and no new host transcendental")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

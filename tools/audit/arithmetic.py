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
  intrinsics feature that exposes them, and portable SIMD (`std::simd`), whose reductions are a
  different tree.
* **Classified, not forbidden** (the S.2 amendment, `docs/speed/byte-neutrality.md`):
  - `std::arch`, `#[target_feature]` and runtime feature detection, only in a file that
    `simd_kernels.json` declares with what it does and the test that holds every one of its paths
    to the scalar one. Inside such a file the fused, horizontal and approximate intrinsics stay
    refused, and so does enabling `fma`: a declared kernel may only do, lane by lane, the scalar
    operations in the scalar order.
  - `#[global_allocator]`, only in `ALLOCATOR_FILES`: an allocator moves addresses, and an address
    reaches no output the suites check (`docs/pointers-that-reach-the-output.md`).
  - one `[profile]` table, `ALLOWED_PROFILES`, in the workspace's `Cargo.toml`, key for key.
* **Forbidden in build configuration.** Any `rustflags`, `RUSTFLAGS`, `target-cpu`, other
  `[profile]` section or key, or `CARGO_PROFILE_*` override, in `.cargo/config*`, any `Cargo.toml`
  or a workflow, and fast-math flags in a `build.rs`. LLVM folds `powf` of a constant differently
  between profiles, which is a byte that moved once already between a debug and a release test
  run, so a profile is admitted only when every suite has run on it.
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
import tomllib
import pathlib
import re
import shutil
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[2]
RATCHET = pathlib.Path(__file__).resolve().parent / "arithmetic.json"
SIMD_RECORD = "tools/audit/simd_kernels.json"

FORBIDDEN_SOURCE = {
    r"\bmul_add\s*\(": "fused multiply-add rounds once where the reference rounds twice",
    r"\b(fadd|fsub|fmul|fdiv|frem)_fast\b": "a fast-math intrinsic licenses reassociation",
    r"\b(fadd|fsub|fmul|fdiv|frem)_algebraic\b|\balgebraic_(add|sub|mul|div|rem)\b":
        "an algebraic float operation licenses reassociation",
    r"feature\s*\(\s*core_intrinsics": "the intrinsics feature exposes the fast-math operations",
    r"\b(std|core)::simd\b|portable_simd": "a SIMD reduction sums in a different tree",
}

# Allowed only in a file `simd_kernels.json` declares.
SIMD_SOURCE = {
    r"\b(std|core)::arch\b|\btarget_feature\b|\bis_\w+_feature_detected\b":
        "explicit SIMD is admitted only in a file simd_kernels.json declares with its test",
}

# Refused even in a declared file: each fuses, reduces across lanes or approximates.
SIMD_IN_KERNEL = {
    r"fmadd|fmsub|fnmadd|fnmsub|\bvfma|\bvfms|\bfma\b": "a fused multiply-add rounds once",
    r"hadd|hsub|_reduce_|\bvaddv|\bvpadd|\bvmaxv|\bvminv|_dp_p":
        "a horizontal operation sums across lanes, in a tree of its own",
    r"_rcp|_rsqrt|\bvrecpe|\bvrsqrte": "an approximate reciprocal is not the scalar division",
}

# `#[global_allocator]`, by file: the shipped binary only, so no library forces it on a caller.
ALLOCATOR_FILES = {"crates/gatk-cli/src/main.rs"}
ALLOCATOR = {r"global_allocator": "a global allocator is admitted only in ALLOCATOR_FILES"}

# The one profile table admitted, in the workspace's own Cargo.toml, exactly. `dist` is the
# shipped binary's: whole-program LTO and one codegen unit, neither of which licenses LLVM to
# reorder or fuse floating point (docs/speed/byte-neutrality.md says how that was shown).
ALLOWED_PROFILES = {"dist": {"inherits": "release", "lto": "fat", "codegen-units": 1}}

# Where the REFERENCE itself rounds once, a fused multiply-add is the port, not a speed change:
# `Java11BLAS.dgemvN`, which Breeze's dense matrix-vector product calls, is a chain of `Math.fma`.
# One line each, by path and text, so a second call site is still a violation.
FUSED_IN_REFERENCE = {
    ("crates/gatk-tools/src/create_read_count_panel_of_normals.rs", "a.mul_add(b, c)"),
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


def scan(
    path: pathlib.Path, rules: dict[str, str], root: pathlib.Path, skip: set[str] = frozenset()
) -> list[str]:
    """Every line of `path` matching a rule. Every line is read: none is skipped as a comment.

    `skip` holds whole stripped lines the caller has already checked another way, such as the
    header of the one admitted profile table.
    """
    found = []
    text = path.read_text(encoding="utf-8", errors="replace")
    for number, line in enumerate(text.splitlines(), 1):
        if line.strip() in skip:
            continue
        for pattern, reason in rules.items():
            if re.search(pattern, line):
                if (str(path.relative_to(root)), line.strip()) in FUSED_IN_REFERENCE:
                    continue
                found.append(f"{path.relative_to(root)}:{number}: {line.strip()}\n    -> {reason}")
    return found


def transcendentals(root: pathlib.Path) -> dict[str, int]:
    counts = {}
    for path in rust_sources(root):
        n = len(TRANSCENDENTAL.findall(path.read_text(encoding="utf-8", errors="replace")))
        if n:
            counts[str(path.relative_to(root))] = n
    return counts


def simd_files(root: pathlib.Path) -> tuple[dict[str, dict], list[str]]:
    """The declared SIMD files, and what is wrong with the record itself."""
    record = root / SIMD_RECORD
    declared = json.loads(record.read_text())["sources"] if record.exists() else {}
    errors = []
    for source, entry in sorted(declared.items()):
        if not (root / source).exists():
            errors.append(f"{SIMD_RECORD}: {source} is declared but does not exist")
        for key in ("what", "test"):
            if not entry.get(key):
                errors.append(f"{SIMD_RECORD}: {source} has no `{key}`")
        test = entry.get("test", "")
        if test and not (root / test).exists():
            errors.append(f"{SIMD_RECORD}: {source} names the test {test}, which does not exist")
    return declared, errors


def profile_violations(root: pathlib.Path) -> tuple[list[str], set[str]]:
    """The workspace manifest's profile tables against `ALLOWED_PROFILES`, key for key.

    Returns the violations and the table headers that passed, which the line scan then skips.
    """
    manifest = root / "Cargo.toml"
    if not manifest.exists():
        return [], set()
    profiles = tomllib.loads(manifest.read_text(encoding="utf-8")).get("profile", {})
    errors = []
    for name, table in sorted(profiles.items()):
        if ALLOWED_PROFILES.get(name) != table:
            errors.append(
                f"Cargo.toml: [profile.{name}] = {table}\n"
                f"    -> only {ALLOWED_PROFILES} is admitted, exactly; see "
                "docs/speed/byte-neutrality.md"
            )
    passed = {f"[profile.{name}]" for name in profiles if ALLOWED_PROFILES.get(name) == profiles[name]}
    return errors, passed


def check(root: pathlib.Path, ratchet: dict[str, int]) -> list[str]:
    declared, violations = simd_files(root)
    for path in rust_sources(root):
        relative = str(path.relative_to(root))
        violations += scan(path, FORBIDDEN_SOURCE, root)
        violations += scan(path, SIMD_IN_KERNEL if relative in declared else SIMD_SOURCE, root)
        if relative not in ALLOCATOR_FILES:
            violations += scan(path, ALLOCATOR, root)
        if path.name == "build.rs":
            violations += scan(path, FORBIDDEN_BUILD_SCRIPT, root)
    profile_errors, admitted = profile_violations(root)
    violations += profile_errors
    for path in config_files(root):
        skip = admitted if path == root / "Cargo.toml" else set()
        violations += scan(path, FORBIDDEN_CONFIG, root, skip)
    for path, count in sorted(transcendentals(root).items()):
        allowed = ratchet.get(path, 0)
        if count > allowed:
            violations.append(
                f"{path}: {count} host transcendental call(s), {allowed} recorded\n"
                "    -> the platform libm is not the JVM's; use jmath, or record the call site "
                "with --record and say why in the change"
            )
    return violations


PLANTED_KERNEL = "crates/planted/src/kernel.rs"
PLANTED_TEST = "crates/planted/tests/kernel.rs"


def build_scratch(scratch: pathlib.Path) -> None:
    """A copy of the tree, plus one declared SIMD kernel of its own and that kernel's test.

    The planted kernel is what the in-kernel plants go into, so the self-test does not depend on
    which real files the record declares.
    """
    scratch.mkdir()
    shutil.copytree(ROOT / "crates", scratch / "crates", ignore=shutil.ignore_patterns("target"))
    shutil.copy(ROOT / "Cargo.toml", scratch / "Cargo.toml")
    record = json.loads((ROOT / SIMD_RECORD).read_text())
    record["sources"][PLANTED_KERNEL] = {"what": "a lane-wise add", "test": PLANTED_TEST}
    (scratch / SIMD_RECORD).parent.mkdir(parents=True)
    (scratch / SIMD_RECORD).write_text(json.dumps(record))
    for path, text in (
        (PLANTED_KERNEL, "use std::arch::x86_64::_mm256_add_pd;\n"),
        (PLANTED_TEST, "#[test]\nfn lanes() {}\n"),
    ):
        (scratch / path).parent.mkdir(parents=True, exist_ok=True)
        (scratch / path).write_text(text)


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
        ("crates/planted/src/lib.rs", '#[target_feature(enable = "avx2")]\nunsafe fn f() {}\n'),
        ("crates/planted/src/lib.rs", 'pub fn f() -> bool { is_x86_feature_detected!("avx2") }\n'),
        ("crates/planted/src/lib.rs", "#[global_allocator]\nstatic A: std::alloc::System = std::alloc::System;\n"),
        # In a declared SIMD file: a fused, a horizontal and an approximate intrinsic, and `fma`.
        ("{simd}", "fn f(a: __m256d) -> __m256d { _mm256_fmadd_pd(a, a, a) }\n"),
        ("{simd}", "fn f(a: __m256d) -> __m256d { _mm256_hadd_pd(a, a) }\n"),
        ("{simd}", "fn f(a: float64x2_t) -> f64 { vaddvq_f64(a) }\n"),
        ("{simd}", "fn f(a: __m256) -> __m256 { _mm256_rcp_ps(a) }\n"),
        ("{simd}", '#[target_feature(enable = "avx2,fma")]\nunsafe fn f() {}\n'),
        # A declared file whose test is gone, and a declared file that is gone.
        ("{simd-test}", None),
        ("{simd}", None),
        ("crates/planted/src/lib.rs", "pub fn f(x: f64) -> f64 { x.exp() }\n"),
        ("crates/planted/build.rs", 'fn main() { cc::Build::new().flag("-ffast-math"); }\n'),
        (".cargo/config.toml", '[build]\nrustflags = ["-C", "target-cpu=native"]\n'),
        ("Cargo.toml", "\n[profile.release]\nlto = true\n"),
        ("Cargo.toml", '\n[profile.dist.package."*"]\nopt-level = 2\n'),
        ("crates/gatk-cli/Cargo.toml", '\n[profile.dist]\ninherits = "release"\nlto = "fat"\n'),
        (".github/workflows/fast.yml", "env:\n  CARGO_PROFILE_DIST_LTO: thin\n"),
        (".github/workflows/fast.yml", "env:\n  RUSTFLAGS: -Ctarget-cpu=native\n"),
    ]
    with tempfile.TemporaryDirectory(prefix="arithmetic-self-test-") as tmp:
        scratch = pathlib.Path(tmp) / "tree"
        for index, (relative, text) in enumerate([(None, None)] + plants):
            if scratch.exists():
                shutil.rmtree(scratch)
            build_scratch(scratch)
            if relative is None:
                # The scratch tree before any plant, with its declared kernel: it must be clean, or
                # every plant below would count as caught for the wrong reason.
                if check(scratch, ratchet):
                    failures.append("the scratch tree is not clean before planting")
                continue
            target = scratch / relative.replace("{simd}", PLANTED_KERNEL).replace(
                "{simd-test}", PLANTED_TEST
            )
            if text is None:
                target.unlink()
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                with open(target, "a", encoding="utf-8") as fh:
                    fh.write(text)
            if not check(scratch, ratchet):
                failures.append(f"not caught: {relative}: {(text or '(deleted)').strip()}")
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

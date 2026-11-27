#!/usr/bin/env python3
"""Prove that a branch changed comments and nothing else.

A commenting pull request claims to leave behaviour alone. The test suite is one witness to that
claim; this script is a stricter one, because it does not depend on what the tests happen to
exercise. It compares every file the branch touched against the base revision with the comments
taken out, and fails if anything else differs.

- A Rust file is compared line by line after dropping blank lines and every line whose first
  non-whitespace characters are ``//`` (so ``//``, ``///`` and ``//!``). A comment trailing code on
  the same line is part of that code line, so editing one is reported: the check would otherwise
  have to parse string literals to know where the code ends.
- A Python file is compared by its abstract syntax tree with every docstring removed, so a
  docstring may be added or rewritten and PEP 8 reformatting is invisible, but a renamed variable
  or a changed constant is not.
- Any other file is ignored, since a commenting branch also edits lists and documentation.

It also refuses an em-dash anywhere in an added line, which the house style forbids.

Usage:
    comments_only.py [BASE]    compare the working tree against BASE (default: origin/main)
"""

from __future__ import annotations

import ast
import subprocess
import sys

EM_DASH = "—"


def git(*args: str) -> str:
    """Run git in the current repository and return its standard output."""
    return subprocess.run(
        ["git", *args], check=True, capture_output=True, text=True
    ).stdout


def base_text(base: str, path: str) -> str | None:
    """The file as it stands at ``base``, or None when the branch created it."""
    try:
        return git("show", f"{base}:{path}")
    except subprocess.CalledProcessError:
        return None


def rust_code(text: str) -> list[str]:
    """The non-blank lines of a Rust source that are not whole-line comments."""
    return [
        line.rstrip()
        for line in text.splitlines()
        if line.strip() and not line.strip().startswith("//")
    ]


def without_docstrings(tree: ast.AST) -> ast.AST:
    """Remove the docstring, if any, from the module and from every class and function."""
    kinds = (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
    for node in ast.walk(tree):
        if not isinstance(node, kinds) or not node.body:
            continue
        first = node.body[0]
        if (
            isinstance(first, ast.Expr)
            and isinstance(first.value, ast.Constant)
            and isinstance(first.value.value, str)
        ):
            # A body may not be empty, so a function that was only a docstring keeps a `pass`.
            node.body = node.body[1:] or [ast.Pass()]
    return tree


def python_code(text: str) -> str:
    """A canonical dump of a Python source's syntax tree without its docstrings."""
    return ast.dump(without_docstrings(ast.parse(text)))


def main() -> int:
    """Compare every changed Rust and Python file against the base, and report the ones that
    differ in more than their comments."""
    base = sys.argv[1] if len(sys.argv) > 1 else "origin/main"
    changed = git("diff", "--name-only", "--diff-filter=AM", base).split()
    failed = False
    for path in changed:
        if not path.endswith((".rs", ".py")):
            continue
        with open(path, encoding="utf-8") as handle:
            new = handle.read()
        old = base_text(base, path)
        if old is None:
            # A new file is new code, which a commenting branch does not add.
            print(f"{path}: new file")
            failed = True
            continue
        if path.endswith(".rs"):
            same = rust_code(old) == rust_code(new)
        else:
            same = python_code(old) == python_code(new)
        if not same:
            print(f"{path}: differs from {base} in more than its comments")
            failed = True
    sources = [path for path in changed if path.endswith((".rs", ".py"))]
    diff = git("diff", base, "--", *sources) if sources else ""
    for line in diff.splitlines():
        if line.startswith("+") and not line.startswith("+++") and EM_DASH in line:
            print(f"em-dash in an added line: {line[1:].strip()}")
            failed = True
    if failed:
        return 1
    print(f"{len(changed)} changed files, comments only")
    return 0


if __name__ == "__main__":
    sys.exit(main())

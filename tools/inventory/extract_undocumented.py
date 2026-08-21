#!/usr/bin/env python3
"""Recover the tools that `gatkdoc` does not document.

`gatk --list` reports 311 runnable tools; the generated `gatkdoc/*.json` covers only 268 of
them. The 43-tool gap is not noise: it includes real tools (`CheckDuplicateMarking`,
`CollectF1R2Counts`, `AddOriginalAlignmentTags`, ...) that simply carry no doc annotation.

Deriving the inventory from `gatkdoc` alone therefore produces a list that is quietly 14%
short, with nothing to indicate it. The CLI is the ground truth for *which tools exist*;
`gatkdoc` is only the richer source for *what arguments they take*. This script closes the gap
by parsing `<Tool> --help` for whatever the docs missed.

Usage:
  ./extract_undocumented.py <gatk-local.jar> <inventory.json> [-o merged.json]
"""

import argparse
import collections
import json
import pathlib
import re
import subprocess
import sys

ANSI = re.compile(r"\x1b\[[0-9;]*m")
# --NAME,-S <Type>   or   --NAME <Type>
ARG_LINE = re.compile(r"^--(?P<name>[A-Za-z0-9_.\-]+)(?:,(?P<short>-[A-Za-z0-9_\-]+))?\s+<(?P<type>[^>]+)>\s*(?P<rest>.*)$")
DEFAULT_VALUE = re.compile(r"Default value:\s*(?P<v>.+?)\.(?:\s|$)")
POSSIBLE_VALUES = re.compile(r"Possible values:\s*\{(?P<v>[^}]*)\}")
LIST_ENTRY = re.compile(r"^ {4}([A-Za-z][A-Za-z0-9]*)\s")


def run(jar, jvm_args, tool_args):
    p = subprocess.run(
        ["java", *jvm_args, "-jar", str(jar), *tool_args],
        capture_output=True, text=True,
    )
    return ANSI.sub("", p.stdout + p.stderr)


def cli_tool_names(jar):
    """The authoritative set of runnable tools, as the CLI itself reports them."""
    out = run(jar, [], ["--list"])
    names = set()
    for line in out.splitlines():
        m = LIST_ENTRY.match(line)
        if m:
            names.add(m.group(1))
    return names


def parse_help(text):
    """Extract arguments from a `--help` rendering.

    Descriptions wrap across indented continuation lines, and `Default value:` /
    `Possible values:` can land on any of them, so each argument's lines are joined before
    the fields are pulled out.
    """
    args, current, buf = [], None, []

    def flush():
        if current is None:
            return
        blob = " ".join(buf)
        dv = DEFAULT_VALUE.search(blob)
        pv = POSSIBLE_VALUES.search(blob)
        default = dv.group("v").strip() if dv else None
        if default in ("null", "none", ""):
            default = None
        members = None
        if pv:
            members = [x.strip() for x in pv.group("v").split(",") if x.strip()]
        current["default"] = default
        current["enum_members"] = members
        args.append(current)

    section_required = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("Required Arguments"):
            section_required = True
            continue
        if stripped.startswith("Optional Arguments") or stripped.startswith("Advanced Arguments"):
            section_required = False
            continue
        m = ARG_LINE.match(stripped)
        if m:
            flush()
            buf = [m.group("rest")]
            current = {
                "name": f"--{m.group('name')}",
                "short": m.group("short"),
                "type": m.group("type").strip(),
                "required": section_required,
                "kind": "required" if section_required else "optional",
                "min": None,
                "max": None,
                "deprecated": False,
            }
        elif current is not None and line.startswith(" ") and stripped:
            buf.append(stripped)
    flush()
    return args


def classes_in_jar(jar):
    """Every class in the jar, as `simple name -> [package.Class]`.

    An undocumented tool's origin cannot be guessed from its name: twelve of the tools
    `gatk --list` reports are Picard's, bundled in the same jar, and calling them GATK-origin
    puts them under the wrong milestone and the wrong archetype for good. The jar knows: a class
    under `picard/` is Picard's and one under `org/broadinstitute/hellbender/` is GATK's.
    """
    import zipfile

    found = collections.defaultdict(list)
    with zipfile.ZipFile(jar) as z:
        for entry in z.namelist():
            if not entry.endswith(".class") or "$" in entry:
                continue
            fqn = entry[: -len(".class")].replace("/", ".")
            found[fqn.rsplit(".", 1)[-1]].append(fqn)
    return found


def resolve_class(by_simple_name, name):
    """The one class a tool's name resolves to, or `(None, None, "gatk")` when it is ambiguous.

    Ambiguity is possible in principle (two packages, one simple name) and is reported rather
    than guessed at: a wrong package is worse than an absent one.
    """
    candidates = by_simple_name.get(name, [])
    tools = [c for c in candidates if c.startswith(("picard.", "org.broadinstitute.hellbender."))]
    if len(tools) != 1:
        return None, None, "gatk"
    fqn = tools[0]
    package, cls = fqn.rsplit(".", 1)
    return package, fqn, ("picard" if fqn.startswith("picard.") else "gatk")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("jar", type=pathlib.Path)
    ap.add_argument("inventory", type=pathlib.Path)
    ap.add_argument("-o", "--out", type=pathlib.Path)
    a = ap.parse_args()

    inv = json.loads(a.inventory.read_text())
    documented = {t["name"]: t for t in inv["tools"]}
    cli = cli_tool_names(a.jar)

    if not cli:
        sys.exit("`--list` returned no tools; refusing to produce a partial inventory")

    phantom = sorted(set(documented) - cli)
    missing = sorted(cli - set(documented))
    print(f"cli={len(cli)} documented={len(documented)} undocumented={len(missing)} phantom={len(phantom)}")
    if phantom:
        print(f"  WARNING: documented but not runnable: {phantom[:5]}")

    by_simple_name = classes_in_jar(a.jar)

    recovered = []
    unresolved = []
    for i, name in enumerate(missing, 1):
        text = run(a.jar, [], [name, "--help"])
        parsed = parse_help(text)
        if not parsed:
            print(f"  [{i}/{len(missing)}] {name}: NO ARGUMENTS PARSED")
        package, cls, origin = resolve_class(by_simple_name, name)
        if cls is None:
            unresolved.append(name)
        recovered.append({
            "name": name,
            "display_name": name,
            "origin": origin,
            "package": package,
            "class": cls,
            "group": "Undocumented",
            "archetype": "unclassified",
            "spark": "Spark" in name,
            "beta": False,
            "experimental": False,
            "deprecated": False,
            "summary": None,
            "branch": f"tool/{origin}-{re.sub(r'[^a-z0-9]+', '', name.lower())}",
            "kind": "tool",
            "documented": False,
            "arguments": parsed,
        })
        print(f"  [{i}/{len(missing)}] {name}: {len(parsed)} arguments")

    if unresolved:
        # Not fatal: the tool still exists and its arguments were read. But a tool whose class
        # the jar cannot name is one whose origin is a guess, so it is printed rather than
        # buried.
        print(f"  WARNING: {len(unresolved)} undocumented tools resolved to no single class: "
              f"{unresolved[:5]}")
    by_origin = collections.Counter(t["origin"] for t in recovered)
    print(f"  undocumented by origin: {dict(by_origin)}")

    for t in inv["tools"]:
        t["documented"] = True
    inv["tools"].extend(recovered)
    inv["counts"]["tools"] = len(inv["tools"])
    inv["counts"]["undocumented_tools"] = len(recovered)
    inv["counts"]["arguments"] = sum(len(t["arguments"]) for t in inv["tools"])
    inv["counts"]["cli_reported_tools"] = len(cli)

    out = a.out or a.inventory
    out.write_text(json.dumps(inv, indent=2))
    print(f"merged: tools={inv['counts']['tools']} arguments={inv['counts']['arguments']} -> {out}")


if __name__ == "__main__":
    main()

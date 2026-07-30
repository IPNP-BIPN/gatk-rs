#!/usr/bin/env python3
"""Generate the per-tool status dashboard from the inventory and the ports' manifests.

The status of the programme lived in three hand-written README tables that had already drifted:
picard-rs listed five byte-identical tools while its tree held suites for forty-four, and nothing
said which of the 311 tools had been started at all. A table that is maintained by hand is a table
that flatters, so this one is generated from two machine-readable sources and cannot say anything
neither of them says:

* `tools/inventory/generated/inventory.json`, the 311-tool ground truth, generated from the pinned
  reference's own CLI;
* each port's `tools/conformance/manifest.json`, which declares every suite, the tools it covers,
  and whether the oracle re-derives its goldens (`oracle-backed`) or nothing does (`unchecked`).

A tool's state is therefore one of:

| state | meaning |
|---|---|
| `oracle-backed` | at least one suite covers it and CI re-derives its goldens every run |
| `unchecked` | a suite covers it, but no CI step has ever re-derived its goldens |
| `not started` | no suite mentions it |

None of these is "byte-identical over its whole argument surface", which is the programme's actual
bar. Coverage per tool (the t-wise arrays, the fuzzer's branch threshold) is not yet measured, so
the dashboard says so in the column rather than leaving it blank and letting a reader assume.

    python3 tools/dashboard/generate.py                 # writes docs/STATUS.md
    python3 tools/dashboard/generate.py --check         # fails if the committed file is stale
"""

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
INVENTORY = REPO / "tools" / "inventory" / "generated" / "inventory.json"
STATUS = REPO / "docs" / "STATUS.md"
VENDORED = Path(__file__).resolve().parent / "ports"

# The sibling ports, as they sit on a maintainer's machine. CI has only this repository checked
# out, so what the dashboard reads is the *vendored* copy under tools/dashboard/ports, refreshed
# from these paths with --refresh. That keeps the generated file reproducible from committed data
# alone: a dashboard that could only be regenerated on one machine is a hand-written table with
# extra steps.
PORTS = [
    ("picard-rs", REPO.parent / "Picard" / "tools" / "conformance" / "manifest.json"),
    ("htsjdk-rs", REPO.parent / "htsjdk" / "tools" / "conformance" / "manifest.json"),
]

# Each port's argument-coverage measurement, written by its own CI from the covering arrays. A
# port with no such file has not run its arrays, and every tool it covers reports `not measured`
# rather than nothing, so a reader cannot mistake silence for coverage.
MEASURED = {
    "picard-rs": REPO.parent / "Picard" / "tools" / "coverage" / "measured.json",
    "htsjdk-rs": REPO.parent / "htsjdk" / "tools" / "coverage" / "measured.json",
}


def summarize(manifest):
    """The subset of a port's manifest the dashboard needs, and nothing else."""
    return {
        "suites": [
            {
                "id": suite["id"],
                "status": suite["status"],
                "tools": suite.get("tools", []),
                "cases": len(suite["cases"]),
            }
            for suite in manifest["suites"]
        ]
    }


def refresh():
    """Re-vendor each port's suite list from the sibling working copy."""
    VENDORED.mkdir(exist_ok=True)
    written, missing = [], []
    for name, path in PORTS:
        if not path.exists():
            missing.append((name, path))
            continue
        with open(path) as fh:
            summary = summarize(json.load(fh))
        summary["source"] = str(path)
        measured = MEASURED.get(name)
        if measured and measured.exists():
            with open(measured) as fh:
                summary["coverage"] = json.load(fh).get("tools", {})
        out = VENDORED / f"{name}.json"
        out.write_text(json.dumps(summary, indent=2) + "\n")
        written.append((name, len(summary["suites"])))
    return written, missing


def load_manifests():
    """Read the vendored summaries. A port with no vendored file is reported, not skipped."""
    found, missing = {}, []
    for name, path in PORTS:
        vendored = VENDORED / f"{name}.json"
        if vendored.exists():
            with open(vendored) as fh:
                found[name] = json.load(fh)
        else:
            missing.append((name, vendored))
    return found, missing


# Strongest claim first. A tool covered by both an oracle-backed suite and a weaker one is
# reported at the stronger state, because that suite really does re-derive its goldens; the weaker
# suites still show in its `suites` column, so nothing is hidden.
RANK = {"oracle-backed": 3, "golden-pending": 2, "unchecked": 1, "not started": 0}


def coverage_cell(entry, manifests):
    """The argument-coverage column for one tool.

    `not measured` is the honest answer for a tool whose arrays have never been run, and it is
    also the answer for a tool that has an array but no port binary to run it against: an array
    run against the reference alone measures the reference, not the port.

    Where there is a measurement, the cell carries the strength, the fraction, and the number of
    distinct outputs the accepted rows produced. The last one is not decoration: an array whose
    accepted rows all produce the same output covers its argument pairs without observing them, so
    a high fraction over one distinct output says nothing about the port.
    """
    measured = manifests.get(entry["port"], {}).get("coverage", {}).get(entry["tool"])
    if not measured:
        return "not measured"
    cell = (
        f"t={measured['t']}, {measured['matched']}/{measured['rows']} rows "
        f"({measured['share']:.0%})"
    )
    if measured["distinct_outputs"] <= 1:
        cell += ", **1 distinct output**"
    return cell


def tool_states(manifests):
    """Map tool name -> {'state', 'suites', 'port', 'cases'}."""
    states = {}
    for port, manifest in manifests.items():
        for suite in manifest["suites"]:
            for tool in suite.get("tools", []):
                entry = states.setdefault(
                    tool,
                    {"state": "unchecked", "suites": [], "port": port, "cases": 0, "tool": tool},
                )
                entry["suites"].append(suite["id"])
                entry["cases"] += suite["cases"]
                if RANK.get(suite["status"], 0) > RANK[entry["state"]]:
                    entry["state"] = suite["status"]
    return states


def render(inventory, manifests, missing):
    tools = inventory["tools"]
    states = tool_states(manifests)

    counted = {"oracle-backed": 0, "golden-pending": 0, "unchecked": 0, "not started": 0}
    rows = []
    for tool in sorted(tools, key=lambda t: (t["origin"], t["name"])):
        entry = states.get(tool["name"])
        state = entry["state"] if entry else "not started"
        counted[state] += 1
        rows.append((tool, state, entry))

    lines = []
    lines.append("# Status")
    lines.append("")
    lines.append(
        "Generated by `tools/dashboard/generate.py` from the tool inventory and the ports' "
        "conformance manifests. Do not edit by hand."
    )
    lines.append("")
    reference = inventory["reference"]
    lines.append(
        "Reference: "
        + ", ".join(
            f"{name} {ref['tag']} (`{ref['sha'][:12]}`)" for name, ref in reference.items()
        )
    )
    lines.append("")

    total = len(tools)
    lines.append("| state | tools | share |")
    lines.append("|---|---:|---:|")
    for state in ("oracle-backed", "golden-pending", "unchecked", "not started"):
        lines.append(f"| {state} | {counted[state]} | {counted[state] / total:.1%} |")
    lines.append(f"| **total** | **{total}** | |")
    lines.append("")
    lines.append(
        "`oracle-backed` means CI re-derives the tool's goldens in the pinned container on every "
        "run and compares them. It does **not** mean the tool is byte-identical over its whole "
        "argument surface: that is what the argument-coverage column is for. A t-wise array "
        "(`tools/coverage/covering.py`, sized in "
        "[what-pairwise-coverage-costs.md](what-pairwise-coverage-costs.md)) is run against the "
        "reference and the port, and the column reports the fraction of its rows on which the two "
        "answered identically, rejections included. `not measured` means the array has never been "
        "run against a port binary, which is still true of most tools here."
    )
    lines.append("")

    if missing:
        lines.append("## Ports not found")
        lines.append("")
        for name, path in missing:
            lines.append(f"- `{name}`: no manifest at `{path}`")
        lines.append("")

    for origin in ("picard", "gatk"):
        subset = [r for r in rows if r[0]["origin"] == origin]
        started = [r for r in subset if r[1] != "not started"]
        lines.append(f"## {origin}-origin tools ({len(started)} of {len(subset)} started)")
        lines.append("")
        if started:
            lines.append("| tool | archetype | state | suites | cases | argument coverage |")
            lines.append("|---|---|---|---|---:|---|")
            for tool, state, entry in started:
                suites = ", ".join(sorted(set(entry["suites"])))
                lines.append(
                    f"| `{tool['name']}` | {tool['archetype']} | {state} | {suites} | "
                    f"{entry['cases']} | {coverage_cell(entry, manifests)} |"
                )
            lines.append("")
        not_started = [r[0]["name"] for r in subset if r[1] == "not started"]
        lines.append(f"<details><summary>{len(not_started)} not started</summary>")
        lines.append("")
        lines.append(", ".join(f"`{n}`" for n in not_started))
        lines.append("")
        lines.append("</details>")
        lines.append("")

    return "\n".join(lines) + "\n"


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true")
    ap.add_argument(
        "--refresh",
        action="store_true",
        help="re-vendor the ports' suite lists from the sibling working copies, then regenerate",
    )
    args = ap.parse_args(argv)

    if args.refresh:
        written, absent = refresh()
        for name, count in written:
            print(f"vendored {name}: {count} suites")
        for name, path in absent:
            print(f"  no manifest for {name} at {path}")

    with open(INVENTORY) as fh:
        inventory = json.load(fh)
    manifests, missing = load_manifests()
    rendered = render(inventory, manifests, missing)

    if args.check:
        current = STATUS.read_text() if STATUS.exists() else ""
        if current != rendered:
            print("docs/STATUS.md is stale: regenerate with tools/dashboard/generate.py")
            return 1
        print("docs/STATUS.md matches the inventory and the manifests")
        return 0

    STATUS.parent.mkdir(exist_ok=True)
    STATUS.write_text(rendered)
    states = tool_states(manifests)
    print(
        f"wrote {STATUS}: {len(inventory['tools'])} tools, "
        f"{len(states)} with a suite, {len(manifests)} port manifest(s) read"
    )
    for name, path in missing:
        print(f"  no manifest for {name} at {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

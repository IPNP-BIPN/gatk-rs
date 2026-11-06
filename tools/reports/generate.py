#!/usr/bin/env python3
"""Generate the per-tool validation report from what the repository commits.

`docs/STATUS.md` says which state a tool is in; this says what that state rests on, so that the 311
green suites become a claim somebody else can check: which inputs, how many records, how many files
compared byte for byte, how much of the argument surface, and what was left out and why.

Every number is read from a committed source, and a number that no source holds says so:

* suites and their tools, from the manifests (the siblings' through the summaries the dashboard
  vendors, which carry no goldens, so their records read `not vendored`);
* records and files, by counting the golden's lines. A golden that labels its rows by `argv` is
  counted per tool; one that does not is counted for the suite, and the column says `suite`;
* argument coverage and the `t` level, from each port's `measured.json`;
* inputs and held-out arguments, from `tools/coverage/fixtures.json`. An argument with no values is
  one no array varies, and the entry's `$why` says why;
* quarantined fields, from each suite's `compare` rules;
* fuzzer coverage, from `tools/fuzz/measured.json`: the reference's JaCoCo probes the covering
  array's rows reach, how many more the fuzzer's mutants reach, and on how many kept cases the
  port diverged.

    python3 tools/reports/generate.py            # writes docs/VALIDATION.md
    python3 tools/reports/generate.py --check    # fails if the committed file is stale
"""

import argparse
import gzip
import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
INVENTORY = REPO / "tools" / "inventory" / "generated" / "inventory.json"
MANIFEST = REPO / "tools" / "conformance" / "manifest.json"
MEASURED = REPO / "tools" / "coverage" / "measured.json"
FUZZED = REPO / "tools" / "fuzz" / "measured.json"
FIXTURES = REPO / "tools" / "coverage" / "fixtures.json"
PORTS = REPO / "tools" / "dashboard" / "ports"
REPORT = REPO / "docs" / "VALIDATION.md"


def golden_rows(path):
    """Lines of a golden, or nothing when it is not there."""
    full = REPO / path
    if not full.is_file():
        return []
    opener = gzip.open if full.suffix == ".gz" else open
    with opener(full, "rt", encoding="utf-8", errors="replace") as handle:
        return [line for line in handle if line.strip() and not line.startswith("#")]


def count_suite(suite):
    """(records, files, {tool: (records, files)}) over a suite's goldens."""
    records = files = 0
    per_tool = {}
    labels = {}
    # A pending suite's golden, if it has one, predates what the suite now runs.
    if suite["status"] == "golden-pending":
        return 0, 0, per_tool
    for case in suite["cases"]:
        if not case.get("golden"):
            continue
        for row in golden_rows(case["golden"]):
            records += 1
            fields = row.split("\t", 3)
            is_file = fields[0] == "file"
            files += is_file
            if fields[0] == "argv" and len(fields) > 2:
                labels[fields[1]] = fields[2].split("\\t", 1)[0].strip()
            tool = labels.get(fields[1]) if len(fields) > 1 else None
            if tool:
                found = per_tool.setdefault(tool, [0, 0])
                found[0] += 1
                found[1] += is_file
    return records, files, per_tool


def quarantined(suite):
    """What a suite's comparison sets aside, by rule."""
    return [rule["rule"] for rule in suite.get("compare", {}).get("rules", [])]


def gather():
    inventory = json.loads(INVENTORY.read_text())
    manifest = json.loads(MANIFEST.read_text())
    measured = json.loads(MEASURED.read_text()).get("tools", {})
    fixtures = json.loads(FIXTURES.read_text()).get("per_tool", {})
    fuzzed = json.loads(FUZZED.read_text()).get("tools", {}) if FUZZED.exists() else {}
    tools = {}
    for entry in inventory["tools"]:
        name = entry["name"] if isinstance(entry, dict) else entry
        tools[name] = {"suites": [], "records": 0, "files": 0, "by_suite": False, "quarantine": set(),
                       "port": None, "coverage": None, "fuzz": fuzzed.get(name)}
    for suite in manifest["suites"]:
        records, files, per_tool = count_suite(suite)
        for name in suite.get("tools", []):
            tool = tools.get(name)
            if tool is None:
                continue
            tool["port"] = "gatk-rs"
            tool["suites"].append((suite["id"], suite["status"]))
            tool["quarantine"].update(quarantined(suite))
            if suite["status"] == "golden-pending":
                continue
            if name in per_tool:
                tool["records"] += per_tool[name][0]
                tool["files"] += per_tool[name][1]
            else:
                # Shared with the suite's other tools, or a golden with no per-run labels.
                tool["records"] += records
                tool["files"] += files
                tool["by_suite"] = True
    for path in sorted(PORTS.glob("*.json")):
        port = json.loads(path.read_text())
        for suite in port["suites"]:
            for name in suite.get("tools", []):
                tool = tools.get(name)
                if tool is None or tool["port"] == "gatk-rs":
                    continue
                tool["port"] = path.stem
                tool["suites"].append((suite["id"], suite["status"]))
        for name, value in port.get("coverage", {}).items():
            if name in tools and tools[name]["port"] == path.stem:
                tools[name]["coverage"] = value
    for name, value in measured.items():
        if name in tools and tools[name]["coverage"] is None:
            tools[name]["coverage"] = value
    return inventory, tools, fixtures


def coverage_text(value):
    if not value:
        return "not measured", "-"
    share = value.get("share")
    rows = value.get("rows")
    matched = value.get("matched")
    return f"{matched}/{rows} ({share:.1%})" if share is not None else f"{matched}/{rows}", str(value.get("t", "-"))


def fuzz_text(value):
    if not value:
        return "not measured"
    gained = value["probes_total"] - value["probes_from_seeds"]
    text = f"{value['probes_total']} probes (+{gained})"
    if "diverged" in value:
        text += f", {value['diverged']}/{value['replayed']} diverge"
    return text


def fixture_summary(entry):
    """Distinct input files an array draws on, and the arguments no array varies."""
    if not entry:
        return [], []
    inputs = set()
    for key, values in entry.items():
        if key.startswith("$"):
            continue
        for value in values:
            if isinstance(value, str) and value.startswith("/work/fixtures/"):
                inputs.add(value.rsplit("/", 1)[1])
    held = sorted(k for k, values in entry.items() if not k.startswith("$") and not values)
    return sorted(inputs), held


def render(inventory, tools, fixtures):
    reference = inventory.get("reference")
    lines = [
        "# Validation report",
        "",
        "Generated by `tools/reports/generate.py` from committed data. Do not edit by hand.",
        "",
        f"Reference: `{json.dumps(reference, sort_keys=True)}`",
        "",
        "How to read it: **records** and **files** count the lines of the goldens CI re-derives from "
        "the pinned container on every run, and the files among them, compared byte for byte. A "
        "count marked `suite` belongs to a suite the tool shares, or whose golden does not label "
        "its runs per tool. **Argument rows** are the rows of the tool's covering array whose "
        "outcome matched the reference, at the **t** level shown. **Fuzzer branches** are the reference's JaCoCo "
        "probes reached (a proxy for branches, not a percentage), how many of them only the "
        "fuzzer's mutants reached, and how many of the cases it kept the port answered "
        "differently, re-run out of process (`tools/fuzz/run_fuzz.py`). A sibling port's goldens are not vendored "
        "here, so its counts read `not vendored`; its own repository holds them.",
        "",
        "Each run's inputs, toolchain and container are in the `provenance` artefact CI keeps "
        "(`tools/provenance/record.py`).",
        "",
        "| tool | port | state | records | files | argument rows | t | fuzzer branches | inputs | held-out arguments | quarantined |",
        "|---|---|---|---:|---:|---:|---:|---|---:|---:|---|",
    ]
    details = []
    for name in sorted(tools):
        tool = tools[name]
        states = {status for _, status in tool["suites"]}
        state = (
            "not started" if not states
            else "oracle-backed" if "oracle-backed" in states
            else "golden-pending" if "golden-pending" in states
            else "unchecked"
        )
        if tool["port"] in (None, "gatk-rs"):
            mark = " suite" if tool["by_suite"] else ""
            known = state not in ("not started", "golden-pending")
            records = f"{tool['records']}{mark}" if known else "-"
            files = f"{tool['files']}{mark}" if known else "-"
        else:
            records = files = "not vendored"
        covered, level = coverage_text(tool["coverage"])
        inputs, held = fixture_summary(fixtures.get(name))
        quarantine = ", ".join(sorted(tool["quarantine"])) or "none"
        lines.append(
            f"| {name} | {tool['port'] or '-'} | {state} | {records} | {files} | {covered} | {level} "
            f"| {fuzz_text(tool['fuzz'])} | {len(inputs)} | {len(held)} | {quarantine} |"
        )
        entry = fixtures.get(name)
        if entry:
            why = " ".join(entry.get("$why", []))
            why = re.sub(r"\s+", " ", why).strip()
            block = [f"### {name}", ""]
            block.append("Suites: " + (", ".join(f"`{i}`" for i, _ in tool["suites"]) or "none") + ".")
            block.append("")
            block.append("Inputs: " + (", ".join(f"`{i}`" for i in inputs) or "none") + ".")
            block.append("")
            block.append("Held out: " + (", ".join(f"`{h}`" for h in held) or "none") + ".")
            if why:
                block.extend(["", why])
            block.append("")
            details.extend(block)
    lines.extend(["", "## What each array draws on and leaves out", ""])
    lines.extend(details)
    return "\n".join(lines).rstrip() + "\n"


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    text = render(*gather())
    if args.check:
        if not REPORT.is_file() or REPORT.read_text() != text:
            print("docs/VALIDATION.md is stale: run python3 tools/reports/generate.py")
            return 1
        print("docs/VALIDATION.md is current")
        return 0
    REPORT.write_text(text)
    print(f"wrote {REPORT}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

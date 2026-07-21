#!/usr/bin/env python3
"""Derive the tool inventory mechanically from the pinned reference's own documentation.

375 tools and ~10,800 arguments are never maintained by hand: a hand-written list of that size
drifts within weeks and there is no way to notice. This reads the `gatkdoc/*.json` that GATK
itself generates, so regenerating is the same thing as re-deriving from the pinned reference.

Emits:
  inventory.json  normalized tools, arguments, types, defaults, enum members, bounds
  summary.md      counts by origin, group and archetype, for the progress dashboard

Usage:
  ./generate_inventory.py <gatkdoc-dir> [-o outdir]
"""

import argparse
import collections
import json
import pathlib
import re
import sys

# Pinned reference versions this inventory is derived from. Recorded in the output so an
# inventory can never be silently attributed to the wrong reference.
REFERENCE = {
    "gatk": {"tag": "4.6.2.0", "sha": "76edc75c26504da94bbaee66584e107e76ee15de"},
    "picard": {"tag": "3.4.0", "sha": "6c3f23bc2e0d229d75e9f9b04200396bcd067526"},
    "htsjdk": {"tag": "4.2.0", "sha": "4cc010022ac038fb30f26e6f9717fabff3e808c1"},
}

# Maps a documented tool group onto the port archetype it belongs to. The archetype is the
# unit of amortization: it is ported once, and members after the first cost only the delta.
GROUP_TO_ARCHETYPE = {
    "Metrics": "metrics-collector",
    "Diagnostics and Quality Control": "reporting-walker",
    "Read Data Manipulation": "record-transform",
    "Variant Annotations": "annotation",
    "Variant Evaluation and Refinement": "variant-walker",
    "Variant Manipulation": "variant-transform",
    "Variant Filtering": "variant-transform",
    "Coverage Analysis": "locus-walker",
    "Reference": "reference-utility",
    "Intervals Manipulation": "interval-utility",
    "Short Variant Discovery": "assembly-caller",
    "Structural Variant Discovery": "sv-caller",
    "Copy Number Variant Discovery": "cnv-segmentation",
    "Methylation-Specific Tools": "record-transform",
    "Genotyping Arrays Manipulation": "array-utility",
    "Base Calling": "base-calling",
    "Metagenomics": "metagenomics",
    "Flow Based Tools": "flow-based",
    "Flow Annotations": "annotation",
    "Read Filters": "read-filter",
    "Other": "unclassified",
}


def split_class(stem):
    """`org_broadinstitute_hellbender_tools_X` -> (package, ClassName)."""
    parts = stem.split("_")
    return ".".join(parts[:-1]), parts[-1]


def parse_bool(v):
    if isinstance(v, bool):
        return v
    return str(v).strip().lower() in ("true", "yes")


def normalize_default(v):
    if v is None:
        return None
    s = str(v).strip()
    return None if s in ("", "null", "NA", "none") else s


def normalize_argument(a):
    options = a.get("options") or []
    enum_members = [o["name"] for o in options if isinstance(o, dict) and "name" in o]
    synonyms = a.get("synonyms")
    return {
        "name": a.get("name"),
        "short": None if synonyms in (None, "NA", "") else synonyms,
        "type": a.get("type"),
        "required": parse_bool(a.get("required")),
        "kind": a.get("kind"),
        "default": normalize_default(a.get("defaultValue")),
        "enum_members": enum_members or None,
        "min": normalize_default(a.get("minValue")),
        "max": normalize_default(a.get("maxValue")),
        "deprecated": parse_bool(a.get("deprecated")),
    }


def branch_name(origin, tool_name):
    slug = re.sub(r"[^a-z0-9]+", "", tool_name.lower())
    return f"tool/{origin}-{slug}"


def classify(stem, name, group, n_args):
    """Separate the four things gatkdoc documents side by side.

    The distinction that matters, and that a naive split gets wrong: the `Metrics` group holds
    *metric definition classes* (`InsertSizeMetrics`, `AlignmentSummaryMetrics`), which are the
    output-file schemas emitted by collector tools, not runnable tools. They have no arguments
    and no CLI. Counting them as tools inflates the tool count by 57 and reports a whole
    archetype as having zero arguments, which is the tell.
    """
    if group == "Read Filters" or "ReadFilter" in name:
        return "read-filter"
    if name.endswith("Annotation") or "_annotator_" in stem:
        return "annotation"
    if group == "Metrics" and n_args == 0:
        return "metrics-definition"
    return "tool"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gatkdoc", type=pathlib.Path)
    ap.add_argument("-o", "--outdir", type=pathlib.Path, default=pathlib.Path("."))
    args = ap.parse_args()

    if not args.gatkdoc.is_dir():
        sys.exit(f"not a directory: {args.gatkdoc}")

    tools, filters, annotations, metric_defs = [], [], [], []
    skipped = []

    for path in sorted(args.gatkdoc.glob("*.json")):
        try:
            d = json.loads(path.read_text())
        except Exception as e:  # a malformed doc must be visible, not silently dropped
            skipped.append((path.name, str(e)))
            continue
        if not isinstance(d, dict) or "name" not in d:
            skipped.append((path.name, "no 'name' key"))
            continue

        stem = path.stem
        package, cls = split_class(stem)
        # gatkdoc renders Picard-sourced tools with a " (Picard)" display suffix. The
        # invocable name is what the CLI dispatches on and what branch names derive from, so
        # strip it and keep it only as an origin signal. Leaving it in silently produces an
        # inventory whose tool names do not exist.
        raw_name = d["name"]
        picard_suffix = raw_name.endswith(" (Picard)")
        name = raw_name[: -len(" (Picard)")] if picard_suffix else raw_name
        origin = "picard" if (picard_suffix or stem.startswith("picard")) else "gatk"
        group = d.get("group") or "Other"
        raw_args = d.get("arguments") or []
        if isinstance(raw_args, dict):  # defensive: some doclet versions group by kind
            raw_args = [a for lst in raw_args.values() if isinstance(lst, list) for a in lst]

        entry = {
            "name": name,
            "display_name": raw_name,
            "origin": origin,
            "package": package,
            "class": f"{package}.{cls}",
            "group": group,
            "archetype": GROUP_TO_ARCHETYPE.get(group, "unclassified"),
            "spark": "Spark" in name,
            "beta": parse_bool(d.get("beta")),
            "experimental": parse_bool(d.get("experimental")),
            "deprecated": parse_bool(d.get("deprecated")),
            "summary": d.get("summary"),
            "branch": branch_name(origin, name),
            "arguments": [normalize_argument(a) for a in raw_args],
        }

        kind = classify(stem, name, group, len(entry["arguments"]))
        entry["kind"] = kind
        {"read-filter": filters, "annotation": annotations,
         "metrics-definition": metric_defs, "tool": tools}[kind].append(entry)

    total_args = sum(len(t["arguments"]) for t in tools)
    inventory = {
        "reference": REFERENCE,
        "generated_from": str(args.gatkdoc),
        "counts": {
            "tools": len(tools),
            "read_filters": len(filters),
            "annotations": len(annotations),
            "metrics_definitions": len(metric_defs),
            "arguments": total_args,
            "skipped": len(skipped),
        },
        "tools": tools,
        "read_filters": filters,
        "annotations": annotations,
        "metrics_definitions": metric_defs,
        "skipped": skipped,
    }

    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "inventory.json").write_text(json.dumps(inventory, indent=2, sort_keys=False))

    write_summary(inventory, args.outdir)

    print(f"tools={len(tools)} filters={len(filters)} annotations={len(annotations)} "
          f"metric_defs={len(metric_defs)} arguments={total_args} skipped={len(skipped)}")
    print(f"wrote {args.outdir / 'inventory.json'} and {args.outdir / 'summary.md'}")


def write_summary(inventory, outdir):
    tools = inventory["tools"]
    filters = inventory["read_filters"]
    annotations = inventory["annotations"]
    metric_defs = inventory["metrics_definitions"]
    skipped = inventory["skipped"]
    total_args = sum(len(t["arguments"]) for t in tools)
    undocumented = [t for t in tools if not t.get("documented", True)]

    by_origin = collections.Counter(t["origin"] for t in tools)
    by_arch = collections.Counter(t["archetype"] for t in tools)
    args_by_arch = collections.Counter()
    for t in tools:
        args_by_arch[t["archetype"]] += len(t["arguments"])
    spark = sum(1 for t in tools if t["spark"])

    lines = [
        "# Tool inventory",
        "",
        "Generated by `tools/inventory/generate_inventory.py`. Do not edit by hand.",
        "",
        f"Reference: GATK {REFERENCE['gatk']['tag']} (`{REFERENCE['gatk']['sha'][:12]}`), "
        f"Picard {REFERENCE['picard']['tag']}, htsjdk {REFERENCE['htsjdk']['tag']}",
        "",
        "| | Count |",
        "|---|---:|",
        f"| Runnable tools (CLI ground truth) | {len(tools)} |",
        f"| of which undocumented in gatkdoc | {len(undocumented)} |",
        f"| of which Spark | {spark} |",
        f"| of which GATK-origin | {by_origin['gatk']} |",
        f"| of which Picard-origin | {by_origin['picard']} |",
        f"| Read filters | {len(filters)} |",
        f"| Annotations | {len(annotations)} |",
        f"| Metric definition classes | {len(metric_defs)} |",
        f"| Arguments (tools only) | {total_args} |",
        f"| Mean arguments per tool | {total_args / max(len(tools), 1):.1f} |",
        "",
        "## By archetype",
        "",
        "The archetype is the unit of amortization: ported once, then each further member",
        "costs only the delta. Largest first, because that is the order that pays.",
        "",
        "| Archetype | Tools | Arguments |",
        "|---|---:|---:|",
    ]
    for arch, n in by_arch.most_common():
        lines.append(f"| `{arch}` | {n} | {args_by_arch[arch]} |")

    lines += ["", "## Largest argument surfaces", "", "| Tool | Origin | Arguments |", "|---|---|---:|"]
    for t in sorted(tools, key=lambda x: -len(x["arguments"]))[:12]:
        lines.append(f"| `{t['name']}` | {t['origin']} | {len(t['arguments'])} |")

    if skipped:
        lines += ["", "## Skipped", ""]
        lines += [f"- `{n}`: {why}" for n, why in skipped[:20]]

    (outdir / "summary.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

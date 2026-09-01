#!/usr/bin/env python3
"""Generate t-wise covering arrays over a tool's arguments, from the generated inventory.

Why this exists: HaplotypeCaller has 174 arguments, so "every parameter" cannot mean every
combination (2^174 of them for the booleans alone). Combinatorial interaction testing is the
established answer: cover every *t-way* interaction at least once, with a number of rows that
grows logarithmically in the number of arguments rather than exponentially. The program commits to
t=2 everywhere and t=3 on the critical path, and this is what produces those rows.

The algorithm is IPOG (In-Parameter-Order-General, Lei et al.), which builds a t-way array for the
first t parameters and then grows it one parameter at a time: horizontally, by extending existing
rows with the value that covers the most still-uncovered tuples, and vertically, by adding rows for
whatever is left. It is deterministic here: ties break on the first candidate in declared order, so
the same inventory always yields the same array and a regenerated array is diffable.

Correctness is not assumed. `--verify` enumerates every t-way tuple of the domain and asserts the
array covers it; that check is exhaustive rather than sampled, and it is what makes a coverage
claim a measurement.

    python3 tools/coverage/covering.py --tool CollectQualityYieldMetrics --t 2
    python3 tools/coverage/covering.py --tool CollectQualityYieldMetrics --t 2 --verify
    python3 tools/coverage/covering.py --tool CollectQualityYieldMetrics --json out.json
"""

import argparse
import itertools
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import domains  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
INVENTORY = REPO / "tools" / "inventory" / "generated" / "inventory.json"


def load_tool(name, inventory_path=None):
    with open(inventory_path or INVENTORY) as fh:
        inventory = json.load(fh)
    for tool in inventory["tools"]:
        if tool["name"] == name:
            return tool
    raise SystemExit(f"no tool {name!r} in the inventory")


def tuples_of(params, domain, t):
    """Every t-way combination of values, as (parameter-tuple, value-index-tuple)."""
    for combo in itertools.combinations(params, t):
        index_ranges = [range(len(domain[p])) for p in combo]
        for indices in itertools.product(*index_ranges):
            yield combo, indices


def _covered_by(row, combo, indices):
    return all(row.get(p) == i for p, i in zip(combo, indices))


def constraints_for(declared, tool_name):
    """The clauses that apply to this tool.

    A clause without `tools` is global. Scoping matters: a queryname-sorted BAM is an invalid input
    for a coordinate-requiring collector and a perfectly ordinary one for SortSam or SamToFastq, so
    forbidding it everywhere would delete a dimension those tools need.
    """
    out = []
    for clause in declared:
        scope = clause.get("tools")
        if scope is None or tool_name in scope:
            out.append(clause)
    return out


def compile_constraints(declared, domain, excluded):
    """Turn declared value constraints into index constraints over the array's parameters.

    A constraint forbids a *combination*: `--FLOW_MODE=true` alone, or a queryname input together
    with `--ASSUME_SORTED=false`. Without them, an array over a tool that refuses such rows spends
    most of its rows being rejected: eight of `CollectQualityYieldMetrics`'s eleven were, and a
    coverage figure that counts rows which could never produce output overstates itself
    (picard-rs decision 0009).

    Three cases for a constraint that names a parameter the array does not vary:

    * the parameter is held at the forbidden value, so the clause is always half-satisfied and the
      remaining keys carry it;
    * the parameter is held at some other value, so the clause can never fire and is dropped;
    * the parameter is not in the tool at all, so the clause is dropped.

    All three are reported rather than applied silently, because a constraint that quietly does
    nothing is indistinguishable from one that works.
    """
    held = {e["argument"]: e.get("held_at") for e in excluded}
    compiled, notes = [], []
    for clause in declared:
        forbid = clause.get("forbid", {})
        indices, drop = {}, None
        for name, value in forbid.items():
            if name in domain:
                match = [i for i, (_, v) in enumerate(domain[name]) if str(v) == str(value)]
                if not match:
                    drop = f"{name}={value} is not in that argument's domain"
                    break
                indices[name] = match[0]
            elif name in held:
                if str(held[name]) != str(value):
                    drop = f"{name} is held at {held[name]}, so the clause can never fire"
                    break
                notes.append(f"{name} is held at the forbidden value {value}; clause still applies")
            else:
                drop = f"{name} is not an argument of this tool"
                break
        if drop:
            notes.append(f"dropped constraint {forbid}: {drop}")
            continue
        if not indices:
            # Every key was held at its forbidden value: the whole tool is forbidden, which is a
            # statement about the fixtures rather than about the array.
            notes.append(f"constraint {forbid} is satisfied by the held values alone")
            continue
        compiled.append({"indices": indices, "why": clause.get("why", "")})
    return compiled, notes


def violates(assignment, constraints):
    """True when a (possibly partial) assignment satisfies every key of some forbidden clause."""
    for clause in constraints:
        indices = clause["indices"]
        if all(assignment.get(p) == i for p, i in indices.items()):
            return True
    return False


def ipog(params, domain, t, constraints=()):
    """Build a t-way covering array. Rows are {parameter: value index}; None means unset."""
    if len(params) < t:
        raise SystemExit(f"{len(params)} parameters cannot support t={t}")

    # Step 1: the exhaustive array over the first t parameters, minus the forbidden rows.
    head, rest = params[:t], params[t:]
    rows = [
        dict(zip(head, indices))
        for indices in itertools.product(*[range(len(domain[p])) for p in head])
        if not violates(dict(zip(head, indices)), constraints)
    ]
    if not rows:
        raise SystemExit(
            f"the constraints forbid every combination of the first {t} arguments "
            f"({', '.join(head)}); reorder or relax them"
        )

    for param in rest:
        covered_params = [p for p in rows[0]]
        # The tuples this new parameter introduces: it, paired with any t-1 of the covered ones.
        pending = set()
        for others in itertools.combinations(covered_params, t - 1):
            for indices in itertools.product(*[range(len(domain[p])) for p in others]):
                for value in range(len(domain[param])):
                    assignment = dict(zip(others + (param,), indices + (value,)))
                    if violates(assignment, constraints):
                        continue
                    pending.add((others + (param,), indices + (value,)))

        # Horizontal growth: extend each existing row with its best value.
        for row in rows:
            if not pending:
                break
            best_value, best_gain, best_hit = None, -1, set()
            for value in range(len(domain[param])):
                if violates({**row, param: value}, constraints):
                    continue
                hit = {
                    (combo, indices)
                    for combo, indices in pending
                    if indices[-1] == value
                    and all(row.get(p) == i for p, i in zip(combo[:-1], indices[:-1]))
                }
                if len(hit) > best_gain:
                    best_value, best_gain, best_hit = value, len(hit), hit
            if best_value is None:
                # Every value of this parameter would make the row forbidden. The row cannot be
                # extended, so it is dropped rather than emitted as a run that will be rejected.
                row["__drop__"] = True
                continue
            row[param] = best_value
            pending -= best_hit
        rows = [r for r in rows if not r.pop("__drop__", False)]

        # Vertical growth: one row per tuple still uncovered, merged where rows do not conflict.
        for combo, indices in sorted(pending, key=lambda x: (tuple(x[0]), x[1])):
            assignment = dict(zip(combo, indices))
            for row in rows:
                if all(row.get(p) in (None, v) for p, v in assignment.items()) and not violates(
                    {**row, **assignment}, constraints
                ):
                    row.update(assignment)
                    break
            else:
                rows.append(dict(assignment))

    # Any parameter left unset in a row takes its first permitted value: an unset parameter would
    # mean the tool runs with the argument absent, which is a different test than the array claims.
    for row in rows:
        for param in params:
            if param in row:
                continue
            for value in range(len(domain[param])):
                if not violates({**row, param: value}, constraints):
                    row[param] = value
                    break
            else:
                row[param] = 0
    return rows


def verify(rows, params, domain, t, constraints=()):
    """Exhaustively check that every permitted t-way tuple appears in some row.

    Forbidden tuples are counted separately and never reported missing: a tuple the tool refuses is
    not coverage that was skipped, it is coverage that does not exist. Reporting it as missing
    would make a correct array look broken; folding it into the denominator would make an
    incomplete one look finished.
    """
    missing, forbidden, total = [], 0, 0
    for combo, indices in tuples_of(params, domain, t):
        if violates(dict(zip(combo, indices)), constraints):
            forbidden += 1
            continue
        total += 1
        if not any(_covered_by(row, combo, indices) for row in rows):
            missing.append((combo, indices))
    return total, missing, forbidden


def as_arguments(row, domain):
    """Render one row as the argument list a runner would pass."""
    return [f"{param}={domain[param][index][1]}" for param, index in sorted(row.items())]


def report(tool, domain, excluded, rows, t, total_tuples, missing, numeric_policy,
           forbidden=0, constraint_notes=()):
    covered = total_tuples - len(missing)
    return {
        "tool": tool["name"],
        "origin": tool["origin"],
        "t": t,
        "numeric_policy": numeric_policy,
        "arguments_total": len(tool["arguments"]),
        "arguments_in_array": len(domain),
        "arguments_excluded": len(excluded),
        "excluded": excluded,
        "rows": len(rows),
        "tuples_total": total_tuples,
        "tuples_covered": covered,
        "tuples_forbidden": forbidden,
        "constraint_notes": list(constraint_notes),
        "coverage": round(covered / total_tuples, 6) if total_tuples else 0.0,
        "array": [
            {
                "row": n,
                "labels": {p: domain[p][i][0] for p, i in sorted(row.items())},
                "arguments": as_arguments(row, domain),
            }
            for n, row in enumerate(rows)
        ],
    }


def summarize_all(inventory_path, t, fixtures, numeric_policy, limit=None, excluded_by_hand=None,
                  declared_constraints=(), per_tool_fixtures=None):
    """Size the whole coverage programme: rows per tool, and the total across the inventory.

    This is the number the plan needs and did not have. A tool's row count is how many oracle runs
    its t-wise claim costs, so the sum is the size of the coverage programme in runs, and the
    excluded-argument column says how much of each tool the claim leaves out.
    """
    with open(inventory_path or INVENTORY) as fh:
        inventory = json.load(fh)

    rows_total = args_in = args_out = 0
    per_tool = []
    for tool in inventory["tools"][:limit]:
        # A tool that declares its own fixture values gets them here, so sizing the programme
        # counts the array a tool will actually run rather than the one its shared paths imply.
        own = domains.for_tool(fixtures, per_tool_fixtures or {}, tool["name"])
        domain, excluded = domains.build(tool, own, numeric_policy, excluded_by_hand)
        params = sorted(domain)
        if len(params) < t:
            per_tool.append(
                {"tool": tool["name"], "rows": None, "in_array": len(params),
                 "excluded": len(excluded), "note": f"fewer than t={t} usable arguments"}
            )
            args_out += len(excluded)
            continue
        constraints, _ = compile_constraints(
            constraints_for(declared_constraints, tool["name"]), domain, excluded
        )
        rows = ipog(params, domain, t, constraints)
        rows_total += len(rows)
        args_in += len(params)
        args_out += len(excluded)
        per_tool.append(
            {"tool": tool["name"], "rows": len(rows), "in_array": len(params),
             "excluded": len(excluded), "note": None}
        )
    return {
        "t": t,
        "numeric_policy": numeric_policy,
        "tools": len(per_tool),
        "tools_sized": sum(1 for e in per_tool if e["rows"] is not None),
        "oracle_runs": rows_total,
        "arguments_in_arrays": args_in,
        "arguments_excluded": args_out,
        "per_tool": per_tool,
    }


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tool")
    ap.add_argument("--all", action="store_true", help="size every tool in the inventory")
    ap.add_argument("--t", type=int, default=2)
    ap.add_argument(
        "--numeric-policy",
        choices=("strict", "perturb"),
        default="strict",
        help="what to do with a numeric argument that has a default and no declared bounds. "
        "strict holds it at its default and excludes it; perturb also offers default +/- one step "
        "(see domains.py for why that is legitimate in a differential test)",
    )
    ap.add_argument("--inventory")
    ap.add_argument(
        "--fixtures",
        help="JSON mapping argument name to the values a repository is willing to pass for it. "
        "Path-typed arguments have no domain without one.",
    )
    ap.add_argument("--verify", action="store_true", help="exhaustively check the array")
    ap.add_argument("--json", help="write the full report here")
    args = ap.parse_args(argv)
    if not args.tool and not args.all:
        ap.error("pass --tool <name> or --all")

    raw_fixtures = json.load(open(args.fixtures)) if args.fixtures else {}
    fixtures, hand_excluded, per_tool_fixtures = domains.load_fixtures(raw_fixtures)
    declared_constraints = (
        raw_fixtures.get("constraints", []) if isinstance(raw_fixtures, dict) else []
    )

    if args.all:
        summary = summarize_all(
            args.inventory, args.t, fixtures, args.numeric_policy,
            excluded_by_hand=hand_excluded, declared_constraints=declared_constraints,
            per_tool_fixtures=per_tool_fixtures,
        )
        print(
            f"t={summary['t']} policy={summary['numeric_policy']} "
            f"tools sized={summary['tools_sized']}/{summary['tools']} "
            f"oracle runs={summary['oracle_runs']} "
            f"arguments in arrays={summary['arguments_in_arrays']} "
            f"excluded={summary['arguments_excluded']}"
        )
        widest = sorted(
            (e for e in summary["per_tool"] if e["rows"]), key=lambda e: -e["rows"]
        )[:10]
        print(f"\n{'tool':44} {'rows':>6} {'in array':>9} {'excluded':>9}")
        for entry in widest:
            print(
                f"{entry['tool']:44} {entry['rows']:>6} {entry['in_array']:>9} "
                f"{entry['excluded']:>9}"
            )
        if args.json:
            with open(args.json, "w") as fh:
                json.dump(summary, fh, indent=2)
            print(f"\nwrote {args.json}")
        return 0

    tool = load_tool(args.tool, args.inventory)
    fixtures = domains.for_tool(fixtures, per_tool_fixtures, tool["name"])
    domain, excluded = domains.build(tool, fixtures, args.numeric_policy, hand_excluded)
    params = sorted(domain)

    if len(params) < args.t:
        print(
            f"{tool['name']}: only {len(params)} of {len(tool['arguments'])} arguments have a "
            f"declared domain, which is fewer than t={args.t}. Supply fixtures."
        )
        for entry in excluded[:10]:
            print(f"  excluded {entry['argument']} ({entry['type']}): {entry['reason']}")
        return 1

    constraints, notes = compile_constraints(
        constraints_for(declared_constraints, tool["name"]), domain, excluded
    )
    rows = ipog(params, domain, args.t, constraints)
    total, missing, forbidden = (0, [], 0)
    if args.verify or args.json:
        total, missing, forbidden = verify(rows, params, domain, args.t, constraints)

    for note in notes:
        print(f"constraint note: {note}")
    print(
        f"{tool['name']}: t={args.t} policy={args.numeric_policy} "
        f"constraints={len(constraints)} rows={len(rows)} "
        f"arguments in array={len(params)}/{len(tool['arguments'])} excluded={len(excluded)}"
    )
    if args.verify:
        print(
            f"tuples={total} covered={total - len(missing)} missing={len(missing)} "
            f"forbidden by constraints={forbidden}"
        )
        for combo, indices in missing[:5]:
            print(f"  MISSING {combo} = {indices}")
    if excluded:
        print("excluded arguments (each one narrows the claim):")
        for entry in excluded[:10]:
            print(f"  {entry['argument']} ({entry['type']}): {entry['reason']}")
        if len(excluded) > 10:
            print(f"  ... and {len(excluded) - 10} more")

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(
                report(
                    tool, domain, excluded, rows, args.t, total, missing, args.numeric_policy,
                    forbidden, notes,
                ),
                fh,
                indent=2,
            )
        print(f"wrote {args.json}")

    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

#!/usr/bin/env python3
"""Turn an argument's declared schema into the set of values a covering array may use.

The rule that governs this file: **a value domain is declared, never guessed.** A covering array
over invented values proves nothing, and the failure is silent, because the array still looks
covered. So every value here traces to something the inventory records (`enum_members`, `min`,
`max`, `default`) or to a fixture the repository owns, and anything that cannot be given a domain
is *excluded by name with a reason* rather than quietly given one value.

The exclusions are as much a part of the coverage claim as the array is: a tool reported as
"t=2 covered" with half its arguments excluded is covered over half a tool. `report()` therefore
carries both numbers.
"""

# Types whose values are file paths. A covering array cannot invent these: a path must exist, be
# readable, and carry content the tool accepts. They take a fixture, or the argument is excluded.
PATH_TYPES = {
    "File",
    "GATKPath",
    "PicardHtsPath",
    "Path",
    "String[]",
    "List[File]",
    "List[GATKPath]",
    "List[PicardHtsPath]",
    "List[Path]",
    "FeatureInput",
}

BOOLEAN_TYPES = {"boolean", "Boolean"}
INT_TYPES = {"int", "Integer", "long", "Long", "short", "Short", "byte", "Byte"}
FLOAT_TYPES = {"double", "Double", "float", "Float"}


def _bound(arg, key):
    raw = arg.get(key)
    if raw in (None, "-Infinity", "Infinity", "NaN"):
        return None
    try:
        return float(raw)
    except ValueError:
        return None


def _within_bounds(arg, value):
    low, high = _bound(arg, "min"), _bound(arg, "max")
    return (low is None or value >= low) and (high is None or value <= high)


def _numeric_domain(arg, numeric_policy):
    """Default, both bounds when finite, and one value just outside each finite bound.

    The out-of-bounds values are deliberate: an argument's rejection path is part of its behaviour,
    and the reference implementation's error text is output like any other. They are marked so a
    runner can choose to expect a non-zero exit.

    `numeric_policy` decides what happens to the common case of a numeric argument with a default
    and no declared bounds, which is 84 of HaplotypeCaller's 174 arguments:

    * `strict` gives it one value, so it is excluded from the array and held at its default. The
      claim is narrower and every value in it came from the inventory.
    * `perturb` also offers the default moved one step in each direction. These values are not
      declared anywhere, which is normally disqualifying, but a *differential* test does not need
      to know the right answer: the oracle defines it for whatever value is passed. What the value
      must be is *accepted*, and a small perturbation of a default is.

    Whichever is used is recorded in the report, because a coverage percentage means something
    different under each.
    """
    values = []
    default = arg.get("default")
    if default not in (None, "null", "[]"):
        values.append(("default", default))
        if numeric_policy == "perturb":
            as_int = arg["type"] in INT_TYPES
            try:
                base = float(default)
            except ValueError:
                base = None
            if base is not None:
                step = 1 if as_int else max(abs(base) * 0.1, 0.1)
                for label, delta in (("default_minus", -step), ("default_plus", step)):
                    moved = base + delta
                    if _within_bounds(arg, moved):
                        values.append((label, str(int(moved) if as_int else moved)))

    for edge, key, outside in (("min", "min", -1), ("max", "max", 1)):
        number = _bound(arg, key)
        if number is None:
            continue
        as_int = arg["type"] in INT_TYPES
        bound = int(number) if as_int else number
        values.append((edge, str(bound)))
        step = 1 if as_int else 1.0
        values.append((f"outside_{edge}", str(bound + outside * step)))

    return values


def domain_of(arg, fixtures, numeric_policy="strict", excluded_by_hand=None):
    """Return (values, exclusion). Exactly one of the two is meaningful.

    `values` is a list of (label, value) pairs; the label is what the coverage report prints, so a
    failing row can be read without decoding the value. `fixtures` maps an argument name to the
    values a repository is willing to pass for it, which is how paths get a domain.

    `excluded_by_hand` maps an argument to the reason a repository refuses to vary it. It exists
    for arguments that are real but not *about the tool*: `--help` and `--version` print usage and
    exit, so a row carrying them measures the argument parser. Refusing them is a judgement, so it
    is written down next to the fixtures rather than hard-coded here.
    """
    name = arg["name"]
    if excluded_by_hand and name in excluded_by_hand:
        return [], f"excluded by the repository: {excluded_by_hand[name]}"
    if name in fixtures:
        # An EMPTY list is a declaration, not an oversight: it says this repository must not give
        # this tool this argument at all. Picard has tools whose arguments are mutually exclusive
        # by another argument's value -- `FilterSamReads` refuses a `READ_LIST_FILE` unless
        # `FILTER` names a read-list filter, `MarkDuplicatesWithMateCigar` refuses `ASSUME_SORTED`
        # beside `ASSUME_SORT_ORDER` -- and a row carries every argument the domain holds, so one
        # such argument refuses every row of the array. Declaring the empty list drops it from the
        # rows and records why, which is the same bargain the exclusions below make.
        if not fixtures[name]:
            return [], "excluded by the repository: no value is declared for this tool"
        return [(f"fixture[{i}]", v) for i, v in enumerate(fixtures[name])], None

    if arg.get("deprecated"):
        return [], "deprecated: exercising it would test a path upstream has withdrawn"

    kind = arg["type"]

    if kind in BOOLEAN_TYPES:
        return [("true", "true"), ("false", "false")], None

    if arg.get("enum_members"):
        return [(m, m) for m in arg["enum_members"]], None

    if kind in INT_TYPES or kind in FLOAT_TYPES:
        values = _numeric_domain(arg, numeric_policy)
        if not values:
            return [], (
                "numeric with no declared bounds and no default: any value would be invented, "
                "and an invented value covers nothing"
            )
        return values, None

    if kind in PATH_TYPES:
        return [], f"{kind} needs a fixture; declare one in the repository's fixtures file"

    if kind in ("String", "List[String]"):
        default = arg.get("default")
        if default in (None, "null", "[]"):
            return [], "free-form string with no default: a fixture must supply the values"
        return [("default", default)], None

    return [], f"unmodelled type {kind}"


def load_fixtures(raw):
    """Accept either {"values": ..., "exclude": ..., "per_tool": ...} or a bare mapping.

    Returns `(values, exclude, per_tool)`. `values` is what every tool sees for an argument name;
    `per_tool` maps a tool name to the values that *replace* them for that tool.

    Why the third: `--INPUT` is one name over tools that do not read the same kind of file at all.
    A shared list of BAMs gives `BedToIntervalList` an array whose every row the reference refuses
    for the same reason, which measures the fixture rather than the tool. The override is per tool
    and not per type because a type is not enough either: two tools taking a `File` can want a BED
    and an interval list.
    """
    if isinstance(raw, dict) and ("values" in raw or "exclude" in raw or "per_tool" in raw):
        per_tool = {
            name: entry
            for name, entry in (raw.get("per_tool") or {}).items()
            if not name.startswith("$")
        }
        return raw.get("values", {}), raw.get("exclude", {}), per_tool
    return raw or {}, {}, {}


def for_tool(values, per_tool, name):
    """The fixture values a named tool sees: the shared ones, with its own overriding.

    An override replaces an argument's whole list rather than extending it, because the point of
    declaring one is that the shared value is *wrong* for this tool, not that it is incomplete.
    """
    own = per_tool.get(name)
    if not own:
        return values
    merged = dict(values)
    merged.update({key: entry for key, entry in own.items() if not key.startswith("$")})
    return merged


def build(tool, fixtures, numeric_policy="strict", excluded_by_hand=None):
    """Split a tool's arguments into a covering-array domain and a list of exclusions."""
    domain, excluded = {}, []
    for arg in tool["arguments"]:
        values, why = domain_of(arg, fixtures, numeric_policy, excluded_by_hand)
        if why is not None or len(values) < 2:
            # A one-value argument cannot participate in an interaction: it is held at that value
            # for every row, which is coverage of one level, not of a pair.
            # `source` is what a runner needs: a fixture value must be passed on the command line
            # (the tool has no default for it), while an argument held at its own default is
            # simply omitted, which is what "default" means. Re-serializing a default is how the
            # first run of CollectAlignmentSummaryMetrics rejected all sixteen rows: the default of
            # a List[String] came back as one token, brackets and commas included.
            source = "none"
            if values:
                source = "fixture" if values[0][0].startswith("fixture[") else "default"
            excluded.append(
                {
                    "argument": arg["name"],
                    "type": arg["type"],
                    "reason": why or f"only one value available ({values[0][0]})",
                    "held_at": values[0][1] if values else None,
                    "source": source,
                }
            )
            continue
        domain[arg["name"]] = values
    return domain, excluded

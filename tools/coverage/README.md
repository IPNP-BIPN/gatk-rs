# Covering arrays

Generates t-wise covering arrays over a tool's arguments, from
`tools/inventory/generated/inventory.json`. See
[docs/what-pairwise-coverage-costs.md](../../docs/what-pairwise-coverage-costs.md) for the measured
sizing of the whole programme and for why the two value policies exist.

```sh
python3 tools/coverage/covering.py --tool HaplotypeCaller --t 2 --verify
python3 tools/coverage/covering.py --tool HaplotypeCaller --t 3 --json hc-t3.json
python3 tools/coverage/covering.py --all --t 2
```

| file | what it does |
|---|---|
| `covering.py` | IPOG generation, exhaustive verification, per-tool and whole-inventory reports |
| `domains.py` | turns an argument's declared schema into its value domain, or excludes it with a reason |

## Fixtures

Path-typed arguments (`--INPUT`, `--REFERENCE_SEQUENCE`, ...) have no domain a generator can
invent: the value must be a file that exists and holds content the tool accepts. Supply them with

```sh
python3 tools/coverage/covering.py --tool CollectQualityYieldMetrics --t 2 \
  --fixtures picard-rs/tools/coverage/fixtures.json
```

where the file maps an argument name to the values that repository is willing to pass:

```json
{
  "--INPUT": ["tests/data/small.bam", "tests/data/unmapped.bam"],
  "--OUTPUT": ["/tmp/out.txt"]
}
```

Fixtures belong to the repository that runs the array, not to the inventory: `gatk-rs` cannot know
which BAMs `picard-rs` keeps. Every argument left without one is reported by name, so the gap is
visible rather than silent.

`--INPUT` is one argument name over tools that do not read the same kind of file, so a single
shared list is wrong for some of them: giving `BedToIntervalList` a BAM produces an array whose
every row the reference refuses for the same reason, which measures the fixture and not the tool.
A tool may therefore declare its own values, which **replace** the shared ones for that tool only:

```json
{
  "values": { "--INPUT": ["/work/fixtures/small.bam"] },
  "per_tool": {
    "BedToIntervalList": { "--INPUT": ["/work/fixtures/targets.bed"] }
  }
}
```

The override is per tool rather than per type because a type is not enough either: two tools
taking a `File` can want a BED and an interval list. It replaces the list rather than extending it,
because the reason to declare one is that the shared value is wrong here, not that it is
incomplete.

## The rule this code exists to enforce

A covering array over invented values proves nothing, and fails silently, because the array still
looks covered. So a value is used only when it traces to the inventory (`enum_members`, `min`,
`max`, `default`) or to a declared fixture, and everything else is **excluded by name with a
reason**. The exclusion list is part of the coverage claim: a tool at "t=2 covered" with half its
arguments excluded is covered over half a tool.

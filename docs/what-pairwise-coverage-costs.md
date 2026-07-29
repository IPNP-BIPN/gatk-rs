# What pairwise coverage actually costs

The plan commits to t-wise covering arrays over every tool's arguments, t=2 everywhere and t=3 on
the critical path, and sizes the programme without ever computing what that means in runs. Now it
is computed: `tools/coverage/covering.py` generates the arrays from the same
`tools/inventory/generated/inventory.json` that gives the 311 tools and 13,130 arguments, and
verifies each one exhaustively.

## The numbers

Every tool in the inventory, t=2, arrays verified tuple by tuple:

| value policy | oracle runs | arguments in the arrays | excluded |
|---|---:|---:|---:|
| `strict` | **19,437** | 6,312 of 13,130 (48%) | 6,818 |
| `perturb` | **21,918** | 8,674 of 13,130 (66%) | 4,456 |

Median tool: **19 rows**. The distribution is very long-tailed: `VCFComparator` alone needs 3,102,
because `--read-filter` has 56 members and `--annotation` has 55, and pairwise over two domains
that size cannot cost less than 56 x 55 = 3,080 rows. The array comes within 1% of that floor,
which is the useful check on the generator: IPOG is not doing anything clever there, and nothing
can.

For scale at the other end:

| tool | arguments | rows (t=2) | rows (t=3) |
|---|---:|---:|---:|
| `CollectQualityYieldMetrics` | 23 | 16 | — |
| `HaplotypeCaller`, strict | 174 | 62 | 325 |
| `HaplotypeCaller`, perturb | 174 | 83 | — |

**174 arguments, 62 runs for every pair.** That is the result that makes the programme's coverage
commitment tractable, and it is why the plan says covering arrays rather than "test the
parameters".

## The two policies, and why there are two

Most arguments carry their value domain in the inventory: a boolean has two values, an enum has its
members, a bounded numeric has its bounds. The awkward case is the numeric with a default and no
declared bounds, which is 84 of HaplotypeCaller's 174 arguments.

* `strict` gives it one value, so it is held at its default and excluded from the array. Every
  value used came from the reference's own documentation.
* `perturb` also offers the default moved one step each way. Those values are not declared
  anywhere. In a differential test that is acceptable, and this is the one place worth being
  precise about why: the test does not need to know the *right* output for a value, because the
  oracle defines it. The value only has to be *accepted*. A small perturbation of a default is.

The policy is recorded in every report, because a coverage percentage means a different thing under
each. Neither is the default answer for the programme; the choice is per tool, and it belongs in
the tool's conformance manifest entry.

## Constraints, and what they were worth

Some combinations are invalid: `CollectQualityYieldMetrics` refuses `FLOW_MODE=true` outright
("obsolete. Flow support now provided by CollectQualityYieldMetricsFlow"), and the coordinate
collectors refuse a queryname-sorted input under either value of `ASSUME_SORTED`. An array that
ignores this spends most of its rows being rejected, and a coverage figure computed over rows that
could never produce output overstates itself.

Constraints are declared in the fixtures file (`forbid`, with a `why`, optionally scoped to
`tools`), compiled into forbidden tuples, and honoured during generation. `--verify` counts them
separately: a forbidden tuple is never reported missing, because coverage that does not exist is
not coverage that was skipped. Measured on the first two tools, running each row against the
oracle:

| tool | before | after |
|---|---|---|
| `CollectQualityYieldMetrics` | 11 rows, **3 accepted**, 3 distinct outputs | 10 rows, **10 accepted**, 4 distinct outputs |
| `CollectAlignmentSummaryMetrics` | 16 rows, **12 accepted**, 9 distinct outputs | 16 rows, **16 accepted**, 11 distinct outputs |

Every row now runs, and both tools produce more distinct outputs than before, which is the point:
rows spent on rejections are rows not spent on the tool. The rejections themselves are behaviour
the port owes, and belong in the conformance manifest as their own cases.

Scoping matters and was learned by measuring. The queryname clause started as two argument-pair
constraints (`with ASSUME_SORTED=false`, `with CREATE_INDEX=true`) and neither shape was right: the
input is invalid for those collectors under any combination, while being an ordinary input for
`SortSam` or `SamToFastq`. A global clause would have deleted a dimension those tools need.

## What is still excluded, and what it would take

Under `strict`, 52% of arguments are outside the arrays. The breakdown is not exotic:

* **paths** (`File`, `GATKPath`, `PicardHtsPath`, and their list forms) have no domain a generator
  can invent. A path must exist and hold content the tool accepts, so it takes a *fixture*, and
  fixtures are owned by the repository that runs the array, not by the inventory. Supplying them is
  the single biggest lever on the excluded count.
* **free-form strings with no default**, same reason.
* **unbounded numerics**, under `strict` only.

Every exclusion is emitted by name with its reason, in the report and on the console. That is the
part that must not be lost: a tool reported as "t=2 covered" with half its arguments excluded is
covered over half a tool, and a coverage number that hides which half is worse than no number.

## Reproducing

```sh
python3 tools/coverage/covering.py --tool HaplotypeCaller --t 2 --verify
python3 tools/coverage/covering.py --tool HaplotypeCaller --t 2 --json hc.json
python3 tools/coverage/covering.py --all --t 2                     # 30s, the table above
python3 tools/coverage/covering.py --all --t 2 --numeric-policy perturb
```

`--verify` is exhaustive rather than sampled: it enumerates every t-way tuple of the domain and
asserts the array covers it. The generator is deterministic, so a regenerated array is diffable and
a change in the arrays means a change in the inventory.

## What this does not yet do

The arrays are rows of argument assignments. Running them against the oracle and the port, and
comparing the outputs, is the next piece (the plan's 0.3), and it needs the fixtures above before
any tool with an `--INPUT` can be exercised at all. Until then these numbers size the work; they do
not perform it.

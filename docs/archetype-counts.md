# Archetype counts: what the inventory says, and one thing the plan got wrong

The program's feasibility rests on archetype amortisation, and therefore on how many tools each
archetype actually holds. This records the measured counts and one correction, so the numbers
that size Phases 2 through 5 come from the generator rather than from a table someone typed.

## The correction

The plan's group table has a row reading **"57 | Metrics | metrics collector"**, presented as
57 tools sharing one shape. That is a conflation, and it is the *same* conflation the inventory
already caught once when the tool count was corrected from 375 to 311.

`tools/inventory/generated/inventory.json` keeps the two apart:

```json
"tools": 311,
"metrics_definitions": 57,
```

**57 is the count of metric definition classes** — `InsertSizeMetrics`, `AlignmentSummaryMetrics`
and friends, which are output-file schemas with no command line. They are not tools and cannot
be ported "one after another" because there is nothing to invoke.

There is also no group named "Metrics" in the inventory at all. The nearest group is
**Diagnostics and Quality Control, 56 tools**, which is a different set: it includes tools that
write no metrics file, and excludes metrics-writing tools filed under other groups.

## What the archetype actually holds

Measured independently in picard-rs by `tools/stratify/stratify.py`, which classifies every
Picard tool that writes a `MetricsFile`:

| | |
|---|---:|
| Picard tools writing a `MetricsFile` | **44** |
| strata, by porting machinery | **11** |
| tools in the largest stratum | **6** |

See `picard-rs/docs/decisions/0001-the-metrics-archetype-is-not-homogeneous.md` for the
stratification and why environment requirements are separated from porting machinery.

## Why this matters to the sizing

"57 collectors, one shape" and "44 tools across 11 strata, the largest holding 6" imply very
different amortisation. The first says a single up-front cost is repaid 56 times; the second
says it is repaid at most 5 times per stratum, and 11 up-front costs must be paid.

That does not settle whether the program is feasible. It does mean the delta must be measured
**within a stratum**, and weighted by that stratum's size, before any phase after Phase 1 is
sized. The delta remains unmeasured.

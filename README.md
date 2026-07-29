# gatk-rs

Native Rust reimplementation of the Broad Institute's
[GATK](https://github.com/broadinstitute/gatk), targeting **byte-identical** output against a
pinned reference build. Work in progress.

> **This is not the official GATK.** It is an independent reimplementation, not affiliated
> with or endorsed by the Broad Institute.

## Reference version

Ported from GATK `4.6.2.0`, commit `76edc75c26504da94bbaee66584e107e76ee15de`, which pins:

| Dependency | Version |
|---|---|
| Picard | 3.4.0 |
| htsjdk | 4.2.0 |
| Barclay | 5.0.0 |
| Intel GKL | 0.8.11 |

All three ports in this program use those exact pins, so they are mutually coherent by
construction.

## Origin

This grows out of [broadinstitute/gatk#9384](https://github.com/broadinstitute/gatk/pull/9384)
("Native Apple Silicon (arm64) support for GATK"), which established that GATK's
floating-point callers are **bio-identical but not bit-identical** across CPU architectures:
Java `Math.log/exp/pow` differ by roughly 1 ULP between architectures, and only `StrictMath`
is portable.

That finding located precisely where the reference implementation stops being a property of
the algorithm and starts being a property of the JVM and the host. This project makes that
boundary explicit, testable, and eliminated.

## Method

Every feature branch ports a **named symbol** of the pinned reference source, read from the
pinned clone and translated. Behavior is never reconstructed from documentation or memory.

```
main                                  shared infrastructure only
└── tool/gatk-<toolname>               one per tool
    └── feat/gatk-<toolname>-<feature> one per ported symbol
```

The tool inventory is **generated**, not hand-written: 375 tools and 10,796 arguments are
derived mechanically from the reference's own machine-readable tool documentation, along with
argument schemas, branch names, and the differential test matrix.

## What exists today

| layer | state |
|---|---|
| tool inventory | generated: 311 tools, 13,130 arguments, from the pinned reference's own CLI |
| status dashboard | generated from the inventory and the ports' manifests: [docs/STATUS.md](docs/STATUS.md) |
| covering arrays | generated and verified per tool: [what pairwise coverage costs](docs/what-pairwise-coverage-costs.md) |
| oracle image | digest-pinned `linux/amd64`, GATK 4.6.2.0, probe asserts the contract during the build |
| `gatk-readfilter` | 55 of the 56 read filters, oracle-backed: 79 instances over 59 records, 4,661 decisions identical to the reference. The exception is `JexlExpressionReadTagValueFilter`, which needs a JEXL expression engine |
| `gatk-engine` | intervals: parsing against a sequence dictionary, union and merge, overlap detection |

The read filters come first because they are stateless, touch no floating point, and every tool
that reads reads runs a chain of them. A wrong filter does not produce a wrong number, it produces
a different set of reads, and every number downstream inherits that.

## Coverage

"Every parameter" is defined operationally, because exhaustive does not exist here:
HaplotypeCaller alone has 174 arguments, implying on the order of 2^174 combinations.

1. **t-wise covering arrays**, generated per tool. `t=2` everywhere, `t=3` on the critical
   path. Combinatorial interaction testing is the established standard for this problem.
2. **Coverage-guided differential fuzzing** against the instrumented Java oracle, steering
   toward reference branches the covering arrays never reach. Divergences are minimized and
   promoted to permanent regression cases.
3. **The reference's own test resources**, which encode the edge cases upstream cared about.

## Bit-identity contract

Goldens come from the pinned reference in a digest-pinned `linux/amd64` container on JDK 17,
produced only on real x86-64 CI. Emulated x86-64 on Apple Silicon does not expose AVX, so GKL
native paths can silently fail to load and yield goldens matching no real machine; the oracle
runner asserts the resolved provider state and fails rather than degrading.

Fields legitimately allowed to vary are canonicalized under explicitly declared rules, and
every comparison records what was compared raw versus canonicalized. Values that cannot be
matched exactly are quarantined with their measured divergence rate and reported as
**bio-identical** rather than **bit-identical**.

## Part of a three-repository program

| Repo | Ports | Depends on |
|---|---|---|
| `htsjdk-rs` | htsjdk 4.2.0 | (none) |
| `picard-rs` | Picard 3.4.0 | `htsjdk-rs` |
| `gatk-rs` | GATK 4.6.2.0 | `picard-rs`, `htsjdk-rs` |

The topology mirrors upstream. GATK's `Main.getPackageList()` returns
`["org.broadinstitute.hellbender", "picard"]`, which is why `gatk MarkDuplicates` runs Picard's
code, and why the dependency runs in that direction here too.

## Relationship to `fulcrumgenomics/riker`

[riker](https://github.com/fulcrumgenomics/riker) is an independent, MIT-licensed Rust
reimplementation of Picard's QC tools from the maintainers of Picard and htsjdk. It overlaps this
program's Picard layer only, not its GATK layer. The distinction that governs the relationship:
**riker targets functional equivalence, this program targets byte equivalence.** riker is the
better tool to *use*; this is a byte-for-byte reproduction of the existing one, bugs included.

That makes riker a source of divergence candidates and never a source to port from. Its `ERRATA`
documents exactly where a careful reimplementer chooses to differ from Picard, and every entry is
pinned as a conformance case in `picard-rs`, measured against the reference rather than trusted.
See `picard-rs` for the two entries pinned so far, and htsjdk-rs decision 0020 for a case where
byte comparison surfaced a behaviour (alignment-block cycle binning is not in riker's errata) that
careful reimplementation did not.

## Commit attribution

Commits are co-authored with the model that wrote them. On 2026-07-21 the history of all three
repositories was rewritten to add that trailer uniformly, at the maintainer's request, changing
every commit SHA. `gatk-rs` pins no dependencies yet, so nothing here was invalidated; the note is
recorded for symmetry with the other two repositories, whose SHA-pinned historical builds are no
longer bit-reproducible as a result. The trade was made deliberately.

## License

Apache License 2.0, matching GATK. See `LICENSE`.

Worth stating because it is easy to get wrong: GATK is **Apache 2.0**, not BSD-3-Clause.
Several third-party descriptions (including the Homebrew formula) say BSD-3-Clause; the
authoritative source is `LICENSE.TXT` in the pinned clone, which is Apache 2.0 with a Broad
preamble. GitHub's API reports `NOASSERTION` for the same reason. The three repositories in
this program therefore do not share one licence: `htsjdk-rs` and `picard-rs` are MIT,
`gatk-rs` is Apache 2.0.

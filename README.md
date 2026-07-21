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

## License

GATK is BSD-3-Clause. License for this port is to be finalized before the repository is made
public.

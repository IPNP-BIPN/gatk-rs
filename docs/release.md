# Release checklist

What has to be true before the first tag of gatk-rs, picard-rs and htsjdk-rs, and where each item
stands. This is the page the tag is decided on: keep it current, and tag only when every row is
green. A tag, a release and a publish are each the maintainer's explicit decision.

States: **green** done and verified; **amber** done but not yet verified, or nearly done; **red**
open. Last updated 2026-11-24.

## Blockers

| # | blocker | owner (session) | PR / issue | state |
|---|---|---|---|---|
| 1 | The chr20 30x best-practices chain passes at scale, port against GATK 4.7.0.0, x86-64 and arm64 | d1-d2-chr20-30x-benchmar; ci-gatk-vs-port-at-scale; fix-two-chr20-divergence | gatk-rs #1465 (draft), #1493, #1479, #1472 | **red**: all open. Known at scale: VariantRecalibrator's `##contig` lines (#1479), BaseRecalibrator on bgzipped known sites (#1472), CollectReadCounts dropping `@SQ` M5/UR (no PR yet) |
| 2 | `main` green in all three repositories, on the commit that is tagged | everyone | CI on `main` | **amber**: the last completed `main` runs are green (gatk-rs `4d2cb248`, picard-rs `2fbd7c48`, htsjdk-rs `cb1f5804`); the runs on the current heads (`43faae30`, `c7921ed9`, `95cb236a`) are queued |
| 3 | `measured.json` frozen from a CI run of the final tree (`tools/coverage/`, `tools/fuzz/`) | test-resources | gatk-rs #1478 (draft) | **red**: to be frozen after #1478's CI run |
| 4 | SAM and CRAM in every tool that reads or writes reads | sam-cram-everywhere | gatk-rs #1481 (merged, read walkers); the rest in progress | **red**: CRAM is still a port limitation in Picard-origin and Spark tools ([known-limits.md](known-limits.md)) |
| 5 | Milestone X closed | x-tools-close; x-spark-ml; b1-run-python-gcnv-htsge | gatk-rs #12; #1504 (X.6 R), #1510 (X.2 ML), #1466 and #1505 (gCNV) | **red**: Spark, ML, Python (#806) and R (#807) boxes open |
| 6 | Milestone GPU closed | gpu-milestone | gatk-rs #33; #1508 | **red**: #1508 open; eligibility guard, two-architecture determinism gate and most kernels unticked |
| 7 | #1454: the held-out arguments back in the arrays | heldout-a-l (working); heldout-m-z | gatk-rs #1454; #1462, #1470, #1471 (merged); #1469 (open) | **red**: issue open, VALIDATION.md still counts 596 held-out arguments over 119 tools |
| 8 | The fuzzer on every tool | fuzz-gatk; fuzz-picard | gatk-rs #1476 (draft); picard-rs #481 (merged) | **red**: gatk-rs's 18 skipped tools not merged |
| 9 | t=3 on the critical path | t-3-flowbased-likelihood; b2-b3-held-out-arguments | gatk-rs #1454, #1501 | **red**: #1501 open |
| 10 | Intel and Graviton measured: faster than Java on CPU and peak RSS, with no trade-off | unassigned | gatk-rs #1493 measures `ubuntu-latest` and `ubuntu-24.04-arm` (Neoverse N2) | **red**: no run on a named Intel part or on AWS Graviton |
| 11 | The reference's own instabilities recorded as limits, not quarantines | quarantines | gatk-rs #1477 | **amber**: open |

## Release engineering

| # | item | PR | state |
|---|---|---|---|
| R1 | Release workflow: `v*` tag → draft release, Linux x86-64, Linux arm64, macOS arm64, sha256, pinned 1.97.1, no `target-cpu=native`, smoke-tested on `qemu64` | gatk-rs #1511; picard-rs (release prep PR) | **amber**: dry run on the PR |
| R2 | One version per workspace (`[workspace.package] version`) | the same PRs, and htsjdk-rs's | **amber** |
| R3 | LICENSE, NOTICE and THIRD-PARTY-NOTICES cover GATK (Apache-2.0), Picard and htsjdk (MIT) and every ported third party; each archive carries the licence of every linked crate | the same PRs | **amber** |
| R4 | README (install, usage, what is byte-identical and how it is proven), CHANGELOG, [known-limits.md](known-limits.md) | the same PRs | **amber** |
| R5 | crates.io | none | **later**: htsjdk-rs, picard-rs, bwa-mem4 and fgkl are git dependencies pinned by revision, which binaries accept and crates.io does not |

## Open questions for the maintainer

- **One version across the three repositories?** picard-rs is at 0.2.0 and htsjdk-rs at 0.4.0,
  both already tagged and released; gatk-rs has never been tagged and is at 0.1.0. Each workspace
  now carries one version; aligning the three numbers is a decision, not a fix.
- **Licence of ported Apache-2.0 code in the MIT repositories.** picard-rs and htsjdk-rs port
  Apache Commons Math (Apache-2.0) into MIT trees. THIRD-PARTY-NOTICES now carries the Apache
  licence and Commons Math's NOTICE; whether those files should say the ported portions remain
  under Apache-2.0 is the maintainer's call.
- **glibc.** The Linux binaries are built on Ubuntu 24.04 and need glibc 2.39 or later. Building on
  22.04 (glibc 2.35) or against musl would reach older hosts.

## Cutting the release (once every blocker is green)

1. Merge the release-prep PRs; set the version in `Cargo.toml` and the date in `CHANGELOG.md`.
2. With the maintainer's OK: `git tag -s vX.Y.Z && git push origin vX.Y.Z`. The workflow builds,
   smoke-tests and attaches the archives to a **draft** release.
3. Check the three archives (sha256, `bin/gatk --version`, THIRD-PARTY-LICENSES.txt), then publish
   the draft by hand.

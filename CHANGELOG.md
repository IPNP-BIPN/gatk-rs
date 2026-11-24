# Changelog

All notable changes to this repository, newest first. Every crate here carries the workspace's one
version (`[workspace.package]` in `Cargo.toml`), because they are one port of one program.

The format is [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); this project does not yet
follow semantic versioning, because 0.x is what an unreleased port is.

## [0.1.0] - unreleased

The first release. Reference: **GATK 4.7.0.0** (`0cde69eed303`), with Picard 3.5.0 and htsjdk
5.0.0 underneath (htsjdk-rs pinned by revision).

### Added

- `gatk`, one binary taking the reference's own command line (`gatk <Tool> <args>`), for Linux
  x86-64, Linux arm64 and macOS arm64, each archive with its sha256 and the licence of every crate
  linked into it (`THIRD-PARTY-LICENSES.txt`).
- All 312 tools of the reference's inventory are oracle-backed: CI re-derives each tool's goldens
  from the pinned GATK container on real x86-64 and compares the port's output byte for byte
  ([docs/VALIDATION.md](docs/VALIDATION.md), [docs/STATUS.md](docs/STATUS.md)).
- What is not byte-identical, and why, is listed in [docs/known-limits.md](docs/known-limits.md).

### Not yet

- The crates are not on crates.io: htsjdk-rs, bwa-mem4 and fgkl are git dependencies
  pinned by revision, which a binary is fine with and crates.io is not. Publishing comes later.

# hdf5-pure, the port's HDF5 reader

Issue #1294, then W6 of the plan to 100%.

**Verdict: hdf5-pure reads every HDF5 file the port is handed. `hdf5_writer.rs` stays in-house.**
hdf5-pure 0.48.0 is a normal dependency of `gatk-tools`, and `gatk_tools::hdf5_reader` is a thin layer
over it. There is no C library anywhere: no `hdf5` or `hdf5-sys` crate, no libhdf5. The three
`PortLimitation` refusals that HDF5 used to cause are gone:

| Tool | What it reads now |
| --- | --- |
| `CreateReadCountPanelOfNormals` | HDF5 read-counts files, beside TSV ones |
| `DetermineGermlineContigPloidy` | HDF5 read-counts files, beside TSV ones |
| `DenoiseReadCounts` | HDF5 read-counts files, and the `--count-panel-of-normals` panel, which turns on the SVD denoising |

`FilterIntervals` read its counts as text and so would have failed on an HDF5 file. It now reads
them through the same function. The VETS tools (`TrainVariantAnnotationsModel`,
`ScoreVariantAnnotations`) read their annotation files through hdf5-pure as well.

## How a counts file is read

`SimpleCountCollection.read` chooses its reader from the file's first 8 bytes (`IOUtils.isHDF5File`),
not from its name, and `gatk_tools::collect_read_counts::read_text` does the same. An HDF5 file is
read from the four paths `HDF5SimpleCountCollection` writes. It is then handed back as the TSV text
`SimpleCountCollection.write` would make of the same collection: the stored dictionary, the
sample's read group, and each count cast with Java's `(int)`. Every runner keeps its one TSV parser,
and an HDF5 input cannot be parsed differently from its TSV twin.

## Why hdf5-pure and not the hand-written reader

Before this, `hdf5_reader.rs` was a reader of about 400 lines, written from the specification. It
handled only the two layouts it had seen: the C library's default (version 0 superblock, symbol
table groups) and the port's own writer's. hdf5-pure reads the whole format, checks each object
header's `lookup3` checksum where the layout has one, and was already pinned and audited here as a
test dependency. The hand-written reader is not deleted. It moved to
`crates/gatk-tools/tests/hdf5_spec_reader/` and remains the second, independent reader in
`hdf5_second_reader.rs`, so the writer is still checked by two readings of the specification that
share no code.

What it was checked against: the reference's own `CollectReadCounts --format HDF5`,
`SimpleCountCollection.writeHDF5` and `CreateReadCountPanelOfNormals` outputs (the covering-array
fixtures `counts.hdf5`, `pon_s*.counts.hdf5`, `pon.hdf5`, `pon_gc.hdf5`, `pon_one.hdf5`). It read every
dataset of every file: fixed-length strings, double matrices, and the chunked matrices
`HDF5Utils.writeChunkedDoubleMatrix` splits into `chunk_<i>` datasets.

## Features

`default-features = false, features = ["std", "checksum"]`. Deflate and chunked layout are left out
because no file the reference writes needs them. The Broad's `HDF5File` creates every dataset
contiguous and unfiltered, and its "chunks" are separate datasets, not HDF5 chunked layout. If a
reference file ever needs a filter, hdf5-pure refuses it with an error, and the error is read as
absent. That fails a covering-array row instead of producing a guess.

## Why the writer stays

- Every HDF5 file the port writes holds only doubles and fixed-length strings. The in-house writer
  is about 330 lines and already measured by covering arrays.
- The reference stamps its object headers with a time, so files are compared as trees, not bytes.
  hdf5-pure's writer would get no closer to the reference's bytes, and a 0.x crate on the byte path
  could change those bytes at any minor release. The reader is off the byte path: a different
  reading of the same tree is a failure, not a drift.

## Pinning and licences

`hdf5-pure = { version = "=0.48.0", ... }` together with `Cargo.lock`. hdf5-pure's own manifest names
its four `-core`, `-filter`, `-format` and `-space` crates with `=`, so the one pin fixes all five.
The rust-version is 1.89, below the pinned toolchain. The crates are MIT OR Apache-2.0.
`hdf5-pure-filter` always compiles a port of the HDF5 library's scale-offset filter, which The HDF
Group licenses under the 3-clause BSD licence. That notice is reproduced in the repository's
`NOTICE`, because the licence asks for it in binary redistributions. ZFP, the other BSD-licensed
code in the crate, is behind a feature this build does not turn on.

Upgrading hdf5-pure is now a change on the read path of four tools, so it goes through the arrays
that read HDF5 (`CreateReadCountPanelOfNormals`, `DenoiseReadCounts`, `FilterIntervals`, the VETS
tools).

# hdf5-pure as a second reader of the port's HDF5

Issue #1294.

**Verdict: adopt as a test-only reader. `hdf5_writer.rs` stays in-house.** hdf5-pure 0.48.0 is a
dev-dependency of `gatk-tools` and reads back every shape the writer produces in `cargo test`, with no
Docker. Nothing in the port's dependencies changes. Replacing the writer remains what the issue said:
an option for the day a tool needs chunked layout, deflate, integer or compound types, or dense group
storage. No tool needs any of those today.

## What was tried

`crates/gatk-tools/tests/hdf5_second_reader.rs` builds a tree with hdf5-pure, starting at the root
(groups, datasets, datatype, shape, and every value compared bit for bit). The tree must be the one
that was written, and it must agree with what `gatk_tools::hdf5_reader` reads from the same bytes.

| Test | What it writes | Result |
| --- | --- | --- |
| `every_shape_the_writer_has` | 48 datasets in 9 groups. Doubles: ±0.0, NaN, ±∞, `MIN_POSITIVE`, the smallest subnormal, `MAX`. Also an empty array, a 2×3 matrix and a 1×n row. Fixed-length strings: `""`, a tab and a newline. A group 6 levels deep, and one with 40 children | Same tree and bits from both readers |
| `collect_read_counts_file` | `collect_read_counts::write_hdf5`, the function the tool calls | Same tree and bits; `/counts/values` is `[1, 3]` |
| `a_header_whose_checksum_is_wrong_is_refused` | The root's link `values` renamed `walues`, checksum left alone | hdf5-pure refuses it with a checksum error |

The third test is the control. It shows that hdf5-pure checks each object header's `lookup3`
checksum, so the files the first two tests open have been checked, not just tolerated. A first
control, which flipped a structural byte, failed with an EOF before reaching the checksum. That is
why the control renames the link instead.

The 3 tests take less than 0.01 s. hdf5-pure is built with `default-features = false,
features = ["std", "checksum"]`, which leaves out deflate. The writer never compresses, and leaving
deflate out means the dev-dependency adds 5 packages to `Cargo.lock` (hdf5-pure and its four
`-core`/`-filter`/`-format`/`-space` crates, all pinned with `=` to the same version) and no other
third-party code. The licence is MIT OR Apache-2.0 and the rust-version is 1.89, below the pinned 1.97.1.

## What it caught

Nothing: the two readers and the written tree agree everywhere. Before this, the port's HDF5 was read
back only by h5py in the covering-array harness, which needs Docker. Now a writer change that
produces a file the specification does not allow fails `cargo test`.

## Why the writer stays

- Every HDF5 file the port writes holds only doubles and fixed-length strings. The in-house writer
  is about 330 lines and already measured by two covering arrays.
- hdf5-pure is a 0.x crate that moves quickly: 0.47.0 when the issue was written, 0.48.0 now. On the
  byte path, each minor release could break the API and change the bytes. Off the byte path, a break
  only means updating a test.
- The reference stamps its object headers with a time, so the files are compared as trees, not as
  bytes. A different writer would not get any closer to the reference's bytes.

## Pinning

`hdf5-pure = { version = "0.48.0", ... }` together with `Cargo.lock`. Upgrading it is a test-only
change. If an upgrade starts refusing a file, read that as a finding about the writer or about
hdf5-pure, not as noise.

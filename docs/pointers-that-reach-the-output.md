# Pointers that reach the output

`BwaMemIndexImageCreator` writes a file whose bytes are not a function of its input. Building one
reference twice in a single process gives two images of the same length that differ in nine bytes,
and those bytes read `70 1a 79 f8 ff 7f 00 00`: an address on the JVM's heap, written into the file
by BWA's own serialiser through JNI.

That is a different animal from the two hazards already written down here:

  * [an unspecified order that reaches the output](an-unspecified-order-that-reaches-the-output.md)
    is deterministic for a given JDK and a given key set; it is unspecified, not unstable;
  * a log4j timestamp is stripped by the dump because it is a clock, and a clock is not an answer.

A pointer cannot be stripped by masking, because which bytes carry one is not knowable from the
file alone, and it cannot be pinned by a golden, because address-space layout randomisation moves
it on every run. A golden holding those bytes would go red on a run that changed nothing.

## What the suite holds instead

  * the SIZE, which is a function of the reference: 1333 bytes for five repeats of eight bases,
    1551 for twenty, 1334 for the same sequence under a longer contig name;
  * the file the tool wrote, which is what `--OUTPUT` and its default decide;
  * the refusal, when the reference cannot be read;
  * and the instability itself, as a `stable` row that reads `false` for every case, so a future
    run that starts producing identical images fails this suite rather than passing it quietly.

## What that means for the port

The port builds the image itself (`crates/gatk-tools/src/bwa_index_image.rs`): the same packing,
the same random replacement of ambiguous bases, the same transform, occurrence counts and sampled
suffix array, laid out the way `bwa_idx2mem` lays them out, with zero where BWA writes an address.

Masking turned out to be possible after all, but not from the bytes alone: it needs the layout.
`bwt_size` and `n_sa` sit at fixed offsets and place the `bntseq_t` block, whose `n_seqs` and
`n_holes` place everything after it, so the address fields (the `sa` pointer of `bwt_t`, the
`anns`, `ambs` and `fp_pac` pointers of `bntseq_t`, and the `name` and `anno` pointers of each
contig) and the structures' padding are known exactly. `tools/coverage/run_array.py` zeroes those
and compares every other byte, which is what the covering array measures. The golden above still
holds the size and the instability, because the conformance dump reads the image without the
layout.

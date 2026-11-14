# What a GenomicsDB query hands the callers

`GenotypeGVCFs`, `GnarlyGenotyper` and `SelectVariants` take `-V gendb://workspace`. GATK 4.7.0.0
reads it with `GenomicsDBFeatureReader` (GenomicsDB 1.5.5), which does not read the workspace in
Java at all: it asks the native library for a stream of VCF **text** and parses that with htsjdk's
`VCFCodec`. So the port reproduces the text, and every caller downstream of it is the code that
already reads a VCF.

The text is three native stages composed:

1. `VariantQueryProcessor::scan_and_operate` walks the array's cells (one per call begin, and a
   copy at each block's END) and cuts them into the intervals over which no sample's call changes,
   with spanning deletions handled;
2. `BroadCombinedGVCFOperator::operate` merges the calls of each interval into one record;
3. htslib's `vcf_format` prints it, floats through `kputd`.

`crates/gatk-tools/src/genomicsdb_query.rs` is those three, and
`crates/gatk-cli/src/runners_genomicsdb_query.rs` the startup around them. Both were written
against the C++ sources of GenomicsDB 1.5.5, htslib and TileDB, and checked against the stream the
oracle produces (a `GenomicsDBQueryStream` read directly from Java, byte for byte, over the whole
array and over a sub-interval) before any caller was run on it.

## How the text differs from a CombineGVCFs output

Measured on the shared corpus's `cg_a`, `cg_b` and `cg_c`: imported, then exported, against
`cg_combined.g.vcf` (the same three through CombineGVCFs).

| | GenomicsDB | CombineGVCFs |
|---|---|---|
| allele order under a spanning deletion | `C *,T,<NON_REF>`: the `*` of the deletion that began earlier comes first | `C T,*,<NON_REF>` |
| ID | the union over every call overlapping the position, so a deletion's ID reaches the positions it spans (`rs102;rs100`) | the record's own |
| ID order | libstdc++'s `unordered_set<string>` iteration order | — |
| allele-specific raw annotations | kept and combined (element-wise sums, or histograms) | dropped by default |
| FORMAT order | `DP` last (`GT:AD:GQ:MIN_DP:PL:SB:DP`) | header order |
| a sample with no call | `.` in every field, GT included | `./.` |
| floats | float32, printed by htslib (`ExcessHet=3.0103`) | the input's text |
| `-L` | each query clips the blocks to its interval | records overlapping the interval |
| FILTER | always `.`: GATK does not ask for it | the inputs' |

Downstream, the callers differ from their CombineGVCFs runs only where the input did: the allele
order at the spanning-deletion site and the empty sample's GT. `GnarlyGenotyper`'s body was
identical.

## Three things that are not in the code's shape

* **Iteration order of unordered containers reaches the output.** The ID column is a
  `std::unordered_set<std::string>` and the histogram annotations (`AS_RAW_*RankSum`) come out of
  an `std::unordered_map` keyed by vid-field index. The port reproduces libstdc++'s (GCC 7)
  hashtable: `_Hash_bytes` (64-bit MurmurHash2, seed `0xc70f6907`), the prime rehash policy, and
  insertion at bucket begin. See also
  [an-unspecified-order-that-reaches-the-output.md](an-unspecified-order-that-reaches-the-output.md).
* **A spanning deletion is rewritten in place.** At the first position after a deletion's begin,
  `handle_deletions` drops the call's INFO, reduces its alleles to `*` (and `<NON_REF>`), picks the
  deletion allele with the lowest homozygous PL, and remaps PL, AD and GT to the reduced alleles.
  The rewritten call is what every later position of the span sees.
* **The arrays are listed by name.** Without `-L`, the reader walks the array directories sorted
  by contig NAME and then start (`ChrArrayFolderComparator`), not in dictionary order.

## The export options GATK sets

`GATKGenomicsDBUtils.createExportConfiguration` asks for GT only when the tool calls genotypes,
never for FILTER, and passes two caps. Each tool overrides `getGenomicsDBOptions`:

| tool | GT | max diploid alternates | max genotype count |
|---|---|---|---|
| GenotypeGVCFs | always | `--genomicsdb-max-alternate-alleles` (50) | `--max-genotype-count` (1024) |
| GnarlyGenotyper | always | `--max-alternate-alleles` + 1 | `--max-genotype-count` |
| SelectVariants | `--call-genotypes` | `--genomicsdb-max-alternate-alleles` | 1024 |

The options refuse a maximum that is not one more than `--max-alternate-alleles` (6 for
SelectVariants, whose genotype arguments are the defaults), with a `BadInput`.

## What the port refuses, in its own words

* `--genomicsdb-use-bcf-codec`: the stream is then BCF read by `BCF2Codec`, whose records were
  measured to differ from the VCF codec's and were not reproduced.
* a workspace that is not on a local file system (`gendb.gs://`, `gendb.s3://`, ...);
* an array with more than one fragment, which is what several import batches write: the native
  read merges fragments, and that merge was not measured;
* a vid map that combines QUAL, an allele-dependent INFO field combined by sum, mean, median or
  concatenation, or a multi-dimensional FORMAT field: GATK's own vid maps have none, and the
  native operator's behaviour on them was not measured;
* a record merging more than 199 distinct IDs, past the part of libstdc++'s prime list the port
  carries.

## How it is checked

* the `genomicsdb-query` conformance suite (`GenomicsDBQueryDump`): the corpus imported in the
  run, then 19 caller command lines, four of them refusals GATK makes (no reference, no workspace,
  and too few alternates for two tools), replayed by `crates/gatk-cli/tests/genomicsdb_query_end_to_end.rs` through the
  port's own `GenomicsDBImport`;
* the covering arrays of the three tools, which take `gendb:///work/fixtures/cg_ws` (the same
  three samples, imported by `MakeFixtures.java`) as one of their `--variant` values. The pair of
  that value with `--genomicsdb-use-bcf-codec=true` is forbidden in `tools/coverage/fixtures.json`,
  since every such row would measure the same refusal;
* unit tests of `kputd`, the libstdc++ iteration order (both measured cases), and the genotype
  index order.

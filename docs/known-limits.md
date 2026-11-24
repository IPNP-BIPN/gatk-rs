# Known limits

What a release of gatk-rs does **not** claim, in one place. Everything else that
[VALIDATION.md](VALIDATION.md) lists as oracle-backed is compared byte for byte against GATK
4.7.0.0 in CI. A limit is one of four kinds:

1. **A port limitation**: the port refuses, by name, something it has not ported. It never answers
   with different bytes instead.
2. **An exclusion of the reference's own tests**: an invocation of GATK's `src/test/resources`
   suite that the harness cannot compare, and why.
3. **An instability of the reference**: GATK itself gives different bytes from run to run or from
   host to host, so there is no single byte to match.
4. **A deliberate divergence**: where GATK 4.7.0.0's htsjdk 5.0.0 loses or corrupts data and a
   later htsjdk release fixes it, the port takes the fix and differs from the reference on purpose.

Counts are as of 2026-11-24 (`origin/main` `43faae30`). Two of the sources below, PR #1477 and
PR #1478, are not merged yet; their lists may still move.

## 1. Port limitations

A refusal is the exception `gatk_rs::PortLimitation`
([main_entry.rs](../crates/gatk-tools/src/main_entry.rs), `PORT_LIMITATION`): the `::` makes it
impossible to mistake for a Java class, and its message ends "This message is the port's own and
not GATK's." The status is non-zero, so a pipeline stops rather than continuing on a guess.

There are 107 construction sites over 20 files, and about 45 more reasons reach it through helpers
(`crates/gatk-cli/src/runners.rs` alone holds 65). By kind:

| kind | what is refused |
|---|---|
| formats | CRAM in and out of Picard-origin and Spark tools (being closed: see the release checklist); BCF in VcfFormatConverter; HDF5 counts and panels of normals; block-compressed and `.bci` SV evidence; Funcotator `vcf`/`simpleXSV` data sources; an indexless `--features-file` |
| engine features | some GATK read filters; `-L unmapped`; more than one `--input` where GATK merges them; block gathering in GatherVcfs; a Spark cluster master; Spark sort sampling over several partitions; PrintFileDiagnostics analyzers |
| second algorithms | `--pedigree` family priors (CalculateGenotypePosteriors, VariantEval, VariantAnnotator); GC-bias correction in DenoiseReadCounts; the SOMATIC splitter and `--keep-original-ac` in LeftAlignAndTrimVariants; parallel collection in CalibrateDragstrModel; FlowBasedHMM and FlowFeatureMapper's GVCF writer; CpxVariantInterpreter; SVCluster `DEFRAGMENT_CNV`; FilterAlignmentArtifacts realignment; the PathSeq Bloom filter |
| PairHMM | only `LOGLESS_CACHING` (the pure-Java order); every other `-pairHMM` value is refused |
| DepthOfCoverage | the TABLE format, `--partition-type` beyond sample, gene coverage, `--include-deletions`, partitioned quantiles |
| R and JavaScript | `--chart`, `--produce-plot`, AnalyzeCovariates `--plots-report-file`; Nashorn JavaScript in FilterVcf |
| unseeded randomness | Mutect2 `--genotype-germline-sites-fraction` below 1; GroundTruthReadsBuilder `--subsampling-ratio` (an unseeded `java.util.Random` has no byte to match) |
| unmeasured hash order | HashMap order in PathSeqBuildReferenceTaxonomy; HashSet layout in ModelSegments |
| JVM error text | IOException, native-loader and Spark task-failure messages are not reproduced |
| remote storage | `gs://`, HDFS, remote GenomicsDB workspaces |
| GenomicsDBImport | 21 options, among them `--bypass-feature-reader`, `--consolidate`, `--genomicsdb-update-workspace-path`, `--merge-contigs-into-num-partitions` |
| Python-backed tools | a missing Python package, `--python-script`, gcnvkernel's host-compiled C |

`grep -rn PORT_LIMITATION crates/` lists every site with its message.

picard-rs and htsjdk-rs have no such type: they refuse through `Failure::Escaped` ("not ported"),
about 52 sites in picard-rs and 9 in htsjdk-rs.

## 2. Exclusions of the reference's own tests

PR #1478 replays GATK 4.7.0.0's integration tests (`src/test/resources`) on both sides: 1,655
invocations extracted over 168 tools, 1,613 of them for ported tools, and 297 tests not extracted,
each with its reason. **189 invocations are excluded**, listed with their reason in
`tools/coverage/test_resources/exclusions.json`:

- **42 AVX runs.** The test asks for `-pairHMM AVX_LOGLESS_CACHING`, which the reference itself
  refuses on the CI runner ("Machine does not support AVX PairHMM"): HaplotypeCaller 37,
  HaplotypeCallerSpark 2, DirichletAlleleDepthAndFraction 2, Mutect2 1.
- **147 port limitations** (section 1), by reason: CRAM input to Spark tools 36; CalibrateDragstrModel
  parallel collection 12; `--chart` 12; `-L unmapped` 12; non-plain PrintSVEvidence 8; GenotypeGVCFs
  deprecated RAW_MQ 7; DepthOfCoverage 6; GenotypeGVCFs `--input-is-somatic` 5; VariantEval
  `--pedigree` 5; Funcotator data sources 5; ReadNameReadFilter 4; AnnotateIntervals tracks 4;
  CollectSVEvidence `.gz` 4; CalculateGenotypePosteriors `--pedigree` 3; Pileup `--metadata` 3;
  GenomicsDBImport 3; and 18 more of one or two each.

## 3. Instabilities of the reference

Where GATK has no single answer, the port matches the one the oracle gives on real x86-64 and
says so; nothing here is a port bug.

| what | why the reference is unstable | what the port claims | source |
|---|---|---|---|
| `Double.toString` / `DecimalFormat` above 2^53 | Java 17's pre-Schubfach converter is not shortest there (JDK-4511638) | 68 values differ, all above 2^53; 0 of 903,121 below it | #115 |
| gCNV (DetermineGermlineContigPloidy, GermlineCNVCaller, PostprocessGermlineCNVCalls) | MKL, `-march` native and ifunc libm: as shipped there is no byte-identical target | byte-identical against a frozen environment with reference BLAS (`gcnv_frozen_env.sh`) | #545, #547, #549, #1505, [gcnv-reproducible-oracle.md](gcnv-reproducible-oracle.md) |
| `Math.exp` (FilterMutectCalls, ActivityProfile) | HotSpot's intrinsic differs by host; aarch64 Java has no exp stub and runs fdlibm | `java_math_exp`: 0 of 2,000,000 differ against arm64 `Math.exp` (glibc: 192,457); aarch64 Java itself differs from the x86-64 goldens on 3 ActivityProfile lines and 1 Mutect line | #1475, #1499, [numeric-functions-a-ported-call-site-reaches.md](numeric-functions-a-ported-call-site-reaches.md) |
| NVScoreVariants after the dense layers | MKL SGEMM sums in an order set by thread count and batch size (4 threads vs 1: 2,047 of 2,840 values, up to 5 ulp) | the VCF and its `%.3f` scores are byte-identical; intermediate floats are not claimed | [nvscorevariants.md](nvscorevariants.md), #1510 |
| CreateReadCountPanelOfNormals rank | the eigensample count is decided on rounding noise in Spark fold order (8, 6 and 7 on one fixture) | the rank is not claimed; see the doc | [a-rank-is-not-a-byte.md](a-rank-is-not-a-byte.md), #1270 |
| HaplotypeCaller, RampedHaplotypeCaller, Mutect2 debug outputs | print `System.identityHashCode` | compared with those masked | #1477 |
| TrainVariantAnnotationsModel, ScoreVariantAnnotations HDF5 | Java HDF5 stamps the clock (12 bytes) | compared with the stamp masked | #1477 |
| AnalyzeSaturationMutagenesis `--write-rejected-reads` | races the async BAM writer (4 of 18 reference runs differ) | held at `false` | #1477 |
| BwaMemIndexImageCreator | writes heap pointers under ASLR (9 bytes) | compared with them masked | #1477, [pointers-that-reach-the-output.md](pointers-that-reach-the-output.md) |
| GenomicsDBImport | names fragments by UUID, thread and clock; writes a `Date` and an absolute path | 133 of 134 files identical | #1477 |

## 4. Deliberate divergences from htsjdk 5.0.0

Taken from htsjdk 5.0.1, 5.0.2 and 6.0.0 because 5.0.0 loses data. Each one differs from GATK
4.7.0.0 and Picard 3.5.0 only on inputs where the reference's output is wrong. The full table, with
upstream issues, is htsjdk-rs `docs/htsjdk-upstream.md` (htsjdk-rs #293).

| where | 5.0.0 | the port | source |
|---|---|---|---|
| CRAM, mapped read with SEQ `*` | writes CIGAR `0M`, drops NM/MD | keeps CIGAR, NM and MD | htsjdk-rs #298; htsjdk #1837, #1843 |
| CRAM tags and names | loses unsigned `B` arrays and bytes 0x80–0xFF | keeps them | htsjdk-rs #298; htsjdk #1840, #1887 |
| CRAM, unpaired read | drops mate flags and PNEXT; drops an RG naming no `@RG` | keeps them | htsjdk-rs #298; htsjdk #1891 |
| CRAM with a mid-file EOF container (samtools cat < 1.13) | stops there | reads to the end | htsjdk-rs #298; htsjdk #1876 |
| BCF 2.1, interior missing value | shifts later values, cuts AD/PL | keeps it in place (a trailing `.` is still padding) | htsjdk-rs #295; htsjdk #1827, #1828 |
| BCF 2.1, Integer INFO list with a null | NullPointerException | writes it | htsjdk-rs #295; htsjdk #1889 |
| Indexed BAM/CRAM, multi-interval query, references out of `@SQ` order | drops records | returns every record | htsjdk-rs #294; htsjdk #1877; gatk-rs #1527 |
| Merged VCF headers, inputs' dictionaries differ | drops contigs | keeps every contig | htsjdk-rs #296; htsjdk #1878; gatk-rs #1527 |

Already htsjdk 6.0.0's behaviour rather than 5.0.0's: an interior empty BGZF block does not end the
file; ReferenceSource is keyed by contig name; `@RG DT` is copied as written. An empty interior SAM
field is refused at every validation stringency (5.0.0 goes by stringency).

## Platforms

- Goldens come only from real x86-64. The port runs, and CI tests it, on x86-64 and arm64; a
  difference between the two is a failure, except where the reference itself differs by
  architecture (the `Math.exp` row above).
- The release binaries are built without `target-cpu=native` and are run on QEMU's baseline
  `qemu64` CPU before they are packaged. Linux binaries need glibc 2.39 or later (Ubuntu 24.04,
  the build host).
- The GPU path (Metal, `--features gatk-engine/gpu`) is not in the release binaries.

/*
 * The argument declarations of three ported tools, taken from the reference.
 *
 * The Barclay parser is measured already, argument by argument and mechanism by mechanism. What is
 * not measured is what any ONE tool declares: the flattened namespace its own fields and its
 * inherited collections produce, and therefore which command lines it accepts. That is the layer
 * Milestone C's per-tool declarations have to reproduce, and it cannot be read off the parser.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE NAMESPACE IS FLAT AND INHERITED: a read walker declares a handful of arguments of its
 *     own and ends up with dozens, the rest coming from the collections its superclasses hold;
 *   - AND WHICH PARSER IS ASKED DECIDES HOW MANY THERE ARE: a parser built straight from the
 *     instance sees 38 where the one the tool hands out sees 70, the gap being the plugin
 *     descriptors' four read-filter arguments and the standard collections the tool adds when it
 *     builds its own. Every one of the 38 is in the 70, so the instance list is a subset and not
 *     a different reading;
 *   - THE ORDER IS THE PARSER'S OWN, subclass first, and it is what the usage prints in;
 *   - A REQUIRED ARGUMENT IS REQUIRED BY ITS DECLARATION AND NOT BY THE TOOL, and the two
 *     archetypes do not mirror each other: a read walker REQUIRES `--input` and has no `--variant`
 *     at all, while a variant walker requires `--variant` and takes `--input` as an OPTIONAL
 *     argument, so `-V` on the first is "not a recognized option" and `-I` on the second parses;
 *   - THE SHORT NAME IS PART OF THE DECLARATION, and several arguments have none;
 *   - A COLLECTION ARGUMENT IS DECLARED AS ONE, which is what lets it be repeated;
 *   - THE DEFAULT IS THE FIELD'S OWN VALUE at construction time, so it is read off an instance and
 *     not off the annotation;
 *   - AN ARGUMENT NAMED TWICE IS AN ERROR ONLY IF IT IS A SCALAR: `--input` is declared as a
 *     collection, so naming it twice parses, while `--output` is a scalar and naming it twice is
 *     refused;
 *   - AN UNKNOWN ARGUMENT IS REFUSED BY THE PARSER AND NOT BY THE TOOL, with a message that names
 *     it;
 *   - AND A MISSING REQUIRED ARGUMENT IS REFUSED THE SAME WAY, before the tool runs at all.
 *
 * Five more, added when the declarations had to carry enough to BUILD a parser rather than only to
 * count one.
 *
 *   - THE TYPE IS THE UNDERLYING FIELD'S, which for a collection is its ELEMENT class: `--input`
 *     is a `List<GATKPath>` and reports `GATKPath`, so the conversion a value goes through is the
 *     element's and the collection is only how many of them there may be;
 *   - PRIMITIVE IS A SEPARATE QUESTION FROM THE CLASS, the class being boxed either way, and it is
 *     the one the null check asks;
 *   - HIDDEN, ADVANCED AND COMMON ARE THREE DIFFERENT FLAGS, and they are what decides which
 *     section of the usage an argument is printed in, or whether it is printed at all;
 *   - A BOUNDED RANGE IS FOUR NULLABLE DOUBLES and not a pair, the recommended range being
 *     declared beside the hard one rather than instead of it;
 *   - AND THE DOCUMENTATION IS PART OF THE DECLARATION: it is the annotation's own string, which
 *     is what the usage text wraps, so it belongs to the argument and not to the renderer.
 *
 * Output:
 *
 *     count\t<tool>\t<how many named arguments it declares>
 *     def\t<tool>\t<index>\t<longName>|<aliases>|<required>|<collection>|<default>|<type>|
 *         <primitive>|<flag>|<hidden>|<advanced>|<common>|<minElements>|<maxElements>|
 *         <minValue>|<maxValue>|<minRecommended>|<maxRecommended>|<mutex>|<plugin>
 *     doc\t<tool>\t<index>\t<the documentation string, escaped>
 *     parse\t<tool>\t<case>\tok|E:<exception class>:<message>
 *
 * Usage: ToolArgumentDeclarationDump
 */

import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.NamedArgumentDefinition;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class ToolArgumentDeclarationDump {

    public static void main(final String[] args) {
        declarations("CountReads",
                new org.broadinstitute.hellbender.tools.CountReads());
        declarations("CountVariants",
                new org.broadinstitute.hellbender.tools.walkers.CountVariants());
        declarations("PrintReads",
                new org.broadinstitute.hellbender.tools.PrintReads());
        // Four more archetypes, so the list is not three walkers of two kinds: a read walker that
        // writes a recalibrated file, a variant walker with a large argument surface of its own, a
        // tool that is no walker at all, and one that takes a list of files rather than one.
        declarations("ApplyBQSR",
                new org.broadinstitute.hellbender.tools.walkers.bqsr.ApplyBQSR());
        declarations("SelectVariants",
                new org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants());
        declarations("IndexFeatureFile",
                new org.broadinstitute.hellbender.tools.IndexFeatureFile());
        // `GatherVcfs` is Picard's tool of that name, which GATK dispatches to; the GATK one is
        // GatherVcfsCloud, and the two declare different arguments. Naming it wrongly here made
        // the generator's cross-check report a `--COMMENT` the parser did not declare, which is
        // Picard's argument and not this tool's.
        declarations("GatherVcfsCloud",
                new org.broadinstitute.hellbender.tools.GatherVcfsCloud());
        // Two more tools that are no walkers, chosen because the port can RUN them: each takes a
        // file and answers, so a declaration for one of them is a command line the binary can be
        // handed end to end rather than only parsed.
        //
        // `CheckTerminatorBlock` and `BuildBamIndex` would have been two more and are not here:
        // they are PICARD's tools, which GATK dispatches to, so their declarations come from
        // Picard's own parser and belong in picard-rs's measurement rather than this one.
        // The two read walkers that share `CountReads`' plumbing exactly: reads in, a number or a
        // report out. They are here because an archetype is the unit of this milestone -- a tool
        // of one archetype shares its argument shape AND its file plumbing with the others, so the
        // second and third cost a fraction of the first.
        declarations("CountBases",
                new org.broadinstitute.hellbender.tools.CountBases());
        declarations("FlagStat",
                new org.broadinstitute.hellbender.tools.FlagStat());
        // Two archetypes neither the read walkers nor the variant walkers reach: a REFERENCE
        // walker, whose traversal is the FASTA rather than a file of records, and an interval
        // utility, which is a `GATKTool` with no traversal at all and writes a DIRECTORY of files.
        declarations("CountBasesInReference",
                new org.broadinstitute.hellbender.tools.walkers.fasta.CountBasesInReference());
        declarations("SplitIntervals",
                new org.broadinstitute.hellbender.tools.walkers.SplitIntervals());
        // A LOCUS walker, whose traversal is the pileup at each position rather than a record or
        // a base, and a second interval utility, which shares `SplitIntervals`' plumbing. Both
        // dumps gain them together: an enum missing from the enums table is an argument the
        // port's parser drops whole, which cost `SplitIntervals` a whole CI round.
        declarations("Pileup",
                new org.broadinstitute.hellbender.tools.walkers.qc.Pileup());
        declarations("PreprocessIntervals",
                new org.broadinstitute.hellbender.tools.copynumber.PreprocessIntervals());
        // A read walker that WRITES a BAM, which is `PrintReads`' plumbing with a filter in front
        // of it, and a `GATKTool` that opens the reads only to read their header. Both dumps
        // together, as ever.
        declarations("PrintDistantMates",
                new org.broadinstitute.hellbender.tools.PrintDistantMates());
        declarations("GetSampleName",
                new org.broadinstitute.hellbender.tools.GetSampleName());
        // A `CommandLineProgram` that is no GATKTool at all -- fifteen arguments, two interval
        // lists in and a verdict out -- and a read walker that rewrites qualities and writes a BAM.
        declarations("CompareIntervalLists",
                new org.broadinstitute.hellbender.tools.CompareIntervalLists());
        declarations("FixMisencodedBaseQualityReads",
                new org.broadinstitute.hellbender.tools.FixMisencodedBaseQualityReads());
        // A second LOCUS walker, which compares the engine's pileup to a samtools one, and a
        // third interval utility, which annotates each interval with its GC content.
        declarations("CheckPileup",
                new org.broadinstitute.hellbender.tools.walkers.qc.CheckPileup());
        declarations("AnnotateIntervals",
                new org.broadinstitute.hellbender.tools.copynumber.AnnotateIntervals());
        // A second REFERENCE walker, which writes a FASTA rather than a number, and a third LOCUS
        // walker, which counts the reads over each interval. `CollectReadCounts` brings an enum of
        // its own -- `--format`, whose two constants are the TSV and the HDF5 it writes -- and an
        // enum missing from the enums table is an argument the port's parser drops whole, so the
        // two dumps gain the pair together as ever.
        declarations("FastaReferenceMaker",
                new org.broadinstitute.hellbender.tools.walkers.fasta.FastaReferenceMaker());
        declarations("CollectReadCounts",
                new org.broadinstitute.hellbender.tools.copynumber.CollectReadCounts());
        // A third LOCUS walker, which is the first tool here whose traversal drives a FEATURE
        // source -- its `-V` is the allele frequencies it summarises the pileup at -- and a read
        // walker that clears one flag bit and writes the reads back, which is `PrintReads`'
        // plumbing with a mutation in the middle.
        declarations("GetPileupSummaries",
                new org.broadinstitute.hellbender.tools.walkers.contamination.GetPileupSummaries());
        declarations("UnmarkDuplicates",
                new org.broadinstitute.hellbender.tools.walkers.UnmarkDuplicates());
        // A `CommandLineProgram` that is no GATKTool -- it reads the table `GetPileupSummaries`
        // writes, so the two together are the first CHAIN of ported tools a command line can run
        // end to end -- and a `GATKTool` with no traversal at all, which opens the reads only to
        // print their header.
        declarations("CalculateContamination",
                new org.broadinstitute.hellbender.tools.walkers.contamination.CalculateContamination());
        declarations("PrintReadsHeader",
                new org.broadinstitute.hellbender.tools.PrintReadsHeader());
        // A third LOCUS walker, which writes TWO files -- a BED of runs and a summary of counts --
        // and a second REFERENCE utility, which writes a shifted FASTA with its dictionary and two
        // interval lists beside it. Neither shape has been declared here before: every tool
        // measured so far writes one file or a directory of them.
        declarations("CallableLoci",
                new org.broadinstitute.hellbender.tools.walkers.coverage.CallableLoci());
        declarations("ShiftFasta",
                new org.broadinstitute.hellbender.tools.walkers.fasta.ShiftFasta());
        // Two more of the record-transform archetype, both already ported and neither runnable
        // until its arguments are declared. `RevertBaseQualityScores` is the one that ABORTS on a
        // read without `OQ` rather than skipping it, and `AddOriginalAlignmentTags` is the first
        // that writes tags rather than changing the read.
        declarations("RevertBaseQualityScores",
                new org.broadinstitute.hellbender.tools.walkers.RevertBaseQualityScores());
        declarations("AddOriginalAlignmentTags",
                new org.broadinstitute.hellbender.tools.AddOriginalAlignmentTags());
        // The next two. `LeftAlignIndels` is a read walker that needs a REFERENCE, which is the
        // first required argument in this dump that is not the reads, and `DumpTabixIndex` is a
        // tool that is no walker at all and whose whole namespace is fourteen arguments.
        declarations("LeftAlignIndels",
                new org.broadinstitute.hellbender.tools.LeftAlignIndels());
        declarations("DumpTabixIndex",
                new org.broadinstitute.hellbender.tools.DumpTabixIndex());
        // The next two. `ReadAnonymizer` is a read walker that needs a reference AND rewrites the
        // bases it reads, and `PrintFileDiagnostics` declares fifteen arguments and is no walker.
        declarations("ReadAnonymizer",
                new org.broadinstitute.hellbender.tools.walkers.ReadAnonymizer());
        declarations("PrintFileDiagnostics",
                new org.broadinstitute.hellbender.tools.PrintFileDiagnostics());
        // The next two, and both write MORE than one output: `SplitReads` writes one file per
        // key of up to three splitters, and `ClipReads` writes a statistics file beside its BAM.
        declarations("SplitReads",
                new org.broadinstitute.hellbender.tools.SplitReads());
        declarations("ClipReads",
                new org.broadinstitute.hellbender.tools.ClipReads());
        // The next two, and both need a REFERENCE: `SplitNCigarReads` splits a read at every `N`
        // and `MethylationTypeCaller` writes a VCF rather than reads.
        declarations("SplitNCigarReads",
                new org.broadinstitute.hellbender.tools.walkers.rnaseq.SplitNCigarReads());
        declarations("MethylationTypeCaller",
                new org.broadinstitute.hellbender.tools.walkers.MethylationTypeCaller());
        // The next two. `BaseRecalibrator` writes a GATKReport rather than reads and takes a
        // FeatureInput of known sites, and `GtfToBed` reads an annotation and writes a BED.
        declarations("BaseRecalibrator",
                new org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator());
        declarations("GtfToBed",
                new org.broadinstitute.hellbender.tools.walkers.conversion.GtfToBed());
        // Two tools of the record-transform archetype that are no WALKERS: both extend `GATKTool`
        // and override `traverse()`, so a second reads source is opened by hand and the engine's
        // filter, transformer and interval machinery never runs. `TransferReadTags` walks two files
        // in lockstep and copies tags from the unmapped one, and `PostProcessReadsForRSEM` reorders
        // a query-name-sorted file into the pairs RSEM will read.
        declarations("TransferReadTags",
                new org.broadinstitute.hellbender.tools.walkers.qc.TransferReadTags());
        declarations("PostProcessReadsForRSEM",
                new org.broadinstitute.hellbender.tools.walkers.qc.PostProcessReadsForRSEM());
        // A `VariantWalker` that writes a TABLE rather than a VCF, and a tool that is no GATK tool
        // at all: `CompareBaseQualities` extends `PicardCommandLineProgram`, so its namespace is
        // Picard's argument set and not the engine's, which nothing declared here has been.
        declarations("VariantsToTable",
                new org.broadinstitute.hellbender.tools.walkers.variantutils.VariantsToTable());
        declarations("CompareBaseQualities",
                new org.broadinstitute.hellbender.tools.validation.CompareBaseQualities());
        // The first two tools here that WRITE a VCF, which is a writer this port has not used from a
        // runner before. `RemoveNearbyIndels` buffers one indel at a time and drops any pair closer
        // than a spacing; `UpdateVCFSequenceDictionary` replaces the header's dictionary and passes
        // every record through, and its `--source-dictionary` is the first argument here that takes
        // a dictionary from any of four file kinds.
        declarations("RemoveNearbyIndels",
                new org.broadinstitute.hellbender.tools.walkers.validation.RemoveNearbyIndels());
        declarations("UpdateVCFSequenceDictionary",
                new org.broadinstitute.hellbender.tools.walkers.variantutils.UpdateVCFSequenceDictionary());
        declarations("PrintBGZFBlockInformation",
                new org.broadinstitute.hellbender.tools.PrintBGZFBlockInformation());
        declarations("CreateHadoopBamSplittingIndex",
                new org.broadinstitute.hellbender.tools.spark.CreateHadoopBamSplittingIndex());
        // Two REFERENCE utilities. `FastaAlternateReferenceMaker` is the maker with a VCF applied,
        // and it declares TWO `FeatureInput`s: the driving variants and `--snp-mask`. The pair
        // matters to the parser only in what it does not do, which is check them against each
        // other: `--snp-mask-priority` without `--snp-mask` is refused by the TOOL, from
        // `onTraversalStart`, so the command line that carries it parses cleanly.
        // `CompareReferences` requires a LIST of references, `--references-to-compare`, where
        // every tool declared here so far takes the one `--reference`.
        declarations("FastaAlternateReferenceMaker",
                new org.broadinstitute.hellbender.tools.walkers.fasta.FastaAlternateReferenceMaker());
        declarations("CompareReferences",
                new org.broadinstitute.hellbender.tools.reference.CompareReferences());
        // Two more REFERENCE utilities, and each brings a shape the dump has not carried.
        // `CheckReferenceCompatibility` takes the reads and a VCF as the things to CHECK a
        // reference against, so its required argument is not the reference at all, and it refuses
        // a command line that names both. `ComposeSTRTableFile` is a DRAGstr tool: it writes a
        // binary table beside a reference and declares the sampling arguments that decide what
        // goes in it.
        declarations("CheckReferenceCompatibility",
                new org.broadinstitute.hellbender.tools.reference.CheckReferenceCompatibility());
        declarations("ComposeSTRTableFile",
                new org.broadinstitute.hellbender.tools.dragstr.ComposeSTRTableFile());
        // A CHAIN, which is why the pair is this one: `CalculateMixingFractions` reads a VCF and
        // the reads and writes a table of one fraction per sample, and
        // `AnnotateVcfWithExpectedAlleleFraction` reads THAT table beside a VCF and writes the
        // expected fraction into each record. The second tool's input is the first tool's output,
        // so the corpus can carry a table the reference itself produced.
        declarations("CalculateMixingFractions",
                new org.broadinstitute.hellbender.tools.walkers.validation.CalculateMixingFractions());
        declarations("AnnotateVcfWithExpectedAlleleFraction",
                new org.broadinstitute.hellbender.tools.walkers.validation.AnnotateVcfWithExpectedAlleleFraction());
        // The third and fourth of the validation walkers, and the pair completes a family. Both
        // read a VCF; `AnnotateVcfWithBamDepth` also reads the reads and writes one Integer INFO
        // field per record, and `CountFalsePositives` writes a TABLE of counts per variant type
        // over a target territory that `-L` decides. The first is the sibling whose default tool
        // header lines DO reach the file, which is what makes the pair with
        // `AnnotateVcfWithExpectedAlleleFraction` worth having declared together.
        declarations("AnnotateVcfWithBamDepth",
                new org.broadinstitute.hellbender.tools.walkers.validation.AnnotateVcfWithBamDepth());
        declarations("CountFalsePositives",
                new org.broadinstitute.hellbender.tools.walkers.validation.CountFalsePositives());
        // Two tools whose answer is a TABLE built from a traversal rather than a transformed file.
        // `EvaluateInfoFieldConcordance` reads one INFO key out of each of TWO VCFs, which makes it
        // the first tool declared here with a second variant input beside the driving one, and
        // `CallCopyRatioSegments` is the first copy-number tool: it reads segments, computes a
        // length-weighted mean and deviation over the copy-neutral ones TWICE, and calls each
        // segment against the second pair.
        declarations("EvaluateInfoFieldConcordance",
                new org.broadinstitute.hellbender.tools.walkers.validation.EvaluateInfoFieldConcordance());
        declarations("CallCopyRatioSegments",
                new org.broadinstitute.hellbender.tools.copynumber.CallCopyRatioSegments());
        // A VCF transformed and a table collected, and both are bigger than anything declared here
        // so far. `VariantFiltration` has SIXTY-ONE arguments and writes its filters into the
        // FILTER column and into a genotype's FT; `CollectAllelicCounts` walks loci over a BAM and
        // counts the reference and the alternate base at each site of an interval list, which is
        // the first copy-number tool here that reads reads.
        declarations("VariantFiltration",
                new org.broadinstitute.hellbender.tools.walkers.filters.VariantFiltration());
        declarations("CollectAllelicCounts",
                new org.broadinstitute.hellbender.tools.copynumber.CollectAllelicCounts());
        // A filter over intervals and a collector over two BAMs. `FilterIntervals` reads the
        // annotated intervals `AnnotateIntervals` writes and the counts `CollectReadCounts` writes,
        // which makes it the second CHAIN declared here; `GetNormalArtifactData` is the first
        // Mutect tool of any kind, and it reads a tumour and a normal at once.
        declarations("FilterIntervals",
                new org.broadinstitute.hellbender.tools.copynumber.FilterIntervals());
        declarations("GetNormalArtifactData",
                new org.broadinstitute.hellbender.tools.walkers.mutect.GetNormalArtifactData());
        // A validator and a normaliser, both over a VCF. `ValidateVariants` writes NOTHING: its
        // whole answer is whether it threw and what it said, and its `--validation-type-to-exclude`
        // brings the `ValidationType` enum, which nothing declared here points at. Its GVCF mode is
        // three checks rather than one, the last of which runs at the end of the traversal.
        // `LeftAlignAndTrimVariants` writes a VCF and is the first declared tool that can SPLIT a
        // record into several: `--split-multi-allelics` turns one multi-allelic line into one per
        // alternate, each trimmed on its own.
        declarations("ValidateVariants",
                new org.broadinstitute.hellbender.tools.walkers.variantutils.ValidateVariants());
        declarations("LeftAlignAndTrimVariants",
                new org.broadinstitute.hellbender.tools.walkers.variantutils.LeftAlignAndTrimVariants());
        // A gatherer and a comparator, and neither is an ordinary walker. `GatherBQSRReports`
        // merges the recalibration tables a scattered run wrote, which makes it the third CHAIN
        // declared here: its input is what `BaseRecalibrator` writes. `Concordance` is the first
        // `AbstractConcordanceWalker` declared here, and it drives TWO variant files at once: the
        // truth and the evaluation, matched variant by variant.
        declarations("GatherBQSRReports",
                new org.broadinstitute.hellbender.tools.walkers.bqsr.GatherBQSRReports());
        declarations("Concordance",
                new org.broadinstitute.hellbender.tools.walkers.validation.Concordance());
        // A coverage walker and a genotype refiner. `DepthOfCoverage` is the first declared tool
        // that writes a DIRECTORY of tables rather than one file, and it partitions its counts by
        // sample, read group and library at once. `CalculateGenotypePosteriors` reads a supporting
        // callset beside its own and rewrites every genotype's likelihoods, which makes it the
        // fourth CHAIN declared here: the corpus's population VCF is the support.
        declarations("DepthOfCoverage",
                new org.broadinstitute.hellbender.tools.walkers.coverage.DepthOfCoverage());
        declarations("CalculateGenotypePosteriors",
                new org.broadinstitute.hellbender.tools.walkers.variantutils.CalculateGenotypePosteriors());
        // A denoiser and a validator. `DenoiseReadCounts` is the next link of the copy-number
        // chain: it reads the counts `CollectReadCounts` writes and standardises them, and its
        // panel of normals is an HDF5 file the corpus does not carry, so the argument that takes
        // one is the measurable refusal. `ValidateBasicSomaticShortMutations` checks a callset
        // against the pileups of a TUMOUR and a NORMAL at once, which makes it the second declared
        // tool that drives two BAMs.
        //
        // `VariantAnnotator` was the first choice for the pair and cannot be declared: it carries
        // no summary in the usage golden or the inventory, and the declaration generator refuses a
        // tool it cannot document.
        declarations("DenoiseReadCounts",
                new org.broadinstitute.hellbender.tools.copynumber.DenoiseReadCounts());
        declarations("ValidateBasicSomaticShortMutations",
                new org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup.ValidateBasicSomaticShortMutations());
        // Two SV-adjacent collectors, and both write a FORMAT this table has not seen. `PrintReadCounts`
        // rewrites a counts file into the SV pipeline's own shape, and `CollectSVEvidence` walks a
        // BAM for the paired-end, split-read and depth evidence that pipeline consumes. Each takes
        // a `--sample-name` of its own, which is the first argument here whose value has to agree
        // with the reads rather than with a file.
        declarations("PrintReadCounts",
                new org.broadinstitute.hellbender.tools.sv.PrintReadCounts());
        declarations("CollectSVEvidence",
                new org.broadinstitute.hellbender.tools.walkers.sv.CollectSVEvidence());
        // Two builders from the PathSeq subsystem, neither of them Spark despite the package they
        // sit in. `PathSeqBuildKmers` writes the host reference's k-mer set, and it carries an
        // argument whose legal values are narrower than its declared bounds: `--kmer-size` is 1 to
        // 31 to the parser and must be ODD to the tool, which a declaration cannot say.
        // `PathSeqBuildReferenceTaxonomy` declares TWO optional catalogue inputs and requires one
        // of them, so "At least one of --refseq-catalog or --genbank-catalog must be specified" is
        // a refusal the parser never raises.
        declarations("PathSeqBuildKmers",
                new org.broadinstitute.hellbender.tools.spark.pathseq.PathSeqBuildKmers());
        declarations("PathSeqBuildReferenceTaxonomy",
                new org.broadinstitute.hellbender.tools.spark.pathseq.PathSeqBuildReferenceTaxonomy());
        // A structural-variant annotator and a shard converter, paired because each declares a
        // shape the table has not held. `SVAnnotate` is a `VariantWalker` whose `--output` is
        // OPTIONAL and whose absence means stdout, and its `--max-breakend-as-cnv-length` is
        // declared `minValue = 0` with a default of `-1`: a default OUTSIDE its own range, which is
        // legal because `isValueOutOfRange` runs on a value the command line SET and never on the
        // default. `ConvertHeaderlessHadoopBamShardToBam` is the opposite shape, a
        // `CommandLineProgram` whose three arguments are all required and all plain files, and the
        // only declared tool that reads a header from one file and records from another.
        declarations("SVAnnotate",
                new org.broadinstitute.hellbender.tools.walkers.sv.SVAnnotate());
        declarations("ConvertHeaderlessHadoopBamShardToBam",
                new org.broadinstitute.hellbender.tools.ConvertHeaderlessHadoopBamShardToBam());
        // Three structural-variant evidence tools, the rest of the family `CollectSVEvidence`
        // writes and `PrintReadCounts` reads. `CondenseDepthEvidence` is a `FeatureWalker` whose
        // driving file is an argument of its own, `--depth-evidence`, rather than the engine's
        // `--feature-file`. `PrintSVEvidence` and `SiteDepthtoBAF` are `MultiFeatureWalker`s: the
        // first takes any number of evidence files of one type and merges them, the second reads
        // site depths and writes B-allele frequencies. All three choose their output's format from
        // its NAME, so `--output` is a `GATKPath` whose extension is part of its value.
        declarations("CondenseDepthEvidence",
                new org.broadinstitute.hellbender.tools.sv.CondenseDepthEvidence());
        declarations("PrintSVEvidence",
                new org.broadinstitute.hellbender.tools.sv.PrintSVEvidence());
        declarations("SiteDepthtoBAF",
                new org.broadinstitute.hellbender.tools.sv.SiteDepthtoBAF());
        // Eight tools whose ports are oracle-backed and whose inputs are tables rather than reads.
        // `GatherTranches` gathers VQSR tranches and declares `--mode`, the
        // `VariantRecalibratorArgumentCollection$Mode` enum the class-name key was made for. The
        // three Mutect gathers each concatenate or sum one kind of table. The four copy-number
        // utilities rewrite segment and region files. All but `GatherTranches` are undocumented.
        declarations("GatherTranches",
                new org.broadinstitute.hellbender.tools.walkers.vqsr.GatherTranches());
        declarations("GatherPileupSummaries",
                new org.broadinstitute.hellbender.tools.walkers.contamination.GatherPileupSummaries());
        declarations("GatherNormalArtifactData",
                new org.broadinstitute.hellbender.tools.walkers.mutect.GatherNormalArtifactData());
        declarations("MergeMutectStats",
                new org.broadinstitute.hellbender.tools.walkers.mutect.MergeMutectStats());
        declarations("CombineSegmentBreakpoints",
                new org.broadinstitute.hellbender.tools.copynumber.utils.CombineSegmentBreakpoints());
        declarations("MergeAnnotatedRegions",
                new org.broadinstitute.hellbender.tools.copynumber.utils.MergeAnnotatedRegions());
        declarations("MergeAnnotatedRegionsByAnnotation",
                new org.broadinstitute.hellbender.tools.copynumber.utils.MergeAnnotatedRegionsByAnnotation());
        declarations("TagGermlineEvents",
                new org.broadinstitute.hellbender.tools.copynumber.utils.TagGermlineEvents());
        // Eleven walkers whose ports are oracle-backed: three variant filters, two Mutect
        // mitochondrial filters, a read-orientation model, a funcotation filter, two concordance
        // tools, an allele-frequency check, a duplicate-set downsampler and the example
        // `MultiFeatureWalker`, whose walk the SV evidence tools already share.
        declarations("FilterVariantTranches",
                new org.broadinstitute.hellbender.tools.walkers.vqsr.FilterVariantTranches());
        declarations("FilterFuncotations",
                new org.broadinstitute.hellbender.tools.funcotator.FilterFuncotations());
        declarations("LearnReadOrientationModel",
                new org.broadinstitute.hellbender.tools.walkers.readorientation.LearnReadOrientationModel());
        declarations("AlleleFrequencyQC",
                new org.broadinstitute.hellbender.tools.walkers.varianteval.AlleleFrequencyQC());
        declarations("CalculateAverageCombinedAnnotations",
                new org.broadinstitute.hellbender.tools.CalculateAverageCombinedAnnotations());
        declarations("DownsampleByDuplicateSet",
                new org.broadinstitute.hellbender.tools.walkers.consensus.DownsampleByDuplicateSet());
        declarations("ExampleMultiFeatureWalker",
                new org.broadinstitute.hellbender.tools.examples.ExampleMultiFeatureWalker());
        declarations("MTLowHeteroplasmyFilterTool",
                new org.broadinstitute.hellbender.tools.walkers.mutect.filtering.MTLowHeteroplasmyFilterTool());
        declarations("MergeMutect2CallsWithMC3",
                new org.broadinstitute.hellbender.tools.walkers.validation.MergeMutect2CallsWithMC3());
        declarations("NuMTFilterTool",
                new org.broadinstitute.hellbender.tools.walkers.mutect.filtering.NuMTFilterTool());
        declarations("ReferenceBlockConcordance",
                new org.broadinstitute.hellbender.tools.walkers.validation.ReferenceBlockConcordance());
        // Thirteen tools whose ports are oracle-backed and whose inputs are VCFs or GVCFs: the GVCF
        // genotyper, combiner and reblocker, VQSR's application, the DRAGstr calibration,
        // an allele-specific read counter, the structural-variant clusterers and their concordance
        // and stratification, and two undocumented ones, a VCF comparator and the annotator.
        declarations("ASEReadCounter",
                new org.broadinstitute.hellbender.tools.walkers.rnaseq.ASEReadCounter());
        declarations("ApplyVQSR",
                new org.broadinstitute.hellbender.tools.walkers.vqsr.ApplyVQSR());
        declarations("CalibrateDragstrModel",
                new org.broadinstitute.hellbender.tools.dragstr.CalibrateDragstrModel());
        declarations("CombineGVCFs",
                new org.broadinstitute.hellbender.tools.walkers.CombineGVCFs());
        // `GenotypeGVCFs`, whose instance-built parser refuses its own definitions: the count row
        // says so, and the tool's own parser, plugins included, is what is declared.
        declarations("GenotypeGVCFs",
                new org.broadinstitute.hellbender.tools.walkers.GenotypeGVCFs());
        declarations("GnarlyGenotyper",
                new org.broadinstitute.hellbender.tools.walkers.gnarlyGenotyper.GnarlyGenotyper());
        declarations("ReblockGVCF",
                new org.broadinstitute.hellbender.tools.walkers.variantutils.ReblockGVCF());
        declarations("GroupedSVCluster",
                new org.broadinstitute.hellbender.tools.walkers.sv.GroupedSVCluster());
        declarations("JointGermlineCNVSegmentation",
                new org.broadinstitute.hellbender.tools.walkers.sv.JointGermlineCNVSegmentation());
        declarations("SVCluster",
                new org.broadinstitute.hellbender.tools.walkers.sv.SVCluster());
        declarations("SVConcordance",
                new org.broadinstitute.hellbender.tools.walkers.sv.SVConcordance());
        declarations("SVStratify",
                new org.broadinstitute.hellbender.tools.walkers.sv.SVStratify());
        declarations("VCFComparator",
                new org.broadinstitute.hellbender.tools.walkers.variantutils.VCFComparator());
        declarations("VariantAnnotator",
                new org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotator());
        // Twenty-seven tools whose ports are oracle-backed and that no earlier lot declared: the flow
        // tools, the BQSR plots, saturation mutagenesis, the BWA index image, the two CRAM
        // utilities, the F1R2 counter, the panels of normals, the scalable-VQSR extractor, the
        // alignment-artefact and Mutect filters, the funcotators and their downloader, the RNA
        // expression counter, the ground-truth tools, two callers, the local assembler, the
        // segment modeller, VariantEval and VariantRecalibrator. `StructuralVariantDiscoverer` is a
        // Spark tool and waits for the Spark argument surface.
        declarations("AddFlowBaseQuality",
                new org.broadinstitute.hellbender.tools.walkers.groundtruth.AddFlowBaseQuality());
        declarations("AddFlowSNVQuality",
                new org.broadinstitute.hellbender.tools.walkers.featuremapping.AddFlowSNVQuality());
        declarations("AnalyzeCovariates",
                new org.broadinstitute.hellbender.tools.walkers.bqsr.AnalyzeCovariates());
        declarations("AnalyzeSaturationMutagenesis",
                new org.broadinstitute.hellbender.tools.AnalyzeSaturationMutagenesis());
        declarations("BwaMemIndexImageCreator",
                new org.broadinstitute.hellbender.tools.BwaMemIndexImageCreator());
        declarations("CRAMIssue8768Detector",
                new org.broadinstitute.hellbender.tools.CRAMIssue8768Detector());
        declarations("CollectF1R2Counts",
                new org.broadinstitute.hellbender.tools.walkers.readorientation.CollectF1R2Counts());
        declarations("CreateReadCountPanelOfNormals",
                new org.broadinstitute.hellbender.tools.copynumber.CreateReadCountPanelOfNormals());
        declarations("CreateSomaticPanelOfNormals",
                new org.broadinstitute.hellbender.tools.walkers.mutect.CreateSomaticPanelOfNormals());
        declarations("ExtractVariantAnnotations",
                new org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.ExtractVariantAnnotations());
        declarations("FilterAlignmentArtifacts",
                new org.broadinstitute.hellbender.tools.walkers.realignmentfilter.FilterAlignmentArtifacts());
        declarations("FilterMutectCalls",
                new org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FilterMutectCalls());
        declarations("FlowFeatureMapper",
                new org.broadinstitute.hellbender.tools.walkers.featuremapping.FlowFeatureMapper());
        declarations("FlowPairHMMAlignReadsToHaplotypes",
                new org.broadinstitute.hellbender.tools.walkers.featuremapping.FlowPairHMMAlignReadsToHaplotypes());
        declarations("FuncotateSegments",
                new org.broadinstitute.hellbender.tools.funcotator.FuncotateSegments());
        declarations("Funcotator",
                new org.broadinstitute.hellbender.tools.funcotator.Funcotator());
        declarations("FuncotatorDataSourceDownloader",
                new org.broadinstitute.hellbender.tools.funcotator.FuncotatorDataSourceDownloader());
        declarations("GeneExpressionEvaluation",
                new org.broadinstitute.hellbender.tools.walkers.rnaseq.GeneExpressionEvaluation());
        declarations("GroundTruthReadsBuilder",
                new org.broadinstitute.hellbender.tools.walkers.groundtruth.GroundTruthReadsBuilder());
        declarations("GroundTruthScorer",
                new org.broadinstitute.hellbender.tools.walkers.groundtruth.GroundTruthScorer());
        declarations("HaplotypeBasedVariantRecaller",
                new org.broadinstitute.hellbender.tools.walkers.variantrecalling.HaplotypeBasedVariantRecaller());
        declarations("LocalAssembler",
                new org.broadinstitute.hellbender.tools.LocalAssembler());
        declarations("ModelSegments",
                new org.broadinstitute.hellbender.tools.copynumber.ModelSegments());
        declarations("RampedHaplotypeCaller",
                new org.broadinstitute.hellbender.tools.walkers.haplotypecaller.RampedHaplotypeCaller());
        declarations("SplitCRAM",
                new org.broadinstitute.hellbender.tools.SplitCRAM());
        declarations("VariantEval",
                new org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEval());
        declarations("VariantRecalibrator",
                new org.broadinstitute.hellbender.tools.walkers.vqsr.VariantRecalibrator());

        // A read walker: its own input, and the arguments it inherits.
        parse("CountReads", "no-arguments", new String[]{});
        parse("CountReads", "input-only", new String[]{"-I", "/dev/null"});
        parse("CountReads", "long-input", new String[]{"--input", "/dev/null"});
        parse("CountReads", "input-twice",
                new String[]{"-I", "/dev/null", "-I", "/dev/null"});
        parse("CountReads", "unknown-argument", new String[]{"--no-such-argument", "1"});
        parse("CountReads", "an-interval", new String[]{"-I", "/dev/null", "-L", "chr1"});
        parse("CountReads", "two-intervals",
                new String[]{"-I", "/dev/null", "-L", "chr1", "-L", "chr2"});
        parse("CountReads", "a-variant-argument",
                new String[]{"-I", "/dev/null", "-V", "/dev/null"});

        // A variant walker, which has no `--input` at all.
        parse("CountVariants", "variant-only", new String[]{"-V", "/dev/null"});
        parse("CountVariants", "an-input", new String[]{"-V", "/dev/null", "-I", "/dev/null"});
        parse("CountVariants", "no-arguments", new String[]{});

        // The output argument a print tool requires on top of its input.
        parse("PrintReads", "input-only", new String[]{"-I", "/dev/null"});
        // A tool that is no walker takes a positional-looking input and nothing else required.
        parse("IndexFeatureFile", "no-arguments", new String[]{});
        parse("IndexFeatureFile", "input-only", new String[]{"-I", "/dev/null"});
        parse("IndexFeatureFile", "an-interval", new String[]{"-I", "/dev/null", "-L", "chr1"});

        parse("PrintReads", "input-and-output",
                new String[]{"-I", "/dev/null", "-O", "/dev/null"});
        // The output is a scalar where the input is a collection, so naming it twice is refused.
        parse("PrintReads", "output-twice", new String[]{
            "-I", "/dev/null", "-O", "/dev/null", "-O", "/dev/null"});

        // POSITIONAL arguments, which no tool declared here has had. `CompareBaseQualities` takes
        // exactly two SAM files that way (`@PositionalArguments(minElements = 2, maxElements = 2)`),
        // so all four counts are measured: none, one, the two it wants, and three. What the parser
        // says for the three that are wrong is the whole reason these cases exist -- a port that
        // guessed the messages would guess four of them.
        parse("CompareBaseQualities", "positional-none", new String[]{});
        parse("CompareBaseQualities", "positional-one", new String[]{"/dev/null"});
        parse("CompareBaseQualities", "positional-two",
            new String[]{"/dev/null", "/dev/null"});
        parse("CompareBaseQualities", "positional-three",
            new String[]{"/dev/null", "/dev/null", "/dev/null"});
        // And the same pair with an argument that is named, so the order of the two kinds is
        // measured rather than assumed.
        parse("CompareBaseQualities", "positional-and-named",
            new String[]{"/dev/null", "-O", "/dev/null", "/dev/null"});
    }

    /**
     * Every named argument the tool declares, in the parser's own order, and how many a parser
     * built straight from the instance would have seen instead.
     *
     * The two are not the same list. A parser constructed from the instance knows nothing about
     * the plugin descriptors or the standard argument collections the tool adds when it builds its
     * own, so it sees 38 arguments where the tool's parser sees 70. The four read-filter arguments
     * are the visible half of that gap; the rest are the common and advanced ones the usage text
     * does not print either.
     */
    static void declarations(final String tool, final Object target) {
        // The parser built over the instance alone holds no plugin arguments, and `GenotypeGVCFs`
        // declares a mutex, `keep-combined-raw-annotations` against
        // `keep-specific-combined-raw-annotation`, whose target only its annotation plugin
        // contributes. That parser refuses the definition outright, so its count is recorded as
        // refused rather than as a number, and no argument is reported as the tool's alone.
        List<NamedArgumentDefinition> instance = null;
        try {
            instance = new CommandLineArgumentParser(target).getNamedArgumentDefinitions();
        } catch (final org.broadinstitute.barclay.argparser.CommandLineException.CommandLineParserInternalException e) {
            instance = null;
        }
        final List<NamedArgumentDefinition> definitions =
                ((CommandLineArgumentParser) ((org.broadinstitute.hellbender.cmdline
                        .CommandLineProgram) target).getCommandLineParser())
                        .getNamedArgumentDefinitions();
        System.out.printf("count\t%s\tinstance=%s tool=%d%n", tool,
                instance == null ? "refused" : String.valueOf(instance.size()),
                definitions.size());
        final List<String> seen = new ArrayList<>();
        if (instance != null) {
            for (final NamedArgumentDefinition definition : instance) {
                seen.add(definition.getLongName());
            }
        }
        for (int i = 0; i < definitions.size(); i++) {
            final NamedArgumentDefinition definition = definitions.get(i);
            final List<String> shorts = new ArrayList<>(definition.getArgumentAliases());
            // The class the parser converts a value to is the UNDERLYING field's, which for a
            // collection is its element class and not the collection's own.
            final String type = definition.getUnderlyingFieldClass().getSimpleName();
            // The same simple name belongs to more than one class: GATK declares two enums called
            // `Mode`, the variant-output filter's and the VQSR one, and a table keyed by the simple
            // name answers with whichever was seen first (IPNP-BIPN/gatk-rs#1179). The BINARY name
            // is what identifies a class, nested ones included, so the declaration carries it and
            // the enum table is joined on it.
            final String typeClass = definition.getUnderlyingFieldClass().getName();
            final boolean primitive = definition.getUnderlyingField().getType().isPrimitive();
            // NOT sorted. The message a mutex violation prints joins this list in the order it
            // holds, and that order is neither alphabetical nor the annotation's: `quantize-quals`
            // declares no mutex at all, and Barclay fills its list in reverse as it walks the
            // declarations, so it reads `static-quantized-quals round-down-quantized` -- the order
            // those two are declared in. Sorting it here made the generated declarations hold the
            // alphabetical order and the port print a message the reference never writes.
            final List<String> mutex = new ArrayList<>(definition.getMutexTargetList());
            final Object plugin = definition.getDescriptorForControllingPlugin();
            System.out.printf(
                    "def\t%s\t%d\t%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%d|%d|%s|%s|%s|%s|%s|%s|%s%n",
                    tool, i,
                    definition.getLongName(),
                    String.join(",", shorts),
                    definition.isOptional() ? "optional" : "required",
                    definition.isCollection() ? "collection" : "scalar",
                    escape(String.valueOf(definition.getDefaultValueAsString())),
                    type,
                    primitive ? "primitive" : "boxed",
                    definition.isFlag() ? "flag" : "valued",
                    definition.isHidden() ? "hidden" : "printed",
                    definition.isAdvanced() ? "advanced" : "plain",
                    definition.isCommon() ? "common" : "own",
                    definition.getMinElements(),
                    definition.getMaxElements(),
                    String.valueOf(definition.getMinValue()),
                    String.valueOf(definition.getMaxValue()),
                    String.valueOf(definition.getMinRecommendedValue()),
                    String.valueOf(definition.getMaxRecommendedValue()),
                    mutex.isEmpty() ? "none" : String.join(",", mutex),
                    plugin == null ? "none" : plugin.getClass().getSimpleName(),
                    typeClass);
            // The documentation is a line of its own: it is prose, and prose carries the pipe the
            // line above uses as its separator.
            System.out.printf("doc\t%s\t%d\t%s%n", tool, i,
                    escape(String.valueOf(definition.getDocString())));
            if (instance != null && !seen.contains(definition.getLongName())) {
                System.out.printf("only-on-the-tool\t%s\t%s%n", tool,
                        definition.getLongName());
            }
        }
    }

    /** The dump's own escaping, since this harness has no shared helper. */
    static String escape(final String text) {
        return text.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
    }

    /** What the parser makes of one command line, before the tool runs at all. */
    static void parse(final String tool, final String label, final String[] argv) {
        final Object target = switch (tool) {
            case "CountReads" -> new org.broadinstitute.hellbender.tools.CountReads();
            case "CountVariants" -> new org.broadinstitute.hellbender.tools.walkers.CountVariants();
            case "IndexFeatureFile" -> new org.broadinstitute.hellbender.tools.IndexFeatureFile();
            case "CompareBaseQualities" ->
                    new org.broadinstitute.hellbender.tools.validation.CompareBaseQualities();
            default -> new org.broadinstitute.hellbender.tools.PrintReads();
        };
        String result;
        try {
            final PrintStream sink = new PrintStream(new ByteArrayOutputStream());
            result = ((org.broadinstitute.hellbender.cmdline.CommandLineProgram) target)
                    .getCommandLineParser().parseArguments(sink, argv) ? "ok" : "not-parsed";
        } catch (final Exception | AssertionError e) {
            result = "E:" + e.getClass().getName() + ":" + e.getMessage();
        }
        System.out.printf("parse\t%s\t%s\t%s%n", tool, label,
                escape(result));
    }
}

/*
 * PathSeqBuildReferenceTaxonomy's output, taken from the reference.
 *
 * The NCBI taxonomy dump and one or two accession catalogs turned into the tree PathSeq scores
 * against, trimmed to the taxa the reference actually holds. The output is a Kryo serialisation, so
 * what is printed here is the tree read back: every node with its parent, rank, name and length,
 * and the accession-to-taxon map.
 *
 * Thirteen behaviours this is built to catch.
 *
 *   - THE MAP THE TOOL WRITES IS KEYED BY THE WHOLE CONTIG NAME, NOT BY THE ACCESSION, whatever
 *     the parameter is called: `addReferenceAccessionToTaxon` is handed the record's name, so an
 *     entry looks up as `ref|NC_BACT.1|` rather than as `NC_BACT.1`, and the accession is only ever
 *     the key the catalog is searched with;
 *   - A REFERENCE NAME IS READ FOR taxid| FIRST AND ref| SECOND: the scan stops at the first taxid
 *     and a name carrying both gives the taxid, while a name carrying neither falls back to the
 *     FIRST WORD of the first bar-delimited token;
 *   - A TAXON ID THAT IS NOT AN INTEGER IS A UserException wherever it is read, the reference name
 *     included;
 *   - THE TWO CATALOG FORMATS ARE JUST TWO COLUMN PAIRS, taxon 0 and accession 2 for RefSeq,
 *     taxon 6 and accession 1 for GenBank;
 *   - THE CATALOG IS SPLIT WITH A LIMIT, so everything past the last column it needs stays in one
 *     field and is never looked at;
 *   - A LINE WITH TOO FEW COLUMNS IS REFUSED WITH A MESSAGE THAT SAYS GenBank WHATEVER THE FORMAT,
 *     and the line number it names counts from one;
 *   - AN EMPTY LINE ENDS THE CATALOG SILENTLY, the read loop stopping on it, so everything after a
 *     blank line is invisible and its accessions are merely reported as not found;
 *   - ONLY "scientific name" ROWS OF names.dmp COUNT, and a taxon named there but absent from the
 *     reference is still added to the map;
 *   - nodes.dmp GIVES EVERY NODE ITS RANK AND PARENT, leaves the root's parent unset, and names a
 *     node it has never seen tax_<id>;
 *   - A NODE MISSING A NAME, A PARENT OR A RANK IS DROPPED FROM THE TREE, and so is one that cannot
 *     be reached from the root;
 *   - THE TREE IS THEN TRIMMED TO THE TAXA THAT HOLD ACCESSIONS AND THEIR ANCESTORS, so a branch of
 *     the taxonomy disappears when nothing in the reference sits on it, and a node's length is the
 *     TOTAL of every contig that landed on it;
 *   - --min-non-virus-contig-length DROPS SHORT CONTIGS FROM THE MAP BUT NOT FROM THE TREE, whose
 *     lengths still count them, and never drops one whose path contains the virus node 10239;
 *   - AND NEITHER CATALOG AT ALL IS A UserException before anything is read.
 *
 * Output:
 *
 *     fixture\t<name>=<the whole file, escaped>
 *     tree\t<label>\t<taxid>\t<name>\t<parent>\t<rank>\t<length>
 *     accession\t<label>\t<accession>=<taxid>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: PathSeqBuildReferenceTaxonomyDump
 */

import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveOutputStream;
import org.broadinstitute.hellbender.tools.spark.pathseq.PSScorer;
import org.broadinstitute.hellbender.tools.spark.pathseq.PSTaxonomyDatabase;
import org.broadinstitute.hellbender.tools.spark.pathseq.PSTree;
import org.broadinstitute.hellbender.tools.spark.pathseq.PathSeqBuildReferenceTaxonomy;

import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.zip.GZIPOutputStream;

public class PathSeqBuildReferenceTaxonomyDump {

    /** One contig of the reference: the name it carries, and how long it is. */
    record Contig(String name, int length) {}

    static Contig contig(final String name, final int length) {
        return new Contig(name, length);
    }

    /**
     * The taxonomy the fixture declares: root, two superkingdoms, and one species under each.
     * The virus node is 10239, which is the number the length filter is exempted by.
     */
    static final String NAMES = String.join("\n",
            "1\t|\troot\t|\t\t|\tscientific name\t|",
            "1\t|\tall\t|\t\t|\tsynonym\t|",
            "2\t|\tBacteria\t|\t\t|\tscientific name\t|",
            "10239\t|\tViruses\t|\t\t|\tscientific name\t|",
            "562\t|\tEscherichia coli\t|\t\t|\tscientific name\t|",
            "11234\t|\tMeasles morbillivirus\t|\t\t|\tscientific name\t|",
            "9606\t|\tHomo sapiens\t|\t\t|\tscientific name\t|",
            "40674\t|\tMammalia\t|\t\t|\tscientific name\t|") + "\n";

    static final String NODES = String.join("\n",
            "1\t|\t1\t|\tno rank\t|",
            "2\t|\t1\t|\tsuperkingdom\t|",
            "10239\t|\t1\t|\tsuperkingdom\t|",
            "562\t|\t2\t|\tspecies\t|",
            "11234\t|\t10239\t|\tspecies\t|",
            "40674\t|\t1\t|\tclass\t|",
            "9606\t|\t40674\t|\tspecies\t|") + "\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("pathseq-taxonomy-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# PathSeqBuildReferenceTaxonomyDump: the NCBI taxonomy trimmed to a reference");

        final Path taxdump = dir.resolve("taxdump.tar.gz");
        writeTaxdump(taxdump);
        System.out.printf("fixture\tnames.dmp=%s%n", ReferenceQueryDump.escape(NAMES));
        System.out.printf("fixture\tnodes.dmp=%s%n", ReferenceQueryDump.escape(NODES));

        // The reference, whose names carry every form the parser knows.
        final List<Contig> contigs = List.of(
                contig("ref|NC_VIRUS.1| a virus", 300),
                contig("ref|NC_BACT.1| a bacterium", 1000),
                contig("ref|NC_SHORT.1| a short bacterium", 100),
                contig("taxid|562| named by its taxon", 800),
                contig("ACC_PLAIN.1 named by its first word", 900),
                contig("gi|9|ref|NC_BOTH.1|taxid|11234|", 700));
        final Path fasta = writeReference(dir, "reference", contigs);

        // The RefSeq catalog: taxon in column 0, accession in column 2.
        final String refseq = String.join("\n",
                "11234\tsomething\tNC_VIRUS.1\tmore\tcolumns\tignored",
                "562\tsomething\tNC_BACT.1\tmore\tcolumns\tignored",
                "562\tsomething\tNC_SHORT.1\tmore\tcolumns\tignored") + "\n";
        final Path refseqPath = writeGz(dir, "refseq.catalog.gz", refseq);
        System.out.printf("fixture\trefseq.catalog=%s%n", ReferenceQueryDump.escape(refseq));

        // The GenBank catalog: accession in column 1, taxon in column 6.
        final String genbank = String.join("\n",
                "a\tACC_PLAIN.1\tc\td\te\tf\t9606\th") + "\n";
        final Path genbankPath = writeGz(dir, "genbank.catalog.gz", genbank);
        System.out.printf("fixture\tgenbank.catalog=%s%n", ReferenceQueryDump.escape(genbank));

        run(dir, "both-catalogs", fasta, taxdump, refseqPath, genbankPath, 0);
        run(dir, "min-length-500", fasta, taxdump, refseqPath, genbankPath, 500);
        run(dir, "refseq-only", fasta, taxdump, refseqPath, null, 0);
        run(dir, "genbank-only", fasta, taxdump, null, genbankPath, 0);
        run(dir, "no-catalog", fasta, taxdump, null, null, 0);

        // A catalog whose blank line hides everything after it.
        final String truncated = "11234\tsomething\tNC_VIRUS.1\n\n562\tsomething\tNC_BACT.1\n";
        final Path truncatedPath = writeGz(dir, "truncated.catalog.gz", truncated);
        System.out.printf("fixture\ttruncated.catalog=%s%n", ReferenceQueryDump.escape(truncated));
        run(dir, "blank-line-truncates", fasta, taxdump, truncatedPath, null, 0);

        // A catalog line with too few columns, refused with a message that says GenBank.
        final String narrow = "11234\tsomething\n";
        final Path narrowPath = writeGz(dir, "narrow.catalog.gz", narrow);
        System.out.printf("fixture\tnarrow.catalog=%s%n", ReferenceQueryDump.escape(narrow));
        run(dir, "narrow-catalog", fasta, taxdump, narrowPath, null, 0);

        // A reference name whose taxon id is not a number.
        final Path badTaxon = writeReference(dir, "bad-taxon",
                List.of(contig("taxid|abc| not a number", 100)));
        run(dir, "bad-taxon-id", badTaxon, taxdump, refseqPath, null, 0);

        // A reference whose accessions are in no catalog at all.
        final Path unknown = writeReference(dir, "unknown",
                List.of(contig("ref|NC_NOWHERE.1| in no catalog", 100)));
        run(dir, "no-relevant-taxa", unknown, taxdump, refseqPath, null, 0);
    }

    static Path writeReference(final Path dir, final String label, final List<Contig> contigs)
            throws Exception {
        final Path fasta = dir.resolve(label + ".fasta");
        final StringBuilder text = new StringBuilder();
        for (final Contig contig : contigs) {
            text.append('>').append(contig.name()).append('\n');
            final StringBuilder bases = new StringBuilder();
            for (int index = 0; index < contig.length(); index++) {
                bases.append("ACGT".charAt(index % 4));
            }
            for (int index = 0; index < bases.length(); index += 60) {
                text.append(bases, index, Math.min(index + 60, bases.length())).append('\n');
            }
        }
        Files.writeString(fasta, text.toString(), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve(label + ".dict")});
        return fasta;
    }

    static Path writeGz(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        try (final OutputStream out = new GZIPOutputStream(Files.newOutputStream(path))) {
            out.write(text.getBytes(StandardCharsets.UTF_8));
        }
        return path;
    }

    static void writeTaxdump(final Path path) throws Exception {
        try (final TarArchiveOutputStream tar =
                     new TarArchiveOutputStream(new GZIPOutputStream(Files.newOutputStream(path)))) {
            addEntry(tar, "names.dmp", NAMES);
            addEntry(tar, "nodes.dmp", NODES);
        }
    }

    static void addEntry(final TarArchiveOutputStream tar, final String name, final String text)
            throws Exception {
        final byte[] bytes = text.getBytes(StandardCharsets.UTF_8);
        final TarArchiveEntry entry = new TarArchiveEntry(name);
        entry.setSize(bytes.length);
        tar.putArchiveEntry(entry);
        tar.write(bytes);
        tar.closeArchiveEntry();
    }

    static void run(final Path dir, final String label, final Path fasta, final Path taxdump,
                    final Path refseq, final Path genbank, final int minLength) throws Exception {
        final Path out = dir.resolve(label + ".db");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-R", fasta.toString(),
                "-O", out.toString(),
                "--tax-dump", taxdump.toString(),
                "--min-non-virus-contig-length", Integer.toString(minLength)));
        if (refseq != null) {
            argv.addAll(Arrays.asList("--refseq-catalog", refseq.toString()));
        }
        if (genbank != null) {
            argv.addAll(Arrays.asList("--genbank-catalog", genbank.toString()));
        }
        try {
            new PathSeqBuildReferenceTaxonomy().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            return;
        }
        final PSTaxonomyDatabase database = PSScorer.readTaxonomyDatabase(out.toString());
        final PSTree tree = database.tree;
        for (final int id : new java.util.TreeSet<>(tree.getNodeIDs())) {
            System.out.printf("tree\t%s\t%d\t%s\t%d\t%s\t%d%n", label, id, tree.getNameOf(id),
                    tree.getParentOf(id), tree.getRankOf(id), tree.getLengthOf(id));
        }
        for (final Map.Entry<String, Integer> entry
                : new TreeMap<>(database.accessionToTaxId).entrySet()) {
            System.out.printf("accession\t%s\t%s=%d%n", label, entry.getKey(), entry.getValue());
        }
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}

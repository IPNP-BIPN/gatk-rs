/*
 * RampedHaplotypeCaller's ramps, taken from the reference.
 *
 * A ramp is a state file: an off ramp halts the caller at one step and writes a zip, an on ramp
 * restarts it from that zip. What is decidable without running the caller is the file's own shape,
 * the two orderings the ramps compare their contents under, and what the comparison refuses.
 *
 * Nine behaviours this is built to catch.
 *
 *   - A ZIP ENTRY IS NAMED `<contig>-<start>-<end>/<name>` when it belongs to a region and by its
 *     bare name when it does not, so the region prefix is a directory made out of coordinates;
 *   - `info.json` IS WRITTEN LAST and with a two-space indent, whatever was added to it;
 *   - AN INFO KEY IS A DOTTED PATH, and every missing level of it is created as an object;
 *   - THE HAPLOTYPE TABLE IS A CSV WITH A FIXED HEADER, whose reference column is 1 or 0 and whose
 *     score is `Double.toString`;
 *   - A REFERENCE HAPLOTYPE SORTS LAST, because the comparator subtracts the reference flags and a
 *     true is one;
 *   - THE HAPLOTYPE COMPARATOR FALLS THROUGH SIX KEYS and compares the score by SIGN rather than by
 *     difference, so two scores a hair apart are still ordered;
 *   - THE READ COMPARATOR STARTS FROM THE STRAND, so a reverse read sorts after a forward one
 *     whatever their positions;
 *   - THE VERIFICATION REFUSES A SIZE MISMATCH BEFORE IT COMPARES ANYTHING, and otherwise names the
 *     index it failed on;
 *   - AND THE BAM INDEX PATH IS BUILT BY `String.replace`, which replaces EVERY `.bam` in the path
 *     rather than the last one.
 *
 * Output:
 *
 *     entry\t<label>=<the zip entry names, comma separated, in write order>
 *     content\t<label>=<one entry's bytes, escaped>
 *     order\t<label>=<the sorted ids, comma separated>
 *     compare\t<label>=<the comparator's sign>
 *     path\t<label>=<a derived path>
 *     name\t<label>=<a derived name>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: RampedHaplotypeCallerDump
 */

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ramps.OffRampBase;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ramps.RampUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

public class RampedHaplotypeCallerDump {

    /** The off ramp's own methods are protected, so a subclass is what reaches them. */
    static class Ramp extends OffRampBase {
        Ramp(final String filename) throws Exception {
            super(filename);
        }

        void put(final htsjdk.samtools.util.Locatable loc, final String name,
                 final byte[] bytes) throws Exception {
            addEntry(loc, name, bytes);
        }

        void putHaplotypes(final htsjdk.samtools.util.Locatable loc,
                           final String name, final List<Haplotype> value) throws Exception {
            addHaplotypes(loc, name, value);
        }

        void putInfo(final String name, final Object value) {
            addInfo(info, name, value);
        }

        String locSuffix(final htsjdk.samtools.util.Locatable loc) {
            return getLocFilenameSuffix(loc);
        }

        String suppName(final String name, final boolean supp) {
            return getReadSuppName(name, supp);
        }

        Path indexPath(final Path bam) {
            return getBamIndexPath(bam);
        }
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 199980))));
        return header;
    }

    static Haplotype haplotype(final String bases, final boolean reference, final int start,
                               final int end, final String cigar, final double score) {
        final Haplotype h = new Haplotype(bases.getBytes(), reference,
                new SimpleInterval("chr1", start, end), TextCigarCodec.decode(cigar));
        h.setScore(score);
        return h;
    }

    static GATKRead read(final SAMFileHeader header, final String name, final int start,
                         final String cigar, final String bases, final boolean reverse) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setMappingQuality(60);
        record.setReadBases(bases.getBytes());
        final byte[] qualities = new byte[bases.length()];
        Arrays.fill(qualities, (byte) 30);
        record.setBaseQualities(qualities);
        record.setReadNegativeStrandFlag(reverse);
        return new SAMRecordToGATKReadAdapter(record);
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("ramped-haplotype-caller-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# RampedHaplotypeCallerDump: the shape of a ramp state file, and the "
                + "orderings it compares under");

        final SAMFileHeader header = header();
        final SimpleInterval loc = new SimpleInterval("chr1", 1000, 2000);

        // The three haplotypes the table is written from: one reference and two not, one of the
        // two carrying a score a hair from the other's.
        final List<Haplotype> haplotypes = new ArrayList<>(List.of(
                haplotype("ACGTACGT", false, 1000, 1007, "8M", -12.5),
                haplotype("ACGTTCGT", true, 1000, 1007, "8M", 0.0),
                haplotype("ACGTACGA", false, 1000, 1007, "8M", -12.5000000001)));

        // A ramp with one region and one root-level entry.
        final Path rampPath = dir.resolve("ramp.zip");
        final Ramp ramp = new Ramp(rampPath.toString());
        ramp.put(loc, "reads.txt", "one\ntwo\n".getBytes());
        ramp.putHaplotypes(loc, "haplotypes.csv", haplotypes);
        ramp.put(null, "root.txt", "at the root\n".getBytes());
        // A dotted key, whose missing levels are created.
        ramp.putInfo("assembly.region.count", 3);
        ramp.putInfo("assembly.region.name", "first");
        ramp.putInfo("flat", true);
        ramp.close();

        // The entry names, in the order they were written.
        final List<String> names = new ArrayList<>();
        try (final ZipFile zip = new ZipFile(rampPath.toFile())) {
            for (final ZipEntry entry : Collections.list(zip.entries())) {
                names.add(entry.getName());
            }
            System.out.printf("entry\tramp=%s%n", String.join(",", names));
            for (final String name : new String[] {"chr1-1000-2000/haplotypes.csv", "info.json"}) {
                final ZipEntry entry = zip.getEntry(name);
                System.out.printf("content\t%s=%s%n", name, ReferenceQueryDump.escape(
                        new String(zip.getInputStream(entry).readAllBytes())));
            }
        }

        // The two orderings.
        final List<Haplotype> sortedHaplotypes = new ArrayList<>(haplotypes);
        sortedHaplotypes.sort(new RampUtils.HaplotypeComparator());
        final List<String> haplotypeOrder = new ArrayList<>();
        for (final Haplotype h : sortedHaplotypes) {
            haplotypeOrder.add(h.getBaseString());
        }
        System.out.printf("order\thaplotypes=%s%n", String.join(",", haplotypeOrder));
        // The score comparison is by SIGN, so a difference of 1e-10 still orders.
        System.out.printf("compare\tscore-hair=%d%n", Integer.signum(
                new RampUtils.HaplotypeComparator().compare(haplotypes.get(0), haplotypes.get(2))));
        System.out.printf("compare\treference-last=%d%n", Integer.signum(
                new RampUtils.HaplotypeComparator().compare(haplotypes.get(1), haplotypes.get(0))));

        final List<GATKRead> reads = new ArrayList<>(List.of(
                read(header, "r-reverse-early", 1000, "8M", "ACGTACGT", true),
                read(header, "r-forward-late", 5000, "8M", "ACGTACGT", false),
                read(header, "r-forward-early", 1000, "8M", "ACGTACGT", false)));
        final List<GATKRead> sortedReads = new ArrayList<>(reads);
        sortedReads.sort(new RampUtils.GATKReadComparator());
        final List<String> readOrder = new ArrayList<>();
        for (final GATKRead r : sortedReads) {
            readOrder.add(r.getName());
        }
        System.out.printf("order\treads=%s%n", String.join(",", readOrder));
        System.out.printf("compare\tstrand-first=%d%n", Integer.signum(
                new RampUtils.GATKReadComparator().compare(reads.get(0), reads.get(1))));

        // The verification, and what it refuses.
        verify("haplotypes-same", () -> RampUtils.compareHaplotypes(haplotypes,
                new ArrayList<>(haplotypes)));
        verify("haplotypes-size", () -> RampUtils.compareHaplotypes(haplotypes,
                haplotypes.subList(0, 2)));
        final List<Haplotype> changed = new ArrayList<>(haplotypes.subList(0, 2));
        changed.add(haplotype("TTTTTTTT", false, 1000, 1007, "8M", -12.5000000001));
        verify("haplotypes-different", () -> RampUtils.compareHaplotypes(haplotypes, changed));
        verify("reads-same", () -> RampUtils.compareReads(reads, new ArrayList<>(reads)));
        verify("reads-size", () -> RampUtils.compareReads(reads, reads.subList(0, 2)));

        // The derived names and paths.
        System.out.printf("name\tloc-suffix=%s%n", ramp.locSuffix(loc));
        System.out.printf("name\tsupp-true=%s%n", ramp.suppName("readname", true));
        System.out.printf("name\tsupp-false=%s%n", ramp.suppName("readname", false));
        // `String.replace` replaces EVERY occurrence, so a path with `.bam` in a directory name
        // has that replaced too.
        System.out.printf("path\tplain=%s%n",
                ramp.indexPath(Path.of("/x/reads.bam")).toString());
        System.out.printf("path\ttwice=%s%n",
                ramp.indexPath(Path.of("/x/run.bam.d/reads.bam")).toString());
        System.out.printf("path\tnone=%s%n",
                ramp.indexPath(Path.of("/x/reads.cram")).toString());
    }

    interface Check {
        void run();
    }

    static void verify(final String label, final Check check) {
        try {
            check.run();
            System.out.printf("error\t%s\tNONE:accepted%n", label);
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(cause.getMessage())));
        }
    }
}

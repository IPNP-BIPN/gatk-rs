//! `LocalAssembler`: the graph it writes and the sequences it reads off it.
//!
//! The de Bruijn assembly itself is not ported. What is ported is everything the two output files
//! are made of: the kmer size the overlaps rest on, a contig's orientation and its reverse
//! complement, the sequence a traversal spells out, and the shape of every line in the GFA and
//! the FASTA.
//!
//! Ported from `org.broadinstitute.hellbender.tools.LocalAssembler`.

/// `LocalAssembler.Kmer.KSIZE`, which must be an odd number under thirty-two.
pub const KMER_SIZE: usize = 31;
/// `LocalAssembler.MIN_THIN_OBS_DEFAULT`: how often a contig must be seen to survive.
pub const MIN_THIN_OBS_DEFAULT: i32 = 4;
/// The GFA version the header names.
pub const GFA_HEADER: &str = "H\tVN:Z:2.0";

/// The reverse complement of a sequence, which is what a contig's other orientation is.
pub fn reverse_complement(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other,
        })
        .collect()
}

/// One step of a traversal: a contig and whether it is taken reverse-complemented.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Step {
    pub contig: String,
    pub reverse_complemented: bool,
}

impl Step {
    /// `Contig.toString`, which appends `RC` to the reverse-complemented orientation and nothing
    /// to the forward one.
    pub fn name(&self) -> String {
        if self.reverse_complemented {
            format!("{}RC", self.contig)
        } else {
            self.contig.clone()
        }
    }

    /// `Contig.toRef`, which the GFA uses instead: a `+` or a `-` rather than a suffix.
    pub fn reference(&self) -> String {
        format!(
            "{}{}",
            self.contig,
            if self.reverse_complemented { "-" } else { "+" }
        )
    }
}

/// `Traversal.toString`: the steps joined by `+`, so a path can read `c1+c3RC+c4RC`.
pub fn traversal_name(steps: &[Step]) -> String {
    steps.iter().map(Step::name).collect::<Vec<_>>().join("+")
}

/// `Traversal.getSequence`: the first contig's leading `KSIZE - 1` bases, then every contig from
/// its `KSIZE - 1`st base on.
///
/// Adjacent contigs therefore overlap by THIRTY bases and each contributes its length less thirty,
/// which is why three contigs of 59, 61 and 60 bases spell out 120 and not 180.
pub fn traversal_sequence(steps: &[Step], sequence_of: impl Fn(&str) -> Option<String>) -> String {
    if steps.is_empty() {
        return String::new();
    }
    let oriented = |step: &Step| {
        let sequence = sequence_of(&step.contig)?;
        Some(if step.reverse_complemented {
            reverse_complement(&sequence)
        } else {
            sequence
        })
    };
    let Some(first) = oriented(&steps[0]) else {
        return String::new();
    };
    let mut sequence = first[..KMER_SIZE - 1].to_string();
    for step in steps {
        let Some(contig) = oriented(step) else {
            return String::new();
        };
        sequence.push_str(&contig[KMER_SIZE - 1..]);
    }
    sequence
}

/// `Traversal.getSequenceLength`, which counts KMERS and adds the overlap back once.
pub fn traversal_sequence_length(kmer_counts: &[usize]) -> usize {
    kmer_counts.iter().sum::<usize>() + KMER_SIZE - 1
}

/// `writeTraversals`' header line: the assembly name, an underscore, `t` and a ONE-based counter,
/// then a space and the traversal's own name.
pub fn fasta_header(assembly_name: Option<&str>, traversal_number: usize, name: &str) -> String {
    match assembly_name {
        Some(assembly) => format!(">{assembly}_t{traversal_number} {name}"),
        None => format!(">t{traversal_number} {name}"),
    }
}

/// `writeContig`'s `S` line: the id, the LENGTH, the sequence and three observation counts.
pub fn gfa_segment(
    id: &str,
    sequence: &str,
    max_observations: i32,
    first_observations: i32,
    last_observations: i32,
) -> String {
    format!(
        "S\t{id}\t{}\t{sequence}\tMO:i:{max_observations}\tFO:i:{first_observations}\tLO:i:{last_observations}",
        sequence.len()
    )
}

/// `writeEdge`'s `E` line, whose two coordinates are the overlap's place in each contig.
///
/// The first contig's overlap runs from its length less `KSIZE - 1` to its end, written with a
/// `$` to say so, and the second's runs from nought to `KSIZE - 1`. The CIGAR is that same length.
pub fn gfa_edge(from: &Step, to: &Step, from_length: usize) -> String {
    format!(
        "E\t*\t{}\t{}\t{}\t{from_length}$\t0\t{}\t{}M",
        from.reference(),
        to.reference(),
        from_length - KMER_SIZE + 1,
        KMER_SIZE - 1,
        KMER_SIZE - 1
    )
}

/// `writePaths`' `O` line: the steps by reference, separated by SPACES rather than by `+`.
pub fn gfa_path(steps: &[Step]) -> String {
    format!(
        "O\t*\t{}",
        steps
            .iter()
            .map(Step::reference)
            .collect::<Vec<_>>()
            .join(" ")
    )
}

/// A traversal parsed back out of a name like `c1+c3RC+c4RC`.
pub fn parse_traversal_name(name: &str) -> Vec<Step> {
    name.split('+')
        .filter(|part| !part.is_empty())
        .map(|part| match part.strip_suffix("RC") {
            Some(contig) => Step {
                contig: contig.to_string(),
                reverse_complemented: true,
            },
            None => Step {
                contig: part.to_string(),
                reverse_complemented: false,
            },
        })
        .collect()
}

/// A traversal parsed back out of a GFA path's references, like `c1+ c3- c4-`.
pub fn parse_gfa_path(references: &str) -> Vec<Step> {
    references
        .split_whitespace()
        .map(|part| Step {
            contig: part[..part.len() - 1].to_string(),
            reverse_complemented: part.ends_with('-'),
        })
        .collect()
}

/// The reverse of a traversal: the steps backwards, each in its other orientation.
pub fn reverse_traversal(steps: &[Step]) -> Vec<Step> {
    steps
        .iter()
        .rev()
        .map(|step| Step {
            contig: step.contig.clone(),
            reverse_complemented: !step.reverse_complemented,
        })
        .collect()
}

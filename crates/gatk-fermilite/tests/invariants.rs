//! What must hold of an assembly whatever the bytes are. Goldens are frozen from x86-64 CI by the
//! repository's owner, never here: the comparison with the bundled library is in
//! `tools/fermilite-conformance`.

use gatk_fermilite::{
    assemble, write_dump, Contig, FmlError, FmlOpt, MAG_F_AGGRESSIVE, MAG_F_POPOPEN,
};

struct Xs(u64);

impl Xs {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn genome(&mut self, n: usize) -> Vec<u8> {
        (0..n)
            .map(|_| b"ACGT"[(self.next() >> 40) as usize & 3])
            .collect()
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            c => *c,
        })
        .collect()
}

fn reads_of(g: &[u8], n: usize, l: usize, r: &mut Xs) -> Vec<(Vec<u8>, Vec<u8>)> {
    (0..n)
        .map(|i| {
            let p = (r.next() >> 20) as usize % (g.len() - l + 1);
            let s = g[p..p + l].to_vec();
            let s = if i & 1 == 0 { s } else { revcomp(&s) };
            let q = vec![b'I'; l];
            (s, q)
        })
        .collect()
}

fn run(flag: Option<i32>, reads: &[(Vec<u8>, Vec<u8>)]) -> Result<Vec<Contig>, FmlError> {
    let refs: Vec<(&[u8], &[u8])> = reads
        .iter()
        .map(|(a, b)| (a.as_slice(), b.as_slice()))
        .collect();
    let opt = match flag {
        Some(f) => FmlOpt::default().with_cleaning_flag(f),
        None => FmlOpt::default(),
    };
    assemble(&opt, &refs)
}

fn check(contigs: &[Contig]) {
    for (i, c) in contigs.iter().enumerate() {
        assert_eq!(c.seq.len(), c.cov.len(), "contig {i}");
        assert!(c.seq.iter().all(|b| b"ACGTN".contains(b)));
        assert!(c.cov.iter().all(|&b| (33..=126).contains(&b)));
        assert!(c.nsr >= 0);
        for o in &c.ovlp {
            let target = &contigs[o.target()];
            assert!(o.overlap_len() > 0);
            // an edge in the graph is there from both of its ends
            let back = target.ovlp.iter().filter(|b| {
                b.target() == i
                    && b.is_rc() != o.is_target_rc()
                    && b.is_target_rc() != o.is_rc()
                    && b.overlap_len() == o.overlap_len()
            });
            assert!(back.count() >= 1, "contig {i} -> {}", o.target());
        }
    }
}

#[test]
fn one_haplotype_makes_one_contig_of_it() {
    let mut r = Xs(0x1234_5678);
    let g = r.genome(1000);
    let reads = reads_of(&g, 400, 100, &mut r);
    let contigs = run(None, &reads).unwrap();
    check(&contigs);
    assert_eq!(contigs.len(), 1);
    let c = &contigs[0];
    let rc = revcomp(&c.seq);
    let has = |h: &[u8], n: &[u8]| h.windows(n.len()).any(|w| w == n);
    assert!(
        has(&g, &c.seq) || has(&g, &rc),
        "not a substring of the genome"
    );
    assert!(c.seq.len() > 900);
}

#[test]
fn a_rerun_is_identical_and_keeps_no_state() {
    let mut r = Xs(99);
    let g = r.genome(600);
    let reads = reads_of(&g, 250, 75, &mut r);
    let h = r.genome(600);
    let other = reads_of(&h, 250, 75, &mut r);
    let first = run(None, &reads).unwrap();
    let _ = run(Some(MAG_F_AGGRESSIVE | MAG_F_POPOPEN), &other);
    assert_eq!(first, run(None, &reads).unwrap());
}

#[test]
fn the_aggressive_flag_is_valid_too() {
    let mut r = Xs(7);
    let g = r.genome(800);
    let mut h = g.clone();
    h[400] = if h[400] == b'A' { b'C' } else { b'A' };
    let mut reads = reads_of(&g, 200, 100, &mut r);
    reads.extend(reads_of(&h, 200, 100, &mut r));
    check(&run(None, &reads).unwrap());
    check(&run(Some(MAG_F_AGGRESSIVE | MAG_F_POPOPEN), &reads).unwrap());
}

#[test]
fn nothing_to_assemble_kills_the_jvm_as_r53_does() {
    assert!(matches!(run(None, &[]), Err(FmlError::Abort(_))));
    let empty = vec![(Vec::new(), Vec::new()); 3];
    assert!(matches!(run(None, &empty), Err(FmlError::Abort(_))));
    let n = vec![(vec![b'N'; 100], vec![b'I'; 100]); 5];
    assert!(matches!(run(None, &n), Err(FmlError::Abort(_))));
}

#[test]
fn a_nul_in_the_quals_is_undefined_and_refused() {
    let mut r = Xs(5);
    let g = r.genome(300);
    let mut reads = reads_of(&g, 50, 75, &mut r);
    reads[3].1[10] = 0;
    assert!(matches!(run(None, &reads), Err(FmlError::Undefined(_))));
}

#[test]
fn a_base_byte_under_the_table_that_is_a_pointer_is_refused() {
    let mut r = Xs(6);
    let g = r.genome(300);
    let mut reads = reads_of(&g, 50, 75, &mut r);
    reads[3].0[10] = 0x90;
    assert!(matches!(run(None, &reads), Err(FmlError::Undefined(_))));
    reads[3].0[10] = 0xE3;
    assert!(matches!(run(None, &reads), Err(FmlError::Undefined(_))));
}

#[test]
fn quals_past_the_buffer_are_javas_exception() {
    let mut r = Xs(8);
    let g = r.genome(300);
    let mut reads = reads_of(&g, 50, 75, &mut r);
    reads[49].1.push(b'I');
    assert!(matches!(
        run(None, &reads),
        Err(FmlError::JavaBufferOverflow)
    ));
}

#[test]
fn the_dump_has_the_harness_format() {
    let contigs = vec![Contig {
        seq: b"ACGT".to_vec(),
        cov: b"!!!!".to_vec(),
        nsr: 2,
        ovlp: Vec::new(),
    }];
    let mut out = Vec::new();
    write_dump(&mut out, 3, "x", 192, 7, &contigs).unwrap();
    assert_eq!(
        out,
        b"A 3 x flag=192 reads=7 contigs=1\nC 0 4 2 ACGT !!!!\n"
    );
}

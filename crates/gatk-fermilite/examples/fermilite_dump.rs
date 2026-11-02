//! The corpus of `tools/fermilite-conformance` through the port, as `FermiLiteDump.java` prints it:
//! `FERMILITE_CORPUS=corpus.bin FERMILITE_DUMP=dump.txt cargo run --release -p gatk-fermilite --example fermilite_dump`.

use gatk_fermilite::{assemble, write_dump, write_dump_abort, FmlError, FmlOpt};
use std::io::{BufWriter, Write};

fn take<'a>(d: &mut &'a [u8], n: usize) -> &'a [u8] {
    let (a, b) = d.split_at(n);
    *d = b;
    a
}

fn take_u32(d: &mut &[u8]) -> u32 {
    u32::from_be_bytes(take(d, 4).try_into().unwrap())
}

fn main() -> std::io::Result<()> {
    let corpus = std::fs::read(std::env::var("FERMILITE_CORPUS").expect("FERMILITE_CORPUS"))?;
    let dump = std::env::var("FERMILITE_DUMP").expect("FERMILITE_DUMP");
    let mut out = BufWriter::new(std::fs::File::create(dump)?);
    let mut d = &corpus[..];
    assert_eq!(take(&mut d, 6), b"FMLC1\n");
    let mut idx = 0u64;
    let mut totals = [0u64; 4];
    while !d.is_empty() {
        let l = take_u32(&mut d) as usize;
        let name = std::str::from_utf8(take(&mut d, l)).unwrap();
        let flag = take_u32(&mut d) as i32;
        let n = take_u32(&mut d) as usize;
        let mut reads = Vec::with_capacity(n);
        for _ in 0..n {
            let l = take_u32(&mut d) as usize;
            let b = take(&mut d, l);
            let l = take_u32(&mut d) as usize;
            let q = take(&mut d, l);
            reads.push((b, q));
        }
        let opt = FmlOpt::default().with_cleaning_flag(flag);
        match assemble(&opt, &reads) {
            Ok(contigs) => {
                totals[0] += 1;
                totals[1] += contigs.len() as u64;
                write_dump(&mut out, idx, name, flag, n, &contigs)?;
            }
            Err(FmlError::Abort(_)) => {
                totals[2] += 1;
                write_dump_abort(&mut out, idx, name, flag, n)?;
            }
            Err(FmlError::JavaBufferOverflow) => {
                writeln!(
                    out,
                    "J {idx} {name} flag={flag} reads={n} BufferOverflowException"
                )?;
            }
            Err(FmlError::Undefined(why)) => {
                totals[3] += 1;
                writeln!(out, "U {idx} {name} flag={flag} reads={n} {why}")?;
            }
        }
        idx += 1;
    }
    out.flush()?;
    eprintln!(
        "{idx} records: {} assembled into {} contigs, {} aborted, {} refused as undefined",
        totals[0], totals[1], totals[2], totals[3]
    );
    Ok(())
}

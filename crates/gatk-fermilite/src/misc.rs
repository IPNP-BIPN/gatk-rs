//! `misc.c`: `fml_assemble` and the stages it strings together.

use crate::bfc::fml_correct_core;
use crate::khash::Hash64;
use crate::mag::{edge_is_del, mag_g_clean, mag_g_merge, Mag};
use crate::rld::Rld;
use crate::rope::{rle_dec1, rle_nptr, Mrope};
use crate::unitig::fml_fmi2mag_core;
use crate::{cvttsd2si, seq_nt6_signed, seq_revcomp6, seq_reverse, Bseq, Contig, FmlOpt, FmlOvlp};

const MAG_MIN_NSR_COEF: f64 = 0.1;

/// `fml_opt_adjust`.
fn fml_opt_adjust(opt: &mut FmlOpt, seqs: &[Bseq]) {
    let mut tot_len = 0u64;
    for s in seqs {
        tot_len += s.l_seq as u64; // compute total length
    }
    let mut log_len = 10;
    while log_len < 32 {
        // compute ceil(log2(tot_len))
        if 1u64 << log_len > tot_len {
            break;
        }
        log_len += 1;
    }
    if opt.ec_k == 0 {
        opt.ec_k = (log_len + 12) / 2;
    }
    if opt.ec_k % 2 == 0 {
        opt.ec_k += 1;
    }
    opt.mag_opt.min_elen = cvttsd2si(tot_len as f64 / seqs.len() as i32 as f64 * 2.5 + 0.499);
}

/// `is_rev_same`.
fn is_rev_same(s: &[u8]) -> bool {
    let l = s.len();
    if l & 1 != 0 {
        return false;
    }
    for i in 0..l >> 1 {
        if s[i] as i32 + s[l - 1 - i] as i32 != 5 {
            return false;
        }
    }
    true
}

/// `fml_fmi_gen`, single-threaded.
fn fml_fmi_gen(seqs: Vec<Bseq>) -> Rld {
    let mut mr = Mrope::new();
    let mut str_: Vec<u8> = Vec::new();
    for mut s in seqs {
        if s.l_seq == 0 {
            continue;
        }
        for c in s.seq.iter_mut() {
            *c = seq_nt6_signed(*c);
        }
        if s.seq.contains(&5) {
            continue;
        }
        if is_rev_same(&s.seq) {
            s.seq.pop();
        }
        seq_reverse(&mut s.seq);
        str_.extend_from_slice(&s.seq);
        str_.push(0);
        seq_revcomp6(&mut s.seq);
        str_.extend_from_slice(&s.seq);
        str_.push(0);
    }
    mr.insert_multi(&str_);
    drop(str_);

    let (mut e, mut di) = Rld::new();
    mr.for_each_block(|block| {
        let mut q = 2;
        let end = 2 + rle_nptr(block);
        while q < end {
            let (c, l) = rle_dec1(block, &mut q);
            e.enc(&mut di, l, c as u8);
        }
    });
    e.enc_finish(&mut di);
    e
}

/// `fml_mag_clean`; `mag_g_trim_open` returns at once with `trim_len` 0.
fn fml_mag_clean(opt: &FmlOpt, g: &mut Mag) {
    let mut o = opt.mag_opt.clone();
    o.min_merge_len = opt.min_merge_len;
    mag_g_merge(g, true, opt.min_merge_len);
    mag_g_clean(g, &o);
}

/// `fml_mag2utg`.
fn fml_mag2utg(g: Mag) -> Vec<Contig> {
    let mut h = Hash64::new();
    let mut j = 0u64;
    for p in &g.v {
        if p.len < 0 {
            continue;
        }
        let (k, _) = h.put(p.k[0]);
        h.vals[k as usize] = j << 1;
        let (k, _) = h.put(p.k[1]);
        h.vals[k as usize] = (j << 1) | 1;
        j += 1;
    }
    let mut utg = Vec::with_capacity(j as usize);
    for p in g.v {
        if p.len < 0 {
            continue;
        }
        let mut ovlp = Vec::new();
        for from in 0..2u32 {
            for s in &p.nei[from as usize] {
                if edge_is_del(s) {
                    continue;
                }
                let k = h.get(s.x);
                c_assert!(k != h.end(), "fml_mag2utg: k != kh_end(h)");
                let val = h.vals[k as usize];
                ovlp.push(FmlOvlp {
                    raw_len: (((s.y as u32) & 0x7fff_ffff) | (from << 31)) as i32,
                    raw_id: (((val >> 1) as u32 & 0x7fff_ffff) | ((val & 1) as u32) << 31) as i32,
                });
            }
        }
        let seq = p.seq.iter().map(|&c| b"$ACGTN"[c as usize]).collect();
        utg.push(Contig {
            seq,
            cov: p.cov,
            nsr: p.nsr,
            ovlp,
        });
    }
    utg
}

/// `fml_assemble`.
pub(crate) fn fml_assemble(opt0: &FmlOpt, mut seqs: Vec<Bseq>) -> Vec<Contig> {
    let mut opt = opt0.clone();
    fml_opt_adjust(&mut opt, &seqs);
    if opt.ec_k >= 0 {
        fml_correct_core(&opt, false, &mut seqs);
    }
    let kcov = fml_correct_core(&opt, true, &mut seqs);
    let e = fml_fmi_gen(seqs);
    let mut g = fml_fmi2mag_core(&e, opt.min_asm_ovlp, opt.min_merge_len);
    drop(e);
    let scaled = kcov as f64 * MAG_MIN_NSR_COEF;
    let o = &mut opt.mag_opt;
    o.min_ensr = if o.min_ensr as f64 > scaled {
        o.min_ensr
    } else {
        cvttsd2si(scaled + 0.499)
    };
    o.min_ensr = if o.min_ensr < opt0.max_cnt {
        o.min_ensr
    } else {
        opt0.max_cnt
    };
    o.min_ensr = if o.min_ensr > opt0.min_cnt {
        o.min_ensr
    } else {
        opt0.min_cnt
    };
    o.min_insr = o.min_ensr - 1;
    fml_mag_clean(&opt, &mut g);
    fml_mag2utg(g)
}

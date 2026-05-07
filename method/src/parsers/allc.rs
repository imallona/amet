//! allc / methylpy format.
//!
//! Columns: chr, pos (1-based), strand, context, mc, cov, methylated_flag.
//!
//! Per-strand records: + strand pos = 1-based position of the C; - strand pos = 1-based
//! position of the C on - strand, which is the G at the same position on + strand.
//!
//! 0-based start of the CpG dinucleotide on + strand:
//!   + strand: pos - 1
//!   - strand: pos - 2
//!
//! Both strands map to the same 0-based start; mc and cov are summed across strands.

use crate::MethCall;
use crate::error::{AmetError, Result};
use crate::io::open_read;
use crate::reference::CpgReference;
use std::collections::BTreeMap;
use std::io::BufRead;
use std::path::Path;

pub fn read(path: &Path, reference: &CpgReference) -> Result<Vec<MethCall>> {
    let reader = open_read(path)?;
    let mut merged: BTreeMap<(u32, u64), (u32, u32)> = BTreeMap::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 7 {
            continue;
        }
        let context = fields[3];
        if !context.starts_with("CG") {
            continue;
        }
        let pos_1based: u64 = match fields[1].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };
        let strand = fields[2];
        let mc: u32 = fields[4].parse().map_err(|_| AmetError::Parse {
            file: path.display().to_string(),
            line: i + 1,
            msg: "mc is not an integer".to_string(),
        })?;
        let cov: u32 = fields[5].parse().map_err(|_| AmetError::Parse {
            file: path.display().to_string(),
            line: i + 1,
            msg: "cov is not an integer".to_string(),
        })?;

        let cpg_start_0based = match strand {
            "+" => {
                if pos_1based < 1 {
                    return Err(AmetError::Parse {
                        file: path.display().to_string(),
                        line: i + 1,
                        msg: format!("pos {} cannot be < 1 on + strand", pos_1based),
                    });
                }
                pos_1based - 1
            }
            "-" => {
                if pos_1based < 2 {
                    return Err(AmetError::Parse {
                        file: path.display().to_string(),
                        line: i + 1,
                        msg: format!("pos {} cannot be < 2 on - strand", pos_1based),
                    });
                }
                pos_1based - 2
            }
            other => {
                return Err(AmetError::Parse {
                    file: path.display().to_string(),
                    line: i + 1,
                    msg: format!("unrecognised strand '{}'", other),
                });
            }
        };

        let chrom_id = match reference.chrom_id(fields[0]) {
            Some(c) => c,
            None => continue,
        };

        let entry = merged.entry((chrom_id, cpg_start_0based)).or_insert((0, 0));
        entry.0 += mc;
        entry.1 += cov;
    }

    Ok(merged
        .into_iter()
        .map(|((chrom_id, pos), (m, t))| MethCall {
            chrom_id,
            pos,
            m,
            t,
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn make_ref() -> CpgReference {
        // CpGs at 0-based starts 9, 19, 29 on chr1.
        // C at 1-based pos 10 = 0-based 9; G at 1-based pos 11 is the same CpG.
        CpgReference {
            chrom_names: vec!["chr1".into()],
            chrom_id_of: [("chr1".into(), 0u32)].into_iter().collect(),
            positions: vec![vec![9, 19, 29]],
        }
    }

    fn write(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        write!(f, "{}", content).unwrap();
        f.flush().unwrap();
        f
    }

    #[test]
    fn plus_strand_offset() {
        let f = write("chr1\t10\t+\tCGN\t1\t2\t1\n");
        let r = make_ref();
        let calls = read(f.path(), &r).unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].pos, 9);
    }

    #[test]
    fn minus_strand_offset() {
        let f = write("chr1\t11\t-\tCGN\t1\t2\t1\n");
        let r = make_ref();
        let calls = read(f.path(), &r).unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].pos, 9);
    }

    #[test]
    fn both_strands_merge_to_same_position() {
        let f = write(
            "chr1\t10\t+\tCGN\t1\t2\t1\n\
             chr1\t11\t-\tCGN\t2\t3\t1\n",
        );
        let r = make_ref();
        let calls = read(f.path(), &r).unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].pos, 9);
        assert_eq!(calls[0].m, 3);
        assert_eq!(calls[0].t, 5);
    }

    #[test]
    fn non_cpg_context_filtered() {
        let f = write(
            "chr1\t10\t+\tCHG\t1\t2\t1\n\
             chr1\t10\t+\tCHH\t1\t2\t1\n\
             chr1\t10\t+\tCGN\t1\t2\t1\n",
        );
        let r = make_ref();
        let calls = read(f.path(), &r).unwrap();
        assert_eq!(calls.len(), 1);
    }

    #[test]
    fn includes_unmethylated_calls() {
        // mc = 0 records (unmethylated CpGs) are required for entropy calculation.
        let f = write("chr1\t10\t+\tCGN\t0\t5\t0\n");
        let r = make_ref();
        let calls = read(f.path(), &r).unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].m, 0);
        assert_eq!(calls[0].t, 5);
    }

    #[test]
    fn unknown_chrom_dropped() {
        let f = write(
            "chrX\t10\t+\tCGN\t1\t1\t1\n\
             chr1\t10\t+\tCGN\t1\t1\t1\n",
        );
        let r = make_ref();
        let calls = read(f.path(), &r).unwrap();
        assert_eq!(calls.len(), 1);
    }

    #[test]
    fn unknown_strand_errors() {
        let f = write("chr1\t10\t.\tCGN\t1\t1\t1\n");
        let r = make_ref();
        assert!(read(f.path(), &r).is_err());
    }

    #[test]
    fn multiple_cpgs_sorted() {
        let f = write(
            "chr1\t30\t+\tCGN\t1\t1\t1\n\
             chr1\t10\t+\tCGN\t1\t1\t1\n\
             chr1\t20\t+\tCGN\t1\t1\t1\n",
        );
        let r = make_ref();
        let calls = read(f.path(), &r).unwrap();
        assert_eq!(calls[0].pos, 9);
        assert_eq!(calls[1].pos, 19);
        assert_eq!(calls[2].pos, 29);
    }

    #[test]
    fn three_strand_records_merge_correctly() {
        let f = write(
            "chr1\t10\t+\tCGN\t1\t1\t1\n\
             chr1\t10\t+\tCGN\t0\t1\t0\n\
             chr1\t11\t-\tCGN\t1\t2\t1\n",
        );
        let r = make_ref();
        let calls = read(f.path(), &r).unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].m, 2);
        assert_eq!(calls[0].t, 4);
    }

    #[test]
    fn cgg_and_cgt_contexts_kept() {
        // CG-prefixed contexts (CGA, CGC, CGG, CGT, CGN) are all CpG dinucleotides.
        let f = write(
            "chr1\t10\t+\tCGA\t1\t1\t1\n\
             chr1\t20\t+\tCGT\t1\t1\t1\n\
             chr1\t30\t+\tCGG\t1\t1\t1\n",
        );
        let r = make_ref();
        let calls = read(f.path(), &r).unwrap();
        assert_eq!(calls.len(), 3);
    }
}

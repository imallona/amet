//! Bian scTrioSeq2 singleC format (Bismark-derived per-cytosine table).
//!
//! Header row present:
//!   #Chr  Pos  Ref  Chain  Total  Met  UnMet  MetRate  Ref_context  Type
//!
//! Pos is 1-based; Chain is `+` or `-`; Type is `CpG` / `CHG` / `CHH`.
//! Per-strand records: + strand pos = 1-based C position; - strand pos = 1-based
//! position of the C on - strand, which is the G at the same position on +
//! strand. Both strands map to the same 0-based start of the CpG dinucleotide
//! on + strand:
//!     + strand: pos - 1
//!     - strand: pos - 2
//!
//! mc = Met (col 6); cov = Total (col 5). CHG / CHH rows are dropped.

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
        if fields.len() < 10 {
            continue;
        }
        if fields[9] != "CpG" {
            continue;
        }
        let pos_1based: u64 = match fields[1].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };
        let strand = fields[3];
        let total: u32 = fields[4].parse().map_err(|_| AmetError::Parse {
            file: path.display().to_string(),
            line: i + 1,
            msg: "Total is not an integer".to_string(),
        })?;
        let met: u32 = fields[5].parse().map_err(|_| AmetError::Parse {
            file: path.display().to_string(),
            line: i + 1,
            msg: "Met is not an integer".to_string(),
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
        entry.0 += met;
        entry.1 += total;
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

    const HEADER: &str = "#Chr\tPos\tRef\tChain\tTotal\tMet\tUnMet\tMetRate\tRef_context\tType\n";

    #[test]
    fn header_skipped() {
        let f = write(HEADER);
        let r = make_ref();
        let calls = read(f.path(), &r).unwrap();
        assert!(calls.is_empty());
    }

    #[test]
    fn plus_strand_offset() {
        let body = "chr1\t10\tC\t+\t2\t1\t1\t0.5\tCGN\tCpG\n";
        let f = write(&format!("{HEADER}{body}"));
        let calls = read(f.path(), &make_ref()).unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].pos, 9);
        assert_eq!(calls[0].m, 1);
        assert_eq!(calls[0].t, 2);
    }

    #[test]
    fn minus_strand_offset() {
        let body = "chr1\t11\tG\t-\t3\t2\t1\t0.66\tCGN\tCpG\n";
        let f = write(&format!("{HEADER}{body}"));
        let calls = read(f.path(), &make_ref()).unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].pos, 9);
        assert_eq!(calls[0].m, 2);
        assert_eq!(calls[0].t, 3);
    }

    #[test]
    fn both_strands_merge() {
        let body = "chr1\t10\tC\t+\t2\t1\t1\t0.5\tCGN\tCpG\n\
                    chr1\t11\tG\t-\t3\t2\t1\t0.66\tCGN\tCpG\n";
        let f = write(&format!("{HEADER}{body}"));
        let calls = read(f.path(), &make_ref()).unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].pos, 9);
        assert_eq!(calls[0].m, 3);
        assert_eq!(calls[0].t, 5);
    }

    #[test]
    fn non_cpg_filtered() {
        let body = "chr1\t10\tC\t+\t1\t1\t0\t1\tCAA\tCHH\n\
                    chr1\t10\tC\t+\t1\t1\t0\t1\tCAG\tCHG\n\
                    chr1\t10\tC\t+\t2\t1\t1\t0.5\tCGN\tCpG\n";
        let f = write(&format!("{HEADER}{body}"));
        let calls = read(f.path(), &make_ref()).unwrap();
        assert_eq!(calls.len(), 1);
    }

    #[test]
    fn unmethylated_kept() {
        let body = "chr1\t10\tC\t+\t5\t0\t5\t0\tCGN\tCpG\n";
        let f = write(&format!("{HEADER}{body}"));
        let calls = read(f.path(), &make_ref()).unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].m, 0);
        assert_eq!(calls[0].t, 5);
    }

    #[test]
    fn unknown_strand_errors() {
        let body = "chr1\t10\tC\t.\t1\t1\t0\t1\tCGN\tCpG\n";
        let f = write(&format!("{HEADER}{body}"));
        assert!(read(f.path(), &make_ref()).is_err());
    }

    #[test]
    fn unknown_chrom_dropped() {
        let body = "chrZZ\t10\tC\t+\t1\t1\t0\t1\tCGN\tCpG\n\
                    chr1\t10\tC\t+\t1\t1\t0\t1\tCGN\tCpG\n";
        let f = write(&format!("{HEADER}{body}"));
        let calls = read(f.path(), &make_ref()).unwrap();
        assert_eq!(calls.len(), 1);
    }

    #[test]
    fn three_strand_records_merge() {
        let body = "chr1\t10\tC\t+\t1\t1\t0\t1\tCGN\tCpG\n\
                    chr1\t10\tC\t+\t1\t0\t1\t0\tCGN\tCpG\n\
                    chr1\t11\tG\t-\t2\t1\t1\t0.5\tCGN\tCpG\n";
        let f = write(&format!("{HEADER}{body}"));
        let calls = read(f.path(), &make_ref()).unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].m, 2);
        assert_eq!(calls[0].t, 4);
    }

    #[test]
    fn multiple_cpgs_sorted() {
        let body = "chr1\t30\tC\t+\t1\t1\t0\t1\tCGN\tCpG\n\
                    chr1\t10\tC\t+\t1\t1\t0\t1\tCGN\tCpG\n\
                    chr1\t20\tC\t+\t1\t1\t0\t1\tCGN\tCpG\n";
        let f = write(&format!("{HEADER}{body}"));
        let calls = read(f.path(), &make_ref()).unwrap();
        assert_eq!(calls[0].pos, 9);
        assert_eq!(calls[1].pos, 19);
        assert_eq!(calls[2].pos, 29);
    }
}

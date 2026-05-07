//! scNMT-seq cpg_level format.
//!
//! Columns (header row present): chr, pos, met_reads, nonmet_reads, rate.
//!
//! pos is the 1-based G position of the strand-collapsed CpG dinucleotide.
//! amet uses the 0-based start of the C on the + strand, so the conversion is pos - 2:
//!
//!   genome:           ... C  G  ...
//!   1-based:              N  N+1
//!   amet pos = (N+1) - 2 = N - 1   (the C, 0-based)

use crate::MethCall;
use crate::error::{AmetError, Result};
use crate::io::open_read;
use crate::reference::CpgReference;
use std::io::BufRead;
use std::path::Path;

pub fn read(path: &Path, reference: &CpgReference) -> Result<Vec<MethCall>> {
    let reader = open_read(path)?;
    let mut calls: Vec<MethCall> = Vec::new();
    let mut lines = reader.lines();

    if let Some(first) = lines.next() {
        let first = first?;
        if let Some(parsed) = try_parse_line(&first, reference, path, 1) {
            calls.push(parsed?);
        }
    }

    for (i, line) in lines.enumerate() {
        let line = line?;
        let line_no = i + 2;
        if let Some(parsed) = try_parse_line(&line, reference, path, line_no) {
            calls.push(parsed?);
        }
    }

    calls.sort_by_key(|c| (c.chrom_id, c.pos));
    Ok(calls)
}

fn try_parse_line(
    line: &str,
    reference: &CpgReference,
    path: &Path,
    line_no: usize,
) -> Option<Result<MethCall>> {
    let trimmed = line.trim();
    if trimmed.is_empty() || trimmed.starts_with('#') {
        return None;
    }
    let fields: Vec<&str> = trimmed.split('\t').collect();
    if fields.len() < 4 {
        return None;
    }
    let pos_1based: u64 = fields[1].parse().ok()?;
    let met: u32 = fields[2].parse().ok()?;
    let nonmet: u32 = fields[3].parse().ok()?;

    if pos_1based < 2 {
        return Some(Err(AmetError::Parse {
            file: path.display().to_string(),
            line: line_no,
            msg: format!(
                "scNMT pos {} cannot be < 2 (would underflow on pos-2 conversion)",
                pos_1based
            ),
        }));
    }
    let pos_0based = pos_1based - 2;

    let chrom_id = reference.chrom_id(fields[0])?;
    Some(Ok(MethCall {
        chrom_id,
        pos: pos_0based,
        m: met,
        t: met + nonmet,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn make_ref() -> CpgReference {
        CpgReference {
            chrom_names: vec!["chr1".into()],
            chrom_id_of: [("chr1".into(), 0u32)].into_iter().collect(),
            positions: vec![vec![8, 18, 28]],
        }
    }

    fn write(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        write!(f, "{}", content).unwrap();
        f.flush().unwrap();
        f
    }

    #[test]
    fn coordinate_offset_pos_minus_2() {
        // 1-based G pos 10, 20, 30 -> 0-based C pos 8, 18, 28.
        let f = write(
            "chr\tpos\tmet_reads\tnonmet_reads\trate\n\
             chr1\t10\t1\t0\t1.0\n\
             chr1\t20\t0\t2\t0.0\n\
             chr1\t30\t1\t1\t0.5\n",
        );
        let r = make_ref();
        let calls = read(f.path(), &r).unwrap();
        assert_eq!(calls.len(), 3);
        assert_eq!(calls[0].pos, 8);
        assert_eq!(calls[1].pos, 18);
        assert_eq!(calls[2].pos, 28);
    }

    #[test]
    fn t_equals_m_plus_nonmet() {
        let f = write("chr\tpos\tmet_reads\tnonmet_reads\trate\nchr1\t10\t3\t7\t0.3\n");
        let r = make_ref();
        let calls = read(f.path(), &r).unwrap();
        assert_eq!(calls[0].m, 3);
        assert_eq!(calls[0].t, 10);
    }

    #[test]
    fn header_skipped_when_first_line_non_numeric() {
        let f = write("chr\tpos\tmet_reads\tnonmet_reads\trate\nchr1\t10\t1\t0\t1.0\n");
        let r = make_ref();
        let calls = read(f.path(), &r).unwrap();
        assert_eq!(calls.len(), 1);
    }

    #[test]
    fn no_header_first_line_kept() {
        let f = write("chr1\t10\t1\t0\t1.0\nchr1\t20\t1\t1\t0.5\n");
        let r = make_ref();
        let calls = read(f.path(), &r).unwrap();
        assert_eq!(calls.len(), 2);
    }

    #[test]
    fn underflow_pos_errors() {
        let f = write("chr\tpos\tmet_reads\tnonmet_reads\trate\nchr1\t1\t1\t0\t1.0\n");
        let r = make_ref();
        assert!(read(f.path(), &r).is_err());
    }

    #[test]
    fn unknown_chrom_dropped() {
        let f = write(
            "chr\tpos\tmet_reads\tnonmet_reads\trate\n\
             chr1\t10\t1\t0\t1.0\n\
             chrX\t99\t1\t0\t1.0\n",
        );
        let r = make_ref();
        let calls = read(f.path(), &r).unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].chrom_id, 0);
    }

    #[test]
    fn output_sorted() {
        let f = write(
            "chr\tpos\tmet_reads\tnonmet_reads\trate\n\
             chr1\t30\t1\t0\t1.0\n\
             chr1\t10\t0\t1\t0.0\n\
             chr1\t20\t1\t1\t0.5\n",
        );
        let r = make_ref();
        let calls = read(f.path(), &r).unwrap();
        assert_eq!(calls[0].pos, 8);
        assert_eq!(calls[1].pos, 18);
        assert_eq!(calls[2].pos, 28);
    }
}

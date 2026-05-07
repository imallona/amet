use crate::error::{AmetError, Result};
use crate::io::open_read;
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::path::{Path, PathBuf};

/// Path of the cached CpG index for a FASTA: `<fasta>.cpg` next to the input.
pub fn cpg_index_path(fasta: &Path) -> PathBuf {
    let mut s = fasta.as_os_str().to_owned();
    s.push(".cpg");
    PathBuf::from(s)
}

/// Ensure a CpG index exists next to the FASTA. Builds it if absent. Returns the index path.
///
/// The index format is the same TSV (chrom\t0-based-pos) consumed by [`crate::reference::read_cpg_reference`],
/// so callers can pass the returned path through unchanged.
pub fn ensure_cpg_index(fasta: &Path) -> Result<PathBuf> {
    let index = cpg_index_path(fasta);
    if index.exists() {
        return Ok(index);
    }
    eprintln!(
        "[amet] CpG index not found for {}; scanning FASTA",
        fasta.display()
    );
    write_cpg_index_from_fasta(fasta, &index)?;
    eprintln!("[amet] wrote CpG index: {}", index.display());
    Ok(index)
}

/// Stream a FASTA and write the CpG index TSV.
fn write_cpg_index_from_fasta(fasta: &Path, index: &Path) -> Result<()> {
    let reader = open_read(fasta)?;
    let out = File::create(index).map_err(AmetError::Io)?;
    let mut writer = BufWriter::new(out);
    scan_fasta(reader, |chrom, pos| {
        writeln!(writer, "{}\t{}", chrom, pos).map_err(AmetError::Io)
    })?;
    writer.flush().map_err(AmetError::Io)?;
    Ok(())
}

/// Stream a FASTA and call `emit(chrom, pos)` for every CpG (0-based start of the C on +).
///
/// Single-pass, O(1) extra memory per chromosome. Tracks the last byte of the previous
/// sequence line so that CpGs spanning a line break are caught.
fn scan_fasta<R, F>(reader: R, mut emit: F) -> Result<()>
where
    R: BufRead,
    F: FnMut(&str, u64) -> Result<()>,
{
    let mut chrom = String::new();
    let mut prev_byte: Option<u8> = None;
    let mut pos: u64 = 0;

    for line in reader.lines() {
        let line = line.map_err(AmetError::Io)?;
        if let Some(stripped) = line.strip_prefix('>') {
            chrom = stripped
                .split_ascii_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            prev_byte = None;
            pos = 0;
            continue;
        }
        if chrom.is_empty() {
            continue;
        }
        let bytes = line.as_bytes();
        if bytes.is_empty() {
            continue;
        }
        if matches!(prev_byte, Some(b'C') | Some(b'c')) && (bytes[0] == b'G' || bytes[0] == b'g') {
            emit(&chrom, pos - 1)?;
        }
        for i in 0..bytes.len().saturating_sub(1) {
            let a = bytes[i];
            let b = bytes[i + 1];
            if (a == b'C' || a == b'c') && (b == b'G' || b == b'g') {
                emit(&chrom, pos + i as u64)?;
            }
        }
        prev_byte = bytes.last().copied();
        pos += bytes.len() as u64;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_fasta(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        write!(f, "{}", content).unwrap();
        f.flush().unwrap();
        f
    }

    fn collect(content: &str) -> Vec<(String, u64)> {
        let mut out = Vec::new();
        scan_fasta(content.as_bytes(), |c, p| {
            out.push((c.to_string(), p));
            Ok(())
        })
        .unwrap();
        out
    }

    #[test]
    fn single_chrom_one_cpg() {
        assert_eq!(collect(">chr1\nAAACGAAA\n"), vec![("chr1".into(), 3)]);
    }

    #[test]
    fn no_cpg() {
        assert_eq!(collect(">chr1\nAAAA\n"), Vec::<(String, u64)>::new());
    }

    #[test]
    fn case_insensitive() {
        assert_eq!(collect(">x\nacgT\n"), vec![("x".into(), 1)]);
        assert_eq!(collect(">x\nAcGt\n"), vec![("x".into(), 1)]);
    }

    #[test]
    fn cpg_spans_line_break() {
        // C at end of line 1, G at start of line 2; CG starts at position 3 of chrom.
        let v = collect(">chr1\nAAAC\nGAAA\n");
        assert_eq!(v, vec![("chr1".into(), 3)]);
    }

    #[test]
    fn multiple_cpgs_one_line() {
        // CG at 0 and 5
        assert_eq!(
            collect(">chr1\nCGAAACG\n"),
            vec![("chr1".into(), 0), ("chr1".into(), 5)]
        );
    }

    #[test]
    fn header_first_token_only() {
        assert_eq!(
            collect(">chr1 GRCm38 dna:chromosome\nACGT\n"),
            vec![("chr1".into(), 1)]
        );
    }

    #[test]
    fn multiple_chromosomes_reset_position() {
        let v = collect(">a\nACG\n>b\nACG\n");
        assert_eq!(v, vec![("a".into(), 1), ("b".into(), 1)]);
    }

    #[test]
    fn cross_chrom_no_cpg_emitted() {
        // Last char of chrom a is 'C', first char of chrom b is 'G'. Must not emit.
        let v = collect(">a\nAAAC\n>b\nGAAA\n");
        assert!(v.is_empty());
    }

    #[test]
    fn writes_cpg_index_file() {
        let fa = write_fasta(">chr1\nACGTACG\n");
        let idx = cpg_index_path(fa.path());
        write_cpg_index_from_fasta(fa.path(), &idx).unwrap();
        let body = std::fs::read_to_string(&idx).unwrap();
        assert_eq!(body, "chr1\t1\nchr1\t5\n");
        std::fs::remove_file(idx).ok();
    }

    #[test]
    fn ensure_skips_existing_index() {
        let fa = write_fasta(">chr1\nACGT\n");
        let idx = cpg_index_path(fa.path());
        // Pre-write a sentinel index that the FASTA scanner would NOT produce.
        std::fs::write(&idx, "sentinel\t999\n").unwrap();
        let returned = ensure_cpg_index(fa.path()).unwrap();
        assert_eq!(returned, idx);
        let body = std::fs::read_to_string(&idx).unwrap();
        assert!(body.starts_with("sentinel"));
        std::fs::remove_file(idx).ok();
    }
}

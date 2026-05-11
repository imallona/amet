use crate::error::{AmetError, Result};
use crate::io::open_read;
use fs2::FileExt;
use std::fs::OpenOptions;
use std::io::{BufRead, BufWriter, Write};
use std::path::{Path, PathBuf};

/// Path of the cached CpG index for a FASTA: `<fasta>.cpg` next to the input.
pub fn cpg_index_path(fasta: &Path) -> PathBuf {
    let mut s = fasta.as_os_str().to_owned();
    s.push(".cpg");
    PathBuf::from(s)
}

fn cpg_index_lock_path(fasta: &Path) -> PathBuf {
    let mut s = fasta.as_os_str().to_owned();
    s.push(".cpg.lock");
    PathBuf::from(s)
}

/// Ensure `<fasta>.cpg` exists and return its path. Safe under concurrent calls:
/// an exclusive flock on `<fasta>.cpg.lock` plus write-to-temp + atomic rename
/// prevents partial files and duplicate scans.
pub fn ensure_cpg_index(fasta: &Path) -> Result<PathBuf> {
    let index = cpg_index_path(fasta);
    if index.exists() {
        return Ok(index);
    }

    let lock_path = cpg_index_lock_path(fasta);
    let lock = OpenOptions::new()
        .create(true)
        .read(true)
        .write(true)
        .truncate(false)
        .open(&lock_path)
        .map_err(AmetError::Io)?;
    lock.lock_exclusive().map_err(AmetError::Io)?;

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

/// Stream a FASTA and write the CpG index TSV via sibling temp + atomic rename.
fn write_cpg_index_from_fasta(fasta: &Path, index: &Path) -> Result<()> {
    let parent = index.parent().unwrap_or(Path::new("."));
    let mut tmp = tempfile::Builder::new()
        .prefix(".cpg.tmp.")
        .tempfile_in(parent)
        .map_err(AmetError::Io)?;
    {
        let reader = open_read(fasta)?;
        let mut writer = BufWriter::new(tmp.as_file_mut());
        scan_fasta(reader, |chrom, pos| {
            writeln!(writer, "{}\t{}", chrom, pos).map_err(AmetError::Io)
        })?;
        writer.flush().map_err(AmetError::Io)?;
    }
    tmp.persist(index).map_err(|e| AmetError::Io(e.error))?;
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

    #[test]
    fn concurrent_ensure_yields_intact_index() {
        let dir = tempfile::tempdir().unwrap();
        let fa_path = dir.path().join("g.fa");
        let mut content = String::new();
        for c in 0..8 {
            content.push_str(&format!(">chr{}\n", c));
            for _ in 0..200 {
                content.push_str("ACGTACGTACGTACGTACGT\n");
            }
        }
        std::fs::write(&fa_path, &content).unwrap();

        let n_threads = 16;
        let handles: Vec<_> = (0..n_threads)
            .map(|_| {
                let p = fa_path.clone();
                std::thread::spawn(move || ensure_cpg_index(&p).unwrap())
            })
            .collect();
        for h in handles {
            let returned = h.join().unwrap();
            assert_eq!(returned, cpg_index_path(&fa_path));
        }

        let idx = cpg_index_path(&fa_path);
        let reference = crate::reference::read_cpg_reference(&idx).unwrap();
        let total_cpgs: usize = reference.positions.iter().map(|v| v.len()).sum();
        assert_eq!(total_cpgs, 8000);
    }
}

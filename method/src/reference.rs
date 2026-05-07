use crate::error::{AmetError, Result};
use crate::io::open_read;
use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

/// All CpG positions in the genome, indexed by chromosome.
///
/// `chrom_id_of` maps a chromosome name to a numeric ID; `chrom_names` is the inverse.
/// `positions[chrom_id]` is a sorted ascending vector of 0-based CpG start positions on that chromosome.
#[derive(Debug, Clone, Default)]
pub struct CpgReference {
    pub chrom_names: Vec<String>,
    pub chrom_id_of: HashMap<String, u32>,
    pub positions: Vec<Vec<u64>>,
}

impl CpgReference {
    pub fn chrom_id(&self, name: &str) -> Option<u32> {
        self.chrom_id_of.get(name).copied()
    }

    pub fn chrom_name(&self, id: u32) -> &str {
        &self.chrom_names[id as usize]
    }

    /// Number of CpGs on a chromosome.
    pub fn n_cpgs(&self, chrom_id: u32) -> usize {
        self.positions[chrom_id as usize].len()
    }
}

/// Read a CpG reference from a TSV file with two columns: chromosome name and 0-based position.
/// Lines may optionally start with a '#' comment marker (skipped). The file must be sorted by
/// (chrom, pos); we verify this and error on violation.
pub fn read_cpg_reference(path: &Path) -> Result<CpgReference> {
    let reader = open_read(path)?;
    let mut chrom_names: Vec<String> = Vec::new();
    let mut chrom_id_of: HashMap<String, u32> = HashMap::new();
    let mut positions: Vec<Vec<u64>> = Vec::new();
    let mut current_chrom: Option<u32> = None;
    let mut last_pos: u64 = 0;

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let mut fields = trimmed.split('\t');
        let chrom = fields.next().ok_or_else(|| AmetError::Parse {
            file: path.display().to_string(),
            line: i + 1,
            msg: "missing chromosome".to_string(),
        })?;
        let pos: u64 = fields
            .next()
            .ok_or_else(|| AmetError::Parse {
                file: path.display().to_string(),
                line: i + 1,
                msg: "missing position".to_string(),
            })?
            .parse()
            .map_err(|_| AmetError::Parse {
                file: path.display().to_string(),
                line: i + 1,
                msg: "position is not an integer".to_string(),
            })?;

        let chrom_id = match chrom_id_of.get(chrom) {
            Some(&id) => id,
            None => {
                let id = chrom_names.len() as u32;
                chrom_names.push(chrom.to_string());
                chrom_id_of.insert(chrom.to_string(), id);
                // Register the other prefix convention so lookups work whether
                // the caller uses `19` or `chr19`. The canonical name in
                // chrom_names (and therefore in output) is whatever the
                // reference used.
                if let Some(stripped) = chrom.strip_prefix("chr") {
                    if !stripped.is_empty() && !chrom_id_of.contains_key(stripped) {
                        chrom_id_of.insert(stripped.to_string(), id);
                    }
                } else {
                    let prefixed = format!("chr{}", chrom);
                    chrom_id_of.entry(prefixed).or_insert(id);
                }
                positions.push(Vec::new());
                id
            }
        };

        if Some(chrom_id) != current_chrom {
            current_chrom = Some(chrom_id);
            last_pos = 0;
        }
        if pos < last_pos {
            return Err(AmetError::Parse {
                file: path.display().to_string(),
                line: i + 1,
                msg: format!(
                    "positions must be sorted ascending within each chromosome (got {} after {})",
                    pos, last_pos
                ),
            });
        }
        last_pos = pos;
        positions[chrom_id as usize].push(pos);
    }

    Ok(CpgReference {
        chrom_names,
        chrom_id_of,
        positions,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        write!(f, "{}", content).unwrap();
        f.flush().unwrap();
        f
    }

    #[test]
    fn parse_simple() {
        let f = write("chr1\t100\nchr1\t200\nchr2\t50\n");
        let r = read_cpg_reference(f.path()).unwrap();
        assert_eq!(r.chrom_names, vec!["chr1", "chr2"]);
        assert_eq!(r.positions[0], vec![100, 200]);
        assert_eq!(r.positions[1], vec![50]);
    }

    #[test]
    fn skip_comments_and_blanks() {
        let f = write("# header\nchr1\t100\n\nchr1\t200\n");
        let r = read_cpg_reference(f.path()).unwrap();
        assert_eq!(r.positions[0], vec![100, 200]);
    }

    #[test]
    fn error_on_unsorted() {
        let f = write("chr1\t200\nchr1\t100\n");
        assert!(read_cpg_reference(f.path()).is_err());
    }

    #[test]
    fn lookup_chrom() {
        let f = write("chr1\t100\nchr2\t50\n");
        let r = read_cpg_reference(f.path()).unwrap();
        assert_eq!(r.chrom_id("chr1"), Some(0));
        assert_eq!(r.chrom_id("chr2"), Some(1));
        assert_eq!(r.chrom_id("chrX"), None);
        assert_eq!(r.chrom_name(0), "chr1");
    }

    #[test]
    fn lookup_handles_chr_prefix_either_way() {
        // Reference uses UCSC naming. Lookups must succeed for Ensembl-style names too.
        let f = write("chr1\t100\nchr19\t200\n");
        let r = read_cpg_reference(f.path()).unwrap();
        assert_eq!(r.chrom_id("chr1"), Some(0));
        assert_eq!(r.chrom_id("1"), Some(0));
        assert_eq!(r.chrom_id("chr19"), Some(1));
        assert_eq!(r.chrom_id("19"), Some(1));
        // Canonical names in output stay UCSC-style.
        assert_eq!(r.chrom_name(0), "chr1");
        assert_eq!(r.chrom_name(1), "chr19");
    }

    #[test]
    fn lookup_handles_ensembl_reference_too() {
        // Reference uses Ensembl naming. Lookups must succeed for UCSC-style names too.
        let f = write("1\t100\n19\t200\n");
        let r = read_cpg_reference(f.path()).unwrap();
        assert_eq!(r.chrom_id("1"), Some(0));
        assert_eq!(r.chrom_id("chr1"), Some(0));
        assert_eq!(r.chrom_id("19"), Some(1));
        assert_eq!(r.chrom_id("chr19"), Some(1));
        assert_eq!(r.chrom_name(0), "1");
    }
}

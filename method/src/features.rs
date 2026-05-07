use crate::error::{AmetError, Result};
use crate::io::open_read;
use crate::reference::CpgReference;
use std::io::BufRead;
use std::path::Path;

/// One BED feature, with its CpG range looked up against the reference.
///
/// `cpg_start_idx..cpg_end_idx` is the half-open index range into
/// `CpgReference.positions[chrom_id]` of CpGs falling inside this feature.
#[derive(Debug, Clone)]
pub struct Feature {
    pub feature_id: String,
    pub chrom_id: u32,
    pub start: u64,
    pub end: u64,
    pub cpg_start_idx: usize,
    pub cpg_end_idx: usize,
}

/// Read a BED file. The 4th column is used as `feature_id` if present, otherwise
/// `chrom:start-end`. Features whose chromosome is not in the reference are skipped
/// with a stderr warning.
pub fn read_features(path: &Path, reference: &CpgReference) -> Result<Vec<Feature>> {
    let reader = open_read(path)?;
    let mut features = Vec::new();
    let mut warned_chroms: std::collections::HashSet<String> = std::collections::HashSet::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') || trimmed.starts_with("track") {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 3 {
            return Err(AmetError::Parse {
                file: path.display().to_string(),
                line: i + 1,
                msg: "BED requires at least 3 columns".to_string(),
            });
        }
        let chrom = fields[0];
        let start: u64 = fields[1].parse().map_err(|_| AmetError::Parse {
            file: path.display().to_string(),
            line: i + 1,
            msg: "start is not an integer".to_string(),
        })?;
        let end: u64 = fields[2].parse().map_err(|_| AmetError::Parse {
            file: path.display().to_string(),
            line: i + 1,
            msg: "end is not an integer".to_string(),
        })?;
        let feature_id = if fields.len() > 3 {
            fields[3].to_string()
        } else {
            format!("{}:{}-{}", chrom, start, end)
        };

        let chrom_id = match reference.chrom_id(chrom) {
            Some(c) => c,
            None => {
                if warned_chroms.insert(chrom.to_string()) {
                    eprintln!(
                        "warning: chromosome {} in BED is not in the CpG reference; skipping its features",
                        chrom
                    );
                }
                continue;
            }
        };

        let positions = &reference.positions[chrom_id as usize];
        let cpg_start_idx = positions.partition_point(|&p| p < start);
        let cpg_end_idx = positions.partition_point(|&p| p < end);

        features.push(Feature {
            feature_id,
            chrom_id,
            start,
            end,
            cpg_start_idx,
            cpg_end_idx,
        });
    }

    Ok(features)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn make_ref() -> CpgReference {
        CpgReference {
            chrom_names: vec!["chr1".to_string()],
            chrom_id_of: [("chr1".to_string(), 0u32)].into_iter().collect(),
            positions: vec![vec![10, 20, 30, 50, 100, 150]],
        }
    }

    fn write(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        write!(f, "{}", content).unwrap();
        f.flush().unwrap();
        f
    }

    #[test]
    fn cpg_range_lookup() {
        let r = make_ref();
        let f = write("chr1\t15\t60\tregion1\n");
        let feats = read_features(f.path(), &r).unwrap();
        assert_eq!(feats.len(), 1);
        // CpGs at 20, 30, 50 fall inside [15, 60).
        assert_eq!(feats[0].cpg_start_idx, 1);
        assert_eq!(feats[0].cpg_end_idx, 4);
    }

    #[test]
    fn empty_feature() {
        let r = make_ref();
        let f = write("chr1\t11\t19\tempty\n");
        let feats = read_features(f.path(), &r).unwrap();
        assert_eq!(feats[0].cpg_start_idx, feats[0].cpg_end_idx);
    }

    #[test]
    fn unknown_chrom_skipped() {
        let r = make_ref();
        let f = write("chrX\t100\t200\tx1\nchr1\t10\t100\tknown\n");
        let feats = read_features(f.path(), &r).unwrap();
        assert_eq!(feats.len(), 1);
        assert_eq!(feats[0].feature_id, "known");
    }

    #[test]
    fn default_feature_id() {
        let r = make_ref();
        let f = write("chr1\t10\t100\n");
        let feats = read_features(f.path(), &r).unwrap();
        assert_eq!(feats[0].feature_id, "chr1:10-100");
    }
}

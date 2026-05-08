//! Per-cell methylation file parsers.
//!
//! All parsers normalise to a vector of MethCall records, one per CpG dinucleotide,
//! with pos = 0-based start of the C on the + strand. Per-strand records are merged
//! and counts summed.

pub mod allc;
pub mod bismark;
pub mod scnmt;

use crate::MethCall;
use crate::error::Result;
use crate::reference::CpgReference;
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CellFormat {
    Allc,
    Scnmt,
    Bismark,
}

impl CellFormat {
    pub fn parse(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "allc" | "methylpy" => Some(Self::Allc),
            "scnmt" | "cpg_level" => Some(Self::Scnmt),
            "bismark" | "singlec" => Some(Self::Bismark),
            _ => None,
        }
    }

    pub fn detect_from_path(path: &Path) -> Self {
        let name = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("")
            .to_lowercase();
        if name.contains(".cpg_level.") || name.contains(".scnmt.") {
            Self::Scnmt
        } else if name.contains(".singlec.") || name.contains(".bismark.") {
            Self::Bismark
        } else {
            // allc is the most general single-cell format; default when nothing matches.
            Self::Allc
        }
    }
}

pub fn read_cell(
    path: &Path,
    format: CellFormat,
    reference: &CpgReference,
) -> Result<Vec<MethCall>> {
    match format {
        CellFormat::Allc => allc::read(path, reference),
        CellFormat::Scnmt => scnmt::read(path, reference),
        CellFormat::Bismark => bismark::read(path, reference),
    }
}

#[cfg(test)]
mod detection_tests {
    use super::*;

    #[test]
    fn detect_allc() {
        assert_eq!(
            CellFormat::detect_from_path(Path::new("cell1.allc.tsv.gz")),
            CellFormat::Allc
        );
    }

    #[test]
    fn detect_scnmt() {
        assert_eq!(
            CellFormat::detect_from_path(Path::new("cell1.cpg_level.tsv.gz")),
            CellFormat::Scnmt
        );
    }

    #[test]
    fn detect_default_is_allc() {
        assert_eq!(
            CellFormat::detect_from_path(Path::new("cell1.tsv.gz")),
            CellFormat::Allc
        );
    }

    #[test]
    fn from_str_aliases() {
        assert_eq!(CellFormat::parse("ALLC"), Some(CellFormat::Allc));
        assert_eq!(CellFormat::parse("cpg_level"), Some(CellFormat::Scnmt));
        assert_eq!(CellFormat::parse("Bismark"), Some(CellFormat::Bismark));
        assert_eq!(CellFormat::parse("singleC"), Some(CellFormat::Bismark));
        assert_eq!(CellFormat::parse("nonsense"), None);
    }

    #[test]
    fn detect_bismark() {
        assert_eq!(
            CellFormat::detect_from_path(Path::new("GSM_xxx.singleC.txt.gz")),
            CellFormat::Bismark
        );
    }
}

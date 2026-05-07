use crate::error::{AmetError, Result};
use crate::io::open_read;
use std::collections::HashMap;
use std::io::BufRead;
use std::path::{Path, PathBuf};

/// One row of the manifest TSV. `extra` carries any non-required columns by name.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CellRow {
    pub cell_id: String,
    pub group: String,
    pub path: PathBuf,
    pub format: Option<String>,
    pub extra: HashMap<String, String>,
}

/// Parse a tab-separated manifest with required columns `cell_id`, `path`,
/// optional column `group` (default "all"), optional column `format`, and any
/// number of extra columns that are passed through verbatim.
///
/// `group_column` chooses which column to read as the group; default "group".
pub fn read_manifest(path: &Path, group_column: &str) -> Result<Vec<CellRow>> {
    let reader = open_read(path)?;
    let mut lines = reader.lines();
    let header = lines
        .next()
        .ok_or_else(|| AmetError::Manifest("manifest is empty".to_string()))??;
    let header_fields: Vec<&str> = header.split('\t').collect();

    let cell_id_col = header_fields
        .iter()
        .position(|&h| h == "cell_id")
        .ok_or_else(|| AmetError::Manifest("missing required column 'cell_id'".to_string()))?;
    let path_col = header_fields
        .iter()
        .position(|&h| h == "path")
        .ok_or_else(|| AmetError::Manifest("missing required column 'path'".to_string()))?;
    let group_col = header_fields.iter().position(|&h| h == group_column);
    let format_col = header_fields.iter().position(|&h| h == "format");

    let mut rows = Vec::new();
    let mut warned_missing_group = false;

    for (i, line) in lines.enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < header_fields.len() {
            return Err(AmetError::Parse {
                file: path.display().to_string(),
                line: i + 2,
                msg: format!(
                    "row has {} fields, header has {}",
                    fields.len(),
                    header_fields.len()
                ),
            });
        }
        let cell_id = fields[cell_id_col].to_string();
        let path = PathBuf::from(fields[path_col]);
        let group = match group_col {
            Some(c) => fields[c].to_string(),
            None => {
                if !warned_missing_group {
                    eprintln!(
                        "warning: no '{}' column in manifest; treating all cells as one group (all)",
                        group_column
                    );
                    warned_missing_group = true;
                }
                "all".to_string()
            }
        };
        let format = format_col.map(|c| fields[c].to_string());
        let mut extra = HashMap::new();
        for (idx, header_name) in header_fields.iter().enumerate() {
            if idx == cell_id_col || idx == path_col {
                continue;
            }
            if Some(idx) == group_col || Some(idx) == format_col {
                continue;
            }
            extra.insert((*header_name).to_string(), fields[idx].to_string());
        }
        rows.push(CellRow {
            cell_id,
            group,
            path,
            format,
            extra,
        });
    }

    if rows.is_empty() {
        return Err(AmetError::Manifest(
            "manifest contains no cells".to_string(),
        ));
    }

    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_manifest(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        write!(f, "{}", content).unwrap();
        f.flush().unwrap();
        f
    }

    #[test]
    fn parse_basic_manifest() {
        let f = write_manifest(
            "cell_id\tgroup\tpath\nA\texcitatory\t/tmp/a.allc.gz\nB\tinhibitory\t/tmp/b.allc.gz\n",
        );
        let rows = read_manifest(f.path(), "group").unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].cell_id, "A");
        assert_eq!(rows[0].group, "excitatory");
        assert_eq!(rows[1].cell_id, "B");
        assert_eq!(rows[1].group, "inhibitory");
    }

    #[test]
    fn parse_with_extra_columns() {
        let f = write_manifest(
            "cell_id\tgroup\tpath\tbatch\tdonor\nA\tex\t/tmp/a\tb1\td1\nB\tex\t/tmp/b\tb2\td1\n",
        );
        let rows = read_manifest(f.path(), "group").unwrap();
        assert_eq!(rows[0].extra.get("batch"), Some(&"b1".to_string()));
        assert_eq!(rows[0].extra.get("donor"), Some(&"d1".to_string()));
    }

    #[test]
    fn parse_with_format_override() {
        let f = write_manifest("cell_id\tgroup\tpath\tformat\nA\tex\t/tmp/a\tallc\n");
        let rows = read_manifest(f.path(), "group").unwrap();
        assert_eq!(rows[0].format, Some("allc".to_string()));
    }

    #[test]
    fn alternative_group_column() {
        let f = write_manifest("cell_id\tcluster\tpath\nA\tneuron\t/tmp/a\n");
        let rows = read_manifest(f.path(), "cluster").unwrap();
        assert_eq!(rows[0].group, "neuron");
    }

    #[test]
    fn missing_group_column_defaults_to_all() {
        let f = write_manifest("cell_id\tpath\nA\t/tmp/a\nB\t/tmp/b\n");
        let rows = read_manifest(f.path(), "group").unwrap();
        assert_eq!(rows[0].group, "all");
        assert_eq!(rows[1].group, "all");
    }

    #[test]
    fn empty_manifest_errors() {
        let f = write_manifest("");
        assert!(read_manifest(f.path(), "group").is_err());
    }

    #[test]
    fn missing_required_column_errors() {
        let f = write_manifest("cell_id\tgroup\nA\tex\n");
        assert!(read_manifest(f.path(), "group").is_err());
    }
}

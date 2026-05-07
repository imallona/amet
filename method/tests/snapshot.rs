//! End-to-end regression test against checked-in golden outputs.
//!
//! Runs amet against the fixture in tests/snapshot/data/ and asserts the
//! resulting cell_feature.tsv and feature.tsv match tests/snapshot/golden/
//! byte-for-byte. The fixture is small (4 cells, 1 feature, 6 CpGs) and
//! purpose-built so that:
//!
//!   - cells A and B (group g1) have mirror-image patterns 000111 vs 111000:
//!     identical i_total, but jsd > 0 within the group.
//!   - cells C and D (group g2) have identical 010101 patterns:
//!     identical i_total, jsd == 0 within the group.
//!
//! If amet's scoring math, parsers, or output format change, this test
//! catches the drift. To regenerate the goldens after an intentional
//! change, run with `UPDATE_SNAPSHOTS=1 cargo test snapshot`.
//!
//! The fixture and goldens live next to this file so they ship with the
//! repo and run identically in CI without external data.

use flate2::read::MultiGzDecoder;
use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::tempdir;

fn binary_path() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_amet"))
}

fn fixture_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/snapshot")
}

fn read_gz(path: &Path) -> String {
    let f = fs::File::open(path).unwrap();
    let mut s = String::new();
    MultiGzDecoder::new(f).read_to_string(&mut s).unwrap();
    s
}

#[test]
fn snapshot_matches_golden() {
    let fix = fixture_root();
    let work = tempdir().unwrap();

    // Copy the FASTA into the workdir so the .cpg sidecar is created in
    // the tempdir, not in the checked-in fixture.
    let fa_src = fix.join("data/tiny.fa");
    let fa = work.path().join("tiny.fa");
    fs::copy(&fa_src, &fa).unwrap();

    let bed = fix.join("data/features.bed");
    let cells_dir = fix.join("data");

    // Manifest with absolute paths to the fixture cell files.
    let manifest_path = work.path().join("cells.tsv");
    let mut manifest = fs::File::create(&manifest_path).unwrap();
    writeln!(manifest, "cell_id\tgroup\tpath").unwrap();
    for (id, group, file) in [
        ("A", "g1", "cellA.allc.tsv"),
        ("B", "g1", "cellB.allc.tsv"),
        ("C", "g2", "cellC.allc.tsv"),
        ("D", "g2", "cellD.allc.tsv"),
    ] {
        writeln!(
            manifest,
            "{}\t{}\t{}",
            id,
            group,
            cells_dir.join(file).display()
        )
        .unwrap();
    }
    drop(manifest);

    let prefix = work.path().join("run");
    let status = Command::new(binary_path())
        .args([
            "--genome",
            fa.to_str().unwrap(),
            "--features",
            bed.to_str().unwrap(),
            "--cells",
            manifest_path.to_str().unwrap(),
            "--output-prefix",
            prefix.to_str().unwrap(),
            "--min-cpgs-per-feature",
            "3",
            "--min-cells-per-group",
            "2",
            "--i-max-lag",
            "3",
            "--threads",
            "1",
        ])
        .status()
        .expect("running amet");
    assert!(status.success(), "amet exited with non-zero status");

    let cf_actual = read_gz(&work.path().join("run.cell_feature.tsv.gz"));
    let feat_actual = read_gz(&work.path().join("run.feature.tsv.gz"));

    let cf_golden_path = fix.join("golden/cell_feature.tsv");
    let feat_golden_path = fix.join("golden/feature.tsv");

    if std::env::var("UPDATE_SNAPSHOTS").is_ok() {
        fs::write(&cf_golden_path, &cf_actual).unwrap();
        fs::write(&feat_golden_path, &feat_actual).unwrap();
        eprintln!("[snapshot] wrote new goldens to {}", fix.display());
        return;
    }

    let cf_golden = fs::read_to_string(&cf_golden_path).unwrap();
    let feat_golden = fs::read_to_string(&feat_golden_path).unwrap();

    assert_eq!(
        cf_actual,
        cf_golden,
        "cell_feature output drifted vs golden at {}.\n\
         To accept the new output, rerun with UPDATE_SNAPSHOTS=1.",
        cf_golden_path.display()
    );
    assert_eq!(
        feat_actual,
        feat_golden,
        "feature output drifted vs golden at {}.\n\
         To accept the new output, rerun with UPDATE_SNAPSHOTS=1.",
        feat_golden_path.display()
    );
}

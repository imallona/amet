//! End-to-end smoke tests that exercise the CLI binary on tiny synthetic inputs.
//!
//! Running `cargo test --test integration` builds the binary if needed.

use std::io::Write;
use std::path::PathBuf;
use std::process::Command;
use tempfile::tempdir;

fn binary_path() -> PathBuf {
    // CARGO_BIN_EXE_amet is set by Cargo when running integration tests.
    PathBuf::from(env!("CARGO_BIN_EXE_amet"))
}

fn write_file(dir: &std::path::Path, name: &str, content: &str) -> PathBuf {
    let p = dir.join(name);
    let mut f = std::fs::File::create(&p).unwrap();
    f.write_all(content.as_bytes()).unwrap();
    p
}

#[test]
fn end_to_end_two_cells_one_feature_allc() {
    let dir = tempdir().unwrap();

    // CpG reference: 6 CpGs at 0-based positions 9, 19, 29, 39, 49, 59 on chr1.
    let cpgs = write_file(
        dir.path(),
        "cpgs.tsv",
        "chr1\t9\nchr1\t19\nchr1\t29\nchr1\t39\nchr1\t49\nchr1\t59\n",
    );

    // Feature covers all 6.
    let bed = write_file(dir.path(), "feat.bed", "chr1\t0\t100\tregion1\n");

    // Cell A: alternating 010101.
    let cell_a = write_file(
        dir.path(),
        "cellA.allc.tsv",
        "chr1\t10\t+\tCGN\t0\t1\t0\n\
         chr1\t20\t+\tCGN\t1\t1\t1\n\
         chr1\t30\t+\tCGN\t0\t1\t0\n\
         chr1\t40\t+\tCGN\t1\t1\t1\n\
         chr1\t50\t+\tCGN\t0\t1\t0\n\
         chr1\t60\t+\tCGN\t1\t1\t1\n",
    );

    // Cell B: shuffled 011010.
    let cell_b = write_file(
        dir.path(),
        "cellB.allc.tsv",
        "chr1\t10\t+\tCGN\t0\t1\t0\n\
         chr1\t20\t+\tCGN\t1\t1\t1\n\
         chr1\t30\t+\tCGN\t1\t1\t1\n\
         chr1\t40\t+\tCGN\t0\t1\t0\n\
         chr1\t50\t+\tCGN\t1\t1\t1\n\
         chr1\t60\t+\tCGN\t0\t1\t0\n",
    );

    let manifest = write_file(
        dir.path(),
        "cells.tsv",
        &format!(
            "cell_id\tgroup\tpath\nA\tg1\t{}\nB\tg1\t{}\n",
            cell_a.display(),
            cell_b.display()
        ),
    );

    let prefix = dir.path().join("run1");
    let status = Command::new(binary_path())
        .args([
            "--cpg-reference",
            cpgs.to_str().unwrap(),
            "--features",
            bed.to_str().unwrap(),
            "--cells",
            manifest.to_str().unwrap(),
            "--output-prefix",
            prefix.to_str().unwrap(),
            "--min-cpgs-per-feature",
            "3",
        ])
        .status()
        .expect("running amet binary");
    assert!(status.success(), "amet exited with non-zero status");

    let cf_path = dir.path().join("run1.cell_feature.tsv.gz");
    let feat_path = dir.path().join("run1.feature.tsv.gz");
    assert!(cf_path.exists(), "cell_feature output missing");
    assert!(feat_path.exists(), "feature output missing");

    // Decompress and inspect the cell_feature output.
    let cf_content = read_gz(&cf_path);
    let lines: Vec<&str> = cf_content.lines().collect();
    assert_eq!(lines.len(), 3, "header + 2 cells, got {}", lines.len());
    let header = lines[0];
    assert!(header.contains("i_total"));
    assert!(header.contains("i_1"));

    // The 010101 cell should have higher i_total than the 011010 cell.
    let row_a: Vec<&str> = lines[1].split('\t').collect();
    let row_b: Vec<&str> = lines[2].split('\t').collect();
    let i_total_col = header.split('\t').position(|h| h == "i_total").unwrap();
    let a_score: f64 = row_a[i_total_col].parse().unwrap();
    let b_score: f64 = row_b[i_total_col].parse().unwrap();
    assert!(
        a_score > b_score,
        "010101 (A) should have higher i_total than 011010 (B); got {} vs {}",
        a_score,
        b_score
    );

    // Feature output should have one row for group g1.
    let feat_content = read_gz(&feat_path);
    let feat_lines: Vec<&str> = feat_content.lines().collect();
    assert_eq!(feat_lines.len(), 2);
    assert!(feat_lines[1].contains("region1"));
    assert!(feat_lines[1].contains("g1"));
}

#[test]
fn iid_cells_have_low_i_total_at_any_p() {
    // Synthesize cells whose calls are independent draws at marginal p ≈ 0.5
    // and verify i_total is small. We approximate by a deterministic balanced sequence
    // with no comethylation: 000111000111 (low lag-1 MI when shuffled, but 010011000110
    // is closer to iid).
    let dir = tempdir().unwrap();

    let cpgs = (0..30)
        .map(|i| format!("chr1\t{}\n", i * 10 + 9))
        .collect::<String>();
    let cpg_path = write_file(dir.path(), "cpgs.tsv", &cpgs);
    let bed = write_file(dir.path(), "feat.bed", "chr1\t0\t1000\tr\n");

    // Pseudo-iid sequence: 30 bits with p = 0.5, no obvious adjacent correlation.
    let bits = [
        1u8, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1,
    ];
    let mut content = String::new();
    for (i, b) in bits.iter().enumerate() {
        let pos = (i * 10 + 10) as u64; // 1-based + strand C position
        content.push_str(&format!("chr1\t{}\t+\tCGN\t{}\t1\t{}\n", pos, b, b));
    }
    let cell = write_file(dir.path(), "cell.allc.tsv", &content);

    let manifest = write_file(
        dir.path(),
        "cells.tsv",
        &format!("cell_id\tgroup\tpath\nA\tg\t{}\n", cell.display()),
    );
    let prefix = dir.path().join("run");
    let status = Command::new(binary_path())
        .args([
            "--cpg-reference",
            cpg_path.to_str().unwrap(),
            "--features",
            bed.to_str().unwrap(),
            "--cells",
            manifest.to_str().unwrap(),
            "--output-prefix",
            prefix.to_str().unwrap(),
            "--min-cpgs-per-feature",
            "3",
            "--i-max-lag",
            "3",
        ])
        .status()
        .unwrap();
    assert!(status.success());

    let cf = read_gz(&dir.path().join("run.cell_feature.tsv.gz"));
    let lines: Vec<&str> = cf.lines().collect();
    let header = lines[0];
    let row: Vec<&str> = lines[1].split('\t').collect();
    let i_total_col = header.split('\t').position(|h| h == "i_total").unwrap();
    let i_total: f64 = row[i_total_col].parse().unwrap();

    // Sanity check rather than strict bound — small finite sample, MM correction adds noise.
    assert!(
        i_total < 1.5,
        "expected pseudo-iid sequence to have low I_total, got {}",
        i_total
    );
}

fn read_gz(path: &std::path::Path) -> String {
    use flate2::read::MultiGzDecoder;
    use std::io::Read;
    let f = std::fs::File::open(path).unwrap();
    let mut s = String::new();
    MultiGzDecoder::new(f).read_to_string(&mut s).unwrap();
    s
}

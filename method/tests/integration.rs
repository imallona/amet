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
    // Rows are looked up by cell_id (column 0): cell_feature row order is
    // cell-interleaved and not guaranteed to follow the manifest.
    let i_total_col = header.split('\t').position(|h| h == "i_total").unwrap();
    let score_of = |cell_id: &str| -> f64 {
        let row = lines[1..]
            .iter()
            .find(|l| l.split('\t').next() == Some(cell_id))
            .unwrap_or_else(|| panic!("no cell_feature row for cell {}", cell_id));
        row.split('\t').nth(i_total_col).unwrap().parse().unwrap()
    };
    let a_score = score_of("A");
    let b_score = score_of("B");
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
    // Synthesize cells whose calls are independent draws at marginal p ~ 0.5
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

    // Sanity check rather than strict bound: small finite sample, MM correction adds noise.
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

#[test]
fn end_to_end_with_genome_flag() {
    // FASTA -> amet derives the CpG index, writes <fasta>.cpg, then runs.
    // Six CpGs on chr1 at positions 9, 19, 29, 39, 49, 59 (0-based).
    let dir = tempdir().unwrap();
    let mut fa = String::from(">chr1\n");
    for i in 0..60 {
        if i == 9 || i == 19 || i == 29 || i == 39 || i == 49 || i == 59 {
            fa.push('C');
        } else if i == 10 || i == 20 || i == 30 || i == 40 || i == 50 || i == 60 {
            fa.push('G');
        } else {
            fa.push('A');
        }
    }
    fa.push('\n');
    let fa_path = write_file(dir.path(), "tiny.fa", &fa);
    let bed = write_file(dir.path(), "feat.bed", "chr1\t0\t100\tregion1\n");

    let cell = write_file(
        dir.path(),
        "cellA.allc.tsv",
        "chr1\t10\t+\tCGN\t0\t1\t0\n\
         chr1\t20\t+\tCGN\t1\t1\t1\n\
         chr1\t30\t+\tCGN\t0\t1\t0\n\
         chr1\t40\t+\tCGN\t1\t1\t1\n\
         chr1\t50\t+\tCGN\t0\t1\t0\n\
         chr1\t60\t+\tCGN\t1\t1\t1\n",
    );
    let manifest = write_file(
        dir.path(),
        "cells.tsv",
        &format!("cell_id\tgroup\tpath\nA\tg1\t{}\n", cell.display()),
    );
    let prefix = dir.path().join("run1");
    let status = Command::new(binary_path())
        .args([
            "--genome",
            fa_path.to_str().unwrap(),
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

    // The CpG index sidecar must exist next to the FASTA.
    let mut cpg_index = fa_path.clone();
    let mut s = cpg_index.into_os_string();
    s.push(".cpg");
    cpg_index = PathBuf::from(s);
    assert!(
        cpg_index.exists(),
        "expected sidecar at {}",
        cpg_index.display()
    );

    let cf = read_gz(&dir.path().join("run1.cell_feature.tsv.gz"));
    let lines: Vec<&str> = cf.lines().collect();
    assert_eq!(lines.len(), 2, "header + 1 cell, got {}", lines.len());
}

#[test]
fn cli_rejects_both_genome_and_cpg_reference() {
    let dir = tempdir().unwrap();
    let bed = write_file(dir.path(), "f.bed", "chr1\t0\t10\tx\n");
    let cells = write_file(dir.path(), "c.tsv", "cell_id\tgroup\tpath\n");
    let cpgs = write_file(dir.path(), "cpgs.tsv", "chr1\t0\n");
    let fa = write_file(dir.path(), "g.fa", ">chr1\nA\n");
    let prefix = dir.path().join("p");
    let status = Command::new(binary_path())
        .args([
            "--genome",
            fa.to_str().unwrap(),
            "--cpg-reference",
            cpgs.to_str().unwrap(),
            "--features",
            bed.to_str().unwrap(),
            "--cells",
            cells.to_str().unwrap(),
            "--output-prefix",
            prefix.to_str().unwrap(),
        ])
        .status()
        .expect("running amet binary");
    assert!(
        !status.success(),
        "expected non-zero exit when both flags given"
    );
}

#[test]
fn mixed_chrom_naming_still_runs() {
    // UCSC-style FASTA (chr1) but Ensembl-style cell + BED (1). Should resolve transparently.
    let dir = tempdir().unwrap();
    let mut fa = String::from(">chr1\n");
    for i in 0..60 {
        let c = if [9u64, 19, 29, 39, 49, 59].contains(&(i as u64)) {
            'C'
        } else if [10u64, 20, 30, 40, 50, 60].contains(&(i as u64)) {
            'G'
        } else {
            'A'
        };
        fa.push(c);
    }
    fa.push('\n');
    let fa_path = write_file(dir.path(), "tiny.fa", &fa);
    // BED uses Ensembl naming.
    let bed = write_file(dir.path(), "feat.bed", "1\t0\t100\tregion1\n");
    // Cell allc uses Ensembl naming.
    let cell = write_file(
        dir.path(),
        "cellA.allc.tsv",
        "1\t10\t+\tCGN\t0\t1\t0\n\
         1\t20\t+\tCGN\t1\t1\t1\n\
         1\t30\t+\tCGN\t0\t1\t0\n\
         1\t40\t+\tCGN\t1\t1\t1\n\
         1\t50\t+\tCGN\t0\t1\t0\n\
         1\t60\t+\tCGN\t1\t1\t1\n",
    );
    let manifest = write_file(
        dir.path(),
        "cells.tsv",
        &format!("cell_id\tgroup\tpath\nA\tg1\t{}\n", cell.display()),
    );
    let prefix = dir.path().join("run1");
    let status = Command::new(binary_path())
        .args([
            "--genome",
            fa_path.to_str().unwrap(),
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
    assert!(
        status.success(),
        "amet should resolve mixed-naming inputs without complaint"
    );
    let cf = read_gz(&dir.path().join("run1.cell_feature.tsv.gz"));
    let lines: Vec<&str> = cf.lines().collect();
    assert_eq!(lines.len(), 2, "header + 1 cell, got {}", lines.len());
}

#[test]
fn multi_features_writes_one_triplet_per_set() {
    // Two BEDs over the same cells, scored in one cell-read pass.
    // Each BED gets its own output triplet keyed by basename.
    let dir = tempdir().unwrap();

    let cpgs = write_file(
        dir.path(),
        "cpgs.tsv",
        "chr1\t9\nchr1\t19\nchr1\t29\nchr1\t39\nchr1\t49\nchr1\t59\n",
    );

    // BED 1: one region named "prom_a" over CpGs 1..3.
    let bed_prom = write_file(dir.path(), "prom.bed", "chr1\t0\t35\tprom_a\n");
    // BED 2: one region named "enh_a" over CpGs 4..6.
    let bed_enh = write_file(dir.path(), "enh.bed", "chr1\t35\t100\tenh_a\n");

    let cell = write_file(
        dir.path(),
        "cell.allc.tsv",
        "chr1\t10\t+\tCGN\t0\t1\t0\n\
         chr1\t20\t+\tCGN\t1\t1\t1\n\
         chr1\t30\t+\tCGN\t0\t1\t0\n\
         chr1\t40\t+\tCGN\t1\t1\t1\n\
         chr1\t50\t+\tCGN\t0\t1\t0\n\
         chr1\t60\t+\tCGN\t1\t1\t1\n",
    );
    let manifest = write_file(
        dir.path(),
        "cells.tsv",
        &format!("cell_id\tgroup\tpath\nA\tg1\t{}\n", cell.display()),
    );

    let prefix = dir.path().join("run");
    let status = Command::new(binary_path())
        .args([
            "--cpg-reference",
            cpgs.to_str().unwrap(),
            "--features",
            bed_prom.to_str().unwrap(),
            "--features",
            bed_enh.to_str().unwrap(),
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

    // Each set must have its own triplet, suffixed with the BED stem.
    let prom_cf = dir.path().join("run.prom.cell_feature.tsv.gz");
    let prom_feat = dir.path().join("run.prom.feature.tsv.gz");
    let prom_pair = dir.path().join("run.prom.pair_counts.tsv.gz");
    let enh_cf = dir.path().join("run.enh.cell_feature.tsv.gz");
    let enh_feat = dir.path().join("run.enh.feature.tsv.gz");
    let enh_pair = dir.path().join("run.enh.pair_counts.tsv.gz");
    for p in [
        &prom_cf, &prom_feat, &prom_pair, &enh_cf, &enh_feat, &enh_pair,
    ] {
        assert!(p.exists(), "expected output {} to exist", p.display());
    }

    // Single-features output paths must NOT exist (no fallback bleed-through).
    assert!(!dir.path().join("run.cell_feature.tsv.gz").exists());
    assert!(!dir.path().join("run.feature.tsv.gz").exists());
    assert!(!dir.path().join("run.pair_counts.tsv.gz").exists());

    // Each set sees only its own feature id.
    let prom_cf_text = read_gz(&prom_cf);
    let enh_cf_text = read_gz(&enh_cf);
    assert!(prom_cf_text.contains("prom_a"));
    assert!(!prom_cf_text.contains("enh_a"));
    assert!(enh_cf_text.contains("enh_a"));
    assert!(!enh_cf_text.contains("prom_a"));
}

#[test]
fn cli_rejects_zero_i_max_lag() {
    // --i-max-lag 0 leaves pair_tables empty, which would later panic on JSD aggregation.
    // The CLI is the right place to reject this; assert the binary exits non-zero.
    let dir = tempdir().unwrap();
    let cpgs = write_file(dir.path(), "cpgs.tsv", "chr1\t9\n");
    let bed = write_file(dir.path(), "f.bed", "chr1\t0\t100\tx\n");
    let cells = write_file(dir.path(), "c.tsv", "cell_id\tgroup\tpath\n");
    let prefix = dir.path().join("p");
    let status = Command::new(binary_path())
        .args([
            "--cpg-reference",
            cpgs.to_str().unwrap(),
            "--features",
            bed.to_str().unwrap(),
            "--cells",
            cells.to_str().unwrap(),
            "--output-prefix",
            prefix.to_str().unwrap(),
            "--i-max-lag",
            "0",
        ])
        .status()
        .expect("running amet binary");
    assert!(
        !status.success(),
        "expected non-zero exit when --i-max-lag is 0"
    );
}

#[test]
fn multi_features_rejects_duplicate_basenames() {
    // Two BEDs whose basenames collapse to the same label must be rejected.
    let dir = tempdir().unwrap();
    let cpgs = write_file(dir.path(), "cpgs.tsv", "chr1\t9\n");
    let sub = dir.path().join("sub");
    std::fs::create_dir_all(&sub).unwrap();
    let bed1 = write_file(dir.path(), "regions.bed", "chr1\t0\t100\tx\n");
    let bed2 = write_file(&sub, "regions.bed", "chr1\t0\t100\ty\n");
    let cells = write_file(dir.path(), "cells.tsv", "cell_id\tgroup\tpath\n");
    let prefix = dir.path().join("run");
    let status = Command::new(binary_path())
        .args([
            "--cpg-reference",
            cpgs.to_str().unwrap(),
            "--features",
            bed1.to_str().unwrap(),
            "--features",
            bed2.to_str().unwrap(),
            "--cells",
            cells.to_str().unwrap(),
            "--output-prefix",
            prefix.to_str().unwrap(),
        ])
        .status()
        .expect("running amet binary");
    assert!(
        !status.success(),
        "expected duplicate-label rejection to exit non-zero"
    );
}

#[test]
fn feature_tsv_multi_feature_multi_group() {
    // Two features in one BED, four cells across two groups. Exercises the
    // per-(feature index, group id) aggregate keying: each feature must get
    // its own feature.tsv row per group, correctly labelled, with no
    // cross-contamination between features.
    let dir = tempdir().unwrap();

    // 12 CpGs at 0-based positions 9, 19, ..., 119 on chr1.
    let cpgs: String = (0..12).map(|i| format!("chr1\t{}\n", i * 10 + 9)).collect();
    let cpg_path = write_file(dir.path(), "cpgs.tsv", &cpgs);

    // featA covers the first 6 CpGs (0-based start < 60), featB the last 6.
    let bed = write_file(
        dir.path(),
        "feat.bed",
        "chr1\t0\t60\tfeatA\nchr1\t60\t130\tfeatB\n",
    );

    // Every cell: featA CpGs unmethylated (m=0), featB CpGs methylated (m=1).
    // The allc + strand C position is the 0-based CpG start + 1, so 10..=120.
    let mut cell_body = String::new();
    for pos in (10..=60).step_by(10) {
        cell_body.push_str(&format!("chr1\t{}\t+\tCGN\t0\t1\t0\n", pos));
    }
    for pos in (70..=120).step_by(10) {
        cell_body.push_str(&format!("chr1\t{}\t+\tCGN\t1\t1\t1\n", pos));
    }
    let cell_paths: Vec<_> = ["A", "B", "C", "D"]
        .iter()
        .map(|id| write_file(dir.path(), &format!("cell{}.allc.tsv", id), &cell_body))
        .collect();

    let manifest = write_file(
        dir.path(),
        "cells.tsv",
        &format!(
            "cell_id\tgroup\tpath\nA\tg1\t{}\nB\tg1\t{}\nC\tg2\t{}\nD\tg2\t{}\n",
            cell_paths[0].display(),
            cell_paths[1].display(),
            cell_paths[2].display(),
            cell_paths[3].display(),
        ),
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
            "--min-cells-per-group",
            "2",
        ])
        .status()
        .expect("running amet binary");
    assert!(status.success(), "amet exited with non-zero status");

    // cell_feature: header + 4 cells x 2 features = 8 data rows.
    let cf = read_gz(&dir.path().join("run.cell_feature.tsv.gz"));
    assert_eq!(
        cf.lines().count(),
        9,
        "expected header + 8 cell_feature rows"
    );

    // feature.tsv: header + (2 features x 2 groups) = 4 data rows, sorted.
    let feat = read_gz(&dir.path().join("run.feature.tsv.gz"));
    let lines: Vec<&str> = feat.lines().collect();
    assert_eq!(
        lines.len(),
        5,
        "expected header + 4 feature rows, got {}",
        lines.len()
    );
    let header: Vec<&str> = lines[0].split('\t').collect();
    let col = |name: &str| header.iter().position(|h| *h == name).unwrap();
    let (fid, grp, ncells, meth) = (
        col("feature_id"),
        col("group"),
        col("n_cells"),
        col("mean_meth_mean"),
    );
    let rows: Vec<Vec<&str>> = lines[1..].iter().map(|l| l.split('\t').collect()).collect();

    // Sorted by (feature_id, group): featA/g1, featA/g2, featB/g1, featB/g2.
    let expect = [
        ("featA", "g1"),
        ("featA", "g2"),
        ("featB", "g1"),
        ("featB", "g2"),
    ];
    for (i, (efid, egrp)) in expect.iter().enumerate() {
        assert_eq!(rows[i][fid], *efid, "row {} feature_id", i);
        assert_eq!(rows[i][grp], *egrp, "row {} group", i);
        assert_eq!(rows[i][ncells], "2", "row {} n_cells", i);
    }

    // featA is all-unmethylated, featB all-methylated. If the feature-index
    // keying were wrong, these per-feature means would be swapped or merged.
    for row in &rows[0..2] {
        let m: f64 = row[meth].parse().unwrap();
        assert!(m < 0.5, "featA mean_meth_mean should be ~0, got {}", m);
    }
    for row in &rows[2..4] {
        let m: f64 = row[meth].parse().unwrap();
        assert!(m > 0.5, "featB mean_meth_mean should be ~1, got {}", m);
    }
}

#[test]
fn cli_rejects_neither_genome_nor_cpg_reference() {
    let dir = tempdir().unwrap();
    let bed = write_file(dir.path(), "f.bed", "chr1\t0\t10\tx\n");
    let cells = write_file(dir.path(), "c.tsv", "cell_id\tgroup\tpath\n");
    let prefix = dir.path().join("p");
    let status = Command::new(binary_path())
        .args([
            "--features",
            bed.to_str().unwrap(),
            "--cells",
            cells.to_str().unwrap(),
            "--output-prefix",
            prefix.to_str().unwrap(),
        ])
        .status()
        .expect("running amet binary");
    assert!(
        !status.success(),
        "expected non-zero exit when neither flag given"
    );
}

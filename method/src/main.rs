use amet::cli::Cli;
use amet::features::{Feature, read_features};
use amet::genome::ensure_cpg_index;
use amet::io::open_write;
use amet::kmer::{PairCounts, build_window, marginal_counts, pair_counts_all_lags};
use amet::manifest::{CellRow, read_manifest};
use amet::parsers::{CellFormat, read_cell};
use amet::reference::read_cpg_reference;
use amet::scores::{i_total::i_total, jsd::JsdAccumulator};
use anyhow::{Context, Result, anyhow};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

fn main() -> Result<()> {
    let cli = Cli::parse_args();

    if cli.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .context("failed to set up thread pool")?;
    }

    let cpg_path: PathBuf = match (&cli.genome, &cli.cpg_reference) {
        (Some(fa), None) => ensure_cpg_index(fa).context("building CpG index from genome")?,
        (None, Some(p)) => p.clone(),
        _ => unreachable!("clap ArgGroup ensures exactly one of --genome and --cpg-reference"),
    };
    if cli.build_cpg_only {
        return Ok(());
    }
    eprintln!("[amet] reading CpG reference: {}", cpg_path.display());
    let reference = read_cpg_reference(&cpg_path).context("reading CpG reference")?;
    let total_cpgs: usize = reference.positions.iter().map(|v| v.len()).sum();
    eprintln!(
        "[amet] reference: {} chromosomes, {} CpGs",
        reference.chrom_names.len(),
        total_cpgs
    );

    let cells_path = cli.cells.as_ref().expect("--cells required");
    let output_prefix = cli
        .output_prefix
        .as_ref()
        .expect("--output-prefix required");

    if cli.features.is_empty() {
        return Err(anyhow!("--features is required"));
    }

    // Resolve per-set labels and check uniqueness up front so we fail before
    // reading feature BEDs and creating per-set outputs.
    let labels: Vec<String> = cli.features.iter().map(|p| features_label(p)).collect();
    {
        let mut seen: HashSet<&str> = HashSet::new();
        for l in &labels {
            if !seen.insert(l.as_str()) {
                return Err(anyhow!(
                    "two --features BEDs resolve to the same label `{}`; rename one of the input files",
                    l
                ));
            }
        }
    }
    let single_set = cli.features.len() == 1;

    let mut sets: Vec<FeatureSet> = Vec::with_capacity(cli.features.len());
    let mut sinks: Vec<Mutex<SetSink>> = Vec::with_capacity(cli.features.len());
    let mut feat_writers: Vec<Box<dyn Write + Send>> = Vec::with_capacity(cli.features.len());
    for (path, label) in cli.features.iter().zip(labels.iter()) {
        eprintln!("[amet] reading features: {}", path.display());
        let features = read_features(path, &reference).context("reading features")?;
        eprintln!("[amet] features ({}): {}", label, features.len());
        let set_prefix = if single_set {
            output_prefix.clone()
        } else {
            with_suffix(output_prefix, &format!(".{}", label))
        };
        let cf_path = with_suffix(&set_prefix, ".cell_feature.tsv.gz");
        let feat_path = with_suffix(&set_prefix, ".feature.tsv.gz");
        let pair_path = with_suffix(&set_prefix, ".pair_counts.tsv.gz");
        let mut cf_writer = open_write(&cf_path).context("opening cell_feature output")?;
        let mut feat_writer = open_write(&feat_path).context("opening feature output")?;
        let mut pair_writer = open_write(&pair_path).context("opening pair_counts output")?;
        write_headers(
            &mut cf_writer,
            &mut feat_writer,
            &mut pair_writer,
            cli.i_max_lag,
        )?;
        sets.push(FeatureSet {
            label: label.clone(),
            features,
            cf_path,
            feat_path,
            pair_path,
        });
        sinks.push(Mutex::new(SetSink {
            cf_writer,
            pair_writer,
            agg: HashMap::new(),
        }));
        feat_writers.push(feat_writer);
    }

    eprintln!("[amet] reading manifest: {}", cells_path.display());
    let manifest = read_manifest(cells_path, &cli.group_column).context("reading manifest")?;
    eprintln!("[amet] cells: {}", manifest.len());

    // Intern group labels: the per-row aggregate key becomes two integers
    // (feature index, group id) instead of two freshly allocated strings.
    // A run has only a handful of distinct groups, so a linear scan is fine.
    let mut group_names: Vec<String> = Vec::new();
    for cell in &manifest {
        if !group_names.iter().any(|g| g == &cell.group) {
            group_names.push(cell.group.clone());
        }
    }

    let i_max_lag = cli.i_max_lag as usize;
    let min_n = cli.min_cpgs_per_feature;

    // Score cells in parallel. Each cell file is read once and scored against
    // every feature set; its cell_feature/pair_counts rows are formatted into a
    // thread-local buffer, then a brief per-set lock flushes the buffer and
    // folds the streaming aggregates. Nothing accumulates across cells, so peak
    // memory stays near `threads * one (cell, set) batch`. cell_feature and
    // pair_counts row order is therefore cell-interleaved; consumers key by
    // cell_id, and feature.tsv below is written in a sorted order.
    manifest.par_iter().try_for_each(|cell| -> Result<()> {
        let format = match cell.format.as_deref() {
            Some(s) => {
                CellFormat::parse(s).unwrap_or_else(|| CellFormat::detect_from_path(&cell.path))
            }
            None => CellFormat::detect_from_path(&cell.path),
        };
        let calls = match read_cell(&cell.path, format, &reference) {
            Ok(c) => c,
            Err(e) => {
                eprintln!("[amet] error reading {}: {}", cell.path.display(), e);
                return Ok(());
            }
        };
        // Every cell's group was interned from this same manifest above.
        let group_id = group_names
            .iter()
            .position(|g| g == &cell.group)
            .expect("cell group interned from manifest") as u32;

        for (set, sink) in sets.iter().zip(sinks.iter()) {
            // Cells are the parallel axis, so features are scored sequentially.
            let rows: Vec<CellFeatureRow> = set
                .features
                .iter()
                .map(|feature| {
                    let window = build_window(
                        feature,
                        &reference,
                        &calls,
                        cli.meth_call_threshold,
                        cli.min_reads_per_cpg,
                    );
                    let mc = marginal_counts(&window);
                    let pair_tables =
                        pair_counts_all_lags(&window, i_max_lag, cli.max_pair_distance);
                    let i_per_lag: Vec<f64> =
                        pair_tables.iter().map(amet::scores::i_total::i_k).collect();
                    CellFeatureRow {
                        feature_id: &feature.feature_id,
                        n_covered: window.n_observed() as u32,
                        n_zeros: mc.counts[0],
                        n_ones: mc.counts[1],
                        mean_meth: window.mean_meth(),
                        i_total_value: i_total(&pair_tables),
                        i_per_lag,
                        pair_tables,
                    }
                })
                .collect();

            // Format outside the lock so the critical section is just the
            // buffered write and the aggregate fold.
            let mut cf_buf: Vec<u8> = Vec::new();
            let mut pair_buf: Vec<u8> = Vec::new();
            for row in &rows {
                write_cell_feature_row(&mut cf_buf, cell, row, min_n, i_max_lag)?;
                write_pair_count_rows(&mut pair_buf, cell, row)?;
            }

            let mut guard = sink.lock().expect("sink mutex poisoned");
            guard.cf_writer.write_all(&cf_buf)?;
            guard.pair_writer.write_all(&pair_buf)?;
            for (feat_idx, row) in rows.iter().enumerate() {
                if row.n_covered >= min_n {
                    let e = guard.agg.entry((feat_idx as u32, group_id)).or_default();
                    e.cov_sum += row.n_covered as u64;
                    e.n_cells += 1;
                    if let Some(m) = row.mean_meth {
                        e.meth.add(m);
                    }
                    e.itotal.add(row.i_total_value);
                    e.jsd.add(&row.pair_tables[0]);
                }
            }
        }
        Ok(())
    })?;

    // Feature-level aggregates per set, written sorted by (feature_id, group).
    // Aggregates are keyed by integer indices; reconstruct the names here and
    // sort lexically so the output order is stable and independent of the
    // HashMap iteration order. The feature index breaks ties so the order is
    // fully determined even if two features share a feature_id.
    for ((set, sink), feat_writer) in sets.iter().zip(sinks).zip(feat_writers.iter_mut()) {
        let agg = sink.into_inner().expect("sink mutex poisoned").agg;
        let w: &mut dyn Write = &mut **feat_writer;
        let mut keys: Vec<(u32, u32)> = agg.keys().copied().collect();
        keys.sort_by(|a, b| {
            let ka = (
                set.features[a.0 as usize].feature_id.as_str(),
                group_names[a.1 as usize].as_str(),
                a.0,
            );
            let kb = (
                set.features[b.0 as usize].feature_id.as_str(),
                group_names[b.1 as usize].as_str(),
                b.0,
            );
            ka.cmp(&kb)
        });
        for key in keys {
            let e = &agg[&key];
            let feature_id = set.features[key.0 as usize].feature_id.as_str();
            let group = group_names[key.1 as usize].as_str();
            let mean_cov = if e.n_cells > 0 {
                e.cov_sum as f64 / e.n_cells as f64
            } else {
                0.0
            };
            let jsd = if e.n_cells >= cli.min_cells_per_group as u64 {
                Some(e.jsd.finish())
            } else {
                None
            };
            write!(
                w,
                "{}\t{}\t{}\t{:.6}\t",
                feature_id, group, e.n_cells, mean_cov
            )?;
            write_opt(w, e.meth.mean())?;
            write!(w, "\t")?;
            write_opt(w, e.meth.var())?;
            write!(w, "\t")?;
            write_opt(w, e.itotal.mean())?;
            write!(w, "\t")?;
            write_opt(w, e.itotal.var())?;
            write!(w, "\t")?;
            write_opt(w, jsd)?;
            writeln!(w)?;
        }

        eprintln!(
            "[amet] done {}: wrote {}, {}, {}",
            set.label,
            set.cf_path.display(),
            set.feat_path.display(),
            set.pair_path.display()
        );
    }
    Ok(())
}

/// Welford online accumulator for mean and sample variance. Used so the
/// feature-level aggregates need not retain every per-cell value.
#[derive(Default)]
struct Welford {
    n: u64,
    mean: f64,
    m2: f64,
}

impl Welford {
    fn add(&mut self, x: f64) {
        self.n += 1;
        let delta = x - self.mean;
        self.mean += delta / self.n as f64;
        let delta2 = x - self.mean;
        self.m2 += delta * delta2;
    }

    fn mean(&self) -> Option<f64> {
        if self.n == 0 { None } else { Some(self.mean) }
    }

    fn var(&self) -> Option<f64> {
        if self.n < 2 {
            None
        } else {
            Some(self.m2 / (self.n - 1) as f64)
        }
    }
}

/// Streaming feature-level aggregate for one (feature, group).
#[derive(Default)]
struct AggEntry {
    cov_sum: u64,
    n_cells: u64,
    meth: Welford,
    itotal: Welford,
    jsd: JsdAccumulator,
}

fn write_cell_feature_row(
    w: &mut dyn Write,
    cell: &CellRow,
    row: &CellFeatureRow<'_>,
    min_n: u32,
    i_max_lag: usize,
) -> std::io::Result<()> {
    write!(
        w,
        "{}\t{}\t{}\t{}\t",
        cell.cell_id, cell.group, row.feature_id, row.n_covered
    )?;
    match row.mean_meth {
        Some(m) => write!(w, "{:.6}", m)?,
        None => write!(w, "NA")?,
    }
    write!(w, "\t{}\t{}", row.n_zeros, row.n_ones)?;
    if row.n_covered >= min_n {
        write!(w, "\t{:.6}", row.i_total_value)?;
        for v in &row.i_per_lag {
            write!(w, "\t{:.6}", v)?;
        }
    } else {
        write!(w, "\tNA")?;
        for _ in 0..i_max_lag {
            write!(w, "\tNA")?;
        }
    }
    writeln!(w)
}

fn write_pair_count_rows(
    w: &mut dyn Write,
    cell: &CellRow,
    row: &CellFeatureRow<'_>,
) -> std::io::Result<()> {
    for (idx, pt) in row.pair_tables.iter().enumerate() {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            cell.cell_id,
            cell.group,
            row.feature_id,
            idx + 1,
            pt.counts[0],
            pt.counts[1],
            pt.counts[2],
            pt.counts[3]
        )?;
    }
    Ok(())
}

struct FeatureSet {
    label: String,
    features: Vec<Feature>,
    cf_path: PathBuf,
    feat_path: PathBuf,
    pair_path: PathBuf,
}

/// Per-set mutable state shared across the parallel cell workers behind a
/// `Mutex`. The feature.tsv writer is not here: it is touched only before and
/// after the parallel section, never concurrently.
struct SetSink {
    cf_writer: Box<dyn Write + Send>,
    pair_writer: Box<dyn Write + Send>,
    /// Keyed by (feature index within the set, interned group id) so the hot
    /// per-row fold does not allocate a fresh key string each time.
    agg: HashMap<(u32, u32), AggEntry>,
}

/// One scored (cell, feature). `feature_id` borrows from the feature set, which
/// outlives the short-lived per-(cell, set) batch this row belongs to.
struct CellFeatureRow<'a> {
    feature_id: &'a str,
    n_covered: u32,
    n_zeros: u32,
    n_ones: u32,
    mean_meth: Option<f64>,
    i_total_value: f64,
    i_per_lag: Vec<f64>,
    pair_tables: Vec<PairCounts>,
}

fn write_headers(
    cf_writer: &mut dyn Write,
    feat_writer: &mut dyn Write,
    pair_writer: &mut dyn Write,
    i_max_lag: u32,
) -> std::io::Result<()> {
    write!(
        cf_writer,
        "cell_id\tgroup\tfeature_id\tn_covered\tmean_meth\tn_zeros\tn_ones\ti_total"
    )?;
    for k in 1..=i_max_lag {
        write!(cf_writer, "\ti_{}", k)?;
    }
    writeln!(cf_writer)?;
    writeln!(
        feat_writer,
        "feature_id\tgroup\tn_cells\tmean_coverage\tmean_meth_mean\tmean_meth_var\ti_total_mean\ti_total_var\tjsd"
    )?;
    writeln!(
        pair_writer,
        "cell_id\tgroup\tfeature_id\tlag\tn00\tn01\tn10\tn11"
    )?;
    Ok(())
}

/// Derive a stable label for a BED path by stripping known BED and compression
/// suffixes if present; otherwise use the file name as-is. Used to disambiguate
/// output paths when multiple --features are supplied.
fn features_label(path: &Path) -> String {
    let raw = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| "features".to_string());
    for suffix in [".bed.gz", ".bed.bgz", ".bed", ".gz", ".bgz"] {
        if let Some(s) = raw.strip_suffix(suffix) {
            return s.to_string();
        }
    }
    raw
}

fn with_suffix(prefix: &Path, suffix: &str) -> PathBuf {
    let mut s = prefix.as_os_str().to_owned();
    s.push(suffix);
    PathBuf::from(s)
}

fn write_opt<W: Write + ?Sized>(w: &mut W, x: Option<f64>) -> std::io::Result<()> {
    match x {
        Some(v) => write!(w, "{:.6}", v),
        None => write!(w, "NA"),
    }
}

#[cfg(test)]
mod tests {
    use super::features_label;
    use std::path::Path;

    #[test]
    fn label_strips_bed_gz() {
        assert_eq!(
            features_label(Path::new("/x/promoters.bed.gz")),
            "promoters"
        );
    }

    #[test]
    fn label_strips_bed() {
        assert_eq!(features_label(Path::new("/x/enhancers.bed")), "enhancers");
    }

    #[test]
    fn label_strips_gz_only() {
        assert_eq!(
            features_label(Path::new("/x/regions.tsv.gz")),
            "regions.tsv"
        );
    }

    #[test]
    fn label_keeps_unknown_extension() {
        assert_eq!(features_label(Path::new("/x/regions.txt")), "regions.txt");
    }

    #[test]
    fn label_handles_multiple_dots() {
        // Common case: dataset-tagged bed names like "mm10.heterochromatin.bed".
        assert_eq!(
            features_label(Path::new("/x/mm10.heterochromatin.bed")),
            "mm10.heterochromatin"
        );
    }

    #[test]
    fn label_strips_bed_bgz() {
        assert_eq!(features_label(Path::new("/x/regions.bed.bgz")), "regions");
    }

    #[test]
    fn label_strips_bgz_only() {
        assert_eq!(
            features_label(Path::new("/x/regions.tsv.bgz")),
            "regions.tsv"
        );
    }
}

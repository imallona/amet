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

    // Resolve per-set labels and check uniqueness up front so we fail before any I/O.
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
    let mut outputs: Vec<SetOutputs> = Vec::with_capacity(cli.features.len());
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
        });
        outputs.push(SetOutputs {
            cf_path,
            feat_path,
            pair_path,
            cf_writer,
            feat_writer,
            pair_writer,
            agg: HashMap::new(),
        });
    }

    eprintln!("[amet] reading manifest: {}", cells_path.display());
    let manifest = read_manifest(cells_path, &cli.group_column).context("reading manifest")?;
    eprintln!("[amet] cells: {}", manifest.len());

    let i_max_lag = cli.i_max_lag as usize;
    let min_n = cli.min_cpgs_per_feature;

    // Stream per cell: each cell file is read once and scored against every
    // feature set. Its cell_feature and pair_counts rows are written out and
    // dropped immediately, so peak memory stays near one (cell, set) batch
    // rather than the full cells x sets x features matrix. Feature-level
    // aggregates are folded incrementally (Welford for mean/variance,
    // JsdAccumulator for jsd), so no per-cell values are retained.
    for cell in &manifest {
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
                continue;
            }
        };

        for (set, out) in sets.iter().zip(outputs.iter_mut()) {
            // Score this cell's features in parallel; collect preserves feature order.
            let rows: Vec<CellFeatureRow> = set
                .features
                .par_iter()
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

            for row in &rows {
                write_cell_feature_row(&mut out.cf_writer, cell, row, min_n, i_max_lag)?;
                write_pair_count_rows(&mut out.pair_writer, cell, row)?;
                if row.n_covered >= min_n {
                    let key = (row.feature_id.to_string(), cell.group.clone());
                    let e = out.agg.entry(key).or_default();
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
    }

    // Feature-level aggregates per set, written sorted by (feature_id, group).
    for (set, out) in sets.iter().zip(outputs.iter_mut()) {
        let mut keys: Vec<_> = out.agg.keys().cloned().collect();
        keys.sort();
        for key in keys {
            let e = &out.agg[&key];
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
                out.feat_writer,
                "{}\t{}\t{}\t{:.6}\t",
                key.0, key.1, e.n_cells, mean_cov
            )?;
            write_opt(&mut out.feat_writer, e.meth.mean())?;
            write!(out.feat_writer, "\t")?;
            write_opt(&mut out.feat_writer, e.meth.var())?;
            write!(out.feat_writer, "\t")?;
            write_opt(&mut out.feat_writer, e.itotal.mean())?;
            write!(out.feat_writer, "\t")?;
            write_opt(&mut out.feat_writer, e.itotal.var())?;
            write!(out.feat_writer, "\t")?;
            write_opt(&mut out.feat_writer, jsd)?;
            writeln!(out.feat_writer)?;
        }

        eprintln!(
            "[amet] done {}: wrote {}, {}, {}",
            set.label,
            out.cf_path.display(),
            out.feat_path.display(),
            out.pair_path.display()
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

/// Streaming feature-level aggregate for one (feature_id, group).
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
}

struct SetOutputs {
    cf_path: PathBuf,
    feat_path: PathBuf,
    pair_path: PathBuf,
    cf_writer: Box<dyn Write>,
    feat_writer: Box<dyn Write>,
    pair_writer: Box<dyn Write>,
    agg: HashMap<(String, String), AggEntry>,
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

/// Derive a stable label for a BED path: strip a compressed-BED extension if any,
/// otherwise use the file name as-is. Used to disambiguate output paths when
/// multiple --features are supplied. The list of suffixes mirrors what
/// `crate::io::open_read` accepts as compressed input.
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

use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(name = "amet")]
#[command(about = "p-decoupled information-theoretic scores for single-cell DNA methylation", long_about = None)]
#[command(version)]
pub struct Cli {
    /// Tab-separated manifest with columns cell_id, group, path, plus optional extras.
    /// `format` column overrides per-cell format auto-detection.
    #[arg(long, value_name = "TSV")]
    pub cells: PathBuf,

    /// BED file of features to score. Features should not overlap.
    #[arg(long, value_name = "BED")]
    pub features: PathBuf,

    /// CpG reference TSV (chrom\tpos), 0-based positions of every CpG to consider.
    #[arg(long, value_name = "TSV")]
    pub cpg_reference: PathBuf,

    /// Output file prefix.
    #[arg(long, value_name = "PREFIX")]
    pub output_prefix: PathBuf,

    /// Manifest column name to use for grouping.
    #[arg(long, default_value = "group")]
    pub group_column: String,

    /// Methylation fraction threshold for binarising m/t to 0/1. Use 0.5 for majority
    /// rule, or e.g. 0.1 to call any position with > 10 percent methylated reads as 1.
    #[arg(long, default_value_t = 0.0)]
    pub meth_call_threshold: f64,

    /// Minimum reads required to consider a CpG observed. Default 1 fits single-cell
    /// data; bulk WGBS users typically set 5-10.
    #[arg(long, default_value_t = 1)]
    pub min_reads_per_cpg: u32,

    /// Minimum covered CpGs required to compute scores for a (cell, feature).
    #[arg(long, default_value_t = 5)]
    pub min_cpgs_per_feature: u32,

    /// Minimum number of cells required to report JSD for a (feature, group).
    #[arg(long, default_value_t = 10)]
    pub min_cells_per_group: u32,

    /// Maximum CpG lag k for the I_total within-cell score: I_total = sum_{k=1..max} I_k.
    #[arg(long, default_value_t = 3)]
    pub i_max_lag: u32,

    /// Number of threads. 0 means all available.
    #[arg(long, default_value_t = 0)]
    pub threads: usize,
}

impl Cli {
    pub fn parse_args() -> Self {
        Self::parse()
    }
}

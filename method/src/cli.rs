use clap::{ArgGroup, Parser};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(name = "amet")]
#[command(about = "p-decoupled information-theoretic scores for single-cell DNA methylation", long_about = None)]
#[command(version)]
#[command(group(ArgGroup::new("cpgs").required(true).args(["genome", "cpg_reference"])))]
pub struct Cli {
    /// Tab-separated manifest with columns cell_id, group, path, plus optional extras.
    /// `format` column overrides per-cell format auto-detection.
    #[arg(long, value_name = "TSV", required_unless_present = "build_cpg_only")]
    pub cells: Option<PathBuf>,

    /// BED file of features to score. Features within a single BED should not overlap.
    /// Pass --features multiple times to score the same cells against several feature
    /// sets in one cell-read pass; each set writes its own output triplet keyed by the
    /// BED basename. With a single --features the output paths are unchanged.
    #[arg(
        long,
        value_name = "BED",
        action = clap::ArgAction::Append,
        required_unless_present = "build_cpg_only"
    )]
    pub features: Vec<PathBuf>,

    /// FASTA of the reference genome. amet derives all CpG positions from it on first
    /// use and caches them to <fasta>.cpg next to the input. Subsequent runs reuse the
    /// cache. Mutually exclusive with --cpg-reference.
    #[arg(long, value_name = "FASTA")]
    pub genome: Option<PathBuf>,

    /// CpG reference TSV (chrom\tpos), 0-based positions of every CpG to consider.
    /// Mutually exclusive with --genome.
    #[arg(long, value_name = "TSV")]
    pub cpg_reference: Option<PathBuf>,

    /// Output file prefix.
    #[arg(
        long,
        value_name = "PREFIX",
        required_unless_present = "build_cpg_only"
    )]
    pub output_prefix: Option<PathBuf>,

    /// Build the <fasta>.cpg index and exit. Requires --genome. Used by snakemake to
    /// materialise the cache once before fanning out concurrent scoring jobs.
    #[arg(long, requires = "genome", conflicts_with = "cpg_reference")]
    pub build_cpg_only: bool,

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
    /// Must be at least 1; lag 1 is required to compute JSD.
    #[arg(long, default_value_t = 3, value_parser = clap::value_parser!(u32).range(1..))]
    pub i_max_lag: u32,

    /// Maximum nucleotide distance allowed between paired CpGs. Pairs whose genomic
    /// distance exceeds this value are not counted. 0 disables the cap.
    #[arg(long, default_value_t = 0)]
    pub max_pair_distance: u64,

    /// Number of threads. 0 means all available.
    #[arg(long, default_value_t = 0)]
    pub threads: usize,
}

impl Cli {
    pub fn parse_args() -> Self {
        Self::parse()
    }
}

# amet

p-decoupled information-theoretic scores for single-cell DNA methylation heterogeneity.

amet quantifies two complementary axes of methylation pattern variability:

- **Within-cell regularity** along consecutive CpGs in one cell, scored by `I_total`, the sum of mutual information across CpG lags 1..k.
- **Across-cell heterogeneity** at a feature within a cell group, scored by `JSD` on per-cell L-mer distributions.

Both scores are decoupled from marginal methylation at the no-signal baseline: a sequence with no comethylation structure scores zero regardless of its methylation level, so no marginal-methylation adjustment is needed.

## Status

Early prototype (v0.1).

## Repository layout

```
amet/
├── method/                    Rust crate (the amet binary and library)
│   ├── Cargo.toml
│   ├── src/                   parsers, scores, CLI, I/O
│   └── tests/                 integration tests
├── workflow/                  Snakemake workflow for simulations and dataset analyses
│   ├── Snakefile
│   ├── config/sim.yaml        simulation parameters
│   ├── envs/r.yml             conda env for the R scripts
│   ├── scripts/               R scripts (data generation, evaluation)
│   └── Rmd/                   reports
├── results/                   gitignored: outputs of running the workflow
├── .github/workflows/         CI definitions
├── README.md
├── LICENSE
└── CLAUDE.md                  design notes for the development environment
```

## Build

```
cd method
cargo build --release
```

The binary lives at `method/target/release/amet`.

## Quick start

```
amet \
  --cpg-reference cpgs.tsv.gz \
  --cells cells.tsv \
  --features features.bed \
  --output-prefix run1
```

`cells.tsv` has columns `cell_id`, `group`, `path` (tab-separated). `features.bed` is a standard BED file. Outputs are written as `run1.cell_feature.tsv.gz` and `run1.feature.tsv.gz`.

See `amet --help` for full options.

## Scores

### Within-cell: `I_total`

For a single cell, with CpG calls along a feature, compute mutual information between CpGs separated by lag k:

```
I_k = H(X_i) + H(X_{i+k}) - H(X_i, X_{i+k})
I_total = sum over k=1..k_max of I_k
```

For an i.i.d. sequence with marginal p, every I_k = 0 regardless of p, so I_total has a p-invariant zero baseline. No adjustment for marginal methylation is needed.

### Across-cell, per group: `JSD`

For each cell, build a 2-mer histogram per feature (4 bins: 00, 01, 10, 11). Compute Jensen-Shannon divergence across the cells in each group. JSD is reported per (feature, group), not pooled.

## Input formats

- `allc` / methylpy.
- `scNMT` cpg_level.

Format is auto-detected from the filename or set explicitly via a `format` column in the manifest.

## Simulator

The `workflow/` directory contains a Snakemake pipeline that generates synthetic single-cell methylation data with known ground truth (Markov-chain single cells, mixtures of Markov chains across cells, fragment-based sparsity), runs amet on it, and produces evaluation metrics. Run with:

```
snakemake -j 4 -s workflow/Snakefile
```

## License

GPL-3.0-or-later.

## Contact

Izaskun Mallona — [izaskun.mallona.work@gmail.com](mailto:izaskun.mallona.work@gmail.com)

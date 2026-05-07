# amet

amet is a tool to quantify within- and across-cells epigenetic heterogeneity using single-cell DNA methylation data.

It produces two complementary scores:

- Within-cell regularity along consecutive CpGs in one cell, scored by `I_total`, the sum of mutual information across CpG lags 1..k.
- Across-cell heterogeneity at a feature within a cell group, scored by `JSD` on per-cell lag-1 2-mer distributions.

A sequence with no comethylation structure scores zero regardless of its methylation level, so no marginal-methylation adjustment is needed.

## Status

Early prototype (v0.1).

## Repository layout

```
amet/
  method/                    Rust crate (the amet binary and library)
    Cargo.toml
    src/                     parsers, scores, CLI, I/O
    tests/                   integration tests
  workflow/                  Snakemake workflow for simulations and dataset analyses
    Snakefile
    config/sim.yaml          simulation parameters
    config/datasets.yaml     dataset paths and prototype subsets
    envs/                    conda envs (rust, bedtools, r-tools, python)
    scripts/                 R / Python / shell helpers
    Rmd/                     reports
  results/                   gitignored: outputs of running the workflow
  .github/workflows/         CI definitions
  README.md
  LICENSE
  AUTHORS
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
  --genome mm10.fa \
  --cells cells.tsv \
  --features features.bed \
  --output-prefix run1
```

Outputs land at `run1.cell_feature.tsv.gz` and `run1.feature.tsv.gz`.

On the first run amet derives every CpG position from the FASTA and writes a sidecar `mm10.fa.cpg` next to the input. Subsequent runs reuse the sidecar. If you already have a CpG list, pass `--cpg-reference cpgs.tsv.gz` instead. Exactly one of `--genome` or `--cpg-reference` is required.

## CLI

| Flag | Required | Default | Description |
|---|---|---|---|
| `--genome` | one of these two | (required) | FASTA of the reference genome. amet derives all CpG positions on first use and caches them to `<fasta>.cpg`. |
| `--cpg-reference` | one of these two | (required) | Tab-separated `chrom\tpos` of every CpG to consider, 0-based. Defines adjacency: any uncovered reference CpG breaks 2-mer pairing across it. |
| `--cells` | yes | (required) | Manifest TSV (see below). |
| `--features` | yes | (required) | Standard BED file of regions to score. Features should not overlap. |
| `--output-prefix` | yes | (required) | Prefix for the two output files. |
| `--group-column` | no | `group` | Manifest column to use as the group label. |
| `--meth-call-threshold` | no | `0.0` | Methylation fraction `m/t` above which a CpG is called methylated. `0.5` for majority rule; `0.1` calls any position with > 10% methylated reads as methylated. |
| `--min-reads-per-cpg` | no | `1` | A CpG is observed only if covered by at least this many reads. Default fits single-cell norms; bulk WGBS users typically set 5-10. |
| `--min-cpgs-per-feature` | no | `5` | A `(cell, feature)` is scored only if at least this many CpGs are covered. Below the threshold, scores are reported as `NA`. |
| `--min-cells-per-group` | no | `10` | A `(feature, group)` reports `jsd` only if at least this many cells pass the per-cell coverage filter. Otherwise `jsd` is `NA`. |
| `--i-max-lag` | no | `3` | Maximum CpG lag k for `I_total = sum_{k=1..max} I_k`. |
| `--threads` | no | `0` (all) | Number of threads. |

### Manifest (`--cells`)

Tab-separated with a header row. Required columns: `cell_id`, `path`. Optional columns: `group` (default `all` with stderr warning), `format` (overrides per-cell format auto-detection). Any extra columns are ignored by amet but preserved in your records.

```
cell_id    group         path                              batch  donor
A12        excitatory    data/ecker/A12.allc.tsv.gz        b1     d1
B07        excitatory    data/ecker/B07.allc.tsv.gz        b1     d1
C03        inhibitory    data/ecker/C03.allc.tsv.gz        b2     d1
```

### Input formats per cell

- `allc` / methylpy: 7+ columns, per-strand or pre-collapsed, `CG`-context filter applied. Default for unrecognised filenames.
- `scNMT` cpg_level: 5 columns with header, `pos` is the 1-based G position.

amet auto-detects from the filename: `*.allc.*` -> allc, `*.cpg_level.*` or `*.scnmt.*` -> scNMT, other -> allc. Override per cell via the `format` column in the manifest (`allc`, `methylpy`, `scnmt`, `cpg_level`).

### CpG reference (`--cpg-reference`)

Two columns, `chrom` and 0-based start of every CpG you want amet to consider (gzip-OK). Must be sorted ascending within each chromosome. Lines starting with `#` are ignored.

```
chr1    9
chr1    19
chr1    127
chr2    34
```

### Features (`--features`)

Standard 3-or-4-column BED. The 4th column, if present, is used as `feature_id`; otherwise amet falls back to `chrom:start-end`. Features are 0-based half-open and must not overlap.

```
chr1    1000    2000    promoter_GENE1
chr1    5000    7000    cgi_chr1_5000
```

### Output: `<prefix>.cell_feature.tsv.gz`

One row per `(cell, feature)`:

```
cell_id  group  feature_id  n_covered  mean_meth  i_total  i_1  i_2  ... i_k
```

`i_k` columns are present up to `--i-max-lag`. When fewer than `--min-cpgs-per-feature` CpGs are covered, score columns are `NA`.

### Output: `<prefix>.feature.tsv.gz`

One row per `(feature, group)`:

```
feature_id  group  n_cells  mean_coverage  jsd
```

JSD is computed across the cells in each group using each cell's lag-1 2-mer distribution. Pooled JSD across groups is not reported; per-group is the only meaningful axis for heterogeneity. If a group has fewer than `--min-cells-per-group` eligible cells, `jsd` is reported as `NA`.

## Scores

### Within-cell: `I_total`

For a single cell, with CpG calls along a feature, compute mutual information between CpGs separated by lag k:

```
I_k = H(X_i) + H(X_{i+k}) - H(X_i, X_{i+k})
I_total = sum over k=1..k_max of I_k
```

For an i.i.d. sequence with marginal p, every I_k = 0 regardless of p, so I_total has a p-invariant zero baseline. No adjustment for marginal methylation is needed.

### Across-cell, per group: `JSD` (lag-1 2-mer)

For each cell, build a lag-1 2-mer histogram per feature (4 bins: 00, 01, 10, 11). Compute Jensen-Shannon divergence across the cells in each group. JSD is reported per (feature, group), not pooled.

## Simulator

The `workflow/` directory contains a Snakemake pipeline that generates synthetic single-cell methylation data with known ground truth (Markov-chain single cells, mixtures of Markov chains across cells, fragment-based sparsity), runs amet on it, and produces evaluation metrics. Run with:

```
snakemake -j 4 -s workflow/Snakefile simulations
```

## License

GPL-3.0-or-later.

## Contact

Izaskun Mallona - [izaskun.mallona.work@gmail.com](mailto:izaskun.mallona.work@gmail.com)

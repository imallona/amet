# amet

amet quantifies within- and across-cells epigenetic heterogeneity from single-cell DNA methylation data.

Two complementary scores:

- Within-cell regularity along consecutive CpGs in one cell, scored by `I_total`, the sum of mutual information across CpG lags 1..k.
- Across-cell heterogeneity at a feature within a cell group, scored by `JSD` on per-cell lag-1 2-mer distributions.

A sequence with no comethylation structure scores zero regardless of its methylation level, so no marginal-methylation adjustment is needed.

## Status

v0.1. The Rust binary is feature-complete for the two scores. The Snakemake workflow runs the simulator-based validation grid and three reference dataset analyses (Argelaguet 2019 mouse gastrulation scNMT-seq, Liu 2021 mouse MOp snmC-seq2, Bian 2018 colorectal scTrio-seq2).

## Repository layout

```
amet/
  Makefile                     entry points for the three dataset analyses + simulations
  method/                      Rust crate
    Cargo.toml
    src/                       parsers, scores, CLI, I/O, CpG-index cache
    tests/                     integration + snapshot tests
  workflow/                    Snakemake pipeline
    Snakefile                  simulations + dispatch to the three dataset rule files
    config/
      sim.yaml                 simulation parameters
      datasets.yaml            dataset paths, prototype subsets, window sizes
    envs/                      conda envs (rust, bedtools, r-tools, python)
    rules/
      common.smk               build_amet, fetch FASTAs, build_cpg_reference
      argelaguet.smk           Argelaguet 2019 gastrulation rules
      crc.smk                  CRC scTrio-seq2 rules
      ecker.smk                Ecker MOp snmC-seq2 rules
    Rmd/                       per-dataset and simulation HTML report sources
    scripts/                   R, Python, and shell helpers (manifest builders,
                               BED stagers, simulator generators, eval scripts)
  results/                     gitignored: workflow outputs (one subdir per task)
  .github/workflows/           CI: cargo fmt, clippy, test
  README.md, LICENSE, AUTHORS
```

## Build

```
cd method
cargo build --release
```

The binary lives at `method/target/release/amet`. CI runs `cargo fmt --check`, `cargo clippy -- -D warnings`, and `cargo test` on every PR.

## Quick start

One-shot scoring of a feature BED across cells in a manifest:

```
amet \
  --genome mm10.fa \
  --cells cells.tsv \
  --features features.bed \
  --output-prefix run1
```

Outputs land at `run1.cell_feature.tsv.gz`, `run1.feature.tsv.gz`, and `run1.pair_counts.tsv.gz`.

On the first run amet derives every CpG position from the FASTA and writes a sidecar `mm10.fa.cpg` next to the input. Subsequent runs reuse the sidecar. The build is protected by an advisory `flock` and atomic rename, so launching many amet jobs in parallel against the same FASTA is safe (no torn cache).

If you already have a CpG list, pass `--cpg-reference cpgs.tsv` instead. Exactly one of `--genome` or `--cpg-reference` is required.

To pre-materialise the CpG cache without running a scoring job:

```
amet --build-cpg-only --genome mm10.fa
```

## Running the dataset analyses

The `Makefile` is the top-level entry point. From the repo root:

```
make argelaguet      # whole-genome Argelaguet (4 Rmds)
make crc             # whole-genome CRC (6 Rmds)
make ecker           # whole-genome Ecker (4 Rmds)
make simulations     # simulation report
make all             # everything above
make dryrun          # snakemake -n for all targets
make unlock          # release a stale snakemake lock
make help            # target list
```

Tunable variables:

| Variable | Default | Description |
|---|---|---|
| `CORES` | 16 | Snakemake `--cores` value |
| `ULIMIT_KB` | 209715200 (200 GB) | Virtual-memory cap; inherited by every amet job |
| `CONDA_ENV` | `snakemake` | Conda env that holds snakemake |
| `CONDA_INIT` | `~/miniconda3/bin/activate` | Conda activation script |

Override on the command line, e.g. `make argelaguet CORES=32`.

### Prototype vs full-run

`workflow/config/datasets.yaml` carries a `prototype:` block. With `prototype.enabled: true` (the default), each dataset's manifest builder caps cells per stratum and restricts to a small list of strata so the runs finish in minutes. Output goes to `results/<dataset>_proto/`.

For full publication runs, set `prototype.enabled: false` and change each dataset's `run_name` (e.g., from `argelaguet_proto` to `argelaguet`) so the full outputs don't overwrite the prototype ones.

### Window sizes

Set per dataset under `<dataset>.window_size` in `datasets.yaml`. Defaults: argelaguet 500 kb (broad gastrulation territories), ecker 10 kb (cortical sub-types), crc 10 kb (tumour windows).

### Sex and mitochondrial chromosomes

Every BED-producing rule filters out X, Y, M, and MT contigs. The CpG reference itself stays complete (it is a property of the genome). Only the regions amet scores are filtered.

### Annotation BEDs

The three dataset rule files expand over a fixed list of annotations defined at the top of each `.smk` file. argelaguet uses the ENCODE ChIP marks (h3k27ac/h3k27me3/h3k9me3/h3k4me1/h3k4me3), Gencode genes and promoters, UCSC RMSK LINEs/SINEs, plus the gastrulation-specific BEDs bundled in the scnmt_gastrulation tarball (E7.5 H3K27ac enhancers, E7.5 H3K4me3, ESC marks). Ecker uses the five ENCODE ChIP marks plus genes/promoters/lines/sines. CRC uses ENCODE ChIP (H3K27me3, H3K9me3, H3K4me3), Segway chromHMM segments, laminB1 LADs, common PMDs and HMDs, UCSC CpG islands, Gencode genes, RMSK repeats, and the patient-specific SCNA tracks.

### Server deployment

The Makefile is designed for a workstation with enough RAM for the whole-genome amet runs (hundreds of GB virtual memory under parallel jobs; the recipes apply `ulimit -v` as a soft cap). It is not designed for laptops.

If you are in the Mark Robinson lab at UZH, `workflow/scripts/internal/setup_barbara_links.sh` populates `results/<dataset>/{cells,raw,features,mm10,hg19}` as symlinks to the pre-staged data tree on `barbara`'s filesystem. See `workflow/scripts/internal/README.md`. Outside that lab, run the per-dataset download rules in each `.smk` instead.

## CLI

| Flag | Default | Description |
|---|---|---|
| `--genome` | (required, mutually exclusive with `--cpg-reference`) | FASTA of the reference genome. amet derives all CpG positions on first use and caches them to `<fasta>.cpg`. |
| `--cpg-reference` | (required, mutually exclusive with `--genome`) | Tab-separated `chrom\tpos` of every CpG, 0-based. Defines adjacency: an uncovered reference CpG breaks 2-mer pairing across it. |
| `--cells` | (required unless `--build-cpg-only`) | Manifest TSV (see below). |
| `--features` | (required unless `--build-cpg-only`) | BED of regions to score. Features must not overlap. |
| `--output-prefix` | (required unless `--build-cpg-only`) | Prefix for the output files. |
| `--build-cpg-only` | off | Only materialise `<fasta>.cpg` and exit. Requires `--genome`. |
| `--group-column` | `group` | Manifest column to use as the group label. |
| `--meth-call-threshold` | `0.0` | Methylation fraction `m/t` above which a CpG is called methylated. `0.5` for majority rule; `0.1` calls any position with > 10% methylated reads as methylated. |
| `--min-reads-per-cpg` | `1` | A CpG is observed only if covered by at least this many reads. Bulk WGBS users typically set 5-10. |
| `--min-cpgs-per-feature` | `5` | A `(cell, feature)` is scored only if at least this many CpGs are covered. Below the threshold, scores are reported as `NA`. |
| `--min-cells-per-group` | `10` | A `(feature, group)` reports `jsd` only if at least this many cells pass the per-cell coverage filter. Otherwise `jsd` is `NA`. |
| `--i-max-lag` | `3` | Maximum CpG lag k for `I_total = sum_{k=1..max} I_k`. |
| `--threads` | `0` (all) | Number of threads. |

### Manifest (`--cells`)

Tab-separated with a header row. Required columns: `cell_id`, `path`. Optional: `group` (default `all`), `format` (overrides per-cell format auto-detection). Extra columns are ignored by amet but preserved for your records.

```
cell_id    group         path                              batch  donor
A12        excitatory    data/A12.allc.tsv.gz              b1     d1
B07        excitatory    data/B07.allc.tsv.gz              b1     d1
C03        inhibitory    data/C03.allc.tsv.gz              b2     d1
```

### Cell formats

| Format | Manifest value(s) | Filename auto-detect | Description |
|---|---|---|---|
| allc / methylpy | `allc`, `methylpy` | default (no marker) | 7+ column allc, per-strand or pre-collapsed. CG context only. |
| scNMT cpg_level | `scnmt`, `cpg_level` | `*.cpg_level.*`, `*.scnmt.*` | 5-column with header; `pos` is the 1-based G position. |
| Bismark singleC | `bismark`, `singlec` | `*.bismark.*`, `*.singlec.*` | 10-column Bismark singleC; CG-only filter applied. |

Override per cell via the `format` column in the manifest.

### CpG reference (`--cpg-reference`)

Two columns, `chrom` and 0-based start of every CpG (gzip-OK). Must be sorted ascending within each chromosome. `#` lines are ignored.

```
chr1    9
chr1    19
chr1    127
chr2    34
```

### Features (`--features`)

3-or-4-column BED. The 4th column, if present, is `feature_id`; otherwise amet falls back to `chrom:start-end`. Features are 0-based half-open and must not overlap.

```
chr1    1000    2000    promoter_GENE1
chr1    5000    7000    cgi_chr1_5000
```

### Outputs

`<prefix>.cell_feature.tsv.gz`. One row per `(cell, feature)`:

```
cell_id  group  feature_id  n_covered  mean_meth  n_zeros  n_ones  i_total  i_1  i_2  ... i_k
```

`i_k` columns are present up to `--i-max-lag`. When fewer than `--min-cpgs-per-feature` CpGs are covered, score columns are `NA`.

`<prefix>.feature.tsv.gz`. One row per `(feature, group)`:

```
feature_id  group  n_cells  mean_coverage  mean_meth_mean  mean_meth_var  i_total_mean  i_total_var  jsd
```

`jsd` is the multi-distribution Jensen-Shannon divergence across the cells in each group, using each cell's lag-1 2-mer distribution. If a group has fewer than `--min-cells-per-group` eligible cells, `jsd` is `NA`. Per-group is the only meaningful axis; pooled JSD is not reported.

`<prefix>.pair_counts.tsv.gz`. One row per `(cell, feature, lag)` with the four 2-mer counts:

```
cell_id  group  feature_id  lag  n00  n01  n10  n11
```

Useful for downstream analyses that want to recompute scores under alternative thresholds, or to diagnose individual cells.

## Scores

### Within-cell: `I_total`

For a single cell, with binarised CpG calls along a feature, compute mutual information between CpGs separated by lag k:

```
I_k     = H(X_i) + H(X_{i+k}) - H(X_i, X_{i+k})
I_total = sum_{k=1..k_max} I_k
```

For an i.i.d. sequence with marginal p, every `I_k = 0` for any p, so `I_total` has a p-invariant zero baseline. No methylation-level adjustment is needed.

`I_total` carries practical methylation dependence in real data (cells with mostly homozygous calls have less information to capture). The zero-baseline property covers the null only; for differential testing across conditions with different marginal methylation, residualise on `mean_meth`.

### Across-cell, per group: `JSD` (lag-1 2-mer)

For each cell, build a lag-1 2-mer histogram per feature (4 bins: 00, 01, 10, 11). Compute multi-distribution Jensen-Shannon divergence across the cells in each group. Reported per `(feature, group)`, never pooled.

## Simulations

`workflow/Snakefile` runs a validation grid that generates synthetic single-cell methylation data with known ground truth (Markov-chain single cells, Markov mixtures across cells, fragment-based sparsity), scores it with amet, and produces an HTML evaluation report covering p-decoupling, lag profile, sparsity, feature length, n-cells, mixture-k, mixture-divergence, wcVI, acVI, feature variability, and consensus perturbation.

Run with:

```
make simulations
```

Outputs go to `results/simulations/`.

## License

GPL-3.0-or-later.

## Contact

Izaskun Mallona - [izaskun.mallona.work@gmail.com](mailto:izaskun.mallona.work@gmail.com)

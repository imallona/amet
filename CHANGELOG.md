# Changelog

Notable changes to amet. Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versioning: [SemVer](https://semver.org/spec/v2.0.0.html).

## [0.1.1] - 2026-05-15

### Added

- `--features` accepts multiple BEDs, scored in one cell-read pass. Each BED
  writes its own triplet `<prefix>.<label>.{cell_feature,feature,pair_counts}.tsv.gz`
  (`<label>` is the BED basename). A single `--features` keeps the old paths.

### Changed

- Streamed scoring: rows written per cell, aggregates updated as it goes. Peak
  memory is bounded by `threads` x one cell-and-feature-set batch, not the full
  cells x feature-sets x features matrix.
- Cells scored in parallel. `cell_feature`/`pair_counts` row order is now
  cell-interleaved (key by `cell_id`); `feature.tsv.gz` stays sorted.
- Workflow: CRC/Ecker/Argelaguet feature rules run amet once per cell-group
  combo. Output names changed from `<annotation>_<combo>.*` to
  `<combo>.<annotation>.*`.
- Makefile: `CORES` 16 -> 40, `ULIMIT_KB` 200 GB -> 100 GB per process.

### Performance

- Pair counting rewritten: compact observed-CpG window, all lags in one sweep.
  Cost scales with observed CpGs, not the feature's total CpG count.

### Fixed

- `--i-max-lag 0` is rejected by the CLI; it no longer panics (lag 1 is needed
  for JSD).
- `features_label` strips `.bgz` and `.bed.bgz`.
- `setup_barbara_links.sh` and `sync_from_barbara.sh` resolved the repo root one
  level too shallow, placing `results/` symlinks under `workflow/`.

## [0.1.0]

Initial release: within-cell `I_total` and across-cell `JSD` scores; allc,
scNMT cpg_level and Bismark singleC parsers; FASTA-derived CpG reference with
sidecar cache; Snakemake workflow for the simulation grid and the CRC, Ecker
and Argelaguet datasets.

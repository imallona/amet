# Internal scripts (Mark Robinson lab, UZH)

These scripts are for the Mark Robinson lab at UZH only. They wire amet up against the pre-staged data tree on `barbara`'s filesystem and reference absolute paths specific to that host. They are not part of the public amet workflow.

Do not run these outside the Mark Robinson lab. The absolute paths under `/home/imallona/src/yamet/` on `barbara` will not exist elsewhere.

| Script | Purpose |
|---|---|
| `setup_barbara_links.sh` | Symlink the dataset-staging directories under `results/<dataset>/` to the existing pre-staged tree. Idempotent. |
| `sync_from_barbara.sh` | rsync subsets of raw inputs from the shared filesystem into local `results/<dataset>/raw/` for laptop smoke tests. |
| `protect_barbara_data.sh` | Set the staged data tree read-only to prevent accidental writes from the workflow. |

## amet workflow entrypoints. Use these to run simulations + the three real
## datasets (Argelaguet, CRC, Ecker) whole-genome on a server with enough RAM.
##
## Not intended for laptops: the whole-genome runs allocate hundreds of GB of
## virtual memory across many parallel amet jobs. The recipes set ulimit -v
## 200 GB as a soft safeguard and let snakemake fan out across CORES cores.
##
## Usage:
##   make argelaguet   # whole-genome Argelaguet (4 Rmds)
##   make crc          # whole-genome CRC (6 Rmds)
##   make ecker        # whole-genome Ecker (4 Rmds)
##   make simulations  # simulations report
##   make all          # everything above
##   make dryrun       # snakemake -n for everything
##   make unlock       # release a stale snakemake lock
##
## Tunable variables (override on the command line):
##   CORES        snakemake --cores value (default 16)
##   ULIMIT_KB    virtual memory cap in KB (default 209715200, i.e. 200 GB)
##   CONDA_ENV    name of the conda env that holds snakemake (default snakemake)
##   CONDA_INIT   path to the conda activation script
##                (default ~/miniconda3/bin/activate)

CORES       ?= 16
ULIMIT_KB   ?= 209715200
CONDA_ENV   ?= snakemake
CONDA_INIT  ?= $(HOME)/miniconda3/bin/activate

WORKFLOW_DIR := workflow

## Standard preamble: activate the snakemake conda env and apply the
## virtual-memory ulimit. snakemake's per-job shells inherit the ulimit, so
## individual amet jobs are bounded by the same cap.
ACTIVATE := source $(CONDA_INIT) && conda activate $(CONDA_ENV) && \
            ulimit -v $(ULIMIT_KB)

SNAKEMAKE := snakemake --use-conda --cores $(CORES) -p

.PHONY: all simulations argelaguet crc ecker dryrun unlock clean help \
        setup-barbara

## Set up symlinks from results/{dataset}/ to a pre-existing data tree on
## barbara so amet does not re-download or re-rsync anything. Run once before
## `make all`. See workflow/scripts/internal/setup_barbara_links.sh for the
## env vars it honors.
setup-barbara:
	bash $(WORKFLOW_DIR)/scripts/internal/setup_barbara_links.sh

all: simulations argelaguet crc ecker

simulations:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && $(SNAKEMAKE) simulations'

argelaguet:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && $(SNAKEMAKE) argelaguet'

crc:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && $(SNAKEMAKE) crc'

ecker:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && $(SNAKEMAKE) ecker'

dryrun:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && \
	  snakemake --cores $(CORES) -n simulations argelaguet crc ecker'

unlock:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && snakemake --unlock'

help:
	@echo "Targets: all simulations argelaguet crc ecker dryrun unlock"
	@echo "Variables: CORES=$(CORES) ULIMIT_KB=$(ULIMIT_KB) CONDA_ENV=$(CONDA_ENV)"

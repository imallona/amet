## amet workflow entrypoints. Use these to run simulations + the three real
## datasets (Argelaguet, CRC, Ecker) on a server with enough RAM.
##
## Not intended for laptops: the runs allocate hundreds of GB of virtual
## memory across many parallel amet jobs. The recipes set ulimit -v
## 200 GB as a soft safeguard and let snakemake fan out across CORES cores.
##
## Usage:
##   make argelaguet                  # proto by default (results/argelaguet_proto/)
##   make crc MODE=full               # full grid (results/crc_full/)
##   make ecker MODE=proto            # explicit proto
##   make simulations                 # simulations report (MODE-agnostic)
##   make all MODE=full               # simulations + 3 datasets in full mode
##   make dryrun MODE=full            # snakemake -n for everything (full)
##   make unlock                      # release a stale snakemake lock
##
## Variables (override on the command line):
##   MODE         proto | full        which dataset config file to load
##                                    (default: proto)
##   CORES        snakemake --cores value (default 16)
##   ULIMIT_KB    virtual memory cap in KB (default 209715200, i.e. 200 GB)
##   CONDA_ENV    name of the conda env that holds snakemake (default snakemake)
##   CONDA_INIT   path to the conda activation script
##                (default ~/miniconda3/bin/activate)

MODE        ?= proto
CORES       ?= 16
ULIMIT_KB   ?= 209715200
CONDA_ENV   ?= snakemake
CONDA_INIT  ?= $(HOME)/miniconda3/bin/activate

WORKFLOW_DIR := workflow
DATASETS_CONFIG := config/datasets_$(MODE).yaml

## Standard preamble: activate the snakemake conda env and apply the
## virtual-memory ulimit. snakemake's per-job shells inherit the ulimit, so
## individual amet jobs are bounded by the same cap.
ACTIVATE := source $(CONDA_INIT) && conda activate $(CONDA_ENV) && \
            ulimit -v $(ULIMIT_KB)

## The trailing `--` stops snakemake's option parsing so subsequent positional
## tokens are unambiguously targets, not extra --configfile values.
SNAKEMAKE := snakemake --use-conda --cores $(CORES) -p \
             --configfile $(DATASETS_CONFIG) --

.PHONY: all simulations argelaguet crc ecker dryrun unlock clean help \
        setup-barbara

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
	  snakemake --cores $(CORES) --configfile $(DATASETS_CONFIG) \
	  -n -- simulations argelaguet crc ecker'

unlock:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && snakemake --unlock'

help:
	@echo "Targets: all simulations argelaguet crc ecker dryrun unlock setup-barbara"
	@echo "Variables: MODE=$(MODE) CORES=$(CORES) ULIMIT_KB=$(ULIMIT_KB) CONDA_ENV=$(CONDA_ENV)"
	@echo "Selected dataset config: $(DATASETS_CONFIG)"

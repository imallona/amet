#!/usr/bin/env bash
## Soft write-protect yamet's read-only data tree on barbara (or anywhere
## else amet runs against an existing yamet checkout). amet only ever reads
## these paths; chmod -R a-w them so a rogue rm/cp/cat-redirect cannot
## clobber yamet's files.
##
## This is a soft barrier (you own the files, you can chmod +w them back),
## not a sudo-level lock. Without sudo we can't use mount --bind -o ro or
## chattr +i. But it does block:
##   - rm <file>                        (refuses on read-only files)
##   - rm -f <file>                     (succeeds, but ONLY because rm only
##                                       needs write on the *parent dir*; -R
##                                       on the parent dir blocks this too)
##   - rm -rf parent/                   (fails because parent has no write
##                                       permission, so its children can't
##                                       be unlinked)
##   - cat foo > <file>                 (refused, file not writable)
##   - any tool's open(O_WRONLY)        (EACCES)
##
## It does not block someone running `chmod +w <path> && rm -rf <path>` by
## hand. If you need stronger protection, ask a sysadmin for `chattr +i`
## (needs root) or move yamet's data onto a read-only NFS export.
##
## Usage:
##   bash workflow/scripts/protect_barbara_data.sh           # protect
##   bash workflow/scripts/protect_barbara_data.sh --revert  # restore +w
##
## Override the yamet root via env: YAMET=/some/path bash ...

set -euo pipefail

YAMET="${YAMET:-$HOME/src/yamet/workflow}"
MODE="protect"
if [[ "${1:-}" == "--revert" ]]; then
    MODE="revert"
fi

paths=(
    "$YAMET/argelaguet/met/cpg_level"
    "$YAMET/argelaguet/sample_metadata.txt"
    "$YAMET/argelaguet/features/genomic_contexts"
    "$YAMET/mm10"
    "$YAMET/data/crc/raw"
    "$YAMET/hg19"
    "$YAMET/data/brain/raw"
)

apply() {
    local p="$1"
    if [[ ! -e "$p" ]]; then
        echo "[protect] skip: $p does not exist"
        return
    fi
    if [[ "$MODE" == "protect" ]]; then
        chmod -R a-w "$p"
        echo "[protect] chmod -R a-w $p"
    else
        chmod -R u+w "$p"
        echo "[protect] chmod -R u+w $p (reverted)"
    fi
}

echo "[protect] yamet root: $YAMET   mode: $MODE"
for p in "${paths[@]}"; do
    apply "$p"
done
echo "[protect] done."

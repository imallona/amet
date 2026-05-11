"""Subset cells.tsv to one (stage, lineage) sanitized pair.

Filter by sanitized stage and lineage, optionally cap at MAX_CELLS per group
with plate-stratified round-robin top-up by file size.

Usage:
    python argelaguet_subset_manifest.py \
        --cells cells.tsv \
        --stage E5-5 --lineage Epiblast \
        --max-cells 20 \
        --out manifests/E5-5_Epiblast.tsv
"""

import argparse
import csv
import os
import re

ap = argparse.ArgumentParser()
ap.add_argument("--cells", required=True)
ap.add_argument("--stage", required=True)
ap.add_argument("--lineage", required=True)
ap.add_argument("--max-cells", type=int, default=20)
ap.add_argument("--out", required=True)
args = ap.parse_args()


def sanitize(x):
    return re.sub(r"[ ._]", "-", str(x))


with open(args.cells) as f:
    reader = csv.DictReader(f, delimiter="\t")
    fieldnames = reader.fieldnames
    rows = list(reader)

stage_col = "stage"
lineage_col = "lineage10x"
plate_col = "plate" if "plate" in fieldnames else None

sub = [
    r for r in rows
    if sanitize(r.get(stage_col)) == args.stage
    and sanitize(r.get(lineage_col)) == args.lineage
]

if sub:
    for r in sub:
        try:
            r["_size"] = os.path.getsize(r["path"])
        except OSError:
            r["_size"] = 0
    if len(sub) > args.max_cells:
        if plate_col:
            plate_groups = {}
            for r in sub:
                p = r.get(plate_col) or "_unknown"
                plate_groups.setdefault(p, []).append(r)
            for p in plate_groups:
                plate_groups[p].sort(key=lambda r: r["_size"], reverse=True)
            picked = []
            order = sorted(plate_groups.keys())
            while len(picked) < args.max_cells and any(plate_groups.values()):
                for p in order:
                    if not plate_groups[p]:
                        continue
                    picked.append(plate_groups[p].pop(0))
                    if len(picked) == args.max_cells:
                        break
            sub = picked
        else:
            sub.sort(key=lambda r: r["_size"], reverse=True)
            sub = sub[: args.max_cells]
    for r in sub:
        r.pop("_size", None)

os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
with open(args.out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    for r in sub:
        w.writerow(r)

print(f"[subset_manifest] {args.stage}/{args.lineage}: {len(sub)} cells -> {args.out}")

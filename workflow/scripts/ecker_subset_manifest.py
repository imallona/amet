"""Subset cells.tsv to one (sub_region, sub_type) combo for Ecker.

Filter to the sanitized (sub_region, sub_type) pair, then keep the top
--max-cells cells ranked by the `size` column (coverage proxy: source TAR
size, written by the manifest builder) with plate-stratified round-robin:
within each plate cells are ranked by size; across plates the max-cells
slots are distributed evenly so a single high-coverage plate cannot
dominate the pick.

Sanitization: replace every space with '-' to match the wildcard values
the smk rules use (e.g. "IT-L23 Cux1" -> "IT-L23-Cux1").

Usage:
    python ecker_subset_manifest.py \
        --cells cells.tsv \
        --sub-region MOp --sub-type "IT-L23-Cux1" \
        --max-cells 50 \
        --out manifests/MOp_IT-L23-Cux1.tsv
"""

import argparse
import csv
import os

ap = argparse.ArgumentParser()
ap.add_argument("--cells", required=True)
ap.add_argument("--sub-region", required=True)
ap.add_argument("--sub-type", required=True)
ap.add_argument("--max-cells", type=int, default=50)
ap.add_argument("--out", required=True)
args = ap.parse_args()


def sanitize(x):
    return str(x).replace(" ", "-")


with open(args.cells) as f:
    reader = csv.DictReader(f, delimiter="\t")
    fieldnames = reader.fieldnames
    rows = list(reader)

sub = [r for r in rows
       if sanitize(r.get("sub_region", "")) == args.sub_region
       and sanitize(r.get("sub_type", "")) == args.sub_type]


def cell_size(row):
    val = row.get("size")
    if val not in (None, ""):
        try:
            return int(val)
        except ValueError:
            pass
    try:
        return os.path.getsize(row["path"])
    except OSError:
        return 0


if len(sub) > args.max_cells:
    plate_col = "plate" if "plate" in fieldnames else None
    if plate_col:
        plate_groups = {}
        for r in sub:
            p = r.get(plate_col) or "_unknown"
            plate_groups.setdefault(p, []).append(r)
        for p in plate_groups:
            plate_groups[p].sort(key=cell_size, reverse=True)
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
        sub.sort(key=cell_size, reverse=True)
        sub = sub[: args.max_cells]
else:
    sub.sort(key=cell_size, reverse=True)

os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
with open(args.out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    for r in sub:
        w.writerow(r)

print(f"[ecker_subset] {args.sub_region}/{args.sub_type}: {len(sub)} cells -> {args.out}")

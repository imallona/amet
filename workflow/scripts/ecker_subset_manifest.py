"""Subset cells.tsv to one (sub_region, sub_type) combo. Mirrors yamet's
get_ecker_harmonized_files: filter by sanitized sub_region and sub_type, cap
at MAX_CELLS by file size.

Sanitization: replace every space with '-' (yamet's _sanitize). Stage
strings have spaces in yamet's metadata (e.g. "IT-L23 Cux1") so the wildcard
value is "IT-L23-Cux1".

Usage:
    python ecker_subset_manifest.py \
        --cells cells.tsv \
        --sub-region MOp --sub-type "IT-L23-Cux1" \
        --max-cells 20 \
        --out manifests/MOp_IT-L23-Cux1.tsv
"""

import argparse
import csv
import os

ap = argparse.ArgumentParser()
ap.add_argument("--cells", required=True)
ap.add_argument("--sub-region", required=True)
ap.add_argument("--sub-type", required=True)
ap.add_argument("--max-cells", type=int, default=20)
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

if len(sub) > args.max_cells:
    for r in sub:
        try:
            r["_size"] = os.path.getsize(r["path"])
        except OSError:
            r["_size"] = 0
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

print(f"[ecker_subset] {args.sub_region}/{args.sub_type}: {len(sub)} cells -> {args.out}")

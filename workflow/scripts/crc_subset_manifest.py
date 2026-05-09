"""Subset cells.tsv to one (patient, location) combo. Mirrors yamet's
get_harmonized_files for CRC: filter by exact patient + location, optionally
cap at MAX_CELLS by file size.

Usage:
    python crc_subset_manifest.py \
        --cells cells.tsv \
        --patient CRC01 --location NC \
        --max-cells 20 \
        --out manifests/CRC01_NC.tsv
"""

import argparse
import csv
import os

ap = argparse.ArgumentParser()
ap.add_argument("--cells", required=True)
ap.add_argument("--patient", required=True)
ap.add_argument("--location", required=True)
ap.add_argument("--max-cells", type=int, default=20)
ap.add_argument("--out", required=True)
args = ap.parse_args()

with open(args.cells) as f:
    reader = csv.DictReader(f, delimiter="\t")
    fieldnames = reader.fieldnames
    rows = list(reader)

sub = [r for r in rows
       if r.get("patient") == args.patient
       and r.get("location") == args.location]

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

print(f"[crc_subset] {args.patient}/{args.location}: {len(sub)} cells -> {args.out}")

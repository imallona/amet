"""Subset cells.tsv to one (patient, location) combo for CRC.

Filter to exact patient and location, then keep the top --max-cells cells
ranked by the `size` column (coverage proxy: singleC.txt.gz size on disk is
monotonic in observed CpGs). CRC has no per-cell plate metadata, so picks
are plain top-N -- no plate stratification.

Usage:
    python crc_subset_manifest.py \
        --cells cells.tsv \
        --patient CRC01 --location NC \
        --max-cells 50 \
        --out manifests/CRC01_NC.tsv
"""

import argparse
import csv
import os

ap = argparse.ArgumentParser()
ap.add_argument("--cells", required=True)
ap.add_argument("--patient", required=True)
ap.add_argument("--location", required=True)
ap.add_argument("--max-cells", type=int, default=50)
ap.add_argument("--out", required=True)
args = ap.parse_args()

with open(args.cells) as f:
    reader = csv.DictReader(f, delimiter="\t")
    fieldnames = reader.fieldnames
    rows = list(reader)

sub = [r for r in rows
       if r.get("patient") == args.patient
       and r.get("location") == args.location]


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


sub.sort(key=cell_size, reverse=True)
if len(sub) > args.max_cells:
    sub = sub[: args.max_cells]

os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
with open(args.out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    for r in sub:
        w.writerow(r)

print(f"[crc_subset] {args.patient}/{args.location}: {len(sub)} cells -> {args.out}")

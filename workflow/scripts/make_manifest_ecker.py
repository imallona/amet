"""Build cells.tsv for Ecker from meta.tsv.gz + the tar files in raw_dir.

Filter to a region (default MOp) and (in prototype mode) to a set of
cell types, capping at cells_per_group cells per group. The manifest emits
one row per cell with `path` pointing to <cells_dir>/<cell_id>.tsv.gz, which
ecker_extract_tar will materialise on demand.
"""

import argparse
import os
import os.path as op
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--meta", required=True)
ap.add_argument("--raw_dir", required=True)
ap.add_argument("--cells_dir", required=True)
ap.add_argument("--region", default="MOp")
ap.add_argument("--proto_cell_types", default="")
ap.add_argument("--cells_per_group", type=int, default=10)
ap.add_argument("--group_col", default="sub_type")
ap.add_argument("--prototype", default="true")
ap.add_argument("--out", required=True)
args = ap.parse_args()

prototype = args.prototype.lower() in ("true", "1", "yes")
proto_cell_types = [s for s in args.proto_cell_types.split(",") if s]

meta = pd.read_csv(args.meta, sep="\t", compression="gzip")
print(f"[manifest] meta rows: {len(meta)}")

# Region filter (e.g. MOp). Match against any column whose name suggests
# region; if none of the columns exist we keep everything.
region_cols = [c for c in ("region", "sub_region") if c in meta.columns]
if region_cols and args.region:
    mask = pd.Series(False, index=meta.index)
    for c in region_cols:
        mask |= meta[c].astype(str).str.contains(args.region, na=False)
    meta = meta[mask]
    print(f"[manifest] after region={args.region}: {len(meta)}")

# Restrict to cells whose tar is present locally.
have_tar = []
for _, row in meta.iterrows():
    bn = str(row.get("basename", ""))
    if bn and op.exists(op.join(args.raw_dir, bn)):
        have_tar.append(True)
    else:
        have_tar.append(False)
meta = meta[have_tar]
print(f"[manifest] after presence-on-disk: {len(meta)}")

if prototype and proto_cell_types and args.group_col in meta.columns:
    meta = meta[meta[args.group_col].isin(proto_cell_types)]
    print(f"[manifest] after proto_cell_types: {len(meta)}")
    meta = meta.groupby(args.group_col, group_keys=False).head(args.cells_per_group)
    print(f"[manifest] after cells_per_group cap: {len(meta)}")

if args.group_col not in meta.columns:
    raise SystemExit(f"group column '{args.group_col}' missing from meta")

os.makedirs(args.cells_dir, exist_ok=True)


def cell_path(basename):
    cell_id = basename.replace(".tsv.tar", "")
    return op.join(args.cells_dir, cell_id + ".tsv.gz")


out = pd.DataFrame({
    "cell_id": meta["basename"].str.replace(".tsv.tar", "", regex=False),
    "group": meta[args.group_col].astype(str),
    "path": meta["basename"].apply(cell_path).apply(op.abspath),
    "format": "allc",
})
for c in ("cell_class", "major_type", "region", "sub_region", "plate"):
    if c in meta.columns:
        out[c] = meta[c].astype(str).values

out.to_csv(args.out, sep="\t", index=False)
print(f"[manifest] wrote {len(out)} cells across {out['group'].nunique()} groups -> {args.out}")

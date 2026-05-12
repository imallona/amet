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

# Region filter (e.g. MOp). Exact-match against sub_region (the cortical
# area) to avoid mis-grabbing cells from MOpUL/MOp2/etc. Atlas-dissection
# slabs (region = 2C/3C/4B/5D) are not used for filtering.
if "sub_region" in meta.columns and args.region:
    meta = meta[meta["sub_region"].astype(str) == args.region]
    print(f"[manifest] after sub_region == {args.region}: {len(meta)}")

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
    # No head() here. The per-combo subset picks top-N by source-TAR size
    # (proxy for cell coverage) and plate-balances the pick.

# Per-cell coverage proxy: the source TAR contains one allc tsv.gz per cell,
# so TAR size is monotonic in observed CpGs. Computed at manifest time
# because per-cell extracted tsv.gz files don't exist yet.
def tar_size(basename):
    try:
        return os.path.getsize(op.join(args.raw_dir, basename))
    except OSError:
        return 0

meta = meta.copy()
meta["size"] = meta["basename"].apply(tar_size)

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
    "size": meta["size"].values,
})
## sub_type and sub_region must both be available as separate keys (the
## wildcard names {sub_region, sub_type} read these).
for c in ("sub_type", "cell_class", "major_type", "region", "sub_region", "plate"):
    if c in meta.columns:
        out[c] = meta[c].astype(str).values

out.to_csv(args.out, sep="\t", index=False)
print(f"[manifest] wrote {len(out)} cells across {out['group'].nunique()} groups -> {args.out}")

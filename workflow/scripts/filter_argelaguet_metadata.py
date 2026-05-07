"""Drop QC-failing cells, drop TET-TKO plates, keep stage / lineage / embryo."""

import argparse
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--meta_in", required=True)
ap.add_argument("--meta_out", required=True)
args = ap.parse_args()

tet_ko_plates = [f"Plate{i}" for i in range(11, 17)]

meta = pd.read_csv(args.meta_in, sep="\t")
n_in = len(meta)
print(f"[filter_metadata] input rows: {n_in}")

meta = meta[meta["pass_metQC"] == True]
print(f"[filter_metadata] after pass_metQC: {len(meta)}")

meta = meta[~meta["plate"].isin(tet_ko_plates)]
print(f"[filter_metadata] after dropping TET-TKO plates: {len(meta)}")

meta = meta.dropna(subset=["id_met", "stage", "lineage10x"])
print(f"[filter_metadata] after dropping NA on required cols: {len(meta)}")

meta.to_csv(args.meta_out, sep="\t", index=False, compression="gzip")
print(f"[filter_metadata] wrote {args.meta_out}")

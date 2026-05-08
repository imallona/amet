"""Compose Ecker single-cell metadata.

Inputs:
  --nemo  : NeMo MOp_Metadata.tsv.gz (per-cell sub_type / cell_class / region).
  --paper : Paper supplementary xlsx with major_type / sub_region / coverage.

Output:
  TSV with columns: cell_id, basename, sub_type, cell_class, major_type,
  region, sub_region, plate. Rows with missing sub_type are dropped.
"""

import argparse
import pandas as pd
import re

ap = argparse.ArgumentParser()
ap.add_argument("--nemo", required=True)
ap.add_argument("--paper", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()

nemo = pd.read_csv(args.nemo, sep="\t", compression="gzip")
print(f"[compose] nemo rows: {len(nemo)}")

# Cell id is in 'CellID' or first column; basename comes from the AllcPath
# tail after stripping .tsv.gz.
cell_col = "CellID" if "CellID" in nemo.columns else nemo.columns[0]
nemo = nemo.rename(columns={cell_col: "cell_id"})
nemo["basename"] = (
    nemo["AllcPath"].fillna("").apply(lambda p: p.rsplit("/", 1)[-1])
).str.replace(r"\.tsv\.gz$", ".tsv.tar", regex=True)


def _read_paper(path):
    # The supplementary xlsx has multiple sheets; we want the cell-level table.
    sheets = pd.read_excel(path, sheet_name=None, engine="openpyxl")
    # Pick the largest sheet with a CellID-like column.
    best = None
    for name, df in sheets.items():
        if not len(df):
            continue
        cands = [c for c in df.columns if re.search("cell.?id|sample", str(c), re.I)]
        if not cands:
            continue
        if best is None or len(df) > len(best[1]):
            best = (cands[0], df)
    if best is None:
        return pd.DataFrame()
    cell_col, df = best
    df = df.rename(columns={cell_col: "cell_id"})
    return df


paper = _read_paper(args.paper)
print(f"[compose] paper rows: {len(paper)}")

merged = nemo.merge(paper, on="cell_id", how="left", suffixes=("", "_paper"))

# Pick group columns. Names vary across the paper sheets so we map a few
# common variants to canonical column names without altering values.
column_aliases = {
    "sub_type": ["SubType", "sub_type", "Subtype"],
    "cell_class": ["CellClass", "cell_class", "MajorType_Cluster"],
    "major_type": ["MajorType", "major_type", "MajorClass"],
    "region": ["Region", "region", "MajorRegion"],
    "sub_region": ["SubRegion", "sub_region"],
    "plate": ["Plate", "plate"],
}
for canonical, aliases in column_aliases.items():
    for alias in aliases:
        if alias in merged.columns and canonical != alias:
            merged[canonical] = merged.get(canonical).fillna(merged[alias]) if canonical in merged.columns else merged[alias]
            break

keep = ["cell_id", "basename"] + [c for c in column_aliases if c in merged.columns]
merged = merged[keep]
merged = merged.dropna(subset=["sub_type"]) if "sub_type" in merged.columns else merged
print(f"[compose] after dropna sub_type: {len(merged)}")

merged.to_csv(args.out, sep="\t", index=False, compression="gzip")
print(f"[compose] wrote {args.out}")

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
    """Parse the Ecker supplementary xlsx.

    Sheet1 starts with a column-key dictionary (rows 0..N) and then the actual
    cell-level table whose first column holds the cell id but is unnamed.
    Locate the header row by scanning for one whose non-NA cell count jumps
    relative to the dictionary preamble; treat the first column of the
    sub-table as cell_id.
    """
    raw = pd.read_excel(path, sheet_name=0, engine="openpyxl",
                        header=None, nrows=200)
    counts = raw.notna().sum(axis=1)
    header_row = None
    for i, c in enumerate(counts):
        if c >= 6 and (i == 0 or counts.iloc[i - 1] <= 2):
            header_row = i
            break
    if header_row is None:
        return pd.DataFrame()
    df = pd.read_excel(path, sheet_name=0, engine="openpyxl",
                       header=header_row)
    first_col = df.columns[0]
    if pd.isna(first_col) or str(first_col).startswith("Unnamed"):
        df = df.rename(columns={first_col: "cell_id"})
    else:
        df.insert(0, "cell_id", df[first_col])
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

"""
Reshape 10x Dcode binarized matrix files from wide to long format.

Input:  wide CSV — one row per cell, epitopes as columns
Output: long CSV — one row per CDR3b + epitope pair (binders only)

Output schema matches VDJdb/McPAS for uniform loading:
  barcode, donor, cdr3b, cdr3a, epitope, HLA, antigen, pathogen,
  dextramer_count, cell_phenotype
"""

import pandas as pd
import re
from pathlib import Path

DONOR_FILES = sorted(
    Path("Databases/10xDcode").glob("vdj_v1_hs_aggregated_donor*_binarized_matrix.csv")
)
OUTPUT = Path("Databases/10xDcode/10xDcode_long.csv")

BINDER_RE = re.compile(r"^(.+)_binder$")
EPITOPE_RE = re.compile(r"^([^_]+)_([^_]+)_(.+?)_([^_]+)$")


def parse_epitope_col(col: str) -> dict | None:
    """Parse 'A0201_GILGFVFTL_Flu-MP_Influenza' into components."""
    m = EPITOPE_RE.match(col)
    if not m:
        return None
    return {
        "HLA": m.group(1),
        "epitope": m.group(2),
        "antigen": m.group(3),
        "pathogen": m.group(4),
    }


def extract_chains(cdr3_str: str) -> tuple[str, str]:
    """Return (cdr3b, cdr3a) from 'TRA:seq;TRB:seq' string."""
    cdr3b, cdr3a = "", ""
    if pd.isna(cdr3_str):
        return cdr3b, cdr3a
    for part in str(cdr3_str).split(";"):
        chain, _, seq = part.partition(":")
        if chain == "TRB" and not cdr3b:
            cdr3b = seq
        elif chain == "TRA" and not cdr3a:
            cdr3a = seq
    return cdr3b, cdr3a


def cell_phenotype(row: pd.Series) -> str:
    cd4 = row.get("CD4", 0) or 0
    cd8 = row.get("CD8a", 0) or 0
    if cd4 > cd8:
        return "CD4"
    if cd8 > cd4:
        return "CD8"
    return "DN"


def reshape_donor(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)

    # identify binder columns and their base epitope column name
    binder_cols = [c for c in df.columns if BINDER_RE.match(c)]
    epitope_base = {c: BINDER_RE.match(c).group(1) for c in binder_cols}

    # parse epitope metadata once
    epitope_meta = {
        base: parse_epitope_col(base)
        for base in epitope_base.values()
    }
    # skip columns that don't match the expected format (e.g. negative controls)
    valid_binder_cols = [
        c for c in binder_cols if epitope_meta[epitope_base[c]] is not None
    ]

    rows = []
    for _, cell in df.iterrows():
        cdr3b, cdr3a = extract_chains(cell.get("cell_clono_cdr3_aa", ""))
        if not cdr3b:
            continue
        phenotype = cell_phenotype(cell)
        for bcol in valid_binder_cols:
            if not cell[bcol]:
                continue
            base = epitope_base[bcol]
            meta = epitope_meta[base]
            rows.append({
                "barcode": cell["barcode"],
                "donor": cell["donor"],
                "cdr3b": cdr3b,
                "cdr3a": cdr3a,
                "epitope": meta["epitope"],
                "HLA": meta["HLA"],
                "antigen": meta["antigen"],
                "pathogen": meta["pathogen"],
                "dextramer_count": cell.get(base, 0),
                "cell_phenotype": phenotype,
            })

    return pd.DataFrame(rows)


def main():
    all_donors = []
    for path in DONOR_FILES:
        print(f"Processing {path.name}...")
        donor_df = reshape_donor(path)
        print(f"  {len(donor_df):,} binder rows")
        all_donors.append(donor_df)

    combined = pd.concat(all_donors, ignore_index=True)
    combined.drop_duplicates(subset=["cdr3b", "epitope", "HLA", "donor"], inplace=True)
    combined.to_csv(OUTPUT, index=False)
    print(f"\nSaved {len(combined):,} rows to {OUTPUT}")
    print(combined.head())


if __name__ == "__main__":
    main()

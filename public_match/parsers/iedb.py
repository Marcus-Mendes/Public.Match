import pandas as pd
from pathlib import Path

IEDB_PATH = Path("Databases/IEDB/iedb.xlsx")


def _coalesce_cdr3(curated: pd.Series, calculated: pd.Series) -> pd.Series:
    """Use curated CDR3 when available, fall back to calculated."""
    return curated.where(curated.notna() & (curated != ""), calculated)


def load(path: Path = IEDB_PATH) -> pd.DataFrame:
    df = pd.read_excel(path, dtype_backend="numpy_nullable")
    df = df[df["Receptor - Type"] == "alphabeta"].copy()

    c1_beta  = df["Chain 1 - Type"] == "beta"
    c1_alpha = df["Chain 1 - Type"] == "alpha"
    c2_beta  = df["Chain 2 - Type"] == "beta"
    c2_alpha = df["Chain 2 - Type"] == "alpha"

    def _cdr3(sub, curated_col, calc_col):
        return _coalesce_cdr3(sub[curated_col].astype(str), sub[calc_col].astype(str)).str.upper().str.strip()

    def _meta(sub):
        return {
            "epitope":   sub["Epitope - Name"].astype(str).str.strip(),
            "antigen":   sub["Epitope - Source Molecule"].astype(str).str.strip(),
            "pathogen":  sub["Epitope - Source Organism"].astype(str).str.strip(),
            "HLA":       sub["Assay - MHC Allele Names"].astype(str).str.strip(),
            "source_db": "IEDB",
        }

    parts = []

    # Chain 1 = beta; alpha partner from Chain 2 where available
    if c1_beta.any():
        sub = df[c1_beta]
        cdr3a = _cdr3(sub, "Chain 2 - CDR3 Curated", "Chain 2 - CDR3 Calculated")
        parts.append(pd.DataFrame({
            "cdr3b": _cdr3(sub, "Chain 1 - CDR3 Curated", "Chain 1 - CDR3 Calculated"),
            "cdr3a": cdr3a.where(c2_alpha[c1_beta].values),
            **_meta(sub),
        }))

    # Chain 2 = beta (and Chain 1 is not beta, to avoid double-counting)
    mask2 = c2_beta & ~c1_beta
    if mask2.any():
        sub = df[mask2]
        cdr3a = _cdr3(sub, "Chain 1 - CDR3 Curated", "Chain 1 - CDR3 Calculated")
        parts.append(pd.DataFrame({
            "cdr3b": _cdr3(sub, "Chain 2 - CDR3 Curated", "Chain 2 - CDR3 Calculated"),
            "cdr3a": cdr3a.where(c1_alpha[mask2].values),
            **_meta(sub),
        }))

    if not parts:
        return pd.DataFrame(columns=["cdr3b", "cdr3a", "epitope", "antigen", "pathogen", "HLA", "source_db"])

    out = pd.concat(parts, ignore_index=True)
    out = out[out["cdr3b"].str.match(r"^[ACDEFGHIKLMNPQRSTVWY]+$", na=False)]
    return out.drop_duplicates(subset=["cdr3b", "epitope"]).reset_index(drop=True)

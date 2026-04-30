import pandas as pd
from pathlib import Path

IEDB_PATH = Path("Databases/IEDB/iedb.xlsx")


def _coalesce_cdr3(curated: pd.Series, calculated: pd.Series) -> pd.Series:
    """Use curated CDR3 when available, fall back to calculated."""
    return curated.where(curated.notna() & (curated != ""), calculated)


def load(path: Path = IEDB_PATH) -> pd.DataFrame:
    df = pd.read_excel(path, dtype_backend="numpy_nullable")

    # keep alpha-beta TCRs only
    df = df[df["Receptor - Type"] == "alphabeta"].copy()

    # rows where chain 1 is beta
    c1_beta = df["Chain 1 - Type"] == "beta"
    # rows where chain 2 is beta (paired receptors)
    c2_beta = df["Chain 2 - Type"] == "beta"

    def extract_beta(mask: pd.Series, curated_col: str, calc_col: str) -> pd.DataFrame:
        sub = df[mask].copy()
        cdr3b = _coalesce_cdr3(sub[curated_col].astype(str), sub[calc_col].astype(str))
        return pd.DataFrame({
            "cdr3b": cdr3b.str.upper().str.strip(),
            "epitope": sub["Epitope - Name"].astype(str).str.strip(),
            "antigen": sub["Epitope - Source Molecule"].astype(str).str.strip(),
            "pathogen": sub["Epitope - Source Organism"].astype(str).str.strip(),
            "HLA": sub["Assay - MHC Allele Names"].astype(str).str.strip(),
            "source_db": "IEDB",
        })

    parts = []
    if c1_beta.any():
        parts.append(extract_beta(c1_beta, "Chain 1 - CDR3 Curated", "Chain 1 - CDR3 Calculated"))
    if c2_beta.any():
        parts.append(extract_beta(c2_beta, "Chain 2 - CDR3 Curated", "Chain 2 - CDR3 Calculated"))

    if not parts:
        return pd.DataFrame(columns=["cdr3b", "epitope", "antigen", "pathogen", "HLA", "source_db"])

    out = pd.concat(parts, ignore_index=True)
    out = out[out["cdr3b"].str.match(r"^[ACDEFGHIKLMNPQRSTVWY]+$", na=False)]
    return out.drop_duplicates(subset=["cdr3b", "epitope"]).reset_index(drop=True)

import pandas as pd
from pathlib import Path

CEDAR_PATH = Path("Databases/CEDAR/cedar.xlsx")


def _coalesce_cdr3(curated: pd.Series, calculated: pd.Series) -> pd.Series:
    return curated.where(curated.notna() & (curated.astype(str) != ""), calculated)


def load(path: Path = CEDAR_PATH) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name="Sheet1")
    df = df[df["Receptor - Type"] == "alphabeta"].copy()

    c1_alpha = df["Chain 1 - Type"] == "alpha"
    c2_beta  = df["Chain 2 - Type"] == "beta"

    def _cdr3a(sub):
        return _coalesce_cdr3(
            sub["Chain 1 - CDR3 Curated"].astype(str),
            sub["Chain 1 - CDR3 Calculated"].astype(str),
        ).str.upper().str.strip()

    def _cdr3b(sub):
        return _coalesce_cdr3(
            sub["Chain 2 - CDR3 Curated"].astype(str),
            sub["Chain 2 - CDR3 Calculated"].astype(str),
        ).str.upper().str.strip()

    def _meta(sub):
        return {
            "epitope":   sub["Epitope - Name"].astype(str).str.strip(),
            "antigen":   sub["Epitope - Source Molecule"].astype(str).str.strip(),
            "pathogen":  sub["Epitope - Source Organism"].astype(str).str.strip(),
            "HLA":       sub["Assay - MHC Allele Names"].astype(str).str.strip(),
            "source_db": "CEDAR",
        }

    parts = []

    # --- Paired rows (Chain 1 = alpha AND Chain 2 = beta) ---
    paired_mask = c1_alpha & c2_beta
    if paired_mask.any():
        sub = df[paired_mask]
        parts.append(pd.DataFrame({
            "cdr3b": _cdr3b(sub),
            "cdr3a": _cdr3a(sub),
            **_meta(sub),
        }))

    # --- Beta-only rows (no alpha chain) ---
    beta_only = ~c1_alpha & c2_beta
    if beta_only.any():
        sub = df[beta_only]
        parts.append(pd.DataFrame({
            "cdr3b": _cdr3b(sub),
            "cdr3a": pd.NA,
            **_meta(sub),
        }))

    # --- Alpha-only rows (no beta chain) ---
    alpha_only = c1_alpha & ~c2_beta
    if alpha_only.any():
        sub = df[alpha_only]
        parts.append(pd.DataFrame({
            "cdr3b": pd.NA,
            "cdr3a": _cdr3a(sub),
            **_meta(sub),
        }))

    if not parts:
        return pd.DataFrame(columns=["cdr3b", "cdr3a", "epitope", "antigen", "pathogen", "HLA", "source_db"])

    out = pd.concat(parts, ignore_index=True)

    aa_pat = r"^[ACDEFGHIKLMNPQRSTVWY]+$"
    cdr3b_ok = out["cdr3b"].str.match(aa_pat, na=True)
    cdr3a_ok = out["cdr3a"].str.match(aa_pat, na=True)
    out = out[cdr3b_ok & cdr3a_ok]

    return out.drop_duplicates(subset=["cdr3b", "cdr3a", "epitope"]).reset_index(drop=True)

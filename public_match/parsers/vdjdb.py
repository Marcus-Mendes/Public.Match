import pandas as pd
from pathlib import Path

VDJDB_PATH = Path("Databases/VDJdb/vdjdb.slim.txt")

_AA_PAT = r"^[ACDEFGHIKLMNPQRSTVWY]+$"


def load(path: Path = VDJDB_PATH, min_score: int = 1) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)

    df = df[df["species"] == "HomoSapiens"].copy()
    df = df[df["vdjdb.score"].fillna(0).astype(int) >= min_score]

    tra = df[df["gene"] == "TRA"].copy()
    trb = df[df["gene"] == "TRB"].copy()

    # complex.ids where BOTH chains are present
    paired_ids = set(tra["complex.id"]) & set(trb["complex.id"])

    def _meta(sub: pd.DataFrame) -> dict:
        return {
            "epitope":   sub["antigen.epitope"].astype(str).str.strip(),
            "antigen":   sub["antigen.gene"].astype(str).str.strip(),
            "pathogen":  sub["antigen.species"].astype(str).str.strip(),
            "HLA":       sub["mhc.a"].astype(str).str.strip(),
            "source_db": "VDJdb",
        }

    parts = []

    # --- Paired: TRB rows joined with their TRA partner on complex.id ---
    trb_p = trb[trb["complex.id"].isin(paired_ids)].copy()
    tra_lookup = (
        tra[tra["complex.id"].isin(paired_ids)][["complex.id", "cdr3"]]
        .rename(columns={"cdr3": "cdr3a"})
        .drop_duplicates("complex.id")        # one alpha per complex
    )
    merged = trb_p.merge(tra_lookup, on="complex.id", how="left")
    parts.append(pd.DataFrame({
        "cdr3b": merged["cdr3"].str.upper().str.strip(),
        "cdr3a": merged["cdr3a"].str.upper().str.strip(),
        **_meta(merged),
    }))

    # --- Unpaired TRB rows (no TRA partner found) ---
    trb_u = trb[~trb["complex.id"].isin(paired_ids)]
    if len(trb_u):
        parts.append(pd.DataFrame({
            "cdr3b": trb_u["cdr3"].str.upper().str.strip(),
            "cdr3a": pd.NA,
            **_meta(trb_u),
        }))

    # --- Unpaired TRA rows (no TRB partner found) ---
    tra_u = tra[~tra["complex.id"].isin(paired_ids)]
    if len(tra_u):
        parts.append(pd.DataFrame({
            "cdr3b": pd.NA,
            "cdr3a": tra_u["cdr3"].str.upper().str.strip(),
            **_meta(tra_u),
        }))

    out = pd.concat(parts, ignore_index=True)

    # validate each chain where present; drop rows with invalid AA sequences
    cdr3b_ok = out["cdr3b"].str.match(_AA_PAT, na=True)
    cdr3a_ok = out["cdr3a"].str.match(_AA_PAT, na=True)
    out = out[cdr3b_ok & cdr3a_ok]

    return out.drop_duplicates(subset=["cdr3b", "cdr3a", "epitope"]).reset_index(drop=True)

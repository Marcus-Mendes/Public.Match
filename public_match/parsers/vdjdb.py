import pandas as pd
from pathlib import Path
import zipfile
import io

VDJDB_PATH = Path("Databases/VDJdb/VDJdb_05302026.zip")


def load(path: Path = VDJDB_PATH, min_score: int = 1) -> pd.DataFrame:
    with zipfile.ZipFile(path) as z:
        tsv_name = next(n for n in z.namelist() if n.endswith(".tsv"))
        with z.open(tsv_name) as f:
            df = pd.read_csv(f, sep="\t", low_memory=False)

    df = df[(df["Gene"] == "TRB") & (df["Species"] == "HomoSapiens")].copy()
    df = df[df["Score"].fillna(0).astype(int) >= min_score]

    out = pd.DataFrame({
        "cdr3b": df["CDR3"].str.upper().str.strip(),
        "epitope": df["Epitope"].str.strip(),
        "antigen": df["Epitope gene"].str.strip(),
        "pathogen": df["Epitope species"].str.strip(),
        "HLA": df["MHC A"].str.strip(),
        "source_db": "VDJdb",
    })

    out = out[out["cdr3b"].str.match(r"^[ACDEFGHIKLMNPQRSTVWY]+$", na=False)]
    return out.dropna(subset=["cdr3b", "epitope"]).drop_duplicates(subset=["cdr3b", "epitope"]).reset_index(drop=True)

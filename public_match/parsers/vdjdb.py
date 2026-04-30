import pandas as pd
from pathlib import Path
import zipfile
import io

VDJDB_DIR = Path("Databases/VDJdb")


def _find_zip(directory: Path) -> Path:
    zips = sorted(directory.glob("VDJdb_*.zip"))
    if not zips:
        raise FileNotFoundError(f"No VDJdb_*.zip found in {directory}")
    return zips[-1]


def load(path: Path = None, min_score: int = 1) -> pd.DataFrame:
    if path is None:
        path = _find_zip(VDJDB_DIR)
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

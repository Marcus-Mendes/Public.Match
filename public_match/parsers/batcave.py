import pandas as pd
from pathlib import Path

BATCAVE_MHCI_PATH = Path("Databases/BATCAVE/TCR_pMHCI_mutational_scan.csv")
BATCAVE_MHCII_PATH = Path("Databases/BATCAVE/TCR_pMHCII_mutational_scan.csv")

# keep only native epitope rows with meaningful activation
_MHCI_ACTIVITY_THRESHOLD = 20.0   # 0–100 scale
_MHCII_ACTIVITY_THRESHOLD = 0.2   # 0–1 scale


def _load_file(path: Path, activity_threshold: float) -> pd.DataFrame:
    df = pd.read_csv(path)

    # filter: human TCRs, native epitope only, positive binders
    df = df[df["tcr_source_organism"].str.lower() == "human"].copy()
    df = df[df["peptide"] == df["index_peptide"]].copy()
    df = df[df["peptide_activity"] >= activity_threshold].copy()
    df = df[df["cdr3b"].notna()].copy()

    return pd.DataFrame({
        "cdr3b":    df["cdr3b"].str.upper().str.strip(),
        "cdr3a":    df["cdr3a"].astype(str).str.upper().str.strip(),
        "epitope":  df["index_peptide"].str.strip(),
        "antigen":  df["peptide_type"].str.strip(),
        "pathogen": df["peptide_type"].str.strip(),
        "HLA":      df["mhc"].str.strip(),
        "source_db": "BATCAVE",
    })


def load(
    mhci_path: Path = BATCAVE_MHCI_PATH,
    mhcii_path: Path = BATCAVE_MHCII_PATH,
) -> pd.DataFrame:
    parts = []
    if mhci_path.exists():
        parts.append(_load_file(mhci_path, _MHCI_ACTIVITY_THRESHOLD))
    if mhcii_path.exists():
        parts.append(_load_file(mhcii_path, _MHCII_ACTIVITY_THRESHOLD))

    if not parts:
        return pd.DataFrame(columns=["cdr3b", "epitope", "antigen", "pathogen", "HLA", "source_db"])

    out = pd.concat(parts, ignore_index=True)
    out = out[out["cdr3b"].str.match(r"^[ACDEFGHIKLMNPQRSTVWY]+$", na=False)]
    return out.drop_duplicates(subset=["cdr3b", "epitope"]).reset_index(drop=True)

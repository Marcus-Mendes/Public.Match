import pandas as pd
from pathlib import Path

_CDR3B_ALIASES = ["cdr3b", "cdr3_beta", "cdr3_b", "junction_aa", "cdr3", "CDR3", "sequence", "seq"]

_COL_ALIASES = {
    "epitope":  ["epitope", "peptide", "antigen_peptide"],
    "antigen":  ["antigen", "antigen_gene", "gene", "antigen_protein"],
    "pathogen": ["pathogen", "organism", "species", "pathogen_species", "source_organism"],
    "HLA":      ["HLA", "mhc", "allele", "hla_allele", "MHC"],
}


def _find_col(columns, aliases):
    cols_lower = {c.lower(): c for c in columns}
    for alias in aliases:
        if alias in columns:
            return alias
        if alias.lower() in cols_lower:
            return cols_lower[alias.lower()]
    return None


def load(path: Path, cdr3_col: str = None, source_name: str = None) -> pd.DataFrame:
    path = Path(path)
    sep = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","
    df = pd.read_csv(path, sep=sep, low_memory=False)

    if source_name is None:
        source_name = path.stem

    # auto-detect first; fall back to the user-supplied hint only when needed
    cdr3b_col = _find_col(df.columns, _CDR3B_ALIASES)
    if cdr3b_col is None and cdr3_col:
        if cdr3_col not in df.columns:
            raise ValueError(
                f"Column '{cdr3_col}' not found in {path.name}. "
                f"Available columns: {list(df.columns)}"
            )
        cdr3b_col = cdr3_col
    if cdr3b_col is None:
        raise ValueError(
            f"Cannot find CDR3β column in {path.name}. "
            f"Tried: {_CDR3B_ALIASES}. "
            f"Use --custom-db-cdr3-col to specify the column name."
        )

    out = pd.DataFrame({
        "cdr3b":     df[cdr3b_col].astype(str).str.upper().str.strip(),
        "source_db": source_name,
    })

    for schema_col, aliases in _COL_ALIASES.items():
        found = _find_col(df.columns, aliases)
        out[schema_col] = df[found].astype(str).str.strip() if found else ""

    out = out[out["cdr3b"].str.match(r"^[ACDEFGHIKLMNPQRSTVWY]+$", na=False)]
    return out.drop_duplicates(subset=["cdr3b", "epitope"]).reset_index(drop=True)

"""
Public.Match CLI

Usage
-----
# FASTA input (original)
python -m public_match --input sequences.fasta --db iedb vdjdb mcpas tenx

# TSV/CSV input
python -m public_match --input repertoire.tsv --seq-col cdr3_beta --name-col barcode

# Custom database
python -m public_match --input sequences.fasta --custom-db my_db.csv
python -m public_match --input sequences.fasta --custom-db my_db.tsv --custom-db-cdr3-col junction_aa
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

from public_match.database import load_databases, ALL_DBS
from public_match.matcher import match

_FASTA_EXTENSIONS = {".fasta", ".fa", ".faa", ".fas"}
_CDR3B_ALIASES = ["cdr3b", "cdr3_beta", "cdr3_b", "junction_aa", "cdr3", "CDR3", "sequence"]
_NAME_ALIASES  = ["name", "id", "cell_id", "barcode", "clone_id", "sample_id"]


def _detect_format(path: Path) -> str:
    if path.suffix.lower() in _FASTA_EXTENSIONS:
        return "fasta"
    if path.suffix.lower() in (".tsv", ".csv", ".txt"):
        return "tabular"
    with open(path) as f:
        first_char = f.read(1)
    return "fasta" if first_char == ">" else "tabular"


def _find_col(columns, aliases):
    cols_lower = {c.lower(): c for c in columns}
    for alias in aliases:
        if alias in columns:
            return alias
        if alias.lower() in cols_lower:
            return cols_lower[alias.lower()]
    return None


def parse_fasta(path: Path) -> dict[str, str]:
    sequences = {}
    current_name = None
    current_seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name:
                    sequences[current_name] = "".join(current_seq).upper()
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_name:
        sequences[current_name] = "".join(current_seq).upper()
    return sequences


def parse_tabular(path: Path, seq_col: str = None, name_col: str = None) -> dict[str, str]:
    sep = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","
    df = pd.read_csv(path, sep=sep)

    # resolve CDR3β column
    if seq_col:
        if seq_col not in df.columns:
            sys.exit(f"Error: --seq-col '{seq_col}' not found. Available columns: {list(df.columns)}")
        cdr3_col = seq_col
    else:
        cdr3_col = _find_col(df.columns, _CDR3B_ALIASES)
        if cdr3_col is None:
            sys.exit(
                f"Error: cannot find CDR3β column in {path.name}. "
                f"Tried: {_CDR3B_ALIASES}. Use --seq-col to specify."
            )

    # resolve name/ID column
    if name_col:
        if name_col not in df.columns:
            sys.exit(f"Error: --name-col '{name_col}' not found. Available columns: {list(df.columns)}")
        id_col = name_col
    else:
        id_col = _find_col(df.columns, _NAME_ALIASES)

    sequences: dict[str, str] = {}
    for i, row in df.iterrows():
        seq = str(row[cdr3_col]).upper().strip()
        if not seq or seq == "NAN":
            continue
        name = str(row[id_col]) if id_col else f"seq_{i + 1}"
        sequences[name] = seq

    return sequences


def main():
    parser = argparse.ArgumentParser(
        prog="public-match",
        description="Match CDR3β sequences against public TCR databases.",
    )
    parser.add_argument("--input", "-i", required=True, type=Path,
                        help="Input file: FASTA (.fasta/.fa) or tabular (.tsv/.csv)")
    parser.add_argument("--output", "-o", default=Path("public_match_results.csv"), type=Path,
                        help="Output CSV file (default: public_match_results.csv)")
    parser.add_argument("--db", nargs="+", default=ALL_DBS, choices=ALL_DBS,
                        help=f"Databases to search (default: all). Choices: {ALL_DBS}")
    parser.add_argument("--method", default="blosum", choices=["exact", "blosum", "edit"],
                        help="Matching method (default: blosum)")
    parser.add_argument("--threshold", type=float, default=0.97,
                        help="Score threshold: min BLOSUM62 score (blosum) or max edit distance (edit). Default: 0.97")

    tsv_group = parser.add_argument_group("tabular input options (TSV/CSV)")
    tsv_group.add_argument("--seq-col", metavar="COL",
                           help="Column containing CDR3β sequences (auto-detected if omitted)")
    tsv_group.add_argument("--name-col", metavar="COL",
                           help="Column to use as sequence ID (auto-detected if omitted)")

    custom_group = parser.add_argument_group("custom database options")
    custom_group.add_argument("--custom-db", nargs="+", type=Path, metavar="PATH",
                              help="One or more custom database files (CSV or TSV)")
    custom_group.add_argument("--custom-db-cdr3-col", metavar="COL",
                              help="CDR3β column name in custom DB files (auto-detected if omitted)")

    args = parser.parse_args()

    if not args.input.exists():
        sys.exit(f"Error: input file not found: {args.input}")

    print("Public.Match — CDR3β sequence search")
    print(f"  Input:     {args.input}")
    print(f"  Databases: {args.db}")
    if args.custom_db:
        print(f"  Custom DB: {[str(p) for p in args.custom_db]}")
    print(f"  Method:    {args.method}  |  Threshold: {args.threshold}")
    print()

    fmt = _detect_format(args.input)
    if fmt == "fasta":
        sequences = parse_fasta(args.input)
    else:
        sequences = parse_tabular(args.input, seq_col=args.seq_col, name_col=args.name_col)

    if not sequences:
        sys.exit("Error: no sequences found in input file.")
    print(f"Loaded {len(sequences)} query sequence(s).\n")

    print("Loading reference databases...")
    reference = load_databases(
        dbs=args.db,
        custom_paths=args.custom_db,
        custom_cdr3_col=args.custom_db_cdr3_col,
    )
    print(f"  Total reference CDR3b entries: {len(reference):,}\n")

    print("Running matches...")
    results = match(
        queries=list(sequences.values()),
        reference=reference,
        method=args.method,
        threshold=args.threshold,
    )

    seq_to_name = {v: k for k, v in sequences.items()}
    results.insert(0, "query_name", results["query_cdr3b"].map(seq_to_name))

    results.to_csv(args.output, index=False)
    print(f"Done. {len(results)} match(es) written to {args.output}")


if __name__ == "__main__":
    main()

"""
Public.Match CLI

Usage
-----
# Beta chain (default)
python -m public_match --input sequences.fasta

# Alpha chain
python -m public_match --input alpha_seqs.fasta --chain alpha

# Paired (TSV with both columns)
python -m public_match --input repertoire.tsv --chain paired

# Paired (two FASTA files, matched by sequence name)
python -m public_match --input beta.fasta --input-alpha alpha.fasta --chain paired

# TSV/CSV input with explicit column names
python -m public_match --input repertoire.tsv --seq-col cdr3_beta --alpha-col cdr3_alpha --chain paired

# Custom database
python -m public_match --input sequences.fasta --custom-db my_db.csv
python -m public_match --input sequences.fasta --custom-db my_db.tsv --custom-db-cdr3-col junction_aa
"""

import argparse
import sys
from pathlib import Path
from typing import Literal

import pandas as pd

from public_match.database import load_databases, ALL_DBS
from public_match.matcher import match

_FASTA_EXTENSIONS = {".fasta", ".fa", ".faa", ".fas"}
_CDR3B_ALIASES = ["cdr3b", "cdr3_beta", "cdr3_b", "junction_aa", "cdr3", "CDR3", "sequence"]
_CDR3A_ALIASES = ["cdr3a", "cdr3_alpha", "cdr3_a", "TRA_CDR3", "cdr3_TRA", "junction_aa_alpha"]
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


def parse_tabular(
    path: Path,
    seq_col: str = None,
    name_col: str = None,
    alpha_col: str = None,
    chain: str = "beta",
) -> dict[str, str | tuple[str, str]]:
    sep = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","
    df = pd.read_csv(path, sep=sep)

    # resolve name/ID column
    if name_col:
        if name_col not in df.columns:
            sys.exit(f"Error: --name-col '{name_col}' not found. Available columns: {list(df.columns)}")
        id_col = name_col
    else:
        id_col = _find_col(df.columns, _NAME_ALIASES)

    def _resolve_col(override, aliases, label):
        if override:
            if override not in df.columns:
                sys.exit(f"Error: {label} column '{override}' not found. Available: {list(df.columns)}")
            return override
        found = _find_col(df.columns, aliases)
        if found is None:
            sys.exit(
                f"Error: cannot find {label} column in {path.name}. "
                f"Tried: {aliases}. Use the appropriate --*-col flag to specify."
            )
        return found

    if chain == "alpha":
        a_col = _resolve_col(alpha_col or seq_col, _CDR3A_ALIASES, "CDR3α")
        sequences: dict = {}
        for i, row in df.iterrows():
            seq = str(row[a_col]).upper().strip()
            if not seq or seq == "NAN":
                continue
            name = str(row[id_col]) if id_col else f"seq_{i + 1}"
            sequences[name] = seq
        return sequences

    if chain == "paired":
        b_col = _resolve_col(seq_col, _CDR3B_ALIASES, "CDR3β")
        a_col = _resolve_col(alpha_col, _CDR3A_ALIASES, "CDR3α")
        sequences = {}
        for i, row in df.iterrows():
            seqb = str(row[b_col]).upper().strip()
            seqa = str(row[a_col]).upper().strip()
            if not seqb or seqb == "NAN" or not seqa or seqa == "NAN":
                continue
            name = str(row[id_col]) if id_col else f"seq_{i + 1}"
            sequences[name] = (seqa, seqb)
        return sequences

    # beta (default)
    b_col = _resolve_col(seq_col, _CDR3B_ALIASES, "CDR3β")
    sequences = {}
    for i, row in df.iterrows():
        seq = str(row[b_col]).upper().strip()
        if not seq or seq == "NAN":
            continue
        name = str(row[id_col]) if id_col else f"seq_{i + 1}"
        sequences[name] = seq
    return sequences


def _parse_paired_fasta(beta_path: Path, alpha_path: Path) -> dict[str, tuple[str, str]]:
    """Merge two FASTA files into paired queries, matched by sequence name."""
    beta = parse_fasta(beta_path)
    alpha = parse_fasta(alpha_path)
    common = set(beta) & set(alpha)
    if not common:
        sys.exit("Error: no matching sequence names between --input and --input-alpha FASTA files.")
    only_beta  = set(beta) - common
    only_alpha = set(alpha) - common
    if only_beta:
        print(f"  Warning: {len(only_beta)} sequence(s) in --input have no alpha match — skipped.")
    if only_alpha:
        print(f"  Warning: {len(only_alpha)} sequence(s) in --input-alpha have no beta match — skipped.")
    return {name: (alpha[name], beta[name]) for name in sorted(common)}


def main():
    parser = argparse.ArgumentParser(
        prog="public-match",
        description="Match CDR3 sequences against public TCR databases.",
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
    parser.add_argument("--chain", default="beta", choices=["beta", "alpha", "paired"],
                        help="Chain mode: beta (default), alpha, or paired (both chains must match)")

    tsv_group = parser.add_argument_group("tabular input options (TSV/CSV)")
    tsv_group.add_argument("--seq-col", metavar="COL",
                           help="Column containing CDR3β sequences (auto-detected if omitted)")
    tsv_group.add_argument("--alpha-col", metavar="COL",
                           help="Column containing CDR3α sequences for --chain alpha/paired (auto-detected if omitted)")
    tsv_group.add_argument("--name-col", metavar="COL",
                           help="Column to use as sequence ID (auto-detected if omitted)")
    tsv_group.add_argument("--input-alpha", type=Path, metavar="PATH",
                           help="Second FASTA file with CDR3α sequences for --chain paired (matched by name to --input)")

    custom_group = parser.add_argument_group("custom database options")
    custom_group.add_argument("--custom-db", nargs="+", type=Path, metavar="PATH",
                              help="One or more custom database files (CSV or TSV)")
    custom_group.add_argument("--custom-db-cdr3-col", metavar="COL",
                              help="CDR3β column name in custom DB files (auto-detected if omitted)")

    args = parser.parse_args()

    if not args.input.exists():
        sys.exit(f"Error: input file not found: {args.input}")
    if args.input_alpha and not args.input_alpha.exists():
        sys.exit(f"Error: --input-alpha file not found: {args.input_alpha}")
    if args.chain == "paired" and args.input_alpha is None and _detect_format(args.input) == "fasta":
        sys.exit("Error: --chain paired with FASTA input requires --input-alpha for the CDR3α sequences.")

    print("Public.Match — CDR3 sequence search")
    print(f"  Input:     {args.input}")
    if args.input_alpha:
        print(f"  Input α:   {args.input_alpha}")
    print(f"  Chain:     {args.chain}")
    print(f"  Databases: {args.db}")
    if args.custom_db:
        print(f"  Custom DB: {[str(p) for p in args.custom_db]}")
    print(f"  Method:    {args.method}  |  Threshold: {args.threshold}")
    print()

    # Parse input sequences
    fmt = _detect_format(args.input)
    if args.chain == "paired" and args.input_alpha:
        sequences = _parse_paired_fasta(args.input, args.input_alpha)
    elif fmt == "fasta":
        sequences = parse_fasta(args.input)
    else:
        sequences = parse_tabular(
            args.input,
            seq_col=args.seq_col,
            name_col=args.name_col,
            alpha_col=args.alpha_col,
            chain=args.chain,
        )

    if not sequences:
        sys.exit("Error: no sequences found in input file.")
    print(f"Loaded {len(sequences)} query sequence(s).\n")

    print("Loading reference databases...")
    reference = load_databases(
        dbs=args.db,
        custom_paths=args.custom_db,
        custom_cdr3_col=args.custom_db_cdr3_col,
        chain=args.chain,
    )
    print(f"  Total reference entries: {len(reference):,}\n")

    print("Running matches...")
    queries = list(sequences.values())
    results = match(
        queries=queries,
        reference=reference,
        method=args.method,
        threshold=args.threshold,
        chain=args.chain,
    )

    # attach query names
    seq_to_name = {v: k for k, v in sequences.items()}
    results.insert(0, "query_name", results.iloc[:, 0].map(seq_to_name))

    results.to_csv(args.output, index=False)
    print(f"Done. {len(results)} match(es) written to {args.output}")


if __name__ == "__main__":
    main()

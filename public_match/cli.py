"""
Public.Match CLI

Usage
-----
python -m public_match \
    --input sequences.fasta \
    --output results.csv \
    --db iedb vdjdb mcpas tenx \
    --method blosum \
    --threshold 0.97
"""

import argparse
import sys
from pathlib import Path

from public_match.database import load_databases, ALL_DBS
from public_match.matcher import match


def parse_fasta(path: Path) -> dict[str, str]:
    """Parse a FASTA file into {name: sequence}."""
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


def main():
    parser = argparse.ArgumentParser(
        prog="public-match",
        description="Match CDR3β sequences against public TCR databases.",
    )
    parser.add_argument("--input", "-i", required=True, type=Path,
                        help="FASTA file with CDR3β sequences")
    parser.add_argument("--output", "-o", default=Path("public_match_results.csv"), type=Path,
                        help="Output CSV file (default: public_match_results.csv)")
    parser.add_argument("--db", nargs="+", default=ALL_DBS, choices=ALL_DBS,
                        help=f"Databases to search (default: all). Choices: {ALL_DBS}")
    parser.add_argument("--method", default="blosum", choices=["exact", "blosum", "edit"],
                        help="Matching method (default: blosum)")
    parser.add_argument("--threshold", type=float, default=0.97,
                        help="Score threshold: min BLOSUM62 score (blosum) or max edit distance (edit). Default: 0.97")

    args = parser.parse_args()

    if not args.input.exists():
        sys.exit(f"Error: input file not found: {args.input}")

    print(f"Public.Match — CDR3β sequence search")
    print(f"  Input:     {args.input}")
    print(f"  Databases: {args.db}")
    print(f"  Method:    {args.method}  |  Threshold: {args.threshold}")
    print()

    sequences = parse_fasta(args.input)
    if not sequences:
        sys.exit("Error: no sequences found in input file.")
    print(f"Loaded {len(sequences)} query sequence(s).\n")

    print("Loading reference databases...")
    reference = load_databases(args.db)
    print(f"  Total reference CDR3b entries: {len(reference):,}\n")

    print("Running matches...")
    results = match(
        queries=list(sequences.values()),
        reference=reference,
        method=args.method,
        threshold=args.threshold,
    )

    # attach query names back to results
    seq_to_name = {v: k for k, v in sequences.items()}
    results.insert(0, "query_name", results["query_cdr3b"].map(seq_to_name))

    results.to_csv(args.output, index=False)
    print(f"Done. {len(results)} match(es) written to {args.output}")


if __name__ == "__main__":
    main()

"""
CDR3b sequence matching engine.

Scoring methods
---------------
exact   — identical sequences only
blosum  — TCRMatch-style normalized BLOSUM62 score (same-length sequences)
          score = BLOSUM62(q, r) / sqrt(BLOSUM62(q,q) * BLOSUM62(r,r))
          threshold typically 0.97
edit    — Levenshtein edit distance, threshold in absolute AA differences
"""

import numpy as np
import pandas as pd
from typing import Literal

# BLOSUM62 matrix as a flat lookup: AA index (0-19) -> score matrix
_AA = "ACDEFGHIKLMNPQRSTVWY"
_AA_IDX = {aa: i for i, aa in enumerate(_AA)}

_BLOSUM62_RAW = np.array([
# A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
[ 4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2],  # A
[ 0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2],  # C
[-2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3],  # D
[-1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2],  # E
[-2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3],  # F
[ 0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3],  # G
[-2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2],  # H
[-1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1],  # I
[-1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2],  # K
[-1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1],  # L
[-1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1],  # M
[-2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2],  # N
[-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3],  # P
[-1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1],  # Q
[-1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2],  # R
[ 1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2],  # S
[ 0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2],  # T
[ 0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1],  # V
[-3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2],  # W
[-2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  7],  # Y
], dtype=np.float32)


def _encode(seq: str) -> np.ndarray | None:
    """Convert AA string to index array; return None if unknown residue."""
    idx = []
    for aa in seq:
        i = _AA_IDX.get(aa)
        if i is None:
            return None
        idx.append(i)
    return np.array(idx, dtype=np.int16)


def _blosum_self(idx: np.ndarray) -> float:
    return float(_BLOSUM62_RAW[idx, idx].sum())


def _blosum_pair(idx_q: np.ndarray, idx_r: np.ndarray) -> float:
    return float(_BLOSUM62_RAW[idx_q, idx_r].sum())


def _levenshtein(a: str, b: str) -> int:
    if len(a) < len(b):
        a, b = b, a
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a):
        curr = [i + 1]
        for j, cb in enumerate(b):
            curr.append(min(prev[j + 1] + 1, curr[j] + 1, prev[j] + (ca != cb)))
        prev = curr
    return prev[-1]


def match(
    queries: list[str],
    reference: pd.DataFrame,
    method: Literal["exact", "blosum", "edit"] = "blosum",
    threshold: float = 0.97,
) -> pd.DataFrame:
    """
    Match query CDR3b sequences against the reference database.

    Parameters
    ----------
    queries     : list of CDR3b amino acid strings
    reference   : DataFrame with columns [cdr3b, epitope, antigen, pathogen, HLA, source_db]
    method      : scoring method
    threshold   : minimum score (blosum) or maximum edit distance (edit)

    Returns
    -------
    DataFrame with columns [query_cdr3b, cdr3b, score, epitope, antigen, pathogen, HLA, source_db]
    """
    results = []
    ref_seqs = reference["cdr3b"].tolist()

    if method == "exact":
        query_set = set(queries)
        mask = reference["cdr3b"].isin(query_set)
        hits = reference[mask].copy()
        hits.insert(0, "query_cdr3b", hits["cdr3b"])
        hits["score"] = 1.0
        return hits.reset_index(drop=True)

    if method == "blosum":
        # group reference by length for fast pre-filtering
        ref_by_len: dict[int, list[tuple[int, np.ndarray]]] = {}
        ref_self: dict[int, list[float]] = {}
        for i, seq in enumerate(ref_seqs):
            idx = _encode(seq)
            if idx is None:
                continue
            L = len(seq)
            ref_by_len.setdefault(L, []).append((i, idx))
            ref_self.setdefault(L, []).append(_blosum_self(idx))

        for query in queries:
            q_idx = _encode(query.upper())
            if q_idx is None:
                continue
            L = len(query)
            if L not in ref_by_len:
                continue
            q_self = _blosum_self(q_idx)
            for (ref_row, r_idx), r_self in zip(ref_by_len[L], ref_self[L]):
                denom = np.sqrt(q_self * r_self)
                if denom == 0:
                    continue
                score = _blosum_pair(q_idx, r_idx) / denom
                if score >= threshold:
                    row = reference.iloc[ref_row].to_dict()
                    row["query_cdr3b"] = query
                    row["score"] = round(float(score), 4)
                    results.append(row)

    elif method == "edit":
        max_dist = int(threshold)
        for query in queries:
            q = query.upper()
            for i, ref_seq in enumerate(ref_seqs):
                if abs(len(q) - len(ref_seq)) > max_dist:
                    continue
                dist = _levenshtein(q, ref_seq)
                if dist <= max_dist:
                    row = reference.iloc[i].to_dict()
                    row["query_cdr3b"] = query
                    row["score"] = dist
                    results.append(row)

    if not results:
        return pd.DataFrame(columns=["query_cdr3b", "cdr3b", "score", "epitope", "antigen", "pathogen", "HLA", "source_db"])

    out = pd.DataFrame(results)
    cols = ["query_cdr3b", "cdr3b", "score", "epitope", "antigen", "pathogen", "HLA", "source_db"]
    return out[cols].sort_values(["query_cdr3b", "score"], ascending=[True, False]).reset_index(drop=True)

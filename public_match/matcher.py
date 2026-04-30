"""
CDR3 sequence matching engine.

Scoring methods
---------------
exact   — identical sequences only
blosum  — TCRMatch-style normalized BLOSUM62 score (same-length sequences)
          score = BLOSUM62(q, r) / sqrt(BLOSUM62(q,q) * BLOSUM62(r,r))
          threshold typically 0.97
edit    — Levenshtein edit distance, threshold in absolute AA differences

Chain modes
-----------
beta    — match CDR3β only (default)
alpha   — match CDR3α only; queries are CDR3α strings
paired  — match both chains; queries are (cdr3a, cdr3b) tuples; both must
          independently meet the threshold; score = mean of the two scores
"""

import numpy as np
import pandas as pd
from typing import Literal

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


def _blosum_score(q_idx, r_idx, q_self, r_self):
    denom = np.sqrt(q_self * r_self)
    if denom == 0:
        return None
    return _blosum_pair(q_idx, r_idx) / denom


def _index_by_len(seqs):
    """Group encoded sequences by length for fast BLOSUM pre-filtering."""
    by_len: dict[int, list[tuple[int, np.ndarray]]] = {}
    self_scores: dict[int, list[float]] = {}
    for i, seq in enumerate(seqs):
        idx = _encode(seq)
        if idx is None:
            continue
        L = len(seq)
        by_len.setdefault(L, []).append((i, idx))
        self_scores.setdefault(L, []).append(_blosum_self(idx))
    return by_len, self_scores


_BETA_COLS  = ["query_cdr3b", "cdr3b", "score", "epitope", "antigen", "pathogen", "HLA", "source_db"]
_ALPHA_COLS = ["query_cdr3a", "cdr3a", "score", "epitope", "antigen", "pathogen", "HLA", "source_db"]
_PAIRED_COLS = ["query_cdr3a", "query_cdr3b", "cdr3a", "cdr3b", "score", "epitope", "antigen", "pathogen", "HLA", "source_db"]


def match(
    queries: list[str] | list[tuple[str, str]],
    reference: pd.DataFrame,
    method: Literal["exact", "blosum", "edit"] = "blosum",
    threshold: float = 0.97,
    chain: Literal["beta", "alpha", "paired"] = "beta",
) -> pd.DataFrame:
    """
    Match query CDR3 sequences against the reference database.

    Parameters
    ----------
    queries   : For chain='beta'/'alpha': list of CDR3 strings.
                For chain='paired': list of (cdr3a, cdr3b) tuples.
    reference : Unified reference DataFrame.
    method    : 'exact', 'blosum', or 'edit'.
    threshold : Min BLOSUM score or max edit distance.
    chain     : 'beta' (default), 'alpha', or 'paired'.

    Returns
    -------
    DataFrame with match results.
    """
    if chain == "beta":
        return _match_single(queries, reference, method, threshold, col="cdr3b")
    if chain == "alpha":
        return _match_single(queries, reference, method, threshold, col="cdr3a")
    # paired
    return _match_paired(queries, reference, method, threshold)


def _match_single(
    queries: list[str],
    reference: pd.DataFrame,
    method: str,
    threshold: float,
    col: str,
) -> pd.DataFrame:
    query_col = f"query_{col}"
    empty = pd.DataFrame(columns=[query_col, col, "score", "epitope", "antigen", "pathogen", "HLA", "source_db"])
    results = []
    ref_seqs = reference[col].tolist()

    if method == "exact":
        query_set = set(q.upper() for q in queries)
        mask = reference[col].isin(query_set)
        hits = reference[mask].copy()
        hits.insert(0, query_col, hits[col])
        hits["score"] = 1.0
        return hits.reset_index(drop=True)[[query_col, col, "score", "epitope", "antigen", "pathogen", "HLA", "source_db"]]

    if method == "blosum":
        by_len, self_scores = _index_by_len(ref_seqs)
        for query in queries:
            q = query.upper()
            q_idx = _encode(q)
            if q_idx is None:
                continue
            L = len(q)
            if L not in by_len:
                continue
            q_self = _blosum_self(q_idx)
            for (ref_row, r_idx), r_self in zip(by_len[L], self_scores[L]):
                score = _blosum_score(q_idx, r_idx, q_self, r_self)
                if score is not None and score >= threshold:
                    row = reference.iloc[ref_row].to_dict()
                    row[query_col] = q
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
                    row[query_col] = q
                    row["score"] = dist
                    results.append(row)

    if not results:
        return empty

    out = pd.DataFrame(results)
    out_cols = [query_col, col, "score", "epitope", "antigen", "pathogen", "HLA", "source_db"]
    return out[out_cols].sort_values([query_col, "score"], ascending=[True, False]).reset_index(drop=True)


def _match_paired(
    queries: list[tuple[str, str]],
    reference: pd.DataFrame,
    method: str,
    threshold: float,
) -> pd.DataFrame:
    empty = pd.DataFrame(columns=_PAIRED_COLS)
    results = []
    ref_a = reference["cdr3a"].tolist()
    ref_b = reference["cdr3b"].tolist()

    if method == "exact":
        for cdr3a_q, cdr3b_q in queries:
            qa, qb = cdr3a_q.upper(), cdr3b_q.upper()
            mask = (reference["cdr3a"] == qa) & (reference["cdr3b"] == qb)
            for _, row in reference[mask].iterrows():
                r = row.to_dict()
                r["query_cdr3a"] = qa
                r["query_cdr3b"] = qb
                r["score"] = 1.0
                results.append(r)

    elif method == "blosum":
        by_len_a, self_a = _index_by_len(ref_a)
        by_len_b, self_b = _index_by_len(ref_b)
        # build a reverse map: ref_row → (La, idx_a, self_a, Lb, idx_b, self_b)
        ref_encoded = {}
        for L, entries in by_len_a.items():
            for pos, (ref_row, idx) in enumerate(entries):
                ref_encoded.setdefault(ref_row, {})["a"] = (L, idx, self_a[L][pos])
        for L, entries in by_len_b.items():
            for pos, (ref_row, idx) in enumerate(entries):
                ref_encoded.setdefault(ref_row, {})["b"] = (L, idx, self_b[L][pos])

        for cdr3a_q, cdr3b_q in queries:
            qa, qb = cdr3a_q.upper(), cdr3b_q.upper()
            qa_idx = _encode(qa)
            qb_idx = _encode(qb)
            if qa_idx is None or qb_idx is None:
                continue
            La, Lb = len(qa), len(qb)
            qa_self = _blosum_self(qa_idx)
            qb_self = _blosum_self(qb_idx)

            for ref_row, chains in ref_encoded.items():
                if "a" not in chains or "b" not in chains:
                    continue
                ra_L, ra_idx, ra_self = chains["a"]
                rb_L, rb_idx, rb_self = chains["b"]
                if ra_L != La or rb_L != Lb:
                    continue
                score_a = _blosum_score(qa_idx, ra_idx, qa_self, ra_self)
                score_b = _blosum_score(qb_idx, rb_idx, qb_self, rb_self)
                if score_a is None or score_b is None:
                    continue
                if score_a >= threshold and score_b >= threshold:
                    row = reference.iloc[ref_row].to_dict()
                    row["query_cdr3a"] = qa
                    row["query_cdr3b"] = qb
                    row["score"] = round((score_a + score_b) / 2, 4)
                    results.append(row)

    elif method == "edit":
        max_dist = int(threshold)
        for cdr3a_q, cdr3b_q in queries:
            qa, qb = cdr3a_q.upper(), cdr3b_q.upper()
            for i, (ra, rb) in enumerate(zip(ref_a, ref_b)):
                if abs(len(qa) - len(ra)) > max_dist or abs(len(qb) - len(rb)) > max_dist:
                    continue
                dist_a = _levenshtein(qa, ra)
                if dist_a > max_dist:
                    continue
                dist_b = _levenshtein(qb, rb)
                if dist_b <= max_dist:
                    row = reference.iloc[i].to_dict()
                    row["query_cdr3a"] = qa
                    row["query_cdr3b"] = qb
                    row["score"] = dist_a + dist_b
                    results.append(row)

    if not results:
        return empty

    out = pd.DataFrame(results)
    return out[_PAIRED_COLS].sort_values(["query_cdr3a", "score"], ascending=[True, False]).reset_index(drop=True)

"""
Epitope-to-CDR3 reverse matching.

Given a list of epitope peptide sequences, find all reference entries
where the epitope column matches (exact, BLOSUM62, or edit distance).
Returns a DataFrame of CDR3α/β sequences that recognise those epitopes.

Reuses scoring helpers from matcher.py; does NOT modify that module.
"""

import pandas as pd
from typing import Literal

from .matcher import _encode, _blosum_self, _blosum_score, _levenshtein, _index_by_len

_OUTPUT_COLS = [
    "query_epitope", "epitope", "score",
    "cdr3a", "cdr3b", "HLA", "antigen", "pathogen", "source_db",
]


def match_by_epitope(
    epitopes: list[str],
    reference: pd.DataFrame,
    method: Literal["exact", "blosum", "edit"] = "exact",
    threshold: float = 0.9,
) -> pd.DataFrame:
    """
    Match query epitope sequences against the reference database.

    Parameters
    ----------
    epitopes  : List of epitope peptide strings to search for.
    reference : Unified reference DataFrame (must contain 'epitope' column).
    method    : 'exact', 'blosum', or 'edit'.
    threshold : Min BLOSUM score (blosum) or max edit distance (edit).

    Returns
    -------
    DataFrame with query_epitope prepended; one row per matching CDR3 entry.
    """
    ref = reference.dropna(subset=["epitope"]).copy()
    ref = ref[ref["epitope"].str.strip() != ""].reset_index(drop=True)

    empty = pd.DataFrame(columns=_OUTPUT_COLS)
    if ref.empty or not epitopes:
        return empty

    for col in _OUTPUT_COLS:
        if col not in ref.columns:
            ref[col] = None

    ref_epitopes = ref["epitope"].tolist()
    results = []

    if method == "exact":
        query_set = {e.upper().strip() for e in epitopes}
        mask = ref["epitope"].str.upper().str.strip().isin(query_set)
        hits = ref[mask].copy()
        if hits.empty:
            return empty
        hits.insert(0, "query_epitope", hits["epitope"])
        hits["score"] = 1.0
        return hits[_OUTPUT_COLS].reset_index(drop=True)

    if method == "blosum":
        by_len, self_scores = _index_by_len(ref_epitopes)
        for epi in epitopes:
            q = epi.upper().strip()
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
                    row = ref.iloc[ref_row].to_dict()
                    row["query_epitope"] = q
                    row["score"] = round(float(score), 4)
                    results.append(row)

    elif method == "edit":
        max_dist = int(threshold)
        for epi in epitopes:
            q = epi.upper().strip()
            for i, ref_epi in enumerate(ref_epitopes):
                if not isinstance(ref_epi, str):
                    continue
                if abs(len(q) - len(ref_epi)) > max_dist:
                    continue
                dist = _levenshtein(q, ref_epi)
                if dist <= max_dist:
                    row = ref.iloc[i].to_dict()
                    row["query_epitope"] = q
                    row["score"] = dist
                    results.append(row)

    if not results:
        return empty

    out = pd.DataFrame(results)
    for col in _OUTPUT_COLS:
        if col not in out.columns:
            out[col] = None

    ascending_score = method == "edit"
    return (
        out[_OUTPUT_COLS]
        .sort_values(["query_epitope", "score"], ascending=[True, ascending_score])
        .reset_index(drop=True)
    )

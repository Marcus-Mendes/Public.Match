"""
Public.Match — Streamlit web UI
Run from the repo root: streamlit run app.py
"""

import streamlit as st
import pandas as pd

from public_match.database import load_databases_cached, build_cache, ALL_DBS, CACHE_PATH
from public_match.matcher import match

st.set_page_config(page_title="Public.Match", page_icon="🧬", layout="wide")

DB_LABELS = {
    "iedb": "IEDB",
    "vdjdb": "VDJdb",
    "mcpas": "McPAS-TCR",
    "tenx": "10x Genomics pMHC",
    "mixtcrpred": "MixTCRpred",
    "batcave": "BATCAVE",
}

_SOURCE_LABELS = {
    "iedb": "IEDB",
    "vdjdb": "VDJdb",
    "mcpas": "McPAS",
    "tenx": "10xDcode",
    "mixtcrpred": "MixTCRpred",
    "batcave": "BATCAVE",
}

EXAMPLE_FASTA = (
    ">cell_001\nCASSLAPGATNEKLFF\n"
    ">cell_002\nGILGFVFTL\n"
    ">cell_003\nCASSLGQTNEKLFF\n"
    ">cell_004\nCASSPGTGASYEQYF\n"
    ">cell_005\nCASSFRGGAFF"
)


def parse_fasta(text: str) -> dict[str, str]:
    sequences: dict[str, str] = {}
    current_name, current_seq = None, []
    is_fasta = ">" in text
    for line in text.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_name:
                sequences[current_name] = "".join(current_seq).upper()
            current_name = line[1:].split()[0]
            current_seq = []
        elif is_fasta:
            current_seq.append(line)
        else:
            seq_id = f"seq_{len(sequences) + 1}"
            sequences[seq_id] = line.upper()
    if is_fasta and current_name:
        sequences[current_name] = "".join(current_seq).upper()
    return sequences


# ── Sidebar ────────────────────────────────────────────────────────────────────

with st.sidebar:
    st.title("⚙️ Settings")

    selected_dbs = st.multiselect(
        "Databases to search",
        options=ALL_DBS,
        default=ALL_DBS,
        format_func=lambda x: DB_LABELS.get(x, x),
    )

    method = st.selectbox(
        "Matching method",
        options=["blosum", "exact", "edit"],
        index=0,
        help=(
            "**BLOSUM62**: TCRMatch-style normalized substitution score (0–1).  \n"
            "**Edit**: Levenshtein distance (integer).  \n"
            "**Exact**: identical sequences only."
        ),
    )

    if method == "blosum":
        threshold = st.slider("Min BLOSUM62 score", 0.80, 1.00, 0.97, 0.01,
                              help="0.97 is the standard TCRMatch threshold")
    elif method == "edit":
        threshold = float(st.slider("Max edit distance (AA)", 0, 5, 1, 1))
    else:
        threshold = 1.0

    st.divider()

    if st.button("Rebuild database cache",
                 help="Re-run after updating source database files"):
        with st.spinner("Building cache from source files…"):
            st.cache_resource.clear()
            build_cache()
        st.success("Cache rebuilt!")
        st.rerun()

    st.divider()
    st.caption("Public.Match · BTC Hackathon 2026")


# ── Header ─────────────────────────────────────────────────────────────────────

st.title("🧬 Public.Match")
st.markdown(
    "Match patient CDR3β sequences against public TCR databases — "
    "IEDB, VDJdb, McPAS-TCR, 10x Genomics, MixTCRpred, and BATCAVE — in one query."
)

st.divider()


# ── Eager database load ────────────────────────────────────────────────────────

@st.cache_resource
def get_all_reference() -> pd.DataFrame:
    return load_databases_cached(ALL_DBS)


status = st.empty()

if not CACHE_PATH.exists():
    status.info(
        "⏳ First launch — building database cache from source files. "
        "This takes 1–2 minutes and only happens once."
    )

with st.spinner("Loading reference databases…"):
    all_reference = get_all_reference()

db_counts = all_reference["source_db"].value_counts()
status.success(
    f"✓ **{len(all_reference):,} reference entries** loaded across "
    f"**{db_counts.shape[0]} databases** — "
    + "  |  ".join(f"{k}: {v:,}" for k, v in db_counts.items())
)

# Filter to user-selected databases
if set(selected_dbs) == set(ALL_DBS):
    reference = all_reference
else:
    labels = [_SOURCE_LABELS[d] for d in selected_dbs]
    reference = all_reference[all_reference["source_db"].isin(labels)].reset_index(drop=True)

st.divider()


# ── Input ──────────────────────────────────────────────────────────────────────

st.subheader("Input sequences")

input_tab, example_tab = st.tabs(["Upload / paste", "Example input"])

sequences: dict[str, str] = {}

with input_tab:
    input_mode = st.radio("Input method", ["Upload FASTA", "Paste sequences"], horizontal=True)

    raw_text = ""
    if input_mode == "Upload FASTA":
        uploaded = st.file_uploader("FASTA file (.fasta / .fa / .txt)", type=["fasta", "fa", "txt"])
        if uploaded:
            raw_text = uploaded.read().decode("utf-8")
    else:
        raw_text = st.text_area(
            "Paste CDR3β sequences (FASTA format or one sequence per line)",
            height=180,
            placeholder=">cell_001\nCASSLAPGATNEKLFF\n>cell_002\nCASSLGQTNEKLFF",
        )

    if raw_text.strip():
        sequences = parse_fasta(raw_text)
        st.success(f"✓ {len(sequences)} sequence(s) loaded")

with example_tab:
    st.code(EXAMPLE_FASTA, language="text")
    if st.button("Load example"):
        sequences = parse_fasta(EXAMPLE_FASTA)
        st.success(f"✓ {len(sequences)} example sequences loaded")

st.divider()


# ── Run ────────────────────────────────────────────────────────────────────────

run_disabled = not sequences or not selected_dbs
st.button("▶  Run Public.Match", type="primary", disabled=run_disabled, key="run_btn")

if not sequences:
    st.info("Upload or paste CDR3β sequences above to get started.")

if st.session_state.get("run_btn"):
    with st.spinner(f"Matching {len(sequences)} sequence(s) against {len(reference):,} reference entries…"):
        results = match(
            queries=list(sequences.values()),
            reference=reference,
            method=method,
            threshold=threshold,
        )

    if not results.empty:
        seq_to_name = {v: k for k, v in sequences.items()}
        results.insert(0, "query_name", results["query_cdr3b"].map(seq_to_name))

    # ── Results ────────────────────────────────────────────────────────────────

    st.subheader("Results")

    if results.empty:
        st.warning(
            "No matches found. Try lowering the BLOSUM62 threshold, "
            "increasing the edit distance, or selecting more databases."
        )
    else:
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Total matches", f"{len(results):,}")
        c2.metric("Queries with hits", results["query_cdr3b"].nunique())
        c3.metric("Unique epitopes", results["epitope"].nunique())
        c4.metric("Databases with hits", results["source_db"].nunique())

        st.markdown("**Hits by database**")
        db_hit_counts = (
            results["source_db"]
            .value_counts()
            .rename_axis("Database")
            .reset_index(name="Hits")
        )
        st.bar_chart(db_hit_counts.set_index("Database"))

        st.markdown("**All matches**")
        st.dataframe(results, use_container_width=True, height=420)

        st.download_button(
            "⬇ Download results CSV",
            data=results.to_csv(index=False),
            file_name="public_match_results.csv",
            mime="text/csv",
        )

"""
Public.Match — Streamlit web UI
Run from the repo root: streamlit run app.py
"""

import streamlit as st
import pandas as pd

from public_match.database import load_databases_cached, build_cache, ALL_DBS, CACHE_PATH
from public_match.matcher import match

st.set_page_config(
    page_title="Public.Match — Break Through Cancer",
    page_icon="🧬",
    layout="wide",
)

# ── BTC Theme ──────────────────────────────────────────────────────────────────

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');

/* ── Global ── */
html, body, .stApp {
    font-family: 'Inter', sans-serif;
    background-color: #f4f7fb;
}

/* ── Sidebar ── */
[data-testid="stSidebar"] {
    background: linear-gradient(180deg, #0b1f3a 0%, #112b52 100%);
    border-right: 1px solid #1a3a6b;
}
[data-testid="stSidebar"] label,
[data-testid="stSidebar"] p,
[data-testid="stSidebar"] span,
[data-testid="stSidebar"] div,
[data-testid="stSidebar"] .stMarkdown {
    color: #ccd9eb !important;
}
[data-testid="stSidebar"] h1,
[data-testid="stSidebar"] h2,
[data-testid="stSidebar"] h3 {
    color: #ffffff !important;
}
[data-testid="stSidebar"] .stSelectbox > div > div,
[data-testid="stSidebar"] .stMultiSelect > div > div {
    background-color: #1a3a6b !important;
    border-color: #2a5298 !important;
    color: #ffffff !important;
}
[data-testid="stSidebar"] hr {
    border-color: #1a3a6b !important;
}

/* ── Main area ── */
[data-testid="stAppViewContainer"] > .main {
    background-color: #f4f7fb;
}
[data-testid="block-container"] {
    padding-top: 1.5rem;
}

/* ── Headings ── */
h1 { color: #0b1f3a; font-weight: 700; }
h2, h3 { color: #1e4d8f; font-weight: 600; }

/* ── Primary button ── */
.stButton > button[kind="primary"] {
    background-color: #1e6ab0;
    color: #ffffff;
    border: none;
    border-radius: 5px;
    font-weight: 600;
    letter-spacing: 0.4px;
    padding: 0.5rem 2.5rem;
    font-size: 1rem;
    transition: background-color 0.2s;
}
.stButton > button[kind="primary"]:hover {
    background-color: #155a9a;
    color: #ffffff;
}
.stButton > button[kind="primary"]:disabled {
    background-color: #8aabc8;
    color: #dce9f5;
}

/* ── Secondary buttons ── */
.stButton > button:not([kind="primary"]) {
    background-color: transparent;
    border: 1.5px solid #1e6ab0;
    color: #1e6ab0;
    border-radius: 5px;
    font-weight: 500;
    transition: all 0.2s;
}
.stButton > button:not([kind="primary"]):hover {
    background-color: #1e6ab0;
    color: #ffffff;
}

/* ── Metrics ── */
[data-testid="stMetric"] {
    background-color: #ffffff;
    border: 1px solid #dde6f0;
    border-top: 4px solid #1e6ab0;
    border-radius: 6px;
    padding: 1rem 1.25rem;
    box-shadow: 0 1px 4px rgba(0,0,0,0.05);
}
[data-testid="stMetricLabel"] { color: #5a7a9a !important; font-size: 0.8rem !important; text-transform: uppercase; letter-spacing: 0.5px; }
[data-testid="stMetricValue"] { color: #0b1f3a !important; font-weight: 700 !important; }

/* ── Alerts ── */
[data-testid="stAlert"] {
    border-radius: 6px;
}

/* ── Tabs ── */
.stTabs [data-baseweb="tab-list"] {
    border-bottom: 2px solid #dde6f0;
}
.stTabs [data-baseweb="tab"] {
    color: #5a7a9a;
    font-weight: 500;
}
.stTabs [aria-selected="true"] {
    color: #1e6ab0 !important;
    border-bottom-color: #1e6ab0 !important;
}

/* ── Dataframe ── */
[data-testid="stDataFrame"] {
    border: 1px solid #dde6f0;
    border-radius: 6px;
    overflow: hidden;
}

/* ── Divider ── */
hr { border-color: #dde6f0 !important; }

/* ── Input widgets ── */
.stTextArea textarea, .stFileUploader {
    border-color: #b0c8e0 !important;
    border-radius: 6px !important;
}
.stTextArea textarea:focus {
    border-color: #1e6ab0 !important;
    box-shadow: 0 0 0 2px rgba(30,106,176,0.15) !important;
}
</style>
""", unsafe_allow_html=True)

# ── Constants ─────────────────────────────────────────────────────────────────

DB_LABELS = {
    "iedb": "IEDB",
    "vdjdb": "VDJdb",
    "mcpas": "McPAS-TCR",
    "tenx": "10x Genomics pMHC",
    "mixtcrpred": "MixTCRpred",
    "batcave": "BATCAVE",
    "neotcr": "NeoTCR",
}

_SOURCE_LABELS = {
    "iedb": "IEDB",
    "vdjdb": "VDJdb",
    "mcpas": "McPAS",
    "tenx": "10xDcode",
    "mixtcrpred": "MixTCRpred",
    "batcave": "BATCAVE",
    "neotcr": "NeoTCR",
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
    st.markdown("## ⚙️ Settings")

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

    if st.button("↺  Rebuild database cache",
                 help="Re-run after updating source database files"):
        with st.spinner("Building cache from source files…"):
            st.cache_resource.clear()
            build_cache()
        st.success("Cache rebuilt!")
        st.rerun()

    st.divider()
    st.markdown(
        "<p style='font-size:0.75rem; color:#6a8aaa;'>Public.Match · BTC Hackathon 2026<br>"
        "Karchin Lab · Johns Hopkins</p>",
        unsafe_allow_html=True,
    )


# ── Header banner ──────────────────────────────────────────────────────────────

st.markdown("""
<div style="
    background: linear-gradient(135deg, #0b1f3a 0%, #1a3a6b 100%);
    border-radius: 10px;
    padding: 2rem 2.5rem;
    margin-bottom: 1.5rem;
    display: flex;
    align-items: center;
    gap: 1.5rem;
">
    <div>
        <div style="display:flex; align-items:center; gap:0.5rem; margin-bottom:0.4rem;">
            <span style="font-size:1.8rem;">🧬</span>
            <span style="font-size:1.8rem; font-weight:700; color:#ffffff; letter-spacing:-0.5px;">Public.Match</span>
        </div>
        <p style="color:#a8c4e0; margin:0; font-size:0.95rem; max-width:680px;">
            Match patient CDR3&beta; sequences against public TCR databases &mdash;
            IEDB, VDJdb, McPAS-TCR, 10x Genomics, MixTCRpred, BATCAVE, and NeoTCR &mdash; in a single query.
        </p>
        <div style="display:flex; gap:6px; margin-top:0.9rem;">
            <span style="background:#f0912a; border-radius:3px; width:10px; height:10px; display:inline-block;"></span>
            <span style="background:#3db5b0; border-radius:3px; width:10px; height:10px; display:inline-block;"></span>
            <span style="background:#8b5cf6; border-radius:3px; width:10px; height:10px; display:inline-block;"></span>
            <span style="background:#d946a8; border-radius:3px; width:10px; height:10px; display:inline-block;"></span>
            <span style="background:#ffffff; border-radius:3px; width:10px; height:10px; display:inline-block;"></span>
            <span style="background:#1e6ab0; border-radius:3px; width:10px; height:10px; display:inline-block;"></span>
        </div>
    </div>
</div>
""", unsafe_allow_html=True)


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
status.markdown(
    f"""<div style="background:#eaf4ee; border:1px solid #b2ddc0; border-left:4px solid #2e9e5b;
    border-radius:6px; padding:0.75rem 1rem; font-size:0.9rem; color:#1a4d2e;">
    ✅ &nbsp;<strong>{len(all_reference):,} reference entries</strong> loaded across
    <strong>{db_counts.shape[0]} databases</strong> &nbsp;—&nbsp;
    {"&nbsp; | &nbsp;".join(f"<strong>{k}</strong>: {v:,}" for k, v in db_counts.items())}
    </div>""",
    unsafe_allow_html=True,
)

# Filter to selected databases
if set(selected_dbs) == set(ALL_DBS):
    reference = all_reference
else:
    labels = [_SOURCE_LABELS[d] for d in selected_dbs]
    reference = all_reference[all_reference["source_db"].isin(labels)].reset_index(drop=True)

st.divider()


# ── Input ──────────────────────────────────────────────────────────────────────

st.subheader("Input sequences")

input_tab, example_tab = st.tabs(["📂  Upload / paste", "💡  Example input"])

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
    if st.button("Load example sequences"):
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

    # ── Results ───────────────────────────────────────────────────────────────

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

        st.markdown("#### Hits by database")
        db_hit_counts = (
            results["source_db"]
            .value_counts()
            .rename_axis("Database")
            .reset_index(name="Hits")
        )
        st.bar_chart(db_hit_counts.set_index("Database"), color="#1e6ab0")

        st.markdown("#### All matches")
        st.dataframe(results, use_container_width=True, height=420)

        st.download_button(
            "⬇  Download results CSV",
            data=results.to_csv(index=False),
            file_name="public_match_results.csv",
            mime="text/csv",
        )

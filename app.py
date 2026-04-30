"""
Public.Match — Streamlit web UI
Run from the repo root: streamlit run app.py
"""

import streamlit as st
import pandas as pd
from pathlib import Path

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

html, body, .stApp {
    font-family: 'Inter', sans-serif;
    background-color: #f4f7fb;
}
[data-testid="stSidebar"] {
    background: linear-gradient(180deg, #0b1f3a 0%, #112b52 100%);
    border-right: 1px solid #1a3a6b;
}
[data-testid="stSidebar"] label,
[data-testid="stSidebar"] p,
[data-testid="stSidebar"] span,
[data-testid="stSidebar"] div,
[data-testid="stSidebar"] .stMarkdown { color: #ccd9eb !important; }
[data-testid="stSidebar"] h1,
[data-testid="stSidebar"] h2,
[data-testid="stSidebar"] h3 { color: #ffffff !important; }
[data-testid="stSidebar"] .stSelectbox > div > div,
[data-testid="stSidebar"] .stMultiSelect > div > div {
    background-color: #1a3a6b !important;
    border-color: #2a5298 !important;
    color: #ffffff !important;
}
[data-testid="stSidebar"] hr { border-color: #1a3a6b !important; }
[data-testid="stAppViewContainer"] > .main { background-color: #f4f7fb; }
[data-testid="block-container"] { padding-top: 1.5rem; }
h1 { color: #0b1f3a; font-weight: 700; }
h2, h3 { color: #1e4d8f; font-weight: 600; }
.stButton > button[kind="primary"] {
    background-color: #1e6ab0; color: #ffffff; border: none;
    border-radius: 5px; font-weight: 600; letter-spacing: 0.4px;
    padding: 0.5rem 2.5rem; font-size: 1rem; transition: background-color 0.2s;
}
.stButton > button[kind="primary"]:hover { background-color: #155a9a; color: #ffffff; }
.stButton > button[kind="primary"]:disabled { background-color: #8aabc8; color: #dce9f5; }
.stButton > button:not([kind="primary"]) {
    background-color: transparent; border: 1.5px solid #1e6ab0;
    color: #1e6ab0; border-radius: 5px; font-weight: 500; transition: all 0.2s;
}
.stButton > button:not([kind="primary"]):hover { background-color: #1e6ab0; color: #ffffff; }
[data-testid="stMetric"] {
    background-color: #ffffff; border: 1px solid #dde6f0;
    border-top: 4px solid #1e6ab0; border-radius: 6px;
    padding: 1rem 1.25rem; box-shadow: 0 1px 4px rgba(0,0,0,0.05);
}
[data-testid="stMetricLabel"] { color: #5a7a9a !important; font-size: 0.8rem !important; text-transform: uppercase; letter-spacing: 0.5px; }
[data-testid="stMetricValue"] { color: #0b1f3a !important; font-weight: 700 !important; }
[data-testid="stAlert"] { border-radius: 6px; }
.stTabs [data-baseweb="tab-list"] { border-bottom: 2px solid #dde6f0; }
.stTabs [data-baseweb="tab"] { color: #5a7a9a; font-weight: 500; }
.stTabs [aria-selected="true"] { color: #1e6ab0 !important; border-bottom-color: #1e6ab0 !important; }
[data-testid="stDataFrame"] { border: 1px solid #dde6f0; border-radius: 6px; overflow: hidden; }
hr { border-color: #dde6f0 !important; }
.stTextArea textarea, .stFileUploader { border-color: #b0c8e0 !important; border-radius: 6px !important; }
.stTextArea textarea:focus { border-color: #1e6ab0 !important; box-shadow: 0 0 0 2px rgba(30,106,176,0.15) !important; }
</style>
""", unsafe_allow_html=True)

# ── Constants ──────────────────────────────────────────────────────────────────

DB_LABELS = {
    "iedb":       "IEDB",
    "vdjdb":      "VDJdb",
    "mcpas":      "McPAS-TCR",
    "tenx":       "10x Genomics pMHC",
    "mixtcrpred": "MixTCRpred",
    "batcave":    "BATCAVE",
    "neotcr":     "NeoTCR",
}

_SOURCE_LABELS = {
    "iedb":       "IEDB",
    "vdjdb":      "VDJdb",
    "mcpas":      "McPAS",
    "tenx":       "10xDcode",
    "mixtcrpred": "MixTCRpred",
    "batcave":    "BATCAVE",
    "neotcr":     "NeoTCR",
}

_CHAIN_LABELS = {
    "beta":   "CDR3β only",
    "alpha":  "CDR3α only",
    "paired": "Paired α + β",
}

def _read_example(filename: str) -> str:
    path = Path(filename)
    return path.read_text().strip() if path.exists() else ""

EXAMPLE_BETA         = _read_example("example_input.fasta")
EXAMPLE_ALPHA        = _read_example("example_input_alpha.fasta")
EXAMPLE_PAIRED_BETA  = _read_example("example_input_paired_beta.fasta")
EXAMPLE_PAIRED_ALPHA = _read_example("example_input_paired_alpha.fasta")

# CDR3b and CDR3a column name aliases for TSV auto-detection
_CDR3B_ALIASES = ["cdr3b", "cdr3_beta", "cdr3_b", "junction_aa", "cdr3", "CDR3", "TRB_CDR3"]
_CDR3A_ALIASES = ["cdr3a", "cdr3_alpha", "cdr3_a", "TRA_CDR3", "cdr3_TRA", "junction_aa_alpha"]
_NAME_ALIASES  = ["name", "id", "cell_id", "barcode", "clone_id", "sample_id"]


def _find_col(columns, aliases):
    cols_lower = {c.lower(): c for c in columns}
    for alias in aliases:
        if alias in columns:
            return alias
        if alias.lower() in cols_lower:
            return cols_lower[alias.lower()]
    return None


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

    chain = st.radio(
        "Chain mode",
        options=["beta", "alpha", "paired"],
        format_func=lambda x: _CHAIN_LABELS[x],
        help=(
            "**CDR3β only**: match beta chain sequences (default, uses fast cache).  \n"
            "**CDR3α only**: match alpha chain sequences.  \n"
            "**Paired α+β**: both chains must independently meet the threshold; "
            "score = mean of the two scores."
        ),
    )

    if chain != "beta":
        st.warning("⚠️ Alpha / paired modes load from source files — slower than beta.", icon=None)

    st.divider()

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
    border-radius: 10px; padding: 2rem 2.5rem; margin-bottom: 1.5rem;
">
    <div style="display:flex; align-items:center; gap:0.5rem; margin-bottom:0.4rem;">
        <span style="font-size:1.8rem;">🧬</span>
        <span style="font-size:1.8rem; font-weight:700; color:#ffffff; letter-spacing:-0.5px;">Public.Match</span>
    </div>
    <p style="color:#a8c4e0; margin:0; font-size:0.95rem; max-width:700px;">
        Match patient CDR3 sequences against public TCR databases &mdash;
        IEDB, VDJdb, McPAS-TCR, 10x Genomics, MixTCRpred, BATCAVE, and NeoTCR &mdash;
        across beta, alpha, or paired chain modes.
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
""", unsafe_allow_html=True)


# ── Reference database load ────────────────────────────────────────────────────

@st.cache_resource
def get_beta_reference() -> pd.DataFrame:
    return load_databases_cached(ALL_DBS, chain="beta")

@st.cache_resource
def get_alpha_reference() -> pd.DataFrame:
    return load_databases_cached(ALL_DBS, chain="alpha")

@st.cache_resource
def get_paired_reference() -> pd.DataFrame:
    return load_databases_cached(ALL_DBS, chain="paired")

_REFERENCE_LOADERS = {
    "beta":   get_beta_reference,
    "alpha":  get_alpha_reference,
    "paired": get_paired_reference,
}

status = st.empty()

if chain == "beta" and not CACHE_PATH.exists():
    status.info("⏳ First launch — building database cache. This takes 1–2 minutes and only happens once.")
elif chain != "beta":
    status.info(f"⏳ Loading source files for {_CHAIN_LABELS[chain]} mode…")

with st.spinner("Loading reference databases…"):
    all_reference = _REFERENCE_LOADERS[chain]()

db_counts = all_reference["source_db"].value_counts()
status.markdown(
    f"""<div style="background:#eaf4ee; border:1px solid #b2ddc0; border-left:4px solid #2e9e5b;
    border-radius:6px; padding:0.75rem 1rem; font-size:0.9rem; color:#1a4d2e;">
    ✅ &nbsp;<strong>{len(all_reference):,} reference entries</strong> loaded &nbsp;·&nbsp;
    mode: <strong>{_CHAIN_LABELS[chain]}</strong> &nbsp;—&nbsp;
    {"&nbsp; | &nbsp;".join(f"<strong>{k}</strong>: {v:,}" for k, v in db_counts.items())}
    </div>""",
    unsafe_allow_html=True,
)

# Filter to selected databases
if set(selected_dbs) != set(ALL_DBS):
    labels = [_SOURCE_LABELS[d] for d in selected_dbs]
    reference = all_reference[all_reference["source_db"].isin(labels)].reset_index(drop=True)
else:
    reference = all_reference

st.divider()


# ── Input ──────────────────────────────────────────────────────────────────────

st.subheader("Input sequences")

# Clear persisted queries when chain mode changes
if st.session_state.get("_last_chain") != chain:
    st.session_state["queries"] = {}
    st.session_state["_last_chain"] = chain

def _save_queries(q: dict):
    st.session_state["queries"] = q

if chain == "beta":
    input_tab, example_tab = st.tabs(["📂  Upload / paste", "💡  Example"])
    with input_tab:
        mode = st.radio("Input method", ["Upload FASTA", "Paste sequences"], horizontal=True)
        raw = ""
        if mode == "Upload FASTA":
            f = st.file_uploader("FASTA (.fasta / .fa / .txt)", type=["fasta", "fa", "txt"])
            if f:
                raw = f.read().decode("utf-8")
        else:
            raw = st.text_area("CDR3β sequences (FASTA or one per line)", height=160,
                               placeholder=">cell_001\nCASSLAPGATNEKLFF")
        if raw.strip():
            _save_queries(parse_fasta(raw))
            st.success(f"✓ {len(st.session_state['queries'])} CDR3β sequence(s) loaded")
    with example_tab:
        st.code(EXAMPLE_BETA, language="text")
        if st.button("Load example"):
            _save_queries(parse_fasta(EXAMPLE_BETA))
            st.success(f"✓ {len(st.session_state['queries'])} example sequences loaded")

elif chain == "alpha":
    input_tab, example_tab = st.tabs(["📂  Upload / paste", "💡  Example"])
    with input_tab:
        mode = st.radio("Input method", ["Upload FASTA", "Paste sequences"], horizontal=True)
        raw = ""
        if mode == "Upload FASTA":
            f = st.file_uploader("FASTA (.fasta / .fa / .txt)", type=["fasta", "fa", "txt"])
            if f:
                raw = f.read().decode("utf-8")
        else:
            raw = st.text_area("CDR3α sequences (FASTA or one per line)", height=160,
                               placeholder=">cell_001\nCAVSANSGTYKYIF")
        if raw.strip():
            _save_queries(parse_fasta(raw))
            st.success(f"✓ {len(st.session_state['queries'])} CDR3α sequence(s) loaded")
    with example_tab:
        st.code(EXAMPLE_ALPHA, language="text")
        if st.button("Load example"):
            _save_queries(parse_fasta(EXAMPLE_ALPHA))
            st.success(f"✓ {len(st.session_state['queries'])} example sequences loaded")

else:  # paired
    tsv_tab, fasta_tab, example_tab = st.tabs(["📋  Upload TSV/CSV", "📂  Two FASTA files", "💡  Example"])

    with tsv_tab:
        st.markdown("Upload a TSV/CSV with **both** CDR3α and CDR3β columns. Column names are auto-detected.")
        tsv_file = st.file_uploader("TSV / CSV file", type=["tsv", "csv", "txt"])
        if tsv_file:
            sep = "\t" if tsv_file.name.endswith(".tsv") or tsv_file.name.endswith(".txt") else ","
            df_in = pd.read_csv(tsv_file, sep=sep)
            b_col = _find_col(df_in.columns, _CDR3B_ALIASES)
            a_col = _find_col(df_in.columns, _CDR3A_ALIASES)
            n_col = _find_col(df_in.columns, _NAME_ALIASES)
            if not b_col or not a_col:
                st.error(f"Could not find CDR3β ({_CDR3B_ALIASES[:3]}…) or CDR3α ({_CDR3A_ALIASES[:3]}…) columns. "
                         f"Columns found: {list(df_in.columns)}")
            else:
                parsed = {}
                for i, row in df_in.iterrows():
                    seqb = str(row[b_col]).upper().strip()
                    seqa = str(row[a_col]).upper().strip()
                    if not seqb or seqb == "NAN" or not seqa or seqa == "NAN":
                        continue
                    name = str(row[n_col]) if n_col else f"seq_{i+1}"
                    parsed[name] = (seqa, seqb)
                _save_queries(parsed)
                st.success(f"✓ {len(st.session_state['queries'])} paired sequence(s) loaded from {tsv_file.name} "
                           f"(β: `{b_col}`, α: `{a_col}`)")

    with fasta_tab:
        st.markdown("Upload two FASTA files — sequences are **matched by name**.")
        col1, col2 = st.columns(2)
        with col1:
            fb = st.file_uploader("CDR3β FASTA", type=["fasta", "fa", "txt"], key="paired_beta")
        with col2:
            fa = st.file_uploader("CDR3α FASTA", type=["fasta", "fa", "txt"], key="paired_alpha")
        if fb and fa:
            beta_seqs  = parse_fasta(fb.read().decode("utf-8"))
            alpha_seqs = parse_fasta(fa.read().decode("utf-8"))
            common = set(beta_seqs) & set(alpha_seqs)
            only_b = set(beta_seqs) - common
            only_a = set(alpha_seqs) - common
            if not common:
                st.error("No matching sequence names between the two files.")
            else:
                _save_queries({name: (alpha_seqs[name], beta_seqs[name]) for name in sorted(common)})
                msg = f"✓ {len(st.session_state['queries'])} paired sequence(s) matched by name."
                if only_b:
                    msg += f" ({len(only_b)} β-only skipped)"
                if only_a:
                    msg += f" ({len(only_a)} α-only skipped)"
                st.success(msg)

    with example_tab:
        st.markdown("Example paired input (alpha + beta, matched by name):")
        col1, col2 = st.columns(2)
        with col1:
            st.markdown("**CDR3β**")
            st.code(EXAMPLE_PAIRED_BETA, language="text")
        with col2:
            st.markdown("**CDR3α**")
            st.code(EXAMPLE_PAIRED_ALPHA, language="text")
        if st.button("Load example"):
            beta_seqs  = parse_fasta(EXAMPLE_PAIRED_BETA)
            alpha_seqs = parse_fasta(EXAMPLE_PAIRED_ALPHA)
            common = set(beta_seqs) & set(alpha_seqs)
            _save_queries({name: (alpha_seqs[name], beta_seqs[name]) for name in sorted(common)})
            st.success(f"✓ {len(st.session_state['queries'])} example paired sequences loaded")

st.divider()


# ── Run ────────────────────────────────────────────────────────────────────────

queries = st.session_state.get("queries", {})

run_disabled = not queries or not selected_dbs
st.button("▶  Run Public.Match", type="primary", disabled=run_disabled, key="run_btn")

if not queries:
    st.info(f"Provide {'paired' if chain == 'paired' else 'CDR3' + ('β' if chain == 'beta' else 'α')} "
            f"sequences above to get started.")

if st.session_state.get("run_btn"):
    query_list = list(queries.values())
    with st.spinner(f"Matching {len(query_list)} sequence(s) against {len(reference):,} reference entries…"):
        results = match(
            queries=query_list,
            reference=reference,
            method=method,
            threshold=threshold,
            chain=chain,
        )

    # Attach query names
    if not results.empty:
        if chain == "paired":
            tuple_to_name = {v: k for k, v in queries.items()}
            results.insert(0, "query_name",
                           list(zip(results["query_cdr3a"], results["query_cdr3b"]))
                           if "query_cdr3a" in results.columns
                           else results.get("query_cdr3b", pd.Series()))
            results["query_name"] = results["query_name"].map(
                lambda t: tuple_to_name.get(t, str(t)) if isinstance(t, tuple) else t
            )
        elif chain == "alpha":
            seq_to_name = {v: k for k, v in queries.items()}
            results.insert(0, "query_name", results["query_cdr3a"].map(seq_to_name))
        else:
            seq_to_name = {v: k for k, v in queries.items()}
            results.insert(0, "query_name", results["query_cdr3b"].map(seq_to_name))

    # ── Results ───────────────────────────────────────────────────────────────

    st.subheader("Results")

    if results.empty:
        st.warning(
            "No matches found. Try lowering the BLOSUM62 threshold, "
            "increasing the edit distance, or selecting more databases."
        )
    else:
        query_col = "query_cdr3b" if chain in ("beta", "paired") else "query_cdr3a"

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Total matches", f"{len(results):,}")
        c2.metric("Queries with hits", results[query_col].nunique())
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

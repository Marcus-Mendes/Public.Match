"""
Public.Match v2 — AlphaFold-inspired UI
Run from repo root: streamlit run app_v2.py
"""

import base64
import streamlit as st
import pandas as pd
from pathlib import Path

from public_match.database import load_databases_cached, build_cache, ALL_DBS, CACHE_PATH
from public_match.matcher import match

st.set_page_config(
    page_title="Public.Match",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# ── Helpers ────────────────────────────────────────────────────────────────────

def _b64(path: str) -> str | None:
    p = Path(path)
    return base64.b64encode(p.read_bytes()).decode() if p.exists() else None

def _read_example(filename: str) -> str:
    p = Path(filename)
    return p.read_text().strip() if p.exists() else ""

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
            sequences[f"seq_{len(sequences)+1}"] = line.upper()
    if is_fasta and current_name:
        sequences[current_name] = "".join(current_seq).upper()
    return sequences

EXAMPLE_BETA        = _read_example("example_input.fasta")
EXAMPLE_ALPHA       = _read_example("example_input_alpha.fasta")
EXAMPLE_PAIRED_BETA = _read_example("example_input_paired_beta.fasta")
EXAMPLE_PAIRED_ALPHA= _read_example("example_input_paired_alpha.fasta")

_CDR3B_ALIASES = ["cdr3b","cdr3_beta","cdr3_b","junction_aa","cdr3","CDR3","TRB_CDR3"]
_CDR3A_ALIASES = ["cdr3a","cdr3_alpha","cdr3_a","TRA_CDR3","cdr3_TRA","junction_aa_alpha"]
_NAME_ALIASES  = ["name","id","cell_id","barcode","clone_id","sample_id"]
_SOURCE_LABELS = {
    "iedb":"IEDB","vdjdb":"VDJdb","mcpas":"McPAS","tenx":"10xDcode",
    "mixtcrpred":"MixTCRpred","batcave":"BATCAVE","neotcr":"NeoTCR","cedar":"CEDAR",
}
DB_DISPLAY = {
    "iedb":"IEDB","vdjdb":"VDJdb","mcpas":"McPAS-TCR","tenx":"10x Genomics",
    "mixtcrpred":"MixTCRpred","batcave":"BATCAVE","neotcr":"NeoTCR","cedar":"CEDAR",
}
SOURCE_COLORS = {
    "IEDB":"#1e6ab0","VDJdb":"#3db5b0","McPAS":"#8b5cf6",
    "10xDcode":"#f0912a","MixTCRpred":"#d946a8","BATCAVE":"#e84545",
    "NeoTCR":"#2e9e5b","CEDAR":"#c0772a",
}

logo_b64 = _b64("webpage_design/btc_logo.png")

# ── Global CSS ─────────────────────────────────────────────────────────────────

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');

/* ── Strip Streamlit chrome ── */
#MainMenu, header[data-testid="stHeader"], footer { display: none !important; }
[data-testid="stSidebar"] { display: none !important; }
[data-testid="collapsedControl"] { display: none !important; }

/* ── Base ── */
html, body, .stApp {
    font-family: 'Inter', sans-serif;
    background: #ffffff;
    color: #1a1a2e;
}
.block-container {
    max-width: 100% !important;
    padding: 0 !important;
    margin: 0 !important;
}
section[data-testid="stVerticalBlock"] > div { gap: 0 !important; }

/* ── Nav ── */
.pm-nav {
    position: sticky; top: 0; z-index: 100;
    background: #ffffff;
    border-bottom: 1px solid #e8edf2;
    padding: 0 3rem;
    display: flex; align-items: center; justify-content: space-between;
    height: 76px;
    box-shadow: 0 1px 6px rgba(0,0,0,0.06);
}
.pm-nav-brand { display: flex; align-items: center; gap: 14px; text-decoration: none; }
.pm-nav-brand .pm-brand-text {
    font-size: 1.65rem; font-weight: 700; letter-spacing: -0.5px;
    font-style: italic; line-height: 1;
}
.pm-brand-public { color: #1e6ab0; }
.pm-brand-match  { color: #4a5568; }
.pm-nav-links { display: flex; gap: 2rem; }
.pm-nav-links a {
    color: #4a5568; font-size: 0.9rem; font-weight: 500;
    text-decoration: none; transition: color 0.2s;
}
.pm-nav-links a:hover { color: #1e6ab0; }

/* ── Hero ── */
.pm-hero {
    background: linear-gradient(160deg, #f0f4ff 0%, #ffffff 55%, #f5fffe 100%);
    padding: 4rem 3rem 3rem;
    text-align: center;
    border-bottom: 1px solid #e8edf2;
}
.pm-hero h1 {
    font-size: 2.6rem; font-weight: 700;
    color: #0b1f3a; margin: 0 0 0.6rem;
    letter-spacing: -0.8px; line-height: 1.2;
}
.pm-hero p {
    font-size: 1.05rem; color: #5a6a7e;
    max-width: 620px; margin: 0 auto 2rem;
    line-height: 1.6;
}
.pm-dots { display: flex; gap: 7px; justify-content: center; margin-bottom: 2rem; }
.pm-dot { width: 10px; height: 10px; border-radius: 3px; }

/* ── Search card ── */
.pm-search-card {
    background: #ffffff;
    border: 1.5px solid #d0dae8;
    border-radius: 16px;
    box-shadow: 0 4px 24px rgba(30,106,176,0.10);
    padding: 2rem 2.5rem 1.5rem;
    max-width: 860px;
    margin: 0 auto;
}

/* ── Chain pills ── */
.pm-chain-row { display: flex; gap: 8px; margin-bottom: 1.2rem; justify-content: center; }
.pm-chain-pill {
    padding: 6px 20px; border-radius: 20px; font-size: 0.85rem; font-weight: 600;
    border: 1.5px solid #d0dae8; cursor: pointer; transition: all 0.15s;
    background: #f7faff; color: #4a5568; user-select: none;
}
.pm-chain-pill.active {
    background: #1e6ab0; color: white; border-color: #1e6ab0;
}
.pm-chain-pill:hover:not(.active) { border-color: #1e6ab0; color: #1e6ab0; }

/* ── DB chips ── */
.pm-db-row { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 1rem; }
.pm-db-chip {
    padding: 4px 14px; border-radius: 20px; font-size: 0.8rem; font-weight: 600;
    border: 1.5px solid; cursor: pointer; transition: all 0.15s; user-select: none;
}

/* ── Streamlit widget overrides ── */
.stTextArea textarea {
    border: none !important; box-shadow: none !important;
    font-family: 'Courier New', monospace !important;
    font-size: 0.9rem !important; resize: vertical !important;
    background: #f8fafc !important;
    border-radius: 8px !important;
    padding: 12px !important;
}
.stTextArea textarea:focus {
    box-shadow: 0 0 0 2px rgba(30,106,176,0.2) !important;
}
.stTextArea label { font-weight: 600 !important; color: #1a2f4e !important; font-size: 0.85rem !important; }

/* ── Search button ── */
.stButton > button[kind="primary"] {
    background: #1e6ab0 !important; color: white !important;
    border: none !important; border-radius: 8px !important;
    font-weight: 600 !important; font-size: 1rem !important;
    padding: 0.65rem 3rem !important; width: 100% !important;
    letter-spacing: 0.2px !important; transition: background 0.2s !important;
    margin-top: 1rem !important;
}
.stButton > button[kind="primary"]:hover { background: #155a9a !important; }
.stButton > button[kind="primary"]:disabled { background: #a8c4de !important; }

/* Secondary buttons */
.stButton > button:not([kind="primary"]) {
    background: transparent !important; color: #1e6ab0 !important;
    border: 1.5px solid #1e6ab0 !important; border-radius: 8px !important;
    font-weight: 500 !important; transition: all 0.2s !important;
}
.stButton > button:not([kind="primary"]):hover {
    background: #1e6ab0 !important; color: white !important;
}

/* ── Selectbox / slider overrides ── */
.stSelectbox label, .stSlider label { font-size: 0.82rem !important; color: #5a6a7e !important; font-weight: 500 !important; }
[data-testid="stSlider"] { padding-top: 0 !important; }

/* ── Status banner ── */
.pm-status {
    background: #f0faf5; border-left: 4px solid #2e9e5b;
    border-radius: 0 8px 8px 0; padding: 0.7rem 1.2rem;
    font-size: 0.88rem; color: #1a4d2e;
    margin: 1.5rem 3rem 0;
}

/* ── Results section ── */
.pm-results-section { padding: 2rem 3rem; }
.pm-results-header { font-size: 1.3rem; font-weight: 700; color: #0b1f3a; margin-bottom: 1.2rem; }

/* ── Metric cards ── */
[data-testid="stMetric"] {
    background: #f8faff !important; border: 1px solid #dde8f5 !important;
    border-top: 3px solid #1e6ab0 !important; border-radius: 10px !important;
    padding: 1rem 1.25rem !important;
}
[data-testid="stMetricLabel"] { color: #5a7a9a !important; font-size: 0.75rem !important; text-transform: uppercase; letter-spacing: 0.6px; }
[data-testid="stMetricValue"] { color: #0b1f3a !important; font-weight: 700 !important; font-size: 1.8rem !important; }

/* ── Source DB badges in table ── */
.db-badge {
    display: inline-block; padding: 2px 10px; border-radius: 12px;
    font-size: 0.75rem; font-weight: 600; color: white;
}

/* ── Divider ── */
hr { border-color: #e8edf2 !important; margin: 0 !important; }

/* ── Footer ── */
.pm-footer {
    background: #0b1f3a; color: #8aabc8;
    padding: 2rem 3rem; margin-top: 3rem;
    display: flex; justify-content: space-between; align-items: center;
    font-size: 0.85rem;
}
.pm-footer a { color: #a8c4e0; text-decoration: none; }
.pm-footer a:hover { color: white; }
.pm-footer-brand { font-weight: 700; color: white; font-size: 1rem; }
</style>
""", unsafe_allow_html=True)

# ── Navigation ─────────────────────────────────────────────────────────────────

logo_html = (
    f'<img src="data:image/png;base64,{logo_b64}" height="52" style="display:block;">'
    if logo_b64 else '🧬'
)
st.markdown(f"""
<div class="pm-nav">
    <div class="pm-nav-brand">
        {logo_html}
        <span class="pm-brand-text">
            <span class="pm-brand-public">Public</span><span class="pm-brand-match">.Match</span>
        </span>
    </div>
    <div class="pm-nav-links">
        <a href="#">About</a>
        <a href="#">Databases</a>
        <a href="https://github.com/Marcus-Mendes/Public.Match" target="_blank">GitHub</a>
    </div>
</div>
""", unsafe_allow_html=True)

# ── Hero ───────────────────────────────────────────────────────────────────────

st.markdown("""
<div class="pm-hero">
    <h1>Search Public TCR Databases</h1>
    <p>
        Match patient CDR3 sequences against IEDB, VDJdb, McPAS-TCR,
        10x Genomics, MixTCRpred, BATCAVE, and NeoTCR in a single query.
    </p>
    <div class="pm-dots">
        <div class="pm-dot" style="background:#f0912a"></div>
        <div class="pm-dot" style="background:#3db5b0"></div>
        <div class="pm-dot" style="background:#8b5cf6"></div>
        <div class="pm-dot" style="background:#d946a8"></div>
        <div class="pm-dot" style="background:#1e6ab0"></div>
        <div class="pm-dot" style="background:#2e9e5b"></div>
        <div class="pm-dot" style="background:#e84545"></div>
    </div>
</div>
""", unsafe_allow_html=True)

# ── Session state defaults ─────────────────────────────────────────────────────

if "chain" not in st.session_state:
    st.session_state["chain"] = "beta"
if "selected_dbs" not in st.session_state:
    st.session_state["selected_dbs"] = list(ALL_DBS)
if "queries" not in st.session_state:
    st.session_state["queries"] = {}
if "_last_chain_v2" not in st.session_state:
    st.session_state["_last_chain_v2"] = "beta"

if st.session_state["_last_chain_v2"] != st.session_state["chain"]:
    st.session_state["queries"] = {}
    st.session_state["_last_chain_v2"] = st.session_state["chain"]

chain = st.session_state["chain"]

# ── Search card ────────────────────────────────────────────────────────────────

_, card_col, _ = st.columns([1, 10, 1])

with card_col:
    st.markdown('<div class="pm-search-card">', unsafe_allow_html=True)

    # Chain mode buttons
    c_b, c_a, c_p = st.columns(3)
    with c_b:
        if st.button("CDR3β only", key="chain_beta",
                     type="primary" if chain == "beta" else "secondary"):
            st.session_state["chain"] = "beta"
            st.rerun()
    with c_a:
        if st.button("CDR3α only", key="chain_alpha",
                     type="primary" if chain == "alpha" else "secondary"):
            st.session_state["chain"] = "alpha"
            st.rerun()
    with c_p:
        if st.button("Paired α + β", key="chain_paired",
                     type="primary" if chain == "paired" else "secondary"):
            st.session_state["chain"] = "paired"
            st.rerun()

    st.divider()

    # Sequence input
    if chain == "beta":
        label = "Enter CDR3β sequences (FASTA format or one per line)"
        placeholder = ">cell_001\nCASSLAPGATNEKLFF\n>cell_002\nCASSLGQTNEKLFF"
        example = EXAMPLE_BETA
    elif chain == "alpha":
        label = "Enter CDR3α sequences (FASTA format or one per line)"
        placeholder = ">cell_001\nCAGPTTSGTYKYIF\n>cell_002\nCVVSASYNTDKLIF"
        example = EXAMPLE_ALPHA
    else:
        label = "Enter CDR3β sequences (FASTA or one per line) — alpha input below"
        placeholder = ">cell_001\nCASSLVSPSEQFF\n>cell_002\nCASSLFTGTNEQFF"
        example = EXAMPLE_PAIRED_BETA

    raw_b = st.text_area(label, height=160, placeholder=placeholder, key="seq_input_b")

    if chain == "paired":
        raw_a = st.text_area(
            "Enter CDR3α sequences (matched by name to CDR3β above)",
            height=120,
            placeholder=">cell_001\nCAGPTTSGTYKYIF\n>cell_002\nCVVSASYNTDKLIF",
            key="seq_input_a",
        )
    else:
        raw_a = ""

    # Process text input into queries immediately
    if chain == "paired" and raw_b.strip() and raw_a.strip():
        beta_seqs = parse_fasta(raw_b)
        alpha_seqs = parse_fasta(raw_a)
        common = set(beta_seqs) & set(alpha_seqs)
        if common:
            st.session_state["queries"] = {
                name: (alpha_seqs[name], beta_seqs[name]) for name in sorted(common)
            }
    elif chain != "paired" and raw_b.strip():
        st.session_state["queries"] = parse_fasta(raw_b)

    # File upload + example row
    up_col, ex_col = st.columns([2, 1])
    with up_col:
        uploaded = st.file_uploader(
            "or upload FASTA" if chain != "paired" else "or upload CDR3β FASTA",
            type=["fasta","fa","txt"], label_visibility="collapsed",
            key="v2_upload_b",
        )
        if uploaded:
            raw = uploaded.read().decode("utf-8")
            st.session_state["queries"] = parse_fasta(raw) if chain != "paired" else st.session_state["queries"]
    with ex_col:
        if st.button("Load example →", key="v2_example"):
            if chain == "paired":
                b = parse_fasta(EXAMPLE_PAIRED_BETA)
                a = parse_fasta(EXAMPLE_PAIRED_ALPHA)
                common = set(b) & set(a)
                st.session_state["queries"] = {name: (a[name], b[name]) for name in sorted(common)}
            elif chain == "alpha":
                st.session_state["queries"] = parse_fasta(EXAMPLE_ALPHA)
            else:
                st.session_state["queries"] = parse_fasta(EXAMPLE_BETA)
            st.rerun()

    queries = st.session_state.get("queries", {})
    if queries:
        n = len(queries)
        chain_label = "paired" if chain == "paired" else ("CDR3β" if chain == "beta" else "CDR3α")
        st.success(f"✓ {n} {chain_label} sequence{'s' if n > 1 else ''} ready")

    st.divider()

    # Database selection
    st.markdown("<p style='font-size:0.85rem;font-weight:600;color:#1a2f4e;margin-bottom:8px;'>Databases</p>",
                unsafe_allow_html=True)
    db_cols = st.columns(len(ALL_DBS))
    selected_dbs = []
    for i, db in enumerate(ALL_DBS):
        with db_cols[i]:
            checked = db in st.session_state["selected_dbs"]
            if st.checkbox(DB_DISPLAY[db], value=checked, key=f"v2_db_{db}"):
                selected_dbs.append(db)
    st.session_state["selected_dbs"] = selected_dbs

    st.divider()

    # Method + threshold
    m_col, t_col = st.columns([1, 2])
    with m_col:
        method = st.selectbox("Method", ["blosum", "exact", "edit"],
                              format_func=lambda x: {"blosum":"BLOSUM62","exact":"Exact","edit":"Edit distance"}[x],
                              key="v2_method")
    with t_col:
        if method == "blosum":
            threshold = st.slider("Min score", 0.80, 1.00, 0.97, 0.01, key="v2_thresh")
        elif method == "edit":
            threshold = float(st.slider("Max edit distance", 0, 5, 1, key="v2_edit"))
        else:
            threshold = 1.0
            st.markdown("<br>", unsafe_allow_html=True)

    run_disabled = not queries or not selected_dbs
    st.button("Search Public Databases", type="primary", disabled=run_disabled, key="v2_run")

    st.markdown('</div>', unsafe_allow_html=True)

# ── Database loading ───────────────────────────────────────────────────────────

@st.cache_resource
def _get_beta_ref():  return load_databases_cached(ALL_DBS, chain="beta")
@st.cache_resource
def _get_alpha_ref(): return load_databases_cached(ALL_DBS, chain="alpha")
@st.cache_resource
def _get_paired_ref():return load_databases_cached(ALL_DBS, chain="paired")

_loaders = {"beta": _get_beta_ref, "alpha": _get_alpha_ref, "paired": _get_paired_ref}

with st.spinner("Loading reference databases…"):
    all_ref = _loaders[chain]()

db_counts = all_ref["source_db"].value_counts()
badges = "".join(
    f'<span style="background:{SOURCE_COLORS.get(k,"#888")};color:white;'
    f'padding:2px 10px;border-radius:12px;font-size:0.78rem;font-weight:600;">'
    f'{k}: {v:,}</span> '
    for k, v in db_counts.items()
)
st.markdown(
    f'<div class="pm-status">✅ &nbsp;<strong>{len(all_ref):,} reference entries</strong> '
    f'loaded &nbsp;—&nbsp; {badges}</div>',
    unsafe_allow_html=True,
)

# Filter to selected dbs
if set(selected_dbs) != set(ALL_DBS):
    labels = [_SOURCE_LABELS[d] for d in selected_dbs]
    reference = all_ref[all_ref["source_db"].isin(labels)].reset_index(drop=True)
else:
    reference = all_ref

# ── Run matching ───────────────────────────────────────────────────────────────

if st.session_state.get("v2_run") and queries:
    query_list = list(queries.values())
    with st.spinner(f"Searching {len(query_list)} sequence(s) across {len(reference):,} entries…"):
        results = match(
            queries=query_list,
            reference=reference,
            method=method,
            threshold=threshold,
            chain=chain,
        )

    if not results.empty:
        seq_to_name = {v: k for k, v in queries.items()}
        if chain == "alpha":
            results.insert(0, "query_name", results["query_cdr3a"].map(seq_to_name))
        elif chain == "paired":
            t2n = {v: k for k, v in queries.items()}
            results.insert(0, "query_name",
                pd.Series(zip(results["query_cdr3a"], results["query_cdr3b"])).map(
                    lambda t: t2n.get(t, str(t))
                ).values)
        else:
            results.insert(0, "query_name", results["query_cdr3b"].map(seq_to_name))

    # ── Results display ──────────────────────────────────────────────────────

    st.markdown('<div class="pm-results-section">', unsafe_allow_html=True)
    st.markdown('<div class="pm-results-header">Results</div>', unsafe_allow_html=True)

    if results.empty:
        st.warning("No matches found. Try lowering the threshold or selecting more databases.")
    else:
        query_col = "query_cdr3a" if chain == "alpha" else "query_cdr3b"

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Total Matches",    f"{len(results):,}")
        c2.metric("Queries with Hits", results[query_col].nunique())
        c3.metric("Unique Epitopes",   results["epitope"].nunique())
        c4.metric("Databases Hit",     results["source_db"].nunique())

        st.markdown("#### Hits by database")
        st.bar_chart(
            results["source_db"].value_counts().rename_axis("Database").reset_index(name="Hits").set_index("Database"),
            color="#1e6ab0",
        )

        # Style the source_db column as colored badges
        def _badge(source):
            col = SOURCE_COLORS.get(source, "#888")
            return f'<span class="db-badge" style="background:{col}">{source}</span>'

        display_df = results.copy()
        display_df["source_db"] = display_df["source_db"].apply(_badge)

        st.markdown("#### All matches")
        st.markdown(
            display_df.to_html(escape=False, index=False,
                               classes="pm-table",
                               border=0),
            unsafe_allow_html=True,
        )
        st.markdown("""<style>
        .pm-table { width:100%; border-collapse:collapse; font-size:0.85rem; margin-top:0.5rem; }
        .pm-table th { background:#f0f4ff; color:#1a2f4e; font-weight:600; padding:8px 12px;
                       border-bottom:2px solid #d0dae8; text-align:left; }
        .pm-table td { padding:7px 12px; border-bottom:1px solid #eef1f6; vertical-align:middle; }
        .pm-table tr:hover td { background:#f7faff; }
        </style>""", unsafe_allow_html=True)

        st.markdown("<br>", unsafe_allow_html=True)
        st.download_button(
            "⬇  Download results CSV",
            data=results.to_csv(index=False),
            file_name="public_match_results.csv",
            mime="text/csv",
        )

    st.markdown('</div>', unsafe_allow_html=True)

# ── Footer ─────────────────────────────────────────────────────────────────────

st.markdown("""
<div class="pm-footer">
    <div>
        <div class="pm-footer-brand">🧬 Public.Match</div>
        <div style="margin-top:4px;">BTC Hackathon 2026 · Karchin Lab · Johns Hopkins</div>
    </div>
    <div style="text-align:right;">
        <a href="https://github.com/Marcus-Mendes/Public.Match" target="_blank">GitHub</a>
        &nbsp;·&nbsp;
        <a href="https://breakthroughcancer.org" target="_blank">Break Through Cancer</a>
    </div>
</div>
""", unsafe_allow_html=True)

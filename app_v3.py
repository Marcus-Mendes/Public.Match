"""
Public.Match v3 — CDR3 Search + Epitope Search
Run from repo root: streamlit run app_v3.py --server.port 8503
"""

import base64
import tempfile
import streamlit as st
import pandas as pd
from pathlib import Path

from public_match.database import load_databases_cached, ALL_DBS
from public_match.matcher import match
from public_match.epitope_matcher import match_by_epitope
from public_match.parsers import custom as custom_parser

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

def parse_epitopes(text: str) -> list[str]:
    """One epitope per line; skip blanks and comment lines."""
    return [
        line.strip().upper()
        for line in text.strip().splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]

EXAMPLE_BETA        = _read_example("example_input.fasta")
EXAMPLE_ALPHA       = _read_example("example_input_alpha.fasta")
EXAMPLE_PAIRED_BETA = _read_example("example_input_paired_beta.fasta")
EXAMPLE_PAIRED_ALPHA= _read_example("example_input_paired_alpha.fasta")
EXAMPLE_EPITOPES    = "GILGFVFTL\nNLVPMVATAV\nKLGGALQAK"  # flu, CMV, EBV

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
.pm-hero-brand {
    font-size: 3.4rem; font-weight: 800; font-style: italic;
    letter-spacing: -1.2px; line-height: 1; margin: 0 0 0.5rem;
}
.pm-hero-subtitle {
    font-size: 1.5rem; font-weight: 600;
    color: #0b1f3a; margin: 0 0 0.8rem;
    letter-spacing: -0.4px; line-height: 1.3;
}
.pm-hero p {
    font-size: 1.05rem; color: #5a6a7e;
    max-width: 620px; margin: 0 auto 2rem;
    line-height: 1.6;
}
.pm-dots { display: flex; gap: 7px; justify-content: center; margin-bottom: 2rem; }
.pm-dot { width: 10px; height: 10px; border-radius: 3px; }

/* ── Search card — style the middle Streamlit column as a card ── */
[data-testid="stHorizontalBlock"]:first-of-type
    > [data-testid="stColumn"]:nth-child(2)
    > div:first-child {
    background: #ffffff;
    border: 1.5px solid #d0dae8;
    border-radius: 16px;
    box-shadow: 0 4px 24px rgba(30,106,176,0.10);
    padding: 2rem 2.5rem 1.5rem !important;
}

/* ── Mode toggle label ── */
.pm-mode-label {
    font-size: 0.75rem; font-weight: 700; color: #8a9ab5;
    text-transform: uppercase; letter-spacing: 0.8px;
    margin-bottom: 6px;
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
        <a href="." >Home</a>
        <a href="?page=about">About</a>
        <a href="?page=databases">Databases</a>
        <a href="https://github.com/Marcus-Mendes/Public.Match" target="_blank">GitHub</a>
    </div>
</div>
""", unsafe_allow_html=True)

# ── Page routing ───────────────────────────────────────────────────────────────

_page = st.query_params.get("page", "home")

if _page == "about":
    readme = Path("README.md").read_text() if Path("README.md").exists() else "_README.md not found._"
    st.markdown('<div style="max-width:860px;margin:3rem auto;padding:0 2rem 4rem;">', unsafe_allow_html=True)
    st.markdown(readme)
    st.markdown('</div>', unsafe_allow_html=True)
    st.stop()

if _page == "databases":
    dbmd = Path("DATABASES.md").read_text() if Path("DATABASES.md").exists() else "_DATABASES.md not found._"
    st.markdown('<div style="max-width:860px;margin:3rem auto;padding:0 2rem 4rem;">', unsafe_allow_html=True)
    st.markdown(dbmd)
    st.markdown('</div>', unsafe_allow_html=True)
    st.stop()

# ── Hero ───────────────────────────────────────────────────────────────────────

st.markdown("""
<div class="pm-hero">
    <h1 class="pm-hero-brand">
        <span class="pm-brand-public">Public</span><span class="pm-brand-match">.Match</span>
    </h1>
    <h2 class="pm-hero-subtitle">Search Public TCR Databases</h2>
    <p>
        Match CDR3 sequences against public databases — or search in reverse
        by epitope to find all known T cells that recognise it.
    </p>
    <div class="pm-dots">
        <div class="pm-dot" style="background:#f0912a"></div>
        <div class="pm-dot" style="background:#3db5b0"></div>
        <div class="pm-dot" style="background:#8b5cf6"></div>
        <div class="pm-dot" style="background:#d946a8"></div>
        <div class="pm-dot" style="background:#1e6ab0"></div>
        <div class="pm-dot" style="background:#2e9e5b"></div>
        <div class="pm-dot" style="background:#e84545"></div>
        <div class="pm-dot" style="background:#c0772a"></div>
    </div>
</div>
""", unsafe_allow_html=True)

# ── Session state defaults ─────────────────────────────────────────────────────

for key, default in [
    ("search_mode", "cdr3"),
    ("chain", "beta"),
    ("selected_dbs", list(ALL_DBS)),
    ("queries", {}),
    ("epitope_queries", []),
    ("_last_chain_v3", "beta"),
]:
    if key not in st.session_state:
        st.session_state[key] = default

if st.session_state["_last_chain_v3"] != st.session_state["chain"]:
    st.session_state["queries"] = {}
    st.session_state["_last_chain_v3"] = st.session_state["chain"]

search_mode = st.session_state["search_mode"]
chain = st.session_state["chain"]

# ── Search card ────────────────────────────────────────────────────────────────

_, card_col, _ = st.columns([1, 10, 1])

with card_col:

    # ── Mode toggle (CDR3 Search vs Epitope Search) ──────────────────────────
    st.markdown('<div class="pm-mode-label">Search mode</div>', unsafe_allow_html=True)
    m_cdr3, m_epi = st.columns(2)
    with m_cdr3:
        if st.button("CDR3 Search", key="mode_cdr3",
                     type="primary" if search_mode == "cdr3" else "secondary"):
            st.session_state["search_mode"] = "cdr3"
            st.rerun()
    with m_epi:
        if st.button("Epitope Search", key="mode_epi",
                     type="primary" if search_mode == "epitope" else "secondary"):
            st.session_state["search_mode"] = "epitope"
            st.rerun()

    st.divider()

    # ── CDR3 Search mode ─────────────────────────────────────────────────────
    if search_mode == "cdr3":

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

        up_col, ex_col = st.columns([2, 1])
        with up_col:
            uploaded = st.file_uploader(
                "or upload FASTA" if chain != "paired" else "or upload CDR3β FASTA",
                type=["fasta","fa","txt"], label_visibility="collapsed",
                key="v3_upload_b",
            )
            if uploaded:
                raw = uploaded.read().decode("utf-8")
                st.session_state["queries"] = parse_fasta(raw) if chain != "paired" else st.session_state["queries"]
        with ex_col:
            if st.button("Load example →", key="v3_example"):
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

    # ── Epitope Search mode ──────────────────────────────────────────────────
    else:
        st.markdown(
            "<p style='font-size:0.9rem;color:#5a6a7e;margin-bottom:4px;'>"
            "Enter epitope peptide sequences (one per line) to find all CDR3s that recognise them.</p>",
            unsafe_allow_html=True,
        )
        raw_epi = st.text_area(
            "Epitope sequences (one per line)",
            height=160,
            placeholder="GILGFVFTL\nNLVPMVATAV\nKLGGALQAK",
            key="epi_input",
        )
        if raw_epi.strip():
            st.session_state["epitope_queries"] = parse_epitopes(raw_epi)

        up_col, ex_col = st.columns([2, 1])
        with up_col:
            uploaded_epi = st.file_uploader(
                "or upload epitope list (.txt)",
                type=["txt","csv"], label_visibility="collapsed",
                key="v3_upload_epi",
            )
            if uploaded_epi:
                st.session_state["epitope_queries"] = parse_epitopes(
                    uploaded_epi.read().decode("utf-8")
                )
        with ex_col:
            if st.button("Load example →", key="v3_epi_example"):
                st.session_state["epitope_queries"] = parse_epitopes(EXAMPLE_EPITOPES)
                st.rerun()

        epitope_queries = st.session_state.get("epitope_queries", [])
        if epitope_queries:
            st.success(f"✓ {len(epitope_queries)} epitope{'s' if len(epitope_queries) > 1 else ''} ready")

    st.divider()

    # Database selection (shared by both modes)
    st.markdown("<p style='font-size:0.85rem;font-weight:600;color:#1a2f4e;margin-bottom:8px;'>Databases</p>",
                unsafe_allow_html=True)
    db_cols = st.columns(len(ALL_DBS))
    selected_dbs = []
    for i, db in enumerate(ALL_DBS):
        with db_cols[i]:
            checked = db in st.session_state["selected_dbs"]
            if st.checkbox(DB_DISPLAY[db], value=checked, key=f"v3_db_{db}"):
                selected_dbs.append(db)
    st.session_state["selected_dbs"] = selected_dbs

    # Custom database (optional)
    with st.expander("＋ Add a custom database (optional)"):
        st.markdown(
            "<p style='font-size:0.85rem;color:#5a6a7e;margin:0 0 8px;'>"
            "Upload a CSV or TSV file. CDR3β/α columns are auto-detected from common names "
            "(cdr3b, cdr3_beta, junction_aa, cdr3a, etc.). "
            "Use the hint below only if auto-detection fails.</p>",
            unsafe_allow_html=True,
        )
        custom_uploads = st.file_uploader(
            "Custom database file(s)",
            type=["csv", "tsv", "txt"],
            accept_multiple_files=True,
            key="v3_custom_db",
            label_visibility="collapsed",
        )
        hint_col, only_col = st.columns([2, 1])
        with hint_col:
            custom_cdr3_hint = st.text_input(
                "CDR3β column name hint",
                placeholder="e.g. my_cdr3_column  (leave blank for auto)",
                key="v3_custom_cdr3_col",
                label_visibility="collapsed",
            )
        with only_col:
            custom_only = st.checkbox(
                "Custom database only",
                help="Skip all built-in databases and search only your uploaded file(s).",
                key="v3_custom_only",
            )

    st.divider()

    # Method + threshold (shared by both modes)
    m_col, t_col = st.columns([1, 2])
    with m_col:
        method = st.selectbox("Method", ["blosum", "exact", "edit"],
                              format_func=lambda x: {"blosum":"BLOSUM62","exact":"Exact","edit":"Edit distance"}[x],
                              key="v3_method")
    with t_col:
        if method == "blosum":
            threshold = st.slider("Min score", 0.80, 1.00, 0.97, 0.01, key="v3_thresh")
        elif method == "edit":
            threshold = float(st.slider("Max edit distance", 0, 5, 1, key="v3_edit"))
        else:
            threshold = 1.0
            st.markdown("<br>", unsafe_allow_html=True)

    custom_uploads = st.session_state.get("v3_custom_db") or []
    custom_only    = st.session_state.get("v3_custom_only", False)
    has_db         = bool(selected_dbs) or bool(custom_uploads)

    if search_mode == "cdr3":
        queries = st.session_state.get("queries", {})
        run_disabled = not queries or not has_db
    else:
        epitope_queries = st.session_state.get("epitope_queries", [])
        run_disabled = not epitope_queries or not has_db

    st.button("Search Public Databases", type="primary", disabled=run_disabled, key="v3_run")

# ── Database loading ───────────────────────────────────────────────────────────

@st.cache_resource
def _get_beta_ref():   return load_databases_cached(ALL_DBS, chain="beta")
@st.cache_resource
def _get_alpha_ref():  return load_databases_cached(ALL_DBS, chain="alpha")
@st.cache_resource
def _get_paired_ref(): return load_databases_cached(ALL_DBS, chain="paired")

_loaders = {"beta": _get_beta_ref, "alpha": _get_alpha_ref, "paired": _get_paired_ref}

# Epitope search always uses the beta reference (broadest coverage)
_ref_chain = chain if search_mode == "cdr3" else "beta"

with st.spinner("Loading reference databases…"):
    all_ref = _loaders[_ref_chain]()

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

# Build reference from built-in DBs (skip if custom-only)
custom_uploads  = st.session_state.get("v3_custom_db") or []
custom_only     = st.session_state.get("v3_custom_only", False)
custom_cdr3_hint = st.session_state.get("v3_custom_cdr3_col", "")

if custom_only:
    reference = pd.DataFrame(columns=all_ref.columns)
elif set(selected_dbs) != set(ALL_DBS):
    labels = [_SOURCE_LABELS[d] for d in selected_dbs]
    reference = all_ref[all_ref["source_db"].isin(labels)].reset_index(drop=True)
else:
    reference = all_ref

# Load and append custom database files
_ref_chain_for_custom = chain if search_mode == "cdr3" else "beta"
custom_load_errors = []
if custom_uploads:
    custom_frames = []
    for uf in custom_uploads:
        try:
            cdf = _load_custom_upload(uf, cdr3_hint=custom_cdr3_hint)
            cdf = _apply_chain_filter(cdf, _ref_chain_for_custom)
            custom_frames.append(cdf)
        except Exception as exc:
            custom_load_errors.append(f"{uf.name}: {exc}")
    if custom_frames:
        custom_combined = pd.concat(custom_frames, ignore_index=True)
        reference = pd.concat([reference, custom_combined], ignore_index=True)

for err in custom_load_errors:
    st.warning(f"⚠️ Could not load custom DB — {err}")

# ── Run matching ───────────────────────────────────────────────────────────────

def _badge(source):
    col = SOURCE_COLORS.get(source, "#888")
    return f'<span class="db-badge" style="background:{col}">{source}</span>'


def _make_aggregate(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    """
    Collapse full results to unique CDR3/epitope pairs.
    Adds db_count (number of distinct databases) and databases (comma-separated list).
    group_cols should not include source_db or score.
    """
    valid_cols = [c for c in group_cols if c in df.columns]
    return (
        df.groupby(valid_cols, dropna=False)
        .agg(
            db_count=("source_db", "nunique"),
            databases=("source_db", lambda x: ", ".join(sorted(x.dropna().unique()))),
            max_score=("score", "max"),
        )
        .reset_index()
        .sort_values("db_count", ascending=False)
        .reset_index(drop=True)
    )


_AA_PAT = r"^[ACDEFGHIKLMNPQRSTVWY]+$"

def _load_custom_upload(uploaded_file, cdr3_hint: str = "") -> pd.DataFrame:
    """Write an uploaded file to a temp path and parse it via custom_parser.load()."""
    suffix = Path(uploaded_file.name).suffix.lower() or ".csv"
    with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp:
        tmp.write(uploaded_file.read())
        tmp_path = Path(tmp.name)
    try:
        return custom_parser.load(
            tmp_path,
            cdr3_col=cdr3_hint.strip() or None,
            source_name=Path(uploaded_file.name).stem,
        )
    finally:
        tmp_path.unlink(missing_ok=True)


def _apply_chain_filter(df: pd.DataFrame, chain: str) -> pd.DataFrame:
    """Apply the same CDR3 chain filter used by load_databases()."""
    if chain in ("beta", "paired"):
        df = df[df["cdr3b"].str.match(_AA_PAT, na=False) & (df["cdr3b"].str.len() >= 8)]
    if chain in ("alpha", "paired"):
        df = df[df["cdr3a"].notna() & (df["cdr3a"].str.len() >= 8)]
    return df.reset_index(drop=True)


TABLE_CSS = """<style>
.pm-table { width:100%; border-collapse:collapse; font-size:0.85rem; margin-top:0.5rem; }
.pm-table th { background:#f0f4ff; color:#1a2f4e; font-weight:600; padding:8px 12px;
               border-bottom:2px solid #d0dae8; text-align:left; }
.pm-table td { padding:7px 12px; border-bottom:1px solid #eef1f6; vertical-align:middle; }
.pm-table tr:hover td { background:#f7faff; }
</style>"""

if st.session_state.get("v3_run"):

    # ── CDR3 Search ───────────────────────────────────────────────────────────
    if search_mode == "cdr3":
        queries = st.session_state.get("queries", {})
        if queries:
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

            st.markdown('<div class="pm-results-section">', unsafe_allow_html=True)
            st.markdown('<div class="pm-results-header">CDR3 Search Results</div>', unsafe_allow_html=True)

            if results.empty:
                st.warning("No matches found. Try lowering the threshold or selecting more databases.")
            else:
                query_col = "query_cdr3a" if chain == "alpha" else "query_cdr3b"
                c1, c2, c3, c4 = st.columns(4)
                c1.metric("Total Matches",     f"{len(results):,}")
                c2.metric("Queries with Hits", results[query_col].nunique())
                c3.metric("Unique Epitopes",   results["epitope"].nunique())
                c4.metric("Databases Hit",     results["source_db"].nunique())

                st.markdown("#### Hits by database")
                st.bar_chart(
                    results["source_db"].value_counts().rename_axis("Database").reset_index(name="Hits").set_index("Database"),
                    color="#1e6ab0",
                )

                display_df = results.copy()
                display_df["source_db"] = display_df["source_db"].apply(_badge)
                st.markdown("#### All matches")
                st.markdown(display_df.to_html(escape=False, index=False, classes="pm-table", border=0), unsafe_allow_html=True)
                st.markdown(TABLE_CSS, unsafe_allow_html=True)

                st.markdown("<br>", unsafe_allow_html=True)
                dl1, dl2 = st.columns(2)
                with dl1:
                    st.download_button(
                        "⬇  Download full results CSV",
                        data=results.to_csv(index=False),
                        file_name="public_match_cdr3_results.csv",
                        mime="text/csv",
                        key="v3_dl_cdr3_full",
                    )
                with dl2:
                    if chain == "alpha":
                        agg_cols = ["query_name", "query_cdr3a", "cdr3a", "epitope"]
                    elif chain == "paired":
                        agg_cols = ["query_name", "query_cdr3a", "query_cdr3b", "cdr3a", "cdr3b", "epitope"]
                    else:
                        agg_cols = ["query_name", "query_cdr3b", "cdr3b", "epitope"]
                    agg_df = _make_aggregate(results, agg_cols)
                    st.download_button(
                        "⬇  Download aggregate CSV",
                        data=agg_df.to_csv(index=False),
                        file_name="public_match_cdr3_aggregate.csv",
                        mime="text/csv",
                        key="v3_dl_cdr3_agg",
                    )

            st.markdown('</div>', unsafe_allow_html=True)

    # ── Epitope Search ────────────────────────────────────────────────────────
    else:
        epitope_queries = st.session_state.get("epitope_queries", [])
        if epitope_queries:
            with st.spinner(f"Searching {len(epitope_queries)} epitope(s) across {len(reference):,} entries…"):
                results = match_by_epitope(
                    epitopes=epitope_queries,
                    reference=reference,
                    method=method,
                    threshold=threshold,
                )

            st.markdown('<div class="pm-results-section">', unsafe_allow_html=True)
            st.markdown('<div class="pm-results-header">Epitope Search Results</div>', unsafe_allow_html=True)

            if results.empty:
                st.warning("No matches found. Try a different epitope, lower the threshold, or select more databases.")
            else:
                c1, c2, c3, c4 = st.columns(4)
                c1.metric("Total Matches",   f"{len(results):,}")
                c2.metric("Unique CDR3β",    results["cdr3b"].nunique())
                c3.metric("Unique CDR3α",    results["cdr3a"].nunique())
                c4.metric("Databases Hit",   results["source_db"].nunique())

                st.markdown("#### Hits by database")
                st.bar_chart(
                    results["source_db"].value_counts().rename_axis("Database").reset_index(name="Hits").set_index("Database"),
                    color="#2e9e5b",
                )

                display_df = results.copy()
                display_df["source_db"] = display_df["source_db"].apply(_badge)
                st.markdown("#### All matching CDR3 sequences")
                st.markdown(display_df.to_html(escape=False, index=False, classes="pm-table", border=0), unsafe_allow_html=True)
                st.markdown(TABLE_CSS, unsafe_allow_html=True)

                st.markdown("<br>", unsafe_allow_html=True)
                dl1, dl2 = st.columns(2)
                with dl1:
                    st.download_button(
                        "⬇  Download full results CSV",
                        data=results.to_csv(index=False),
                        file_name="public_match_epitope_results.csv",
                        mime="text/csv",
                        key="v3_dl_epi_full",
                    )
                with dl2:
                    agg_df = _make_aggregate(results, ["query_epitope", "epitope", "cdr3a", "cdr3b"])
                    st.download_button(
                        "⬇  Download aggregate CSV",
                        data=agg_df.to_csv(index=False),
                        file_name="public_match_epitope_aggregate.csv",
                        mime="text/csv",
                        key="v3_dl_epi_agg",
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

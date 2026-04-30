# OTS — Observed TCR Space (Oxford/OPIG)

Source: https://opig.stats.ox.ac.uk/webapps/ots  
Publication: Cell Reports 2024 — "The Observed T Cell Receptor Space database enables paired-chain repertoire mining, coherence analysis, and language modeling"  
Entries: 5.35M redundant / **1.63M non-redundant** paired αβ TCR sequences; updated quarterly

## Why OTS matters

OTS is the largest curated collection of **full-length, consistently processed, paired αβ TCR sequences** from 50+ public studies. Unlike IEDB/VDJdb it is not filtered to known-epitope TCRs — it represents the entire public repertoire, which makes it the best resource for:
- Identifying *public* TCRs (sequences shared across individuals)
- Distinguishing antigen-specific from bystander T cells by publicness scoring
- Paired αβ language model embeddings

## Download (manual — browser required)

OTS uses a JavaScript-rendered interface that cannot be fully automated. To download:

1. Go to https://opig.stats.ox.ac.uk/webapps/ots/ots_paired/
2. Submit the search form **with no filters** to select all sequences
3. Click the **"here"** link next to "A shell-script with the commands to download all the data-units can be downloaded"
4. Run the downloaded `bulk_download.sh` script — it issues one `curl` per study (~50 studies)
5. Place all downloaded files in this folder

## Automated download helper

Once you have `bulk_download.sh` from step 3 above, run:

```bash
cd Databases/OTS
bash bulk_download.sh
```

Each downloaded file is a gzipped FASTA or CSV per study. Use `ots_parse.py` (to be added) to merge into a unified CSV.

## Expected columns (per-study CSV format)

| Column | Description |
|---|---|
| `cdr3_alpha` | CDR3 alpha amino acid sequence |
| `cdr3_beta` | CDR3 beta amino acid sequence |
| `v_alpha` | TRAV gene |
| `j_alpha` | TRAJ gene |
| `v_beta` | TRBV gene |
| `j_beta` | TRBJ gene |
| `species` | Human / Mouse |
| `disease` | Disease/condition label |
| `study` | Source study ID |

## Unified schema mapping

```
cdr3_alpha  ← cdr3_alpha
cdr3_beta   ← cdr3_beta
epitope     ← (not available — OTS does not include epitope annotations)
mhc         ← (not available)
source      ← "OTS"
```

**Note:** OTS has no epitope labels. Use it for publicness scoring — if a query TCR matches an OTS sequence, it is a known public TCR, even without confirmed epitope specificity.

## Data size warning

The full OTS database is ~several GB uncompressed. Do not commit the data files to git.
Add `Databases/OTS/*.csv` and `Databases/OTS/*.gz` to `.gitignore`.

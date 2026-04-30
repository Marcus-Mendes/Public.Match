# Database Reference

This document describes the reference databases included in `Databases/` and how they are used by Public.Match.

---

## IEDB — Immune Epitope Database

**File:** `Databases/IEDB/iedb.xlsx`
**Source:** https://www.iedb.org (receptor export)
**Entries:** ~217,000
**Format:** Excel (.xlsx), manually exported from the IEDB receptor search

### Key columns

| Column | Description |
|---|---|
| `Chain 1 - CDR3 Curated` / `Chain 1 - CDR3 Calculated` | CDR3 sequence for chain 1 (TRA or TRB) |
| `Chain 1 - Type` | Chain type (Alpha / Beta) |
| `Chain 2 - CDR3 Curated` / `Chain 2 - CDR3 Calculated` | CDR3 sequence for chain 2 |
| `Epitope - Name` | Epitope peptide sequence |
| `Epitope - Source Organism` | Pathogen or tissue of origin |
| `Assay - MHC Allele Names` | HLA restriction |
| `Receptor - Type` | TCR / BCR |

### Notes
- Prefer `CDR3 Curated` over `CDR3 Calculated` when available; fall back to Calculated.
- Filter to `Receptor - Type == TCR` and chain type `Beta` for CDR3β matching.

---

## VDJdb

**File:** `Databases/VDJdb/VDJdb_05302026.zip` → `VDJdb_05302026.tsv`
**Source:** https://vdjdb.cdr3.net
**Entries:** ~210,000
**Format:** ZIP-compressed TSV

### Key columns

| Column | Description |
|---|---|
| `CDR3` | CDR3 amino acid sequence |
| `Gene` | TRA or TRB |
| `V` / `J` | V and J gene calls |
| `Species` | HomoSapiens / MusMusculus |
| `MHC A` / `MHC B` | HLA alleles |
| `MHC class` | MHCI or MHCII |
| `Epitope` | Epitope peptide |
| `Epitope gene` | Antigen gene name |
| `Epitope species` | Pathogen |
| `Score` | VDJdb confidence score (0–3; use ≥1 recommended) |

### Notes
- Filter to `Gene == TRB` for CDR3β matching.
- Recommended to filter `Score >= 1` to exclude low-confidence entries.

---

## McPAS-TCR

**File:** `Databases/McPAS/McPAS-TCR.csv`
**Source:** http://friedmanlab.weizmann.ac.il/McPAS-TCR
**Entries:** ~40,700
**Format:** CSV (latin-1 encoded)

### Key columns

| Column | Description |
|---|---|
| `CDR3.beta.aa` | CDR3β amino acid sequence |
| `CDR3.alpha.aa` | CDR3α amino acid sequence |
| `TRBV` / `TRBJ` | V and J gene calls |
| `Epitope.peptide` | Epitope peptide |
| `MHC` | HLA restriction |
| `Pathology` | Disease/pathology category |
| `Category` | Infectious disease / Cancer / Autoimmune |
| `Species` | Human / Mouse |

### Notes
- Filter `Species == Human` for human repertoire analysis.
- Many entries lack an epitope peptide; use `Pathology` for disease-level annotation.

---

## 10x Genomics Dcode (Human Donors)

**Files:** `Databases/10xDcode/vdj_v1_hs_aggregated_donor{1–4}_binarized_matrix.csv`
**Source:** 10x Genomics public dataset (4 healthy donors)
**Entries:** ~46k–78k cells per donor
**Format:** CSV

### Key columns

| Column | Description |
|---|---|
| `cell_clono_cdr3_aa` | Paired CDR3 sequences: `TRA:…;TRB:…` (semicolon-separated) |
| `cell_clono_cdr3_nt` | Nucleotide CDR3 sequences |
| `CD3`, `CD4`, `CD8a`, … | Surface marker counts (CITE-seq) |
| `{HLA}_{peptide}_{antigen}_{pathogen}` | Dextramer UMI count per epitope (e.g. `A0201_GILGFVFTL_Flu-MP_Influenza`) |
| `{HLA}_{peptide}_{antigen}_{pathogen}_binder` | Boolean — True if cell is a confirmed binder for that epitope |

### Epitope column format

```
A0201_GILGFVFTL_Flu-MP_Influenza
│     │         │      └─ Pathogen
│     │         └─ Antigen name
│     └─ Epitope peptide
└─ HLA allele
```

### Reshaped (long) format

Run `scripts/reshape_10x_dcode.py` to convert all four donor files into a single long-format table:

**Output:** `Databases/10xDcode/10xDcode_long.csv` (~18,800 rows, one per CDR3β–epitope binder pair)

| Column | Description |
|---|---|
| `cdr3b` | CDR3β amino acid sequence |
| `cdr3a` | CDR3α amino acid sequence (first TRA chain) |
| `epitope` | Epitope peptide |
| `HLA` | HLA allele (e.g. `A0201`) |
| `antigen` | Antigen name (e.g. `Flu-MP`) |
| `pathogen` | Pathogen (e.g. `Influenza`, `CMV`, `Cancer`) |
| `dextramer_count` | Raw dextramer UMI count |
| `cell_phenotype` | `CD4`, `CD8`, or `DN` (based on surface markers) |
| `barcode` | Cell barcode |
| `donor` | Donor ID (donor1–donor4) |

### Notes
- Only confirmed binder cells (`_binder == True`) are included in the long format.
- Donors 1–4 cover multiple HLA supertypes; useful for cross-HLA public TCR discovery.
- Duplicates are dropped on `(cdr3b, epitope, HLA, donor)`.

---

## Database Summary

| Database | Entries | CDR3β col | Epitope col | HLA col |
|---|---|---|---|---|
| IEDB | ~217k | `Chain N - CDR3 Curated` | `Epitope - Name` | `Assay - MHC Allele Names` |
| VDJdb | ~210k | `CDR3` (where Gene=TRB) | `Epitope` | `MHC A` |
| McPAS | ~41k | `CDR3.beta.aa` | `Epitope.peptide` | `MHC` |
| 10x Dcode | ~230k cells | parsed from `cell_clono_cdr3_aa` | column name | column name |

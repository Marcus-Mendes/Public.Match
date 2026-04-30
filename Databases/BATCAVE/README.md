# BATCAVE — Comprehensive TCR–Epitope Mutational Scan Database

Source: https://github.com/meyer-lab-cshl/BATMAN-paper  
Publication: Cell Systems 2025 — "T cell receptor cross-reactivity prediction improved by a comprehensive mutational scan database"  
PMID: 38370810  
Entries: 19,145 (MHC I) + 5,730 (MHC II) = ~24,875 TCR–peptide activity measurements

## What makes BATCAVE unique

Unlike IEDB/VDJdb/McPAS which record only native epitope–TCR pairs, BATCAVE contains **complete single amino acid mutational scans**: each of 151 TCRs is tested against every single-residue variant of 25 immunogenic epitopes. This makes it the gold standard for studying TCR cross-reactivity and binding tolerance.

## Files

| File | Entries | MHC class | Description |
|---|---|---|---|
| `TCR_pMHCI_mutational_scan.csv` | 19,145 | I | Human and mouse TCRs vs. MHC class I epitope variants |
| `TCR_pMHCII_mutational_scan.csv` | 5,730 | II | Human TCRs vs. MHC class II epitope variants |

## Columns (both files share the same schema)

| Column | Description |
|---|---|
| `tcr` | TCR name/identifier |
| `va` | Alpha V gene (short form) |
| `vb` | Beta V gene (short form) |
| `cdr3a` | CDR3 alpha amino acid sequence |
| `cdr3b` | CDR3 beta amino acid sequence |
| `trav` | TRAV gene segment |
| `traj` | TRAJ gene segment |
| `trbv` | TRBV gene segment |
| `trbd` | TRBD gene segment |
| `trbj` | TRBJ gene segment |
| `assay` | Experimental readout (e.g. `TNF Secretion`, `TScan II`) |
| `tcr_source_organism` | `Human` or `Mouse` |
| `index_peptide` | The native/reference epitope for this TCR |
| `mhc` | HLA/MHC restriction |
| `pmid` | PubMed ID |
| `peptide_type` | `autoimmune`, `viral`, `cancer` |
| `peptide` | Actual peptide tested (may be single-mutant of `index_peptide`) |
| `peptide_activity` | TCR activation readout (0–100 for class I; 0–1 for class II) |

## Unified schema mapping

```
cdr3_alpha   ← cdr3a
cdr3_beta    ← cdr3b
epitope      ← index_peptide   (native epitope; use peptide for mutant-aware matching)
mhc          ← mhc
v_alpha      ← trav
j_alpha      ← traj
v_beta       ← trbv
j_beta       ← trbj
```

## Refresh

```bash
BASE="https://github.com/meyer-lab-cshl/BATMAN-paper/raw/main/results_batman/tcr_epitope_datasets/mutational_scan_datasets/database"
curl -L "${BASE}/TCR_pMHCI_mutational_scan_database.xlsx" -o pMHCI.xlsx
curl -L "${BASE}/TCR_pMHCII_mutational_scan_database.xlsx" -o pMHCII.xlsx
python3 -c "import openpyxl, csv, sys; [csv.writer(open(f'{n}.csv','w')).writerows(openpyxl.load_workbook(f'{n}.xlsx',read_only=True).active.iter_rows(values_only=True)) for n in ['pMHCI','pMHCII']]"
```

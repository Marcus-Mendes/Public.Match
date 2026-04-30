# 10x Genomics — Chromium Single Cell Immune Profiling (pMHC Dextramer)

Source: 10x Genomics published dataset — 4-donor human PBMC  
Reference: 10x Genomics Chromium Single Cell V(D)J + Feature Barcode  
Entries: ~189,515 cells across 4 donors

## Files

| File | Donor | Cells |
|---|---|---|
| `vdj_v1_hs_aggregated_donor1_binarized_matrix.csv` | donor1 | ~46,526 |
| `vdj_v1_hs_aggregated_donor2_binarized_matrix.csv` | donor2 | ~77,854 |
| `vdj_v1_hs_aggregated_donor3_binarized_matrix.csv` | donor3 | ~37,824 |
| `vdj_v1_hs_aggregated_donor4_binarized_matrix.csv` | donor4 | ~27,308 |

## Key columns

| Column | Description |
|---|---|
| `barcode` | Cell barcode |
| `donor` | Donor ID |
| `cell_clono_cdr3_aa` | CDR3 amino acid sequences — format: `TRA:seq;TRB:seq` (may have multiple TRA) |
| `cell_clono_cdr3_nt` | CDR3 nucleotide sequences (same format) |
| `CD3`, `CD8a`, `CD4`, ... | Cell surface marker counts (CITE-seq) |
| `A0201_GILGFVFTL_Flu-MP_Influenza` | pMHC tetramer UMI count — format: `HLA_epitope_antigen_pathology` |
| `A0201_GILGFVFTL_Flu-MP_Influenza_binder` | Binary binder call (1.0 = binder) |

## Epitopes covered (examples)

Covers 55 pMHC specificities including:
- CMV pp65, IE-1 (HLA-A*01:01, A*02:01, A*24:02, B*35:01, B*07:02)
- EBV BMLF1, EBNA-3A/3B, LMP1/2 (multiple HLAs)
- Influenza MP (HLA-A*02:01)
- HIV Gag (HLA-A*02:01)
- Cancer antigens: MART-1, gp100, NY-ESO-1, MAGE-A1/3, Tyrosinase (HLA-A*02:01)
- Negative controls (NC)

## Parsing CDR3

```python
# cell_clono_cdr3_aa = "TRA:CAASXXX;TRA:CAAYYYY;TRB:CASSZZZ"
chains = cell_clono_cdr3_aa.split(';')
alpha = [c[4:] for c in chains if c.startswith('TRA:')]
beta  = [c[4:] for c in chains if c.startswith('TRB:')]
```

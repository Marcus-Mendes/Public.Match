# IEDB — Immune Epitope Database

Source: https://www.iedb.org  
Downloaded: 2026-04-30  
TCR entries: 226,280 rows (223,186 with at least one CDR3)

## Files

| File | Description | Rows |
|---|---|---|
| `tcr_full_v3_slim.csv` | TCR receptor export — CDR3 + epitope + MHC columns only | 226,280 |
| `epitope_table_export_*.csv` | Full epitope metadata table (no CDR3; use for epitope lookups) | ~75,000 |

## tcr_full_v3_slim.csv columns

| Column | Description |
|---|---|
| `receptor_id` | IEDB receptor ID |
| `receptor_type` | e.g. `alphabeta` |
| `epitope` | Epitope peptide name |
| `source_molecule` | Antigen source molecule |
| `source_organism` | Source organism (e.g. `Homo sapiens`) |
| `mhc_allele` | MHC restriction (e.g. `HLA-A*02:01`) |
| `chain1_type` | `alpha` or `beta` |
| `chain1_v` | Chain 1 V gene (curated) |
| `chain1_j` | Chain 1 J gene (curated) |
| `chain1_cdr3` | Chain 1 CDR3 amino acid sequence |
| `chain2_type` | `alpha` or `beta` |
| `chain2_v` | Chain 2 V gene (curated) |
| `chain2_j` | Chain 2 J gene (curated) |
| `chain2_cdr3` | Chain 2 CDR3 amino acid sequence |

Note: chain1 = alpha for `alphabeta` receptors; chain2 = beta.

## Refresh

```bash
curl -L "https://www.iedb.org/downloader.php?file_name=doc/receptor_full_v3.zip" \
  -o receptor_full_v3.zip
unzip -p receptor_full_v3.zip tcr_full_v3.csv | python3 iedb_slim.py > tcr_full_v3_slim.csv
```

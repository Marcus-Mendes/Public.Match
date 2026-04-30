# McPAS-TCR — Pathology-Associated TCR Database

Source: https://friedmanlab.weizmann.ac.il/McPAS-TCR  
Downloaded: 2026-04-30  
Entries: ~40,779 rows (manually curated)

## File

`McPAS-TCR.csv` — comma-separated, downloaded directly from the Weizmann Institute portal.

## Key columns

| Column | Description |
|---|---|
| `CDR3.alpha.aa` | Alpha chain CDR3 amino acid sequence |
| `CDR3.beta.aa` | Beta chain CDR3 amino acid sequence |
| `Epitope.peptide` | Epitope peptide sequence |
| `MHC` | MHC restriction |
| `Category` | Disease category (Pathogens, Autoimmune, Cancer, Other) |
| `Pathology` | Specific disease/condition |
| `TRAV` | Alpha V gene |
| `TRAJ` | Alpha J gene |
| `TRBV` | Beta V gene |
| `TRBD` | Beta D gene |
| `TRBJ` | Beta J gene |
| `Species` | Human or Mouse |
| `Antigen.protein` | Full antigen protein name |
| `T.Cell.Type` | CD4/CD8 |

## Unified schema mapping

```
cdr3_alpha  ← CDR3.alpha.aa
cdr3_beta   ← CDR3.beta.aa
epitope     ← Epitope.peptide
mhc         ← MHC
v_alpha     ← TRAV
v_beta      ← TRBV
```

## Refresh

Visit https://friedmanlab.weizmann.ac.il/McPAS-TCR and use "Download CSV" on the full database.

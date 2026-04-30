# MixTCRpred Training Dataset

Source: https://github.com/GfellerLab/MixTCRpred  
Publication: Nature Communications 2024 — "Deep learning predictions of TCR-epitope interactions reveal epitope-specific chains in dual alpha T cells"  
Zenodo: https://zenodo.org/records/10806391  
Entries: 17,715 αβ TCR–epitope pairs across 146 pMHC specificities

## File

`full_training_set_146pmhc.csv` — the complete training dataset used to build MixTCRpred.

## Columns

| Column | Description |
|---|---|
| `epitope` | Epitope peptide sequence |
| `cdr3_TRA` | CDR3 alpha amino acid sequence |
| `cdr3_TRB` | CDR3 beta amino acid sequence |
| `TRAV` | Alpha V gene |
| `TRAJ` | Alpha J gene |
| `TRBV` | Beta V gene |
| `TRBJ` | Beta J gene |
| `MHC` | MHC restriction (e.g. `HLA-A*02:01`, `H2-IAb`) |
| `MHC_class` | `MHCI` or `MHCII` |
| `species` | `HomoSapiens` or `MusMusculus` |

## Unified schema mapping

```
cdr3_alpha  ← cdr3_TRA
cdr3_beta   ← cdr3_TRB
epitope     ← epitope
mhc         ← MHC
v_alpha     ← TRAV
j_alpha     ← TRAJ
v_beta      ← TRBV
j_beta      ← TRBJ
```

## Coverage

- 146 distinct pMHC specificities
- Both MHC class I and class II
- Human and mouse TCRs
- Includes viral (flu, EBV, CMV), cancer (NY-ESO-1, MART-1), and autoimmune antigens

## Refresh

```bash
curl -L "https://raw.githubusercontent.com/GfellerLab/MixTCRpred/main/full_training_set_146pmhc.csv" \
  -o full_training_set_146pmhc.csv
```

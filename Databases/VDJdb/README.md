# VDJdb — 2025-12-29 release

Source: https://vdjdb.cdr3.net  
Release: 2025-12-29 (`vdjdb-2025-12-29.zip`)  
Entries: 145,408 rows (α and β chains listed separately; linked by `complex.id`)

## File

`vdjdb.slim.txt` — tab-separated, extracted from the official release zip.

## Key columns

| Column | Description |
|---|---|
| `gene` | `TRA` (alpha) or `TRB` (beta) |
| `cdr3` | CDR3 amino acid sequence |
| `v.segm` | V gene segment |
| `j.segm` | J gene segment |
| `species` | e.g. `HomoSapiens` |
| `antigen.epitope` | Epitope peptide sequence |
| `antigen.gene` | Antigen gene name |
| `antigen.species` | Pathogen/source species |
| `mhc.a` | MHC alpha chain (e.g. `HLA-A*02:01`) |
| `mhc.b` | MHC beta chain |
| `mhc.class` | `MHCI` or `MHCII` |
| `complex.id` | Links TRA/TRB rows from the same TCR clone |
| `vdjdb.score` | Confidence score (0–3) |

## Unified schema mapping

```
cdr3_alpha  ← rows where gene == "TRA", column cdr3
cdr3_beta   ← rows where gene == "TRB", column cdr3
epitope     ← antigen.epitope
mhc         ← mhc.a
v_alpha     ← v.segm (TRA rows)
v_beta      ← v.segm (TRB rows)
```

## Refresh

```bash
VER=$(date +%Y-%m-%d)
curl -L "https://github.com/antigenomics/vdjdb-db/releases/latest/download/vdjdb-${VER}.zip" \
  -o vdjdb-latest.zip
unzip -p vdjdb-latest.zip vdjdb.slim.txt > vdjdb.slim.txt
```

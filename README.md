# Public.Match

**Hackathon Project** | Spatial Immunology & TCR Repertoire Analysis

---

## The Problem

T-cell receptor (TCR) repertoires from cancer patients contain valuable signals — but identifying which TCRs are "public" (shared across individuals and described in curated databases) requires tools that are tightly coupled to a single database. **TCRMatch**, for example, only queries IEDB. Researchers working with VDJdb, McPAS-TCR, or custom epitope databases have no unified solution.

## Our Idea

**Public.Match** is a generalized TCR public-sequence matching tool that:

- Accepts patient TCR repertoires (CDR3α/β sequences) as input
- Searches across multiple curated databases — **IEDB**, **VDJdb**, **McPAS-TCR**, **10x Genomics pMHC**, and others
- Returns matched public TCRs with epitope annotations, HLA restrictions, and match scores
- Is built with **Claude Code**, using AI-assisted development to rapidly prototype and extend the tool during the hackathon

Think of it as TCRMatch — but database-agnostic.

## Why It Matters

| Today | With Public.Match |
|---|---|
| Run TCRMatch → IEDB only | Single query → IEDB + VDJdb + McPAS + 10x |
| Manual format conversion per DB | Unified input/output schema |
| No cross-database deduplication | Merged, ranked hits across all sources |

Identifying public TCRs that recognize known epitopes helps distinguish **antigen-specific** from **bystander** T cells — a key step in spatial immunology pipelines like our own [Soleil](https://github.com/Marcus-Mendes) engagement scoring framework.

## Approach

1. **Unify database schemas** — normalize all sources into a common `cdr3_alpha / cdr3_beta / epitope / mhc / v_gene / j_gene / source` format
2. **Extend or wrap TCRMatch** — reuse its edit-distance / GLIPH2-inspired scoring logic against any database
3. **Build a CLI** — `public-match --input repertoire.tsv --db iedb vdjdb mcpas 10x --score 0.97`
4. **Claude Code as co-developer** — use Claude Code to accelerate implementation, handle format parsing edge cases, and generate test cases

## Databases

| Database | Entries | CDR3α | CDR3β | Epitope | HLA | Folder |
|---|---|---|---|---|---|---|
| [IEDB](https://www.iedb.org/) | 226,280 TCR records | ✓ | ✓ | ✓ | ✓ | `Databases/IEDB/` |
| [VDJdb](https://vdjdb.cdr3.net/) | 145,408 chain records | ✓ | ✓ | ✓ | ✓ | `Databases/VDJdb/` |
| [McPAS-TCR](https://friedmanlab.weizmann.ac.il/McPAS-TCR/) | 40,779 | ✓ | ✓ | ✓ | ✓ | `Databases/McPAS/` |
| [10x Genomics pMHC](https://www.10xgenomics.com/) | 189,515 cells / 4 donors | ✓ | ✓ | ✓ (55 pMHC) | ✓ | `Databases/10xDcode/` |
| [MixTCRpred](https://github.com/GfellerLab/MixTCRpred) | 17,715 αβ pairs | ✓ | ✓ | ✓ (146 pMHC) | ✓ | `Databases/MixTCRpred/` |
| [BATCAVE](https://github.com/meyer-lab-cshl/BATMAN-paper) | 24,875 TCR–peptide measurements | ✓ | ✓ | ✓ (mutational scan) | ✓ | `Databases/BATCAVE/` |
| [OTS](https://opig.stats.ox.ac.uk/webapps/ots/) | 1.63M non-redundant paired αβ | ✓ | ✓ | — (publicness only) | — | `Databases/OTS/` (manual download) |

### Unified schema

All databases are mapped to a common record format:

```
cdr3_alpha    CDR3 amino acid sequence of the alpha chain
cdr3_beta     CDR3 amino acid sequence of the beta chain
epitope       Epitope peptide sequence
mhc           MHC/HLA restriction (e.g. HLA-A*02:01)
v_alpha       TRAV gene
j_alpha       TRAJ gene
v_beta        TRBV gene
j_beta        TRBJ gene
source        Database of origin (IEDB / VDJdb / McPAS / 10x)
```

### Databases on the roadmap

| Database | Entries | Notes |
|---|---|---|
| [TCRdb 2.0](https://guolab.wchscu.cn/TCRdb2/) | ~700M sequences | Broad clinical coverage; no epitope labels |
| [STCRDab](http://opig.stats.ox.ac.uk/webapps/stcrdab/) | ~1,000 | 3D structural data from PDB |
| [PIRD](https://db.cngb.org/pird/) | large | Pan Immune Repertoire Database; China National GeneBank |
| [ePytope-TCR datasets](https://www.cell.com/cell-genomics/fulltext/S2666-979X(25)00202-2) | 21 datasets / 762 epitopes | 2025 benchmarking collection on Zenodo |

## Hackathon Deliverable

A working CLI prototype that takes a CDR3 repertoire file and returns matched public TCRs from all four databases, with a unified output format and match score.

---

*Built at Hackathon · 2026 · with Claude Code*

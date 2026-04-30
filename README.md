# Public.Match

**Hackathon Project** | TCR Repertoire Analysis

---

## The Problem

T-cell receptor (TCR) repertoires from cancer patients contain valuable signals — but identifying which TCRs are "public" (shared across individuals and described in curated databases) requires tools that are tightly coupled to a single database. **TCRMatch**, for example, only queries IEDB. Researchers working with VDJdb, McPAS-TCR, or custom epitope databases have no unified solution.

## Our Idea

**Public.Match** is a generalized TCR public-sequence matching tool that:

- Accepts patient TCR repertoires (CDR3α/β sequences) as input
- Searches across multiple curated databases — **IEDB**, **VDJdb**, and others
- Returns matched public TCRs with epitope annotations, HLA restrictions, and match scores
- Is built with **Claude Code**, using AI-assisted development to rapidly prototype and extend the tool during the hackathon

Think of it as TCRMatch — but database-agnostic.

## Why It Matters

| Today | With Public.Match |
|---|---|
| Run TCRMatch → IEDB only | Single query → IEDB + VDJdb + custom DBs |
| Manual format conversion per DB | Unified input/output schema |
| No cross-database deduplication | Merged, ranked hits across all sources |

Identifying public TCRs that recognize known epitopes helps distinguish **antigen-specific** from **bystander** T cells — a key step in spatial immunology pipelines like our own [Soleil](https://github.com/Marcus-Mendes) engagement scoring framework.

## Approach

1. **Unify database schemas** — normalize IEDB, VDJdb (and optionally McPAS-TCR) into a common CDR3 + epitope + HLA format
2. **Extend or wrap TCRMatch** — reuse its edit-distance / GLIPH2-inspired scoring logic against any database
3. **Build a CLI** — `public-match --input repertoire.tsv --db iedb vdjdb --score 0.97`
4. **Claude Code as co-developer** — use Claude Code to accelerate implementation, handle format parsing edge cases, and generate test cases

## Databases Targeted

| Database | Antigen Coverage | CDR3β entries |
|---|---|---|
| [IEDB](https://www.iedb.org/) | Broad (viral, tumor, autoimmune) | ~500k |
| [VDJdb](https://vdjdb.cdr3.net/) | Viral-focused (CMV, EBV, flu) | ~90k |
| McPAS-TCR | Pathology-associated | ~5k |

## Hackathon Deliverable

A working CLI prototype that takes a CDR3 repertoire file and returns matched public TCRs from at least two databases, with a unified output format and match score.

---

*Built at [Hackathon Name] · 2026*

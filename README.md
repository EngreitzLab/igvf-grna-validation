# igvf-grna-validation

Tools for validating and correcting IGVF Per-Guide Metadata files before submission to the IGVF Data Portal.

## Background

The [IGVF Consortium](https://igvf.org/) requires gRNA library metadata to be submitted in a standardized TSV format defined in:

> `input/FINAL_FINAL_Format_for_Per_Guide_Metadata_Finalized_11_26_25.pdf`

This repository contains a validator that checks files against that spec and per-file correction scripts that fix common issues.

---

## Scripts

### `validate_grna_file.py` — Validator

Reads a Per-Guide Metadata TSV (or TSV.gz) and prints a structured diagnostic report. **Does not modify any files.**

**Usage:**
```bash
python3 validate_grna_file.py <input_file.tsv.gz>
python3 validate_grna_file.py <input_file.tsv.gz> --json-out problems.json
```

**Requirements:**
- Python 3, pandas, requests
- GENCODE v43 GTF (`input/IGVFFI9573KOZR.gtf.gz`) — downloaded automatically on first run if absent

**Checks performed:**

| Check | Field(s) | Description |
|---|---|---|
| S1 | structure | Required columns present; optional columns noted; unknown columns flagged |
| S2 | guide_id, spacer | guide_id values are unique; each guide_id maps to exactly one spacer |
| T1 | type | Values are from the spec enum: `targeting`, `safe-targeting`, `non-targeting`, `positive control`, `negative control`, `variant` |
| T2 | targeting | Values are `True` or `False` (title case) |
| T3 | targeting | `non-targeting` and `safe-targeting` rows must have `targeting=False` |
| T4 | strand | Non-empty and `+`/`-` for targeting rows |
| T5 | guide_chr/start/end | Non-empty for targeting rows |
| E1 | genomic_element | Non-empty for targeting rows |
| E2 | genomic_element | Values are from the spec enum: `promoter`, `enhancer`, `insulator`, `silencer`, `distal element`, `splice site`, `gene`, `variant` |
| E2b | genomic_element | Warns when `enhancer` is used; recommends `distal element` |
| E3 | genomic_element | Must be empty for `non-targeting` rows |
| N1 | intended_target_name | Non-empty for targeting rows |
| N2 | intended_target_name | Must be an Ensembl ID (`ENSG…`) for `promoter`/`gene`/`splice site` elements; must be a coord string (`chr:start-end`) for `distal element`/`enhancer`/`insulator`/`silencer` |
| N3 | intended_target_name | Must be empty for `non-targeting` rows |
| C1 | intended_target_chr/start/end | Non-empty for targeting rows |
| C2 | intended_target_chr/start/end | Target window must contain the guide position (per-row check) |
| C3 | intended_target_start/end | Target window must span all guides in a guide-prefix group (cross-row check) |
| P1 | putative_target_genes | Required for `positive control` rows with a distal-type element |
| D1 | description | Warns if column is absent or all-empty |

**Output format:**

Issues are grouped by field and labeled `✗` (error) or `!` (warning). A summary line reports total error and warning category counts. Example:

```
=== IGVF Per-Guide Metadata Validation Report ===
File: output/IGVFFI1207NRVS.tsv.gz
Rows: 220  |  targeting: 130  |  positive control: 60  |  non-targeting: 20  |  safe-targeting: 10

── targeting ─────────────────────────────────────────────────────────────────
  ✗  220 rows have invalid targeting value (must be 'True' or 'False'): {'TRUE': 130, ...}

── Fields with no issues ─────────────────────────────────────────────────────
  ✓  guide_id
  ✓  spacer
  ...

══ SUMMARY ═══════════════════════════════════════════════════════════════════
  1 error category, 0 warning categories found
```

The optional `--json-out` flag writes a machine-readable report for use with `fix_interactive.py`.

---

### `fix_interactive.py` — Interactive fixer

Reads a JSON problem report produced by `validate_grna_file.py --json-out` and offers each fixable issue interactively. Shared fix functions are also imported by the per-file scripts below.

```bash
python3 fix_interactive.py <input_file> --problems <report.json> [--output <output.tsv.gz>]
```

### Per-file fix scripts

Each file with non-trivial issues has a dedicated script that applies targeted corrections and writes a corrected TSV.gz to `output/`.

| Script | File | Notes |
|---|---|---|
| `fix_IGVFFI9754AGFB.py` | 37,637-row screen | Patch script applied on top of `fix_grna_file.py` output |
| `fix_IGVFFI1207NRVS.py` | 220-row GATA1 enhancer screen | ENSG lookup for 5 TSS control genes |
| `fix_IGVFFI4290JQVQ.py` | 296-row EC enhancer screen | ENSG lookup for 5 gene-element targets |
| `fix_IGVFFI0580WJFK.py` | 16,201-row WTC11 random screen | TSS intersection against GENCODE v43 for 83 promoter windows |

---

## Column spec (summary)

| Column | Required | Notes |
|---|---|---|
| guide_id | yes | Unique per row |
| spacer | yes | Guide spacer sequence |
| targeting | yes | `True` or `False` (title case) |
| type | yes | `targeting`, `safe-targeting`, `non-targeting`, `positive control`, `negative control`, `variant` |
| guide_chr | yes | |
| guide_start | yes | |
| guide_end | yes | |
| strand | yes | `+` or `-` |
| pam | yes | |
| genomic_element | yes (if targeting=True) | `promoter`, `distal element`, `enhancer`, `insulator`, `silencer`, `splice site`, `gene`, `variant` |
| intended_target_name | yes (if targeting=True) | ENSG ID for promoter/gene/splice site; `chr:start-end` for distal/enhancer/insulator/silencer |
| intended_target_chr | yes (if targeting=True) | Must match guide_chr |
| intended_target_start | yes (if targeting=True) | Must be ≤ guide_start |
| intended_target_end | yes (if targeting=True) | Must be ≥ guide_end |
| putative_target_genes | optional | Required for positive control rows with distal-type element |
| reporter | optional | |
| imperfect | optional | |
| description | optional | Recommended: human-readable label for the guide group |

---

## Reference data

- **GENCODE v43 GTF** (`input/IGVFFI9573KOZR.gtf.gz`): used for gene name → Ensembl ID lookup and TSS position intersection. Downloaded automatically from the IGVF Data Portal on first run.
- **Spec PDF** (`input/FINAL_FINAL_Format_for_Per_Guide_Metadata_Finalized_11_26_25.pdf`): authoritative format specification.

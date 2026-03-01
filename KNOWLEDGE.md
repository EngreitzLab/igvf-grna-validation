# IGVF Per-Guide Metadata: Knowledge Reference

This file captures all domain knowledge needed to fix future IGVF Per-Guide Metadata
(gRNA) files. It is a living reference document — update it as new patterns are found.

---

## 1. Spec Overview

**File format:** Tab-separated values (TSV), gzip-compressed, with a header row.
**14 required columns + 4 optional columns** (order unspecified):

Required:
```
guide_id, spacer, targeting, type, guide_chr, guide_start, guide_end, strand, pam,
genomic_element, intended_target_name, intended_target_chr, intended_target_start,
intended_target_end
```

Optional (marked `(-)` in spec — may be absent or empty):
```
putative_target_genes, reporter, imperfect, description
```

**Source spec PDF:** `input/FINAL_FINAL_Format_for_Per_Guide_Metadata_Finalized_11_26_25.pdf`

---

## 2. Column-by-Column Rules

### `type`
- Required for all rows.
- Must be one of (full spec enum): `"targeting"`, `"safe-targeting"`, `"non-targeting"`,
  `"positive control"`, `"negative control"`, `"variant"`.
- Determines how all other columns are validated.

### `targeting`
- Required for all rows.
- Must be `"True"` or `"False"` (title case — Python boolean string representation).
- If `type ∈ {"non-targeting", "safe-targeting"}`, `targeting` **must** be `"False"`.
- If `type == "targeting"`, `targeting` should be `"True"` (warn if not).
- Common failure: input files use `"TRUE"` / `"FALSE"` (uppercase) or `"true"` / `"false"` (lowercase) — must be title case.
- Element fields (`guide_chr`, `guide_start`, `guide_end`, `strand`, `pam`, `genomic_element`,
  `intended_target_name`, `intended_target_chr`, `intended_target_start`, `intended_target_end`)
  are **required** when `targeting == "True"`, and **not required** (can be empty) when `targeting == "False"`.

### `genomic_element`
- Required when `targeting == "True"`; should be empty when `targeting == "False"`.
- Full valid enum (per spec): `"variant"`, `"promoter"`, `"enhancer"`, `"insulator"`,
  `"silencer"`, `"distal element"`, `"splice site"`, `"gene"`.
- **Preferred value for non-promoter regulatory elements**: `"distal element"` (use this
  unless a more specific classification like `"enhancer"` is firmly established).

#### Classification rules (applied in order):
1. `"TSS"` in `description` (case-sensitive substring match) → `"promoter"`
2. `guide_id` starts with `"9p21_add_"` AND no TSS in description → `"distal element"`
3. All other targeting rows → `"promoter"`

### `intended_target_name`
- Required when `targeting == "True"`; should be empty when `targeting == "False"`.
- Value depends on `genomic_element`:
  - `"promoter"`, `"splice site"`, `"gene"`: Ensembl gene ID matching `^ENSG\d+$`
  - `"enhancer"`, `"insulator"`, `"silencer"`, `"distal element"`: coordinate string `^chr\S+:\d+-\d+$`
  - `"variant"`: normalized SPDI id

### `intended_target_chr`
- Required when `targeting == "True"`; should be empty when `targeting == "False"`.
- Value: `guide_chr` of the target group (first value in group, should be uniform).

### `intended_target_start`
- Required **only** for `type == "targeting"` rows.
- Value: `min(guide_start)` across all targeting guides with the same `description`.

### `intended_target_end`
- Required **only** for `type == "targeting"` rows.
- Value: `max(guide_end)` across all targeting guides with the same `description`.

### `putative_target_genes`, `reporter`, `imperfect`
- No automated transformations applied in first file — leave as-is unless spec changes.

### `description`
- Free-text field; holds target name / locus information.
- Used as grouping key for coordinate bounds computation.
- Used as source for gene symbol extraction.
- **Do not modify** this column.

---

## 3. Gene Symbol Extraction Algorithm

Used to derive `intended_target_name` for promoter guides.

```python
def extract_gene_candidates(desc: str) -> list:
    """Return ordered list of candidate gene symbols from a description string."""
    cleaned = re.sub(r"_TSS\d*", "", desc)    # strip _TSS, _TSS2, _TSS123 …
    if "_and_" in cleaned:
        parts = cleaned.split("_and_")
        # priority: first token, then last, then any middle tokens
        candidates = [parts[0], parts[-1]] + parts[1:-1]
    else:
        candidates = [cleaned]
    # Also try stripping _alt suffix from each candidate
    extras = [c[:-4] for c in candidates if c.endswith("_alt")]
    return candidates + extras
```

Then look up each candidate in the GENCODE v43 gene map:
```python
for candidate in candidates:
    if candidate in gene_map:
        return gene_map[candidate]
```

---

## 4. Resolution Strategy Hierarchy for `intended_target_name` Failures

When the direct gene symbol lookup fails, apply these passes in order:

### Pass 1: Mouse gene capitalization
- Condition: `fallback[0].isupper() and fallback[1].islower()` (e.g., `"Cmip"`, `"Crip2"`)
- Fix: `fallback.upper()` → look up in gene_map (e.g., `"CMIP"` → ENSG ID)
- In first file: resolved 30 descriptions

### Pass 2: Enhancer descriptions → coordinate string
- Condition: `description.startswith("Enhancer")`
- Fix: use `"chr:start-end"` coordinate string derived from guide coordinate bounds
- This is valid because enhancer loci are "distal element"-like even if genomic_element was set to "promoter" — but check carefully; in first file these were classified correctly as `"distal element"` and this pass fixed the name accordingly
- In first file: resolved 60 descriptions

### Pass 3: Ensembl REST API (renamed/obsolete human genes)
- Uses two endpoints:
  - Batch POST: `https://rest.ensembl.org/lookup/symbol/homo_sapiens` with `{"symbols": [...]}`
  - Xref GET: `https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{symbol}?object_type=gene`
- Rate limit: sleep 70ms between xref GET requests (≤15/sec)
- In first file: resolved 75 descriptions

### Pass 4 (fallback): Gene symbol placeholder
- Leave the gene symbol as the value for `intended_target_name`
- Flag for manual review
- In first file: 4 remaining (see §6 below)

---

## 5. Reference Files

### GENCODE v43 GTF
- IGVF reference file ID: `IGVFFI9573KOZR`
- Download URL: `https://api.data.igvf.org/reference-files/IGVFFI9573KOZR/@@download/IGVFFI9573KOZR.gtf.gz`
- Cached locally at: `input/IGVFFI9573KOZR.gtf.gz`
- Parser: extract `gene` feature rows; parse `gene_id` (strip version suffix `.N`) and `gene_name` attributes

```python
def load_gene_map(gtf_path: str) -> dict:
    """Return {gene_name: ensembl_id_no_version} from GTF gene features."""
    gene_map = {}
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = fields[8]
            gid = re.search(r'gene_id "([^"]+)"', attrs)
            gnm = re.search(r'gene_name "([^"]+)"', attrs)
            if gid and gnm:
                ensembl_id = gid.group(1).split(".")[0]
                gene_map[gnm.group(1)] = ensembl_id
    return gene_map
```

### Ensembl REST API endpoints
- Batch symbol lookup: `POST https://rest.ensembl.org/lookup/symbol/homo_sapiens`
  - Body: `{"symbols": ["GENE1", "GENE2", ...]}`
  - Headers: `Content-Type: application/json`, `Accept: application/json`
- Xref/symbol lookup: `GET https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{symbol}?object_type=gene`
  - Filter response: `entry["type"] == "gene"`, return `entry["id"]`

---

## 6. Known Failure Patterns from First File (IGVFFI9754AGFB)

### Unresolvable gene symbols (manual review required)
| Gene Symbol | Reason |
|---|---|
| `AKAP2` | Renamed in Ensembl; no current xref match |
| `LOC101927497` | NCBI LOC entry; no standard Ensembl ID |
| `LOC391322` | NCBI LOC entry; no standard Ensembl ID |
| `SMIM11B` | Possibly obsolete or merged gene |

### Pattern: `LOC\d+` genes
- These are NCBI LocusLink provisional gene names.
- They typically have no Ensembl ID and cannot be resolved automatically.
- Flag these in the validation report; they require manual cross-reference via NCBI Gene.

### Pattern: Mouse gene names
- Identified by: first char uppercase, second char lowercase (e.g., `Cmip`, `Crip2`, `Slc7a8`)
- Fix: uppercase the entire symbol, then look up in GENCODE v43 (human gene name)

### Pattern: Enhancer descriptions
- Identified by: description starts with `"Enhancer"`
- Fix: use coordinate string `"chr:start-end"` from guide bounds

---

## 7. Template Scripts

The following scripts in `260301-gRNAFileCorrection/` serve as templates for future files:

- **`fix_grna_file.py`** — applies all automated fixes (hardcoded to first file's paths)
- **`validate_grna_file.py`** — CLI validator, reads any `.tsv.gz` and reports issues

To adapt for a new file:
1. Copy `fix_grna_file.py`, update `INPUT_FILE` and `OUTPUT_FILE` constants.
2. Check if new file has the same `guide_id` prefix patterns for `genomic_element` classification (the `9p21_add_` prefix is specific to the first file's distal element guides).
3. Run `validate_grna_file.py <new_file>` first to triage which issues need fixing.

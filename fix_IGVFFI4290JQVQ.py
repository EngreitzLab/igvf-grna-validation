#!/usr/bin/env python3
"""
Fix script for IGVFFI4290JQVQ (296-row CRISPRi EC enhancer screen).

Data breakdown:
  64  targeting gene-element rows  (SOX2_TSS x14, KDR_TSS x10, HAND1_TSS x10,
                                    HEY2_TSS x15, DLL4_TSS x15)
  192 targeting enhancer rows      (per-window coords already correct)
  20  non-targeting controls
  20  safe-targeting controls

Issues fixed:
  1. targeting case: TRUE/FALSE → True/False
  2. description: absent column → fill from intended_target_name (per-row label)
  3. intended_target_chr/start/end for gene rows: empty → span of guide coords
     (grouped by intended_target_name, e.g. SOX2_TSS)
  4. genomic_element='gene' → 'promoter' for gene rows
  5. intended_target_name for gene rows: 'GENE_TSS' → ENSG ID
  6. genomic_element='enhancer' → 'distal element'

NOTE: All 296 rows share one guide_id prefix (221204_EC_Enhancer_Screen), so
fix_coords_recompute (which groups by prefix) cannot be used — it would collapse
all rows into a single giant span.  Gene-row coords are instead grouped by
intended_target_name.

Usage:
    python3 fix_IGVFFI4290JQVQ.py
"""

import os
import re
import sys

import pandas as pd

# ── Imports from shared modules ────────────────────────────────────────────────
from fix_interactive import fix_targeting_case, COLUMN_ORDER, _is_nan_like
from validate_grna_file import load_gene_map, GTF_PATH, GTF_URL, download_gtf

INPUT_FILE  = "input/IGVFFI4290JQVQ.tsv"
OUTPUT_FILE = "output/IGVFFI4290JQVQ.tsv.gz"


# ── Fix helpers ────────────────────────────────────────────────────────────────

def fill_description_from_intended_target_name(df: pd.DataFrame) -> pd.DataFrame:
    """Fill description column from intended_target_name where description is empty.

    Must run BEFORE fix_gene_rows(), which overwrites intended_target_name with
    ENSG IDs.  This preserves human-readable labels (e.g. 'SOX2_TSS',
    'chr3:181756680-181756981') as the description.

    Does NOT use fix_description_from_guide_id — all 296 rows share the same
    guide_id prefix and it carries no group-level information.
    """
    if "description" not in df.columns:
        df["description"] = ""

    mask = _is_nan_like(df["description"]) & ~_is_nan_like(df["intended_target_name"])
    n = mask.sum()
    df.loc[mask, "description"] = df.loc[mask, "intended_target_name"]
    print(f"  description filled for {n} rows from intended_target_name")
    return df


def fix_coords_gene_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Compute intended_target_chr/start/end for gene rows.

    Groups by intended_target_name (e.g. 'SOX2_TSS') and assigns the min
    guide_start and max guide_end as the spanning target coordinates.
    Enhancer rows are untouched — their per-window coords are already correct.
    """
    is_true = df["targeting"].str.upper() == "TRUE"
    mask = is_true & (df["genomic_element"] == "gene")

    tgt = df.loc[mask].copy()
    tgt["_gs"] = pd.to_numeric(tgt["guide_start"], errors="coerce")
    tgt["_ge"] = pd.to_numeric(tgt["guide_end"],   errors="coerce")

    grp = tgt.groupby("intended_target_name").agg(
        _chr=("guide_chr", "first"),
        _start=("_gs", "min"),
        _end=("_ge", "max"),
    )

    for itn, row in grp.iterrows():
        pmask = mask & (df["intended_target_name"] == itn)
        df.loc[pmask, "intended_target_chr"]   = row["_chr"]
        df.loc[pmask, "intended_target_start"] = int(row["_start"])
        df.loc[pmask, "intended_target_end"]   = int(row["_end"])

    print(f"  Recomputed coords for {len(grp)} gene-row groups:")
    for itn, row in grp.iterrows():
        print(f"    {itn}: {row['_chr']}:{int(row['_start'])}-{int(row['_end'])}")

    return df


def fix_gene_rows(df: pd.DataFrame, gene_map: dict) -> pd.DataFrame:
    """Fix genomic_element and intended_target_name for gene-element rows.

    - genomic_element: 'gene' → 'promoter'
    - intended_target_name: 'GENE_TSS' → ENSG ID  (strip '_TSS' suffix for lookup)
    """
    is_true = df["targeting"].str.upper() == "TRUE"
    mask = is_true & (df["genomic_element"] == "gene")
    n_rows = mask.sum()

    # Fix genomic_element
    df.loc[mask, "genomic_element"] = "promoter"
    print(f"  genomic_element: 'gene' → 'promoter' for {n_rows} rows")

    # Fix intended_target_name
    failures = []
    for idx in df.index[mask]:
        itn = df.at[idx, "intended_target_name"]
        # e.g. 'SOX2_TSS' → 'SOX2'
        symbol = re.sub(r"_TSS$", "", itn)
        ensg = gene_map.get(symbol)
        if ensg:
            df.at[idx, "intended_target_name"] = ensg
        else:
            failures.append((idx, itn, symbol))

    if failures:
        print(f"  WARNING: {len(failures)} gene rows with no ENSG lookup:")
        for idx, itn, sym in failures[:10]:
            print(f"    intended_target_name={itn!r}  symbol={sym!r}")
    else:
        print(f"  All {n_rows} gene rows got ENSG IDs.")

    return df


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    # Load gene map
    if not os.path.exists(GTF_PATH):
        download_gtf(GTF_URL, GTF_PATH)
    print("Parsing GENCODE v43 GTF …", file=sys.stderr)
    gene_map = load_gene_map(GTF_PATH)
    print(f"  Loaded {len(gene_map):,} gene_name → ENSEMBL mappings", file=sys.stderr)

    # Load input
    print(f"Loading {INPUT_FILE} …")
    df = pd.read_csv(INPUT_FILE, sep="\t", dtype=str).fillna("")
    print(f"  {len(df):,} rows, {len(df.columns)} columns")

    # Step 1: fix targeting case
    print("\nStep 1: fix targeting case (TRUE/FALSE → True/False)")
    df, msg = fix_targeting_case(df)
    print(f"  {msg}")

    # Step 2: fill description from intended_target_name (BEFORE step 4 overwrites it)
    print("\nStep 2: fill description from intended_target_name")
    df = fill_description_from_intended_target_name(df)

    # Step 3: compute coords for gene rows (BEFORE step 4 renames intended_target_name)
    print("\nStep 3: compute intended_target_chr/start/end for gene rows")
    df = fix_coords_gene_rows(df)

    # Step 4: fix genomic_element + intended_target_name for gene rows
    print("\nStep 4: fix genomic_element ('gene'→'promoter') and intended_target_name (symbol→ENSG)")
    df = fix_gene_rows(df, gene_map)

    # Step 5: convert enhancer → distal element
    print("\nStep 5: convert genomic_element 'enhancer' → 'distal element'")
    n = (df["genomic_element"] == "enhancer").sum()
    df.loc[df["genomic_element"] == "enhancer", "genomic_element"] = "distal element"
    print(f"  {n} rows updated")

    # Reorder columns to canonical order
    extra = [c for c in df.columns if c not in COLUMN_ORDER]
    if extra:
        print(f"\n  Note: extra columns kept at end: {extra}")
    cols_out = [c for c in COLUMN_ORDER if c in df.columns] + extra
    df = df[cols_out]

    # Write output
    os.makedirs(os.path.dirname(OUTPUT_FILE) or ".", exist_ok=True)
    df.to_csv(OUTPUT_FILE, sep="\t", index=False, compression="gzip")
    size = os.path.getsize(OUTPUT_FILE)
    print(f"\nWritten: {OUTPUT_FILE} ({size:,} bytes, {len(df):,} rows)")


if __name__ == "__main__":
    main()

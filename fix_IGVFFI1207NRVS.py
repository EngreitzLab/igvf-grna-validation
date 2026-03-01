#!/usr/bin/env python3
"""
Fix script for IGVFFI1207NRVS (220-row CRISPRi screen — GATA1 enhancers + TSS controls).

Issues to fix:
  1. targeting case: TRUE/FALSE → True/False
  2. intended_target_name for 60 promoter (positive control) rows: coord strings → ENSG IDs
  3. intended_target_chr/start/end for promoter rows: per-guide coords → spanning coords
  4. description: absent column → fill from guide_id prefix

Enhancer rows (130 GATA1 tiling guides) are left untouched — their per-window coords
and intended_target_name are already correct.

Usage:
    python3 fix_IGVFFI1207NRVS.py
"""

import os
import re
import sys

import pandas as pd

# ── Imports from shared modules ────────────────────────────────────────────────
from fix_interactive import fix_targeting_case, fix_description_from_guide_id, COLUMN_ORDER
from validate_grna_file import load_gene_map, GTF_PATH, GTF_URL, download_gtf

INPUT_FILE  = "input/IGVFFI1207NRVS.tsv"
OUTPUT_FILE = "output/IGVFFI1207NRVS.tsv.gz"


# ── Fix helpers ────────────────────────────────────────────────────────────────

def fix_intended_target_name_promoter(df: pd.DataFrame, gene_map: dict) -> pd.DataFrame:
    """Replace coord-string intended_target_name with ENSG IDs for promoter rows.

    Mask: targeting==True AND genomic_element=='promoter'
    Gene symbol extracted from guide_id by stripping trailing _N then 'TSS.' prefix.
    """
    is_true = df["targeting"].str.upper() == "TRUE"
    mask = is_true & (df["genomic_element"] == "promoter")
    n_rows = mask.sum()

    failures = []
    for idx in df.index[mask]:
        guide_id = df.at[idx, "guide_id"]
        # TSS.CDKN2A_1  →  strip _1  →  TSS.CDKN2A  →  strip TSS.  →  CDKN2A
        prefix = re.sub(r"_\d+$", "", guide_id)
        symbol = re.sub(r"^TSS\.", "", prefix)

        ensg = gene_map.get(symbol)
        if ensg:
            df.at[idx, "intended_target_name"] = ensg
        else:
            failures.append((idx, guide_id, symbol))

    if failures:
        print(f"  WARNING: {len(failures)} promoter rows with no ENSG lookup:")
        for idx, gid, sym in failures[:10]:
            print(f"    guide_id={gid!r}  symbol={sym!r}")
    else:
        print(f"  All {n_rows} promoter rows got ENSG IDs.")

    return df


def fix_coords_promoter_only(df: pd.DataFrame) -> pd.DataFrame:
    """Recompute intended_target_chr/start/end for promoter rows only.

    Groups by guide_id prefix (strip trailing _N); assigns min(guide_start) /
    max(guide_end) for all rows in that group.  Enhancer rows are untouched.
    """
    is_true = df["targeting"].str.upper() == "TRUE"
    mask = is_true & (df["genomic_element"] == "promoter")

    tgt = df.loc[mask].copy()
    tgt["_prefix"] = tgt["guide_id"].str.replace(r"_\d+$", "", regex=True)
    tgt["_gs"] = pd.to_numeric(tgt["guide_start"], errors="coerce")
    tgt["_ge"] = pd.to_numeric(tgt["guide_end"],   errors="coerce")

    grp = tgt.groupby("_prefix").agg(
        _chr=("guide_chr", "first"),
        _start=("_gs", "min"),
        _end=("_ge", "max"),
    )

    prefix_col = df["guide_id"].str.replace(r"_\d+$", "", regex=True)
    for prefix, row in grp.iterrows():
        pmask = mask & (prefix_col == prefix)
        df.loc[pmask, "intended_target_chr"]   = row["_chr"]
        df.loc[pmask, "intended_target_start"] = int(row["_start"])
        df.loc[pmask, "intended_target_end"]   = int(row["_end"])

    print(f"  Recomputed coords for {len(grp)} promoter-prefix groups.")
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

    # Step 2: fix intended_target_name for promoter rows
    print("\nStep 2: fix intended_target_name for promoter rows (coord → ENSG)")
    df = fix_intended_target_name_promoter(df, gene_map)

    # Step 3: recompute coords for promoter rows only
    print("\nStep 3: recompute intended_target_chr/start/end for promoter rows")
    df = fix_coords_promoter_only(df)

    # Step 4: convert enhancer → distal element
    print("\nStep 4: convert genomic_element 'enhancer' → 'distal element'")
    n = (df["genomic_element"] == "enhancer").sum()
    df.loc[df["genomic_element"] == "enhancer", "genomic_element"] = "distal element"
    print(f"  {n} rows updated")

    # Step 5: fill description from guide_id
    print("\nStep 5: fill description column from guide_id prefix")
    df, msg = fix_description_from_guide_id(df)
    print(f"  {msg}")

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

#!/usr/bin/env python3
"""
Patch script for IGVFFI9754AGFB — corrects two issues in the existing output:

  1. targeting case: 'true'/'false' → 'True'/'False'
     The original fix_grna_file.py used str.lower(), producing lowercase.
     Spec requires title case.

  2. genomic_element for Enhancer-description rows with coord-string intended_target_name:
     Rows with description starting with 'Enhancer' fell through to the
     default 'promoter' classification in the original script. Phase 2 then
     correctly resolved their intended_target_name to a coordinate string,
     but left genomic_element='promoter' — a mismatch.
     Fix: 'promoter' → 'distal element' ONLY for rows where intended_target_name
     is already a coord string (chr:start-end).  Enhancer_*_and_GENE descriptions
     that resolved to ENSG IDs are dual-targeting (enhancer + TSS) and correctly
     remain as 'promoter'.

Known remaining issues (pre-existing, cannot be auto-resolved):
  - 46 rows for AKAP2, SMIM11B, LOC391322, LOC101927497 still use gene
    symbols as intended_target_name placeholder.

Input:  output/IGVFFI9754AGFB.tsv.gz  (existing broken output)
Output: output/IGVFFI9754AGFB.tsv.gz  (overwrite in place)

Usage:
    python3 fix_IGVFFI9754AGFB.py
"""

import os
import pandas as pd

FILE = "output/IGVFFI9754AGFB.tsv.gz"

print(f"Loading {FILE} …")
df = pd.read_csv(FILE, sep="\t", dtype=str).fillna("")
print(f"  {len(df):,} rows, {len(df.columns)} columns")

# ── Fix 1: targeting case ─────────────────────────────────────────────────────
print("\nFix 1: targeting case ('true'/'false' → 'True'/'False')")
mapping = {"true": "True", "false": "False", "True": "True", "False": "False"}
before = df["targeting"].value_counts().to_dict()
df["targeting"] = df["targeting"].map(mapping).fillna(df["targeting"])
after = df["targeting"].value_counts().to_dict()
print(f"  Before: {before}")
print(f"  After:  {after}")

# ── Fix 2: genomic_element for Enhancer-description rows ─────────────────────
print("\nFix 2: genomic_element 'promoter' → 'distal element' for Enhancer-description rows")
# Only convert rows whose intended_target_name is already a coord string.
# Enhancer_*_and_GENE descriptions that resolved to ENSG IDs are dual-targeting
# (enhancer + TSS gene) and should remain 'promoter'.
is_coord_itn = df["intended_target_name"].str.match(r"^chr\S+:\d+-\d+$")
mask = (df["genomic_element"] == "promoter") & df["description"].str.startswith("Enhancer") & is_coord_itn
n = mask.sum()
df.loc[mask, "genomic_element"] = "distal element"
print(f"  {n} rows updated")
print(f"  genomic_element value_counts now:")
print(df["genomic_element"].value_counts().to_string())

# ── Write output ──────────────────────────────────────────────────────────────
df.to_csv(FILE, sep="\t", index=False, compression="gzip")
size = os.path.getsize(FILE)
print(f"\nOverwritten: {FILE} ({size:,} bytes, {len(df):,} rows)")

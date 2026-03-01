#!/usr/bin/env python3
"""
Fix IGVFFI4634ZBZN and IGVFFI4575UMXX.

Fixes applied to both files:
  1. targeting: "TRUE" → "True", "FALSE" → "False" (title case per spec)
  2. genomic_element: empty → "promoter" for all targeting=True rows
  3. Add empty 'description' column if absent (optional per spec, avoids validator warning)
"""

import gzip
import os
import pandas as pd

FILES = [
    ("input/IGVFFI4634ZBZN.tsv", "output/IGVFFI4634ZBZN.tsv.gz"),
    ("input/IGVFFI4575UMXX.tsv", "output/IGVFFI4575UMXX.tsv.gz"),
]

# Files that need coordinate recomputation (target coords set to guide min/max per target group)
FILES_NEEDING_COORD_RECOMPUTE = {"IGVFFI4575UMXX"}


def recompute_target_coords(df):
    """
    For each targeting=True row, set intended_target_chr/start/end to the
    min/max of guide positions across all rows sharing the same intended_target_name.
    Overwrites existing values. Returns modified df.
    """
    is_true = df["targeting"] == "True"
    tgt = df.loc[is_true].copy()
    tgt["_gs"] = pd.to_numeric(tgt["guide_start"], errors="coerce")
    tgt["_ge"] = pd.to_numeric(tgt["guide_end"],   errors="coerce")

    grp = tgt.groupby("intended_target_name").agg(
        _chr  =("guide_chr",   "first"),
        _start=("_gs", "min"),
        _end  =("_ge", "max"),
    )

    for itn, row in grp.iterrows():
        mask = is_true & (df["intended_target_name"] == itn)
        df.loc[mask, "intended_target_chr"]   = row["_chr"]
        df.loc[mask, "intended_target_start"] = int(row["_start"])
        df.loc[mask, "intended_target_end"]   = int(row["_end"])
    return df, len(grp)

os.makedirs("output", exist_ok=True)

COLUMN_ORDER = [
    "guide_id", "spacer", "targeting", "type",
    "guide_chr", "guide_start", "guide_end", "strand", "pam",
    "genomic_element", "intended_target_name",
    "intended_target_chr", "intended_target_start", "intended_target_end",
    "putative_target_genes", "reporter", "imperfect", "description",
]

for input_path, output_path in FILES:
    file_id = os.path.splitext(os.path.basename(input_path))[0]
    print(f"\n{'='*60}")
    print(f"Input:  {input_path}")
    print(f"Output: {output_path}")

    df = pd.read_csv(input_path, sep="\t", dtype=str).fillna("")
    print(f"  Loaded {len(df):,} rows, {len(df.columns)} columns")

    # ── Fix 1: targeting case ─────────────────────────────────────────────────
    before = df["targeting"].value_counts().to_dict()
    df["targeting"] = df["targeting"].map(
        {"TRUE": "True", "True": "True", "FALSE": "False", "False": "False"}
    ).fillna(df["targeting"])
    after = df["targeting"].value_counts().to_dict()
    print(f"  targeting: {before} → {after}")

    # ── Fix 2: genomic_element → "promoter" for targeting=True rows ──────────
    is_targeting_true = df["targeting"] == "True"
    ge_empty = is_targeting_true & (df["genomic_element"] == "")
    n_ge = ge_empty.sum()
    df.loc[ge_empty, "genomic_element"] = "promoter"
    print(f"  genomic_element: filled 'promoter' for {n_ge:,} targeting=True rows")

    # ── Fix 3: add description column if absent ───────────────────────────────
    if "description" not in df.columns:
        df["description"] = ""
        print("  description: added empty column")

    # ── Fix 4 (IGVFFI4575UMXX only): recompute target coords from guide bounds ─
    if file_id in FILES_NEEDING_COORD_RECOMPUTE:
        # Capture pre-fix coords for worst-case reporting
        is_true = df["targeting"] == "True"
        tdf_before = df.loc[is_true].copy()
        tdf_before["_gs"] = pd.to_numeric(tdf_before["guide_start"], errors="coerce")
        tdf_before["_ge"] = pd.to_numeric(tdf_before["guide_end"],   errors="coerce")
        tdf_before["_ts"] = pd.to_numeric(tdf_before["intended_target_start"], errors="coerce")
        tdf_before["_te"] = pd.to_numeric(tdf_before["intended_target_end"],   errors="coerce")
        worst_start_over = (tdf_before["_ts"] - tdf_before["_gs"]).max()
        worst_end_under  = (tdf_before["_ge"] - tdf_before["_te"]).max()

        df, n_groups = recompute_target_coords(df)
        print(f"  recompute_target_coords: updated {n_groups:,} target groups; "
              f"worst pre-fix start overshoot={int(worst_start_over):,} bp, "
              f"end undershoot={int(worst_end_under):,} bp")

    # ── Reorder columns to canonical order ────────────────────────────────────
    extra = [c for c in df.columns if c not in COLUMN_ORDER]
    if extra:
        print(f"  WARNING: unexpected extra columns kept at end: {extra}")
    cols_out = [c for c in COLUMN_ORDER if c in df.columns] + extra
    df = df[cols_out]

    # ── Write output ──────────────────────────────────────────────────────────
    df.to_csv(output_path, sep="\t", index=False, compression="gzip")
    print(f"  Written: {output_path} ({os.path.getsize(output_path):,} bytes, {len(df):,} rows)")

print("\nDone.")

#!/usr/bin/env python3
"""
Universal interactive fixer for IGVF Per-Guide Metadata TSV files.

Reads a JSON problem report produced by:
    python3 validate_grna_file.py <file> --json-out <report.json>

Then walks through each fixable issue, prompts the user to apply or skip,
and writes the corrected file.

Usage:
    python3 fix_interactive.py <input_file> --problems <problems.json> [--output <output_file>]

If --output is not given, defaults to output/<basename>.tsv.gz
"""

import argparse
import json
import os
import re
import sys

import pandas as pd


# ── Constants ─────────────────────────────────────────────────────────────────

COLUMN_ORDER = [
    "guide_id", "spacer", "targeting", "type",
    "guide_chr", "guide_start", "guide_end", "strand", "pam",
    "genomic_element", "intended_target_name",
    "intended_target_chr", "intended_target_start", "intended_target_end",
    "putative_target_genes", "reporter", "imperfect", "description",
]

_NAN_STRINGS = {"", "nan", "NaN", "NAN", "na", "NA", "N/A", "n/a", "none", "None", "null", "NULL"}

def _is_nan_like(series: "pd.Series") -> "pd.Series":
    return series.isin(_NAN_STRINGS)


# ── Fix functions ─────────────────────────────────────────────────────────────

def fix_targeting_case(df):
    """Normalize TRUE/FALSE → True/False."""
    mapping = {"TRUE": "True", "True": "True", "FALSE": "False", "False": "False"}
    df["targeting"] = df["targeting"].map(mapping).fillna(df["targeting"])
    return df, "targeting normalized to title case"


def fix_targeting_consistency(df):
    """Set targeting=False for non-targeting/safe-targeting rows."""
    mask = df["type"].isin({"non-targeting", "safe-targeting"})
    n = mask.sum()
    df.loc[mask, "targeting"] = "False"
    return df, f"targeting set to 'False' for {n:,} non-targeting/safe-targeting rows"


def fix_genomic_element_promoter(df):
    """Fill empty genomic_element with 'promoter' for targeting=True rows."""
    is_true = df["targeting"].str.upper() == "TRUE"
    mask = is_true & _is_nan_like(df["genomic_element"])
    n = mask.sum()
    df.loc[mask, "genomic_element"] = "promoter"
    return df, f"genomic_element set to 'promoter' for {n:,} rows"


def fix_coords_recompute(df):
    """Recompute intended_target_chr/start/end from guide bounds per target group."""
    is_true = df["targeting"].str.upper() == "TRUE"
    tgt = df.loc[is_true].copy()
    tgt["_gs"] = pd.to_numeric(tgt["guide_start"], errors="coerce")
    tgt["_ge"] = pd.to_numeric(tgt["guide_end"], errors="coerce")
    grp = tgt.groupby("intended_target_name").agg(
        _chr=("guide_chr", "first"), _start=("_gs", "min"), _end=("_ge", "max")
    )
    for itn, row in grp.iterrows():
        mask = is_true & (df["intended_target_name"] == itn)
        df.loc[mask, "intended_target_chr"]   = row["_chr"]
        df.loc[mask, "intended_target_start"] = int(row["_start"])
        df.loc[mask, "intended_target_end"]   = int(row["_end"])
    return df, f"recomputed coordinates for {len(grp):,} target groups"


def fix_description_from_guide_id(df):
    """Fill description with guide group name (guide_id stripped of trailing _N)."""
    group_name = df["guide_id"].str.replace(r"_\d+$", "", regex=True)
    if "description" not in df.columns:
        df["description"] = group_name
    else:
        df["description"] = df["description"].where(
            ~_is_nan_like(df["description"]), group_name
        )
    n = (~_is_nan_like(df["description"])).sum()
    return df, f"description filled for {n:,} rows from guide_id"


# ── Registry ──────────────────────────────────────────────────────────────────

FIX_FUNCTIONS = {
    "fix_targeting_case":           fix_targeting_case,
    "fix_targeting_consistency":    fix_targeting_consistency,
    "fix_genomic_element_promoter": fix_genomic_element_promoter,
    "fix_coords_recompute":         fix_coords_recompute,
    "fix_description_from_guide_id": fix_description_from_guide_id,
}

FIX_DESCRIPTIONS = {
    "fix_targeting_case":           "Map TRUE→True / FALSE→False (title case per spec)",
    "fix_targeting_consistency":    "Set targeting='False' for non-targeting/safe-targeting rows",
    "fix_genomic_element_promoter": "Fill empty genomic_element with 'promoter' for targeting=True rows",
    "fix_coords_recompute":         "Recompute intended_target_chr/start/end from guide position min/max",
    "fix_description_from_guide_id": "Fill description with guide group name (guide_id stripped of _N suffix)",
}


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Interactively fix IGVF Per-Guide Metadata TSV files.")
    parser.add_argument("input_file", help="Input TSV or TSV.gz file")
    parser.add_argument("--problems", required=True, metavar="FILE",
                        help="JSON problem report from validate_grna_file.py --json-out")
    parser.add_argument("--output", metavar="FILE",
                        help="Output path (default: output/<basename>.tsv.gz)")
    args = parser.parse_args()

    # Determine output path
    if args.output:
        output_path = args.output
    else:
        basename = os.path.basename(args.input_file)
        # Strip .tsv or .tsv.gz
        stem = re.sub(r"\.tsv(\.gz)?$", "", basename)
        output_path = os.path.join("output", stem + ".tsv.gz")

    # Load problem report
    with open(args.problems) as f:
        report = json.load(f)
    issues = report["issues"]

    print(f"\n{'═'*70}")
    print(f"  Interactive Fixer")
    print(f"  Input:    {args.input_file}")
    print(f"  Problems: {args.problems}")
    print(f"  Output:   {output_path}")
    print(f"{'═'*70}\n")

    # Load data
    print(f"Loading {args.input_file} …", file=sys.stderr)
    df = pd.read_csv(args.input_file, sep="\t", dtype=str).fillna("")
    print(f"  {len(df):,} rows, {len(df.columns)} columns\n", file=sys.stderr)

    # Partition issues
    fixable   = [iss for iss in issues if iss.get("fix_type") in FIX_FUNCTIONS]
    unfixable = [iss for iss in issues if iss.get("fix_type") not in FIX_FUNCTIONS]

    # Deduplicate fixable by fix_type so the same automated fix isn't offered twice
    seen_fix_types = set()
    unique_fixable = []
    for iss in fixable:
        ft = iss["fix_type"]
        if ft not in seen_fix_types:
            seen_fix_types.add(ft)
            unique_fixable.append(iss)
    fixable = unique_fixable

    # Show unfixable issues (informational only)
    if unfixable:
        print(f"── Issues with no automated fix ({len(unfixable)}) " + "─" * 30)
        for iss in unfixable:
            tag = "ERROR" if iss["severity"] == "error" else "WARN "
            print(f"  [{tag}] {iss['field']}: {iss['message']}")
            print(f"         → No automated fix available. Manual review required.")
        print()

    if not fixable:
        print("No automatically fixable issues found. Nothing to do.\n")
        sys.exit(0)

    # Interactive loop
    applied_fixes = []
    apply_all = False

    for i, iss in enumerate(fixable, 1):
        print(f"══ Issue {i}/{len(fixable)} " + "═" * 50)
        print(f"  Field:    {iss['field']}")
        print(f"  Severity: {iss['severity'].upper()}")
        print(f"  Problem:  {iss['message']}")
        print(f"  Fix:      {FIX_DESCRIPTIONS[iss['fix_type']]}")

        if apply_all:
            choice = "y"
            print(f"  [auto-applying all]")
        else:
            choice = input("  Apply this fix? [y]es / [n]o / [a]ll / [q]uit: ").strip().lower()

        if choice == "q":
            print("\nQuitting without saving.")
            sys.exit(0)
        elif choice == "a":
            apply_all = True
            choice = "y"

        if choice == "y":
            df, summary = FIX_FUNCTIONS[iss["fix_type"]](df)
            print(f"  ✓ Applied: {summary}")
            applied_fixes.append(iss["fix_type"])
        else:
            print(f"  – Skipped.")
        print()

    if not applied_fixes:
        print("No fixes were applied. Output file not written.\n")
        sys.exit(0)

    # Reorder columns to canonical order
    extra = [c for c in df.columns if c not in COLUMN_ORDER]
    if extra:
        print(f"  Note: extra columns kept at end: {extra}")
    cols_out = [c for c in COLUMN_ORDER if c in df.columns] + extra
    df = df[cols_out]

    # Write output
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False, compression="gzip")
    size = os.path.getsize(output_path)
    print(f"{'═'*70}")
    print(f"  Written: {output_path} ({size:,} bytes, {len(df):,} rows)")
    print(f"  Applied fixes: {', '.join(applied_fixes)}")
    print(f"{'═'*70}\n")


if __name__ == "__main__":
    main()

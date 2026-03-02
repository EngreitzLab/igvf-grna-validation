#!/usr/bin/env python3
"""
Fix script for IGVFFI0580WJFK (16,201-row WTC11 CRISPRi random screen).

Data breakdown:
  15,401 targeting rows
     400 safe-targeting controls
     400 non-targeting controls

  Of the 15,401 targeting rows:
    1,231  genomic_element='promoter'  (83 unique windows, need ENSG lookup)
   13,493  genomic_element='distal_element'  (underscore → space)
      677  genomic_element empty + intended_target_name empty  (flagged, not fixed)

Issues fixed:
  1. targeting case: TRUE/FALSE → True/False
  2. genomic_element: 'distal_element' (underscore) → 'distal element' (space)
  3. description: absent column → copy from intended_target_name (where non-empty)
  4. intended_target_name for promoter rows: coord strings → ENSG IDs
     via GENCODE v43 TSS intersection

Issue flagged but NOT fixed:
  5. 677 targeting rows with empty genomic_element + empty intended_target_name

Usage:
    python3 fix_IGVFFI0580WJFK.py
"""

import bisect
import gzip
import os
import re
import sys

import pandas as pd

# ── Imports from shared modules ────────────────────────────────────────────────
from fix_interactive import fix_targeting_case, COLUMN_ORDER, _is_nan_like
from validate_grna_file import GTF_PATH, GTF_URL, download_gtf

INPUT_FILE  = "input/IGVFFI0580WJFK.tsv"
OUTPUT_FILE = "output/IGVFFI0580WJFK.tsv.gz"


# ── TSS map ───────────────────────────────────────────────────────────────────

def load_tss_map(gtf_path: str) -> dict:
    """Build {chr_str: sorted [(tss0, ensg_no_version, gene_name), ...]} from GTF.

    Uses 'gene' features only (one TSS per gene, avoids multi-transcript redundancy).
    TSS is computed in 0-based coordinates:
      + strand: tss0 = gene_start - 1   (GTF is 1-based inclusive)
      - strand: tss0 = gene_end   - 1
    """
    tss_map: dict[str, list] = {}
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chrom  = fields[0]
            start1 = int(fields[3])   # 1-based inclusive start
            end1   = int(fields[4])   # 1-based inclusive end
            strand = fields[6]
            attrs  = fields[8]

            gid_m = re.search(r'gene_id "([^"]+)"', attrs)
            gnm_m = re.search(r'gene_name "([^"]+)"', attrs)
            if not (gid_m and gnm_m):
                continue

            ensg = gid_m.group(1).split(".")[0]
            name = gnm_m.group(1)
            tss0 = (start1 - 1) if strand == "+" else (end1 - 1)

            if chrom not in tss_map:
                tss_map[chrom] = []
            tss_map[chrom].append((tss0, ensg, name))

    # Sort each chromosome list by tss0 for bisect
    for chrom in tss_map:
        tss_map[chrom].sort(key=lambda x: x[0])

    return tss_map


# ── Promoter ENSG lookup via TSS intersection ──────────────────────────────────

def lookup_promoter_ensg(df: pd.DataFrame, tss_map: dict) -> pd.DataFrame:
    """Replace coord-string intended_target_name with ENSG IDs for promoter rows.

    Mask: genomic_element == 'promoter'
    For each unique window (coord string):
      - 0 TSS overlaps  → FAILURE (leave coord string, warn)
      - 1 TSS overlap   → assign ENSG ID
      - >1 TSS overlaps → AMBIGUOUS: assign ENSG of gene with TSS closest to
                          window center; print all matches for review

    Coordinate system:
      intended_target_name uses 0-based half-open BED: 'chrN:w_start-w_end'
      TSS overlap: w_start <= tss0 < w_end
    """
    mask_promoter = df["genomic_element"] == "promoter"
    n_rows = mask_promoter.sum()

    # Collect unique windows
    windows = df.loc[mask_promoter, "intended_target_name"].unique()
    print(f"  {n_rows} promoter rows, {len(windows)} unique windows")

    resolved   = {}   # window → ensg
    failures   = []
    ambiguous  = []

    coord_re = re.compile(r"^(chr\S+):(\d+)-(\d+)$")

    for win in windows:
        m = coord_re.match(win)
        if not m:
            failures.append((win, "could not parse coord string"))
            continue

        chrom, w_start, w_end = m.group(1), int(m.group(2)), int(m.group(3))
        center = (w_start + w_end) / 2.0

        chrom_list = tss_map.get(chrom, [])
        if not chrom_list:
            failures.append((win, f"chromosome {chrom!r} not in GTF"))
            continue

        # Use bisect to find candidates in [w_start, w_end)
        keys = [t[0] for t in chrom_list]
        lo = bisect.bisect_left(keys, w_start)
        hi = bisect.bisect_left(keys, w_end)
        hits = chrom_list[lo:hi]   # all TSSs with w_start <= tss0 < w_end

        if len(hits) == 0:
            failures.append((win, "no TSS in window"))
        elif len(hits) == 1:
            resolved[win] = hits[0][1]   # ensg
        else:
            # Ambiguous: pick closest to center
            best = min(hits, key=lambda t: abs(t[0] - center))
            resolved[win] = best[1]
            ambiguous.append((win, hits, best))

    # Print summary
    n_res  = len(resolved)
    n_fail = len(failures)
    n_amb  = len(ambiguous)
    print(f"  Resolved: {n_res}   Failures: {n_fail}   Ambiguous (assigned closest): {n_amb}")

    if failures:
        print(f"\n  FAILURE windows ({n_fail}) — coord string left unchanged:")
        for win, reason in failures:
            print(f"    {win}  [{reason}]")

    if ambiguous:
        print(f"\n  AMBIGUOUS windows ({n_amb}) — assigned gene with TSS closest to window center:")
        for win, hits, best in ambiguous:
            all_str = ", ".join(f"{nm} ({ensg}, tss0={tss0})" for tss0, ensg, nm in hits)
            print(f"    {win}")
            print(f"      All hits: {all_str}")
            print(f"      Assigned: {best[2]} ({best[1]})")

    # Apply resolved mappings
    for win, ensg in resolved.items():
        row_mask = mask_promoter & (df["intended_target_name"] == win)
        df.loc[row_mask, "intended_target_name"] = ensg

    return df


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    # Load TSS map
    if not os.path.exists(GTF_PATH):
        download_gtf(GTF_URL, GTF_PATH)
    print("Parsing GENCODE v43 GTF for TSS positions …", file=sys.stderr)
    tss_map = load_tss_map(GTF_PATH)
    n_genes = sum(len(v) for v in tss_map.values())
    print(f"  Loaded {n_genes:,} gene TSSs across {len(tss_map):,} chromosomes",
          file=sys.stderr)

    # Load input
    print(f"\nLoading {INPUT_FILE} …")
    df = pd.read_csv(INPUT_FILE, sep="\t", dtype=str).fillna("")
    print(f"  {len(df):,} rows, {len(df.columns)} columns")

    # ── Step 1: fix targeting case ─────────────────────────────────────────────
    print("\nStep 1: fix targeting case (TRUE/FALSE → True/False)")
    df, msg = fix_targeting_case(df)
    print(f"  {msg}")

    # ── Step 2: fix genomic_element underscore → space ─────────────────────────
    print("\nStep 2: fix genomic_element 'distal_element' → 'distal element'")
    n_de = (df["genomic_element"] == "distal_element").sum()
    df["genomic_element"] = df["genomic_element"].str.replace(
        "distal_element", "distal element", regex=False
    )
    print(f"  {n_de} rows updated")

    # ── Step 3: fill description from intended_target_name ────────────────────
    # Must run BEFORE step 4, which overwrites intended_target_name with ENSG IDs.
    print("\nStep 3: fill description column from intended_target_name")
    if "description" not in df.columns:
        df["description"] = ""
    mask_desc = _is_nan_like(df["description"]) & ~_is_nan_like(df["intended_target_name"])
    n_desc = mask_desc.sum()
    df.loc[mask_desc, "description"] = df.loc[mask_desc, "intended_target_name"]
    print(f"  {n_desc} rows filled")

    # ── Step 4: TSS intersection for promoter rows ────────────────────────────
    print("\nStep 4: TSS intersection for promoter rows (coord string → ENSG ID)")
    df = lookup_promoter_ensg(df, tss_map)

    # ── Step 5: flag the 677 empty-target rows ────────────────────────────────
    is_true = df["targeting"].str.upper() == "TRUE"
    empty_mask = (
        is_true
        & _is_nan_like(df["genomic_element"])
        & _is_nan_like(df["intended_target_name"])
    )
    n_empty = empty_mask.sum()
    if n_empty:
        empty_rows = df.loc[empty_mask]
        chrs  = empty_rows["guide_chr"].value_counts().to_dict()
        gids  = empty_rows["guide_id"].tolist()
        print(f"\nStep 5: WARNING — {n_empty} targeting rows with empty genomic_element "
              f"AND empty intended_target_name (cannot determine target; left as-is)")
        print(f"  Chromosomes: {chrs}")
        print(f"  guide_id range: {gids[0]} … {gids[-1]}")
    else:
        print("\nStep 5: No empty-target targeting rows found.")

    # ── Reorder columns and write output ──────────────────────────────────────
    extra = [c for c in df.columns if c not in COLUMN_ORDER]
    if extra:
        print(f"\n  Note: extra columns kept at end: {extra}")
    cols_out = [c for c in COLUMN_ORDER if c in df.columns] + extra
    df = df[cols_out]

    os.makedirs(os.path.dirname(OUTPUT_FILE) or ".", exist_ok=True)
    df.to_csv(OUTPUT_FILE, sep="\t", index=False, compression="gzip")
    size = os.path.getsize(OUTPUT_FILE)
    print(f"\nWritten: {OUTPUT_FILE} ({size:,} bytes, {len(df):,} rows)")


if __name__ == "__main__":
    main()

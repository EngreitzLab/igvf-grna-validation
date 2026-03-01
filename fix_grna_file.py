#!/usr/bin/env python3
"""
Fix gRNA metadata file to IGVF Per-Guide Metadata Submission Format spec.
Reads: input/IGVFFI9754AGFB.tsv.gz
Writes: output/IGVFFI9754AGFB.tsv.gz
"""

import gzip
import os
import re
import sys
import time

import pandas as pd
import requests


# ── Constants ────────────────────────────────────────────────────────────────

INPUT_FILE  = "input/IGVFFI9754AGFB.tsv.gz"
OUTPUT_FILE = "output/IGVFFI9754AGFB.tsv.gz"
GTF_PATH    = "input/IGVFFI9573KOZR.gtf.gz"
GTF_URL     = ("https://api.data.igvf.org/reference-files/IGVFFI9573KOZR/"
               "@@download/IGVFFI9573KOZR.gtf.gz")


# ── Step 1: Download & parse GENCODE v43 GTF ─────────────────────────────────

def download_gtf(url: str, dest: str) -> None:
    print(f"Downloading GENCODE GTF to {dest} …")
    r = requests.get(url, stream=True, timeout=120)
    r.raise_for_status()
    with open(dest, "wb") as fh:
        for chunk in r.iter_content(chunk_size=65536):
            fh.write(chunk)
    print("Download complete.")


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
                ensembl_id = gid.group(1).split(".")[0]   # drop version suffix
                gene_map[gnm.group(1)] = ensembl_id
    return gene_map


if not os.path.exists(GTF_PATH):
    download_gtf(GTF_URL, GTF_PATH)

print("Parsing GENCODE v43 GTF …")
gene_map = load_gene_map(GTF_PATH)
print(f"  Loaded {len(gene_map):,} gene_name → ENSEMBL mappings")


# ── Step 2: Load input file ───────────────────────────────────────────────────

print(f"Reading {INPUT_FILE} …")
df = pd.read_csv(INPUT_FILE, sep="\t", dtype=str).fillna("")
print(f"  {len(df):,} rows loaded")


# ── Step 3: Compute per-target coordinate bounds ──────────────────────────────
#
#  Group only targeting rows by their exact `description` value.
#  target_chr  = first guide_chr in group (should be uniform)
#  target_start = min(guide_start) across group
#  target_end   = max(guide_end)   across group

targeting_mask = df["type"] == "targeting"
tdf = df[targeting_mask].copy()
tdf["_gs"] = pd.to_numeric(tdf["guide_start"], errors="coerce")
tdf["_ge"] = pd.to_numeric(tdf["guide_end"],   errors="coerce")

coord_bounds = (
    tdf.groupby("description", sort=False)
       .agg(target_chr=("guide_chr", "first"),
            target_start=("_gs", "min"),
            target_end=("_ge", "max"))
       .reset_index()
)
coord_bounds["target_start"] = coord_bounds["target_start"].astype(int)
coord_bounds["target_end"]   = coord_bounds["target_end"].astype(int)
coord_bounds_dict = {
    row["description"]: {
        "chr":   row["target_chr"],
        "start": row["target_start"],
        "end":   row["target_end"],
    }
    for _, row in coord_bounds.iterrows()
}


# ── Step 4a: Fix `targeting` boolean (all rows) ───────────────────────────────

df["targeting"] = df["targeting"].str.lower()                    # TRUE→true, FALSE→false
df.loc[df["type"] == "safe-targeting", "targeting"] = "false"    # safe-targeting must be false


# ── Step 4b: Set `genomic_element` ───────────────────────────────────────────
#
#  Rules (only for type=="targeting"; others left blank):
#   • description contains "TSS" (case-sensitive) → "promoter"
#   • guide_id starts with "9p21_add_" AND no TSS   → "distal element"
#   • all other targeting guides                    → "promoter"

def get_genomic_element(row):
    if row["type"] != "targeting":
        return ""
    if "TSS" in row["description"]:
        return "promoter"
    if row["guide_id"].startswith("9p21_add_"):
        return "distal element"
    return "promoter"

df["genomic_element"] = df.apply(get_genomic_element, axis=1)


# ── Step 4c: Set `intended_target_chr/start/end` ─────────────────────────────

def fill_target_coords(df_in, coord_bounds_dict):
    chr_col   = df_in["description"].map(lambda d: coord_bounds_dict.get(d, {}).get("chr",   ""))
    start_col = df_in["description"].map(lambda d: coord_bounds_dict.get(d, {}).get("start", ""))
    end_col   = df_in["description"].map(lambda d: coord_bounds_dict.get(d, {}).get("end",   ""))
    # Convert integers to strings (non-targeting rows get "" from the dict miss)
    start_col = start_col.map(lambda v: str(v) if v != "" else "")
    end_col   = end_col.map(lambda v: str(v) if v != "" else "")
    return chr_col, start_col, end_col

chr_col, start_col, end_col = fill_target_coords(df[targeting_mask], coord_bounds_dict)
df["intended_target_chr"]   = ""
df["intended_target_start"] = ""
df["intended_target_end"]   = ""
df.loc[targeting_mask, "intended_target_chr"]   = chr_col.values
df.loc[targeting_mask, "intended_target_start"] = start_col.values
df.loc[targeting_mask, "intended_target_end"]   = end_col.values


# ── Step 4d: Set `intended_target_name` ──────────────────────────────────────
#
#  Promoter guides  → ENSEMBL gene ID (fall back to gene symbol on miss)
#  Distal-element   → "chr:start-end" coordinate string
#  Non-targeting / safe-targeting → ""
#
#  Gene symbol extraction from description:
#   1. Remove all _TSS\d* occurrences (handles TSS2, TSS, etc.)
#   2. Split on _and_; try candidates in order: first, last, middle tokens
#   3. For each candidate also try stripping a trailing _alt suffix
#   4. Return first token found in gene_map (GENCODE v43)

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


lookup_failures = {}   # desc → fallback symbol used

def get_intended_target_name(row):
    if row["type"] != "targeting":
        return ""
    ge   = row["genomic_element"]
    desc = row["description"]
    if ge == "distal element":
        bounds = coord_bounds_dict.get(desc)
        if bounds:
            return f"{bounds['chr']}:{bounds['start']}-{bounds['end']}"
        return ""
    # promoter
    candidates = extract_gene_candidates(desc)
    for candidate in candidates:
        if candidate in gene_map:
            return gene_map[candidate]
    # All candidates failed → report and fall back to first candidate
    fallback = candidates[0] if candidates else desc
    lookup_failures[desc] = fallback
    return fallback

df["intended_target_name"] = df.apply(get_intended_target_name, axis=1)


# ── Step 4e: Phase 2 – resolve remaining lookup_failures ─────────────────────

# Pass 1: Mouse gene capitalization (e.g. Cmip → CMIP → ENSG…)
mouse_resolved = {}
for desc, fallback in list(lookup_failures.items()):
    if len(fallback) > 1 and fallback[0].isupper() and fallback[1].islower():
        upper = fallback.upper()
        if upper in gene_map:
            mouse_resolved[desc] = gene_map[upper]
            del lookup_failures[desc]
print(f"Mouse gene resolution: {len(mouse_resolved)} descriptions resolved")

# Pass 2: Enhancer descriptions → coordinate string from guide bounds
enhancer_resolved = {}
for desc, fallback in list(lookup_failures.items()):
    if desc.startswith("Enhancer"):
        bounds = coord_bounds_dict.get(desc)
        enhancer_resolved[desc] = (
            f"{bounds['chr']}:{bounds['start']}-{bounds['end']}" if bounds else fallback
        )
        del lookup_failures[desc]
print(f"Enhancer description resolution: {len(enhancer_resolved)} descriptions resolved")

# Pass 3: Ensembl REST API for renamed human genes
def ensembl_batch_lookup(symbols):
    url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"
    r = requests.post(url, json={"symbols": symbols},
                      headers={"Content-Type": "application/json",
                               "Accept": "application/json"}, timeout=60)
    if r.status_code == 200:
        return {sym: info["id"] for sym, info in r.json().items()
                if info and isinstance(info, dict) and "id" in info}
    return {}

def ensembl_xref_lookup(symbol):
    url = f"https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{symbol}"
    r = requests.get(url, headers={"Accept": "application/json"},
                     params={"object_type": "gene"}, timeout=30)
    if r.status_code == 200:
        for entry in r.json():
            if entry.get("type") == "gene":
                return entry["id"]
    return None

unique_syms = list(set(lookup_failures.values()))
print(f"Querying Ensembl REST API for {len(unique_syms)} unique symbols …")
batch = ensembl_batch_lookup(unique_syms)
print(f"  Batch POST resolved: {len(batch)}")
xref = {}
for sym in [s for s in unique_syms if s not in batch]:
    result = ensembl_xref_lookup(sym)
    if result:
        xref[sym] = result
    time.sleep(0.07)
print(f"  Xref GET resolved:   {len(xref)}")
ensembl_sym_map = {**batch, **xref}

ensembl_resolved = {}
final_failures = {}
for desc, fallback in lookup_failures.items():
    if fallback in ensembl_sym_map:
        ensembl_resolved[desc] = ensembl_sym_map[fallback]
    else:
        final_failures[desc] = fallback
print(f"Ensembl API resolution: {len(ensembl_resolved)} descriptions resolved")

# Apply all Phase 2 resolutions to the dataframe
for resolutions in [mouse_resolved, enhancer_resolved, ensembl_resolved]:
    for desc, name in resolutions.items():
        mask = (df["description"] == desc) & (df["type"] == "targeting")
        df.loc[mask, "intended_target_name"] = name


# ── Step 5: Write output ──────────────────────────────────────────────────────

os.makedirs("output", exist_ok=True)
df.to_csv(OUTPUT_FILE, sep="\t", index=False, compression="gzip")
print(f"\nOutput written → {OUTPUT_FILE}")


# ── Summary ───────────────────────────────────────────────────────────────────

print("\n=== Summary of changes ===")
print(f"  Total rows:                {len(df):,}")
print(f"  targeting = true:          {(df['targeting'] == 'true').sum():,}")
print(f"  targeting = false:         {(df['targeting'] == 'false').sum():,}")
print(f"  genomic_element=promoter:  {(df['genomic_element'] == 'promoter').sum():,}")
print(f"  genomic_element=distal:    {(df['genomic_element'] == 'distal element').sum():,}")
tm = targeting_mask
print(f"  targeting rows with intended_target_name:  "
      f"{(tm & (df['intended_target_name'] != '')).sum():,}")
print(f"  targeting rows with intended_target_chr:   "
      f"{(tm & (df['intended_target_chr'] != '')).sum():,}")

print(f"\n  Phase 2 resolution summary:")
print(f"    Mouse gene capitalization:  {len(mouse_resolved):3d} descriptions resolved")
print(f"    Enhancer coordinate string: {len(enhancer_resolved):3d} descriptions resolved")
print(f"    Ensembl REST API:           {len(ensembl_resolved):3d} descriptions resolved")

if final_failures:
    print(f"\n=== Remaining lookup failures ({len(final_failures)} unique descriptions) ===")
    print("  These rows still use the gene symbol as a placeholder.")
    print("  Manual review required:\n")
    for desc, fallback in sorted(final_failures.items()):
        print(f"  WARNING  description={desc!r:60s}  fallback={fallback!r}")
else:
    print("\n  All gene symbols resolved to ENSEMBL IDs successfully.")

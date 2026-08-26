#!/usr/bin/env python3
"""
Correct the 300-gene promoter guide library (300genes_guide_metadata_v43_promoter_fixed.tsv).

Three corrections, in order (each depends on the previous):

  1. LIFT   — three guides sit on alt-haplotype contigs (validator check T6). Each spacer
              has exactly one exact match on the primary assembly, established by
              verify_guide_coords_against_genome.py, so replace the alt placement with the
              primary one. Lifted coordinates are written spacer-only (0-based half-open,
              PAM excluded) per spec.

  2. WINDOW — every one of the library's 301 element windows equals the exact span of the
              guides in its `description` group. Lifting moves two CSNK2B guides onto
              chr6, so that group's window is recomputed under the same rule.

  3. ENSG   — spec p.3: for genomic_element == "promoter", intended_target_name must be the
              ENSEMBL id of the regulated gene. The file carries gene symbols. Map each to
              its GENCODE v43 gene id, resolving symbol collisions to the primary-assembly
              copy (every collision here is a primary gene shadowed by *_alt duplicates).

Usage:
    python3 fix_300genes_promoter_library.py <input.tsv> \
        --gtf input/IGVFFI9573KOZR.gtf.gz \
        --coord-verification problems/coord_verification.json \
        --out output/300genes_guide_metadata_v43_corrected.tsv
"""

import argparse
import gzip
import json
import re
import sys

import pandas as pd

PRIMARY_CONTIG_RE = re.compile(r"^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$")


def load_primary_gene_ids(gtf_path: str) -> dict:
    """Map gene_name → unversioned gene_id, keeping only primary-assembly gene records.

    GENCODE lists a separate gene record per alt haplotype, so a symbol like CSNK2B maps to
    seven ids; only one lives on the primary assembly. Restricting to primary contigs makes
    the mapping unique without needing per-row coordinate arbitration.
    """
    mapping, dropped = {}, {}
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if f[2] != "gene":
                continue
            name = re.search(r'gene_name "([^"]+)"', f[8])
            gid = re.search(r'gene_id "([^"]+)"', f[8])
            if not (name and gid):
                continue
            name, gid = name.group(1), gid.group(1).split(".")[0]
            if PRIMARY_CONTIG_RE.match(f[0]):
                mapping.setdefault(name, set()).add(gid)
            else:
                dropped.setdefault(name, set()).add(gid)

    ambiguous = {k: sorted(v) for k, v in mapping.items() if len(v) > 1}
    return {k: next(iter(v)) for k, v in mapping.items() if len(v) == 1}, ambiguous, dropped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("guides")
    ap.add_argument("--gtf", default="input/IGVFFI9573KOZR.gtf.gz")
    ap.add_argument("--coord-verification", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--changelog", default=None)
    args = ap.parse_args()

    df = pd.read_csv(args.guides, sep="\t", dtype=str)
    original = df.copy()
    log = []

    def note(msg):
        log.append(msg)
        print(msg)

    note(f"Read {len(df):,} rows from {args.guides}")

    # ── 1. LIFT alt-contig guides to primary assembly ─────────────────────────
    ver = json.load(open(args.coord_verification))
    lift_hits, alt_spacers = ver["lift_hits"], ver["alt_spacers"]

    note("\n── 1. Lift alt-contig guides to primary assembly ──")
    lifted_ids = []
    for gid, hits in lift_hits.items():
        if len(hits) != 1:
            note(f"  SKIP {gid}: {len(hits)} primary matches (need exactly 1)")
            continue
        h = hits[0]
        idx = df.index[df.guide_id == gid]
        assert len(idx) == 1, f"{gid} matched {len(idx)} rows"
        i = idx[0]
        row = df.loc[i]
        # the spacer's orientation at the primary locus must agree with the recorded strand
        assert row.strand == h["strand"], (
            f"{gid}: recorded strand {row.strand} != primary-match strand {h['strand']}")
        assert h["end"] - h["start"] == len(alt_spacers[gid]), "lift span != spacer length"

        note(f"  {gid}: {row.guide_chr}:{row.guide_start}-{row.guide_end} "
             f"→ {h['contig']}:{h['start']}-{h['end']} ({h['strand']}) "
             f"[spacer-only, {h['end'] - h['start']} bp]")
        df.loc[i, "guide_chr"]   = h["contig"]
        df.loc[i, "guide_start"] = str(h["start"])
        df.loc[i, "guide_end"]   = str(h["end"])
        # intended_target_chr follows the guide, but only where it was already populated
        if pd.notna(row.intended_target_chr) and str(row.intended_target_chr).strip():
            df.loc[i, "intended_target_chr"] = h["contig"]
        lifted_ids.append(gid)

    # The two CSNK2B guides join the existing CSNK2B promoter group: their description must
    # match so they share an element, and the alt-contig provenance moves into this log.
    csnk2b_lifted = [g for g in lifted_ids if g.startswith("CSNK2B_")]
    if csnk2b_lifted:
        target_desc = "CSNK2B_P1P2"
        for gid in csnk2b_lifted:
            i = df.index[df.guide_id == gid][0]
            note(f"  {gid}: description {df.loc[i, 'description']!r} → {target_desc!r} "
                 f"(joins the CSNK2B promoter element)")
            df.loc[i, "description"] = target_desc
            df.loc[i, "intended_target_name"] = "CSNK2B"

    # ── 2. Recompute element windows for groups touched by the lift ───────────
    note("\n── 2. Recompute element windows (window == span of guides in the group) ──")
    is_targeting = df.type == "targeting"
    tdf = df[is_targeting].copy()
    tdf["gstart"] = tdf.guide_start.astype(int)
    tdf["gend"] = tdf.guide_end.astype(int)

    touched_groups = set(
        df.loc[df.guide_id.isin(lifted_ids) & is_targeting, "description"])
    for desc in sorted(touched_groups):
        sub = tdf[tdf.description == desc]
        chrs = set(sub.guide_chr)
        assert len(chrs) == 1, f"{desc} spans multiple contigs: {chrs}"
        new_start, new_end = int(sub.gstart.min()), int(sub.gend.max())
        mask = (df.description == desc) & is_targeting
        old = (df.loc[mask, "intended_target_start"].iloc[0],
               df.loc[mask, "intended_target_end"].iloc[0])
        note(f"  {desc} ({mask.sum()} guides on {chrs.pop()}): "
             f"{old[0]}-{old[1]} → {new_start}-{new_end}")
        df.loc[mask, "intended_target_start"] = str(new_start)
        df.loc[mask, "intended_target_end"]   = str(new_end)

    # Assert the window rule now holds for every group, not just the touched ones
    chk = df[df.type == "targeting"].copy()
    chk["gstart"] = chk.guide_start.astype(int)
    chk["gend"] = chk.guide_end.astype(int)
    agg = chk.groupby("description").agg(
        ws=("intended_target_start", "first"), we=("intended_target_end", "first"),
        gmin=("gstart", "min"), gmax=("gend", "max"))
    violations = agg[(agg.ws.astype(int) != agg.gmin) | (agg.we.astype(int) != agg.gmax)]
    note(f"  window rule holds for {len(agg) - len(violations)}/{len(agg)} groups")
    if len(violations):
        note(f"  VIOLATIONS:\n{violations.to_string()}")

    # ── 3. Gene symbol → GENCODE v43 ENSG for promoter rows ──────────────────
    note("\n── 3. intended_target_name: gene symbol → GENCODE v43 ENSG ──")
    gene_map, ambiguous, alt_only = load_primary_gene_ids(args.gtf)
    note(f"  loaded {len(gene_map):,} unique primary-assembly gene_name → ENSG mappings")
    if ambiguous:
        note(f"  {len(ambiguous)} symbols still ambiguous on the primary assembly "
             f"(not used unless referenced): {list(ambiguous)[:5]}")

    needs_ensg = (df.type == "targeting") & (df.genomic_element == "promoter")
    symbols = sorted(set(df.loc[needs_ensg, "intended_target_name"]))
    unresolved = [s for s in symbols if s not in gene_map]
    collided = [s for s in symbols if s in ambiguous]
    note(f"  {needs_ensg.sum():,} promoter rows across {len(symbols)} distinct symbols")
    note(f"  resolved: {len(symbols) - len(unresolved)}   unresolved: {len(unresolved)}   "
         f"primary-assembly collisions: {len(collided)}")
    if unresolved:
        note(f"  UNRESOLVED (left unchanged): {unresolved}")
    if collided:
        note(f"  COLLIDED (left unchanged): {collided}")

    shadowed = [s for s in symbols if s in alt_only and s in gene_map]
    note(f"  {len(shadowed)} symbols had *_alt duplicate gene records, resolved to the "
         f"primary copy: {shadowed[:12]}")

    resolvable = needs_ensg & df.intended_target_name.isin(gene_map)
    df.loc[resolvable, "intended_target_name"] = (
        df.loc[resolvable, "intended_target_name"].map(gene_map))
    note(f"  rewrote {resolvable.sum():,} rows to ENSG ids")

    # ── Write ────────────────────────────────────────────────────────────────
    changed = (original.fillna("__NA__") != df.fillna("__NA__")).any(axis=1).sum()
    note(f"\n{changed:,} of {len(df):,} rows changed")
    for col in df.columns:
        n = (original[col].fillna("__NA__") != df[col].fillna("__NA__")).sum()
        if n:
            note(f"  {col}: {n:,} cells changed")

    df.to_csv(args.out, sep="\t", index=False)
    note(f"\nWrote {args.out}")

    if args.changelog:
        with open(args.changelog, "w") as fh:
            fh.write("\n".join(log) + "\n")
        print(f"Wrote {args.changelog}", file=sys.stderr)


if __name__ == "__main__":
    main()

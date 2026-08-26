#!/usr/bin/env python3
"""
Verify every positioned guide's spacer against the GRCh38 reference, and locate
primary-assembly matches for guides that sit on non-primary contigs (check T6).

Streams the gzipped reference one contig at a time so peak memory stays at a single
chromosome, except for the contigs named in --lift-to, which are held in memory to
search for the alt-contig spacers.

Usage:
    python3 verify_guide_coords_against_genome.py <guides.tsv> \
        --fasta input/GRCh38_ucsc_IGVFFI6815WBWB.fasta.gz \
        --lift-to chr6 chr17 \
        --out-json coord_verification.json
"""

import argparse
import gzip
import json
import re
import sys

import pandas as pd

PRIMARY_CONTIG_RE = re.compile(r"^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$")
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def iter_fasta(path: str, wanted: set):
    """Yield (contig_name, uppercase_sequence) for contigs in *wanted*, streaming."""
    opener = gzip.open if path.endswith(".gz") else open
    name, chunks, keep = None, [], False
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if keep and name:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                keep = name in wanted
                chunks = []
            elif keep:
                chunks.append(line.strip())
    if keep and name:
        yield name, "".join(chunks).upper()


def classify_match(seq: str, row, spacer: str) -> dict:
    """Compare the reference sequence at the row's coordinates to the spacer.

    Tries the spec convention (0-based half-open, spacer only) first, then the
    variants the file might actually use, and reports which one fits.
    """
    start, end = int(row.guide_start), int(row.guide_end)
    slen = len(spacer)
    rc = reverse_complement(spacer)

    candidates = {
        # (label, slice) — spec: 0-based, half-open, PAM excluded
        "exact":              (start, end),
        # guide_end includes the 3 bp PAM
        "pam_included":       (start, end - 3),
        # 1-based start (off-by-one)
        "one_based":          (start - 1, start - 1 + slen),
        # window is spacer+PAM on the minus strand, so the spacer is at the far end
        "pam_included_left":  (start + 3, end),
    }
    for label, (s, e) in candidates.items():
        if s < 0 or e > len(seq) or e - s != slen:
            continue
        ref = seq[s:e]
        if ref == spacer:
            return {"convention": label, "strand_of_spacer": "+", "ref": ref}
        if ref == rc:
            return {"convention": label, "strand_of_spacer": "-", "ref": ref}

    return {"convention": None, "strand_of_spacer": None,
            "ref": seq[start:end] if end <= len(seq) else ""}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("guides")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--lift-to", nargs="*", default=["chr6", "chr17"],
                    help="primary contigs to hold in memory and search for alt spacers")
    ap.add_argument("--out-json", default=None)
    args = ap.parse_args()

    df = pd.read_csv(args.guides, sep="\t", dtype=str)
    positioned = df.dropna(subset=["guide_chr", "guide_start", "guide_end"]).copy()
    print(f"{len(positioned):,} positioned guides across "
          f"{positioned.guide_chr.nunique()} contigs", file=sys.stderr)

    non_primary = positioned[~positioned.guide_chr.str.match(PRIMARY_CONTIG_RE)]
    alt_spacers = {r.guide_id: r.spacer for r in non_primary.itertuples()}
    print(f"{len(alt_spacers)} guides on non-primary contigs: "
          f"{list(alt_spacers)}", file=sys.stderr)

    wanted = set(positioned.guide_chr) | set(args.lift_to)
    verdicts, lift_hits = {}, {gid: [] for gid in alt_spacers}

    for contig, seq in iter_fasta(args.fasta, wanted):
        rows = positioned[positioned.guide_chr == contig]
        for r in rows.itertuples():
            verdicts[r.guide_id] = classify_match(seq, r, r.spacer)
        print(f"  {contig}: {len(seq):>12,} bp  |  {len(rows):>4} guides checked",
              file=sys.stderr)

        # search this contig for every alt-contig spacer, if it is a lift target
        if contig in args.lift_to:
            for gid, spacer in alt_spacers.items():
                rc = reverse_complement(spacer)
                for strand, probe in (("+", spacer), ("-", rc)):
                    pos = seq.find(probe)
                    while pos != -1:
                        lift_hits[gid].append({
                            "contig": contig, "start": pos, "end": pos + len(probe),
                            "strand": strand,
                        })
                        pos = seq.find(probe, pos + 1)

    # ── Report ────────────────────────────────────────────────────────────────
    conv = pd.Series([v["convention"] for v in verdicts.values()]).value_counts(dropna=False)
    print("\n=== Spacer-vs-reference agreement (all positioned guides) ===")
    print(conv.to_string())

    unmatched = [g for g, v in verdicts.items() if v["convention"] is None]
    if unmatched:
        print(f"\n{len(unmatched)} guides whose spacer matches the reference at NO "
              f"tested convention: {unmatched[:10]}")

    print("\n=== Primary-assembly matches for non-primary-contig guides ===")
    for gid, hits in lift_hits.items():
        print(f"\n{gid}  spacer={alt_spacers[gid]}  "
              f"(alt placement verdict: {verdicts.get(gid, {}).get('convention')})")
        if not hits:
            print("  no exact match on any lift-target contig")
        for h in hits:
            print(f"  {h['contig']}:{h['start']}-{h['end']} ({h['strand']})")

    if args.out_json:
        with open(args.out_json, "w") as fh:
            json.dump({"verdicts": verdicts, "lift_hits": lift_hits,
                       "alt_spacers": alt_spacers}, fh, indent=2)
        print(f"\nWrote {args.out_json}")


if __name__ == "__main__":
    main()

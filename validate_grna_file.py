#!/usr/bin/env python3
"""
Validate an IGVF Per-Guide Metadata TSV.gz file against the spec.
Prints a structured diagnostic report. Does NOT modify any files.

Usage:
    python3 validate_grna_file.py <input_file.tsv.gz>
"""

import gzip
import os
import re
import sys

import pandas as pd
import requests


# ── Constants ─────────────────────────────────────────────────────────────────

GTF_PATH = "input/IGVFFI9573KOZR.gtf.gz"
GTF_URL  = ("https://api.data.igvf.org/reference-files/IGVFFI9573KOZR/"
             "@@download/IGVFFI9573KOZR.gtf.gz")

REQUIRED_COLUMNS = [
    "guide_id", "spacer", "targeting", "type", "guide_chr", "guide_start",
    "guide_end", "strand", "pam", "genomic_element", "intended_target_name",
    "intended_target_chr", "intended_target_start", "intended_target_end",
    "putative_target_genes", "reporter", "imperfect", "description",
]

VALID_TYPES = {"targeting", "safe-targeting", "non-targeting"}
VALID_GE    = {"promoter", "distal element"}
ENSG_RE     = re.compile(r"^ENSG\d+$")
COORD_RE    = re.compile(r"^chr\S+:\d+-\d+$")


# ── GTF / gene map ────────────────────────────────────────────────────────────

def download_gtf(url: str, dest: str) -> None:
    print(f"Downloading GENCODE GTF to {dest} …", file=sys.stderr)
    r = requests.get(url, stream=True, timeout=120)
    r.raise_for_status()
    with open(dest, "wb") as fh:
        for chunk in r.iter_content(chunk_size=65536):
            fh.write(chunk)
    print("Download complete.", file=sys.stderr)


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


# ── Placeholder classification ────────────────────────────────────────────────

def classify_placeholder(desc: str, fallback: str) -> str:
    """Return a human-readable fix strategy for a gene symbol placeholder."""
    # Enhancer check must come before mouse-gene check (both match capital + lowercase)
    if desc.startswith("Enhancer"):
        return "enhancer description (use coordinate string)"
    if re.match(r"^LOC\d+$", fallback):
        return "NCBI LOC entry (no Ensembl ID — manual review)"
    if len(fallback) > 1 and fallback[0].isupper() and fallback[1].islower():
        return "mouse gene name (try uppercasing)"
    return "renamed/obsolete human gene (try Ensembl REST API)"


def extract_gene_candidates(desc: str) -> list:
    """Return ordered list of candidate gene symbols from a description string."""
    cleaned = re.sub(r"_TSS\d*", "", desc)
    if "_and_" in cleaned:
        parts = cleaned.split("_and_")
        candidates = [parts[0], parts[-1]] + parts[1:-1]
    else:
        candidates = [cleaned]
    extras = [c[:-4] for c in candidates if c.endswith("_alt")]
    return candidates + extras


# ── Issue accumulator ─────────────────────────────────────────────────────────

class Issue:
    def __init__(self, field: str, severity: str, message: str,
                 count: int = 0, examples: list = None):
        self.field    = field
        self.severity = severity   # "error" or "warn"
        self.message  = message
        self.count    = count
        self.examples = examples or []

    def symbol(self):
        return "✗" if self.severity == "error" else "!"


# ── Main validation ───────────────────────────────────────────────────────────

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 validate_grna_file.py <input_file.tsv.gz>")
        sys.exit(1)

    input_file = sys.argv[1]

    # Load gene map
    if not os.path.exists(GTF_PATH):
        download_gtf(GTF_URL, GTF_PATH)
    print("Parsing GENCODE v43 GTF …", file=sys.stderr)
    gene_map = load_gene_map(GTF_PATH)
    print(f"  Loaded {len(gene_map):,} gene_name → ENSEMBL mappings", file=sys.stderr)

    # Load input
    print(f"Reading {input_file} …", file=sys.stderr)
    df = pd.read_csv(input_file, sep="\t", dtype=str).fillna("")
    print(f"  {len(df):,} rows loaded", file=sys.stderr)

    issues = []

    # ── Structural: required columns ──────────────────────────────────────────
    missing_cols = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    extra_cols   = [c for c in df.columns if c not in REQUIRED_COLUMNS]

    if missing_cols:
        issues.append(Issue(
            "structure", "error",
            f"Missing required columns: {missing_cols}",
            count=len(missing_cols),
        ))
    if extra_cols:
        issues.append(Issue(
            "structure", "warn",
            f"Extra columns (not in spec): {extra_cols}",
            count=len(extra_cols),
        ))

    # Abort further checks if key columns are missing
    essential = {"type", "targeting", "genomic_element", "intended_target_name",
                 "intended_target_chr", "intended_target_start", "intended_target_end",
                 "guide_chr", "guide_start", "guide_end", "guide_id", "description"}
    if essential & set(missing_cols):
        _print_report(input_file, df, issues)
        return

    # Convenience masks
    mask_tgt  = df["type"] == "targeting"
    mask_safe = df["type"] == "safe-targeting"
    mask_non  = df["type"] == "non-targeting"
    mask_other = ~(mask_tgt | mask_safe | mask_non)

    # ── type field ────────────────────────────────────────────────────────────
    bad_type = mask_other
    if bad_type.any():
        examples = df.loc[bad_type, "type"].unique().tolist()[:5]
        issues.append(Issue(
            "type", "error",
            f"{bad_type.sum():,} rows have invalid type values",
            count=int(bad_type.sum()), examples=examples,
        ))

    # ── targeting field ───────────────────────────────────────────────────────
    bad_case = df["targeting"].isin({"TRUE", "FALSE", "True", "False", "1", "0"})
    if bad_case.any():
        examples = df.loc[bad_case, "targeting"].value_counts().to_dict()
        issues.append(Issue(
            "targeting", "error",
            f"{bad_case.sum():,} rows have non-lowercase targeting value (e.g. TRUE/FALSE)",
            count=int(bad_case.sum()), examples=[str(examples)],
        ))

    bad_targeting_vals = ~df["targeting"].isin({"true", "false", "TRUE", "FALSE",
                                                 "True", "False", "1", "0"})
    if bad_targeting_vals.any():
        examples = df.loc[bad_targeting_vals, "targeting"].unique().tolist()[:5]
        issues.append(Issue(
            "targeting", "error",
            f"{bad_targeting_vals.sum():,} rows have unrecognized targeting value",
            count=int(bad_targeting_vals.sum()), examples=examples,
        ))

    safe_true = mask_safe & (df["targeting"].str.lower() == "true")
    if safe_true.any():
        issues.append(Issue(
            "targeting", "error",
            f"{safe_true.sum():,} safe-targeting rows have targeting='true' (must be 'false')",
            count=int(safe_true.sum()),
        ))

    # ── genomic_element (targeting rows) ──────────────────────────────────────
    ge_empty = mask_tgt & (df["genomic_element"] == "")
    if ge_empty.any():
        issues.append(Issue(
            "genomic_element", "error",
            f"{ge_empty.sum():,} targeting rows have empty genomic_element",
            count=int(ge_empty.sum()),
        ))

    ge_invalid = mask_tgt & (~df["genomic_element"].isin(VALID_GE)) & (df["genomic_element"] != "")
    if ge_invalid.any():
        examples = df.loc[ge_invalid, "genomic_element"].unique().tolist()[:5]
        issues.append(Issue(
            "genomic_element", "error",
            f"{ge_invalid.sum():,} targeting rows have invalid genomic_element",
            count=int(ge_invalid.sum()), examples=examples,
        ))

    # genomic_element must be empty for non-targeting/safe-targeting
    ge_nonempty_nontgt = (~mask_tgt) & (df["genomic_element"] != "")
    if ge_nonempty_nontgt.any():
        issues.append(Issue(
            "genomic_element", "error",
            f"{ge_nonempty_nontgt.sum():,} non-targeting/safe-targeting rows have non-empty genomic_element",
            count=int(ge_nonempty_nontgt.sum()),
        ))

    # ── intended_target_name ──────────────────────────────────────────────────
    itn_empty = mask_tgt & (df["intended_target_name"] == "")
    if itn_empty.any():
        issues.append(Issue(
            "intended_target_name", "error",
            f"{itn_empty.sum():,} targeting rows have empty intended_target_name",
            count=int(itn_empty.sum()),
        ))

    # Format checks for non-empty targeting rows
    mask_tgt_filled = mask_tgt & (df["intended_target_name"] != "")
    if mask_tgt_filled.any():
        sub = df[mask_tgt_filled]

        # Promoter guides must be ENSG*
        promoter_rows = sub[sub["genomic_element"] == "promoter"]
        if not promoter_rows.empty:
            not_ensg = promoter_rows[
                ~promoter_rows["intended_target_name"].str.match(r"^ENSG\d+$")
            ]
            if not not_ensg.empty:
                # Split into: (a) coord strings in promoter rows (ge mismatch), (b) symbol placeholders
                is_coord = not_ensg["intended_target_name"].str.match(r"^chr\S+:\d+-\d+$")
                ge_mismatch = not_ensg[is_coord]
                placeholders = not_ensg[~is_coord]

                if not ge_mismatch.empty:
                    n_desc = len(ge_mismatch["description"].unique())
                    ex = ge_mismatch["description"].unique().tolist()[:3]
                    issues.append(Issue(
                        "intended_target_name", "error",
                        f"{ge_mismatch.shape[0]:,} targeting rows have genomic_element='promoter' but "
                        f"coord-string intended_target_name ({n_desc} unique descriptions — "
                        f"fix_grna_file.py resolved these as enhancers but left genomic_element as 'promoter'):\n"
                        f"       e.g. {', '.join(ex)}",
                        count=int(ge_mismatch.shape[0]),
                    ))

                if not placeholders.empty:
                    placeholder_classes: dict[str, list] = {}
                    for _, row in placeholders.drop_duplicates("description").iterrows():
                        candidates = extract_gene_candidates(row["description"])
                        fallback   = candidates[0] if candidates else row["description"]
                        cat = classify_placeholder(row["description"], fallback)
                        placeholder_classes.setdefault(cat, []).append(fallback)

                    n_total    = len(placeholders["description"].unique())
                    class_lines = []
                    manual_review = []
                    for cat, syms in sorted(placeholder_classes.items()):
                        class_lines.append(f"       - {cat}: {len(syms):4d}  (e.g. {', '.join(syms[:3])})")
                        if "manual" in cat:
                            manual_review.extend(syms)

                    issues.append(Issue(
                        "intended_target_name", "error",
                        f"{placeholders.shape[0]:,} targeting promoter rows have gene-symbol placeholder "
                        f"(not ENSG* or coord) ({n_total} unique descriptions):\n" + "\n".join(class_lines),
                        count=int(placeholders.shape[0]),
                        examples=manual_review[:10],
                    ))

        # Distal element guides must be coord string
        distal_rows = sub[sub["genomic_element"] == "distal element"]
        if not distal_rows.empty:
            not_coord = distal_rows[
                ~distal_rows["intended_target_name"].str.match(r"^chr\S+:\d+-\d+$")
            ]
            if not not_coord.empty:
                examples = not_coord["intended_target_name"].unique().tolist()[:5]
                issues.append(Issue(
                    "intended_target_name", "error",
                    f"{not_coord.shape[0]:,} targeting distal-element rows have malformed intended_target_name (expected chr:start-end)",
                    count=int(not_coord.shape[0]), examples=examples,
                ))

    # intended_target_name must be empty for non-targeting/safe-targeting
    itn_nonempty_nontgt = (~mask_tgt) & (df["intended_target_name"] != "")
    if itn_nonempty_nontgt.any():
        issues.append(Issue(
            "intended_target_name", "error",
            f"{itn_nonempty_nontgt.sum():,} non-targeting/safe-targeting rows have non-empty intended_target_name",
            count=int(itn_nonempty_nontgt.sum()),
        ))

    # ── intended_target_chr/start/end ─────────────────────────────────────────
    for col in ["intended_target_chr", "intended_target_start", "intended_target_end"]:
        coord_empty = mask_tgt & (df[col] == "")
        if coord_empty.any():
            issues.append(Issue(
                col, "error",
                f"{coord_empty.sum():,} targeting rows have empty {col}",
                count=int(coord_empty.sum()),
            ))
        coord_nonempty_nontgt = (~mask_tgt) & (df[col] != "")
        if coord_nonempty_nontgt.any():
            issues.append(Issue(
                col, "error",
                f"{coord_nonempty_nontgt.sum():,} non-targeting/safe-targeting rows have non-empty {col}",
                count=int(coord_nonempty_nontgt.sum()),
            ))

    # Coordinate consistency: chr matches, start ≤ min(guide_start), end ≥ max(guide_end)
    if mask_tgt.any():
        tdf = df[mask_tgt].copy()
        tdf["_gs"] = pd.to_numeric(tdf["guide_start"], errors="coerce")
        tdf["_ge"] = pd.to_numeric(tdf["guide_end"],   errors="coerce")
        tdf["_ts"] = pd.to_numeric(tdf["intended_target_start"], errors="coerce")
        tdf["_te"] = pd.to_numeric(tdf["intended_target_end"],   errors="coerce")

        # Chr mismatch
        chr_mismatch = tdf["intended_target_chr"] != tdf["guide_chr"]
        chr_mismatch &= tdf["intended_target_chr"] != ""
        if chr_mismatch.any():
            issues.append(Issue(
                "intended_target_chr", "error",
                f"{chr_mismatch.sum():,} targeting rows have intended_target_chr ≠ guide_chr",
                count=int(chr_mismatch.sum()),
            ))

        # Start > guide_start (target start should be ≤ any guide start in the group)
        start_too_large = tdf["_ts"] > tdf["_gs"]
        start_too_large = start_too_large.fillna(False)
        if start_too_large.any():
            issues.append(Issue(
                "intended_target_start", "warn",
                f"{start_too_large.sum():,} targeting rows have intended_target_start > guide_start",
                count=int(start_too_large.sum()),
            ))

        # End < guide_end
        end_too_small = tdf["_te"] < tdf["_ge"]
        end_too_small = end_too_small.fillna(False)
        if end_too_small.any():
            issues.append(Issue(
                "intended_target_end", "warn",
                f"{end_too_small.sum():,} targeting rows have intended_target_end < guide_end",
                count=int(end_too_small.sum()),
            ))

    _print_report(input_file, df, issues)


def _print_report(input_file: str, df: pd.DataFrame, issues: list) -> None:
    mask_tgt  = df["type"] == "targeting"   if "type" in df.columns else pd.Series(False, index=df.index)
    mask_safe = df["type"] == "safe-targeting" if "type" in df.columns else pd.Series(False, index=df.index)
    mask_non  = df["type"] == "non-targeting"  if "type" in df.columns else pd.Series(False, index=df.index)

    print()
    print("=== IGVF Per-Guide Metadata Validation Report ===")
    print(f"File: {input_file}")
    print(f"Rows: {len(df):,}  |  "
          f"targeting: {mask_tgt.sum():,}  |  "
          f"safe-targeting: {mask_safe.sum():,}  |  "
          f"non-targeting: {mask_non.sum():,}")
    print()

    # Group issues by field
    fields_seen = []
    field_issues: dict[str, list] = {}
    for iss in issues:
        if iss.field not in field_issues:
            fields_seen.append(iss.field)
            field_issues[iss.field] = []
        field_issues[iss.field].append(iss)

    # Special: structural section
    if "structure" in field_issues:
        print(f"── Structural {'─' * 63}")
        for iss in field_issues["structure"]:
            print(f"  {iss.symbol()}  {iss.message}")
        print()

    all_fields = [f for f in fields_seen if f != "structure"]

    # For each remaining field, print a section
    for field in all_fields:
        header = f"── {field} "
        print(header + "─" * max(1, 78 - len(header)))
        field_ok = True
        for iss in field_issues[field]:
            field_ok = False
            print(f"  {iss.symbol()}  {iss.message}")
            if iss.examples:
                # Don't print example list if it's embedded in message already
                pass
        if field_ok:
            print(f"  ✓  OK")
        print()

    # Print ✓ for fields with no issues
    clean_fields = []
    for col in ["structure", "type", "targeting", "genomic_element",
                "intended_target_name", "intended_target_chr",
                "intended_target_start", "intended_target_end"]:
        if col not in field_issues:
            clean_fields.append(col)

    if clean_fields:
        header = "── Fields with no issues "
        print(header + "─" * max(1, 78 - len(header)))
        for col in clean_fields:
            print(f"  ✓  {col}")
        print()

    # Summary
    print("══ SUMMARY " + "═" * 67)
    n_errors = sum(1 for iss in issues if iss.severity == "error")
    n_warns  = sum(1 for iss in issues if iss.severity == "warn")
    if issues:
        print(f"  {n_errors} error categor{'y' if n_errors == 1 else 'ies'}, "
              f"{n_warns} warning categor{'y' if n_warns == 1 else 'ies'} found")
    else:
        print("  ✓  No issues found — file passes all validation checks")
    if n_errors > 0:
        print("  Suggested fix: run fix_grna_file.py (handles all automated corrections)")
    # Count manual-review placeholders from issues
    for iss in issues:
        if "manual" in iss.message.lower() and iss.examples:
            print(f"  Manual review required for: {', '.join(iss.examples)}")
    print()


if __name__ == "__main__":
    main()

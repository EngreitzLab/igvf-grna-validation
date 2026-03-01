#!/usr/bin/env python3
"""
Validate an IGVF Per-Guide Metadata TSV.gz file against the spec.
Prints a structured diagnostic report. Does NOT modify any files.

Spec: FINAL_FINAL_Format_for_Per_Guide_Metadata_Finalized_11_26_25.pdf

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

COORD_SPAN_WARN_BP = 50   # overshoot ≤ this → warn; > this → error

GTF_PATH = "input/IGVFFI9573KOZR.gtf.gz"
GTF_URL  = ("https://api.data.igvf.org/reference-files/IGVFFI9573KOZR/"
             "@@download/IGVFFI9573KOZR.gtf.gz")

# Required columns (non-optional per spec — no leading "-" in spec)
REQUIRED_COLUMNS = [
    "guide_id", "spacer", "targeting", "type",
    "guide_chr", "guide_start", "guide_end", "strand", "pam",
    "genomic_element", "intended_target_name",
    "intended_target_chr", "intended_target_start", "intended_target_end",
]

# Optional columns (marked with "(-)" in spec)
OPTIONAL_COLUMNS = [
    "putative_target_genes", "reporter", "imperfect", "description",
]

ALL_KNOWN_COLUMNS = REQUIRED_COLUMNS + OPTIONAL_COLUMNS

# Valid type values per spec (full enum)
VALID_TYPES = {
    "targeting", "safe-targeting", "non-targeting",
    "positive control", "negative control", "variant",
}

# Per spec: if type ∈ these, targeting must be False
TYPES_REQUIRE_FALSE = {"non-targeting", "safe-targeting"}

# Valid targeting values per spec (boolean, title case)
VALID_TARGETING = {"True", "False"}

# Valid genomic_element values per spec (full enum)
VALID_GE = {
    "variant", "promoter", "enhancer", "insulator", "silencer",
    "distal element", "splice site", "gene",
}

# genomic_element types whose intended_target_name must be an ENSEMBL ID
GE_NEEDS_ENSG = {"promoter", "splice site", "gene"}

# genomic_element types whose intended_target_name must be a coord string
GE_NEEDS_COORD = {"enhancer", "insulator", "silencer", "distal element"}

# genomic_element types where putative_target_genes is required for positive controls
GE_NEEDS_PTG = {"enhancer", "insulator", "silencer", "distal element"}

ENSG_RE   = re.compile(r"^ENSG\d+$")
COORD_RE  = re.compile(r"^chr\S+:\d+-\d+$")
VALID_STRAND = {"+", "-"}


# ── GTF / gene map ────────────────────────────────────────────────────────────

def download_gtf(url: str, dest: str) -> None:
    print(f"Downloading GENCODE GTF to {dest} …", file=sys.stderr)
    r = requests.get(url, stream=True, timeout=120)
    r.raise_for_status()
    with open(dest, "wb") as fh:
        for chunk in r.iter_content(chunk_size=65536):
            fh.write(chunk)
    print("Download complete.", file=sys.stderr)


# Strings considered equivalent to NaN/missing for optional fields
_NAN_STRINGS = {"", "nan", "NaN", "NAN", "na", "NA", "N/A", "n/a", "none", "None", "null", "NULL"}

def _is_nan_like(series: "pd.Series") -> "pd.Series":
    """Return boolean mask: True where the value is empty or a NaN-equivalent string."""
    return series.isin(_NAN_STRINGS)


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

    # Load gene map (used only for context; not needed for format validation)
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

    # ── S1: Required and optional columns ─────────────────────────────────────
    missing_required = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    missing_optional = [c for c in OPTIONAL_COLUMNS if c not in df.columns]
    extra_cols       = [c for c in df.columns if c not in ALL_KNOWN_COLUMNS]

    if missing_required:
        issues.append(Issue(
            "structure", "error",
            f"Missing required columns: {missing_required}",
            count=len(missing_required),
        ))
    if missing_optional:
        issues.append(Issue(
            "structure", "warn",
            f"Missing optional columns (will be treated as empty): {missing_optional}",
            count=len(missing_optional),
        ))
    if extra_cols:
        issues.append(Issue(
            "structure", "warn",
            f"Extra columns (not in spec): {extra_cols}",
            count=len(extra_cols),
        ))

    # Abort further checks if any essential required column is missing
    essential = {
        "type", "targeting", "genomic_element", "intended_target_name",
        "intended_target_chr", "intended_target_start", "intended_target_end",
        "guide_chr", "guide_start", "guide_end", "guide_id", "spacer",
    }
    if essential & set(missing_required):
        _print_report(input_file, df, issues)
        return

    # Backfill missing optional columns as empty strings
    for col in OPTIONAL_COLUMNS:
        if col not in df.columns:
            df[col] = ""

    # ── S2: guide_id uniqueness and 1-to-1 mapping to spacer ─────────────────
    dup_guids = df.duplicated("guide_id").sum()
    if dup_guids:
        examples = df.loc[df.duplicated("guide_id", keep=False), "guide_id"].unique().tolist()[:5]
        issues.append(Issue(
            "guide_id", "error",
            f"{dup_guids:,} duplicate guide_id values (spec requires unique guide_ids)",
            count=int(dup_guids), examples=examples,
        ))

    multi_spacer = (df.groupby("guide_id")["spacer"].nunique() > 1).sum()
    if multi_spacer:
        issues.append(Issue(
            "spacer", "error",
            f"{multi_spacer:,} guide_ids map to more than one spacer (spec requires 1-to-1 mapping)",
            count=int(multi_spacer),
        ))

    # ── T1: type field — must be from the spec enum ───────────────────────────
    bad_type = ~df["type"].isin(VALID_TYPES)
    if bad_type.any():
        val_counts = df.loc[bad_type, "type"].value_counts().to_dict()
        issues.append(Issue(
            "type", "error",
            f"{bad_type.sum():,} rows have invalid type value "
            f"(valid: {sorted(VALID_TYPES)}): {val_counts}",
            count=int(bad_type.sum()),
        ))

    # Type masks used throughout
    mask_type = {t: (df["type"] == t) for t in VALID_TYPES}
    mask_requires_false = mask_type["non-targeting"] | mask_type["safe-targeting"]

    # ── T2: targeting field — must be "True" or "False" (title case) ──────────
    bad_targeting = ~df["targeting"].isin(VALID_TARGETING)
    if bad_targeting.any():
        val_counts = df.loc[bad_targeting, "targeting"].value_counts().to_dict()
        issues.append(Issue(
            "targeting", "error",
            f"{bad_targeting.sum():,} rows have invalid targeting value "
            f"(must be 'True' or 'False', title case): {val_counts}",
            count=int(bad_targeting.sum()),
        ))

    # ── T3: targeting consistency with type ───────────────────────────────────
    # Spec: "if type ∈ {non-targeting, safe-targeting}, targeting should be False"
    nontgt_true = mask_requires_false & (df["targeting"].str.upper() == "TRUE")
    if nontgt_true.any():
        by_type = df.loc[nontgt_true, "type"].value_counts().to_dict()
        issues.append(Issue(
            "targeting", "error",
            f"{nontgt_true.sum():,} non-targeting/safe-targeting rows have targeting='True' "
            f"(spec requires 'False'): {by_type}",
            count=int(nontgt_true.sum()),
        ))

    # Warn: type="targeting" rows should have targeting="True"
    tgt_type_false = mask_type["targeting"] & (df["targeting"].str.upper() == "FALSE")
    if tgt_type_false.any():
        issues.append(Issue(
            "targeting", "warn",
            f"{tgt_type_false.sum():,} rows with type='targeting' have targeting='False' "
            f"(expected 'True')",
            count=int(tgt_type_false.sum()),
        ))

    # Primary masks for field-presence checks.
    # Use case-insensitive matching so element checks work even when the
    # targeting case is wrong (case error is already reported above).
    is_true  = df["targeting"].str.upper() == "TRUE"   # element fields required
    is_false = df["targeting"].str.upper() == "FALSE"  # element fields not required

    # For the "must be empty/NaN" hygiene check: only enforce on non-targeting
    # rows, NOT safe-targeting.  Safe-targeting guides target real genomic loci
    # and may legitimately have element fields populated.
    must_be_empty = is_false & (df["type"] == "non-targeting")

    # ── T4: strand — must be "+" or "-" for targeting=True rows ──────────────
    if "strand" in df.columns:
        strand_empty = is_true & _is_nan_like(df["strand"])
        if strand_empty.any():
            issues.append(Issue(
                "strand", "error",
                f"{strand_empty.sum():,} targeting=True rows have empty strand",
                count=int(strand_empty.sum()),
            ))

        bad_strand = is_true & ~_is_nan_like(df["strand"]) & (~df["strand"].isin(VALID_STRAND))
        if bad_strand.any():
            examples = df.loc[bad_strand, "strand"].unique().tolist()[:5]
            issues.append(Issue(
                "strand", "error",
                f"{bad_strand.sum():,} targeting=True rows have invalid strand value "
                f"(must be '+' or '-'): {examples}",
                count=int(bad_strand.sum()),
            ))

    # ── T5: guide_chr / guide_start / guide_end required for targeting=True ───
    for col in ["guide_chr", "guide_start", "guide_end"]:
        if col in df.columns:
            empty = is_true & _is_nan_like(df[col])
            if empty.any():
                issues.append(Issue(
                    col, "error",
                    f"{empty.sum():,} targeting=True rows have empty/NaN {col}",
                    count=int(empty.sum()),
                ))

    # ── E1: genomic_element must be non-empty for targeting=True rows ─────────
    ge_empty = is_true & _is_nan_like(df["genomic_element"])
    if ge_empty.any():
        issues.append(Issue(
            "genomic_element", "error",
            f"{ge_empty.sum():,} targeting=True rows have empty/NaN genomic_element",
            count=int(ge_empty.sum()),
        ))

    # ── E2: genomic_element must be a valid enum value ────────────────────────
    ge_invalid = is_true & (~df["genomic_element"].isin(VALID_GE)) & ~_is_nan_like(df["genomic_element"])
    if ge_invalid.any():
        examples = df.loc[ge_invalid, "genomic_element"].unique().tolist()[:8]
        issues.append(Issue(
            "genomic_element", "error",
            f"{ge_invalid.sum():,} targeting=True rows have invalid genomic_element "
            f"(valid: {sorted(VALID_GE)}): {examples}",
            count=int(ge_invalid.sum()),
        ))

    # ── E2b: warn when "enhancer" is used — prefer "distal element" ───────────
    ge_enhancer = is_true & (df["genomic_element"] == "enhancer")
    if ge_enhancer.any():
        issues.append(Issue(
            "genomic_element", "warn",
            f"{ge_enhancer.sum():,} targeting=True rows use genomic_element='enhancer'; "
            f"prefer 'distal element' unless a more specific classification is known",
            count=int(ge_enhancer.sum()),
        ))

    # ── E3: genomic_element must be empty/NaN for non-targeting rows ─────────
    ge_nonempty_nontgt = must_be_empty & ~_is_nan_like(df["genomic_element"])
    if ge_nonempty_nontgt.any():
        examples = df.loc[ge_nonempty_nontgt, "genomic_element"].unique().tolist()[:5]
        issues.append(Issue(
            "genomic_element", "error",
            f"{ge_nonempty_nontgt.sum():,} non-targeting rows have non-empty genomic_element "
            f"(must be empty or NaN for type='non-targeting'): e.g. {examples}",
            count=int(ge_nonempty_nontgt.sum()),
        ))

    # ── N1: intended_target_name — required for targeting=True ───────────────
    itn_empty = is_true & _is_nan_like(df["intended_target_name"])
    if itn_empty.any():
        issues.append(Issue(
            "intended_target_name", "error",
            f"{itn_empty.sum():,} targeting=True rows have empty/NaN intended_target_name",
            count=int(itn_empty.sum()),
        ))

    # ── N2: intended_target_name format per genomic_element ──────────────────
    # Only check rows with valid ge and non-empty intended_target_name
    mask_check_itn = is_true & ~_is_nan_like(df["intended_target_name"]) & df["genomic_element"].isin(VALID_GE)
    if mask_check_itn.any():
        sub = df[mask_check_itn]

        for ge_type in GE_NEEDS_ENSG:
            ge_rows = sub[sub["genomic_element"] == ge_type]
            if not ge_rows.empty:
                not_ensg = ge_rows[~ge_rows["intended_target_name"].str.match(r"^ENSG\d+$")]
                if not not_ensg.empty:
                    examples = not_ensg["intended_target_name"].unique().tolist()[:5]
                    issues.append(Issue(
                        "intended_target_name", "error",
                        f"{not_ensg.shape[0]:,} targeting rows with genomic_element='{ge_type}' "
                        f"have non-ENSG intended_target_name (expected ENSG ID): {examples}",
                        count=int(not_ensg.shape[0]),
                    ))

        for ge_type in GE_NEEDS_COORD:
            ge_rows = sub[sub["genomic_element"] == ge_type]
            if not ge_rows.empty:
                not_coord = ge_rows[~ge_rows["intended_target_name"].str.match(r"^chr\S+:\d+-\d+$")]
                if not not_coord.empty:
                    examples = not_coord["intended_target_name"].unique().tolist()[:5]
                    issues.append(Issue(
                        "intended_target_name", "error",
                        f"{not_coord.shape[0]:,} targeting rows with genomic_element='{ge_type}' "
                        f"have malformed intended_target_name "
                        f"(expected chr:start-end): {examples}",
                        count=int(not_coord.shape[0]),
                    ))

    # ── N3: intended_target_name must be empty/NaN for non-targeting rows ────
    itn_nonempty_nontgt = must_be_empty & ~_is_nan_like(df["intended_target_name"])
    if itn_nonempty_nontgt.any():
        issues.append(Issue(
            "intended_target_name", "error",
            f"{itn_nonempty_nontgt.sum():,} non-targeting rows have non-empty "
            f"intended_target_name (must be empty or NaN for type='non-targeting')",
            count=int(itn_nonempty_nontgt.sum()),
        ))

    # ── C1: intended_target_chr/start/end required for targeting=True ─────────
    for col in ["intended_target_chr", "intended_target_start", "intended_target_end"]:
        coord_empty = is_true & _is_nan_like(df[col])
        if coord_empty.any():
            issues.append(Issue(
                col, "error",
                f"{coord_empty.sum():,} targeting=True rows have empty/NaN {col}",
                count=int(coord_empty.sum()),
            ))

        coord_nonempty_nontgt = must_be_empty & ~_is_nan_like(df[col])
        if coord_nonempty_nontgt.any():
            issues.append(Issue(
                col, "error",
                f"{coord_nonempty_nontgt.sum():,} non-targeting rows have non-empty {col} "
                f"(must be empty or NaN for type='non-targeting')",
                count=int(coord_nonempty_nontgt.sum()),
            ))

    # ── C2: coordinate consistency for targeting=True rows ────────────────────
    if is_true.any():
        tdf = df[is_true].copy()
        tdf["_gs"] = pd.to_numeric(tdf["guide_start"], errors="coerce")
        tdf["_ge"] = pd.to_numeric(tdf["guide_end"],   errors="coerce")
        tdf["_ts"] = pd.to_numeric(tdf["intended_target_start"], errors="coerce")
        tdf["_te"] = pd.to_numeric(tdf["intended_target_end"],   errors="coerce")

        # intended_target_chr must match guide_chr
        chr_mismatch = (
            ~_is_nan_like(tdf["intended_target_chr"]) &
            ~_is_nan_like(tdf["guide_chr"]) &
            (tdf["intended_target_chr"] != tdf["guide_chr"])
        )
        if chr_mismatch.any():
            issues.append(Issue(
                "intended_target_chr", "error",
                f"{chr_mismatch.sum():,} targeting rows have intended_target_chr ≠ guide_chr",
                count=int(chr_mismatch.sum()),
            ))

        # intended_target_start should be ≤ guide_start (target spans its guides)
        start_over = (tdf["_ts"] - tdf["_gs"]).clip(lower=0).fillna(0)  # >0 means target_start > guide_start
        start_small = (start_over > 0) & (start_over <= COORD_SPAN_WARN_BP)
        start_large = start_over > COORD_SPAN_WARN_BP
        if start_small.any():
            issues.append(Issue(
                "intended_target_start", "warn",
                f"{start_small.sum():,} targeting rows have intended_target_start > guide_start "
                f"by ≤ {COORD_SPAN_WARN_BP} bp (minor boundary effect)",
                count=int(start_small.sum()),
            ))
        if start_large.any():
            worst = (
                tdf.loc[start_large, ["intended_target_name", "_ts", "_gs"]]
                .assign(over=start_over[start_large])
                .nlargest(3, "over")[["intended_target_name", "over"]]
            )
            worst_str = ", ".join(
                f"{r['intended_target_name']} overshoot={int(r['over'])} bp"
                for _, r in worst.iterrows()
            )
            issues.append(Issue(
                "intended_target_start", "error",
                f"{start_large.sum():,} targeting rows have intended_target_start > guide_start "
                f"by > {COORD_SPAN_WARN_BP} bp — target window does not span its guides; "
                f"worst: [{worst_str}]",
                count=int(start_large.sum()),
            ))

        # intended_target_end should be ≥ guide_end
        end_under = (tdf["_ge"] - tdf["_te"]).clip(lower=0).fillna(0)   # >0 means target_end < guide_end
        end_small = (end_under > 0) & (end_under <= COORD_SPAN_WARN_BP)
        end_large = end_under > COORD_SPAN_WARN_BP
        if end_small.any():
            issues.append(Issue(
                "intended_target_end", "warn",
                f"{end_small.sum():,} targeting rows have intended_target_end < guide_end "
                f"by ≤ {COORD_SPAN_WARN_BP} bp (minor boundary effect)",
                count=int(end_small.sum()),
            ))
        if end_large.any():
            worst = (
                tdf.loc[end_large, ["intended_target_name", "_te", "_ge"]]
                .assign(under=end_under[end_large])
                .nlargest(3, "under")[["intended_target_name", "under"]]
            )
            worst_str = ", ".join(
                f"{r['intended_target_name']} overshoot={int(r['under'])} bp"
                for _, r in worst.iterrows()
            )
            issues.append(Issue(
                "intended_target_end", "error",
                f"{end_large.sum():,} targeting rows have intended_target_end < guide_end "
                f"by > {COORD_SPAN_WARN_BP} bp — target window does not span its guides; "
                f"worst: [{worst_str}]",
                count=int(end_large.sum()),
            ))

    # ── P1: putative_target_genes — required for positive controls with
    #         distal-type genomic_element ────────────────────────────────────
    pc_distal = mask_type["positive control"] & df["genomic_element"].isin(GE_NEEDS_PTG)
    ptg_empty_pc = pc_distal & (df["putative_target_genes"] == "")
    if ptg_empty_pc.any():
        issues.append(Issue(
            "putative_target_genes", "error",
            f"{ptg_empty_pc.sum():,} positive control rows with distal-type genomic_element "
            f"have empty putative_target_genes (required per spec)",
            count=int(ptg_empty_pc.sum()),
        ))

    # Warn: non-empty putative_target_genes values don't look like ENSG IDs
    def _ptg_looks_valid(val: str) -> bool:
        val = val.strip().lstrip("[").rstrip("]")
        parts = [p.strip().strip('"\'') for p in re.split(r'[,;|]\s*', val) if p.strip().strip('"\'')]
        return all(ENSG_RE.match(p) for p in parts) if parts else True

    ptg_filled = df[df["putative_target_genes"] != ""]
    if not ptg_filled.empty:
        bad_ptg = ptg_filled[~ptg_filled["putative_target_genes"].apply(_ptg_looks_valid)]
        if not bad_ptg.empty:
            examples = bad_ptg["putative_target_genes"].unique().tolist()[:5]
            issues.append(Issue(
                "putative_target_genes", "warn",
                f"{len(bad_ptg):,} rows have putative_target_genes that don't appear to be "
                f"ENSG IDs (expected ENSG* or empty): {examples}",
                count=len(bad_ptg),
            ))

    _print_report(input_file, df, issues)


def _print_report(input_file: str, df: pd.DataFrame, issues: list) -> None:
    def _count(col, val):
        return (df[col] == val).sum() if col in df.columns else 0

    print()
    print("=== IGVF Per-Guide Metadata Validation Report ===")
    print(f"File: {input_file}")

    type_parts = [f"Rows: {len(df):,}"]
    for t in ["targeting", "safe-targeting", "non-targeting",
              "positive control", "negative control", "variant"]:
        n = _count("type", t)
        if n:
            type_parts.append(f"{t}: {n:,}")
    print("  |  ".join(type_parts))
    print()

    # Group issues by field, preserving insertion order
    fields_seen: list = []
    field_issues: dict = {}
    for iss in issues:
        if iss.field not in field_issues:
            fields_seen.append(iss.field)
            field_issues[iss.field] = []
        field_issues[iss.field].append(iss)

    if "structure" in field_issues:
        print(f"── Structural {'─' * 63}")
        for iss in field_issues["structure"]:
            print(f"  {iss.symbol()}  {iss.message}")
        print()

    for field in [f for f in fields_seen if f != "structure"]:
        header = f"── {field} "
        print(header + "─" * max(1, 78 - len(header)))
        for iss in field_issues[field]:
            print(f"  {iss.symbol()}  {iss.message}")
        print()

    # Summarise clean fields
    all_checked = [
        "structure", "guide_id", "spacer", "targeting", "type",
        "strand", "guide_chr", "guide_start", "guide_end",
        "genomic_element", "intended_target_name",
        "intended_target_chr", "intended_target_start", "intended_target_end",
        "putative_target_genes",
    ]
    clean = [c for c in all_checked if c not in field_issues]
    if clean:
        header = "── Fields with no issues "
        print(header + "─" * max(1, 78 - len(header)))
        for col in clean:
            print(f"  ✓  {col}")
        print()

    print("══ SUMMARY " + "═" * 67)
    n_err  = sum(1 for iss in issues if iss.severity == "error")
    n_warn = sum(1 for iss in issues if iss.severity == "warn")
    if not issues:
        print("  ✓  No issues found — file passes all validation checks")
    else:
        print(f"  {n_err} error categor{'y' if n_err == 1 else 'ies'}, "
              f"{n_warn} warning categor{'y' if n_warn == 1 else 'ies'} found")
        if n_err > 0:
            print("  Suggested fix: run fix_grna_file.py (handles all automated corrections)")
    print()


if __name__ == "__main__":
    main()

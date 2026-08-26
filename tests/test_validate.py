"""
Unit tests for validate_grna_file.validate_df().

Each test targets a specific spec check (labeled S1, T1, E1, N1, C1, P1, D1, etc.).
Tests use the row-factory helpers from conftest.py; no GTF file is needed.
"""

import pytest
import pandas as pd

from conftest import (
    make_row, make_nt_row, make_safe_row, make_pc_row,
    df_from, errors, warnings, has_error, has_warning,
)
from validate_grna_file import validate_df


# ── Helpers ───────────────────────────────────────────────────────────────────

def clean(issues):
    """Return True if there are no error-severity issues."""
    return not errors(issues)


# ══════════════════════════════════════════════════════════════════════════════
# S1 — Column structure
# ══════════════════════════════════════════════════════════════════════════════

def test_s1_missing_required_column_is_error():
    df = df_from(make_row())
    df = df.drop(columns=["genomic_element"])
    issues = validate_df(df)
    assert has_error(issues, "structure", "genomic_element")


def test_s1_missing_optional_column_is_warning():
    df = df_from(make_row())
    df = df.drop(columns=["description"])
    issues = validate_df(df)
    assert has_warning(issues, "structure", "description")
    assert not errors(issues, "structure")


def test_s1_extra_column_is_warning():
    df = df_from(make_row())
    df["custom_score"] = "1.0"
    issues = validate_df(df)
    assert has_warning(issues, "structure", "custom_score")


def test_s1_missing_essential_column_aborts_early():
    # When an essential column (e.g. targeting) is missing, downstream checks
    # should not run — only S1 issues should be returned.
    df = df_from(make_row())
    df = df.drop(columns=["targeting"])
    issues = validate_df(df)
    fields = {i.field for i in issues}
    assert "structure" in fields
    # No downstream field errors (targeting, genomic_element, etc.)
    assert "targeting" not in fields
    assert "genomic_element" not in fields


# ══════════════════════════════════════════════════════════════════════════════
# S2 — guide_id uniqueness and spacer mapping
# ══════════════════════════════════════════════════════════════════════════════

def test_s2_duplicate_guide_id_is_error():
    df = df_from(make_row(guide_id="g1"), make_row(guide_id="g1", spacer="AAAAAAAAAA"))
    issues = validate_df(df)
    assert has_error(issues, "guide_id", "duplicate")


def test_s2_unique_guide_ids_clean():
    df = df_from(make_row(guide_id="g1"), make_row(guide_id="g2"))
    issues = validate_df(df)
    assert not errors(issues, "guide_id")


def test_s2_same_guide_id_different_spacers_is_error():
    df = df_from(
        make_row(guide_id="g1", spacer="ACGTACGTACGT"),
        make_row(guide_id="g1", spacer="TTTTTTTTTTTT"),
    )
    issues = validate_df(df)
    assert has_error(issues, "spacer", "more than one spacer")


def test_s2_same_spacer_different_guide_ids_clean():
    # Two rows with different guide_ids but same spacer is fine
    df = df_from(
        make_row(guide_id="g1", spacer="ACGTACGTACGT"),
        make_row(guide_id="g2", spacer="ACGTACGTACGT"),
    )
    issues = validate_df(df)
    assert not errors(issues, "spacer")


# ══════════════════════════════════════════════════════════════════════════════
# T1 — type enum
# ══════════════════════════════════════════════════════════════════════════════

def test_t1_invalid_type_is_error():
    df = df_from(make_row(type="tss-targeting"))
    issues = validate_df(df)
    assert has_error(issues, "type")


def test_t1_all_valid_types_clean():
    rows = [
        make_row(guide_id="g1", type="targeting", targeting="True",
                 genomic_element="promoter", intended_target_name="ENSG00000000001"),
        make_nt_row(guide_id="g2", type="non-targeting"),
        make_safe_row(guide_id="g3", type="safe-targeting"),
        make_row(guide_id="g4", type="positive control", targeting="True",
                 genomic_element="promoter", intended_target_name="ENSG00000000004"),
        make_row(guide_id="g5", type="negative control", targeting="False",
                 genomic_element="", intended_target_name="",
                 intended_target_chr="", intended_target_start="", intended_target_end=""),
        make_row(guide_id="g6", type="variant", targeting="True",
                 genomic_element="variant", intended_target_name="ENSG00000000006"),
    ]
    df = df_from(*rows)
    issues = validate_df(df)
    assert not errors(issues, "type")


# ══════════════════════════════════════════════════════════════════════════════
# T2 — targeting field case
# ══════════════════════════════════════════════════════════════════════════════

def test_t2_lowercase_true_is_error():
    df = df_from(make_row(targeting="true"))
    issues = validate_df(df)
    assert has_error(issues, "targeting", "title case")


def test_t2_uppercase_true_is_error():
    df = df_from(make_row(targeting="TRUE"))
    issues = validate_df(df)
    assert has_error(issues, "targeting", "title case")


def test_t2_lowercase_false_is_error():
    df = df_from(make_nt_row(targeting="false"))
    issues = validate_df(df)
    assert has_error(issues, "targeting", "title case")


def test_t2_title_case_clean():
    df = df_from(make_row(targeting="True"), make_nt_row(targeting="False"))
    issues = validate_df(df)
    assert not errors(issues, "targeting")


# ══════════════════════════════════════════════════════════════════════════════
# T3 — targeting / type consistency
# ══════════════════════════════════════════════════════════════════════════════

def test_t3_non_targeting_with_true_is_error():
    df = df_from(make_nt_row(targeting="True"))
    issues = validate_df(df)
    assert has_error(issues, "targeting", "non-targeting")


def test_t3_safe_targeting_with_true_is_error():
    df = df_from(make_safe_row(targeting="True"))
    issues = validate_df(df)
    assert has_error(issues, "targeting", "safe-targeting")


def test_t3_type_targeting_with_false_is_warning():
    df = df_from(make_row(type="targeting", targeting="False"))
    issues = validate_df(df)
    assert has_warning(issues, "targeting", "type='targeting'")


def test_t3_consistency_clean():
    df = df_from(make_row(targeting="True"), make_nt_row(targeting="False"))
    issues = validate_df(df)
    # No T3 errors (T2 may pass too)
    t3_errs = [i for i in errors(issues, "targeting")
               if "non-targeting" in i.message or "safe-targeting" in i.message]
    assert not t3_errs


# ══════════════════════════════════════════════════════════════════════════════
# T4 — strand
# ══════════════════════════════════════════════════════════════════════════════

def test_t4_empty_strand_for_targeting_is_error():
    df = df_from(make_row(strand=""))
    issues = validate_df(df)
    assert has_error(issues, "strand", "empty")


def test_t4_invalid_strand_value_is_error():
    df = df_from(make_row(strand="."))
    issues = validate_df(df)
    assert has_error(issues, "strand", "invalid")


def test_t4_valid_strands_clean():
    df = df_from(make_row(guide_id="g1", strand="+"), make_row(guide_id="g2", strand="-"))
    issues = validate_df(df)
    assert not errors(issues, "strand")


def test_t4_empty_strand_for_non_targeting_allowed():
    # non-targeting rows don't need a strand
    df = df_from(make_nt_row(strand=""))
    issues = validate_df(df)
    assert not errors(issues, "strand")


# ══════════════════════════════════════════════════════════════════════════════
# T5 / T5b — guide coordinates
# ══════════════════════════════════════════════════════════════════════════════

def test_t5_empty_guide_chr_for_targeting_is_error():
    df = df_from(make_row(guide_chr=""))
    issues = validate_df(df)
    assert has_error(issues, "guide_chr", "empty")


def test_t5_empty_guide_start_for_targeting_is_error():
    df = df_from(make_row(guide_start=""))
    issues = validate_df(df)
    assert has_error(issues, "guide_start", "empty")


def test_t5b_non_numeric_guide_start_is_error():
    df = df_from(make_row(guide_start="not_a_number"))
    issues = validate_df(df)
    assert has_error(issues, "guide_start", "non-numeric")


def test_t5b_non_numeric_guide_end_is_error():
    df = df_from(make_row(guide_end="chr1:100"))
    issues = validate_df(df)
    assert has_error(issues, "guide_end", "non-numeric")


def test_t5b_whitespace_guide_start_treated_as_empty():
    # Bug 5 regression: whitespace-only values must be flagged as empty
    df = df_from(make_row(guide_start="   "))
    issues = validate_df(df)
    assert has_error(issues, "guide_start", "empty")


# ══════════════════════════════════════════════════════════════════════════════
# T5c — guide window must span exactly the spacer (PAM excluded)
#
# Spec, guide_end: "should not include the PAM sequence".
#
# Dimensions:
#   window_width:  == len(spacer)  |  == len(spacer)+len(pam)  |  neither
#   strand:        "+" (PAM to the right)  |  "-" (PAM to the left)
#   pam_length:    3 (NGG)  |  6 (NNGRRT)   — must not be hardcoded
#   row_kind:      targeting (has coords)  |  safe-targeting (coords, targeting=False)
#                  |  non-targeting (no coords, must be skipped)
# ══════════════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("strand", ["+", "-"])
def test_t5c_window_equal_to_spacer_is_clean(strand):
    df = df_from(make_row(spacer="A" * 20, guide_start="1000", guide_end="1020",
                          strand=strand))
    issues = validate_df(df)
    for field in ("guide_start", "guide_end"):
        assert not has_error(issues, field, "PAM"), (
            f"Expected no PAM error for a spacer-width window (20 bp spacer, "
            f"[1000,1020), strand {strand}); got a {field} error. Check the "
            f"span == spacer_len comparison in T5c.")


@pytest.mark.parametrize("strand,flagged,clean", [
    ("+", "guide_end",   "guide_start"),   # PAM is 3' → right of the spacer in genome coords
    ("-", "guide_start", "guide_end"),     # PAM is 3' → left of the spacer in genome coords
])
def test_t5c_pam_in_window_is_flagged_on_the_strand_appropriate_end(strand, flagged, clean):
    df = df_from(make_row(spacer="A" * 20, guide_start="1000", guide_end="1023",
                          strand=strand, pam="NGG"))
    issues = validate_df(df)
    assert has_error(issues, flagged, "include the PAM"), (
        f"A 23 bp window around a 20 bp spacer on the '{strand}' strand puts the 3 bp PAM "
        f"in {flagged}; expected an error on {flagged} but found none. "
        f"Got: {[(i.field, i.message) for i in errors(issues)]}")
    assert not has_error(issues, clean, "include the PAM"), (
        f"The PAM is not in {clean} for a '{strand}'-strand guide, but T5c flagged it — "
        f"the strand→end mapping is inverted.")


def test_t5c_uses_actual_pam_length_not_hardcoded_three():
    # 6 bp PAM (SaCas9 NNGRRT): window is spacer+6, so the length must come from the column
    df = df_from(make_row(spacer="A" * 20, guide_start="1000", guide_end="1026",
                          strand="+", pam="NNGRRT"))
    issues = validate_df(df)
    assert has_error(issues, "guide_end", "include the PAM"), (
        "A 26 bp window around a 20 bp spacer with a 6 bp PAM should be recognised as "
        "PAM-inclusive; T5c appears to hardcode a 3 bp PAM. "
        f"Got: {[(i.field, i.message) for i in errors(issues)]}")


def test_t5c_span_mismatch_unrelated_to_pam_is_error():
    # off by 7 — neither len(spacer) nor len(spacer)+len(pam), so it is not a PAM problem
    df = df_from(make_row(spacer="A" * 20, guide_start="1000", guide_end="1027",
                          strand="+", pam="NGG"))
    issues = validate_df(df)
    assert has_error(issues, "guide_end", "not the PAM length"), (
        "A 27 bp window around a 20 bp spacer with a 3 bp PAM matches neither the spacer "
        "width nor spacer+PAM, and must still be reported as a coordinate error rather "
        f"than passing silently. Got: {[(i.field, i.message) for i in errors(issues)]}")
    assert not has_error(issues, "guide_end", "include the PAM"), (
        "An unexplained 7 bp overhang must not be reported as a PAM-inclusion error.")


def test_t5c_applies_to_safe_targeting_rows():
    # safe-targeting rows carry coordinates with targeting=False, so the check must not be
    # gated on targeting=True — this is where the 300-gene library's 60 such rows live
    df = df_from(make_safe_row(spacer="C" * 19, guide_start="5000", guide_end="5022",
                               strand="+", pam="NGG"))
    issues = validate_df(df)
    assert has_error(issues, "guide_end", "include the PAM"), (
        "A safe-targeting row (targeting=False) with a 22 bp window around a 19 bp spacer "
        "must still be flagged; T5c looks to be gated on targeting=True. "
        f"Got: {[(i.field, i.message) for i in errors(issues)]}")


def test_t5c_skips_rows_without_coordinates():
    df = df_from(make_nt_row())
    issues = validate_df(df)
    for field in ("guide_start", "guide_end"):
        assert not has_error(issues, field, "PAM"), (
            f"A non-targeting row has no coordinates, so T5c must skip it; got a {field} "
            f"error instead. Check the coords_ok guard.")


def test_t5c_ignores_whitespace_around_spacer():
    df = df_from(make_row(spacer="  " + "A" * 20 + " ", guide_start="1000",
                          guide_end="1020"))
    issues = validate_df(df)
    assert not has_error(issues, "guide_end", "PAM"), (
        "A padded 20 bp spacer measured 23 characters, so the window looked PAM-inclusive; "
        "T5c must strip whitespace before taking len(spacer).")


# ══════════════════════════════════════════════════════════════════════════════
# T6 — non-primary-assembly contigs
#
# Dimensions:
#   contig_class:  primary (chr1..22,X,Y,M)  |  alt  |  unplaced  |  unlocalized
#                  |  patch fix  |  decoy
#   which_field:   guide_chr  |  intended_target_chr  |  both
#   row_kind:      targeting  |  safe-targeting (targeting=False, has coords)
#                  |  non-targeting (no coords)
# ══════════════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("contig", [
    "chr6_GL000252v2_alt",     # alt haplotype
    "chr17_KI270858v1_alt",
    "chr1_KI270762v1_alt",
    "chr15_KI270905v1_alt",
    "chrUn_GL000195v1",        # unplaced
    "chr1_KI270706v1_random",  # unlocalized
    "chr2_KN538363v1_fix",     # patch scaffold
    "chrEBV",                  # decoy
])
def test_t6_non_primary_guide_chr_is_warning(contig):
    df = df_from(make_row(guide_chr=contig, intended_target_chr=contig))
    issues = validate_df(df)
    assert has_warning(issues, "guide_chr", "non-primary-assembly"), (
        f"{contig!r} is not a primary-assembly contig and must raise a guide_chr warning. "
        f"Got warnings: {[(i.field, i.message) for i in warnings(issues)]}")


@pytest.mark.parametrize("contig", [
    "chr1", "chr9", "chr10", "chr22", "chrX", "chrY", "chrM",
])
def test_t6_primary_contigs_are_clean(contig):
    df = df_from(make_row(guide_chr=contig, intended_target_chr=contig))
    issues = validate_df(df)
    for field in ("guide_chr", "intended_target_chr"):
        assert not has_warning(issues, field, "non-primary-assembly"), (
            f"{contig!r} IS a primary-assembly contig but T6 flagged {field} — "
            f"PRIMARY_CONTIG_RE is too strict.")


def test_t6_flags_intended_target_chr_independently():
    # guide on primary, intended window on an alt contig: only the window is wrong
    df = df_from(make_row(guide_chr="chr6", intended_target_chr="chr6_GL000252v2_alt"))
    issues = validate_df(df)
    assert has_warning(issues, "intended_target_chr", "non-primary-assembly"), (
        "intended_target_chr on an alt contig must be flagged even when guide_chr is "
        f"primary. Got warnings: {[(i.field, i.message) for i in warnings(issues)]}")
    assert not has_warning(issues, "guide_chr", "non-primary-assembly"), (
        "guide_chr is chr6 (primary) and must not be flagged; T6 is not evaluating the "
        "two columns independently.")


def test_t6_applies_to_non_targeting_rows_too():
    # safe-targeting rows carry guide positions but have targeting=False, so the check
    # must not be gated on targeting=True (SAFE_TARGETING_1656 regression)
    df = df_from(make_safe_row(guide_chr="chr17_KI270858v1_alt", intended_target_chr=""))
    issues = validate_df(df)
    assert has_warning(issues, "guide_chr", "non-primary-assembly"), (
        "A safe-targeting row (targeting=False) on an alt contig must be flagged; T6 looks "
        f"to be gated on targeting=True. Got: {[(i.field, i.message) for i in warnings(issues)]}")


def test_t6_empty_coords_are_not_flagged():
    df = df_from(make_nt_row())
    issues = validate_df(df)
    for field in ("guide_chr", "intended_target_chr"):
        assert not has_warning(issues, field, "non-primary-assembly"), (
            f"A non-targeting row has an empty {field}, which must not be treated as a "
            f"non-primary contig. Check the _is_nan_like guard in T6.")


def test_t6_is_warning_not_error():
    # the spec does not forbid alt contigs, so a legitimately alt library must still pass
    df = df_from(make_row(guide_chr="chr6_GL000252v2_alt",
                          intended_target_chr="chr6_GL000252v2_alt"))
    issues = validate_df(df)
    assert not errors(issues, "guide_chr"), (
        "T6 must be warn-severity: the spec permits alt contigs, so flagging them as "
        f"errors would fail valid libraries. Got errors: "
        f"{[i.message for i in errors(issues, 'guide_chr')]}")


# ══════════════════════════════════════════════════════════════════════════════
# E1-E3 — genomic_element
# ══════════════════════════════════════════════════════════════════════════════

def test_e1_empty_genomic_element_for_targeting_is_error():
    df = df_from(make_row(genomic_element=""))
    issues = validate_df(df)
    assert has_error(issues, "genomic_element", "empty")


def test_e2_invalid_genomic_element_is_error():
    df = df_from(make_row(genomic_element="TSS"))
    issues = validate_df(df)
    assert has_error(issues, "genomic_element", "invalid")


def test_e2_distal_element_underscore_is_error():
    # "distal_element" (with underscore) is NOT a valid enum value
    df = df_from(make_row(genomic_element="distal_element",
                          intended_target_name="chr1:900-1100"))
    issues = validate_df(df)
    assert has_error(issues, "genomic_element", "invalid")


def test_e2b_enhancer_is_warning():
    df = df_from(make_row(genomic_element="enhancer",
                          intended_target_name="chr1:900-1100"))
    issues = validate_df(df)
    assert has_warning(issues, "genomic_element", "distal element")


def test_e2_all_valid_genomic_elements_clean():
    valid_elements = [
        ("promoter",      "ENSG00000000001"),
        ("gene",          "ENSG00000000002"),
        ("splice site",   "ENSG00000000003"),
        ("distal element","chr1:900-1100"),
        ("enhancer",      "chr1:900-1200"),   # valid enum; will warn but not error
        ("insulator",     "chr1:900-1300"),
        ("silencer",      "chr1:900-1400"),
        ("variant",       "ENSG00000000007"),
    ]
    rows = [
        make_row(guide_id=f"g{i}", genomic_element=ge, intended_target_name=itn)
        for i, (ge, itn) in enumerate(valid_elements)
    ]
    df = df_from(*rows)
    issues = validate_df(df)
    assert not has_error(issues, "genomic_element", "invalid")


def test_e3_non_targeting_row_with_genomic_element_is_error():
    df = df_from(make_nt_row(genomic_element="promoter"))
    issues = validate_df(df)
    assert has_error(issues, "genomic_element", "non-targeting")


def test_e3_safe_targeting_with_genomic_element_allowed():
    # safe-targeting guides target real loci; element fields are permitted
    df = df_from(make_safe_row(genomic_element="promoter"))
    issues = validate_df(df)
    assert not has_error(issues, "genomic_element", "non-targeting")


# ══════════════════════════════════════════════════════════════════════════════
# N1-N3 — intended_target_name
# ══════════════════════════════════════════════════════════════════════════════

def test_n1_empty_itn_for_targeting_is_error():
    df = df_from(make_row(intended_target_name=""))
    issues = validate_df(df)
    assert has_error(issues, "intended_target_name", "empty")


def test_n2_promoter_requires_ensg():
    df = df_from(make_row(genomic_element="promoter",
                          intended_target_name="chr1:900-1100"))
    issues = validate_df(df)
    assert has_error(issues, "intended_target_name", "ENSG")


def test_n2_promoter_valid_ensg_clean():
    df = df_from(make_row(genomic_element="promoter",
                          intended_target_name="ENSG00000000001"))
    issues = validate_df(df)
    assert not errors(issues, "intended_target_name")


def test_n2_gene_requires_ensg():
    df = df_from(make_row(genomic_element="gene",
                          intended_target_name="TP53"))
    issues = validate_df(df)
    assert has_error(issues, "intended_target_name", "ENSG")


def test_n2_splice_site_requires_ensg():
    df = df_from(make_row(genomic_element="splice site",
                          intended_target_name="chr1:900-1100"))
    issues = validate_df(df)
    assert has_error(issues, "intended_target_name", "ENSG")


def test_n2_distal_element_requires_coord():
    df = df_from(make_row(genomic_element="distal element",
                          intended_target_name="ENSG00000000001"))
    issues = validate_df(df)
    assert has_error(issues, "intended_target_name", "chr:start-end")


def test_n2_distal_element_valid_coord_clean():
    df = df_from(make_row(genomic_element="distal element",
                          intended_target_name="chr1:900-1100"))
    issues = validate_df(df)
    assert not errors(issues, "intended_target_name")


def test_n2_insulator_requires_coord():
    df = df_from(make_row(genomic_element="insulator",
                          intended_target_name="ENSG00000000001"))
    issues = validate_df(df)
    assert has_error(issues, "intended_target_name", "chr:start-end")


def test_n2_silencer_requires_coord():
    df = df_from(make_row(genomic_element="silencer",
                          intended_target_name="ENSG00000000001"))
    issues = validate_df(df)
    assert has_error(issues, "intended_target_name", "chr:start-end")


def test_n2_variant_element_no_format_restriction():
    # genomic_element='variant' is in neither GE_NEEDS_ENSG nor GE_NEEDS_COORD
    df = df_from(make_row(genomic_element="variant",
                          intended_target_name="rs12345"))
    issues = validate_df(df)
    assert not errors(issues, "intended_target_name")


def test_n3_non_targeting_with_itn_is_error():
    df = df_from(make_nt_row(intended_target_name="ENSG00000000001"))
    issues = validate_df(df)
    assert has_error(issues, "intended_target_name", "non-targeting")


def test_n3_safe_targeting_with_itn_allowed():
    df = df_from(make_safe_row(intended_target_name="ENSG00000000002"))
    issues = validate_df(df)
    assert not has_error(issues, "intended_target_name", "non-targeting")


# ══════════════════════════════════════════════════════════════════════════════
# C1 / C1b — target coordinates presence and type
# ══════════════════════════════════════════════════════════════════════════════

def test_c1_empty_target_chr_is_error():
    df = df_from(make_row(intended_target_chr=""))
    issues = validate_df(df)
    assert has_error(issues, "intended_target_chr", "empty")


def test_c1_empty_target_start_is_error():
    df = df_from(make_row(intended_target_start=""))
    issues = validate_df(df)
    assert has_error(issues, "intended_target_start", "empty")


def test_c1_empty_target_end_is_error():
    df = df_from(make_row(intended_target_end=""))
    issues = validate_df(df)
    assert has_error(issues, "intended_target_end", "empty")


def test_c1_non_targeting_with_target_chr_is_error():
    df = df_from(make_nt_row(intended_target_chr="chr1"))
    issues = validate_df(df)
    assert has_error(issues, "intended_target_chr", "non-targeting")


def test_c1_non_targeting_with_target_start_is_error():
    df = df_from(make_nt_row(intended_target_start="100"))
    issues = validate_df(df)
    assert has_error(issues, "intended_target_start", "non-targeting")


def test_c1b_non_numeric_target_start_is_error():
    # Bug 2 regression: non-numeric values must be flagged
    df = df_from(make_row(intended_target_start="start"))
    issues = validate_df(df)
    assert has_error(issues, "intended_target_start", "non-numeric")


def test_c1b_non_numeric_target_end_is_error():
    df = df_from(make_row(intended_target_end="end"))
    issues = validate_df(df)
    assert has_error(issues, "intended_target_end", "non-numeric")


def test_c1b_whitespace_target_start_treated_as_empty():
    df = df_from(make_row(intended_target_start="  "))
    issues = validate_df(df)
    assert has_error(issues, "intended_target_start", "empty")


# ══════════════════════════════════════════════════════════════════════════════
# C2 — coordinate span consistency (per-row)
# ══════════════════════════════════════════════════════════════════════════════

def test_c2_target_chr_mismatch_is_error():
    df = df_from(make_row(guide_chr="chr1", intended_target_chr="chr2"))
    issues = validate_df(df)
    assert has_error(issues, "intended_target_chr", "≠")


def test_c2_target_start_greater_than_guide_start_large_is_error():
    # intended_target_start > guide_start by more than 50 bp → error
    df = df_from(make_row(guide_start="1000", intended_target_start="1200",
                          intended_target_end="1500"))
    issues = validate_df(df)
    assert has_error(issues, "intended_target_start", "> guide_start")


def test_c2_target_start_greater_than_guide_start_small_is_warning():
    # overshoot ≤ 50 bp → warning only
    df = df_from(make_row(guide_start="1000", intended_target_start="1030",
                          intended_target_end="1500"))
    issues = validate_df(df)
    assert has_warning(issues, "intended_target_start", "≤ 50 bp")
    assert not has_error(issues, "intended_target_start", "> guide_start")


def test_c2_target_end_less_than_guide_end_large_is_error():
    df = df_from(make_row(guide_end="1100", intended_target_end="900"))
    issues = validate_df(df)
    assert has_error(issues, "intended_target_end", "< guide_end")


def test_c2_target_end_less_than_guide_end_small_is_warning():
    df = df_from(make_row(guide_end="1020", intended_target_end="1000"))
    issues = validate_df(df)
    assert has_warning(issues, "intended_target_end", "≤ 50 bp")
    assert not has_error(issues, "intended_target_end", "< guide_end")


def test_c2_coords_span_guides_clean():
    df = df_from(make_row(
        guide_start="1000", guide_end="1020",
        intended_target_start="900", intended_target_end="1100",
    ))
    issues = validate_df(df)
    assert not errors(issues, "intended_target_start")
    assert not errors(issues, "intended_target_end")
    assert not errors(issues, "intended_target_chr")


# ══════════════════════════════════════════════════════════════════════════════
# C3 — per-group coordinate consistency
# ══════════════════════════════════════════════════════════════════════════════

def test_c3_per_guide_coords_flagged():
    # Both guides in the same prefix group (GENE1_TSS) have coords that only
    # cover their own position — neither spans the full group.
    rows = [
        make_row(guide_id="GENE1_TSS_1", guide_start="1000", guide_end="1020",
                 intended_target_start="1000", intended_target_end="1020"),
        make_row(guide_id="GENE1_TSS_2", guide_start="1100", guide_end="1120",
                 intended_target_start="1100", intended_target_end="1120"),
    ]
    df = df_from(*rows)
    issues = validate_df(df)
    assert has_error(issues, "intended_target_start", "span all guides")


def test_c3_per_group_coords_clean():
    # Both guides share coords that span the entire group.
    rows = [
        make_row(guide_id="GENE1_TSS_1", guide_start="1000", guide_end="1020",
                 intended_target_start="900", intended_target_end="1200"),
        make_row(guide_id="GENE1_TSS_2", guide_start="1100", guide_end="1120",
                 intended_target_start="900", intended_target_end="1200"),
    ]
    df = df_from(*rows)
    issues = validate_df(df)
    assert not has_error(issues, "intended_target_start", "span all guides")


def test_c3_uniform_prefix_multichrom_not_flagged():
    # Random-screen naming: after stripping _N the guide_id collapses to one prefix
    # ("Lib_Crop") spanning multiple chromosomes — the prefix does NOT encode elements,
    # so C3 must not emit a spurious span error. Regression for the released
    # *_Random_Screen_Crop_<N> guide files (K562 IGVFFI7272YLVI / WTC11 IGVFFI0580WJFK).
    rows = [
        make_row(guide_id="Lib_Crop_1", guide_chr="chr1", guide_start="1000",
                 guide_end="1020", intended_target_start="1000", intended_target_end="1020"),
        make_row(guide_id="Lib_Crop_2", guide_chr="chr5", guide_start="900000",
                 guide_end="900020", intended_target_start="900000", intended_target_end="900020"),
    ]
    issues = validate_df(df_from(*rows))
    assert not has_error(issues, "intended_target_start", "span all guides")


def test_c3_uniform_prefix_huge_span_not_flagged():
    # Same prefix, one chromosome, but guides > C3_MAX_ELEMENT_BP apart → not one element.
    rows = [
        make_row(guide_id="Lib_Crop_1", guide_chr="chr1", guide_start="1000",
                 guide_end="1020", intended_target_start="1000", intended_target_end="1020"),
        make_row(guide_id="Lib_Crop_2", guide_chr="chr1", guide_start="5000000",
                 guide_end="5000020", intended_target_start="5000000", intended_target_end="5000020"),
    ]
    issues = validate_df(df_from(*rows))
    assert not has_error(issues, "intended_target_start", "span all guides")


def test_c3_edge_overhang_within_tolerance_clean():
    # A designed element window with an edge guide overhanging by ≤ 50 bp is not an error.
    rows = [
        make_row(guide_id="GENE1_TSS_1", guide_start="895", guide_end="915",
                 intended_target_start="900", intended_target_end="1400"),
        make_row(guide_id="GENE1_TSS_2", guide_start="1300", guide_end="1320",
                 intended_target_start="900", intended_target_end="1400"),
    ]
    issues = validate_df(df_from(*rows))
    assert not has_error(issues, "intended_target_start", "span all guides")


# ══════════════════════════════════════════════════════════════════════════════
# P1 — putative_target_genes
# ══════════════════════════════════════════════════════════════════════════════

def test_p1_positive_control_distal_empty_ptg_is_error():
    df = df_from(make_pc_row(putative_target_genes=""))
    issues = validate_df(df)
    assert has_error(issues, "putative_target_genes", "positive control")


def test_p1_positive_control_distal_valid_ptg_clean():
    df = df_from(make_pc_row(putative_target_genes='["ENSG00000000003"]'))
    issues = validate_df(df)
    assert not errors(issues, "putative_target_genes")


def test_p1_ptg_bare_value_not_array_is_error():
    # spec types putative_target_genes as string[]; a bare ENSG id must be an array literal
    df = df_from(make_pc_row(putative_target_genes="ENSG00000000003"))
    issues = validate_df(df)
    assert has_error(issues, "putative_target_genes", "array")


def test_p1_positive_control_promoter_no_ptg_required():
    # For promoter (not a distal-type element), putative_target_genes is not required
    df = df_from(make_row(type="positive control", genomic_element="promoter",
                          intended_target_name="ENSG00000000001",
                          putative_target_genes=""))
    issues = validate_df(df)
    assert not has_error(issues, "putative_target_genes", "positive control")


def test_p1_ptg_multiple_ensg_array_clean():
    df = df_from(make_pc_row(putative_target_genes='["ENSG00000000003", "ENSG00000000004"]'))
    issues = validate_df(df)
    assert not errors(issues, "putative_target_genes")


def test_p1_ptg_non_ensg_value_warns():
    df = df_from(make_pc_row(putative_target_genes='["TP53"]'))
    issues = validate_df(df)
    assert has_warning(issues, "putative_target_genes", "ENSG")


def test_p1_targeting_row_no_ptg_not_required():
    df = df_from(make_row(type="targeting", putative_target_genes=""))
    issues = validate_df(df)
    assert not errors(issues, "putative_target_genes")


# ══════════════════════════════════════════════════════════════════════════════
# D1 — description
# ══════════════════════════════════════════════════════════════════════════════

def test_d1_absent_description_warns():
    # When description is absent, validate_df() backfills it as empty strings
    # before the D1 check runs, so the warning says "all empty" (not "absent").
    # The S1 structural warning also fires for the missing optional column.
    df = df_from(make_row())
    df = df.drop(columns=["description"])
    issues = validate_df(df)
    assert has_warning(issues, "description")
    assert has_warning(issues, "structure", "description")


def test_d1_all_empty_description_warns():
    df = df_from(make_row(description=""))
    issues = validate_df(df)
    assert has_warning(issues, "description", "all empty")


def test_d1_filled_description_clean():
    df = df_from(make_row(description="GENE1_TSS"))
    issues = validate_df(df)
    assert not warnings(issues, "description")


# ══════════════════════════════════════════════════════════════════════════════
# _is_nan_like — whitespace / NaN-string handling (bug 5 regression)
# ══════════════════════════════════════════════════════════════════════════════

def test_nan_strings_treated_as_empty():
    """Fields containing 'NA', 'None', 'null', etc. should be treated as missing."""
    for nan_val in ["NA", "na", "N/A", "None", "null", "NaN", "nan"]:
        df = df_from(make_row(intended_target_name=nan_val))
        issues = validate_df(df)
        assert has_error(issues, "intended_target_name", "empty"), \
            f"'{nan_val}' should be treated as empty for intended_target_name"


def test_whitespace_only_treated_as_empty():
    """Whitespace-only values must be treated as missing (bug 5)."""
    df = df_from(make_row(intended_target_name="   "))
    issues = validate_df(df)
    assert has_error(issues, "intended_target_name", "empty")


# ══════════════════════════════════════════════════════════════════════════════
# Full clean file — no issues expected (except optional-column warnings)
# ══════════════════════════════════════════════════════════════════════════════

def test_fully_valid_file_no_errors():
    """A correctly-formed file should produce zero error-severity issues."""
    rows = [
        make_row(guide_id="g1", description="GENE1_TSS"),
        make_row(guide_id="g2", genomic_element="distal element",
                 intended_target_name="chr1:900-1100", description="ENH1"),
        make_nt_row(guide_id="nt1"),
        make_safe_row(guide_id="st1"),
        make_pc_row(guide_id="pc1"),
    ]
    df = df_from(*rows)
    issues = validate_df(df)
    assert not errors(issues), f"Expected no errors; got: {[i.message for i in errors(issues)]}"


def test_fully_valid_file_enhancer_warning_only():
    """'enhancer' is valid but should produce a warning, not an error."""
    df = df_from(make_row(genomic_element="enhancer",
                          intended_target_name="chr1:900-1100"))
    issues = validate_df(df)
    assert not errors(issues, "genomic_element")
    assert has_warning(issues, "genomic_element", "distal element")

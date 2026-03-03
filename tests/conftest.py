"""Shared test helpers for validate_grna_file tests."""

import sys
import os

# Ensure the project root is on the path so validate_grna_file can be imported
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pandas as pd
import pytest

from validate_grna_file import validate_df, Issue


# ── Row factories ──────────────────────────────────────────────────────────────

def make_row(**overrides) -> dict:
    """Return a fully valid targeting row dict, with keyword overrides applied."""
    base = {
        "guide_id":               "g1",
        "spacer":                 "ACGTACGTACGT",
        "targeting":              "True",
        "type":                   "targeting",
        "guide_chr":              "chr1",
        "guide_start":            "1000",
        "guide_end":              "1020",
        "strand":                 "+",
        "pam":                    "NGG",
        "genomic_element":        "promoter",
        "intended_target_name":   "ENSG00000000001",
        "intended_target_chr":    "chr1",
        "intended_target_start":  "900",
        "intended_target_end":    "1100",
        "description":            "GENE1_TSS",
    }
    base.update(overrides)
    return base


def make_nt_row(**overrides) -> dict:
    """Return a fully valid non-targeting row dict, with keyword overrides applied."""
    base = {
        "guide_id":               "nt1",
        "spacer":                 "TTTTTTTTTTTT",
        "targeting":              "False",
        "type":                   "non-targeting",
        "guide_chr":              "",
        "guide_start":            "",
        "guide_end":              "",
        "strand":                 "",
        "pam":                    "NGG",
        "genomic_element":        "",
        "intended_target_name":   "",
        "intended_target_chr":    "",
        "intended_target_start":  "",
        "intended_target_end":    "",
        "description":            "NT_control",
    }
    base.update(overrides)
    return base


def make_safe_row(**overrides) -> dict:
    """Return a fully valid safe-targeting row dict, with keyword overrides applied."""
    base = {
        "guide_id":               "st1",
        "spacer":                 "CCCCCCCCCCCC",
        "targeting":              "False",
        "type":                   "safe-targeting",
        "guide_chr":              "chr2",
        "guide_start":            "5000",
        "guide_end":              "5020",
        "strand":                 "-",
        "pam":                    "NGG",
        "genomic_element":        "promoter",
        "intended_target_name":   "ENSG00000000002",
        "intended_target_chr":    "chr2",
        "intended_target_start":  "4900",
        "intended_target_end":    "5100",
        "description":            "SAFE_ctrl",
    }
    base.update(overrides)
    return base


def make_pc_row(**overrides) -> dict:
    """Return a valid positive-control row targeting a distal element."""
    base = {
        "guide_id":               "pc1",
        "spacer":                 "GGGGGGGGGGGG",
        "targeting":              "True",
        "type":                   "positive control",
        "guide_chr":              "chr3",
        "guide_start":            "2000",
        "guide_end":              "2020",
        "strand":                 "+",
        "pam":                    "NGG",
        "genomic_element":        "distal element",
        "intended_target_name":   "chr3:1900-2100",
        "intended_target_chr":    "chr3",
        "intended_target_start":  "1900",
        "intended_target_end":    "2100",
        "putative_target_genes":  "ENSG00000000003",
        "description":            "ENH_ctrl",
    }
    base.update(overrides)
    return base


def df_from(*rows) -> pd.DataFrame:
    """Build a DataFrame from one or more row dicts, filling missing values with ''."""
    return pd.DataFrame(list(rows)).fillna("")


# ── Issue query helpers ────────────────────────────────────────────────────────

def errors(issues: list, field: str = None) -> list:
    """Return error-severity issues, optionally filtered to a specific field."""
    return [i for i in issues if i.severity == "error"
            and (field is None or i.field == field)]


def warnings(issues: list, field: str = None) -> list:
    """Return warn-severity issues, optionally filtered to a specific field."""
    return [i for i in issues if i.severity == "warn"
            and (field is None or i.field == field)]


def has_error(issues: list, field: str, fragment: str = "") -> bool:
    """True if there is an error for *field* whose message contains *fragment*."""
    return any(
        fragment.lower() in i.message.lower()
        for i in issues
        if i.severity == "error" and i.field == field
    )


def has_warning(issues: list, field: str, fragment: str = "") -> bool:
    """True if there is a warning for *field* whose message contains *fragment*."""
    return any(
        fragment.lower() in i.message.lower()
        for i in issues
        if i.severity == "warn" and i.field == field
    )

#!/usr/bin/env python3
"""
Run the vendored IGVF-DACC guide RNA frictionless check on one or more files.

Usage:
    python3 run_external_check.py output/IGVFFI0580WJFK.tsv.gz [...]

Exit code 0 if all files pass; non-zero otherwise.

The check is vendored from:
  https://github.com/IGVF-DACC/checkfiles/blob/dev/src/checkfiles/guide_rna_sequences_check.py
Run `make update-external` to fetch the latest version.
"""

import sys
import os

# Make the vendored module importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "external"))

try:
    from frictionless import validate
    from guide_rna_sequences_check import GuideRnaSequencesCheck
except ImportError as e:
    print(f"ERROR: {e}")
    print("Install dependencies with: pip install frictionless")
    sys.exit(2)


def check_file(filepath: str) -> bool:
    """Run the IGVF-DACC check on *filepath*. Returns True if valid."""
    print(f"\n{'='*60}")
    print(f"External check: {filepath}")
    print(f"{'='*60}")

    report = validate(filepath, checks=[GuideRnaSequencesCheck()])

    if report.valid:
        stats = report.stats
        print(f"  PASSED  ({stats.get('rows', '?')} rows, 0 errors)")
        return True

    # Print errors grouped by field
    error_rows = []
    for task in report.tasks:
        for err in task.errors:
            error_rows.append(err)

    print(f"  FAILED  ({len(error_rows)} error(s))")
    # Show up to 20 errors, then summarise
    shown = error_rows[:20]
    for err in shown:
        row_num = getattr(err, 'row_number', '?')
        field = getattr(err, 'field_name', '?')
        note = getattr(err, 'note', str(err))
        print(f"    row {row_num}, field={field}: {note}")
    if len(error_rows) > 20:
        print(f"    ... and {len(error_rows) - 20} more errors")

    return False


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    files = sys.argv[1:]
    results = {f: check_file(f) for f in files}

    print(f"\n{'='*60}")
    print("Summary:")
    all_pass = True
    for f, ok in results.items():
        status = "PASS" if ok else "FAIL"
        print(f"  {status}  {f}")
        if not ok:
            all_pass = False

    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()

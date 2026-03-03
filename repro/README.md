# Bug Repro Fixtures (Failing Cases Only)

Each file below is intentionally tiny and crafted to trigger one specific bug from the review.

## 1) `fix_grna_file.py` lowercases `targeting`

Fixture: `repro/fix_grna_targeting_case.tsv`

Repro command:

```bash
python3 - <<'PY'
import csv

with open("repro/fix_grna_targeting_case.tsv", newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

# same logic as fix_grna_file.py lines 106-107
for r in rows:
    r["targeting"] = r["targeting"].lower()
    if r["type"] == "safe-targeting":
        r["targeting"] = "false"

for r in rows:
    print(r["guide_id"], r["type"], r["targeting"])
PY
```

Expected bad behavior: output is `true/false` (lowercase), not spec-required `True/False`.

## 2) Validator misses non-numeric coordinates in C2

Fixture: `repro/validate_nonnumeric.tsv`

Setup once (offline GTF stub so validator does not download):

```bash
mkdir -p input
printf 'chr1\tsrc\tgene\t1\t1000\t.\t+\t.\tgene_id "ENSG000001.1"; gene_name "GENE1";\n' | gzip -c > input/IGVFFI9573KOZR.gtf.gz
```

Repro command:

```bash
python3 validate_grna_file.py repro/validate_nonnumeric.tsv
```

Expected bad behavior: no explicit error that `guide_start` / `intended_target_start` are non-numeric.

## 3) Interactive `fix_coords_recompute` crashes on NaN numeric conversion

Fixtures: `repro/fix_coords_nan.tsv`, `repro/problems_fix_coords.json`

Repro command:

```bash
printf 'a\n' | python3 fix_interactive.py repro/fix_coords_nan.tsv --problems repro/problems_fix_coords.json --output repro/out.tsv.gz
```

Expected bad behavior: crashes with `ValueError: cannot convert float NaN to integer`.

## 4) `fix_easy_files` NaN-to-int conversion crashes

Repro command:

```bash
python3 - <<'PY'
# same failing cast pattern used by fix_easy_files after numeric coercion
x = float("nan")
print(int(x))
PY
```

Expected bad behavior: same `ValueError: cannot convert float NaN to integer`.

## 5) Whitespace-only values are not treated as empty

Fixture: `repro/validate_whitespace.tsv`

Repro command:

```bash
python3 validate_grna_file.py repro/validate_whitespace.tsv
```

Expected bad behavior: whitespace in `intended_target_start/end` is not flagged as empty in C1; downstream numeric checks also miss it.

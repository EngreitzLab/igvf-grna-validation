# CLAUDE.md — igvf-grna-validation

Agent guide for this repo. Human overview in [README.md](README.md); domain reference in
[KNOWLEDGE.md](KNOWLEDGE.md); the authoritative format spec is the PDF in `input/`.

## What this repo is

A **validator plus a set of one-way fixers** for IGVF Per-Guide Metadata (gRNA library) TSVs,
run before submission to the IGVF Data Portal. Two kinds of code, and the split matters:

| Kind | Files | Contract |
|---|---|---|
| **Checkers** | `validate_grna_file.py`, `verify_guide_coords_against_genome.py`, `run_external_check.py` | **Never write to the input.** Read, report, exit. |
| **Fixers** | `fixers/*.py` | Read one input, write a corrected copy to `output/`. Never edit in place. |

Everything else is data: `input/` (files being fixed, plus the spec PDF and downloaded
genome data — see debt #5), `output/` (corrected results), `problems/` (machine-readable
reports), `repro/` (small fixtures reproducing past bugs), `external/` (vendored DACC
checker), `tests/`.

---

## Hard rules

These are the ones that have actually caused damage. Read them before running anything.

1. **Run every script from the repo root.** Scripts resolve `input/…` and `output/…`
   relative to the working directory, so `python3 fixers/fix_X.py` works and
   `cd fixers && python3 fix_X.py` silently fails or writes to the wrong place.

2. **Never `import` a fixer to inspect it.** `fixers/fix_grna_file.py`,
   `fix_easy_files.py` and `fix_IGVFFI9754AGFB.py` have **no `if __name__ == "__main__"`
   guard** — importing them executes them and overwrites tracked files in `output/`. To
   inspect one, read it or `ast.parse` it. To check that its imports resolve, write a
   throwaway probe next to it rather than importing the module.

3. **`git status` after running any fixer.** Fixers overwrite `output/*.tsv.gz`, which is
   git-tracked. Regenerating an unrelated output is a silent, committable side effect.
   Revert what you did not mean to change.

4. **Never commit genome reference data.** The GENCODE GTF (~50 MB) and the GRCh38 FASTA
   (~983 MB) are downloaded on demand and today live in `input/`, each with its own
   `.gitignore` line. Before adding a third reference file, do debt #5 (move all of them to
   a wholesale-ignored `reference/`) rather than adding a fourth ignore line.

5. **Keep `validate_grna_file.py` genome-free.** It may use the GENCODE GTF (auto-downloaded,
   small) but must never require the FASTA. Any check that needs sequence goes in
   `verify_guide_coords_against_genome.py`. This keeps the common path fast and dependency-light.

6. **Prefer a TSV-only check over a genome-backed one.** If a property is derivable from the
   file alone, it belongs in the validator. Check T5c is the worked example: PAM-in-window
   looks like it needs the reference, but `guide_end - guide_start` vs `len(spacer)` decides
   it — so it is a validator check, and the genome was used only to *prove the rule* once.

---

## Layout

```
validate_grna_file.py                 the validator — spec checks, GTF only, no FASTA
verify_guide_coords_against_genome.py genome-backed verification (needs the GRCh38 FASTA)
run_external_check.py                 wrapper around the vendored IGVF-DACC checker
fixers/                               one-way correction scripts (see below)
tests/                                pytest suite for validate_df(); conftest.py row factories
input/                                files being corrected, the spec PDF, and the
                                      downloaded GTF + FASTA (both gitignored) — see debt #5
output/                               corrected results (git-tracked) + changelogs
problems/                             --json-out reports, read by fixers/fix_interactive.py
repro/                                minimal TSVs reproducing past validator bugs
external/                             vendored IGVF-DACC checker (`make update-external`)
```

`input/` is 993 MB, essentially all of it re-downloadable reference data.

---

## Script conventions

**Naming.** Say what the script *does*, not which file it was written for.
`fix_300genes_promoter_library.py` is the pattern to follow;
`fix_IGVFFI0580WJFK.py` is not — an accession tells the reader nothing. When adding a fixer,
name it for the correction or the library, and record the accession in the docstring.

**Interface.** New fixers take `argparse` arguments — input path, output path, and any
reference paths — rather than module-level `INPUT_FILE`/`OUTPUT_FILE` constants. Constants
make a script single-use and unrunnable on a second file.

**Structure.** Every script gets a `main()` and an `if __name__ == "__main__":` guard, so it
can be imported without side effects.

**Logging.** A fixer prints what it changed, per step, with before → after values, and
supports `--changelog <path>` to persist that log. A silent fixer is unreviewable.

**Self-verification.** End each correction step with an `assert` that the invariant now holds
(`fixers/fix_300genes_promoter_library.py` asserts every row spans exactly its spacer after
the PAM trim). Do not rely on re-running the validator to catch your own bug.

**Shared helpers** live in `fixers/fix_interactive.py` and are imported by the per-file
fixers. Fixers that also need validator internals must bootstrap the repo root:

```python
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from validate_grna_file import load_gene_map, GTF_PATH, GTF_URL, download_gtf
```

---

## Adding a validator check

1. **Anchor it in the spec.** Quote the sentence from the PDF in a comment. If the spec does
   not require it, the check is a `warn`, never an `error` — an over-strict error fails valid
   libraries. T6 (non-primary contigs) is `warn` for exactly this reason.
2. **Label it** in the existing scheme (`S*` structure, `T*` guide fields, `E*`
   genomic_element, `N*` intended_target_name, `C*` coordinates, `P*` putative_target_genes,
   `D*` description). Place it in numeric order in `validate_df()`.
3. **Do not gate on `targeting == True` reflexively.** Safe-targeting rows carry coordinates
   with `targeting=False`; several checks must cover them. This was a real bug source — see
   `test_t6_applies_to_non_targeting_rows_too`.
4. **Add it to the checks table in [README.md](README.md).**
5. **Test it** in `tests/test_validate.py`: partition the input space, one case per
   subdomain plus boundaries, and give every assertion a message saying what was expected,
   what was received, and what to check. Follow the strategy comment blocks above the T5c
   and T6 sections.
6. **Run the suite** — `make test` or `pytest tests/ -q`. If a new check makes a
   pre-existing test fail, first ask whether the *fixture* was wrong: T5c exposed test rows
   with 12 bp spacers and 20 bp coordinate spans, which had never been coherent.

---

## Environment

Local runs use the parent project's pixi environment:

```bash
pixi run --manifest-path ../../pixi.toml python validate_grna_file.py <file>
```

Plain `python3` also works given pandas + requests. `requirements.txt` is currently wrong
(see debt #4) — trust the imports, not that file. Never use system pip or brew.

---

## Known debt

Numbered so it can be worked off incrementally. Do not fix these opportunistically inside an
unrelated change; each deserves its own commit.

1. **3 fixers lack `__main__` guards** — `fixers/fix_grna_file.py`, `fix_easy_files.py`,
   `fix_IGVFFI9754AGFB.py`. Wrap their top-level bodies in `main()`. This is the highest-value
   fix: it removes the import-executes-and-overwrites hazard in rule 2.
2. **6 of 8 fixers hardcode input/output paths** as module constants. Convert to `argparse`.
3. **`_is_nan_like` is defined twice** — in `validate_grna_file.py` and in
   `fixers/fix_interactive.py`. Two copies of the NaN-string vocabulary will drift. Promote one
   to a shared module and drop the `_` prefix while moving it (a leading underscore on a
   cross-module import signals private but protects nothing).
4. **`requirements.txt` is inaccurate** — it lists `frictionless`, which nothing imports, and
   omits `requests`, which `validate_grna_file.py` needs.
5. **`input/` conflates three things**: the spec PDF, auto-downloaded reference data (993 MB
   of the directory), and actual input TSVs. Move the reference data into a new `reference/`,
   update `GTF_PATH` in `validate_grna_file.py` and the `--fasta` default in
   `verify_guide_coords_against_genome.py`, then replace the two per-file `.gitignore` lines
   with a single `reference/`. Rules 1 and 4 and the layout block above all need updating in
   the same commit.
6. **Accession-named fixers** (`fix_IGVFFI*.py`) should be renamed for what they correct, per
   the naming convention above.
7. **`verify_guide_coords_against_genome.py` does not auto-download its FASTA**, unlike the
   validator's GTF. Add the same download-on-first-run path so the script is self-sufficient.

# PhyloHPC Master Reorganization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reorganize the `phylohpc` repo per good code practices — fix the broken step-4 visualization, consolidate config, untrack build artifacts, group scattered scripts, and delete dead code — without breaking either pipeline stack or the test suite.

**Architecture:** All work happens on a new branch `reorg-cleanup` off `master`; `master` stays untouched as the archive. Changes are mechanical (git mv / git rm / reference updates) except Task 1, which reconstructs a deleted script under TDD. Each task ends green (`pytest tests/`) with no dangling references.

**Tech Stack:** Nextflow (`.nf`, step1/2 primary), Snakemake (`.smk`, step4 authoritative), Python 3.10 (`workflow/*.py`, pytest), R (analysis scripts), git submodule `phylogeny`.

## Global Constraints

- **NEVER touch the `phylogeny/` submodule** or its contents — the parent repo only stores its gitlink commit.
- **`pytest tests/` must pass (or improve)** after every task. Baseline: `tests/test_hog_hierarchy_report.py` currently fails at collection (missing module) — Task 1 fixes it; all other tests currently pass and must stay passing.
- **Two pipeline stacks are both live** — Nextflow owns step1/2 + SLURM; `step4_ancestry.smk` is the *authoritative* step-4 impl. Do not remove either stack. `step4.ancestry.nf` stays (documented-incomplete, out of scope).
- **Use `git mv` / `git rm`** for tracked files (preserve history), never raw `rm` on tracked paths.
- **Every reference update must be grep-verified**: after moving/deleting a path, `grep -rn '<oldname>'` across `README.md docs/ instructions/ config/ *.nf *.smk workflow/ .github/ *.R` must return no live references (historical dev logs `WHATIDID.md` / `TODOs.md` are exempt — leave them as-is).
- **Do not modify `phylogeny/`, `data/` input files, `img/` legitimate assets, `docs/step2_report.html`** (published gh-pages demo), or the `.nf`/`.smk` pipeline logic beyond the specific reference/path edits each task names.
- Commit after each task with a conventional-commit message. Do not merge to `master` until the final review passes.

## Non-Goals (explicitly out of scope)

- Consolidating `check_job.py` vs `check_job.v2.py` (test-coupled rename; noted as Minor follow-up).
- Resolving the `nextflow.config` params ↔ `config/*.yaml` duplication (design change, not a move; noted for follow-up).
- Deleting `report_step1.py` or `visualize_ancestry.py` (both have passing tests; keep).
- Fixing `step4.ancestry.nf` incompleteness (documented, separate work).

---

### Task 0: Create the working branch

**Files:** none (git branch only)

- [ ] **Step 1: Confirm clean tree and branch off master**

Run:
```bash
cd /home/grygoriyzolotarov/Documents/projects/phylohpc
git status --short            # expect empty
git checkout -b reorg-cleanup
git branch --show-current     # expect: reorg-cleanup
```
Expected: on branch `reorg-cleanup`, working tree clean.

- [ ] **Step 2: Record baseline test state**

Run: `python -m pytest tests/ -q 2>&1 | tail -20`
Expected: all tests pass EXCEPT `tests/test_hog_hierarchy_report.py` which errors at collection with `ModuleNotFoundError: No module named 'visualize_hog_hierarchy'` (or `ImportError`). Record the pass count; Task 1 must restore this test, all others stay green.

---

### Task 1: Recreate `workflow/visualize_hog_hierarchy.py` (fix broken step-4 viz)

**Context:** Commit `990531a` deleted this file but left every caller referencing it: `step4.ancestry.nf:156`, `step4_ancestry.smk:137,160`, `workflow/build_hog_report.py:36`, `tests/test_hog_hierarchy_report.py:12`. The callers were rewired to a **newer** interface than the last-committed version — the recovered blob (`git show 990531a~1:workflow/visualize_hog_hierarchy.py`) is a useful starting scaffold (Sankey HTML generator, ~622 lines, has `load_tsvs`, `build_data`, `HTML_TEMPLATE`, `main`) but LACKS the symbols/args the current callers need. The test + callers are the spec.

**Files:**
- Create: `workflow/visualize_hog_hierarchy.py`
- Reference (read, do not modify): `tests/test_hog_hierarchy_report.py`, `workflow/build_hog_report.py`, `workflow/link_hog_levels.py`, `step4_ancestry.smk`, `step4.ancestry.nf`

**Interface the recreated module MUST expose (the contract):**
- `build_data(stats: pd.DataFrame, links: pd.DataFrame, levels: list[str]) -> dict` — per-HG Sankey data structure (from the recovered blob).
- `og_short_py(og_id: str) -> str` — shorten an OG id, **preserving a `:like:<label>` suffix** (per `test_og_short_py_preserves_like_label`). Required by `tests/test_hog_hierarchy_report.py` import.
- `parse_newick_tree(...)` — parse a pruned newick into the structure `build_hog_report.py` consumes (imported at `build_hog_report.py:36` alongside `HTML_TEMPLATE, build_data`).
- `HTML_TEMPLATE` — module-level HTML string.
- `main()` with an argparse CLI accepting **`--links` (nargs+), `--stats` (nargs+), `--trees` (nargs+), `--levels`, `--output`** — note `--trees` is new vs the recovered blob (callers pass `--trees ${all_trees}`).

- [ ] **Step 1: Read the spec surfaces**

Run:
```bash
sed -n '1,120p' tests/test_hog_hierarchy_report.py
sed -n '1,60p' workflow/build_hog_report.py
git show 990531a~1:workflow/visualize_hog_hierarchy.py > /tmp/vhh_recovered.py; wc -l /tmp/vhh_recovered.py
grep -n 'og_short\|parse_newick\|build_data\|HTML_TEMPLATE\|--trees\|process_hg' tests/test_hog_hierarchy_report.py workflow/build_hog_report.py step4_ancestry.smk step4.ancestry.nf
```
Read all of `tests/test_hog_hierarchy_report.py` — it is the authoritative behavioral spec.

- [ ] **Step 2: Run the failing test to confirm the gap**

Run: `python -m pytest tests/test_hog_hierarchy_report.py -q`
Expected: FAIL/ERROR — `ImportError: cannot import name 'og_short_py'` or `ModuleNotFoundError: No module named 'visualize_hog_hierarchy'`.

- [ ] **Step 3: Reconstruct the module**

Start from `/tmp/vhh_recovered.py` as scaffold. Add/adjust to satisfy the contract above: implement `og_short_py` (preserving `:like:` suffix), `parse_newick_tree` (matching how `build_hog_report.py` calls it — read that file), add the `--trees` CLI argument and thread trees into `build_data`/`HTML_TEMPLATE` as the callers expect. Keep it a self-contained stdlib + pandas script (`#!/usr/bin/env python3`). Do NOT invent features the callers/tests don't use (YAGNI).

- [ ] **Step 4: Run the full hog-hierarchy test to green**

Run: `python -m pytest tests/test_hog_hierarchy_report.py -v`
Expected: PASS (all tests in the file).

- [ ] **Step 5: Confirm no regressions and caller references resolve**

Run:
```bash
python -m pytest tests/ -q 2>&1 | tail -5
python -c "import ast,sys; ast.parse(open('workflow/visualize_hog_hierarchy.py').read()); print('syntax ok')"
python workflow/visualize_hog_hierarchy.py --help
```
Expected: full suite green (Task-0 baseline pass count + the previously-erroring file now passing), `--help` lists `--links --stats --trees --levels --output`.

- [ ] **Step 6: Commit**

```bash
git add workflow/visualize_hog_hierarchy.py
git commit -m "fix(step4): restore missing visualize_hog_hierarchy.py to caller contract"
```

---

### Task 2: Consolidate step-4 config into `config/`

**Context:** `config/` holds `step1.yaml`, `step2.yaml` (Snakemake configs), but the step-4 config sits at root as `config_ancestry.yaml`, and `configs/` (plural) holds only the dead `configs/config.txt` (read solely by archived `workflow/old/*.sh`). `ancestry_ids.txt` is a step-4 input referenced by the step-4 config.

**Files:**
- Rename: `config_ancestry.yaml` → `config/step4.yaml`
- Rename: `ancestry_ids.txt` → `config/ancestry_ids.txt`
- Delete: `configs/config.txt` (and empty `configs/` dir)
- Modify: `config/step4.yaml` (`ids:` value), `step4_ancestry.smk:15,19`, `README.md` (L302, L323, L336, L349 area), `docs/manual.md` (refs to `configs/config.txt`)

- [ ] **Step 1: Move the files (preserve history)**

```bash
git mv config_ancestry.yaml config/step4.yaml
git mv ancestry_ids.txt config/ancestry_ids.txt
git rm configs/config.txt
```

- [ ] **Step 2: Update the moved config's internal ref**

In `config/step4.yaml`, change the `ids:` line from `ids:             ancestry_ids.txt` to `ids:             config/ancestry_ids.txt`. (Other paths in the file — `data/...`, `results/...` — are relative to repo-root CWD and stay unchanged.)

- [ ] **Step 3: Point the Snakemake file at the new config**

In `step4_ancestry.smk`:
- Line 15: `configfile: "config_ancestry.yaml"` → `configfile: "config/step4.yaml"`
- Line 19: default `config.get("ids", "ancestry_ids.txt")` → `config.get("ids", "config/ancestry_ids.txt")`
- Line 9 (comment/example in header): `ids=ancestry_ids.txt` → `ids=config/ancestry_ids.txt`

- [ ] **Step 4: Update docs references**

In `README.md`: replace `ancestry_ids.txt` with `config/ancestry_ids.txt` at the generation line (`... > ancestry_ids.txt`), the `--ids ancestry_ids.txt` nextflow example, and the `ids=ancestry_ids.txt` snakemake `--config` example. In `docs/manual.md`: remove or update any line referencing `configs/config.txt` (it is deleted).

- [ ] **Step 5: Verify config loads and no dangling refs**

```bash
python -c "import yaml; d=yaml.safe_load(open('config/step4.yaml')); print(d['ids']); assert d['ids']=='config/ancestry_ids.txt'"
grep -rn 'config_ancestry.yaml\|configs/config.txt' README.md docs/ instructions/ *.smk *.nf workflow/ && echo 'DANGLING REFS FOUND' || echo 'clean'
grep -rn '[^/]ancestry_ids.txt' README.md *.smk | grep -v 'config/ancestry_ids.txt' && echo 'CHECK root ancestry_ids refs' || echo 'ids refs ok'
python -m pytest tests/ -q 2>&1 | tail -3
```
Expected: yaml loads with correct `ids`, no dangling refs, tests still green.

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "refactor(config): move step4 config + ancestry_ids into config/, drop dead configs/config.txt"
```

---

### Task 3: Untrack build artifacts, clean junk, rewrite `.gitignore`

**Context:** `.snakemake/` (Snakemake state) and `downstream/*.png` (regenerable stats output) are force-tracked build artifacts. `img/phylo/` holds junk: `btag094.pdf` (stray 1.4 MB paper figure) and extensionless `Bolinf`/`Bolmic` (byte-identical download-error dupes next to real `.png`s). `.gitignore` has duplicate lines (`*.json`, `*.tab` twice; explicit `species_list` after `*_list`) and over-broad patterns (`*.tsv`, `*.json`, `species*`, `*html`, `slurm*`, `ids*`, `*_list`, `*.pdf`) that force legit inputs to be `-f`-added and silently drop new ones.

**Files:**
- Untrack: `.snakemake/` (all), `downstream/` (all)
- Delete: `img/phylo/btag094.pdf`, `img/phylo/Bolinf`, `img/phylo/Bolmic`
- Rewrite: `.gitignore`

- [ ] **Step 1: Confirm nothing references `downstream/`**

Run: `grep -rn 'downstream/' README.md docs/ instructions/ *.R *.nf *.smk workflow/ && echo 'REF FOUND — stop and report' || echo 'downstream/ unreferenced, safe to untrack'`
Expected: `downstream/` unreferenced. If a ref is found, stop and report (do not delete).

- [ ] **Step 2: Untrack artifacts and delete junk**

```bash
git rm -r --cached .snakemake
git rm -r downstream
git rm img/phylo/btag094.pdf img/phylo/Bolinf img/phylo/Bolmic
```

- [ ] **Step 3: Rewrite `.gitignore` cleanly**

Replace the entire file with (deduplicated; drops over-broad `*.tsv`/`*.json`/`species*`/`*html`/`slurm*`/`ids*`/`*_list`/`*.pdf`/`*.tab` so tracked inputs like `data/*.tsv`, `workflow/models/*.json`, `config/species_list`, `docs/*.html` stay tracked — outputs stay ignored via their directories):

```gitignore
# Nextflow / Snakemake run state
.nextflow*
work/
work_*/
.snakemake/

# Pipeline outputs
results/
downstream/
reports/
report-*.html
timeline*.html
trace*.txt
logs/
hmms/
pep2hg/
demo/

# Scratch / local inputs
tmp/
test/
data/input*

# Python
__pycache__/
*.pyc

# Editor / OS
*.log
```

- [ ] **Step 4: Verify legit tracked files are not newly ignored**

```bash
git check-ignore data/species_info.tsv workflow/models/models.json config/step4.yaml docs/step2_report.html && echo 'BAD: legit file ignored' || echo 'legit inputs safe'
git status --short
python -m pytest tests/ -q 2>&1 | tail -3
```
Expected: none of the legit files are ignored; `git status` shows the removals staged; tests green.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "chore: untrack .snakemake/ + downstream/, delete img junk, dedupe/tighten .gitignore"
```

---

### Task 4: Group R analysis scripts under `R/`, relocate the documented tool

**Context:** Seven loose R files sit at repo root with inconsistent casing. `train.R` is the one documented, argparse-driven tool (README L218, pairs with `workflow/predict_resources.py`) → move next to its Python sibling. The rest are analysis/scratch → group under `R/`. `_export_models.R` (superseded by `train.R`) and `workflow/predict_resources.R` (superseded by `workflow/predict_resources.py`, no caller, no test) are dead. `resources.R` and `generax_stats.R` contain `source('helper.R')`.

**Files:**
- Move: `resources.R`, `downstream_stats.R`, `generax_stats.R`, `helper.R` → `R/`; `runtime.r` → `R/runtime.R` (fix casing)
- Move: `train.R` → `workflow/train.R`
- Delete: `_export_models.R`, `workflow/predict_resources.R`
- Modify: `R/resources.R` + `R/generax_stats.R` (`source()` path), `README.md` (L218 `Rscript train.R`, L396 script-name mentions, L581-area if present), `docs/manual.md` (train.R path L296, script table L577-581)

- [ ] **Step 1: Create `R/` and move**

```bash
mkdir -p R
git mv resources.R R/resources.R
git mv downstream_stats.R R/downstream_stats.R
git mv generax_stats.R R/generax_stats.R
git mv helper.R R/helper.R
git mv runtime.r R/runtime.R
git mv train.R workflow/train.R
git rm _export_models.R workflow/predict_resources.R
```

- [ ] **Step 2: Fix internal `source()` paths (scripts run from repo root)**

In `R/resources.R` and `R/generax_stats.R`, change `source('helper.R')` → `source('R/helper.R')`.

- [ ] **Step 3: Update docs to the new paths**

In `README.md`: `Rscript train.R \` → `Rscript workflow/train.R \` (L218 area); the monitoring-section mentions of `downstream_stats.R`, `resources.R`, `generax_stats.R` → `R/downstream_stats.R`, `R/resources.R`, `R/generax_stats.R`. In `docs/manual.md`: update the `train.R` invocation path to `workflow/train.R` and the script-table entries for the moved/removed R scripts (remove `_export_models.R` and `workflow/predict_resources.R` rows).

- [ ] **Step 4: Verify no dangling refs**

```bash
grep -rn "source('helper.R')\|source(\"helper.R\")" R/ && echo 'BAD source path' || echo 'source paths ok'
grep -rn '[^/]train\.R\|_export_models\|predict_resources\.R\b' README.md docs/ instructions/ *.nf *.smk workflow/ R/ | grep -v 'workflow/train.R\|workflow/predict_resources.py' && echo 'CHECK refs' || echo 'refs clean'
python -m pytest tests/ -q 2>&1 | tail -3
```
Expected: `source('R/helper.R')` in both scripts, no dangling `train.R`/`_export_models`/`predict_resources.R` refs, tests green (predict_resources.py test unaffected).

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "refactor(R): group analysis scripts under R/, move train.R to workflow/, drop superseded R scripts"
```

---

### Task 5: Relocate loose root inputs and root doc

**Context:** `species_list`, `sps_annotate` (both analysis-set inputs) and `Gene_families.md` (documentation) sit at repo root. `species_list` is read by `config/step1.yaml` and referenced across README. `sps_annotate` has a duplicate entry (`Axidam` twice). `Gene_families.md` documents `data/gene_families_searchinfo.csv` fields → belongs in `docs/`. (Runs from Task 3's clean `.gitignore` so new paths aren't ignored.)

**Files:**
- Move: `species_list` → `config/species_list`; `sps_annotate` → `config/sps_annotate`; `Gene_families.md` → `docs/Gene_families.md`
- Modify: `config/step1.yaml:5` (`species_list:` value), `README.md` (species_list mentions L109/L116/L119/L126/L127, `Gene_families.md` link L108)

- [ ] **Step 1: Dedupe then move**

```bash
# remove the duplicate Axidam line in sps_annotate (keep first occurrence), preserving order
awk '!seen[$0]++' sps_annotate > sps_annotate.tmp && mv sps_annotate.tmp sps_annotate
git mv species_list config/species_list
git mv sps_annotate config/sps_annotate
git mv Gene_families.md docs/Gene_families.md
```

- [ ] **Step 2: Update `config/step1.yaml`**

Change `species_list:  species_list` → `species_list:  config/species_list`.

- [ ] **Step 3: Update README references**

In `README.md`: the input-table row and the `check_tree.py ... species_list ...` / `prepare_fasta.sh species_list ...` command lines → `config/species_list`. The `Gene_families.md` link in the input table → `docs/Gene_families.md`.

- [ ] **Step 4: Verify**

```bash
python -c "import yaml; d=yaml.safe_load(open('config/step1.yaml')); assert d['species_list']=='config/species_list', d['species_list']; print('ok')"
grep -rn '[^/]species_list\b' README.md config/ *.nf *.smk | grep -v 'config/species_list' && echo 'CHECK species_list refs' || echo 'species_list refs clean'
grep -rn 'Gene_families.md' README.md docs/ | grep -v 'docs/Gene_families.md' && echo 'CHECK Gene_families refs' || echo 'ok'
git grep -c '^Axidam' config/sps_annotate  # expect 1
python -m pytest tests/ -q 2>&1 | tail -3
```
Expected: yaml ok, no dangling refs, single `Axidam`, tests green.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "refactor: move species_list/sps_annotate to config/, Gene_families.md to docs/, dedupe sps_annotate"
```

---

### Task 6: Delete dead code and a mis-saved file

**Context:** `workflow/old/` (26 legacy bash/py files, zero live refs) and `workflow/s04_possvm.sh` (zero refs, oldest file) are dead. Four docs-only one-off scripts (`get_tips.py`, `strip_len.py`, `strip_support.py`, `prune_tree.py`) have no tests and no pipeline caller — only `docs/manual.md` mentions. `instructions/Design.md` is a mis-saved copy of the `frontend-design` Claude skill (nothing to do with phylohpc).

**Files:**
- Delete: `workflow/old/` (recursive), `workflow/s04_possvm.sh`, `workflow/get_tips.py`, `workflow/strip_len.py`, `workflow/strip_support.py`, `workflow/prune_tree.py`, `instructions/Design.md`
- Modify: `docs/manual.md` (remove the one-off-script usage lines)

- [ ] **Step 1: Re-confirm the one-offs are unreferenced by live code/tests**

```bash
for s in get_tips strip_len strip_support prune_tree s04_possvm; do
  echo "== $s =="; grep -rln "$s" tests/ *.nf *.smk workflow/*.py config/ 2>/dev/null | grep -v "workflow/$s"
done
grep -rln 'workflow/old/' *.nf *.smk workflow/*.py tests/ config/ README.md
```
Expected: no live (non-docs, non-self) references. If any test/pipeline reference appears, stop and report.

- [ ] **Step 2: Delete**

```bash
git rm -r workflow/old
git rm workflow/s04_possvm.sh workflow/get_tips.py workflow/strip_len.py workflow/strip_support.py workflow/prune_tree.py
git rm instructions/Design.md
```

- [ ] **Step 3: Scrub docs references**

In `docs/manual.md`, remove the usage lines/sections that invoke `get_tips.py`, `strip_len.py`, `strip_support.py`, `prune_tree.py` (they no longer exist). Leave surrounding valid content intact.

- [ ] **Step 4: Verify no dangling refs and tests green**

```bash
grep -rn 'workflow/old\|s04_possvm\|get_tips\|strip_len\|strip_support\|prune_tree\|instructions/Design' README.md docs/ instructions/ *.nf *.smk workflow/ config/ && echo 'DANGLING' || echo 'clean'
python -m pytest tests/ -q 2>&1 | tail -3
```
Expected: no dangling refs, tests green.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "chore: delete dead workflow/old + orphan one-off scripts + mis-saved instructions/Design.md"
```

---

## Post-Plan: Final review & integration

After all tasks: dispatch the whole-branch code review (`reorg-cleanup` vs `master`), address Critical/Important findings, then use `superpowers:finishing-a-development-branch` to merge `reorg-cleanup` into `master`. The handout branch (further stripping of `R/`, `TODOs.md`, `WHATIDID.md`, `instructions/`) is a follow-up cut from the reorganized `master`, not part of this plan.

## Self-Review Notes

- **Coverage:** Broken viz (T1), config sprawl (T2), artifact/gitignore hygiene (T3), R-script layout (T4), loose root inputs (T5), dead code (T6) — all exploration findings mapped. Snakemake retained (Global Constraints), phylogeny submodule untouched (Global Constraints).
- **Ordering dependency:** T2 before T5 (both touch config-relative inputs); T3 (clean `.gitignore`) before T5 (so `config/species_list` etc. aren't ignored). Sequential execution enforces this.
- **Test safety:** every task ends with `pytest tests/`; only truly test-free/caller-free files are deleted. `report_step1.py`, `visualize_ancestry.py`, `predict_resources.py`, `check_job*.py` (tested) are all kept.

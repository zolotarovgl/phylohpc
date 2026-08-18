# `dev/phylohpc` — the development checkout

**This is the clean clone. Pipeline changes are made HERE, committed and pushed; analysis
directories pull them.** Created 2026-08-18.

## Why it exists

Every existing checkout was an *analysis* directory with the pipeline inside it, so code,
config and outputs were entangled and no clone was clean enough to develop in:

| checkout | HEAD | state |
|---|---|---|
| `projects/2022_scRNAseq_Mlei/Phylogenies_v3` | `45b57fd` (= origin/master) | 5 modified config/data files; carries the Mlei run + `work/` + `results/` |
| `projects/2025_phylogeny/phylohpc` | `45b57fd` (= origin/master) | **9,293** dirty paths, plus `results/` and `work/` |
| `projects/2026_chromo/phylohpc` | `06097ad` on branch `claude/pipeline-cleanup-fixes-observability` | **9 commits of pipeline work that are NOT in master**, 51 dirty paths, plus `results/` |

That last row is the point: real pipeline development (report refactor, SLURM log routing,
step2 param fixes, failure diagnostics, 17 added chromatin families) is sitting on a branch
inside a chromatin analysis dir. It is pushed to origin — `origin/claude/pipeline-cleanup-fixes-observability`,
0 ahead / 0 behind — but not merged, so no other project sees it.

## Rules for this dir

- **No analysis runs here.** No `results/`, no `work/`, no project data. If it needs a
  `--outdir`, it is not a dev task.
- Branch from `master`, push, and merge on GitHub. Analysis dirs then `git pull`.
- The `phylogeny` submodule is pinned at `05d88cd` (a descendant of `origin/blastology`, not a
  tag); POSSVM at `1ee6113`. Bump deliberately, never as a side effect.

## First jobs, in order

1. **Merge `claude/pipeline-cleanup-fixes-observability` into master** — 9 commits, already
   reviewed enough to be pushed, currently invisible to every other project.
2. **The species-tree bug that killed the last GeneRax run** (2026-08-14, 475 tasks, exit 10).
   GeneRax's own log: `[Error] Cannot parse species tree file (species_tree.newick)`. The file
   staged into the task dir was `data/species_tree.newick` (md5 `95521216…`, **4 polytomies**),
   while `config/step2.yaml` asks for `data/species_tree.binary.newick` (md5 `4ae36437…`). So
   `params.SPECIES_TREE` and the YAML disagree, and GeneRax gets the tree its own docs say it
   cannot take. It is a path bug, not MPI: `Found generax` / `Found mpirun` both succeeded and
   `prterun` only relayed the child's exit code.
   → the fix is a `nextflow.config` / YAML precedence fix + a **fail-fast check** (run
   `workflow/check_tree.py` on `params.SPECIES_TREE` before the DAG, so this cannot recur).
3. **Stop `errorStrategy` from mislabelling failures.** `step2.nf` reports exit 10 as
   "Family parsing error — ignored", which sent the whole diagnosis toward the family files;
   the true message was in `.command.err` on `/no_backup` and went unread for four days.
   Print the first lines of the tool's log on failure instead of a guessed label.

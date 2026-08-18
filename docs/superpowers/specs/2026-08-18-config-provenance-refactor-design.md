# phylohpc: config, provenance and failure-reporting refactor

**Status:** design approved 2026-08-18 (Grygoriy). Branch `refactor/config-provenance`.
**Scope:** the configuration / provenance / reporting layer. **The DAG and the processes are
not touched.**

## 1. The problem, with evidence

On 2026-08-14 a `step2` run in `Phylogenies_v3` lost **all 475 GeneRax tasks**. Every one
exited 10 within ~65 s. GeneRax's own log:

```
[00:00:00] Filtering invalid families...
[Error] Cannot parse species tree file (species_tree.newick)
```

The staged tree was `data/species_tree.newick` (md5 `95521216`, 53 internal nodes, **4
multifurcating**) while `config/step2.yaml` asked for `data/species_tree.binary.newick`
(md5 `4ae36437`). GeneRax requires a strictly binary tree. Four compounding causes:

1. **The GeneRax channel read `params.SPECIES_TREE`** (uppercase, set only in
   `nextflow.config`), while the documented switch is `--species_tree` (lowercase). The
   documented flag could not affect GeneRax at all.
2. **`config/step2.yaml` was never read by the Nextflow path.** It is the *Snakemake*
   config (`docs/manual.md:81`), there is no `-params-file` anywhere in the repo, and the
   run was `nextflow run … -profile slurm,precise`. Every edit to it configured nothing.
3. **The failure was mislabelled.** `errorStrategy` maps exit 10 to
   `"Family parsing error — ignored"` — a guess. It sent the diagnosis at the family files
   for four days; the real message was one line in `.command.err` on `/no_backup`.
4. **The run reported success** with `results/generax/` empty, and the bad state propagated
   into `gene_lists/` and from there into Fig 3a panel labels.

None of this is a DAG or a process bug. It is entirely the configuration and reporting
layer, which is why that is the whole scope here.

Two more measurements that shape the design:

- **Three params exist in two casings**: `species_tree`/`SPECIES_TREE`,
  `outdir`/`OUTDIR`, `refnames`/`REFNAMES`. 41 params in `nextflow.config`, 29 distinct
  `params.*` referenced across the `.nf` files.
- **Two implementations exist** for the same pipeline: `step1.nf`/`step2.nf` (694 lines)
  and `step1.smk`/`step2.smk` (327 lines), reading *different* config sources.

## 2. Goals / non-goals

**Goals**
- Exactly one place a run's science parameters live, and it is read.
- A run is reproducible from `run.yaml` + a pipeline revision.
- A failure names itself, is counted, and cannot be mistaken for success.
- One entry point whose pre-flight checks cannot be bypassed.

**Non-goals (deliberately deferred)**
- Splitting `step2.nf` (562 lines) into per-stage modules. Next pass.
- Renaming `step1`/`step2`/`step4` to `search`/`phylogeny`/`ancestral`.
- JSON-schema validation of `run.yaml` (`nf-schema`). Worth doing once the parameter set
  has been consolidated, not before.
- Fixing GeneRax's MPI build. Separate issue; the species-tree bug is not MPI.

## 3. Design

### 3.1 Parameter model

- **`nextflow.config` keeps infrastructure only**: executors, queues, profiles, retry
  policy, per-family resource models. No science parameters.
- **`run.yaml` holds every science parameter** and is passed with `-params-file`.
- **One canonical name per parameter**, lower_snake_case. The uppercase twins are removed;
  during a deprecation window `step2.nf` accepts an uppercase key only when the lowercase
  one is absent, and logs a warning naming the file to fix.
- Nothing may read a parameter that is not in the canonical list.

### 3.2 Provenance

The pipeline is run **by revision from GitHub**:

```bash
nextflow run zolotarovgl/phylohpc -r v1.0 -params-file run.yaml -profile slurm
```

Nextflow caches the repo and stamps the revision into the run report, so
`run.yaml` + revision is the complete record. Analysis directories no longer clone the
pipeline; development happens in `dev/phylohpc` and is tagged.

Every run prints a **resolved-configuration block** before the DAG is built: params file,
pipeline revision, output dir, work dir, species tree (with binary/not), family count,
`run_generax` on/off. This block is the direct antidote to "which file configured this?".

### 3.3 Entry point

One script, `phylohpc`, with subcommands sharing a single flag parser and one pre-flight:

| command | does |
|---|---|
| `phylohpc run <step> -p run.yaml` | foreground, for an interactive `salloc`/`srun` session |
| `phylohpc submit <step> -p run.yaml` | `sbatch`, same checks, same logging |
| `phylohpc rerun-failed <step> -p run.yaml` | re-runs only the failed families (§3.5) |
| `phylohpc report` | replaces `make_report.sh` |
| `phylohpc init` | writes a commented `run.yaml` skeleton — convenience only |

`submit_nf.sh` and `run_interactive.sh` are absorbed into it.

### 3.4 Pre-flight

Runs in the wrapper (fails in ~1 s) **and** inside the workflow (so it applies however the
pipeline is launched):

- every input path in `run.yaml` exists;
- when `run_generax` is on: the species tree parses and is **strictly binary**, reported as
  a count of multifurcating nodes plus the `check_tree.py --random-resolve` command that
  fixes it; `generax` and `mpirun` are on `PATH` (warn, since the bioconda build is
  non-MPI);
- the family table is readable and non-empty.

*Already implemented on this branch* (`af89c31`) for the species tree; verified against the
exact tree that killed the 2026-08-14 run — 4 multifurcating nodes detected.

### 3.5 Failure reporting and rerun

- **No invented labels.** On task failure the pipeline prints the first lines of that
  task's own stderr / tool log.
- **`reports/failures.tsv`**: `family, exit_code, first_error, workdir`.
- **Tolerance**: the run FAILS if failures exceed **5 % of families or 20 families,
  whichever is smaller**. `--strict` fails on the first failure. Below tolerance the run
  ends `COMPLETE-WITH-FAILURES`, never a silent success.
- **An expected output directory that ends up empty is an error**, independent of exit
  codes — the specific hole through which the 2026-08-14 run reported success.
- **`rerun-failed`** is thin because the hook exists: `step2.nf` already accepts
  `params.ids` (one HG id per line). It reads `reports/failures.tsv`, writes
  `reports/rerun.ids`, and relaunches with `--ids reports/rerun.ids -resume`; successful
  families stay cached.

### 3.6 Removals

- `step1.smk`, `step2.smk`, `workflow/step4_ancestry.smk` — the Snakemake implementation.
  Nextflow is the supported path (decided 2026-08-18); two implementations reading
  different configs is what let a config file become decorative.
- `submit_nf.sh`, `run_interactive.sh`, `make_report.sh` — absorbed into `phylohpc`.

## 4. Migration for existing analysis directories

`Phylogenies_v3`, `2025_phylogeny`, `2026_chromo` each hold a clone plus local edits.
Per directory: write `run.yaml` from what that project actually ran (recoverable from
`.nextflow.log`), delete the clone's pipeline files, keep `data/`, `results/`, `run.yaml`.
⚠ `2026_chromo` also holds the 9 cleanup commits — those are already merged on this
branch, so nothing is lost when its clone goes.

## 5. Verification

The refactor is config-layer only, so the test is that **the same inputs give the same
outputs**:

1. Re-run `step2` for `Phylogenies_v3` on this branch with `run_generax: false` and
   compare POSSVM output against the existing `results/possvm_prev/` — must be identical.
2. Re-run with `run_generax: true` and a **binary** species tree; the previous run produced
   0 trees, so any success is progress, and the count of families completed is the metric.
3. Negative tests, each of which must fail fast with a named cause: multifurcating species
   tree; missing family table; a typo'd key in `run.yaml`; `generax` absent from `PATH`.
4. `rerun-failed` on a run with a deliberately broken family: exactly that family re-runs.

## 6. Risks

- **`-r <revision>` pins a run to a tag**; fixes do not reach it until re-tagged and re-run.
  Accepted: it is the price of the run being reproducible at all.
- **Deleting the `.smk` files** removes the local/small-run path some projects may use.
  Mitigated by `phylohpc run`, which runs the same DAG locally.
- **A tolerance threshold can hide a systematic failure** below the cut. Mitigated by the
  empty-output check and by `failures.tsv` always being written.
- **The merge on this branch is untested.** Nothing has been run since; a `-stub` or dry
  run on the cluster is the first task in the plan.

# PhyloHPC — Pipeline Improvement Plan (v2, post-review)

Status: revised after a three-angle critique (scientific/reproducibility,
workflow-engineering, pragmatism). Changes from v1 are listed at the bottom.
Scope: the orchestration layer (`step1/2/4` `.nf`/`.smk`, `nextflow.config`,
`config/*.yaml`, `workflow/` helpers). The `phylogeny` submodule and the
scientific methods are fixed inputs except where their *invocation* lives in this
layer.

## Goals
1. Protect result integrity (no silent data loss into ancestral reconstruction).
2. Make runs reproducible and portable off the CRG HPC.
3. Stop the two engines diverging — preferably by not maintaining two.
4. Reduce maintenance burden and onboarding friction.

## Guiding principles
- Behaviour-preserving refactors first; change scientific outputs only deliberately and documented.
- Fail loud, not silent. Validate positive outputs, not just exit codes.
- One source of truth for any value both engines consume.
- Every change ships with a test or an explicit "why untestable".

---

## Phase 0 — Decide the engine strategy (DO FIRST)
The single highest-leverage decision, not an afterthought. Today the engines
already disagree (`max_n` 2000 vs 3000; step4 tree source + POSSVM flags; resource
model exists only in `.nf`). Maintaining both means doing W1/W3/W5 twice and
building W6 to police the duplication.

**Recommendation: make Nextflow primary; freeze or delete the Snakemake port.**
README and all real runs are Nextflow; SLURM, the resource-prediction model, and
GeneRax are Nextflow-native. If a collaborator genuinely needs Snakemake, keep it
but drop the cross-engine conformance test (W6b). This decision gates W5/W6 effort
(roughly halves it) and converts much of W8 into a deletion win.

**Acceptance.** A written decision in this doc; if "NF primary", `.smk` files moved
to `archive/` or clearly marked unsupported.

---

## Phase A — Protect results

### W1 — Failure manifest + positive output validation (HIGH)
**Problem.** `errorStrategy 'ignore'` on `ALN`/`PHY`/`PVM`/`GR_watcher` (`step2.nf:107,162,237,365`)
drops failed HGs with no record; you cannot tell "OG absent because the family
failed" from "absent biologically" — this threatens ancestral-content conclusions.
Crash-logging alone is insufficient: a POSSVM run can exit 0 yet emit an
empty/degenerate `ortholog_groups.csv`, and GeneRax "Exit 10" is ignore-by-design.
**Change.**
- Capture failures reliably: a task killed by `errorStrategy 'ignore'` emits **no**
  output, so do NOT rely on `collectFile` of a failed task's output. Either (a) wrap
  the tool in `set +e`, decide success in-script, and always exit 0 writing either
  the real output or a declared **optional** `${id}.failed` marker; or (b) harvest
  failures in `workflow.onComplete{}` from the trace records. Prefer (a) — it keeps
  the signal in the dataflow.
- `collectFile` markers → `results/failures.tsv` (id, step, exit, stderr tail),
  distinguishing **ignored-by-design** (e.g. GeneRax exit 10) from genuine failure.
- Positive validation pass: assert every retained HG produced a **non-empty** OG CSV
  and that the PAM column count reconciles against the POSSVM outputs.
- `workflow.onComplete`: print `N/M families failed at <step>`; exit non-zero if the
  failure rate exceeds `params.max_fail_frac` (default 0.1).
**Effort.** ~1 day. **Risk.** Low.
**Acceptance.** A deliberately-broken HG and an empty-OG HG both surface in
`failures.tsv`; the run prints counts; exceeding the threshold exits non-zero.

### W2 — Trustworthy resume (HIGH, behaviour-changing)
**Problem.** `ALN`/`PHY`/`GR_watcher` hand-roll caching by checking
`file("${OUTDIR}/.../${id}...").exists()` and symlinking (`step2.nf:120-124,178-187,403-413`),
plus a **third** cache: the IQ-TREE `.ckp.gz` checkpoint branch (`step2.nf:189-203`).
This fights `-resume`, copies files onto themselves via `publishDir mode:'copy'`,
has no staleness check, and records fake provenance.
**Change.** Replace with **`storeDir`** (the correct tool: it persists outside
`work/`, surviving `rm -rf work` on scratch — which bare `-resume` does not). Remove
**all three** hand-rolled caches, including the ckp branch. Document the `storeDir`
staleness trade-off (keys on output filename; a changed input under the same `${id}`
reuses the stale result). Note `cache 'lenient'` is currently set only in the
`slurm` profile (`nextflow.config:19`), not `local`.
**Effort.** ~1.5 days (incl. a real cluster `-resume`/`storeDir` test). **Risk.**
Medium-High — could silently recompute a large run if mis-set. Do not ship without
a cluster test.
**Acceptance.** Re-running after `rm -rf work` reuses stored outputs; no symlink
self-copies; `trace` shows real execution; a changed input invalidates the entry.

### W6a — `nextflow -preview` smoke test in CI (HIGH, cheap)
**Problem.** Nothing parses/executes the `.nf`; `test_step2_smk.py` even reimplements
smk logic in Python. That is how the `--logfile` bug and the null
`file(params.SPECIES_TREE)` crash could ship.
**Change.** Add `nextflow run <step>.nf -preview -profile test,local` (builds the
DAG, validates args/params, no execution — needs only Java). Wire into
`.github/workflows`. (`-stub` is NOT usable: no process defines a `stub:` block.)
**Effort.** ~0.5 day. **Risk.** Low.
**Acceptance.** CI fails on an invalid `.nf` arg or a null required param.

---

## Phase B — Correctness & reproducibility

### W-step4 — Reconcile step4 (HIGH correctness; promoted out of v1's W5)
**Problem.** Not param drift — a method divergence. NF (`step4.ancestry.nf:181`)
consumes raw `${gene_trees_dir}/${id}.treefile` (default `results/gene_trees`) via
the `main.py possvm` wrapper. SMK (`step4_ancestry.smk:25,70,82-91`) consumes
`${id}.generax.tree` from `results/possvm` via `possvm.py` directly with
`-min_support_transfer 0 -method lpa`. Different input trees AND different POSSVM
parameters → different orthogroups, and the SMK path contradicts the README's
"raw gene trees … to avoid circularity" (`README.md:308`).
**Change.** Decide the intended method (raw trees per the README), make both engines
match exactly (or, post-Phase-0, keep only NF). Assert "step4 consumes raw gene
trees" as a tested invariant.
**Effort.** ~0.5 day. **Risk.** Medium (changes step4 outputs deliberately).
**Acceptance.** Both engines (or the surviving one) read raw trees with identical
POSSVM flags; a test pins the tree-source path.

### W3 — Reproducibility, descoped to the high-value half (HIGH)
**Problem.** IQ-TREE2 runs `-bb` (UFBoot) with **no `-seed`** *and* `-nt AUTO`
(`phylogeny/helper/functions.py:134`); FastTree (the default `fast` profile, which
feeds step4) is unseeded (`functions.py:170`); MAFFT/DIAMOND/MCL are thread-sensitive,
so even **HG identity** — the unit of everything downstream — is not reproducible.
Only the submodule commit is pinned at runtime.
**Change (cheap, do now):**
- Add `params.seed` (default 42); thread `-seed` AND a **fixed `-nt N`** (not AUTO)
  into the IQ-TREE call. Requires a small `phylogeny` submodule change + re-pin.
- Seed the FastTree path or explicitly document it as non-reproducible.
- Record tool + DB versions to `results/versions.yml`; **checksum `Pfam-A.hmm`** and
  record its release; version `genefam.csv` (its threshold column changes results).
- Assert the submodule commit at startup.
**Change (optional, separate effort — NOT 1.5 days):** a Singularity/Docker profile.
This is a *recurring* tax (image rebuilds on every tool bump; Singularity + OpenMPI +
GeneRax on SLURM is fiddly). `module load` + the pinned `environment.yaml` is
"good enough" reproducibility for one HPC; treat the container as a later, optional
workstream.
**Effort.** ~1 day for the cheap half. **Risk.** Medium (submodule change).
**Acceptance.** Same-input, **same-thread-count** runs reproduce HG membership,
alignments, trees, and POSSVM OG CSVs — state the thread-count caveat honestly.

### W-possvm — POSSVM invocation correctness (MEDIUM)
Even with the POSSVM *algorithm* out of scope, its invocation lives here:
- The wrapper silently rescales `min_support_transfer` (÷100 if it exceeds the max
  observed support; `functions.py:202-205`), so FastTree (supports 0–1) and IQ-TREE
  (0–100) trees are thresholded on different effective criteria. **Log when this
  fires**; have W1's manifest flag it.
- README advertises clade-scoping via `--ignoretips` (`README.md:29,308`) but step4
  passes `--outgroup` (`step4.ancestry.nf:95`), which only roots; OGs are filtered
  post-hoc in `link_hog_levels.py`. Reconcile README↔code and decide whether OGs must
  be scoped at definition time.

### W4 — Portability (MEDIUM)
**Problem.** Hardcoded CRG paths/queues: `pfam_db` (`step1.nf:11`, `nextflow.config:75`),
`prepare_db_dir` (`config/step1.yaml`), and `submit_nf.sh:4-5`'s SBATCH
`-p genoa64 --qos=pipelines`, README workdirs under `/no_backup/asebe/...`.
**Change.** Required params with no default; a `crg` profile for the site values;
parametrise the SBATCH header in `submit_nf.sh`. Document required overrides.
**Effort.** ~0.5 day. **Risk.** Low.

### W5 — Param hygiene & shared scientific config (MEDIUM)
**Problem.** Mixed casing + the fallback ladder (`step2.nf:4-23`); `step2.nf` eagerly
calls `file(params.SPECIES_TREE)` at parse time (`:477-478`) → null-crash without
`-profile`; `step1.nf:25-34` doesn't skip `#` lines in `genefam.csv` (`.smk` does);
`max_n` 2000 vs 3000.
**Change.** One casing convention; real in-file defaults in `step2.nf`; skip comments
in genefam parsing. For shared values, note the engines **cannot** read one native
file: use a single **YAML** as source of truth, loaded by Nextflow via
`-params-file` and Snakemake via `configfile`, scoped to **scientific scalars only**
(`max_n`, `node_names`, POSSVM flags, tree source). The Groovy-parsed `resources.tsv`
model (`step2.nf:38-79`) has no Snakemake equivalent — another reason Phase 0 matters.
**Effort.** ~1 day. **Risk.** Medium — do after W6a so the smoke test catches regressions.

---

## Phase C — Lock it in

### W6b — Cross-engine conformance (MEDIUM; only if keeping both engines)
Run `.nf` and `.smk` on a tiny fixture and diff OG CSVs — but **canonicalize first**
(sort by gene id; normalize arbitrary OG relabeling) or the diff flaps. Hidden cost
the v1 understated: the "3-HG fixture" runs MAFFT→tree→POSSVM(→GeneRax) and needs
real tools + curated FASTA/tree data + wall-time; GitHub runners have no SLURM, so it
runs `local`-executor only, likely nightly/self-hosted. **If Phase 0 drops Snakemake,
cut this entirely.**

---

## Not recommended

### W7 — nf-core modularization (CUT)
A 3–5 day (realistically 2–3×) rewrite to `modules/`/`subworkflows/`/`nf-validation`
serves future contributors who don't exist for a one-scientist repo, risks
destabilizing a pipeline that produces publication results, and **doubles** the
two-engine tax. Pursue modularization only opportunistically, only post-Phase-0, and
only single-engine.

### W8 — Remaining cleanup (LOW, anytime)
Merge identical `PVM`/`PVM_PREV`; remove `step4_ancestry.smk`'s dead second
implementation (`rule visualize_hierarchy`) — but confirm `build_report` subsumes the
`og_links/og_stats` artifacts first; make `REPORT`'s `species_info.tsv` a param.
(`generax.nf`, `workflow/old/`, `predict_resources.R` already archived.)

---

## Sequencing
0. **Phase 0 decision** (engine strategy).
1. **Phase A:** W1 (corrected mechanism) → W6a → W2.
2. **Phase B:** W-step4 → W3 (cheap half) → W-possvm → W4 → W5.
3. **Phase C:** W6b (only if dual-engine). W8 anytime. W7 not recommended.

## Open decisions
1. Engine strategy (Phase 0) — **resolve first**.
2. Container tech (Singularity vs Docker) — only if/when the optional container half of W3 is funded.
3. The 4 `diamond`-strategy chromatin families (PCRing/Kelch/WD40/zf-C2H2): support diamond search, add as `GA`, or leave out (a real methodological change needing a versioning story).

## Out of scope
The scientific methods themselves (POSSVM/GeneRax/PastML algorithms, F81, MPPA) —
but NOT their invocation, which W-possvm and W3 own.

---

## Changes from v1 (after subagent review)
- Added **Phase 0** (engine strategy) as the first, highest-leverage decision.
- **W1**: corrected the capture mechanism (a `errorStrategy 'ignore'` task emits no
  output; use `set +e`+optional marker or `onComplete` trace) and added positive
  output validation + ignored-by-design vs failed distinction.
- Promoted the **step4 divergence** from W5 to its own HIGH correctness workstream.
- **W3**: split into a cheap high-value half (seed + fixed `-nt` + versions/checksums
  + submodule assertion) and an optional container; added FastTree/MAFFT/DIAMOND/MCL
  and HG-identity determinism; Pfam/genefam versioning.
- **W6**: split into W6a (cheap `-preview` smoke test, promoted to Phase A) and W6b
  (expensive conformance, conditional on keeping both engines); `-stub` not usable.
- **W2**: target `storeDir` (not bare `-resume`); remove the ckp branch too; flag the
  profile-specific `cache 'lenient'` and the cluster-test requirement.
- Added **W-possvm** (support-threshold rescaling logging; `--ignoretips` vs
  `--outgroup` README↔code reconciliation).
- **W4**: added `submit_nf.sh` SBATCH header.
- **W5**: clarified shared config via `-params-file`/`configfile`, scoped to scalars.
- **W7**: cut (was LOW/optional).

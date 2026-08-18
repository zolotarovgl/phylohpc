# phylohpc config/provenance/failure-reporting refactor — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make a phylohpc run configured from exactly one file that is actually read, reproducible from that file plus a pipeline revision, and impossible to mistake a mass failure for success.

**Architecture:** `nextflow.config` keeps infrastructure only; all science parameters move to a `run.yaml` passed with `-params-file`, under one canonical lower_snake_case name each. A single `phylohpc` wrapper provides `run` / `submit` / `rerun-failed` / `report` / `init`, sharing one pre-flight that also exists inside the workflow. Task failures print their own stderr, are collected into `reports/failures.tsv`, and fail the run past a tolerance. The DAG and the processes are not touched.

**Tech Stack:** Nextflow DSL2 (Groovy), bash, Python 3 (validators), SLURM. No new runtime dependencies.

**Spec:** `docs/superpowers/specs/2026-08-18-config-provenance-refactor-design.md`

## Global Constraints

- Branch: **`refactor/config-provenance`**. `master` stays a clean mirror of `origin/master`; never commit there.
- **Do not modify** the DAG shape, process bodies, or the science of any step. Config/report layer only.
- Canonical parameter names are **lower_snake_case**. Uppercase twins (`SPECIES_TREE`, `OUTDIR`, `REFNAMES`, `IQTREE2_MODEL`, `REFSPECIES`, `SUBS_MODEL`, `MAX_SPR`, `NCPU_GENERAX`, `TREE_METHOD`, `MAFFT_OPT`) are accepted for one deprecation window **only when the lowercase key is absent**, and must log a warning naming the file to fix.
- The repo has **no test framework** (the pytest suite was deleted by the merged cleanup branch). Tests in this plan are bash scripts under `tests/`, run by `tests/run_tests.sh`, exiting non-zero on failure. Do not reintroduce pytest.
- Test fixtures use the committed 7-species dataset (`data/input.fasta`, species Aque Cowc Dmel Mlei Mmus Nvec Tadh) and small hand-written newick files. **No test may require SLURM or GeneRax.**
- GeneRax requires a **strictly binary** species tree. Never "fix" a multifurcating tree silently; fail with the `workflow/check_tree.py --random-resolve` instruction.
- Every commit message states what was measured, not what was intended.

---

### Task 1: Test harness

**Files:**
- Create: `tests/run_tests.sh`
- Create: `tests/lib.sh`
- Create: `tests/data/binary.newick`
- Create: `tests/data/multifurcating.newick`
- Create: `tests/test_fixtures.sh`

**Interfaces:**
- Consumes: nothing.
- Produces: `tests/lib.sh` exporting `assert_eq <expected> <actual> <msg>`, `assert_contains <haystack> <needle> <msg>`, `assert_fails <cmd...>`, `assert_ok <cmd...>`; `tests/run_tests.sh` running every `tests/test_*.sh` and reporting a pass/fail count. Later tasks add `tests/test_<topic>.sh` files that source `tests/lib.sh`.

- [ ] **Step 1: Write the failing test**

`tests/test_fixtures.sh`:
```bash
#!/usr/bin/env bash
set -uo pipefail
source "$(dirname "$0")/lib.sh"

# the two fixture trees must exist and differ in exactly the property we test on
assert_ok  test -s "$(dirname "$0")/data/binary.newick"
assert_ok  test -s "$(dirname "$0")/data/multifurcating.newick"
assert_eq  "0" "$(python3 "$(dirname "$0")/../workflow/count_multifurcations.py" "$(dirname "$0")/data/binary.newick")" \
           "binary fixture has no multifurcations"
assert_eq  "1" "$(python3 "$(dirname "$0")/../workflow/count_multifurcations.py" "$(dirname "$0")/data/multifurcating.newick")" \
           "multifurcating fixture has exactly one"
```

- [ ] **Step 2: Run it to verify it fails**

Run: `bash tests/test_fixtures.sh`
Expected: FAIL — `tests/lib.sh` does not exist yet.

- [ ] **Step 3: Write the harness and fixtures**

`tests/lib.sh`:
```bash
#!/usr/bin/env bash
# Minimal assertions. Every test script sources this and exits non-zero on failure.
_FAILED=0
assert_eq()       { if [[ "$1" != "$2" ]]; then echo "  FAIL: ${3:-assert_eq}: expected '$1', got '$2'"; _FAILED=1; else echo "  ok: ${3:-assert_eq}"; fi; }
assert_contains() { if [[ "$1" != *"$2"* ]]; then echo "  FAIL: ${3:-assert_contains}: '$2' not found"; _FAILED=1; else echo "  ok: ${3:-assert_contains}"; fi; }
assert_ok()       { if "$@" >/dev/null 2>&1; then echo "  ok: $*"; else echo "  FAIL: expected success: $*"; _FAILED=1; fi; }
assert_fails()    { if "$@" >/dev/null 2>&1; then echo "  FAIL: expected failure: $*"; _FAILED=1; else echo "  ok (failed as expected): $*"; fi; }
trap '[[ $_FAILED -eq 0 ]] || exit 1' EXIT
```

`tests/data/binary.newick`:
```
((A:1,B:1):1,(C:1,D:1):1);
```

`tests/data/multifurcating.newick`:
```
((A:1,B:1,C:1):1,D:1);
```

`tests/run_tests.sh`:
```bash
#!/usr/bin/env bash
# Run every tests/test_*.sh. Exit non-zero if any fails.
cd "$(dirname "$0")/.."
fail=0
for t in tests/test_*.sh; do
  echo "== $t"
  bash "$t" || fail=1
done
[[ $fail -eq 0 ]] && echo "ALL TESTS PASSED" || echo "TESTS FAILED"
exit $fail
```

`workflow/count_multifurcations.py`:
```python
#!/usr/bin/env python3
"""Print the number of multifurcating (>2 children) nodes in a newick file.

Used by the pre-flight in both the wrapper and step2.nf: GeneRax requires a
strictly binary species tree and its own error message ("Cannot parse species
tree file") does not say why.
"""
import sys

def count_multifurcations(newick: str) -> int:
    kids, bad = [], 0
    for ch in newick:
        if ch == '(':
            kids.append(1)
        elif ch == ',' and kids:
            kids[-1] += 1
        elif ch == ')':
            n = kids.pop() if kids else 0
            if n > 2:
                bad += 1
    return bad

if __name__ == '__main__':
    print(count_multifurcations(open(sys.argv[1]).read().strip()))
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `chmod +x tests/run_tests.sh && bash tests/run_tests.sh`
Expected: `ALL TESTS PASSED`.

- [ ] **Step 5: Verify against the real failure**

Run: `python3 workflow/count_multifurcations.py ~/nobackup/nextflow/phylohpc/work_step2/0c/d870ef*/species_tree.newick`
Expected: `4` — the tree that killed the 2026-08-14 run. If this prints 0, the parser is wrong; stop and fix it before continuing.

- [ ] **Step 6: Commit**

```bash
git add tests/ workflow/count_multifurcations.py
git commit -m "test: bash harness plus a multifurcation counter verified on the tree that killed the 2026-08-14 run (4 nodes)"
```

---

### Task 2: Canonical parameter names

**Files:**
- Modify: `step2.nf:1-30` (the params resolution block), `step2.nf:529-560` (species tree block from Task 0's fix)
- Modify: `step1.nf:1-20`
- Modify: `nextflow.config` (params + both profile blocks)
- Create: `workflow/params.nf`
- Create: `tests/test_params.sh`

**Interfaces:**
- Consumes: `tests/lib.sh` from Task 1.
- Produces: `workflow/params.nf` defining `def canonical(name, upper)` — returns `params[name]` when present, else `params[upper]` with a deprecation warning, else the default. Later tasks call `canonical('species_tree','SPECIES_TREE')` etc.

- [ ] **Step 1: Write the failing test**

`tests/test_params.sh`:
```bash
#!/usr/bin/env bash
set -uo pipefail
source "$(dirname "$0")/lib.sh"
cd "$(dirname "$0")/.."

# no uppercase param may be READ anywhere except through the deprecation shim
hits=$(grep -nE 'params\.(SPECIES_TREE|OUTDIR|REFNAMES|REFSPECIES|SUBS_MODEL|MAX_SPR|NCPU_GENERAX|TREE_METHOD|MAFFT_OPT|IQTREE2_MODEL)' step1.nf step2.nf | grep -v 'workflow/params.nf' | wc -l)
assert_eq "0" "$hits" "no direct uppercase params.* reads outside the shim"

# every param the .nf files read must exist in nextflow.config or be defined with a default
assert_ok grep -q 'species_tree' nextflow.config
assert_ok grep -q 'def canonical' workflow/params.nf
```

- [ ] **Step 2: Run it to verify it fails**

Run: `bash tests/test_params.sh`
Expected: FAIL — `workflow/params.nf` missing and uppercase reads still present.

- [ ] **Step 3: Write the shim**

`workflow/params.nf`:
```groovy
// One canonical name per parameter, lower_snake_case.
//
// WHY THIS EXISTS: params.SPECIES_TREE (uppercase, nextflow.config only) and
// params.species_tree (lowercase, the documented --species_tree flag) both existed and
// the GeneRax channel read the uppercase one. The documented flag therefore could not
// affect GeneRax, and on 2026-08-14 that handed it a multifurcating tree and killed all
// 475 tasks. Read every parameter through canonical() so a name can never fork again.
def canonical(String name, String upper, Object fallback = null) {
    if( params.containsKey(name) && params[name] != null )
        return params[name]
    if( params.containsKey(upper) && params[upper] != null ) {
        log.warn "DEPRECATED parameter '${upper}' — rename it to '${name}' in your params file / nextflow.config"
        return params[upper]
    }
    return fallback
}
```

- [ ] **Step 4: Replace every uppercase read**

In `step2.nf`, add `includeConfig`-style inclusion at the top (`evaluate(new File("${projectDir}/workflow/params.nf"))` is not valid DSL2; instead paste `canonical()` into `step2.nf` directly above the params block and keep `workflow/params.nf` as the single source copied by both step files — note this duplication in a comment). Then replace, one at a time:

| was | becomes |
|---|---|
| `params.SPECIES_TREE` | `canonical('species_tree','SPECIES_TREE',"${projectDir}/data/species_tree.newick")` |
| `params.OUTDIR` | `canonical('outdir','OUTDIR',"${projectDir}/results")` |
| `params.REFNAMES` | `canonical('refnames','REFNAMES',"${projectDir}/data/Mmus_gene_names.csv")` |
| `params.REFSPECIES` | `canonical('refspecies','REFSPECIES','Mmus')` |
| `params.SUBS_MODEL` | `canonical('subs_model','SUBS_MODEL','LG')` |
| `params.MAX_SPR` | `canonical('max_spr','MAX_SPR',3)` |
| `params.NCPU_GENERAX` | `canonical('ncpu_generax','NCPU_GENERAX',2)` |
| `params.TREE_METHOD` | `canonical('tree_method','TREE_METHOD','iqtree2')` |
| `params.MAFFT_OPT` | `canonical('mafft_opt','MAFFT_OPT','')` |
| `params.IQTREE2_MODEL` | `canonical('iqtree2_model','IQTREE2_MODEL','TEST')` |

In `nextflow.config`, rename the keys in the `params` block and in both `fast` and `precise` profiles to the lowercase names.

- [ ] **Step 5: Run tests**

Run: `bash tests/run_tests.sh`
Expected: `ALL TESTS PASSED`.

- [ ] **Step 6: Commit**

```bash
git add step1.nf step2.nf nextflow.config workflow/params.nf tests/test_params.sh
git commit -m "refactor(params): one canonical lower_snake_case name per parameter, uppercase accepted only via a warning shim"
```

---

### Task 3: run.yaml is the only science config, and the run says what it resolved

**Files:**
- Create: `config/run.example.yaml`
- Modify: `step2.nf` (add the resolved-configuration banner after the params block)
- Modify: `nextflow.config` (strip science params, keep infrastructure)
- Create: `tests/test_banner.sh`

**Interfaces:**
- Consumes: `canonical()` from Task 2.
- Produces: a `printResolvedConfig(step)` function in `step2.nf` printing one `key: value` per line under a `== resolved configuration` header; `config/run.example.yaml` documenting every key.

- [ ] **Step 1: Write the failing test**

`tests/test_banner.sh`:
```bash
#!/usr/bin/env bash
set -uo pipefail
source "$(dirname "$0")/lib.sh"
cd "$(dirname "$0")/.."

assert_ok test -f config/run.example.yaml
# the example must mention every canonical key the workflow reads
for k in species_tree outdir refnames refspecies subs_model max_spr ncpu_generax tree_method mafft_opt iqtree2_model run_generax ids family_info; do
  assert_ok grep -q "^${k}:" config/run.example.yaml
done
assert_ok grep -q 'printResolvedConfig' step2.nf
# nextflow.config must no longer carry science parameters
assert_eq "0" "$(grep -cE '^\s+(species_tree|refnames|subs_model|max_spr)\s*=' nextflow.config)" "science params removed from nextflow.config"
```

- [ ] **Step 2: Run it to verify it fails**

Run: `bash tests/test_banner.sh`
Expected: FAIL — no `config/run.example.yaml`.

- [ ] **Step 3: Write `config/run.example.yaml`**

```yaml
# phylohpc run configuration -- THE file that configures a run.
#   nextflow run zolotarovgl/phylohpc -r v1.0 -params-file run.yaml -profile slurm
# Copy this next to your data, edit, and COMMIT IT: run.yaml + the pipeline revision
# printed at start-up are the complete record of what was run.

# ---- inputs
infasta:       data/input.fasta            # proteomes, one multi-fasta
family_info:   genefam.csv                 # family table: name, HMMs, inflation, min hits, cutoff, category, class
ids:           ""                          # optional: file with one HG id per line, restricts the run
refnames:      data/Mmus_gene_names.csv    # gene-name table for the POSSVM name transfer
refspecies:    Mmus                        # POSSVM reference species -- EVERY name comes from this species

# ---- outputs
outdir:        results

# ---- species tree
# MUST be strictly binary when run_generax is true. Resolve polytomies first:
#   python workflow/check_tree.py in.newick ids out.newick --random-resolve --seed 1
# and note that randomly resolved nodes carry no evidence: never read a clade claim off one.
species_tree:  data/species_tree.newick

# ---- methods
tree_method:   iqtree2                     # iqtree2 | fasttree
iqtree2_model: TEST
mafft_opt:     "--maxiterate 1000 --localpair"
subs_model:    LG

# ---- reconciliation (optional)
run_generax:   false                       # needs an MPI GeneRax on PATH; the bioconda build is non-MPI
max_spr:       3
ncpu_generax:  2

# ---- failure policy
max_failed_frac:  0.05                     # fail the run above this fraction of failed families
max_failed_count: 20                       # ... or this many, whichever is smaller
```

- [ ] **Step 4: Add the banner to `step2.nf`**

Insert immediately after the params resolution block:
```groovy
// Print what actually got resolved, BEFORE the DAG is built. The 2026-08-14 failure was
// invisible partly because nothing ever stated which file configured the run.
def printResolvedConfig(String step) {
    log.info """
== resolved configuration (${step})
   params file   : ${workflow.configFiles.join(', ')}
   revision      : ${workflow.revision ?: 'local checkout'} (${workflow.commitId ?: 'no commit id'})
   outdir        : ${canonical('outdir','OUTDIR')}
   work dir      : ${workflow.workDir}
   species tree  : ${canonical('species_tree','SPECIES_TREE')}
   run_generax   : ${params.run_generax}
   tree method   : ${canonical('tree_method','TREE_METHOD')}
   family info   : ${params.family_info}
   ids restrict  : ${params.ids ?: '(all clusters)'}
"""
}
printResolvedConfig('step2')
```

- [ ] **Step 5: Strip science params from `nextflow.config`**

Remove `iqtree2_model`, `refspecies`, `refnames`, `species_tree`, `subs_model`, `max_spr`, `ncpu_generax`, `tree_method`, `mafft_opt` from the `params` block and from the `fast`/`precise` profiles; leave executor, queue, cpus, memory, retry and resource-model settings. Leave `outdir` with its default (a run without an outdir must still work).

- [ ] **Step 6: Run tests and commit**

Run: `bash tests/run_tests.sh` → `ALL TESTS PASSED`
```bash
git add config/run.example.yaml step2.nf nextflow.config tests/test_banner.sh
git commit -m "feat(config): run.yaml is the only science config; every run prints its resolved configuration and revision"
```

---

### Task 4: Pre-flight, in the workflow and in the wrapper

**Files:**
- Modify: `step2.nf` (extend the existing species-tree guard into a single `preflightChecks()`)
- Create: `workflow/preflight.sh`
- Create: `tests/test_preflight.sh`

**Interfaces:**
- Consumes: `workflow/count_multifurcations.py` (Task 1), `canonical()` (Task 2).
- Produces: `workflow/preflight.sh` usable as `bash workflow/preflight.sh <params.yaml> <step>`, exiting non-zero with a named cause; `preflightChecks()` in `step2.nf` doing the same before the DAG.

- [ ] **Step 1: Write the failing test**

`tests/test_preflight.sh`:
```bash
#!/usr/bin/env bash
set -uo pipefail
source "$(dirname "$0")/lib.sh"
cd "$(dirname "$0")/.."
tmp=$(mktemp -d)

cat > "$tmp/ok.yaml" <<Y
infasta: data/input.fasta
family_info: data/gene_families_searchinfo.csv
refnames: data/Mmus_gene_names.csv
species_tree: tests/data/binary.newick
run_generax: false
outdir: $tmp/results
Y
assert_ok bash workflow/preflight.sh "$tmp/ok.yaml" step2

cat > "$tmp/bad_tree.yaml" <<Y
infasta: data/input.fasta
family_info: data/gene_families_searchinfo.csv
refnames: data/Mmus_gene_names.csv
species_tree: tests/data/multifurcating.newick
run_generax: true
outdir: $tmp/results
Y
assert_fails bash workflow/preflight.sh "$tmp/bad_tree.yaml" step2
out=$(bash workflow/preflight.sh "$tmp/bad_tree.yaml" step2 2>&1 || true)
assert_contains "$out" "multifurcating" "names the actual problem"
assert_contains "$out" "check_tree.py"  "tells you how to fix it"

cat > "$tmp/missing.yaml" <<Y
infasta: does/not/exist.fasta
family_info: data/gene_families_searchinfo.csv
refnames: data/Mmus_gene_names.csv
species_tree: tests/data/binary.newick
run_generax: false
outdir: $tmp/results
Y
assert_fails bash workflow/preflight.sh "$tmp/missing.yaml" step2
rm -rf "$tmp"
```

- [ ] **Step 2: Run it to verify it fails**

Run: `bash tests/test_preflight.sh`
Expected: FAIL — `workflow/preflight.sh` missing.

- [ ] **Step 3: Write `workflow/preflight.sh`**

```bash
#!/usr/bin/env bash
# Pre-flight for a phylohpc run. Reads the params YAML, checks every input, and fails
# with a NAMED cause in about a second.
#
# It exists because the 2026-08-14 run staged a multifurcating species tree, GeneRax said
# only "Cannot parse species tree file", the errorStrategy relabelled that as a family
# parsing error, and the run reported success with an empty output directory.
set -uo pipefail
YAML="${1:?usage: preflight.sh <params.yaml> <step>}"
STEP="${2:-step2}"
cd "$(dirname "$0")/.."

get() { python3 -c "import sys,yaml;print(yaml.safe_load(open(sys.argv[1])).get(sys.argv[2],'') or '')" "$YAML" "$1"; }

fail=0
say()  { printf '   %s\n' "$*"; }
err()  { printf 'PREFLIGHT ERROR: %s\n' "$*" >&2; fail=1; }

[[ -f "$YAML" ]] || { err "params file not found: $YAML"; exit 1; }
say "params file : $YAML"

for key in infasta family_info refnames; do
  v=$(get "$key")
  [[ -n "$v" ]] || { err "'$key' is not set in $YAML"; continue; }
  [[ -e "$v" ]] || err "'$key' points at a missing file: $v"
done

tree=$(get species_tree)
generax=$(get run_generax)
if [[ "$generax" == "True" || "$generax" == "true" || "$generax" == "1" ]]; then
  [[ -e "$tree" ]] || err "run_generax is on but species_tree is missing: $tree"
  if [[ -e "$tree" ]]; then
    n=$(python3 workflow/count_multifurcations.py "$tree")
    if [[ "$n" != "0" ]]; then
      err "species tree '$tree' has $n multifurcating node(s); GeneRax needs a strictly binary tree.
  Resolve it:  python workflow/check_tree.py '$tree' <ids> <out.newick> --random-resolve --seed 1
  then set species_tree to that file. (This is what killed the 2026-08-14 run.)"
    else
      say "species tree: $tree (strictly binary)"
    fi
  fi
  command -v generax >/dev/null || say "WARNING: generax not on PATH -- the bioconda build is non-MPI; load an MPI build"
  command -v mpirun  >/dev/null || say "WARNING: mpirun not on PATH"
fi

[[ $fail -eq 0 ]] || exit 1
say "preflight OK ($STEP)"
```

- [ ] **Step 4: Extend the in-workflow guard**

In `step2.nf`, wrap the existing species-tree check plus input existence into `preflightChecks()`, called right after `printResolvedConfig('step2')`. Keep the wording identical to `preflight.sh` so the same failure reads the same either way.

- [ ] **Step 5: Run tests and commit**

Run: `bash tests/run_tests.sh` → `ALL TESTS PASSED`
```bash
git add workflow/preflight.sh step2.nf tests/test_preflight.sh
git commit -m "feat(preflight): one named-cause check in both the wrapper and the workflow; multifurcating trees fail in ~1s"
```

---

### Task 5: Failure accounting

**Files:**
- Modify: `step2.nf:416-440` (the GeneRax `errorStrategy`), `step2.nf:43-48` (`workflow.onComplete`)
- Create: `tests/test_failure_policy.sh`

**Interfaces:**
- Consumes: `canonical()` (Task 2).
- Produces: `reports/failures.tsv` with header `family\texit_code\tfirst_error\tworkdir`; `params.max_failed_frac` (0.05) and `params.max_failed_count` (20) honoured; `--strict` fails on first failure.

- [ ] **Step 1: Write the failing test**

`tests/test_failure_policy.sh`:
```bash
#!/usr/bin/env bash
set -uo pipefail
source "$(dirname "$0")/lib.sh"
cd "$(dirname "$0")/.."

# the guessed label must be gone
assert_eq "0" "$(grep -c 'Family parsing error' step2.nf)" "no invented failure label remains"
# the policy must be present and configurable
assert_ok grep -q 'max_failed_frac'  step2.nf
assert_ok grep -q 'max_failed_count' step2.nf
assert_ok grep -q 'failures.tsv'     step2.nf
assert_ok grep -q 'strict'           step2.nf
```

- [ ] **Step 2: Run it to verify it fails**

Run: `bash tests/test_failure_policy.sh`
Expected: FAIL — `Family parsing error` is still in `step2.nf:420`.

- [ ] **Step 3: Replace the guessed label with accounting**

In `step2.nf`, above the processes:
```groovy
// Failure accounting. A task NEVER gets a guessed diagnosis: errReport() already prints
// the tail of its own .command.err, and every failure is recorded here so the run cannot
// end green with an empty output directory (2026-08-14: 475 of 475 GeneRax tasks failed
// and the run reported success).
@groovy.transform.Field
def FAILURES = java.util.Collections.synchronizedList([])

def recordFailure(task) {
    def errf = task.workDir?.resolve('.command.err')
    def first = (errf && errf.exists())
        ? (errf.text.readLines().find { it.trim() && !it.startsWith('+') } ?: '(empty stderr)')
        : '(.command.err not found)'
    FAILURES << [task.tag ?: task.index, task.exitStatus, first.take(200), task.workDir]
}
```

Replace the GeneRax `errorStrategy` body:
```groovy
    errorStrategy {
        errReport(task)
        def max_attempts = 5
        if( task.attempt <= max_attempts && task.exitStatus != 10 )
            return 'retry'
        recordFailure(task)
        return params.strict ? 'terminate' : 'ignore'
    }
```
(Exit 10 is not retried — it was deterministic across all 475 tasks — but it is now recorded, not labelled.)

In `workflow.onComplete`, replace the body with:
```groovy
workflow.onComplete {
    def outdir = canonical('outdir','OUTDIR')
    def n_fail = FAILURES.size()
    def n_fam  = file("${outdir}/clusters").exists() ? file("${outdir}/clusters").list().findAll{ it.endsWith('.fasta') }.size() : 0
    if( n_fail ) {
        def rep = file("reports/failures.tsv")
        rep.parent.mkdirs()
        rep.text = "family\texit_code\tfirst_error\tworkdir\n" +
                   FAILURES.collect{ it.join('\t') }.join('\n') + '\n'
        log.warn "${n_fail} task(s) failed -- see reports/failures.tsv"
    }
    // an expected output directory that ends up EMPTY is an error in its own right
    if( params.run_generax ) {
        def gr = file("${outdir}/generax")
        if( !gr.exists() || gr.list().size() == 0 )
            log.error "run_generax was on but ${outdir}/generax is EMPTY -- treat this run as FAILED"
    }
    def frac  = n_fam ? (n_fail / (double) n_fam) : 0
    def limit = Math.min(params.max_failed_count as int, Math.ceil(n_fam * (params.max_failed_frac as double)) as int)
    if( n_fail > limit )
        log.error "step2 FAILED: ${n_fail} of ${n_fam} families failed (${String.format('%.1f', frac*100)}%), over the tolerance of ${limit}"
    else
        log.info "step2 ${workflow.success ? 'completed' : 'FAILED'} after ${workflow.duration}; ${n_fail} failed family(ies) within tolerance"
}
```

Add the defaults near the top of `step2.nf`:
```groovy
params.strict           = params.containsKey('strict') ? params.strict : false
params.max_failed_frac  = params.containsKey('max_failed_frac')  ? params.max_failed_frac  : 0.05
params.max_failed_count = params.containsKey('max_failed_count') ? params.max_failed_count : 20
```

- [ ] **Step 4: Run tests**

Run: `bash tests/run_tests.sh`
Expected: `ALL TESTS PASSED`.

- [ ] **Step 5: Commit**

```bash
git add step2.nf tests/test_failure_policy.sh
git commit -m "feat(failures): record instead of label; failures.tsv, a tolerance, --strict, and an empty-output check"
```

---

### Task 6: The `phylohpc` entry point

**Files:**
- Create: `phylohpc`
- Delete: `submit_nf.sh`, `run_interactive.sh`, `make_report.sh`
- Create: `tests/test_cli.sh`

**Interfaces:**
- Consumes: `workflow/preflight.sh` (Task 4), `reports/failures.tsv` (Task 5), `params.ids` (already in `step2.nf:132`).
- Produces: `phylohpc {run|submit|rerun-failed|report|init} <step> -p <params.yaml> [-- extra nextflow args]`.

- [ ] **Step 1: Write the failing test**

`tests/test_cli.sh`:
```bash
#!/usr/bin/env bash
set -uo pipefail
source "$(dirname "$0")/lib.sh"
cd "$(dirname "$0")/.."

assert_ok test -x phylohpc
out=$(./phylohpc 2>&1 || true);        assert_contains "$out" "usage"        "bare invocation prints usage"
out=$(./phylohpc bogus 2>&1 || true);  assert_contains "$out" "unknown"      "unknown subcommand is rejected"
assert_fails ./phylohpc run step2 -p does_not_exist.yaml
# init writes a params file that preflight accepts
tmp=$(mktemp -d); ( cd "$tmp" && "$OLDPWD/phylohpc" init >/dev/null 2>&1 ); assert_ok test -f "$tmp/run.yaml"; rm -rf "$tmp"
# rerun-failed needs a failures.tsv and says so
out=$(./phylohpc rerun-failed step2 -p config/run.example.yaml 2>&1 || true)
assert_contains "$out" "failures.tsv" "rerun-failed explains what it needs"
```

- [ ] **Step 2: Run it to verify it fails**

Run: `bash tests/test_cli.sh`
Expected: FAIL — `phylohpc` does not exist.

- [ ] **Step 3: Write `phylohpc`**

```bash
#!/usr/bin/env bash
# phylohpc -- the one entry point.
#
#   phylohpc run          step2 -p run.yaml      # foreground (interactive salloc/srun)
#   phylohpc submit       step2 -p run.yaml      # sbatch
#   phylohpc rerun-failed step2 -p run.yaml      # only the families in reports/failures.tsv
#   phylohpc report                              # build the HTML report
#   phylohpc init                                # write a commented run.yaml here
#
# Every mode runs the same pre-flight; there is no second front door that skips it.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

usage() { sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-1}"; }
[[ $# -ge 1 ]] || usage

CMD="$1"; shift || true
STEP=""; PARAMS=""; PROFILE="${PROFILE:-local,precise}"; WORKDIR="${WORKDIR:-work/}"; EXTRA=()
case "$CMD" in run|submit|rerun-failed) STEP="${1:?need a step, e.g. step2}"; shift || true ;; esac
while [[ $# -gt 0 ]]; do
  case "$1" in
    -p|--params)  PARAMS="$2"; shift 2 ;;
    --profile)    PROFILE="$2"; shift 2 ;;
    -w|--workdir) WORKDIR="$2"; shift 2 ;;
    --) shift; EXTRA+=("$@"); break ;;
    *)  EXTRA+=("$1"); shift ;;
  esac
done

run_nf() {
  local step="$1"; shift
  bash "$HERE/workflow/preflight.sh" "$PARAMS" "$step"
  nextflow run "$HERE/$step.nf" -ansi-log false -profile "$PROFILE" -resume \
      -w "$WORKDIR" -params-file "$PARAMS" \
      -with-report "reports/report.$step.html" -with-trace "reports/trace.$step.txt" "$@"
}

case "$CMD" in
  init)
    [[ -f run.yaml ]] && { echo "run.yaml already exists here"; exit 1; }
    cp "$HERE/config/run.example.yaml" run.yaml
    echo "wrote run.yaml -- edit it, commit it: it is the record of what you ran"
    ;;
  report)
    python3 "$HERE/workflow/report_step2.py" "${EXTRA[@]}"
    ;;
  run)
    [[ -n "$PARAMS" ]] || usage
    mkdir -p reports logs
    run_nf "$STEP" "${EXTRA[@]}" 2>&1 | tee "logs/$STEP.$(date +%Y%m%d-%H%M%S).log"
    ;;
  submit)
    [[ -n "$PARAMS" ]] || usage
    mkdir -p reports logs
    bash "$HERE/workflow/preflight.sh" "$PARAMS" "$STEP"   # fail before queueing
    sbatch --job-name "$STEP" --output "reports/slurm.$STEP.out" \
           --wrap "cd $PWD && $HERE/phylohpc run $STEP -p $PARAMS --profile ${PROFILE/local/slurm} -w $WORKDIR ${EXTRA[*]:-}"
    ;;
  rerun-failed)
    [[ -f reports/failures.tsv ]] || { echo "no reports/failures.tsv here -- run the step first"; exit 1; }
    mkdir -p reports
    tail -n +2 reports/failures.tsv | cut -f1 | sort -u > reports/rerun.ids
    echo "re-running $(wc -l < reports/rerun.ids) failed families (reports/rerun.ids)"
    run_nf "$STEP" --ids reports/rerun.ids "${EXTRA[@]}"
    ;;
  *) echo "unknown subcommand: $CMD"; usage ;;
esac
```

- [ ] **Step 4: Delete the old entry points**

```bash
git rm submit_nf.sh run_interactive.sh make_report.sh
```

- [ ] **Step 5: Run tests**

Run: `chmod +x phylohpc && bash tests/run_tests.sh`
Expected: `ALL TESTS PASSED`.

- [ ] **Step 6: Commit**

```bash
git add phylohpc tests/test_cli.sh
git commit -m "feat(cli): one phylohpc entry point (run/submit/rerun-failed/report/init); pre-flight cannot be bypassed"
```

---

### Task 7: Remove the Snakemake implementation and update the docs

**Files:**
- Delete: `step1.smk`, `step2.smk`, `workflow/step4_ancestry.smk`
- Modify: `README.md`, `docs/manual.md`
- Create: `tests/test_docs.sh`

**Interfaces:**
- Consumes: everything above.
- Produces: documentation whose commands match the CLI.

- [ ] **Step 1: Write the failing test**

`tests/test_docs.sh`:
```bash
#!/usr/bin/env bash
set -uo pipefail
source "$(dirname "$0")/lib.sh"
cd "$(dirname "$0")/.."

assert_eq "0" "$(ls *.smk workflow/*.smk 2>/dev/null | wc -l)" "no Snakemake files remain"
assert_eq "0" "$(grep -rc 'snakemake' README.md | cut -d: -f1 | paste -sd+ | bc 2>/dev/null || echo 0)" "README no longer documents snakemake"
assert_ok grep -q 'phylohpc run'    README.md
assert_ok grep -q 'phylohpc submit' README.md
assert_ok grep -q 'params-file'     README.md
assert_eq "0" "$(grep -c 'submit_nf.sh' README.md)" "no references to the deleted wrapper"
```

- [ ] **Step 2: Run it to verify it fails**

Run: `bash tests/test_docs.sh`
Expected: FAIL — the `.smk` files exist.

- [ ] **Step 3: Delete and rewrite**

```bash
git rm step1.smk step2.smk workflow/step4_ancestry.smk
```
Rewrite the README's Run section to:
````markdown
### Run

```bash
phylohpc init                          # writes a commented run.yaml
$EDITOR run.yaml                       # point it at your data; COMMIT IT
phylohpc submit step1 -p run.yaml      # search + cluster
phylohpc submit step2 -p run.yaml      # align + trees + [generax] + possvm
phylohpc rerun-failed step2 -p run.yaml
```

Or straight from GitHub, pinned to a revision — `run.yaml` plus that revision is the
complete record of a run:

```bash
nextflow run zolotarovgl/phylohpc -r v1.0 -params-file run.yaml -profile slurm
```
````
Delete the Snakemake sections from `docs/manual.md`, including the line describing
`config/step2.yaml` as the Snakemake config.

- [ ] **Step 4: Run tests and commit**

Run: `bash tests/run_tests.sh` → `ALL TESTS PASSED`
```bash
git add -A README.md docs/manual.md tests/test_docs.sh
git commit -m "refactor: drop the Snakemake implementation; docs describe the phylohpc CLI and params-file only"
```

---

### Task 8: Verification against the previous run (REQUIRES THE CLUSTER)

**Files:**
- Create: `docs/verification-2026-08-refactor.md`
- Create: `~/ant/gzolotarov/projects/2022_scRNAseq_Mlei/Phylogenies_v3/run.yaml`

**Interfaces:**
- Consumes: the whole CLI.
- Produces: a written comparison against `results/possvm_prev/`.

⚠ **This task cannot be done from the workstation** — it needs an interactive HPC session. An agent should prepare the `run.yaml` and the comparison script and then STOP, leaving the run to Grygoriy.

- [ ] **Step 1: Write `run.yaml` for `Phylogenies_v3`**

Reconstruct from that directory's `.nextflow.log` (grep for the `nextflow run` line and the resolved params) and from `config/step2.yaml`. Set `run_generax: false` for the comparison run and `species_tree: data/species_tree.binary.newick`.

- [ ] **Step 2: Write the comparison script**

`docs/verification-2026-08-refactor.md` records the command and the expected result:
```bash
phylohpc run step2 -p run.yaml -w work_verify/ --outdir results_verify
diff <(sort results/possvm_prev/*.ortholog_groups.csv | md5sum) \
     <(sort results_verify/possvm/*.ortholog_groups.csv | md5sum)
```
Expected: identical. A difference means the refactor changed behaviour and must be
explained before anything is tagged.

- [ ] **Step 3: Negative tests on the cluster**

Each must fail fast with a named cause: multifurcating `species_tree`; missing `family_info`; `generax` absent with `run_generax: true`; a `run.yaml` key that does not exist.

- [ ] **Step 4: Record the outcome and commit**

Write the measured result — not the intent — into `docs/verification-2026-08-refactor.md`, then:
```bash
git add docs/verification-2026-08-refactor.md
git commit -m "docs: verification of the refactor against results/possvm_prev"
```

- [ ] **Step 5: Tag**

Only after Step 4 passes:
```bash
git tag -a v1.0 -m "config/provenance/failure-reporting refactor"
git push origin refactor/config-provenance --tags
```

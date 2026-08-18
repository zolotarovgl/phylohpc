#!/usr/bin/env bash
set -uo pipefail
source "$(dirname "$0")/lib.sh"
cd "$(dirname "$0")/.."

# no uppercase param may be READ anywhere except through the deprecation shim
# (canonical() itself references the uppercase names as string literals and via
# params[upper] -- not as params.UPPER -- so it does not trip this grep)
hits=$(grep -nE 'params\.(SPECIES_TREE|OUTDIR|REFNAMES|REFSPECIES|SUBS_MODEL|MAX_SPR|NCPU_GENERAX|TREE_METHOD|MAFFT_OPT|IQTREE2_MODEL)' step1.nf step2.nf | wc -l)
assert_eq "0" "$hits" "no direct uppercase params.* reads outside the shim"

# every param the .nf files read must exist in nextflow.config or be defined with a default
assert_ok grep -q 'species_tree' nextflow.config
# canonical() lives in step2.nf directly -- evaluate("workflow/params.nf") was tried and
# measured NOT to put the function into DSL2 script scope ("Unknown method invocation
# `canonical`"), so there is no workflow/params.nf; grep the file where it actually lives.
assert_ok grep -q 'def canonical' step2.nf

# the uppercase keys must never be DEFAULTED anywhere in nextflow.config, in any
# profile -- that is what makes containsKey(upper) a reliable "user passed a legacy
# flag" signal for canonical()'s hard-error branch (round 2, 2026-08-18)
cfg_upper_hits=$(grep -nE '^\s*(SPECIES_TREE|OUTDIR|REFNAMES|REFSPECIES|SUBS_MODEL|MAX_SPR|NCPU_GENERAX|TREE_METHOD|MAFFT_OPT|IQTREE2_MODEL)\s*=' nextflow.config | wc -l)
assert_eq "0" "$cfg_upper_hits" "nextflow.config never defaults an uppercase key (any profile)"

# ── Runtime cases: actually exercise canonical()'s shim through step2.nf -preview ──
# (a few seconds each -- these run the real DAG build/params-resolution code, not a grep)
_NF="nextflow run step2.nf -preview"
_TMPDIR=$(mktemp -d)
trap '[[ $_FAILED -eq 0 ]] || exit 1; rm -rf "$_TMPDIR"' EXIT

# 1. uppercase-only: legacy SPECIES_TREE must now be a HARD ERROR (round 2), never a
#    silent pass-through and never a mere warning. The error must name both spellings.
cat > "$_TMPDIR/upper_only.yaml" << 'YAML'
run_generax: true
SPECIES_TREE: tests/data/binary.newick
YAML
out=$($_NF -params-file "$_TMPDIR/upper_only.yaml" 2>&1)
if echo "$out" | grep -q 'step2 completed OK'; then
    echo "  FAIL: uppercase-only: run passed -- legacy SPECIES_TREE was not refused"
    _FAILED=1
else
    echo "  ok: uppercase-only: run correctly failed"
fi
assert_contains "$out" "SPECIES_TREE" "uppercase-only: error names the legacy key"
assert_contains "$out" "species_tree" "uppercase-only: error names the canonical key"

# 2. both keys set: ambiguity is refused too, unconditionally -- lowercase does NOT
#    silently win any more (round 2 supersedes round 1's precedence rule).
cat > "$_TMPDIR/both.yaml" << 'YAML'
run_generax: true
species_tree: tests/data/binary.newick
SPECIES_TREE: tests/data/multifurcating.newick
YAML
out=$($_NF -params-file "$_TMPDIR/both.yaml" 2>&1)
if echo "$out" | grep -q 'step2 completed OK'; then
    echo "  FAIL: both-set: run passed -- ambiguous SPECIES_TREE+species_tree was not refused"
    _FAILED=1
else
    echo "  ok: both-set: run correctly failed (ambiguity refused)"
fi
assert_contains "$out" "SPECIES_TREE" "both-set: error names the legacy key"

# 3. lowercase only: must succeed, no legacy-key error.
cat > "$_TMPDIR/lower_only.yaml" << 'YAML'
run_generax: true
species_tree: tests/data/binary.newick
YAML
out=$($_NF -params-file "$_TMPDIR/lower_only.yaml" 2>&1)
assert_contains "$out" "step2 completed OK" "lowercase-only: run succeeded on the binary tree"

# 4. default applies: neither set -- must resolve the projectDir default and succeed.
cat > "$_TMPDIR/neither.yaml" << 'YAML'
run_generax: true
YAML
out=$($_NF -params-file "$_TMPDIR/neither.yaml" 2>&1)
assert_contains "$out" "step2 completed OK" "default-applies: run succeeded on the projectDir default tree"

# 5. profile-driven method tuning still works: fast -> fasttree, precise -> iqtree2,
#    exercised via the log.info line added in step2.nf for exactly this purpose.
out_fast=$($_NF -profile fast -params-file "$_TMPDIR/neither.yaml" 2>&1)
assert_contains "$out_fast" "tree_method: fasttree" "profile fast resolves tree_method=fasttree"
out_precise=$($_NF -profile precise -params-file "$_TMPDIR/neither.yaml" 2>&1)
assert_contains "$out_precise" "tree_method: iqtree2" "profile precise resolves tree_method=iqtree2"

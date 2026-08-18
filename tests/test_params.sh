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

# ── Runtime cases: actually exercise canonical()'s shim through step2.nf -preview ──
# (a few seconds each -- these run the real DAG build/params-resolution code, not a grep)
_NF="nextflow run step2.nf -preview"
_TMPDIR=$(mktemp -d)
trap '[[ $_FAILED -eq 0 ]] || exit 1; rm -rf "$_TMPDIR"' EXIT

# 1. uppercase-only: legacy SPECIES_TREE must actually be USED (not silently dropped)
#    and must warn. Uses the multifurcating fixture so a pass-through failure is loud.
cat > "$_TMPDIR/upper_only.yaml" << 'YAML'
run_generax: true
SPECIES_TREE: tests/data/multifurcating.newick
YAML
out=$($_NF -params-file "$_TMPDIR/upper_only.yaml" 2>&1)
assert_contains "$out" "DEPRECATED parameter 'SPECIES_TREE'" "uppercase-only: warns"
assert_contains "$out" "multifurcating node(s)" "uppercase-only: legacy value was actually used (run failed on it)"
if echo "$out" | grep -q 'step2 completed OK'; then
    echo "  FAIL: uppercase-only: run passed -- SPECIES_TREE was silently ignored"
    _FAILED=1
else
    echo "  ok: uppercase-only: run correctly failed on the legacy tree"
fi

# 2. lowercase wins: both set, opposite trees -- must succeed (lowercase binary tree
#    used) and must NOT warn about SPECIES_TREE.
cat > "$_TMPDIR/both.yaml" << 'YAML'
run_generax: true
species_tree: tests/data/binary.newick
SPECIES_TREE: tests/data/multifurcating.newick
YAML
out=$($_NF -params-file "$_TMPDIR/both.yaml" 2>&1)
if echo "$out" | grep -q "DEPRECATED parameter 'SPECIES_TREE'"; then
    echo "  FAIL: lowercase-wins: warned even though lowercase was set"
    _FAILED=1
else
    echo "  ok: lowercase-wins: no deprecation warning"
fi
assert_contains "$out" "step2 completed OK" "lowercase-wins: run succeeded on the binary tree"

# 3. default applies: neither set -- must resolve the projectDir default and succeed,
#    with no warning.
cat > "$_TMPDIR/neither.yaml" << 'YAML'
run_generax: true
YAML
out=$($_NF -params-file "$_TMPDIR/neither.yaml" 2>&1)
if echo "$out" | grep -q "DEPRECATED parameter"; then
    echo "  FAIL: default-applies: warned with neither uppercase nor lowercase set"
    _FAILED=1
else
    echo "  ok: default-applies: no deprecation warning"
fi
assert_contains "$out" "step2 completed OK" "default-applies: run succeeded on the projectDir default tree"

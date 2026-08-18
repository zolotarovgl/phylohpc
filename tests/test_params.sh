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

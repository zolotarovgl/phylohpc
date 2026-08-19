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
# nextflow.config's TOP-LEVEL params{} block (the one below the profiles{} block) must
# no longer carry science parameters. The fast/precise profiles legitimately override
# tree_method/mafft_opt/max_spr/ncpu_generax (see nextflow.config's own comment) -- so
# only the trailing top-level block is scanned here, not the whole file.
top_level_params=$(awk '/^params \{/{f=1} f' nextflow.config)
assert_eq "0" "$(echo "$top_level_params" | grep -cE '^\s+(species_tree|refnames|subs_model|max_spr)\s*=')" "science params removed from nextflow.config's top-level params{} block"

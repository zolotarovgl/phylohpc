#!/usr/bin/env bash
set -uo pipefail
source "$(dirname "$0")/lib.sh"
cd "$(dirname "$0")/.."

assert_eq "0" "$(ls *.smk workflow/*.smk 2>/dev/null | wc -l)" "no Snakemake files remain"
assert_eq "0" "$(grep -c 'snakemake' README.md || true)" "README no longer documents snakemake"
assert_ok grep -q 'phylohpc run'    README.md
assert_ok grep -q 'phylohpc submit' README.md
assert_ok grep -q 'params-file'     README.md
assert_eq "0" "$(grep -c 'submit_nf.sh' README.md)" "no references to the deleted wrapper"

# Every wrapper deleted in this refactor must leave no dangling reference anywhere
# under docs/ or in a tracked *.md at the repo root. Historical design docs that
# describe the deletion itself (docs/superpowers/plans|specs/, and the pre-existing
# docs/pipeline_improvement_plan.md analysis) are exempt -- they are the record of
# the decision, not live usage instructions.
DELETED_NAMES=(submit_nf.sh run_interactive.sh make_report.sh '\.smk')
DOC_FILES=$(git ls-files '*.md' 'docs/**/*.md' \
    | grep -v '^docs/superpowers/' \
    | grep -v '^docs/pipeline_improvement_plan\.md$' \
    | sort -u)
for name in "${DELETED_NAMES[@]}"; do
    hits=$(grep -l "$name" $DOC_FILES 2>/dev/null | wc -l)
    assert_eq "0" "$hits" "no live doc references '$name'"
done

# CRITICAL 3 (2026-08-19 whole-branch review): config/step1.yaml and config/step2.yaml
# were still tracked, self-describing as "Snakemake configuration for step*.smk", while
# the .smk files were deleted -- exactly "a file nothing reads" from the 2026-08-14
# incident, just relocated. Guard against it recurring: no TRACKED file may self-describe
# as a Snakemake config, and no tracked file may reference a deleted .smk by name.
# Exclude this test itself (it necessarily contains the grep patterns it's checking
# for) and archive/ (retired code kept for provenance, never rerun -- see archive/README.md).
ALL_TRACKED=$(git ls-files | grep -v '^tests/test_docs\.sh$' | grep -v '^archive/')
snakemake_config_hits=$(grep -lE '^\s*#.*Snakemake configuration for' $ALL_TRACKED 2>/dev/null | wc -l)
assert_eq "0" "$snakemake_config_hits" "no tracked file self-describes as a Snakemake config"

smk_ref_hits=$(grep -lE '\.smk\b' $ALL_TRACKED 2>/dev/null \
    | grep -v '^docs/superpowers/' \
    | grep -v '^docs/pipeline_improvement_plan\.md$' \
    | wc -l)
assert_eq "0" "$smk_ref_hits" "no live tracked file references a .smk file"

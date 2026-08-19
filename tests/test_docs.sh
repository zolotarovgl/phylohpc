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

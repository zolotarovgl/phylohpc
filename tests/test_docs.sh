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

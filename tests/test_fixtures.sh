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

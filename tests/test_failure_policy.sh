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

# ── behavioural tests of the actual tolerance/writer logic (workflow/failure_policy.py),
# which step2.nf shells out to -- this is not string-grepping, it runs the real function.
PY="python3 workflow/failure_policy.py"

# tolerance_limit / decide table: n_failed, n_families, max_frac, max_count, strict -> expect
check_decide() {
    got=$($PY decide "$1" "$2" "$3" "$4" "$5")
    assert_eq "$6" "$got" "decide(failed=$1 fam=$2 frac=$3 count=$4 strict=$5)"
}

# 0 failures -> always ok
check_decide 0   500 0.05 20 false ok
# 1 of 500 -> well within both frac (25) and count (20) -> ok
check_decide 1   500 0.05 20 false ok
# 20 of 500 -> frac limit is 25, count limit is 20, min=20 -> at the boundary -> ok
check_decide 20  500 0.05 20 false ok
# 21 of 500 -> exceeds the count limit (20) even though frac limit (25) would allow it -> fail
check_decide 21  500 0.05 20 false fail
# 5% boundary exactly, small n where frac limit binds under count limit: 25 fam, 5%=1.25->ceil 2, count 20 -> min=2
check_decide 2   25  0.05 20 false ok
check_decide 3   25  0.05 20 false fail
# all failed (500 of 500) -> fail regardless of tolerance
check_decide 500 500 0.05 20 false fail
# --strict: even 1 failure fails the run, no matter the tolerance
check_decide 1   500 0.05 20 true  fail
check_decide 0   500 0.05 20 true  ok

# ── TSV writer: fabricated input -> exact 5-column TSV with header. The `process`
# column (added 2026-08-19) records WHICH stage failed -- see CRITICAL 2, failures used
# to be recorded from GR_watcher only, so a non-GeneRax failure was invisible.
tmp=$(mktemp -d)
printf 'HG0001\t10\tCannot parse species tree file\t/work/ab/cd1234\tGR_watcher\nHG0002\t137\tOut of memory\t/work/ef/gh5678\tALN\n' \
    | $PY write-tsv "$tmp/failures.tsv"
expected=$'family\texit_code\tfirst_error\tworkdir\tprocess\nHG0001\t10\tCannot parse species tree file\t/work/ab/cd1234\tGR_watcher\nHG0002\t137\tOut of memory\t/work/ef/gh5678\tALN'
got=$(cat "$tmp/failures.tsv")
assert_eq "$expected" "$got" "failures.tsv exact 5-column content with header"

# malformed row (wrong field count) must be rejected, not silently written
if printf 'HG0003\t10\tonly-three-fields\n' | $PY write-tsv "$tmp/bad.tsv" >/dev/null 2>&1; then
    echo "  FAIL: malformed row was silently accepted"
    _FAILED=1
else
    echo "  ok: malformed row rejected"
fi

rm -rf "$tmp"

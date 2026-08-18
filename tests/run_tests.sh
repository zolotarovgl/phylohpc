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

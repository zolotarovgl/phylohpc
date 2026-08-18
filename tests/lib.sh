#!/usr/bin/env bash
# Minimal assertions. Every test script sources this and exits non-zero on failure.
_FAILED=0
assert_eq()       { if [[ "$1" != "$2" ]]; then echo "  FAIL: ${3:-assert_eq}: expected '$1', got '$2'"; _FAILED=1; else echo "  ok: ${3:-assert_eq}"; fi; }
assert_contains() { if [[ "$1" != *"$2"* ]]; then echo "  FAIL: ${3:-assert_contains}: '$2' not found"; _FAILED=1; else echo "  ok: ${3:-assert_contains}"; fi; }
assert_ok()       { if "$@" >/dev/null 2>&1; then echo "  ok: $*"; else echo "  FAIL: expected success: $*"; _FAILED=1; fi; }
assert_fails()    { if "$@" >/dev/null 2>&1; then echo "  FAIL: expected failure: $*"; _FAILED=1; else echo "  ok (failed as expected): $*"; fi; }
trap '[[ $_FAILED -eq 0 ]] || exit 1' EXIT

#!/usr/bin/env bash
set -uo pipefail
source "$(dirname "$0")/lib.sh"
cd "$(dirname "$0")/.."

assert_ok test -x phylohpc
out=$(./phylohpc 2>&1 || true);        assert_contains "$out" "usage"        "bare invocation prints usage"
out=$(./phylohpc bogus 2>&1 || true);  assert_contains "$out" "unknown"      "unknown subcommand is rejected"
assert_fails ./phylohpc run step2 -p does_not_exist.yaml
# init writes a params file that preflight accepts
tmp=$(mktemp -d); ( cd "$tmp" && "$OLDPWD/phylohpc" init >/dev/null 2>&1 ); assert_ok test -f "$tmp/run.yaml"; rm -rf "$tmp"
# rerun-failed needs a failures.tsv and says so
out=$(./phylohpc rerun-failed step2 -p config/run.example.yaml 2>&1 || true)
assert_contains "$out" "failures.tsv" "rerun-failed explains what it needs"

#!/usr/bin/env bash
# Real-execution test of the FAILURES accumulation pattern step2.nf uses: a top-level
# `def FAILURES = ...list`, appended to INLINE inside an errorStrategy closure, read from
# workflow.onComplete. Runs a real nextflow pipeline (no -preview) with a process that
# actually exits 10, so this exercises the real Groovy scoping bug the reviewer found
# (a separate top-level recordFailure(task) method throws MissingPropertyException when
# called from errorStrategy) and confirms the inline-append fix actually works, in both
# places, for real.
set -uo pipefail
source "$(dirname "$0")/lib.sh"
cd "$(dirname "$0")/.."

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

cat > "$tmp/fixture.nf" <<'NF'
nextflow.enable.dsl=2

def FAILURES = java.util.Collections.synchronizedList([])

process FAILER {
    tag 'boom'
    errorStrategy {
        // EXACT same pattern as step2.nf's GR_watcher: append happens INLINE in the
        // errorStrategy closure body, not via a separate top-level method.
        FAILURES << ['boom', task.exitStatus, 'deliberate exit 10 for the test', task.workDir]
        return 'ignore'
    }
    script:
    """
    exit 10
    """
}

workflow.onComplete {
    def rep = file("failures.tsv")
    rep.text = "family\texit_code\tfirst_error\tworkdir\n" +
               FAILURES.collect{ it.join('\t') }.join('\n') + (FAILURES.size() ? '\n' : '')
    log.info "FAILURES.size()=${FAILURES.size()}"
}

workflow {
    FAILER()
}
NF

cd "$tmp"
out=$(nextflow run fixture.nf 2>&1)
rc=$?
cd - >/dev/null

echo "$out" | tail -20

assert_eq "0" "$(echo "$out" | grep -c 'MissingPropertyException')" "no MissingPropertyException (the scoping bug) in the real run"
assert_eq "1" "$(echo "$out" | grep -c 'FAILURES.size()=1')" "onComplete saw exactly 1 recorded failure"
assert_ok test -f "$tmp/failures.tsv"
n_data_rows=$(tail -n +2 "$tmp/failures.tsv" | grep -c '.')
assert_eq "1" "$n_data_rows" "failures.tsv has exactly one data row"
assert_contains "$(cat "$tmp/failures.tsv")" "boom	10	deliberate exit 10 for the test" "failures.tsv row content"

rm -rf "$tmp/work" "$tmp/.nextflow"* 2>/dev/null

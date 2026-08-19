#!/usr/bin/env bash
# Real-execution test (no -preview) for CRITICAL 2 (2026-08-19 whole-branch review):
# FAILURES was only ever appended from GR_watcher (GeneRax) -- an ALN/PHY/PVM failure
# 'ignore'd with nothing recorded, n_fail stayed 0, and the run reported success with an
# empty output directory. This proves TWO different (fictional) process names both land in
# FAILURES with a `process` column, and that the REAL workflow/failure_policy.py writer
# (shelled out to, not reimplemented) produces the correct 5-column failures.tsv from them
# -- end to end, the same way step2.nf's own onComplete does it.
#
# Same self-contained-fixture technique as tests/test_failure_recording.sh (which this
# extends): a real nextflow run, not a grep. Deliberately does NOT exercise retries --
# reproduced separately, this sandbox's nextflow build (25.10.4) has a retry bug where a
# multi-statement dynamic errorStrategy closure never actually retries regardless of what
# it returns, which is an unrelated engine limitation, not a CRITICAL 2 concern (CRITICAL 2
# is about WHICH processes get recorded, not the retry ladder).
set -uo pipefail
source "$(dirname "$0")/lib.sh"
cd "$(dirname "$0")/.."

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

cat > "$tmp/fixture.nf" <<'NF'
nextflow.enable.dsl=2

def FAILURES = java.util.Collections.synchronizedList([])

// Stands in for ALN -- a non-GeneRax process that used to fail silently, uncounted.
process ALN_LIKE {
    tag 'fam1'
    errorStrategy {
        FAILURES << ['fam1', task.exitStatus, 'deliberate exit 3 for the test', task.workDir, 'ALN']
        return 'ignore'
    }
    script:
    """
    exit 3
    """
}

// A second, different process -- proves the column records WHICH stage failed, not just
// THAT one failed.
process PVM_LIKE {
    tag 'fam2'
    errorStrategy {
        FAILURES << ['fam2', task.exitStatus, 'deliberate exit 4 for the test', task.workDir, 'PVM']
        return 'ignore'
    }
    script:
    """
    exit 4
    """
}

workflow.onComplete {
    def rep = file("reports/failures.tsv")
    rep.parent.mkdirs()
    def rowLines = FAILURES.collect{ it.join('\t') }.join('\n') + '\n'
    def writeProc = ["python3", "workflow/failure_policy.py", "write-tsv", rep.toString()].execute()
    writeProc.getOutputStream().withWriter { it << rowLines }
    writeProc.waitFor()
    log.info "writeProc.exitValue()=${writeProc.exitValue()}"
    log.info "FAILURES.size()=${FAILURES.size()}"
}

workflow {
    ALN_LIKE()
    PVM_LIKE()
}
NF

cp -r workflow "$tmp/workflow"

cd "$tmp"
out=$(nextflow run fixture.nf 2>&1)
rc=$?
cd - >/dev/null

echo "$out" | tail -20

assert_eq "1" "$(echo "$out" | grep -c 'writeProc.exitValue()=0')" "failure_policy.py write-tsv succeeded"
assert_eq "1" "$(echo "$out" | grep -c 'FAILURES.size()=2')" "onComplete saw exactly 2 recorded failures, from two different processes"
assert_ok test -f "$tmp/reports/failures.tsv"

got_header=$(head -1 "$tmp/reports/failures.tsv")
assert_eq "family	exit_code	first_error	workdir	process" "$got_header" "failures.tsv has the 5-column header including process"

aln_row=$(grep '^fam1' "$tmp/reports/failures.tsv")
pvm_row=$(grep '^fam2' "$tmp/reports/failures.tsv")
assert_contains "$aln_row" "	ALN" "fam1's row records process=ALN (not GeneRax)"
assert_contains "$pvm_row" "	PVM" "fam2's row records process=PVM (not GeneRax)"

rm -rf "$tmp/work" "$tmp/.nextflow"* 2>/dev/null

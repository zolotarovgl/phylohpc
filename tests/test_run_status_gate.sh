#!/usr/bin/env bash
# Real-execution test (no -preview) for CRITICAL 1 (2026-08-19 whole-branch review): a
# FAILED step2 run still exited 0, because everything workflow.onComplete could do was
# log.error, which does not touch nextflow's own exit code. step2.nf now writes
# reports/run_status.txt ("ok"/"fail") as the last act of onComplete, and `phylohpc` reads
# it after nextflow returns.
#
# This exercises the REAL phylohpc CLI logic (run_nf / check_run_status / the tee pipeline
# and its exit-code propagation) against a real, non-preview nextflow run -- but the run is
# of a small self-contained fixture *.nf standing in for step2.nf, not the full pipeline:
# reproduced separately (isolation confirmed with a two-line onComplete, same technique as
# test_failure_recording.sh) is a genuine RETRY BUG in this sandbox's nextflow build
# (25.10.4): a dynamic errorStrategy closure with more than one statement never actually
# retries, unconditionally on the first failed attempt, regardless of what it returns --
# so exercising step2.nf's real ALN/PHY/GR_watcher retry ladders end-to-end here is not
# possible without either retrying zero times (not what production does) or hitting an
# unrelated engine limitation that has nothing to do with this review's findings. A fixture
# with maxRetries not in play (the same restriction test_failure_recording.sh already
# works under) isolates exactly the phylohpc-side logic CRITICAL 1 changed.
set -uo pipefail
source "$(dirname "$0")/lib.sh"
cd "$(dirname "$0")/.."
HERE="$PWD"

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

mkdir -p "$tmp/proj/workflow" "$tmp/proj/config"
cp "$HERE/phylohpc" "$tmp/proj/phylohpc"
chmod +x "$tmp/proj/phylohpc"

# A no-op preflight -- this test is about the run_status.txt gate, not preflight itself.
cat > "$tmp/proj/workflow/preflight.sh" <<'SH'
#!/usr/bin/env bash
echo "   preflight OK (stub)"
SH
chmod +x "$tmp/proj/workflow/preflight.sh"

# A fixture stepX.nf that writes reports/run_status.txt exactly the way step2.nf's
# onComplete does, driven by a param so both outcomes are exercised.
cat > "$tmp/proj/stepfail.nf" <<'NF'
nextflow.enable.dsl=2
params.verdict = 'ok'
process NOOP {
    script:
    """
    true
    """
}
workflow.onComplete {
    def f = file("reports/run_status.txt")
    f.parent.mkdirs()
    f.text = params.verdict + "\n"
}
workflow { NOOP() }
NF

cd "$tmp/proj"

# 1. verdict=fail -> phylohpc run must exit non-zero and reports/run_status.txt says fail.
cat > "$tmp/proj/fail.yaml" <<'Y'
verdict: fail
Y
# run_nf hardcodes "$HERE/$step.nf" -- point step at our fixture by naming it stepfail.
sed -i 's/step2/stepfail/' phylohpc 2>/dev/null || true
out_fail=$(./phylohpc run stepfail -p fail.yaml 2>&1)
rc_fail=$?
echo "$out_fail" | tail -10
if [[ $rc_fail -eq 0 ]]; then
    echo "  FAIL: phylohpc run exited 0 for verdict=fail"
    _FAILED=1
else
    echo "  ok (rc=$rc_fail): phylohpc run exited non-zero for verdict=fail"
fi
assert_ok test -f reports/run_status.txt
assert_eq "fail" "$(tr -d '[:space:]' < reports/run_status.txt)" "run_status.txt says fail"

# 2. verdict=ok -> phylohpc run must exit 0.
rm -f reports/run_status.txt
cat > "$tmp/proj/ok.yaml" <<'Y'
verdict: ok
Y
out_ok=$(./phylohpc run stepfail -p ok.yaml -w work2 2>&1)
rc_ok=$?
if [[ $rc_ok -eq 0 ]]; then
    echo "  ok: phylohpc run exited 0 for verdict=ok"
else
    echo "  FAIL: phylohpc run exited $rc_ok for verdict=ok"
    _FAILED=1
fi

# 3. nextflow ITSELF exits 0 (a real successful run) but never writes run_status.txt
# (no onComplete write at all -- the "nextflow died before reaching onComplete" case,
# and more generally "the file is simply missing"). This must ALSO be treated as
# failure, never as a silent success by omission -- exercises check_run_status's own
# absence branch specifically, as opposed to test 1's explicit verdict=fail.
rm -f reports/run_status.txt
cat > "$tmp/proj/crash.nf" <<'NF'
nextflow.enable.dsl=2
process NOOP2 {
    script:
    """
    true
    """
}
// deliberately no onComplete -- reports/run_status.txt is never written, even though
// this run succeeds from nextflow's own point of view.
workflow { NOOP2() }
NF
sed -i 's/stepfail/crash/' phylohpc
out_crash=$(./phylohpc run crash -p ok.yaml -w work3 2>&1)
rc_crash=$?
assert_contains "$out_crash" "did not reach onComplete" "phylohpc names the real cause"
if [[ $rc_crash -eq 0 ]]; then
    echo "  FAIL: phylohpc run exited 0 when run_status.txt was never written"
    _FAILED=1
else
    echo "  ok (rc=$rc_crash): phylohpc run exited non-zero when run_status.txt is absent"
fi
assert_eq "0" "$(test -f reports/run_status.txt && echo 1 || echo 0)" "run_status.txt correctly absent"

cd "$HERE" >/dev/null

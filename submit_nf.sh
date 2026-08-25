#!/usr/bin/env bash
#SBATCH --no-requeue
# 1GB was not enough for the driver: the 2026-08-24 run was OOM-killed at 20:02:34 with
# ~1500 families in flight ("Detected 1 oom_kill event in StepId=27863248.batch"), taking
# the whole workflow down mid-GeneRax.
#SBATCH --mem=8GB
#SBATCH -p genoa64
#SBATCH --qos=pipelines

set -euo pipefail

_term() {
    kill -s SIGTERM "$pid" 2>/dev/null || true
    wait "$pid" 2>/dev/null || true
}
trap _term TERM

module load Java

# Pin the engine. The 2026-08-16 run used 25.04.2; nextflow then self-updated to 26.04.6,
# which changed config validation and put the resume cache for 1455 ALN / 476 PHY / 473
# PVM_PREV tasks at risk. Override with NXF_VER=... on the command line if you mean to.
export NXF_VER="${NXF_VER:-25.04.2}"

# Keep the JVM inside the cgroup above, so a heap overrun surfaces as a Java OOM we can
# read rather than a silent SIGKILL of the whole driver.
export NXF_OPTS="${NXF_OPTS:--Xms1g -Xmx6g}"

# iqtree2 is NOT in the `tfs` conda env (which supplies generax, mpirun and mafft), so every
# PHY/ALN task in the 2026-08-24 run died instantly with "iqtree2 not found in PATH" -- 1308
# PHY and 159 ALN failures. No single env has all three; `phylo` carries iqtree2. Appended,
# not prepended, so tfs keeps priority for everything it does provide.
_IQTREE_ENV="${IQTREE_ENV:-$HOME/miniconda3/envs/phylo/bin}"
if [[ -x "$_IQTREE_ENV/iqtree2" ]]; then
    export PATH="$PATH:$_IQTREE_ENV"
else
    echo "WARNING: no iqtree2 at $_IQTREE_ENV -- PHY will fail" >&2
fi
command -v iqtree2 >/dev/null || echo "WARNING: iqtree2 still not on PATH" >&2

mkdir -p reports

NF_EXTRA_OPTS=()
PASSTHRU=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --report)
            # Optional filename argument
            if [[ $# -ge 2 && "$2" != -* ]]; then
                NF_EXTRA_OPTS+=("-with-report" "$2")
                shift 2
            else
                NF_EXTRA_OPTS+=("-with-report" "reports/report.html")
                shift 1
            fi
            ;;
        --trace)
            if [[ $# -ge 2 && "$2" != -* ]]; then
                NF_EXTRA_OPTS+=("-with-trace" "$2")
                shift 2
            else
                NF_EXTRA_OPTS+=("-with-trace" "reports/trace.txt")
                shift 1
            fi
            ;;
        --timeline)
            if [[ $# -ge 2 && "$2" != -* ]]; then
                NF_EXTRA_OPTS+=("-with-timeline" "$2")
                shift 2
            else
                NF_EXTRA_OPTS+=("-with-timeline" "reports/timeline.html")
                shift 1
            fi
            ;;
        *)
            PASSTHRU+=("$1")
            shift
            ;;
    esac
done

mkdir -p reports
nextflow run -ansi-log false \
    "${NF_EXTRA_OPTS[@]}" \
    "${PASSTHRU[@]}" & pid=$!

wait $pid
exit 0
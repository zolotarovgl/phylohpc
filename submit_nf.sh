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

# Pre-flight the environment. The 2026-08-24 run was launched from the `tfs` conda env,
# which has generax + mpirun but NO iqtree2 -- so 1308 PHY and 159 ALN tasks died within
# milliseconds on "iqtree2 not found in PATH" while GeneRax ran happily, and the workflow
# reported nothing wrong. Refuse to start rather than burn a night on that again.
#
# The intended environment (workflow/environment.yaml) is:
#     conda activate phylo     # iqtree2, generax, mafft, hmmer, diamond, mcl, POSSVM deps
#     module load OpenMPI      # phylo's openmpi record is the conda-forge EXTERNAL stub
#                              # (0 files) -- mpirun is expected to come from the cluster
_missing=()
for _t in iqtree2 mafft generax mpirun; do
    command -v "$_t" >/dev/null 2>&1 || _missing+=("$_t")
done
if (( ${#_missing[@]} )); then
    echo "ERROR: not on PATH: ${_missing[*]}" >&2
    echo "       run 'conda activate phylo' and 'module load OpenMPI' before submitting;" >&2
    echo "       see workflow/environment.yaml. Set SKIP_ENV_CHECK=1 to override." >&2
    [[ -n "${SKIP_ENV_CHECK:-}" ]] || exit 78
fi
echo "env: iqtree2=$(command -v iqtree2) generax=$(command -v generax) mpirun=$(command -v mpirun)"

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
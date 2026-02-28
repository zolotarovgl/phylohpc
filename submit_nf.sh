#!/usr/bin/env bash
#SBATCH --no-requeue
#SBATCH --mem=1GB
#SBATCH -p genoa64
#SBATCH --qos=pipelines

set -euo pipefail

_term() {
    kill -s SIGTERM "$pid" 2>/dev/null || true
    wait "$pid" 2>/dev/null || true
}
trap _term TERM

module load Java

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
nextflow run -resume -ansi-log false \
    "${NF_EXTRA_OPTS[@]}" \
    "${PASSTHRU[@]}" & pid=$!

wait $pid
exit 0
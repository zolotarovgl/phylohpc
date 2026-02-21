#!/usr/bin/env bash
#SBATCH --no-requeue
#SBATCH --mem=6G
#SBATCH -p genoa64
#SBATCH --qos=pipelines

set -e
set -u
set -o pipefail

_term() {
	kill -s SIGTERM $pid
	wait $pid
}
trap _term TERM

module load Java

export NXF_JVM_ARGS="-Xms2g -Xmx5g"

nextflow run -ansi-log false "$@" & pid=$!
wait $pid
exit 0

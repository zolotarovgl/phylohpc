#!/usr/bin/env bash
#SBATCH --no-requeue
#SBATCH --mem=1GB
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

#export NXF_JVM_ARGS="-Xms2g -Xmx5g"
#export NXF_DEBUG=2
#export NXF_OPTS='-Dnextflow.trace.flush=true'

nextflow run -resume -ansi-log false "$@" \
	-with-report report.html \
	-with-trace trace.txt \
	-with-timeline timeline.html & pid=$!
	
wait $pid
exit 0

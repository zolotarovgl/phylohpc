#!/usr/bin/env bash
if [ $# -ne 2 ]; then
    echo "usage: $0 config.txt output.tsv" >&2
    exit 1
fi
CONFIG=$1
OUTFILE=$2
source workflow/functions.sh
read_config "$CONFIG"
> "$OUTFILE"
for PREF in $(ls ${CLUSTER_DIR}/*_cluster.tsv | xargs -n1 basename | sed 's/_cluster.tsv//g'); do
    awk -v PREF="$PREF" '{print $2"\t"PREF"."$1}' ${CLUSTER_DIR}/${PREF}_cluster.tsv >> "$OUTFILE"
done


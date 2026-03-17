#!/bin/bash
set -e

check_file() {
    local file="$1"
    if [[ -f "$file" && -s "$file" ]]; then
        STATUS=1
    else
        STATUS=0
    fi
    echo $STATUS
}

source workflow/functions.sh
CONFIG=configs/config.txt

read_config "$CONFIG"

if [[ $# -eq 1 ]]; then
    IFS="." read -r PREF FAMILY HG <<< "$1"
elif [[ $# -eq 2 ]]; then
    PREF=$1
    FAMILY=$2
else
    echo -e "Check the status of the PREF FAMILY\nUsage: $0 PREF FAMILY\n   or: $0 PREF.FAMILY"
    exit 1
fi


FILEPREF=${PREF}.${FAMILY}
STATUS_S=$(check_file $SEARCH_DIR/${PREF}.${FAMILY}.genes.list)
STATUS_C=$(check_file $CLUSTER_DIR/${PREF}.${FAMILY}_cluster.tsv)

echo -e "#prefix\tsearch\tclustering"
echo -e "$FILEPREF\t$STATUS_S\t$STATUS_C"

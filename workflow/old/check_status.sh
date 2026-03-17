#!/bin/bash
#set -u -x
set -e
trap 'echo "Undefined variable at line $LINENO: $BASH_COMMAND"' ERR

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
elif [[ $# -eq 3 ]]; then
    PREF=$1
    FAMILY=$2
    HG=$3
else
    echo -e "Check the status of the PREF FAMILY HG\nUsage: $0 PREF FAMILY HG\n   or: $0 PREF.FAMILY.HG"
    exit 1
fi

FILEPREF=${PREF}.${FAMILY}.${HG}
CLU_FILE=${CLUSTER_DIR}/${FILEPREF}.fasta
ALN_FILE=$ALIGN_DIR/${FILEPREF}.aln.fasta
TREE_FILE=$TREE_DIR/${FILEPREF}.treefile
POSSVM_FILE=$TREE_DIR/${FILEPREF}.treefile.ortholog_groups.csv

STATUS_C=$(check_file "$CLU_FILE")
STATUS_A=$(check_file "$ALN_FILE")
STATUS_T=$(check_file "$TREE_FILE")
STATUS_P=$(check_file "$POSSVM_FILE")
# check if the phylogeny has finished 

if [[ -f "${TREE_DIR}/${FILEPREF}.log" ]]; then
    STATUS_T=$(grep -c "Consensus tree written" "${TREE_DIR}/${FILEPREF}.log" || true)
else
    STATUS_T=0
fi

echo -e "#prefix\tcluster\talign\ttree\tpossvm"
echo -e "$FILEPREF\t$STATUS_C\t$STATUS_A\t$STATUS_T\t$STATUS_P"

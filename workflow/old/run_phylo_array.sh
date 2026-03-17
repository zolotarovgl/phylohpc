#!/bin/bash
# Load configuration file
CONFIG_FILE="configs/config.txt"

if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Configuration file $CONFIG_FILE not found!"
    exit 1
fi

# Read variables from the configuration file
while IFS='=' read -r key value; do
    if [[ ! -z "$key" && ! "$key" =~ ^# ]]; then  # Ignore empty lines and comments
        eval "${key}='${value}'"
    fi
done < "$CONFIG_FILE"

###################################
# Parse provided arguments - need to replace with globally available ones!
###################################

if [ "$#" -lt 7 ]; then
    echo 'Provide all 5 arguments!'
    echo -e "Array script\nUsage: $0 <INPUT_LIST> [OUTDIR] [REFNAMES] [REFSPECIES] [NCPU] [TREE_METHOD] [MAFFT_OPT]"
    exit 1
fi


INPUT_LIST=$1
OUTDIR=$2
REFNAMES=$3
REFSPECIES=$4
NCPU=$5
TREE_METHOD=$6
MAFFT_OPT=$7

# get file prefix from prefix file
id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INPUT_LIST")
echo $id
# this is the prefix for which the phylogeny should be ran
# here $id contains ${PREFIX} ${FAMILY} ${CLUSTER}
PREFIX=$(echo $id | cut -f 1 -d .)
FAMILY=$(echo $id | cut -f 2 -d .)
HG=$(echo $id | cut -f 3 -d .)
echo "${PREFIX}.${FAMILY}.${HG}"

#############################################
bash workflow/run_phylo.sh ${CLUSTER_DIR}/${PREFIX}.${FAMILY}.${HG}.fasta ${PREFIX}.${FAMILY}.${HG} $OUTDIR $REFNAMES $REFSPECIES $NCPU $TREE_METHOD "$MAFFT_OPT"


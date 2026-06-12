#!/bin/bash
set -euo pipefail
source workflow/functions.sh
CONFIG_FILE="configs/config.txt"
SCRIPT_NAME="run_generax_array.sh"
read_config $CONFIG_FILE

if [ "$#" -lt 1 ]; then
    echo 'Provide arguments!'
    echo -e "Array script\nUsage: $0 <INPUT_LIST> "
    exit 1
fi

####################################################
INPUT_LIST=$1
echo "${SCRIPT_NAME}: INPUT_LIST: ${INPUT_LIST}"
# For each job in array:
ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INPUT_LIST")
# here $id contains ${PREF} ${FAMILY} ${CLUSTER}
#####################################################
# Array input
echo $ID
PREF=$(echo $ID | cut -f 1 -d .)
FAMILY=$(echo $ID| cut -f 2 -d .)
HG=$(echo $ID | cut -f 3 -d .)
echo "${SCRIPT_NAME}: ${SLURM_ARRAY_TASK_ID}: ${PREF}.${FAMILY}.${HG}"
# Inputs
ALIGNMENT=${ALIGN_DIR}/${PREF}.${FAMILY}.${HG}.aln.fasta
GENETREE=${TREE_DIR}/${PREF}.${FAMILY}.${HG}.treefile
IQTREEFILE=${TREE_DIR}/${PREF}.${FAMILY}.${HG}.iqtree
echo "${SCRIPT_NAME}: INPUTS: ${ALIGNMENT} ${GENETREE} ${IQTREEFILE}" 

#######################################################
check_file_exists $ALIGNMENT
check_file_exists $GENETREE
check_file_exists $IQTREEFILE
check_file_exists $SPECIES_TREE
#########################################################
# Outputs
mkdir -p $GENERAX_DIR
OUTDIR=$GENERAX_DIR/${PREF}.${FAMILY}.${HG}
echo "bash workflow/run_generax.sh --alignment $ALIGNMENT --genetree $GENETREE --iqtreefile $IQTREEFILE --speciestree $SPECIES_TREE --prefix ${PREF}.${FAMILY}.${HG} --maxspr $GENERAX_MAXSPR --outdir $OUTDIR --ncpu $NCPU --tmpdir $GENERAX_TMP"
bash workflow/run_generax.sh --alignment $ALIGNMENT --genetree $GENETREE --iqtreefile $IQTREEFILE --speciestree $SPECIES_TREE --prefix ${PREF}.${FAMILY}.${HG} --maxspr $GENERAX_MAXSPR --outdir $GENERAX_DIR --ncpu $NCPU --tmpdir $GENERAX_TMP

#!/bin/bash
if [ "$#" -lt 8 ]; then
    echo -e "Run simple phylogeny\nUsage: $0 <INPUT_FASTA> <OUT_PREFIX> <OUTDIR> <REFNAMES> <REFSPECIES> <NCPU> <TREE_METHOD> <MAFFT_OPT>"
    exit 1
fi

# Required arguments
INPUT_FASTA=$1
OUT_PREFIX=$2
OUTDIR=$3
REFNAMES=$4
REFSPECIES=$5
NCPU=$6
TREE_METHOD=$7
MAFFT_OPT=$8


# Define directories and files
ALN_DIR=${OUTDIR}/align
TREE_DIR=${OUTDIR}/gene_trees
ALN_FILE=${ALN_DIR}/${OUT_PREFIX}.aln.fasta
TREE_PREFIX=${TREE_DIR}/${OUT_PREFIX}
TREE_FILE=${TREE_PREFIX}.treefile
POSSVM_OUT=${TREE_FILE}.ortholog_groups.csv

############################
# Logging prefix:
LOGGING="$0: ${OUT_PREFIX}:"
# Print configuration
echo "Configuration:"
echo "OUTDIR: $OUTDIR"
echo "ALN_DIR: $ALN_DIR"
echo "TREE_DIR: $TREE_DIR"
echo "REFNAMES: $REFNAMES"
echo "NCPU: $NCPU"
echo "TREE_METHOD: $TREE_METHOD"
echo "MAFFT_OPT: $MAFFT_OPT"


echo "INPUT_FASTA: $INPUT_FASTA"
echo "OUT_PREFIX: $OUT_PREFIX"
echo "ALN_FILE: $ALN_FILE"
echo "TREE_FILE: $TREE_FILE"
echo "POSSVM_OUT: $POSSVM_OUT"
########################

if [ -f "$POSSVM_OUT" ]; then
    echo "Skipping POSSVM: $POSSVM already exists."
else
if [ -f "$TREE_FILE" ]; then
    echo "Skipping phylogeny: $TREE_FILE already exists."
else
if [ -f "$ALN_FILE" ]; then
    echo "Skipping alignment: $ALN_FILE already exists."
else
    mkdir -p $ALN_DIR
    echo "${LOGGING} Running alignment..."
    echo "${LOGGING} python phylogeny/main.py align -f $INPUT_FASTA -o $ALN_FILE -c $NCPU -m ${MAFFT_OPT}"
    python phylogeny/main.py align -f $INPUT_FASTA -o $ALN_FILE -c $NCPU -m "$MAFFT_OPT"
    echo "Alignment done."
fi
    echo "Running phylogeny..."
    mkdir -p $TREE_DIR
    echo "${LOGGING} python phylogeny/main.py phylogeny -f $ALN_FILE --outprefix $TREE_PREFIX -c $NCPU --method $TREE_METHOD"
    python phylogeny/main.py phylogeny -f $ALN_FILE --outprefix $TREE_PREFIX -c $NCPU --method $TREE_METHOD
    [ ! -f "$TREE_FILE" ] && echo "File not found: $TREE_FILE" && exit 1    
    echo "Tree done."
fi
    echo "${LOGGING} python phylogeny/main.py possvm -t $TREE_FILE --refsps $REFSPECIES -r $REFNAMES -o ${OUT_PREFIX}"
    python phylogeny/main.py possvm -t $TREE_FILE --refsps $REFSPECIES -r $REFNAMES -o ${OUT_PREFIX}"." 
fi
echo "${LOGGING} Phylogeny done."

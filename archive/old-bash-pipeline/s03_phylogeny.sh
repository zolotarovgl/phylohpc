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

# Argument parsing
if [ "$#" -ne 5 ]; then
    echo "Command: $0 $@"
    echo -e "Launch 2 arrays for big and small groups.\nUsage: $0 <PREFIX; e.g. tfs> <FAMILY> <MIN_SEQ> <MAX_SEQ> <NCPU>"
    exit 1
fi

PREFIX=$1
FAMILY=$2
MIN_SEQ=$3 # minimal allowed number of sequences
MAX_SEQ=$4 # maximum number of sequences to send to the long job
NCPU=$5

# Create output directories
mkdir -p $ALIGN_DIR
mkdir -p $TREE_DIR

# Your script logic here
# Example: echo commands for testing
echo "Launching jobs with the following parameters:"
echo "Prefix: $PREFIX, Family: $FAMILY, Min Seq: $MIN_SEQ, Max Seq: $MAX_SEQ, NCPU: $NCPU"

#######################################################################################
# RUN
#######################################################################################

##################
# 1. Filter clusters
##################
echo "Cluster filtering"
mkdir -p tmp/clustering
cat ${CLUSTER_DIR}/${PREFIX}.${FAMILY}_cluster.tsv | cut -f 1 | sort | uniq -c | awk -v MAX_SEQ=${MAX_SEQ} '$1>=MAX_SEQ {print $2}' > tmp/clustering/${PREFIX}.${FAMILY}.clusters_big
cat ${CLUSTER_DIR}/${PREFIX}.${FAMILY}_cluster.tsv | cut -f 1 | sort | uniq -c | awk -v MAX_SEQ=${MAX_SEQ} -v MIN_SEQ=${MIN_SEQ} '$1<MAX_SEQ&&$1>MIN_SEQ {print $2}' > tmp/clustering/${PREFIX}.${FAMILY}.clusters_small


# now, prepare fasta files for phylogenies 
for ID in $(cat tmp/clustering/${PREFIX}.${FAMILY}.clusters_{small,big}); do
if [ ! -e "${SEARCH_DIR}/${PREFIX}.${FAMILY}.domains.fasta.fai" ]; then
samtools faidx ${SEARCH_DIR}/${PREFIX}.${FAMILY}.domains.fasta
fi
xargs samtools faidx ${SEARCH_DIR}/${PREFIX}.${FAMILY}.domains.fasta < <(awk -v ID=$ID '$1==ID {print $2}' ${CLUSTER_DIR}/${PREFIX}.${FAMILY}_cluster.tsv) > results/clusters/${PREFIX}.${FAMILY}.${ID}.fasta
done

echo "done."
####################
# 2. Launch arrays
####################
# This script is the main center of handling big and small groups and adjusting the parameters
mkdir -p logs/arrays
# create the list of prefixes for array execution => launch an array job 
if [ -s tmp/clustering/${PREFIX}.${FAMILY}.clusters_big ]; then
# create prefix file 
awk -v PREF=$PREFIX -v FAMILY=$FAMILY '{print PREF"."FAMILY"."$1}' tmp/clustering/${PREFIX}.${FAMILY}.clusters_big > tmp/clustering/${PREFIX}.${FAMILY}.clusters_big.prefs 
# run job array on a file prefix 
N=$(wc -l < tmp/clustering/${PREFIX}.${FAMILY}.clusters_big.prefs)
echo "Submitting ${N} big jobs for ${PREFIX}.${FAMILY}"
# Usage: run_phylo_array.sh <INPUT_LIST> [OUTDIR] [REFNAMES] [REFSPECIES] [NCPU] [TREE_METHOD] 
sbatch --job-name=${PROJECT}.s3.${PREFIX}.${FAMILY}.big \
       --output=logs/arrays/s3.${PREFIX}.${FAMILY}.big_%A_%a.out \
       --error=logs/arrays/s3.${PREFIX}.${FAMILY}.big_%A_%a.err \
       --array=1-$N \
       --qos $QOS_S3_LONG \
       --cpus-per-task=$NCPU \
       --time=${TIME_S3_LONG} \
       --mem=${MEM_LONG} \
       workflow/run_phylo_array.sh tmp/clustering/${PREFIX}.${FAMILY}.clusters_big.prefs results/ $REFNAMES $REFSPECIES $NCPU $TREE_METHOD_BIG "$MAFFT_OPT_BIG"
fi
if [ -s tmp/clustering/${PREFIX}.${FAMILY}.clusters_small ]; then
# create prefix file 
awk -v PREF=$PREFIX -v FAMILY=$FAMILY '{print PREF"."FAMILY"."$1}' tmp/clustering/${PREFIX}.${FAMILY}.clusters_small > tmp/clustering/${PREFIX}.${FAMILY}.clusters_small.prefs 
# run job array on a file prefix 
N=$(wc -l < tmp/clustering/${PREFIX}.${FAMILY}.clusters_small.prefs)
echo "Submitting ${N} small jobs for ${PREFIX}.${FAMILY}"
sbatch --job-name=${PROJECT}.s3.${PREFIX}.${FAMILY}.small \
       --output=logs/arrays/s3.${PREFIX}.${FAMILY}.small_%A_%a.out \
       --error=logs/arrays/s3.${PREFIX}.${FAMILY}.small_%A_%a.err \
       --array=1-$N \
       --qos $QOS_S3_SHORT \
       --cpus-per-task=$NCPU \
       --time=${TIME_S3_SHORT} \
       --mem=${MEM_SHORT} \
       workflow/run_phylo_array.sh tmp/clustering/${PREFIX}.${FAMILY}.clusters_small.prefs results/ $REFNAMES $REFSPECIES $NCPU $TREE_METHOD_SMALL "$MAFFT_OPT_SMALL"
fi
echo "$0 done."
exit 0

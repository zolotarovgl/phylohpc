#!/bin/bash
set -e

# Ensure all arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <INPUT_FASTA> <OUTPUT_FILE> <MAX_N>"
    exit 1
fi

INPUT_FASTA=$1
OUTPUT_FILE=$2
MAX_N=$3

# Infer the output prefix for the per-cluster fasta files 


echo "$(date '+%Y-%m-%d %H:%M:%S') $0 start"
source workflow/functions.sh
read_config configs/config.txt


echo "Clustering: ${INPUT_FASTA} => ${OUTPUT_FILE}"
echo "python phylogeny/main.py cluster -f $INPUT_FASTA --out_file $OUTPUT_FILE -c $S2_NCPU -m $MAX_N"
python phylogeny/main.py cluster -f $INPUT_FASTA --out_file $OUTPUT_FILE -c $S2_NCPU -m $MAX_N -i $S2_INFLATION
###############################
echo "$(date '+%Y-%m-%d %H:%M:%S') $0 splitting fasta..."
samtools faidx $INPUT_FASTA
OUTPREFIX=$(echo $OUTPUT_FILE | sed -E 's/_cluster.tsv$//g')
for ID in $(cut -f 1 $OUTPUT_FILE | sort | uniq); do 
xargs samtools faidx $INPUT_FASTA < <(awk -v ID=${ID} '$1==ID { print $2}' $OUTPUT_FILE) > ${OUTPREFIX}.${ID}.fasta
echo ${OUTPREFIX}.${ID}.fasta
done 
###############################
echo "$(date '+%Y-%m-%d %H:%M:%S') $0 done"
exit 0  

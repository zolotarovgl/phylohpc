#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --output=/users/asebe/gzolotarov/projects/2024_hpc_phylogeny/pipeline_new/logs/%x.%j.out
#SBATCH --error=/users/asebe/gzolotarov/projects/2024_hpc_phylogeny/pipeline_new/logs/%x.%j.err

#SBATCH --nodes=1 #always 1
#SBATCH --ntasks=1 #always 1
#SBATCH --cpus-per-task=1 #change depending on your needs


# Ensure all arguments are provided
if [ "$#" -ne 3 ]; then
    echo -e "Run simple POSSVM\nUsage: $0 <PREFIX> <FAMILY> <REFERENCE_NAMES>"
    exit 1
fi


PREFIX=$1
FAMILY=$2
REFNAMES=$3

for TREEFILE in $(ls results/gene_trees/${PREFIX}.${FAMILY}.HG*.treefile); do 
echo $TREEFILE
if [ ! -e "$TREEFILE.ortholog_groups.newick" ]; then
python phylogeny/main.py possvm -t $TREEFILE -r $REFNAMES
fi
done

echo "$0 done."
exit 0

#!/bin/bash
set -e

source workflow/functions.sh

if [[ $# -lt 2 ]]; then
    echo -e "Submit hmmsearch and clustering for a family\nUsage: $0 CONFIG PREF FAMILY\n   or: $0 CONFIG PREF.FAMILY"
    exit 1
fi

CONFIG=$1
read_config "$CONFIG"
shift

if [[ $# -eq 1 ]]; then
    IFS="." read -r PREF FAMILY <<< "$1"
elif [[ $# -eq 2 ]]; then
    PREF=$1
    FAMILY=$2
else
    echo -e "Usage: $0 CONFIG PREF FAMILY\n   or: $0 CONFIG PREF.FAMILY"
    exit 1
fi

FILEPREF=${PREF}.${FAMILY}
echo -e "Prefix: ${PREF}; Family: ${FAMILY}"


# s01 HMM search 
job_name=${PROJECT}.s1.${PREF}.${FAMILY}
job_out=$(sbatch --time=${TIME_S1} --qos=${QOS_S1} --mem=${MEM_S1} --job-name=$job_name --output=$LOG_DIR/s1/%x.%j.out --error=$LOG_DIR/s1/%x.%j.out workflow/s01_search.sh ${FAMILY}  ${GENEFAM_INFO} ${INFASTA} $SEARCH_DIR)
jid=$(echo $job_out | awk '{print $4}')
jid1=$jid
echo "Submitted ${jid}: ${job_name}"
# s02 Clustering 
job_name=${PROJECT}.s2.${PREF}.${FAMILY}
job_out=$(sbatch --qos=${QOS_S2} --time=${TIME_S2} --dependency=afterok:$jid --mem=${MEM_S2} --job-name=$job_name --output=$LOG_DIR/s2/%x.%j.out --error=$LOG_DIR/s2/%x.%j.out workflow/s02_cluster.sh ${SEARCH_DIR}/${PREF}.${FAMILY}.domains.fasta ${CLUSTER_DIR}/${PREF}.${FAMILY}_cluster.tsv $MAX_N)
jid=$(echo $job_out | awk '{print $4}')
echo "Submitted ${jid}: ${job_name}"

jid2=$jid

#scancel $jid1
#scancel $jid2

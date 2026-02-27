source workflow/functions.sh
read_config configs/config.txt

mkdir -p $SEARCH_DIR
mkdir -p $CLUSTER_DIR


PREF=tfs
FAMILY=Forkhead

# s01 HMM search 
JOB_NAME=${PROJECT}.s1.search.${PREF}.${FAMILY}
LOG_OUTPUT=$LOG_DIR/s1/%x.%j.out
LOG_ERROR=$LOG_DIR/s1/%x.%j.out
job_out=$(sbatch \
    --time=${TIME_S1} \
    --qos=${QOS_S1} \
    --mem=${MEM_S1} \
    --job-name=${JOB_NAME} \
    --output=${LOG_OUTPUT} \
    --error=${LOG_ERROR} \
    workflow/s01_search.sh \
    ${FAMILY} ${GENEFAM_INFO} ${INFASTA} ${SEARCH_DIR})
jid=$(echo "$job_out" | awk '{print $4}')
echo "Submitted ${jid}: ${JOB_NAME}"

# s02 Clustering
JOB_NAME=${PROJECT}.s2.cluster.${PREF}.${FAMILY}
LOG_OUTPUT=$LOG_DIR/s2/%x.%j.out
LOG_ERROR=$LOG_DIR/s2/%x.%j.out
job_out=$(sbatch \
    --qos=${QOS_S2} \
    --time=${TIME_S2} \
    --dependency=afterok:${jid} \
    --mem=${MEM_S2} \
    --job-name=${JOB_NAME} \
    --output="${LOG_OUTPUT}" \
    --error="${LOG_ERROR}" \
    workflow/s02_cluster.sh \
    "${SEARCH_DIR}/${PREF}.${FAMILY}.domains.fasta" \
    "${CLUSTER_DIR}/${PREF}.${FAMILY}_cluster.tsv" \
    ${MAX_N})
jid=$(echo "$job_out" | awk '{print $4}')
echo "Submitted ${jid}: ${JOB_NAME}"



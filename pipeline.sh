#rm -rf logs
#rm -rf results
#rm -rf tmp/


#################################################

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


printf "%-20s: %s\n" "NCPU" "$NCPU"
printf "%-20s: %s\n" "QOS_S1" "$QOS_S1"
printf "%-20s: %s\n" "QOS_S2" "$QOS_S2"
printf "%-20s: %s\n" "MEM_S1" "$MEM_S1"
printf "%-20s: %s\n" "MEM_S2" "$MEM_S2"
printf "%-20s: %s\n" "TIME_S1" "$TIME_S1"
printf "%-20s: %s\n" "TIME_S2" "$TIME_S2"
printf "%-20s: %s\n" "TIME_S3" "$TIME_S3"
printf "%-20s: %s\n" "TIME_S3_SHORT" "$TIME_S3_SHORT"
printf "%-20s: %s\n" "TIME_S3_LONG" "$TIME_S3_LONG"

printf "%-20s: %s\n" "MIN_SEQ" "$MIN_SEQ"
printf "%-20s: %s\n" "MAX_SEQ" "$MAX_SEQ"
printf "%-20s: %s\n" "MAX_N" "$MAX_N"
printf "%-20s: %s\n" "GENEFAM_INFO" "$GENEFAM_INFO"
printf "%-20s: %s\n" "REFNAMES" "$REFNAMES"
printf "%-20s: %s\n" "INFASTA" "$INFASTA"
printf "%-20s: %s\n" "PROJECT" "$PROJECT"

echo -e "\nDirectories:"
printf "%-20s: %s\n" "OUTDIR" "$OUTDIR"
printf "%-20s: %s\n" "LOG_DIR" "$LOG_DIR"
printf "%-20s: %s\n" "SEARCH_DIR" "$SEARCH_DIR"
printf "%-20s: %s\n" "CLUSTER_DIR" "$CLUSTER_DIR"
printf "%-20s: %s\n" "TREE_DIR" "$TREE_DIR"
printf "%-20s: %s\n" "ALIGN_DIR" "$ALIGN_DIR"


################################################
mkdir -p $OUTDIR
mkdir -p $LOG_DIR
mkdir -p $SEARCH_DIR
mkdir -p $CLUSTER_DIR
mkdir -p $TREE_DIR
mkdir -p $ALIGN_DIR
echo -e '\n\n\n'


for L in $(grep -v '#' $GENEFAM_INFO  | awk '{print $NF"."$1}'); do

    PREF=$(echo $L | cut -f 1 -d .)
    FAMILY=$(echo $L | cut -f 2 -d .)
    echo "${PREF}-${FAMILY}"
    
    OUTFILE=${SEARCH_DIR}/${PREF}.${FAMILY}.domains.csv
    
    if [ ! -f $OUTFILE ]; then
        # s01 HMM search 
        job_out=$(sbatch --time=${TIME_S1} --qos=${QOS_S1} --mem=${MEM_S1} --job-name=${PROJECT}s1.${PREF}.${FAMILY} --output=$LOG_DIR/s1/%x.%j.out --error=$LOG_DIR/s1/%x.%j.out workflow/s01_search.sh ${FAMILY}  ${GENEFAM_INFO} ${INFASTA} $SEARCH_DIR)
        jid=$(echo $job_out | awk '{print $4}')
        echo "Submitted ${jid}: ${PROJECT}s1.${PREF}.${FAMILY}"
        # s02 Clustering 
        job_out=$(sbatch --qos=${QOS_S2} --time=${TIME_S2} --dependency=afterok:$jid --mem=${MEM_S2} --job-name=${PROJECT}s2.${PREF}.${FAMILY} --output=$LOG_DIR/s2/%x.%j.out --error=$LOG_DIR/s2/%x.%j.out workflow/s02_cluster.sh ${SEARCH_DIR}/${PREF}.${FAMILY}.domains.fasta ${CLUSTER_DIR}/${PREF}.${FAMILY}_cluster.tsv $MAX_N)
        jid=$(echo $job_out | awk '{print $4}')
        echo "Submitted ${jid}: ${PROJECT}s2.${PREF}.${FAMILY}"
        # launch a job (launcher) that will launch 2 arrays
        job_out=$(sbatch \
            --time=${TIME_S3} \
            --dependency=afterok:$jid \
            --qos=vshort \
            --job-name=${PROJECT}.s3.launcher.${PREF}.${FAMILY} \
            --nodes=1 \
            --ntasks=1 \
            --cpus-per-task=1\
            --output=$LOG_DIR/s3/%x.%j.out \
            --error=$LOG_DIR/s3/%x.%j.err \
            workflow/s03_phylogeny.sh $PREF $FAMILY $MIN_SEQ $MAX_SEQ $NCPU)
        jid=$(echo $job_out | awk '{print $4}')
        echo "Submitted ${jid}: ${PROJECT}s3.${PREF}.${FAMILY}"
    else
        echo "Found ${OUTFILE}. Skipping ...."
    fi   
done

# Parhyale TFs phylogenies

```bash
# Collect genomic taxonomic info
cat /users/asebe/xgraubove/genomes/genome_sources_2025-03-20.tsv | grep -f species_list | awk -F '\t' '{print $2}' > list
# get crustacean tree
cat ~/genomes/crustacea/genome_info.tab  | cut -f 2 >> list
echo "Parhyale hawaiensis" >> list
~/species_tree/ncbi_tree.py --input list > species_tree.newick
~/species_tree/tax_info.py --input list > tax_info.tab
```

```
13:39:06.665 [ERRO] taxonomy data not found, please download and uncompress ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz, and copy "names.dmp", "nodes.dmp", "delnodes.dmp", and "merged.dmp" to /users/asebe/gzolotarov/.taxonkit
```

Prepare the input
```bash
bash workflow/prepare_fasta.sh species_list data/input.fasta
# add Parhaw- prefixed TFs 
cat ../data/List_all_TFs_PFAM.fa | bioawk -c fastx '{print ">Parhaw_"$name"\n"$seq}' > Parhaw_TFs.fasta
#cat Parhaw_TFs.fasta >> data/input.fasta
bioawk -c fastx '{print ">Parhaw_"$name"\n"$seq}' ../data/Parhyale_hawaiensis_longest_prot_2_short_name.fasta > Parhaw_all.fasta
cat ~/genomes/crustacea/data/*_long.pep.fasta >> data/input.fasta
cat Parhaw_all.fasta >> data/input.fasta
samtools faidx data/input.fasta
```


# Summarize 

Compare the list of input TFs to the list of the proteins picked up by the pipeline

```bash
cat Parhaw_TFs.fasta  | bioawk -c fastx '{print $name}' | sort | uniq > ids.tfs
N_PICK=$(cat results/search/*domains.csv  | grep -w -f ids.tfs | wc -l | awk '{print $1}')
N_INPUT=$(wc -l ids.tfs | awk '{print $1}')
cat results/search/*list | grep Parhaw | sort | uniq > ids.pick
N=$(wc -l ids.pick | awk '{print $1}') 
N_MISSING=$(comm ids.tfs ids.pick -23 | wc -l | awk '{print $1}')
N_NEW=$(comm ids.tfs ids.pick -13  | wc -l | awk '{print $1}')
echo -e "# picked by the pipeline: $N\n# from input: ${N_PICK} / ${N_INPUT}\n# missing: $N_MISSING\n# new $N_NEW"

```


```bash
# get pep2cluster mapping 
./workflow/get_gene2cluster.sh configs/config.txt pep2cluster.tab
```



# Add Mmus possvm run:
```bash
cat hg_status.tab | grep 1111 | cut -f 1 > tmp/torun
cat results/gene_trees/*groups.csv  | grep Parhaw | cut -f 2 | sort | uniq  | cut -f 1 -d : | sed -E 's/\.[0-9]+$//g' > tmp/torun
for ID in $(cat tmp/torun); do 
PREF=${ID%%.*}
FAMILY=${ID#*.};FAMILY=${FAMILY%%.*}
HG=${ID##*.}
echo $PREF $FAMILY $HG
TREE_FILE=${TREE_DIR}/${PREF}.${FAMILY}.${HG}.treefile
mkdir -p tmp/trees/
cp $TREE_FILE tmp/trees/
REFSPECIES=Dmel
REFNAMES=data/${REFSPECIES}_gene_names.csv
python phylogeny/main.py possvm -t tmp/trees/${ID}.treefile --refsps $REFSPECIES -r $REFNAMES -o ${PREF}.${FAMILY}.${HG}"."
mv tmp/trees/${ID}.treefile.ortholog_groups.csv tmp/trees/${ID}.${REFSPECIES}.treefile.ortholog_groups.csv
mv tmp/trees/${ID}.treefile.ortholog_groups.newick tmp/trees/${ID}.${REFSPECIES}.treefile.ortholog_groups.newick
mv tmp/trees/${ID}.treefile.pairs_orthologs.csv tmp/trees/${ID}.${REFSPECIES}.treefile.pairs_orthologs.csv
REFSPECIES=Mmus
REFNAMES=data/${REFSPECIES}_gene_names.csv
python phylogeny/main.py possvm -t tmp/trees/${ID}.treefile --refsps $REFSPECIES -r $REFNAMES -o ${PREF}.${FAMILY}.${HG}"."
mv tmp/trees/${ID}.treefile.ortholog_groups.csv tmp/trees/${ID}.${REFSPECIES}.treefile.ortholog_groups.csv
mv tmp/trees/${ID}.treefile.ortholog_groups.newick tmp/trees/${ID}.${REFSPECIES}.treefile.ortholog_groups.newick
mv tmp/trees/${ID}.treefile.pairs_orthologs.csv tmp/trees/${ID}.${REFSPECIES}.treefile.pairs_orthologs.csv
done 

```






# Useful commands

Check the status of the HGs:
```bash
for FILEPREF in $(basename -a -s .fasta results/clusters/*.fasta); do  ./workflow/check_status.sh $FILEPREF; done | grep -v '#' | awk '{print $1"\t"$2$3$4$5}' > hg_status.tab
```


```bash
# Check the status of a single HG
#./workflow/check_status.sh tfs HMGbox_Sox HG1
get_hg_job_status() {
    local ID="$1"
    local line
    line=$(list_all_jobs | awk -v id="$ID" '$2 ~ (id "$") {print $0}' | tail -n 1)
    if [[ -n "$line" ]]; then
        #echo $line
      echo "$line" | awk '{print $4}'
    else
        echo "NOTFOUND"
    fi
}

# submit missing alignment jobs
#python submit.py configs/config.txt --cpus 8 --mem 1G --mode align --pref tfs --family Homeodomains --hg HG13 -n 

# Filter HGs by size - get the list of all decent HGs
grep -c '>' results/clusters/*.fasta \
  | awk -F ':' '$2>=5{print $1}' \
  | xargs -n1 basename -s .fasta > tmp/all_ids

for FILEPREF in $(cat tmp/all_ids); do  ./workflow/check_status.sh $FILEPREF; done | grep -v '#' | awk '{print $1"\t"$2$3$4$5}' > hg_status.tab

for ID in $(cat tmp/all_ids); do JOB_STATUS=$(get_hg_job_status $ID); echo -e $ID"\t"$JOB_STATUS; done > job_status.tab
awk 'FNR==NR{d[$1]=$2;next}{print $1"\t"$2"\t"d[$1]}' job_status.tab hg_status.tab | grep -v -E 'RUNNING|PENDING' > hg_status.filt.tab
cat hg_status.filt.tab | cut -f 2 | sort | uniq -c | sort -rn

######################################################
# Submit alignment for 1000
######################################################
cat hg_status.filt.tab | grep 1000 | cut -f 1 > tmp/torun



while IFS="." read -r PREF FAMILY HG; do
    echo "$PREF $FAMILY $HG"
    python submit.py configs/config.txt --cpus 8 --mem 3G --mode align --pref $PREF --family $FAMILY --hg $HG 
done < tmp/torun


# submit if NOT RUNNING already...
# issue 
######################################################
# Submit phylogeny for 1100
######################################################
cat hg_status.filt.tab | grep 1100 | cut -f 1 > tmp/torun

while IFS="." read -r PREF FAMILY HG; do
    echo "$PREF $FAMILY $HG"
    python submit.py configs/config.txt --mode phylogeny --cpus 8 --mem 1G --time "12:00:00"  --pref $PREF --family $FAMILY --hg $HG 
done < tmp/torun
######################################################
# Possum for 1110
# The reference species is Dmel everywhere 
######################################################

read_config configs/config.txt
for ID in $(./get_hg_status.sh  | grep 1110 | cut -f 1); do 
PREF=${ID%%.*}
FAMILY=${ID#*.};FAMILY=${FAMILY%%.*}
HG=${ID##*.}
echo $PREF $FAMILY $HG
TREE_FILE=${TREE_DIR}/${PREF}.${FAMILY}.${HG}.treefile
python phylogeny/main.py possvm -t $TREE_FILE --refsps $REFSPECIES -r $REFNAMES -o ${PREF}.${FAMILY}.${HG}"."
done 


#####################################################
# Resubmit alignment and phylogeny job
#####################################################

./get_hg_status.sh  | grep 1110

PREF=tfs
FAMILY=bZIP
HG=HG4
read_config configs/config.txt
TREE_FILE=${TREE_DIR}/${PREF}.${FAMILY}.${HG}.treefile
python phylogeny/main.py possvm -t $TREE_FILE --refsps $REFSPECIES -r $REFNAMES -o ${PREF}.${FAMILY}.${HG}"."




# launch a job (launcher) that will launch 2 arrays
job_out=$(sbatch \
    --time=${TIME_S3} \
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


```


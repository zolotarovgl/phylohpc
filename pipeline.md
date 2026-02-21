

## Step1 - prepare the inputs   


```bash
bash workflow/prepare_fasta.sh species_list_small data/input.fasta
```



```bash  
mkdir -p info # a folder to store JSON info files 
python check_families.py configs/config.txt --genefam genefam.csv --json info/families.json  --output family_status.tab 
# Now run with --resubmit to re-submit the missing families 
python check_families.py configs/config.txt --genefam genefam.csv --json info/families.json  --output family_status.tab --resubmit
```

```bash
python submit_hg.py --pref tfs --family Forkhead --hg HG1  --mode align configs/config.txt
```



### Nextflow pipeline 

 - [ ] run on the correct node 

```
	DO NOT EXECUTE NEXTFLOW PIPELINES ON LOGIN NODES
	FOLLOW INSTRUCTIONS HERE
	https://dokuwiki.hpc.crg.es/doku.php?id=sit:nextflow_on_new_cluster

```



Nextflow pipeline: 
`ids.txt` - contains the list of the `PREF.FAMILY.HG` to run the analysis for


```bash
 grep -l '>Clacla' results/clusters/*fasta | xargs -n1 basename | sed 's/.fasta//g' | sort | uniq > ids.txt 
```

```bash
module load Java
```

```bash
sbatch submit_nf.sh main.nf -profile slurm -w /no_backup/asebe/gzolotarov/work/
```



To clean: 

```bash
nextflow clean -failed -f
```



---

For example, to submit alignment and phylogeny jobs:

```bash
# Select HGs to run run only for HGs with Clacla
grep -l '>Clacla' results/clusters/*fasta | xargs -n1 basename | sed 's/.fasta//g' | sort | uniq > ids
for ID in $(cat ids); do ./workflow/check_status.sh $ID; done | grep -v '#' | awk '{print $1"\t"$2$3$4$5}' > hg_status.tab
cat hg_status.tab  | grep 1000 | cut -f 1  > tmp/torun


# submit alignment jobs:
TIME=00:30:00
MEM=500M
NCPU=4
while IFS="." read -r PREF FAMILY HG; do
    python submit_hg.py configs/config.txt --time $TIME --cpus $NCPU --mem $MEM --mode align --pref $PREF --family $FAMILY --hg $HG --mafft ""  --json aln.info.json
done < tmp/torun

./workflow/get_hg_status.sh > hg_status.tab
cat hg_status.tab | grep 1100 | cut -f 1  | grep -w -f ids > tmp/torun

# remove gap-only sequences
for ID in $(cat tmp/torun); do
echo $ID
clean_fasta results/align/${ID}.aln.fasta tmp_fa 20; mv tmp_fa results/align/${ID}.aln.fasta
done

TIME=1-00:00:00
MEM=500M
NCPU=6
while IFS="." read -r PREF FAMILY HG; do
    python submit_hg.py configs/config.txt --time $TIME --cpus $NCPU --mem $MEM --mode phylogeny --pref $PREF --family $FAMILY --hg $HG  --json phy.info.json
done < tmp/torun
```
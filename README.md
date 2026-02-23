# TODOs

- [ ] `step1` - search and clustering pipeline   
- [ ] allow the phylogeny script to rerun the iqtree if it finds the outputs? 
- [ ] limit the iqtree run times 
- [ ] generax + 2nd possvm 
- [ ] mafft oom errors (code 1 instesad of 137)


# Pipeline 

0. prepare the input data - `data/input.fasta`    
1. search and clustering $\rightarrow$ homology groups. Filtering $\rightarrow$ `ids.txt`  
2. Alignment and phylogenies per homology group   
3. Gather the annotations per species  


## step 0

```bash
bash workflow/prepare_fasta.sh species_list data/input.fasta
```

## DEV step 1

Run the pythong wrapper for the step1 

```bash
mkdir -p info # a folder to store JSON info files 
python check_families.py configs/config.txt --genefam genefam.csv --json info/families.json  --output family_status.tab --resubmit
mkdir -p tmp
grep -l '>Clacla' results/clusters/*fasta | xargs -n1 basename | sed 's/.fasta//g' > tmp/hg_soi
for f in results/clusters/*fasta; do      printf "%s\t%s\n" "$(basename "$f" | sed 's/.fasta//g')" "$(grep -c '>' "$f")"; done > tmp/hg_count
cat tmp/hg_count | awk 'FNR==NR{d[$1]=$2;next}d[$1]>=30{print $1}' - tmp/hg_soi > ids.txt
```

## step 2

Once you have step1 results and the list of homology groups to process `ids.txt`, run the pipeline: 

```bash
mamba activate phylo 
WORKDIR=/no_backup/asebe/gzolotarov/work/
sbatch submit_nf.sh step2.nf -profile slurm -w $WORKDIR
```


## step3   

Gather the annotations per species of interest:  

```bash
# no splitting by the group 
mkdir -p results/annotations
for SP in Clacla Corcan Osclob Axidam Halduj Spolac; do
	echo $SP
	python workflow/gather_anno.py --id ${SP} --search-dir results/search/ --tree-dir results/possvm/ > results/annotations/${SP}.tab
done

# with splittinb by the functional groups 
for SP in Clacla Corcan Osclob Axidam Halduj Spolac; do
for PREF in tfs neu ion; do
	echo "${SP} ${PREF}"
	python workflow/gather_anno.py --prefix ${PREF} --id ${SP} --search-dir results/search/ --tree-dir results/possvm/ > results/annotations/${SP}.${PREF}.tab
done
done
# still very few proteins classified 
```




# Resource usage 

`dowstream_stats.R` - explores and plots resource usage. 


# TODOs

- [ ] `step1` - search and clustering pipeline   
- [ ] job duration and memory prediction  
- [ ] proper job re-submission rules   
- [ ] allow the phylogeny script to rerun the iqtree if it finds the outputs? 
- [ ] limit the iqtree run times 
- [ ] generax + 2nd possvm 
- [ ] mafft oom errors (code 1 instesad of 137)
- [x] `PHY` job time extension 

Testing PHY job re-submission for timeouts.
The problem, is that during each execution, the nextflow will start over - one needs to be able to predict how long the jobs will take. 


# Pipeline 

0. prepare the input data - `data/input.fasta`    
1. search and clustering $\rightarrow$ homology groups. Select which HGs to run the pipeline for $\rightarrow$ `ids.txt`  
2. Alignment and phylogenies per homology group   
3. Gather the annotations per species  


# Prepare input data 

```bash
bash workflow/prepare_fasta.sh species_list data/input.fasta
```

# Step1 

## Interactive  
```bash
# Interactive session
module load Java 
mamba activate phylo 
WORKDIR=work_step1
nextflow run -profile local -w $WORKDIR -resume step1.nf --genefam_info genefam.csv --infasta data/input.fasta 
```


## SLURM  
```bash
# SLURM submssions
WORKDIR=work_step1
mkdir -p $WORKDIR
mkdir -p reports
sbatch --time=01:00:00 -J step1 submit_nf.sh step1.nf -profile slurm -w $WORKDIR --report reports/report.step1.html --trace reports/trace.step1.txt --timeline reports/timeline.step1.html 
```

Note: Use `-profile slurm`  to run using the SLURM scheduler instead of locally. Use interactive jobs if `-profile local` unless you want Emyr coming to your desk!  

__Note:__ the re-clustering of the fig files does not work here!   


# Filter homology groups   

Filter and get the list of homology groups to run the alignment for: 
```bash
python select_hgs.py --out ids.txt --soi Mmus --min_seqs 10 --min_sps 3
```
`--soi` - will keep only HG ids with `Mmus` sequnce  

# Step 2

## Resource-informed submission:

Predict resources to generate `resources.tsv`. All IDs not in this table will get default values.


```bash
Rscript predict_resources.R \
  --ids_fn ids.txt \
  --cluster_dir results/clusters \
  --models_rds workflow/models/models.rds \
  --outfile resources.tsv \
  --min_mem 100 --min_time 5 --max_mem 10000 --max_time 2880 --increase 0.1 
```

```bash
WORKDIR=work_step2
sbatch -J step2 submit_nf.sh step2.nf -profile slurm -w $WORKDIR --report reports/report.step2.html --trace reports/trace.step2.txt --timeline reports/timeline.step2.html 
```


# Gather annotations 

Gather the annotations per species of interest:  

```bash
# no splitting by the group 
SP=Nvec
mkdir -p results/annotations 
python workflow/gather_anno.py --id ${SP} --search-dir results/search/ --tree-dir results/possvm/ > results/annotations/${SP}.tab

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


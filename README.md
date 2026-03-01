# Pipeline 

0. prepare the input data - `data/input.fasta`    
1. search and clustering $\rightarrow$ homology groups. Select which HGs to run the pipeline for $\rightarrow$ `ids.txt`  
2. Alignment, phylogeny, possvm     
- optional: GeneRax + POSSVM
3. Gather the annotations per species  


# Prepare input data 
From a species list, get the proteomes from Xavi's database. You can as well just use custom proteomes concatenated into `data/input.fasta`  
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
nextflow run -profile local -w $WORKDIR -resume step1.nf --genefam_info genefam.csv --infasta data/input.fasta -with-report reports/report.step1.html -with-trace reports/trace.step1.html
```


## SLURM  
```bash
# SLURM submssions
WORKDIR=work_step1
mkdir -p $WORKDIR
mkdir -p reports
sbatch --time=01:00:00 -J step1 submit_nf.sh step1.nf -profile slurm -w $WORKDIR --report reports/report.step1.html --trace reports/trace.step1.txt --timeline reports/timeline.step1.html 
```

__Note__: Use `-profile slurm`  to run using the SLURM scheduler instead of locally. Use interactive jobs if `-profile local` unless you want Emyr coming to your desk!  


## Homology groups filtering  

Filter and get the list of homology groups to run the alignment for: 
```bash
python select_hgs.py --out ids.txt --soi Mmus --min_seqs 5 --min_sps 3
```
* `--soi` - will keep only HG ids with `Mmus` sequences (it makes sense to filter by the reference species)  

Output: `ids.txt` file with selected homology groups. 

# Step 2

Predict resources to generate `resources.tsv`. All IDs not in this table will get default values.

```bash
# models.json should contain the mem and time models for each job stored as coefficients - see _export_models.r which will expor workflow/models/models.rds
python workflow/predict_resources.py --ids_fn ids.txt --cluster_dir results/clusters --models_json workflow/models/models.json --outfile resources.tsv --min_mem 100 --min_time 5 --max_mem 10000 --max_time 2880 --increase 0.1 
```

Submit with predicted resources   
```bash
# Interactive session
WORKDIR=work_step2
nextflow run -profile local -w $WORKDIR -resume step2.nf --genefam_info genefam.csv --infasta data/input.fasta -with-report reports/report.step2.html -with-trace reports/trace.step2.html

# SLURM
WORKDIR=work_step2
sbatch -J step2 submit_nf.sh step2.nf -profile slurm -w $WORKDIR --report reports/report.step2.html --trace reports/trace.step2.txt --timeline reports/timeline.step2.html 
```

# Step 3 - GeneRax   


```bash

# list HGs with trees:
find results/gene_trees -type f -name "*.treefile" \
  -exec basename {} .treefile \; > ids_generax


# INTERACTIVE
module load OpenMPI
module load Java 
mamba activate phylo
WORKDIR=work_generax
IDS=ids_generax
nextflow run -profile local -w $WORKDIR -resume generax.nf --ids $IDS 
# SLURM
WORKDIR=work_generax
mkdir -p $WORKDIR
IDS=ids_generax
sbatch -J generax submit_nf.sh generax.nf --ids $IDS -profile slurm -w $WORKDIR --report reports/report.generax.html --trace reports/trace.generax.txt --timeline reports/timeline.generax.html 
```




----
----
----

# Dev

# Gather annotations 

After the phylogenies are done, you can gather the annotations per species 
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

# TODOs

- [ ] proper environment with `openmpi` for generax  
- [ ] generax resource prediction   
- [ ] quantile regression for resource prediction   
- [ ] `phylo` environment with `Rscript` support  
- [ ] allow the phylogeny script to rerun the iqtree if it finds the outputs? 
- [ ] mafft oom errors (code 1 instesad of 137) - proper handling 
- [x] better subclustering logic in `phylogeny/`  
- [x] generax family error handling - raises exit 10  
- [x] Clustering: proper subclustering - local and global model  
- [x] generax
- [x] `step1` - search and clustering pipeline   
- [x] `PHY` job time extension 
- [x] job duration and memory prediction  
- [x] proper job re-submission rules   

# Pipeline 

0. prepare the input data - `data/input.fasta`    
1. search and clustering $\rightarrow$ homology groups. Select which HGs to run the pipeline for $\rightarrow$ `ids.txt`  
2. Alignment, phylogeny, possvm     
- optional: GeneRax + POSSVM
3. Gather the annotations per species  


# Profiles   
There are 2 main flavours to run the pipeline:  
- `fast` - for testing, automatic `mafft`, `fasttree` and `max_spr` of the generax is set to 2 (suboptimal) 
- `precise` - LINSI mode, IQTREE2 with model testing, max_spr up to 7   

Execution:  
- `local` - use the interactive HPC session or run on your machine   
- `slurm` - configured to be used on the CRG HPC system   

To run a pipeline, combine the flavour and the executor. For instance `-profile local,fast` will run the fast pipeline locally. For submitting the jobs via slurm, you have to use a species sbatch script `submit.nf` (vis the commands below).  


# Prepare inputs   
If you are planning to run GeneRax, you have to make sure i) your species tree contains all the prefixes present in the input fasta file; ii) strictly binary (expected by GeneRax) - i.e. no polytomies are present. To check the tree:   

```bash
python workflow/check_tree.py  data/species_tree.full.newick species_list data/species_tree.newick 
```

From a species list, get the proteomes from Xavi's database. You can as well just use custom proteomes concatenated into `data/input.fasta`  
```bash
bash workflow/prepare_fasta.sh species_list data/input.fasta
```



# Step1 

Interactive: 

```bash
module load Java 
mamba activate phylo 


WORKDIR=/no_backup/asebe/gzolotarov/nextflow/phylohpc/work_step1
nextflow run -profile local -w $WORKDIR -resume step1.nf --genefam_info genefam.csv --infasta data/input.fasta -with-report reports/report.step1.html -with-trace reports/trace.step1.html
```

SLURM:  
```bash
module load Java 
mamba activate phylo 

WORKDIR=/no_backup/asebe/gzolotarov/nextflow/phylohpc/work_step1
sbatch --time=01:00:00 -J step1 submit_nf.sh step1.nf -profile slurm -w $WORKDIR --report reports/report.step1.html --trace reports/trace.step1.txt --timeline reports/timeline.step1.html 
```

__Note__: Use `-profile slurm`  to run using the SLURM scheduler instead of locally. Use interactive jobs if `-profile local` unless you want Emyr coming to your desk!  


## Homology groups filtering  

Filter and get the list of homology groups to run the following steps for: 
```bash
python select_hgs.py --out ids.txt --soi Mmus --min_seqs 30 --min_sps 10
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
module load OpenMPI
module load Java 
mamba activate phylo
WORKDIR=/no_backup/asebe/gzolotarov/nextflow/phylohpc/work_step2
PROFILE=local,fast
nextflow run -resume -profile $PROFILE -w $WORKDIR  step2.nf --run_generax --genefam_info genefam.csv --infasta data/input.fasta -with-report reports/report.step2.html -with-trace reports/trace.step2.html

# SLURM
PROFILE=slurm,fast
sbatch -J step2 submit_nf.sh step2.nf -profile $PROFILE --run_generax -w $WORKDIR --report reports/report.step2.html --trace reports/trace.step2.txt --timeline reports/timeline.step2.html 
```

`--run_generax` - use this flag to run `GeneRax` prior to `POSSVM`. 


# Runtimes   

How long does the slurm execution take with the linsi, fasttree and small SPR value? 

# Gather annotations 

After the phylogenies are done, you can gather the annotations per species 
Gather the annotations per species of interest:  

```bash
# no splitting by the group 
SP=Nvec
TREEDIR=results/possvm/ # use possvm if no generax available, or results/generax 
python workflow/gather_annotations.py --search-dir results/search/ --tree-dir $TREEDIR --id Nvec
```




# Resource usage 

`dowstream_stats.R` - explores and plots resource usage. 

---

# TODOs

- [ ] 2 execution profiles - fast and precise
- [ ] generax: missing species in the tree - tree checks!  
- [ ] generax speedup - does increasing the number of cores make a difference?  
- [ ] generax resource prediction   
- [ ] quantile regression for resource prediction   
- [ ] `phylo` environment with `Rscript` support  
- [ ] allow the phylogeny script to rerun the iqtree if it finds the outputs? 
- [ ] mafft oom errors (code 1 instesad of 137) - proper handling 
- [x] `generax.nf` - proper OOM and OOT handling  
- [x] proper environment with `openmpi` for generax  
- [x] `step2.nf` - make sure the processes are correctly cached and not rerun   
- [x] generax caching issue   
- [x] better subclustering logic in `phylogeny/`  
- [x] generax family error handling - raises exit 10  
- [x] Clustering: proper subclustering - local and global model  
- [x] generax
- [x] `step1` - search and clustering pipeline   
- [x] `PHY` job time extension 
- [x] job duration and memory prediction  
- [x] proper job re-submission rules   

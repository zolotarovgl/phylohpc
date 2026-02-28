Nextflow pipeline for the step1 

TODOs:
- [ ] separate sbatch script to submit the step1 

# Prepare input data 

```bash
bash workflow/prepare_fasta.sh species_list data/input.fasta
```

# Step1 

```bash
# Interactive session
module load Java 
mamba activate phylo 
nextflow run -profile local -resume step1.nf --genefam_info genefam.csv --infasta data/input.fasta 
```

```bash
# SLURM submssions
WORKDIR=work_step1
mkdir -p $WORKDIR
mkdir -p reports
sbatch --time=01:00:00 -J step1 submit_nf.sh step1.nf -profile slurm -w $WORKDIR --report reports/report.step1.html --trace reports/trace.step1.txt --timeline reports/timeline.step1.html 
```

Note: Use `-profile slurm`  to run using the SLURM scheduler instead of locally. Use interactive jobs if `-profile local` unless you want Emyr coming to your desk!  


# Filter homology groups   

Filter and get the list of homology groups to run the alignment for: 
```bash
python select_hgs.py --out ids.txt --soi Mmus --min_seqs 10 --min_sps 3
```
`--soi` - will keep only HG ids with `Mmus` sequnce  
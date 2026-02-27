Nextflow pipeline for the step1 

```bash
bash workflow/prepare_fasta.sh species_list data/input.fasta
```

```bash
module load Java 
mamba activate phylo 
nextflow run -profile local -resume step1.nf --genefam_info genefam.csv --infasta data/input.fasta 
```

Filter and get the list of homology groups to run the alignment for: 

```bash
t.b.a
```
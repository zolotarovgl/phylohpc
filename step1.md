Nextflow pipeline for the step1 


```bash
module load Java 
mamaba activate phylo 
nextflow run -resume step1.nf --genefam_info genefam.csv --infasta data/input.fasta 
```

Filter and get the list of homology groups to run the alignment for: 

```bash
t.b.a
```
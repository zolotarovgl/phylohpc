# TODOs

- [ ] allow the phylogeny script to rerun the iqtree if it finds the outputs? 


# Pipeline 

## step 0

```bash
bash workflow/prepare_fasta.sh species_list data/input.fasta
```

## step 1

```bash
mkdir -p info # a folder to store JSON info files 
python check_families.py configs/config.txt --genefam genefam.csv --json info/families.json  --output family_status.tab --resubmit

mkdir -p tmp
grep -l '>Clacla' results/clusters/*fasta | xargs -n1 basename | sed 's/.fasta//g' > tmp/hg_soi
for f in results/clusters/*fasta; do      printf "%s\t%s\n" "$(basename "$f" | sed 's/.fasta//g')" "$(grep -c '>' "$f")"; done > tmp/hg_count
cat tmp/hg_count | awk 'FNR==NR{d[$1]=$2;next}d[$1]>=30{print $1}' - tmp/hg_soi > ids.txt
```

## step 2

```bash
mamba activate phylo 
WORKDIR=/no_backup/asebe/gzolotarov/work/
sbatch submit_nf.sh step2.nf -profile slurm -w $WORKDIR
```



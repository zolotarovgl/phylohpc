# PhyloHPC — Gene Family Phylogenomics Pipeline

A Nextflow pipeline for large-scale gene family phylogenomics: from raw proteomes to ortholog annotations, gene tree–species tree reconciliation, and clade-level ancestral gene content reconstruction.

![](img/metazoan_tree.TFevol.svg)

---


## Overview

```
Multi-species proteomes
        │
        ▼
  step1.nf ── PFAM domain search (HMMER) - **genefam.csv**
              MCL clustering
              → Homology Groups (HGs)
        │
        ▼
  step2.nf ── Multiple sequence alignment (MAFFT)
              Gene tree inference (IQ-TREE2 / FastTree)
              Ortholog prediction (POSSVM)
              [optional] Gene tree–species tree reconciliation (GeneRax)
        │
        ▼
  workflow/gather_annotations.py ── Per-species ortholog annotation tables
        │
        ▼
  workflow/step4.ancestry.nf ── Clade-scoped ortholog calling (POSSVM + --ignoretips)
                        Presence/absence matrix
                        Ancestral state reconstruction (PastML, MPPA + F81)
                        Interactive visualisation (D3.js HTML)
```

**Documentation:** full user manual → [`docs/manual.md`](docs/manual.md) · Step 4 ancestral reconstruction → [`docs/ancestry.md`](docs/ancestry.md)

---

## Clone 


```bash
git clone --recurse-submodules https://github.com/zolotarovgl/phylohpc.git
```


## Requirements

**Environment:** `workflow/environment.yaml` is a minimal, portable spec — only
the direct tools, loosely pinned, so it solves across machines. For local macOS
work on Apple Silicon, use the Rosetta/x86_64 environment in
`workflow/environment.macos-x86_64.yaml`.

```bash
mamba env create -f workflow/environment.yaml
```


```bash
mamba activate phylo
```

On Apple Silicon Macs:

```bash
conda create -n phylo-macos
conda activate phylo-macos
conda config --env --set subdir osx-64
mamba env update -n phylo-macos -f workflow/environment.macos-x86_64.yaml
```

If `mamba` is not installed, use:

```bash
conda env update -n phylo-macos -f workflow/environment.macos-x86_64.yaml
```

The Nextflow config sets `PYTHONNOUSERSITE=1` for local and slurm runs so
user-level `~/.local` Python packages do not override the environment.
The macOS environment uses Python 3.10 because the current `clipkit` build
requires it for the alignment trimming step in `step2.nf`.

Key tools: HMMER, MAFFT, IQ-TREE2, FastTree, MCL, GeneRax, POSSVM, PastML, ete3, Python 3.10.

**On the CRG HPC**, load the required modules before running any pipeline:
```bash
mamba activate phylo
module load Java 
module load OpenMPI       
```

**GeneRax reconciliation (`--run_generax`) is optional and not in the conda env.**
It needs an MPI-compiled GeneRax + `mpirun` (the bioconda build is non-MPI and
crashes). Use it on the CRG HPC via `module load OpenMPI` with an MPI GeneRax, or
run locally without `--run_generax` (POSSVM runs on the raw trees). POSSVM ships
with the `phylogeny/` submodule cloned above.

## Input data

| File | Description |
|---|---|
| `data/input.fasta` | Concatenated proteomes for all species. Sequence IDs must start with a 3–6 character species prefix followed by `_` (e.g. `Mmus_ENSMUSP00000001`). |
| `data/species_tree.full.newick` | Full species tree with **named internal nodes** (required for GeneRax and ancestral reconstruction). |
| `data/species_tree.newick` | Strictly binary species tree (no polytomies) used by GeneRax. Generated from the full tree — see below. |
| `data/Mmus_gene_names.csv` | Reference species gene name mapping (`protein_id<TAB>gene_name`). |
| `genefam.csv` | Gene family definitions (HMM profiles, clustering parameters — see `docs/Gene_families.md`). |
| `config/species_list` | One species prefix per line; defines the analysis set. |



## Test data

### Fetch proteomes from a sequence database
```bash
DB_PATH=~/ant/xgraubove/genomes/data/ # path to Xavi's database 
bash workflow/prepare_fasta.sh config/species_list data/input.fasta $DB_PATH
```

### Prepare the binary species tree

Note: `--random-resolve` is for the cases where the tree contains polytomies
```bash
python workflow/check_tree.py --random-resolve data/species_tree.full.newick config/species_list data/species_tree.newick
```

This validates the full tree, prunes it to the species in `config/species_list`, and resolves any polytomies, writing a GeneRax-compatible binary tree to `data/species_tree.newick`.


### Test genefam.csv 

```bash
cat data/gene_families_searchinfo.csv | grep -E 'Forkhead|T-box' > genefam.csv
```




---

## Execution profiles

Profiles are combined on the command line with `-profile <executor>,<flavour>`.

**Executors**

| Profile | Description |
|---|---|
| `local` | Run interactively on the current machine / HPC interactive session. |
| `slurm` | Submit each process as a SLURM job on the CRG HPC (`genoa64` queue). |

**Flavours**

| Profile | Alignment | Tree method | GeneRax SPR |
|---|---|---|---|
| `fast` | MAFFT default | FastTree | 2 |
| `precise` | MAFFT L-INS-i | IQ-TREE2 (model test) | 3 |

Use `fast` for testing and small families; `precise` for publication-quality trees.

---

## Step 1 — PFAM search and homology group clustering

Searches the input proteome with HMMER against the defined PFAM profiles, then clusters domain-containing sequences into Homology Groups (HGs) using MCL.

```bash
WORKDIR=/no_backup/asebe/gzolotarov/nextflow/phylohpc/work_step1/
WORKDIR=work
PROFILE=local,fast

# NB: do not forget to set the path to pfam_db Pfam-A.hmm file!
# Interactive
nextflow run step1.nf \
    -profile $PROFILE \
    -resume \
    -w $WORKDIR \
	--pfam_db /home/grygoriyzolotarov/ant/xgraubove/data/pfam/Pfam-A.hmm \
    --genefam_info genefam.csv \
    --infasta data/input.fasta \
    -with-report reports/report.step1.html \
    -with-trace  reports/trace.step1.txt

# SLURM
PROFILE=slurm,fast
sbatch --time=01:00:00 -J step1 submit_nf.sh step1.nf \
    -resume -profile $PROFILE \
    -w $WORKDIR \
    --report  reports/report.step1.html \
    --trace   reports/trace.step1.txt \
    --timeline reports/timeline.step1.html
```

**Outputs:** `results/search/` (domain sequences, gene lists) and `results/clusters/` (per-HG FASTA files).


### Select HGs for downstream analysis

```bash
# Inspect cluster sizes
python workflow/get_seqstat.py results/clusters/*.fasta | sort -k2 -n

# Select HGs that contain the reference species, have ≥20 sequences
# across ≥5 species, and write their IDs to ids.txt
python workflow/select_hgs.py \
    --out ids.txt \
    --soi Mmus \
    --min_seqs 20 \
    --min_sps 5
```

`ids.txt` controls which HGs are processed in all subsequent steps.

---

## Step 2 — Gene tree inference and ortholog prediction

For each HG in `ids.txt`: align sequences (MAFFT), infer a gene tree (FastTree or IQ-TREE2), and predict ortholog groups with POSSVM.

Optionally, reconcile gene trees with the species tree using GeneRax before POSSVM.

### Optional: resource prediction

Before submitting large runs, predict memory and runtime requirements from a previous trace file:

```bash
python workflow/get_seqstat.py results/clusters/*.fasta > seq_stat.tab

Rscript workflow/train.R \
    --trace     reports/trace.step2.txt \
    --seq_stats seq_stat.tab \
    --outfile   workflow/models/models.json \
    --plotfile  workflow/models/models.pdf \
    --tau 0.95

python workflow/predict_resources.py \
    --ids_fn      ids.txt \
    --cluster_dir results/clusters \
    --models_json workflow/models/models.json \
    --defaults_json workflow/models/defaults.json \
    --outfile     resources.tsv \
    --max_mem 100000 --max_time 2880 --increase 0.5
```

### Run

```bash
# Interactive
WORKDIR=/no_backup/asebe/gzolotarov/nextflow/phylohpc/work_step2/
WORKDIR=work/
PROFILE=local,fast
nextflow run step2.nf \
    -profile $PROFILE \
    -resume -w $WORKDIR \
    --run_generax \
    --family_info genefam.csv \
    --species_tree data/species_tree.full.newick \
    -with-report reports/report.step2.html \
    -with-trace  reports/trace.step2.txt

# SLURM
PROFILE=slurm,precise
sbatch -J step2 -o reports/slurm.step2.out submit_nf.sh step2.nf \
    -profile $PROFILE \
    -resume -w $WORKDIR \
    --run_generax \
    --family_info genefam.csv \
    --species_tree data/species_tree.full.newick \
    --report   reports/report.step2.html \
    --trace    reports/trace.step2.txt \
    --timeline reports/timeline.step2.html
```

`--run_generax` enables gene tree–species tree reconciliation with GeneRax before POSSVM runs. Without this flag, POSSVM runs directly on the IQ-TREE2/FastTree trees.

**Outputs:**
- `results/align/` — trimmed alignments
- `results/gene_trees/` — raw gene trees (`.treefile`)
- `results/possvm/` — POSSVM ortholog groups (when `--run_generax`)
- `results/possvm_prev/` — POSSVM on raw trees (when `--run_generax`)
- `results/generax/` — reconciled gene trees (when `--run_generax`)
- `results/report_step2.html` - an interactive html report 
---


## Create an interactive report 

The script now takes as input results dir and figures out much of the rest
```bash
python workflow/report_step2.py --results_dir results --species_tree data/species_tree.full.newick --family_info genefam.csv --output report2.html
```

---

## Monitoring resource usage

```bash
# Get SLURM stats for completed jobs from a trace file
cat reports/trace.step2.txt | grep -E "COMPLETED|CACHED" | cut -f3 | grep -v native > job_ids
python workflow/check_job.v2.py -f job_ids > job_stats.tab
```

`R/downstream_stats.R` and `R/resources.R` plot CPU/memory efficiency across jobs. `R/generax_stats.R` analyses GeneRax runtime scaling with SPR radius.



## Misc

Download phylopics
```bash
comm <(cut -f 1 data/species_info.tsv | sort | uniq ) <(basename -a img/phylo/*.png | sed 's/.png//g' | sort | uniq) -32 | awk -F '\t' 'FNR==NR{d[$1]=$2;next}{print $1"\t"d[$1]}' data/species_info.tsv  - > tmp/phylopic_missing
wc -l tmp/phylopic_missing


for PREF in $(cat tmp/phylopic_missing | cut -f 1 ); do
SPECIES_NAME=$(awk -F '\t' -v ID=$PREF '$1==ID{print $2}' data/species_info.tsv )
echo $PREF
echo $SPECIES_NAME
python workflow/download_phylopic.py "${SPECIES_NAME}" -o img/phylo/${PREF}.png
done

```
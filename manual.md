# phylohpc — User Manual

A Nextflow pipeline for phylogenomic analysis on HPC clusters.
It takes a multi-species proteome and a gene-family definition file and produces:
- Sequence homology groups (HGs) from PFAM domain search + MCL clustering
- Multiple-sequence alignments (MAFFT)
- Gene trees (FastTree / IQ-TREE 2)
- Ortholog groups via POSSVM
- Reconciled gene trees via GeneRax *(optional)*
- Per-species ortholog annotations

---

## Table of Contents

1. [Installation](#installation)
2. [Repository layout](#repository-layout)
3. [Input files](#input-files)
4. [Quick-start](#quick-start)
5. [Step 1 — Search & clustering](#step-1--search--clustering-step1nf)
6. [Filtering homology groups](#filtering-homology-groups)
7. [Step 2 — Phylogenies](#step-2--phylogenies-step2nf)
   - [Resource prediction (optional)](#resource-prediction-optional)
   - [Running the pipeline](#running-the-pipeline)
8. [GeneRax-only pipeline](#generax-only-pipeline-generaxnf)
9. [Gathering annotations](#gathering-annotations)
10. [Monitoring jobs](#monitoring-jobs)
11. [Tree utilities](#tree-utilities)
12. [Post-processing / analysis](#post-processing--analysis)
13. [Parameter reference](#parameter-reference)
14. [Profiles reference](#profiles-reference)
15. [Output structure](#output-structure)

---

## Installation

### Prerequisites

| Tool | Version tested |
|---|---|
| Java | ≥ 11 |
| Nextflow | ≥ 23 |
| OpenMPI | any (required for GeneRax) |
| conda / mamba | any |

### Environment setup

```bash
mamba env create -f workflow/environment.yaml
mamba activate phylo
```

The `phylo` environment includes all bioinformatics tools used by the pipeline
(HMMER, DIAMOND, MCL, MAFFT, IQ-TREE 2, FastTree, GeneRax, POSSVM, ETE3,
NumPy, pandas, …).

---

## Repository layout

```
phylohpc/
├── step1.nf              # PFAM search + MCL clustering pipeline
├── step2.nf              # Alignment, phylogeny, POSSVM, GeneRax pipeline
├── generax.nf            # GeneRax-only re-run pipeline
├── nextflow.config       # Default parameters and execution profiles
├── submit_nf.sh          # SLURM wrapper for launching Nextflow itself
├── data/
│   ├── input.fasta           # Multi-species proteome (user-provided)
│   ├── species_tree.newick   # Pruned, strictly binary species tree
│   ├── species_tree.full.newick  # Full species tree (before pruning)
│   └── Mmus_gene_names.csv   # Reference species gene name table
├── configs/
│   └── config.txt        # Legacy config (no longer read by pipelines)
├── workflow/
│   ├── remove_gaponly.py     # Gap-only column removal (used by step2.nf)
│   ├── select_hgs.py         # Filter HGs → ids.txt
│   ├── get_seqstat.py        # FASTA sequence statistics
│   ├── predict_resources.py  # Predict SLURM memory/time per HG
│   ├── gather_annotations.py # Gather per-species ortholog annotations
│   ├── check_job.py          # Inspect individual SLURM job stats
│   ├── check_job.v2.py       # Inspect a batch of SLURM job stats
│   ├── check_tree.py         # Prune & validate species tree
│   ├── helper.py             # Shared Python utilities
│   ├── prune_tree.py         # Prune tree to a tip subset
│   ├── strip_len.py          # Strip branch lengths from Newick
│   ├── strip_support.py      # Strip bootstrap support from Newick
│   ├── get_tips.py           # List all tip labels in a tree
│   └── models/
│       ├── models.json       # Fitted resource-prediction coefficients
│       └── defaults.json     # Per-job default and large-family resources
├── phylogeny/
│   └── main.py               # Core computation engine (called by NF processes)
├── tests/                    # Pytest test suite
├── docs/
│   └── manual.md             # This file
└── *.R                       # Post-run analysis scripts
```

---

## Input files

### `data/input.fasta`

Concatenated protein sequences for all species. Headers must follow the
`<species_prefix>_<gene_id>` convention, e.g.:

```
>Mmus_ENSMUSG00000001
MKTIIALSYIFCLVFA...
>Hsap_ENSG00000139618
MPIGSKERPTFFEIFKTRCNKADLTHSQISDFHTYQITSFSTQNLQR...
```

All downstream tools derive the species prefix from everything before the
first `_` in the header.

### `genefam.csv` (or `genefam.csv` path via `--genefam_info`)

Tab-separated, seven columns, one row per gene family:

| Col | Field | Example |
|---|---|---|
| 1 | Family name | `TF_GATA` |
| 2 | HMM profile(s), comma-separated | `PF00320` |
| 3 | MCL inflation parameter | `1.4` |
| 4 | Minimum sequences per cluster | `4` |
| 5 | HMMER score threshold | `ga` (gathering cutoff) or `tc` |
| 6 | Group label | `TFs` |
| 7 | Output prefix | `tfs` |

### `data/species_tree.newick`

Strictly binary Newick species tree, with tip labels matching the species
prefixes used in `input.fasta`. See [Tree utilities](#tree-utilities) for
how to prepare this from a full reference tree.

---

## Quick-start

```bash
# 1. Activate environment and load modules
mamba activate phylo
module load Java OpenMPI

# 2. Run PFAM search & clustering
nextflow run step1.nf -profile local,fast \
  --genefam_info genefam.csv \
  --infasta data/input.fasta

# 3. Filter homology groups
python workflow/select_hgs.py --out ids.txt --soi Mmus --min_seqs 20 --min_sps 5

# 4. Run alignment, phylogeny, POSSVM (+ GeneRax)
nextflow run step2.nf -profile local,precise --run_generax

# 5. Gather annotations
python workflow/gather_annotations.py \
  --search-dir results/search \
  --tree-dir results/possvm results/possvm_prev \
  --id sps_annotate --outdir results/annotations --split-prefix
```

---

## Step 1 — Search & clustering (`step1.nf`)

Runs an HMMER profile search against the input proteome for each gene family,
then clusters the domain-extracted sequences with DIAMOND + MCL.

### Interactive (local) run

```bash
module load Java
mamba activate phylo

WORKDIR=/path/to/work/step1
nextflow run step1.nf \
  -profile local \
  -w $WORKDIR \
  -resume \
  --genefam_info genefam.csv \
  --infasta data/input.fasta \
  -with-report reports/report.step1.html \
  -with-trace  reports/trace.step1.txt
```

### SLURM batch run

```bash
WORKDIR=/path/to/work/step1
sbatch -J step1 -o reports/slurm.step1.out submit_nf.sh step1.nf \
  -profile slurm \
  -resume \
  -w $WORKDIR \
  --report   reports/report.step1.html \
  --trace    reports/trace.step1.txt \
  --timeline reports/timeline.step1.html
```

> **Important:** Always use `-profile slurm` when submitting via `sbatch`.
> Running without a SLURM profile in a batch job will try to execute all
> processes locally inside the submission node.

### Key parameters

| Parameter | Default | Description |
|---|---|---|
| `--genefam_info` | `genefam.csv` | Gene family definition file |
| `--infasta` | `data/input.fasta` | Multi-species protein FASTA |
| `--pfam_db` | *(Pfam-A.hmm path)* | Path to Pfam-A HMM database |
| `--domain_expand` | `30` | Residues added each side of HMM hit |
| `--s1_ncpu` | `4` | CPUs per hmmsearch job |
| `--s2_ncpu` | `4` | CPUs per clustering job |
| `--s2_inflation` | `1.1` | MCL inflation parameter |
| `--max_n` | `2000` | Max sequences; larger clusters are sub-clustered |
| `--search_dir` | `results/search` | Output directory for search results |
| `--cluster_dir` | `results/clusters` | Output directory for cluster FASTA files |

### Outputs

```
results/
├── search/
│   ├── <prefix>.<family>.domains.fasta   # Domain-extracted sequences
│   ├── <prefix>.<family>.domains.csv     # Domain hit table
│   └── <prefix>.<family>.genes.list      # Gene IDs with hits
└── clusters/
    ├── <prefix>.<family>_cluster.tsv     # Gene → cluster mapping
    └── <prefix>.<family>.<HG>.fasta      # Per-HG sequence FASTA
```

---

## Filtering homology groups

After Step 1, choose which HGs to carry forward based on size, species
composition, and presence of a species of interest (SOI).

```bash
python workflow/select_hgs.py \
  --cluster_dir results/clusters \
  --out ids.txt \
  --soi Mmus \
  --min_seqs 20 \
  --min_sps 5
```

| Flag | Default | Description |
|---|---|---|
| `--cluster_dir` / `-d` | `results/clusters` | Directory of HG FASTA files |
| `--out` | stdout | Output file, one HG id per line |
| `--soi` | *(none)* | Keep only HGs containing this species prefix |
| `--min_seqs` | `1` | Minimum number of sequences |
| `--min_sps` | `1` | Minimum number of distinct species |

The resulting `ids.txt` is the primary input to `step2.nf`.

To inspect sequence counts and lengths before filtering:

```bash
python workflow/get_seqstat.py results/clusters/*.fasta \
  | grep -wf ids.txt \
  | sort -k2 -n
```

Output columns: `HG_id  n_sequences  median_length`

---

## Step 2 — Phylogenies (`step2.nf`)

Runs alignment (MAFFT), tree inference (FastTree or IQ-TREE 2), ortholog
prediction (POSSVM), and optionally gene-tree / species-tree reconciliation
(GeneRax).

### Resource prediction (optional)

Per-HG resource requests can be predicted from sequence statistics using a
pre-fitted quantile regression model, avoiding both under- and over-allocation.

**1. Gather sequence statistics**

```bash
python workflow/get_seqstat.py results/clusters/*.fasta > seq_stat.tab
```

**2. Train models** *(skip if using the bundled `workflow/models/`)*

Requires a Nextflow trace file from a previous run:

```bash
Rscript train.R \
  --trace    reports/trace.step2.txt \
  --seq_stats seq_stat.tab \
  --outfile  workflow/models/models.json \
  --plotfile workflow/models/models.pdf \
  --tau      0.95
```

`--tau` sets the quantile: 0.95 means resources will cover 95 % of jobs
without running out. Higher values are more conservative but waste more
allocation.

**3. Predict resources for `ids.txt`**

```bash
python workflow/predict_resources.py \
  --ids_fn       ids.txt \
  --cluster_dir  results/clusters \
  --models_json  workflow/models/models.json \
  --defaults_json workflow/models/defaults.json \
  --outfile      resources.tsv \
  --max_mem      100000 \
  --max_time     2880 \
  --increase     0.5
```

| Flag | Description |
|---|---|
| `--max_mem` | Hard cap on predicted memory (MB) |
| `--max_time` | Hard cap on predicted time (minutes) |
| `--increase` | Safety margin added on top of predictions (0.5 = +50 %) |

The `resources.tsv` file is automatically picked up by `step2.nf` if it
exists at `${projectDir}/resources.tsv`.

**`workflow/models/defaults.json` structure**

```json
{
  "ALN": {
    "mem": 500, "time": 30,
    "large_nseq_threshold": 1000,
    "large_mem": 50000, "large_time": 360
  },
  "PHY": { "mem": 500, "time": 30, ... },
  "PVM": { "mem": 500, "time": 5 },
  "GR_watcher": { "mem": 2048, "time": 60, ... }
}
```

Jobs with `nseq ≥ large_nseq_threshold` receive the `large_mem`/`large_time`
values regardless of model predictions.

### Running the pipeline

**Interactive (local)**

```bash
module load Java OpenMPI
mamba activate phylo

WORKDIR=/path/to/work/step2
PROFILE=local,precise   # or: local,fast

nextflow run step2.nf \
  -profile $PROFILE \
  -w $WORKDIR \
  -resume \
  --run_generax \
  -with-report  reports/report.step2.html \
  -with-trace   reports/trace.step2.txt \
  -with-timeline reports/timeline.step2.html
```

**SLURM batch**

```bash
WORKDIR=/path/to/work/step2
PROFILE=slurm,precise

sbatch -J step2 -o reports/slurm.step2.out submit_nf.sh step2.nf \
  -profile $PROFILE \
  -resume \
  --run_generax \
  -w $WORKDIR \
  --report   reports/report.step2.html \
  --trace    reports/trace.step2.txt \
  --timeline reports/timeline.step2.html
```

### Key parameters

| Parameter | Default | Description |
|---|---|---|
| `--ids` | `ids.txt` | File listing HG ids to process |
| `--resources_tsv` | `resources.tsv` | Per-HG resource predictions (optional) |
| `--run_generax` | `false` | Enable GeneRax reconciliation + second POSSVM |
| `--OUTDIR` | `results/` | Root output directory |
| `--REFSPECIES` | `Mmus` | Reference species prefix for POSSVM |
| `--REFNAMES` | `data/Mmus_gene_names.csv` | Gene name table for reference species |
| `--SPECIES_TREE` | `data/species_tree.newick` | Species tree for GeneRax |
| `--MAFFT_OPT` | profile-dependent | MAFFT alignment options |
| `--TREE_METHOD` | profile-dependent | `fasttree` or `iqtree2` |
| `--IQTREE2_MODEL` | `TEST` | IQ-TREE 2 substitution model |
| `--SUBS_MODEL` | `LG` | Substitution model for GeneRax |
| `--MAX_SPR` | profile-dependent | GeneRax max SPR radius |
| `--NCPU_GENERAX` | profile-dependent | CPUs per GeneRax job |
| `--tag_prefix` | *(none)* | Prefix added to Nextflow job tags |

### Workflow logic

```
ids.txt
  └─ ALN (MAFFT + gap-only removal)
       └─ PHY (FastTree / IQ-TREE 2)
            ├─ [run_generax=false] PVM (POSSVM on gene tree)
            └─ [run_generax=true]
                 ├─ PVM_PREV  (POSSVM on original gene tree → results/possvm_prev/)
                 └─ GR_watcher (GeneRax reconciliation)
                      └─ PVM  (POSSVM on reconciled tree → results/possvm/)
```

Existing outputs are automatically reused via symlinks (caching without
re-computation even across `-resume` boundaries).

### Outputs

```
results/
├── align/           # *.aln.fasta — trimmed MSAs
├── gene_trees/      # *.treefile  — raw gene trees
├── possvm/          # *.ortholog_groups.csv, *.newick, *.pairs_orthologs.csv
├── possvm_prev/     # Same but from pre-GeneRax trees (when --run_generax)
└── generax/         # *.generax.tree, *.generax.log (when --run_generax)
```

---

## GeneRax-only pipeline (`generax.nf`)

Use when you already have alignments and gene trees and only want to
(re-)run GeneRax + POSSVM — for example, after changing `--MAX_SPR` or the
substitution model.

```bash
nextflow run generax.nf \
  -profile local,precise \
  -resume \
  --ids ids.txt \
  --ALIGN_DIR results/align \
  --TREE_DIR  results/gene_trees
```

| Parameter | Default | Description |
|---|---|---|
| `--ids` | `ids.txt` | HG ids to process |
| `--ALIGN_DIR` | `results/align` | Directory with `*.aln.fasta` files |
| `--TREE_DIR` | `results/gene_trees` | Directory with `*.treefile` files |
| `--SPECIES_TREE` | *(from nextflow.config)* | Species tree for GeneRax |
| `--REFSPECIES` | *(from nextflow.config)* | Reference species for POSSVM |
| `--REFNAMES` | *(from nextflow.config)* | Gene name table for reference species |

---

## Gathering annotations

After the pipeline, collect per-species ortholog annotations from the POSSVM
output files.

```bash
python workflow/gather_annotations.py \
  --search-dir results/search \
  --tree-dir   results/possvm results/possvm_prev \
  --id         sps_annotate \
  --outdir     results/annotations \
  --split-prefix
```

`--id` accepts either:
- A single species prefix: `--id Hsap`
- A file of prefixes, one per line: `--id sps_annotate`

`--tree-dir` accepts multiple directories in **priority order** — a gene
annotation found in the first directory is kept and later directories are
not consulted for that gene.

| Flag | Description |
|---|---|
| `--search-dir` | Search results directory |
| `--tree-dir` | POSSVM output dir(s), space-separated, highest priority first |
| `--id` | Species prefix or file of prefixes |
| `--outfile` | Output TSV (single species, no `--split-prefix`) |
| `--outdir` | Output directory (multiple species or `--split-prefix`) |
| `--prefix` | Restrict to gene families with this prefix, e.g. `tfs` |
| `--split-prefix` | Write one file per family prefix per species |

**Output format** (TSV, one gene per line):

```
Hsap_BRCA2    OG0001    homolog    1.00
Hsap_TP53     OG0002    ortholog   1.00
Hsap_GENE99   pfam.HG042:Unclassified
```

Column 2 is either an ortholog-group ID (from POSSVM) or
`<HG_id>:Unclassified` / `<family>:Unclassified` if no ortholog group was
assigned.

---

## Monitoring jobs

### Check a single SLURM job

```bash
python workflow/check_job.py <jobid> [<jobid> ...]
# or comma-separated:
python workflow/check_job.py 123456,123457,123458
```

Output columns: `JobID  State  ExitCode  Elapsed  Timelimit  ReqMem  MaxRSS`

### Check all jobs from a trace file

```bash
cat reports/trace.step2.txt \
  | grep -E 'COMPLETED|CACHED' \
  | cut -f3 \
  | grep -v native > job_ids

python workflow/check_job.v2.py -f job_ids > job_stats.tab
```

This aggregates `MaxRSS` across all sub-tasks of array jobs, giving the
true peak memory for each submission.

---

## Tree utilities

### Prepare the species tree

Before running GeneRax, ensure the species tree:
1. Contains exactly the species prefixes present in `input.fasta`
2. Is strictly binary (no polytomies)

```bash
python workflow/check_tree.py \
  data/species_tree.full.newick \
  species_list \
  data/species_tree.newick
```

- `species_list` — one species prefix per line
- Output is written to the third argument and validated as strictly binary

### Other tree utilities

```bash
# List all tip labels
python workflow/get_tips.py tree.newick

# Prune to a subset of tips
python workflow/prune_tree.py input.newick tips.txt output.newick

# Strip branch lengths from a Newick
python workflow/strip_len.py tree.newick > tree.topology.newick

# Strip support values from a Newick
python workflow/strip_support.py tree.newick > tree.nosupport.newick
```

---

## Post-processing / analysis

These R scripts are run **after** the pipeline and are not part of the
Nextflow workflow:

| Script | Purpose |
|---|---|
| `resources.R` | Plot memory/time usage scaling with sequence count |
| `downstream_stats.R` | Explore annotation completeness across species |
| `generax_stats.R` | GeneRax log-likelihood gain vs SPR radius and runtime breakdown |
| `train.R` | Fit quantile regression models from trace + sequence stats |
| `_export_models.R` | Export fitted R models to `workflow/models/models.json` |

---

## Parameter reference

### Global defaults (`nextflow.config`)

These values apply to all pipelines unless overridden on the command line or
by a profile.

| Parameter | Value | Used by |
|---|---|---|
| `genefam_info` | `genefam.csv` | step1.nf |
| `infasta` | `data/input.fasta` | step1.nf |
| `pfam_db` | *(Pfam-A.hmm path)* | step1.nf |
| `domain_expand` | `30` | step1.nf |
| `s1_ncpu` | `4` | step1.nf |
| `s2_ncpu` | `4` | step1.nf |
| `s2_inflation` | `1.1` | step1.nf |
| `max_n` | `2000` | step1.nf |
| `search_dir` | `results/search` | step1.nf |
| `cluster_dir` | `results/clusters` | step1.nf |
| `REFSPECIES` | `Mmus` | step2.nf, generax.nf |
| `REFNAMES` | `data/Mmus_gene_names.csv` | step2.nf, generax.nf |
| `SPECIES_TREE` | `data/species_tree.newick` | step2.nf, generax.nf |
| `IQTREE2_MODEL` | `TEST` | step2.nf |
| `SUBS_MODEL` | `LG` | step2.nf, generax.nf |

---

## Profiles reference

Profiles are combined with `,` on the command line, e.g. `-profile slurm,precise`.

### Execution profiles

| Profile | Executor | Notes |
|---|---|---|
| `local` | Local process | Use in interactive HPC sessions |
| `slurm` | SLURM (`genoa64` queue, `normal` QOS) | Use in `sbatch` jobs |

### Computation flavour profiles

| Profile | MAFFT | Tree method | GeneRax `max_spr` | GeneRax CPUs |
|---|---|---|---|---|
| `fast` | default | FastTree | 2 | 8 |
| `precise` | L-INS-i (`--maxiterate 1000 --localpair`) | IQ-TREE 2 | 3 | 2 |

---

## Output structure

```
results/
├── search/
│   ├── <pfx>.<family>.domains.fasta    # Extracted domain sequences
│   ├── <pfx>.<family>.domains.csv      # HMMER domain hit table
│   └── <pfx>.<family>.genes.list       # All genes with HMM hits
├── clusters/
│   ├── <pfx>.<family>_cluster.tsv      # Gene → HG mapping
│   └── <pfx>.<family>.<HG>.fasta       # Per-HG FASTA
├── align/
│   └── <HG>.aln.fasta                  # Trimmed MSA
├── gene_trees/
│   └── <HG>.treefile                   # Gene tree (Newick)
├── possvm/
│   ├── <HG>.*.ortholog_groups.csv      # Ortholog group table
│   ├── <HG>.*.ortholog_groups.newick   # Labelled gene tree
│   └── <HG>.*.pairs_orthologs.csv      # Pairwise ortholog pairs
├── possvm_prev/                        # Same, from pre-GeneRax trees
├── generax/
│   ├── <HG>.generax.tree               # Reconciled gene tree
│   └── <HG>.generax.log               # GeneRax log
└── annotations/
    └── <species>.<prefix>.tsv          # Per-species annotation tables
```

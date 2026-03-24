# PhyloHPC — Gene Family Phylogenomics Pipeline

A Nextflow pipeline for large-scale gene family phylogenomics: from raw proteomes to ortholog annotations, gene tree–species tree reconciliation, and clade-level ancestral gene content reconstruction.

![](img/metazoan_tree.TFevol.svg)

---

## Overview

```
Multi-species proteomes
        │
        ▼
  step1.nf ── PFAM domain search (HMMER)
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
  gather_annotations.py ── Per-species ortholog annotation tables
        │
        ▼
  step4.ancestry.nf ── Clade-scoped ortholog calling (POSSVM + --ignoretips)
                        Presence/absence matrix
                        Ancestral state reconstruction (PastML, MPPA + F81)
                        Interactive visualisation (D3.js HTML)
```

---

## Requirements

**Environment:** `workflow/environment.yaml` is the original Linux/HPC export.
For local macOS work on Apple Silicon, use the Rosetta/x86_64 environment in
`workflow/environment.macos-x86_64.yaml`.

```bash
mamba env create -f workflow/environment.yaml
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
module load Java          # for Nextflow
module load OpenMPI       # for GeneRax only
```

---

## Input data

| File | Description |
|---|---|
| `data/input.fasta` | Concatenated proteomes for all species. Sequence IDs must start with a 3–6 character species prefix followed by `_` (e.g. `Mmus_ENSMUSP00000001`). |
| `data/species_tree.full.newick` | Full species tree with **named internal nodes** (required for GeneRax and ancestral reconstruction). |
| `data/species_tree.newick` | Strictly binary species tree (no polytomies) used by GeneRax. Generated from the full tree — see below. |
| `data/Mmus_gene_names.csv` | Reference species gene name mapping (`protein_id<TAB>gene_name`). |
| `genefam.csv` | Gene family definitions (HMM profiles, clustering parameters — see `Gene_families.md`). |
| `species_list` | One species prefix per line; defines the analysis set. |

### Prepare the binary species tree

Note: `--random-resolve` is for the cases where the tree contains polytomies

```bash
python workflow/check_tree.py --random-resolve data/species_tree.full.newick species_list data/species_tree.newick
```

This validates the full tree, prunes it to the species in `species_list`, and resolves any polytomies, writing a GeneRax-compatible binary tree to `data/species_tree.newick`.

### Fetch proteomes from a sequence database

Here, use the database from Xavi. 

```bash
bash workflow/prepare_fasta.sh species_list data/input.fasta
bash workflow/prepare_fasta.sh species_list data/input.fasta /path/to/proteome_db
```

Or concatenate custom proteomes into `data/input.fasta` directly.

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


# Interactive
nextflow run step1.nf \
    -profile local,fast \
    -resume \
    -w $WORKDIR \
    --genefam_info genefam.csv \
    --infasta data/input.fasta \
    -with-report reports/report.step1.html \
    -with-trace  reports/trace.step1.txt

# SLURM
sbatch --time=01:00:00 -J step1 submit_nf.sh step1.nf \
    -resume -profile slurm,fast \
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

Rscript train.R \
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

---


## Create an interactive report 

```bash
python workflow/report_step2.py --possvm_dir results/possvm --search_dir results/search/ --cluster_dir results/clusters --family_info genefam.csv --species_tree data/species_tree.full.newick --output report2.html
```

## Step 3 — Gather per-species annotations

Collect POSSVM ortholog assignments into per-species annotation tables. Uses the GeneRax-reconciled POSSVM output when available, falling back to the raw-tree output otherwise.

```bash
python workflow/gather_annotations.py \
    --search-dir results/search/ \
    --tree-dir   results/possvm/ results/possvm_prev \
    --id         sps_annotate \
    --outdir     results/annotations/ \
    --split-prefix
```

`sps_annotate` is a file listing the species prefixes to annotate. Outputs one TSV per species (or per species × gene family prefix with `--split-prefix`), with columns: gene ID, OG name, reference OG, reference gene name.

---

## Step 4 — Ancestral gene content reconstruction

Determines which orthogroups were likely present in the last common ancestor (LCA) of a named clade, using clade-scoped ortholog calling and probabilistic ancestral state reconstruction.

### What it does

1. **`EXTRACT_CLADE`** — finds the named node (e.g. `Bilateria`) in the full species tree with named internal nodes, and outputs the in-clade species list, the out-of-clade ignore list for POSSVM, and the pruned subtree.

2. **`PVM_CLADE`** — re-runs POSSVM on the **raw gene trees** (not GeneRax-reconciled, to avoid circularity) with `--ignoretips` pointing at out-of-clade species, scoping ortholog group definitions to the target clade.

3. **`BUILD_PAM`** — aggregates all per-HG POSSVM outputs into a binary species × OG presence/absence matrix (PAM).

4. **`ANCESTRAL_RECON`** — runs [PastML](https://pastml.pasteur.fr) (MPPA + F81 model) on the PAM against the pruned species tree. The F81 model allows asymmetric gain/loss rates via estimated stationary frequencies; MPPA integrates over model parameter uncertainty for more reliable probability estimates. Outputs `P(present)` at every internal node.

5. **`VISUALIZE`** — generates a self-contained interactive HTML file for exploring the results (see below).

### Run

```bash
cat ids.txt | grep -E 'RFX' > ancestry_ids.txt
# Single clade
nextflow run step4.ancestry.nf \
    -profile local,fast \
    --node_names Bilateria \
    --ids ancestry_ids.txt

# Multiple clades in one run
nextflow run step4.ancestry.nf -profile local,fast --node_names "Metazoa,Bilateria,Deuterostomia,Chordata" --ids ancestry_ids.txt
```

`--node_names` must match named internal nodes in `data/species_tree.full.newick`.

**Key parameters**

| Parameter | Default | Description |
|---|---|---|
| `--node_names` | — | Comma-separated clade node name(s) **(required)** |
| `--ids` | `ids.txt` | HG IDs to process |
| `--SPECIES_TREE` | `data/species_tree.full.newick` | Full tree with named internal nodes |
| `--REFSPECIES` | `Mmus` | Reference species for POSSVM |
| `--gene_trees_dir` | `results/gene_trees` | Directory of raw `.treefile` outputs |
| `--min_presence` | `2` | Minimum species required to retain an OG in the PAM |

**Outputs** (under `results/ancestry/{node}/`):

| File | Description |
|---|---|
| `{node}.in_species.txt` | Species within the clade |
| `{node}.pruned.tree` | Pruned species subtree used for ASR |
| `{node}.pam.tsv` | Binary presence/absence matrix (species × OGs) |
| `{node}.ancestral_states.tsv` | Per-OG summary: `P_at_root`, gain/loss support, species counts |
| `{node}.node_probs.tsv` | `P(present)` at every internal node for every OG |
| `{node}.html` | Interactive visualisation (open in any browser) |

**`ancestral_states.tsv` columns**

| Column | Description |
|---|---|
| `og` | Orthogroup identifier |
| `n_present` | Number of in-clade species carrying the OG |
| `n_total` | Total in-clade species with data |
| `P_at_root` | Marginal posterior probability of presence at the LCA |
| `support` | Classification: `present` (≥0.9), `likely_present` (≥0.5), `likely_absent` (≥0.1), `absent` |

The visualisation can also be generated independently from existing outputs:

```bash
python workflow/visualize_ancestry.py \
    --tree       results/ancestry/Bilateria/Bilateria.pruned.tree \
    --node_probs results/ancestry/Bilateria/Bilateria.node_probs.tsv \
    --states     results/ancestry/Bilateria/Bilateria.ancestral_states.tsv \
    --pam        results/ancestry/Bilateria/Bilateria.pam.tsv \
    --output     Bilateria.html \
    --node       Bilateria
```

---

## Monitoring resource usage

```bash
# Get SLURM stats for completed jobs from a trace file
cat reports/trace.step2.txt | grep -E "COMPLETED|CACHED" | cut -f3 | grep -v native > job_ids
python workflow/check_job.v2.py -f job_ids > job_stats.tab
```

`downstream_stats.R` and `resources.R` plot CPU/memory efficiency across jobs. `generax_stats.R` analyses GeneRax runtime scaling with SPR radius.

---

## Repository layout

```
step1.nf                        PFAM search + clustering
step2.nf                        Alignment, gene trees, POSSVM, GeneRax
generax.nf                      Standalone GeneRax re-run
step4.ancestry.nf               Ancestral gene content reconstruction
nextflow.config                 Profiles and default parameters
workflow/
  environment.yaml              Conda environment
  extract_clade.py              Extract in/out-of-clade species from named node
  build_pam.py                  Build presence/absence matrix from POSSVM output
  ancestral_reconstruction.py   PastML wrapper + output parsing
  visualize_ancestry.py         Self-contained D3.js HTML generator
  gather_annotations.py         Per-species ortholog annotation tables
  select_hgs.py                 Filter HGs for downstream analysis
  get_seqstat.py                Sequence statistics per HG
  predict_resources.py          Resource prediction from trained models
  check_tree.py                 Validate and prune species tree
  check_job.v2.py               Query SLURM job statistics
  models/                       Pre-trained resource prediction models
data/
  species_tree.full.newick      Full tree with named internal nodes
  species_tree.newick           Binary tree for GeneRax
  Mmus_gene_names.csv           Reference gene name mapping
```

---

## TODOs

- [ ] HG phylogenetic profiles (presence/absence across all species)
- [ ] GeneRax: assess whether sharing information across families improves reconciliation
- [ ] Resource efficiency reports
- [ ] Re-clustering: prevent DIAMOND reruns; evaluate MMseqs2
- [ ] Allow IQ-TREE2 to resume from existing checkpoint files
- [ ] Proper handling of MAFFT OOM errors (exit code 1 vs. 137)

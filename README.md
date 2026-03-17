# PhyloHPC — Gene Family Phylogenomics Pipeline

Gene family phylogenomics: raw proteomes → ortholog annotation → gene tree/species tree reconciliation → ancestral gene content reconstruction.

Full documentation and diagrams: see the **[project mural](docs/)**.

---

## Install

```bash
mamba env create -f workflow/environment.yaml
mamba activate phylo

# HPC only
module load Java OpenMPI
```

---

## Input

| File | Description |
|---|---|
| `data/input.fasta` | Concatenated proteomes. IDs must start with a species prefix + `_` (e.g. `Mmus_ENSMUSP00000001`). |
| `data/species_tree.full.newick` | Full species tree with **named internal nodes**. |
| `data/species_tree.newick` | Binary tree for GeneRax — generate with `python workflow/check_tree.py data/species_tree.full.newick species_list data/species_tree.newick`. |
| `ids.txt` | One HG ID per line — controls which homology groups are processed. |

---

## Step 1 — PFAM search + clustering

```bash
nextflow run step1.nf -profile local,fast -resume -w /tmp/work_step1
```

Select HGs for downstream analysis:

```bash
python workflow/select_hgs.py --out ids.txt --soi Mmus --min_seqs 20 --min_sps 5
```

**Outputs:** `results/search/`, `results/clusters/`

---

## Step 2 — Alignment, gene trees, orthologs

```bash
# fast (FastTree, local)
nextflow run step2.nf -profile local,fast -resume -w /tmp/work_step2

# with GeneRax reconciliation
nextflow run step2.nf -profile local,precise -resume --run_generax -w /tmp/work_step2

# Snakemake (easier for local testing)
snakemake -s step2.smk -j4
snakemake -s step2.smk -j4 --config run_generax=1
snakemake -s step2.smk -n   # dry-run
```

**Outputs:**

| Path | Description |
|---|---|
| `results/align/` | Trimmed alignments |
| `results/gene_trees/` | Raw gene trees (`.treefile`) |
| `results/possvm/` | POSSVM ortholog groups |
| `results/possvm_prev/` | POSSVM on original IQ-TREE2 trees (only when `--run_generax`) |
| `results/generax/` | GeneRax-reconciled trees (only when `--run_generax`) |
| **`results/report_step2.html`** | **Interactive gene tree browser** — open in any browser |

### Step 2 report (`report_step2.html`)

Self-contained HTML, no server needed. Features:
- Gene tree viewer with zoom/pan, collapse-to-OG, expand-all
- Species / clade colour modes; per-species heatmap across all HGs
- **GeneRax vs Original toggle** — appears for each HG when `--run_generax` was used; switches between the GeneRax-reconciled and the raw IQ-TREE2 tree without reloading
- Tooltip: species, OG membership, subtree leaf counts

---

## Step 3 — Per-species annotation tables

```bash
python workflow/gather_annotations.py \
    --search-dir results/search/ \
    --tree-dir   results/possvm/ results/possvm_prev/ \
    --id         sps_annotate \
    --outdir     results/annotations/ \
    --split-prefix
```

---

## Step 4 — Ancestral gene content reconstruction

```bash
nextflow run step4.ancestry.nf -profile local,fast --node_names Bilateria --ids ids.txt
```

**Output: `results/ancestry/Bilateria/Bilateria.html`** — interactive D3.js ancestral state viewer.

Key files under `results/ancestry/{node}/`:

| File | Description |
|---|---|
| `{node}.pam.tsv` | Binary species × OG presence/absence matrix |
| `{node}.ancestral_states.tsv` | `P_at_root`, gain/loss support, species counts per OG |
| `{node}.node_probs.tsv` | `P(present)` at every internal node for every OG |
| **`{node}.html`** | **Interactive visualisation** |

---

## Repository layout

```
step1.nf                 PFAM search + clustering
step2.nf                 Alignment, gene trees, POSSVM, GeneRax (Nextflow)
step2.smk                Same pipeline as a Snakefile (easier for local testing)
config/step2.yaml        Snakemake defaults
step4.ancestry.nf        Ancestral gene content reconstruction
nextflow.config          Nextflow profiles and defaults
workflow/
  report_step2.py        Step 2 interactive HTML report generator
  visualize_ancestry.py  Step 4 interactive HTML report generator
  gather_annotations.py  Per-species ortholog annotation tables
  environment.yaml       Conda environment
  models/                Pre-trained resource prediction models
data/
  species_tree.full.newick
  species_tree.newick
  Mmus_gene_names.csv
```

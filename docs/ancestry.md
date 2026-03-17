# Ancestral Gene Content Reconstruction — `step4.ancestry.nf`

This document explains what the ancestral reconstruction pipeline does, how to run it, and how to interpret its outputs and the interactive visualisation.

---

## What problem does it solve?

After running `step2.nf` you have, for each Homology Group, a set of ortholog groups (OGs) — clusters of genes inferred to descend from a single ancestral copy. What you don't yet know is: **was a given OG already present in the last common ancestor (LCA) of a particular clade, or did it arise later within that clade?**

`step4.ancestry.nf` answers this by:

1. Calling OGs scoped only to the species within the target clade (ignoring outgroups during ortholog assignment)
2. Building a presence/absence matrix across those species
3. Fitting a probabilistic gain/loss model to each OG using the clade's species tree
4. Computing the posterior probability that each OG was present at the LCA

---

## Quick start

```bash
mamba activate phylo

nextflow run step4.ancestry.nf \
    -profile local,precise \
    --node_names Bilateria \
    --ids ids.txt
```

For multiple clades in one run:

```bash
nextflow run step4.ancestry.nf \
    -profile slurm,precise \
    --node_names "Bilateria,Metazoa,Deuterostomia,Vertebrata" \
    --ids ids.txt
```

`--node_names` must match named internal nodes in `data/species_tree.full.newick`. To see which named nodes are available:

```python
from ete3 import Tree
t = Tree("data/species_tree.full.newick", format=1)
print([n.name for n in t.traverse() if not n.is_leaf() and n.name])
```

---

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--node_names` | — | Comma-separated clade node name(s) **(required)** |
| `--ids` | `ids.txt` | HG IDs to analyse |
| `--SPECIES_TREE` | `data/species_tree.full.newick` | Species tree with named internal nodes |
| `--REFSPECIES` | `Mmus` | Reference species for POSSVM ortholog naming |
| `--REFNAMES` | `data/Mmus_gene_names.csv` | Reference gene name mapping |
| `--gene_trees_dir` | `results/gene_trees` | Location of raw `.treefile` outputs from step2 |
| `--min_presence` | `2` | Minimum number of clade species that must carry an OG for it to be included in the reconstruction |
| `--OUTDIR` | `results` | Root output directory |

---

## How it works — step by step

### Step 1 — `EXTRACT_CLADE`

Reads the full species tree and finds the node named by `--node_names`. Outputs:

- **`{node}.in_species.txt`** — the species within the clade
- **`{node}.ignore_species.txt`** — all other species in the full tree (passed to POSSVM as `--ignoretips`)
- **`{node}.pruned.tree`** — the subtree rooted at that node, used for ancestral reconstruction

### Step 2 — `PVM_CLADE`

For every HG in `ids.txt`, re-runs POSSVM on the **raw gene tree** (not the GeneRax-reconciled tree) with the out-of-clade species ignored.

**Why raw trees?** GeneRax reconciliation pulls gene tree topology towards the species tree to minimise duplication/loss cost. Using reconciled trees for ASR would introduce a circularity — the reconciled tree already encodes assumptions about which ancestral nodes carry the gene. The raw IQ-TREE/FastTree topology is independent evidence.

**Why `--ignoretips`?** Ortholog group boundaries are defined by the reference species. If a group of cnidarian sequences would normally pull an OG definition outward (because the reference gene has a cnidarian outparalog), ignoring those out-of-clade species keeps the OG definitions biologically meaningful within the target clade.

### Step 3 — `BUILD_PAM`

Aggregates all per-HG POSSVM CSVs into a single binary **species × OG** presence/absence matrix. A cell is `1` if the species has at least one gene assigned to that OG, `0` otherwise. OGs present in fewer than `--min_presence` species are dropped.

### Step 4 — `ANCESTRAL_RECON`

Runs [PastML](https://pastml.pasteur.fr) (Ishikawa et al. 2019, *Mol Biol Evol*) on the PAM using the **MPPA + F81** method.

**The model:** F81 is a continuous-time Markov model on binary characters (0=absent, 1=present). It has two effective rate parameters: a gain rate and a loss rate, expressed through the stationary frequencies of the two states. MPPA (Marginal Posterior Probability Approximation) computes `P(state | all tip data)` at every internal node by integrating over model parameter uncertainty rather than plugging in point estimates — this makes it more conservative and reliable than standard ML marginal reconstruction.

**Output at each internal node:** `P(present)` — the probability that the OG was present in the ancestral lineage at that point in the tree.

### Step 5 — `VISUALIZE`

Generates a self-contained HTML file for interactive exploration (see below).

---

## Output files

All outputs land in `results/ancestry/{node}/`.

| File | Description |
|---|---|
| `{node}.in_species.txt` | Species within the clade |
| `{node}.ignore_species.txt` | Species used as POSSVM ignoretips |
| `{node}.pruned.tree` | Pruned species subtree (newick) |
| `{node}.pam.tsv` | Binary presence/absence matrix (species × OGs) |
| `{node}.ancestral_states.tsv` | Per-OG summary at the LCA |
| `{node}.node_probs.tsv` | `P(present)` at every internal node for every OG |
| `{node}.html` | Interactive D3.js visualisation |

### `ancestral_states.tsv` — column reference

| Column | Description |
|---|---|
| `og` | Orthogroup identifier (e.g. `tfs.HG001.OG1`) |
| `n_present` | Number of in-clade species that carry the OG |
| `n_total` | Total in-clade species with data |
| `P_at_root` | Marginal posterior P(present) at the LCA |
| `support` | Qualitative label (see table below) |

**Support labels:**

| Label | `P_at_root` | Interpretation |
|---|---|---|
| `present` | ≥ 0.9 | Strong evidence the OG was present at the LCA |
| `likely_present` | 0.5 – 0.9 | Moderate support for ancestral presence |
| `likely_absent` | 0.1 – 0.5 | Moderate support for ancestral absence |
| `absent` | < 0.1 | Strong evidence the OG was absent at the LCA |
| `uncertain` | NaN | Reconstruction failed or insufficient data |

### `node_probs.tsv` — column reference

| Column | Description |
|---|---|
| `og` | Orthogroup identifier |
| `node` | Internal node name from the species tree |
| `P_present` | `P(present)` at this node for this OG |

This table lets you trace the exact branch where a gain or loss most likely occurred: look for the deepest node where `P_present` exceeds 0.5, and the shallowest node where it drops below 0.5.

---

## Interactive visualisation

Open `{node}.html` in any browser. No server or internet connection required.

### Layout

```
┌──────────────────┬────────────────────────────────────────────┐
│  Sidebar         │  Tree panel                                │
│                  │                                            │
│  [Search box]    │  Root ──── InternalNode ──── Species       │
│                  │              (coloured by P(present))      │
│  ▾ tfs           │                                            │
│    ▾ HG001       │  [Colour legend: red ── white ── blue]     │
│      OG1  4/30   │                                            │
│      OG2  12/30  │                                            │
│  ▾ chr           │                                            │
│    ...           │                                            │
└──────────────────┴────────────────────────────────────────────┘
```

### Controls

| Action | Effect |
|---|---|
| Click an OG in the sidebar | Recolours the tree by `P(present)` for that OG |
| Search box | Filters the OG list in real time (matches anywhere in the OG name) |
| Click an internal node | Collapses / expands that clade |
| Scroll wheel | Zoom in / out |
| Click and drag | Pan the tree |
| Hover over any node | Tooltip: node name, `P(present)`, clade size |

### Colour scale

The diverging **red → white → blue** scale maps directly to `P(present)`:

- **Deep red** — effectively absent (P ≈ 0)
- **White** — uncertain (P ≈ 0.5)
- **Deep blue** — effectively present (P ≈ 1)

Leaf nodes are coloured by their **observed state** from the PAM (red = absent, blue = present), so you can immediately see whether the ancestral call is driven by a coherent phylogenetic signal or by a patchy species distribution.

---

## Running on specific gene families only

By default `--ids ids.txt` includes all HGs selected in step 2. To analyse only transcription factors:

```bash
grep "^tfs\." ids.txt > ids_tfs.txt

nextflow run step4.ancestry.nf \
    -profile local,precise \
    --node_names Bilateria \
    --ids ids_tfs.txt
```

---

## Re-generating the visualisation from existing outputs

If you want to tweak the HTML after the pipeline has run, call the script directly without re-running the full pipeline:

```bash
python workflow/visualize_ancestry.py \
    --tree       results/ancestry/Bilateria/Bilateria.pruned.tree \
    --node_probs results/ancestry/Bilateria/Bilateria.node_probs.tsv \
    --states     results/ancestry/Bilateria/Bilateria.ancestral_states.tsv \
    --pam        results/ancestry/Bilateria/Bilateria.pam.tsv \
    --output     Bilateria.html \
    --node       Bilateria
```

`--pam` is optional. Without it, leaf nodes are shown in a neutral grey rather than coloured by observed state.

---

## Comparing results across clades

If you ran step4 on multiple clades, you can join the `ancestral_states.tsv` files to find OGs whose ancestral status differs between nested clades — for example, OGs ancestral in Metazoa but not in Bilateria (suggesting they were lost in the bilaterian lineage), or OGs ancestral in Bilateria but absent in Vertebrata.

```python
import pandas as pd

clades = ["Metazoa", "Bilateria", "Deuterostomia", "Vertebrata"]
dfs = []
for c in clades:
    df = pd.read_csv(f"results/ancestry/{c}/{c}.ancestral_states.tsv", sep="\t")
    df["clade"] = c
    dfs.append(df[["og", "clade", "P_at_root", "support"]])

combined = pd.concat(dfs).pivot(index="og", columns="clade", values="support")
print(combined)
```

---

## Reference

Ishikawa SA, Zhukova A, Iwasaki W, Gascuel O (2019).
A fast likelihood method to reconstruct and visualize ancestral scenarios.
*Molecular Biology and Evolution* 36(9): 2069–2085.
doi: [10.1093/molbev/msz131](https://doi.org/10.1093/molbev/msz131)

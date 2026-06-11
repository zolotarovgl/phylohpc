# Pipeline stream instructions

## Scope
- Improve the main pipeline itself.
- Focus on `step1.nf`, `step2.nf`, `workflow/gather_annotations.py`, `workflow/select_hgs.py`, `workflow/predict_resources.py`, and supporting scripts directly used by those stages.
- `step1.smk` is a simple Snakemake port of step 1 for exploratory runs; it should closely mirror the working sibling Snakemake implementation and use the same default inputs and output paths as `step1.nf`, while still preserving the `step1.nf` output naming scheme for search outputs and per-HG FASTAs.
- Preserve the end-to-end contract from proteomes to HGs, trees, POSSVM outputs, annotations, and report generation.
- Prefer improvements that make local and SLURM runs more robust, outputs more consistent, and downstream interpretation easier.

## Stream boundaries
- Start with:
  - `step1.nf`
  - `step1.smk`
  - `step2.nf`
  - `workflow/gather_annotations.py`
  - `workflow/select_hgs.py`
  - `workflow/predict_resources.py`
  - `workflow/report_step2.py` when the change affects step-2 report generation or its inputs
  - related tests under `tests/`
- Do not broaden into `report` or `hogs` unless:
  - the user explicitly asks
  - the bug clearly crosses the stream boundary
  - a dependency in another stream is the direct cause of the issue

## Pipeline layout
- Inputs:
  - `data/input.fasta` — concatenated proteomes
  - `species_list` — species prefixes included in the run
  - `genefam.csv` — family definitions for step 1 and step 2
  - `data/species_tree.full.newick` — named internal nodes; canonical full species tree
  - `data/species_tree.newick` — binary/pruned tree used in some downstream contexts
  - `data/Mmus_gene_names.csv` — reference gene-name mapping for POSSVM naming
- `step1.nf`:
  - `SEARCH` runs PFAM/HMMER domain search through `phylogeny/main.py hmmsearch`
  - `CLUSTER` runs MCL clustering through `phylogeny/main.py cluster`
  - main outputs: `results/search/` and `results/clusters/`
- HG selection:
  - `workflow/select_hgs.py` selects the HGs that proceed downstream
  - `ids.txt` is the main step-2 input list
- `step2.nf`:
  - `ALN` creates trimmed MSAs in `results/align/`
  - `PHY` creates per-HG gene trees in `results/gene_trees/`
    - the per-HG phylogeny log should be retained as `results/gene_trees/{id}.log`
  - optional `GR_watcher` runs GeneRax and writes reconciled trees to `results/generax/`
  - `PVM` and `PVM_PREV` run POSSVM and write ortholog outputs to `results/possvm/` and `results/possvm_prev/`
  - `REPORT` builds `results/report_step2.html`
- `workflow/gather_annotations.py`:
  - combines search results, cluster membership, and POSSVM outputs into per-species annotation tables
  - must preserve the source dimension when both `results/possvm/` and `results/possvm_prev/` exist, or deliberately choose one source explicitly

## Desired functionality
- Keep the old-manual workflow recognizable even when implementation details change.
- Preserve the main contract from proteomes to clustered HGs, trimmed alignments, gene trees, POSSVM orthology calls, annotation tables, and the step-2 report.
- `gather_annotations.py` is not optional glue; it is part of the expected pipeline behavior and should remain compatible with current outputs.
- `resources.tsv` is an optimization aid, not a biological output.

## Source-awareness rule
- Treat each HG as potentially having two parallel annotation sources:
  - GeneRax-backed POSSVM outputs in `results/possvm/`
  - non-GeneRax or IQ-TREE2-backed POSSVM outputs in `results/possvm_prev/`
- Do not assume orthogroup names, gene-to-OG assignments, reference orthologs, OG support, or reference support are identical between those sources for the same HG.
- Any utility that displays, stores, filters, joins, exports, caches, or navigates by orthogroup names or reference-ortholog metadata must keep the source dimension intact unless the caller explicitly requests a merged summary.
- The step-1 HMMER wrapper should preserve the family-wide unmerged parsed hit table as `results/search/{pref}.{family}.domains_ummerged.csv`; this is a clearer retained copy of the combined `*.domtable.csv.tmp` content and must not replace or alter the existing downstream files.
- The `phylogeny` CLI wrapper should default its log to `<outprefix>.log` when no explicit `--logfile` is provided; pipeline callers should pass and publish that log explicitly rather than discarding IQ-TREE/FastTree stdout-stderr to `/dev/null`.

## Key files
- `step1.nf`
- `step1.smk`
- `step2.nf`
- `workflow/gather_annotations.py`
- `workflow/select_hgs.py`
- `workflow/predict_resources.py`
- `workflow/report_step2.py`
- `ids.txt`
- `genefam.csv`
- `data/species_tree.full.newick`

## Regeneration
- After editing `step1.smk`, dry-run it against a config or run it into the non-conflicting exploratory directories before treating the workflow as valid:
```bash
snakemake -s step1.smk --cores 4
```
- The default `config/step1.yaml` should mirror `step1.nf`: `data/input.fasta`, `genefam.csv`, `results/search`, and `results/clusters`. Use `--config search_dir=... cluster_dir=...` only when you explicitly want a separate exploratory output root.
- After editing `workflow/report_step2.py`, regenerate the step-2 report with the same CLI shape used by the pipeline:
```bash
python workflow/report_step2.py \
    --possvm_dir results/possvm \
    --possvm_prev_dir results/possvm_prev \
    --search_dir results/search \
    --cluster_dir results/clusters \
    --family_info data/gene_families_searchinfo.csv \
    --species_tree data/species_tree.full.newick \
    --output results/report_step2.html
```
- If `results/possvm_prev/` is absent for the run you are testing, omit `--possvm_prev_dir`.

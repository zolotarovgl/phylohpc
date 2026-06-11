# hOGs stream instructions

## Scope
- Test, stabilize, and improve step 4 hierarchical orthogroup inference and reporting.
- Prioritize correctness of clade extraction, clade-level POSSVM outputs, cross-level OG linking, output layout, and `hog_hierarchy.html`.

## Active implementation
- Treat `step4_ancestry.smk` as the active implementation for step 4.
- Ignore `step4.ancestry.nf` unless the user explicitly asks about it.
- Prefer Snakemake here because it preserves the full output layout, not only declared outputs, while the Nextflow step-4 implementation is still incomplete.
- If `README.md`, `step4.ancestry.nf`, and the active step-4 code disagree, follow `step4_ancestry.smk` and the Python scripts it calls.

## Stream boundaries
- Start with:
  - `step4_ancestry.smk`
  - `workflow/build_hog_report.py`
  - `workflow/visualize_hog_hierarchy.py`
  - `workflow/link_hog_levels.py`
  - `workflow/extract_clade.py`
  - `config_ancestry.yaml`
  - `ancestry_ids.txt`
  - related tests and `results/ancestry/` outputs
- Do not broaden into `pipeline` or `report` unless:
  - the user explicitly asks
  - the bug clearly crosses the stream boundary
  - a dependency in another stream is the direct cause of the issue

## Desired functionality
- The step-4 goal is hierarchical orthogroup inference across nested clades, not standalone per-clade OG tables.
- For each HG, step 4 should:
  - extract in-clade and out-of-clade species for every requested node
  - rerun POSSVM separately for each clade level
  - filter out-of-clade genes before cross-level comparisons
  - link broader-level OGs to narrower-level OGs by shared gene membership
  - summarize splits and retentions across consecutive levels in a single report
- When changing step-4 logic, verify that downstream report-building still works on the expected file names and directory structure.
- The preferred final artifact is `results/ancestry/hog_hierarchy.html`, a self-contained HTML report for exploring hOG relationships across levels.

## Output and pipeline notes
- Preferred step-4 output layout:
  - `results/ancestry/{node}/{node}.in_species.txt`
  - `results/ancestry/{node}/{node}.ignore_species.txt`
  - `results/ancestry/{node}/{node}.pruned.tree`
  - `results/ancestry/{node}/possvm/{hg}.ortholog_groups.csv`
  - `results/ancestry/hog_hierarchy.html`
- `ancestry_ids.txt` is the HG subset used for step 4 and hOG reporting.
- POSSVM is called with `-method lpa` (label propagation); support values do not affect OG inference.
- GeneRax trees (`{hg}.generax.tree`) can be fed into POSSVM and carry posterior probability support values on internal nodes.
- Species tree node names must match exactly; `"Filosoa"` does not exist.

## Key files
- `step4_ancestry.smk`
- `workflow/build_hog_report.py`
- `workflow/visualize_hog_hierarchy.py`
- `workflow/link_hog_levels.py`
- `workflow/extract_clade.py`
- `config_ancestry.yaml`
- `ancestry_ids.txt`
- `data/species_tree.full.newick`

## Regeneration
- Always run after editing `workflow/visualize_hog_hierarchy.py` or `workflow/build_hog_report.py`:
```bash
python workflow/build_hog_report.py \
    --ancestry_dir results/ancestry \
    --nodes "Metazoa,Bilateria,Euarchontoglires" \
    --ids ancestry_ids.txt \
    --output results/ancestry/hog_hierarchy.html
```

## hOG hierarchy HTML — architecture
- Single self-contained HTML file with embedded D3.js (v7, CDN)
- Data embedded as `%%DATA_JSON%%` and replaced at build time
- Main views:
  1. Sankey or metro diagram for hOG flow across taxonomic levels
  2. OG detail popup (`#popup-backdrop`) with species tree and presence or absence
  3. Gene tree popup (`#gtpopup-backdrop`) with full gene tree and OG membership boxes

## OG detail popup — species tree (`renderPopupTree`)
- Shows the full species tree for the OG's taxonomic level (`TREES[node.level]`)
- Uses a cladogram layout with leaves aligned to the right
- Internal-node click collapses to a triangle showing `[nPresent/nTotal]`
- Collapse state lives in `popupState.collapsed`
- `#popup-rh-slider` controls row height via `popupState.rowH`
- Presence and absence are computed bottom-up from `popupState.presentSpecies`
- MRCA uses a Dollo-parsimony interpretation and is drawn as a blue circle
- Secondary losses are colored red
- The side panel shows species counts, MRCA clade name, and loss nodes
- Do not introduce species-color rendering into this popup tree; keep present black, absent grey, loss red, MRCA blue

## Gene tree popup (`renderGeneTree`)
- Parses POSSVM newick from `meta.gene_tree` with `parseNewickJS()`
- Supports cladogram or branch-length layout
- Draws OG membership boxes behind the tree plus bar-panel columns to the right
- Level-color legend toggles visibility through `gtState.hiddenLevels`
- Column-header abbreviations can isolate a single taxonomic level
- OG hover uses the shared tooltip layer and must render above the popup

## Patterns to follow when modifying
- `workflow/visualize_hog_hierarchy.py` contains the `HTML_TEMPLATE`; JS and CSS are embedded there.
- Follow the existing tree layout and collapse patterns already used in the file instead of introducing a parallel rendering model.
- Keep the report self-contained.
